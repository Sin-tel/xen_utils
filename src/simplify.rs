//! Simplifying intervals.

use diophantine::{Matrix, cvp_exact, lll};

use crate::Error;
use crate::notation::Notation;
use crate::primes::Subgroup;
use crate::util::{LLL_DELTA, subtract};

/// How far around the best so far to look, in each direction along each reduced
/// basis vector. The search steps by this much until nothing nearer is better,
/// so the radius bounds one step rather than the whole walk.
const SEARCH_RADIUS: i64 = 2;

/// Rewrites a just interval as the simplest one the temperament makes equal to
/// it.
///
/// The intervals a temperament calls the same are the cosets of its comma
/// lattice, so simplifying is a closest vector problem.
///
/// "Simplest" here means the Wilson norm [`sopfr`].
#[derive(Debug, Clone)]
pub struct Simplifier {
    subgroup: Subgroup,
    lattice: Matrix<i64>,
    weights: Matrix<f64>,
}

impl Simplifier {
    /// Builds the simplifier of the temperament `notation` notates.
    ///
    /// # Errors
    /// Returns [`Error::InvalidDimensions`] if the comma lattice cannot be
    /// computed or reduced.
    pub fn new(notation: &Notation) -> Result<Self, Error> {
        let subgroup = notation.subgroup().clone();
        let weights = subgroup.weights();

        let commas = notation.temperament().comma_basis()?;
        let lattice = lll(&commas, LLL_DELTA, &weights)?;

        Ok(Simplifier {
            subgroup,
            lattice,
            weights,
        })
    }

    /// The reduced basis of the comma lattice being searched.
    pub fn lattice(&self) -> &Matrix<i64> {
        &self.lattice
    }

    /// The simplest just interval the temperament makes equal to `interval`.
    ///
    /// # Errors
    /// Returns [`Error::InvalidDimensions`] if `interval` does not have one
    /// entry per basis element of the subgroup.
    pub fn simplify(&self, interval: &[i64]) -> Result<Vec<i64>, Error> {
        Ok(self
            .simplify_within(interval, SEARCH_RADIUS)?
            .swap_remove(0))
    }

    /// Returns a sorted list of `count` just intervals the temperament makes equal to
    /// `interval`, simplest first.
    ///
    /// # Errors
    /// Returns [`Error::InvalidDimensions`] if `interval` does not have one
    /// entry per basis element of the subgroup.
    pub fn candidates(&self, interval: &[i64], count: usize) -> Result<Matrix<i64>, Error> {
        let mut found = self.simplify_within(interval, SEARCH_RADIUS)?;
        found.truncate(count);
        Ok(found)
    }

    /// [`simplify`](Self::simplify), stepping by a ball of the given radius
    /// rather than the one [`SEARCH_RADIUS`] fixes. This is what the radius is
    /// chosen by.
    fn simplify_within(&self, interval: &[i64], radius: i64) -> Result<Matrix<i64>, Error> {
        if interval.len() != self.subgroup.dim() {
            return Err(Error::InvalidDimensions(format!(
                "interval has {} entries, expected {}",
                interval.len(),
                self.subgroup.dim()
            )));
        }
        if self.lattice.is_empty() {
            return Ok(vec![interval.to_vec()]);
        }

        let seed = subtract(
            interval,
            &cvp_exact(interval, &self.lattice, &self.weights)?,
        );
        Ok(self.search(&seed, radius))
    }

    /// Steps from `seed` to the nearest simpler interval until there is none,
    /// looking `radius` along each reduced basis vector at a time, and returns
    /// everything it passed over, simplest first.
    ///
    /// Ties under [`sopfr`] are common, so we use the quadratic norm and
    /// the coordinates themselves to break ties.
    fn search(&self, seed: &[i64], radius: i64) -> Matrix<i64> {
        // This runs for every point of a ball of `width.pow(lattice.len())`
        // points, possibly several balls deep, so it is written to allocate once
        // per point.

        let width = (2 * radius + 1) as usize;
        let mut best = seed.to_vec();
        let mut seen = Vec::new();

        // Every round strictly lowers the rank, which cannot go on for ever, so
        // this terminates whatever the radius.
        loop {
            let mut improved = best.clone();
            let mut improved_rank = self.rank(&improved);

            for point in 0..width.pow(self.lattice.len() as u32) {
                let mut candidate = best.clone();
                let mut remaining = point;
                for row in &self.lattice {
                    let digit = (remaining % width) as i64 - radius;
                    remaining /= width;
                    if digit != 0 {
                        for (c, &r) in candidate.iter_mut().zip(row) {
                            *c -= digit * r;
                        }
                    }
                }

                let rank = self.rank(&candidate);
                if rank < improved_rank {
                    improved.clone_from(&candidate);
                    improved_rank = rank;
                }
                seen.push((rank, candidate));
            }
            if improved == best {
                break;
            }
            best = improved;
        }

        seen.sort_unstable();
        seen.dedup();
        seen.into_iter().map(|(_, interval)| interval).collect()
    }

    /// What the search minimizes: the Wilson norm, and the quadratic norm
    /// standing in for it as a tie-break.
    fn rank(&self, interval: &[i64]) -> (i64, i64) {
        let primes = self.subgroup.basis();
        let squared = interval
            .iter()
            .zip(primes)
            .map(|(&e, &p)| i64::from(p) * i64::from(p) * e * e)
            .sum();
        (sopfr(interval, primes), squared)
    }
}

/// The Wilson norm of an interval: the sum of the prime factors of its
/// numerator and denominator, with repetition. `81/80` is `3*4 + 2*4 + 5 = 25`.
pub fn sopfr(interval: &[i64], primes: &[u32]) -> i64 {
    interval
        .iter()
        .zip(primes)
        .map(|(&exponent, &prime)| i64::from(prime) * exponent.abs())
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::temperament::Temperament;

    /// A notation of the equal temperament of `divisions` over `subgroup`, and
    /// its simplifier.
    ///
    /// The largest of the run rather than the recommended one, since these tests
    /// stack a notation's first accidental and the recommendation is free to
    /// keep none. Simplifying does not depend on the notation either way.
    fn simplifier(divisions: i64, subgroup: &str) -> (Notation, Simplifier) {
        let subgroup: Subgroup = subgroup.parse().unwrap();
        let temperament = Temperament::equal(divisions, &subgroup).unwrap();
        let options = Notation::options(&temperament).unwrap();
        let notation = options.last().unwrap().clone();
        let simplifier = Simplifier::new(&notation).unwrap();
        (notation, simplifier)
    }

    /// The just interval of `ups` of the first accidental stacked on `C`.
    fn stack(notation: &Notation, ups: i64) -> Vec<i64> {
        let mut spelling = vec![0; notation.rank()];
        spelling[2] = ups;
        notation.to_just(&spelling).unwrap()
    }

    /// What `interval` simplifies to, as a ratio.
    fn simplified(notation: &Notation, simplifier: &Simplifier, interval: &[i64]) -> (u64, u64) {
        let simplified = simplifier.simplify(interval).unwrap();
        notation.subgroup().to_ratio(&simplified).unwrap()
    }

    #[test]
    fn the_wilson_norm_is_the_sum_of_prime_factors() {
        let primes = [2, 3, 5, 7, 11];
        // 81/80 is 3+3+3+3 over 2+2+2+2+5, and 45/44 is 3+3+5 over 2+2+11.
        assert_eq!(sopfr(&[-4, 4, -1, 0, 0], &primes), 25);
        assert_eq!(sopfr(&[-2, 2, 1, 0, -1], &primes), 26);
        assert_eq!(sopfr(&[1, 0, 0, 0, 0], &primes), 2);
        assert_eq!(sopfr(&[0; 5], &primes), 0);
    }

    #[test]
    fn a_stack_of_ups_simplifies_to_the_octave() {
        // The accidental of 41et is worth one step, so 41 of them are an octave,
        // though as a just interval they are 81/80 forty one times over.
        let (notation, simplifier) = simplifier(41, "2.3.5.7.11");
        let stacked = stack(&notation, 41);
        assert_eq!(stacked, vec![-164, 164, -41, 0, 0]);

        let simplified = simplifier.simplify(&stacked).unwrap();
        assert_eq!(simplified, vec![1, 0, 0, 0, 0]);
        assert_eq!(notation.note(&notation.spell(&simplified).unwrap()), "C6");
    }

    #[test]
    fn forty_one_ups_are_the_interval_table_of_41et() {
        // The stack reaches every step of 41et, and what it simplifies to is
        // what those steps are usually called.
        let (notation, simplifier) = simplifier(41, "2.3.5.7.11");
        for (ups, ratio) in [
            (0, (1, 1)),
            (8, (8, 7)),
            (9, (7, 6)),
            (11, (6, 5)),
            (13, (5, 4)),
            (17, (4, 3)),
            (19, (11, 8)),
            (24, (3, 2)),
            (28, (8, 5)),
            (30, (5, 3)),
            (33, (7, 4)),
            (37, (15, 8)),
        ] {
            assert_eq!(
                simplified(&notation, &simplifier, &stack(&notation, ups)),
                ratio,
                "{ups} ups"
            );
        }
    }

    #[test]
    fn the_l1_norm_is_what_decides() {
        // A step of 41et is 64/63 or 81/80, both of norm 25, and not 45/44,
        // which the quadratic norm alone prefers at 182 against 233.
        let (notation, simplifier) = simplifier(41, "2.3.5.7.11");
        assert_eq!(
            simplified(&notation, &simplifier, &stack(&notation, 1)),
            (64, 63)
        );

        // Under the quadratic norm the seed really is the undecimal comma, so
        // it is the walk that is doing this, not the reduction.
        let seed = simplifier.simplify_within(&stack(&notation, 1), 0).unwrap();
        assert_eq!(notation.subgroup().to_ratio(&seed[0]).unwrap(), (45, 44));
    }

    #[test]
    fn a_radius_of_two_is_enough_and_a_radius_of_one_is_not() {
        // A radius of one settles in a hollow two steps wide. Fifteen octaves
        // down, 9et reaches 1/39366, of norm 29, by stepping two at once; a
        // radius of one stops at 1/34560, of norm 30.
        {
            let (notation, simplifier) = simplifier(9, "2.3.5.7");
            let far = [-15, 0, 0, 0];
            let primes = notation.subgroup().basis();
            let ratio = |interval: &[i64]| notation.subgroup().to_ratio(interval).unwrap();
            assert_eq!(
                ratio(&simplifier.simplify_within(&far, 1).unwrap()[0]),
                (1, 34560)
            );
            assert_eq!(
                ratio(&simplifier.simplify_within(&far, 2).unwrap()[0]),
                (1, 39366)
            );
            assert_eq!(sopfr(&[-8, 3, 1, 0], primes), 30);
            assert_eq!(sopfr(&[-1, 9, 0, 0], primes), 29);
        }

        for (divisions, subgroup) in [
            (15, "2.3.5"),
            (22, "2.3.5"),
            (19, "2.3.5.7"),
            (31, "2.3.5.7"),
            (41, "2.3.5.7.11"),
            (72, "2.3.5.7.11"),
        ] {
            let (notation, simplifier) = simplifier(divisions, subgroup);
            assert!(
                notation.rank() > 2,
                "{divisions}et over {subgroup} has no accidental"
            );
            for ups in -30..=60 {
                let stacked = stack(&notation, ups);
                assert_eq!(
                    simplifier.simplify_within(&stacked, 3).unwrap()[0],
                    simplifier.simplify_within(&stacked, 2).unwrap()[0],
                    "{divisions}et over {subgroup}, {ups} ups: a wider step moved the answer"
                );
            }
        }
    }

    #[test]
    fn a_second_round_can_still_improve_on_the_first() {
        // `search`'s round loop keeps re-centring on the best point found so
        // far, rather than searching one ball around the seed and stopping -
        // which is what it could afford to do if the first ball always settled
        // the question, the way `Notation::respell`'s single ball does.
        //
        // Augmented's accidental is a comma that moves quickly along the fifth
        // chain, so six of them lands somewhere the first ball's local optimum,
        // `[-31, 24, -3]`, is not the simplest interval: a second, genuinely
        // improving round - not just the confirming pass after it - trades the
        // fifth for two more fifths and drops the accidental's prime entirely,
        // one lighter under sopfr. Stopping after one ball would silently
        // return the heavier answer instead.
        let subgroup: Subgroup = "2.3.5".parse().unwrap();
        let comma = subgroup.factorize(128, 125).unwrap();
        let temperament = Temperament::from_commas(&[comma], &subgroup).unwrap();
        let notation = Notation::from_temperament(&temperament).unwrap();
        let simplifier = Simplifier::new(&notation).unwrap();

        let stacked = stack(&notation, 6);
        let all = simplifier.candidates(&stacked, usize::MAX).unwrap();

        assert_eq!(
            subgroup.to_ratio(&all[0]).unwrap(),
            (282_429_536_481, 274_877_906_944)
        );

        // The first ball's own local optimum is still among everything the
        // walk passed over, but no longer the answer: it ranks behind the
        // point the second round found.
        let one_ball = vec![-31i64, 24, -3];
        let one_ball_rank = all
            .iter()
            .position(|candidate| candidate == &one_ball)
            .expect("the first ball's local optimum is one of the points the walk visited");
        assert!(
            one_ball_rank > 0,
            "the walk should have found something better than the first ball's optimum"
        );
    }

    #[test]
    fn simplifying_keeps_the_pitch_and_settles() {
        for (divisions, subgroup) in [(41, "2.3.5.7.11"), (31, "2.3.5.7"), (22, "2.3.5")] {
            let (notation, simplifier) = simplifier(divisions, subgroup);
            assert!(
                notation.rank() > 2,
                "{divisions}et over {subgroup} has no accidental"
            );
            let temperament = notation.temperament();
            for ups in -20..=60 {
                let stacked = stack(&notation, ups);
                let simplified = simplifier.simplify(&stacked).unwrap();

                assert_eq!(
                    temperament.map(&simplified).unwrap(),
                    temperament.map(&stacked).unwrap()
                );
                assert!(
                    sopfr(&simplified, notation.subgroup().basis())
                        <= sopfr(&stacked, notation.subgroup().basis())
                );
                assert_eq!(simplifier.simplify(&simplified).unwrap(), simplified);
            }
        }
    }

    #[test]
    fn every_notation_of_a_temperament_simplifies_alike() {
        // The comma lattice is the temperament's, so the notation only decides
        // how the answer is spelled.
        let subgroup: Subgroup = "2.3.5.7.11".parse().unwrap();
        let temperament = Temperament::equal(41, &subgroup).unwrap();
        let options = Notation::options(&temperament).unwrap();
        assert_eq!(options.len(), 3);

        let interval = subgroup.factorize(45, 32).unwrap();
        let simplest = Simplifier::new(&options[0])
            .unwrap()
            .simplify(&interval)
            .unwrap();
        for option in &options {
            let simplifier = Simplifier::new(option).unwrap();
            assert_eq!(simplifier.simplify(&interval).unwrap(), simplest);
        }
    }

    #[test]
    fn the_candidates_are_ranked_readings_of_one_pitch() {
        // Fifteen steps of 41et, where the simplest reading is not the most
        // convenient spelling: 9/7 costs three marks and means something, and
        // the runner up is spelled no better and means much less.
        let (notation, simplifier) = simplifier(41, "2.3.5.7.11");
        let subgroup = notation.subgroup();
        let target = stack(&notation, 15);

        let candidates = simplifier.candidates(&target, 3).unwrap();
        let ratios: Vec<(u64, u64)> = candidates
            .iter()
            .map(|candidate| subgroup.to_ratio(candidate).unwrap())
            .collect();
        assert_eq!(ratios, vec![(9, 7), (32, 25), (35, 27)]);

        // The first is what simplifying gives, the rest are ranked behind it,
        // and every one of them is the same pitch to the temperament.
        assert_eq!(candidates[0], simplifier.simplify(&target).unwrap());
        let pitch = notation.temperament().map(&target).unwrap();
        let mut norms = Vec::new();
        for candidate in simplifier.candidates(&target, usize::MAX).unwrap() {
            assert_eq!(notation.temperament().map(&candidate).unwrap(), pitch);
            norms.push(sopfr(&candidate, subgroup.basis()));
        }
        assert!(norms.windows(2).all(|pair| pair[0] <= pair[1]));

        // The walk sees a few hundred of them, so asking for a handful is free
        // and asking for more than there are gives what there is.
        assert_eq!(norms.len(), 625);
        assert_eq!(simplifier.candidates(&target, 900).unwrap().len(), 625);
        assert_eq!(simplifier.candidates(&target, 0).unwrap().len(), 0);
    }

    #[test]
    fn just_intonation_has_nothing_to_simplify() {
        let subgroup: Subgroup = "2.3.5.7".parse().unwrap();
        let notation = Notation::from_ji(&subgroup).unwrap();
        let simplifier = Simplifier::new(&notation).unwrap();

        assert!(simplifier.lattice().is_empty());
        let interval = subgroup.factorize(225, 224).unwrap();
        assert_eq!(simplifier.simplify(&interval).unwrap(), interval);
    }

    #[test]
    fn the_coordinates_have_to_fit_the_subgroup() {
        let (_, simplifier) = simplifier(41, "2.3.5.7.11");
        assert!(simplifier.simplify(&[0, 0]).is_err());
    }
}
