//! The search behind [`Notation::options`].
//!
//! The whole of the public surface is that one function, so this module is
//! tested through it rather than on its own.
//!
//! The shape of the search is: derive one accidental per prime beyond 3, sort
//! them into the ones every notation keeps, the ones the run takes on one at a
//! time, and the ones no notation keeps, then fix one comma per dropped
//! accidental and read the run off the [`Plan`].

use diophantine::{Matrix, cvp_exact, eye, kernel_right, lll, solve_diophantine, transpose};

use crate::Error;
use crate::notation::{Notation, accidental, fifth_chain};
use crate::primes::Weighting;
use crate::temperament::Temperament;
use crate::util::{column, combination, difference, first_column, is_zero, select, subtract};

/// The derivation of the notations one temperament offers.
///
/// Everything the search does hangs off the accidentals and what each is worth
/// to the temperament, so those are worked out once here and the rest of the
/// search is methods on them. Accidentals are passed around as indices into
/// [`Self::accidentals`] throughout.
pub(crate) struct Search<'a> {
    temperament: &'a Temperament,
    /// One accidental per prime beyond 3, in order of the primes.
    accidentals: Matrix<i64>,
    /// What each accidental is worth to the temperament.
    images: Matrix<i64>,
    /// The accidentals the temperament does not temper out. One it does temper
    /// out would raise by nothing, so it is never worth keeping.
    useful: Vec<usize>,
}

/// How the search sorts the accidentals, which is the whole of the decision
/// being made: the run follows from this and the commas.
struct Plan {
    /// Kept by every notation in the run.
    necessary: Vec<usize>,
    /// Taken on one at a time in this order, each splitting apart spellings the
    /// notation before it wrote alike.
    optional: Vec<usize>,
    /// Every accidental some notation drops, in the order their commas are
    /// fixed: the optional ones first, then the ones no notation keeps.
    dropped: Vec<usize>,
}

impl Plan {
    /// The accidentals kept by the notation that has taken on `extras` of the
    /// optional ones, in order of the primes.
    fn kept(&self, extras: usize) -> Vec<usize> {
        let mut kept: Vec<usize> = self
            .necessary
            .iter()
            .chain(&self.optional[..extras])
            .copied()
            .collect();
        kept.sort_unstable();
        kept
    }

    /// The accidentals still kept wherever the `position`th dropped accidental
    /// is dropped, and so the ones its comma may be built from: the necessary
    /// ones, and the optional ones the run has already taken on.
    ///
    /// One that no notation keeps is dropped even where every optional
    /// accidental is present, so it may use them all.
    fn available(&self, position: usize) -> Vec<usize> {
        let taken = &self.optional[..position.min(self.optional.len())];
        self.necessary.iter().chain(taken).copied().collect()
    }
}

impl<'a> Search<'a> {
    /// Derives the accidentals of the temperament's subgroup and what each one
    /// is worth to it.
    ///
    /// # Errors
    /// Returns [`Error::Unsupported`] if some prime has no accidental.
    pub(crate) fn new(temperament: &'a Temperament) -> Result<Self, Error> {
        let subgroup = temperament.subgroup();
        let accidentals: Matrix<i64> = (2..subgroup.dim())
            .map(|index| accidental(subgroup, index))
            .collect::<Result<_, Error>>()?;
        let images = temperament.map_all(&accidentals)?;
        let useful = (0..accidentals.len())
            .filter(|&index| !is_zero(&images[index]))
            .collect();
        Ok(Search {
            temperament,
            accidentals,
            images,
            useful,
        })
    }

    /// Every notation the temperament offers, smallest first.
    pub(crate) fn run(&self) -> Result<Vec<Notation>, Error> {
        let plan = self.plan()?;
        let commas = self.commas(&plan)?;
        (0..=plan.optional.len())
            .map(|extras| self.notation(&plan, &commas, extras))
            .collect()
    }

    /// Sorts the accidentals into the ones every notation keeps, the ones the
    /// run takes on one at a time, and the ones no notation keeps.
    fn plan(&self) -> Result<Plan, Error> {
        let candidates = self.candidates();
        let necessary = self.necessary(&candidates)?;

        // The rest, in the order the run takes them on, passing over any worth
        // exactly what one already kept is worth: that would only be the
        // accidental already there over again. `81/80` and `64/63` being one
        // interval is a property of 41et, not a second symbol to read.
        let mut optional: Vec<usize> = Vec::new();
        for index in candidates {
            if necessary.contains(&index) {
                continue;
            }
            let mut present = necessary.iter().chain(&optional);
            if !present.any(|&kept| same_worth(&self.images[kept], &self.images[index])) {
                optional.push(index);
            }
        }

        // Everything else no notation keeps: the ones passed over, the ones
        // tempered out, and the ones an equal temperament may not reach for.
        let never_kept = (0..self.accidentals.len())
            .filter(|index| !necessary.contains(index) && !optional.contains(index));
        let dropped = optional.iter().copied().chain(never_kept).collect();

        Ok(Plan {
            necessary,
            optional,
            dropped,
        })
    }

    /// The accidentals a notation may use, in the order the run takes them on:
    /// by prime, except that an equal temperament puts the one worth a single
    /// step first.
    ///
    /// An equal temperament notated with accidentals wants the finest of them to
    /// be worth one step, which is the rule ups and downs is built on: with no
    /// symbol for a single step, single steps can only be reached by walking the
    /// fifth chain, which is not how anyone writes such a temperament. So if no
    /// accidental is worth a single step, an equal temperament gets none at all,
    /// and is left with the fifth chain if that reaches every note and with no
    /// notation if it does not. 25et over `2.3.5` is the plain case: its chain
    /// closes after five notes and its syntonic comma is worth two steps. The
    /// answer there is a subgroup whose accidental does fit, not a coarser
    /// accidental.
    ///
    /// This is the one place rank 1 is singled out, and it has to be: above rank
    /// 1 no accidental can generate the tempered lattice by itself, so there is
    /// nothing for "worth one step" to generalise to.
    fn candidates(&self) -> Vec<usize> {
        if self.temperament.rank() != 1 {
            return self.useful.clone();
        }
        let step = self
            .useful
            .iter()
            .copied()
            .find(|&index| self.images[index][0].abs() == 1);
        let Some(step) = step else {
            return Vec::new();
        };
        std::iter::once(step)
            .chain(self.useful.iter().copied().filter(|&index| index != step))
            .collect()
    }

    /// The smallest subset of `candidates` that makes a notation possible at
    /// all, which every notation in the run therefore keeps.
    ///
    /// A notation of the same rank as the temperament keeps `rank - 2` of them,
    /// so that is the smallest this can come to; where no such subset spans, the
    /// notation is forced to be larger. Keeping every candidate spans unless the
    /// candidates have been cut down, which only happens for an equal
    /// temperament with no accidental worth a single step.
    fn necessary(&self, candidates: &[usize]) -> Result<Vec<usize>, Error> {
        for size in 0..=candidates.len() {
            // Subsets in lexicographic order, so that the accidentals offered
            // first are the ones kept where there is a choice.
            for subset in subsets(candidates, size) {
                if self.spans(&subset)? {
                    return Ok(subset);
                }
            }
        }
        Err(self.nothing_spans())
    }

    /// Whether the octave, the fifth and the accidentals at `keep` generate the
    /// whole tempered lattice. Where they do not, some interval the temperament
    /// distinguishes has no spelling.
    fn spans(&self, keep: &[usize]) -> Result<bool, Error> {
        let images = self.temperament.map_all(&self.generators(keep))?;
        let identity: Matrix<i64> = eye(self.temperament.rank());
        Ok(solve_diophantine(&transpose(&images), &identity).is_ok())
    }

    /// The notational generators of a notation keeping the accidentals at
    /// `keep`: the octave, the fifth, then those accidentals.
    fn generators(&self, keep: &[usize]) -> Matrix<i64> {
        let mut generators = fifth_chain(self.temperament.dim());
        generators.extend(select(&self.accidentals, keep));
        generators
    }

    /// Why no subset of the candidates spans.
    fn nothing_spans(&self) -> Error {
        let subgroup = self.temperament.subgroup();
        if self.temperament.rank() == 1 {
            return Error::Unsupported(format!(
                "the fifth chain of this equal temperament does not reach every note, and no accidental of {subgroup} is worth a single step of it"
            ));
        }
        Error::Unsupported(format!(
            "the octave, the fifth and the accidentals of {subgroup} do not generate this rank {} temperament",
            self.temperament.rank()
        ))
    }

    /// One comma per dropped accidental, in the order [`Plan::dropped`] gives
    /// them.
    ///
    /// Fixing these once and for all is what makes taking an accidental on only
    /// ever remove a comma, so that a smaller notation's kernel always contains
    /// a larger one's and any spelling can be simplified onto a smaller
    /// notation's.
    ///
    /// It is also what makes unimodularity automatic for every subset. Nothing
    /// is ever built from an accidental that no notation keeps, so substituting
    /// the commas back is triangular and terminates at the fifth chain. That
    /// last clause matters: letting a never kept accidental stand in for another
    /// one only stacks one substitution on top of another, and it is why 41et
    /// writes `7/4` as `vBb` - `64/63` reaching for `81/80`, which the notation
    /// does keep - rather than as a detour along the chain.
    fn commas(&self, plan: &Plan) -> Result<Matrix<i64>, Error> {
        plan.dropped
            .iter()
            .enumerate()
            .map(|(position, &index)| self.comma(&plan.available(position), index))
            .collect()
    }

    /// The notational comma of the accidental at `index`: the difference between
    /// it and the replacement a notation without it must use.
    ///
    /// The replacement is built from the octave, the fifth and the accidentals
    /// at `available`, and has to be worth what the accidental is worth, so that
    /// the two are the same pitch. An accidental the temperament tempers out is
    /// replaced by nothing at all, since the temperament already calls it a
    /// unison.
    ///
    /// Which replacement is usually still a choice, settled in two steps. A
    /// plain stack of the accidentals still available keeps the nominal and the
    /// sharps that just intonation gives the prime and changes only the number
    /// of accidentals, so where there is such a stack it is the one wanted. This
    /// is what writes `33/32` as two syntonic commas in 41et and as one septimal
    /// comma in 31et. Failing that, see [`simplest_comma`](Self::simplest_comma).
    ///
    /// The stack has to come first: using the simplest comma everywhere spells
    /// 41et's `7/4` as `vvA#` instead of `vBb` and 31et's `11/8` as `vvGb`
    /// instead of `^F`, because a simpler comma is not the same as the intended
    /// one.
    fn comma(&self, available: &[usize], index: usize) -> Result<Vec<i64>, Error> {
        // The temperament already spells this accidental as a unison, so there
        // is nothing for it to be replaced by.
        if !self.useful.contains(&index) {
            return Ok(self.accidentals[index].clone());
        }

        // The accidentals still available, lowest prime first, so that a stack
        // reaches for the lower primes where it has a choice.
        let mut order: Vec<usize> = available.iter().copied().filter(|&o| o != index).collect();
        order.sort_unstable();
        let stack = select(&self.accidentals, &order);

        let target = self.temperament.map(&self.accidentals[index])?;
        let images = self.temperament.map_all(&stack)?;
        if let Some(counts) = preferred_solution(&target, &images)? {
            return Ok(difference(&self.accidentals[index], &counts, &stack));
        }

        // Failing that, the same generators with the fifth chain added.
        let mut every = stack;
        every.extend(fifth_chain(self.temperament.dim()));
        self.simplest_comma(&every, &target, index)
    }

    /// The simplest comma replacing the accidental at `index` by a combination
    /// of `generators` worth `target`.
    ///
    /// How far along the fifth chain to walk is a choice. The replacements that
    /// work differ by the commas the notation could temper out, so they form a
    /// coset of that lattice, and the element of smallest [`Weighting`] norm is
    /// taken. Preferring a short walk instead would run the accidentals away:
    /// asked to replace a one step accidental where only a two step one is
    /// available, that answers with a stack of 25 of them and an octave off
    /// rather than take five fifths.
    fn simplest_comma(
        &self,
        generators: &Matrix<i64>,
        target: &[i64],
        index: usize,
    ) -> Result<Vec<i64>, Error> {
        let columns = transpose(&self.temperament.map_all(generators)?);
        let counts = solve_diophantine(&columns, &column(target)).map_err(|_| {
            Error::Unsupported(format!(
                "the accidental {:?} of {} cannot be replaced, though it was not kept",
                self.accidentals[index],
                self.temperament.subgroup()
            ))
        })?;
        let rough = difference(&self.accidentals[index], &first_column(&counts), generators);

        // The commas the notation could temper out, which is what any two
        // replacements differ by.
        let freedom = kernel_right(&columns)?;
        if freedom.is_empty() || freedom[0].is_empty() {
            return Ok(rough);
        }
        let lattice: Matrix<i64> = transpose(&freedom)
            .iter()
            .map(|counts| combination(counts, generators, self.temperament.dim()))
            .collect();

        let weights = self.temperament.subgroup().weights(Weighting::default());
        let reduced = lll(&lattice, 0.99, &weights)?;
        Ok(subtract(&rough, &cvp_exact(&rough, &reduced, &weights)?))
    }

    /// The notation that has taken on `extras` of the optional accidentals: it
    /// keeps those and the necessary ones, and its kernel is spanned by the
    /// commas of everything else.
    fn notation(
        &self,
        plan: &Plan,
        commas: &Matrix<i64>,
        extras: usize,
    ) -> Result<Notation, Error> {
        let kept = plan.kept(extras);
        let kernel = plan
            .dropped
            .iter()
            .zip(commas)
            .filter(|(index, _)| !kept.contains(index))
            .map(|(_, comma)| comma.clone())
            .collect();
        Notation::assemble(
            self.temperament.subgroup(),
            select(&self.accidentals, &kept),
            kernel,
        )
    }
}

/// Whether two accidentals are worth the same to the temperament, so that one of
/// them is the other over again. The two may point in opposite directions.
fn same_worth(one: &[i64], other: &[i64]) -> bool {
    one == other || one.iter().zip(other).all(|(a, b)| *a == -b)
}

/// The subsets of `items` of size `size`, in lexicographic order.
fn subsets(items: &[usize], size: usize) -> Vec<Vec<usize>> {
    if size == 0 {
        return vec![Vec::new()];
    }
    let mut result = Vec::new();
    for (position, &item) in items.iter().enumerate() {
        for mut rest in subsets(&items[position + 1..], size - 1) {
            let mut subset = vec![item];
            subset.append(&mut rest);
            result.push(subset);
        }
    }
    result
}

/// The greatest common divisor of `a` and `b`, which is zero only if both are.
fn gcd(a: i64, b: i64) -> i64 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

/// Writes `target` as an integer combination of `generators`, given in order of
/// preference, most preferred first. `None` if it cannot be written as one at
/// all.
///
/// A combination is usually not unique, and the preference order picks one out.
/// The last generator is used only if `target` cannot be reached without it,
/// and then by as little as possible; then the one before it, and so on. Where
/// the amount is not pinned down either, the counts that work form an
/// arithmetic progression, and the one nearest none is taken.
fn preferred_solution(target: &[i64], generators: &Matrix<i64>) -> Result<Option<Vec<i64>>, Error> {
    let Some((last, rest)) = generators.split_last() else {
        return Ok(is_zero(target).then(Vec::new));
    };
    let rest = rest.to_vec();

    // Reaching the target without the least preferred generator at all.
    if let Some(mut counts) = preferred_solution(target, &rest)? {
        counts.push(0);
        return Ok(Some(counts));
    }

    let columns = transpose(generators);
    let Ok(solution) = solve_diophantine(&columns, &column(target)) else {
        return Ok(None);
    };

    // The counts of the last generator that work differ by whatever multiple of
    // it the others can make up for, so they run in steps of `period`.
    let kernel = kernel_right(&columns)?;
    let period = if kernel.len() == generators.len() {
        kernel[generators.len() - 1].iter().copied().fold(0, gcd)
    } else {
        0
    };

    let mut count = solution[generators.len() - 1][0];
    if period != 0 {
        count = count.rem_euclid(period);
        if 2 * count > period {
            count -= period;
        }
    }

    let reduced = subtract(
        target,
        &combination(&[count], &vec![last.clone()], target.len()),
    );
    let mut counts =
        preferred_solution(&reduced, &rest)?.expect("the remainder is reachable by construction");
    counts.push(count);
    Ok(Some(counts))
}
