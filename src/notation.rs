//! Notation systems as linear maps.

use diophantine::{Matrix, cvp_exact, eye, kernel_left, lll, solve_diophantine, transpose};

use crate::Error;
use crate::primes::Subgroup;
use crate::search::Search;
use crate::temperament::Temperament;
use crate::util::{LLL_DELTA, column, combination, first_column};

/// The nominals, in order of the fifth chain.
/// Nominals for note at `f` fifths spells `NOMINALS[(f + 1) mod 7]`.
const NOMINALS: [char; 7] = ['F', 'C', 'G', 'D', 'A', 'E', 'B'];

/// The octave number of the note with notation coordinates `(0, 0, ..)`.
const CENTRE_OCTAVE: i64 = 5;

/// Raising and lowering symbols for each accidental.
/// These exist only so that [`Notation::note`] can print something legible;
/// real microtonal accidentals are not in unicode, so anything that has to
/// look right should render the notation coordinates itself.
///
/// Only consulted when a notation keeps more than one accidental - see
/// [`accidental_symbol`]. A prime beyond `19` in such a notation has no entry
/// here and so no symbol; that is the ceiling this stops at.
const PRIME_SYMBOLS: [(u32, char, char); 6] = [
    (5, '^', 'v'),
    (7, '>', '<'),
    (11, 't', 'd'),
    (13, '*', '%'),
    (17, '/', '\\'),
    (19, ')', '('),
];

/// The largest interval, in cents, that counts as an accidental: half an
/// apotome, half of the sharp `2187/2048 = 3^7 / 2^11`, or 56.8 cents.
///
/// Bounding accidentals here is what lets every prime be written without
/// augmented or diminished intervals.
///
/// The value is `600 * (7 * log2(3) - 11)`.
const MAX_ACCIDENTAL_CENTS: f64 = 56.842_503_028_855_52;

/// A notation system: a set of symbols, and what each one maps to in the temperament.
///
/// Notation coordinates are counts of notational generators. The first two are
/// always the octave `2/1` and the fifth `3/2`, which together give the
/// nominals and the sharps and flats. Each remaining coordinate counts one accidental,
/// which raises or lowers by a small interval that is not a sharp.
///
/// **A notation is its generators.** Where they sit in the temperament fixes
/// everything else: the map [`pitch`](Self::pitch) from a written note down to
/// what it sounds, and the kernel of that map, the
/// [enharmonics](Self::enharmonics) - the written notes the temperament calls
/// one pitch.
///
/// A just interval is spelled by asking the temperament what pitch it is and then
/// asking which of that pitch's spellings reads best, which is [`spell`](Self::spell), and the
/// alternatives it passed over are [`respell`](Self::respell).
#[derive(Debug, Clone)]
pub struct Notation {
    /// The octave, the fifth, then the accidentals, as prime interval vectors.
    generators: Matrix<i64>,
    /// What each generator maps to in the temperament, one per row.
    /// This is the map from notation coordinates to pitches.
    images: Matrix<i64>,
    /// A reduced basis of the written notes that map to unison, in notation coordinates.
    enharmonics: Matrix<i64>,
    /// The quadratic form standing in for [`spelling_cost`].
    weights: Matrix<f64>,
    temperament: Temperament,
}

impl Notation {
    /// Builds the notation of `subgroup` as just intonation, with an accidental
    /// for every prime beyond 3.
    ///
    /// Spelling is a bijection here: nothing is tempered out, so no two written
    /// notes are one pitch and there are no enharmonics.
    ///
    /// # Errors
    /// Returns [`Error::Unsupported`] if the subgroup has more primes beyond 3
    /// than there are accidental symbols, or if some prime has no accidental.
    pub fn from_ji(subgroup: &Subgroup) -> Result<Self, Error> {
        let accidentals = (2..subgroup.dim())
            .map(|index| accidental(subgroup, index))
            .collect::<Result<Matrix<i64>, Error>>()?;
        Notation::from_accidentals(&Temperament::just(subgroup)?, &accidentals)
    }

    /// Builds the recommended notation of `temperament`: the smallest of
    /// the run [`options`](Self::options) gives that
    /// [keeps every nominal](Self::keeps_nominals).
    ///
    /// # Errors
    /// Returns [`Error::Unsupported`] if the notation needs more accidentals
    /// than there are symbols, if some prime has no accidental, or if no
    /// notation can be derived at all.
    pub fn from_temperament(temperament: &Temperament) -> Result<Self, Error> {
        let options = Notation::options(temperament)?;
        for option in &options {
            if option.keeps_nominals()? {
                return Ok(option.clone());
            }
        }
        // A temperament that has to walk the fifth chain to reach a prime can
        // end up with no notation keeping its nominals.
        // Those only have a single notation, so there is nothing to choose.
        Ok(options
            .first()
            .expect("there is always at least one option")
            .clone())
    }

    /// Every notation `temperament` offers, smallest first.
    ///
    /// The first is the smallest notation there is, and the last keeps every
    /// accidental that is useful.
    ///
    /// An accidental tempered out by the temperament is always dropped, and so are
    /// accidentals that map to the same pitch as earlier accidentals.
    /// Example: in 41 equal temperament, `81/80` and `64/63`, so the second is
    /// dropped.
    /// Of the rest, the smallest set that lets the notation reach every pitch at all is
    /// *necessary* and is always kept; the remainder are *optional*, and the
    /// options enable them one at a time.
    ///
    /// # Errors
    /// Returns [`Error::Unsupported`] if some prime has no accidental, if a
    /// notation in the run needs more accidentals than there are symbols, or if
    /// no notation exists at all - an equal temperament whose fifth chain does
    /// not reach every note and which has no accidental worth a single step of
    /// it has none.
    pub fn options(temperament: &Temperament) -> Result<Vec<Self>, Error> {
        Search::new(temperament)?.run()
    }

    /// Builds the notation of `temperament` with `accidentals` as its extra
    /// generators.
    ///
    /// An accidental is an interval vector over `temperament`'s subgroup. The
    /// octave and fifth are always present; this list supplies every generator
    /// beyond those two. The vectors need not be the accidentals derived by
    /// [`from_ji`](Self::from_ji).
    ///
    /// # Errors
    /// Returns [`Error::InvalidDimensions`] if an accidental is not an interval
    /// of the temperament's subgroup, and [`Error::Unsupported`] if it cannot
    /// be named, if an accidental beyond the first has no symbol, or if the
    /// generators do not reach every pitch of the temperament.
    pub fn from_accidentals(
        temperament: &Temperament,
        accidentals: &[Vec<i64>],
    ) -> Result<Self, Error> {
        let subgroup = temperament.subgroup();
        for (index, accidental) in accidentals.iter().enumerate() {
            if accidental.len() != subgroup.dim() {
                return Err(Error::InvalidDimensions(format!(
                    "accidental {index} has {} entries, expected {} for {subgroup}",
                    accidental.len(),
                    subgroup.dim()
                )));
            }
            if accidental[2..].iter().all(|&exponent| exponent == 0) {
                return Err(Error::Unsupported(format!(
                    "accidental {index} of {subgroup} has no prime beyond 2 and 3 to name"
                )));
            }
        }
        // A single accidental never needs a name of its own - see
        // `accidental_symbol` - so only two or more accidentals require every
        // one of them to have a prime in `PRIME_SYMBOLS`.
        if accidentals.len() > 1 {
            for a in accidentals {
                let prime = accidental_prime(subgroup, a);
                if !PRIME_SYMBOLS.iter().any(|&(p, ..)| p == prime) {
                    return Err(Error::Unsupported(format!(
                        "notation over {subgroup} keeps more than one accidental, and prime {prime} has no symbol"
                    )));
                }
            }
        }

        let mut generators = fifth_chain(subgroup.dim());
        generators.extend_from_slice(accidentals);
        let images = temperament.map_all(&generators)?;
        if !spans(&images, temperament.rank()) {
            return Err(Error::Unsupported(format!(
                "the octave, the fifth and these accidentals of {subgroup} do not reach every pitch of this rank {} temperament",
                temperament.rank()
            )));
        }

        // The written notes that map to unison.
        // Reduced so the CVP search has a good basis to work with.
        let weights = spelling_cost_l2(generators.len());
        let enharmonics = kernel_left(&images)?;
        let enharmonics = lll(&enharmonics, LLL_DELTA, &weights).unwrap_or(enharmonics);

        Ok(Notation {
            generators,
            images,
            enharmonics,
            weights,
            temperament: temperament.clone(),
        })
    }

    /// The notational generators as prime interval vectors, in the order their
    /// notation coordinates count them: the octave, the fifth, then the
    /// accidentals. Every accidental is an ascending interval.
    pub fn generators(&self) -> &Matrix<i64> {
        &self.generators
    }

    /// A reduced basis of the enharmonic lattice, in notation coordinates: the
    /// written notes the temperament calls a unison.
    ///
    /// This is the freedom the notation has over the temperament, so it has
    /// `rank() - temperament.rank()` members and is empty exactly when spelling
    /// is a bijection. Every notation of an equal temperament has one, since a
    /// notation is free on its octave, its fifth and its accidentals and so can
    /// never close the circle of fifths.
    pub fn enharmonics(&self) -> &Matrix<i64> {
        &self.enharmonics
    }

    /// The temperament being notated. [`from_ji`](Self::from_ji) notates just
    /// intonation itself, which tempers nothing out.
    pub fn temperament(&self) -> &Temperament {
        &self.temperament
    }

    /// The just intonation subgroup of the temperament being notated.
    pub fn subgroup(&self) -> &Subgroup {
        self.temperament.subgroup()
    }

    /// The number of notation coordinates.
    pub fn rank(&self) -> usize {
        self.generators.len()
    }

    /// The dimension of the subgroup being notated, i.e. the length of the interval
    /// vectors this notation spells.
    pub fn dim(&self) -> usize {
        self.subgroup().dim()
    }

    /// Maps notation coordinates to a pitch in the temperament.
    ///
    /// # Errors
    /// Returns [`Error::InvalidDimensions`] if `coordinates` does not have one
    /// entry per notation coordinate.
    pub fn pitch(&self, coordinates: &[i64]) -> Result<Vec<i64>, Error> {
        self.check(coordinates)?;
        Ok(combination(
            coordinates,
            &self.images,
            self.temperament.rank(),
        ))
    }

    /// The simplest way to write a just interval.
    ///
    /// # Errors
    /// Returns [`Error::InvalidDimensions`] if `interval` does not have one
    /// entry per basis element of the subgroup, and [`Error::Unsupported`] if
    /// the temperament puts it out of the notation's reach.
    pub fn spell(&self, interval: &[i64]) -> Result<Vec<i64>, Error> {
        if interval.len() != self.dim() {
            return Err(Error::InvalidDimensions(format!(
                "interval has {} entries, expected {}",
                interval.len(),
                self.dim()
            )));
        }
        let target = self.temperament.map(interval)?;
        let solution = solve_diophantine(&transpose(&self.images), &column(&target))
            .map_err(|_| Error::Unsupported(format!("{interval:?} cannot be written")))?;
        let seed = first_column(&solution);
        Ok(self
            .respell(&seed, 1)?
            .pop()
            .expect("respell always returns at least the note it was given"))
    }

    /// Returns a list of `count` different ways to write the same pitch.
    /// Sorted by simplicity, best first.
    ///
    /// # Errors
    /// Returns [`Error::InvalidDimensions`] if `coordinates` does not have one
    /// entry per notation coordinate.
    pub fn respell(&self, coordinates: &[i64], count: usize) -> Result<Vec<Vec<i64>>, Error> {
        self.check(coordinates)?;
        if self.enharmonics.is_empty() {
            return Ok(vec![coordinates.to_vec()]);
        }

        // `solve_diophantine` hands back an arbitrary vector, so reduce first.
        let mut centre = coordinates.to_vec();
        centre[1] -= NOMINAL_CENTRE;
        let seed: Vec<i64> = match cvp_exact(&centre, &self.enharmonics, &self.weights) {
            Ok(near) => coordinates.iter().zip(near).map(|(a, b)| a - b).collect(),
            Err(_) => coordinates.to_vec(),
        };

        let mut found = Vec::new();
        let mut candidate = vec![0; self.rank()];
        let mut steps = vec![-RESPELL_WIDTH; self.enharmonics.len()];
        loop {
            for slot in 0..self.rank() {
                candidate[slot] = seed[slot]
                    + steps
                        .iter()
                        .zip(&self.enharmonics)
                        .map(|(&step, row)| step * row[slot])
                        .sum::<i64>();
            }
            found.push(candidate.clone());

            let mut place = 0;
            while place < steps.len() && steps[place] == RESPELL_WIDTH {
                steps[place] = -RESPELL_WIDTH;
                place += 1;
            }
            if place == steps.len() {
                break;
            }
            steps[place] += 1;
        }

        // Compared rather than keyed, so that ranking a candidate does not clone
        // its coordinates to break the tie with.
        // TODO: This tie-break doesn't quite work.
        found.sort_unstable_by(|a, b| {
            spelling_cost(a)
                .cmp(&spelling_cost(b))
                .then_with(|| spelling_break_ties(a).cmp(&spelling_break_ties(b)))
        });
        found.dedup();
        found.truncate(count.max(1));
        Ok(found)
    }

    /// Whether every prime is written on the nominal just intonation gives it.
    ///
    /// # Errors
    /// Returns [`Error::Unsupported`] if some prime has no accidental.
    pub fn keeps_nominals(&self) -> Result<bool, Error> {
        let nominals = NOMINALS.len() as i64;
        for index in 2..self.dim() {
            let accidental = accidental(self.subgroup(), index)?;
            let mut prime = vec![0i64; self.dim()];
            prime[index] = 1;

            let wanted = just_nominal(&accidental, index);
            if self.spell(&prime)?[1].rem_euclid(nominals) != wanted {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// Converts notation coordinates back to just interval.
    ///
    /// This gives one just interval of the many the pitch could be interpreted as.
    /// Query [`Simplifier`](crate::Simplifier) to find simpler ones.
    ///
    /// # Errors
    /// Returns [`Error::InvalidDimensions`] if `coordinates` does not have one
    /// entry per notation coordinate.
    pub fn to_just(&self, coordinates: &[i64]) -> Result<Vec<i64>, Error> {
        self.check(coordinates)?;
        Ok(combination(coordinates, &self.generators, self.dim()))
    }

    /// Checks that `coordinates` has one entry per notation coordinate.
    fn check(&self, coordinates: &[i64]) -> Result<(), Error> {
        if coordinates.len() != self.rank() {
            return Err(Error::InvalidDimensions(format!(
                "interval has {} entries, expected {}",
                coordinates.len(),
                self.rank()
            )));
        }
        Ok(())
    }

    /// Writes notation coordinates as a note in scientific pitch notation,
    /// such as `C5`, `Eb4`, `vE5`.
    ///
    /// The coordinates are read as an interval up from `C5`, so `(0, 0, ..)`
    /// is `C5` itself and `(0, 1, 0, ..)`, a fifth up, is `G5`. Accidentals
    /// come before the nominal and sharps and flats after it. The symbols are
    /// placeholders for debugging, since real microtonal accidentals are not in
    /// unicode; see [`accidental_symbol`] for how one is picked.
    ///
    /// # Panics
    /// Panics if `interval` does not have one entry per notation coordinate.
    pub fn note(&self, interval: &[i64]) -> String {
        assert_eq!(
            interval.len(),
            self.rank(),
            "interval must have one entry per notation coordinate"
        );
        let (octaves, fifths) = (interval[0], interval[1]);

        let accidentals = &self.generators()[2..];
        let mut name = String::new();
        for (&count, generator) in interval[2..].iter().zip(accidentals) {
            let (up, down) = accidental_symbol(self.subgroup(), generator, accidentals.len());
            let symbol = if count < 0 { down } else { up };
            name.extend(std::iter::repeat_n(symbol, count.unsigned_abs() as usize));
        }

        // The fifth chain runs F C G D A E B and then wraps round to F#, so
        // both the nominal and the number of sharps come from `fifths + 1`.
        name.push(NOMINALS[(fifths + 1).rem_euclid(7) as usize]);
        let sharps = (fifths + 1).div_euclid(7);
        name.push_str(&if sharps < 0 { "b" } else { "#" }.repeat(sharps.unsigned_abs() as usize));

        // An octave is seven nominals and a fifth is four, which fixes how far
        // up the staff the note sits and so which octave it lands in.
        let steps = 7 * octaves + 4 * fifths;
        name.push_str(&(CENTRE_OCTAVE + steps.div_euclid(7)).to_string());

        name
    }
}

/// The octave `2/1` and the fifth `3/2` as interval vectors over a subgroup of
/// `dim` primes. Every notation is generated by these two and its accidentals.
pub(crate) fn fifth_chain(dim: usize) -> Matrix<i64> {
    let mut octave = vec![0i64; dim];
    octave[0] = 1;
    let mut fifth = vec![0i64; dim];
    (fifth[0], fifth[1]) = (-1, 1);
    vec![octave, fifth]
}

/// The nominal just intonation gives the prime at `index`, as a position on the
/// fifth chain modulo the seven nominals.
///
/// Just intonation spells a prime as its own accidental on top of a stretch of
/// the fifth chain. The accidental `a` has exponent `s = +-1` on the prime
/// itself, so `p = s * a - s * a[1] fifths - .. octaves`, which puts the prime
/// `-s * a[1]` fifths along. The syntonic comma has `s = -1` and four threes, so
/// just intonation spells 5 four fifths up, on `E`.
fn just_nominal(accidental: &[i64], index: usize) -> i64 {
    (-accidental[index] * accidental[1]).rem_euclid(NOMINALS.len() as i64)
}

/// How far up and down the fifth chain to look for an accidental.
const MAX_FIFTH_OFFSET: i64 = 12;

/// Chooses the accidental for the prime at `index` of `subgroup`: the smallest
/// detour along the fifth chain that lands within [`MAX_ACCIDENTAL_CENTS`] of
/// the prime.
///
/// The candidates are the prime shifted by 0, 1, -1, 2, -2, .. fifths, each
/// reduced by octaves. The first one small enough wins, returned as an
/// ascending interval.
/// This gives `81/80` for 5, `64/63` for 7 and `33/32` for 11.
///
/// # Errors
/// Returns [`Error::Unsupported`] if no offset within [`MAX_FIFTH_OFFSET`]
/// gives a small enough interval.
pub(crate) fn accidental(subgroup: &Subgroup, index: usize) -> Result<Vec<i64>, Error> {
    let offsets = std::iter::once(0).chain((1..=MAX_FIFTH_OFFSET).flat_map(|k| [k, -k]));
    for offset in offsets {
        let mut interval = vec![0i64; subgroup.dim()];
        interval[index] = 1;
        interval[1] = offset;

        // Octave reduction: whichever power of two lands closest to unison.
        interval[0] = -(subgroup.to_cents(&interval) / 1200.0).round_ties_even() as i64;

        if subgroup.to_cents(&interval).abs() < MAX_ACCIDENTAL_CENTS {
            return Ok(subgroup.ascending(&interval));
        }
    }
    Err(Error::Unsupported(format!(
        "no accidental within {MAX_FIFTH_OFFSET} fifths of {} in {subgroup}",
        subgroup.basis()[index]
    )))
}

/// The prime `generator` raises or lowers: the one basis element beyond the
/// octave and the fifth where it is nonzero.
///
/// Every accidental has support `{2, 3, p}` by construction, so this is always
/// exactly one index.
///
/// # Panics
/// Panics if `generator` has no such index, which cannot happen for a
/// generator [`accidental`] produced.
fn accidental_prime(subgroup: &Subgroup, generator: &[i64]) -> u32 {
    let index = generator[2..]
        .iter()
        .position(|&e| e != 0)
        .expect("an accidental has support on its own prime beyond 2 and 3");
    subgroup.basis()[index + 2]
}

/// The raising and lowering symbols for one of a notation's accidentals, given
/// how many accidentals it keeps in total.
///
/// A notation with only one accidental is using it as ups and downs does, for
/// a generic small step, so it gets the generic `^`/`v` - unless it is the
/// quartertone `33/32`, which already has its own ASCII shorthand, `t`/`d`.
/// Once a second accidental is in play there is no such generic reading left,
/// so each one gets its own fixed symbol from [`PRIME_SYMBOLS`], keyed by
/// prime rather than by position: the same prime then always prints the same
/// way, whichever other accidentals it shares the notation with.
///
/// # Panics
/// Panics if `total > 1` and `generator`'s prime has no entry in
/// [`PRIME_SYMBOLS`]; [`Notation::from_accidentals`] checks this ahead of time.
fn accidental_symbol(subgroup: &Subgroup, generator: &[i64], total: usize) -> (char, char) {
    let prime = accidental_prime(subgroup, generator);
    if total == 1 {
        return if prime == 11 { ('t', 'd') } else { ('^', 'v') };
    }
    PRIME_SYMBOLS
        .iter()
        .find(|&&(p, ..)| p == prime)
        .map(|&(_, up, down)| (up, down))
        .expect("Notation::from_accidentals checks every accidental beyond the first has a symbol")
}

/// The middle of the seven naturals, as a fifth coordinate.
///
/// `F C G D A E B` are the fifth coordinates `-1` to `5`, so `D` at `2` is the
/// middle of them.
const NOMINAL_CENTRE: i64 = 2;

/// How far either way [`Notation::respell`] walks each enharmonic.
///
/// The answers wanted are the first few and the basis is reduced, so this only
/// has to be wide enough that nothing better lies outside it.
const RESPELL_WIDTH: i64 = 1;

/// A sharp is worth two accidental marks.
const COST_FIFTH: i64 = 2;
const COST_MARK: i64 = 7;

/// What a written note costs to read.
///
/// The fifths are counted from [`NOMINAL_CENTRE`] (D) rather than from `C`, since
/// measuring from `C` makes the flat side cheaper.
/// The octave does not appear at all.
fn spelling_cost(coordinates: &[i64]) -> i64 {
    let marks: i64 = coordinates[2..].iter().map(|c| c.abs()).sum();
    COST_MARK * marks + COST_FIFTH * (coordinates[1] - NOMINAL_CENTRE).abs()
}

/// Tie-breaker for `spelling_cost`.
fn spelling_break_ties(coordinates: &[i64]) -> i64 {
    (coordinates[1] - NOMINAL_CENTRE).abs()
}

/// Whether generators with these `images` reach every pitch of a rank `rank`
/// temperament. Where they do not, some pitch cannot be written at all.
pub(crate) fn spans(images: &Matrix<i64>, rank: usize) -> bool {
    let identity: Matrix<i64> = eye(rank);
    solve_diophantine(&transpose(images), &identity).is_ok()
}

/// A quadratic stand-in for [`spelling_cost`], for the lattice algorithms, which want
/// a form rather than a count.
fn spelling_cost_l2(rank: usize) -> Matrix<f64> {
    (0..rank)
        .map(|row| {
            (0..rank)
                .map(|col| match (row == col, row < 2) {
                    (false, _) => 0.0,
                    (true, true) => (COST_FIFTH * COST_FIFTH) as f64,
                    (true, false) => (COST_MARK * COST_MARK) as f64,
                })
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn notation(subgroup: &str) -> Notation {
        Notation::from_ji(&subgroup.parse::<Subgroup>().unwrap()).unwrap()
    }

    /// The notation of the temperament of `subgroup` that tempers out `commas`.
    fn tempered(subgroup: &str, commas: &[(u64, u64)]) -> Notation {
        let subgroup: Subgroup = subgroup.parse().unwrap();
        let commas: Vec<Vec<i64>> = commas
            .iter()
            .map(|&(num, den)| subgroup.factorize(num, den).unwrap())
            .collect();
        let temperament = Temperament::from_commas(&commas, &subgroup).unwrap();
        Notation::from_temperament(&temperament).unwrap()
    }

    /// Every notation the equal temperament of `divisions` over `subgroup`
    /// offers.
    fn et_options(divisions: i64, subgroup: &str) -> Vec<Notation> {
        let subgroup: Subgroup = subgroup.parse().unwrap();
        let temperament = Temperament::equal(divisions, &subgroup).unwrap();
        Notation::options(&temperament).unwrap()
    }

    /// Every notation the temperament of `subgroup` tempering out `commas`
    /// offers.
    fn options_of(subgroup: &str, commas: &[(u64, u64)]) -> Vec<Notation> {
        let subgroup: Subgroup = subgroup.parse().unwrap();
        let commas: Vec<Vec<i64>> = commas
            .iter()
            .map(|&(num, den)| subgroup.factorize(num, den).unwrap())
            .collect();
        let temperament = Temperament::from_commas(&commas, &subgroup).unwrap();
        Notation::options(&temperament).unwrap()
    }

    /// The ranks of a run of notations.
    fn ranks(options: &[Notation]) -> Vec<usize> {
        options.iter().map(Notation::rank).collect()
    }

    /// The accidentals of a notation, as ratios.
    fn accidental_ratios(n: &Notation) -> Vec<(u64, u64)> {
        n.generators()[2..]
            .iter()
            .map(|a| n.subgroup().to_ratio(a).unwrap())
            .collect()
    }

    /// The note name of a ratio, read as an interval up from C5.
    fn note_of(n: &Notation, num: u64, den: u64) -> String {
        let interval = n.subgroup().factorize(num, den).unwrap();
        n.note(&n.spell(&interval).unwrap())
    }

    #[test]
    fn pythagorean_mapping() {
        let n = notation("2.3");
        // 2/1 is an octave and 3/1 is an octave plus a fifth.
        assert_eq!(n.spell(&[1, 0]).unwrap(), vec![1, 0]);
        assert_eq!(n.spell(&[0, 1]).unwrap(), vec![1, 1]);
        assert_eq!(n.rank(), 2);
        assert_eq!(n.dim(), 2);
        assert!(accidental_ratios(&n).is_empty());
    }

    #[test]
    fn generators_map_to_unit_vectors() {
        for subgroup in ["2.3", "2.3.5", "2.3.7", "2.3.11", "2.3.5.7.11"] {
            let n = notation(subgroup);
            for (index, generator) in n.generators().iter().enumerate() {
                let mut unit = vec![0; n.rank()];
                unit[index] = 1;
                assert_eq!(n.spell(generator).unwrap(), unit);
            }
        }
    }

    #[test]
    fn derived_accidentals() {
        // The three the fifth-chain derivation produces.
        assert_eq!(accidental_ratios(&notation("2.3.5")), vec![(81, 80)]);
        assert_eq!(accidental_ratios(&notation("2.3.7")), vec![(64, 63)]);
        assert_eq!(accidental_ratios(&notation("2.3.11")), vec![(33, 32)]);
        assert_eq!(
            accidental_ratios(&notation("2.3.5.7.11")),
            vec![(81, 80), (64, 63), (33, 32)]
        );
        // Beyond 11 the same rule keeps giving the usual choices.
        assert_eq!(accidental_ratios(&notation("2.3.13")), vec![(1053, 1024)]);
        assert_eq!(accidental_ratios(&notation("2.3.19")), vec![(513, 512)]);
    }

    #[test]
    fn accidentals_are_small_and_ascending() {
        let n = notation("2.3.5.7.11.13");
        for a in &n.generators()[2..] {
            let cents = n.subgroup().to_cents(a);
            assert!(cents > 0.0 && cents < MAX_ACCIDENTAL_CENTS);
        }
    }

    #[test]
    fn five_limit_mapping() {
        let n = notation("2.3.5");
        // 5 is four fifths up, lowered by a syntonic comma.
        assert_eq!(n.spell(&[0, 0, 1]).unwrap(), vec![0, 4, -1]);
        assert_eq!(
            n.generators(),
            &vec![vec![1, 0, 0], vec![-1, 1, 0], vec![-4, 4, -1]]
        );
    }

    #[test]
    fn spelling_checks_dimensions() {
        let n = notation("2.3");
        assert!(n.spell(&[1, 0, 0]).is_err());
        assert!(n.spell(&[1]).is_err());
    }

    #[test]
    fn notes_along_the_fifth_chain() {
        let n = notation("2.3");
        // The nominals, as fifths up and down from C5 with no octave shift.
        assert_eq!(n.note(&[0, 0]), "C5");
        assert_eq!(n.note(&[0, 1]), "G5");
        assert_eq!(n.note(&[0, 2]), "D6");
        assert_eq!(n.note(&[0, -1]), "F4");
        assert_eq!(n.note(&[0, -2]), "Bb3");
        // Seven fifths is a sharp, and it wraps back round to F.
        assert_eq!(n.note(&[0, 6]), "F#8");
        assert_eq!(n.note(&[0, 13]), "F##12");
        assert_eq!(n.note(&[0, -8]), "Fb0");
    }

    #[test]
    fn notes_across_octaves() {
        let n = notation("2.3");
        assert_eq!(n.note(&[1, 0]), "C6");
        assert_eq!(n.note(&[-1, 0]), "C4");
        assert_eq!(n.note(&[3, 0]), "C8");
        // B4 is five fifths down three octaves, just under C5.
        assert_eq!(n.note(&[-3, 5]), "B4");
        // B5 is the same chain with one octave less taken off.
        assert_eq!(n.note(&[-2, 5]), "B5");
    }

    #[test]
    fn notes_from_ratios() {
        let n = notation("2.3");
        assert_eq!(note_of(&n, 1, 1), "C5");
        assert_eq!(note_of(&n, 3, 2), "G5");
        assert_eq!(note_of(&n, 9, 8), "D5");
        assert_eq!(note_of(&n, 4, 3), "F5");
        assert_eq!(note_of(&n, 16, 9), "Bb5");
        assert_eq!(note_of(&n, 2, 1), "C6");
        assert_eq!(note_of(&n, 1, 2), "C4");
        // The pythagorean major third, 81/64.
        assert_eq!(note_of(&n, 81, 64), "E5");
        // A sharp above C5.
        assert_eq!(note_of(&n, 2187, 2048), "C#5");
    }

    #[test]
    fn accidentals_match_fjs() {
        let n = notation("2.3.5.7.11.13.17.19");

        let fjs_accidentals = vec![
            (81, 80),
            (64, 63),
            (33, 32),
            (1053, 1024),
            (4131, 4096),
            (513, 512),
        ];
        for (i, a) in n.generators()[2..].iter().enumerate() {
            let ratio = n.subgroup().to_ratio(a).unwrap();
            assert_eq!(ratio, fjs_accidentals[i])
        }
    }

    #[test]
    fn notes_with_accidentals() {
        // The just third, seventh and eleventh are each a pythagorean interval
        // bent by one accidental.
        assert_eq!(note_of(&notation("2.3.5"), 5, 4), "vE5");
        assert_eq!(note_of(&notation("2.3.7"), 7, 4), "vBb5");
        // 33/32 alone is the quartertone, so it gets t/d rather than ^/v.
        assert_eq!(note_of(&notation("2.3.11"), 11, 8), "tF5");

        // Over the full subgroup they keep those spellings, and each symbol is
        // now its prime's own: 11 keeps t/d even alongside the others.
        let n = notation("2.3.5.7.11");
        assert_eq!(note_of(&n, 5, 4), "vE5");
        assert_eq!(note_of(&n, 7, 4), "<Bb5");
        assert_eq!(note_of(&n, 11, 8), "tF5");
        // 25/16 stacks two syntonic commas; 35/32 is a whole tone bent by both.
        assert_eq!(note_of(&n, 25, 16), "vvG#5");
        assert_eq!(note_of(&n, 35, 32), "v<D5");
    }

    #[test]
    fn tempering_out_an_accidental_drops_it() {
        // Meantone tempers out 81/80 and archytas 64/63, which are exactly the
        // accidentals for 5 and 7, so both are notated as pythagorean is.
        let pythagorean = notation("2.3");
        for (subgroup, comma) in [("2.3.5", (81, 80)), ("2.3.7", (64, 63))] {
            let n = tempered(subgroup, &[comma]);
            assert_eq!(n.rank(), 2);
            assert!(accidental_ratios(&n).is_empty());
            // The comma is a unison, so it is written as one.
            let comma = n.subgroup().factorize(comma.0, comma.1).unwrap();
            assert_eq!(n.spell(&comma).unwrap(), vec![0; n.rank()]);
            // The octave and the fifth are still what they always were.
            assert_eq!(n.generators()[0][..2], pythagorean.generators()[0][..]);
            assert_eq!(n.generators()[1][..2], pythagorean.generators()[1][..]);
        }
    }

    #[test]
    fn tempered_primes_land_on_nominals() {
        // The meantone third is a plain E and the archytas seventh a plain Bb,
        // where just intonation needs an accidental on each.
        let meantone = tempered("2.3.5", &[(81, 80)]);
        assert_eq!(meantone.spell(&[0, 0, 1]).unwrap(), vec![0, 4]);
        assert_eq!(note_of(&meantone, 5, 4), "E5");
        assert_eq!(note_of(&meantone, 81, 64), "E5");

        let archytas = tempered("2.3.7", &[(64, 63)]);
        assert_eq!(archytas.spell(&[0, 0, 1]).unwrap(), vec![4, -2]);
        assert_eq!(note_of(&archytas, 7, 4), "Bb5");
        assert_eq!(note_of(&archytas, 16, 9), "Bb5");
    }

    #[test]
    fn only_tempered_accidentals_are_dropped() {
        // Marvel tempers out 225/224, which is neither accidental, so both stay
        // and the notation spells as the just intonation one does. The two are
        // not equal, since each carries the temperament it notates.
        let marvel = tempered("2.3.5.7", &[(225, 224)]);
        assert_eq!(marvel.rank(), 4);
        assert_eq!(marvel.generators(), notation("2.3.5.7").generators());
        assert_ne!(marvel.temperament(), notation("2.3.5.7").temperament());

        // Septimal meantone tempers out 81/80 as well, dropping just that one.
        let n = tempered("2.3.5.7", &[(81, 80), (225, 224)]);
        assert_eq!(accidental_ratios(&n), vec![(64, 63)]);
        assert_eq!(note_of(&n, 5, 4), "E5");
        assert_eq!(note_of(&n, 7, 4), "vBb5");
    }

    #[test]
    fn equal_temperaments_keep_the_fifth_chain() {
        // 12et tempers out 81/80, so it is notated as meantone is: rank 2 over
        // a rank 1 temperament, which is what leaves C# and Db to differ.
        let subgroup = Subgroup::p_limit(5);
        let n = Notation::from_temperament(&Temperament::equal(12, &subgroup).unwrap()).unwrap();
        // The notations differ only in the temperament they carry.
        let meantone = tempered("2.3.5", &[(81, 80)]);
        assert_eq!(n.generators(), meantone.generators());
        assert_eq!(n.generators(), meantone.generators());
        assert_eq!(note_of(&n, 5, 4), "E5");
    }

    #[test]
    fn the_smallest_notation_drops_an_accidental_that_is_not_tempered_out() {
        // Septimal meantone has two notations. Keeping the accidental for 7
        // writes the harmonic seventh as a lowered Bb; the smallest writes the
        // same pitch with sharps and flats alone, as A#.
        let options = options_of("2.3.5.7", &[(81, 80), (225, 224)]);
        assert_eq!(ranks(&options), vec![2, 3]);

        assert!(accidental_ratios(&options[0]).is_empty());
        assert_eq!(note_of(&options[0], 7, 4), "A#5");
        assert_eq!(note_of(&options[0], 5, 4), "E5");
        assert_eq!(note_of(&options[1], 7, 4), "vBb5");
    }

    #[test]
    fn the_smallest_notation_has_the_rank_of_the_temperament() {
        // Schismatic reaches 5 along the fifth chain too, eight fifths down,
        // so its just third is written as a diminished fourth.
        let schismatic = options_of("2.3.5", &[(32805, 32768)]);
        assert_eq!(schismatic[0].rank(), 2);
        assert_eq!(note_of(&schismatic[0], 5, 4), "Fb5");

        // Marvel is rank 3, so its smallest notation keeps one accidental. Only
        // the one for 5 works; keeping the one for 7 leaves 5 unreachable.
        let marvel = options_of("2.3.5.7", &[(225, 224)]);
        assert_eq!(marvel[0].rank(), 3);
        assert_eq!(accidental_ratios(&marvel[0]), vec![(81, 80)]);
        assert_eq!(note_of(&marvel[0], 5, 4), "vE5");
        assert_eq!(note_of(&marvel[0], 7, 4), "vvA#5");

        // Where the smallest already reaches the rank of the temperament and
        // keeps its nominals, it is what is recommended too.
        for (subgroup, comma) in [("2.3.5", (81, 80)), ("2.3.7", (64, 63))] {
            assert_eq!(
                options_of(subgroup, &[comma])[0].generators(),
                tempered(subgroup, &[comma]).generators()
            );
        }
    }

    #[test]
    fn equal_temperament_offers_the_step_when_the_chain_reaches() {
        // The fifth chain of 22et reaches every note, so the fifth chain alone
        // is a notation and the step accidental is only an offer. Its major
        // third at nine fifths up is a defining property, not a mistake.
        let options = et_options(22, "2.3.5.7");
        assert_eq!(ranks(&options), vec![2, 3]);
        assert!(accidental_ratios(&options[0]).is_empty());
        assert_eq!(note_of(&options[0], 5, 4), "D#5");
        assert_eq!(note_of(&options[0], 7, 4), "Bb5");
        // Its step is the syntonic comma, which puts 5 back on E.
        assert_eq!(accidental_ratios(&options[1]), vec![(81, 80)]);
        assert_eq!(note_of(&options[1], 5, 4), "vE5");
        assert_eq!(note_of(&options[1], 7, 4), "Bb5");

        // 31et is meantone, so 64/63 is its step instead and the two swap over.
        let options = et_options(31, "2.3.5.7");
        assert_eq!(ranks(&options), vec![2, 3]);
        assert_eq!(note_of(&options[0], 5, 4), "E5");
        assert_eq!(note_of(&options[0], 7, 4), "A#5");
        assert_eq!(accidental_ratios(&options[1]), vec![(64, 63)]);
        assert_eq!(note_of(&options[1], 7, 4), "vBb5");
    }

    #[test]
    fn equal_temperament_stacks_further_accidentals_on_the_step() {
        // 41et: the fifth chain, then one accidental worth a step, then a
        // second for 11. The 11 is worth two steps, so with only the step it is
        // written as two of them. 64/63 is worth a step as well, so it is the
        // step over again and is never offered as a second accidental - but it
        // is still what the seventh is spelled with.
        let options = et_options(41, "2.3.5.7.11");
        assert_eq!(ranks(&options), vec![2, 3, 4]);
        assert_eq!(accidental_ratios(&options[1]), vec![(81, 80)]);
        assert_eq!(note_of(&options[1], 7, 4), "vBb5");
        assert_eq!(note_of(&options[1], 11, 8), "^^F5");
        assert_eq!(accidental_ratios(&options[2]), vec![(81, 80), (33, 32)]);
        // 11 keeps its own t/d now that it shares the notation with 5.
        assert_eq!(note_of(&options[2], 11, 8), "tF5");

        // 72et needs its step, since its fifth chain closes early, and has two
        // further accidentals on top: 64/63 is worth two steps and 33/32 three.
        let options = et_options(72, "2.3.5.7.11");
        assert_eq!(ranks(&options), vec![3, 4, 5]);
        assert_eq!(accidental_ratios(&options[0]), vec![(81, 80)]);
        assert_eq!(note_of(&options[0], 7, 4), "vvBb5");
        assert_eq!(note_of(&options[0], 11, 8), "^^^F5");
        assert_eq!(
            accidental_ratios(&options[2]),
            vec![(81, 80), (64, 63), (33, 32)]
        );
    }

    #[test]
    fn an_accidental_worth_what_another_is_worth_is_not_offered() {
        // 41et spells 5/4 and 7/4 with the same accidental, since it tempers out
        // the difference between 81/80 and 64/63. That is a property of 41et, not
        // an excuse for a second symbol, so there is no notation with both.
        for options in [et_options(41, "2.3.5.7"), et_options(41, "2.3.5.7.11")] {
            for n in &options {
                assert!(!accidental_ratios(n).contains(&(64, 63)));
            }
        }

        // Nor at higher rank: pele tempers out 5120/5103, so its accidentals for
        // 5 and 7 are one interval and it has a single notation.
        let options = options_of("2.3.5.7", &[(5120, 5103)]);
        assert_eq!(ranks(&options), vec![3]);
        assert_eq!(accidental_ratios(&options[0]), vec![(81, 80)]);
        assert_eq!(note_of(&options[0], 7, 4), "vBb5");

        // 31et is the same story one prime up: 33/32 is worth what 64/63 is.
        let options = et_options(31, "2.3.5.7.11");
        assert_eq!(ranks(&options), vec![2, 3]);
        assert_eq!(accidental_ratios(&options[1]), vec![(64, 63)]);
        assert_eq!(note_of(&options[1], 11, 8), "^F5");
    }

    #[test]
    fn equal_temperament_with_nothing_to_offer() {
        // A meantone equal temperament tempers out the accidental for 5, so
        // nothing is left to offer. Its fifth chain reaches every note, so the
        // fifth chain alone is all it gets - and all it needs, the enharmonics
        // of 12et being the point of it.
        for divisions in [5, 7, 12, 19] {
            let options = et_options(divisions, "2.3.5");
            assert_eq!(ranks(&options), vec![2]);
            assert_eq!(note_of(&options[0], 5, 4), "E5");
        }

        // 15et and 34et need their accidental and have nothing beyond it.
        for divisions in [15, 34] {
            let options = et_options(divisions, "2.3.5");
            assert_eq!(ranks(&options), vec![3]);
            assert_eq!(accidental_ratios(&options[0]), vec![(81, 80)]);
            assert_eq!(note_of(&options[0], 5, 4), "vE5");
        }
    }

    #[test]
    fn equal_temperament_whose_accidental_is_worth_more_than_a_step() {
        // The fifth chain of 25et closes after five notes, so it needs an
        // accidental, and its syntonic comma is worth two steps rather than one.
        // An accidental is one step or it is nothing - that is the rule ups and
        // downs is built on - so these have no notation.
        for divisions in [25, 51, 54] {
            let subgroup: Subgroup = "2.3.5".parse().unwrap();
            let temperament = Temperament::equal(divisions, &subgroup).unwrap();
            assert!(Notation::options(&temperament).is_err());
        }

        // The answer is a subgroup whose accidental does fit. 24et over 2.3.5 is
        // contorted and saturates to 12et, but over 2.3.5.11 its quartertone is
        // 33/32 and worth exactly one step.
        let options = et_options(24, "2.3.5.11");
        assert_eq!(ranks(&options), vec![3]);
        assert_eq!(accidental_ratios(&options[0]), vec![(33, 32)]);
        assert_eq!(note_of(&options[0], 5, 4), "E5");
        // The sole accidental is the quartertone itself, so it is t/d.
        assert_eq!(note_of(&options[0], 11, 8), "tF5");
    }

    #[test]
    fn equal_temperament_gets_single_step_accidental() {
        // A constructed example that maps 81/80 to two steps and 64/63 to one, and does not have a single circle of fifths.
        let subgroup: Subgroup = "2.3.5.7".parse().unwrap();
        let temperament =
            Temperament::from_mapping(&vec![vec![25, 40, 60, 69]], &subgroup).unwrap();
        let options = Notation::options(&temperament).unwrap();
        assert_eq!(ranks(&options), vec![3]);
        assert_eq!(accidental_ratios(&options[0]), vec![(64, 63)]);
    }

    #[test]
    fn no_rank_two_notation_when_the_fifth_chain_does_not_reach() {
        // Blackwood's fifth chain closes after five notes, and porcupine
        // reaches 5 only through its own comma. Neither has a rank 2 notation,
        // so neither has anything smaller than the one it started with.
        for comma in [(256, 243), (250, 243)] {
            let options = options_of("2.3.5", &[comma]);
            assert_eq!(ranks(&options), vec![3]);
            assert_eq!(note_of(&options[0], 5, 4), "vE5");
        }
    }

    #[test]
    fn rank_two_temperaments_offer_the_whole_run() {
        // Huygens keeps the meantone spelling of 5, so its accidental for 5 is
        // tempered out and only 7 and 11 are left to offer. All three notations
        // spell 5 and 7 the same way; the 11 is what each of them adds.
        let options = options_of("2.3.5.7.11", &[(81, 80), (126, 125), (99, 98)]);
        assert_eq!(ranks(&options), vec![2, 3, 4]);
        assert_eq!(accidental_ratios(&options[1]), vec![(64, 63)]);
        assert_eq!(accidental_ratios(&options[2]), vec![(64, 63), (33, 32)]);
        for n in &options {
            assert_eq!(note_of(n, 5, 4), "E5");
        }
        assert_eq!(note_of(&options[0], 7, 4), "A#5");
        assert_eq!(note_of(&options[1], 7, 4), "vBb5");
        // 11 keeps t/d once 7 is kept alongside it too.
        assert_eq!(note_of(&options[2], 11, 8), "tF5");
    }

    #[test]
    fn a_dropped_accidental_becomes_a_stack_of_the_lower_ones() {
        // 11-limit marvel keeps an accidental for 5 and may drop the other two.
        // 64/63 needs the fifth chain, but 33/32 is worth exactly what the two
        // below it are worth together, so it becomes one of each.
        let options = options_of("2.3.5.7.11", &[(225, 224), (385, 384)]);
        assert_eq!(ranks(&options), vec![3, 4, 5]);
        assert_eq!(note_of(&options[1], 11, 8), "^>F5");
        // Dropping 64/63 too substitutes its own spelling into that one.
        assert_eq!(note_of(&options[0], 11, 8), "^^^Gbb5");
    }

    #[test]
    fn a_stack_is_the_shortest_one_that_works() {
        // 72et is worth one step for 81/80, two for 64/63 and three for 33/32,
        // so a notation keeping the first two may write the third either as
        // three of the first or as one of each. Two marks beat three.
        let options = et_options(72, "2.3.5.7.11");
        assert_eq!(ranks(&options), vec![3, 4, 5]);
        assert_eq!(note_of(&options[1], 11, 8), "^>F5");
        // With only 81/80 kept there is nothing to choose and it is three.
        assert_eq!(note_of(&options[0], 11, 8), "^^^F5");
    }

    #[test]
    fn the_recommendation_is_the_smallest_that_keeps_the_nominals() {
        // 41et has three notations. The smallest spells every prime off its
        // nominal, and the largest only turns the two marks of the middle one
        // into one of another kind, so the middle one is what is wanted.
        let subgroup: Subgroup = "2.3.5.7.11".parse().unwrap();
        let temperament = Temperament::equal(41, &subgroup).unwrap();
        let options = Notation::options(&temperament).unwrap();
        assert_eq!(ranks(&options), vec![2, 3, 4]);
        assert_eq!(
            Notation::from_temperament(&temperament)
                .unwrap()
                .generators(),
            options[1].generators()
        );

        // Where in the run it lands is not fixed: septimal meantone, marvel and
        // schismatic want the second of two, huygens the third of three.
        for (subgroup, commas, wanted) in [
            ("2.3.5.7", &[(81, 80), (225, 224)][..], 1),
            // Huygens wants the second of three: its rank 3 notation already
            // spells every prime on its own nominal, so the fourth accidental
            // is a symbol bought for nothing.
            ("2.3.5.7.11", &[(81, 80), (126, 125), (99, 98)][..], 1),
            ("2.3.5.7", &[(225, 224)][..], 1),
            ("2.3.5", &[(32805, 32768)][..], 1),
        ] {
            assert_eq!(
                tempered(subgroup, commas).generators(),
                options_of(subgroup, commas)[wanted].generators()
            );
        }
    }

    #[test]
    fn a_sharp_is_the_same_nominal() {
        // Flattone writes 11/8 as F#, where just intonation has tF. Seven fifths
        // leave the letter alone, so that is the same nominal and the smaller
        // notation is kept rather than taking on an accidental for 11.
        let flattone = &[(45, 44), (81, 80)][..];
        let options = options_of("2.3.5.11", flattone);
        assert_eq!(ranks(&options), vec![2, 3]);
        assert_eq!(note_of(&options[0], 11, 8), "F#5");
        // The larger notation keeps an accidental for 11 and then does not use
        // it: one sharp is cheaper to read than one mark, and both are on F.
        assert_eq!(note_of(&options[1], 11, 8), "F#5");
        assert!(options[0].keeps_nominals().unwrap());
        assert_eq!(
            tempered("2.3.5.11", flattone).generators(),
            options[0].generators()
        );
    }

    #[test]
    fn a_prime_off_its_nominal_is_refused() {
        // The spellings the rule is there to reject, each the smallest notation
        // of its temperament: 5 on F rather than E, 7 on A rather than B, and 11
        // on G rather than F.
        for (subgroup, commas) in [
            ("2.3.5", &[(32805, 32768)][..]),
            ("2.3.5.7", &[(81, 80), (225, 224)][..]),
            ("2.3.5.7.11", &[(81, 80), (126, 125), (99, 98)][..]),
        ] {
            assert!(!options_of(subgroup, commas)[0].keeps_nominals().unwrap());
        }
        // Just intonation keeps every nominal, by definition.
        for subgroup in ["2.3", "2.3.5", "2.3.5.7.11", "2.3.5.7.11.13"] {
            assert!(notation(subgroup).keeps_nominals().unwrap());
        }
    }

    #[test]
    fn a_temperament_accepts_each_requested_11_limit_prefix() {
        let subgroup: Subgroup = "2.3.5.7.11".parse().unwrap();
        let temperament = Temperament::equal(41, &subgroup).unwrap();
        let accidentals =
            [(81, 80), (64, 63), (33, 32)].map(|(num, den)| subgroup.factorize(num, den).unwrap());

        let notations: Vec<Notation> = (1..=accidentals.len())
            .map(|count| Notation::from_accidentals(&temperament, &accidentals[..count]).unwrap())
            .collect();

        assert_eq!(
            notations.iter().map(accidental_ratios).collect::<Vec<_>>(),
            vec![
                vec![(81, 80)],
                vec![(81, 80), (64, 63)],
                vec![(81, 80), (64, 63), (33, 32)],
            ]
        );
        assert_eq!(note_of(&notations[0], 11, 8), "^^F5");
        assert_eq!(note_of(&notations[2], 11, 8), "tF5");
    }

    #[test]
    fn an_accidental_list_must_reach_every_pitch() {
        let subgroup: Subgroup = "2.3.5".parse().unwrap();
        let temperament = Temperament::equal(25, &subgroup).unwrap();
        assert!(Notation::from_accidentals(&temperament, &[]).is_err());

        let syntonic = subgroup.factorize(81, 80).unwrap();
        let notation = Notation::from_accidentals(&temperament, &[syntonic]).unwrap();
        assert_eq!(notation.rank(), 3);
        assert_eq!(note_of(&notation, 5, 4), "vE5");
    }

    #[test]
    fn too_many_accidentals() {
        // 5, 7, 11, 13, 17 and 19 all have symbols; 23 does not.
        assert!(Notation::from_ji(&Subgroup::p_limit(19)).is_ok());
        assert!(Notation::from_ji(&Subgroup::p_limit(23)).is_err());
    }

    #[test]
    fn the_circle_of_fifths_never_closes() {
        // An equal temperament must temper out the pythagorean comma, and no
        // notation can: a notation is free on its octave, its fifth and its
        // accidentals. That difference is the enharmonic lattice, and in 12et
        // it is exactly what leaves C# and Db to differ.
        let subgroup = Subgroup::p_limit(5);
        let t = Temperament::equal(12, &subgroup).unwrap();
        let n = Notation::from_accidentals(&t, &[]).unwrap();
        assert_eq!(n.enharmonics().len(), 1);
        // Twelve fifths less seven octaves, which is the pythagorean comma.
        assert_eq!(n.enharmonics()[0], vec![-7, 12]);
        let comma = n.to_just(&n.enharmonics()[0]).unwrap();
        assert_eq!(subgroup.to_ratio(&comma).unwrap(), (531441, 524288));
    }

    #[test]
    fn enharmonics_are_the_freedom_over_the_temperament() {
        // 22et: the fifth chain alone leaves one enharmonic, and the notation
        // with an accidental leaves two - its octave being twenty two steps and
        // its fifth thirteen, which is the whole of ups and downs in 22et.
        let subgroup: Subgroup = "2.3.5.7".parse().unwrap();
        let t = Temperament::equal(22, &subgroup).unwrap();
        let syntonic = subgroup.factorize(81, 80).unwrap();
        let bare = Notation::from_accidentals(&t, &[]).unwrap();
        let raised = Notation::from_accidentals(&t, &[syntonic]).unwrap();

        assert_eq!(bare.enharmonics().len(), 1);
        assert_eq!(raised.enharmonics().len(), 2);
        for e in raised.enharmonics() {
            assert_eq!(raised.pitch(e).unwrap(), vec![0]);
        }

        // Every requested notation has one enharmonic per extra coordinate over
        // the temperament, and each is a non-unison it tempers out.
        for (divisions, subgroup) in [(12, "2.3.5"), (22, "2.3.5.7"), (41, "2.3.5.7.11")] {
            let subgroup: Subgroup = subgroup.parse().unwrap();
            let t = Temperament::equal(divisions, &subgroup).unwrap();
            let accidentals = Notation::from_ji(&subgroup).unwrap().generators()[2..].to_vec();
            for count in 0..=accidentals.len() {
                let Ok(n) = Notation::from_accidentals(&t, &accidentals[..count]) else {
                    continue;
                };
                assert_eq!(n.enharmonics().len(), n.rank() - t.rank());
                for e in n.enharmonics() {
                    // Worth nothing as a pitch, but not the unison on the page.
                    assert_eq!(n.pitch(e).unwrap(), vec![0; t.rank()]);
                    assert!(e.iter().any(|&x| x != 0));
                }
            }
        }
    }

    #[test]
    fn a_bijection_has_no_enharmonics() {
        // Where the notation has the rank of the temperament its kernel is the
        // temperament's, so there is nothing left over to choose.
        let subgroup: Subgroup = "2.3.5.7".parse().unwrap();
        let commas = [vec![-4, 4, -1, 0], vec![1, 2, -3, 1]];
        let t = Temperament::from_commas(&commas, &subgroup).unwrap();
        let n = Notation::from_accidentals(&t, &[]).unwrap();
        assert_eq!(n.rank(), t.rank());
        assert!(n.enharmonics().is_empty());
    }

    #[test]
    fn spelling_cost_agrees_l2() {
        let rank = 5;
        let weights = spelling_cost_l2(5);

        // Octaves are free only for the l1 cost
        assert_eq!(spelling_cost(&vec![5, 2, 0, 0, 0]), 0);

        for i in 1..rank {
            // Relative to D5
            let mut interval = vec![-1, 2, 0, 0, 0];
            interval[i] += 1;

            let w_l1 = spelling_cost(&interval);
            let w_l2 = weights[i][i];

            assert_eq!((w_l1 * w_l1) as f64, w_l2);
        }
    }
}
