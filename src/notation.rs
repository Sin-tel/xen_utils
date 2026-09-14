//! Notation systems as linear maps.

use diophantine::{Matrix, integer_inverse, kernel_left, transpose};

use crate::Error;
use crate::primes::Subgroup;
use crate::search::Search;
use crate::temperament::Temperament;
use crate::util::combination;

/// The nominals, in order of the fifth chain. `F` is one fifth below `C`, so
/// the fifth coordinate `f` picks out `NOMINALS[(f + 1) mod 7]`.
const NOMINALS: [char; 7] = ['F', 'C', 'G', 'D', 'A', 'E', 'B'];

/// The octave number of the note with notation coordinates `(0, 0, ..)`, i.e.
/// of the `C` that the fifth chain is centred on.
const CENTRE_OCTAVE: i64 = 5;

/// Raising and lowering symbols for each accidental, in order. These exist
/// only so that [`Notation::note`] can print something legible; real
/// microtonal accidentals are not in unicode, so anything that has to look
/// right should render the notation coordinates itself.
const ACCIDENTAL_SYMBOLS: [(char, char); 4] = [('^', 'v'), ('>', '<'), ('+', '-'), ('*', '%')];

/// How far up and down the fifth chain to look for an accidental.
const MAX_FIFTH_OFFSET: i64 = 64;

/// The largest interval, in cents, that counts as an accidental: half an
/// apotome, half of the sharp `2187/2048 = 3^7 / 2^11`, or 56.8 cents.
///
/// Anything wider than this is closer to the neighbouring point of the fifth
/// chain than to this one, so it belongs there as a sharp or a flat instead.
/// Bounding accidentals here is what lets every prime be written directly,
/// with no augmented or diminished interval needed.
///
/// The value is `600 * (7 * log2(3) - 11)`.
const MAX_ACCIDENTAL_CENTS: f64 = 56.842_503_028_855_52;

/// A notation system: a linear map from the interval vectors of a just
/// intonation [`Subgroup`] to notation coordinates.
///
/// Notation coordinates are counts of notational generators. The first two are
/// always the octave `2/1` and the fifth `3/2`, which together give the
/// nominals and the sharps and flats - a sharp is seven fifths less four
/// octaves, `2187/2048`. Each remaining coordinate counts one accidental,
/// which raises or lowers by a small interval that is not a sharp.
///
/// The map is built so that each generator gets the corresponding unit vector:
/// `2/1` maps to `(1, 0, ..)`, `3/2` to `(0, 1, 0, ..)` and the first
/// accidental to `(0, 0, 1, ..)`. What is left over is the kernel, spanned by
/// the notational [commas](Self::commas): the intervals the notation spells
/// identically. [`from_ji`](Self::from_ji) has no commas at all, so spelling
/// is a bijection; [`from_temperament`](Self::from_temperament) drops the
/// accidentals a temperament makes redundant, and those become its commas.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Notation {
    mapping: Matrix<i64>,
    generators: Matrix<i64>,
    commas: Matrix<i64>,
    subgroup: Subgroup,
}

impl Notation {
    /// Builds the notation of `subgroup` as just intonation, with an accidental
    /// for every prime beyond 3 and so no commas.
    ///
    /// # Errors
    /// Returns [`Error::Unsupported`] if the subgroup has more primes beyond 3
    /// than there are accidental symbols, or if some prime has no accidental.
    pub fn from_ji(subgroup: &Subgroup) -> Result<Self, Error> {
        let accidentals = (2..subgroup.dim())
            .map(|index| accidental(subgroup, index))
            .collect::<Result<Matrix<i64>, Error>>()?;
        Notation::assemble(subgroup, accidentals, Vec::new())
    }

    /// Builds the notation of `temperament` worth recommending: the smallest of
    /// the run [`options`](Self::options) gives that
    /// [keeps every nominal](Self::keeps_nominals).
    ///
    /// Each further accidental costs another symbol to read, so the smallest
    /// notation is the one wanted - but not at the price of moving a prime onto
    /// another letter. 41et is the plain case. Its smallest notation writes
    /// `5/4` as `Fb5` and `11/8` as `D###5`, which are the right pitches spelled
    /// off their nominals; the next one up keeps `81/80` and writes `vE5` and
    /// `^^F5`, and that is the one this returns. Taking the third accidental as
    /// well would only turn `^^F5` into `>F5`.
    ///
    /// Septimal meantone is the same story one accidental down: `A#5` for the
    /// harmonic seventh is off its nominal, so the notation that writes `vBb5`
    /// is the one returned.
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
        // A temperament whose accidentals are all worth more than a step of it
        // has to walk the fifth chain to replace them, and can end up with no
        // notation that keeps its nominals - 13et over `2.3.5` is the smallest
        // such. Every one of those offers a single notation, though, so there is
        // nothing to choose between and this returns the only one there is.
        Ok(options
            .last()
            .expect("there is always at least one option")
            .clone())
    }

    /// Every notation `temperament` offers, smallest first.
    ///
    /// The first is the smallest notation there is, and the last keeps every
    /// accidental worth keeping; which of them to use is a matter of taste, so
    /// the choice is left open. Enabling an accidental splits apart spellings
    /// the one before it wrote alike, at the cost of another symbol to read.
    ///
    /// Each accidental falls into one of three classes. One the temperament
    /// tempers out would raise by nothing and is always dropped, and so is one
    /// worth exactly what an accidental already kept is worth, since it would
    /// only be that accidental over again - `81/80` and `64/63` being the same
    /// interval is a property of 41et, not a second symbol to read. Of the rest,
    /// the smallest set that makes the notation span the subgroup at all is
    /// *necessary* and is always kept; the remainder are *optional*, and the
    /// run enables them one at a time.
    ///
    /// The first notation therefore has the rank of the temperament wherever one
    /// of that rank exists, and the last keeps an accidental for every prime the
    /// temperament does not spell on the fifth chain already.
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

    /// Builds the notation with the given accidentals and commas, deriving the
    /// mapping from them.
    pub(crate) fn assemble(
        subgroup: &Subgroup,
        accidentals: Matrix<i64>,
        commas: Matrix<i64>,
    ) -> Result<Self, Error> {
        if accidentals.len() > ACCIDENTAL_SYMBOLS.len() {
            return Err(Error::Unsupported(format!(
                "notation over {subgroup} needs {} accidentals, but only {} are named",
                accidentals.len(),
                ACCIDENTAL_SYMBOLS.len()
            )));
        }

        // The generators, as interval vectors: the octave, the fifth, then the
        // accidentals, in the order their notation coordinates count them.
        let mut generators = fifth_chain(subgroup.dim());
        generators.extend(accidentals);
        let rank = generators.len();

        // The mapping is whatever sends each generator to its unit vector and
        // each comma to zero, i.e. `mapping * transpose(square) = [I | 0]`.
        // Inverting needs the commas to fill the generators out to a square
        // matrix, and only succeeds if the two together span the subgroup.
        let mut square = generators.clone();
        square.extend(commas.iter().cloned());
        let inverse = integer_inverse(&square)?;
        let mapping = transpose(&inverse).into_iter().take(rank).collect();

        Ok(Notation {
            mapping,
            generators,
            commas,
            subgroup: subgroup.clone(),
        })
    }

    /// The mapping matrix (rows = notation coordinates, columns = basis
    /// elements of [`Self::subgroup`]).
    pub fn mapping(&self) -> &Matrix<i64> {
        &self.mapping
    }

    /// The notational generators as interval vectors, in the order their
    /// notation coordinates count them: the octave, the fifth, then the
    /// accidentals. Every accidental is an ascending interval.
    pub fn generators(&self) -> &Matrix<i64> {
        &self.generators
    }

    /// A basis for the kernel: the commas this notation spells as a unison.
    /// Empty when spelling is a bijection.
    pub fn commas(&self) -> &Matrix<i64> {
        &self.commas
    }

    /// A basis for the enharmonic lattice: the intervals `temperament` calls a
    /// unison but this notation spells apart.
    ///
    /// This is the freedom the notation has over the temperament, so it has
    /// `rank() - temperament.rank()` members and is empty exactly when the two
    /// ranks agree and spelling is a bijection. Every notation of an equal
    /// temperament has one, since a notation is free on its octave, its fifth
    /// and its accidentals by construction and so can never close the circle of
    /// fifths: the enharmonic of 12et is the pythagorean comma, which is what
    /// leaves `C#` and `Db` to differ.
    ///
    /// The basis is whatever falls out of the kernel computation, so it is not
    /// reduced to small intervals; [`to_interval`](Self::to_interval) puts each
    /// one back in notation coordinates, where the `[0, 1, -13]` of 22et reads
    /// as its fifth being thirteen of its accidental.
    ///
    /// # Errors
    /// Returns [`Error::InvalidSubgroup`] if `temperament` is over a different
    /// subgroup than this notation.
    pub fn enharmonics(&self, temperament: &Temperament) -> Result<Matrix<i64>, Error> {
        if temperament.subgroup() != &self.subgroup {
            return Err(Error::InvalidSubgroup(format!(
                "this notation is over {}, but the temperament is over {}",
                self.subgroup,
                temperament.subgroup()
            )));
        }
        // The combinations of the generators the temperament sends to nothing.
        let images = temperament.map_all(&self.generators)?;
        Ok(kernel_left(&images)?
            .iter()
            .map(|counts| combination(counts, &self.generators, self.dim()))
            .collect())
    }

    /// Whether every prime is spelled on the nominal just intonation gives it.
    ///
    /// A notation that drops an accidental has to respell the prime it belonged
    /// to. Where the replacement is a stack of the accidentals that are left,
    /// the nominal is kept and only the number of marks changes; where the fifth
    /// chain has to be walked instead, the prime can land on another letter.
    /// 41et writing `5/4` as `Fb5` and `11/8` as `D###5` is what that looks
    /// like, against the `vE5` and `^^F5` of the notation above it.
    ///
    /// Sharps and flats do not count, since seven fifths leave the letter alone:
    /// flattone writes `11/8` as `F#5` where just intonation has `^F5`, and that
    /// is the same nominal.
    ///
    /// # Errors
    /// Returns [`Error::Unsupported`] if some prime has no accidental.
    pub fn keeps_nominals(&self) -> Result<bool, Error> {
        for index in 2..self.dim() {
            let mut prime = vec![0i64; self.dim()];
            prime[index] = 1;
            let fifths = self.to_interval(&prime)?[1];
            let just = just_nominal(&accidental(&self.subgroup, index)?, index);
            if fifths.rem_euclid(NOMINALS.len() as i64) != just {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// The just intonation subgroup being notated.
    pub fn subgroup(&self) -> &Subgroup {
        &self.subgroup
    }

    /// The number of notation coordinates.
    pub fn rank(&self) -> usize {
        self.mapping.len()
    }

    /// The rank of the subgroup being notated, i.e. the length of the interval
    /// vectors this notation maps.
    pub fn dim(&self) -> usize {
        self.subgroup.dim()
    }

    /// Rewrites a just interval in notation coordinates.
    ///
    /// # Errors
    /// Returns [`Error::InvalidDimensions`] if `interval` does not have one
    /// entry per basis element of the subgroup.
    pub fn to_interval(&self, interval: &[i64]) -> Result<Vec<i64>, Error> {
        if interval.len() != self.dim() {
            return Err(Error::InvalidDimensions(format!(
                "interval has {} entries, expected {}",
                interval.len(),
                self.dim()
            )));
        }
        Ok(self
            .mapping
            .iter()
            .map(|row| row.iter().zip(interval).map(|(a, b)| a * b).sum())
            .collect())
    }

    /// Writes notation coordinates as a note in scientific pitch notation,
    /// such as `C5`, `Eb4` or `vE5`.
    ///
    /// The coordinates are read as an interval up from `C5`, so `(0, 0, ..)`
    /// is `C5` itself and `(0, 1, 0, ..)`, a fifth up, is `G5`. Accidentals
    /// come before the nominal and sharps and flats after it. The symbols are
    /// placeholders for debugging, since real microtonal accidentals are not in
    /// unicode.
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

        let mut name = String::new();
        for (&count, &(up, down)) in interval[2..].iter().zip(&ACCIDENTAL_SYMBOLS) {
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

/// Chooses the accidental for the prime at `index` of `subgroup`: the smallest
/// detour along the fifth chain that lands within [`MAX_ACCIDENTAL_CENTS`] of
/// the prime.
///
/// The candidates are the prime shifted by 0, 1, -1, 2, -2, .. fifths, each
/// reduced by octaves to within a tritone of unison. The first one small
/// enough wins, returned as an ascending interval. This gives `81/80` for 5,
/// `64/63` for 7 and `33/32` for 11.
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
        let temperament = Temperament::et(divisions, &subgroup).unwrap();
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
        n.note(&n.to_interval(&interval).unwrap())
    }

    #[test]
    fn pythagorean_mapping() {
        let n = notation("2.3");
        // 2/1 is an octave and 3/1 is an octave plus a fifth.
        assert_eq!(n.mapping(), &vec![vec![1, 1], vec![0, 1]]);
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
                assert_eq!(n.to_interval(generator).unwrap(), unit);
            }
        }
    }

    #[test]
    fn derived_accidentals() {
        // The three the fifth-chain search is meant to produce.
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
        assert_eq!(n.to_interval(&[0, 0, 1]).unwrap(), vec![0, 4, -1]);
        assert_eq!(
            n.mapping(),
            &vec![vec![1, 1, 0], vec![0, 1, 4], vec![0, 0, -1]]
        );
    }

    #[test]
    fn to_interval_checks_dimensions() {
        let n = notation("2.3");
        assert!(n.to_interval(&[1, 0, 0]).is_err());
        assert!(n.to_interval(&[1]).is_err());
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
    fn notes_with_accidentals() {
        // The just third, seventh and eleventh are each a pythagorean interval
        // bent by one accidental.
        assert_eq!(note_of(&notation("2.3.5"), 5, 4), "vE5");
        assert_eq!(note_of(&notation("2.3.7"), 7, 4), "vBb5");
        assert_eq!(note_of(&notation("2.3.11"), 11, 8), "^F5");

        // Over the full subgroup they keep those spellings, and the symbols
        // are assigned in order of the primes.
        let n = notation("2.3.5.7.11");
        assert_eq!(note_of(&n, 5, 4), "vE5");
        assert_eq!(note_of(&n, 7, 4), "<Bb5");
        assert_eq!(note_of(&n, 11, 8), "+F5");
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
            assert_eq!(n.commas().len(), 1);
            assert_eq!(n.subgroup().to_ratio(&n.commas()[0]).unwrap(), comma);
            // The octave and the fifth are still mapped as they always were.
            assert_eq!(n.mapping()[0][..2], pythagorean.mapping()[0][..]);
            assert_eq!(n.mapping()[1][..2], pythagorean.mapping()[1][..]);
        }
    }

    #[test]
    fn tempered_primes_land_on_nominals() {
        // The meantone third is a plain E and the archytas seventh a plain Bb,
        // where just intonation needs an accidental on each.
        let meantone = tempered("2.3.5", &[(81, 80)]);
        assert_eq!(meantone.to_interval(&[0, 0, 1]).unwrap(), vec![0, 4]);
        assert_eq!(note_of(&meantone, 5, 4), "E5");
        assert_eq!(note_of(&meantone, 81, 64), "E5");

        let archytas = tempered("2.3.7", &[(64, 63)]);
        assert_eq!(archytas.to_interval(&[0, 0, 1]).unwrap(), vec![4, -2]);
        assert_eq!(note_of(&archytas, 7, 4), "Bb5");
        assert_eq!(note_of(&archytas, 16, 9), "Bb5");
    }

    #[test]
    fn only_tempered_accidentals_are_dropped() {
        // Marvel tempers out 225/224, which is neither accidental, so both stay
        // and the notation is the just intonation one.
        let marvel = tempered("2.3.5.7", &[(225, 224)]);
        assert_eq!(marvel.rank(), 4);
        assert_eq!(marvel, notation("2.3.5.7"));

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
        let n = Notation::from_temperament(&Temperament::et(12, &subgroup).unwrap()).unwrap();
        assert_eq!(n, tempered("2.3.5", &[(81, 80)]));
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
                options_of(subgroup, &[comma])[0],
                tempered(subgroup, &[comma])
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
        assert_eq!(note_of(&options[2], 11, 8), ">F5");

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
            let temperament = Temperament::et(divisions, &subgroup).unwrap();
            assert!(Notation::options(&temperament).is_err());
        }

        // The answer is a subgroup whose accidental does fit. 24et over 2.3.5 is
        // contorted and saturates to 12et, but over 2.3.5.11 its quartertone is
        // 33/32 and worth exactly one step.
        let options = et_options(24, "2.3.5.11");
        assert_eq!(ranks(&options), vec![3]);
        assert_eq!(accidental_ratios(&options[0]), vec![(33, 32)]);
        assert_eq!(note_of(&options[0], 5, 4), "E5");
        assert_eq!(note_of(&options[0], 11, 8), "^F5");
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
        assert_eq!(note_of(&options[2], 11, 8), ">F5");
    }

    #[test]
    fn kernels_nest_along_the_run() {
        // Every comma of a larger notation is a comma of every smaller one, so
        // a larger notation's spelling can always be simplified onto a smaller
        // one's. This is what fixing one comma per accidental buys.
        for (subgroup, commas) in [
            ("2.3.5.7.11", &[(81, 80), (126, 125), (99, 98)][..]),
            ("2.3.5.7.11", &[(225, 224), (385, 384)][..]),
            ("2.3.5.7", &[(225, 224)][..]),
        ] {
            let options = options_of(subgroup, commas);
            for pair in options.windows(2) {
                for comma in pair[1].commas() {
                    assert_eq!(pair[0].to_interval(comma).unwrap(), vec![0; pair[0].rank()]);
                }
            }
        }
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
    fn the_recommendation_is_the_smallest_that_keeps_the_nominals() {
        // 41et has three notations. The smallest spells every prime off its
        // nominal, and the largest only turns the two marks of the middle one
        // into one of another kind, so the middle one is what is wanted.
        let subgroup: Subgroup = "2.3.5.7.11".parse().unwrap();
        let temperament = Temperament::et(41, &subgroup).unwrap();
        let options = Notation::options(&temperament).unwrap();
        assert_eq!(ranks(&options), vec![2, 3, 4]);
        assert_eq!(
            Notation::from_temperament(&temperament).unwrap(),
            options[1]
        );

        // Where in the run it lands is not fixed: septimal meantone, marvel and
        // schismatic want the second of two, huygens the third of three.
        for (subgroup, commas, wanted) in [
            ("2.3.5.7", &[(81, 80), (225, 224)][..], 1),
            ("2.3.5.7.11", &[(81, 80), (126, 125), (99, 98)][..], 2),
            ("2.3.5.7", &[(225, 224)][..], 1),
            ("2.3.5", &[(32805, 32768)][..], 1),
        ] {
            assert_eq!(
                tempered(subgroup, commas),
                options_of(subgroup, commas)[wanted]
            );
        }
    }

    #[test]
    fn a_sharp_is_the_same_nominal() {
        // Flattone writes 11/8 as F#, where just intonation has ^F. Seven fifths
        // leave the letter alone, so that is the same nominal and the smaller
        // notation is kept rather than taking on an accidental for 11.
        let flattone = &[(45, 44), (81, 80)][..];
        let options = options_of("2.3.5.11", flattone);
        assert_eq!(ranks(&options), vec![2, 3]);
        assert_eq!(note_of(&options[0], 11, 8), "F#5");
        assert_eq!(note_of(&options[1], 11, 8), "^F5");
        assert!(options[0].keeps_nominals().unwrap());
        assert_eq!(tempered("2.3.5.11", flattone), options[0]);
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
    fn a_run_with_a_choice_always_keeps_its_nominals_somewhere() {
        // The recommendation falls back to the last notation where nothing keeps
        // its nominals, which happens only where the fifth chain has to be walked
        // to replace every accidental - 13et over 2.3.5 is the smallest such.
        // Every one of those offers a single notation, so the fallback never
        // decides anything.
        for subgroup in ["2.3.5", "2.3.5.7", "2.3.5.7.11"] {
            let subgroup: Subgroup = subgroup.parse().unwrap();
            for divisions in 5..=72 {
                let t = Temperament::et(divisions, &subgroup).unwrap();
                let Ok(options) = Notation::options(&t) else {
                    continue;
                };
                assert!(
                    options.len() == 1 || options.iter().any(|n| n.keeps_nominals().unwrap()),
                    "{divisions}et over {subgroup} offers a choice but keeps no nominals"
                );
            }
        }
    }

    #[test]
    fn too_many_accidentals() {
        assert!(Notation::from_ji(&Subgroup::p_limit(17)).is_err());
    }

    #[test]
    fn the_circle_of_fifths_never_closes() {
        // An equal temperament must temper out the pythagorean comma, and no
        // notation can: a notation is free on its octave, its fifth and its
        // accidentals. That difference is the enharmonic lattice, and in 12et
        // it is exactly what leaves C# and Db to differ.
        let subgroup = Subgroup::p_limit(5);
        let t = Temperament::et(12, &subgroup).unwrap();
        let options = Notation::options(&t).unwrap();
        assert_eq!(ranks(&options), vec![2]);

        let n = &options[0];
        let enharmonics = n.enharmonics(&t).unwrap();
        assert_eq!(enharmonics.len(), 1);
        assert_eq!(
            subgroup.to_ratio(&enharmonics[0]).unwrap(),
            (531441, 524288)
        );
        // Twelve fifths less seven octaves, in notation coordinates.
        assert_eq!(n.to_interval(&enharmonics[0]).unwrap(), vec![-7, 12]);
    }

    #[test]
    fn enharmonics_are_the_freedom_over_the_temperament() {
        // 22et: the fifth chain alone leaves one enharmonic, and the notation
        // with an accidental leaves two - its octave being twenty two steps and
        // its fifth thirteen, which is the whole of ups and downs in 22et.
        let subgroup: Subgroup = "2.3.5.7".parse().unwrap();
        let t = Temperament::et(22, &subgroup).unwrap();
        let options = Notation::options(&t).unwrap();

        assert_eq!(options[0].enharmonics(&t).unwrap().len(), 1);
        let coordinates: Vec<Vec<i64>> = options[1]
            .enharmonics(&t)
            .unwrap()
            .iter()
            .map(|e| options[1].to_interval(e).unwrap())
            .collect();
        assert_eq!(coordinates, vec![vec![0, 1, -13], vec![1, 0, -22]]);

        // Every notation in the run, of any temperament: an enharmonic is
        // tempered out, is not spelled as a unison, and the lattice has exactly
        // the rank the notation has over the temperament.
        for (divisions, subgroup) in [(12, "2.3.5"), (22, "2.3.5.7"), (41, "2.3.5.7.11")] {
            let subgroup: Subgroup = subgroup.parse().unwrap();
            let t = Temperament::et(divisions, &subgroup).unwrap();
            for n in &Notation::options(&t).unwrap() {
                let enharmonics = n.enharmonics(&t).unwrap();
                assert_eq!(enharmonics.len(), n.rank() - t.rank());
                for e in &enharmonics {
                    assert_eq!(t.map(e).unwrap(), vec![0; t.rank()]);
                    assert!(n.to_interval(e).unwrap().iter().any(|&x| x != 0));
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
        let n = &Notation::options(&t).unwrap()[0];
        assert_eq!(n.rank(), t.rank());
        assert!(n.enharmonics(&t).unwrap().is_empty());
    }

    #[test]
    fn enharmonics_need_the_matching_subgroup() {
        let n = notation("2.3.5");
        let other: Subgroup = "2.3.7".parse().unwrap();
        let t = Temperament::et(12, &other).unwrap();
        assert!(n.enharmonics(&t).is_err());
    }
}
