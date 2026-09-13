//! Notation systems as linear maps.

use diophantine::{Matrix, integer_det, integer_inverse, transpose};

use crate::Error;
use crate::primes::Subgroup;
use crate::temperament::Temperament;

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

	/// Builds the notation of `temperament`, over the subgroup it tempers.
	///
	/// Both settings start from the just intonation notation and drop
	/// accidentals the temperament has made redundant. They differ in how far
	/// they go, and there is no one right answer: a temperament may well have
	/// more than one notation worth using.
	///
	/// With `minimal` false, only accidentals the temperament tempers out are
	/// dropped. Such an accidental raises and lowers by nothing, so the prime
	/// it was written for is already sitting on the fifth chain. Meantone and
	/// archytas come out notated exactly as pythagorean is. Anything else the
	/// temperament does is left alone, which is always a correct notation but
	/// not always the intended one.
	///
	/// With `minimal` true, as many accidentals are dropped as the temperament
	/// will allow: the result has the same rank as the temperament, so
	/// spelling is a bijection and there are no enharmonics at all. This gets
	/// septimal meantone written with nothing but sharps and flats, where the
	/// other setting keeps an accidental for 7 - both spell the same pitch,
	/// `A#5` one way and `vBb5` the other. It also fails far more readily,
	/// since most temperaments have no such notation.
	///
	/// # Errors
	/// Returns [`Error::Unsupported`] if the notation needs more accidentals
	/// than there are symbols, if some prime has no accidental, or if
	/// `minimal` was asked for and no set of accidentals gives it.
	pub fn from_temperament(temperament: &Temperament, minimal: bool) -> Result<Self, Error> {
		let subgroup = temperament.subgroup();
		let accidentals = (2..subgroup.dim())
			.map(|index| accidental(subgroup, index))
			.collect::<Result<Matrix<i64>, Error>>()?;

		if minimal {
			return Notation::minimal(temperament, &accidentals);
		}

		let (kept, commas) = accidentals.into_iter().partition(|a| {
			temperament
				.map(a)
				.is_ok_and(|steps| steps.iter().any(|&s| s != 0))
		});
		Notation::assemble(subgroup, kept, commas)
	}

	/// Builds the notation of `temperament` that has no enharmonics: its
	/// kernel is the whole comma lattice, so it has the same rank as the
	/// temperament and spelling is a bijection.
	///
	/// The octave and the fifth are always generators, so such a notation
	/// keeps exactly `rank - 2` accidentals, and which ones is a real choice.
	/// A set works when the generators and the commas together still span the
	/// subgroup, i.e. when the square matrix of the two is unimodular, which
	/// is the same as saying the chosen generators generate the temperament.
	/// Sets keeping the accidentals of the lower primes are tried first.
	fn minimal(temperament: &Temperament, accidentals: &Matrix<i64>) -> Result<Self, Error> {
		let subgroup = temperament.subgroup();
		let rank = temperament.rank();
		let keep = rank
			.checked_sub(2)
			.filter(|&keep| keep <= accidentals.len());
		let keep = keep.ok_or_else(|| {
			Error::Unsupported(format!(
				"a rank {rank} notation over {subgroup} cannot be built from the octave, \
				 the fifth and {} accidentals",
				accidentals.len()
			))
		})?;
		// The reduced basis, so that the commas come out small and ascending
		// like the accidentals do. Any basis of the same lattice would give the
		// same notation.
		let commas = temperament.reduced_comma_basis()?;

		// Subsets of the accidentals, smallest indices first, so that the
		// accidentals of the lower primes are the ones kept where there is a
		// choice.
		for mask in 0u32..(1 << accidentals.len()) {
			if mask.count_ones() as usize != keep {
				continue;
			}
			let kept: Matrix<i64> = accidentals
				.iter()
				.enumerate()
				.filter(|(index, _)| mask & (1 << index) != 0)
				.map(|(_, a)| a.clone())
				.collect();

			let mut square = fifth_chain(subgroup.dim());
			square.extend(kept.iter().chain(&commas).cloned());
			if integer_det(&square)?.abs() == 1 {
				return Notation::assemble(subgroup, kept, commas);
			}
		}

		Err(Error::Unsupported(format!(
			"no {keep} of the accidentals of {subgroup} generate this rank {rank} temperament"
		)))
	}

	/// Builds the notation with the given accidentals and commas, deriving the
	/// mapping from them.
	fn assemble(
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
	/// placeholders for debugging; see [`ACCIDENTAL_SYMBOLS`].
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
fn fifth_chain(dim: usize) -> Matrix<i64> {
	let mut octave = vec![0i64; dim];
	octave[0] = 1;
	let mut fifth = vec![0i64; dim];
	(fifth[0], fifth[1]) = (-1, 1);
	vec![octave, fifth]
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
fn accidental(subgroup: &Subgroup, index: usize) -> Result<Vec<i64>, Error> {
	let offsets = std::iter::once(0).chain((1..=MAX_FIFTH_OFFSET).flat_map(|k| [k, -k]));
	for offset in offsets {
		let mut interval = vec![0i64; subgroup.dim()];
		interval[index] = 1;
		interval[1] = offset;

		// Octave reduction: whichever power of two lands closest to unison.
		interval[0] = -(subgroup.to_cents(&interval) / 1200.0).round_ties_even() as i64;

		let cents = subgroup.to_cents(&interval);
		if cents.abs() < MAX_ACCIDENTAL_CENTS {
			if cents < 0.0 {
				for exponent in &mut interval {
					*exponent = -*exponent;
				}
			}
			return Ok(interval);
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
	fn tempered(subgroup: &str, commas: &[(u64, u64)], minimal: bool) -> Notation {
		let subgroup: Subgroup = subgroup.parse().unwrap();
		let commas: Vec<Vec<i64>> = commas
			.iter()
			.map(|&(num, den)| subgroup.factorize(num, den).unwrap())
			.collect();
		let temperament = Temperament::from_commas(&commas, &subgroup).unwrap();
		Notation::from_temperament(&temperament, minimal).unwrap()
	}

	/// Whether the minimal notation of that temperament exists.
	fn has_minimal(subgroup: &str, commas: &[(u64, u64)]) -> bool {
		let subgroup: Subgroup = subgroup.parse().unwrap();
		let commas: Vec<Vec<i64>> = commas
			.iter()
			.map(|&(num, den)| subgroup.factorize(num, den).unwrap())
			.collect();
		let temperament = Temperament::from_commas(&commas, &subgroup).unwrap();
		Notation::from_temperament(&temperament, true).is_ok()
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
			let n = tempered(subgroup, &[comma], false);
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
		let meantone = tempered("2.3.5", &[(81, 80)], false);
		assert_eq!(meantone.to_interval(&[0, 0, 1]).unwrap(), vec![0, 4]);
		assert_eq!(note_of(&meantone, 5, 4), "E5");
		assert_eq!(note_of(&meantone, 81, 64), "E5");

		let archytas = tempered("2.3.7", &[(64, 63)], false);
		assert_eq!(archytas.to_interval(&[0, 0, 1]).unwrap(), vec![4, -2]);
		assert_eq!(note_of(&archytas, 7, 4), "Bb5");
		assert_eq!(note_of(&archytas, 16, 9), "Bb5");
	}

	#[test]
	fn only_tempered_accidentals_are_dropped() {
		// Marvel tempers out 225/224, which is neither accidental, so both stay
		// and the notation is the just intonation one.
		let marvel = tempered("2.3.5.7", &[(225, 224)], false);
		assert_eq!(marvel.rank(), 4);
		assert_eq!(marvel, notation("2.3.5.7"));

		// Septimal meantone tempers out 81/80 as well, dropping just that one.
		let n = tempered("2.3.5.7", &[(81, 80), (225, 224)], false);
		assert_eq!(accidental_ratios(&n), vec![(64, 63)]);
		assert_eq!(note_of(&n, 5, 4), "E5");
		assert_eq!(note_of(&n, 7, 4), "vBb5");
	}

	#[test]
	fn equal_temperaments_keep_the_fifth_chain() {
		// 12et tempers out 81/80, so it is notated as meantone is: rank 2 over
		// a rank 1 temperament, which is what leaves C# and Db to differ.
		let subgroup = Subgroup::p_limit(5);
		let n =
			Notation::from_temperament(&Temperament::et(12, &subgroup).unwrap(), false).unwrap();
		assert_eq!(n, tempered("2.3.5", &[(81, 80)], false));
		assert_eq!(note_of(&n, 5, 4), "E5");
	}

	#[test]
	fn minimal_drops_an_accidental_that_is_not_tempered_out() {
		// Septimal meantone has two notations. Keeping the accidental for 7
		// writes the harmonic seventh as a lowered Bb; the minimal one writes
		// the same pitch with sharps and flats alone, as A#.
		let commas = [(81, 80), (225, 224)];
		let plain = tempered("2.3.5.7", &commas, false);
		let minimal = tempered("2.3.5.7", &commas, true);

		assert_eq!(plain.rank(), 3);
		assert_eq!(note_of(&plain, 7, 4), "vBb5");

		assert_eq!(minimal.rank(), 2);
		assert!(accidental_ratios(&minimal).is_empty());
		assert_eq!(note_of(&minimal, 7, 4), "A#5");
		assert_eq!(note_of(&minimal, 5, 4), "E5");
	}

	#[test]
	fn minimal_has_the_rank_of_the_temperament() {
		// Schismatic reaches 5 along the fifth chain too, eight fifths down,
		// so its just third is written as a diminished fourth.
		let schismatic = tempered("2.3.5", &[(32805, 32768)], true);
		assert_eq!(schismatic.rank(), 2);
		assert_eq!(note_of(&schismatic, 5, 4), "Fb5");

		// Marvel is rank 3, so its minimal notation keeps one accidental. Only
		// the one for 5 works; keeping the one for 7 leaves 5 unreachable.
		let marvel = tempered("2.3.5.7", &[(225, 224)], true);
		assert_eq!(marvel.rank(), 3);
		assert_eq!(accidental_ratios(&marvel), vec![(81, 80)]);
		assert_eq!(note_of(&marvel, 5, 4), "vE5");
		assert_eq!(note_of(&marvel, 7, 4), "vvA#5");

		// Where rule one already reaches the rank of the temperament, both
		// settings agree.
		for (subgroup, comma) in [("2.3.5", (81, 80)), ("2.3.7", (64, 63))] {
			assert_eq!(
				tempered(subgroup, &[comma], true),
				tempered(subgroup, &[comma], false)
			);
		}
	}

	#[test]
	fn minimal_fails_when_the_fifth_chain_does_not_reach() {
		// An equal temperament is rank 1, and a notation needs the octave and
		// the fifth at the very least. 12et is exactly the case where the
		// enharmonics are the point.
		let subgroup = Subgroup::p_limit(5);
		let et = Temperament::et(12, &subgroup).unwrap();
		assert!(Notation::from_temperament(&et, true).is_err());
		assert!(Notation::from_temperament(&et, false).is_ok());

		// Blackwood is rank 2, but its fifth chain closes after five notes, so
		// the octave and the fifth alone do not span it.
		assert!(!has_minimal("2.3.5", &[(256, 243)]));

		// Porcupine reaches 5 only through its own comma, not the fifth chain.
		assert!(!has_minimal("2.3.5", &[(250, 243)]));
	}

	#[test]
	fn too_many_accidentals() {
		assert!(Notation::from_ji(&Subgroup::p_limit(17)).is_err());
	}
}
