//! Notation systems as linear maps.

use diophantine::{Matrix, integer_inverse, transpose};

use crate::Error;
use crate::primes::Subgroup;

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
/// accidental to `(0, 0, 1, ..)`. Here every prime beyond 3 gets an accidental
/// of its own, so the map is a bijection: every interval of the subgroup has
/// exactly one spelling, and no two intervals are spelled alike.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Notation {
	mapping: Matrix<i64>,
	generators: Matrix<i64>,
	subgroup: Subgroup,
}

impl Notation {
	/// Builds the notation over `subgroup`, choosing an accidental for each
	/// prime beyond 3.
	///
	/// # Errors
	/// Returns [`Error::Unsupported`] if the subgroup has more primes beyond 3
	/// than there are accidental symbols, or if some prime has no accidental.
	pub fn new(subgroup: &Subgroup) -> Result<Self, Error> {
		let dim = subgroup.dim();
		if dim - 2 > ACCIDENTAL_SYMBOLS.len() {
			return Err(Error::Unsupported(format!(
				"notation over {subgroup} needs {} accidentals, but only {} are named",
				dim - 2,
				ACCIDENTAL_SYMBOLS.len()
			)));
		}

		// Rows are the notational generators as interval vectors: the octave,
		// the fifth, then one accidental per remaining prime.
		let mut generators = vec![vec![0i64; dim]; dim];
		generators[0][0] = 1;
		generators[1][0] = -1;
		generators[1][1] = 1;
		for (index, row) in generators.iter_mut().enumerate().skip(2) {
			*row = accidental(subgroup, index)?;
		}

		// The mapping is whatever sends each generator to its unit vector, i.e.
		// `mapping * transpose(generators) = I`. Inverting is only this direct
		// because the generators here span the whole subgroup; a notation that
		// spelled two intervals alike would need its commas to fill out the
		// square matrix first.
		let mapping = transpose(&integer_inverse(&generators)?);

		Ok(Notation {
			mapping,
			generators,
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
		Notation::new(&subgroup.parse::<Subgroup>().unwrap()).unwrap()
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
	fn too_many_accidentals() {
		assert!(Notation::new(&Subgroup::p_limit(17)).is_err());
	}
}
