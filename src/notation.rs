//! Notation systems as linear maps.

use diophantine::{
	Matrix, cvp_exact, eye, integer_inverse, kernel_right, lll, solve_diophantine, transpose,
};

use crate::Error;
use crate::primes::{Subgroup, Weighting};
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
	/// A temperament may have more than one notation worth using, so this picks
	/// an end of the run [`options`](Self::options) gives: the smallest
	/// notation with `minimal` true, and the one keeping every accidental
	/// worth keeping with it false.
	///
	/// Septimal meantone is the plain case of the two differing. The smallest
	/// notation writes the harmonic seventh with sharps and flats alone, as
	/// `A#5`; the largest keeps an accidental for 7 and writes the same pitch
	/// `vBb5`.
	///
	/// # Errors
	/// Returns [`Error::Unsupported`] if the notation needs more accidentals
	/// than there are symbols, if some prime has no accidental, or if no
	/// notation can be derived at all.
	pub fn from_temperament(temperament: &Temperament, minimal: bool) -> Result<Self, Error> {
		let options = Notation::options(temperament)?;
		let option = if minimal {
			options.first()
		} else {
			options.last()
		};
		Ok(option.expect("there is always at least one option").clone())
	}

	/// Every notation `temperament` offers, smallest first.
	///
	/// The first is the smallest notation there is, and the last keeps every
	/// accidental worth keeping; which of them to use is a matter of taste, so
	/// the choice is left open. Enabling an accidental splits apart spellings
	/// the one before it wrote alike, at the cost of another symbol to read.
	///
	/// Each accidental falls into one of three classes. One the temperament
	/// tempers out would raise by nothing and is always dropped. Of the rest,
	/// the smallest set that makes the notation span the subgroup at all is
	/// *necessary* and is always kept; the remainder are *optional*, and the
	/// run enables them one at a time in order of prime.
	///
	/// The first notation therefore has the rank of the temperament wherever one
	/// of that rank exists, and the last keeps an accidental for every prime the
	/// temperament does not spell on the fifth chain already.
	///
	/// # Errors
	/// Returns [`Error::Unsupported`] if some prime has no accidental, or if a
	/// notation in the run needs more accidentals than there are symbols.
	pub fn options(temperament: &Temperament) -> Result<Vec<Self>, Error> {
		let subgroup = temperament.subgroup();
		let accidentals: Matrix<i64> = (2..subgroup.dim())
			.map(|index| accidental(subgroup, index))
			.collect::<Result<_, Error>>()?;

		let images: Matrix<i64> = accidentals
			.iter()
			.map(|a| temperament.map(a))
			.collect::<Result<_, Error>>()?;
		// An accidental the temperament tempers out raises by nothing, so it is
		// never worth keeping.
		let useful: Vec<usize> = (0..accidentals.len())
			.filter(|&index| images[index].iter().any(|&step| step != 0))
			.collect();

		let necessary = Notation::necessary(temperament, &accidentals, &useful)?;
		let dropped: Vec<usize> = (0..accidentals.len())
			.filter(|index| !necessary.contains(index))
			.collect();
		let optional: Vec<usize> = dropped
			.iter()
			.copied()
			.filter(|index| useful.contains(index))
			.collect();

		// One comma per accidental that may be dropped, fixed once and for all
		// so that enabling an accidental only ever removes a comma. A smaller
		// notation's kernel therefore always contains a larger one's.
		let commas: Matrix<i64> = dropped
			.iter()
			.map(|&index| Notation::comma(temperament, &accidentals, &necessary, &useful, index))
			.collect::<Result<_, Error>>()?;

		(0..=optional.len())
			.map(|extras| {
				let keep: Vec<usize> = necessary
					.iter()
					.chain(&optional[..extras])
					.copied()
					.collect();
				let kept = (0..accidentals.len())
					.filter(|index| keep.contains(index))
					.map(|index| accidentals[index].clone())
					.collect();
				let commas = dropped
					.iter()
					.zip(&commas)
					.filter(|(index, _)| !keep.contains(index))
					.map(|(_, comma)| comma.clone())
					.collect();
				Notation::assemble(subgroup, kept, commas)
			})
			.collect()
	}

	/// Which accidentals every notation of `temperament` must keep, as indices
	/// into `accidentals`: the smallest subset of the `useful` ones that makes
	/// a notation possible at all, preferring the lower primes.
	///
	/// The octave, the fifth and the kept accidentals have to generate the
	/// temperament, since otherwise some interval the temperament distinguishes
	/// has no spelling. Taking `keep` of them works exactly when their images
	/// generate the whole of the tempered lattice.
	///
	/// A notation of the same rank as the temperament keeps `rank - 2` of them,
	/// so that is the smallest this can come to; where no such subset spans,
	/// the notation is forced to be larger. Keeping every useful accidental
	/// always spans, so the search always finds something.
	fn necessary(
		temperament: &Temperament,
		accidentals: &Matrix<i64>,
		useful: &[usize],
	) -> Result<Vec<usize>, Error> {
		let rank = temperament.rank();
		let identity: Matrix<i64> = eye(rank);

		for keep in 0..=useful.len() {
			// Subsets of `keep` indices, in lexicographic order, so that the
			// accidentals of the lower primes are the ones kept where there is
			// a choice.
			for subset in subsets(useful, keep) {
				let mut generators = fifth_chain(temperament.dim());
				generators.extend(subset.iter().map(|&index| accidentals[index].clone()));
				let images: Matrix<i64> = generators
					.iter()
					.map(|g| temperament.map(g))
					.collect::<Result<_, Error>>()?;
				if solve_diophantine(&transpose(&images), &identity).is_ok() {
					return Ok(subset);
				}
			}
		}

		Err(Error::Unsupported(format!(
			"the octave, the fifth and the accidentals of {} do not generate this rank {rank} \
			 temperament",
			temperament.subgroup()
		)))
	}

	/// The notational comma of the accidental at `index`: the difference
	/// between it and the replacement a notation without it must use.
	///
	/// The replacement is built from the octave, the fifth, the necessary
	/// accidentals and the accidentals of the lower primes - the ones still
	/// available once the higher ones have been dropped in turn - and has to be
	/// worth what the accidental is worth, so that the two are the same pitch.
	/// An accidental the temperament tempers out is replaced by nothing at all,
	/// since the temperament already calls it a unison.
	///
	/// Which replacement is usually still a choice, settled in two steps.
	///
	/// A plain stack of the accidentals still available keeps the nominal and
	/// the sharps that just intonation gives the prime and changes only the
	/// number of accidentals, so where there is such a stack it is the one
	/// wanted. This is what writes `33/32` as two syntonic commas in 41et and
	/// as one septimal comma in 31et.
	///
	/// Failing that the fifth chain has to be walked, and how far is a choice
	/// again. The replacements that work differ by the commas the notation could
	/// temper out, so they form a coset of that lattice, and the one whose comma
	/// is simplest under [`Weighting`] is taken. Preferring a short walk instead
	/// would run the accidentals away: 51et writes `64/63` as a stack of 25
	/// syntonic commas rather than take five fifths.
	fn comma(
		temperament: &Temperament,
		accidentals: &Matrix<i64>,
		necessary: &[usize],
		useful: &[usize],
		index: usize,
	) -> Result<Vec<i64>, Error> {
		// The temperament already spells this accidental as a unison.
		if !useful.contains(&index) {
			return Ok(accidentals[index].clone());
		}

		// The accidentals still available, lowest prime first. One the
		// temperament tempers out raises by nothing, so it is no use here.
		let stack: Matrix<i64> = useful
			.iter()
			.copied()
			.filter(|&other| other != index && (other < index || necessary.contains(&other)))
			.map(|other| accidentals[other].clone())
			.collect();
		let target = temperament.map(&accidentals[index])?;

		let images = Notation::images(temperament, &stack)?;
		if let Some(counts) = preferred_solution(&target, &images)? {
			return Ok(difference(&accidentals[index], &counts, &stack));
		}

		// The same generators with the octave and the fifth added.
		let mut every = stack;
		every.extend(fifth_chain(temperament.dim()));
		let images = Notation::images(temperament, &every)?;
		let columns = transpose(&images);

		let wanted: Matrix<i64> = target.iter().map(|&step| vec![step]).collect();
		let counts = solve_diophantine(&columns, &wanted).map_err(|_| {
			Error::Unsupported(format!(
				"the accidental {:?} of {} cannot be replaced, though it was not kept",
				accidentals[index],
				temperament.subgroup()
			))
		})?;
		let rough = difference(
			&accidentals[index],
			&counts.iter().map(|row| row[0]).collect::<Vec<i64>>(),
			&every,
		);

		let freedom = kernel_right(&columns)?;
		if freedom.is_empty() || freedom[0].is_empty() {
			return Ok(rough);
		}
		let zero = vec![0i64; temperament.dim()];
		let lattice: Matrix<i64> = transpose(&freedom)
			.iter()
			.map(|counts| {
				difference(&zero, counts, &every)
					.iter()
					.map(|x| -x)
					.collect()
			})
			.collect();

		let weights = temperament.subgroup().weights(Weighting::default());
		let reduced = lll(&lattice, 0.99, &weights)?;
		let closest = cvp_exact(&rough, &reduced, &weights)?;
		Ok(rough
			.iter()
			.zip(&closest)
			.map(|(value, near)| value - near)
			.collect())
	}

	/// The images of `generators` under `temperament`.
	fn images(temperament: &Temperament, generators: &Matrix<i64>) -> Result<Matrix<i64>, Error> {
		generators.iter().map(|g| temperament.map(g)).collect()
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

/// The greatest common divisor of `a` and `b`, which is zero only if both are.
fn gcd(a: i64, b: i64) -> i64 {
	let (mut a, mut b) = (a.abs(), b.abs());
	while b != 0 {
		(a, b) = (b, a % b);
	}
	a
}

/// `interval` less the combination of `generators` given by `counts`.
fn difference(interval: &[i64], counts: &[i64], generators: &Matrix<i64>) -> Vec<i64> {
	(0..interval.len())
		.map(|prime| {
			interval[prime]
				- counts
					.iter()
					.zip(generators)
					.map(|(count, g)| count * g[prime])
					.sum::<i64>()
		})
		.collect()
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
		return Ok(target.iter().all(|&x| x == 0).then(Vec::new));
	};
	let rest = rest.to_vec();

	// Reaching the target without the least preferred generator at all.
	if let Some(mut counts) = preferred_solution(target, &rest)? {
		counts.push(0);
		return Ok(Some(counts));
	}

	let columns = transpose(generators);
	let wanted: Matrix<i64> = target.iter().map(|&x| vec![x]).collect();
	let Ok(solution) = solve_diophantine(&columns, &wanted) else {
		return Ok(None);
	};

	// The counts of the last generator that work differ by whatever multiple of
	// it the others can make up for, so they run in steps of `period`.
	let kernel = diophantine::kernel_right(&columns)?;
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

	let reduced: Vec<i64> = target
		.iter()
		.zip(last)
		.map(|(&t, &l)| t - count * l)
		.collect();
	let mut counts =
		preferred_solution(&reduced, &rest)?.expect("the remainder is reachable by construction");
	counts.push(count);
	Ok(Some(counts))
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
		// 41et: the fifth chain, then one accidental worth a step, then 64/63,
		// which is worth a step as well, then 33/32, which is worth two and so
		// is written as two of the first while it is dropped.
		let options = et_options(41, "2.3.5.7.11");
		assert_eq!(ranks(&options), vec![2, 3, 4, 5]);
		assert_eq!(accidental_ratios(&options[1]), vec![(81, 80)]);
		assert_eq!(note_of(&options[1], 7, 4), "vBb5");
		assert_eq!(note_of(&options[1], 11, 8), "^^F5");
		assert_eq!(
			accidental_ratios(&options[3]),
			vec![(81, 80), (64, 63), (33, 32)]
		);
		assert_eq!(note_of(&options[3], 11, 8), "+F5");

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
		// accidental, and its syntonic comma is worth two steps rather than
		// one. Two steps is still enough to reach what the chain misses, since
		// two and the fifteen of the fifth share no factor with 25.
		for divisions in [25, 51, 54] {
			let options = et_options(divisions, "2.3.5");
			assert_eq!(ranks(&options), vec![3]);
			assert_eq!(accidental_ratios(&options[0]), vec![(81, 80)]);
			assert_eq!(note_of(&options[0], 5, 4), "vE5");
		}
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
	fn the_two_settings_are_the_ends_of_the_run() {
		// 41et has four notations; the flag picks the outer two.
		let subgroup: Subgroup = "2.3.5.7.11".parse().unwrap();
		let temperament = Temperament::et(41, &subgroup).unwrap();
		let options = Notation::options(&temperament).unwrap();
		assert_eq!(options.len(), 4);
		assert_eq!(
			Notation::from_temperament(&temperament, true).unwrap(),
			options[0]
		);
		assert_eq!(
			Notation::from_temperament(&temperament, false).unwrap(),
			*options.last().unwrap()
		);
	}

	#[test]
	fn too_many_accidentals() {
		assert!(Notation::from_ji(&Subgroup::p_limit(17)).is_err());
	}
}
