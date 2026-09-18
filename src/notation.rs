//! Notation systems as linear maps.

use diophantine::{Matrix, cvp_exact, eye, kernel_left, lll, solve_diophantine, transpose};

use crate::Error;
use crate::notation_options::NotationOptions;
use crate::primes::Subgroup;
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
const PRIME_SYMBOLS: [(u32, char, char); 6] = [
    (5, '^', 'v'),
    (7, '>', '<'),
    (11, 't', 'd'),
    (13, '*', '%'),
    (17, '/', '\\'),
    (19, ')', '('),
];

/// Arbitrary symbols for primes beyond [`PRIME_SYMBOLS`], handed out in order.
/// Past these, pairs of Greek letters and then CJK ideographs, which never run out.
const FALLBACK_SYMBOLS: [(char, char); 5] =
    [('!', '?'), ('[', ']'), ('{', '}'), ('+', '~'), ('\'', ',')];

/// Characters [`Notation::note`] uses for something else, so no accidental may.
const RESERVED_SYMBOLS: &str = "FCGDAEB#b-0123456789";

/// The largest interval, in cents, that counts as an accidental:
/// half an apotome, about 56.8 cents.
///
/// This lets every prime be written without augmented or diminished intervals.
/// 1200*log2(sqrt(2187/2048))
const MAX_ACCIDENTAL_CENTS: f64 = 56.842_503_028_855_52;

/// An accidental, raising or lowering a note by a small interval.
#[derive(Debug, Clone)]
pub struct Accidental {
    /// The interval vector of the accidental.
    pub vector: Vec<i64>,
    /// The raising and lowering symbols.
    pub symbols: (char, char),
}

/// A notation system: a set of symbols, and what each one maps to in the temperament.
///
/// Notation coordinates are counts of notational generators. The first two are
/// always the octave `2/1` and the fifth `3/2`, which together give the
/// nominals and the sharps and flats. Other coordinates count accidentals,
/// which raises or lowers by a small interval that is not a sharp.
#[derive(Debug, Clone)]
pub struct Notation {
    /// The octave, the fifth, then the accidentals, as prime interval vectors.
    generators: Matrix<i64>,
    /// The accidental symbols.
    accidentals: Vec<(char, char)>,
    /// What each generator maps to in the temperament, one per row.
    /// This is the map from notation coordinates to tempered intervals.
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
    /// Returns [`Error::Unsupported`] if some prime has no accidental.
    pub fn from_ji(subgroup: &Subgroup) -> Result<Self, Error> {
        let accidentals = derive_accidentals(subgroup)?;
        Notation::from_accidentals(&Temperament::from_ji(subgroup)?, &accidentals)
    }

    /// Builds the recommended notation of `temperament`: the smallest of
    /// the run [`options`](Self::options) gives that
    /// [keeps every nominal](Self::keeps_nominals).
    ///
    /// # Errors
    /// Returns [`Error::Unsupported`] if some prime has no accidental, or if no
    /// notation can be derived at all.
    pub fn from_temperament(temperament: &Temperament) -> Result<Self, Error> {
        Self::from_temperament_with(temperament, &derive_accidentals(temperament.subgroup())?)
    }

    /// Builds the recommended notation of `temperament` given a set of `accidentals`.
    pub fn from_temperament_with(
        temperament: &Temperament,
        accidentals: &[Accidental],
    ) -> Result<Self, Error> {
        let options = Notation::options_with(temperament, accidentals)?;
        for option in &options {
            if option.keeps_nominals()? {
                return Ok(option.clone());
            }
        }
        Ok(options
            .first()
            .expect("there is always at least one option")
            .clone())
    }

    /// The best notation `temperament` offers of each size, smallest first.
    ///
    /// Every subset of the accidentals is a candidate. An accidental tempered
    /// out is never kept, nor one worth what an earlier accidental is worth,
    /// up to direction: in 41et `64/63` is worth what `81/80` is. A subset
    /// must reach every pitch, and an equal temperament keeping any accidental
    /// must keep one worth a single step.
    ///
    /// Of each size the best is the one with fewest primes off their nominals
    /// (see [`keeps_nominals`](Self::keeps_nominals)), then the lowest total
    /// [`nominal_costs`](Self::nominal_costs), then those costs compared one
    /// prime at a time, lowest first. The best of one size need not keep what
    /// the best of the size below keeps.
    ///
    /// # Errors
    /// Returns [`Error::Unsupported`] if some prime has no accidental, or if
    /// no notation exists at all - an equal temperament whose fifth chain does
    /// not reach every note and which has no accidental worth a single step of
    /// it has none.
    pub fn options(temperament: &Temperament) -> Result<Vec<Self>, Error> {
        Self::options_with(temperament, &derive_accidentals(temperament.subgroup())?)
    }

    /// [`options`](Self::options) over the given `accidentals`.
    pub fn options_with(
        temperament: &Temperament,
        accidentals: &[Accidental],
    ) -> Result<Vec<Self>, Error> {
        NotationOptions::new(temperament, accidentals)?.search()
    }

    /// The best notation of `temperament` keeping exactly `count` of
    /// `accidentals`, ranked as [`options`](Self::options) ranks them.
    ///
    /// # Errors
    /// Returns [`Error::Unsupported`] if no `count` of them make a notation.
    pub fn with_count(
        temperament: &Temperament,
        accidentals: &[Accidental],
        count: usize,
    ) -> Result<Self, Error> {
        NotationOptions::new(temperament, accidentals)?.with_count(count)
    }

    /// Builds the notation of `temperament` with `accidentals` as its extra
    /// generators.
    ///
    /// This is internal to crate since it doesn't do any checks.
    ///
    /// # Errors
    /// Returns [`Error::InvalidDimensions`] if an accidental is not an interval
    /// of the temperament's subgroup, and [`Error::Unsupported`] if two
    /// symbols collide, or if the generators do not reach every pitch of the
    /// temperament.
    pub(crate) fn from_accidentals(
        temperament: &Temperament,
        accidentals: &[Accidental],
    ) -> Result<Self, Error> {
        let subgroup = temperament.subgroup();
        check_symbols(accidentals)?;

        let mut generators = fifth_chain(subgroup.dim());
        generators.extend(accidentals.iter().map(|a| a.vector.clone()));
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
            accidentals: accidentals.iter().map(|a| a.symbols).collect(),
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

    /// The temperament being notated.
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
        self.respell_within(coordinates, count, RESPELL_WIDTH)
    }

    /// [`respell`](Self::respell), walking `width` enharmonics either way.
    #[doc(hidden)]
    pub fn respell_within(
        &self,
        coordinates: &[i64],
        count: usize,
        width: i64,
    ) -> Result<Vec<Vec<i64>>, Error> {
        self.check(coordinates)?;
        Ok(cheapest(
            coordinates,
            &self.enharmonics,
            &self.weights,
            width,
            count,
        ))
    }

    /// The cheapest way to write each prime beyond 3 on its just nominal: at
    /// the degree just intonation gives it, or `None` where no spelling of its
    /// pitch has that degree.
    ///
    /// The degree is exact rather than modulo the seven nominals: a note ten
    /// sharps up an octave lower is not on the nominal, though its letter is.
    ///
    /// This asks for the prime on its letter, whatever [`spell`](Self::spell)
    /// would choose: 41et's largest notation spells `7/4` as `tA`, but writes
    /// it on its nominal as `vBb`.
    ///
    /// # Errors
    /// Returns [`Error::Unsupported`] if some prime has no accidental.
    pub fn nominal_spellings(&self) -> Result<Vec<Option<Vec<i64>>>, Error> {
        // The notation coordinates and their degree side by side, so that one
        // solve finds a spelling of the right pitch at the right degree.
        let with_degree: Matrix<i64> = self
            .images
            .iter()
            .enumerate()
            .map(|(index, image)| {
                let mut row = image.clone();
                row.push(GENERATOR_DEGREES.get(index).copied().unwrap_or(0));
                row
            })
            .collect();
        // The enharmonics that keep the degree, reduced for the search.
        let level = kernel_left(&with_degree)?;
        let level = lll(&level, LLL_DELTA, &self.weights).unwrap_or(level);

        let mut spellings = Vec::new();
        for (index, nominal) in (2..self.dim()).zip(just_nominals(self.subgroup())?) {
            let mut prime = vec![0i64; self.dim()];
            prime[index] = 1;
            let mut target = self.temperament.map(&prime)?;
            target.push(nominal.degree);
            let spelling = solve_diophantine(&transpose(&with_degree), &column(&target))
                .ok()
                .map(|solution| {
                    let seed = first_column(&solution);
                    cheapest(&seed, &level, &self.weights, RESPELL_WIDTH, 1).remove(0)
                });
            spellings.push(spelling);
        }
        Ok(spellings)
    }

    /// What [`nominal_spellings`](Self::nominal_spellings) cost to read.
    ///
    /// # Errors
    /// Returns [`Error::Unsupported`] if some prime has no accidental.
    pub fn nominal_costs(&self) -> Result<Vec<Option<i64>>, Error> {
        Ok(self
            .nominal_spellings()?
            .iter()
            .map(|spelling| spelling.as_deref().map(spelling_cost))
            .collect())
    }

    /// Whether every prime beyond 3 can be written on the nominal just
    /// intonation gives it, at a cost no worse than either just intonation
    /// spends on it or this notation spends on its cheapest spelling of it.
    ///
    /// # Errors
    /// Returns [`Error::Unsupported`] if some prime has no accidental.
    pub fn keeps_nominals(&self) -> Result<bool, Error> {
        Ok(self.nominal_verdict()?.failures == 0)
    }

    /// [`nominal_costs`](Self::nominal_costs), and how many primes fail
    /// [`keeps_nominals`](Self::keeps_nominals).
    pub(crate) fn nominal_verdict(&self) -> Result<NominalVerdict, Error> {
        let costs = self.nominal_costs()?;
        let nominals = just_nominals(self.subgroup())?;
        let mut failures = 0;
        for ((index, cost), nominal) in (2..self.dim()).zip(&costs).zip(nominals) {
            let mut prime = vec![0i64; self.dim()];
            prime[index] = 1;
            let best = spelling_cost(&self.spell(&prime)?);
            if cost.is_none_or(|cost| cost > best.max(nominal.cost)) {
                failures += 1;
            }
        }
        Ok(NominalVerdict { costs, failures })
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

        let mut name = String::new();
        for (&count, (up, down)) in interval[2..].iter().zip(&self.accidentals) {
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

/// How well a notation keeps its nominals.
pub(crate) struct NominalVerdict {
    /// [`Notation::nominal_costs`].
    pub(crate) costs: Vec<Option<i64>>,
    /// How many primes fail [`Notation::keeps_nominals`].
    pub(crate) failures: usize,
}

/// Where just intonation writes a prime: its degree, and what that costs.
#[derive(Debug, Clone, Copy)]
pub(crate) struct JustNominal {
    /// Letters up from the unison, seven to the octave.
    pub(crate) degree: i64,
    /// [`spelling_cost`] of the just spelling.
    pub(crate) cost: i64,
}

/// The degree of the octave and the fifth. Every accidental has degree zero.
const GENERATOR_DEGREES: [i64; 2] = [7, 4];

/// The degree of notation coordinates: how many letters up they sit.
pub(crate) fn degree(coordinates: &[i64]) -> i64 {
    GENERATOR_DEGREES
        .iter()
        .zip(coordinates)
        .map(|(d, c)| d * c)
        .sum()
}

/// Where just intonation writes each prime beyond 3.
///
/// Just intonation writes a prime as its own accidental on top of a stretch of
/// the fifth chain. The accidental `a` has exponent `s = +-1` on the prime, so
/// the prime is `s * a` less `s * a[1]` threes and `s * a[0]` twos: the fifth
/// coordinate `-s * a[1]`, the octave coordinate `-s * (a[0] + a[1])` and one
/// mark. Together with 7 for the octave and 11 for the three, these degrees
/// are a linear map from interval vectors to degrees, `(7, 11, 16, 20, 24, ..)`,
/// whose kernel holds the apotome and every derived accidental.
pub(crate) fn just_nominals(subgroup: &Subgroup) -> Result<Vec<JustNominal>, Error> {
    (2..subgroup.dim())
        .map(|index| {
            let a = derive_accidental_vector(subgroup, index)?;
            let s = a[index];
            let spelling = [-s * (a[0] + a[1]), -s * a[1], 1];
            Ok(JustNominal {
                degree: degree(&spelling),
                cost: spelling_cost(&spelling),
            })
        })
        .collect()
}

/// Derives the default accidentals for a subgroup.
pub fn derive_accidentals(subgroup: &Subgroup) -> Result<Vec<Accidental>, Error> {
    let mut result = Vec::new();
    let mut fallbacks = 0;
    for index in 2..subgroup.dim() {
        let vector = derive_accidental_vector(subgroup, index)?;
        let prime = subgroup.basis()[index];
        let symbols = match PRIME_SYMBOLS.iter().find(|&&(p, ..)| p == prime) {
            Some(&(_, up, down)) => (up, down),
            None => {
                fallbacks += 1;
                fallback_symbols(fallbacks - 1)
            }
        };
        result.push(Accidental { vector, symbols });
    }
    Ok(result)
}

/// The `n`th pair of arbitrary symbols for a prime with no symbol of its own.
fn fallback_symbols(n: usize) -> (char, char) {
    const GREEK_PAIRS: usize = 12;
    let pair = |start: u32, k: usize| {
        let up = char::from_u32(start + 2 * k as u32).expect("a valid code point");
        let down = char::from_u32(start + 2 * k as u32 + 1).expect("a valid code point");
        (up, down)
    };
    match n.checked_sub(FALLBACK_SYMBOLS.len()) {
        None => FALLBACK_SYMBOLS[n],
        Some(k) if k < GREEK_PAIRS => pair('α' as u32, k),
        Some(k) => pair(0x4E00, k - GREEK_PAIRS),
    }
}

/// Checks that every symbol of `accidentals` is distinct, and that none of
/// them is a character [`Notation::note`] already uses.
fn check_symbols(accidentals: &[Accidental]) -> Result<(), Error> {
    let symbols: Vec<char> = accidentals
        .iter()
        .flat_map(|a| [a.symbols.0, a.symbols.1])
        .collect();
    for (index, &symbol) in symbols.iter().enumerate() {
        if RESERVED_SYMBOLS.contains(symbol) {
            return Err(Error::Unsupported(format!(
                "accidental symbol '{symbol}' is already used for nominals, sharps, flats or octaves"
            )));
        }
        if symbols[..index].contains(&symbol) {
            return Err(Error::Unsupported(format!(
                "accidental symbol '{symbol}' is used more than once"
            )));
        }
    }
    Ok(())
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
pub(crate) fn derive_accidental_vector(
    subgroup: &Subgroup,
    index: usize,
) -> Result<Vec<i64>, Error> {
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

/// The middle of the seven naturals, as a fifth coordinate.
///
/// `F C G D A E B` are the fifth coordinates `-1` to `5`, so `D` at `2` is the
/// middle of them.
const NOMINAL_CENTRE: i64 = 2;

/// How far Notation::respell searches for.
const RESPELL_WIDTH: i64 = 2;

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

/// The `count` cheapest spellings in the coset of `coordinates` modulo
/// `lattice`, best first: a closest vector under `weights`, then a walk of
/// `width` steps either way along each member of `lattice` under
/// [`spelling_cost`] itself.
fn cheapest(
    coordinates: &[i64],
    lattice: &Matrix<i64>,
    weights: &Matrix<f64>,
    width: i64,
    count: usize,
) -> Vec<Vec<i64>> {
    if lattice.is_empty() {
        return vec![coordinates.to_vec()];
    }
    let rank = coordinates.len();

    // `solve_diophantine` hands back an arbitrary vector, so reduce first.
    let mut centre = coordinates.to_vec();
    centre[1] -= NOMINAL_CENTRE;
    let seed: Vec<i64> = match cvp_exact(&centre, lattice, weights) {
        Ok(near) => coordinates.iter().zip(near).map(|(a, b)| a - b).collect(),
        Err(_) => coordinates.to_vec(),
    };

    let mut found = Vec::new();
    let mut candidate = vec![0; rank];
    let mut steps = vec![-width; lattice.len()];
    loop {
        for slot in 0..rank {
            candidate[slot] = seed[slot]
                + steps
                    .iter()
                    .zip(lattice)
                    .map(|(&step, row)| step * row[slot])
                    .sum::<i64>();
        }
        found.push(candidate.clone());

        let mut place = 0;
        while place < steps.len() && steps[place] == width {
            steps[place] = -width;
            place += 1;
        }
        if place == steps.len() {
            break;
        }
        steps[place] += 1;
    }

    // Compared rather than keyed, so that ranking a candidate does not clone
    // its coordinates to break the tie with.
    found.sort_unstable_by(|a, b| {
        spelling_cost(a)
            .cmp(&spelling_cost(b))
            .then_with(|| spelling_break_ties(a).cmp(&spelling_break_ties(b)))
            .then_with(|| a.cmp(b))
    });
    found.dedup();
    found.truncate(count.max(1));
    found
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
                .map(|col| match (row == col, row) {
                    (false, _) => 0.0,
                    // Register costs nothing. The form is then only
                    // semidefinite, but no enharmonic is a stack of octaves,
                    // so it is definite on the enharmonic lattice.
                    (true, 0) => 0.0,
                    (true, 1) => (COST_FIFTH * COST_FIFTH) as f64,
                    (true, _) => (COST_MARK * COST_MARK) as f64,
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

    /// Every notation the equal temperament of `divisions` over `subgroup` offers.
    fn et_options(divisions: i64, subgroup: &str) -> Vec<Notation> {
        let subgroup: Subgroup = subgroup.parse().unwrap();
        let temperament = Temperament::equal(divisions, &subgroup).unwrap();
        Notation::options(&temperament).unwrap()
    }

    /// Every notation the temperament of `subgroup` tempering out `commas` offers.
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
        assert_eq!(n.note(&[-3, 5]), "B4");
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
        assert_eq!(note_of(&n, 81, 64), "E5");
        assert_eq!(note_of(&n, 2187, 2048), "C#5");
    }

    #[test]
    fn accidentals_match_fjs() {
        let n = notation("2.3.5.7.11.13.17.19");

        let fjs_accidentals = [
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
        assert_eq!(note_of(&notation("2.3.5"), 5, 4), "vE5");
        assert_eq!(note_of(&notation("2.3.7"), 7, 4), "<Bb5");
        assert_eq!(note_of(&notation("2.3.11"), 11, 8), "tF5");

        let n = notation("2.3.5.7.11");
        assert_eq!(note_of(&n, 5, 4), "vE5");
        assert_eq!(note_of(&n, 7, 4), "<Bb5");
        assert_eq!(note_of(&n, 11, 8), "tF5");
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
        assert_eq!(note_of(&n, 7, 4), "<Bb5");
    }

    #[test]
    fn notation_12et_is_meantone() {
        // 12et tempers out 81/80, so it has the same notation as meantone.
        let subgroup = Subgroup::p_limit(5);
        let n = Notation::from_temperament(&Temperament::equal(12, &subgroup).unwrap()).unwrap();
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
        assert_eq!(note_of(&options[0], 5, 4), "E5");
        assert_eq!(note_of(&options[1], 5, 4), "E5");
        assert_eq!(note_of(&options[0], 7, 4), "A#5");
        assert_eq!(note_of(&options[1], 7, 4), "<Bb5");
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
        assert_eq!(note_of(&options[1], 7, 4), "<Bb5");
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
        assert_eq!(note_of(&options[1], 11, 8), ">F5");
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
        assert_eq!(note_of(&options[1], 7, 4), "<Bb5");
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
    fn spelling_the_spelled_settles() {
        // The best spelling of a pitch is already `spell`'s fixed point:
        // reading it back as a just interval and spelling again should not
        // move it - the same property `Simplifier::simplify` settles into on
        // its own answer.
        for (divisions, subgroup) in [(41, "2.3.5.7.11"), (31, "2.3.5.7"), (22, "2.3.5")] {
            let notation = et_options(divisions, subgroup).pop().unwrap();
            for ups in -20..=60 {
                let mut spelling = vec![0; notation.rank()];
                spelling[2] = ups;
                let interval = notation.to_just(&spelling).unwrap();

                let spelled = notation.spell(&interval).unwrap();
                let round_trip = notation
                    .spell(&notation.to_just(&spelled).unwrap())
                    .unwrap();
                assert_eq!(
                    round_trip, spelled,
                    "{divisions}et over {subgroup}, {ups} ups"
                );
            }
        }
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
        let accidentals = derive_accidentals(&subgroup).unwrap();

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

        let syntonic = derive_accidentals(&subgroup).unwrap()[0].clone();
        let notation = Notation::from_accidentals(&temperament, &[syntonic]).unwrap();
        assert_eq!(notation.rank(), 3);
        assert_eq!(note_of(&notation, 5, 4), "vE5");
    }

    /// Every octave-free interval with exponents in -2..=2 beyond the octave.
    fn small_intervals(dim: usize) -> Vec<Vec<i64>> {
        (0..5i64.pow(dim as u32 - 1))
            .map(|mut code| {
                let mut interval = vec![0];
                for _ in 1..dim {
                    interval.push(code % 5 - 2);
                    code /= 5;
                }
                interval
            })
            .collect()
    }

    #[test]
    fn register_does_not_change_the_spelling() {
        // Octaves cost nothing to write, so moving a pitch by octaves moves
        // only its octave coordinate. The seed once charged for them and got
        // this wrong far from C5.
        for n in et_options(22, "2.3.5.7")
            .into_iter()
            .chain(options_of("2.3.5.7.11", &[(441, 440), (896, 891)]))
        {
            for interval in small_intervals(n.dim()) {
                let spelling = n.spell(&interval).unwrap();
                for octaves in [-12, -3, 4, 9] {
                    let mut moved = interval.clone();
                    moved[0] += octaves;
                    let mut expected = spelling.clone();
                    expected[0] += octaves;
                    assert_eq!(n.spell(&moved).unwrap(), expected);
                }
            }
        }
    }

    #[test]
    fn respell_walks_wide_enough() {
        // Pele's rank 4 notation has pitches whose best spelling a walk of one
        // enharmonic either way misses, `E##` found as `vdF##`. Two finds them all.
        let options = options_of("2.3.5.7.11", &[(441, 440), (896, 891)]);
        let n = &options[1];
        let mut missed_at_one = 0;
        for interval in small_intervals(n.dim()) {
            let spelling = n.spell(&interval).unwrap();
            let wide = n.respell_within(&spelling, 3, 5).unwrap();
            assert_eq!(n.respell(&spelling, 3).unwrap(), wide);
            if n.respell_within(&spelling, 1, 1).unwrap()[0] != wide[0] {
                missed_at_one += 1;
            }
        }
        assert!(missed_at_one > 0);
    }

    #[test]
    fn just_nominals_are_the_diatonic_degrees() {
        // (7, 11, 16, 20, 24, 26): third, seventh, fourth and sixth, as just
        // intonation writes them, one mark each on vE, <Bb, tF and *Ab.
        let subgroup: Subgroup = "2.3.5.7.11.13".parse().unwrap();
        let nominals = just_nominals(&subgroup).unwrap();
        let degrees: Vec<i64> = nominals.iter().map(|n| n.degree).collect();
        let costs: Vec<i64> = nominals.iter().map(|n| n.cost).collect();
        assert_eq!(degrees, vec![16, 20, 24, 26]);
        assert_eq!(costs, vec![11, 15, 13, 19]);
        // Every derived accidental has degree zero under it.
        let map: Vec<i64> = [7, 11].into_iter().chain(degrees).collect();
        for a in derive_accidentals(&subgroup).unwrap() {
            let dot: i64 = a.vector.iter().zip(&map).map(|(x, d)| x * d).sum();
            assert_eq!(dot, 0);
        }
    }

    #[test]
    fn nominal_costs_look_past_the_cheapest_spelling() {
        // Schismatic's fifth chain puts 5 on F, and nowhere else.
        let options = options_of("2.3.5", &[(32805, 32768)]);
        assert_eq!(options[0].nominal_costs().unwrap(), vec![None]);
        assert_eq!(options[1].nominal_costs().unwrap(), vec![Some(11)]);

        // Flattone writes 11 as F#, which is on F.
        let options = options_of("2.3.5.11", &[(45, 44), (81, 80)]);
        assert_eq!(options[0].nominal_costs().unwrap(), vec![Some(4), Some(8)]);

        // 41et's largest notation writes 7/4 as tA, but vBb is there too, at
        // what just intonation spends on it, so the nominals are kept.
        let options = et_options(41, "2.3.5.7.11");
        assert_eq!(note_of(&options[2], 7, 4), "tA5");
        assert_eq!(options[2].nominal_costs().unwrap()[1], Some(15));
        assert!(options[2].keeps_nominals().unwrap());
    }

    #[test]
    fn a_cheaper_spelling_elsewhere_does_not_hide_the_nominal() {
        // Semaphore and pele each have a spelling cheaper than the one on the
        // nominal, which the old rule, reading only `spell`, took as a miss.
        let semaphore = options_of("2.3.7", &[(49, 48)]);
        assert!(semaphore[0].keeps_nominals().unwrap());

        let pele = tempered("2.3.5.7.11", &[(441, 440), (896, 891)]);
        assert_eq!(accidental_ratios(&pele), vec![(81, 80), (33, 32)]);
        assert!(pele.keeps_nominals().unwrap());
    }

    #[test]
    fn a_tie_goes_to_the_lower_primes() {
        // Miracle's two accidental notations tie three ways: whichever pair is
        // kept, one of 5, 7 and 11 needs two marks. The simpler 5 and then 7
        // are preferred, so it is 11 that takes them.
        let miracle = &[(225, 224), (1029, 1024), (385, 384)][..];
        let options = options_of("2.3.5.7.11", miracle);
        assert_eq!(accidental_ratios(&options[1]), vec![(81, 80), (64, 63)]);
        assert_eq!(note_of(&options[1], 11, 8), "^>F5");
        assert_eq!(
            tempered("2.3.5.7.11", miracle).generators(),
            options[1].generators()
        );
    }

    #[test]
    fn a_notation_of_a_given_size() {
        let subgroup: Subgroup = "2.3.5.7.11".parse().unwrap();
        let commas = [(225, 224), (1029, 1024), (385, 384)]
            .map(|(n, d)| subgroup.factorize(n, d).unwrap());
        let miracle = Temperament::from_commas(&commas, &subgroup).unwrap();
        let accidentals = derive_accidentals(&subgroup).unwrap();

        // The fifth chain alone does not reach every pitch of miracle.
        assert!(Notation::with_count(&miracle, &accidentals, 0).is_err());
        let one = Notation::with_count(&miracle, &accidentals, 1).unwrap();
        assert_eq!(accidental_ratios(&one), vec![(81, 80)]);
        let three = Notation::with_count(&miracle, &accidentals, 3).unwrap();
        assert_eq!(three.rank(), 5);
        assert!(Notation::with_count(&miracle, &accidentals, 4).is_err());

        // Only the images matter: in 41et 49/48 is a step as good as 81/80.
        let subgroup: Subgroup = "2.3.5.7".parse().unwrap();
        let t = Temperament::equal(41, &subgroup).unwrap();
        let septimal = Accidental {
            vector: subgroup.factorize(49, 48).unwrap(),
            symbols: ('^', 'v'),
        };
        let n = Notation::with_count(&t, &[septimal], 1).unwrap();
        let syntonic = Notation::with_count(&t, &derive_accidentals(&subgroup).unwrap()[..1], 1).unwrap();
        assert_eq!(n.enharmonics(), syntonic.enharmonics());
        assert_eq!(n.nominal_costs().unwrap(), syntonic.nominal_costs().unwrap());
    }

    #[test]
    fn primes_beyond_the_table_get_arbitrary_symbols() {
        // Past 19 the symbols are arbitrary, but every one is still distinct:
        // enough primes to run through the ASCII and Greek pools and into CJK.
        let n = Notation::from_ji(&Subgroup::p_limit(199)).unwrap();
        assert_eq!(n.accidentals.len(), 44);
        assert!(check_symbols(&derive_accidentals(n.subgroup()).unwrap()).is_ok());
        assert_eq!(note_of(&notation("2.3.23"), 23, 16), "!Gb5");
    }

    #[test]
    fn colliding_symbols_are_refused() {
        let subgroup: Subgroup = "2.3.5.7".parse().unwrap();
        let temperament = Temperament::from_ji(&subgroup).unwrap();
        let mut accidentals = derive_accidentals(&subgroup).unwrap();
        assert!(Notation::from_accidentals(&temperament, &accidentals).is_ok());

        // The same symbol for two accidentals.
        accidentals[1].symbols = ('^', '<');
        assert!(Notation::from_accidentals(&temperament, &accidentals).is_err());
        // Raising and lowering with one symbol.
        accidentals[1].symbols = ('>', '>');
        assert!(Notation::from_accidentals(&temperament, &accidentals).is_err());
        // A flat.
        accidentals[1].symbols = ('>', 'b');
        assert!(Notation::from_accidentals(&temperament, &accidentals).is_err());
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
        let bare = Notation::from_accidentals(&t, &[]).unwrap();

        let syntonic = derive_accidentals(&subgroup).unwrap()[0].clone();
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
            let accidentals = derive_accidentals(&subgroup).unwrap();
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
        assert_eq!(spelling_cost(&[5, 2, 0, 0, 0]), 0);

        for i in 1..rank {
            // Relative to D5
            let mut interval = vec![-1, 2, 0, 0, 0];
            interval[i] += 1;

            let w_l1 = spelling_cost(&interval);
            let w_l2 = weights[i][i];

            assert_eq!((w_l1 * w_l1) as f64, w_l2);
        }
    }

    #[test]
    fn custom_accidental_25edo() {
        let subgroup: Subgroup = "2.3.5".parse().unwrap();
        let temperament = Temperament::equal(25, &subgroup).unwrap();

        // Fails on defaults because 81/80 does not reach the whole temperament
        // or is not worth a step (25edo over 2.3.5 is known to have no default notation).
        assert!(Notation::options(&temperament).is_err());

        // But we can supply 25/24 as a custom accidental.
        let custom_acc = Accidental {
            vector: subgroup.factorize(25, 24).unwrap(),
            symbols: ('^', 'v'),
        };

        let options = Notation::options_with(&temperament, &[custom_acc]).unwrap();
        assert_eq!(ranks(&options), vec![3]);
        assert_eq!(note_of(&options[0], 5, 4), "vvE5");
    }

    #[test]
    fn custom_accidental_41edo() {
        let subgroup: Subgroup = "2.3.5.7".parse().unwrap();
        let temperament = Temperament::equal(41, &subgroup).unwrap();

        // We can notate using a single accidental for 49/48 or 50/49.
        let acc_49_48 = Accidental {
            vector: subgroup.factorize(49, 48).unwrap(),
            symbols: ('^', 'v'),
        };

        let acc_50_49 = Accidental {
            vector: subgroup.factorize(50, 49).unwrap(),
            symbols: ('^', 'v'),
        };

        let options_49_48 =
            Notation::options_with(&temperament, std::slice::from_ref(&acc_49_48)).unwrap();
        assert_eq!(options_49_48.last().unwrap().rank(), 3);

        let options_50_49 =
            Notation::options_with(&temperament, std::slice::from_ref(&acc_50_49)).unwrap();
        assert_eq!(options_50_49.last().unwrap().rank(), 3);
    }
}
