//! The search behind [`Notation::options`].
//!
//! The whole of the public surface is that one function, so this module is
//! tested through it rather than on its own.
//!
//! The search decides **which accidentals a notation keeps**, and that is all it
//! decides. A notation is its generators: where they sit in the temperament
//! fixes the pitch of every written note and which written notes are the same
//! pitch, and there is no map from just intonation to choose. So the shape of
//! the search is: derive one accidental per prime beyond 3, work out what each
//! is worth to the temperament, and check them one at a time to see which can be
//! dropped.

use diophantine::Matrix;

use crate::Error;
use crate::notation::{Notation, accidental, fifth_chain, spans};
use crate::temperament::Temperament;
use crate::util::{is_zero, select};

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
/// being made: the run follows from this and nothing else.
struct Plan {
    /// Kept by every notation in the run.
    necessary: Vec<usize>,
    /// Taken on one at a time in this order, each splitting apart pitches the
    /// notation before it wrote alike.
    optional: Vec<usize>,
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
        (0..=plan.optional.len())
            .map(|extras| {
                Notation::build(
                    self.temperament,
                    select(&self.accidentals, &plan.kept(extras)),
                )
            })
            .collect()
    }

    /// Sorts the accidentals into the ones every notation keeps and the ones the
    /// run takes on one at a time. Everything else no notation keeps.
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

        Ok(Plan {
            necessary,
            optional,
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
    /// 1 no accidental can reach every pitch by itself, so there is nothing for
    /// "worth one step" to generalise to.
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
    /// so that is the smallest this can come to; where no such subset reaches
    /// every pitch, the notation is forced to be larger. Keeping every candidate
    /// works unless the candidates have been cut down, which only happens for an
    /// equal temperament with no accidental worth a single step.
    fn necessary(&self, candidates: &[usize]) -> Result<Vec<usize>, Error> {
        for size in 0..=candidates.len() {
            // Subsets in lexicographic order, so that the accidentals offered
            // first are the ones kept where there is a choice.
            for subset in subsets(candidates, size) {
                if self.reaches(&subset)? {
                    return Ok(subset);
                }
            }
        }
        Err(self.nothing_reaches())
    }

    /// Whether the octave, the fifth and the accidentals at `keep` reach every
    /// pitch of the temperament. Where they do not, some pitch has no spelling.
    fn reaches(&self, keep: &[usize]) -> Result<bool, Error> {
        let mut generators = fifth_chain(self.temperament.dim());
        generators.extend(select(&self.accidentals, keep));
        let images = self.temperament.map_all(&generators)?;
        Ok(spans(&images, self.temperament.rank()))
    }

    /// Why no subset of the candidates reaches every pitch.
    fn nothing_reaches(&self) -> Error {
        let subgroup = self.temperament.subgroup();
        if self.temperament.rank() == 1 {
            return Error::Unsupported(format!(
                "the fifth chain of this equal temperament does not reach every note, and no accidental of {subgroup} is worth a single step of it"
            ));
        }
        Error::Unsupported(format!(
            "the octave, the fifth and the accidentals of {subgroup} do not reach every pitch of this rank {} temperament",
            self.temperament.rank()
        ))
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
