//! The search behind [`Notation::options`].
//!
//! The whole of the public surface is that one function, so this module is
//! tested through it rather than on its own.

use diophantine::Matrix;

use crate::Error;
use crate::notation::{Notation, accidental, fifth_chain, spans};
use crate::temperament::Temperament;
use crate::util::{is_zero, select};

/// The derivation of the notations one temperament offers.
///
/// Everything the search does hangs off the accidentals and what each is maps
/// to in the temperament, so those are worked out once here and the rest of the
/// search is methods on them. Accidentals are passed around as indices into
/// [`Self::accidentals`] throughout.
pub(crate) struct Search<'a> {
    temperament: &'a Temperament,
    /// One accidental per prime beyond 3, in order of the primes.
    accidentals: Matrix<i64>,
    /// What each accidental maps to in the temperament.
    images: Matrix<i64>,
    /// The accidentals the temperament does not temper out. One it does temper
    /// out would raise by nothing, so it is never worth keeping.
    useful: Vec<usize>,
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
        let candidates = self.candidates();
        let necessary = self.necessary(&candidates)?;

        // All other ones are optional, unless they map to the same interval
        // as one before it.
        let mut optional: Vec<usize> = Vec::new();
        for index in candidates {
            if necessary.contains(&index) {
                continue;
            }
            let mut present = necessary.iter().chain(&optional);
            if !present.any(|&kept| equal_up_to_sign(&self.images[kept], &self.images[index])) {
                optional.push(index);
            }
        }

        (0..=optional.len())
            .map(|extras| {
                // The accidentals kept by the notation that has taken on `extras` of the optional ones, in order.
                let mut keep: Vec<usize> = necessary
                    .iter()
                    .chain(&optional[..extras])
                    .copied()
                    .collect();
                keep.sort_unstable();
                Notation::from_accidentals(self.temperament, &select(&self.accidentals, &keep))
            })
            .collect()
    }

    /// The accidentals a notation may use, in order.
    ///
    /// For an equal temperament, we first look for an accidental that maps to
    /// one step, since that one is preferred over all others.
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

    /// The smallest subset of `candidates` that makes a notation possible at all.
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
                "the fifth chain of this equal temperament does not reach every note, and no accidental of {subgroup} maps to a single step of it"
            ));
        }
        Error::Unsupported(format!(
            "the octave, the fifth and the accidentals of {subgroup} do not reach every pitch of this rank {} temperament",
            self.temperament.rank()
        ))
    }
}

/// Two accidentals do the same thing if they agree up to sign.
fn equal_up_to_sign(one: &[i64], other: &[i64]) -> bool {
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
