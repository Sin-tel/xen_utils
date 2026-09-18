use diophantine::Matrix;

use crate::Error;
use crate::notation::{Accidental, Notation, fifth_chain, spans};
use crate::temperament::Temperament;
use crate::util::{is_zero, select};

pub(crate) struct NotationOptions<'a> {
    temperament: &'a Temperament,
    accidentals: Vec<Accidental>,
    /// What each accidental maps to in the temperament.
    images: Matrix<i64>,
}

impl<'a> NotationOptions<'a> {
    pub(crate) fn new(
        temperament: &'a Temperament,
        accidentals: &[Accidental],
    ) -> Result<Self, Error> {
        let derived: Matrix<i64> = accidentals.iter().map(|a| a.vector.clone()).collect();
        let derived_images = temperament.map_all(&derived)?;
        let useful: Vec<usize> = (0..derived.len())
            .filter(|&index| !is_zero(&derived_images[index]))
            .collect();

        let filtered_accidentals: Vec<Accidental> =
            useful.iter().map(|&i| accidentals[i].clone()).collect();

        Ok(NotationOptions {
            temperament,
            accidentals: filtered_accidentals,
            images: select(&derived_images, &useful),
        })
    }

    /// Every notation the temperament offers, smallest first.
    pub(crate) fn search(&self) -> Result<Vec<Notation>, Error> {
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
                let keep_accs: Vec<Accidental> =
                    keep.iter().map(|&i| self.accidentals[i].clone()).collect();
                Notation::from_accidentals(self.temperament, &keep_accs)
            })
            .collect()
    }

    fn candidates(&self) -> Vec<usize> {
        if self.temperament.rank() != 1 {
            return (0..self.accidentals.len()).collect();
        }
        // For an equal temperament, we first look for an accidental that maps to
        // one step, since that one is preferred over all others.
        let step = (0..self.accidentals.len()).find(|&index| self.images[index][0].abs() == 1);
        match step {
            Some(step) => std::iter::once(step)
                .chain((0..self.accidentals.len()).filter(|&index| index != step))
                .collect(),
            None => Vec::new(),
        }
    }

    /// The smallest subset of `candidates` that makes a notation possible at all.
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

    /// Whether the octave, the fifth and the accidentals at `keep` reach every
    /// pitch of the temperament.
    fn spans(&self, keep: &[usize]) -> Result<bool, Error> {
        let mut generators = fifth_chain(self.temperament.dim());
        generators.extend(keep.iter().map(|&i| self.accidentals[i].vector.clone()));
        let images = self.temperament.map_all(&generators)?;
        Ok(spans(&images, self.temperament.rank()))
    }

    /// Why no subset of the candidates reaches every pitch.
    fn nothing_spans(&self) -> Error {
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

fn equal_up_to_sign(one: &[i64], other: &[i64]) -> bool {
    one == other || one.iter().zip(other).all(|(a, b)| *a == -b)
}

/// Subsets of `items` of size `size`, in lexicographic order.
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
