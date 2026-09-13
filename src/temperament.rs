//! Regular temperaments as integer linear maps.

use diophantine::{Matrix, kernel_left, kernel_right, lll, saturation, transpose};

use crate::Error;
use crate::primes::{Subgroup, Weighting};

/// A regular temperament: a linear map from the interval vectors of a just
/// intonation subgroup to a free abelian group of lower rank.
///
/// Internally this is stored as a mapping matrix (rows are generators,
/// columns are basis elements of the subgroup), always kept in canonical
/// form: saturated (i.e. defactored - no contorsion) and reduced to Hermite
/// normal form. Two mapping matrices describe the same temperament exactly
/// when they have the same canonical form, regardless of choice of
/// generators.
///
/// The [`Subgroup`] the mapping is over is carried along with it, since a
/// mapping matrix means nothing without knowing which rationals its columns
/// refer to.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Temperament {
    mapping: Matrix<i64>,
    subgroup: Subgroup,
}

impl Temperament {
    /// Builds a temperament from a mapping matrix over `subgroup` (rows =
    /// generators, columns = basis elements).
    ///
    /// The mapping is canonicalized, so the result does not depend on which
    /// basis of generators the input happened to use.
    pub fn from_mapping(mapping: &Matrix<i64>, subgroup: &Subgroup) -> Result<Self, Error> {
        if mapping.is_empty() {
            return Err(Error::InvalidDimensions("mapping must be non-empty".into()));
        }
        if mapping.iter().any(|row| row.len() != subgroup.dim()) {
            return Err(Error::InvalidDimensions(format!(
                "mapping rows must have {} entries, one per basis element of {subgroup}",
                subgroup.dim()
            )));
        }
        let mapping = saturation(mapping)?;
        Ok(Temperament {
            mapping,
            subgroup: subgroup.clone(),
        })
    }

    /// Builds the temperament over `subgroup` obtained by tempering out the
    /// given commas.
    pub fn from_commas(commas: &[Vec<i64>], subgroup: &Subgroup) -> Result<Self, Error> {
        if commas.is_empty() {
            return Err(Error::InvalidDimensions("need at least one comma".into()));
        }
        if commas.iter().any(|c| c.len() != subgroup.dim()) {
            return Err(Error::InvalidDimensions(format!(
                "commas must have {} entries, one per basis element of {subgroup}",
                subgroup.dim()
            )));
        }

        // kernel_left wants the commas as columns.
        let comma_matrix = transpose(&commas.to_vec());
        let mapping = kernel_left(&comma_matrix)?;
        if mapping.is_empty() {
            return Err(Error::InvalidDimensions(
                "commas span the whole space; resulting temperament has rank 0".into(),
            ));
        }
        Temperament::from_mapping(&mapping, subgroup)
    }

    /// Builds the equal temperament of `edo` divisions of the octave over
    /// `subgroup`: each basis element `p` is mapped to the number of steps of
    /// `1/edo` octaves that best approximates it, found by scaling `log2(p)` by
    /// `edo` and rounding.
    pub fn et(edo: i64, subgroup: &Subgroup) -> Result<Self, Error> {
        let map: Vec<i64> = subgroup
            .basis()
            .iter()
            .map(|&p| (edo as f64 * f64::from(p).log2()).round_ties_even() as i64)
            .collect();
        Temperament::from_mapping(&vec![map], subgroup)
    }

    /// The canonical mapping matrix (rows = generators, columns = basis
    /// elements of [`Self::subgroup`]).
    pub fn mapping(&self) -> &Matrix<i64> {
        &self.mapping
    }

    /// The just intonation subgroup being tempered.
    pub fn subgroup(&self) -> &Subgroup {
        &self.subgroup
    }

    /// The rank of the temperament (number of independent generators).
    pub fn rank(&self) -> usize {
        self.mapping.len()
    }

    /// The rank of the subgroup being tempered, i.e. the length of the
    /// interval vectors this temperament maps.
    pub fn dim(&self) -> usize {
        self.subgroup.dim()
    }

    /// A basis for the lattice of commas tempered out by this temperament.
    ///
    /// This is whatever basis falls out of the kernel computation and is
    /// usually not musically sensible on its own; see [`Self::reduced_comma_basis`].
    pub fn comma_basis(&self) -> Result<Vec<Vec<i64>>, Error> {
        let columns = kernel_right(&self.mapping)?;
        Ok(transpose(&columns))
    }

    /// A basis for the comma lattice, reduced (via LLL under [`Weighting::Wilson`])
    /// to small, musically sensible commas. Each returned comma is normalized to
    /// be greater than unison.
    pub fn reduced_comma_basis(&self) -> Result<Vec<Vec<i64>>, Error> {
        let commas = self.comma_basis()?;
        if commas.is_empty() {
            return Ok(commas);
        }

        let weights = self.subgroup.weights(Weighting::Wilson);
        let reduced = lll(&commas, 0.99, &weights)?;
        Ok(reduced
            .into_iter()
            .map(|comma| {
                if self.subgroup.to_cents(&comma) < 0.0 {
                    comma.iter().map(|x| -x).collect()
                } else {
                    comma
                }
            })
            .collect())
    }

    /// Applies the mapping to an interval, returning its tempered representation
    /// as a vector of generator counts.
    pub fn map(&self, interval: &[i64]) -> Result<Vec<i64>, Error> {
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
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The reduced comma basis as ratios, in the order LLL returns them: the
    /// Lovasz condition fixes that order, so it is part of what is being tested.
    fn comma_ratios(t: &Temperament) -> Vec<(u64, u64)> {
        t.reduced_comma_basis()
            .unwrap()
            .iter()
            .map(|comma| t.subgroup().to_ratio(comma).unwrap())
            .collect()
    }

    #[test]
    fn et_12edo_5limit() {
        let s = Subgroup::p_limit(5);
        // 12edo: 2 -> 12 steps, 3 -> 19 steps, 5 -> 28 steps.
        let t = Temperament::et(12, &s).unwrap();
        assert_eq!(t.mapping(), &vec![vec![12, 19, 28]]);
    }

    #[test]
    fn et_over_subgroup() {
        // 2.3.7 is not a prime limit: 7 is mapped, 5 is not part of the group.
        let s: Subgroup = "2.3.7".parse().unwrap();
        let t = Temperament::et(12, &s).unwrap();
        assert_eq!(t.mapping(), &vec![vec![12, 19, 34]]);
        assert_eq!(t.subgroup(), &s);
        assert_eq!(t.dim(), 3);
        assert_eq!(t.rank(), 1);
    }

    #[test]
    fn mapping_must_match_subgroup() {
        let s = Subgroup::p_limit(5);
        assert!(Temperament::from_mapping(&vec![vec![12, 19]], &s).is_err());
        assert!(Temperament::from_commas(&[vec![-4, 4]], &s).is_err());
    }

    #[test]
    fn reduced_commas_of_12et() {
        let s = Subgroup::p_limit(5);
        let t = Temperament::et(12, &s).unwrap();
        let commas = t.reduced_comma_basis().unwrap();
        // Rank 1 over a rank 3 group, so the comma lattice has rank 2.
        assert_eq!(commas.len(), 2);
        for comma in &commas {
            // Every returned comma is tempered out, and is an ascending interval.
            assert_eq!(t.map(comma).unwrap(), vec![0]);
            assert!(s.to_cents(comma) > 0.0);
            // Reduction should find commas far smaller than an octave.
            assert!(s.to_cents(comma) < 100.0);
        }
    }

    #[test]
    fn reduced_commas_of_miracle() {
        let s = Subgroup::p_limit(11);
        let m31 = Temperament::et(31, &s).unwrap();
        let m41 = Temperament::et(41, &s).unwrap();
        assert_eq!(m31.mapping(), &vec![vec![31, 49, 72, 87, 107]]);
        assert_eq!(m41.mapping(), &vec![vec![41, 65, 95, 115, 142]]);

        // The join of 31et and 41et is 11-limit miracle: stacking the two maps
        // gives a rank 2 temperament tempering out everything both temper out.
        let mapping = vec![m31.mapping()[0].clone(), m41.mapping()[0].clone()];
        let miracle = Temperament::from_mapping(&mapping, &s).unwrap();
        assert_eq!(miracle.rank(), 2);

        // Rank 2 over a rank 5 group, so the comma lattice has rank 3.
        assert_eq!(
            comma_ratios(&miracle),
            vec![(225, 224), (385, 384), (441, 440)]
        );
    }

    #[test]
    fn meantone_from_comma() {
        let s = Subgroup::p_limit(5);
        // Tempering out 81/80 leaves rank 2 meantone.
        let t = Temperament::from_commas(&[vec![-4, 4, -1]], &s).unwrap();
        assert_eq!(t.rank(), 2);
        assert_eq!(t.map(&[-4, 4, -1]).unwrap(), vec![0, 0]);
    }
}
