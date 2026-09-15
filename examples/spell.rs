//! Spells a few just intervals that a temperament calls one pitch.
//!
//! The point is what a notation can and cannot decide. Asked for a pitch it can
//! decide nothing: the spellings of one pitch are a coset of the enharmonic
//! lattice and the notation has no opinion about which member to show. Asked
//! for an interval it decides everything, and the answer carries the harmonic
//! reading the temperament threw away.

use xen_utils::{Notation, Simplifier, Subgroup, Temperament};

fn main() {
    for (divisions, subgroup, ratios) in [
        (
            12,
            "2.3.5",
            vec![(16, 15), (25, 24), (135, 128), (9, 8), (10, 9)],
        ),
        (
            41,
            "2.3.5.7.11",
            vec![(14, 11), (81, 64), (11, 7), (128, 81)],
        ),
    ] {
        let subgroup: Subgroup = subgroup.parse().unwrap();
        let t = Temperament::et(divisions, &subgroup).unwrap();
        let n = Notation::from_temperament(&t).unwrap();
        let simplifier = Simplifier::new(&Notation::options(&t).unwrap()[0]).unwrap();

        println!(
            "\n{divisions}et over {subgroup}, notation of rank {}",
            n.rank()
        );
        println!(
            "{:>10}  {:>5}  {:>8}  {:>10}",
            "interval", "steps", "spelled", "simplest"
        );
        for (num, den) in ratios {
            let interval = subgroup.factorize(num, den).unwrap();
            let steps = t.map(&interval).unwrap()[0];
            let simplified = simplifier.simplify(&interval).unwrap();
            let (sn, sd) = subgroup.to_ratio(&simplified).unwrap();
            println!(
                "{:>10}  {steps:>5}  {:>8}  {:>10}",
                format!("{num}/{den}"),
                n.note(&n.to_notation(&interval).unwrap()),
                format!("{sn}/{sd}")
            );
        }
    }
}
