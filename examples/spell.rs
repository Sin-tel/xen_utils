//! Spells a list of just intervals in the recommended notation.
//!
//! Run it, change the derivation, run it again and diff: the intervals whose
//! spelling moves are exactly the ones the comma choice reaches. Two notations
//! built on the same generators agree on everything those generators span and
//! can only disagree beyond it.

use xen_utils::{Notation, Simplifier, Subgroup, Temperament};

fn main() {
    for (divisions, subgroup, rank, ratios) in [
        (
            12,
            "2.3.5",
            2,
            vec![(16, 15), (25, 24), (135, 128), (9, 8), (10, 9)],
        ),
        (
            41,
            "2.3.5.7.11",
            3,
            vec![
                // the 5-limit, which the octave, the fifth and 81/80 span
                (9, 8),
                (5, 4),
                (6, 5),
                (45, 32),
                (5, 3),
                (15, 8),
                (25, 16),
                (81, 64),
                // and beyond it, where the commas decide
                (7, 4),
                (7, 6),
                (7, 5),
                (11, 8),
                (11, 9),
                (14, 11),
                (11, 7),
            ],
        ),
    ] {
        let subgroup: Subgroup = subgroup.parse().unwrap();
        let t = Temperament::et(divisions, &subgroup).unwrap();
        let options = Notation::options(&t).unwrap();
        // Pinned by rank, so that two derivations are compared at the same
        // notation rather than at whichever each recommends.
        let n = options.iter().find(|n| n.rank() == rank).unwrap();
        let simplifier = Simplifier::new(&options[0]).unwrap();

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
