//! The candidate spellings of a pitch in an equal temperament.

use xen_utils::{Notation, Simplifier, Subgroup, Temperament};

const DIVISIONS: i64 = 41;
const SUBGROUP: &str = "2.3.5.7.11";
const SHOWN: usize = 3;

fn main() {
    for (divisions, subgroup) in [(12, "2.3.5"), (DIVISIONS, SUBGROUP)] {
        let subgroup: Subgroup = subgroup.parse().unwrap();
        let temperament = Temperament::equal(divisions, &subgroup).unwrap();
        let notation = Notation::from_temperament(&temperament).unwrap();
        let simplifier = Simplifier::new(&temperament).unwrap();

        println!(
            "\n{divisions}et over {subgroup}, notation of rank {}",
            notation.rank()
        );
        println!(
            "enharmonics: {}",
            notation
                .enharmonics()
                .iter()
                .map(|e| format!("{e:?}"))
                .collect::<Vec<_>>()
                .join("  ")
        );
        println!("\n{:>4}  {:>9}  ways to write it", "step", "reading");

        for steps in 0..divisions + 1 {
            let interval = temperament.map_inverse(&[steps]).unwrap();
            let seed = notation.spell(&interval).unwrap();
            let reading = simplifier.simplify(&interval).unwrap();
            let (num, den) = subgroup.to_ratio(&reading).unwrap();

            let ways: Vec<String> = notation
                .respell(&seed, SHOWN)
                .unwrap()
                .iter()
                .map(|c| format!("{:10}", notation.note(c)))
                .collect();
            println!(
                "{steps:4}  {:>9}  {}",
                format!("{num}/{den}"),
                ways.join("")
            );
        }
    }
}
