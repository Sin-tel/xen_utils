//! The candidate **spellings** of a pitch, for an end user to cycle through.
//!
//! One pitch is any number of written notes, and nothing in the temperament
//! prefers one of them: they differ by an enharmonic, which is what the
//! temperament calls a unison and the notation still writes apart. So this is a
//! question the notation has to answer with a ranking of its own, and the
//! ranking is what a note costs to write - seven half-apotomes for an
//! accidental mark and two for a fifth away from the middle of the naturals.
//!
//! This is one of the two questions a pitch has; `cargo run --example
//! candidates` is the other. They are separate: which just interval is meant
//! and how the note is written are decided by different things, and 12et
//! writing `D` is not an opinion about `9/8` against `10/9`.

use xen_utils::{Notation, Simplifier, Subgroup, Temperament};

// TODO: `reached` breaks on valid ETs.
const DIVISIONS: i64 = 41;
const SUBGROUP: &str = "2.3.5.7.11";
const SHOWN: usize = 5;

fn main() {
    for (divisions, subgroup) in [(12, "2.3.5"), (DIVISIONS, SUBGROUP)] {
        let subgroup: Subgroup = subgroup.parse().unwrap();
        let temperament = Temperament::equal(divisions, &subgroup).unwrap();
        let options = Notation::options(&temperament).unwrap();
        let notation = Notation::from_temperament(&temperament).unwrap();
        let simplifier = Simplifier::new(&options[0]).unwrap();

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
            let interval = reached(&subgroup, &temperament, steps, divisions);
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

/// Some just interval worth `step` steps: the shortest stack of fifths that
/// reaches it, octave reduced so that the notes come out in one register. Any one will do: every spelling of it is in the same
/// coset, and the simplifier reduces every interval worth the same alike.
fn reached(subgroup: &Subgroup, t: &Temperament, step: i64, divisions: i64) -> Vec<i64> {
    let mut fifth = vec![0i64; subgroup.dim()];
    (fifth[0], fifth[1]) = (-1, 1);
    let per_fifth = t.map(&fifth).unwrap()[0];
    for count in (0..=divisions).flat_map(|k| [k, -k]) {
        let remainder = step - count * per_fifth;
        if remainder.rem_euclid(divisions) == 0 {
            let mut interval = vec![0i64; subgroup.dim()];
            (interval[0], interval[1]) = (-count, count);
            interval[0] -= (t.map(&interval).unwrap()[0] - step) / divisions;
            return interval;
        }
    }
    unreachable!("the fifth chain reaches every step of a notatable equal temperament")
}
