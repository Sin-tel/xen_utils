//! Rough timings of the operations expected to run on a hot path: `simplify`,
//! `simplifications`, `spell`, `spellings` and `to_interval`. Building a `Notation` or a
//! `Simplifier` is not measured, since that happens once and is reused.

use std::time::Instant;

use xen_utils::{Notation, Simplifier, Subgroup, Temperament};

/// Runs `f` `iterations` times and prints the average time per call.
fn time(name: &str, iterations: u32, mut f: impl FnMut()) {
    let start = Instant::now();
    for _ in 0..iterations {
        f();
    }
    let elapsed = start.elapsed();
    println!(
        "{name:32} {:>10.1} ns/call  ({iterations} calls in {elapsed:?})",
        elapsed.as_nanos() as f64 / f64::from(iterations)
    );
}

fn main() {
    for (label, divisions, subgroup) in [
        ("41et  2.3.5.7.11", 41, "2.3.5.7.11"),
        ("72et  2.3.5.7.11", 72, "2.3.5.7.11"),
        ("41et  2.3.5.7.11.13", 41, "2.3.5.7.11.13"),
    ] {
        println!("\n== {label} ==");
        let subgroup: Subgroup = subgroup.parse().unwrap();
        let temperament = Temperament::equal(divisions, &subgroup).unwrap();
        let notation = Notation::from_temperament(&temperament).unwrap();
        let simplifier = Simplifier::new(&temperament).unwrap();
        println!(
            "notation length {}, comma lattice rank {}",
            notation.len(),
            simplifier.lattice().len()
        );

        // A spread of intervals, not just the unison, so the walk does real work.
        let intervals: Vec<Vec<i64>> = (0..divisions)
            .map(|ups| {
                let mut spelling = vec![0; notation.len()];
                spelling[2] = ups;
                notation.to_interval(&spelling).unwrap()
            })
            .collect();
        let tempered: Vec<Vec<i64>> = intervals
            .iter()
            .map(|interval| temperament.temper(interval).unwrap())
            .collect();
        let spellings: Vec<Vec<i64>> = tempered
            .iter()
            .map(|t| notation.spell(t).unwrap())
            .collect();

        let mut i = 0;
        time("Notation::spell", 10_000, || {
            i = (i + 1) % tempered.len();
            std::hint::black_box(notation.spell(&tempered[i]).unwrap());
        });

        let mut i = 0;
        time("Notation::spellings(4)", 10_000, || {
            i = (i + 1) % tempered.len();
            std::hint::black_box(notation.spellings(&tempered[i], 4).unwrap());
        });

        let mut i = 0;
        time("Notation::to_interval", 1_000_000, || {
            i = (i + 1) % spellings.len();
            std::hint::black_box(notation.to_interval(&spellings[i]).unwrap());
        });

        let mut i = 0;
        time("Simplifier::simplify", 10_000, || {
            i = (i + 1) % tempered.len();
            std::hint::black_box(simplifier.simplify(&tempered[i]).unwrap());
        });

        let mut i = 0;
        time("Simplifier::simplifications(8)", 10_000, || {
            i = (i + 1) % tempered.len();
            std::hint::black_box(simplifier.simplifications(&tempered[i], 8).unwrap());
        });
    }
}
