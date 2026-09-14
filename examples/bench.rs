//! Rough timings of the operations expected to run on a hot path: `simplify`,
//! `candidates`, `to_notation` and `to_just`. Building a `Notation` or a
//! `Simplifier` is not measured, since that happens once and is reused.
//!
//! No harness, no statistics - just wall time over enough iterations that
//! noise averages out. Good enough to tell "negligible" from "worth fixing".

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
        let temperament = Temperament::et(divisions, &subgroup).unwrap();
        let notation = Notation::from_temperament(&temperament).unwrap();
        let simplifier = Simplifier::new(&notation).unwrap();
        println!(
            "notation rank {}, comma lattice rank {}",
            notation.rank(),
            simplifier.lattice().len()
        );

        // A spread of intervals, not just the unison, so the walk does real work.
        let intervals: Vec<Vec<i64>> = (0..divisions)
            .map(|ups| {
                let mut spelling = vec![0; notation.rank()];
                spelling[2] = ups;
                notation.to_just(&spelling).unwrap()
            })
            .collect();
        let coordinates: Vec<Vec<i64>> = intervals
            .iter()
            .map(|interval| notation.to_notation(interval).unwrap())
            .collect();

        let mut i = 0;
        time("Notation::to_notation", 1_000_000, || {
            i = (i + 1) % intervals.len();
            std::hint::black_box(notation.to_notation(&intervals[i]).unwrap());
        });

        let mut i = 0;
        time("Notation::to_just", 1_000_000, || {
            i = (i + 1) % coordinates.len();
            std::hint::black_box(notation.to_just(&coordinates[i]).unwrap());
        });

        let mut i = 0;
        time("Simplifier::simplify", 10_000, || {
            i = (i + 1) % intervals.len();
            std::hint::black_box(simplifier.simplify(&intervals[i]).unwrap());
        });

        let mut i = 0;
        time("Simplifier::candidates(8)", 10_000, || {
            i = (i + 1) % intervals.len();
            std::hint::black_box(simplifier.candidates(&intervals[i], 8).unwrap());
        });
    }
}
