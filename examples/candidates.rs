//! The candidate readings of an interval, for an end user to cycle through.
//!
//! One tempered pitch is any number of just intervals, and which one it "is"
//! depends on what is being played. The simplifier ranks them, so the answer it
//! gives is the first of a list rather than the only one there is.
//!
//! These four steps of 41et are the ones where the simplest reading is not the
//! most convenient spelling, so they are where having the rest of the list
//! matters: `14/11` and `11/7` are the ones worth reading, and cost three marks
//! apiece, while the runners up are spelled plainly and mean much less.

use xen_utils::{Notation, Simplifier, Subgroup, Temperament, simplify::sopfr};

const DIVISIONS: i64 = 41;
const STEPS: [i64; 4] = [14, 15, 20, 27];
const SHOWN: usize = 5;

fn main() {
    let subgroup: Subgroup = "2.3.5.7.11".parse().unwrap();
    let temperament = Temperament::et(DIVISIONS, &subgroup).unwrap();
    let notation = Notation::from_temperament(&temperament).unwrap();
    let simplifier = Simplifier::new(&notation).unwrap();

    for steps in STEPS {
        // The stack of ups that reaches this step, which is the awkward
        // spelling the simplifier is being asked to improve on.
        let mut spelling = vec![0; notation.rank()];
        spelling[2] = steps;
        let stacked = notation.to_just(&spelling).unwrap();

        let tempered = 1200.0 * steps as f64 / DIVISIONS as f64;
        let searched = simplifier.candidates(&stacked, usize::MAX).unwrap().len();
        println!("{steps} steps of {DIVISIONS}et, {tempered:.1}c - {searched} intervals searched");

        for candidate in simplifier.candidates(&stacked, SHOWN).unwrap() {
            let (numerator, denominator) = subgroup.to_ratio(&candidate).unwrap();
            println!(
                "    {:>7}  {:7} norm {:3}  {:+6.1}c",
                format!("{numerator}/{denominator}"),
                notation.note(&notation.to_notation(&candidate).unwrap()),
                sopfr(&candidate, subgroup.basis()),
                subgroup.to_cents(&candidate) - tempered,
            );
        }
        println!();
    }
}
