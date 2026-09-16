//! The candidate **readings** of a pitch, for an end user to cycle through.
//!
//! One tempered pitch is any number of just intervals, and which one it "is"
//! depends on what is being played. The simplifier ranks them, so the answer it
//! gives is the first of a list rather than the only one there is.
//!
//! This is one of the two questions a pitch has; `cargo run --example spellings`
//! is the other. They are separate, and the spelling does not appear here for a
//! reason: every reading on this list is the same pitch, so the notation writes
//! them all the same way. Which just interval is meant and how the note is
//! written are decided by different things and neither constrains the other.

use xen_utils::{Notation, Simplifier, Subgroup, Temperament, simplify::sopfr};

const DIVISIONS: i64 = 41;
const STEPS: [i64; 4] = [14, 15, 20, 27];
const SHOWN: usize = 5;

fn main() {
    let subgroup: Subgroup = "2.3.5.7.11".parse().unwrap();
    let temperament = Temperament::equal(DIVISIONS, &subgroup).unwrap();
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
        let written = notation.note(&notation.spell(&stacked).unwrap());
        println!(
            "{steps} steps of {DIVISIONS}et, {tempered:.1}c, written {written}, {searched} searched"
        );

        for candidate in simplifier.candidates(&stacked, SHOWN).unwrap() {
            let (numerator, denominator) = subgroup.to_ratio(&candidate).unwrap();
            println!(
                "    {:>7}  norm {:3}  {:+6.1}c",
                format!("{numerator}/{denominator}"),
                sopfr(&candidate, subgroup.basis()),
                subgroup.to_cents(&candidate) - tempered,
            );
        }
        println!();
    }
}
