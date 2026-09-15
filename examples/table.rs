//! Spells the whole interval table of an equal temperament, one line per pitch,
//! under every notation in the run.
//!
//! ```text
//! cargo run --example table -- 41 2.3.5.7.11
//! ```
//!
//! This exists to judge a notation the way a reader meets it: not by how it
//! spells a prime, but by how it spells the interval the simplifier actually
//! hands back for each pitch. Each spelling is followed by the marks and the
//! sharps or flats it costs, and each notation by the totals, so that two
//! derivations can be compared over the whole table rather than at one prime.
//!
//! It also prints the enharmonic lattice in notation coordinates. Two notations
//! keeping the same generators have the same one - the map from notation
//! coordinates to pitches depends only on where the generators go - so any two
//! spellings of a pitch differ by an element of it, whichever notation produced
//! them.

use xen_utils::{Notation, Simplifier, Subgroup, Temperament};

fn main() {
    let mut args = std::env::args().skip(1);
    let divisions: i64 = args.next().and_then(|a| a.parse().ok()).unwrap_or(41);
    let subgroup: Subgroup = args
        .next()
        .unwrap_or_else(|| "2.3.5.7.11".to_string())
        .parse()
        .expect("a subgroup such as 2.3.5.7.11");

    let temperament = Temperament::et(divisions, &subgroup).expect("an equal temperament");
    let options = Notation::options(&temperament).expect("a notation");
    let recommended = Notation::from_temperament(&temperament).expect("the run is not empty");
    let simplifier = Simplifier::new(&options[0]).expect("a simplifier");

    println!("{divisions}et over {subgroup}\n");

    for n in &options {
        let mark = if n.mapping() == recommended.mapping() {
            " (recommended)"
        } else {
            ""
        };
        println!(
            "[{}] keeps {}{mark}",
            n.rank(),
            if n.rank() == 2 {
                "nothing".to_string()
            } else {
                n.generators()[2..]
                    .iter()
                    .map(|a| ratio(&subgroup, a))
                    .collect::<Vec<_>>()
                    .join(" ")
            }
        );
        println!("    enharmonics, in notation coordinates:");
        for e in n.enharmonics().expect("an enharmonic lattice") {
            println!("        {:?}", n.to_notation(&e).expect("an enharmonic"));
        }
    }

    let header: Vec<String> = options
        .iter()
        .map(|n| format!("{:14}", format!("[{}]", n.rank())))
        .collect();
    println!("\n{:>4}  {:>9}  {}", "step", "simplest", header.join(""));

    let mut totals = vec![(0i64, 0i64); options.len()];
    for step in 0..=divisions {
        let target = scaled(&subgroup, &temperament, step, divisions);
        let simplified = simplifier.simplify(&target).expect("a simplification");

        let cells: Vec<String> = options
            .iter()
            .enumerate()
            .map(|(index, n)| {
                let coordinates = n.to_notation(&simplified).expect("a just interval");
                let marks: i64 = coordinates[2..].iter().map(|c| c.abs()).sum();
                let sharps = (coordinates[1] + 1).div_euclid(7).abs();
                totals[index].0 += marks;
                totals[index].1 += sharps;
                format!(
                    "{:14}",
                    format!("{} {marks}/{sharps}", n.note(&coordinates))
                )
            })
            .collect();
        println!(
            "{step:4}  {:>9}  {}",
            ratio(&subgroup, &simplified),
            cells.join("")
        );
    }

    let summary: Vec<String> = totals
        .iter()
        .map(|(marks, sharps)| format!("{:14}", format!("{marks} marks {sharps} sh")))
        .collect();
    println!("{:4}  {:>9}  {}", "", "totals", summary.join(""));
}

/// Some just interval worth `step` steps, for the simplifier to reduce. Any one
/// will do, since simplifying is a property of the temperament and every
/// interval worth the same reduces to the same place, so this takes the first
/// the fifth chain reaches.
fn scaled(subgroup: &Subgroup, t: &Temperament, step: i64, divisions: i64) -> Vec<i64> {
    // A generator of one step: the temperament maps the octave to `divisions`,
    // so any interval mapping to 1 will do, and the fifth chain reaches one
    // wherever a notation exists at all.
    let mut octave = vec![0i64; subgroup.dim()];
    octave[0] = 1;
    let mut fifth = vec![0i64; subgroup.dim()];
    fifth[0] = -1;
    fifth[1] = 1;
    let per_fifth = t.map(&fifth).expect("a fifth")[0];

    // Solve `a * divisions + b * per_fifth = step` for small `b`.
    for b in -divisions..=divisions {
        let remainder = step - b * per_fifth;
        if remainder.rem_euclid(divisions) == 0 {
            let a = remainder / divisions;
            return (0..subgroup.dim())
                .map(|i| a * octave[i] + b * fifth[i])
                .collect();
        }
    }
    unreachable!("the fifth chain reaches every step of a notatable equal temperament")
}

fn ratio(subgroup: &Subgroup, interval: &[i64]) -> String {
    match subgroup.to_ratio(interval) {
        Ok((num, den)) => format!("{num}/{den}"),
        Err(_) => format!("{interval:?}"),
    }
}
