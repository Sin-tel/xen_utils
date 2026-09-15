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

use diophantine::{Matrix, hnf, lll};
use xen_utils::{Notation, Simplifier, Subgroup, Temperament, Weighting};

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
        println!("    kernel - the commas it writes away, so the pairs of just");
        println!("    intervals it cannot tell apart:");
        for comma in reduced(&subgroup, n.commas()) {
            println!("        {:14} {:?}", ratio(&subgroup, &comma), comma);
        }
        println!("    the same kernel in hermite normal form:");
        for row in normal(n.commas()) {
            println!("        {row:?}");
        }
        println!("    enharmonics, in notation coordinates, hermite normal form:");
        for row in normal(&enharmonics(n)) {
            println!("        {row:?}");
        }
    }

    let header: Vec<String> = options
        .iter()
        .map(|n| format!("{:14}", format!("[{}]", n.rank())))
        .collect();
    println!("\n{:>4}  {:>9}  {}", "step", "simplest", header.join(""));

    let mut totals = vec![(0i64, 0i64); options.len()];
    let mut coset = (0i64, 0i64);
    let mut symbols = (0i64, 0i64);
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
        // The same pitch, spelled by searching the coset in the recommended
        // notation rather than by mapping the interval.
        let seed = recommended
            .to_notation(&simplified)
            .expect("a just interval");
        let lattice = enharmonics(&recommended);
        let by_marks = best(&seed, &lattice, |marks, sharps| (marks, sharps));
        let by_symbols = best(&seed, &lattice, |marks, sharps| (marks + sharps, marks));
        coset.0 += cost(&by_marks).0;
        coset.1 += cost(&by_marks).1;
        symbols.0 += cost(&by_symbols).0;
        symbols.1 += cost(&by_symbols).1;

        println!(
            "{step:4}  {:>9}  {}{:14}{:14}",
            ratio(&subgroup, &simplified),
            cells.join(""),
            format!(
                "{} {}/{}",
                recommended.note(&by_marks),
                cost(&by_marks).0,
                cost(&by_marks).1
            ),
            format!(
                "{} {}/{}",
                recommended.note(&by_symbols),
                cost(&by_symbols).0,
                cost(&by_symbols).1
            )
        );
    }

    let summary: Vec<String> = totals
        .iter()
        .map(|(marks, sharps)| format!("{:14}", format!("{marks} marks {sharps} sh")))
        .collect();
    println!(
        "{:4}  {:>9}  {}{:14}{:14}",
        "",
        "totals",
        summary.join(""),
        format!("{} marks {} sh", coset.0, coset.1),
        format!("{} marks {} sh", symbols.0, symbols.1)
    );
}

/// What a spelling costs to write: its accidental marks, and its sharps or
/// flats. Seven fifths are a sharp, and a fifth coordinate from -1 to 5 is the
/// range that needs none.
fn cost(coordinates: &[i64]) -> (i64, i64) {
    (
        coordinates[2..].iter().map(|c| c.abs()).sum(),
        (coordinates[1] + 1).div_euclid(7).abs(),
    )
}

/// The cheapest spelling of the pitch `seed` spells, searched over its coset of
/// the enharmonic lattice. Ties are settled by the coordinates so that the
/// answer does not depend on the order the box is walked in.
fn best<K: Ord>(seed: &[i64], lattice: &Matrix<i64>, key: impl Fn(i64, i64) -> K) -> Vec<i64> {
    const WIDTH: i64 = 6;
    let rank = lattice.len();
    let mut chosen = seed.to_vec();
    let mut steps = vec![-WIDTH; rank];
    loop {
        let candidate: Vec<i64> = (0..seed.len())
            .map(|slot| {
                seed[slot]
                    + steps
                        .iter()
                        .zip(lattice)
                        .map(|(&step, row)| step * row[slot])
                        .sum::<i64>()
            })
            .collect();
        let (marks, sharps) = cost(&candidate);
        let (best_marks, best_sharps) = cost(&chosen);
        if (key(marks, sharps), candidate.clone()) < (key(best_marks, best_sharps), chosen.clone())
        {
            chosen = candidate;
        }
        let mut place = 0;
        while place < rank && steps[place] == WIDTH {
            steps[place] = -WIDTH;
            place += 1;
        }
        if place == rank {
            return chosen;
        }
        steps[place] += 1;
    }
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

/// A kernel basis reduced to small ascending intervals, which is the readable
/// form: the stored basis is one comma per dropped accidental and says nothing.
fn reduced(subgroup: &Subgroup, commas: &Matrix<i64>) -> Matrix<i64> {
    if commas.is_empty() {
        return Vec::new();
    }
    let weights = subgroup.weights(Weighting::Wilson);
    lll(commas, 0.99, &weights)
        .expect("a kernel basis")
        .iter()
        .map(|c| subgroup.ascending(c))
        .collect()
}

/// A lattice in hermite normal form, which is the only way to tell two of them
/// apart: a basis is not canonical and two different ones can span the same
/// lattice.
fn normal(lattice: &Matrix<i64>) -> Matrix<i64> {
    if lattice.is_empty() {
        return Vec::new();
    }
    hnf(lattice).expect("a lattice basis")
}

/// The enharmonic lattice in notation coordinates.
fn enharmonics(n: &Notation) -> Matrix<i64> {
    n.enharmonics()
        .expect("an enharmonic lattice")
        .iter()
        .map(|e| n.to_notation(e).expect("an enharmonic"))
        .collect()
}

fn ratio(subgroup: &Subgroup, interval: &[i64]) -> String {
    match subgroup.to_ratio(interval) {
        Ok((num, den)) => format!("{num}/{den}"),
        Err(_) => format!("{interval:?}"),
    }
}
