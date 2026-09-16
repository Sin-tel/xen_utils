//! The other ways to write a note: given a spelling, everything else the
//! temperament calls the same pitch.
//!
//! ```text
//! cargo run --example respell -- 41 2.3.5.7.11
//! ```
//!
//! The spellings of one pitch are a coset of the enharmonic lattice, and that
//! lattice is `kernel_left` of the generator images - so it depends on the
//! generators and the temperament and on nothing else. The notation's own
//! kernel never enters, which is why this question can be asked of the symbol
//! system alone and why two notations built on the same generators answer it
//! identically.
//!
//! Ranking spellings is not the dead end that ranking intervals by notation
//! coordinates would be. There is no prime on the fifth axis, so no norm there
//! means anything about pitch - but the number of symbols on the page is not a
//! statement about pitch. It is a count, and counting is exactly what is wanted
//! here.
//!
//! Which count, though, matters. Two are printed. **Marks before sharps** reads
//! the accidentals first and breaks ties on sharps; **symbols** counts both
//! alike. Marks before sharps runs away: 41et's fifth chain reaches every pitch,
//! so somewhere out along it there is always a spelling with no marks at all and
//! six sharps, and preferring no marks at any price finds it. The run ends by
//! walking a wider box and reporting how much of each ordering survived, which
//! is the test that says so.

use diophantine::{Matrix, lll};
use xen_utils::{Notation, Simplifier, Subgroup, Temperament};

/// How far either way to walk each enharmonic. The answers wanted are the first
/// few, so this only has to be wide enough that nothing better lies outside.
const WIDTH: i64 = 4;

/// How many alternatives to show.
const SHOWN: usize = 6;

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
    let notation = Notation::from_temperament(&temperament).expect("the run is not empty");
    let simplifier = Simplifier::new(&options[0]).expect("a simplifier");

    println!(
        "{divisions}et over {subgroup}, notation of rank {} keeping {}",
        notation.rank(),
        notation.generators()[2..]
            .iter()
            .map(|a| ratio(&subgroup, a))
            .collect::<Vec<_>>()
            .join(" ")
    );
    let lattice = enharmonics(&notation);
    println!("enharmonics, in notation coordinates:");
    for row in &lattice {
        println!("    {row:?}");
    }
    println!(
        "\n{:>4}  {:>9}  {:>12}  ways to write it, simplest first",
        "step", "simplest", "ranked by"
    );

    let mut stable = [(0, 0), (0, 0)];
    for step in 0..divisions {
        let target = reached(&subgroup, &temperament, step, divisions);
        let simplified = simplifier.simplify(&target).expect("a simplification");
        let seed = notation.to_notation(&simplified).expect("a just interval");

        for (which, label) in ["marks", "symbols"].into_iter().enumerate() {
            // The same ranking over a wider box. If the two disagree near the
            // top, the ordering is not converging and the count is wrong.
            let near = ranked(&seed, &lattice, label, WIDTH);
            let far = ranked(&seed, &lattice, label, WIDTH + 3);
            stable[which].0 += usize::from(near[..2] == far[..2]);
            stable[which].1 += 1;

            let ways: Vec<String> = near
                .iter()
                .take(SHOWN)
                .map(|c| {
                    let (marks, sharps) = cost(c);
                    format!("{:11}", format!("{} {marks}/{sharps}", notation.note(c)))
                })
                .collect();
            println!(
                "{:>4}  {:>9}  {label:>12}  {}",
                if label == "marks" {
                    step.to_string()
                } else {
                    String::new()
                },
                if label == "marks" {
                    ratio(&subgroup, &simplified)
                } else {
                    String::new()
                },
                ways.join("")
            );
        }
    }

    println!();
    for (which, label) in ["marks", "symbols"].into_iter().enumerate() {
        let (same, total) = stable[which];
        println!("{label:>12}: the first two agree with a wider box on {same} of {total} pitches");
    }
}

/// The coset of `seed`, sorted by `label`.
///
/// "marks" reads the accidentals first and the sharps only as a tie-break;
/// "symbols" counts everything on the page alike. Ties fall back to the
/// coordinates, so the order never depends on how the box was walked.
fn ranked(seed: &[i64], lattice: &Matrix<i64>, label: &str, width: i64) -> Vec<Vec<i64>> {
    let mut coset = walk(seed, lattice, width);
    coset.sort_by_key(|c| {
        let (marks, sharps) = cost(c);
        let rank = if label == "marks" {
            (marks, sharps)
        } else {
            (marks + sharps, marks)
        };
        (rank, c.clone())
    });
    coset
}

/// What a spelling costs to write: its accidental marks, and its sharps or
/// flats. Seven fifths are a sharp, and a fifth coordinate from -1 to 5 needs
/// none.
fn cost(coordinates: &[i64]) -> (i64, i64) {
    (
        coordinates[2..].iter().map(|c| c.abs()).sum(),
        (coordinates[1] + 1).div_euclid(7).abs(),
    )
}

/// Every spelling within `width` enharmonics of `seed`.
fn walk(seed: &[i64], lattice: &Matrix<i64>, width: i64) -> Vec<Vec<i64>> {
    let mut found = Vec::new();
    let mut steps = vec![-width; lattice.len()];
    loop {
        found.push(
            (0..seed.len())
                .map(|slot| {
                    seed[slot]
                        + steps
                            .iter()
                            .zip(lattice)
                            .map(|(&step, row)| step * row[slot])
                            .sum::<i64>()
                })
                .collect(),
        );
        let mut place = 0;
        while place < lattice.len() && steps[place] == width {
            steps[place] = -width;
            place += 1;
        }
        if place == lattice.len() {
            return found;
        }
        steps[place] += 1;
    }
}

/// The enharmonic lattice in notation coordinates, reduced.
///
/// Reduced because the box is walked around it and a box around a long basis
/// reaches nothing useful. 41et's lattice comes back as "the octave is 41 ups"
/// and "the fifth is 24 ups", and four steps either way along those two never
/// reach `vB#`, which is one mark and one sharp, while happily reaching
/// `vvvvvvvD`, which is seven marks. Reduced, the basis is the limma and the
/// pythagorean comma and the near spellings are near.
fn enharmonics(n: &Notation) -> Matrix<i64> {
    let lattice: Matrix<i64> = n
        .enharmonics()
        .expect("an enharmonic lattice")
        .iter()
        .map(|e| n.to_notation(e).expect("an enharmonic"))
        .collect();
    if lattice.is_empty() {
        return lattice;
    }
    // Unit weights: these coordinates count symbols, not primes.
    let size = lattice[0].len();
    let weights = (0..size)
        .map(|row| {
            (0..size)
                .map(|col| f64::from(u8::from(row == col)))
                .collect()
        })
        .collect();
    lll(&lattice, 0.99, &weights).unwrap_or(lattice)
}

/// Some just interval worth `step` steps, for the simplifier to reduce.
fn reached(subgroup: &Subgroup, t: &Temperament, step: i64, divisions: i64) -> Vec<i64> {
    let mut fifth = vec![0i64; subgroup.dim()];
    fifth[0] = -1;
    fifth[1] = 1;
    let per_fifth = t.map(&fifth).expect("a fifth")[0];
    for count in -divisions..=divisions {
        let remainder = step - count * per_fifth;
        if remainder.rem_euclid(divisions) == 0 {
            let octaves = remainder / divisions;
            let mut interval = vec![0i64; subgroup.dim()];
            interval[0] = octaves - count;
            interval[1] = count;
            return interval;
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
