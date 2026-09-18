//! How wide `spellings` has to walk: its answers at widths 0, 1 and 2 against
//! a much wider walk, over every notation of every temperament in the data
//! files and every equal temperament up to 72 over a few subgroups.

use std::error::Error;

use xen_utils::{Notation, Subgroup, Temperament};

const TEMPERAMENTS: &str = include_str!("../data/temperaments.txt");
const TEMPERAMENTS_BIG: &str = include_str!("../data/temperaments_big.txt");

const WIDTHS: [i64; 3] = [0, 1, 2];
const COUNTS: [usize; 3] = [1, 3, 5];
/// Examples printed per width and count.
const SHOWN: usize = 4;

#[derive(Default)]
struct Tally {
    calls: usize,
    /// Per width, per count: how many calls differ from the reference.
    differ: [[usize; COUNTS.len()]; WIDTHS.len()],
    /// Calls where the reference width was too small to compare that width.
    skipped: [usize; WIDTHS.len()],
    shown: [[usize; COUNTS.len()]; WIDTHS.len()],
    /// Calls where the reference itself moved when widened by one.
    reference_moved: usize,
}

fn main() {
    let mut tally = Tally::default();
    let mut seen = std::collections::HashSet::new();

    for data in [TEMPERAMENTS, TEMPERAMENTS_BIG] {
        for line in data.lines() {
            let line = line.split('#').next().unwrap_or_default().trim();
            let mut fields = line.split('|').map(str::trim);
            let (Some(name), Some(subgroup), Some(definition)) =
                (fields.next(), fields.next(), fields.next())
            else {
                continue;
            };
            if !seen.insert((subgroup.to_string(), definition.to_string())) {
                continue;
            }
            match parse(subgroup, definition) {
                Ok(t) => sweep(&format!("{name} over {subgroup}"), &t, &mut tally),
                Err(error) => eprintln!("{name}: {error}"),
            }
        }
    }
    for subgroup in ["2.3.5", "2.3.5.7", "2.3.5.7.11", "2.3.5.7.11.13"] {
        let subgroup: Subgroup = subgroup.parse().unwrap();
        for divisions in 5..=72 {
            let Ok(t) = Temperament::equal(divisions, &subgroup) else {
                continue;
            };
            sweep(&format!("{divisions}et over {subgroup}"), &t, &mut tally);
        }
    }

    println!("\n{} spellings calls", tally.calls);
    println!("reference moved when widened: {}", tally.reference_moved);
    for (w, width) in WIDTHS.iter().enumerate() {
        let counts: Vec<String> = COUNTS
            .iter()
            .enumerate()
            .map(|(c, count)| format!("top {count}: {}", tally.differ[w][c]))
            .collect();
        println!(
            "width {width}: differs  {}   (not compared: {})",
            counts.join("  "),
            tally.skipped[w]
        );
    }
}

/// The reference width for a lattice of `enharmonics` members: as wide as
/// is affordable.
fn reference_width(enharmonics: usize) -> i64 {
    match enharmonics {
        0..=2 => 8,
        3 => 5,
        4 => 4,
        _ => 3,
    }
}

fn sweep(name: &str, t: &Temperament, tally: &mut Tally) {
    let Ok(options) = Notation::options(t) else {
        return;
    };
    let max_count = *COUNTS.iter().max().unwrap();

    for (index, n) in options.iter().enumerate() {
        let lattice = n.enharmonics();
        if lattice.is_empty() {
            continue;
        }
        let reference = reference_width(lattice.len());

        for interval in intervals(n.dim()) {
            let tempered = t.temper(&interval).unwrap();
            tally.calls += 1;
            let wanted = n.spellings_within(&tempered, max_count, reference).unwrap();

            // Is the reference itself settled? Only where it is cheap to check.
            if lattice.len() <= 3 {
                let wider = n
                    .spellings_within(&tempered, max_count, reference + 1)
                    .unwrap();
                if wider != wanted {
                    tally.reference_moved += 1;
                }
            }

            for (w, &width) in WIDTHS.iter().enumerate() {
                if width >= reference {
                    tally.skipped[w] += 1;
                    continue;
                }
                let got = n.spellings_within(&tempered, max_count, width).unwrap();
                for (c, &count) in COUNTS.iter().enumerate() {
                    let upto = count.min(wanted.len());
                    if got.len() < upto || got[..upto] != wanted[..upto] {
                        tally.differ[w][c] += 1;
                        if tally.shown[w][c] < SHOWN {
                            tally.shown[w][c] += 1;
                            let show = |list: &[Vec<i64>]| {
                                list.iter()
                                    .take(count)
                                    .map(|c| n.note(c))
                                    .collect::<Vec<_>>()
                                    .join(" ")
                            };
                            println!(
                                "width {width} top {count}: {name} [{}] option {index}, {tempered:?}: got [{}] wanted [{}]",
                                n.len(),
                                show(&got),
                                show(&wanted)
                            );
                        }
                    }
                }
            }
        }
    }
}

/// Every interval with exponents in -2..=2 on the primes beyond 2, or -1..=1
/// once there are more than four of them, octave-free.
fn intervals(dim: usize) -> Vec<Vec<i64>> {
    let bound = if dim > 5 { 1 } else { 2 };
    let mut result = Vec::new();
    let mut exponents = vec![-bound; dim - 1];
    loop {
        let mut interval = vec![0];
        interval.extend(&exponents);
        result.push(interval);
        let mut place = 0;
        while place < exponents.len() && exponents[place] == bound {
            exponents[place] = -bound;
            place += 1;
        }
        if place == exponents.len() {
            return result;
        }
        exponents[place] += 1;
    }
}

fn parse(subgroup: &str, definition: &str) -> Result<Temperament, Box<dyn Error>> {
    let subgroup: Subgroup = subgroup.parse()?;
    if let Some(divisions) = definition.strip_prefix("et ") {
        return Ok(Temperament::equal(divisions.trim().parse()?, &subgroup)?);
    }
    let commas = definition
        .split([',', ' '])
        .filter(|field| !field.is_empty())
        .map(|ratio| {
            let (num, den) = ratio.split_once('/').ok_or("not a ratio")?;
            Ok(subgroup.factorize(num.parse()?, den.parse()?)?)
        })
        .collect::<Result<Vec<Vec<i64>>, Box<dyn Error>>>()?;
    Ok(Temperament::from_commas(&commas, &subgroup)?)
}
