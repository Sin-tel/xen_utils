//! Walks one temperament's derivation a step at a time, laying out the choice
//! at each point and marking the one taken.
//!
//! ```text
//! cargo run --example trace              # 72et, the default
//! cargo run --example trace -- miracle   # any name in data/temperaments_big.txt
//! ```
//!
//! The question this is built to answer is what the options are when a dropped
//! accidental is replaced by a stack of the ones still kept, and how that
//! choice reaches the spellings. So each stage prints every stack that works,
//! not only the winner.
//!
//! Nothing here re-implements the decision. The stacks are enumerated - which
//! is just solving a linear system and adding the kernel - and the winner is
//! read back out of `Notation::options`, so this cannot drift from what the
//! library does. It asserts that it found the library's answer among the
//! options it listed.

use std::error::Error;

use diophantine::{Matrix, kernel_right, lll, transpose};
use xen_utils::{Notation, Subgroup, Temperament};

const TEMPERAMENTS: &str = include_str!("../data/temperaments_big.txt");
const DEFAULT: &str = "72et";

/// How far either way to walk each relation when listing the stacks that work.
const WIDTH: i64 = 5;

fn main() {
    let wanted = std::env::args()
        .nth(1)
        .unwrap_or_else(|| DEFAULT.to_string());
    let mut found = false;

    for line in TEMPERAMENTS.lines() {
        let line = line.split('#').next().unwrap_or_default().trim();
        if line.is_empty() {
            continue;
        }
        let mut fields = line.split('|').map(str::trim);
        let (Some(name), Some(subgroup), Some(definition)) =
            (fields.next(), fields.next(), fields.next())
        else {
            continue;
        };
        if name != wanted {
            continue;
        }
        found = true;
        match parse(subgroup, definition) {
            Ok(t) => trace(name, &t),
            Err(error) => eprintln!("{name}: {error}"),
        }
    }

    if !found {
        eprintln!("no temperament named {wanted} in data/temperaments_big.txt");
    }
}

fn trace(name: &str, t: &Temperament) {
    let subgroup = t.subgroup();
    println!("== {name} over {subgroup}, rank {} ==\n", t.rank());

    let ji = match Notation::from_ji(subgroup) {
        Ok(ji) => ji,
        Err(error) => return println!("no just intonation notation: {error}"),
    };
    let accidentals: Matrix<i64> = ji.generators()[2..].to_vec();

    println!("the accidentals, and what each is worth to the temperament");
    for a in &accidentals {
        println!(
            "    {:8}  {}",
            ratio(subgroup, a),
            worth(t, &t.map(a).expect("an accidental of this subgroup"))
        );
    }

    let options = match Notation::options(t) {
        Ok(options) => options,
        Err(error) => return println!("\nno notation: {error}"),
    };

    println!("\nthe run");
    for n in &options {
        println!(
            "    [{}]  keeps {}",
            n.rank(),
            names(subgroup, &n.generators()[2..].to_vec())
        );
    }

    // Each notation in the run keeps one accidental more than the one before,
    // and the comma of that accidental is the one that drops out of the kernel
    // when it is taken on. So the stages are the consecutive pairs, plus a last
    // one for anything no notation ever keeps.
    for pair in 0..options.len() {
        let kept = &options[pair].generators()[2..].to_vec();
        let next: Matrix<i64> = match options.get(pair + 1) {
            Some(n) => n.generators()[2..].to_vec(),
            None => Vec::new(),
        };
        let taken: Vec<Vec<i64>> = next.iter().filter(|a| !kept.contains(a)).cloned().collect();
        let never: Vec<Vec<i64>> = accidentals
            .iter()
            .filter(|a| !kept.contains(a) && !next.contains(a))
            .cloned()
            .collect();

        // The accidental this stage is about: the one the next notation takes
        // on, or - at the end of the run - everything no notation ever keeps.
        let dropped = if taken.is_empty() { never } else { taken };
        for accidental in &dropped {
            stage(t, &ji, &options[pair], kept, accidental);
        }
    }

    spellings(t, &options);
}

/// One dropped accidental: what it has to be worth, what the available
/// accidentals can make of that, and which replacement was taken.
///
/// The options are listed as a box around the one taken, since any two
/// replacements worth the same differ by a relation among the generators. That
/// is the point of the listing - these are the neighbours the answer beat, not
/// a re-derivation of how it was found.
fn stage(t: &Temperament, ji: &Notation, notation: &Notation, kept: &Matrix<i64>, dropped: &[i64]) {
    let subgroup = t.subgroup();
    println!(
        "
replacing {} where the notation keeps {}",
        ratio(subgroup, dropped),
        names(subgroup, kept)
    );

    let target = t.map(dropped).expect("an accidental of this subgroup");
    println!("    it is worth {}", worth(t, &target));

    // The comma the library settled on: the one that leaves the kernel when
    // this accidental is taken on.
    let index = index_of(ji, dropped).expect("an accidental of this notation");
    let comma = notation
        .commas()
        .iter()
        .find(|c| {
            let coordinates = ji.to_notation(c).expect("a comma of this subgroup");
            coordinates[2 + index] != 0
                && ji.generators()[2..].iter().enumerate().all(|(other, a)| {
                    other == index || kept.contains(a) || coordinates[2 + other] == 0
                })
        })
        .expect("a comma for every dropped accidental");
    let chosen = ji.to_notation(comma).expect("a comma of this subgroup");
    println!(
        "    the comma taken is {:12} {:?} in notation coordinates",
        ratio(subgroup, &subgroup.ascending(comma)),
        chosen
    );

    // Step 1 is a stack of the kept accidentals alone, which in the notation
    // basis is exactly a comma with no octave and no fifth to it. Step 2 is
    // everything else, and it adds the fifth chain to the generators.
    let step_two = chosen[0] != 0 || chosen[1] != 0;
    let mut generators = kept.clone();
    if step_two {
        println!("    no stack of them is worth that, so the fifth chain is walked too");
        generators.extend(fifth_chain(subgroup.dim()));
    }

    // The count of each generator is its notation coordinate, up to the sign
    // the comma happens to be written with.
    let sign = if chosen[2 + index] > 0 { -1 } else { 1 };
    let mut taken: Vec<i64> = kept
        .iter()
        .map(|a| sign * chosen[2 + index_of(ji, a).expect("a kept accidental")])
        .collect();
    if step_two {
        taken.push(sign * chosen[0]);
        taken.push(sign * chosen[1]);
    }

    let images = t.map_all(&generators).expect("generators of this subgroup");
    let relations = relations(&images);
    if relations.is_empty() {
        println!("    they have no relation among them, so the replacement is forced");
    } else if !step_two {
        println!("    they do have relations, so there is a choice:");
        for relation in &relations {
            println!("        {}", relation_line(subgroup, &generators, relation));
        }
    }

    let mut options: Vec<(i64, i64, i64, Vec<i64>)> = Vec::new();
    let mut steps = vec![-WIDTH; relations.len()];
    loop {
        let counts: Vec<i64> = combination(&steps, &relations, taken.len())
            .iter()
            .zip(&taken)
            .map(|(a, b)| a + b)
            .collect();
        let comma = difference(dropped, &counts, &generators);
        assert!(
            is_zero(&t.map(&comma).expect("an interval of this subgroup")),
            "a listed replacement is not worth what the accidental is worth"
        );
        let marks = counts[..kept.len()].iter().map(|c| c.abs()).sum();
        let fifths = if step_two { counts[kept.len() + 1] } else { 0 };
        options.push((sopfr(subgroup, &comma), marks, fifths, counts));

        let mut place = 0;
        while place < steps.len() && steps[place] == WIDTH {
            steps[place] = -WIDTH;
            place += 1;
        }
        if place == steps.len() {
            break;
        }
        steps[place] += 1;
    }

    // Each list is sorted the way the code that picks from it sorts: step 1 by
    // the marks it costs, step 2 by the norm of the comma.
    if step_two {
        options.sort_unstable_by_key(|(sopfr, marks, _, counts)| (*sopfr, *marks, counts.clone()));
        println!("    the replacements that work, simplest comma first");
    } else {
        options.sort_unstable_by_key(|(sopfr, marks, _, counts)| (*marks, counts.clone(), *sopfr));
        println!("    the stacks that work, fewest marks first");
    }
    println!(
        "        {:38}  {:>5}  {:>6}  {:>5}",
        "", "marks", "fifths", "sopfr"
    );

    for (place, (sopfr, marks, fifths, counts)) in options.iter().enumerate() {
        let is_taken = *counts == taken;
        if place >= 8 && !is_taken {
            continue;
        }
        let description = match stack_line(subgroup, kept, &counts[..kept.len()]) {
            line if line == "nothing" && step_two => "the fifth chain alone".to_string(),
            line => line,
        };
        println!(
            "        {description:38}  {marks:>5}  {fifths:>6}  {sopfr:>5}{}",
            if is_taken { "  <- taken" } else { "" }
        );
    }
}

/// The octave and the fifth, as interval vectors.
fn fifth_chain(dim: usize) -> Matrix<i64> {
    let mut octave = vec![0; dim];
    octave[0] = 1;
    let mut fifth = vec![0; dim];
    fifth[0] = -1;
    fifth[1] = 1;
    vec![octave, fifth]
}

/// `vector` less the combination `counts` makes of `generators`.
fn difference(vector: &[i64], counts: &[i64], generators: &Matrix<i64>) -> Vec<i64> {
    let made = combination(counts, generators, vector.len());
    vector.iter().zip(made).map(|(a, b)| a - b).collect()
}

/// The sum of the prime factors of numerator and denominator, with
/// multiplicity - the Wilson norm, and an `L1` norm on the prime exponents.
fn sopfr(subgroup: &Subgroup, interval: &[i64]) -> i64 {
    interval
        .iter()
        .zip(subgroup.basis())
        .map(|(e, &p)| e.abs() * i64::from(p))
        .sum()
}

/// How the choice reaches the page: each prime, spelled in every notation of
/// the run, smallest notation last so that the substitutions read downwards.
fn spellings(t: &Temperament, options: &[Notation]) {
    let subgroup = t.subgroup();
    println!("\nhow that reaches the spelling");
    for index in 2..subgroup.dim() {
        let mut harmonic = vec![0i64; subgroup.dim()];
        harmonic[index] = 1;
        harmonic[0] = -(subgroup.to_cents(&harmonic) / 1200.0).floor() as i64;
        let line: Vec<String> = options
            .iter()
            .rev()
            .map(|n| {
                let coordinates = n.to_notation(&harmonic).expect("a just interval");
                format!("[{}] {:8}", n.rank(), n.note(&coordinates))
            })
            .collect();
        println!("    {:8}  {}", ratio(subgroup, &harmonic), line.join("  "));
    }
    println!(
        "\n    a smaller notation substitutes its own commas into a larger one's\n    \
         spelling, which is why the marks compound downwards."
    );
}

/// The relations among `images`: the stacks of them worth nothing, which is
/// what any two stacks worth the same differ by.
///
/// Reduced, because the box walked around them is only meaningful if they are
/// short - `simplest_comma` reduces for the same reason. Without this, orwell's
/// relation comes back as 31 syntonic commas against 7 octaves and a box of
/// four misses `225/224` entirely.
fn relations(images: &Matrix<i64>) -> Matrix<i64> {
    if images.is_empty() {
        return Vec::new();
    }
    let kernel: Matrix<i64> = transpose(&kernel_right(&transpose(images)).unwrap_or_default())
        .into_iter()
        .filter(|r| !is_zero(r))
        .collect();
    if kernel.is_empty() {
        return kernel;
    }
    let weights = unit_weights(kernel[0].len());
    lll(&kernel, 0.99, &weights).unwrap_or(kernel)
}

/// An identity weight matrix. These coordinates count generators rather than
/// primes, so there is nothing to weight them by.
fn unit_weights(size: usize) -> Matrix<f64> {
    (0..size)
        .map(|row| {
            (0..size)
                .map(|col| f64::from(u8::from(row == col)))
                .collect()
        })
        .collect()
}

/// A stack written out, such as `81/80 + 64/63` or `3 x 81/80`.
fn stack_line(subgroup: &Subgroup, kept: &Matrix<i64>, counts: &[i64]) -> String {
    let terms: Vec<String> = counts
        .iter()
        .zip(kept)
        .filter(|&(&count, _)| count != 0)
        .map(|(&count, a)| {
            let name = ratio(subgroup, a);
            match count {
                1 => format!("+ {name}"),
                -1 => format!("- {name}"),
                c if c > 0 => format!("+ {c} x {name}"),
                c => format!("- {} x {name}", -c),
            }
        })
        .collect();
    if terms.is_empty() {
        return "nothing".to_string();
    }
    terms.join(" ").trim_start_matches("+ ").to_string()
}

/// A relation written as an equation, such as `2 x 81/80 = 64/63`.
fn relation_line(subgroup: &Subgroup, kept: &Matrix<i64>, relation: &[i64]) -> String {
    let left: Vec<i64> = relation.iter().map(|&c| c.max(0)).collect();
    let right: Vec<i64> = relation.iter().map(|&c| (-c).max(0)).collect();
    format!(
        "{} = {}",
        stack_line(subgroup, kept, &left),
        stack_line(subgroup, kept, &right)
    )
}

/// What an image is worth, in steps for an equal temperament and as the vector
/// itself for anything else.
fn worth(t: &Temperament, image: &[i64]) -> String {
    if is_zero(image) {
        return "nothing - the temperament tempers it out".to_string();
    }
    if t.rank() == 1 {
        let steps = image[0];
        return format!("{steps} step{}", if steps.abs() == 1 { "" } else { "s" });
    }
    format!("{image:?}")
}

/// The index of `accidental` among the notation's accidentals.
fn index_of(ji: &Notation, accidental: &[i64]) -> Option<usize> {
    ji.generators()[2..].iter().position(|a| a == accidental)
}

fn names(subgroup: &Subgroup, accidentals: &Matrix<i64>) -> String {
    if accidentals.is_empty() {
        return "nothing".to_string();
    }
    accidentals
        .iter()
        .map(|a| ratio(subgroup, a))
        .collect::<Vec<_>>()
        .join(" ")
}

fn ratio(subgroup: &Subgroup, interval: &[i64]) -> String {
    match subgroup.to_ratio(interval) {
        Ok((num, den)) => format!("{num}/{den}"),
        Err(_) => format!("{interval:?}"),
    }
}

fn is_zero(vector: &[i64]) -> bool {
    vector.iter().all(|&e| e == 0)
}

fn combination(counts: &[i64], generators: &Matrix<i64>, dim: usize) -> Vec<i64> {
    let mut total = vec![0; dim];
    for (&count, generator) in counts.iter().zip(generators) {
        for (slot, &e) in total.iter_mut().zip(generator) {
            *slot += count * e;
        }
    }
    total
}

fn parse(subgroup: &str, definition: &str) -> Result<Temperament, Box<dyn Error>> {
    let subgroup: Subgroup = subgroup.parse()?;
    if let Some(divisions) = definition.strip_prefix("et ") {
        return Ok(Temperament::et(divisions.trim().parse()?, &subgroup)?);
    }
    let commas = definition
        .split([',', ' '])
        .filter(|field| !field.is_empty())
        .map(|r| {
            let (num, den) = r.split_once('/').ok_or("not a ratio")?;
            Ok(subgroup.factorize(num.parse()?, den.parse()?)?)
        })
        .collect::<Result<Vec<Vec<i64>>, Box<dyn Error>>>()?;
    Ok(Temperament::from_commas(&commas, &subgroup)?)
}
