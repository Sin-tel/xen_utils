//! Prints a fingerprint of every notation the library derives, over the whole
//! temperament list and every equal temperament to 99 across several subgroups.
//!
//! This exists to be diffed. Redirect it to a file, change the derivation, run
//! it again, and the diff is exactly the set of temperaments whose spelling
//! moved - which is the only way to tell whether a simpler rule is the same
//! rule. The fingerprint is the mapping, the accidentals kept, and how each
//! prime beyond 3 comes out, since those are what a caller can observe.
//!
//! The equal temperament sweep is here for coverage, not for judgement. A
//! temperament can be arbitrarily bad and plenty of these are; sweeping them
//! catches a panic or a rule that moved, and says nothing about whether an
//! answer is sensible. For that, read the named list - and the rank 2 and rank
//! 3 entries in it especially, which is what anyone actually notates.

use std::error::Error;

use xen_utils::{Notation, Subgroup, Temperament};

const TEMPERAMENTS: &str = include_str!("../data/temperaments_big.txt");
const SUBGROUPS: [&str; 6] = [
    "2.3.5",
    "2.3.7",
    "2.3.11",
    "2.3.5.7",
    "2.3.5.11",
    "2.3.5.7.11",
];

fn main() {
    for (number, line) in TEMPERAMENTS.lines().enumerate() {
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
        match parse(subgroup, definition) {
            Ok(t) => show(name, &t),
            Err(error) => eprintln!("line {}: {error}", number + 1),
        }
    }

    for subgroup in SUBGROUPS {
        let subgroup: Subgroup = subgroup.parse().unwrap();
        for divisions in 5..=99 {
            match Temperament::et(divisions, &subgroup) {
                Ok(t) => show(&format!("{divisions}et"), &t),
                Err(_) => println!("{divisions}et {subgroup}: contorted"),
            }
        }
    }
}

fn show(name: &str, t: &Temperament) {
    let subgroup = t.subgroup();
    let options = match Notation::options(t) {
        Ok(options) => options,
        Err(error) => {
            println!("{name} {subgroup}: {error}");
            return;
        }
    };
    let recommended = Notation::from_temperament(t).expect("the run is not empty");
    for n in &options {
        let mark = if n.generators() == recommended.generators() {
            "->"
        } else {
            "  "
        };
        let accidentals: Vec<String> = n.generators()[2..]
            .iter()
            .map(|a| ratio(subgroup, a))
            .collect();
        let spelling: Vec<String> = (2..subgroup.dim())
            .map(|index| {
                let mut harmonic = vec![0i64; subgroup.dim()];
                harmonic[index] = 1;
                harmonic[0] = -(subgroup.to_cents(&harmonic) / 1200.0).floor() as i64;
                let coordinates = n.spell(&harmonic).unwrap();
                let marks: i64 = coordinates[2..].iter().map(|c| c.abs()).sum();
                format!(
                    "{}={} ({marks})",
                    ratio(subgroup, &harmonic),
                    n.note(&coordinates)
                )
            })
            .collect();
        println!(
            "{name} {subgroup} {mark} [{}] {} | {} | {:?}",
            n.rank(),
            accidentals.join(" "),
            spelling.join(" "),
            n.generators()
        );
    }
}

/// An interval as a ratio, or as its exponents when it is too big for one.
fn ratio(subgroup: &Subgroup, interval: &[i64]) -> String {
    match subgroup.to_ratio(interval) {
        Ok((num, den)) => format!("{num}/{den}"),
        Err(_) => format!("{interval:?}"),
    }
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
