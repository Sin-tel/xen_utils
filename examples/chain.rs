//! Collects the octave and fifth part of every comma the search picks, over the
//! temperament list and every equal temperament to 99 across six subgroups.
//!
//! The question is whether step 2 ever walks the fifth chain to anywhere except
//! a point where the chain nearly closes - the limma at five fifths, the apotome
//! at seven, the pythagorean comma at twelve, and so on up the convergents. If
//! it does not, then step 2 is the same operation `notation::accidental`
//! performs one level down: that one corrects a prime onto the fifth chain with
//! a small interval, and this one corrects an accidental.
//!
//! Only the recommended notation of each temperament is counted. The bottom of
//! a long run has had every accidental substituted away and walks wherever it
//! has to; what a caller actually gets handed is the question.

use std::collections::BTreeMap;
use std::error::Error;

use xen_utils::{Notation, Subgroup, Temperament};

const TEMPERAMENTS: &str = include_str!("../data/temperaments.txt");
const SUBGROUPS: [&str; 6] = [
    "2.3.5",
    "2.3.7",
    "2.3.11",
    "2.3.5.7",
    "2.3.5.11",
    "2.3.5.7.11",
];

fn main() {
    // Keyed by the fifth count, since that is what names the point on the chain.
    let mut seen: BTreeMap<i64, (i64, usize, String)> = BTreeMap::new();

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
            Ok(t) => collect(name, &t, &mut seen),
            Err(error) => eprintln!("line {}: {error}", number + 1),
        }
    }

    for subgroup in SUBGROUPS {
        let subgroup: Subgroup = subgroup.parse().unwrap();
        for divisions in 5..=99 {
            if let Ok(t) = Temperament::et(divisions, &subgroup) {
                collect(&format!("{divisions}et"), &t, &mut seen);
            }
        }
    }

    let total: usize = seen.values().map(|(_, count, _)| count).sum();
    println!(
        "{:>6}  {:>7}  {:>7}  {:>6}  {:>6}  {:>7}  {}",
        "fifths", "octaves", "cents", "letter", "sharps", "commas", "first seen"
    );
    for (fifths, (octaves, count, example)) in &seen {
        let cents = (f64::from(3).log2() - 1.0) * 1200.0 * *fifths as f64 + 1200.0 * *octaves as f64;
        // Seven fifths are a sharp and leave the letter alone, so the letter
        // moves by the fifth count modulo seven and the sharps are the rest.
        let letter = fifths.rem_euclid(7);
        let sharps = fifths.div_euclid(7);
        println!(
            "{fifths:6}  {octaves:7}  {cents:7.1}  {letter:6}  {sharps:6}  {count:7}  {example}"
        );
    }
    println!("
{total} commas over the recommended notations");
}

/// Records the octave and fifth part of each comma of each notation.
fn collect(name: &str, t: &Temperament, seen: &mut BTreeMap<i64, (i64, usize, String)>) {
    let subgroup = t.subgroup();
    let Ok(ji) = Notation::from_ji(subgroup) else {
        return;
    };
    let Ok(options) = Notation::options(t) else {
        return;
    };
    let Ok(recommended) = Notation::from_temperament(t) else {
        return;
    };
    for n in &options {
        if n.mapping() != recommended.mapping() {
            continue;
        }
        for comma in n.commas() {
            let notational = ji.to_notation(comma).unwrap();
            // Normalise the direction, so that a comma and its inverse are one
            // point on the chain rather than two.
            let (octaves, fifths) = if notational[1] < 0 || (notational[1] == 0 && notational[0] < 0)
            {
                (-notational[0], -notational[1])
            } else {
                (notational[0], notational[1])
            };
            let entry = seen
                .entry(fifths)
                .or_insert_with(|| (octaves, 0, format!("{name} {subgroup} [{}]", n.rank())));
            entry.1 += 1;
            assert_eq!(
                entry.0, octaves,
                "{fifths} fifths turned up with two different octave counts"
            );
        }
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
