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
//! Every notation in every run is counted, not only the recommended one: a
//! recommendation that keeps all its accidentals has no commas at all, so
//! counting those alone throws most of rank 2 away and leaves nothing to read.
//!
//! Rank 1 and the higher ranks are counted apart. An equal temperament sweep is
//! a coverage test rather than a judgement - most of what it turns up is
//! nobody's notation - so the rank 2 and up column is the one that says whether
//! a walk is a reasonable thing to have done. That column is thin, which is the
//! honest state of the evidence: the named list is short on rank 2 and rank 3
//! entries, and lengthening it is the way to fill it in.

use std::collections::BTreeMap;
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
    // Keyed by the fifth count, since that is what names the point on the chain.
    let mut seen: BTreeMap<i64, Row> = BTreeMap::new();

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

    // for subgroup in SUBGROUPS {
    //     let subgroup: Subgroup = subgroup.parse().unwrap();
    //     for divisions in 5..=99 {
    //         if let Ok(t) = Temperament::et(divisions, &subgroup) {
    //             collect(&format!("{divisions}et"), &t, &mut seen);
    //         }
    //     }
    // }

    let equal: usize = seen.values().map(|row| row.equal).sum();
    let higher: usize = seen.values().map(|row| row.higher).sum();
    println!(
        "{:>6}  {:>7}  {:>7}  {:>6}  {:>6}  {:>7}  {:>9}  first seen above rank 1",
        "fifths", "octaves", "cents", "letter", "sharps", "rank 1", "rank 2 up",
    );
    for (fifths, row) in &seen {
        let cents =
            (f64::from(3).log2() - 1.0) * 1200.0 * *fifths as f64 + 1200.0 * row.octaves as f64;
        // Seven fifths are a sharp and leave the letter alone, so the letter
        // moves by the fifth count modulo seven and the sharps are the rest.
        let letter = fifths.rem_euclid(7);
        let sharps = fifths.div_euclid(7);
        println!(
            "{fifths:6}  {:7}  {cents:7.1}  {letter:6}  {sharps:6}  {:7}  {:9}  {}",
            row.octaves, row.equal, row.higher, row.example
        );
    }
    println!(
        "
{equal} commas from equal temperaments, {higher} from rank 2 and up"
    );
}

/// One point on the chain of fifths: its octave count, how many commas landed
/// on it from a rank 1 temperament and from a higher one, and where it first
/// turned up above rank 1.
#[derive(Default)]
struct Row {
    octaves: i64,
    equal: usize,
    higher: usize,
    example: String,
}

/// Records the octave and fifth part of each comma of each notation.
fn collect(name: &str, t: &Temperament, seen: &mut BTreeMap<i64, Row>) {
    let subgroup = t.subgroup();
    let Ok(ji) = Notation::from_ji(subgroup) else {
        return;
    };
    let Ok(options) = Notation::options(t) else {
        return;
    };
    for n in &options {
        for comma in n.commas() {
            let notational = ji.to_notation(comma).unwrap();
            // Normalise the direction, so that a comma and its inverse are one
            // point on the chain rather than two.
            let (octaves, fifths) =
                if notational[1] < 0 || (notational[1] == 0 && notational[0] < 0) {
                    (-notational[0], -notational[1])
                } else {
                    (notational[0], notational[1])
                };
            let entry = seen.entry(fifths).or_default();
            entry.octaves = octaves;
            if t.rank() == 1 {
                entry.equal += 1;
            } else {
                entry.higher += 1;
                if entry.example.is_empty() {
                    entry.example = format!("{name} {subgroup} [{}]", n.rank());
                }
            }
            assert_eq!(
                entry.octaves, octaves,
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
