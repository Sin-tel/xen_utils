//! Lists the notational commas of every notation each temperament offers, next
//! to how it spells the primes. The commas are the choice being made: the
//! spelling is only one representative of its enharmonic class, but the kernel
//! is the notation.

use std::error::Error;

use diophantine::{Matrix, lll};
use xen_utils::{Notation, Subgroup, Temperament, Weighting};

const TEMPERAMENTS: &str = include_str!("../data/temperaments.txt");

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
}

fn show(name: &str, t: &Temperament) {
    let subgroup = t.subgroup();
    println!("{name} ({subgroup}), rank {}", t.rank());
    let options = match Notation::options(t) {
        Ok(options) => options,
        Err(error) => {
            println!("    {error}");
            return;
        }
    };
    for n in &options {
        let accidentals: Vec<String> = n.generators()[2..]
            .iter()
            .map(|a| ratio(subgroup, a))
            .collect();
        let commas: Vec<String> = reduced(subgroup, n.commas())
            .iter()
            .map(|c| ratio(subgroup, c))
            .collect();
        let spelling: Vec<String> = (2..subgroup.dim())
            .map(|index| {
                let mut harmonic = vec![0i64; subgroup.dim()];
                harmonic[index] = 1;
                harmonic[0] = -(subgroup.to_cents(&harmonic) / 1200.0).floor() as i64;
                format!(
                    "{} {}",
                    ratio(subgroup, &harmonic),
                    n.note(&n.to_interval(&harmonic).unwrap())
                )
            })
            .collect();
        println!(
            "    [{}] {:22} {:28} commas {}",
            n.rank(),
            accidentals.join(" "),
            spelling.join("  "),
            commas.join(" ")
        );
    }
}

/// A kernel basis reduced to small, ascending commas, the way
/// `Temperament::reduced_comma_basis` does it. The stored basis has one comma
/// per dropped accidental, which is the right thing to build the notation from
/// but not the right thing to read.
fn reduced(subgroup: &Subgroup, commas: &Matrix<i64>) -> Matrix<i64> {
    if commas.is_empty() {
        return Vec::new();
    }
    let weights = subgroup.weights(Weighting::Wilson);
    let basis = lll(commas, 0.99, &weights).expect("a kernel basis over this subgroup");
    basis
        .into_iter()
        .map(|comma| {
            if subgroup.to_cents(&comma) < 0.0 {
                comma.iter().map(|x| -x).collect()
            } else {
                comma
            }
        })
        .collect()
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
