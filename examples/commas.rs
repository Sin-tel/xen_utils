//! Lists the notational commas of every notation each temperament offers, next
//! to how it spells the primes. The commas are the choice being made: the
//! spelling is only one representative of its enharmonic class, but the kernel
//! is the notation.
//!
//! The enharmonics are the other half of that: the intervals the temperament
//! calls a unison but the notation still spells apart. Together the two account
//! for every comma the temperament tempers out, and the enharmonics are what a
//! spelling may be reduced by. An equal temperament always has some, since no
//! notation can close the circle of fifths - in 12et the enharmonic is the
//! pythagorean comma, which is what leaves `C#` and `Db` to differ.
//!
//! `->` marks the notation `Notation::from_temperament` recommends: the smallest
//! one that spells every prime on the nominal just intonation gives it.

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
    let recommended = Notation::from_temperament(t).expect("the run is not empty");
    for n in &options {
        let accidentals: Vec<String> = n.generators()[2..]
            .iter()
            .map(|a| ratio(subgroup, a))
            .collect();
        let commas: Vec<String> = reduced(subgroup, n.commas())
            .iter()
            .map(|c| ratio(subgroup, c))
            .collect();
        let enharmonics: Vec<String> = reduced(subgroup, &n.enharmonics().unwrap())
            .iter()
            .map(|e| ratio(subgroup, e))
            .collect();
        let spelling: Vec<String> = (2..subgroup.dim())
            .map(|index| {
                let mut harmonic = vec![0i64; subgroup.dim()];
                harmonic[index] = 1;
                harmonic[0] = -(subgroup.to_cents(&harmonic) / 1200.0).floor() as i64;
                format!(
                    "{} {}",
                    ratio(subgroup, &harmonic),
                    n.note(&n.to_notation(&harmonic).unwrap())
                )
            })
            .collect();
        let line = format!(
            "{} [{}] {:22} {:28} commas {:26} enharmonics {}",
            if n.mapping() == recommended.mapping() { " ->" } else { "   " },
            n.rank(),
            accidentals.join(" "),
            spelling.join("  "),
            commas.join(" "),
            enharmonics.join(" ")
        );
        println!("{}", line.trim_end());
    }
}

/// A lattice basis reduced to small, ascending intervals, the way
/// `Temperament::reduced_comma_basis` does it. Neither basis the library hands
/// back is reduced: the commas are one per dropped accidental and the
/// enharmonics are whatever the kernel computation gave, which is the right
/// thing to build on but not the right thing to read.
fn reduced(subgroup: &Subgroup, commas: &Matrix<i64>) -> Matrix<i64> {
    if commas.is_empty() {
        return Vec::new();
    }
    let weights = subgroup.weights(Weighting::Wilson);
    let basis = lll(commas, 0.99, &weights).expect("a kernel basis over this subgroup");
    basis
        .iter()
        .map(|comma| subgroup.ascending(comma))
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
