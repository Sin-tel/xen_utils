//! Everything a temperament tempers out, written in the notation basis.
//!
//! The just intonation notation - octave, fifth, one accidental per prime
//! beyond 3 - is a change of basis and nothing more: its generator matrix is
//! triangular with `+-1` down the diagonal, since each accidental has exponent
//! `+-1` on its own prime and support `{2, 3, p}`. So `to_notation` and
//! `to_just` are mutually inverse there, and any comma may be read in either
//! basis without losing anything.
//!
//! What that buys is legibility. `5120/5103` is `[10, -6, 1, -1]` over the
//! primes and `[0, 0, -1, 1]` in the notation basis, which says `81/80` and
//! `64/63` are one interval and says it on sight. This prints, for each
//! temperament:
//!
//! * the comma lattice, reduced, in both bases;
//! * the same lattice in Hermite normal form over the notation basis, which is
//!   canonical - it depends on the lattice and the column order, not on which
//!   basis went in;
//! * the commas `Notation::options` actually chose, in the notation basis,
//!   one line per notation in the run.
//!
//! The columns are ordered `[octave, fifth, a_1, a_2, ...]` throughout, with
//! `a_i` the accidental for the `i`th prime beyond 3.

use std::error::Error;

use diophantine::{Matrix, hnf};
use xen_utils::{Notation, Subgroup, Temperament};

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
    println!("\n{name} ({subgroup}), rank {}", t.rank());

    let ji = match Notation::from_ji(subgroup) {
        Ok(ji) => ji,
        Err(error) => {
            println!("    no just intonation notation: {error}");
            return;
        }
    };

    // The change of basis is only a change of basis if it is invertible, so
    // check that rather than assume it.
    for (index, generator) in ji.generators().iter().enumerate() {
        let mut unit = vec![0; ji.rank()];
        unit[index] = 1;
        assert_eq!(ji.to_notation(generator).unwrap(), unit);
        assert_eq!(ji.to_just(&unit).unwrap(), *generator);
    }

    println!(
        "    columns  [octave, fifth, {}]",
        ji.generators()[2..]
            .iter()
            .map(|a| ratio(subgroup, a))
            .collect::<Vec<_>>()
            .join(", ")
    );

    let commas = t.reduced_comma_basis().expect("a comma basis");
    println!("    comma lattice, reduced:");
    for comma in &commas {
        println!(
            "        {:14} {:22} {:?}",
            ratio(subgroup, comma),
            format!("{:?}", comma),
            ji.to_notation(comma).unwrap()
        );
    }

    let notational: Matrix<i64> = commas.iter().map(|c| ji.to_notation(c).unwrap()).collect();
    if !notational.is_empty() {
        println!("    the same lattice in hermite normal form:");
        for row in hnf(&notational).expect("a lattice basis") {
            println!("        {row:?}");
        }
    }

    let options = match Notation::options(t) {
        Ok(options) => options,
        Err(error) => {
            println!("    no notation: {error}");
            return;
        }
    };
    println!("    the commas each notation in the run chose:");
    for n in &options {
        let chosen: Vec<String> = n
            .commas()
            .iter()
            .map(|c| format!("{:?}", ji.to_notation(c).unwrap()))
            .collect();
        println!("        [{}]  {}", n.rank(), chosen.join("  "));
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
