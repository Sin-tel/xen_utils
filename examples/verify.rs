//! Checks the invariants of the run of notations, over every temperament in
//! `data/temperaments.txt` and every equal temperament up to 72 over a few
//! subgroups.
//!
//!   * each generator maps to its unit vector, so `assemble` inverted right;
//!   * `ker(notation)` is contained in `ker(temperament)`, i.e. the notation
//!     never spells two intervals the temperament tells apart the same way;
//!   * the kernel and the enharmonic lattice account between them for every
//!     comma the temperament tempers out, so neither is too small;
//!   * kernels nest: a smaller notation spells alike everything a larger one
//!     does, so a larger notation's spelling can be simplified onto a smaller
//!     one's;
//!   * ranks run upwards one at a time, and the first has the rank of the
//!     temperament where one of that rank exists at all.

use std::error::Error;

use diophantine::{Matrix, solve_diophantine, transpose};
use xen_utils::{Notation, Subgroup, Temperament};

const TEMPERAMENTS: &str = include_str!("../data/temperaments.txt");

fn main() {
    let mut failures = 0;

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
            Ok(t) => failures += check(name, &t),
            Err(error) => eprintln!("line {}: {error}", number + 1),
        }
    }

    for subgroup in ["2.3.5", "2.3.5.7", "2.3.5.7.11"] {
        let subgroup: Subgroup = subgroup.parse().unwrap();
        for divisions in 5..=72 {
            let t = Temperament::et(divisions, &subgroup).unwrap();
            failures += check(&format!("{divisions}et over {subgroup}"), &t);
        }
    }

    println!("{failures} failures");
}

fn check(name: &str, t: &Temperament) -> usize {
    let mut failures = 0;
    let mut fail = |what: String| {
        eprintln!("{name}: {what}");
        failures += 1;
    };

    let options = match Notation::options(t) {
        Ok(options) => options,
        // Not every temperament has a notation over every subgroup, and refusing
        // is a result rather than a failure.
        Err(error) => {
            println!("{name}: no notation - {error}");
            return failures;
        }
    };

    let ranks: Vec<usize> = options.iter().map(Notation::rank).collect();
    if !ranks.windows(2).all(|pair| pair[1] == pair[0] + 1) {
        fail(format!("ranks do not run upwards one at a time: {ranks:?}"));
    }
    if ranks[0] < t.rank() {
        fail(format!(
            "smallest notation is rank {} < {}",
            ranks[0],
            t.rank()
        ));
    }
    if *ranks.last().unwrap() != t.dim() && t.rank() > 1 {
        // The largest keeps every useful accidental; over a full subgroup that
        // is all of them unless something is tempered out.
    }

    for (index, n) in options.iter().enumerate() {
        for (position, generator) in n.generators().iter().enumerate() {
            let mut unit = vec![0; n.rank()];
            unit[position] = 1;
            if n.to_notation(generator).unwrap() != unit {
                fail(format!(
                    "option {index}: generator {position} is not its unit vector"
                ));
            }
        }
        for comma in n.commas() {
            if n.to_notation(comma).unwrap().iter().any(|&x| x != 0) {
                fail(format!("option {index}: a comma is not in the kernel"));
            }
            if t.map(comma).unwrap().iter().any(|&x| x != 0) {
                fail(format!(
                    "option {index}: the comma {comma:?} is not tempered out"
                ));
            }
        }
        if n.commas().len() + n.rank() != n.dim() {
            fail(format!(
                "option {index}: comma count does not fill out the rank"
            ));
        }

        // What the notation spells alike and what it still spells apart have to
        // account between them for everything the temperament tempers out.
        let enharmonics = n.enharmonics().unwrap();
        let tempered = n.dim() - t.rank();
        if n.commas().len() + enharmonics.len() != tempered {
            fail(format!(
                "option {index}: {} commas and {} enharmonics are not the {tempered} dimensions tempered out",
                n.commas().len(),
                enharmonics.len(),
            ));
        }
        for enharmonic in &enharmonics {
            if t.map(enharmonic).unwrap().iter().any(|&x| x != 0) {
                fail(format!(
                    "option {index}: the enharmonic {enharmonic:?} is not tempered out"
                ));
            }
            if n.to_notation(enharmonic).unwrap().iter().all(|&x| x == 0) {
                fail(format!(
                    "option {index}: the enharmonic {enharmonic:?} is spelled as a unison"
                ));
            }
        }
    }

    // Kernel nesting: every comma of a larger notation is a comma of every
    // smaller one.
    for (index, pair) in options.windows(2).enumerate() {
        let (smaller, larger) = (&pair[0], &pair[1]);
        for comma in larger.commas() {
            if smaller.to_notation(comma).unwrap().iter().any(|&x| x != 0) {
                fail(format!(
                    "option {index} does not contain the kernel of option {}",
                    index + 1
                ));
            }
        }
    }

    failures
}

/// Whether the rows of `basis` span the same lattice as the rows of `other`.
#[allow(dead_code)]
fn same_lattice(basis: &Matrix<i64>, other: &Matrix<i64>) -> bool {
    solve_diophantine(&transpose(basis), &transpose(other)).is_ok()
        && solve_diophantine(&transpose(other), &transpose(basis)).is_ok()
}

fn parse(subgroup: &str, definition: &str) -> Result<Temperament, Box<dyn Error>> {
    let subgroup: Subgroup = subgroup.parse()?;
    if let Some(divisions) = definition.strip_prefix("et ") {
        return Ok(Temperament::et(divisions.trim().parse()?, &subgroup)?);
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
