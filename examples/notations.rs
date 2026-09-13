//! Derives a notation for every temperament in `data/temperaments.txt` and
//! prints how each one spells the primes it tempers.
//!
//! Two notations are derived for each: the plain one, which only drops the
//! accidentals the temperament tempers out, and the minimal one, which has no
//! enharmonics at all and often does not exist. Each cell gives the rank of
//! the notation followed by its spelling of every prime beyond 3, octave
//! reduced. Run with `cargo run --example notations`.

use std::error::Error;

use xen_utils::{Notation, Subgroup, Temperament};

const TEMPERAMENTS: &str = include_str!("../data/temperaments.txt");

const HEADER: [&str; 5] = [
	"temperament",
	"subgroup",
	"rank",
	"plain notation",
	"minimal notation",
];

fn main() {
	let mut rows = vec![HEADER.map(String::from)];

	for (number, line) in TEMPERAMENTS.lines().enumerate() {
		let line = line.split('#').next().unwrap_or_default().trim();
		if line.is_empty() {
			continue;
		}
		match row(line) {
			Ok(row) => rows.push(row),
			Err(error) => eprintln!("line {}: {error}", number + 1),
		}
	}

	print_table(&rows);
}

/// Parses one line of the list and derives both notations for it.
fn row(line: &str) -> Result<[String; 5], Box<dyn Error>> {
	let mut fields = line.split('|').map(str::trim);
	let (Some(name), Some(subgroup), Some(definition), None) =
		(fields.next(), fields.next(), fields.next(), fields.next())
	else {
		return Err(format!("expected `name | subgroup | commas`, got {line:?}").into());
	};

	let subgroup: Subgroup = subgroup.parse()?;
	let temperament = parse_temperament(definition, &subgroup)?;

	Ok([
		name.to_string(),
		subgroup.to_string(),
		temperament.rank().to_string(),
		describe(&temperament, false),
		describe(&temperament, true),
	])
}

/// Parses the third field: `et N`, or a list of commas.
fn parse_temperament(definition: &str, subgroup: &Subgroup) -> Result<Temperament, Box<dyn Error>> {
	if let Some(divisions) = definition.strip_prefix("et ") {
		return Ok(Temperament::et(divisions.trim().parse()?, subgroup)?);
	}

	let commas = definition
		.split([',', ' '])
		.filter(|field| !field.is_empty())
		.map(|ratio| {
			let (num, den) = ratio
				.split_once('/')
				.ok_or_else(|| format!("{ratio:?} is not a ratio like `81/80`"))?;
			Ok(subgroup.factorize(num.parse()?, den.parse()?)?)
		})
		.collect::<Result<Vec<Vec<i64>>, Box<dyn Error>>>()?;

	Ok(Temperament::from_commas(&commas, subgroup)?)
}

/// The rank of the derived notation and how it spells each prime beyond 3, or
/// a dash if there is no such notation.
fn describe(temperament: &Temperament, minimal: bool) -> String {
	let Ok(notation) = Notation::from_temperament(temperament, minimal) else {
		return "-".to_string();
	};

	let subgroup = notation.subgroup();
	let mut cell = format!("[{}]", notation.rank());
	for index in 2..subgroup.dim() {
		let harmonic = octave_reduce(subgroup, index);
		let (num, den) = subgroup.to_ratio(&harmonic).expect("a single prime fits");
		let coordinates = notation
			.to_interval(&harmonic)
			.expect("built over subgroup");
		cell.push_str(&format!(" {num}/{den} {}", notation.note(&coordinates)));
	}
	cell
}

/// The prime at `index`, brought into the octave above the unison.
fn octave_reduce(subgroup: &Subgroup, index: usize) -> Vec<i64> {
	let mut interval = vec![0i64; subgroup.dim()];
	interval[index] = 1;
	interval[0] = -(subgroup.to_cents(&interval) / 1200.0).floor() as i64;
	interval
}

/// Prints rows padded to a common width, with a rule under the header.
fn print_table(rows: &[[String; 5]]) {
	let mut widths = [0; 5];
	for row in rows {
		for (width, cell) in widths.iter_mut().zip(row) {
			*width = (*width).max(cell.chars().count());
		}
	}

	for (number, row) in rows.iter().enumerate() {
		let line: Vec<String> = row
			.iter()
			.zip(widths)
			.map(|(cell, width)| format!("{cell:width$}"))
			.collect();
		println!("{}", line.join("  ").trim_end());

		if number == 0 {
			let rule: Vec<String> = widths.iter().map(|&width| "-".repeat(width)).collect();
			println!("{}", rule.join("  "));
		}
	}
}
