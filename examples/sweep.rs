//! Sweeps equal temperaments and looks for notations where a spelling has run
//! away: many accidental marks, or many sharps and flats.

use xen_utils::{Notation, Subgroup, Temperament};

fn main() {
    for subgroup in ["2.3.5", "2.3.5.7", "2.3.5.7.11"] {
        let subgroup: Subgroup = subgroup.parse().unwrap();
        for divisions in 5..=99 {
            let Ok(t) = Temperament::et(divisions, &subgroup) else {
                println!("{divisions}et over {subgroup}: contorted");
                continue;
            };
            let Ok(options) = Notation::options(&t) else {
                println!("{divisions}et over {subgroup}: no notation");
                continue;
            };
            for n in &options {
                for index in 2..subgroup.dim() {
                    let mut harmonic = vec![0i64; subgroup.dim()];
                    harmonic[index] = 1;
                    harmonic[0] = -(subgroup.to_cents(&harmonic) / 1200.0).floor() as i64;
                    let coordinates = n.to_notation(&harmonic).unwrap();
                    let marks: i64 = coordinates[2..].iter().map(|c| c.abs()).sum();
                    let sharps = (coordinates[1] + 1).div_euclid(7).abs();
                    if marks > 4 || sharps > 3 {
                        println!(
                            "{divisions}et over {subgroup} [{}]: {} -> {} ({marks} marks, {sharps} sharps)",
                            n.rank(),
                            subgroup.basis()[index],
                            n.note(&coordinates),
                        );
                    }
                }
            }
        }
    }
}
