//! ProBound / MotifCentral models. A MotifCentral fit is one JSON object; a
//! ProBound fit is JSON lines, the last line the final model. Each binding
//! mode of non-zero `size` is one motif:
//!
//! - `coefficients.bindingModes[i].mononucleotide[n * x + letter]`: the
//!   coefficient of `letter` at position `x`;
//! - `...dinucleotide[d - 1][n * n * x + n * a + b]`: of `a` at `x` and `b`
//!   at `x + d`;
//!
//! with `n` letters in `modelSettings.letterOrder`. A window's affinity is
//! `exp(sum of its coefficients)`, so each coefficient `beta` becomes a
//! potential `exp(beta)`. Activity, position bias, modifications and
//! interactions between modes are not part of one mode's model and are
//! left out. IDs are the file name, with `.<mode>` when several modes bind
//! specifically.

use std::path::Path;

use omicsweave_core::Result;
use omicsweave_ffm::{Alphabet, Factor, FactorKind, Motif, Support};
use serde_json::Value;

use super::{Entry, Parsed, bad, potentials, stem};

/// The model: the whole text when it is one JSON value, else its last
/// non-empty line (JSON lines).
fn model(path: &Path, text: &str) -> Result<Value> {
    if let Ok(v) = serde_json::from_str::<Value>(text) {
        return Ok(v);
    }
    let last = text
        .lines()
        .rev()
        .find(|l| !l.trim().is_empty())
        .ok_or_else(|| bad(path, 0, "empty file"))?;
    serde_json::from_str(last).map_err(|e| bad(path, 0, format!("not a ProBound JSON model: {e}")))
}

fn alphabet(path: &Path, settings: &Value) -> Result<Alphabet> {
    let letters = settings["letterOrder"]
        .as_str()
        .ok_or_else(|| bad(path, 0, "modelSettings.letterOrder is missing"))?;
    if letters == "ACGT" {
        return Ok(Alphabet::dna());
    }
    let mut alphabet = Alphabet::new(format!("ProBound {letters}"), letters.as_bytes())
        .map_err(|e| bad(path, 0, e))?;
    // `letterComplement` pairs letters: "C-G,A-T".
    if let Some(pairs) = settings["letterComplement"].as_str() {
        let mut partner = vec![None; letters.len()];
        for pair in pairs.split(',') {
            if let Some((a, b)) = pair.trim().split_once('-') {
                let (Some(i), Some(j)) = (letters.find(a), letters.find(b)) else {
                    continue;
                };
                partner[i] = Some(letters.as_bytes()[j]);
                partner[j] = Some(letters.as_bytes()[i]);
            }
        }
        if let Some(p) = partner.into_iter().collect::<Option<Vec<u8>>>() {
            alphabet = alphabet.with_complement(&p).map_err(|e| bad(path, 0, e))?;
        }
    }
    Ok(alphabet)
}

fn floats(path: &Path, v: &Value, what: &str) -> Result<Vec<f64>> {
    v.as_array()
        .ok_or_else(|| bad(path, 0, format!("{what} is not an array")))?
        .iter()
        .map(|x| {
            x.as_f64().ok_or_else(|| {
                bad(
                    path,
                    0,
                    format!("{what} holds a value that is not a number"),
                )
            })
        })
        .collect()
}

pub(crate) fn parse(path: &Path, text: &str) -> Result<Parsed> {
    let model = model(path, text)?;
    let settings = &model["modelSettings"];
    let alphabet = alphabet(path, settings)?;
    let n = alphabet.len();
    let modes = settings["bindingModes"]
        .as_array()
        .ok_or_else(|| bad(path, 0, "modelSettings.bindingModes is missing"))?;
    let coefficients = model["coefficients"]["bindingModes"]
        .as_array()
        .ok_or_else(|| bad(path, 0, "coefficients.bindingModes is missing"))?;
    let specific: Vec<usize> = (0..modes.len())
        .filter(|&i| modes[i]["size"].as_u64().unwrap_or(0) > 0)
        .collect();
    let mut entries = Vec::new();
    let base = stem(path);
    for &i in &specific {
        let what = |f: &str| format!("binding mode {i} {f}");
        let k = u32::try_from(modes[i]["size"].as_u64().unwrap_or(0))
            .map_err(|_| bad(path, 0, what("size exceeds the 32-bit position limit")))?
            as usize;
        let c = coefficients
            .get(i)
            .ok_or_else(|| bad(path, 0, what("has no coefficients")))?;
        let mono = floats(path, &c["mononucleotide"], &what("mononucleotide"))?;
        if Some(mono.len()) != n.checked_mul(k) {
            return Err(bad(
                path,
                0,
                what(&format!(
                    "has {} mononucleotide values, need {n} x {k}",
                    mono.len()
                )),
            ));
        }
        let mut factors: Vec<Factor> = mono
            .chunks(n)
            .enumerate()
            .map(|(x, t)| Factor {
                kind: FactorKind::Potential,
                positions: vec![x as u32],
                table: potentials(t.to_vec()),
            })
            .collect();
        let di: &[Value] = match &c["dinucleotide"] {
            Value::Null => &[],
            v => v
                .as_array()
                .ok_or_else(|| bad(path, 0, what("dinucleotide is not an array")))?
                .as_slice(),
        };
        for (slot, distance) in di.iter().enumerate() {
            let d = slot + 1;
            let values = floats(path, distance, &what(&format!("dinucleotide distance {d}")))?;
            if values.is_empty() {
                continue;
            }
            if d >= k || Some(values.len()) != (n * n).checked_mul(k - d) {
                return Err(bad(
                    path,
                    0,
                    what(&format!(
                        "distance {d} has {} values, need {n}^2 x {}",
                        values.len(),
                        k.saturating_sub(d)
                    )),
                ));
            }
            for (x, t) in values.chunks(n * n).enumerate() {
                factors.push(Factor {
                    kind: FactorKind::Potential,
                    positions: vec![x as u32, (x + d) as u32],
                    table: potentials(t.to_vec()),
                });
            }
        }
        entries.push(Entry {
            id: if specific.len() == 1 {
                base.clone()
            } else {
                format!("{base}.{i}")
            },
            name: String::new(),
            motif: Motif {
                width: k as u32,
                support: Support::All,
                sites: None,
                factors,
                gaps: None,
            },
        });
    }
    if entries.is_empty() {
        return Err(bad(path, 0, "no binding mode of non-zero size"));
    }
    Ok(Parsed {
        alphabet,
        entries,
        skipped: Vec::new(),
    })
}
