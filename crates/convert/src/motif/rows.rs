//! Motifs written one row per position under a `>` header: counts or
//! frequencies (HOCOMOCO `.pcm` / `.pfm`, Cluster-Buster), log weights
//! (HOCOMOCO `.pwm`), their dinucleotide forms (`.dpcm`, `.dpwm`: one row
//! of sixteen values per pair of adjacent positions, `AA AC ... TT`), and
//! HOMER `.motif` files (probabilities under a `>consensus name threshold`
//! header).

use std::path::Path;

use omicsweave_core::Result;
use omicsweave_ffm::{Factor, FactorKind, Motif, Support};

use super::{Entry, Parsed, bad, count_sites, header, numbers, per_position, potentials};
use crate::format::Format;

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Kind {
    Counts,
    Weights,
    DiCounts,
    DiWeights,
}

impl Kind {
    fn columns(self) -> usize {
        match self {
            Kind::Counts | Kind::Weights => 4,
            Kind::DiCounts | Kind::DiWeights => 16,
        }
    }
}

/// A `>` record: header line number, ID, name, rows.
type Record<'a> = (usize, &'a str, String, String, Vec<Vec<f64>>);

fn records<'a>(path: &Path, text: &'a str) -> Result<Vec<Record<'a>>> {
    let mut out: Vec<Record<'a>> = Vec::new();
    for (k, raw) in text.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if line.starts_with('>') {
            let (id, name) = header(line);
            out.push((k + 1, line, id, name, Vec::new()));
            continue;
        }
        let Some(rec) = out.last_mut() else {
            return Err(bad(path, k + 1, "matrix row before a `>` header"));
        };
        let row = numbers(line).ok_or_else(|| bad(path, k + 1, "row values are not numbers"))?;
        rec.4.push(row);
    }
    Ok(out)
}

pub(crate) fn parse(path: &Path, text: &str, kind: Kind) -> Result<Parsed> {
    let mut entries = Vec::new();
    for (line, _, id, name, rows) in records(path, text)? {
        let cols = kind.columns();
        if rows.is_empty() || rows.iter().any(|r| r.len() != cols) {
            return Err(bad(
                path,
                line,
                format!("motif {id:?}: every row needs {cols} values"),
            ));
        }
        let motif = match kind {
            Kind::Counts => per_position(rows, None),
            Kind::Weights => Motif {
                width: rows.len() as u32,
                support: Support::All,
                sites: None,
                factors: rows
                    .into_iter()
                    .enumerate()
                    .map(|(p, w)| Factor {
                        kind: FactorKind::Potential,
                        positions: vec![p as u32],
                        table: potentials(w),
                    })
                    .collect(),
                gaps: None,
            },
            Kind::DiCounts | Kind::DiWeights => {
                let counts = kind == Kind::DiCounts;
                let sites = if counts { count_sites(&rows) } else { None };
                Motif {
                    width: rows.len() as u32 + 1,
                    support: Support::All,
                    sites,
                    factors: rows
                        .into_iter()
                        .enumerate()
                        .map(|(p, t)| Factor {
                            kind: if counts {
                                FactorKind::Joint
                            } else {
                                FactorKind::Potential
                            },
                            positions: vec![p as u32, p as u32 + 1],
                            table: if counts { t } else { potentials(t) },
                        })
                        .collect(),
                    gaps: None,
                }
            }
        };
        entries.push(Entry { id, name, motif });
    }
    Ok(Parsed::dna(entries))
}

/// HOMER: `>consensus<TAB>name<TAB>log-odds threshold[...]`, then rows of
/// four probabilities. The ID and the name are the name column.
pub(crate) fn parse_homer(path: &Path, text: &str) -> Result<Parsed> {
    let mut entries = Vec::new();
    for (line, head, _, _, rows) in records(path, text)? {
        let fields: Vec<&str> = head.trim_start_matches('>').split('\t').collect();
        let Some(name) = fields.get(1).map(|s| s.trim()).filter(|s| !s.is_empty()) else {
            return Err(bad(path, line, "HOMER header without a name column"));
        };
        if rows.is_empty() || rows.iter().any(|r| r.len() != 4) {
            return Err(bad(
                path,
                line,
                format!("motif {name:?}: every row needs 4 values"),
            ));
        }
        entries.push(Entry {
            id: name.to_string(),
            name: name.to_string(),
            motif: per_position(rows, None),
        });
    }
    Ok(Parsed::dna(entries))
}

/// What `>`-headed text is.
pub(crate) enum Headed {
    /// Not numeric matrices (FASTA).
    Text,
    Motif(Format),
    /// Numeric matrices whose layout cannot be told.
    Ambiguous,
}

/// Tells the `>`-headed numeric formats apart.
pub(crate) fn sniff_headed(text: &str, complete: bool) -> Headed {
    match sniff_shape(text, complete) {
        Some(Some(f)) => Headed::Motif(f),
        Some(None) => Headed::Ambiguous,
        None => Headed::Text,
    }
}

/// `None` when the text is not numeric matrices; `Some(None)` when it is
/// but the layout is ambiguous.
fn sniff_shape(text: &str, complete: bool) -> Option<Option<Format>> {
    let mut motifs: Vec<(Vec<&str>, Vec<usize>, bool)> = Vec::new();
    let mut lines: Vec<&str> = text.lines().collect();
    // The last line of a sniffed prefix may be cut.
    if !complete && !text.ends_with('\n') {
        lines.pop();
    }
    for raw in lines {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some(head) = line.strip_prefix('>') {
            motifs.push((head.split('\t').collect(), Vec::new(), false));
            continue;
        }
        let rec = motifs.last_mut()?;
        if line.contains('[') && line.contains(']') {
            return Some(Some(Format::Jaspar));
        }
        let row = numbers(line)?;
        if row.is_empty() {
            return None;
        }
        rec.1.push(row.len());
        rec.2 |= row.iter().any(|&v| v < 0.0);
    }
    // Drop a last motif the prefix may have cut.
    if !complete && motifs.len() > 1 {
        motifs.pop();
    }
    if motifs.is_empty() || motifs.iter().any(|m| m.1.is_empty()) {
        return None;
    }
    let negative = motifs.iter().any(|m| m.2);
    let homer = motifs.iter().all(|m| {
        m.0.len() >= 3 && m.0[2].trim().parse::<f64>().is_ok() && m.1.iter().all(|&c| c == 4)
    });
    if homer {
        return Some(Some(Format::Homer));
    }
    let all_cols = |c: usize| motifs.iter().all(|m| m.1.iter().all(|&x| x == c));
    if all_cols(16) {
        return Some(Some(if negative { Format::Dpwm } else { Format::Dpcm }));
    }
    let rows_of_four = all_cols(4);
    let four_rows = motifs
        .iter()
        .all(|m| m.1.len() == 4 && m.1.iter().all(|&c| c == m.1[0]));
    Some(match (rows_of_four, four_rows) {
        (true, false) => Some(if negative { Format::Pwm } else { Format::Pcm }),
        (false, true) if !negative => Some(Format::JasparPfm),
        _ => None,
    })
}
