//! GLAM2 text output (`glam2.txt`). Each alignment opens with
//! `Score: S  Columns: W  Sequences: N` and ends with its column table:
//!
//! ```text
//!  A  C  D ... Y Del Ins Score
//!  0  0  0 ...  3   0      5.92      <- column: residue counts, deletions
//!                       0 -0.102     <- junction: residues inserted
//!  0  0  3 ...  7   0      13.8
//! ```
//!
//! An alignment becomes a gapped profile: one joint factor of residue
//! counts per column, the deletion count of each column and the insertion
//! count of each junction, with `sites` the aligned sequences; a column
//! deleted in every sequence keeps its all-zero counts. The table
//! header names the alphabet (`a c g t` for DNA, the amino acids for
//! protein, else the letters of the alphabet file GLAM2 ran with). Scores
//! are derived from the counts and are not kept. IDs are the file name,
//! with `.<k>` (from 1, as GLAM2 numbers them) when the file holds several
//! alignments.

use std::path::Path;

use omicsweave_core::Result;
use omicsweave_ffm::{Alphabet, Factor, FactorKind, Gaps, Motif, Support};

use super::{Entry, Parsed, bad, stem};

/// One alignment's column table.
struct Table {
    line: usize,
    width: usize,
    sequences: f64,
    letters: Vec<u8>,
    counts: Vec<Vec<f64>>,
    deletions: Vec<f64>,
    insertions: Vec<f64>,
}

/// `Score: S  Columns: W  Sequences: N`: the width and aligned sequences.
fn score_line(path: &Path, k: usize, line: &str) -> Result<(usize, f64)> {
    let words: Vec<&str> = line.split_whitespace().collect();
    let field = |name: &str| {
        words
            .iter()
            .position(|w| *w == name)
            .and_then(|i| words.get(i + 1))
            .and_then(|v| v.parse::<u64>().ok())
            .ok_or_else(|| bad(path, k, format!("no `{name}` count on the score line")))
    };
    let width = u32::try_from(field("Columns:")?)
        .map_err(|_| bad(path, k, "column count exceeds the 32-bit position limit"))?;
    Ok((width as usize, field("Sequences:")? as f64))
}

fn counts(path: &Path, k: usize, words: &[&str]) -> Result<Vec<f64>> {
    words
        .iter()
        .map(|w| {
            w.parse::<u64>()
                .map(|v| v as f64)
                .map_err(|_| bad(path, k, format!("count {w:?} is not a whole number")))
        })
        .collect()
}

fn alphabet(path: &Path, line: usize, letters: &[u8]) -> Result<Alphabet> {
    let upper = letters.to_ascii_uppercase();
    if upper == b"ACGT" {
        return Ok(Alphabet::dna());
    }
    let protein = Alphabet::protein();
    if upper == protein.symbols() {
        return Ok(protein);
    }
    Alphabet::new(
        format!("GLAM2 {}", String::from_utf8_lossy(letters)),
        letters,
    )
    .map_err(|e| bad(path, line, e))
}

impl Table {
    /// The profile, or why it cannot be one (`Err` inside: skipped).
    fn motif(self, path: &Path) -> Result<std::result::Result<Motif, String>> {
        if self.width == 0 {
            return Ok(Err("the alignment has no columns".into()));
        }
        if self.counts.len() != self.width {
            return Err(bad(
                path,
                self.line,
                format!(
                    "the column table has {} columns, the score line {}",
                    self.counts.len(),
                    self.width
                ),
            ));
        }
        Ok(Ok(Motif {
            width: self.width as u32,
            support: Support::All,
            sites: (self.sequences > 0.0).then_some(self.sequences),
            factors: self
                .counts
                .into_iter()
                .enumerate()
                .map(|(p, table)| Factor {
                    kind: FactorKind::Joint,
                    positions: vec![p as u32],
                    table,
                })
                .collect(),
            gaps: Some(Gaps {
                deletions: self.deletions,
                insertions: self.insertions,
            }),
        }))
    }
}

pub(crate) fn parse(path: &Path, text: &str) -> Result<Parsed> {
    let mut tables: Vec<Table> = Vec::new();
    // Inside an alignment's column table, still reading rows.
    let mut reading = false;
    for (k, raw) in text.lines().enumerate() {
        let k = k + 1;
        let words: Vec<&str> = raw.split_whitespace().collect();
        if words.first() == Some(&"Score:") {
            let (width, sequences) = score_line(path, k, raw)?;
            tables.push(Table {
                line: k,
                width,
                sequences,
                letters: Vec::new(),
                counts: Vec::new(),
                deletions: Vec::new(),
                insertions: Vec::new(),
            });
            reading = false;
            continue;
        }
        let Some(t) = tables.last_mut() else {
            continue;
        };
        if words.ends_with(&["Del", "Ins", "Score"]) && t.letters.is_empty() {
            let letters = &words[..words.len() - 3];
            if letters.is_empty() || letters.iter().any(|w| w.len() != 1) {
                return Err(bad(path, k, "the column table header names no letters"));
            }
            t.letters = letters.iter().map(|w| w.as_bytes()[0]).collect();
            reading = true;
            continue;
        }
        if !reading {
            continue;
        }
        let n = t.letters.len();
        if words.is_empty() || t.counts.len() == t.width {
            reading = false;
        } else if words.len() == n + 2 {
            // Residue counts, deletions, then the column score.
            let c = counts(path, k, &words[..=n])?;
            if t.counts.len() != t.insertions.len() {
                return Err(bad(path, k, "two columns without a junction between"));
            }
            t.deletions.push(c[n]);
            t.counts.push(c[..n].to_vec());
        } else if words.len() == 2 {
            // Inserted residues, then the junction score.
            if t.counts.len() != t.insertions.len() + 1 {
                return Err(bad(path, k, "a junction not between two columns"));
            }
            t.insertions.push(counts(path, k, &words[..1])?[0]);
        } else {
            return Err(bad(
                path,
                k,
                format!(
                    "a table row of {} fields, need {} (column) or 2 (junction)",
                    words.len(),
                    n + 2
                ),
            ));
        }
    }
    let Some(first) = tables.first() else {
        return Err(bad(path, 0, "no alignment (`Score:` line)"));
    };
    if first.letters.is_empty() {
        return Err(bad(path, first.line, "the alignment has no column table"));
    }
    let letters = first.letters.clone();
    let alphabet = alphabet(path, first.line, &letters)?;
    let base = stem(path);
    let several = tables.len() > 1;
    let mut entries = Vec::with_capacity(tables.len());
    let mut skipped = Vec::new();
    for (i, t) in tables.into_iter().enumerate() {
        if t.letters != letters {
            return Err(bad(path, t.line, "alignments name different alphabets"));
        }
        let id = if several {
            format!("{base}.{}", i + 1)
        } else {
            base.clone()
        };
        match t.motif(path)? {
            Ok(motif) => entries.push(Entry {
                id,
                name: String::new(),
                motif,
            }),
            Err(why) => skipped.push(format!("{}: motif {id:?}: {why}; skipped", path.display())),
        }
    }
    Ok(Parsed {
        alphabet,
        entries,
        skipped,
    })
}

/// Whether text is GLAM2 output, by its banner.
pub(crate) fn is_glam2(first_line: &str) -> bool {
    first_line.starts_with("GLAM2: Gapped Local Alignment of Motifs")
}
