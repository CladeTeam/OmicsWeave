//! TRANSFAC matrix records: two-letter line codes, a `P0` (or `PO`) header
//! naming the columns, numbered position rows (an optional consensus letter
//! at the end), records ended by `//`. The ID is `AC` (else `ID`). The name
//! is the `ID` line when `AC` is the ID (JASPAR's TF name, HOCOMOCO's model
//! name), else `NA`.

use std::path::Path;

use omicsweave_core::Result;
use omicsweave_ffm::Alphabet;

use super::{Entry, Parsed, bad, per_position};

/// Whether text looks like TRANSFAC: line codes and a `P0` / `PO` header.
pub(crate) fn is_transfac(text: &str) -> bool {
    let mut coded = false;
    for line in text
        .lines()
        .map(str::trim_end)
        .filter(|l| !l.trim().is_empty())
    {
        let code = line.get(..2).unwrap_or(line);
        if code == "P0" || code == "PO" {
            return coded || line.len() > 2;
        }
        if !(code.chars().all(|c| c.is_ascii_uppercase()) && code.len() == 2) && line != "//" {
            return false;
        }
        coded = true;
    }
    false
}

#[derive(Default)]
struct Record {
    line: usize,
    ac: Option<String>,
    id: Option<String>,
    na: Option<String>,
    columns: Vec<char>,
    rows: Vec<Vec<f64>>,
}

pub(crate) fn parse(path: &Path, text: &str) -> Result<Parsed> {
    let mut entries = Vec::new();
    let mut rna = None::<bool>;
    let mut rec = Record::default();
    let mut in_matrix = false;
    for (k, raw) in text.lines().enumerate() {
        let n = k + 1;
        let line = raw.trim_end();
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with("//") {
            finish(path, std::mem::take(&mut rec), &mut entries, &mut rna)?;
            in_matrix = false;
            continue;
        }
        let code = line.get(..2).unwrap_or(line);
        let value = line.get(2..).unwrap_or("").trim();
        if rec.line == 0 {
            rec.line = n;
        }
        match code {
            "AC" => rec.ac = Some(value.to_string()),
            "ID" => rec.id = Some(value.to_string()),
            "NA" => rec.na = Some(value.to_string()),
            "P0" | "PO" => {
                rec.columns = value
                    .split_whitespace()
                    .map(|w| w.chars().next().unwrap_or(' ').to_ascii_uppercase())
                    .collect();
                in_matrix = true;
            }
            "XX" => in_matrix = false,
            _ if in_matrix && code.chars().all(|c| c.is_ascii_digit()) => {
                let words: Vec<&str> = line.split_whitespace().collect();
                let w = rec.columns.len();
                if words.len() < 1 + w {
                    return Err(bad(
                        path,
                        n,
                        format!("{} columns, need {w}", words.len() - 1),
                    ));
                }
                let row = words[1..=w]
                    .iter()
                    .map(|v| v.parse::<f64>())
                    .collect::<std::result::Result<Vec<_>, _>>()
                    .map_err(|_| bad(path, n, "position values are not numbers"))?;
                rec.rows.push(row);
            }
            _ => {}
        }
    }
    if !rec.rows.is_empty() {
        finish(path, rec, &mut entries, &mut rna)?;
    }
    Ok(Parsed {
        alphabet: if rna == Some(true) {
            Alphabet::rna()
        } else {
            Alphabet::dna()
        },
        entries,
        skipped: Vec::new(),
    })
}

fn finish(path: &Path, rec: Record, out: &mut Vec<Entry>, rna: &mut Option<bool>) -> Result<()> {
    if rec.rows.is_empty() {
        return Ok(());
    }
    let (id, name) = match (rec.ac, rec.id, rec.na) {
        (Some(ac), id, na) => (ac, id.or(na).unwrap_or_default()),
        (None, Some(id), na) => (id, na.unwrap_or_default()),
        (None, None, _) => return Err(bad(path, rec.line, "matrix record without AC or ID")),
    };
    // Columns A, C, G, T (or U) in any order.
    let is_rna = rec.columns.contains(&'U');
    if rna.replace(is_rna).is_some_and(|r| r != is_rna) {
        return Err(bad(path, rec.line, "records mix DNA and RNA columns"));
    }
    let order = if is_rna {
        ['A', 'C', 'G', 'U']
    } else {
        ['A', 'C', 'G', 'T']
    };
    let index: Vec<usize> = order
        .iter()
        .map(|c| rec.columns.iter().position(|x| x == c))
        .collect::<Option<_>>()
        .filter(|_| rec.columns.len() == 4)
        .ok_or_else(|| {
            bad(
                path,
                rec.line,
                format!(
                    "matrix {id:?}: P0 columns {:?} are not A C G T",
                    rec.columns
                ),
            )
        })?;
    let rows = rec
        .rows
        .into_iter()
        .map(|r| index.iter().map(|&j| r[j]).collect())
        .collect();
    out.push(Entry {
        id,
        name,
        motif: per_position(rows, None),
    });
    Ok(())
}
