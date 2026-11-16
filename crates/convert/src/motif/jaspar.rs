//! JASPAR count matrices: `.jaspar` (`>ID name` then `A [ ... ]` rows) and
//! raw PFM (`>ID name` then four rows of numbers, A, C, G, T). Both store
//! counts; a motif's sites are its largest column sum.

use std::path::Path;

use omicsweave_core::Result;

use super::{Entry, Parsed, bad, header, numbers, per_position, stem};

/// A motif being collected: its header line and symbol rows.
struct Pending {
    line: usize,
    id: String,
    name: String,
    rows: Vec<(Option<char>, Vec<f64>)>,
}

fn finish(path: &Path, p: Pending, out: &mut Vec<Entry>) -> Result<()> {
    if p.rows.len() != 4 {
        return Err(bad(
            path,
            p.line,
            format!(
                "motif {:?} has {} rows, need 4 (A, C, G, T)",
                p.id,
                p.rows.len()
            ),
        ));
    }
    let width = p.rows[0].1.len();
    if width == 0 || p.rows.iter().any(|(_, r)| r.len() != width) {
        return Err(bad(
            path,
            p.line,
            format!("motif {:?} rows differ in length", p.id),
        ));
    }
    // Rows are A, C, G, T unless their letters say otherwise.
    let mut by_symbol: [Option<&Vec<f64>>; 4] = [None; 4];
    for (k, (letter, row)) in p.rows.iter().enumerate() {
        let slot = match letter {
            None => k,
            Some(c) => match c.to_ascii_uppercase() {
                'A' => 0,
                'C' => 1,
                'G' => 2,
                'T' | 'U' => 3,
                other => {
                    return Err(bad(
                        path,
                        p.line,
                        format!("motif {:?}: unknown row letter {other:?}", p.id),
                    ));
                }
            },
        };
        if by_symbol[slot].replace(row).is_some() {
            return Err(bad(
                path,
                p.line,
                format!("motif {:?}: a row repeats", p.id),
            ));
        }
    }
    let by_symbol = by_symbol.map(|r| r.expect("four unique rows"));
    let columns: Vec<Vec<f64>> = (0..width)
        .map(|j| {
            by_symbol
                .iter()
                .map(|r| {
                    // SAFETY: all four rows were checked to have `width`
                    // cells, and the outer loop has j < width.
                    unsafe { *r.get_unchecked(j) }
                })
                .collect()
        })
        .collect();
    let sites = columns
        .iter()
        .map(|c| c.iter().sum::<f64>())
        .fold(0.0f64, f64::max);
    out.push(Entry {
        id: p.id,
        name: p.name,
        motif: per_position(columns, (sites > 0.0).then_some(sites)),
    });
    Ok(())
}

pub(crate) fn parse_brackets(path: &Path, text: &str) -> Result<Parsed> {
    let mut entries = Vec::new();
    let mut pending: Option<Pending> = None;
    for (k, raw) in text.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('>') {
            if let Some(p) = pending.take() {
                finish(path, p, &mut entries)?;
            }
            let (id, name) = header(line);
            pending = Some(Pending {
                line: k + 1,
                id,
                name,
                rows: Vec::new(),
            });
            continue;
        }
        let p = pending
            .as_mut()
            .ok_or_else(|| bad(path, k + 1, "matrix row before a `>` header"))?;
        let (Some(open), Some(close)) = (line.find('['), line.rfind(']')) else {
            return Err(bad(path, k + 1, "expected a row such as `A [ 1 2 3 ]`"));
        };
        if close <= open {
            return Err(bad(path, k + 1, "expected `[` before `]`"));
        }
        let letter = line[..open].trim();
        let mut letters = letter.chars();
        let letter = match (letters.next(), letters.next()) {
            (Some(c), None) => c,
            _ => return Err(bad(path, k + 1, format!("row letter {letter:?}"))),
        };
        let values = numbers(&line[open + 1..close])
            .ok_or_else(|| bad(path, k + 1, "row values are not numbers"))?;
        p.rows.push((Some(letter), values));
    }
    if let Some(p) = pending {
        finish(path, p, &mut entries)?;
    }
    Ok(Parsed::dna(entries))
}

pub(crate) fn parse_raw(path: &Path, text: &str) -> Result<Parsed> {
    let mut entries = Vec::new();
    let mut pending: Option<Pending> = None;
    for (k, raw) in text.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if line.starts_with('>') {
            if let Some(p) = pending.take() {
                finish(path, p, &mut entries)?;
            }
            let (id, name) = header(line);
            pending = Some(Pending {
                line: k + 1,
                id,
                name,
                rows: Vec::new(),
            });
            continue;
        }
        // A file of one motif may have no header: its ID is the file name.
        let p = pending.get_or_insert_with(|| Pending {
            line: k + 1,
            id: stem(path),
            name: String::new(),
            rows: Vec::new(),
        });
        let values = numbers(line).ok_or_else(|| bad(path, k + 1, "row values are not numbers"))?;
        p.rows.push((None, values));
    }
    if let Some(p) = pending {
        finish(path, p, &mut entries)?;
    }
    Ok(Parsed::dna(entries))
}
