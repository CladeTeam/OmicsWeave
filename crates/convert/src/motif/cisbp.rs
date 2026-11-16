//! CIS-BP PWM files: one motif per file, a `Pos A C G T` (or `U`) header,
//! then a numbered row of probabilities per position. The ID is the file
//! name; a file with no rows (CIS-BP ships such shells) is passed over.

use std::path::Path;

use omicsweave_core::Result;
use omicsweave_ffm::Alphabet;

use super::{Entry, Parsed, bad, per_position, stem};

pub(crate) fn parse(path: &Path, text: &str) -> Result<Parsed> {
    let mut lines = text
        .lines()
        .enumerate()
        .filter(|(_, l)| !l.trim().is_empty());
    let Some((_, head)) = lines.next() else {
        return Err(bad(path, 0, "empty file"));
    };
    let columns: Vec<&str> = head.split_whitespace().collect();
    if columns.first() != Some(&"Pos") {
        return Err(bad(path, 1, "expected a `Pos A C G T` header"));
    }
    let letters = &columns[1..];
    let rna = letters.contains(&"U");
    let order: [&str; 4] = if rna {
        ["A", "C", "G", "U"]
    } else {
        ["A", "C", "G", "T"]
    };
    let index: Vec<usize> = order
        .iter()
        .map(|c| letters.iter().position(|x| x == c))
        .collect::<Option<_>>()
        .filter(|_| letters.len() == 4)
        .ok_or_else(|| bad(path, 1, format!("columns {letters:?} are not A C G T / U")))?;
    let mut rows = Vec::new();
    for (k, line) in lines {
        let words: Vec<&str> = line.split_whitespace().collect();
        if words.len() != 5 {
            return Err(bad(path, k + 1, format!("{} fields, need 5", words.len())));
        }
        let values = words[1..]
            .iter()
            .map(|w| w.parse::<f64>())
            .collect::<std::result::Result<Vec<_>, _>>()
            .map_err(|_| bad(path, k + 1, "values are not numbers"))?;
        rows.push(index.iter().map(|&j| values[j]).collect());
    }
    let alphabet = if rna {
        Alphabet::rna()
    } else {
        Alphabet::dna()
    };
    let id = stem(path);
    if rows.is_empty() {
        return Ok(Parsed {
            alphabet,
            entries: Vec::new(),
            skipped: vec![format!(
                "{}: motif {id:?} has no positions; skipped",
                path.display()
            )],
        });
    }
    Ok(Parsed {
        alphabet,
        entries: vec![Entry {
            id,
            name: String::new(),
            motif: per_position(rows, None),
        }],
        skipped: Vec::new(),
    })
}
