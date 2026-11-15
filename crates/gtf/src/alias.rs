//! Chromosome names: renaming at read, and matching an annotation's
//! seqids to another file's reference names (a BAM header, a FASTA index).

use std::collections::HashMap;

use omicsweave_core::{Error, Result};

/// Name pairs, used in both directions.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Aliases {
    to: HashMap<String, String>,
}

impl Aliases {
    pub fn new() -> Aliases {
        Aliases::default()
    }

    /// Reads `from` as `to` (and, when matching, `to` as `from`).
    pub fn insert(&mut self, from: impl Into<String>, to: impl Into<String>) {
        self.to.insert(from.into(), to.into());
    }

    pub fn from_pairs<A: Into<String>, B: Into<String>>(
        pairs: impl IntoIterator<Item = (A, B)>,
    ) -> Aliases {
        let mut a = Aliases::new();
        for (from, to) in pairs {
            a.insert(from, to);
        }
        a
    }

    /// Two whitespace-separated columns per line (`from to`), as in a
    /// UCSC `chromAlias` table reduced to two columns; `#` lines and blank
    /// lines are skipped.
    pub fn parse(text: &str) -> Result<Aliases> {
        let mut a = Aliases::new();
        for (i, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut cols = line.split_whitespace();
            match (cols.next(), cols.next()) {
                (Some(from), Some(to)) => a.insert(from, to),
                _ => {
                    return Err(Error::Format {
                        line: i as u64 + 1,
                        detail: "alias line needs two names".into(),
                    });
                }
            }
        }
        Ok(a)
    }

    /// The name `name` is read as.
    pub fn rename<'a>(&'a self, name: &'a str) -> &'a str {
        self.to.get(name).map_or(name, String::as_str)
    }

    pub fn len(&self) -> usize {
        self.to.len()
    }

    pub fn is_empty(&self) -> bool {
        self.to.is_empty()
    }

    fn reverse(&self, name: &str) -> impl Iterator<Item = &str> {
        self.to
            .iter()
            .filter(move |(_, to)| to.as_str() == name)
            .map(|(from, _)| from.as_str())
    }
}

/// UCSC ↔ Ensembl spellings of `name`: with or without `chr`, and the
/// mitochondrion as `chrM`, `MT` or `M`.
pub fn builtin_alternatives(name: &str) -> Vec<String> {
    let mut out = Vec::with_capacity(3);
    match name {
        "chrM" | "chrMT" => out.extend(["MT".to_string(), "M".to_string(), "chrM".to_string()]),
        "MT" | "M" => out.extend(["chrM".to_string(), "MT".to_string(), "M".to_string()]),
        _ => {}
    }
    if let Some(rest) = name.strip_prefix("chr") {
        out.push(rest.to_string());
    } else {
        out.push(format!("chr{name}"));
    }
    out.retain(|n| n != name);
    out
}

/// For each of `targets`, the index of the name in `names` it refers to:
/// the same name, else an alias of it (either direction), else a built-in
/// UCSC ↔ Ensembl spelling ([`builtin_alternatives`]).
pub fn match_names(
    names: &[&str],
    targets: &[&str],
    aliases: Option<&Aliases>,
) -> Vec<Option<u32>> {
    let by_name: HashMap<&str, u32> = names
        .iter()
        .enumerate()
        .rev()
        .map(|(i, &n)| (n, i as u32))
        .collect();
    targets
        .iter()
        .map(|&t| {
            if let Some(&i) = by_name.get(t) {
                return Some(i);
            }
            if let Some(a) = aliases {
                if let Some(&i) = by_name.get(a.rename(t)) {
                    return Some(i);
                }
                for from in a.reverse(t) {
                    if let Some(&i) = by_name.get(from) {
                        return Some(i);
                    }
                }
            }
            builtin_alternatives(t)
                .iter()
                .find_map(|alt| by_name.get(alt.as_str()).copied())
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_ucsc_and_ensembl_names() {
        let names = ["1", "2", "MT", "X"];
        let targets = ["chr1", "2", "chrM", "chrX", "chrUn"];
        assert_eq!(
            match_names(&names, &targets, None),
            vec![Some(0), Some(1), Some(2), Some(3), None]
        );
        let a = Aliases::parse("# from to\nchrUn\tX\n").unwrap();
        assert_eq!(match_names(&names, &["chrUn"], Some(&a)), vec![Some(3)]);
    }
}
