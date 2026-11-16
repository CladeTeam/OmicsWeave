//! MEME motif format (versions 4 and 5): `MOTIF id [alt]` records with a
//! `letter-probability matrix`. The background, strands, `log-odds matrix`
//! and `URL` lines are not part of a model and are not read.

use std::path::Path;

use omicsweave_core::Result;
use omicsweave_ffm::Alphabet;

use super::{Entry, Parsed, bad, numbers, per_position};

const PROTEIN: &str = "ACDEFGHIKLMNPQRSTVWY";

pub(crate) fn parse(path: &Path, text: &str) -> Result<Parsed> {
    let lines: Vec<&str> = text.lines().collect();
    let mut alphabet: Option<Alphabet> = None;
    let mut entries = Vec::new();
    let mut skipped = Vec::new();
    let mut i = 0;
    let mut seen_version = false;
    let mut background = String::new();
    while i < lines.len() {
        let line = lines[i].trim();
        let n = i + 1;
        i += 1;
        if line.is_empty() {
            continue;
        }
        if line.starts_with("MEME version") {
            seen_version = true;
        } else if let Some(rest) = line.strip_prefix("ALPHABET=") {
            alphabet = Some(short_alphabet(path, n, rest.trim())?);
        } else if line.starts_with("Background letter frequencies") {
            // Letter / frequency pairs, possibly over several lines.
            while let Some(l) = lines.get(i).map(|l| l.trim()).filter(|l| !l.is_empty()) {
                if l.starts_with("MOTIF") {
                    break;
                }
                background.push_str(l);
                background.push(' ');
                i += 1;
            }
        } else if line.starts_with("ALPHABET") {
            let (a, next) = custom_alphabet(path, &lines, i - 1)?;
            alphabet = Some(a);
            i = next;
        } else if let Some(rest) = line.strip_prefix("MOTIF") {
            let mut words = rest.split_whitespace();
            let id = words
                .next()
                .ok_or_else(|| bad(path, n, "MOTIF without an identifier"))?
                .to_string();
            let name = words.next().unwrap_or("").to_string();
            // Without an alphabet line MEME guesses it; so do we, from the
            // matrix width and the background letters.
            let alphabet = alphabet.get_or_insert_with(|| guess(&lines[i..], &background));
            match matrix(path, &lines, &mut i, alphabet.len())? {
                Some((rows, nsites)) if !rows.is_empty() => entries.push(Entry {
                    id,
                    name,
                    motif: per_position(rows, nsites),
                }),
                _ => skipped.push(format!(
                    "{}: motif {id:?} has no letter-probability matrix; skipped",
                    path.display()
                )),
            }
        }
    }
    if !seen_version {
        return Err(bad(
            path,
            0,
            "not a MEME motif file (no `MEME version` line)",
        ));
    }
    Ok(Parsed {
        alphabet: alphabet.unwrap_or_else(Alphabet::dna),
        entries,
        skipped,
    })
}

/// DNA unless the next matrix has 20 columns (protein) or the background
/// names `U` but not `T` (RNA).
fn guess(rest: &[&str], background: &str) -> Alphabet {
    let alength = rest
        .iter()
        .map(|l| l.trim())
        .take_while(|l| !l.starts_with("MOTIF"))
        .find_map(|l| l.strip_prefix("letter-probability matrix"))
        .and_then(|h| {
            key_values(h.trim_start_matches(':'))
                .into_iter()
                .find(|(k, _)| *k == "alength")
                .and_then(|(_, v)| v.parse::<usize>().ok())
        });
    let letters: Vec<&str> = background.split_whitespace().step_by(2).collect();
    if alength == Some(20) || letters.len() == 20 {
        Alphabet::protein()
    } else if letters.contains(&"U") && !letters.contains(&"T") {
        Alphabet::rna()
    } else {
        Alphabet::dna()
    }
}

/// `ALPHABET= ACGT`: the letters of one of MEME's standard alphabets.
fn short_alphabet(path: &Path, line: usize, letters: &str) -> Result<Alphabet> {
    match letters {
        "ACGT" => Ok(Alphabet::dna()),
        "ACGU" => Ok(Alphabet::rna()),
        PROTEIN => Ok(Alphabet::protein()),
        other => Err(bad(
            path,
            line,
            format!("unknown alphabet {other:?}; define it with an ALPHABET block"),
        )),
    }
}

/// An `ALPHABET "name" [type] ... END ALPHABET` block starting at line
/// index `start`; returns the alphabet and the index after the block.
fn custom_alphabet(path: &Path, lines: &[&str], start: usize) -> Result<(Alphabet, usize)> {
    let head = lines[start].trim();
    let name = head
        .split('"')
        .nth(1)
        .unwrap_or_else(|| head.trim_start_matches("ALPHABET").trim())
        .to_string();
    let mut core: Vec<u8> = Vec::new();
    let mut pairs: Vec<(u8, u8)> = Vec::new();
    let mut aliases: Vec<(u8, u8)> = Vec::new();
    let mut i = start + 1;
    loop {
        let Some(raw) = lines.get(i) else {
            return Err(bad(path, start + 1, "ALPHABET block without END ALPHABET"));
        };
        let line = strip_comment(raw).trim();
        i += 1;
        if line.is_empty() {
            continue;
        }
        if line.starts_with("END ALPHABET") {
            break;
        }
        let n = i;
        let symbol = |word: &str| -> Result<u8> {
            match word.as_bytes() {
                [b] => Ok(*b),
                _ => Err(bad(
                    path,
                    n,
                    format!("symbol {word:?} is not one character"),
                )),
            }
        };
        let words = symbol_words(line);
        let Some(first) = words.first() else {
            return Err(bad(path, n, "alphabet line without a symbol"));
        };
        if let Some(eq) = words.iter().position(|w| w == "=") {
            let sym = symbol(first)?;
            let targets = words.get(eq + 1).map(String::as_bytes).unwrap_or(&[]);
            // One core symbol: an alias. Several: an ambiguous symbol, which
            // is invalid input in FFM.
            if let [t] = targets {
                aliases.push((sym, *t));
            }
        } else if let Some(t) = words.iter().position(|w| w == "~") {
            let a = symbol(first)?;
            let b = symbol(
                words
                    .get(t + 1)
                    .ok_or_else(|| bad(path, n, "complement without a partner"))?,
            )?;
            core.extend([a, b]);
            pairs.push((a, b));
        } else {
            core.push(symbol(first)?);
        }
    }
    let order = meme_order(&core);
    let mut alphabet = Alphabet::new(name, &order).map_err(|e| bad(path, start + 1, e))?;
    if !pairs.is_empty() {
        let partner = |s: u8| {
            pairs.iter().find_map(|&(a, b)| {
                if a == s {
                    Some(b)
                } else if b == s {
                    Some(a)
                } else {
                    None
                }
            })
        };
        if let Some(partners) = order
            .iter()
            .map(|&s| partner(s))
            .collect::<Option<Vec<u8>>>()
        {
            alphabet = alphabet
                .with_complement(&partners)
                .map_err(|e| bad(path, start + 1, e))?;
        }
    }
    for (alias, target) in aliases {
        if !order.contains(&alias) {
            alphabet = alphabet
                .with_alias(alias, target)
                .map_err(|e| bad(path, start + 1, e))?;
        }
    }
    Ok((alphabet, i))
}

/// A line up to a `#` outside quotes.
fn strip_comment(line: &str) -> &str {
    let mut quoted = false;
    for (k, c) in line.char_indices() {
        match c {
            '"' => quoted = !quoted,
            '#' if !quoted => return &line[..k],
            _ => {}
        }
    }
    line
}

/// The words of an alphabet line, with quoted names and colors dropped and
/// `~` / `=` as words of their own.
fn symbol_words(line: &str) -> Vec<String> {
    let mut words = Vec::new();
    let mut quoted = false;
    let mut cur = String::new();
    for c in line.chars() {
        match c {
            '"' => quoted = !quoted,
            _ if quoted => {}
            '~' | '=' => {
                if !cur.is_empty() {
                    words.push(std::mem::take(&mut cur));
                }
                words.push(c.to_string());
            }
            c if c.is_whitespace() => {
                if !cur.is_empty() {
                    words.push(std::mem::take(&mut cur));
                }
            }
            c => cur.push(c),
        }
    }
    if !cur.is_empty() {
        words.push(cur);
    }
    // A color (six hex digits) follows a symbol or its name.
    let mut out = Vec::new();
    for (k, w) in words.iter().enumerate() {
        let is_color = k > 0
            && w.len() == 6
            && w.chars().all(|c| c.is_ascii_hexdigit())
            && !matches!(words[k - 1].as_str(), "=");
        if !is_color {
            out.push(w.clone());
        }
    }
    out
}

/// The order of core symbols in MEME matrices.
fn meme_order(core: &[u8]) -> Vec<u8> {
    let mut order = core.to_vec();
    // MEME orders matrix columns uppercase A-Z, lowercase a-z, digits, then
    // `*`, `-`, `.` (byte order), whatever order the block defines them in.
    order.sort_by_key(|&c| {
        let class = if c.is_ascii_uppercase() {
            0
        } else if c.is_ascii_lowercase() {
            1
        } else if c.is_ascii_digit() {
            2
        } else {
            3
        };
        (class, c)
    });
    order
}

/// Reads the matrix of the current motif: skips to its `letter-probability
/// matrix` line (stopping at the next `MOTIF`), then reads `w` rows (or rows
/// until a non-numeric line when `w=` is absent). Returns rows and `nsites`.
type Matrix = (Vec<Vec<f64>>, Option<f64>);

fn matrix(path: &Path, lines: &[&str], i: &mut usize, alength: usize) -> Result<Option<Matrix>> {
    while *i < lines.len() {
        let line = lines[*i].trim();
        if line.starts_with("MOTIF") {
            return Ok(None);
        }
        *i += 1;
        let Some(rest) = line.strip_prefix("letter-probability matrix") else {
            continue;
        };
        let header_line = *i;
        let fields = key_values(rest.trim_start_matches(':'));
        let get = |k: &str| fields.iter().find(|(key, _)| *key == k).map(|(_, v)| *v);
        let parse = |k: &str| -> Result<Option<f64>> {
            get(k)
                .map(|v| {
                    v.parse::<f64>()
                        .map_err(|_| bad(path, header_line, format!("{k}= {v:?} is not a number")))
                })
                .transpose()
        };
        if let Some(a) = parse("alength")?
            && a != alength as f64
        {
            return Err(bad(
                path,
                header_line,
                format!("alength= {a} but the alphabet has {alength} symbols"),
            ));
        }
        let w = get("w")
            .map(|w| {
                w.parse::<u32>()
                    .map(|n| n as usize)
                    .map_err(|_| bad(path, header_line, "w must be a non-negative 32-bit integer"))
            })
            .transpose()?;
        let nsites = parse("nsites")?.filter(|&n| n > 0.0);
        let mut rows = Vec::new();
        while *i < lines.len() && w.is_none_or(|w| rows.len() < w) {
            let line = lines[*i].trim();
            if line.is_empty() && w.is_some() {
                *i += 1;
                continue;
            }
            let Some(row) = numbers(line).filter(|r| !r.is_empty()) else {
                break;
            };
            if row.len() != alength {
                return Err(bad(
                    path,
                    *i + 1,
                    format!("{} values, the alphabet has {alength}", row.len()),
                ));
            }
            rows.push(row);
            *i += 1;
        }
        if let Some(w) = w
            && rows.len() != w
        {
            return Err(bad(
                path,
                header_line,
                format!("w= {w} but {} rows", rows.len()),
            ));
        }
        return Ok(Some((rows, nsites)));
    }
    Ok(None)
}

/// `key= value` pairs (MEME writes `alength= 4 w= 11 ...`).
fn key_values(text: &str) -> Vec<(&str, &str)> {
    let words: Vec<&str> = text.split_whitespace().collect();
    let mut out = Vec::new();
    let mut k = 0;
    while k < words.len() {
        let w = words[k];
        if let Some(key) = w.strip_suffix('=') {
            if let Some(v) = words.get(k + 1) {
                out.push((key, *v));
            }
            k += 2;
        } else if let Some((key, v)) = w.split_once('=') {
            out.push((key, v));
            k += 1;
        } else {
            k += 1;
        }
    }
    out
}
