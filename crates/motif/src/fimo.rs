//! FIMO's `--text` output for [`HitEntry`] rows.

use std::io::{self, Write};

use omicsweave_core::{Strand, Strings};
use omicsweave_ffm::Alphabet;

use crate::reduce::HitEntry;

/// Names and bytes the text output prints beside each hit.
pub struct FimoText<'a> {
    /// The scanned alphabet; reverse-strand matches print its symbols.
    pub alphabet: &'a Alphabet,
    /// Per motif: its ID and its alternate ID (may be empty).
    pub motif_ids: &'a [&'a str],
    pub motif_alt_ids: &'a [&'a str],
    /// Per sequence: its name and the scanned bytes.
    pub sequence_names: Strings<'a>,
    pub sequences: Strings<'a>,
    /// Per sequence: added to positions (0 when absent), e.g. a region start.
    pub origins: Option<&'a [u64]>,
}

/// Writes `hits` as FIMO 5.5.9 `--text --no-qvalue` does: a header before
/// the first row, 1-based inclusive coordinates, `%g` scores, `%.3g`
/// p-values, an empty q-value column and the matched bases: as given on the
/// forward strand, as the alphabet's complement symbols on the reverse.
pub fn write_text<W: Write>(mut out: W, hits: &[HitEntry], text: &FimoText<'_>) -> io::Result<()> {
    if hits.is_empty() {
        return Ok(());
    }
    out.write_all(
        b"motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n",
    )?;
    let lut = text.alphabet.lut();
    let symbols = text.alphabet.symbols();
    let complement = text.alphabet.complement();
    let mut matched = Vec::new();
    for h in hits {
        let m = h.motif as usize;
        let s = h.sequence as usize;
        let origin = text.origins.map_or(0, |o| o[s]);
        let window = &text.sequences.get(s)[h.start as usize..h.end as usize];
        matched.clear();
        let strand = match (h.strand, complement) {
            (Strand::Reverse, Some(comp)) => {
                matched.extend(window.iter().rev().map(|&b| match lut[b as usize] {
                    0xFF => b,
                    c => symbols[comp[c as usize] as usize],
                }));
                '-'
            }
            _ => {
                matched.extend_from_slice(window);
                '+'
            }
        };
        out.write_all(text.motif_ids[m].as_bytes())?;
        out.write_all(b"\t")?;
        out.write_all(text.motif_alt_ids.get(m).copied().unwrap_or("").as_bytes())?;
        out.write_all(b"\t")?;
        out.write_all(text.sequence_names.get(s))?;
        write!(
            out,
            "\t{}\t{}\t{}\t{}\t{}\t\t",
            origin + h.start as u64 + 1,
            origin + h.end as u64,
            strand,
            format_g(h.score, 6),
            format_g(h.pvalue, 3),
        )?;
        out.write_all(&matched)?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

/// C's `%.{precision}g`.
pub fn format_g(x: f64, precision: usize) -> String {
    let p = precision.max(1);
    if x.is_nan() {
        return if x.is_sign_negative() { "-nan" } else { "nan" }.into();
    }
    if x.is_infinite() {
        return if x < 0.0 { "-inf" } else { "inf" }.into();
    }
    if x == 0.0 {
        return if x.is_sign_negative() { "-0" } else { "0" }.into();
    }
    let sci = format!("{:.*e}", p - 1, x);
    let (mantissa, exp) = sci.split_once('e').expect("exponent");
    let exp: i32 = exp.parse().expect("exponent digits");
    if exp < -4 || exp >= p as i32 {
        let mantissa = strip_zeros(mantissa);
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{mantissa}e{sign}{:02}", exp.abs())
    } else {
        let fixed = format!("{:.*}", (p as i32 - 1 - exp) as usize, x);
        strip_zeros(&fixed).to_string()
    }
}

/// C's `%.{precision}e`.
pub fn format_e(x: f64, precision: usize) -> String {
    if !x.is_finite() {
        return format_g(x, 1);
    }
    let sci = format!("{:.*e}", precision, x);
    let (mantissa, exp) = sci.split_once('e').expect("exponent");
    let exp: i32 = exp.parse().expect("exponent digits");
    let sign = if exp < 0 { '-' } else { '+' };
    format!("{mantissa}e{sign}{:02}", exp.abs())
}

fn strip_zeros(s: &str) -> &str {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g_formatting_matches_c() {
        assert_eq!(format_g(17.5455, 6), "17.5455");
        assert_eq!(format_g(-3.12, 6), "-3.12");
        assert_eq!(format_g(100000.0, 6), "100000");
        assert_eq!(format_g(1234567.0, 6), "1.23457e+06");
        assert_eq!(format_g(1.23e-5, 3), "1.23e-05");
        assert_eq!(format_g(0.0001, 3), "0.0001");
        assert_eq!(format_g(9.996e-5, 3), "0.0001");
        assert_eq!(format_g(2.5e-10, 3), "2.5e-10");
        assert_eq!(format_g(0.5, 3), "0.5");
        assert_eq!(format_e(0.00105, 2), "1.05e-03");
        assert_eq!(format_e(0.0, 2), "0.00e+00");
        assert_eq!(format_e(123456.0, 2), "1.23e+05");
    }
}
