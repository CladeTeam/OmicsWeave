//! Decode the external oracle output without consensus implementation helpers.

/// The sequence lines of a FASTA file, joined.
pub fn fasta_sequence(text: &str) -> Vec<u8> {
    text.lines()
        .filter(|l| !l.starts_with('>'))
        .flat_map(|l| l.trim_end().bytes())
        .collect()
}
