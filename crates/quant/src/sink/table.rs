//! Small tables: per-barcode metrics as CSV, histograms, mismatch tables,
//! M-bias and per-variant counts as TSV.

use std::io::{BufWriter, Write};
use std::path::Path;

use omicsweave_core::Result;

use super::{push_u64, write_gz};
use crate::accum::Delivery;
use crate::model::{BarcodeMetrics, Histogram, Mbias, MismatchTable, VariantTable};

/// `barcode,<column>,...`, a barcode per line in row order.
pub fn write_metrics(m: &BarcodeMetrics, path: &Path, d: &mut Delivery<'_>) -> Result<()> {
    let mut out = BufWriter::with_capacity(1 << 20, d.outputs.create(path, "csv")?);
    let mut line = b"barcode".to_vec();
    for (name, _) in &m.columns {
        line.push(b',');
        line.extend_from_slice(name.as_bytes());
    }
    line.push(b'\n');
    out.write_all(&line)?;
    for (r, b) in m.barcodes.names.iter().enumerate() {
        line.clear();
        line.extend_from_slice(b);
        for (_, v) in &m.columns {
            line.push(b',');
            push_u64(&mut line, v[r]);
        }
        line.push(b'\n');
        out.write_all(&line)?;
    }
    out.flush()?;
    Ok(())
}

/// `<x axis>\t<y axis>\treads`, a nonzero cell per line, ascending (BGZF).
pub fn write_histogram(h: &Histogram, path: &Path, d: &mut Delivery<'_>) -> Result<()> {
    let mut text = format!("{}\t{}\treads\n", h.axes[0], h.axes[1]).into_bytes();
    for &(x, y, n) in &h.cells {
        push_u64(&mut text, x as u64);
        text.push(b'\t');
        push_u64(&mut text, y as u64);
        text.push(b'\t');
        push_u64(&mut text, n);
        text.push(b'\n');
    }
    write_gz(d, path, "histogram", &text)
}

/// `mate\tcycle\treference\tA\tC\tG\tT`: per mate (1, 2), cycle
/// (0-based) and transcript-sense reference base, the read bases over it.
pub fn write_mismatches(t: &MismatchTable, path: &Path, d: &mut Delivery<'_>) -> Result<()> {
    let mut out = BufWriter::with_capacity(1 << 20, d.outputs.create(path, "tsv")?);
    out.write_all(b"mate\tcycle\treference\tA\tC\tG\tT\n")?;
    let mut line = Vec::new();
    for (m, cycles) in t.cycles.iter().enumerate() {
        for (c, counts) in cycles.iter().enumerate() {
            for (r, row) in counts.iter().enumerate() {
                if row.iter().all(|&v| v == 0) {
                    continue;
                }
                line.clear();
                push_u64(&mut line, m as u64 + 1);
                line.push(b'\t');
                push_u64(&mut line, c as u64);
                line.push(b'\t');
                line.push(b"ACGT"[r]);
                for &v in row {
                    line.push(b'\t');
                    push_u64(&mut line, v);
                }
                line.push(b'\n');
                out.write_all(&line)?;
            }
        }
    }
    out.flush()?;
    Ok(())
}

/// `context\tmate\tposition\tmethylated\tunmethylated`: per context,
/// mate (1, 2) and 1-based position in the read as sequenced.
pub fn write_mbias(m: &Mbias, path: &Path, d: &mut Delivery<'_>) -> Result<()> {
    let mut out = BufWriter::with_capacity(1 << 20, d.outputs.create(path, "tsv")?);
    out.write_all(b"context\tmate\tposition\tmethylated\tunmethylated\n")?;
    let mut line = Vec::new();
    for (c, name) in crate::accum::methylation::CONTEXTS.iter().enumerate() {
        for mate in 0..2 {
            for (p, v) in m.counts[mate * 3 + c].iter().enumerate() {
                line.clear();
                line.extend_from_slice(name.as_bytes());
                line.push(b'\t');
                push_u64(&mut line, mate as u64 + 1);
                line.push(b'\t');
                push_u64(&mut line, p as u64 + 1);
                line.push(b'\t');
                push_u64(&mut line, v[0]);
                line.push(b'\t');
                push_u64(&mut line, v[1]);
                line.push(b'\n');
                out.write_all(&line)?;
            }
        }
    }
    out.flush()?;
    Ok(())
}

/// `chrom\tpos\tref\talt\tref_count\talt_count\tother_count\talt_fraction`
/// (BGZF): a site per line in panel order, 1-based; the fraction is
/// `alt / (ref + alt)`, empty without either.
pub fn write_variants(t: &VariantTable, path: &Path, d: &mut Delivery<'_>) -> Result<()> {
    let mut g = super::GzText::create_unindexed(d, path, "variants")?;
    g.write(b"chrom\tpos\tref\talt\tref_count\talt_count\tother_count\talt_fraction\n")?;
    let mut line = Vec::with_capacity(1 << 16);
    for (s, c) in t.sites.iter().zip(&t.counts) {
        line.extend_from_slice(t.contigs[s.contig as usize].as_bytes());
        line.push(b'\t');
        push_u64(&mut line, s.pos as u64 + 1);
        line.push(b'\t');
        line.push(s.reference);
        line.push(b'\t');
        line.push(s.alt);
        for v in c {
            line.push(b'\t');
            push_u64(&mut line, *v);
        }
        line.push(b'\t');
        let depth = c[0] + c[1];
        if depth > 0 {
            line.extend_from_slice(format!("{:.6}", c[1] as f64 / depth as f64).as_bytes());
        }
        line.push(b'\n');
        if line.len() >= 1 << 16 {
            g.write(&line)?;
            line.clear();
        }
    }
    g.write(&line)?;
    g.close()
}
