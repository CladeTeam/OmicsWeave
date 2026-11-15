//! [`CountMatrix`] → a 10x-style directory `<name>/`:
//!
//! - `barcodes.tsv.gz`: one row key per line, in row order;
//! - `features.tsv.gz`: `id`, `name`, `type`, `chrom`, `start`, `end`
//!   (0-based half-open), one feature per line, in column order (the six
//!   columns of Cell Ranger ARC);
//! - the matrix, rows × features: `matrix.cpk` (CellPack, `u32`, one record
//!   per row, appended in order), or `matrix.mtx.gz` (Matrix Market,
//!   features × barcodes as Cell Ranger writes it, 1-based), or
//!   `matrix.csv.gz` (dense, a feature per line, a column per row key);
//!   further layers are `<layer>.cpk` / `<layer>.mtx.gz` / `<layer>.csv.gz`;
//! - with `table`, `counts.tsv`: featureCounts' layout (`Geneid`, `Chr`,
//!   `Start` (1-based), `End`, `Strand`, `Length`, a column per row key).
//!
//! The matrix is read once: every requested file is fed row by row.

use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use cellpack::CpkWriter;
use cellpack::format::cpk::{Dtype, SourceLayout};
use cellpack::scans::{
    CsrRef, DEFAULT_LANE_MIN, DEFAULT_MAX_LANES, DEFAULT_MODE_SHARE, DEFAULT_MU_INDICES,
    DEFAULT_MU_VALUES, EncodeOptions, Indices, LayoutChoice, Thresholds, ValueChoice,
};
use cellpack::store::CpkOptions;
use omicsweave_bgzf::{Writer, deflate_block_at};
use omicsweave_convert::human;
use omicsweave_core::{Error, Result};

use super::{GzText, TEXT_LEVEL, push_u64, write_gz};
use crate::accum::Delivery;
use crate::model::{CountMatrix, FeatureList, MatrixRow};
use crate::task::MatrixFormat;

/// What a matrix sink wrote.
pub struct Written {
    pub nonzero: u64,
    /// Paths relative to the output directory.
    pub files: Vec<String>,
    /// With `column_sums`, each layer's total per column.
    pub column_sums: Vec<Vec<u64>>,
}

fn cpk_error(e: cellpack::Error) -> Error {
    Error::invalid("cellpack", e.to_string())
}

/// CellPack's default encoding.
pub fn cpk_options() -> CpkOptions {
    CpkOptions::new(EncodeOptions {
        renumber: true,
        thresholds: Thresholds {
            indices: DEFAULT_MU_INDICES,
            values: DEFAULT_MU_VALUES,
        },
        max_lanes: DEFAULT_MAX_LANES,
        lane_min: DEFAULT_LANE_MIN,
        column_modes: Some(DEFAULT_MODE_SHARE),
        layout: LayoutChoice::Auto,
        values: ValueChoice::Auto,
    })
}

/// Rows staged for one `.cpk` layer.
struct CpkLayer {
    writer: CpkWriter,
    indptr: Vec<u64>,
    indices: Vec<u32>,
    data: Vec<u8>,
    n_cols: u32,
}

impl CpkLayer {
    fn row(&mut self, cols: &[u32], values: &[u32]) {
        for (&c, &v) in cols.iter().zip(values) {
            if v != 0 {
                self.indices.push(c);
                self.data.extend_from_slice(&v.to_le_bytes());
            }
        }
        self.indptr.push(self.indices.len() as u64);
    }

    fn flush(&mut self) -> Result<()> {
        if self.indptr.len() > 1 {
            self.writer
                .append(CsrRef {
                    n_genes: self.n_cols,
                    value_width: 4,
                    indptr: &self.indptr,
                    indices: Indices::U32(&self.indices),
                    data: &self.data,
                })
                .map_err(cpk_error)?;
        }
        self.indptr.clear();
        self.indptr.push(0);
        self.indices.clear();
        self.data.clear();
        Ok(())
    }

    fn staged(&self) -> usize {
        self.indptr.len() - 1
    }
}

/// One `.mtx.gz` layer: entries go to a temporary BGZF body, which follows
/// the header once the count of entries is known.
struct MtxLayer {
    body: Writer<BufWriter<File>>,
    path: PathBuf,
    nnz: u64,
    line: Vec<u8>,
}

/// Writes `m` as `<name>/` under the delivery's directory.
pub fn write(m: CountMatrix, name: &str, table: bool, d: &mut Delivery<'_>) -> Result<Written> {
    write_with(m, name, table, false, d)
}

/// [`write`], also summing every layer per column when `column_sums`.
pub fn write_with(
    mut m: CountMatrix,
    name: &str,
    table: bool,
    column_sums: bool,
    d: &mut Delivery<'_>,
) -> Result<Written> {
    if table && m.features.table().is_none() {
        return Err(Error::invalid(
            "matrix",
            format!("{name}: a featureCounts table needs feature coordinates"),
        ));
    }
    let n_rows = m.rows.len();
    let n_cols = m.features.len();
    let mut files = Vec::new();
    let rel = |f: &str| format!("{name}/{f}");

    let mut text = Vec::new();
    for b in &m.rows.names {
        text.extend_from_slice(b);
        text.push(b'\n');
    }
    let p = d.path(&rel("barcodes.tsv.gz"))?;
    write_gz(d, &p, "barcodes", &text)?;
    files.push(rel("barcodes.tsv.gz"));
    let p = d.path(&rel("features.tsv.gz"))?;
    features_file(&*m.features, &p, d)?;
    files.push(rel("features.tsv.gz"));

    let layers = m.layers.clone();
    let named = m.named_layers;
    let file_of = |l: usize, ext: &str| {
        if l == 0 && !named {
            format!("matrix.{ext}")
        } else {
            format!("{}.{ext}", layers[l])
        }
    };
    let dense_bytes = (n_rows as u64)
        .saturating_mul(n_cols as u64)
        .saturating_mul(4 * layers.len() as u64);
    let dense_needed = table || d.matrix == MatrixFormat::Csv;
    if dense_needed && dense_bytes > d.buffer.max(64 << 20) {
        return Err(Error::invalid(
            "matrix",
            format!(
                "{name}: a dense {n_rows} × {n_cols} table needs {}, more than the budget allows; \
                 write cpk or mtx",
                human(dense_bytes)
            ),
        ));
    }
    let mut dense: Vec<Vec<u32>> = if dense_needed {
        (0..layers.len())
            .map(|_| vec![0u32; n_rows * n_cols])
            .collect()
    } else {
        Vec::new()
    };

    let mut cpk: Vec<CpkLayer> = Vec::new();
    let mut mtx: Vec<MtxLayer> = Vec::new();
    match d.matrix {
        MatrixFormat::Cpk => {
            if n_cols == 0 {
                return Err(Error::invalid(
                    "matrix",
                    format!("{name}: no features to count over"),
                ));
            }
            for l in 0..layers.len() {
                let path = d.path(&rel(&file_of(l, "cpk")))?;
                drop(d.outputs.create(&path, "cpk")?);
                let tmp = d
                    .outputs
                    .temporary(&path)
                    .expect("created above")
                    .to_owned();
                let mut writer = CpkWriter::create(
                    &tmp,
                    n_cols as u32,
                    Dtype::U32,
                    SourceLayout::Csr,
                    cpk_options(),
                )
                .map_err(cpk_error)?;
                writer.set_threads(d.pool.threads());
                cpk.push(CpkLayer {
                    writer,
                    indptr: vec![0],
                    indices: Vec::new(),
                    data: Vec::new(),
                    n_cols: n_cols as u32,
                });
                files.push(rel(&file_of(l, "cpk")));
            }
        }
        MatrixFormat::Mtx => {
            for l in 0..layers.len() {
                let path = d.tmp.join(format!(
                    ".omicsweave-mtx-{}-{}-{l}",
                    std::process::id(),
                    name.replace('/', "_")
                ));
                let f = File::create(&path)?;
                mtx.push(MtxLayer {
                    body: Writer::new(
                        BufWriter::with_capacity(1 << 20, f),
                        TEXT_LEVEL,
                        Box::new(d.pool.ordered(4 * d.pool.threads())),
                    )?,
                    path,
                    nnz: 0,
                    line: Vec::new(),
                });
            }
        }
        MatrixFormat::Csv => {}
    }

    let mut sums: Vec<Vec<u64>> = if column_sums {
        (0..layers.len()).map(|_| vec![0u64; n_cols]).collect()
    } else {
        Vec::new()
    };
    let mut row = MatrixRow::default();
    let mut next_row = 0usize;
    let mut nonzero = 0u64;
    while m.next_row(&mut row)? {
        let r = row.row as usize;
        let n = row.cols.len();
        nonzero += n as u64;
        for (l, s) in sums.iter_mut().enumerate() {
            for (i, &c) in row.cols.iter().enumerate() {
                s[c as usize] += row.values[l * n + i] as u64;
            }
        }
        // Rows without any nonzero cell are empty records.
        for layer in &mut cpk {
            for _ in next_row..r {
                layer.row(&[], &[]);
                if layer.staged() >= 4096 {
                    layer.flush()?;
                }
            }
        }
        for (l, layer) in cpk.iter_mut().enumerate() {
            layer.row(&row.cols, &row.values[l * n..(l + 1) * n]);
            if layer.staged() >= 4096 || layer.indices.len() >= 1 << 22 {
                layer.flush()?;
            }
        }
        next_row = r + 1;
        for (l, layer) in mtx.iter_mut().enumerate() {
            for (i, &c) in row.cols.iter().enumerate() {
                let v = row.values[l * n + i];
                if v == 0 {
                    continue;
                }
                layer.line.clear();
                push_u64(&mut layer.line, c as u64 + 1);
                layer.line.push(b' ');
                push_u64(&mut layer.line, r as u64 + 1);
                layer.line.push(b' ');
                push_u64(&mut layer.line, v as u64);
                layer.line.push(b'\n');
                layer.body.write(&layer.line)?;
                layer.nnz += 1;
            }
        }
        for (l, t) in dense.iter_mut().enumerate() {
            for (i, &c) in row.cols.iter().enumerate() {
                t[c as usize * n_rows + r] = row.values[l * n + i];
            }
        }
    }
    for mut layer in cpk {
        for _ in next_row..n_rows {
            layer.row(&[], &[]);
            if layer.staged() >= 4096 {
                layer.flush()?;
            }
        }
        layer.flush()?;
        layer.writer.finish().map_err(cpk_error)?;
    }
    for (l, layer) in mtx.into_iter().enumerate() {
        let file = file_of(l, "mtx.gz");
        let path = d.path(&rel(&file))?;
        mtx_file(layer, &path, n_rows, n_cols, d)?;
        files.push(rel(&file));
    }
    if d.matrix == MatrixFormat::Csv {
        for (l, t) in dense.iter().enumerate() {
            let file = file_of(l, "csv.gz");
            let path = d.path(&rel(&file))?;
            dense_csv(t, &m, &path, d)?;
            files.push(rel(&file));
        }
    }
    if table {
        let path = d.path(&rel("counts.tsv"))?;
        counts_tsv(&dense[0], &m, &path, d)?;
        files.push(rel("counts.tsv"));
    }
    Ok(Written {
        nonzero,
        files,
        column_sums: sums,
    })
}

/// `features.tsv.gz`, a feature per line, streamed in chunks.
fn features_file(f: &dyn FeatureList, path: &Path, d: &mut Delivery<'_>) -> Result<()> {
    let mut g = GzText::create(d, path, "features")?;
    let mut text = Vec::with_capacity(1 << 16);
    for i in 0..f.len() {
        f.line(i, &mut text);
        text.push(b'\n');
        if text.len() >= 1 << 16 {
            g.write(&text)?;
            text.clear();
        }
    }
    g.write(&text)?;
    g.close()
}

/// The header block, then the body's blocks (BGZF members concatenate
/// into one gzip stream).
fn mtx_file(
    layer: MtxLayer,
    path: &Path,
    n_rows: usize,
    n_cols: usize,
    d: &mut Delivery<'_>,
) -> Result<()> {
    let MtxLayer {
        body,
        path: tmp,
        nnz,
        ..
    } = layer;
    let done = body.finish()?;
    done.out.into_inner().map_err(|e| e.into_error())?;
    let header =
        format!("%%MatrixMarket matrix coordinate integer general\n%\n{n_cols} {n_rows} {nnz}\n");
    let mut head = Vec::new();
    deflate_block_at(TEXT_LEVEL, header.as_bytes(), &mut head)?;
    let mut out = BufWriter::with_capacity(1 << 20, d.outputs.create(path, "mtx")?);
    out.write_all(&head)?;
    let mut src = File::open(&tmp)?;
    let mut buf = vec![0u8; 1 << 20];
    loop {
        let n = src.read(&mut buf)?;
        if n == 0 {
            break;
        }
        out.write_all(&buf[..n])?;
    }
    out.flush()?;
    drop(src);
    let _ = std::fs::remove_file(&tmp);
    Ok(())
}

/// Dense `feature × row` CSV.
fn dense_csv(t: &[u32], m: &CountMatrix, path: &Path, d: &mut Delivery<'_>) -> Result<()> {
    let n_rows = m.rows.len();
    let mut g = GzText::create(d, path, "csv")?;
    let mut line = b"feature".to_vec();
    for r in &m.rows.names {
        line.push(b',');
        line.extend_from_slice(r);
    }
    line.push(b'\n');
    g.write(&line)?;
    for c in 0..m.features.len() {
        line.clear();
        m.features.id(c, &mut line);
        for &v in &t[c * n_rows..(c + 1) * n_rows] {
            line.push(b',');
            push_u64(&mut line, v as u64);
        }
        line.push(b'\n');
        g.write(&line)?;
    }
    g.close()
}

/// featureCounts' table: annotation columns, then a count column per row
/// key.
fn counts_tsv(t: &[u32], m: &CountMatrix, path: &Path, d: &mut Delivery<'_>) -> Result<()> {
    let n_rows = m.rows.len();
    let f = m.features.table().expect("checked by write_with");
    let mut out = BufWriter::with_capacity(1 << 20, d.outputs.create(path, "tsv")?);
    let mut line = b"Geneid\tChr\tStart\tEnd\tStrand\tLength".to_vec();
    for r in &m.rows.names {
        line.push(b'\t');
        line.extend_from_slice(r);
    }
    line.push(b'\n');
    out.write_all(&line)?;
    for c in 0..f.len() {
        line.clear();
        line.extend_from_slice(f.id[c].as_bytes());
        line.push(b'\t');
        line.extend_from_slice(f.chrom[c].as_bytes());
        line.push(b'\t');
        push_u64(&mut line, f.start[c] as u64 + 1);
        line.push(b'\t');
        push_u64(&mut line, f.end[c] as u64);
        line.push(b'\t');
        line.push(f.strand[c] as u8);
        line.push(b'\t');
        push_u64(&mut line, f.length[c]);
        for &v in &t[c * n_rows..(c + 1) * n_rows] {
            line.push(b'\t');
            push_u64(&mut line, v as u64);
        }
        line.push(b'\n');
        out.write_all(&line)?;
    }
    out.flush()?;
    Ok(())
}
