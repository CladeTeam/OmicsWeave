//! [`FragmentTable`] → 10x `fragments.tsv.gz` (BGZF; `chrom`, `start`,
//! `end`, `barcode`, `count`, 0-based half-open) and its tabix index
//! `fragments.tsv.gz.tbi`, built as `tabix -p bed` builds it; the indexed
//! BED-like writer ([`IndexedBed`]) is shared with the per-cytosine
//! methylation table.

use std::path::{Path, PathBuf};

use omicsweave_bgzf::{BinningIndex, IndexFormat, Writer};
use omicsweave_core::exec::Inline;
use omicsweave_core::{Error, Result};

use super::{GzText, push_u64};
use crate::accum::Delivery;
use crate::model::FragmentTable;

/// `tbx_conf_bed`: generic format with UCSC (0-based) coordinates,
/// sequence in column 1, begin 2, end 3, `#` comments, no skipped lines.
const TBX_UCSC: i32 = 0x10000;
const BED_CONF: [i32; 6] = [TBX_UCSC, 1, 2, 3, b'#' as i32, 0];

/// Writes the fragments and their index; returns the index's path.
pub fn write(table: FragmentTable, path: &Path, d: &mut Delivery<'_>) -> Result<PathBuf> {
    let FragmentTable {
        chroms,
        barcodes,
        fragments,
    } = table;
    let mut w = IndexedBed::create(d, path, "fragments", &chroms)?;
    let mut line = Vec::with_capacity(128);
    fragments.drain(|f| {
        line.clear();
        line.extend_from_slice(chroms[f.tid as usize].0.as_bytes());
        line.push(b'\t');
        push_u64(&mut line, f.start as u64);
        line.push(b'\t');
        push_u64(&mut line, f.end as u64);
        line.push(b'\t');
        line.extend_from_slice(&barcodes.names[f.row as usize]);
        line.push(b'\t');
        push_u64(&mut line, f.count as u64);
        line.push(b'\n');
        w.push(f.tid, f.start, f.end, &line)
    })?;
    w.finish(path, d)
}

/// BGZF text of BED-like lines in coordinate order, indexed as `tabix -p
/// bed` indexes it.
pub(crate) struct IndexedBed<'c> {
    g: GzText,
    chroms: &'c [(String, u32)],
    /// Tabix numbers sequences in order of appearance; the index is built
    /// over every reference and trimmed to those that appeared.
    seen: Vec<Option<u32>>,
    names: Vec<&'c str>,
    index: BinningIndex,
}

impl<'c> IndexedBed<'c> {
    pub fn create(
        d: &mut Delivery<'_>,
        path: &Path,
        kind: &'static str,
        chroms: &'c [(String, u32)],
    ) -> Result<IndexedBed<'c>> {
        let g = GzText::create(d, path, kind)?;
        let index = BinningIndex::try_new(IndexFormat::Bai, chroms.len(), g.tell())?;
        Ok(IndexedBed {
            g,
            chroms,
            seen: vec![None; chroms.len()],
            names: Vec::new(),
            index,
        })
    }

    /// Writes `line` (with its newline), an interval `[start, end)` of
    /// reference `tid`; lines come sorted by reference and start.
    pub fn push(&mut self, tid: u32, start: u32, end: u32, line: &[u8]) -> Result<()> {
        let t = tid as usize;
        let id = match self.seen.get(t).copied() {
            Some(Some(x)) => x,
            Some(None) => {
                let x = self.names.len() as u32;
                self.names.push(&self.chroms[t].0);
                self.seen[t] = Some(x);
                x
            }
            None => {
                return Err(Error::invalid("tabix", "a line on an unknown reference"));
            }
        };
        self.g.write(line)?;
        self.index
            .push(id as i32, start as i64, end as i64, self.g.tell(), true)
    }

    /// Finishes the text and writes `<path>.tbi`; returns its path.
    pub fn finish(self, path: &Path, d: &mut Delivery<'_>) -> Result<PathBuf> {
        let IndexedBed {
            mut g,
            chroms,
            names,
            index,
            ..
        } = self;
        g.flush_block()?;
        let end = g.tell();
        let done = g.finish()?;
        done.out.into_inner().map_err(|e| e.into_error())?;
        let bai = index.finish(end, &done.blocks)?;
        let tbi = tabix(&bai, &names, chroms.len())?;
        let mut tbi_path = path.as_os_str().to_owned();
        tbi_path.push(".tbi");
        let tbi_path = PathBuf::from(tbi_path);
        d.outputs.write(&tbi_path, "tbi", &tbi)?;
        Ok(tbi_path)
    }
}

/// A tabix index from the BAI body of the same bins over `n_refs`
/// references, of which the first `names.len()` appeared: magic, reference
/// count, the configuration and names, then the BAI's references (the
/// empty trailing ones dropped) and unplaced count; BGZF-compressed.
fn tabix(bai: &[u8], names: &[&str], n_refs: usize) -> Result<Vec<u8>> {
    // Each reference without records is a zero bin count and a zero
    // linear-index length.
    let empty = 8 * (n_refs - names.len());
    let refs_end = bai.len().checked_sub(8 + empty);
    let Some(refs_end) = refs_end.filter(|&e| e >= 8 && &bai[..4] == b"BAI\x01") else {
        return Err(Error::invalid("tabix", "unexpected index layout"));
    };
    if bai[refs_end..bai.len() - 8].iter().any(|&b| b != 0) {
        return Err(Error::invalid("tabix", "references out of order"));
    }
    let mut body = Vec::with_capacity(bai.len() + 64);
    body.extend_from_slice(b"TBI\x01");
    body.extend_from_slice(&(names.len() as i32).to_le_bytes());
    for v in BED_CONF {
        body.extend_from_slice(&v.to_le_bytes());
    }
    let l_nm: usize = names.iter().map(|n| n.len() + 1).sum();
    body.extend_from_slice(&(l_nm as i32).to_le_bytes());
    for n in names {
        body.extend_from_slice(n.as_bytes());
        body.push(0);
    }
    body.extend_from_slice(&bai[8..refs_end]);
    body.extend_from_slice(&bai[bai.len() - 8..]);
    let mut w = Writer::new(
        Vec::new(),
        omicsweave_bgzf::DEFAULT_LEVEL,
        Box::new(Inline::new()),
    )?;
    w.write(&body)?;
    Ok(w.finish()?.out)
}
