//! [`FragmentTable`] → 10x `fragments.tsv.gz` (BGZF; `chrom`, `start`,
//! `end`, `barcode`, `count`, 0-based half-open) and its tabix index
//! `fragments.tsv.gz.tbi` (CSI for references longer than 2^29 bases),
//! built as `tabix -p bed` builds it; the indexed
//! BED-like writer ([`IndexedBed`]) is shared with the per-cytosine
//! methylation table.

use std::io::Read;
use std::path::{Path, PathBuf};

use omicsweave_bgzf::{BinningIndex, IndexFormat, Writer, try_csi_levels};
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
    format: IndexFormat,
}

impl<'c> IndexedBed<'c> {
    pub fn create(
        d: &mut Delivery<'_>,
        path: &Path,
        kind: &'static str,
        chroms: &'c [(String, u32)],
    ) -> Result<IndexedBed<'c>> {
        i32::try_from(chroms.len()).map_err(|_| Error::invalid("tabix", "too many references"))?;
        let g = GzText::create(d, path, kind)?;
        let max_len = chroms.iter().map(|c| u64::from(c.1)).max().unwrap_or(0);
        let format = if max_len <= 1 << 29 {
            IndexFormat::Bai
        } else {
            let (min_shift, levels) = try_csi_levels(max_len, 14)?;
            IndexFormat::Csi { min_shift, levels }
        };
        let index = BinningIndex::try_new(format, chroms.len(), g.tell())?;
        Ok(IndexedBed {
            g,
            chroms,
            seen: vec![None; chroms.len()],
            names: Vec::new(),
            index,
            format,
        })
    }

    /// Writes `line` (with its newline), an interval `[start, end)` of
    /// reference `tid`; lines come sorted by reference and start.
    pub fn push(&mut self, tid: u32, start: u32, end: u32, line: &[u8]) -> Result<()> {
        let t = tid as usize;
        let Some((_, len)) = self.chroms.get(t) else {
            return Err(Error::invalid("tabix", "a line on an unknown reference"));
        };
        if start >= end || end > *len {
            return Err(Error::invalid(
                "tabix",
                "an interval is empty or outside its reference",
            ));
        }
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

    /// Finishes the text and writes `<path>.tbi` or `<path>.csi`.
    pub fn finish(self, path: &Path, d: &mut Delivery<'_>) -> Result<PathBuf> {
        let IndexedBed {
            mut g,
            chroms,
            names,
            index,
            format,
            ..
        } = self;
        g.flush_block()?;
        let end = g.tell();
        let done = g.finish()?;
        done.out.into_inner().map_err(|e| e.into_error())?;
        let raw = index.finish(end, &done.blocks)?;
        let kind = match format {
            IndexFormat::Bai => "tbi",
            IndexFormat::Csi { .. } => "csi",
        };
        let bytes = tabix(&raw, &names, chroms.len(), format)?;
        let mut index_path = path.as_os_str().to_owned();
        index_path.push(format!(".{kind}"));
        let index_path = PathBuf::from(index_path);
        d.outputs.write(&index_path, kind, &bytes)?;
        Ok(index_path)
    }
}

/// Adds the BED configuration and observed sequence names to the binning
/// index, dropping unused trailing references. TBI stores the configuration
/// in its header; CSI stores it in its auxiliary section. Both are BGZF.
fn tabix(raw: &[u8], names: &[&str], n_refs: usize, format: IndexFormat) -> Result<Vec<u8>> {
    let n_names = i32::try_from(names.len())
        .map_err(|_| Error::invalid("tabix", "too many sequence names"))?;
    let missing = n_refs
        .checked_sub(names.len())
        .ok_or_else(|| Error::invalid("tabix", "more names than references"))?;
    let n_refs =
        i32::try_from(n_refs).map_err(|_| Error::invalid("tabix", "too many references"))?;
    let l_nm = names
        .iter()
        .try_fold(0usize, |len, name| {
            len.checked_add(name.len())?.checked_add(1)
        })
        .ok_or_else(|| Error::invalid("tabix", "sequence names exceed index capacity"))?;
    let aux_len = (BED_CONF.len() * 4 + 4)
        .checked_add(l_nm)
        .and_then(|len| i32::try_from(len).ok())
        .ok_or_else(|| Error::invalid("tabix", "sequence names exceed index capacity"))?;
    let l_nm = i32::try_from(l_nm)
        .map_err(|_| Error::invalid("tabix", "sequence names exceed index capacity"))?;
    let mut aux = Vec::with_capacity(aux_len as usize);
    for v in BED_CONF {
        aux.extend_from_slice(&v.to_le_bytes());
    }
    aux.extend_from_slice(&l_nm.to_le_bytes());
    for n in names {
        aux.extend_from_slice(n.as_bytes());
        aux.push(0);
    }
    let mut csi = Vec::new();
    let (index, header, empty_ref) = match format {
        IndexFormat::Bai => (raw, 8, 8),
        IndexFormat::Csi { .. } => {
            flate2::read::MultiGzDecoder::new(raw).read_to_end(&mut csi)?;
            (&csi[..], 20, 4)
        }
    };
    // Empty references contain only a zero bin count, plus a zero
    // linear-index length in BAI.
    let empty = missing
        .checked_mul(empty_ref)
        .and_then(|n| n.checked_add(8))
        .ok_or_else(|| Error::invalid("tabix", "trailing references exceed index capacity"))?;
    let refs_end = index.len().checked_sub(empty);
    let magic = if header == 8 { b"BAI\x01" } else { b"CSI\x01" };
    let Some(refs_end) = refs_end.filter(|&e| e >= header && &index[..4] == magic) else {
        return Err(Error::invalid("tabix", "unexpected index layout"));
    };
    let word = |at: usize| u32::from_le_bytes(index[at..at + 4].try_into().expect("index word"));
    if word(header - 4) != n_refs as u32 {
        return Err(Error::invalid("tabix", "index reference count differs"));
    }
    if let IndexFormat::Csi { min_shift, levels } = format
        && (word(4) != min_shift || word(8) != levels || word(12) != 0)
    {
        return Err(Error::invalid(
            "tabix",
            "unexpected CSI configuration or auxiliary data",
        ));
    }
    if index[refs_end..index.len() - 8].iter().any(|&b| b != 0) {
        return Err(Error::invalid("tabix", "references out of order"));
    }
    let mut w = Writer::new_streaming(
        Vec::new(),
        omicsweave_bgzf::DEFAULT_LEVEL,
        Box::new(Inline::new()),
    )?;
    if header == 8 {
        w.write(b"TBI\x01")?;
        w.write(&n_names.to_le_bytes())?;
        w.write(&aux)?;
    } else {
        w.write(&index[..12])?;
        w.write(&aux_len.to_le_bytes())?;
        w.write(&aux)?;
        w.write(&n_names.to_le_bytes())?;
    }
    w.write(&index[header..refs_end])?;
    w.write(&index[index.len() - 8..])?;
    w.finish_streaming()
}

#[cfg(test)]
mod tests {
    use super::*;
    use htslib_sys as h;
    use omicsweave_convert::stage::Outputs;
    use omicsweave_core::exec::Pool;
    use std::ffi::CString;

    #[test]
    fn tabix_conversion_checks_header_counts_and_csi_auxiliary_layout() {
        let mut bai = b"BAI\x01".to_vec();
        bai.extend_from_slice(&0u32.to_le_bytes());
        bai.extend_from_slice(&0u64.to_le_bytes());
        assert!(tabix(&bai, &[], 0, IndexFormat::Bai).is_ok());
        assert!(tabix(&bai[..7], &[], 0, IndexFormat::Bai).is_err());
        assert!(tabix(&bai, &["extra"], 0, IndexFormat::Bai).is_err());
        assert!(tabix(&bai, &[], usize::MAX, IndexFormat::Bai).is_err());
        assert!(tabix(&bai, &[], 1, IndexFormat::Bai).is_err());
        bai[4..8].copy_from_slice(&1u32.to_le_bytes());
        assert!(tabix(&bai, &[], 0, IndexFormat::Bai).is_err());

        for (min_shift, levels, aux_len, n_refs) in [
            (15u32, 5u32, 0u32, 0u32),
            (14, 6, 0, 0),
            (14, 5, 4, 0),
            (14, 5, 0, 1),
        ] {
            let mut body = b"CSI\x01".to_vec();
            for word in [min_shift, levels, aux_len, n_refs] {
                body.extend_from_slice(&word.to_le_bytes());
            }
            body.extend_from_slice(&0u64.to_le_bytes());
            let mut writer = Writer::new_streaming(
                Vec::new(),
                omicsweave_bgzf::DEFAULT_LEVEL,
                Box::new(Inline::new()),
            )
            .unwrap();
            writer.write(&body).unwrap();
            let raw = writer.finish_streaming().unwrap();
            assert!(
                tabix(
                    &raw,
                    &[],
                    0,
                    IndexFormat::Csi {
                        min_shift: 14,
                        levels: 5
                    }
                )
                .is_err()
            );
        }
    }

    #[test]
    fn invalid_intervals_do_not_write_lines_or_claim_references() {
        let dir = tempfile::tempdir().unwrap();
        let pool = Pool::new(1);
        let mut outputs = Outputs::new();
        let mut d = Delivery {
            outputs: &mut outputs,
            dir: dir.path(),
            pool: &pool,
            matrix: crate::task::MatrixFormat::Mtx,
            tmp: dir.path(),
            buffer: 4 << 20,
        };
        let path = dir.path().join("bounded.tsv.gz");
        let chroms = vec![("c".to_owned(), 1000)];
        let mut writer = IndexedBed::create(&mut d, &path, "fragments", &chroms).unwrap();
        for (tid, start, end) in [(1, 0, 1), (0, 100, 100), (0, 1000, 1001), (0, 20, 10)] {
            assert!(writer.push(tid, start, end, b"invalid\n").is_err());
        }
        assert!(writer.names.is_empty());
        writer.push(0, 12, 13, b"c\t12\t13\n").unwrap();
        writer.finish(&path, &mut d).unwrap();
        outputs.commit().unwrap();
        assert_eq!(query(&path, "c:1-1000"), ["c\t12\t13"]);
    }

    fn query(path: &Path, region: &str) -> Vec<String> {
        let path = CString::new(path.to_str().unwrap()).unwrap();
        let region = CString::new(region).unwrap();
        // SAFETY: all htslib handles and the malloc'd line buffer are released
        // below; the strings remain alive for the calls using their pointers.
        unsafe {
            let tbx = h::tbx_index_load(path.as_ptr());
            assert!(!tbx.is_null(), "htslib must load the index");
            let fp = h::hts_open(path.as_ptr(), c"r".as_ptr());
            assert!(!fp.is_null());
            let itr = h::tbx_itr_querys1(tbx, region.as_ptr());
            assert!(!itr.is_null());
            let mut line = h::kstring_t::default();
            let mut lines = Vec::new();
            loop {
                let ret = h::tbx_itr_next1(fp, tbx, itr, &mut line);
                if ret < 0 {
                    assert_eq!(ret, -1, "htslib region query failed");
                    break;
                }
                lines.push(
                    String::from_utf8(std::slice::from_raw_parts(line.s.cast(), line.l).to_vec())
                        .unwrap(),
                );
            }
            h::free(line.s.cast());
            h::ow_hts_itr_destroy(itr);
            h::tbx_destroy(tbx);
            assert_eq!(h::hts_close(fp), 0);
            lines
        }
    }

    #[test]
    fn bed_indices_query_sparse_records_at_tbi_boundary_and_long_reference_end() {
        for len in [1000u32, 1 << 29, (1 << 29) + 1, 830_000_000, u32::MAX] {
            let dir = tempfile::tempdir().unwrap();
            let pool = Pool::new(2);
            let mut outputs = Outputs::new();
            let mut d = Delivery {
                outputs: &mut outputs,
                dir: dir.path(),
                pool: &pool,
                matrix: crate::task::MatrixFormat::Mtx,
                tmp: dir.path(),
                buffer: 4 << 20,
            };
            // Absent references before, between and after observed references
            // ensure the name dictionary follows the written records.
            let chroms = vec![
                ("absent-first".to_owned(), len),
                ("long".to_owned(), len),
                ("absent-middle".to_owned(), len),
                ("short".to_owned(), 1000),
                ("absent-last".to_owned(), len),
            ];
            let path = dir.path().join("fragments.tsv.gz");
            let mut w = IndexedBed::create(&mut d, &path, "fragments", &chroms).unwrap();
            let first = "long\t10\t20\tBC\t1";
            let last = format!("long\t{}\t{len}\tBC\t2", len - 1);
            let short = "short\t12\t13\tBC\t3";
            for (tid, start, end, line) in [
                (1, 10, 20, first),
                (1, len - 1, len, last.as_str()),
                (3, 12, 13, short),
            ] {
                w.push(tid, start, end, format!("{line}\n").as_bytes())
                    .unwrap();
                w.g.flush_block().unwrap();
            }
            let index = w.finish(&path, &mut d).unwrap();
            let kind = if len <= 1 << 29 { "tbi" } else { "csi" };
            assert_eq!(index.extension().unwrap(), kind);
            let files = outputs.commit().unwrap();
            assert_eq!(files[1], (kind, index.clone()));
            let other = if kind == "tbi" { "csi" } else { "tbi" };
            assert!(
                !path
                    .with_file_name(format!("fragments.tsv.gz.{other}"))
                    .exists()
            );
            assert_eq!(query(&path, "long:11-20"), [first]);
            assert_eq!(query(&path, &format!("long:{len}-{len}")), [last]);
            assert_eq!(query(&path, "long:21-30"), Vec::<String>::new());
            assert_eq!(query(&path, "short:13-13"), [short]);
        }
    }

    #[test]
    fn empty_bed_indices_load_without_phantom_names() {
        for len in [1000u32, 830_000_000] {
            let dir = tempfile::tempdir().unwrap();
            let pool = Pool::new(1);
            let mut outputs = Outputs::new();
            let mut d = Delivery {
                outputs: &mut outputs,
                dir: dir.path(),
                pool: &pool,
                matrix: crate::task::MatrixFormat::Mtx,
                tmp: dir.path(),
                buffer: 4 << 20,
            };
            let path = dir.path().join("empty.tsv.gz");
            let chroms = vec![("long".to_owned(), len)];
            IndexedBed::create(&mut d, &path, "fragments", &chroms)
                .unwrap()
                .finish(&path, &mut d)
                .unwrap();
            outputs.commit().unwrap();
            let path = CString::new(path.to_str().unwrap()).unwrap();
            // SAFETY: the loaded handle is checked and released once.
            unsafe {
                let tbx = h::tbx_index_load(path.as_ptr());
                assert!(!tbx.is_null());
                h::tbx_destroy(tbx);
            }
        }
    }
}
