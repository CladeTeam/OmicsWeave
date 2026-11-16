//! In-memory reference and indexed BCF fixtures.

use std::collections::HashMap;
use std::ffi::CString;
use std::path::Path;

use omicsweave_bgzf::{BinIndex, locate_index};
use omicsweave_consensus::{Consensus, ConsensusOptions, RecordFilter};
use omicsweave_core::{Error, FileSource, Locus, Result, SequenceRead, orient};
use omicsweave_vcf::BcfReader;

/// Chromosome sequences held in memory.
pub struct Genome(pub HashMap<String, Vec<u8>>);

impl SequenceRead for Genome {
    fn chrom_len(&self, chrom: &str) -> Result<u32> {
        self.0
            .get(chrom)
            .map(|s| s.len() as u32)
            .ok_or_else(|| Error::UnknownChrom(chrom.to_owned()))
    }

    fn sequence(&mut self, locus: &Locus<'_>, out: &mut [u8]) -> Result<()> {
        let seq = self
            .0
            .get(locus.chrom)
            .ok_or_else(|| Error::UnknownChrom(locus.chrom.to_owned()))?;
        let mut at = 0;
        for span in locus.spans() {
            let n = span.len() as usize;
            out[at..at + n].copy_from_slice(&seq[span.start as usize..span.end as usize]);
            at += n;
        }
        orient::orient_bases(locus.strand, out);
        Ok(())
    }
}

/// Rewrites a VCF as an indexed BCF with htslib (bcftools needs an index).
pub fn to_indexed_bcf(vcf: &Path, bcf: &Path) {
    let c = |s: &str| CString::new(s).unwrap();
    // SAFETY: handles come from htslib and are released below.
    unsafe {
        let fp = htslib_sys::hts_open(c(vcf.to_str().unwrap()).as_ptr(), c("r").as_ptr());
        let out = htslib_sys::hts_open(c(bcf.to_str().unwrap()).as_ptr(), c("wb").as_ptr());
        assert!(!fp.is_null() && !out.is_null());
        let hdr = htslib_sys::bcf_hdr_read(fp);
        assert_eq!(htslib_sys::bcf_hdr_write(out, hdr), 0);
        let rec = htslib_sys::bcf_init();
        while htslib_sys::bcf_read(fp, hdr, rec) == 0 {
            assert_eq!(htslib_sys::bcf_write(out, hdr, rec), 0);
        }
        htslib_sys::bcf_destroy(rec);
        htslib_sys::bcf_hdr_destroy(hdr);
        htslib_sys::hts_close(out);
        htslib_sys::hts_close(fp);
        let path = c(bcf.to_str().unwrap());
        assert_eq!(
            htslib_sys::bcf_index_build3(path.as_ptr(), std::ptr::null(), 14, 0),
            0
        );
    }
}

/// Rewrites raw (uncompressed) BCF bytes as BGZF with htslib and indexes
/// them, without re-encoding the records.
pub fn to_indexed_raw_bcf(raw: &[u8], bcf: &Path) {
    let c = |s: &str| CString::new(s).unwrap();
    let path = c(bcf.to_str().unwrap());
    // SAFETY: the handle comes from htslib and is closed below.
    unsafe {
        let fp = htslib_sys::bgzf_open(path.as_ptr(), c("w").as_ptr());
        assert!(!fp.is_null());
        assert_eq!(
            htslib_sys::bgzf_write(fp, raw.as_ptr().cast(), raw.len()),
            raw.len() as isize
        );
        assert_eq!(htslib_sys::bgzf_close(fp), 0);
        assert_eq!(
            htslib_sys::bcf_index_build3(path.as_ptr(), std::ptr::null(), 14, 0),
            0
        );
    }
}

/// An indexed BCF and the CSI beside it.
pub fn open_bcf(bcf: &Path) -> BcfReader<FileSource> {
    let csi = locate_index(bcf).expect("CSI beside the BCF");
    BcfReader::open(
        FileSource::open(bcf).unwrap(),
        BinIndex::read_file(csi).unwrap(),
    )
    .unwrap()
}

/// A builder over the indexed BCF at `bcf`.
pub fn consensus(
    bcf: &Path,
    filter: Option<RecordFilter<'_>>,
    options: ConsensusOptions,
) -> Consensus {
    Consensus::new(open_bcf(bcf), filter, options).unwrap()
}

/// Writes `vcf` as `<dir>/<name>.vcf`, then as an indexed `<dir>/<name>.bcf`,
/// whose path it returns.
pub fn bcf_of(dir: &Path, name: &str, vcf: &[u8]) -> std::path::PathBuf {
    let text = dir.join(format!("{name}.vcf"));
    let bcf = dir.join(format!("{name}.bcf"));
    std::fs::write(&text, vcf).unwrap();
    to_indexed_bcf(&text, &bcf);
    bcf
}
