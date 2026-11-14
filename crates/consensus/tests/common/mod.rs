//! bcftools consensus as the reference, run in process; an in-memory
//! reference sequence; random fixtures.

#![allow(dead_code)]

use std::collections::HashMap;
use std::ffi::{CString, c_char};
use std::path::Path;
use std::sync::Mutex;

use omicsweave_core::{Error, Region, Result, SequenceRead};

/// splitmix64, for deterministic fixtures from a proptest-chosen seed.
pub struct Rng(pub u64);

impl Rng {
    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    pub fn range(&mut self, lo: u64, hi: u64) -> u64 {
        lo + self.next() % (hi - lo + 1)
    }

    pub fn chance(&mut self, percent: u64) -> bool {
        self.next() % 100 < percent
    }

    pub fn pick<'a, T>(&mut self, items: &'a [T]) -> &'a T {
        &items[self.range(0, items.len() as u64 - 1) as usize]
    }
}

/// Chromosome sequences held in memory.
pub struct Genome(pub HashMap<String, Vec<u8>>);

impl SequenceRead for Genome {
    fn chrom_len(&self, chrom: &str) -> Result<u32> {
        self.0
            .get(chrom)
            .map(|s| s.len() as u32)
            .ok_or_else(|| Error::UnknownChrom(chrom.to_owned()))
    }

    fn sequence(&mut self, region: Region<'_>, out: &mut [u8]) -> Result<()> {
        let seq = self
            .0
            .get(region.chrom)
            .ok_or_else(|| Error::UnknownChrom(region.chrom.to_owned()))?;
        out.copy_from_slice(&seq[region.start as usize..region.end as usize]);
        Ok(())
    }
}

/// Rewrites a VCF as an indexed BCF with htslib (bcftools needs an index).
pub fn to_indexed_bcf(vcf: &Path, bcf: &Path) {
    let c = |s: &str| CString::new(s).unwrap();
    // SAFETY: handles come from htslib and are released below.
    unsafe {
        htslib_sys::hts_set_log_level(htslib_sys::HTS_LOG_OFF);
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

static BCFTOOLS: Mutex<()> = Mutex::new(());

/// How a `bcftools consensus` run ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Outcome {
    /// Exited with this status.
    Done(i32),
    /// Reported an error.
    Failed,
    /// Killed by this signal (bcftools' own memory errors).
    Crashed(i32),
}

/// Runs `bcftools consensus` with `args` in a forked child, so that its
/// process-global state (getopt, error exits) and any memory it corrupts
/// stay out of the test process.
pub fn bcftools_consensus(args: &[String]) -> Outcome {
    let _one = BCFTOOLS.lock().unwrap_or_else(|e| e.into_inner());
    let owned: Vec<CString> = std::iter::once("consensus".to_string())
        .chain(args.iter().cloned())
        .map(|a| CString::new(a).unwrap())
        .collect();
    let mut argv: Vec<*mut c_char> = owned.iter().map(|a| a.as_ptr() as *mut c_char).collect();
    argv.push(std::ptr::null_mut());
    // SAFETY: the child only calls into bcftools on memory prepared before
    // the fork and leaves with _exit; the parent waits for it.
    unsafe {
        let pid = libc::fork();
        assert!(pid >= 0, "fork failed");
        if pid == 0 {
            // bcftools runs with htslib's default log level.
            htslib_sys::hts_set_log_level(htslib_sys::HTS_LOG_WARNING);
            let status = bcftools_sys::ow_bcftools_consensus(owned.len() as i32, argv.as_mut_ptr());
            libc::_exit(if status < 0 { 255 } else { status });
        }
        let mut wstatus = 0;
        assert_eq!(libc::waitpid(pid, &mut wstatus, 0), pid);
        if libc::WIFSIGNALED(wstatus) {
            Outcome::Crashed(libc::WTERMSIG(wstatus))
        } else if libc::WEXITSTATUS(wstatus) == 255 {
            Outcome::Failed
        } else {
            Outcome::Done(libc::WEXITSTATUS(wstatus))
        }
    }
}

/// The sequence lines of a FASTA file, joined.
pub fn fasta_sequence(text: &str) -> Vec<u8> {
    text.lines()
        .filter(|l| !l.starts_with('>'))
        .flat_map(|l| l.trim_end().bytes())
        .collect()
}
