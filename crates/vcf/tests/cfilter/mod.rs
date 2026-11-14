//! bcftools' filter.c as the reference: an expression compiled with
//! `filter_init` and tested on every record of a file with
//! `filter_test(filter, rec, NULL)`, in a forked child so that its exits and
//! crashes end only the child. And the same with this crate's [`Filter`].

#![allow(dead_code)]

use std::ffi::CString;
use std::path::Path;

use omicsweave_core::Error;
use omicsweave_vcf::{Filter, Record, VariantReader};

/// What became of one record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Outcome {
    /// The reader rejected the record.
    Rejected,
    Pass(bool),
    /// An error filter.c reports (`error()`), or [`Error::Invalid`].
    Error,
    /// The C process ended otherwise (exit, assertion, signal), or this
    /// crate reports behaviour bcftools leaves undefined.
    Undefined,
}

/// An expression's results over a file: `None` when it does not compile
/// (`Some(Undefined)`-like compile outcomes are reported as `Err`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Run {
    InitError,
    InitUndefined,
    Records(Vec<Outcome>),
}

/// filter.c on `path`.
pub fn c_filter(path: &Path, expr: &str) -> Run {
    let path = CString::new(path.to_str().unwrap()).unwrap();
    let Ok(expr) = CString::new(expr) else {
        return Run::InitUndefined;
    };
    let mut fds = [0i32; 2];
    let mut out = Vec::new();
    let status;
    // SAFETY: the child only calls into htslib/bcftools and async-safe libc
    // functions before _exit; the parent reads the pipe and reaps the child.
    // The test binary runs this from a single test thread.
    unsafe {
        assert_eq!(libc::pipe(fds.as_mut_ptr()), 0);
        let pid = libc::fork();
        assert!(pid >= 0, "fork failed");
        if pid == 0 {
            libc::close(fds[0]);
            let null = libc::open(c"/dev/null".as_ptr(), libc::O_WRONLY);
            libc::dup2(null, 2);
            // A computation that does not end is reported as undefined.
            libc::alarm(20);
            let put = |c: u8| {
                libc::write(fds[1], (&c as *const u8).cast(), 1);
            };
            htslib_sys::hts_set_log_level(htslib_sys::HTS_LOG_OFF);
            let fp = htslib_sys::hts_open(path.as_ptr(), c"r".as_ptr());
            let hdr = htslib_sys::bcf_hdr_read(fp);
            let filter = bcftools_sys::ow_filter_init(hdr, expr.as_ptr());
            if filter.is_null() {
                put(b'I');
                libc::_exit(0);
            }
            put(b'K');
            loop {
                // Every record gets its own object: filter.c caches per
                // record address.
                let rec = htslib_sys::bcf_init();
                let r = htslib_sys::bcf_read(fp, hdr, rec);
                if r == -1 {
                    break;
                }
                if r < -1 {
                    put(b'R');
                    continue;
                }
                match bcftools_sys::ow_filter_test(filter, rec) {
                    1 => put(b'1'),
                    0 => put(b'0'),
                    _ => {
                        put(b'E');
                        libc::_exit(0);
                    }
                }
            }
            put(b'.');
            libc::_exit(0);
        }
        libc::close(fds[1]);
        let mut buf = [0u8; 4096];
        loop {
            let n = libc::read(fds[0], buf.as_mut_ptr().cast(), buf.len());
            if n <= 0 {
                break;
            }
            out.extend_from_slice(&buf[..n as usize]);
        }
        libc::close(fds[0]);
        let mut wstatus = 0;
        assert_eq!(libc::waitpid(pid, &mut wstatus, 0), pid);
        status = libc::WIFEXITED(wstatus) && libc::WEXITSTATUS(wstatus) == 0;
    }
    match out.first() {
        Some(b'I') => return Run::InitError,
        Some(b'K') => {}
        _ => return Run::InitUndefined,
    }
    let mut records: Vec<Outcome> = out[1..]
        .iter()
        .filter(|&&c| c != b'.')
        .map(|&c| match c {
            b'R' => Outcome::Rejected,
            b'1' => Outcome::Pass(true),
            b'0' => Outcome::Pass(false),
            _ => Outcome::Error,
        })
        .collect();
    let finished = out.last() == Some(&b'.') || out.last() == Some(&b'E');
    if !status || !finished {
        records.push(Outcome::Undefined);
    }
    Run::Records(records)
}

/// This crate on the same bytes.
pub fn ours(bytes: &[u8], expr: &str) -> Run {
    let mut reader = VariantReader::new(bytes).unwrap();
    let mut filter = match Filter::new(reader.header(), expr) {
        Ok(f) => f,
        Err(Error::Invalid { .. }) => return Run::InitError,
        Err(_) => return Run::InitUndefined,
    };
    let mut rec = Record::default();
    let mut out = Vec::new();
    loop {
        match reader.read(&mut rec) {
            Ok(false) => break,
            Err(_) => out.push(Outcome::Rejected),
            Ok(true) => match filter.test(reader.header(), &rec) {
                Ok(b) => out.push(Outcome::Pass(b)),
                Err(Error::Invalid { .. }) => {
                    out.push(Outcome::Error);
                    break;
                }
                Err(_) => {
                    out.push(Outcome::Undefined);
                    break;
                }
            },
        }
    }
    Run::Records(out)
}

/// Counts of how comparisons ended.
#[derive(Debug, Default, Clone)]
pub struct Tally {
    pub init_errors: usize,
    pub init_undefined: usize,
    pub records: usize,
    pub passed: usize,
    pub runtime_errors: usize,
    pub undefined: usize,
}

/// Compares the two runs record by record up to the first error or
/// undefined outcome on either side; `Err` describes a difference.
pub fn compare(c: &Run, o: &Run, tally: &mut Tally) -> Result<(), String> {
    match (c, o) {
        (Run::InitUndefined, _) | (_, Run::InitUndefined) => {
            tally.init_undefined += 1;
            Ok(())
        }
        (Run::InitError, Run::InitError) => {
            tally.init_errors += 1;
            Ok(())
        }
        (Run::Records(c), Run::Records(o)) => {
            for (i, (a, b)) in c.iter().zip(o).enumerate() {
                match (a, b) {
                    (Outcome::Undefined, _) | (_, Outcome::Undefined) => {
                        tally.undefined += 1;
                        return Ok(());
                    }
                    (Outcome::Error, Outcome::Error) => {
                        tally.runtime_errors += 1;
                        return Ok(());
                    }
                    (a, b) if a == b => {
                        tally.records += 1;
                        if *a == Outcome::Pass(true) {
                            tally.passed += 1;
                        }
                    }
                    (a, b) => return Err(format!("record {i}: filter.c {a:?}, ours {b:?}")),
                }
            }
            if c.len() != o.len() {
                return Err(format!(
                    "record count: filter.c {}, ours {}",
                    c.len(),
                    o.len()
                ));
            }
            Ok(())
        }
        (c, o) => Err(format!("compile: filter.c {c:?}, ours {o:?}")),
    }
}
