//! bcftools' filter.c as the reference: an expression compiled with
//! `filter_init` and tested on every record of a file with
//! `filter_test(filter, rec, NULL)` in a fresh helper process, compared with
//! this crate's [`Filter`].

#![allow(dead_code)]

use std::ffi::CString;
use std::path::Path;

use omicsweave_core::Error;
use omicsweave_vcf::{Filter, Record, VariantReader};

/// What became of one record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Outcome {
    Rejected,
    Pass(bool),
    Error,
    /// The C oracle exited, crashed or timed out.
    Undefined,
    /// A Rust implementation failure, including its concrete error.
    ImplementationFailure(String),
}

/// An expression's compile result or its per-record outcomes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Run {
    InitError,
    InitUndefined,
    InitFailure(String),
    Records(Vec<Outcome>),
}

/// Run filter.c in a fresh executable, so its global C state, aborts and
/// resource timeout cannot affect concurrently running Rust tests.
pub fn c_filter(path: &Path, expr: &str) -> Run {
    let directory = tempfile::tempdir().unwrap();
    let protocol = directory.path().join("outcomes");
    let diagnostics = directory.path().join("stderr");
    let mut child = std::process::Command::new(std::env::current_exe().unwrap())
        .args([
            "--ignored",
            "--exact",
            "filter_oracle::oracle_process",
            "--nocapture",
            "--test-threads=1",
        ])
        .env("OMICSWEAVE_FILTER_ORACLE_INPUT", path)
        .env("OMICSWEAVE_FILTER_ORACLE_EXPRESSION", expr)
        .env("OMICSWEAVE_FILTER_ORACLE_OUTPUT", &protocol)
        .stdout(std::process::Stdio::null())
        .stderr(std::fs::File::create(&diagnostics).unwrap())
        .spawn()
        .expect("start isolated bcftools oracle");
    let deadline = std::time::Instant::now() + std::time::Duration::from_secs(25);
    let status = loop {
        if let Some(status) = child.try_wait().expect("poll bcftools oracle") {
            break status;
        }
        if std::time::Instant::now() >= deadline {
            let _ = child.kill();
            break child.wait().expect("reap timed-out bcftools oracle");
        }
        std::thread::sleep(std::time::Duration::from_millis(2));
    };
    // B proves the exact entry point ran; . is written only after its complete
    // result. A renamed helper cannot pass through libtest's zero-test exit.
    let out = std::fs::read(&protocol).unwrap_or_else(|error| {
        panic!(
            "oracle helper did not run: {error}; status={status:?}; stderr={}",
            std::fs::read_to_string(&diagnostics).unwrap_or_default()
        )
    });
    assert_eq!(out.first(), Some(&b'B'), "oracle helper handshake missing");
    let out = &out[1..];
    let status = status.success();
    assert!(
        out.iter()
            .all(|byte| matches!(byte, b'I' | b'K' | b'R' | b'0' | b'1' | b'E' | b'.')),
        "corrupt oracle protocol: {out:?}"
    );
    assert!(
        out.iter()
            .position(|&byte| byte == b'.')
            .is_none_or(|at| at + 1 == out.len()),
        "oracle completion marker before end: {out:?}"
    );
    assert!(
        out.first().is_none_or(|byte| matches!(byte, b'I' | b'K')),
        "invalid oracle start: {out:?}"
    );
    assert!(
        out.get(1..)
            .is_none_or(|body| !body.contains(&b'K') && !body.contains(&b'I')),
        "duplicate oracle initialization: {out:?}"
    );
    assert!(
        out.last() != Some(&b'.') || status,
        "oracle failed after reporting a complete result"
    );
    if out.first() == Some(&b'I') {
        assert_eq!(out, b"I.", "incomplete oracle compile-error response");
    }
    if let Some(at) = out.iter().position(|&byte| byte == b'E') {
        assert_eq!(
            &out[at..],
            b"E.",
            "incomplete oracle runtime-error response"
        );
    }
    match out.first() {
        Some(b'I') if status && out.last() == Some(&b'.') => return Run::InitError,
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
    let finished = out.last() == Some(&b'.');
    if !status || !finished {
        records.push(Outcome::Undefined);
    }
    Run::Records(records)
}

/// The only C filter invocation in the child. `--exact` and a fresh process
/// keep bcftools' non-reentrant error handler local to this one expression.
#[test]
#[ignore = "subprocess entry point; invoked by c_filter"]
fn oracle_process() {
    use std::os::fd::AsRawFd;
    let Some(protocol) = std::env::var_os("OMICSWEAVE_FILTER_ORACLE_OUTPUT") else {
        return;
    };
    let output = std::fs::File::create(protocol).unwrap();
    let fd = output.as_raw_fd();
    let path = CString::new(std::env::var("OMICSWEAVE_FILTER_ORACLE_INPUT").unwrap()).unwrap();
    let expression =
        CString::new(std::env::var("OMICSWEAVE_FILTER_ORACLE_EXPRESSION").unwrap()).unwrap();
    // SAFETY: each pointer is checked before use; all C strings outlive the
    // calls. Each record has a distinct address, as filter.c caches by address.
    // The process owns these C allocations until _exit; no other test runs in it.
    unsafe {
        let put = |byte: u8| {
            assert_eq!(libc::write(fd, (&byte as *const u8).cast(), 1), 1);
        };
        std::panic::set_hook(Box::new(move |_| {
            let marker = b'!';
            libc::write(fd, (&marker as *const u8).cast(), 1);
        }));
        put(b'B');
        libc::alarm(20);
        let input = htslib_sys::hts_open(path.as_ptr(), c"r".as_ptr());
        assert!(!input.is_null(), "open oracle input");
        let header = htslib_sys::bcf_hdr_read(input);
        assert!(!header.is_null(), "read oracle header");
        let filter = bcftools_sys::ow_filter_init(header, expression.as_ptr());
        if filter.is_null() {
            put(b'I');
            put(b'.');
            libc::_exit(0);
        }
        put(b'K');
        loop {
            let record = htslib_sys::bcf_init();
            assert!(!record.is_null());
            let status = htslib_sys::bcf_read(input, header, record);
            if status == -1 {
                break;
            }
            if status < -1 {
                put(b'R');
                continue;
            }
            match bcftools_sys::ow_filter_test(filter, record) {
                1 => put(b'1'),
                0 => put(b'0'),
                _ => {
                    put(b'E');
                    put(b'.');
                    libc::_exit(0);
                }
            }
        }
        put(b'.');
        libc::_exit(0);
    }
}

/// This crate on the same bytes. No implementation error classifies the C oracle.
pub fn ours(bytes: &[u8], expr: &str) -> Run {
    let mut reader = VariantReader::new(bytes).unwrap();
    let mut filter = match Filter::new(reader.header(), expr) {
        Ok(f) => f,
        Err(Error::Invalid { .. }) => return Run::InitError,
        Err(error) => return Run::InitFailure(format!("{error:?}")),
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
                Err(error) => {
                    out.push(Outcome::ImplementationFailure(format!("{error:?}")));
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
    pub oracle_init_undefined: usize,
    pub records: usize,
    pub passed: usize,
    pub runtime_errors: usize,
    pub oracle_runtime_undefined: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Comparison {
    Compared,
    OracleUndefined,
}

/// Only an undefined reference can end a comparison without deciding it.
/// Only an undefined C oracle can end a comparison without deciding it.
/// Every implementation failure against a defined oracle is an error.
pub fn compare(c: &Run, o: &Run, tally: &mut Tally) -> Result<Comparison, String> {
    match (c, o) {
        (Run::InitUndefined, _) => {
            tally.oracle_init_undefined += 1;
            Ok(Comparison::OracleUndefined)
        }
        (Run::InitError, Run::InitError) => {
            tally.init_errors += 1;
            Ok(Comparison::Compared)
        }
        (Run::Records(c), Run::Records(o)) => {
            for (i, (a, b)) in c.iter().zip(o).enumerate() {
                match (a, b) {
                    (Outcome::Undefined, _) => {
                        tally.oracle_runtime_undefined += 1;
                        return Ok(Comparison::OracleUndefined);
                    }
                    (Outcome::Error, Outcome::Error) => {
                        tally.runtime_errors += 1;
                        return Ok(Comparison::Compared);
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
            if c.get(o.len()) == Some(&Outcome::Undefined) {
                tally.oracle_runtime_undefined += 1;
                return Ok(Comparison::OracleUndefined);
            }
            if c.len() != o.len() {
                return Err(format!(
                    "record count: filter.c {}, ours {}",
                    c.len(),
                    o.len()
                ));
            }
            Ok(Comparison::Compared)
        }
        (c, o) => Err(format!("compile: filter.c {c:?}, ours {o:?}")),
    }
}
