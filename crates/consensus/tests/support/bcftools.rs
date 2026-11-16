//! A fresh executable isolates C global state and memory failures from Rust tests.

use std::ffi::{CString, c_char};
use std::io::Write as _;
use std::os::unix::process::ExitStatusExt as _;
use std::process::{Command, Stdio};

/// How a `bcftools consensus` invocation ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Outcome {
    Done(i32),
    Failed,
    Crashed(i32),
}

pub fn bcftools_consensus(args: &[String]) -> Outcome {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("arguments");
    let protocol = dir.path().join("outcome");
    let diagnostics = dir.path().join("stderr");
    let mut bytes = Vec::new();
    for arg in args {
        assert!(!arg.as_bytes().contains(&0), "NUL in oracle argument");
        bytes.extend_from_slice(arg.as_bytes());
        bytes.push(0);
    }
    std::fs::write(&input, bytes).unwrap();
    // libtest omits the crate name from fully qualified test names.
    let module = module_path!().split_once("::").unwrap().1;
    let mut child = Command::new(std::env::current_exe().unwrap())
        .args([
            "--ignored",
            "--exact",
            &format!("{module}::oracle_process"),
            "--nocapture",
            "--test-threads=1",
        ])
        .env("OMICSWEAVE_CONSENSUS_ORACLE_INPUT", &input)
        .env("OMICSWEAVE_CONSENSUS_ORACLE_OUTPUT", &protocol)
        .stdout(Stdio::null())
        .stderr(std::fs::File::create(&diagnostics).unwrap())
        .spawn()
        .expect("start isolated bcftools consensus");
    let deadline = std::time::Instant::now() + std::time::Duration::from_secs(30);
    let status = loop {
        if let Some(status) = child.try_wait().expect("poll oracle") {
            break status;
        }
        if std::time::Instant::now() >= deadline {
            let _ = child.kill();
            let status = child.wait().expect("reap timed-out oracle");
            panic!("bcftools consensus timed out: {status:?}; args={args:?}");
        }
        std::thread::sleep(std::time::Duration::from_millis(2));
    };
    let out = std::fs::read(&protocol).unwrap_or_else(|error| {
        panic!(
            "oracle helper did not run: {error}; status={status:?}; args={args:?}; stderr={}",
            std::fs::read_to_string(&diagnostics).unwrap_or_default()
        )
    });
    assert_eq!(out.first(), Some(&b'B'), "oracle handshake: {out:?}");
    assert!(
        !out.contains(&b'!'),
        "oracle helper panicked: args={args:?}"
    );
    if let Some(signal) = status.signal() {
        assert_eq!(out, b"B", "oracle crashed after reporting a result");
        assert_ne!(signal, libc::SIGALRM, "bcftools consensus timed out");
        return Outcome::Crashed(signal);
    }
    assert!(status.success(), "oracle helper failed: {status:?}");
    let result = std::str::from_utf8(&out[1..])
        .unwrap()
        .strip_suffix('.')
        .expect("oracle completion marker; an empty --exact match must fail");
    let code: i32 = result.parse().expect("oracle status");
    if code < 0 {
        Outcome::Failed
    } else {
        Outcome::Done(code)
    }
}

#[test]
#[ignore = "subprocess entry point; invoked by bcftools_consensus"]
fn oracle_process() {
    use std::os::fd::AsRawFd as _;
    let Some(protocol) = std::env::var_os("OMICSWEAVE_CONSENSUS_ORACLE_OUTPUT") else {
        return;
    };
    let mut output = std::fs::File::create(protocol).unwrap();
    let fd = output.as_raw_fd();
    let bytes =
        std::fs::read(std::env::var_os("OMICSWEAVE_CONSENSUS_ORACLE_INPUT").unwrap()).unwrap();
    let owned: Vec<CString> = std::iter::once(CString::new("consensus").unwrap())
        .chain(
            bytes
                .strip_suffix(&[0])
                .unwrap_or(&bytes)
                .split(|&byte| byte == 0)
                .map(|arg| CString::new(arg).unwrap()),
        )
        .collect();
    let mut argv: Vec<*mut c_char> = owned.iter().map(|arg| arg.as_ptr().cast_mut()).collect();
    argv.push(std::ptr::null_mut());
    std::panic::set_hook(Box::new(move |_| {
        // SAFETY: the child owns the file descriptor until it exits.
        unsafe {
            libc::write(fd, b"!".as_ptr().cast(), 1);
        }
    }));
    output.write_all(b"B").unwrap();
    // SAFETY: all C strings outlive the invocation, argv is NUL terminated,
    // and this process runs only this test. The alarm bounds a stuck C oracle.
    let code = unsafe {
        libc::alarm(25);
        bcftools_sys::ow_bcftools_consensus(owned.len() as i32, argv.as_mut_ptr())
    };
    write!(output, "{code}.").unwrap();
    output.flush().unwrap();
    // SAFETY: C process-global state is intentionally discarded at exit.
    unsafe {
        libc::_exit(0);
    }
}
