//! Independent htslib byte encoding and comparison helpers.
use omicsweave_vcf::{BcfEncoder, BcfScratch, BcfSpan, Header, bcf_header_bytes, read_vcf_header};
use std::ffi::CString;
use std::io::Read;
use std::path::Path;

/// One side's output: the header bytes (`None` when the header is
/// refused) and, per data line, the record bytes and span (`None` when
/// the record is refused).
#[derive(Debug, PartialEq, Eq)]
struct Encoded {
    header: Option<Vec<u8>>,
    records: Vec<Option<(Vec<u8>, BcfSpan)>>,
}

fn data_lines(text: &[u8]) -> Vec<&[u8]> {
    let mut input = text;
    if read_vcf_header(&mut input).is_err() {
        return Vec::new();
    }
    let mut lines: Vec<&[u8]> = input.split(|&b| b == b'\n').collect();
    if lines.last() == Some(&&b""[..]) {
        lines.pop();
    }
    lines
}

fn rid(record: &[u8]) -> u32 {
    u32::from_le_bytes(record[8..12].try_into().unwrap())
}

/// htslib's encoding of the VCF at `path`. With `complete`, a first pass
/// registers undeclared names and both header and records use the header
/// it leaves.
fn htslib(path: &Path, complete: bool) -> Encoded {
    let c = |s: &str| CString::new(s).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.bcf");
    let mut accepted = Vec::new();
    // SAFETY: handles come from htslib, are checked before use and are
    // released before returning.
    unsafe {
        let vcf = c(path.to_str().unwrap());
        let mut fp = htslib_sys::hts_open(vcf.as_ptr(), c("r").as_ptr());
        assert!(!fp.is_null());
        let hdr = htslib_sys::bcf_hdr_read(fp);
        if hdr.is_null() {
            htslib_sys::hts_close(fp);
            return Encoded {
                header: None,
                records: Vec::new(),
            };
        }
        let rec = htslib_sys::bcf_init();
        if complete {
            while htslib_sys::bcf_read(fp, hdr, rec) != -1 {}
            htslib_sys::hts_close(fp);
            fp = htslib_sys::hts_open(vcf.as_ptr(), c("r").as_ptr());
            htslib_sys::bcf_hdr_destroy(htslib_sys::bcf_hdr_read(fp));
        }
        let out = htslib_sys::hts_open(c(out_path.to_str().unwrap()).as_ptr(), c("wbu").as_ptr());
        assert!(!out.is_null());
        assert_eq!(htslib_sys::bcf_hdr_write(out, hdr), 0);
        loop {
            let r = htslib_sys::bcf_read(fp, hdr, rec);
            if r == -1 {
                break;
            }
            let span = (htslib_sys::ow_bcf_pos(rec), htslib_sys::ow_bcf_rlen(rec));
            accepted.push((r == 0 && htslib_sys::bcf_write(out, hdr, rec) == 0).then_some(span));
        }
        htslib_sys::bcf_destroy(rec);
        htslib_sys::hts_close(out);
        htslib_sys::bcf_hdr_destroy(hdr);
        htslib_sys::hts_close(fp);
    }
    // `wbu` is BGZF at level 0 for BCF.
    let mut bytes = Vec::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(&out_path).unwrap())
        .read_to_end(&mut bytes)
        .unwrap();
    let l_text = u32::from_le_bytes(bytes[5..9].try_into().unwrap()) as usize;
    let header = bytes[..9 + l_text].to_vec();
    let mut at = 9 + l_text;
    let records = accepted
        .into_iter()
        .map(|span| {
            let (pos, rlen) = span?;
            let l_shared = u32::from_le_bytes(bytes[at..at + 4].try_into().unwrap()) as usize;
            let l_indiv = u32::from_le_bytes(bytes[at + 4..at + 8].try_into().unwrap()) as usize;
            let record = bytes[at..at + 8 + l_shared + l_indiv].to_vec();
            at += record.len();
            let span = BcfSpan {
                rid: rid(&record),
                beg: pos,
                end: pos + rlen,
            };
            Some((record, span))
        })
        .collect();
    assert_eq!(at, bytes.len(), "htslib wrote records it did not report");
    Encoded {
        header: Some(header),
        records,
    }
}

/// This crate's encoding of `text`, against `header_text` when given (as
/// VCF text) or else the header of `text` itself.
fn ours(text: &[u8], header_text: Option<&[u8]>) -> Encoded {
    let mut input = text;
    let (header, consumed) = match read_vcf_header(&mut input) {
        Ok(h) => h,
        Err(_) => {
            return Encoded {
                header: None,
                records: Vec::new(),
            };
        }
    };
    let header = match header_text {
        Some(t) => Header::from_vcf_text(t).unwrap(),
        None => header,
    };
    let encoder = BcfEncoder::new(&header);
    let mut scratch = BcfScratch::default();
    let mut out = vec![0xAA];
    let records = data_lines(text)
        .into_iter()
        .enumerate()
        .map(|(i, line)| {
            out.truncate(1);
            match encoder.encode_line(&mut scratch, line, consumed + 1 + i as u64, &mut out) {
                Ok(span) => Some((out[1..].to_vec(), span)),
                Err(e) => {
                    assert!(
                        matches!(e, omicsweave_core::Error::Format { line, .. }
                            if line == consumed + 1 + i as u64),
                        "{e:?}"
                    );
                    assert_eq!(out, [0xAA], "a refused record leaves the output as it was");
                    None
                }
            }
        })
        .collect();
    Encoded {
        header: Some(bcf_header_bytes(&header).unwrap()),
        records,
    }
}

fn hex(b: &[u8]) -> String {
    b.iter()
        .map(|x| format!("{x:02x}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn compare(ours: &Encoded, want: &Encoded, text: &str) {
    // Long fixtures are cut in messages; the line at fault is named.
    let text = text.get(..4000).unwrap_or(text);
    match (&ours.header, &want.header) {
        (Some(a), Some(b)) => assert!(
            a == b,
            "header bytes differ\nours:   {:?}\nhtslib: {:?}\ntext:\n{text}",
            String::from_utf8_lossy(a),
            String::from_utf8_lossy(b)
        ),
        (a, b) => assert_eq!(
            a.is_some(),
            b.is_some(),
            "header accepted by one side only\n{text}"
        ),
    }
    let lines = data_lines(text.as_bytes());
    for (i, (a, b)) in ours.records.iter().zip(&want.records).enumerate() {
        let line = String::from_utf8_lossy(lines.get(i).copied().unwrap_or(b"?"));
        match (a, b) {
            (Some((ra, sa)), Some((rb, sb))) => {
                assert!(
                    ra == rb,
                    "record {i} bytes differ: {line:?}\nours:   {}\nhtslib: {}\ntext:\n{text}",
                    hex(ra),
                    hex(rb)
                );
                assert_eq!(sa, sb, "record {i} span: {line:?}");
            }
            (a, b) => assert_eq!(
                a.is_some(),
                b.is_some(),
                "record {i} accepted by one side only: {line:?}\ntext:\n{text}"
            ),
        }
    }
    assert_eq!(
        ours.records.len(),
        want.records.len(),
        "record count\n{text}"
    );
}

/// htslib's BCF header bytes for header text given to `bcf_hdr_parse`;
/// `None` when it refuses the text.
fn htslib_header(text: &[u8]) -> Option<Vec<u8>> {
    let c = |s: &str| CString::new(s).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("h.bcf");
    let mut buf = text.to_vec();
    buf.push(0);
    // SAFETY: `buf` is NUL-terminated and outlives the call; the header
    // and file handles are released before returning.
    unsafe {
        let hdr = htslib_sys::bcf_hdr_init(c("r").as_ptr());
        assert!(!hdr.is_null());
        if htslib_sys::bcf_hdr_parse(hdr, buf.as_mut_ptr().cast()) != 0 {
            htslib_sys::bcf_hdr_destroy(hdr);
            return None;
        }
        let out = htslib_sys::hts_open(c(path.to_str().unwrap()).as_ptr(), c("wbu").as_ptr());
        assert!(!out.is_null());
        assert_eq!(htslib_sys::bcf_hdr_write(out, hdr), 0);
        htslib_sys::hts_close(out);
        htslib_sys::bcf_hdr_destroy(hdr);
    }
    let mut bytes = Vec::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(&path).unwrap())
        .read_to_end(&mut bytes)
        .unwrap();
    Some(bytes)
}

/// Header text parsed both ways must give the same BCF header, or be
/// refused by both.
pub fn check_header(text: &str) {
    let ours = Header::from_vcf_text(text.as_bytes())
        .ok()
        .map(|h| bcf_header_bytes(&h).unwrap());
    let want = htslib_header(text.as_bytes());
    match (&ours, &want) {
        (Some(a), Some(b)) => assert!(
            a == b,
            "header bytes differ\nours:   {:?}\nhtslib: {:?}\ntext:\n{text}",
            String::from_utf8_lossy(a),
            String::from_utf8_lossy(b)
        ),
        _ => assert_eq!(
            ours.is_some(),
            want.is_some(),
            "header accepted by one side only\n{text}"
        ),
    }
}

/// The header part of a VCF: up to the end of the first line that starts
/// with `#` but not `##`.
pub fn header_part(text: &str) -> &str {
    let mut end = 0;
    for line in text.split_inclusive('\n') {
        end += line.len();
        if line.starts_with('#') && !line.starts_with("##") {
            break;
        }
    }
    &text[..end]
}

/// Encodes `text` both ways against its own header.
pub fn check(text: &str) {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.vcf");
    std::fs::write(&path, text).unwrap();
    compare(&ours(text.as_bytes(), None), &htslib(&path, false), text);
}

/// Encodes `text` both ways against the header htslib completes by
/// reading all records once.
pub fn check_completed(text: &str) {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.vcf");
    std::fs::write(&path, text).unwrap();
    let want = htslib(&path, true);
    let Some(header) = &want.header else {
        panic!("htslib refused the header:\n{text}");
    };
    // The BCF header text, read back as VCF text (its IDX are kept).
    let completed = &header[9..header.len() - 1];
    compare(&ours(text.as_bytes(), Some(completed)), &want, text);
}
