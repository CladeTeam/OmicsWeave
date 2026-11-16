use super::reference::Encoded;
use omicsweave_sam::{BamEncoder, BamScratch, SamHeader, bam_header_bytes};

/// The data lines after the header, as a line splitter hands them over:
/// without their `\n`, and no empty line after a final `\n`.
pub fn data_lines(rest: &[u8]) -> Vec<&[u8]> {
    if rest.is_empty() {
        return Vec::new();
    }
    let body = rest.strip_suffix(b"\n").unwrap_or(rest);
    body.split(|&c| c == b'\n').collect()
}

/// This crate's encoding, appending every record to one buffer as a
/// worker does.
pub fn ours(sam: &[u8]) -> Encoded {
    let mut input = sam;
    let Ok(header) = SamHeader::read(&mut input) else {
        return Encoded {
            header: None,
            records: Vec::new(),
        };
    };
    assert_eq!(header.raw(), &sam[..sam.len() - input.len()]);
    let encoder = BamEncoder::new(&header);
    let mut scratch = BamScratch::new();
    let mut all = b"prefix".to_vec();
    let mut records = Vec::new();
    for (i, line) in data_lines(input).into_iter().enumerate() {
        let before = all.len();
        let line_no = header.lines() + 1 + i as u64;
        match encoder.encode_line(&mut scratch, line, line_no, &mut all) {
            Ok(span) => {
                let rec = all[before..].to_vec();
                assert_eq!(
                    u32::from_le_bytes(rec[..4].try_into().unwrap()) as usize,
                    rec.len() - 4
                );
                records.push(Some((rec, span)));
            }
            Err(e) => {
                assert_eq!(all.len(), before, "a failed line left bytes behind");
                assert!(
                    matches!(e, omicsweave_core::Error::Format { line, .. } if line == line_no)
                );
                records.push(None);
            }
        }
    }
    assert!(all.starts_with(b"prefix"));
    Encoded {
        header: Some(bam_header_bytes(&header)),
        records,
    }
}
