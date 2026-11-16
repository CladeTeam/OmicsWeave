use omicsweave_bam::{BamHeader, Scanner};
use omicsweave_core::{Error, exec::Inline};

/// Encode a binary BAM header independently of SAM parsing and reference lookup.
fn bytes(text: &[u8], refs: &[(&[u8], u32)]) -> Vec<u8> {
    let mut out = b"BAM\x01".to_vec();
    out.extend_from_slice(&(text.len() as u32).to_le_bytes());
    out.extend_from_slice(text);
    out.extend_from_slice(&(refs.len() as i32).to_le_bytes());
    for &(name, len) in refs {
        out.extend_from_slice(&(name.len() as i32).to_le_bytes());
        out.extend_from_slice(name);
        out.extend_from_slice(&len.to_le_bytes());
    }
    out
}

fn read(raw: &[u8]) -> Result<BamHeader, Error> {
    let mut writer = omicsweave_bgzf::Writer::new(Vec::new(), 0, Box::new(Inline::new())).unwrap();
    writer.write(raw).unwrap();
    let compressed = writer.finish().unwrap().out;
    Scanner::new(&compressed[..], Box::new(Inline::new())).header()
}

#[test]
fn binary_reference_order_lengths_and_first_duplicate_lookup_are_preserved() {
    let text = b"@HD\tVN:1.6\n\0\0";
    let header = read(&bytes(
        text,
        &[(b"chr2\0", 200), (b"chr1\0padding", 123), (b"chr2\0", 999)],
    ))
    .unwrap();
    assert_eq!(header.text(), text);
    assert_eq!(
        header
            .references()
            .iter()
            .map(|r| (r.name.as_str(), r.len))
            .collect::<Vec<_>>(),
        [("chr2", 200), ("chr1", 123), ("chr2", 999)]
    );
    assert_eq!(header.tid("chr2"), Some(0));
    assert_eq!(header.tid("chr1"), Some(1));
    assert_eq!(header.chrom_len("chr2").unwrap(), 200);
    assert_eq!(header.tid("absent"), None);
    assert!(
        matches!(header.chrom_len("absent"), Err(Error::UnknownChrom(name)) if name == "absent")
    );
}

#[test]
fn malformed_and_truncated_headers_are_rejected() {
    let valid = bytes(b"@HD\n", &[(b"chr1\0", 100)]);
    for cut in 0..valid.len() {
        assert!(
            read(&valid[..cut]).is_err(),
            "header cut={cut}/{}",
            valid.len()
        );
    }
    let mut bad_magic = valid.clone();
    bad_magic[0] = b'X';
    assert!(matches!(
        read(&bad_magic),
        Err(Error::BadMagic {
            what: "BAM header",
            ..
        })
    ));
    let mut negative_refs = bytes(b"", &[]);
    negative_refs[8..12].copy_from_slice(&(-1i32).to_le_bytes());
    assert!(matches!(read(&negative_refs), Err(Error::Invalid { .. })));
    for name in [&b""[..], &b"\xff\0"[..]] {
        assert!(
            matches!(
                read(&bytes(b"", &[(name, 100)])),
                Err(Error::Invalid { .. })
            ),
            "reference name={name:?}"
        );
    }
    let empty = read(&bytes(b"", &[])).unwrap();
    assert!(empty.text().is_empty());
    assert!(empty.references().is_empty());
}
