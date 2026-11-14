use omicsweave_core::{Error, FileSource, ReadAt};

const BYTES: &[u8] = include_bytes!("fixtures/read-at.txt");

fn compare_reads(source: &(impl ReadAt + ?Sized)) {
    assert_eq!(source.len(), BYTES.len() as u64);
    assert!(!source.is_empty());
    for offset in 0..=BYTES.len() + 1 {
        for len in 0..=BYTES.len() + 1 {
            let mut guarded = vec![0x55; len + 2];
            let result = source.read_exact_at(&mut guarded[1..1 + len], offset as u64);
            if offset <= BYTES.len() && len <= BYTES.len() - offset {
                result.unwrap();
                assert_eq!(&guarded[1..1 + len], &BYTES[offset..offset + len]);
            } else {
                assert!(matches!(result, Err(Error::OutOfRange { .. })));
                assert!(guarded.iter().all(|&byte| byte == 0x55));
            }
            assert_eq!(guarded[0], 0x55);
            assert_eq!(guarded[len + 1], 0x55);
        }
    }
    for len in [0, 1, 2] {
        let mut buf = vec![0x55; len];
        assert!(matches!(
            source.read_exact_at(&mut buf, u64::MAX),
            Err(Error::OutOfRange { .. })
        ));
        assert!(buf.iter().all(|&byte| byte == 0x55));
    }
}

#[test]
fn memory_and_file_sources_share_read_boundaries() {
    compare_reads(BYTES);
    compare_reads(&BYTES);
    compare_reads(&BYTES.to_vec());
    let file = FileSource::open(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/read-at.txt"
    ))
    .unwrap();
    compare_reads(&file);
}

#[test]
fn empty_source_accepts_only_empty_read_at_zero() {
    let source: &[u8] = &[];
    assert!(source.is_empty());
    source.read_exact_at(&mut [], 0).unwrap();
    assert!(source.read_exact_at(&mut [], 1).is_err());
    assert!(source.read_exact_at(&mut [0], 0).is_err());
}
