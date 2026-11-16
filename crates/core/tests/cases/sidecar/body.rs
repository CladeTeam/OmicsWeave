use omicsweave_core::Error;
use omicsweave_core::sidecar::{BodyReader, BodyWriter};

#[test]
fn oversized_body_read_reports_error_without_advancing() {
    let mut reader = BodyReader::new(&[3, 4, 5], "test");
    assert_eq!(reader.u8().unwrap(), 3);
    assert!(matches!(
        reader.bytes(usize::MAX),
        Err(Error::Truncated { have: 3, .. })
    ));
    assert_eq!(reader.remaining(), 2);
    assert_eq!(reader.bytes(2).unwrap(), &[4, 5]);
    reader.finish().unwrap();
}

#[test]
fn word_encoding_is_little_endian_and_reads_require_complete_fields() {
    let mut writer = BodyWriter(Vec::new());
    writer.u8(0x12);
    writer.u16(0x3456);
    writer.u32(0x789a_bcde);
    writer.u64(0x0123_4567_89ab_cdef);
    assert_eq!(
        writer.0,
        [
            0x12, 0x56, 0x34, 0xde, 0xbc, 0x9a, 0x78, 0xef, 0xcd, 0xab, 0x89, 0x67, 0x45, 0x23,
            0x01
        ]
    );
    let mut reader = BodyReader::new(&writer.0, "words");
    assert!(reader.finish().is_err());
    assert_eq!(reader.u8().unwrap(), 0x12);
    assert_eq!(reader.u16().unwrap(), 0x3456);
    assert_eq!(reader.u32().unwrap(), 0x789a_bcde);
    assert_eq!(reader.u64().unwrap(), 0x0123_4567_89ab_cdef);
    assert_eq!(reader.remaining(), 0);
    reader.finish().unwrap();
    assert!(matches!(reader.u8(), Err(Error::Truncated { .. })));
    assert_eq!(reader.remaining(), 0);

    for width in [2, 4, 8] {
        for len in 0..width {
            let bytes = vec![0xff; len];
            let mut reader = BodyReader::new(&bytes, "partial word");
            let result = match width {
                2 => reader.u16().map(u64::from),
                4 => reader.u32().map(u64::from),
                8 => reader.u64(),
                _ => unreachable!(),
            };
            assert!(
                matches!(result, Err(Error::Truncated { .. })),
                "width={width} len={len}"
            );
            assert_eq!(reader.remaining(), len, "width={width} len={len}");
            assert_eq!(reader.bytes(len).unwrap(), bytes);
            reader.finish().unwrap();
        }
    }
}

#[test]
fn record_counts_are_bounded_by_remaining_bytes_before_allocation() {
    for (count, width, available, valid) in [
        (0u64, 4usize, 0usize, true),
        (2, 4, 8, true),
        (2, 4, 7, false),
        (u64::MAX, 8, 0, false),
        (2, usize::MAX, 0, false),
    ] {
        let mut bytes = count.to_le_bytes().to_vec();
        bytes.resize(8 + available, 0);
        let mut reader = BodyReader::new(&bytes, "records");
        let result = reader.count(width);
        assert_eq!(
            result.is_ok(),
            valid,
            "count={count} width={width} available={available}"
        );
        if valid {
            assert_eq!(result.unwrap(), count as usize);
        }
        assert_eq!(reader.remaining(), available);
    }
}
