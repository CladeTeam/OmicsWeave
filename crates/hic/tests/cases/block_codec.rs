use omicsweave_hic::{Record, decode_block, encode_block};
use proptest::prelude::*;

#[test]
fn encoder_matches_literal_dense_and_sparse_wire_layouts() {
    for version in [8, 9] {
        for float in [false, true] {
            for dense in [false, true] {
                let records = if dense {
                    vec![(10, 20), (11, 20), (10, 21), (11, 21)]
                } else {
                    vec![(10, 20), (100, 40)]
                };
                let records: Vec<_> = records
                    .into_iter()
                    .enumerate()
                    .map(|(i, (x, y))| Record {
                        x,
                        y,
                        value: (i + 1) as f32 + if float { 0.5 } else { 0.0 },
                    })
                    .collect();
                let mut expected = (records.len() as i32).to_le_bytes().to_vec();
                expected.extend_from_slice(&10i32.to_le_bytes());
                expected.extend_from_slice(&20i32.to_le_bytes());
                expected.push(u8::from(float));
                if version == 9 {
                    expected.extend_from_slice(&[0, 0]);
                }
                if dense {
                    // Four values occupy a 2-by-2 rectangle, with no holes.
                    expected.extend_from_slice(&[2, 4, 0, 0, 0, 2, 0]);
                    for record in &records {
                        if float {
                            expected.extend_from_slice(&record.value.to_le_bytes());
                        } else {
                            expected.extend_from_slice(&(record.value as i16).to_le_bytes());
                        }
                    }
                } else {
                    // Two rows, each with one x offset and one value.
                    expected.extend_from_slice(&[1, 2, 0]);
                    for (record, offsets) in records.iter().zip([[0, 0], [90, 20]]) {
                        for field in [offsets[1], 1, offsets[0]] {
                            expected.extend_from_slice(&(field as i16).to_le_bytes());
                        }
                        if float {
                            expected.extend_from_slice(&record.value.to_le_bytes());
                        } else {
                            expected.extend_from_slice(&(record.value as i16).to_le_bytes());
                        }
                    }
                }
                let mut encoded = Vec::new();
                encode_block(version, &records, &mut encoded).unwrap();
                assert_eq!(
                    encoded, expected,
                    "v{version}, float={float}, dense={dense}"
                );
                let mut decoded = Vec::new();
                decode_block(version, &expected, |x, y, value| {
                    decoded.push(Record {
                        x: x as u32,
                        y: y as u32,
                        value,
                    });
                })
                .unwrap();
                assert_eq!(decoded, records);
            }
        }
    }
}

#[test]
fn sparse_extreme_coordinates_never_size_an_ineligible_dense_rectangle() {
    let records = [
        Record {
            x: 0,
            y: 0,
            value: 0.5,
        },
        Record {
            x: i32::MAX as u32,
            y: i32::MAX as u32,
            value: 0.5,
        },
    ];
    let mut encoded = Vec::new();
    encode_block(9, &records, &mut encoded).unwrap();
    let mut expected = Vec::new();
    for field in [2i32, 0, 0] {
        expected.extend_from_slice(&field.to_le_bytes());
    }
    expected.extend_from_slice(&[1, 1, 1, 1]);
    expected.extend_from_slice(&2i32.to_le_bytes());
    for coordinate in [0i32, i32::MAX] {
        for field in [coordinate, 1, coordinate] {
            expected.extend_from_slice(&field.to_le_bytes());
        }
        expected.extend_from_slice(&0.5f32.to_le_bytes());
    }
    assert_eq!(encoded, expected);
    assert_eq!(encoded.len(), 52);
    let mut decoded = Vec::new();
    decode_block(9, &encoded, |x, y, value| {
        decoded.push(Record {
            x: x as u32,
            y: y as u32,
            value,
        });
    })
    .unwrap();
    assert_eq!(decoded, records);
}

fn put_int(bytes: &mut Vec<u8>, wide: bool, value: i32) {
    if wide {
        bytes.extend_from_slice(&value.to_le_bytes());
    } else {
        bytes.extend_from_slice(&(value as i16).to_le_bytes());
    }
}

#[test]
fn sparse_coordinates_are_checked_for_every_field_layout() {
    for wide_x in [false, true] {
        for wide_y in [false, true] {
            for float in [false, true] {
                for bad_x in [false, true] {
                    let mut bytes = 1i32.to_le_bytes().to_vec();
                    let (x0, y0) = if bad_x { (i32::MAX, 0) } else { (0, i32::MIN) };
                    bytes.extend_from_slice(&x0.to_le_bytes());
                    bytes.extend_from_slice(&y0.to_le_bytes());
                    bytes.extend_from_slice(&[
                        u8::from(float),
                        u8::from(wide_x),
                        u8::from(wide_y),
                        1,
                    ]);
                    put_int(&mut bytes, wide_y, 1);
                    put_int(&mut bytes, wide_y, if bad_x { 0 } else { -1 });
                    put_int(&mut bytes, wide_x, 1);
                    put_int(&mut bytes, wide_x, i32::from(bad_x));
                    if float {
                        bytes.extend_from_slice(&3.0f32.to_le_bytes());
                    } else {
                        bytes.extend_from_slice(&3i16.to_le_bytes());
                    }
                    let mut callbacks = 0;
                    assert!(decode_block(9, &bytes, |_, _, _| callbacks += 1).is_err());
                    assert_eq!(callbacks, 0);
                }
            }
        }
    }
}

#[test]
fn codecs_cover_narrow_wide_and_unaligned_records() {
    for version in [8, 9] {
        for wide in [false, true] {
            for short in [false, true] {
                if wide && version < 9 {
                    continue;
                }
                let records: Vec<_> = (0..64u32)
                    .flat_map(|y| {
                        (0..8u32).map(move |x| Record {
                            x: if wide { x * 50_000 } else { x },
                            y: if wide { y * 40_000 } else { y },
                            value: if short { 32767.0 } else { (x + y) as f32 / 7.0 },
                        })
                    })
                    .collect();
                let mut bytes = Vec::new();
                encode_block(version, &records, &mut bytes).unwrap();
                let mut aligned = vec![0xAA];
                aligned.extend(&bytes);
                let mut actual = Vec::new();
                decode_block(version, &aligned[1..], |x, y, value| {
                    actual.push(Record {
                        x: x as u32,
                        y: y as u32,
                        value,
                    })
                })
                .unwrap();
                assert_eq!(actual, records);
                for end in 0..bytes.len() {
                    assert!(decode_block(version, &bytes[..end], |_, _, _| {}).is_err());
                }
            }
        }
    }
    assert!(decode_block(6, &(-1i32).to_le_bytes(), |_, _, _| {}).is_err());
    for records in [
        vec![
            Record {
                x: 0,
                y: 1,
                value: 1.0,
            },
            Record {
                x: 0,
                y: 0,
                value: 2.0,
            },
        ],
        vec![Record {
            x: 0,
            y: 0,
            value: f32::INFINITY,
        }],
    ] {
        assert!(encode_block(9, &records, &mut Vec::new()).is_err());
    }
}

#[test]
fn inflater_grows_reuses_and_rejects_corruption() {
    use omicsweave_hic::{Deflater, Inflater};
    let text = vec![b'A'; 1 << 18];
    let mut encoder = Deflater::new(9);
    let mut raw = Vec::new();
    encoder.deflate(&text, &mut raw).unwrap();
    let mut decoder = Inflater::new();
    let mut out = Vec::new();
    decoder.inflate(&raw, &mut out).unwrap();
    assert_eq!(out, text);
    for end in [0, 1, raw.len() / 2, raw.len() - 1] {
        assert!(decoder.inflate(&raw[..end], &mut out).is_err());
    }
    let mut broken = raw.clone();
    *broken.last_mut().unwrap() ^= 1;
    assert!(decoder.inflate(&broken, &mut out).is_err());
    decoder.inflate(&raw, &mut out).unwrap();
    assert_eq!(out, text);
    encoder.deflate(b"xy", &mut raw).unwrap();
    decoder.inflate(&raw, &mut out).unwrap();
    assert_eq!(out, b"xy");
}

#[test]
fn blocks_of_other_versions_are_unsupported() {
    let record = [Record {
        x: 0,
        y: 0,
        value: 1.0,
    }];
    for version in [6, 7, 10] {
        let mut out = Vec::new();
        assert!(matches!(
            encode_block(version, &record, &mut out),
            Err(omicsweave_core::Error::Unsupported(_))
        ));
        assert!(matches!(
            decode_block(version, &[0; 16], |_, _, _| {}),
            Err(omicsweave_core::Error::Unsupported(_))
        ));
    }
}

proptest! {
    #![proptest_config(ProptestConfig { failure_persistence: Some(Box::new(proptest::test_runner::FileFailurePersistence::Direct(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/cases/block_codec.proptest-regressions")))), ..ProptestConfig::default() })]
    #[test]
    fn arbitrary_block_bytes_do_not_panic(bytes in prop::collection::vec(any::<u8>(),0..512), version in 6i32..=9) {
        let _ = decode_block(version,&bytes,|_,_,_|{});
    }
}
