use omicsweave_vcf::{BT_INT16, BT_INT32, BT_NULL, FormatField, Typed};

#[test]
fn public_typed_views_reject_overflow_without_large_arrays() {
    let bytes = [42, 0, 0, 0, 9, 0, 0, 0];
    let typed = Typed {
        ty: BT_INT32,
        bytes: &bytes,
    };
    assert_eq!(typed.raw(0), 42);
    assert_eq!(typed.raw(1), 9);
    for i in [2, 1usize << (usize::BITS - 2), usize::MAX] {
        assert!(std::panic::catch_unwind(|| typed.raw(i)).is_err());
        assert!(std::panic::catch_unwind(|| typed.value(i)).is_err());
    }
    let f = FormatField {
        key: 0,
        n: 1,
        values: Typed {
            ty: BT_INT16,
            bytes: &bytes,
        },
    };
    assert_eq!(f.size(), 2);
    assert_eq!(f.sample(1).bytes, &bytes[2..4]);
    for i in [4, 1usize << (usize::BITS - 1), usize::MAX] {
        assert!(std::panic::catch_unwind(|| f.sample(i)).is_err());
    }
    let huge = FormatField {
        n: 1usize << (usize::BITS - 1),
        ..f
    };
    assert!(std::panic::catch_unwind(|| huge.size()).is_err());
    let one_byte = FormatField {
        values: Typed {
            ty: omicsweave_vcf::BT_INT8,
            bytes: &bytes,
        },
        ..f
    };
    assert!(std::panic::catch_unwind(|| one_byte.sample(usize::MAX)).is_err());
    for zero in [
        FormatField { n: 0, ..f },
        FormatField {
            values: Typed {
                ty: BT_NULL,
                bytes: &[],
            },
            ..f
        },
    ] {
        assert_eq!(zero.size(), 0);
        assert!(zero.sample(usize::MAX).bytes.is_empty());
    }
    for (ty, width) in [(2, 2), (3, 4), (4, 8), (5, 4)] {
        let bytes = [1; 8];
        let typed = Typed {
            ty,
            bytes: &bytes[..width],
        };
        let index = usize::MAX / width + 1;
        assert!(std::panic::catch_unwind(|| typed.raw(index)).is_err());
        let field = FormatField {
            key: 0,
            n: 1,
            values: typed,
        };
        assert!(std::panic::catch_unwind(|| field.sample(index)).is_err());
        let field = FormatField { n: index, ..field };
        assert!(std::panic::catch_unwind(|| field.size()).is_err());
    }
}

#[test]
fn wide_format_keeps_first_duplicate_values_and_input_order() {
    for count in [8, 32, 33, 128, 255] {
        for duplicate in [false, true] {
            let mut text = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n".to_owned();
            for i in 0..count {
                text.push_str(&format!("##FORMAT=<ID=F{i},Number=1,Type=Integer>\n"));
            }
            text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\nchr1\t1\t.\tA\tC\t.\t.\t.\t");
            let names: Vec<_> = (0..count)
                .map(|i| {
                    if duplicate {
                        (count - 1 - i) % 7
                    } else {
                        count - 1 - i
                    }
                })
                .collect();
            for (i, name) in names.iter().enumerate() {
                if i != 0 {
                    text.push(':');
                }
                text.push_str(&format!("F{name}"));
            }
            for sample in 0..3 {
                text.push('\t');
                for i in 0..count {
                    if i != 0 {
                        text.push(':');
                    }
                    text.push_str(&(i + sample * 1000).to_string());
                }
            }
            text.push('\n');
            let mut reader = omicsweave_vcf::VariantReader::new(text.as_bytes()).unwrap();
            let mut rec = omicsweave_vcf::Record::default();
            assert!(reader.read(&mut rec).unwrap());
            let h = reader.header();
            let unique = if duplicate { count.min(7) } else { count };
            let expected_keys: Vec<_> = names[..unique]
                .iter()
                .map(|name| h.key_id(&format!("F{name}")).unwrap())
                .collect();
            assert_eq!(
                rec.format_fields().map(|f| f.key).collect::<Vec<_>>(),
                expected_keys
            );
            let mut out = Vec::new();
            for (i, name) in names[..unique].iter().enumerate() {
                assert_eq!(
                    rec.format_int32(h, &format!("F{name}"), &mut out).unwrap(),
                    3
                );
                assert_eq!(out, [i as i32, i as i32 + 1000, i as i32 + 2000]);
            }
        }
    }
}

#[test]
fn format_dedup_handles_bitset_last_bit_and_sparse_fallback() {
    for high_id in [511, 512, 4096] {
        let mut text = format!(
            "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n\
             ##FORMAT=<ID=F0,Number=1,Type=Integer,IDX={high_id}>\n"
        );
        for i in 1..33 {
            text.push_str(&format!(
                "##FORMAT=<ID=F{i},Number=1,Type=Integer,IDX={i}>\n"
            ));
        }
        text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t1\t.\tA\tC\t.\t.\t.\tF0");
        for i in 1..33 {
            text.push_str(&format!(":F{i}"));
        }
        text.push_str(":F0\t0");
        for i in 1..34 {
            text.push_str(&format!(":{i}"));
        }
        text.push('\n');
        let mut reader = omicsweave_vcf::VariantReader::new(text.as_bytes()).unwrap();
        let mut rec = omicsweave_vcf::Record::default();
        assert!(reader.read(&mut rec).unwrap());
        let keys: Vec<_> = std::iter::once(high_id).chain(1..33).collect();
        assert_eq!(rec.format_fields().map(|f| f.key).collect::<Vec<_>>(), keys);
        let mut out = Vec::new();
        for i in 0..33 {
            assert_eq!(
                rec.format_int32(reader.header(), &format!("F{i}"), &mut out)
                    .unwrap(),
                1
            );
            assert_eq!(out, [i]);
        }
    }
}
