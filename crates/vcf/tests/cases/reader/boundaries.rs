use crate::support::variants;

use std::io::{Read, Write};

use omicsweave_vcf::{Header, Record, VariantReader};

const COLUMNS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

fn vector(ty: u8, n: usize, bytes: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    if n < 15 {
        out.push((n as u8) << 4 | ty);
    } else {
        out.extend_from_slice(&[0xf0 | ty, 0x13]);
        out.extend_from_slice(&(n as i32).to_le_bytes());
    }
    out.extend_from_slice(bytes);
    out
}

fn bcf(
    version: &str,
    samples: usize,
    rlen: i32,
    alleles: &[&[u8]],
    info: &[(u8, Vec<u8>)],
    formats: &[(u8, Vec<u8>)],
) -> Vec<u8> {
    let mut header = format!(
        "##fileformat=VCFv{version}\n##contig=<ID=chr1,IDX=0>\n\
         ##FILTER=<ID=PASS,Description=\"Pass\",IDX=0>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\",IDX=1>\n\
         ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\",IDX=2>\n\
         ##INFO=<ID=END,Number=1,Type=Integer,Description=\"END\",IDX=3>\n\
         ##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"SVLEN\",IDX=4>\n\
         ##FORMAT=<ID=LEN,Number=1,Type=Integer,Description=\"LEN\",IDX=5>\n{COLUMNS}"
    );
    if samples > 0 {
        header.push_str("\tFORMAT");
        for i in 0..samples {
            header.push_str(&format!("\tS{i}"));
        }
    }
    header.push_str("\n\0");
    let mut out = b"BCF\x02\x02".to_vec();
    out.extend_from_slice(&(header.len() as u32).to_le_bytes());
    out.extend_from_slice(header.as_bytes());
    let mut shared = vector(7, 300, &vec![b'i'; 300]);
    for allele in alleles {
        shared.extend_from_slice(&vector(7, allele.len(), allele));
    }
    shared.push(0); // FILTER=.
    for (key, value) in info {
        shared.extend_from_slice(&[0x11, *key]);
        shared.extend_from_slice(value);
    }
    let mut indiv = Vec::new();
    for (key, value) in formats {
        indiv.extend_from_slice(&[0x11, *key]);
        indiv.extend_from_slice(value);
    }
    out.extend_from_slice(&(shared.len() as u32 + 24).to_le_bytes());
    out.extend_from_slice(&(indiv.len() as u32).to_le_bytes());
    out.extend_from_slice(&0i32.to_le_bytes()); // rid
    out.extend_from_slice(&9u32.to_le_bytes()); // pos
    out.extend_from_slice(&rlen.to_le_bytes());
    out.extend_from_slice(&0f32.to_le_bytes());
    out.extend_from_slice(&(info.len() as u16).to_le_bytes());
    out.extend_from_slice(&(alleles.len() as u16).to_le_bytes());
    out.extend_from_slice(&((samples as u32) | ((formats.len() as u32) << 24)).to_le_bytes());
    out.extend_from_slice(&shared);
    out.extend_from_slice(&indiv);
    out
}

fn oracle(bytes: &[u8]) {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("input.bcf");
    std::fs::write(&path, bytes).unwrap();
    assert_eq!(variants::ours(bytes), variants::htslib(&path));
}

#[test]
fn bcf_integer_widths_extended_vectors_and_sentinels_match_htslib() {
    for version in ["4.3", "4.4"] {
        for (ty, width, missing, end) in [
            (1, 1, i8::MIN as i32, i8::MIN as i32 + 1),
            (2, 2, i16::MIN as i32, i16::MIN as i32 + 1),
            (3, 4, i32::MIN, i32::MIN + 1),
        ] {
            for ploidy in [1, 2, 3, 20] {
                let mut values = Vec::new();
                for sample in 0..4 {
                    for j in 0..ploidy {
                        let v = match sample {
                            0 => {
                                if j == 0 {
                                    2
                                } else {
                                    5
                                }
                            }
                            1 => {
                                if j == 0 {
                                    missing
                                } else {
                                    5
                                }
                            }
                            2 => {
                                if j == 0 {
                                    end
                                } else {
                                    4
                                }
                            }
                            _ => {
                                if j == 0 {
                                    0
                                } else {
                                    end
                                }
                            }
                        };
                        values.extend_from_slice(&v.to_le_bytes()[..width]);
                    }
                }
                let gt = vector(ty, ploidy, &values);
                oracle(&bcf(version, 4, 1, &[b"A", b"C"], &[], &[(1, gt)]));
            }
        }
    }
}

#[test]
fn negative_bcf_rlen_uses_end_svlen_and_format_len() {
    let alleles: Vec<&[u8]> = vec![
        b"A",
        b"C",
        b"C",
        b"C",
        b"C",
        b"C",
        b"C",
        b"<DEL>",
        b"<DEL>",
        b"<NON_REF>",
    ];
    for rlen in [-1, i32::MIN] {
        let info = vec![
            (3, vector(3, 1, &30i32.to_le_bytes())),
            (
                4,
                vector(
                    2,
                    9,
                    &[5i16, 7, 11, 13, 17, 19, 23, 999, 31]
                        .into_iter()
                        .flat_map(i16::to_le_bytes)
                        .collect::<Vec<_>>(),
                ),
            ),
        ];
        let len = vector(
            3,
            1,
            &[40i32, 80]
                .into_iter()
                .flat_map(i32::to_le_bytes)
                .collect::<Vec<_>>(),
        );
        oracle(&bcf("4.5", 2, rlen, &alleles, &info, &[(5, len)]));
    }
    // htslib repairs rlen before unpacking: string lengths and symbolic
    // recognition therefore inspect encoded bytes, including NUL padding.
    oracle(&bcf(
        "4.5",
        1,
        -1,
        &[b"A\0\0", b"<DEL>\0", b"<NON_REF>\0"],
        &[(
            4,
            vector(
                2,
                2,
                &[100i16, 200]
                    .into_iter()
                    .flat_map(i16::to_le_bytes)
                    .collect::<Vec<_>>(),
            ),
        )],
        &[(5, vector(2, 1, &300i16.to_le_bytes()))],
    ));
}

#[test]
fn empty_indiv_and_zero_sample_compatibility_matches_htslib() {
    let mut empty = bcf(
        "4.4",
        1,
        1,
        &[b"A", b"C"],
        &[],
        &[(1, vector(1, 2, &[2, 4]))],
    );
    let header_len = u32::from_le_bytes(empty[5..9].try_into().unwrap()) as usize;
    let record = header_len + 9;
    let shared_len = u32::from_le_bytes(empty[record..record + 4].try_into().unwrap()) as usize;
    empty[record + 4..record + 8].copy_from_slice(&0u32.to_le_bytes());
    empty.truncate(record + 8 + shared_len);
    // htslib's bcf_read1_core deliberately repairs legacy bcf_subset output
    // by discarding n_fmt when the individual block or sample count is zero.
    oracle(&empty);
    let recs = read_all(&empty[..]);
    assert!(recs[0].genotypes().is_none());
    oracle(&bcf(
        "4.4",
        0,
        1,
        &[b"A", b"C"],
        &[],
        &[(1, vector(1, 2, &[2, 4]))],
    ));
}

struct Short<'a>(&'a [u8]);

impl Read for Short<'_> {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        let n = out.len().min(self.0.len()).min(1);
        out[..n].copy_from_slice(&self.0[..n]);
        self.0 = &self.0[n..];
        Ok(n)
    }
}

fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut z = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    z.write_all(bytes).unwrap();
    z.finish().unwrap()
}

fn read_all(source: impl Read) -> Vec<Record> {
    let mut reader = VariantReader::new(source).unwrap();
    let mut rec = Record::default();
    let mut records = Vec::new();
    while reader.read(&mut rec).unwrap() {
        records.push(rec.clone());
    }
    records
}

#[test]
fn short_reads_and_concatenated_gzip_members_preserve_detection_and_records() {
    let text = format!("##fileformat=VCFv4.3\n{COLUMNS}\nchr1\t10\t.\tA\tC\t.\t.\t.\n");
    let binary = bcf("4.3", 0, 1, &[b"A", b"C"], &[], &[]);
    for bytes in [text.as_bytes(), &binary] {
        let want = read_all(bytes);
        assert_eq!(read_all(Short(bytes)), want);
        let mut compressed = Vec::new();
        for part in bytes.chunks(7) {
            compressed.extend_from_slice(&gzip(part));
        }
        assert_eq!(read_all(Short(&compressed)), want);
    }
}

#[test]
fn malformed_record_clears_output_and_bcf_errors_end_stream() {
    let mut bytes = bcf(
        "4.4",
        1,
        1,
        &[b"A", b"C"],
        &[],
        &[
            (1, vector(1, 2, &[2, 4])),
            (2, vector(3, 1, &5i32.to_le_bytes())),
        ],
    );
    bytes.pop(); // Truncation after a fully valid GT.
    let mut reader = VariantReader::new(&bytes[..]).unwrap();
    let mut rec = Record::default();
    assert!(reader.read(&mut rec).is_err());
    assert_eq!(rec, Record::default());
    assert!(!reader.read(&mut rec).unwrap());

    let text = format!(
        "##fileformat=VCFv4.4\n##FORMAT=<ID=GT,Number=1,Type=String>\n{COLUMNS}\tFORMAT\tS\n\
        chr1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/1\n\
        chr1\t11\t.\tA\tC\t.\t.\t.\tGT\t0/q\n\
        chr1\t12\t.\tA\tC\t.\t.\t.\tGT\t1\n"
    );
    let mut reader = VariantReader::new(text.as_bytes()).unwrap();
    assert!(reader.read(&mut rec).unwrap());
    assert!(reader.read(&mut rec).is_err());
    assert_eq!(rec, Record::default());
    assert!(reader.read(&mut rec).unwrap());
    assert_eq!(rec.pos, 11);
}

#[test]
fn invalid_binary_descriptors_and_ignored_fields_are_rejected() {
    let cases = [
        // Lengths must be scalar integers, with no recursive extension.
        bcf("4.4", 0, 1, &[b"A"], &[(3, vec![0xf3, 0xf1, 0x11, 1])], &[]),
        // A null INFO value may not claim a nonzero element count.
        bcf("4.4", 0, 1, &[b"A"], &[(3, vec![0x10])], &[]),
        // A malformed FORMAT after a valid GT must not be ignored.
        bcf(
            "4.4",
            1,
            1,
            &[b"A"],
            &[],
            &[(1, vector(1, 1, &[2])), (2, vec![0x16])],
        ),
        bcf("4.4", 0, 1, &[], &[], &[]),
        bcf("4.4", 0, 1, &[b"A"], &[(127, vector(1, 1, &[1]))], &[]),
    ];
    for bytes in cases {
        let mut reader = VariantReader::new(&bytes[..]).unwrap();
        assert!(reader.read(&mut Record::default()).is_err());
    }
}

#[test]
fn huge_truncated_bcf_lengths_do_not_allocate_claimed_sizes() {
    let mut bytes = b"BCF\x02\x02".to_vec();
    bytes.extend_from_slice(&u32::MAX.to_le_bytes());
    assert!(VariantReader::new(&bytes[..]).is_err());

    let mut bytes = bcf("4.4", 0, 1, &[b"A"], &[], &[]);
    let header_len = u32::from_le_bytes(bytes[5..9].try_into().unwrap()) as usize;
    bytes[header_len + 9..header_len + 13].copy_from_slice(&u32::MAX.to_le_bytes());
    let mut reader = VariantReader::new(&bytes[..]).unwrap();
    assert!(reader.read(&mut Record::default()).is_err());
}

#[test]
fn gzip_crc_and_truncated_trailers_are_errors() {
    let text = format!("##fileformat=VCFv4.3\n{COLUMNS}\nchr1\t10\t.\tA\tC\t.\t.\t.\n");
    let valid = gzip(text.as_bytes());
    for n in [1, 4, 8] {
        let bytes = &valid[..valid.len() - n];
        let mut reader = VariantReader::new(bytes).unwrap();
        let mut rec = Record::default();
        let result = loop {
            match reader.read(&mut rec) {
                Ok(true) => {}
                result => break result,
            }
        };
        assert!(result.is_err());
        assert_eq!(rec, Record::default());
        assert!(!reader.read(&mut rec).unwrap());
    }
    let mut corrupt = valid;
    let crc_at = corrupt.len() - 8;
    corrupt[crc_at] ^= 1;
    let mut reader = VariantReader::new(&corrupt[..]).unwrap();
    let mut rec = Record::default();
    let result = loop {
        match reader.read(&mut rec) {
            Ok(true) => {}
            result => break result,
        }
    };
    assert!(result.is_err());
    assert_eq!(rec, Record::default());
    assert!(!reader.read(&mut rec).unwrap());
}

#[test]
fn metadata_boundaries_are_checked_without_arithmetic_panics() {
    let huge_version = format!("##fileformat=VCFv2147483647.2147483647\n{COLUMNS}\n");
    assert_eq!(Header::parse(&huge_version).unwrap().version, 4_002_000);
    assert!(Header::parse("#not-a-column-header\n").is_err());
    for meta in [
        "##contig=<ID=chr1,IDX=4294967295>",
        "##FILTER=<ID=A,IDX=0>\n##FILTER=<ID=B,IDX=0>",
        "##contig=<ID=chr1,IDX=0>\n##contig=<ID=chr2,IDX=0>",
    ] {
        assert!(Header::parse(&format!("{meta}\n{COLUMNS}\n")).is_err());
    }
    let sparse = Header::parse(&format!("##contig=<ID=chr1,IDX=100>\n{COLUMNS}\n")).unwrap();
    assert_eq!(sparse.contig_id("chr1"), Some(100));
}

#[test]
fn invalid_utf8_identifiers_error_without_aliasing_and_next_line_is_readable() {
    let mut bytes = format!(
        "##fileformat=VCFv4.4\n##FORMAT=<ID=GT,Number=1,Type=String>\n{COLUMNS}\tFORMAT\tS\n"
    )
    .into_bytes();
    let malformed: [&[u8]; 6] = [
        b"chr\xff\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\n",
        b"chr\xfe\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\n",
        b"chr1\t1\t.\tA\tC\t.\t.\tX\xff=1\tGT\t0/1\n",
        b"chr1\t1\t.\tA\tC\t.\t.\tX\xfe=1\tGT\t0/1\n",
        b"chr1\t1\t.\tA\tC\t.\t.\t.\tX\xff\t1\n",
        b"chr1\t1\t.\tA\tC\t.\t.\t.\tX\xfe\t1\n",
    ];
    for line in malformed {
        bytes.extend_from_slice(line);
        bytes.extend_from_slice(b"chr1\t2\t.\tA\tC\t.\t.\t.\tGT\t1/1\n");
    }
    let mut reader = VariantReader::new(&bytes[..]).unwrap();
    let mut rec = Record::default();
    for _ in malformed {
        assert!(matches!(
            reader.read(&mut rec),
            Err(omicsweave_core::Error::Format { .. })
        ));
        assert_eq!(rec, Record::default());
        assert!(reader.read(&mut rec).unwrap());
        assert_eq!(rec.pos, 1);
        assert!(reader.header().contig_id("chr\u{fffd}").is_none());
        assert!(reader.header().key_id("X\u{fffd}").is_none());
    }
    assert!(!reader.read(&mut rec).unwrap());
}

#[test]
fn vcf_and_bcf_headers_reject_invalid_utf8() {
    let mut text = b"##fileformat=VCFv4.4\n##source=invalid-\xff\n".to_vec();
    text.extend_from_slice(COLUMNS.as_bytes());
    text.push(b'\n');
    assert!(VariantReader::new(&text[..]).is_err());

    let mut binary = bcf("4.4", 0, 1, &[b"A", b"C"], &[], &[]);
    let contig_at = binary.windows(4).position(|w| w == b"chr1").unwrap();
    binary[contig_at] = 0xff;
    assert!(VariantReader::new(&binary[..]).is_err());
}

#[test]
fn unicode_identifiers_and_raw_vcf_value_bytes_are_preserved() {
    let mut bytes = format!(
        "##fileformat=VCFv4.4\n##contig=<ID=染色体α>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String>\n{COLUMNS}\tFORMAT\t样本β\n\
         染色体α\t2\t"
    )
    .into_bytes();
    bytes.extend_from_slice(b"id\xfc\tA\xff\tC\xfe\t.\t.\t");
    bytes.extend_from_slice("注释=".as_bytes());
    bytes.extend_from_slice(b"\xfd\tGT:");
    bytes.extend_from_slice("格式".as_bytes());
    bytes.extend_from_slice(b"\t0/1:\xfe\n");
    let mut reader = VariantReader::new(&bytes[..]).unwrap();
    assert_eq!(reader.header().samples, ["样本β"]);
    let mut rec = Record::default();
    assert!(reader.read(&mut rec).unwrap());
    assert_eq!(reader.header().contigs[rec.rid as usize].name, "染色体α");
    assert!(reader.header().key_id("注释").is_some());
    assert!(reader.header().key_id("格式").is_some());
    assert_eq!(rec.allele(0), b"A\xff");
    assert_eq!(rec.allele(1), b"C\xfe");
    assert_eq!(rec.genotypes().unwrap().values, [2, 4]);
    assert!(!reader.read(&mut rec).unwrap());
}

#[test]
fn unicode_bcf_header_and_raw_alleles_are_preserved() {
    let original = bcf("4.4", 0, 2, &[b"A\xff", b"C\xfe"], &[], &[]);
    let old_len = u32::from_le_bytes(original[5..9].try_into().unwrap()) as usize;
    let header = std::str::from_utf8(&original[9..9 + old_len])
        .unwrap()
        .replace("chr1", "染色体α");
    let mut bytes = original[..5].to_vec();
    bytes.extend_from_slice(&(header.len() as u32).to_le_bytes());
    bytes.extend_from_slice(header.as_bytes());
    bytes.extend_from_slice(&original[9 + old_len..]);
    let mut reader = VariantReader::new(&bytes[..]).unwrap();
    assert_eq!(reader.header().contigs[0].name, "染色体α");
    let mut rec = Record::default();
    assert!(reader.read(&mut rec).unwrap());
    assert_eq!(rec.allele(0), b"A\xff");
    assert_eq!(rec.allele(1), b"C\xfe");
}
