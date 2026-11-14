//! Malformed input and invalid requests yield errors, never panics.

use omicsweave_fasta::{
    BgzfWriter, BlockHeader, Error, Fai, FaiEntry, Gzi, Inflater, MAX_BLOCK_LEN, Parsed, Record,
    TwoBitHeader,
};

fn entry() -> FaiEntry {
    FaiEntry {
        name: "s".into(),
        len: 10,
        offset: 3,
        line_bases: 4,
        line_width: 5,
    }
}

#[test]
fn fai_rejects_ranges_outside_the_sequence() {
    let e = entry();
    assert!(matches!(
        e.text_range(5, 11),
        Err(Error::RegionOutOfRange { .. })
    ));
    assert!(matches!(
        e.text_range(6, 5),
        Err(Error::InvalidRegion { .. })
    ));
    assert_eq!(e.text_range(4, 9).unwrap(), 8..14);
    assert_eq!(e.text_range(7, 7).unwrap(), 3..3);
}

#[test]
fn fai_decode_checks_the_line_layout() {
    let e = entry();
    let text = b">s\nACGT\nACGT\nAC\n";
    let mut out = [0u8; 10];
    e.decode(0, &text[3..], &mut out).unwrap();
    assert_eq!(&out, b"ACGTACGTAC");

    let shifted = b">s\nACGTA\nCGT\nAC\n";
    assert!(matches!(
        e.decode(0, &shifted[3..], &mut out),
        Err(Error::Invalid { .. })
    ));
    assert!(matches!(
        e.decode(0, &text[3..8], &mut out),
        Err(Error::Truncated { .. })
    ));
}

#[test]
fn fai_parse_rejects_bad_rows() {
    assert!(Fai::parse(b"s\t10\t3\t4\n").is_err());
    assert!(Fai::parse(b"s\t10\t3\t4\t5\t9\n").is_err());
    assert!(Fai::parse(b"s\tten\t3\t4\t5\n").is_err());
    assert!(Fai::parse(b"s\t10\t3\t5\t4\n").is_err());
    assert_eq!(
        Fai::parse(b"s\t10\t3\t4\t5\n").unwrap().get("s"),
        Some(&entry())
    );
}

#[test]
fn bgzf_detects_corruption() {
    let mut w = BgzfWriter::new(Vec::new());
    std::io::Write::write_all(&mut w, b"ACGTACGTACGT").unwrap();
    let (mut bytes, _) = w.finish().unwrap();
    let mut out = vec![0u8; MAX_BLOCK_LEN];
    let mut inflater = Inflater::new();
    assert_eq!(inflater.inflate_block(&bytes, &mut out).unwrap().1, 12);

    let header = BlockHeader::parse(&bytes).unwrap();
    bytes[header.block_len - 8] ^= 1;
    assert!(matches!(
        inflater.inflate_block(&bytes, &mut out),
        Err(Error::Inflate(_))
    ));
    assert!(matches!(
        inflater.inflate_block(&bytes[..10], &mut out),
        Err(Error::Truncated { .. })
    ));
    assert!(BlockHeader::parse(b"not a gzip block").is_err());
}

#[test]
fn gzi_rejects_unordered_and_short_input() {
    let mut b = 2u64.to_le_bytes().to_vec();
    for v in [100u64, 1000, 50, 2000] {
        b.extend_from_slice(&v.to_le_bytes());
    }
    assert!(matches!(Gzi::parse(&b), Err(Error::Invalid { .. })));
    assert!(matches!(Gzi::parse(&b[..20]), Err(Error::Truncated { .. })));
}

#[test]
fn gzi_bounds_each_block_before_reads_are_planned() {
    let encode = |points: &[(u64, u64)]| {
        let mut bytes = (points.len() as u64).to_le_bytes().to_vec();
        for &(compressed, text) in points {
            bytes.extend(compressed.to_le_bytes());
            bytes.extend(text.to_le_bytes());
        }
        bytes
    };
    let limit = MAX_BLOCK_LEN as u64;
    for points in [vec![(1, limit + 1)], vec![(28, 0), (29, limit + 1)]] {
        assert!(matches!(
            Gzi::parse(&encode(&points)),
            Err(Error::Invalid { .. })
        ));
    }
    let points = [(28, 0), (29 + limit, limit), (57 + limit, limit)];
    let index = Gzi::parse(&encode(&points)).unwrap();
    assert_eq!(&index.starts()[1..], &points);
    assert_eq!(Gzi::parse(&index.to_bytes()).unwrap(), index);
}

#[test]
fn twobit_header_and_records() {
    let mut h = 0x1A41_2743u32.to_be_bytes().to_vec();
    h.extend_from_slice(&[0u8; 12]);
    assert!(matches!(TwoBitHeader::parse(&h), Err(Error::BigEndian)));

    // 8 bases, one N block that runs past the end.
    let mut r = Vec::new();
    for v in [8u32, 1, 6, 4, 0, 0] {
        r.extend_from_slice(&v.to_le_bytes());
    }
    assert_eq!(Record::parse(&r[..4], 0).unwrap(), Parsed::Need(8));
    assert_eq!(Record::parse(&r[..8], 0).unwrap(), Parsed::Need(20));
    assert_eq!(Record::parse(&r[..20], 0).unwrap(), Parsed::Need(24));
    assert!(matches!(Record::parse(&r, 0), Err(Error::Invalid { .. })));
}

#[test]
fn readers_check_keys_and_outputs() {
    use omicsweave_core::{Region, SequenceRead};
    use omicsweave_fasta::FastaReader;

    let text = b">s\nACGT\nACGT\nAC\n";
    let fai = Fai::parse(b"s\t10\t3\t4\t5\n").unwrap();
    let mut reader = FastaReader::new(&text[..], fai, None);
    let mut out = [0u8; 4];
    reader.sequence(Region::new("s", 2, 6), &mut out).unwrap();
    assert_eq!(&out, b"GTAC");
    assert_eq!(reader.chrom_len("t"), Err(Error::UnknownChrom("t".into())));
    assert!(matches!(
        reader.sequence(Region::new("s", 8, 12), &mut out),
        Err(Error::RegionOutOfRange { .. })
    ));
    assert!(matches!(
        reader.sequence(Region::new("s", 0, 3), &mut out),
        Err(Error::OutputSize { .. })
    ));
}
