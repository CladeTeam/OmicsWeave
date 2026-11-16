use crate::support::{Hts, fixture, open, query_all, read_all};
use omicsweave_bam::Record;
use omicsweave_core::{Region, exec::Inline};

#[test]
fn records_with_long_cigars_use_their_cg_tag() {
    let f = fixture(41, 4_000, 30, 6);
    let all = read_all(&f.bam);
    assert!(
        all.iter().any(|r| r.cigar.len() > 65_535),
        "fixture has a long CIGAR"
    );
    let raw = std::fs::read(&f.bam).unwrap();
    let mut scan = omicsweave_bam::Scanner::new(&raw[..], Box::new(Inline::new()));
    let header = scan.header().unwrap();
    let n_refs = header.references().len();
    let mut i = 0;
    while let Some(range) = scan.next_record().unwrap() {
        let rec = Record::parse(scan.slice(range), n_refs).unwrap();
        let want = &all[i];
        assert_eq!(rec.n_cigar(), want.cigar.len());
        assert_eq!(rec.end(), want.end);
        let (a, b) = rec.aux();
        assert_eq!([a, b].concat(), want.aux);
        i += 1;
    }
    assert_eq!(i, all.len());
}

#[test]
fn long_cigar_and_empty_sequence_outputs_match_htslib() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("long.sam");
    let bam = dir.path().join("long.bam");
    let bai = dir.path().join("long.bam.bai");
    let cigar = "1M1I".repeat(33_000) + "1M";
    let sequence = "AC".repeat(33_000) + "G";
    std::fs::write(
        &sam,
        format!(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\nlong\t0\tchr1\t1\t60\t{cigar}\t*\t0\t0\t{sequence}\t*\tNM:i:33000\nmissing\t0\tchr1\t2\t60\t1M\t*\t0\t0\t*\t*\n"
        ),
    )
    .unwrap();
    crate::support::bam_of(&sam, &bam, 0);
    crate::support::index(&bam, &bai, 0);
    let hts = Hts::open(&bam, &bai);
    let expected = hts.query(0, 0, 100_000);
    assert_eq!(expected[0].cigar.len(), 66_001);
    assert!(expected[1].seq.is_empty());
    let mut reader = open(&bam, &bai);
    assert_eq!(
        query_all(&mut reader, Region::new("chr1", 0, 100_000)),
        expected
    );

    let raw = std::fs::read(&bam).unwrap();
    let mut scan = omicsweave_bam::Scanner::new(&raw[..], Box::new(Inline::new()));
    scan.header().unwrap();
    let range = scan.next_record().unwrap().unwrap();
    let data = scan.slice(range);
    for cut in [0, 31, 32, 33, data.len() / 2, data.len() - 1] {
        assert!(Record::parse(&data[..cut], 1).is_err(), "cut {cut}");
    }
}
