//! `BedTable` queries against a brute-force scan of the generated records,
//! for plain text and every gzip form.

use std::ffi::CString;
use std::io::Write;

use omicsweave_bed::BedTable;
use omicsweave_core::{
    Error, IntervalEntry, IntervalRead, IntervalSize, IntervalsOut, Region, StringColumn, Strings,
};
use proptest::prelude::*;

#[test]
fn chromosome_names_are_lossless_utf8() {
    assert!(BedTable::load(b"\xff\t0\t10\n\xfe\t20\t30\n").is_err());
    let table = BedTable::load("染色体\t0\t10\n".as_bytes()).unwrap();
    assert_eq!(table.chroms(), &["染色体"]);
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn range(&mut self, lo: u64, hi: u64) -> u64 {
        lo + self.next() % (hi - lo + 1)
    }
    fn chance(&mut self, percent: u64) -> bool {
        self.next() % 100 < percent
    }
}

#[derive(Debug, Clone)]
struct Rec {
    chrom: String,
    start: u32,
    end: u32,
    /// Fields after the third, tab-joined.
    rest: String,
}

fn random_records(rng: &mut Rng) -> Vec<Rec> {
    let mut recs = Vec::new();
    for _ in 0..rng.range(1, 600) {
        let start = rng.range(0, 50_000) as u32;
        let width = if rng.chance(3) {
            0
        } else {
            rng.range(1, 3_000) as u32
        };
        let fields = rng.range(0, 4);
        let rest: Vec<String> = (0..fields)
            .map(|f| format!("f{f}_{}", rng.next() % 1000))
            .collect();
        recs.push(Rec {
            chrom: format!("chr{}", rng.range(1, 4)),
            start,
            end: start + width,
            rest: rest.join("\t"),
        });
    }
    if rng.chance(50) {
        recs.sort_by(|a, b| (&a.chrom, a.start).cmp(&(&b.chrom, b.start)));
    }
    recs
}

fn render(recs: &[Rec], rng: &mut Rng) -> Vec<u8> {
    let eol: &[u8] = if rng.chance(20) { b"\r\n" } else { b"\n" };
    let spaces = rng.chance(20);
    let mut out = Vec::new();
    if rng.chance(50) {
        out.extend_from_slice(b"track name=t");
        out.extend_from_slice(eol);
    }
    for (i, r) in recs.iter().enumerate() {
        if rng.chance(5) {
            out.extend_from_slice(b"# note");
            out.extend_from_slice(eol);
        }
        if rng.chance(5) {
            out.extend_from_slice(eol);
        }
        let sep = if spaces { " " } else { "\t" };
        let mut line = format!("{}{sep}{}{sep}{}", r.chrom, r.start, r.end);
        if !r.rest.is_empty() {
            line.push_str(sep);
            line.push_str(&if spaces {
                r.rest.replace('\t', " ")
            } else {
                r.rest.clone()
            });
        }
        out.extend_from_slice(line.as_bytes());
        if i + 1 < recs.len() || rng.chance(90) {
            out.extend_from_slice(eol);
        }
    }
    out
}

fn gzip(text: &[u8], level: u32) -> Vec<u8> {
    let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(level));
    enc.write_all(text).unwrap();
    enc.finish().unwrap()
}

/// `text` as bgzip, written by htslib.
fn bgzip(text: &[u8]) -> Vec<u8> {
    use htslib_sys as ffi;
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.bed.gz");
    let p = CString::new(path.to_str().unwrap()).unwrap();
    // SAFETY: valid path and buffer; the handle is closed before reading.
    unsafe {
        let fp = ffi::bgzf_open(p.as_ptr(), c"w".as_ptr());
        assert!(!fp.is_null());
        assert_eq!(
            ffi::bgzf_write(fp, text.as_ptr().cast(), text.len()),
            text.len() as isize
        );
        assert_eq!(ffi::bgzf_close(fp), 0);
    }
    std::fs::read(path).unwrap()
}

/// A query's output, owned.
#[derive(Debug, PartialEq)]
struct Answer {
    entries: Vec<IntervalEntry>,
    annotations: Vec<Vec<u8>>,
}

fn query(table: &mut BedTable, region: Region<'_>) -> Answer {
    let size = table.query(region).unwrap();
    write_answer(table, size)
}

fn write_answer(table: &BedTable, size: IntervalSize) -> Answer {
    let mut entries = vec![IntervalEntry::default(); size.entries];
    let mut offsets = vec![0i64; size.entries + 1];
    let mut data = vec![0u8; size.annotation_bytes];
    table
        .write(IntervalsOut {
            entries: &mut entries,
            annotations: StringColumn {
                offsets: &mut offsets,
                data: &mut data,
            },
        })
        .unwrap();
    let strings = Strings {
        offsets: &offsets,
        data: &data,
    };
    let annotations = entries
        .iter()
        .map(|e| strings.get(e.annotation as usize).to_vec())
        .collect();
    Answer {
        entries,
        annotations,
    }
}

fn brute_force(recs: &[Rec], region: Region<'_>) -> Answer {
    let mut hits: Vec<(usize, &Rec)> = recs
        .iter()
        .enumerate()
        .filter(|(_, r)| {
            !region.is_empty() && r.chrom == region.chrom && region.overlaps(r.start, r.end)
        })
        .collect();
    hits.sort_by_key(|&(i, r)| (r.start, i));
    Answer {
        entries: hits
            .iter()
            .enumerate()
            .map(|(row, (_, r))| IntervalEntry {
                start: r.start,
                end: r.end,
                annotation: row as u32,
            })
            .collect(),
        annotations: hits
            .iter()
            .map(|(_, r)| r.rest.as_bytes().to_vec())
            .collect(),
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 64, ..ProptestConfig::default() })]

    #[test]
    fn queries_match_a_scan(seed in any::<u64>()) {
        let mut rng = Rng(seed);
        let recs = random_records(&mut rng);
        let text = render(&recs, &mut rng);
        let mut members = gzip(&text[..text.len() / 2], 6);
        members.extend_from_slice(&gzip(&text[text.len() / 2..], 1));
        let forms = [text.clone(), gzip(&text, rng.range(0, 9) as u32), members, bgzip(&text)];

        let mut regions: Vec<(String, u32, u32)> = (1..=4).map(|c| (format!("chr{c}"), 0, 60_000)).collect();
        regions.push(("chrUn".into(), 0, 100));
        for _ in 0..40 {
            let start = rng.range(0, 55_000) as u32;
            regions.push((format!("chr{}", rng.range(1, 4)), start, start + rng.range(0, 8_000) as u32));
        }

        for raw in &forms {
            let mut table = BedTable::load(raw).unwrap();
            prop_assert_eq!(table.len(), recs.len());
            for (chrom, start, end) in &regions {
                let region = Region::new(chrom, *start, *end);
                prop_assert_eq!(query(&mut table, region), brute_force(&recs, region), "{}:{}-{}", chrom, start, end);
            }
        }
    }
}

#[test]
fn rejects_bad_lines_and_requests() {
    let line = |text: &[u8]| match BedTable::load(text) {
        Err(Error::Format { line, .. }) => line,
        other => panic!("expected a format error, got {other:?}"),
    };
    assert_eq!(line(b"chr1\t1\t2\nchr1\t5\t3\n"), 2);
    assert_eq!(line(b"# c\nchr1\t1\n"), 2);
    assert_eq!(line(b"chr1\tx\t2\n"), 1);
    assert_eq!(line(b"chr1 1\n"), 1);
    assert!(BedTable::load(b"\x1f\x8bnot gzip").is_err());

    let mut table = BedTable::load(b"chr1\t10\t20\tname\n").unwrap();
    assert!(matches!(
        table.query(Region::new("chr1", 5, 2)),
        Err(Error::InvalidRegion { .. })
    ));
    let mut entries = vec![IntervalEntry::default(); 2];
    let (mut offsets, mut data) = (vec![0i64; 3], vec![0u8; 4]);
    // A failed query stages nothing; writing then is an error.
    assert_eq!(
        table.write(out(&mut entries, &mut offsets, &mut data)),
        Err(Error::NotQueried)
    );
    // A buffer sized for another result is rejected.
    assert_eq!(table.query(Region::new("chr1", 0, 100)).unwrap().entries, 1);
    assert!(matches!(
        table.write(out(&mut entries, &mut offsets, &mut data)),
        Err(Error::OutputSize { .. })
    ));
}

fn out<'a>(
    entries: &'a mut [IntervalEntry],
    offsets: &'a mut [i64],
    data: &'a mut [u8],
) -> IntervalsOut<'a> {
    IntervalsOut {
        entries,
        annotations: StringColumn { offsets, data },
    }
}

#[test]
fn write_copies_the_latest_query() {
    let mut table = BedTable::load(b"chr1\t10\t20\ta\nchr1\t30\t40\tbb\n").unwrap();
    table.query(Region::new("chr1", 0, 15)).unwrap();
    let size = table.query(Region::new("chr1", 25, 35)).unwrap();
    assert_eq!((size.entries, size.annotation_bytes), (1, 2));
    let mut entries = vec![IntervalEntry::default(); 1];
    let (mut offsets, mut data) = (vec![0i64; 2], vec![0u8; 2]);
    let out = IntervalsOut {
        entries: &mut entries,
        annotations: StringColumn {
            offsets: &mut offsets,
            data: &mut data,
        },
    };
    table.write(out).unwrap();
    assert_eq!(
        (entries[0].start, entries[0].end, &data[..]),
        (30, 40, &b"bb"[..])
    );
}

#[test]
fn nested_intervals_keep_stable_order_and_annotations() {
    let mut recs = vec![Rec {
        chrom: "chr1".into(),
        start: 0,
        end: 100_000,
        rest: "outer".into(),
    }];
    for i in 0..4096 {
        recs.push(Rec {
            chrom: "chr1".into(),
            start: i * 10,
            end: i * 10 + 3,
            rest: format!("r{i}"),
        });
    }
    recs.extend([
        Rec {
            chrom: "chr1".into(),
            start: 40_000,
            end: 40_020,
            rest: "tie".into(),
        },
        Rec {
            chrom: "chr1".into(),
            start: 40_001,
            end: 40_001,
            rest: "empty".into(),
        },
    ]);
    let mut text = Vec::new();
    for r in &recs {
        writeln!(text, "{}\t{}\t{}\t{}", r.chrom, r.start, r.end, r.rest).unwrap();
    }
    let mut table = BedTable::load(&text).unwrap();
    for (s, e) in [
        (0, 0),
        (0, 1),
        (100, 100),
        (39_999, 40_010),
        (40_001, 40_002),
        (100_000, 100_001),
    ] {
        let region = Region::new("chr1", s, e);
        assert_eq!(query(&mut table, region), brute_force(&recs, region));
    }
    assert!(table.query(Region::new("chr1", 3, 2)).is_err());
    let mut entries = [];
    let mut offsets = [0];
    let mut data = [];
    assert_eq!(
        table.write(IntervalsOut {
            entries: &mut entries,
            annotations: StringColumn {
                offsets: &mut offsets,
                data: &mut data
            }
        }),
        Err(Error::NotQueried)
    );
}

#[test]
fn large_unsorted_groups_match_stable_scan_across_coordinate_bytes() {
    for mask in [0x0000_00FFu32, 0xFF00_FF00, 0xFFFF_FFFF, 0x8000_FFFF] {
        let mut recs = Vec::new();
        let mut text = Vec::new();
        for i in 0..16_384u32 {
            let mut x = (i / 4).wrapping_mul(0x9E37_79B9);
            x ^= x >> 16;
            let start = x & mask;
            let record = Rec {
                chrom: format!("chr{}", i % 2),
                start,
                end: start.saturating_add(i % 23),
                rest: format!("row{i}"),
            };
            writeln!(
                text,
                "{}\t{}\t{}\t{}",
                record.chrom, record.start, record.end, record.rest
            )
            .unwrap();
            recs.push(record);
        }
        let mut table = BedTable::load(&text).unwrap();
        for chrom in ["chr0", "chr1"] {
            for (start, end) in [(0, u32::MAX), (0, 128), (0x8000_0000, 0x9000_0000)] {
                let region = Region::new(chrom, start, end);
                assert_eq!(query(&mut table, region), brute_force(&recs, region));
            }
        }
    }
}

#[test]
fn staged_spans_and_tree_hits_match_scan_through_reuse() {
    // Rows are interleaved, duplicate starts retain file order, annotation
    // lengths differ, and chrB has nonmonotone ends.
    let recs: Vec<_> = (0..512u32)
        .flat_map(|i| {
            [
                Rec {
                    chrom: "chrA".into(),
                    start: i / 3 * 10,
                    end: i / 3 * 10 + 5,
                    rest: "a".repeat((i % 13) as usize),
                },
                Rec {
                    chrom: "chrB".into(),
                    start: i * 10,
                    end: if i == 0 { 100_000 } else { i * 10 + 5 },
                    rest: format!("b{i}"),
                },
                Rec {
                    chrom: "chrC".into(),
                    start: 5120 - i * 10,
                    end: 5125 - i * 10,
                    rest: format!("c{i}"),
                },
            ]
        })
        .collect();
    let mut text = Vec::new();
    for r in &recs {
        writeln!(text, "{}\t{}\t{}\t{}", r.chrom, r.start, r.end, r.rest).unwrap();
    }
    let mut table = BedTable::load(&text).unwrap();
    for region in [
        Region::new("chrA", 0, 10_000),
        Region::new("chrB", 4900, 4901),
        Region::new("chrA", 130, 141),
        Region::new("chrC", 0, 10_000),
        Region::new("chrB", 0, 10_000),
        Region::new("unknown", 0, 1),
        Region::new("chrA", 0, 0),
        Region::new("chrA", 0, 10_000),
    ] {
        assert_eq!(query(&mut table, region), brute_force(&recs, region));
        let cloned = table.clone();
        let expected = brute_force(&recs, region);
        let size = IntervalSize {
            entries: expected.entries.len(),
            annotation_bytes: expected.annotations.iter().map(Vec::len).sum(),
        };
        assert_eq!(write_answer(&cloned, size), expected);
        assert!(table.query(Region::new("chrA", 3, 2)).is_err());
        let (mut entries, mut offsets, mut data) = ([], [99], []);
        assert_eq!(
            table.write(out(&mut entries, &mut offsets, &mut data)),
            Err(Error::NotQueried)
        );
        assert_eq!(offsets, [99]);
    }
}

#[test]
fn empty_and_default_tables_write_an_empty_column() {
    for mut table in [BedTable::default(), BedTable::load(b"").unwrap()] {
        assert_eq!(
            query(&mut table, Region::new("absent", 0, 1)),
            Answer {
                entries: vec![],
                annotations: vec![]
            }
        );
        assert_eq!(
            query(&mut table, Region::new("absent", u32::MAX, u32::MAX)),
            Answer {
                entries: vec![],
                annotations: vec![]
            }
        );
    }
}

#[test]
fn contiguous_rows_write_exactly_with_output_guards() {
    let recs: Vec<_> = (0..64u32)
        .map(|i| Rec {
            chrom: "chr1".into(),
            start: i * 10,
            end: i * 10 + 5,
            rest: "x".repeat((i % 11) as usize),
        })
        .collect();
    let mut text = Vec::new();
    for r in &recs {
        writeln!(text, "{}\t{}\t{}\t{}", r.chrom, r.start, r.end, r.rest).unwrap();
    }
    let mut table = BedTable::load(&text).unwrap();
    for region in [
        Region::new("chr1", 0, 640),
        Region::new("chr1", 205, 211),
        Region::new("chr1", 0, 0),
    ] {
        let expected = brute_force(&recs, region);
        let size = table.query(region).unwrap();
        let guard = IntervalEntry {
            start: u32::MAX,
            end: u32::MAX,
            annotation: u32::MAX,
        };
        let mut entries = vec![guard; size.entries + 2];
        let mut offsets = vec![-1; size.entries + 3];
        let mut data = vec![255; size.annotation_bytes + 2];
        table
            .write(out(
                &mut entries[1..size.entries + 1],
                &mut offsets[1..size.entries + 2],
                &mut data[1..size.annotation_bytes + 1],
            ))
            .unwrap();
        assert_eq!(entries[0], guard);
        assert_eq!(entries[size.entries + 1], guard);
        assert_eq!(offsets[0], -1);
        assert_eq!(offsets[size.entries + 2], -1);
        assert_eq!(data[0], 255);
        assert_eq!(data[size.annotation_bytes + 1], 255);
        assert_eq!(&entries[1..size.entries + 1], expected.entries);
        let strings = Strings {
            offsets: &offsets[1..size.entries + 2],
            data: &data[1..size.annotation_bytes + 1],
        };
        for (row, annotation) in expected.annotations.iter().enumerate() {
            assert_eq!(strings.get(row), annotation);
        }
    }
}

#[test]
fn gzip_preserves_long_annotations_across_member_boundaries_and_rejects_corruption() {
    let note = vec![b'x'; (1 << 18) + 17];
    let mut text = "染色体\t0\t10\t".as_bytes().to_vec();
    text.extend_from_slice(&note);
    text.extend_from_slice(b"\r\nchr2\t1\t2\ty\n");
    for split in [1, 5, text.len() - 15, text.len() - 14] {
        let mut raw = gzip(&text[..split], 6);
        raw.extend_from_slice(&gzip(&text[split..], 6));
        let mut table = BedTable::load(&raw).unwrap();
        let answer = query(&mut table, Region::new("染色体", 0, 10));
        assert_eq!(answer.annotations.as_slice(), std::slice::from_ref(&note));
    }
    let raw = gzip(&text, 6);
    for missing in [1, 4, 8] {
        assert!(BedTable::load(&raw[..raw.len() - missing]).is_err());
    }
    let mut bad_crc = raw;
    let trailer = bad_crc.len() - 8;
    bad_crc[trailer] ^= 1;
    assert!(BedTable::load(&bad_crc).is_err());
}
