//! Explicit bytes oracle for mixed BED syntax and half-open boundary queries.
use std::io::Write;

use omicsweave_bed::BedTable;
use omicsweave_core::{IntervalEntry, IntervalRead, IntervalsOut, Region, StringColumn};

fn check(table: &mut BedTable, region: Region<'_>, expected: &[(u32, u32, &[u8])]) {
    let size = table.query(region).unwrap();
    assert_eq!(size.entries, expected.len());
    assert_eq!(
        size.annotation_bytes,
        expected.iter().map(|r| r.2.len()).sum::<usize>()
    );
    let mut entries = vec![IntervalEntry::default(); size.entries];
    let mut offsets = vec![0; size.entries + 1];
    let mut data = vec![0; size.annotation_bytes];
    table
        .write(IntervalsOut {
            entries: &mut entries,
            annotations: StringColumn {
                offsets: &mut offsets,
                data: &mut data,
            },
        })
        .unwrap();
    let mut byte = 0;
    for (row, &(start, end, note)) in expected.iter().enumerate() {
        assert_eq!(
            entries[row],
            IntervalEntry {
                start,
                end,
                annotation: row as u32
            }
        );
        assert_eq!(offsets[row], byte as i64);
        assert_eq!(&data[byte..byte + note.len()], note);
        byte += note.len();
    }
    assert_eq!(offsets[expected.len()], byte as i64);
}

#[test]
fn mixed_syntax_keeps_full_annotation_bytes_and_stable_duplicate_rows() {
    // A tab selects tab syntax for the entire line; whitespace fields cannot
    // supply missing tab-delimited coordinates.
    assert!(BedTable::load(b"c 16 17 after\tc\t18\t19\tignored").is_err());
    let text = b"track name=x\r\n#comment\nbrowser position c:1-20\nc 10 12 first second\nc\t0\t20\touter\nc\t10\t10\tzero\nd\t1\t2\tother\nc\t10\t12\t\xff\0\t\x80\t\nc\t10\t12\t\xff\0\t\x80\t\nc\t12\t15\r\nc 16 17 after";
    let mut encoded = Vec::new();
    // The member boundary cuts through a tab-delimited binary annotation.
    let split = text.iter().position(|&b| b == 0xff).unwrap() + 2;
    for member in [&text[..split], &text[split..]] {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        encoder.write_all(member).unwrap();
        encoded.extend_from_slice(&encoder.finish().unwrap());
    }
    for raw in [text.as_slice(), encoded.as_slice()] {
        let mut table = BedTable::load(raw).unwrap();
        assert_eq!(table.chroms(), &["c", "d"]);
        assert_eq!(table.len(), 8);
        assert_eq!(table.annotations().get(0), b"first\tsecond");
        assert_eq!(table.annotations().get(4), b"\xff\0\t\x80\t");
        assert_eq!(table.annotations().get(5), b"\xff\0\t\x80\t");
        let duplicate: &[u8] = b"\xff\0\t\x80\t";
        check(
            &mut table,
            Region::new("c", 9, 11),
            &[
                (0, 20, b"outer"),
                (10, 12, b"first\tsecond"),
                (10, 10, b"zero"),
                (10, 12, duplicate),
                (10, 12, duplicate),
            ],
        );
        check(
            &mut table,
            Region::new("c", 10, 11),
            &[
                (0, 20, b"outer"),
                (10, 12, b"first\tsecond"),
                (10, 12, duplicate),
                (10, 12, duplicate),
            ],
        );
        check(
            &mut table,
            Region::new("c", 12, 13),
            &[(0, 20, b"outer"), (12, 15, b"")],
        );
        check(
            &mut table,
            Region::new("c", 16, 17),
            &[(0, 20, b"outer"), (16, 17, b"after")],
        );
        check(&mut table, Region::new("c", 10, 10), &[]);
        check(&mut table, Region::new("missing", 0, 20), &[]);
    }
}
