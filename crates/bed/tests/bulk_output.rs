//! Independent expected byte concatenation for partial BED query outputs.
use omicsweave_bed::BedTable;
use omicsweave_core::{
    IntervalEntry, IntervalRead, IntervalSize, IntervalsOut, Region, StringColumn,
};
use std::io::Write;

#[test]
fn partial_spans_preserve_large_empty_and_binary_annotations() {
    let notes: Vec<Vec<u8>> = (0..512usize)
        .map(|row| {
            let len = [0, 1, 15, 256, 2048][row % 5];
            (0..len)
                .map(|at| [0, b'\t', 0xff, b'x', b' '][(row + at) % 5])
                .collect()
        })
        .collect();
    let mut text = b"prefix\t0\t1\toutside-before\n".to_vec();
    for (row, note) in notes.iter().enumerate() {
        write!(text, "chr1\t{}\t{}\t", row * 10, row * 10 + 5).unwrap();
        text.extend_from_slice(note);
        text.push(b'\n');
    }
    text.extend_from_slice(b"tail\t0\t1\toutside-after\n");
    let mut table = BedTable::load(&text).unwrap();
    for region in [
        Region::new("chr1", 0, 5120),
        Region::new("chr1", 10, 100),
        Region::new("chr1", 1250, 3815),
        Region::new("chr1", 2560, 2560),
        Region::new("chr1", 17, 18),
        Region::new("chr1", 0, 5120),
    ] {
        let mut expected_entries = Vec::new();
        let mut expected_offsets = vec![0i64];
        let mut expected_data = Vec::new();
        for (row, note) in notes.iter().enumerate() {
            let start = row as u32 * 10;
            let end = start + 5;
            if !region.is_empty() && start < region.end && end > region.start {
                expected_entries.push(IntervalEntry {
                    start,
                    end,
                    annotation: expected_entries.len() as u32,
                });
                expected_data.extend_from_slice(note);
                expected_offsets.push(expected_data.len() as i64);
            }
        }
        let size = table.query(region).unwrap();
        assert_eq!(
            size,
            IntervalSize {
                entries: expected_entries.len(),
                annotation_bytes: expected_data.len(),
            }
        );
        let guard = IntervalEntry {
            start: 333,
            end: 444,
            annotation: 555,
        };
        let mut entries = vec![guard; size.entries + 2];
        let mut offsets = vec![-77; size.entries + 3];
        let mut data = vec![0xa5; size.annotation_bytes + 2];
        for _ in 0..2 {
            table
                .write(IntervalsOut {
                    entries: &mut entries[1..size.entries + 1],
                    annotations: StringColumn {
                        offsets: &mut offsets[1..size.entries + 2],
                        data: &mut data[1..size.annotation_bytes + 1],
                    },
                })
                .unwrap();
            assert_eq!(&entries[1..size.entries + 1], expected_entries);
            assert_eq!(&offsets[1..size.entries + 2], expected_offsets);
            assert_eq!(&data[1..size.annotation_bytes + 1], expected_data);
            assert_eq!((entries[0], entries[size.entries + 1]), (guard, guard));
            assert_eq!((offsets[0], offsets[size.entries + 2]), (-77, -77));
            assert_eq!((data[0], data[size.annotation_bytes + 1]), (0xa5, 0xa5));
        }
    }
}

#[test]
fn interior_reordering_keeps_annotations_paired_with_stable_entries() {
    let mut table = BedTable::load(
        b"c\t0\t5\ta\nc\t30\t35\tbbb\nc\t10\t15\tcc\nc\t20\t25\tdddd\nc\t40\t45\te\n",
    )
    .unwrap();
    // Sorted endpoints have the same file row ids, while the middle rows move.
    let size = table.query(Region::new("c", 0, 100)).unwrap();
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
    assert_eq!(
        entries.iter().map(|e| e.start).collect::<Vec<_>>(),
        [0, 10, 20, 30, 40]
    );
    assert_eq!(offsets, [0, 1, 3, 7, 10, 11]);
    assert_eq!(data, b"accddddbbbe");
}
