use omicsweave_core::{
    Error, IntervalEntry, IntervalSize, IntervalWriter, IntervalsOut, StringColumn, Strings,
};

#[test]
fn interval_writer_roundtrips_bytes_with_guarded_outputs() {
    let notes: [&[u8]; 3] = [b"a\tbc", b"", b"\xff\x00"];
    let guard = IntervalEntry {
        start: u32::MAX,
        end: u32::MAX,
        annotation: u32::MAX,
    };
    let mut entries = [guard; 5];
    let mut offsets = [-99; 6];
    let mut data = [0x55; 8];
    let mut out = IntervalsOut {
        entries: &mut entries[1..4],
        annotations: StringColumn {
            offsets: &mut offsets[1..5],
            data: &mut data[1..7],
        },
    };
    out.check(IntervalSize {
        entries: notes.len(),
        annotation_bytes: notes.iter().map(|note| note.len()).sum(),
    })
    .unwrap();
    let mut writer = IntervalWriter::new(&mut out);
    for (row, note) in notes.iter().enumerate() {
        writer.push(row as u32, row as u32 + 2, note);
        assert_eq!(writer.rows(), row + 1);
    }
    assert_eq!(entries[0], guard);
    assert_eq!(entries[4], guard);
    assert_eq!(offsets[0], -99);
    assert_eq!(offsets[5], -99);
    assert_eq!(data[0], 0x55);
    assert_eq!(data[7], 0x55);
    let strings = Strings {
        offsets: &offsets[1..5],
        data: &data[1..7],
    };
    assert_eq!(strings.len(), notes.len());
    for (row, note) in notes.iter().enumerate() {
        assert_eq!(entries[row + 1].annotation as usize, row);
        assert_eq!(strings.get(row), *note);
    }
}

#[test]
fn malformed_string_rows_panic_without_accepting_truncated_offsets() {
    for offsets in [
        vec![-1, 0],
        vec![0, -1],
        vec![2, 1],
        vec![0, 4],
        vec![1_i64 << 32, (1_i64 << 32) + 1],
        vec![],
        vec![0],
    ] {
        let strings = Strings {
            offsets: &offsets,
            data: b"abc",
        };
        assert!(std::panic::catch_unwind(|| strings.get(0)).is_err());
    }
    let strings = Strings {
        offsets: &[0, 0],
        data: b"",
    };
    assert!(std::panic::catch_unwind(|| strings.get(usize::MAX)).is_err());
}

#[test]
#[cfg(target_pointer_width = "64")]
fn interval_sizes_reject_unrepresentable_output_fields() {
    let mut entries = [];
    let mut offsets = [0];
    let mut data = [];
    let out = IntervalsOut {
        entries: &mut entries,
        annotations: StringColumn {
            offsets: &mut offsets,
            data: &mut data,
        },
    };
    for size in [
        IntervalSize {
            entries: u32::MAX as usize + 2,
            annotation_bytes: 0,
        },
        IntervalSize {
            entries: 0,
            annotation_bytes: i64::MAX as usize + 1,
        },
    ] {
        assert!(matches!(out.check(size), Err(Error::Invalid { .. })));
    }
    // The last row index, rather than the row count, is stored as u32.
    assert!(matches!(
        out.check(IntervalSize {
            entries: u32::MAX as usize + 1,
            annotation_bytes: 0,
        }),
        Err(Error::OutputSize { .. })
    ));
}
