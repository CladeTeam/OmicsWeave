//! Caller-owned interval rows and variable-length output columns.

use std::panic::{AssertUnwindSafe, catch_unwind};

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

const ENTRY_GUARD: IntervalEntry = IntervalEntry {
    start: 0xfedc_ba98,
    end: 0x7654_3210,
    annotation: 0xabcd_ef01,
};

#[test]
fn safe_writer_handles_arbitrary_output_capacities_without_crossing_guards() {
    for entry_count in 0..=1 {
        for offset_count in 0..=2 {
            for byte_count in 0..=2 {
                for annotation in [b"".as_slice(), b"\x00\xff".as_slice()] {
                    let mut entries = [ENTRY_GUARD; 3];
                    let mut offsets = [-71; 4];
                    let mut bytes = [0x93; 4];
                    let result = catch_unwind(AssertUnwindSafe(|| {
                        let mut out = IntervalsOut {
                            entries: &mut entries[1..1 + entry_count],
                            annotations: StringColumn {
                                offsets: &mut offsets[1..1 + offset_count],
                                data: &mut bytes[1..1 + byte_count],
                            },
                        };
                        // Public safe entry points must remain memory safe even
                        // when callers have not validated these capacities.
                        IntervalWriter::new(&mut out).push(u32::MAX, 0, annotation);
                    }));
                    let fits =
                        entry_count >= 1 && offset_count >= 2 && byte_count >= annotation.len();
                    assert_eq!(result.is_ok(), fits);
                    assert_eq!(entries[0], ENTRY_GUARD);
                    assert!(entries[1 + entry_count..].iter().all(|e| *e == ENTRY_GUARD));
                    assert_eq!(offsets[0], -71);
                    assert!(offsets[1 + offset_count..].iter().all(|&o| o == -71));
                    assert_eq!(bytes[0], 0x93);
                    assert!(bytes[1 + byte_count..].iter().all(|&b| b == 0x93));
                    if fits {
                        assert_eq!(entries[1].start, u32::MAX);
                        assert_eq!(entries[1].end, 0);
                        assert_eq!(entries[1].annotation, 0);
                        assert_eq!(&offsets[1..3], &[0, annotation.len() as i64]);
                        assert_eq!(&bytes[1..1 + annotation.len()], annotation);
                    }
                }
            }
        }
    }
}

#[test]
fn safe_writer_rejects_an_extra_row_with_empty_annotation() {
    let mut entries = [ENTRY_GUARD; 3];
    let mut offsets = [-71; 4];
    let mut bytes = [0x93; 3];
    let result = catch_unwind(AssertUnwindSafe(|| {
        let mut out = IntervalsOut {
            entries: &mut entries[1..2],
            annotations: StringColumn {
                offsets: &mut offsets[1..3],
                data: &mut bytes[1..2],
            },
        };
        let mut writer = IntervalWriter::new(&mut out);
        writer.push(3, 5, b"x");
        assert_eq!(writer.rows(), 1);
        writer.push(9, 9, b"");
    }));
    assert!(result.is_err());
    assert_eq!(entries[0], ENTRY_GUARD);
    assert_eq!(entries[2], ENTRY_GUARD);
    assert_eq!(offsets[0], -71);
    assert_eq!(offsets[3], -71);
    assert_eq!(bytes[0], 0x93);
    assert_eq!(bytes[2], 0x93);
}

#[test]
fn interval_buffers_require_exact_sizes() {
    let mut entries = [];
    let mut offsets = [99];
    let mut data = [];
    let out = IntervalsOut {
        entries: &mut entries,
        annotations: StringColumn {
            offsets: &mut offsets,
            data: &mut data,
        },
    };
    assert!(out.check(IntervalSize::default()).is_ok());
    assert!(
        out.check(IntervalSize {
            entries: usize::MAX,
            annotation_bytes: 0
        })
        .is_err()
    );
    assert!(
        out.check(IntervalSize {
            entries: 0,
            annotation_bytes: 1
        })
        .is_err()
    );
}
