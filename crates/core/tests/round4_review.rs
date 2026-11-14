use std::panic::{AssertUnwindSafe, catch_unwind};

use omicsweave_core::{IntervalEntry, IntervalWriter, IntervalsOut, StringColumn};

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
