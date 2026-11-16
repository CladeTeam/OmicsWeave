use std::collections::BTreeMap;
use std::ffi::{CStr, CString, c_char, c_int, c_void};
use std::io::Read;

use htslib_sys as h;
use omicsweave_bgzf::{BinIndex, BinningIndex, Blocks, IndexFormat, Logical, Reader, Writer};
use omicsweave_core::exec::Inline;

unsafe extern "C" {
    fn hts_idx_init(n: c_int, fmt: c_int, offset: u64, shift: c_int, levels: c_int) -> *mut c_void;
    fn hts_idx_push(
        idx: *mut c_void,
        tid: c_int,
        beg: i64,
        end: i64,
        offset: u64,
        mapped: c_int,
    ) -> c_int;
    fn hts_idx_finish(idx: *mut c_void, end: u64) -> c_int;
    fn hts_idx_save_as(
        idx: *const c_void,
        name: *const c_char,
        index: *const c_char,
        fmt: c_int,
    ) -> c_int;
    fn hts_idx_destroy(idx: *mut c_void);
}

struct Index(*mut c_void);
impl Drop for Index {
    fn drop(&mut self) {
        // SAFETY: the index is uniquely owned and destroyed exactly once.
        unsafe { hts_idx_destroy(self.0) }
    }
}

type Record = (i32, i64, i64, bool);
type Bin = (u64, Vec<(u64, u64)>);
type Canonical = (i32, i32, Vec<u8>, Vec<BTreeMap<u32, Bin>>, u64);

fn canonical(raw: &[u8]) -> Canonical {
    let mut text = Vec::new();
    Reader::new(raw, Box::new(Inline::new()))
        .read_to_end(&mut text)
        .unwrap();
    let mut at = 0;
    let mut take = |n: usize| {
        let b = &text[at..at + n];
        at += n;
        b
    };
    assert_eq!(take(4), b"CSI\x01");
    let shift = i32::from_le_bytes(take(4).try_into().unwrap());
    let levels = i32::from_le_bytes(take(4).try_into().unwrap());
    let aux = u32::from_le_bytes(take(4).try_into().unwrap()) as usize;
    let aux = take(aux).to_vec();
    let count = u32::from_le_bytes(take(4).try_into().unwrap());
    let mut refs = Vec::new();
    for _ in 0..count {
        let bins = u32::from_le_bytes(take(4).try_into().unwrap());
        let mut reference = BTreeMap::new();
        for _ in 0..bins {
            let id = u32::from_le_bytes(take(4).try_into().unwrap());
            let loff = u64::from_le_bytes(take(8).try_into().unwrap());
            let count = u32::from_le_bytes(take(4).try_into().unwrap());
            let chunks = (0..count)
                .map(|_| {
                    let a = u64::from_le_bytes(take(8).try_into().unwrap());
                    let b = u64::from_le_bytes(take(8).try_into().unwrap());
                    (a, b)
                })
                .collect();
            assert!(reference.insert(id, (loff, chunks)).is_none());
        }
        refs.push(reference);
    }
    let unmapped = u64::from_le_bytes(take(8).try_into().unwrap());
    assert_eq!(at, text.len());
    (shift, levels, aux, refs, unmapped)
}

fn compare(records: &[Record], shift: u32, levels: u32, refs: usize) {
    let blocks = Blocks {
        starts: (0..=records.len())
            .map(|i| (i as u64 * 100_000, i as u64))
            .collect(),
        end: (records.len() as u64 + 1) * 100_000,
        text_len: records.len() as u64,
    };
    let mut ours = BinningIndex::try_new(
        IndexFormat::Csi {
            min_shift: shift,
            levels,
        },
        refs,
        Logical::new(0, 0),
    )
    .unwrap();
    // SAFETY: parameters are within htslib's index geometry and the new
    // opaque handle is checked before use, then uniquely owned by Index.
    let other =
        Index(unsafe { hts_idx_init(refs as c_int, 0, 0, shift as c_int, levels as c_int) });
    assert!(!other.0.is_null());
    for (i, &(tid, beg, end, mapped)) in records.iter().enumerate() {
        let at = Logical::new((i + 1) as u64, 0);
        ours.push(tid, beg, end, at, mapped).unwrap();
        // SAFETY: the index is live; values and physical offsets describe
        // the same sorted records passed to the Rust builder.
        assert_eq!(
            unsafe { hts_idx_push(other.0, tid, beg, end, blocks.resolve(at), mapped as c_int) },
            0
        );
    }
    let end = Logical::new(records.len() as u64, 0);
    let raw = ours.finish(end, &blocks).unwrap();
    // SAFETY: the index is live and final offset follows its last record.
    assert_eq!(unsafe { hts_idx_finish(other.0, blocks.resolve(end)) }, 0);
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("oracle.csi");
    let cpath = CString::new(path.as_os_str().as_encoded_bytes()).unwrap();
    // SAFETY: C strings outlive this call and the live index is not mutated.
    assert_eq!(
        unsafe { hts_idx_save_as(other.0, c"unused".as_ptr(), cpath.as_ptr(), 0) },
        0
    );
    let expected = std::fs::read(path).unwrap();
    assert_eq!(canonical(&raw), canonical(&expected));
    let ours = BinIndex::parse(&raw).unwrap();
    let other = BinIndex::parse(&expected).unwrap();
    for tid in 0..refs {
        for &(s, e) in &[
            (0, 1),
            (0, 1_000_000),
            (1, 90_000),
            (65_535, 100_001),
            (300_000, 500_000),
        ] {
            assert_eq!(ours.chunks(tid, s, e), other.chunks(tid, s, e));
        }
    }
}

#[test]
fn sparse_csi_loff_chunks_and_counts_match_htslib() {
    for seed in 0..8u64 {
        let mut rng = seed + 1;
        let mut records = Vec::new();
        for tid in 0..3 {
            let mut pos = 0;
            for i in 0usize..120 {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                pos += ((rng >> 32) % 4000) as i64;
                let width = if i % 23 == 0 {
                    150_000
                } else {
                    (rng % 300) + 1
                } as i64;
                records.push((tid, pos, pos + width, !i.is_multiple_of(13)));
            }
        }
        records.push((-1, -1, 0, false));
        compare(&records, 0, 7, 3);
        compare(&records, 14, 5, 3);
    }
    compare(
        &[(1, 70_000, 70_010, true), (1, 200_000, 200_001, false)],
        0,
        7,
        3,
    );
}

#[test]
fn csi_empty_intervals_at_and_between_window_boundaries_match_htslib() {
    for mapped in [false, true] {
        for records in [
            vec![(0, 70_000, 70_000, mapped)],
            vec![(0, 65_536, 65_536, mapped), (0, 70_000, 70_000, mapped)],
            vec![(0, 70_000, 70_000, mapped), (0, 100_000, 100_000, mapped)],
            vec![(0, 10, 20, true), (0, 70_000, 70_000, mapped)],
            vec![(0, 0, 100_000, true), (0, 70_000, 70_000, mapped)],
            vec![(0, 0, 80_000, true), (0, 100_000, 100_000, mapped)],
            vec![
                (0, 0, 80_000, true),
                (0, 100_000, 100_000, mapped),
                (0, 120_000, 120_002, true),
            ],
            vec![(0, 140_000, 140_000, mapped)],
            vec![(0, 140_001, 140_001, mapped)],
            vec![(0, 1 << 21, 1 << 21, mapped)],
        ] {
            compare(&records, 0, 7, 1);
            compare(&records, 1, 7, 1);
        }
    }
}

#[test]
fn sparse_empty_interval_serialization_keeps_boundary_windows_unset() {
    let blocks = Blocks {
        starts: vec![(0, 0)],
        end: 100,
        text_len: 2,
    };
    for mapped in [false, true] {
        for (shift, pos, boundary) in [(0, 70_000, true), (1, 140_000, true), (1, 140_001, false)] {
            let mut index = BinningIndex::try_new(
                IndexFormat::Csi {
                    min_shift: shift,
                    levels: 7,
                },
                1,
                Logical(1),
            )
            .unwrap();
            index.push(0, pos, pos, Logical(2), mapped).unwrap();
            let raw = index.finish(Logical(2), &blocks).unwrap();
            let (_, _, _, refs, _) = canonical(&raw);
            let meta = ((1u64 << 24) - 1) / 7 + 1;
            let mut data_bins = 0;
            for (&bin, (loff, _)) in &refs[0] {
                if bin as u64 != meta {
                    assert_eq!(*loff, if boundary { u64::MAX } else { 1 });
                    data_bins += 1;
                }
            }
            assert!(data_bins > 0);
            let parsed = BinIndex::parse(&raw).unwrap();
            assert_eq!(parsed.counts(0), Some((mapped as u64, (!mapped) as u64)));
        }
    }
}

#[test]
fn billion_window_csi_finishes_and_queries_without_dense_materialization() {
    let blocks = Blocks {
        starts: vec![(0, 0)],
        end: 100,
        text_len: 3,
    };
    let mut index = BinningIndex::try_new(
        IndexFormat::Csi {
            min_shift: 0,
            levels: 10,
        },
        1,
        Logical(0),
    )
    .unwrap();
    index.push(0, 0, 1_000_000_000, Logical(1), true).unwrap();
    index
        .push(0, 1_000_000_001, 1_000_000_002, Logical(2), true)
        .unwrap();
    let raw = index.finish(Logical(2), &blocks).unwrap();
    assert!(raw.len() < 1000);
    let parsed = BinIndex::parse(&raw).unwrap();
    assert_eq!(parsed.counts(0), Some((2, 0)));
    assert!(!parsed.chunks(0, 500_000_000, 500_000_001).is_empty());
    assert!(!parsed.chunks(0, 1_000_000_001, 1_000_000_002).is_empty());
}

#[test]
fn sparse_csi_and_dense_bai_drive_htslib_record_iteration() {
    for format in [
        IndexFormat::Bai,
        IndexFormat::Csi {
            min_shift: 14,
            levels: 5,
        },
        IndexFormat::Csi {
            min_shift: 0,
            levels: 7,
        },
    ] {
        let dir = tempfile::tempdir().unwrap();
        let bam = dir.path().join("records.bam");
        let index = dir.path().join("records.idx");
        let mut writer = Writer::new(Vec::new(), 0, Box::new(Inline::new())).unwrap();
        let text = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1500000\n";
        let mut header = b"BAM\x01".to_vec();
        header.extend_from_slice(&(text.len() as i32).to_le_bytes());
        header.extend_from_slice(text);
        header.extend_from_slice(&1i32.to_le_bytes());
        header.extend_from_slice(&5i32.to_le_bytes());
        header.extend_from_slice(b"chr1\0");
        header.extend_from_slice(&1_500_000i32.to_le_bytes());
        writer.write(&header).unwrap();
        let mut builder = BinningIndex::try_new(format, 1, writer.tell()).unwrap();
        let mut records = Vec::new();
        for i in 0..3000i32 {
            let pos = i * 97;
            let long = i % 73 == 0;
            let mapped = i % 29 != 0 || long;
            let name = format!("r{i}\0");
            let mut body = vec![0u8; 32];
            body[4..8].copy_from_slice(&pos.to_le_bytes());
            body[8..12].copy_from_slice(&((60u32 << 8) | name.len() as u32).to_le_bytes());
            body[12..14].copy_from_slice(&(if long { 3u16 } else { 1 }).to_le_bytes());
            body[14..16].copy_from_slice(&(if mapped { 0u16 } else { 4 }).to_le_bytes());
            body[16..20].copy_from_slice(&10i32.to_le_bytes());
            body[20..24].copy_from_slice(&(-1i32).to_le_bytes());
            body[24..28].copy_from_slice(&(-1i32).to_le_bytes());
            body.extend_from_slice(name.as_bytes());
            let cigar = if long {
                vec![(5u32, 0), (149_990, 3), (5, 0)]
            } else {
                vec![(10, 0)]
            };
            for (len, op) in cigar {
                body.extend_from_slice(&((len << 4) | op).to_le_bytes());
            }
            body.extend_from_slice(&[0x11; 5]);
            body.extend_from_slice(&[30; 10]);
            let mut raw = (body.len() as u32).to_le_bytes().to_vec();
            raw.extend_from_slice(&body);
            writer.flush_try(raw.len()).unwrap();
            builder.amend_last(writer.tell());
            writer.write(&raw).unwrap();
            let end = pos
                + if !mapped {
                    1
                } else if long {
                    150_000
                } else {
                    10
                };
            builder
                .push(0, pos as i64, end as i64, writer.tell(), mapped)
                .unwrap();
            records.push((pos as i64, end as i64, name[..name.len() - 1].to_string()));
        }
        let end = writer.tell();
        let finished = writer.finish().unwrap();
        std::fs::write(&bam, finished.out).unwrap();
        std::fs::write(&index, builder.finish(end, &finished.blocks).unwrap()).unwrap();
        let bam_c = CString::new(bam.as_os_str().as_encoded_bytes()).unwrap();
        let index_c = CString::new(index.as_os_str().as_encoded_bytes()).unwrap();
        // SAFETY: all C handles are checked live and kept until every
        // iterator finishes; paths remain valid for all corresponding calls.
        unsafe {
            let file = h::hts_open(bam_c.as_ptr(), c"r".as_ptr());
            assert!(!file.is_null());
            let header = h::sam_hdr_read(file);
            assert!(!header.is_null());
            let index = h::ow_sam_index_load(file, bam_c.as_ptr(), index_c.as_ptr());
            assert!(!index.is_null());
            let record = h::bam_init1();
            assert!(!record.is_null());
            for start in (0..500_000i64).step_by(9973) {
                for width in [1, 100, 60_000] {
                    let end = start + width;
                    let iterator = h::ow_sam_itr_queryi(index, 0, start, end);
                    assert!(!iterator.is_null());
                    let mut actual = Vec::new();
                    loop {
                        let code = h::ow_sam_itr_next(file, iterator, record);
                        if code < 0 {
                            assert_eq!(code, -1);
                            break;
                        }
                        actual.push(
                            CStr::from_ptr(h::ow_bam_qname(record))
                                .to_str()
                                .unwrap()
                                .to_string(),
                        );
                    }
                    let expected: Vec<_> = records
                        .iter()
                        .filter(|(s, e, _)| *s < end && start < *e)
                        .map(|(_, _, n)| n.clone())
                        .collect();
                    assert_eq!(actual, expected, "{format:?} {start}..{end}");
                    h::ow_hts_itr_destroy(iterator);
                }
            }
            h::bam_destroy1(record);
            h::ow_hts_idx_destroy(index);
            h::sam_hdr_destroy(header);
            assert_eq!(h::hts_close(file), 0);
        }
    }
}
