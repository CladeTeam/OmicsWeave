use omicsweave_core::{ContactQuery, ContactRead, ContactUnit, Error, Region};
use omicsweave_hic::{
    BlockEntry, Chrom, Depth, HicHeader, HicReader, HicWriter, Record, ZoomData, decode_block,
    encode_block,
};
use proptest::prelude::*;
use std::collections::BTreeSet;
use std::io::Cursor;

fn zoom(bb: u32, columns: u32, intra: bool) -> ZoomData {
    ZoomData {
        chr1: 0,
        chr2: u32::from(!intra),
        unit: ContactUnit::Bp,
        bin_size: 1,
        sum_counts: 0.0,
        block_bin_count: bb,
        block_column_count: columns,
        blocks: (0..2048)
            .filter(|i| i % 3 != 0)
            .map(|number| BlockEntry {
                number,
                position: 0,
                size: 1,
            })
            .collect(),
    }
}

proptest! {
    #[test]
    fn block_selection_matches_candidate_enumeration(
        bb in 1u32..16, columns in 1u32..32, x in 0u32..128, y in 0u32..128,
        w in 0u32..64, h in 0u32..64, intra in any::<bool>(), v9 in any::<bool>(), kind in 0u32..3,
    ) {
        let depth = match kind { 0 => Depth::Log(2), 1 => Depth::Log(3), _ => Depth::Linear(4) };
        let version = if v9 {9} else {8};
        let z = zoom(bb,columns,intra);
        let mut numbers = BTreeSet::new();
        if w>0 && h>0 {
            let (x1,x2,y1,y2)=(x as i64,(x+w-1) as i64,y as i64,(y+h-1) as i64);
            let bb=bb as i64;
            let rects = if v9 && intra {
                let (a,b)=(depth.of(x1-y2,bb),depth.of(x2-y1,bb));
                vec![(if x1<=y2 && y1<=x2 {0} else {a.min(b)}, a.max(b)+1,(x1+y1)/2/bb,(x2+y2)/2/bb+1)]
            } else {
                let (cl,ch,rl,rh)=(x1/bb,(x2+1)/bb,y1/bb,(y2+1)/bb);
                let mut rects=vec![(rl,rh,cl,ch)];
                if intra { rects.push((cl,ch,rl,rh)); }
                rects
            };
            for (rl,rh,cl,ch) in rects { for r in rl..=rh { for c in cl..=ch { numbers.insert((r*columns as i64+c) as u32); } } }
        }
        let want:Vec<_>=z.blocks.iter().copied().filter(|b| numbers.contains(&b.number)).collect();
        prop_assert_eq!(z.blocks_for(version,depth,x..x+w,y..y+h), want);
    }

    #[test]
    fn arbitrary_block_bytes_do_not_panic(bytes in prop::collection::vec(any::<u8>(),0..512), version in 6i32..=9) {
        let _ = decode_block(version,&bytes,|_,_,_|{});
    }
}

#[test]
fn huge_queries_and_invalid_grids_are_bounded() {
    let z = zoom(1, 1, false);
    assert_eq!(
        z.blocks_for(8, Depth::Log(2), 0..u32::MAX, 0..u32::MAX),
        z.blocks
    );
    for (bb, cols) in [(0, 1), (1, 0)] {
        assert!(
            zoom(bb, cols, true)
                .blocks_for(9, Depth::Log(2), 0..1, 0..1)
                .is_empty()
        );
    }
}

#[test]
fn codecs_cover_narrow_wide_and_unaligned_records() {
    for version in [8, 9] {
        for wide in [false, true] {
            for short in [false, true] {
                if wide && version < 9 {
                    continue;
                }
                let records: Vec<_> = (0..64u32)
                    .flat_map(|y| {
                        (0..8u32).map(move |x| Record {
                            x: if wide { x * 50_000 } else { x },
                            y: if wide { y * 40_000 } else { y },
                            value: if short { 32767.0 } else { (x + y) as f32 / 7.0 },
                        })
                    })
                    .collect();
                let mut bytes = Vec::new();
                encode_block(version, &records, &mut bytes).unwrap();
                let mut aligned = vec![0xAA];
                aligned.extend(&bytes);
                let mut actual = Vec::new();
                decode_block(version, &aligned[1..], |x, y, value| {
                    actual.push(Record {
                        x: x as u32,
                        y: y as u32,
                        value,
                    })
                })
                .unwrap();
                assert_eq!(actual, records);
                for end in 0..bytes.len() {
                    assert!(decode_block(version, &bytes[..end], |_, _, _| {}).is_err());
                }
            }
        }
    }
    assert!(decode_block(6, &(-1i32).to_le_bytes(), |_, _, _| {}).is_err());
    for records in [
        vec![
            Record {
                x: 0,
                y: 1,
                value: 1.0,
            },
            Record {
                x: 0,
                y: 0,
                value: 2.0,
            },
        ],
        vec![Record {
            x: 0,
            y: 0,
            value: f32::INFINITY,
        }],
    ] {
        assert!(encode_block(9, &records, &mut Vec::new()).is_err());
    }
}

#[test]
fn reader_checks_zero_resolution_and_writer_checks_duplicates() {
    let chrom = Chrom {
        name: "s".into(),
        len: 100,
    };
    let duplicate = HicHeader::new("g", vec![chrom.clone(), chrom.clone()], vec![1]);
    assert!(HicWriter::create(Cursor::new(Vec::new()), duplicate, 1).is_err());
    let mut writer = HicWriter::create(
        Cursor::new(Vec::new()),
        HicHeader::new("g", vec![chrom], vec![1]),
        1,
    )
    .unwrap();
    let overflow = [Record {
        x: 0,
        y: 0,
        value: f32::MAX,
    }; 2];
    assert!(
        writer
            .write_matrix(0, 0, ContactUnit::Bp, 1, &overflow)
            .is_err()
    );
    writer
        .write_matrix(
            0,
            0,
            ContactUnit::Bp,
            1,
            &[Record {
                x: 0,
                y: 0,
                value: 1.0,
            }],
        )
        .unwrap();
    let bytes = writer.finish().unwrap().into_inner();
    let mut reader = HicReader::open(&bytes[..]).unwrap();
    let region = Region::new("s", 0, 10);
    assert!(matches!(
        reader.matrix(ContactQuery::new(region, region, 0), &mut []),
        Err(Error::Invalid { .. })
    ));
    let mut invalid = bytes.clone();
    let pos = invalid.windows(2).position(|w| w == b"s\0").unwrap();
    invalid[pos] = 0xFF;
    assert!(HicReader::open(&invalid[..]).is_err());
    // Repeated chromosome names in a source must not silently overwrite IDs.
    let header = HicHeader::new(
        "g",
        vec![
            Chrom {
                name: "a".into(),
                len: 100,
            },
            Chrom {
                name: "b".into(),
                len: 100,
            },
        ],
        vec![1],
    );
    let writer = HicWriter::create(Cursor::new(Vec::new()), header, 1).unwrap();
    let mut duplicate = writer.finish().unwrap().into_inner();
    let pos = duplicate.windows(2).position(|w| w == b"b\0").unwrap();
    duplicate[pos] = b'a';
    assert!(HicReader::open(&duplicate[..]).is_err());
}

#[test]
fn merged_floats_follow_contact_input_order() {
    let header = HicHeader::new(
        "g",
        vec![Chrom {
            name: "s".into(),
            len: 7,
        }],
        vec![1],
    );
    let mut records = Vec::new();
    let mut want = [[0.0f32; 7]; 7];
    for i in 0..1024usize {
        let (x, y) = (i % 7, (i * 3 / 7) % 7);
        let value = [100_000_000.0, 1.0, -100_000_000.0, 3.0][i % 4];
        records.push(Record {
            x: x as u32,
            y: y as u32,
            value,
        });
        let (a, b) = (x.min(y), x.max(y));
        want[a][b] += value;
    }
    let mut writer = HicWriter::create(Cursor::new(Vec::new()), header, 1).unwrap();
    writer
        .write_matrix(0, 0, ContactUnit::Bp, 1, &records)
        .unwrap();
    let bytes = writer.finish().unwrap().into_inner();
    let mut reader = HicReader::open(&bytes[..]).unwrap();
    let region = Region::new("s", 0, 7);
    let mut out = [0.0; 49];
    reader
        .matrix(ContactQuery::new(region, region, 1), &mut out)
        .unwrap();
    for x in 0..7 {
        for y in 0..7 {
            assert_eq!(out[x * 7 + y].to_bits(), want[x.min(y)][x.max(y)].to_bits());
        }
    }
}

#[test]
fn inflater_grows_reuses_and_rejects_corruption() {
    use omicsweave_hic::{Deflater, Inflater};
    let text = vec![b'A'; 1 << 18];
    let mut encoder = Deflater::new(9);
    let mut raw = Vec::new();
    encoder.deflate(&text, &mut raw).unwrap();
    let mut decoder = Inflater::new();
    let mut out = Vec::new();
    decoder.inflate(&raw, &mut out).unwrap();
    assert_eq!(out, text);
    for end in [0, 1, raw.len() / 2, raw.len() - 1] {
        assert!(decoder.inflate(&raw[..end], &mut out).is_err());
    }
    let mut broken = raw.clone();
    *broken.last_mut().unwrap() ^= 1;
    assert!(decoder.inflate(&broken, &mut out).is_err());
    decoder.inflate(&raw, &mut out).unwrap();
    assert_eq!(out, text);
    encoder.deflate(b"xy", &mut raw).unwrap();
    decoder.inflate(&raw, &mut out).unwrap();
    assert_eq!(out, b"xy");
}

#[test]
fn adjacent_compressed_blocks_share_one_positional_read() {
    use omicsweave_core::{ReadAt, Result};
    use std::{cell::Cell, rc::Rc};
    struct Counted<'a> {
        bytes: &'a [u8],
        reads: Rc<Cell<usize>>,
    }
    impl ReadAt for Counted<'_> {
        fn len(&self) -> u64 {
            self.bytes.len() as u64
        }
        fn read_exact_at(&self, out: &mut [u8], offset: u64) -> Result<()> {
            self.reads.set(self.reads.get() + 1);
            self.bytes.read_exact_at(out, offset)
        }
    }
    let header = HicHeader::new(
        "g",
        vec![Chrom {
            name: "s".into(),
            len: 2_000_000,
        }],
        vec![1000],
    );
    let mut writer = HicWriter::create(Cursor::new(Vec::new()), header, 1).unwrap();
    let records: Vec<_> = (0..2000u32)
        .map(|x| Record {
            x,
            y: x,
            value: 1.0,
        })
        .collect();
    writer
        .write_matrix(0, 0, ContactUnit::Bp, 1000, &records)
        .unwrap();
    let bytes = writer.finish().unwrap().into_inner();
    let reads = Rc::new(Cell::new(0));
    let mut reader = HicReader::open(Counted {
        bytes: &bytes,
        reads: reads.clone(),
    })
    .unwrap();
    let blocks = reader
        .zoom(0, 0, ContactUnit::Bp, 1000)
        .unwrap()
        .unwrap()
        .blocks
        .len();
    assert!(blocks > 1);
    reads.set(0);
    let region = Region::new("s", 0, 2_000_000);
    let n = reader
        .query_entries(ContactQuery::new(region, region, 1000))
        .unwrap();
    assert_eq!(n, 2000);
    assert_eq!(reads.get(), 1);
}

#[test]
fn blocks_of_other_versions_are_unsupported() {
    let record = [Record {
        x: 0,
        y: 0,
        value: 1.0,
    }];
    for version in [6, 7, 10] {
        let mut out = Vec::new();
        assert!(matches!(
            encode_block(version, &record, &mut out),
            Err(omicsweave_core::Error::Unsupported(_))
        ));
        assert!(matches!(
            decode_block(version, &[0; 16], |_, _, _| {}),
            Err(omicsweave_core::Error::Unsupported(_))
        ));
    }
}
