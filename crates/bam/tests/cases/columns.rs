use crate::support::{Hts, Rec, fixture, open};
use omicsweave_core::{
    AlignmentFields, AlignmentRead, AlignmentSize, AlignmentsOut, CigarColumn, Error,
    ModificationColumn, Region, StringColumn,
};

#[test]
fn fields_are_optional_and_sizes_exact() {
    let f = fixture(7, 3_000, 50, 6);
    let mut r = open(&f.bam, &f.bai);
    let region = Region::new("chr1", 0, 3_000_000);
    let all = r.query(region, AlignmentFields::ALL).unwrap();
    let only = r.query(region, AlignmentFields::NONE).unwrap();
    assert_eq!(
        only,
        AlignmentSize {
            records: all.records,
            ..Default::default()
        }
    );
    let mut entries = vec![Default::default(); only.records];
    r.write(AlignmentsOut {
        entries: &mut entries,
        names: None,
        seqs: None,
        quals: None,
        cigars: None,
        aux: None,
        modifications: None,
    })
    .unwrap();
    let names = r.query(region, AlignmentFields::NAME).unwrap();
    let mut offsets = vec![0i64; names.records + 1];
    let mut data = vec![0u8; names.name_bytes];
    let mut entries = vec![Default::default(); names.records];
    // A buffer for a field the query did not ask for is refused.
    let mut so = vec![0i64; names.records + 1];
    let err = r
        .write(AlignmentsOut {
            entries: &mut entries,
            names: Some(StringColumn {
                offsets: &mut offsets,
                data: &mut data,
            }),
            seqs: Some(StringColumn {
                offsets: &mut so,
                data: &mut [],
            }),
            quals: None,
            cigars: None,
            aux: None,
            modifications: None,
        })
        .unwrap_err();
    assert!(matches!(err, Error::Invalid { .. }), "{err}");
    let mut short = vec![Default::default(); names.records.saturating_sub(1)];
    let err = r
        .write(AlignmentsOut {
            entries: &mut short,
            names: Some(StringColumn {
                offsets: &mut offsets,
                data: &mut data,
            }),
            seqs: None,
            quals: None,
            cigars: None,
            aux: None,
            modifications: None,
        })
        .unwrap_err();
    assert!(matches!(err, Error::OutputSize { .. }), "{err}");
}

#[test]
fn selected_columns_match_htslib_with_guarded_outputs() {
    let f = fixture(101, 24, 50, 0);
    let hts = Hts::open(&f.bam, &f.bai);
    let mut reader = open(&f.bam, &f.bai);
    for region in [
        Region::new("chr1", 0, 3_000_000),
        Region::new("empty", 0, 50_000),
    ] {
        let expected = hts.query(
            hts.tid(region.chrom),
            region.start as i64,
            region.end as i64,
        );
        for bits in 0..=AlignmentFields::ALL.bits() {
            let fields = AlignmentFields::from_bits(bits).unwrap();
            let size = reader.query(region, fields).unwrap();
            assert_eq!(size.records, expected.len());
            let n = size.records + 1;
            let mut entries = vec![omicsweave_core::AlignmentEntry::default(); size.records + 2];
            entries[0].start = -777;
            entries[size.records + 1].start = -777;
            let (mut no, mut nd) = (vec![-777; n + 2], vec![0x5a; size.name_bytes + 2]);
            let (mut so, mut sd) = (vec![-777; n + 2], vec![0x5a; size.seq_bytes + 2]);
            let (mut qo, mut qd) = (vec![-777; n + 2], vec![0x5a; size.seq_bytes + 2]);
            let (mut ao, mut ad) = (vec![-777; n + 2], vec![0x5a; size.aux_bytes + 2]);
            let (mut co, mut cd) = (vec![-777; n + 2], vec![0xdeadbeef; size.cigar_ops + 2]);
            // The random records carry no MM tags.
            assert_eq!(size.modifications, 0);
            let mut mo = vec![-777; n + 2];
            reader
                .write(AlignmentsOut {
                    entries: &mut entries[1..size.records + 1],
                    names: fields
                        .contains(AlignmentFields::NAME)
                        .then_some(StringColumn {
                            offsets: &mut no[1..n + 1],
                            data: &mut nd[1..size.name_bytes + 1],
                        }),
                    seqs: fields
                        .contains(AlignmentFields::SEQ)
                        .then_some(StringColumn {
                            offsets: &mut so[1..n + 1],
                            data: &mut sd[1..size.seq_bytes + 1],
                        }),
                    quals: fields
                        .contains(AlignmentFields::QUAL)
                        .then_some(StringColumn {
                            offsets: &mut qo[1..n + 1],
                            data: &mut qd[1..size.seq_bytes + 1],
                        }),
                    aux: fields
                        .contains(AlignmentFields::AUX)
                        .then_some(StringColumn {
                            offsets: &mut ao[1..n + 1],
                            data: &mut ad[1..size.aux_bytes + 1],
                        }),
                    cigars: fields
                        .contains(AlignmentFields::CIGAR)
                        .then_some(CigarColumn {
                            offsets: &mut co[1..n + 1],
                            ops: &mut cd[1..size.cigar_ops + 1],
                        }),
                    modifications: fields.contains(AlignmentFields::MODIFICATIONS).then_some(
                        ModificationColumn {
                            offsets: &mut mo[1..n + 1],
                            entries: &mut [],
                        },
                    ),
                })
                .unwrap();
            assert_eq!((mo[0], mo[n + 1]), (-777, -777));
            if fields.contains(AlignmentFields::MODIFICATIONS) {
                assert!(mo[1..n + 1].iter().all(|&o| o == 0));
            }
            for (flag, offsets, data, part) in [
                (
                    AlignmentFields::NAME,
                    &no,
                    &nd,
                    (|r: &Rec| &r.name) as fn(&Rec) -> &Vec<u8>,
                ),
                (
                    AlignmentFields::SEQ,
                    &so,
                    &sd,
                    (|r: &Rec| &r.seq) as fn(&Rec) -> &Vec<u8>,
                ),
                (
                    AlignmentFields::QUAL,
                    &qo,
                    &qd,
                    (|r: &Rec| &r.qual) as fn(&Rec) -> &Vec<u8>,
                ),
                (
                    AlignmentFields::AUX,
                    &ao,
                    &ad,
                    (|r: &Rec| &r.aux) as fn(&Rec) -> &Vec<u8>,
                ),
            ] {
                assert_eq!((offsets[0], offsets[n + 1]), (-777, -777));
                assert_eq!((data[0], data[data.len() - 1]), (0x5a, 0x5a));
                if fields.contains(flag) {
                    let mut at = 0;
                    assert_eq!(offsets[1], 0);
                    for (i, record) in expected.iter().enumerate() {
                        let bytes = part(record);
                        assert_eq!(&data[1 + at..1 + at + bytes.len()], bytes);
                        at += bytes.len();
                        assert_eq!(offsets[i + 2], at as i64);
                    }
                    assert_eq!(at, data.len() - 2);
                }
            }
            assert_eq!((co[0], co[n + 1]), (-777, -777));
            assert_eq!((cd[0], cd[cd.len() - 1]), (0xdeadbeef, 0xdeadbeef));
            if fields.contains(AlignmentFields::CIGAR) {
                let want: Vec<_> = expected
                    .iter()
                    .flat_map(|r| r.cigar.iter().copied())
                    .collect();
                assert_eq!(&cd[1..cd.len() - 1], want);
                let mut at = 0;
                assert_eq!(co[1], 0);
                for (i, rec) in expected.iter().enumerate() {
                    at += rec.cigar.len();
                    assert_eq!(co[i + 2], at as i64);
                }
            }
            assert_eq!(
                (entries[0].start, entries[size.records + 1].start),
                (-777, -777)
            );
            for (entry, record) in entries[1..size.records + 1].iter().zip(&expected) {
                assert_eq!((entry.start, entry.end), (record.pos, record.end));
            }
        }
    }
}
