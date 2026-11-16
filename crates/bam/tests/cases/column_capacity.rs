use crate::support::{index, modifications::write_bam, open};
use omicsweave_core::{
    AlignmentEntry, AlignmentFields, AlignmentRead, AlignmentSize, AlignmentsOut, CigarColumn,
    Error, ModificationColumn, ModificationEntry, Region, StringColumn,
};

/// One selected column: all destinations carry sentinels so validation must finish before copying.
struct Buffers {
    entries: Vec<AlignmentEntry>,
    offsets: Vec<i64>,
    bytes: Vec<u8>,
    ops: Vec<u32>,
    mods: Vec<ModificationEntry>,
}

impl Buffers {
    fn new(lengths: [usize; 3]) -> Self {
        Self {
            entries: vec![
                AlignmentEntry {
                    start: -777,
                    ..Default::default()
                };
                lengths[0]
            ],
            offsets: vec![-777; lengths[1]],
            bytes: vec![0x5a; lengths[2]],
            ops: vec![0xdeadbeef; lengths[2]],
            mods: vec![
                ModificationEntry {
                    qpos: -777,
                    ..Default::default()
                };
                lengths[2]
            ],
        }
    }

    fn out(&mut self, field: AlignmentFields) -> AlignmentsOut<'_> {
        let mut out = AlignmentsOut {
            entries: &mut self.entries,
            names: None,
            seqs: None,
            quals: None,
            cigars: None,
            aux: None,
            modifications: None,
        };
        if field == AlignmentFields::CIGAR {
            out.cigars = Some(CigarColumn {
                offsets: &mut self.offsets,
                ops: &mut self.ops,
            });
        } else if field == AlignmentFields::MODIFICATIONS {
            out.modifications = Some(ModificationColumn {
                offsets: &mut self.offsets,
                entries: &mut self.mods,
            });
        } else {
            let column = StringColumn {
                offsets: &mut self.offsets,
                data: &mut self.bytes,
            };
            if field == AlignmentFields::NAME {
                out.names = Some(column);
            } else if field == AlignmentFields::SEQ {
                out.seqs = Some(column);
            } else if field == AlignmentFields::QUAL {
                out.quals = Some(column);
            } else if field == AlignmentFields::AUX {
                out.aux = Some(column);
            } else {
                unreachable!();
            }
        }
        out
    }

    fn assert_untouched(&self) {
        assert!(self.entries.iter().all(|e| e.start == -777));
        assert!(self.offsets.iter().all(|&v| v == -777));
        assert!(self.bytes.iter().all(|&v| v == 0x5a));
        assert!(self.ops.iter().all(|&v| v == 0xdeadbeef));
        assert!(self.mods.iter().all(|m| m.qpos == -777));
    }
}

fn column_len(field: AlignmentFields, size: AlignmentSize) -> usize {
    if field == AlignmentFields::NAME {
        size.name_bytes
    } else if field == AlignmentFields::SEQ || field == AlignmentFields::QUAL {
        size.seq_bytes
    } else if field == AlignmentFields::CIGAR {
        size.cigar_ops
    } else if field == AlignmentFields::AUX {
        size.aux_bytes
    } else {
        size.modifications
    }
}

#[test]
fn each_column_requires_exact_capacity_and_failed_writes_keep_staged_results() {
    let dir = tempfile::tempdir().unwrap();
    let bam = write_bam(
        dir.path(),
        "capacities",
        &["read\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACCG\tIIII\tMM:Z:C+m,0,0;\tML:B:C,100,200".into()],
    );
    let bai = dir.path().join("capacities.bai");
    index(&bam, &bai, 0);
    let mut reader = open(&bam, &bai);
    for field in [
        AlignmentFields::NAME,
        AlignmentFields::SEQ,
        AlignmentFields::QUAL,
        AlignmentFields::CIGAR,
        AlignmentFields::AUX,
        AlignmentFields::MODIFICATIONS,
    ] {
        let size = reader.query(Region::new("chr1", 0, 100), field).unwrap();
        let exact = [size.records, size.records + 1, column_len(field, size)];
        assert!(exact.iter().all(|&n| n > 0));
        for dimension in 0..3 {
            for delta in [-1isize, 1] {
                let mut lengths = exact;
                lengths[dimension] = lengths[dimension].checked_add_signed(delta).unwrap();
                let mut buffers = Buffers::new(lengths);
                let err = reader.write(buffers.out(field)).unwrap_err();
                assert!(
                    matches!(err, Error::OutputSize { have, need, .. } if have == lengths[dimension] && need == exact[dimension]),
                    "field={field:?}, dimension={dimension}, delta={delta}: {err}"
                );
                buffers.assert_untouched();
                let mut retry = Buffers::new(exact);
                reader.write(retry.out(field)).unwrap();
                assert_eq!(
                    retry.entries[0].start, 9,
                    "field={field:?}, dimension={dimension}, delta={delta}"
                );
                assert_eq!(retry.offsets, [0, exact[2] as i64]);
            }
        }
        let mut omitted = vec![
            AlignmentEntry {
                start: -777,
                ..Default::default()
            };
            size.records
        ];
        assert!(
            matches!(
                reader.write(AlignmentsOut {
                    entries: &mut omitted,
                    names: None,
                    seqs: None,
                    quals: None,
                    cigars: None,
                    aux: None,
                    modifications: None,
                }),
                Err(Error::Invalid { .. })
            ),
            "missing {field:?}"
        );
        assert_eq!(omitted[0].start, -777);
        reader
            .query(Region::new("chr1", 0, 100), AlignmentFields::NONE)
            .unwrap();
        let mut extra = Buffers::new(exact);
        assert!(
            matches!(reader.write(extra.out(field)), Err(Error::Invalid { .. })),
            "unexpected {field:?}"
        );
        extra.assert_untouched();
    }
}
