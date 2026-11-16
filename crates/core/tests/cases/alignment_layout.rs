//! Exact capacities and field selection for alignment output columns.

use omicsweave_core::{
    AlignmentEntry, AlignmentFields, AlignmentSize, AlignmentsOut, CigarColumn, ModificationColumn,
    ModificationEntry, StringColumn,
};

struct Buffers {
    entries: Vec<AlignmentEntry>,
    offsets: [Vec<i64>; 6],
    names: Vec<u8>,
    seqs: Vec<u8>,
    quals: Vec<u8>,
    cigars: Vec<u32>,
    aux: Vec<u8>,
    modifications: Vec<ModificationEntry>,
}

impl Buffers {
    fn new(size: AlignmentSize) -> Self {
        Self {
            entries: vec![AlignmentEntry::default(); size.records],
            offsets: std::array::from_fn(|_| vec![0; size.records + 1]),
            names: vec![0; size.name_bytes],
            seqs: vec![0; size.seq_bytes],
            quals: vec![0; size.seq_bytes],
            cigars: vec![0; size.cigar_ops],
            aux: vec![0; size.aux_bytes],
            modifications: vec![ModificationEntry::default(); size.modifications],
        }
    }

    fn output(&mut self, fields: AlignmentFields) -> AlignmentsOut<'_> {
        let [
            name_offsets,
            seq_offsets,
            qual_offsets,
            cigar_offsets,
            aux_offsets,
            mod_offsets,
        ] = &mut self.offsets;
        AlignmentsOut {
            entries: &mut self.entries,
            names: fields
                .contains(AlignmentFields::NAME)
                .then_some(StringColumn {
                    offsets: name_offsets,
                    data: &mut self.names,
                }),
            seqs: fields
                .contains(AlignmentFields::SEQ)
                .then_some(StringColumn {
                    offsets: seq_offsets,
                    data: &mut self.seqs,
                }),
            quals: fields
                .contains(AlignmentFields::QUAL)
                .then_some(StringColumn {
                    offsets: qual_offsets,
                    data: &mut self.quals,
                }),
            cigars: fields
                .contains(AlignmentFields::CIGAR)
                .then_some(CigarColumn {
                    offsets: cigar_offsets,
                    ops: &mut self.cigars,
                }),
            aux: fields
                .contains(AlignmentFields::AUX)
                .then_some(StringColumn {
                    offsets: aux_offsets,
                    data: &mut self.aux,
                }),
            modifications: fields.contains(AlignmentFields::MODIFICATIONS).then_some(
                ModificationColumn {
                    offsets: mod_offsets,
                    entries: &mut self.modifications,
                },
            ),
        }
    }
}

#[test]
fn every_field_selection_requires_exactly_its_queried_columns() {
    for size in [
        AlignmentSize::default(),
        AlignmentSize {
            records: 2,
            name_bytes: 5,
            seq_bytes: 7,
            cigar_ops: 3,
            aux_bytes: 11,
            modifications: 4,
        },
    ] {
        for bits in 0..=AlignmentFields::ALL.bits() {
            let fields = AlignmentFields::from_bits(bits).unwrap();
            let mut buffers = Buffers::new(size);
            buffers.output(fields).check(fields, size).unwrap();
            for bit in [
                AlignmentFields::NAME,
                AlignmentFields::SEQ,
                AlignmentFields::QUAL,
                AlignmentFields::CIGAR,
                AlignmentFields::AUX,
                AlignmentFields::MODIFICATIONS,
            ] {
                let mismatched = AlignmentFields::from_bits(bits ^ bit.bits()).unwrap();
                assert!(
                    buffers.output(mismatched).check(fields, size).is_err(),
                    "fields={bits:#x} mismatched column={:#x} size={size:?}",
                    bit.bits(),
                );
            }
        }
    }
}

#[test]
fn each_column_rejects_short_and_oversized_offsets_and_payloads() {
    let size = AlignmentSize {
        records: 2,
        name_bytes: 3,
        seq_bytes: 3,
        cigar_ops: 3,
        aux_bytes: 3,
        modifications: 3,
    };
    for capacity in [2, 4] {
        for column in 0..6 {
            let mut buffers = Buffers::new(size);
            buffers.offsets[column].resize(capacity, 0);
            assert!(
                buffers
                    .output(AlignmentFields::ALL)
                    .check(AlignmentFields::ALL, size)
                    .is_err(),
                "column={column} offset capacity={capacity}",
            );
            let mut buffers = Buffers::new(size);
            match column {
                0 => buffers.names.resize(capacity, 0),
                1 => buffers.seqs.resize(capacity, 0),
                2 => buffers.quals.resize(capacity, 0),
                3 => buffers.cigars.resize(capacity, 0),
                4 => buffers.aux.resize(capacity, 0),
                5 => buffers
                    .modifications
                    .resize(capacity, ModificationEntry::default()),
                _ => unreachable!(),
            }
            assert!(
                buffers
                    .output(AlignmentFields::ALL)
                    .check(AlignmentFields::ALL, size)
                    .is_err(),
                "column={column} payload capacity={capacity}",
            );
        }
    }
}
