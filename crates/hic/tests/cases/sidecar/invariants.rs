use omicsweave_core::ContactUnit;
use omicsweave_core::sidecar::SourceIdentity;
use omicsweave_hic::{BlockSize, HicSidecar, MatrixBlocks};

#[test]
fn sidecar_checks_index_invariants_even_with_a_valid_container_checksum() {
    let identity = SourceIdentity {
        len: 100,
        fingerprint: 7,
    };
    let valid = HicSidecar {
        version: 9,
        matrices: vec![MatrixBlocks {
            chr1: 0,
            chr2: 1,
            unit: ContactUnit::Bp,
            bin_size: 1,
            block_bin_count: 16,
            block_column_count: 8,
            blocks: [1, 2]
                .into_iter()
                .map(|number| BlockSize {
                    number,
                    position: 80,
                    size: 20,
                    inflated: 64,
                    records: 2,
                })
                .collect(),
        }],
    };
    assert_eq!(
        HicSidecar::parse(&valid.to_bytes(identity), identity).unwrap(),
        valid
    );
    for case in 0..8 {
        let mut invalid = valid.clone();
        let matrix = &mut invalid.matrices[0];
        match case {
            0 => matrix.blocks.reverse(),
            1 => matrix.blocks[1].number = matrix.blocks[0].number,
            2 => matrix.blocks[0].position = 81,
            3 => matrix.blocks[0].position = u64::MAX - 10,
            4 => matrix.chr1 = 2,
            5 => matrix.bin_size = 0,
            6 => matrix.block_bin_count = 0,
            7 => matrix.block_column_count = 0,
            _ => unreachable!(),
        }
        assert!(
            HicSidecar::parse(&invalid.to_bytes(identity), identity).is_err(),
            "case {case}"
        );
    }
}
