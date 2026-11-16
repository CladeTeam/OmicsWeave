use omicsweave_bbi::{BBI_SIDECAR_VERSION, BbiSidecar, BlockInfo, Kind, LevelBlocks};
use omicsweave_core::sidecar::{self, SidecarKind, SourceIdentity};
#[test]
fn sidecar_v1_encodes_the_maximum_level_count_and_rejects_overflow() {
    let mut sidecar = BbiSidecar {
        kind: Kind::BigWig,
        uncompress_buf_size: 0,
        levels: vec![
            LevelBlocks {
                reduction: 0,
                blocks: Vec::new()
            };
            u16::MAX as usize
        ],
    };
    let body = sidecar.body();
    assert_eq!(u16::from_le_bytes(body[2..4].try_into().unwrap()), u16::MAX);
    assert_eq!(body.len(), 8 + u16::MAX as usize * 16);
    let source = SourceIdentity {
        len: 0,
        fingerprint: 1,
    };
    let sealed = sidecar::seal(SidecarKind::Bbi, BBI_SIDECAR_VERSION, source, &body);
    assert_eq!(BbiSidecar::parse(&sealed, source).unwrap(), sidecar);
    sidecar.levels.push(LevelBlocks {
        reduction: 1,
        blocks: Vec::new(),
    });
    assert!(std::panic::catch_unwind(|| sidecar.body()).is_err());
    assert!(std::panic::catch_unwind(|| sidecar.to_bytes(source)).is_err());
}

#[test]
fn caller_constructed_block_range_rejects_overflow_in_every_profile() {
    let block = BlockInfo {
        start: (0, 0),
        end: (0, 1),
        offset: u64::MAX,
        size: 1,
        inflated: 1,
        items: 1,
    };
    assert!(std::panic::catch_unwind(|| block.range()).is_err());
    let block = BlockInfo { size: 0, ..block };
    assert_eq!(block.range(), u64::MAX..u64::MAX);
}

mod source_identity {
    use crate::support::Opened;
    use crate::support::libbigwig::Writer;
    use omicsweave_bbi::{BbiSidecar, Kind};
    use omicsweave_core::Error;
    use omicsweave_core::sidecar::sidecar_path;
    #[test]
    fn sidecar_of_a_libbigwig_file_and_its_path() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("lib.bw");
        {
            let mut w = Writer::create(&path, 6, &[("chr1", 100_000), ("chr2", 50_000)]).unwrap();
            let items: Vec<(u32, u32, f32)> =
                (0..5_000).map(|i| (i * 20, i * 20 + 7, i as f32)).collect();
            assert!(w.add_bedgraph("chr1", &items));
            assert!(w.add_fixedstep("chr2", 100, 5, 10, &vec![1.5; 2_000]));
        }
        let built = BbiSidecar::build_file(&path).unwrap();
        assert_eq!(BbiSidecar::load(&path).unwrap(), built);
        assert!(sidecar_path(&path).ends_with("lib.bw.sidecar"));

        let mut ours = Opened::open(&path);
        assert_eq!(built.kind, Kind::BigWig);
        assert_eq!(built.levels.len(), 1 + ours.zooms.len());
        let total: u32 = built.levels[0].blocks.iter().map(|b| b.items).sum();
        assert_eq!(total, 5_000 + 2_000);
        for (i, level) in built.levels.iter().enumerate() {
            for b in &level.blocks {
                assert!(b.inflated <= built.uncompress_buf_size);
            }
            if i > 0 {
                let records: u32 = level.blocks.iter().map(|b| b.items).sum();
                let found: usize = (0..2)
                    .map(|c| ours.summaries(i - 1, c, 0, u32::MAX).len())
                    .sum();
                assert_eq!(records as usize, found);
            }
        }

        // Any change to the file makes the saved sidecar stale.
        let mut bytes = std::fs::read(&path).unwrap();
        bytes.push(0);
        std::fs::write(&path, &bytes).unwrap();
        assert_eq!(BbiSidecar::load(&path), Err(Error::StaleSidecar));
    }
}
