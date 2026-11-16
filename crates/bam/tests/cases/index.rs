use crate::support::{fixture, parse_index};
use omicsweave_bam::{BinIndex, IndexChoice, IndexKind, build_index};
use omicsweave_core::exec::{Inline, Pool};

#[test]
fn built_indexes_match_samtools_index() {
    for (seed, level) in [(11u64, 6), (12, 0), (13, 9)] {
        let f = fixture(seed, 8_000, 20, level);
        let theirs_bai = parse_index(&std::fs::read(&f.bai).unwrap());
        let theirs_csi = parse_index(&std::fs::read(&f.csi).unwrap());
        let theirs_csi12 = parse_index(&std::fs::read(&f.csi12).unwrap());
        let raw = std::fs::read(&f.bam).unwrap();
        let pool = Pool::new(4);
        for (choice, kind, want) in [
            (IndexChoice::Bai, IndexKind::Bai, &theirs_bai),
            (IndexChoice::Auto, IndexKind::Bai, &theirs_bai),
            (
                IndexChoice::Csi { min_shift: 14 },
                IndexKind::Csi,
                &theirs_csi,
            ),
            (
                IndexChoice::Csi { min_shift: 12 },
                IndexKind::Csi,
                &theirs_csi12,
            ),
        ] {
            let (k, bytes) = build_index(&raw[..], Box::new(Inline::new()), choice).unwrap();
            assert_eq!(k, kind);
            assert_eq!(&parse_index(&bytes), want, "{choice:?}");
            let (_, threaded) = build_index(&raw[..], Box::new(pool.ordered(8)), choice).unwrap();
            assert_eq!(threaded, bytes);
            // Our index reads back as htslib's does.
            BinIndex::parse(&bytes).unwrap();
        }
    }
}
