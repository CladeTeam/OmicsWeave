use crate::support::{Hts, Rng, fixture, open, query_all, regions};
use crate::support::{parse_index, with_empty_blocks};
use omicsweave_bam::{IndexChoice, build_index};
use omicsweave_core::Region;
use omicsweave_core::exec::Inline;

#[test]
fn queries_match_htslib_through_bai_and_csi() {
    for (seed, per_ref, gap, level) in [(1u64, 6_000, 30, 6), (2, 2_500, 400, 1), (3, 12_000, 5, 0)]
    {
        let f = fixture(seed, per_ref, gap, level);
        let mut rng = Rng::new(seed + 100);
        let regions = regions(&mut rng, &f.refs, 120);
        let blocks = crate::support::blocks(&std::fs::read(&f.bam).unwrap()).len();
        assert!(blocks > 50, "{blocks} blocks");
        for idx in [&f.bai, &f.csi, &f.csi12] {
            let hts = Hts::open(&f.bam, idx);
            let mut ours = open(&f.bam, idx);
            let (mut total, mut nonempty) = (0, 0);
            for &(chrom, start, end) in &regions {
                let want = hts.query(hts.tid(chrom), start as i64, end as i64);
                let got = query_all(&mut ours, Region::new(chrom, start, end));
                assert_eq!(
                    got.len(),
                    want.len(),
                    "seed={seed}, per_ref={per_ref}, gap={gap}, level={level}, index={idx:?}, region={chrom}:{start}-{end}"
                );
                assert_eq!(
                    got, want,
                    "seed={seed}, per_ref={per_ref}, gap={gap}, level={level}, index={idx:?}, region={chrom}:{start}-{end}"
                );
                total += got.len();
                nonempty += !got.is_empty() as usize;
            }
            assert!(
                total > 20_000 && nonempty > 60,
                "seed={seed}, per_ref={per_ref}, gap={gap}, level={level}, index={idx:?}: {total} records, {nonempty} regions"
            );
        }
    }
}

#[test]
fn empty_blocks_inside_the_file() {
    let f = fixture(21, 6_000, 25, 6);
    let spliced = f.bam.with_file_name("spliced.bam");
    with_empty_blocks(&f.bam, &spliced, 3);
    let bai = f.bam.with_file_name("spliced.bam.bai");
    crate::support::index(&spliced, &bai, 0);
    let raw = std::fs::read(&spliced).unwrap();
    let (_, ours) = build_index(&raw[..], Box::new(Inline::new()), IndexChoice::Bai).unwrap();
    assert_eq!(
        parse_index(&ours),
        parse_index(&std::fs::read(&bai).unwrap())
    );
    let hts = Hts::open(&spliced, &bai);
    let mut r = open(&spliced, &bai);
    let mut rng = Rng::new(5);
    for (chrom, start, end) in regions(&mut rng, &f.refs, 60) {
        let want = hts.query(hts.tid(chrom), start as i64, end as i64);
        assert_eq!(query_all(&mut r, Region::new(chrom, start, end)), want);
    }
}
