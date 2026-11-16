//! Public in-memory and sidecar index construction over htslib-written BCF.

use crate::indexed_fixture::fixture;
use omicsweave_bgzf::{BinIndex, IndexKind};
use omicsweave_core::exec::Inline;
use omicsweave_core::{FileSource, Region};
use omicsweave_vcf::{BcfReader, RecordQuery, build_index, build_index_file};

const HEADER: &str = "##fileformat=VCFv4.4\n\
    ##contig=<ID=chr1,length=2000000>\n\
    ##contig=<ID=chr2,length=2000000>\n\
    ##INFO=<ID=END,Number=1,Type=Integer>\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

fn text() -> String {
    format!(
        "{HEADER}\
        chr1\t1\tfirst\tA\tC\t.\tPASS\t.\n\
        chr1\t100\tdeletion\tA\t<DEL>\t.\tPASS\tEND=120\n\
        chr1\t70001\tlast\tAT\tA\t.\tPASS\t.\n\
        chr2\t5\tother\tA\tG\t.\tPASS\t.\n"
    )
}

#[test]
fn built_indexes_find_boundary_overlaps_and_match_htslib_metadata() {
    for min_shift in [10, 14, 17] {
        let fixture = fixture(&text(), min_shift as i32);
        let bcf = std::fs::read(&fixture.bcf).unwrap();
        let expected = BinIndex::read_file(&fixture.csi).unwrap();
        let bytes = build_index(&bcf[..], Box::new(Inline::new()), min_shift).unwrap();
        let built = BinIndex::parse(&bytes).unwrap();
        assert_eq!(built.kind(), IndexKind::Csi);
        assert_eq!(
            (built.min_shift(), built.levels(), built.n_refs()),
            (expected.min_shift(), expected.levels(), 2)
        );
        for rid in 0..2 {
            assert_eq!(built.counts(rid), expected.counts(rid));
        }
        let mut reader = BcfReader::open(FileSource::open(&fixture.bcf).unwrap(), built).unwrap();
        for (chrom, start, end, positions) in [
            ("chr1", 0, 1, vec![0]),
            ("chr1", 98, 99, vec![]),
            ("chr1", 119, 120, vec![99]),
            ("chr1", 120, 121, vec![]),
            ("chr1", 70001, 70002, vec![70000]),
            ("chr2", 4, 5, vec![4]),
        ] {
            let table = reader
                .records(
                    Region::new(chrom, start, end),
                    &[],
                    &RecordQuery::default(),
                    None,
                )
                .unwrap();
            assert_eq!(
                table.genotypes.start, positions,
                "min_shift={min_shift}, {chrom}:{start}-{end}"
            );
        }
        // The file entry point replaces the selected sidecar atomically.
        for threads in [1, 2] {
            let path = build_index_file(&fixture.bcf, min_shift, threads).unwrap();
            assert_eq!(path, fixture.csi);
            let on_disk = BinIndex::read_file(path).unwrap();
            assert_eq!(on_disk.counts(0), Some((3, 0)));
            assert_eq!(on_disk.counts(1), Some((1, 0)));
        }
    }
}

#[test]
fn empty_bcf_builds_a_queryable_index_and_bad_input_preserves_existing_sidecar() {
    let fixture = fixture(HEADER, 14);
    let bcf = std::fs::read(&fixture.bcf).unwrap();
    let index = build_index(&bcf[..], Box::new(Inline::new()), 14).unwrap();
    let mut reader = BcfReader::open(
        FileSource::open(&fixture.bcf).unwrap(),
        BinIndex::parse(&index).unwrap(),
    )
    .unwrap();
    let empty = reader
        .records(
            Region::new("chr1", 0, 100),
            &[],
            &RecordQuery::default(),
            None,
        )
        .unwrap();
    assert!(empty.genotypes.is_empty());
    assert_eq!(reader.index().n_refs(), 2);
    for input in [&b"not a BGZF BCF"[..], &bcf[..10]] {
        assert!(build_index(input, Box::new(Inline::new()), 14).is_err());
    }
    assert!(build_index(&bcf[..], Box::new(Inline::new()), u32::MAX).is_err());
    let saved = std::fs::read(&fixture.csi).unwrap();
    std::fs::write(&fixture.bcf, &bcf[..10]).unwrap();
    assert!(build_index_file(&fixture.bcf, 14, 1).is_err());
    assert_eq!(std::fs::read(&fixture.csi).unwrap(), saved);
}
