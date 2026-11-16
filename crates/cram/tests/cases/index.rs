use crate::support::*;
use omicsweave_core::FileSource;
use omicsweave_cram::CramIndex;

#[test]
fn indexes_match_samtools_index() {
    let f = fixture(
        4,
        standard_refs(),
        SamSpec {
            per_ref: 2_000,
            gap: 30,
            long_cigars: true,
            unplaced: 30,
        },
    );
    for (name, opts) in [
        ("a", &["version=3.0"][..]),
        (
            "b",
            &[
                "version=3.1",
                "seqs_per_slice=200",
                "slices_per_container=4",
            ],
        ),
        ("c", &["version=3.1", "no_ref=1"]),
    ] {
        let (cram, crai) = f.cram(name, opts, name != "c");
        let want = gunzip(&std::fs::read(&crai).unwrap());
        let built = omicsweave_cram::build_index(FileSource::open(&cram).unwrap()).unwrap();
        assert_eq!(gunzip(&built), want, "{name}");
        let header = omicsweave_cram::CramHeader::read(&FileSource::open(&cram).unwrap()).unwrap();
        let parsed = CramIndex::parse(&built, header.references().len()).unwrap();
        assert_eq!(
            parsed.entries().len(),
            want.split(|&c| c == b'\n')
                .filter(|l| !l.is_empty())
                .count()
        );
    }
}
