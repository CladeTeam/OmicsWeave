use crate::support::*;
use omicsweave_core::{Error, FileSource};
use omicsweave_cram::{CramHeader, Scanner};
use std::path::Path;

#[test]
fn queries_match_htslib_across_versions_and_codecs() {
    let f = fixture(
        1,
        standard_refs(),
        SamSpec {
            per_ref: 3_000,
            gap: 50,
            long_cigars: true,
            unplaced: 20,
        },
    );
    let configs: [(&str, &[&str]); 9] = [
        ("v30", &["version=3.0"]),
        ("v31", &["version=3.1"]),
        ("gzip", &["version=3.0", "use_rans=0", "level=1"]),
        (
            "bz2lzma",
            &[
                "version=3.0",
                "use_rans=0",
                "use_bzip2=1",
                "use_lzma=1",
                "level=9",
            ],
        ),
        (
            "arith",
            &[
                "version=3.1",
                "use_arith=1",
                "use_fqz=1",
                "use_tok=1",
                "level=9",
            ],
        ),
        ("archive", &["version=3.1", "archive", "level=9"]),
        (
            "small",
            &[
                "version=3.1",
                "seqs_per_slice=300",
                "slices_per_container=3",
            ],
        ),
        ("fast", &["version=3.1", "fast", "use_tok=0"]),
        (
            "fqz",
            &["version=3.1", "use_rans=0", "use_fqz=1", "level=1"],
        ),
    ];
    let mut methods = std::collections::BTreeSet::new();
    for (i, (name, opts)) in configs.iter().enumerate() {
        eprintln!(
            "CRAM config={name}, options={opts:?}, query_seed={}",
            100 + i
        );
        let (cram, crai) = f.cram(name, opts, true);
        let covered: std::collections::BTreeSet<_> = block_methods(&cram).into_iter().collect();
        eprintln!("CRAM config={name}, block_methods={covered:?}");
        methods.extend(covered);
        let n = check_queries(&f, &cram, &crai, true, 100 + i as u64);
        assert!(n > 5_000, "{name}: {n} records");
    }
    // Every block compression method was decoded at least once.
    assert_eq!(
        methods.into_iter().collect::<Vec<_>>(),
        (0..=8).collect::<Vec<u8>>(),
        "methods covered"
    );
}

#[test]
fn file_definition_rejects_truncation_magic_and_unsupported_versions() {
    let mut definition = [0u8; 26];
    definition[..6].copy_from_slice(b"CRAM\x03\x01");
    for length in 0..definition.len() {
        assert!(
            matches!(
                CramHeader::read(&definition[..length]),
                Err(Error::Truncated { need: 26, have, .. }) if have == length as u64
            ),
            "definition length={length}"
        );
    }
    for (major, minor) in [(2, 1), (3, 2), (4, 0)] {
        definition[4] = major;
        definition[5] = minor;
        assert!(
            matches!(
                CramHeader::read(&definition[..]),
                Err(Error::Unsupported(_))
            ),
            "version={major}.{minor}"
        );
    }
    definition[..4].copy_from_slice(b"BAM\x01");
    assert!(matches!(
        CramHeader::read(&definition[..]),
        Err(Error::BadMagic { .. })
    ));
    definition[..6].copy_from_slice(b"CRAM\x03\x01");
    assert!(
        matches!(
            CramHeader::read(&definition[..]),
            Err(Error::Invalid { .. })
        ),
        "a complete definition still needs a header container"
    );
}

/// Compression methods of every block of a file.
fn block_methods(cram: &Path) -> Vec<u8> {
    let source = FileSource::open(cram).unwrap();
    let header = omicsweave_cram::CramHeader::read(&source).unwrap();
    let mut scan = Scanner::new(source, &header, true);
    let mut out = Vec::new();
    while let Some(c) = scan.next_container().unwrap() {
        out.extend_from_slice(c.block_methods());
    }
    out
}
