//! CRAM → BAM against htslib: the uncompressed BAM stream is the one
//! htslib writes reading the CRAM, whatever the budget, and the index holds
//! what `samtools index` builds for the output.

use crate::cram_support as cram_common;
use crate::support as common;

use crate::htslib as hts;
use common::rng::Rng;
use common::storage::inflate;
use cram_common::{SamSpec, bam_of_cram, cram_of, random_sam, reference, write_fasta};
use omicsweave_convert::{
    BamIndex, Budget, Conversion, CramToBam, Format, Index, Options, convert,
};

/// Inline and threaded, roomy and with two slices in flight.
fn budgets() -> Vec<Budget> {
    vec![
        Budget::new(1, 1 << 30),
        Budget::new(6, 1 << 30),
        Budget::new(3, 64 << 20),
    ]
}

#[test]
fn bam_matches_htslib_with_its_index() {
    let dir = tempfile::tempdir().unwrap();
    let refs: Vec<(String, u64)> = vec![
        ("chr1".into(), 400_000),
        ("chr2".into(), 120_000),
        ("chrM".into(), 16_569),
    ];
    let mut rng = Rng::new(11);
    let seqs = reference(&mut rng, &refs);
    let fasta = dir.path().join("ref.fa");
    write_fasta(&fasta, &refs, &seqs);
    let sam = dir.path().join("in.sam");
    let spec = SamSpec {
        per_ref: 2_500,
        gap: 40,
        long_cigars: true,
        unplaced: 25,
    };
    std::fs::write(&sam, random_sam(&mut rng, &refs, &seqs, &spec)).unwrap();
    for (name, opts, with_ref) in [
        ("v31", &["version=3.1", "seqs_per_slice=700"][..], true),
        ("v30", &["version=3.0", "lossy_names=1"], true),
        ("noref", &["version=3.1", "no_ref=1"], false),
    ] {
        let cram = dir.path().join(format!("{name}.cram"));
        cram_of(&sam, &cram, with_ref.then_some(fasta.as_path()), opts);
        assert_eq!(Format::detect(&cram).unwrap(), Format::Cram);
        let theirs_bam = dir.path().join(format!("{name}.theirs.bam"));
        bam_of_cram(&cram, with_ref.then_some(fasta.as_path()), &theirs_bam, 6);
        let want = inflate(&theirs_bam);
        let mut first: Option<Vec<u8>> = None;
        for budget in budgets() {
            let dst = dir.path().join(format!("{name}.bam"));
            let mut c = Conversion::between(Format::Cram, Format::Bam, &cram, &dst).unwrap();
            c.configure(&Options {
                index: Some(Index::Auto),
                reference: with_ref.then(|| fasta.clone()),
                ..Options::default()
            })
            .unwrap();
            let report = convert(c, budget.clone(), None).unwrap();
            assert_eq!(report.files[0], (Format::Bam, dst.clone()));
            assert_eq!(inflate(&dst), want, "{name} {budget:?}");
            let ours = dir.path().join(format!("{name}.bam.bai"));
            let theirs = dir.path().join("theirs.bai");
            hts::bam_index(&dst, &theirs, 0);
            assert_eq!(hts::parse_index(&ours), hts::parse_index(&theirs), "{name}");
            let bytes = std::fs::read(&dst).unwrap();
            match &first {
                None => first = Some(bytes),
                Some(f) => assert_eq!(f, &bytes, "{name} {budget:?}"),
            }
        }
        // The .crai written by htslib is told by content.
        let crai = dir.path().join(format!("{name}.cram.crai"));
        cram_common::crai_of(&cram, &crai);
        assert_eq!(Format::detect(&crai).unwrap(), Format::Crai);
    }
}

#[test]
fn reference_is_an_option_of_cram_only() {
    let mut c = Conversion::CramToBam(CramToBam::new("a.cram", "b.bam"));
    assert!(
        c.configure(&Options {
            sort: true,
            ..Options::default()
        })
        .is_err()
    );
    let mut s = Conversion::infer("a.sam", "b.bam").unwrap();
    assert!(
        s.configure(&Options {
            reference: Some("ref.fa".into()),
            ..Options::default()
        })
        .is_err()
    );
    c.configure(&Options {
        reference: Some("ref.fa".into()),
        index: Some(Index::Csi),
        ..Options::default()
    })
    .unwrap();
    let Conversion::CramToBam(c) = c else {
        unreachable!()
    };
    assert_eq!(c.index, BamIndex::Csi { min_shift: 14 });
    assert_eq!(
        Format::from_extension("x.cram".as_ref()).unwrap(),
        Format::Cram
    );
    assert_eq!(
        Format::from_extension("x.cram.crai".as_ref()).unwrap(),
        Format::Crai
    );
}
