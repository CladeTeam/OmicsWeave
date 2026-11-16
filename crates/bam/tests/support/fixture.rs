use super::{Rng, SamSpec, bam_of, index, random_sam};
use std::path::PathBuf;

pub struct Fixture {
    _dir: tempfile::TempDir,
    pub bam: PathBuf,
    pub bai: PathBuf,
    pub csi: PathBuf,
    pub csi12: PathBuf,
    pub refs: Vec<(&'static str, u64)>,
}

pub fn fixture(seed: u64, per_ref: usize, gap: u64, level: u32) -> Fixture {
    let dir = tempfile::tempdir().unwrap();
    let refs = vec![
        ("chr1", 3_000_000),
        ("chr2", 1_200_000),
        ("chrM", 16_569),
        ("empty", 50_000),
    ];
    let mut rng = Rng::new(seed);
    let mut spec_refs = refs.clone();
    spec_refs.pop();
    let text = random_sam(
        &mut rng,
        &SamSpec {
            refs: spec_refs,
            per_ref,
            gap,
            long_cigars: true,
        },
    );
    // The empty reference stays in the header without records.
    let mut header_fix = String::from_utf8(text).unwrap();
    let rg = header_fix.find("@RG\t").unwrap();
    header_fix.insert_str(rg, "@SQ\tSN:empty\tLN:50000\n");
    let sam = dir.path().join("in.sam");
    std::fs::write(&sam, header_fix).unwrap();
    let bam = dir.path().join("in.bam");
    bam_of(&sam, &bam, level);
    let bai = dir.path().join("in.bam.bai");
    index(&bam, &bai, 0);
    let csi = dir.path().join("in.bam.csi");
    index(&bam, &csi, 14);
    let csi12 = dir.path().join("in12.csi");
    index(&bam, &csi12, 12);
    Fixture {
        _dir: dir,
        bam,
        bai,
        csi,
        csi12,
        refs,
    }
}
