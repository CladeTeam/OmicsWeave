//! Reference and file fixtures shared by the behavioral cases.

use super::{Rng, SamSpec, crai_of, cram_of, path_cstring, random_sam};
use htslib_sys as h;
use omicsweave_core::FileSource;
use omicsweave_fasta::{Fai, FastaReader};
use std::path::{Path, PathBuf};

/// Random reference sequences: mostly ACGT, with N runs and soft-masked
/// (lowercase) stretches.
pub fn reference(rng: &mut Rng, refs: &[(String, u64)]) -> Vec<Vec<u8>> {
    refs.iter()
        .map(|(_, len)| {
            let mut s = Vec::with_capacity(*len as usize);
            while (s.len() as u64) < *len {
                let run = rng.range(1, 400) as usize;
                let kind = rng.range(0, 100);
                for _ in 0..run {
                    let b = b"ACGT"[rng.range(0, 4) as usize];
                    s.push(match kind {
                        0..=1 => b'N',
                        2..=12 => b.to_ascii_lowercase(),
                        _ => b,
                    });
                }
            }
            s.truncate(*len as usize);
            s
        })
        .collect()
}

/// Writes FASTA (60 bases a line) and its `.fai` (htslib's).
pub fn write_fasta(path: &Path, refs: &[(String, u64)], seqs: &[Vec<u8>]) {
    let mut text = Vec::new();
    for ((name, _), s) in refs.iter().zip(seqs) {
        text.extend_from_slice(format!(">{name}\n").as_bytes());
        for line in s.chunks(60) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    std::fs::write(path, text).unwrap();
    unsafe {
        assert_eq!(
            h::fai_build3(
                path_cstring(path).as_ptr(),
                std::ptr::null(),
                std::ptr::null()
            ),
            0
        );
    }
}

/// Our reader of the FASTA written by [`write_fasta`].
pub fn fasta_reader(path: &Path) -> FastaReader<FileSource> {
    let mut fai = path.as_os_str().to_owned();
    fai.push(".fai");
    let fai = Fai::parse(&std::fs::read(PathBuf::from(fai)).unwrap()).unwrap();
    FastaReader::new(FileSource::open(path).unwrap(), fai, None)
}

pub struct Fixture {
    pub dir: tempfile::TempDir,
    pub refs: Vec<(String, u64)>,
    pub fasta: PathBuf,
    pub sam: PathBuf,
}

pub fn fixture(seed: u64, refs: Vec<(String, u64)>, spec: SamSpec) -> Fixture {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(seed);
    let seqs = reference(&mut rng, &refs);
    let fasta = dir.path().join("ref.fa");
    write_fasta(&fasta, &refs, &seqs);
    let text = random_sam(&mut rng, &refs, &seqs, &spec);
    let sam = dir.path().join("in.sam");
    std::fs::write(&sam, text).unwrap();
    Fixture {
        dir,
        refs,
        fasta,
        sam,
    }
}

pub fn standard_refs() -> Vec<(String, u64)> {
    vec![
        ("chr1".into(), 600_000),
        ("chr2".into(), 250_000),
        ("chrM".into(), 16_569),
        ("empty".into(), 30_000),
    ]
}

impl Fixture {
    /// The CRAM written with `options` and its htslib index.
    pub fn cram(&self, name: &str, options: &[&str], with_ref: bool) -> (PathBuf, PathBuf) {
        let cram = self.dir.path().join(format!("{name}.cram"));
        let reference = with_ref.then_some(self.fasta.as_path());
        cram_of(&self.sam, &cram, reference, options);
        let crai = self.dir.path().join(format!("{name}.cram.crai"));
        crai_of(&cram, &crai);
        (cram, crai)
    }
}
