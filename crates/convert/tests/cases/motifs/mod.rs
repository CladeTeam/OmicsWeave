//! Motif collections → FFM: every source format keeps its values, IDs and
//! names; content tells the formats apart (and FASTA from motifs); several
//! form one collection; malformed or conflicting sources are errors.

use std::path::{Path, PathBuf};

use omicsweave_convert::{Budget, Conversion, Format, convert, detect_sources};
use omicsweave_ffm::{FactorKind, FfmReader};

/// A factor: kind, positions, table.
type FactorRead = (FactorKind, Vec<u32>, Vec<f64>);

/// A motif: ID, name, width, sites, factors.
type MotifRead = (String, String, u32, Option<f64>, Vec<FactorRead>);

struct Read {
    alphabet: Vec<u8>,
    motifs: Vec<MotifRead>,
}

fn read(path: &Path) -> Read {
    let r = FfmReader::open(std::fs::read(path).unwrap()).unwrap();
    let motifs = (0..r.len())
        .map(|i| {
            let m = r.motif(i).unwrap();
            let v = m.view();
            (
                r.id(i).unwrap().to_string(),
                r.annotation(i).unwrap().get(0).unwrap().to_string(),
                v.width(),
                v.sites(),
                v.factors()
                    .map(|f| (f.kind, f.positions.to_vec(), f.table.to_vec()))
                    .collect(),
            )
        })
        .collect();
    Read {
        alphabet: r.alphabet().symbols().to_vec(),
        motifs,
    }
}

fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p
}

/// Converts `sources` (detected) into `out.ffm` and reads it back.
fn run(dir: &Path, sources: &[PathBuf]) -> (Read, Vec<String>) {
    let from = detect_sources(sources).unwrap();
    let dst = dir.join("out.ffm");
    let mut c = Conversion::between(from, Format::Ffm, sources[0].clone(), &dst).unwrap();
    c.add_sources(sources[1..].iter().cloned()).unwrap();
    let report = convert(c, Budget::new(2, 1 << 30), None).unwrap();
    (read(&dst), report.warnings)
}

const MEME: &str = "MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.3 C 0.2 G 0.2 T 0.3

MOTIF MA0002.1 RUNX1
letter-probability matrix: alength= 4 w= 2 nsites= 26 E= 0
 0.25 0.25 0.25 0.25
 0.1 0.2 0.3 0.4
URL http://example.org/MA0002.1

MOTIF lone
letter-probability matrix: alength= 4 w= 1
1 0 0 0
";

mod collections;
mod matrices;
mod potentials;
