//! FASTA → bgzip FASTA and FASTA → 2bit: the same bytes as the format
//! crates' own writers, for every input storage and budget.

use crate::support as common;

use std::io::Cursor;

use common::budgets;
use common::rng::Rng;
use common::storage::{STORAGES, Storage, write_input};
use omicsweave_convert::{Conversion, FastaToBgzipFasta, FastaToTwoBit, convert};
use omicsweave_fasta::{FastaWriter, TwoBitWriter, compress_fasta, encode_twobit};

fn random_records(rng: &mut Rng, max_len: u64) -> Vec<(String, Vec<u8>)> {
    (0..rng.range(1, 7))
        .map(|i| {
            let len = rng.range(1, max_len) as usize;
            let alphabet: &[u8] = if rng.chance(50) {
                b"ACGTNacgtn"
            } else {
                b"ACGTRYKM"
            };
            let mut seq = Vec::with_capacity(len);
            while seq.len() < len {
                let c = *rng.pick(alphabet);
                let run = rng.range(1, 200) as usize;
                seq.extend(std::iter::repeat_n(c, run.min(len - seq.len())));
            }
            (format!("seq{i}_{}", rng.next() % 1000), seq)
        })
        .collect()
}

fn fasta_text(records: &[(String, Vec<u8>)], width: usize) -> Vec<u8> {
    let mut w = FastaWriter::new(Vec::new(), width);
    for (n, s) in records {
        w.add(n, s).unwrap();
    }
    w.finish().unwrap().0
}

#[test]
fn bgzip_matches_compress_fasta() {
    let dir = tempfile::tempdir().unwrap();
    for seed in 0..4u64 {
        let mut rng = Rng::new(seed);
        let records = random_records(&mut rng, 300_000);
        let text = fasta_text(&records, [0, 60, 80][seed as usize % 3]);
        let (want, want_fai, want_gzi) = compress_fasta(&text[..], Vec::new(), 4).unwrap();
        for storage in STORAGES {
            let src = dir.path().join(format!("in{seed}.fa"));
            write_input(&src, &text, storage);
            for budget in budgets() {
                let dst = dir.path().join("out.fa.gz");
                let mut conversion = FastaToBgzipFasta::new(&src, &dst);
                conversion.level = 4;
                let report = convert(
                    Conversion::FastaToBgzipFasta(conversion),
                    budget.clone(),
                    None,
                )
                .unwrap();
                assert_eq!(report.records, records.len() as u64);
                assert_eq!(std::fs::read(&dst).unwrap(), want, "{storage:?} {budget:?}");
                assert_eq!(
                    std::fs::read(dir.path().join("out.fa.gz.fai")).unwrap(),
                    want_fai.to_bytes()
                );
                assert_eq!(
                    std::fs::read(dir.path().join("out.fa.gz.gzi")).unwrap(),
                    want_gzi.to_bytes()
                );
            }
        }
    }
}

fn expected_twobit(records: &[(String, Vec<u8>)]) -> Vec<u8> {
    let names = records.iter().map(|(n, _)| n.clone()).collect();
    let mut w = TwoBitWriter::create(Cursor::new(Vec::new()), names, 0).unwrap();
    for (n, s) in records {
        let (rec, packed) = encode_twobit(s).unwrap();
        w.add(n, &rec, &mut &packed[..]).unwrap();
    }
    w.finish().unwrap().into_inner()
}

#[test]
fn twobit_matches_the_writer() {
    let dir = tempfile::tempdir().unwrap();
    for seed in 0..4u64 {
        let mut rng = Rng::new(100 + seed);
        let records = random_records(&mut rng, 9_000_000);
        let want = expected_twobit(&records);
        let text = fasta_text(&records, 70);
        for storage in STORAGES {
            for with_fai in [false, true] {
                if with_fai && storage == Storage::Gzip {
                    continue;
                }
                let src = dir.path().join(format!("in{seed}{storage:?}{with_fai}.fa"));
                write_input(&src, &text, storage);
                if with_fai {
                    omicsweave_fasta::build_index_files(&src).unwrap();
                }
                for budget in budgets() {
                    let dst = dir.path().join("out.2bit");
                    let report = convert(
                        Conversion::FastaToTwoBit(FastaToTwoBit::new(&src, &dst)),
                        budget.clone(),
                        None,
                    )
                    .unwrap();
                    assert_eq!(std::fs::read(&dst).unwrap(), want, "{storage:?} {budget:?}");
                    // Under the tight budget the packed-base buffer is under
                    // 1 MiB: a sequence of 4 Mbases (1 MiB packed) spills.
                    let longest = records.iter().map(|(_, s)| s.len()).max().unwrap();
                    if budget.memory < 64 << 20 && longest > 4_000_000 {
                        assert!(
                            report.spilled > 0,
                            "large sequences spill under a tight budget"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn twobit_keeps_an_empty_final_record_with_an_unterminated_header() {
    let dir = tempfile::tempdir().unwrap();
    let records = vec![
        ("a".to_owned(), b"ACGT".to_vec()),
        ("b".to_owned(), Vec::new()),
    ];
    let want = expected_twobit(&records);
    for storage in STORAGES {
        let src = dir.path().join("in.fa");
        write_input(&src, b">a\nACGT\n>b", storage);
        let dst = dir.path().join("out.2bit");
        let report = convert(
            Conversion::FastaToTwoBit(FastaToTwoBit::new(&src, &dst)),
            budgets()[0].clone(),
            None,
        )
        .unwrap();
        assert_eq!(report.records, 2);
        assert_eq!(std::fs::read(&dst).unwrap(), want, "{storage:?}");
    }
}

#[test]
fn a_stale_fai_is_reported() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("a.fa");
    std::fs::write(&src, b">a\nACGT\n>b\nAC\n").unwrap();
    std::fs::write(dir.path().join("a.fa.fai"), b"a\t4\t3\t4\t5\n").unwrap();
    let initial_entries = std::fs::read_dir(dir.path()).unwrap().count();
    let err = convert(
        Conversion::FastaToTwoBit(FastaToTwoBit::new(&src, dir.path().join("a.2bit"))),
        budgets()[0].clone(),
        None,
    )
    .unwrap_err();
    assert!(err.to_string().contains("out of date"), "{err}");
    assert!(!dir.path().join("a.2bit").exists());
    assert_eq!(
        std::fs::read_dir(dir.path()).unwrap().count(),
        initial_entries
    );
}
