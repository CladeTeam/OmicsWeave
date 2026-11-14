//! BED → bigBed: the same file and sidecar as `BigBedWriter` given the
//! records sorted by (chromosome in order of first appearance, start, input
//! order), whether the sort stayed in memory or spilled.

mod common;

use std::io::Cursor;

use common::{Rng, STORAGES, budgets, write_input};
use omicsweave_bbi::{BbiSidecar, BedSchema, BigBedWriter, Chrom, WriteOptions};
use omicsweave_convert::{BedToBigBed, Budget, ChromSizes, Conversion, convert};

struct Rec {
    chrom: usize,
    start: u32,
    end: u32,
    rest: String,
}

fn random_bed(rng: &mut Rng, n: usize, extra: usize) -> (Vec<Rec>, Vec<(String, u32)>) {
    let chroms: Vec<(String, u32)> = (0..rng.range(1, 6))
        .map(|i| {
            (
                format!("chr{}", [7, 1, 12, 3, 9][i as usize]),
                rng.range(50_000, 2_000_000) as u32,
            )
        })
        .collect();
    let recs = (0..n)
        .map(|i| {
            let chrom = rng.range(0, chroms.len() as u64) as usize;
            let len = chroms[chrom].1;
            let start = rng.range(0, len as u64 - 10) as u32;
            let end = (start + rng.range(0, 5_000) as u32).min(len);
            let rest = (0..extra)
                .map(|k| match k {
                    0 => format!("peak{i}"),
                    1 => rng.range(0, 1000).to_string(),
                    2 => (*rng.pick(&["+", "-", "."])).to_string(),
                    _ => format!("x{}", rng.next() % 97),
                })
                .collect::<Vec<_>>()
                .join("\t");
            Rec {
                chrom,
                start,
                end,
                rest,
            }
        })
        .collect();
    (recs, chroms)
}

fn bed_text(recs: &[Rec], chroms: &[(String, u32)], extra: usize) -> Vec<u8> {
    let mut out = String::from("track name=test\n# comment\n");
    for (i, r) in recs.iter().enumerate() {
        out.push_str(&format!("{}\t{}\t{}", chroms[r.chrom].0, r.start, r.end));
        if extra > 0 {
            out.push('\t');
            out.push_str(&r.rest);
        }
        out.push_str(if i % 5 == 0 { "\r\n" } else { "\n" });
        if i % 1000 == 7 {
            out.push('\n');
        }
    }
    out.into_bytes()
}

fn expected(
    recs: &[Rec],
    chroms: &[(String, u32)],
    extra: usize,
    sizes: bool,
) -> (Vec<u8>, BbiSidecar) {
    let mut order: Vec<usize> = Vec::new();
    for r in recs {
        if !order.contains(&r.chrom) {
            order.push(r.chrom);
        }
    }
    let id = |c: usize| order.iter().position(|&o| o == c).unwrap() as u32;
    let mut sorted: Vec<(u32, u32, usize)> = recs
        .iter()
        .enumerate()
        .map(|(i, r)| (id(r.chrom), r.start, i))
        .collect();
    sorted.sort();
    let out_chroms: Vec<Chrom> = order
        .iter()
        .map(|&c| Chrom {
            name: chroms[c].0.clone(),
            len: if sizes {
                chroms[c].1
            } else {
                recs.iter()
                    .filter(|r| r.chrom == c)
                    .map(|r| r.end)
                    .max()
                    .unwrap()
            },
        })
        .collect();
    let fields = 3 + extra as u16;
    let schema = BedSchema {
        field_count: fields,
        defined_field_count: fields.min(12),
        auto_sql: None,
    };
    let mut w = BigBedWriter::create(
        Cursor::new(Vec::new()),
        out_chroms,
        schema,
        WriteOptions::default(),
    )
    .unwrap();
    for (c, _, i) in sorted {
        let r = &recs[i];
        w.add(c, r.start, r.end, r.rest.as_bytes()).unwrap();
    }
    let (out, sidecar) = w.finish().unwrap();
    (out.into_inner(), sidecar)
}

#[test]
fn bigbed_matches_the_writer_in_memory_and_spilled() {
    let dir = tempfile::tempdir().unwrap();
    for seed in 0..3u64 {
        let mut rng = Rng::new(seed);
        let extra = [0, 3, 5][seed as usize];
        let (recs, chroms) = random_bed(&mut rng, 60_000, extra);
        let text = bed_text(&recs, &chroms, extra);
        let sizes_path = dir.path().join("chrom.sizes");
        let sizes_text: String = chroms.iter().map(|(n, l)| format!("{n}\t{l}\n")).collect();
        std::fs::write(&sizes_path, sizes_text).unwrap();
        for with_sizes in [false, true] {
            let (want, want_sidecar) = expected(&recs, &chroms, extra, with_sizes);
            for storage in STORAGES {
                let src = dir.path().join("in.bed");
                write_input(&src, &text, storage);
                let mut all = budgets();
                all.push(Budget::new(4, 8 << 20));
                for budget in all {
                    let dst = dir.path().join("out.bigbed");
                    let sizes = if with_sizes {
                        ChromSizes::File(sizes_path.clone())
                    } else {
                        ChromSizes::FromData
                    };
                    let report = convert(
                        Conversion::BedToBigBed(BedToBigBed::new(&src, &dst, sizes)),
                        budget.clone(),
                        None,
                    )
                    .unwrap();
                    assert_eq!(report.records, recs.len() as u64);
                    assert_eq!(std::fs::read(&dst).unwrap(), want, "{storage:?} {budget:?}");
                    assert_eq!(BbiSidecar::load(&dst).unwrap(), want_sidecar);
                    if budget.memory <= 8 << 20 {
                        assert!(report.spilled > 0, "a tight budget spills");
                    }
                }
            }
        }
    }
}

#[test]
fn records_past_the_chromosome_end_are_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("a.bed");
    std::fs::write(&src, b"chr1\t0\t10\nchr1\t5\t120\n").unwrap();
    let sizes = dir.path().join("sizes");
    std::fs::write(&sizes, b"chr1\t100\n").unwrap();
    let dst = dir.path().join("a.bb");
    let err = convert(
        Conversion::BedToBigBed(BedToBigBed::new(&src, &dst, ChromSizes::File(sizes))),
        Budget::new(2, 1 << 28),
        None,
    )
    .unwrap_err();
    assert!(err.to_string().contains("line 2"), "{err}");
    assert!(!dst.exists());
}
