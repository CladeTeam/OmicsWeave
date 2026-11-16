//! 4DN pairs and Juicer short text → .hic: the file and sidecar
//! `HicBuilder` writes for the same contacts, whatever the storage, thread
//! count or spilling; rows of unknown chromosomes skipped with a warning;
//! the counts straw reads back.

use crate::support as common;

use std::collections::BTreeMap;
use std::ffi::CString;
use std::fmt::Write as _;
use std::io::Cursor;
use std::path::Path;

use common::rng::Rng;
use common::storage::{STORAGES, Storage, write_input};
use omicsweave_convert::{
    Budget, Conversion, Format, JUICER_RESOLUTIONS, Options, PairsDialect, PairsToHic, convert,
};
use omicsweave_hic::{BuildOptions, Chrom, HicBuilder, HicHeader, HicSidecar};

const CHROMS: [(&str, u64); 3] = [("chr1", 900_000), ("chr2", 610_000), ("chrX", 333_333)];

/// One row: chromosome indexes (`None`: not in the list) and positions as
/// written, with its count.
#[derive(Clone)]
struct Row {
    names: (String, String),
    chroms: (Option<u32>, Option<u32>),
    pos: (u32, u32),
    count: f32,
}

fn random_rows(rng: &mut Rng, n: usize, unknown: bool) -> Vec<Row> {
    let mut rows = Vec::with_capacity(n);
    for _ in 0..n {
        let end = |rng: &mut Rng| -> (String, Option<u32>, u32) {
            if unknown && rng.chance(3) {
                let name = *rng.pick(&["chrM", "!", "chrUn_1"]);
                return (name.into(), None, rng.range(0, 10_000) as u32);
            }
            let c = rng.range(0, CHROMS.len() as u64) as u32;
            // Some positions run past the chromosome's end, which counts there.
            let len = CHROMS[c as usize].1;
            (
                CHROMS[c as usize].0.into(),
                Some(c),
                rng.range(1, len + 2_000) as u32,
            )
        };
        let (n1, c1, p1) = end(rng);
        // Near-diagonal pairs, as Hi-C has them, fold into shared cells.
        let (n2, c2, p2) = if let Some(c) = c1
            && rng.chance(70)
        {
            let len = CHROMS[c as usize].1 as u32;
            let p = (p1 as i64 + rng.range(0, 60_000) as i64 - 30_000).clamp(1, len as i64);
            (n1.clone(), c1, p as u32)
        } else {
            end(rng)
        };
        rows.push(Row {
            names: (n1, n2),
            chroms: (c1, c2),
            pos: (p1, p2),
            count: 1.0,
        });
    }
    rows
}

fn header(resolutions: &[u32], genome: &str) -> HicHeader {
    let mut r = resolutions.to_vec();
    r.sort_unstable_by(|a, b| b.cmp(a));
    HicHeader::new(
        genome,
        CHROMS
            .iter()
            .map(|&(n, l)| Chrom {
                name: n.into(),
                len: l,
            })
            .collect(),
        r,
    )
}

/// What `HicBuilder` writes for the known rows.
fn expected(rows: &[Row], resolutions: &[u32], genome: &str) -> (Vec<u8>, HicSidecar) {
    let mut b = HicBuilder::new(header(resolutions, genome)).unwrap();
    for r in rows {
        if let (Some(a), Some(c)) = r.chroms {
            b.add(a, r.pos.0, c, r.pos.1, r.count).unwrap();
        }
    }
    let (out, sidecar) = b
        .finish_with_sidecar(Cursor::new(Vec::new()), &BuildOptions::default())
        .unwrap();
    (out.into_inner(), sidecar)
}

fn pairs_text(rows: &[Row], genome: &str) -> Vec<u8> {
    let mut t = String::from("## pairs format v1.0\n#sorted: none\n#shape: upper triangle\n");
    writeln!(t, "#genome_assembly: {genome}").unwrap();
    for (n, l) in CHROMS {
        writeln!(t, "#chromsize: {n} {l}").unwrap();
    }
    // Columns in an order of their own, with extra ones.
    t.push_str("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 pair_type mapq1 mapq2\n");
    for (i, r) in rows.iter().enumerate() {
        writeln!(
            t,
            "r{i}\t{}\t{}\t{}\t{}\t+\t-\tUU\t0\t60",
            r.names.0, r.pos.0, r.names.1, r.pos.1
        )
        .unwrap();
    }
    t.into_bytes()
}

fn short_text(rows: &[Row], scored: bool) -> Vec<u8> {
    let mut t = String::new();
    for r in rows {
        write!(
            t,
            "0 {} {} 0 16 {} {} 1",
            r.names.0, r.pos.0, r.names.1, r.pos.1
        )
        .unwrap();
        if scored {
            write!(t, " {}", r.count).unwrap();
        }
        t.push('\n');
    }
    t.into_bytes()
}

fn hicpro_text(rows: &[Row]) -> Vec<u8> {
    let mut text = String::new();
    for (i, row) in rows.iter().enumerate() {
        write!(
            text,
            "r{i}\t{}\t{}\t+\t{}\t{}\t-\t0\tfrag_0\tfrag_0\t0\t0",
            row.names.0, row.pos.0, row.names.1, row.pos.1
        )
        .unwrap();
        if i % 2 == 0 {
            text.push_str("\t1");
        }
        text.push_str("\t\n");
    }
    text.into_bytes()
}

#[test]
fn hicpro_positions_counts_unknown_chromosomes_and_storage() {
    let dir = tempfile::tempdir().unwrap();
    let mut rows = random_rows(&mut Rng::new(46), 20_000, true);
    // Exact bin boundaries distinguish direct use from subtracting one.
    for pos in [(5_000, 10_000), (4_999, 5_000), (1, 1)] {
        rows.push(Row {
            names: ("chr1".into(), "chr2".into()),
            chroms: (Some(0), Some(1)),
            pos,
            count: 1.0,
        });
    }
    let resolutions = vec![5_000, 20_000];
    let options = Options {
        chrom_sizes: Some(sizes_file(dir.path())),
        resolutions: Some(resolutions.clone()),
        ..Options::default()
    };
    let want = expected(&rows, &resolutions, "unknown").0;
    for storage in STORAGES {
        let src = dir.path().join("input.validPairs");
        write_input(&src, &hicpro_text(&rows), storage);
        let from = Format::detect(&src).unwrap();
        assert_eq!(
            from,
            if storage == Storage::Plain {
                Format::HicProPairs
            } else {
                Format::GzipHicProPairs
            }
        );
        for budget in budgets() {
            let dst = dir.path().join("output.hic");
            let report = run(&src, &dst, from, &options, budget).unwrap();
            assert_eq!(std::fs::read(&dst).unwrap(), want);
            assert_eq!(
                report.records,
                rows.iter()
                    .filter(|r| r.chroms.0.is_some() && r.chroms.1.is_some())
                    .count() as u64
            );
            assert_eq!(report.warnings.len(), 1);
            for chrom in ["chrM", "chrUn_1", "!"] {
                assert!(report.warnings[0].contains(chrom));
            }
        }
    }
    for name in ["x.allValidPairs", "x.validPairs"] {
        assert_eq!(
            Format::from_extension(Path::new(name)).unwrap(),
            Format::HicProPairs
        );
        assert_eq!(
            Format::from_extension(Path::new(&format!("{name}.gz"))).unwrap(),
            Format::GzipHicProPairs
        );
    }
    let src = dir.path().join("input.validPairs");
    let error = run(
        &src,
        &dir.path().join("no-sizes.hic"),
        Format::HicProPairs,
        &Options::default(),
        Budget::new(1, 64 << 20),
    )
    .unwrap_err();
    assert!(error.to_string().contains("chromosome sizes"));
}

#[test]
fn hicpro_rejects_invalid_columns_and_strands() {
    let dir = tempfile::tempdir().unwrap();
    let options = Options {
        chrom_sizes: Some(sizes_file(dir.path())),
        resolutions: Some(vec![5_000]),
        ..Options::default()
    };
    for line in [
        "r chr1 1 + chr2 2 - 0 f1 f2 30",
        "r chr1 1 0 chr2 2 - 0 f1 f2 30 30",
        "r chr1 1 + chr2 bad - 0 f1 f2 30 30",
        "r chr1 1 + chr2 2 - 0 f1 f2 30 30 1 extra",
    ] {
        let src = dir.path().join("bad.validPairs");
        std::fs::write(&src, format!("{line}\n")).unwrap();
        let error = run(
            &src,
            &dir.path().join("bad.hic"),
            Format::HicProPairs,
            &options,
            Budget::new(1, 64 << 20),
        )
        .unwrap_err();
        assert!(error.to_string().contains("HiC-Pro rows"), "{error}");
    }
}

#[test]
fn hicpro_full_width_positions_clamp_before_binning() {
    let dir = tempfile::tempdir().unwrap();
    let text = "r1\tchr1\t18446744073709551615\t+\tchr2\t18446744073709551615\t-\t0\tf1\tf2\t0\t0\n\
                r2\tchr2\t00018446744073709551615\t-\tchr1\t18446744073709551615\t+\t0\tf2\tf1\t0\t0\tallele\t\r\n\
                r3\tchr1\t0\t+\tchr2\t18446744073709551615\t-\t0\tf1\tf2\t0\t0";
    let end1 = CHROMS[0].1 as u32;
    let end2 = CHROMS[1].1 as u32;
    let rows: Vec<_> = [(end1, end2), (end1, end2), (0, end2)]
        .into_iter()
        .map(|pos| Row {
            names: ("chr1".into(), "chr2".into()),
            chroms: (Some(0), Some(1)),
            pos,
            count: 1.0,
        })
        .collect();
    let resolutions = [5_000, 20_000];
    let want = expected(&rows, &resolutions, "unknown").0;
    let options = Options {
        chrom_sizes: Some(sizes_file(dir.path())),
        resolutions: Some(resolutions.to_vec()),
        ..Options::default()
    };
    for storage in STORAGES {
        let src = dir.path().join("boundary.validPairs");
        write_input(&src, text.as_bytes(), storage);
        let from = Format::detect(&src).unwrap();
        assert_eq!(
            from,
            if storage == Storage::Plain {
                Format::HicProPairs
            } else {
                Format::GzipHicProPairs
            }
        );
        let dst = dir.path().join("boundary.hic");
        let report = run(&src, &dst, from, &options, Budget::new(2, 64 << 20)).unwrap();
        assert_eq!(report.records, 3);
        assert!(report.warnings.is_empty());
        assert_eq!(std::fs::read(&dst).unwrap(), want);
    }
}

#[test]
fn hicpro_long_chromosome_end_is_readable_by_straw() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("long.allValidPairs");
    let dst = dir.path().join("long.hic");
    let sizes = dir.path().join("chrom.sizes");
    std::fs::write(&sizes, "chrA\t830000000\nchrB\t780000000\n").unwrap();
    std::fs::write(
        &src,
        "r1\tchrA\t829999999\t+\tchrA\t830000000\t-\t0\tf1\tf2\t0\t0\n\
         r2\tchrB\t779999999\t-\tchrA\t829999999\t+\t0\tf1\tf2\t0\t0\n",
    )
    .unwrap();
    let options = Options {
        chrom_sizes: Some(sizes),
        resolutions: Some(vec![5_000, 25_000]),
        ..Options::default()
    };
    let report = run(
        &src,
        &dst,
        Format::HicProPairs,
        &options,
        Budget::new(2, 32 << 20),
    )
    .unwrap();
    assert_eq!(report.records, 2);
    for res in [5_000, 25_000] {
        assert_eq!(
            straw(
                &dst,
                "chrA:829000000:830000001",
                "chrA:829000000:830000001",
                res
            ),
            vec![((829_999_999 / res) * res, 830_000_000, 1.0)]
        );
        assert_eq!(
            straw(
                &dst,
                "chrA:829000000:830000001",
                "chrB:779000000:780000001",
                res
            ),
            vec![((829_999_999 / res) * res, (779_999_999 / res) * res, 1.0)]
        );
    }
}

fn sizes_file(dir: &Path) -> std::path::PathBuf {
    let p = dir.join("chrom.sizes");
    let text: String = CHROMS.iter().map(|(n, l)| format!("{n}\t{l}\n")).collect();
    std::fs::write(&p, text).unwrap();
    p
}

fn run(
    src: &Path,
    dst: &Path,
    from: Format,
    options: &Options,
    budget: Budget,
) -> omicsweave_core::Result<omicsweave_convert::Report> {
    let mut c = Conversion::between(from, Format::Hic, src, dst)?;
    c.configure(options)?;
    convert(c, budget, None)
}

fn budgets() -> Vec<Budget> {
    vec![
        Budget::new(1, 1 << 30),
        Budget::new(6, 1 << 30),
        Budget::new(3, 12 << 20),
    ]
}

#[test]
fn leading_whitespace_and_split_scored_first_rows_preserve_contacts() {
    let dir = tempfile::tempdir().unwrap();
    let row = Row {
        names: ("chr1".into(), "chr1".into()),
        chroms: (Some(0), Some(0)),
        pos: (1000, 2000),
        count: 2.5,
    };
    let resolutions = [5000];
    let expected = expected(std::slice::from_ref(&row), &resolutions, "unknown").0;
    let sizes = sizes_file(dir.path());
    let text = format!(
        " \t\r\n#comment\n0 {}chr1 1000 0 16 chr1 2000 1 2.5\n",
        " ".repeat((1 << 20) + 17)
    );
    let options = Options {
        chrom_sizes: Some(sizes),
        resolutions: Some(resolutions.to_vec()),
        ..Options::default()
    };
    for storage in STORAGES {
        let src = dir.path().join("input.short");
        write_input(&src, text.as_bytes(), storage);
        for budget in [Budget::new(1, 1 << 28), Budget::new(3, 12 << 20)] {
            let dst = dir.path().join("output.hic");
            let report = run(&src, &dst, Format::JuicerShort, &options, budget).unwrap();
            assert_eq!(report.records, 1);
            assert_eq!(std::fs::read(&dst).unwrap(), expected);
        }
    }
    let source = format!(
        " \t\r\n{}",
        String::from_utf8(pairs_text(
            std::slice::from_ref(&Row { count: 1.0, ..row }),
            "unknown"
        ))
        .unwrap()
    );
    let src = dir.path().join("input.pairs");
    std::fs::write(&src, source).unwrap();
    let dst = dir.path().join("pairs.hic");
    run(
        &src,
        &dst,
        Format::Pairs,
        &Options {
            resolutions: Some(resolutions.to_vec()),
            ..Options::default()
        },
        Budget::new(1, 1 << 28),
    )
    .unwrap();
    let integer = Row {
        names: ("chr1".into(), "chr1".into()),
        chroms: (Some(0), Some(0)),
        pos: (1000, 2000),
        count: 1.0,
    };
    assert_eq!(
        std::fs::read(dst).unwrap(),
        self::expected(&[integer], &resolutions, "unknown").0
    );
}

#[test]
fn pairs_match_the_builder_in_memory_and_spilled() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(11);
    let rows = random_rows(&mut rng, 150_000, true);
    let resolutions = vec![5_000, 25_000, 100_000, 1_000_000];
    let (want, want_sidecar) = expected(&rows, &resolutions, "hg-test");
    let text = pairs_text(&rows, "hg-test");
    let known = rows
        .iter()
        .filter(|r| r.chroms.0.is_some() && r.chroms.1.is_some())
        .count() as u64;
    let options = Options {
        resolutions: Some(resolutions.clone()),
        sidecar: true,
        ..Options::default()
    };
    for storage in STORAGES {
        let src = dir.path().join("in.pairs");
        write_input(&src, &text, storage);
        let from = Format::detect(&src).unwrap();
        assert_eq!(
            from,
            if storage == Storage::Plain {
                Format::Pairs
            } else {
                Format::GzipPairs
            }
        );
        for budget in budgets() {
            let dst = dir.path().join("out.hic");
            let report = run(&src, &dst, from, &options, budget.clone()).unwrap();
            let got = std::fs::read(&dst).unwrap();
            assert!(
                got == want,
                "{storage:?} {budget:?}: {} vs {} bytes, first difference at {:?}",
                got.len(),
                want.len(),
                got.iter().zip(&want).position(|(a, b)| a != b)
            );
            assert_eq!(HicSidecar::load(&dst).unwrap(), want_sidecar);
            assert_eq!(report.records, known);
            if budget.memory <= 12 << 20 {
                assert!(report.spilled > 0, "a tight budget spills");
            }
            let [warning] = report.warnings.as_slice() else {
                panic!("one warning: {:?}", report.warnings);
            };
            let skipped = rows.len() as u64 - known;
            assert!(
                warning.starts_with(&format!("skipped {skipped} of {} rows", rows.len())),
                "{warning}"
            );
            for name in ["chrM (", "! (", "chrUn_1 ("] {
                assert!(warning.contains(name), "{warning}");
            }
        }
    }
}

#[test]
fn juicer_default_resolutions_and_chrom_sizes_override() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(5);
    let rows = random_rows(&mut rng, 20_000, false);
    let src = dir.path().join("in.pairs.gz");
    write_input(&src, &pairs_text(&rows, "hg-test"), Storage::Bgzf);
    let dst = dir.path().join("out.hic");
    let report = run(
        &src,
        &dst,
        Format::GzipPairs,
        &Options::default(),
        Budget::new(4, 1 << 30),
    )
    .unwrap();
    assert!(report.warnings.is_empty());
    assert_eq!(
        std::fs::read(&dst).unwrap(),
        expected(&rows, &JUICER_RESOLUTIONS, "hg-test").0
    );
    // Chromosome sizes given replace the header's; the genome stays.
    let options = Options {
        chrom_sizes: Some(sizes_file(dir.path())),
        resolutions: Some(vec![10_000, 50_000]),
        ..Options::default()
    };
    run(
        &src,
        &dst,
        Format::GzipPairs,
        &options,
        Budget::new(2, 1 << 30),
    )
    .unwrap();
    assert_eq!(
        std::fs::read(&dst).unwrap(),
        expected(&rows, &[10_000, 50_000], "hg-test").0
    );
}

#[test]
fn short_text_with_and_without_scores() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(3);
    let mut rows = random_rows(&mut rng, 60_000, true);
    let sizes = sizes_file(dir.path());
    let resolutions = vec![5_000, 20_000, 500_000];
    let options = Options {
        chrom_sizes: Some(sizes.clone()),
        resolutions: Some(resolutions.clone()),
        ..Options::default()
    };
    for scored in [false, true] {
        if scored {
            for r in &mut rows {
                r.count = rng.range(0, 6) as f32;
            }
        }
        let (want, _) = expected(&rows, &resolutions, "unknown");
        let text = short_text(&rows, scored);
        for storage in [Storage::Plain, Storage::Gzip] {
            let src = dir.path().join("in.txt");
            write_input(&src, &text, storage);
            let from = Format::detect(&src).unwrap();
            assert_eq!(
                from,
                if storage == Storage::Plain {
                    Format::JuicerShort
                } else {
                    Format::GzipJuicerShort
                }
            );
            for budget in budgets() {
                let dst = dir.path().join("out.hic");
                let report = run(&src, &dst, from, &options, budget.clone()).unwrap();
                assert!(
                    std::fs::read(&dst).unwrap() == want,
                    "scored {scored} {storage:?} {budget:?}"
                );
                assert_eq!(report.warnings.len(), 1);
            }
        }
    }
    // Short text has no header: chromosome sizes are required.
    let src = dir.path().join("in.txt");
    let err = run(
        &src,
        &dir.path().join("x.hic"),
        Format::JuicerShort,
        &Options::default(),
        Budget::new(1, 1 << 30),
    )
    .unwrap_err();
    assert!(err.to_string().contains("chromosome sizes"), "{err}");
}

#[test]
fn fractional_scores_do_not_depend_on_threads_or_spilling() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(8);
    let mut rows = random_rows(&mut rng, 80_000, false);
    for r in &mut rows {
        r.count = (rng.range(1, 1_000_000) as f32) / 997.0;
    }
    let src = dir.path().join("in.txt");
    write_input(&src, &short_text(&rows, true), Storage::Plain);
    let options = Options {
        chrom_sizes: Some(sizes_file(dir.path())),
        resolutions: Some(vec![5_000, 10_000, 250_000]),
        ..Options::default()
    };
    let mut first: Option<Vec<u8>> = None;
    for budget in budgets() {
        let dst = dir.path().join("out.hic");
        run(&src, &dst, Format::JuicerShort, &options, budget).unwrap();
        let got = std::fs::read(&dst).unwrap();
        match &first {
            None => first = Some(got),
            Some(f) => crate::support::assertions::bytes_eq(&got, f),
        }
    }
}

/// straw's raw observed records: (x, y) bin starts and counts.
fn straw(path: &Path, loc1: &str, loc2: &str, res: u32) -> Vec<(u32, u32, f32)> {
    let c = |s: &str| CString::new(s).unwrap();
    let (t, n, f, l1, l2, u) = (
        c("observed"),
        c("NONE"),
        c(path.to_str().unwrap()),
        c(loc1),
        c(loc2),
        c("BP"),
    );
    // SAFETY: valid C strings; the result is released after copying.
    unsafe {
        let r = straw_sys::straw_query(
            t.as_ptr(),
            n.as_ptr(),
            f.as_ptr(),
            l1.as_ptr(),
            l2.as_ptr(),
            u.as_ptr(),
            res as i32,
        );
        assert!(!r.is_null(), "straw failed");
        let len = (*r).n as usize;
        let out = (0..len)
            .map(|i| {
                (
                    *(*r).bin_x.add(i) as u32,
                    *(*r).bin_y.add(i) as u32,
                    *(*r).counts.add(i),
                )
            })
            .collect();
        straw_sys::straw_free(r);
        out
    }
}

#[test]
fn straw_reads_back_the_counts() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(21);
    let rows = random_rows(&mut rng, 40_000, false);
    let src = dir.path().join("in.pairs");
    write_input(&src, &pairs_text(&rows, "hg-test"), Storage::Plain);
    let dst = dir.path().join("out.hic");
    let options = Options {
        resolutions: Some(vec![10_000, 50_000]),
        ..Options::default()
    };
    run(&src, &dst, Format::Pairs, &options, Budget::new(3, 1 << 30)).unwrap();
    for res in [10_000u32, 50_000] {
        for (a, b) in [(0u32, 0u32), (0, 1), (1, 2)] {
            let mut want: BTreeMap<(u32, u32), f32> = BTreeMap::new();
            for r in &rows {
                let (Some(c1), Some(c2)) = r.chroms else {
                    continue;
                };
                let (mut p, mut q) = (r.pos.0, r.pos.1);
                let (c1, c2) = if c1 <= c2 {
                    (c1, c2)
                } else {
                    std::mem::swap(&mut p, &mut q);
                    (c2, c1)
                };
                if (c1, c2) != (a, b) {
                    continue;
                }
                let clamp = |c: u32, p: u32| p.min(CHROMS[c as usize].1 as u32);
                let (mut x, mut y) = (clamp(c1, p) / res, clamp(c2, q) / res);
                if c1 == c2 && x > y {
                    std::mem::swap(&mut x, &mut y);
                }
                *want.entry((x * res, y * res)).or_default() += r.count;
            }
            let (n1, l1) = CHROMS[a as usize];
            let (n2, l2) = CHROMS[b as usize];
            let mut got: Vec<(u32, u32, f32)> =
                straw(&dst, &format!("{n1}:0:{l1}"), &format!("{n2}:0:{l2}"), res);
            got.sort_by_key(|&(x, y, _)| (x, y));
            let want: Vec<(u32, u32, f32)> =
                want.into_iter().map(|((x, y), v)| (x, y, v)).collect();
            assert_eq!(got, want, "{n1} x {n2} at {res}");
        }
    }
}

#[test]
fn bad_resolutions_and_headers_are_errors() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(2);
    let rows = random_rows(&mut rng, 100, false);
    let src = dir.path().join("in.pairs");
    write_input(&src, &pairs_text(&rows, "g"), Storage::Plain);
    let dst = dir.path().join("out.hic");
    let options = Options {
        resolutions: Some(vec![5_000, 25_000, 12_000]),
        ..Options::default()
    };
    let err = run(&src, &dst, Format::Pairs, &options, Budget::new(1, 1 << 30)).unwrap_err();
    assert!(
        err.to_string()
            .contains("resolution 12000 is not a multiple of the finest, 5000"),
        "{err}"
    );
    // No #chromsize lines and no chromosome sizes.
    let text: Vec<u8> = pairs_text(&rows, "g")
        .split(|&c| c == b'\n')
        .filter(|l| !l.starts_with(b"#chromsize"))
        .flat_map(|l| l.iter().copied().chain([b'\n']))
        .collect();
    write_input(&src, &text, Storage::Plain);
    let err = run(
        &src,
        &dst,
        Format::Pairs,
        &Options::default(),
        Budget::new(1, 1 << 30),
    )
    .unwrap_err();
    assert!(err.to_string().contains("#chromsize"), "{err}");
    assert!(!dst.exists());
    // The conversion is chosen by format.
    let c = Conversion::infer(&src, &dst).unwrap();
    assert!(matches!(
        c,
        Conversion::PairsToHic(PairsToHic {
            dialect: PairsDialect::Pairs,
            ..
        })
    ));
}
