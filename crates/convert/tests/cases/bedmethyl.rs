//! bedMethyl → bigBed and bigWig: modkit (tab- and space-separated counts)
//! and ENCODE dialects; every field back from the bigBed, the selected
//! column back from the bigWig; selection errors; the same files whatever
//! the thread count or spilling.

use crate::bbi_support as bbifix;
use crate::support as common;

use std::fmt::Write as _;
use std::path::Path;

use bbifix::Opened;
use common::budgets;
use common::rng::Rng;
use common::storage::{STORAGES, Storage, write_input};
use omicsweave_convert::{Budget, Conversion, Format, Options, convert};

/// One item: chromosome index, start, and its fields after the third.
#[derive(Clone)]
struct Item {
    chrom: usize,
    start: u32,
    fields: Vec<String>,
}

impl Item {
    fn code(&self) -> &str {
        &self.fields[0]
    }
    fn strand(&self) -> &str {
        &self.fields[2]
    }
}

const CHROMS: [(&str, u32); 2] = [("chr1", 400_000), ("chr2", 250_000)];

/// modkit items: at each sampled C, `m` and `h` calls on `+`, and on the
/// `-` strand one base on, in shuffled order.
fn modkit_items(rng: &mut Rng, sites: usize) -> Vec<Item> {
    let mut items = Vec::new();
    for _ in 0..sites {
        let chrom = rng.range(0, 2) as usize;
        // Even positions: a site's `-` call never meets another's `+` call.
        let pos = rng.range(0, CHROMS[chrom].1 as u64 / 2 - 1) as u32 * 2;
        for (strand, at) in [("+", pos), ("-", pos + 1)] {
            for code in ["m", "h"] {
                let valid = rng.range(1, 60);
                let nmod = rng.range(0, valid + 1);
                let canonical = valid - nmod;
                let percent = format!("{:.2}", nmod as f64 / valid as f64 * 100.0);
                let mut f = vec![
                    code.to_owned(),
                    valid.to_string(),
                    strand.to_owned(),
                    at.to_string(),
                    (at + 1).to_string(),
                    "255,0,0".to_owned(),
                    valid.to_string(),
                    percent,
                    nmod.to_string(),
                    canonical.to_string(),
                ];
                for _ in 0..5 {
                    f.push(rng.range(0, 4).to_string());
                }
                items.push(Item {
                    chrom,
                    start: at,
                    fields: f,
                });
            }
        }
    }
    dedup(shuffle(rng, items))
}

/// ENCODE items, one per strand and position.
fn encode_items(rng: &mut Rng, sites: usize) -> Vec<Item> {
    let mut items = Vec::new();
    for _ in 0..sites {
        let chrom = rng.range(0, 2) as usize;
        let at = rng.range(0, CHROMS[chrom].1 as u64 - 1) as u32;
        let coverage = rng.range(1, 2_000);
        items.push(Item {
            chrom,
            start: at,
            fields: vec![
                "Sample".into(),
                coverage.min(1000).to_string(),
                if rng.chance(50) { "+" } else { "-" }.into(),
                at.to_string(),
                (at + 1).to_string(),
                "0,255,0".into(),
                coverage.to_string(),
                rng.range(0, 101).to_string(),
            ],
        });
    }
    dedup(shuffle(rng, items))
}

fn shuffle(rng: &mut Rng, mut items: Vec<Item>) -> Vec<Item> {
    for i in (1..items.len()).rev() {
        let j = rng.range(0, i as u64 + 1) as usize;
        items.swap(i, j);
    }
    items
}

/// One item per (chromosome, start, strand, code).
fn dedup(items: Vec<Item>) -> Vec<Item> {
    let mut seen = std::collections::HashSet::new();
    items
        .into_iter()
        .filter(|i| seen.insert((i.chrom, i.start, i.strand().to_owned(), i.code().to_owned())))
        .collect()
}

/// The text of `items`; `spaces` writes modkit's counts space-separated.
fn text(items: &[Item], spaces: bool) -> Vec<u8> {
    let mut t = String::from("track name=calls\n");
    for i in items {
        write!(t, "{}\t{}\t{}", CHROMS[i.chrom].0, i.start, i.start + 1).unwrap();
        for (k, f) in i.fields.iter().enumerate() {
            let sep = if spaces && k >= 7 { ' ' } else { '\t' };
            write!(t, "{sep}{f}").unwrap();
        }
        t.push('\n');
    }
    t.into_bytes()
}

/// Items in bigBed order: chromosomes by first appearance, then start,
/// ties in input order.
fn sorted(items: &[Item]) -> (Vec<usize>, Vec<Vec<&Item>>) {
    let mut order = Vec::new();
    for i in items {
        if !order.contains(&i.chrom) {
            order.push(i.chrom);
        }
    }
    let by = order
        .iter()
        .map(|&c| {
            let mut v: Vec<&Item> = items.iter().filter(|i| i.chrom == c).collect();
            v.sort_by_key(|i| i.start);
            v
        })
        .collect();
    (order, by)
}

fn run(
    src: &Path,
    dst: &Path,
    to: Format,
    options: &Options,
    budget: Budget,
) -> omicsweave_core::Result<omicsweave_convert::Report> {
    let mut c = Conversion::between(Format::detect(src)?, to, src, dst)?;
    c.configure(options)?;
    convert(c, budget, None)
}

fn check_bigbed(path: &Path, items: &[Item], auto_sql_has: &str) {
    let mut bb = Opened::open(path);
    assert!(bb.auto_sql.as_deref().unwrap().contains(auto_sql_has));
    let (order, by) = sorted(items);
    for (&c, want) in order.iter().zip(&by) {
        let (name, len) = CHROMS[c];
        let got = bb.entries(name, 0, len);
        assert_eq!(got.len(), want.len(), "{name}");
        for (g, w) in got.iter().zip(want) {
            assert_eq!((g.0, g.1), (w.start, w.start + 1));
            assert_eq!(String::from_utf8_lossy(&g.2), w.fields.join("\t"));
        }
    }
}

#[test]
fn modkit_to_bigbed_keeps_every_field() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(1);
    let items = modkit_items(&mut rng, 15_000);
    let mut first: Option<Vec<u8>> = None;
    for spaces in [false, true] {
        for storage in STORAGES {
            let src = dir.path().join("calls.bed");
            write_input(&src, &text(&items, spaces), storage);
            assert_eq!(
                Format::detect(&src).unwrap(),
                if storage == Storage::Plain {
                    Format::BedMethyl
                } else {
                    Format::GzipBedMethyl
                }
            );
            let mut all = budgets();
            all.push(Budget::new(4, 8 << 20));
            for budget in all {
                let dst = dir.path().join("calls.bb");
                let report = run(
                    &src,
                    &dst,
                    Format::BigBed,
                    &Options::default(),
                    budget.clone(),
                )
                .unwrap();
                assert_eq!(report.records, items.len() as u64);
                if budget.memory <= 8 << 20 {
                    assert!(report.spilled > 0, "a tight budget spills");
                }
                let got = std::fs::read(&dst).unwrap();
                match &first {
                    None => {
                        check_bigbed(&dst, &items, "Nvalid_cov");
                        first = Some(got);
                    }
                    Some(f) => assert!(&got == f, "{spaces} {storage:?} {budget:?}"),
                }
            }
        }
    }
}

#[test]
fn encode_to_bigbed_and_bigwig() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(2);
    let items = encode_items(&mut rng, 20_000);
    let src = dir.path().join("wgbs.bed.gz");
    write_input(&src, &text(&items, false), Storage::Bgzf);
    let bb = dir.path().join("wgbs.bb");
    run(
        &src,
        &bb,
        Format::BigBed,
        &Options::default(),
        Budget::new(3, 1 << 30),
    )
    .unwrap();
    check_bigbed(&bb, &items, "percentMeth");

    // Both strands at one position overlap until a strand is chosen.
    let plus: Vec<Item> = dedup_positions(items.iter().filter(|i| i.strand() == "+"));
    let mut both = plus.clone();
    let mut other = plus[0].clone();
    other.fields[2] = "-".into();
    both.push(other);
    let src = dir.path().join("both.bed");
    write_input(&src, &text(&both, false), Storage::Plain);
    let bw = dir.path().join("wgbs.bw");
    let err = run(
        &src,
        &bw,
        Format::BigWig,
        &Options::default(),
        Budget::new(2, 1 << 30),
    )
    .unwrap_err();
    assert!(err.to_string().contains("choose a strand"), "{err}");
    for (value, field) in [("coverage", 6), ("percent", 7), ("score", 1)] {
        let options = Options {
            value: Some(value.into()),
            strand: Some("+".into()),
            ..Options::default()
        };
        run(&src, &bw, Format::BigWig, &options, Budget::new(2, 1 << 30)).unwrap();
        check_bigwig(&bw, &plus, field);
    }
    let options = Options {
        value: Some("nmod".into()),
        ..Options::default()
    };
    let err = run(&src, &bw, Format::BigWig, &options, Budget::new(2, 1 << 30)).unwrap_err();
    assert!(err.to_string().contains("modkit column"), "{err}");
}

/// One item per (chromosome, start), the first.
fn dedup_positions<'a>(items: impl Iterator<Item = &'a Item>) -> Vec<Item> {
    let mut seen = std::collections::HashSet::new();
    items
        .filter(|i| seen.insert((i.chrom, i.start)))
        .cloned()
        .collect()
}

/// The bigWig holds `items` with the value of their `field` (after the
/// third).
fn check_bigwig(path: &Path, items: &[Item], field: usize) {
    let mut bw = Opened::open(path);
    let (order, by) = sorted(items);
    for (&c, want) in order.iter().zip(&by) {
        let (name, len) = CHROMS[c];
        let want: Vec<(u32, u32, f32)> = want
            .iter()
            .map(|i| {
                let v: f64 = i.fields[field].parse().unwrap();
                (i.start, i.start + 1, v as f32)
            })
            .collect();
        assert_eq!(bw.intervals(name, 0, len), want, "{name}");
    }
}

#[test]
fn modkit_to_bigwig_selects_a_code_and_a_strand() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(3);
    let items = modkit_items(&mut rng, 12_000);
    let src = dir.path().join("calls.bedmethyl");
    write_input(&src, &text(&items, true), Storage::Plain);
    let bw = dir.path().join("calls.bw");

    let err = run(
        &src,
        &bw,
        Format::BigWig,
        &Options::default(),
        Budget::new(2, 1 << 30),
    )
    .unwrap_err();
    assert!(
        err.to_string().contains("2 modification codes (h, m)"),
        "{err}"
    );
    let options = Options {
        mod_code: Some("a".into()),
        ..Options::default()
    };
    let err = run(&src, &bw, Format::BigWig, &options, Budget::new(2, 1 << 30)).unwrap_err();
    assert!(
        err.to_string().contains("no item has modification code"),
        "{err}"
    );

    for (code, strand, value, field) in [
        ("m", None, "percent", 7),
        ("h", Some("+"), "valid_coverage", 6),
        ("m", Some("-"), "nmod", 8),
        ("h", None, "nnocall", 14),
    ] {
        let keep: Vec<Item> = items
            .iter()
            .filter(|i| i.code() == code && strand.is_none_or(|s| i.strand() == s))
            .cloned()
            .collect();
        let options = Options {
            value: Some(value.into()),
            mod_code: Some(code.into()),
            strand: strand.map(str::to_owned),
            ..Options::default()
        };
        let mut first: Option<Vec<u8>> = None;
        let mut all = budgets();
        all.push(Budget::new(4, 8 << 20));
        for budget in all {
            let report = run(&src, &bw, Format::BigWig, &options, budget.clone()).unwrap();
            assert_eq!(report.records, keep.len() as u64);
            let got = std::fs::read(&bw).unwrap();
            match &first {
                None => {
                    check_bigwig(&bw, &keep, field);
                    first = Some(got);
                }
                Some(f) => assert!(&got == f, "{code} {strand:?} {budget:?}"),
            }
        }
    }
}

#[test]
fn bad_lines_and_options_are_errors() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("x.bedmethyl");
    let good = "chr1\t10\t11\tm\t5\t+\t10\t11\t255,0,0\t5\t60.00\t3\t2\t0\t0\t1\t0\t0\n";
    for (bad, line) in [
        (
            "chr1\t10\t11\tm\t5\t+\t10\t11\t255,0,0\t5\tx\t3\t2\t0\t0\t1\t0\t0\n",
            2,
        ),
        (
            "chr1\t10\t11\tm\t5\t*\t10\t11\t255,0,0\t5\t60\t3\t2\t0\t0\t1\t0\t0\n",
            2,
        ),
        ("chr1\t10\t11\t.\t7\t-\t10\t11\t0,255,0\t7\t43\n", 2),
    ] {
        std::fs::write(&src, format!("{good}{bad}")).unwrap();
        let c = Conversion::between(
            Format::BedMethyl,
            Format::BigBed,
            &src,
            dir.path().join("x.bb"),
        )
        .unwrap();
        let err = convert(c, Budget::new(1, 1 << 30), None).unwrap_err();
        let omicsweave_core::Error::Format { line: got, .. } = err else {
            panic!("{err}");
        };
        assert_eq!(got, line, "{bad:?}");
    }
    let mut c = Conversion::between(
        Format::BedMethyl,
        Format::BigWig,
        &src,
        dir.path().join("x.bw"),
    )
    .unwrap();
    for options in [
        Options {
            value: Some("fraction".into()),
            ..Options::default()
        },
        Options {
            strand: Some("both".into()),
            ..Options::default()
        },
    ] {
        assert!(c.configure(&options).is_err());
    }
    let mut c = Conversion::between(
        Format::BedMethyl,
        Format::BigBed,
        &src,
        dir.path().join("x.bb"),
    )
    .unwrap();
    let options = Options {
        mod_code: Some("m".into()),
        ..Options::default()
    };
    assert!(c.configure(&options).is_err(), "mod is for bigWig only");
}
