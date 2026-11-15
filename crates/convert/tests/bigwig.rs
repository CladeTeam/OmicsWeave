//! bedGraph and wig → bigWig: the same file and sidecar as `BigWigWriter`
//! given the items sorted by (chromosome in order of first item, start),
//! whether the sort stayed in memory or spilled; and the items, values and
//! zoom summaries libBigWig and the crate read back are those the text
//! describes.

#[path = "../../bbi/tests/common/mod.rs"]
mod bbifix;
mod common;

use std::fmt::Write as _;
use std::io::Cursor;
use std::path::Path;

use bbifix::Opened;
use bbifix::libbigwig::Reader;
use common::{Rng, STORAGES, Storage, budgets, write_input};
use omicsweave_bbi::{BbiSidecar, BigWigWriter, Chrom, WriteOptions};
use omicsweave_convert::{
    BedGraphToBigWig, Budget, ChromSizes, Conversion, Format, WigToBigWig, convert,
};
use omicsweave_core::Error;

/// One item the text describes: chromosome index, start, end, value.
type Item = (usize, u32, u32, f32);

/// A source text and what it describes.
struct Model {
    chroms: Vec<(String, u32)>,
    text: String,
    /// In text order.
    items: Vec<Item>,
}

fn chroms(rng: &mut Rng) -> Vec<(String, u32)> {
    (0..rng.range(1, 5))
        .map(|i| {
            (
                format!("chr{}", [5, 2, 11, 1][i as usize]),
                rng.range(50_000_000, 100_000_000) as u32,
            )
        })
        .collect()
}

/// A value and how it is written; the text parses back to the value.
fn value(rng: &mut Rng) -> (f32, String) {
    let v = match rng.range(0, 10) {
        0 => rng.range(0, 50) as f32,
        1 => 0.0,
        _ => (rng.range(0, 2_000_000) as f32 - 1_000_000.0) / 1000.0,
    };
    let text = if rng.chance(10) {
        format!("{v:e}")
    } else {
        format!("{v}")
    };
    assert_eq!(text.parse::<f64>().unwrap() as f32, v);
    (v, text)
}

fn noise(rng: &mut Rng, out: &mut String) {
    match rng.range(0, 40) {
        0 => out.push_str("# comment\n"),
        1 => out.push('\n'),
        2 => out.push_str("   \t\n"),
        _ => {}
    }
}

fn eol(rng: &mut Rng) -> &'static str {
    if rng.chance(10) { "\r\n" } else { "\n" }
}

fn sep(rng: &mut Rng) -> &'static str {
    [" ", "\t", "  "][rng.range(0, 3) as usize]
}

/// Shuffled bedGraph lines of disjoint items.
fn random_bedgraph(rng: &mut Rng, n: usize) -> Model {
    let chroms = chroms(rng);
    let mut lines = Vec::new();
    let mut pos = vec![0u32; chroms.len()];
    for _ in 0..n {
        let c = rng.range(0, chroms.len() as u64) as usize;
        let start = pos[c] + rng.range(0, 40) as u32;
        let end = start + rng.range(1, 60) as u32;
        if end > chroms[c].1 {
            continue;
        }
        pos[c] = end;
        let (v, t) = value(rng);
        lines.push(((c, start, end, v), t));
    }
    for i in (1..lines.len()).rev() {
        let j = rng.range(0, i as u64 + 1) as usize;
        lines.swap(i, j);
    }
    let mut text = String::from("browser position chr1:1-100\ntrack type=bedGraph name=t\n");
    let mut items = Vec::new();
    for ((c, s, e, v), t) in lines {
        noise(rng, &mut text);
        let (a, b, d) = (sep(rng), sep(rng), sep(rng));
        let _ = write!(text, "{}{a}{s}{b}{e}{d}{t}{}", chroms[c].0, eol(rng));
        items.push((c, s, e, v));
    }
    Model {
        chroms,
        text,
        items,
    }
}

/// Wig sections of every kind, disjoint, in shuffled order: `fixedStep`
/// (span at most the step), `variableStep` (lines sometimes shuffled) and
/// bedGraph lines. Sections are long enough to cross batch boundaries.
fn random_wig(rng: &mut Rng, sections: usize) -> Model {
    let chroms = chroms(rng);
    let mut pos = vec![0u32; chroms.len()];
    let mut blocks: Vec<(String, Vec<Item>)> = Vec::new();
    for _ in 0..sections {
        let c = rng.range(0, chroms.len() as u64) as usize;
        let name = &chroms[c].0;
        let n = if rng.chance(25) {
            rng.range(2_000, 10_000)
        } else {
            rng.range(1, 300)
        };
        let mut text = String::new();
        let mut items = Vec::new();
        let start = pos[c] + rng.range(0, 1000) as u32;
        let mut at = start;
        match rng.range(0, 3) {
            0 => {
                let span = rng.range(1, 30) as u32;
                let step = span + rng.range(0, 20) as u32;
                let _ = write!(text, "fixedStep chrom={name} start={}", start + 1);
                if step != 1 || rng.chance(50) {
                    let _ = write!(text, " step={step}");
                }
                if span != 1 || rng.chance(50) {
                    let _ = write!(text, " span={span}");
                }
                text.push_str(eol(rng));
                for _ in 0..n {
                    let (v, t) = value(rng);
                    noise(rng, &mut text);
                    let _ = write!(text, "{t}{}", eol(rng));
                    items.push((c, at, at + span, v));
                    at += step;
                }
            }
            1 => {
                let span = rng.range(1, 30) as u32;
                let _ = write!(text, "variableStep chrom={name}");
                if span != 1 || rng.chance(50) {
                    let _ = write!(text, "\tspan={span}");
                }
                text.push_str(eol(rng));
                let mut lines = Vec::new();
                for _ in 0..n {
                    let (v, t) = value(rng);
                    lines.push(((c, at, at + span, v), t));
                    at += span + rng.range(0, 40) as u32;
                }
                if rng.chance(30) {
                    for i in (1..lines.len()).rev() {
                        let j = rng.range(0, i as u64 + 1) as usize;
                        lines.swap(i, j);
                    }
                }
                for (it, t) in lines {
                    noise(rng, &mut text);
                    let _ = write!(text, "{}{}{t}{}", it.1 + 1, sep(rng), eol(rng));
                    items.push(it);
                }
            }
            _ => {
                for _ in 0..n {
                    let end = at + rng.range(1, 50) as u32;
                    let (v, t) = value(rng);
                    noise(rng, &mut text);
                    let _ = write!(text, "{name}\t{at}\t{end}\t{t}{}", eol(rng));
                    items.push((c, at, end, v));
                    at = end + rng.range(0, 30) as u32;
                }
            }
        }
        pos[c] = items.iter().map(|i| i.2).max().unwrap_or(at);
        blocks.push((text, items));
    }
    for i in (1..blocks.len()).rev() {
        let j = rng.range(0, i as u64 + 1) as usize;
        blocks.swap(i, j);
    }
    let mut text = String::from("track type=wiggle_0 name=test\n");
    let mut items = Vec::new();
    for (t, i) in blocks {
        text.push_str(&t);
        items.extend(i);
    }
    Model {
        chroms,
        text,
        items,
    }
}

/// Chromosomes in order of first item, with their lengths.
fn order(m: &Model, sizes: bool) -> Vec<(usize, Chrom)> {
    let mut out: Vec<(usize, Chrom)> = Vec::new();
    for &(c, _, e, _) in &m.items {
        match out.iter_mut().find(|(o, _)| *o == c) {
            Some((_, chrom)) => chrom.len = chrom.len.max(e),
            None => out.push((
                c,
                Chrom {
                    name: m.chroms[c].0.clone(),
                    len: e,
                },
            )),
        }
    }
    if sizes {
        for (c, chrom) in &mut out {
            chrom.len = m.chroms[*c].1;
        }
    }
    out
}

/// (start, end, value) of each output chromosome's items, by id.
type ByChrom = Vec<Vec<(u32, u32, f32)>>;

/// Items by output chromosome id, sorted by start.
fn sorted(m: &Model, sizes: bool) -> (Vec<Chrom>, ByChrom) {
    let ord = order(m, sizes);
    let mut by = vec![Vec::new(); ord.len()];
    for &(c, s, e, v) in &m.items {
        let id = ord.iter().position(|(o, _)| *o == c).unwrap();
        by[id].push((s, e, v));
    }
    for items in &mut by {
        items.sort_by_key(|i| i.0);
    }
    (ord.into_iter().map(|(_, c)| c).collect(), by)
}

fn expected(m: &Model, sizes: bool) -> (Vec<u8>, BbiSidecar) {
    let (chroms, by) = sorted(m, sizes);
    let mut w =
        BigWigWriter::create(Cursor::new(Vec::new()), chroms, WriteOptions::default()).unwrap();
    for (id, items) in by.iter().enumerate() {
        for &(s, e, v) in items {
            w.add(id as u32, s, e, v).unwrap();
        }
    }
    let (out, sidecar) = w.finish().unwrap();
    (out.into_inner(), sidecar)
}

fn conversion(wig: bool, src: &Path, dst: &Path, sizes: ChromSizes) -> Conversion {
    if wig {
        Conversion::WigToBigWig(WigToBigWig::new(src, dst, sizes))
    } else {
        Conversion::BedGraphToBigWig(BedGraphToBigWig::new(src, dst, sizes))
    }
}

/// What libBigWig and the crate read from `path` is what `m` describes:
/// chromosomes, items, per-base values, the total summary and every zoom
/// level's summaries.
fn check_read_back(path: &Path, m: &Model, sizes: bool, rng: &mut Rng) {
    let (chroms, by) = sorted(m, sizes);
    let mut theirs = Reader::open(path).expect("libBigWig opens the file");
    let mut ours = Opened::open(path);
    let names: Vec<(String, u32)> = chroms.iter().map(|c| (c.name.clone(), c.len)).collect();
    assert_eq!(theirs.chroms(), names);
    let h = theirs.header();
    assert_eq!(h.kind, 0);
    let covered: u64 = by.iter().flatten().map(|&(s, e, _)| (e - s) as u64).sum();
    let all = || by.iter().flatten().map(|i| i.2);
    assert_eq!(h.bases_covered, covered);
    assert_eq!(h.min, all().fold(f32::INFINITY, f32::min) as f64);
    assert_eq!(h.max, all().fold(f32::NEG_INFINITY, f32::max) as f64);
    let sum: f64 = by
        .iter()
        .flatten()
        .map(|&(s, e, v)| v as f64 * (e - s) as f64)
        .sum();
    assert!(
        (h.sum - sum).abs() <= 1e-6 * sum.abs().max(1.0),
        "{} {sum}",
        h.sum
    );
    assert!(h.zoom_levels > 0);
    assert_eq!(
        theirs.zooms(),
        ours.zooms
            .iter()
            .map(|z| (z.reduction, z.data_offset, z.index_offset))
            .collect::<Vec<_>>()
    );

    for (id, (chrom, items)) in chroms.iter().zip(&by).enumerate() {
        assert_eq!(&theirs.intervals(&chrom.name, 0, chrom.len).unwrap(), items);
        assert_eq!(&ours.intervals(&chrom.name, 0, chrom.len), items);
        for _ in 0..20 {
            let start = rng.range(0, chrom.len as u64) as u32;
            let end = (start + rng.range(1, 20_000) as u32).min(chrom.len);
            let mut want = vec![f32::NAN; (end - start) as usize];
            for &(s, e, v) in items {
                for p in s.max(start)..e.min(end) {
                    want[(p - start) as usize] = v;
                }
            }
            let got = theirs.values(&chrom.name, start, end).unwrap();
            assert!(
                got.iter().zip(&want).all(|(a, b)| bbifix::same(*a, *b)),
                "libBigWig values {}:{start}-{end}",
                chrom.name
            );
            let got = ours.values(&chrom.name, start, end);
            assert!(
                got.iter().zip(&want).all(|(a, b)| bbifix::same(*a, *b)),
                "values {}:{start}-{end}",
                chrom.name
            );
        }
        for zoom in 0..ours.zooms.len() {
            let mut valid = 0u64;
            for r in ours.summaries(zoom, id as u32, 0, chrom.len) {
                let first = items.partition_point(|&(_, e, _)| e <= r.start);
                let (mut n, mut min, mut max) = (0u32, f32::INFINITY, f32::NEG_INFINITY);
                for &(s, e, v) in items[first..].iter().take_while(|i| i.0 < r.end) {
                    n += e.min(r.end) - s.max(r.start);
                    min = min.min(v);
                    max = max.max(v);
                }
                assert_eq!(r.valid, n, "zoom {zoom} covered bases of {r:?}");
                if n > 0 {
                    assert_eq!((r.min, r.max), (min, max), "zoom {zoom} range of {r:?}");
                }
                valid += r.valid as u64;
            }
            let here: u64 = items.iter().map(|&(s, e, _)| (e - s) as u64).sum();
            assert_eq!(valid, here, "zoom {zoom} covered bases of {}", chrom.name);
        }
    }
}

fn write_sizes(dir: &Path, m: &Model) -> std::path::PathBuf {
    let path = dir.join("chrom.sizes");
    let text: String = m
        .chroms
        .iter()
        .map(|(n, l)| format!("{n}\t{l}\n"))
        .collect();
    std::fs::write(&path, text).unwrap();
    path
}

/// Every storage and budget gives the writer's bytes; a tight budget
/// spills.
fn check_conversions(m: &Model, wig: bool, seed: u64) {
    let dir = tempfile::tempdir().unwrap();
    let sizes_path = write_sizes(dir.path(), m);
    let mut rng = Rng::new(seed);
    for with_sizes in [false, true] {
        let (want, want_sidecar) = expected(m, with_sizes);
        for storage in STORAGES {
            let src = dir.path().join("in.txt");
            write_input(&src, m.text.as_bytes(), storage);
            let mut all = budgets();
            all.push(Budget::new(4, 8 << 20));
            if with_sizes && storage != Storage::Plain {
                continue;
            }
            for budget in all {
                let dst = dir.path().join("out.bw");
                let sizes = if with_sizes {
                    ChromSizes::File(sizes_path.clone())
                } else {
                    ChromSizes::FromData
                };
                let report =
                    convert(conversion(wig, &src, &dst, sizes), budget.clone(), None).unwrap();
                assert_eq!(report.records, m.items.len() as u64);
                assert_eq!(report.files[0], (Format::BigWig, dst.clone()));
                assert!(
                    std::fs::read(&dst).unwrap() == want,
                    "{storage:?} {budget:?} sizes {with_sizes}"
                );
                assert_eq!(BbiSidecar::load(&dst).unwrap(), want_sidecar);
                if budget.memory <= 8 << 20 {
                    assert!(report.spilled > 0, "a tight budget spills");
                }
            }
        }
        let dst = dir.path().join("out.bw");
        check_read_back(&dst, m, with_sizes, &mut rng);
    }
}

#[test]
fn bedgraph_matches_the_writer_and_reads_back() {
    for seed in 0..2u64 {
        let mut rng = Rng::new(seed);
        let m = random_bedgraph(&mut rng, 200_000);
        check_conversions(&m, false, seed);
    }
}

#[test]
fn wig_matches_the_writer_and_reads_back() {
    for seed in 0..2u64 {
        let mut rng = Rng::new(100 + seed);
        let m = random_wig(&mut rng, 160);
        assert!(m.text.len() > 1 << 20, "sections cross batches");
        check_conversions(&m, true, seed);
    }
}

#[test]
fn bedgraph_lines_in_wig_end_a_step_section() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("a.wig");
    std::fs::write(
        &src,
        "fixedStep chrom=c1 start=101 step=10 span=5\n1\n2\nc2 0 4 7\nc1 300 310 -1\n\
         variableStep chrom=c2 span=3\n# c\n11 0.5\n6 0.25\n",
    )
    .unwrap();
    let dst = dir.path().join("a.bw");
    convert(
        Conversion::infer(&src, &dst).unwrap(),
        Budget::new(2, 64 << 20),
        None,
    )
    .unwrap();
    let mut r = Reader::open(&dst).unwrap();
    assert_eq!(
        r.chroms(),
        [("c1".to_string(), 310), ("c2".to_string(), 13)]
    );
    assert_eq!(
        r.intervals("c1", 0, 310).unwrap(),
        [(100, 105, 1.0), (110, 115, 2.0), (300, 310, -1.0)]
    );
    assert_eq!(
        r.intervals("c2", 0, 13).unwrap(),
        [(0, 4, 7.0), (5, 8, 0.25), (10, 13, 0.5)]
    );
}

fn line_of(e: Error) -> u64 {
    match e {
        Error::Format { line, .. } => line,
        e => panic!("not a format error: {e}"),
    }
}

#[test]
fn overlaps_and_out_of_range_items_are_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let sizes = dir.path().join("sizes");
    std::fs::write(&sizes, "chr1\t1000\nchr2\t50\n").unwrap();
    let dst = dir.path().join("x.bw");
    for (wig, text, line) in [
        (false, "chr1\t0\t10\t1\nchr2\t0\t5\t1\nchr1\t5\t20\t2\n", 3),
        (false, "chr1\t500\t600\t1\nchr1\t10\t520\t2\n", 1),
        (false, "chr2\t40\t60\t1\n", 1),
        (
            true,
            "fixedStep chrom=chr1 start=1 step=2 span=3\n1\n2\n",
            3,
        ),
        (true, "variableStep chrom=chr2 span=5\n48 1\n", 2),
    ] {
        let src = dir.path().join("x.txt");
        std::fs::write(&src, text).unwrap();
        let e = convert(
            conversion(wig, &src, &dst, ChromSizes::File(sizes.clone())),
            Budget::new(2, 64 << 20),
            None,
        )
        .unwrap_err();
        assert_eq!(line_of(e), line, "{text:?}");
        assert!(!dst.exists());
    }
    let src = dir.path().join("x.txt");
    std::fs::write(&src, "chrX\t0\t5\t1\n").unwrap();
    let e = convert(
        conversion(false, &src, &dst, ChromSizes::File(sizes)),
        Budget::new(2, 64 << 20),
        None,
    )
    .unwrap_err();
    assert!(e.to_string().contains("chrX"), "{e}");
}

#[test]
fn signal_texts_are_told_apart_by_content() {
    let dir = tempfile::tempdir().unwrap();
    let cases: [(&str, Format, Format); 6] = [
        (
            "chr1\t0\t10\t1.5\nchr1 10 20 -2\n",
            Format::BedGraph,
            Format::GzipBedGraph,
        ),
        (
            "track type=bedGraph\nchr1\t0\t10\tx\n",
            Format::BedGraph,
            Format::GzipBedGraph,
        ),
        (
            "# x\nfixedStep chrom=chr1 start=1 step=1\n1\n",
            Format::Wig,
            Format::GzipWig,
        ),
        (
            "track type=wiggle_0\nchr1 0 1 2\n",
            Format::Wig,
            Format::GzipWig,
        ),
        ("chr1\t0\t10\tpeak1\n", Format::Bed, Format::GzipBed),
        ("chr1\t0\t10\t1\t0\t+\n", Format::Bed, Format::GzipBed),
    ];
    for (text, plain, gzip) in cases {
        for storage in STORAGES {
            let path = dir.path().join("f");
            write_input(&path, text.as_bytes(), storage);
            let want = if storage == Storage::Plain {
                plain
            } else {
                gzip
            };
            assert_eq!(Format::detect(&path).unwrap(), want, "{text:?} {storage:?}");
        }
    }
    for (name, format) in [
        ("a.bedGraph", Format::BedGraph),
        ("a.bdg.gz", Format::GzipBedGraph),
        ("a.wig", Format::Wig),
        ("a.wig.gz", Format::GzipWig),
    ] {
        assert_eq!(Format::from_extension(Path::new(name)).unwrap(), format);
    }
    let routes: Vec<_> = Conversion::routes().collect();
    for route in [
        (Format::BedGraph, Format::BigWig),
        (Format::GzipWig, Format::BigWig),
        (Format::BedGraph, Format::BigBed),
    ] {
        assert!(routes.contains(&route), "{route:?}");
    }
}
