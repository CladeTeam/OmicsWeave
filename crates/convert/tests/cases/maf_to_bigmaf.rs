//! MAF → bigMaf: blocks already in kent's written form come back from the
//! bigBed unchanged, `;`-joined; chromosome lengths are the reference's
//! `srcSize`; the same file whatever the storage, thread count or
//! spilling; the reference is told or inferred; overlaps, a misplaced
//! reference and contradicting lengths are errors.

use crate::bbi_support as bbifix;
use crate::support as common;

use std::fmt::Write as _;
use std::path::Path;

use bbifix::Opened;
use common::budgets;
use common::rng::Rng;
use common::storage::{STORAGES, Storage, write_input};
use omicsweave_convert::{BIGMAF_AUTOSQL, Budget, Conversion, Format, Options, convert};

const CHROMS: [(&str, u32); 2] = [("chr1", 3_000_000), ("chr2.alt", 1_500_000)];
const OTHERS: [(&str, u32); 3] = [("panTro4.chr1", 9_000_000), ("mm10.chrX", 70_000), ("x", 5)];

/// One block: reference chromosome, start, size, and its lines as kent
/// writes them.
struct Block {
    chrom: usize,
    start: u32,
    size: u32,
    lines: Vec<String>,
}

fn digits(x: u32) -> usize {
    x.to_string().len()
}

/// A block of `size` reference bases at `start`, with a few other rows;
/// columns padded as `mafWriteDelimiter` pads them.
fn block(rng: &mut Rng, chrom: usize, start: u32, size: u32) -> Block {
    let mut text = String::new();
    for _ in 0..size {
        text.push(*rng.pick(&['A', 'C', 'g', 't']));
        if rng.chance(10) {
            text.push('-');
        }
    }
    let cols = text.len();
    let (name, len) = CHROMS[chrom];
    // (src, start, size, strand, srcSize, text, i line)
    let mut rows = vec![(
        format!("hg38.{name}"),
        start,
        size,
        '+',
        len,
        Some(text),
        None,
    )];
    for &(src, src_len) in &OTHERS {
        if !rng.chance(60) {
            continue;
        }
        if rng.chance(20) {
            let at = rng.range(0, u64::from(src_len)) as u32;
            rows.push((
                src.into(),
                at,
                rng.range(0, 40) as u32,
                '-',
                src_len,
                None,
                None,
            ));
            continue;
        }
        let row: String = (0..cols)
            .map(|_| *rng.pick(&['A', 'C', 'G', 'T', '-', 'n']))
            .collect();
        let bases = row.bytes().filter(|&c| c != b'-').count() as u32;
        if bases > src_len {
            continue;
        }
        let at = rng.range(0, u64::from(src_len - bases) + 1) as u32;
        // kent writes a row of no bases that has an `i` line as an `e` line.
        let info = (bases > 0 && rng.chance(50)).then(|| format!("C 0 I {}", rng.range(0, 99)));
        rows.push((src.into(), at, bases, '+', src_len, Some(row), info));
    }
    let w_src = rows.iter().map(|r| r.0.len()).max().unwrap();
    let w_start = rows.iter().map(|r| digits(r.1)).max().unwrap();
    let w_size = rows
        .iter()
        .map(|r| digits(if r.5.is_some() { r.2 } else { 0 }))
        .max()
        .unwrap();
    let w_len = rows.iter().map(|r| digits(r.4)).max().unwrap();
    let mut lines = vec![format!(
        "a score={}.{:06}",
        rng.range(0, 99_999),
        rng.range(0, 999_999)
    )];
    for (src, at, n, strand, src_len, text, info) in rows {
        match text {
            Some(t) => {
                lines.push(format!(
                    "s {src:<w_src$} {at:>w_start$} {n:>w_size$} {strand} {src_len:>w_len$} {t}"
                ));
                if let Some(i) = info {
                    lines.push(format!("i {src:<w_src$} {i}"));
                }
            }
            None => lines.push(format!(
                "e {src:<w_src$} {at:>w_start$} {n:>w_size$} {strand} {src_len:>w_len$} I"
            )),
        }
    }
    Block {
        chrom,
        start,
        size,
        lines,
    }
}

/// Non-overlapping blocks on both chromosomes, in shuffled order.
fn blocks(rng: &mut Rng, n: usize) -> Vec<Block> {
    let mut out = Vec::new();
    for (chrom, &(_, len)) in CHROMS.iter().enumerate() {
        let mut at = rng.range(0, 1_000) as u32;
        for _ in 0..n / 2 {
            let size = rng.range(1, 120) as u32;
            if at + size > len {
                break;
            }
            out.push(block(rng, chrom, at, size));
            at += size + rng.range(0, 30) as u32;
        }
    }
    for i in (1..out.len()).rev() {
        let j = rng.range(0, i as u64 + 1) as usize;
        out.swap(i, j);
    }
    out
}

fn text(blocks: &[Block]) -> Vec<u8> {
    let mut t = String::from("##maf version=1 scoring=roast.v3.3\n# multiz\n\n");
    for b in blocks {
        for l in &b.lines {
            writeln!(t, "{l}").unwrap();
        }
        t.push('\n');
    }
    t.into_bytes()
}

fn run(
    src: &Path,
    dst: &Path,
    options: &Options,
    budget: Budget,
) -> omicsweave_core::Result<omicsweave_convert::Report> {
    let mut c = Conversion::between(Format::detect(src)?, Format::BigBed, src, dst)?;
    c.configure(options)?;
    convert(c, budget, None)
}

fn check(path: &Path, blocks: &[Block]) {
    let mut bb = Opened::open(path);
    assert_eq!(bb.auto_sql.as_deref(), Some(BIGMAF_AUTOSQL));
    for (c, &(name, len)) in CHROMS.iter().enumerate() {
        let id = bb.chroms.id(name).unwrap();
        assert_eq!(bb.chroms.get(id).unwrap().len, len, "{name}: srcSize");
        let mut want: Vec<&Block> = blocks.iter().filter(|b| b.chrom == c).collect();
        want.sort_by_key(|b| b.start);
        let got = bb.entries(name, 0, len);
        assert_eq!(got.len(), want.len(), "{name}");
        for (g, w) in got.iter().zip(want) {
            assert_eq!((g.0, g.1), (w.start, w.start + w.size));
            let mut block = w.lines.join(";");
            block.push(';');
            assert_eq!(String::from_utf8_lossy(&g.2), block);
        }
    }
}

#[test]
fn blocks_round_trip_in_every_storage_and_budget() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(7);
    let blocks = blocks(&mut rng, 30_000);
    let mut first: Option<Vec<u8>> = None;
    for storage in STORAGES {
        let src = dir.path().join("multiz.maf");
        write_input(&src, &text(&blocks), storage);
        let want = if storage == Storage::Plain {
            Format::Maf
        } else {
            Format::GzipMaf
        };
        assert_eq!(Format::detect(&src).unwrap(), want);
        let mut all = budgets();
        all.push(Budget::new(4, 8 << 20));
        for budget in all {
            let dst = dir.path().join("multiz.bigmaf");
            let report = run(&src, &dst, &Options::default(), budget.clone()).unwrap();
            assert_eq!(report.records, blocks.len() as u64);
            if budget.memory <= 8 << 20 {
                assert!(report.spilled > 0, "a tight budget spills");
            }
            let got = std::fs::read(&dst).unwrap();
            match &first {
                None => {
                    check(&dst, &blocks);
                    first = Some(got);
                }
                Some(f) => assert!(&got == f, "{storage:?} {budget:?}"),
            }
        }
    }
}

#[test]
fn the_reference_is_told_or_inferred() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(8);
    let blocks = blocks(&mut rng, 200);
    let src = dir.path().join("x.maf.gz");
    write_input(&src, &text(&blocks), Storage::Gzip);
    let dst = dir.path().join("x.bb");
    let told = Options {
        ref_species: Some("hg38".into()),
        ..Options::default()
    };
    run(&src, &dst, &told, Budget::new(2, 1 << 30)).unwrap();
    check(&dst, &blocks);

    let wrong = Options {
        ref_species: Some("mm10".into()),
        ..Options::default()
    };
    let err = run(&src, &dst, &wrong, Budget::new(2, 1 << 30)).unwrap_err();
    assert!(
        err.to_string().contains("not with the reference mm10"),
        "{err}"
    );

    let sizes = dir.path().join("sizes");
    std::fs::write(&sizes, "chr1\t3000000\nchr2.alt\t1500000\nchr3\t9\n").unwrap();
    let with_sizes = Options {
        chrom_sizes: Some(sizes.clone()),
        ..Options::default()
    };
    run(&src, &dst, &with_sizes, Budget::new(2, 1 << 30)).unwrap();
    check(&dst, &blocks);
    std::fs::write(&sizes, "chr1\t3000001\nchr2.alt\t1500000\n").unwrap();
    let err = run(&src, &dst, &with_sizes, Budget::new(2, 1 << 30)).unwrap_err();
    assert!(err.to_string().contains("in the chromosome sizes"), "{err}");

    let err = Conversion::between(Format::Bed, Format::BigBed, &src, &dst)
        .unwrap()
        .configure(&told)
        .unwrap_err();
    assert!(err.to_string().contains("ref_species"), "{err}");
}

#[test]
fn reference_inference_skips_other_paragraphs() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("x.maf");
    let dst = dir.path().join("x.bb");
    std::fs::write(
        &src,
        "##maf version=1\n\
         x ignored paragraph\n\
         a\n\
         s wrong.chr1 0 1 + 5 A\n\
         \n\
         a\n\
         s hg38.chr1 1 1 + 5 A\n",
    )
    .unwrap();
    let report = run(&src, &dst, &Options::default(), Budget::new(1, 1 << 30)).unwrap();
    assert_eq!(report.records, 1);
}

#[test]
fn bad_files_are_errors() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("x.maf");
    let dst = dir.path().join("x.bb");
    let a = "a score=1\ns hg38.chr1 10 5 + 100 ACGTA\n\n";
    for (body, says) in [
        (
            format!("{a}a score=2\ns hg38.chr1 12 4 + 100 ACGT\n"),
            "overlap",
        ),
        (
            format!("{a}a score=2\ns hg38.chr1 50 1 + 101 A\n"),
            "101 long here but 100",
        ),
        (
            format!("{a}a score=2\ns hg38.chr1 50 1 + 100 A\ns mm10.c 1 1 + 5 AA\n"),
            "line 7: alignment text of 2 columns",
        ),
    ] {
        std::fs::write(&src, format!("##maf version=1\n{body}")).unwrap();
        let err = run(&src, &dst, &Options::default(), Budget::new(1, 1 << 30)).unwrap_err();
        assert!(err.to_string().contains(says), "{says}: {err}");
    }
    std::fs::write(&src, "##maf\na\ns hg38.chr1 1 1 + 5 A\n").unwrap();
    let c = Conversion::between(Format::Maf, Format::BigBed, &src, &dst).unwrap();
    let err = convert(c, Budget::new(1, 1 << 30), None).unwrap_err();
    assert!(err.to_string().contains("no version"), "{err}");
}
