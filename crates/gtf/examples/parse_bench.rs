//! Parses a synthetic GENCODE-scale GTF and reports time and memory.
//!
//! ```text
//! cargo run --release -p omicsweave-gtf --example parse_bench -- \
//!     [--dir DIR] [--genes N] [--threads N] [--attributes all|none] [--gzip]
//! ```
//!
//! The file (`DIR/synthetic.gtf`, default: the system temp dir) is written
//! once, deterministically: about 35 lines per gene, GENCODE attributes;
//! the default 47,000 genes give about 1.6 million lines. Then it is read, and the
//! interval indexes, derived views and a sweep over sorted positions are
//! timed. Peak RSS is the process high-water mark (`VmHWM`).

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

use omicsweave_gtf::{Annotation, AttributePolicy, Flank, Level, ReadOptions};

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }
}

const CHROMS: [(&str, u32); 25] = [
    ("chr1", 248956422),
    ("chr2", 242193529),
    ("chr3", 198295559),
    ("chr4", 190214555),
    ("chr5", 181538259),
    ("chr6", 170805979),
    ("chr7", 159345973),
    ("chr8", 145138636),
    ("chr9", 138394717),
    ("chr10", 133797422),
    ("chr11", 135086622),
    ("chr12", 133275309),
    ("chr13", 114364328),
    ("chr14", 107043718),
    ("chr15", 101991189),
    ("chr16", 90338345),
    ("chr17", 83257441),
    ("chr18", 80373285),
    ("chr19", 58617616),
    ("chr20", 64444167),
    ("chr21", 46709983),
    ("chr22", 50818468),
    ("chrX", 156040895),
    ("chrY", 57227415),
    ("chrM", 16569),
];

fn generate(path: &PathBuf, genes: usize) -> std::io::Result<u64> {
    let mut rng = Rng(42);
    let mut w = BufWriter::with_capacity(1 << 20, File::create(path)?);
    writeln!(w, "##description: synthetic GENCODE-style annotation")?;
    writeln!(w, "##provider: synthetic")?;
    writeln!(w, "##format: gtf")?;
    let total: u64 = CHROMS[..24].iter().map(|c| c.1 as u64).sum();
    let mut lines = 3u64;
    let types = ["protein_coding", "lncRNA", "processed_pseudogene", "miRNA"];
    for (c, &(chrom, len)) in CHROMS.iter().enumerate() {
        let n = if c == 24 {
            13
        } else {
            (genes as u64 * len as u64 / total) as usize
        };
        let mut pos = 10_000u32;
        let step = (len as u64 / (n as u64 + 1)) as u32;
        for gi in 0..n {
            pos = pos.saturating_add(step / 2 + rng.below(step as u64) as u32);
            if pos as u64 + 3_000_000 > len as u64 {
                break;
            }
            let strand = if rng.below(2) == 0 { '+' } else { '-' };
            let gid = format!("ENSG{:011}.{}", c * 100_000 + gi, 1 + rng.below(9));
            let gname = format!("GENE{c}_{gi}");
            let gtype = types[rng.below(4) as usize];
            let coding = gtype == "protein_coding";
            let n_tx = 1 + rng.below(8) as usize;
            let mut txs = Vec::new();
            let (mut gs, mut ge) = (u32::MAX, 0u32);
            for t in 0..n_tx {
                let n_ex = 1 + rng.below(9) as usize;
                let mut at = pos + rng.below(2000) as u32;
                let mut exons = Vec::new();
                for _ in 0..n_ex {
                    let l = 50 + rng.below(400) as u32;
                    exons.push((at, at + l));
                    at += l + 200 + rng.below(20_000) as u32;
                }
                gs = gs.min(exons[0].0);
                ge = ge.max(exons.last().unwrap().1);
                txs.push((t, exons));
            }
            let gattr = format!(
                "gene_id \"{gid}\"; gene_type \"{gtype}\"; gene_name \"{gname}\"; level 2;"
            );
            writeln!(
                w,
                "{chrom}\tHAVANA\tgene\t{}\t{ge}\t.\t{strand}\t.\t{gattr} hgnc_id \"HGNC:{gi}\"; havana_gene \"OTTHUMG{gi:011}.1\";",
                gs + 1
            )?;
            lines += 1;
            for (t, exons) in txs {
                let tid = format!("ENST{:011}.{}", (c * 100_000 + gi) * 10 + t, 1 + t % 3);
                let tname = format!("{gname}-{}", 201 + t);
                let tattr = format!(
                    "gene_id \"{gid}\"; transcript_id \"{tid}\"; gene_type \"{gtype}\"; gene_name \"{gname}\"; transcript_type \"{gtype}\"; transcript_name \"{tname}\"; level 2;"
                );
                let tags = if t == 0 {
                    " tag \"basic\"; tag \"Ensembl_canonical\"; tag \"MANE_Select\";"
                } else {
                    " tag \"basic\";"
                };
                let (ts, te) = (exons[0].0, exons.last().unwrap().1);
                writeln!(
                    w,
                    "{chrom}\tHAVANA\ttranscript\t{}\t{te}\t.\t{strand}\t.\t{tattr} transcript_support_level \"1\";{tags} havana_transcript \"OTTHUMT{t:011}.1\";",
                    ts + 1
                )?;
                lines += 1;
                let n_ex = exons.len();
                for (k, &(s, e)) in exons.iter().enumerate() {
                    let number = if strand == '+' { k + 1 } else { n_ex - k };
                    let common = format!(
                        "{tattr} exon_number {number}; exon_id \"ENSE{:011}.1\";{tags}",
                        (c * 100_000 + gi) * 100 + t * 10 + k
                    );
                    writeln!(
                        w,
                        "{chrom}\tHAVANA\texon\t{}\t{e}\t.\t{strand}\t.\t{common}",
                        s + 1
                    )?;
                    lines += 1;
                    if coding && n_ex > 2 && k > 0 && k + 1 < n_ex {
                        writeln!(
                            w,
                            "{chrom}\tHAVANA\tCDS\t{}\t{e}\t.\t{strand}\t{}\t{common} protein_id \"ENSP{gi:011}.1\";",
                            s + 1,
                            k % 3
                        )?;
                        lines += 1;
                    } else if coding && n_ex > 2 {
                        writeln!(
                            w,
                            "{chrom}\tHAVANA\tUTR\t{}\t{e}\t.\t{strand}\t.\t{common}",
                            s + 1
                        )?;
                        lines += 1;
                    }
                }
                if coding && n_ex > 2 {
                    let (s, _) = exons[1];
                    let (_, e) = exons[n_ex - 2];
                    let (start, stop) = if strand == '+' {
                        ((s, s + 3), (e - 3, e))
                    } else {
                        ((e - 3, e), (s, s + 3))
                    };
                    for (kind, (cs, ce)) in [("start_codon", start), ("stop_codon", stop)] {
                        writeln!(
                            w,
                            "{chrom}\tHAVANA\t{kind}\t{}\t{ce}\t.\t{strand}\t0\t{tattr}",
                            cs + 1
                        )?;
                    }
                    lines += 2;
                }
            }
        }
    }
    w.flush()?;
    Ok(lines)
}

fn status(key: &str) -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with(key))
                .and_then(|l| l.split_whitespace().nth(1)?.parse().ok())
        })
        .unwrap_or(0)
}

fn mib(bytes: u64) -> f64 {
    bytes as f64 / (1 << 20) as f64
}

fn main() {
    let mut dir = std::env::temp_dir();
    let mut genes = 47_000;
    let mut threads = 0;
    let mut policy = AttributePolicy::All;
    let mut gzip = false;
    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--dir" => dir = PathBuf::from(args.next().expect("--dir DIR")),
            "--genes" => genes = args.next().expect("--genes N").parse().expect("number"),
            "--threads" => threads = args.next().expect("--threads N").parse().expect("number"),
            "--attributes" => {
                policy = match args.next().as_deref() {
                    Some("none") => AttributePolicy::None,
                    _ => AttributePolicy::All,
                }
            }
            "--gzip" => gzip = true,
            other => panic!("unknown argument {other}"),
        }
    }
    let path = dir.join(format!("synthetic-{genes}.gtf"));
    if !path.exists() {
        let t = Instant::now();
        let lines = generate(&path, genes).expect("write synthetic GTF");
        eprintln!("wrote {lines} lines in {:.2?}", t.elapsed());
    }
    let mut source = path.clone();
    if gzip {
        source = dir.join(format!("synthetic-{genes}.gtf.gz"));
        if !source.exists() {
            let mut enc = flate2::write::GzEncoder::new(
                File::create(&source).unwrap(),
                flate2::Compression::new(6),
            );
            std::io::copy(&mut File::open(&path).unwrap(), &mut enc).unwrap();
            enc.finish().unwrap();
        }
    }
    let size = std::fs::metadata(&source).unwrap().len();
    let rss0 = status("VmRSS:") * 1024;
    let opts = ReadOptions {
        threads,
        attributes: policy.clone(),
        ..ReadOptions::default()
    };
    let t = Instant::now();
    let ann = Annotation::read(&source, &opts).expect("read");
    let parse = t.elapsed();
    let peak = status("VmHWM:") * 1024;
    let rss = status("VmRSS:") * 1024;
    let r = ann.report();
    println!(
        "file {} ({:.1} MiB), threads {}, attributes {:?}",
        source.display(),
        mib(size),
        if threads == 0 {
            "auto".to_string()
        } else {
            threads.to_string()
        },
        policy
    );
    println!(
        "records {} genes {} transcripts {} features {} problems {}",
        r.records,
        r.genes,
        r.transcripts,
        r.features,
        r.problems.len()
    );
    println!(
        "parse {:.3} s ({:.0} MiB/s, {:.2} M lines/s)",
        parse.as_secs_f64(),
        mib(size) / parse.as_secs_f64(),
        r.records as f64 / parse.as_secs_f64() / 1e6
    );
    println!(
        "model heap {:.1} MiB; RSS {:.1} MiB after read (from {:.1}); peak RSS {:.1} MiB",
        mib(ann.heap_bytes() as u64),
        mib(rss),
        mib(rss0),
        mib(peak)
    );
    let t = Instant::now();
    for level in [Level::Gene, Level::Transcript, Level::Exon] {
        ann.interval_index(level);
    }
    println!(
        "indexes (gene, transcript, exon) {:.3} s",
        t.elapsed().as_secs_f64()
    );
    let t = Instant::now();
    let introns = ann.introns().len();
    let exonic = ann.exonic_regions().len();
    let junctions = ann.junctions().len();
    let bins = ann.exon_bins(false).set.len();
    let promoters = ann
        .derive()
        .promoters(Level::Transcript, Flank::new(2000, 200))
        .len();
    println!(
        "derived views {:.3} s: introns {introns}, exonic regions {exonic}, junctions {junctions}, exon bins {bins}, promoters {promoters}",
        t.elapsed().as_secs_f64()
    );
    // A sorted scan: 100-base reads every 25 bases on chr1, spliced into
    // two blocks.
    let s = ann.seqid("chr1").unwrap();
    let mut cursor = ann.cursor(Level::Exon, s);
    let mut hits = Vec::new();
    let (mut queries, mut found) = (0u64, 0u64);
    let t = Instant::now();
    let mut p = 0u32;
    while p < 248_000_000 {
        cursor.advance(p);
        hits.clear();
        cursor.overlaps(p, p + 50, &mut hits);
        cursor.overlaps(p + 1050, p + 1100, &mut hits);
        found += hits.len() as u64;
        queries += 2;
        p += 25;
    }
    let el = t.elapsed().as_secs_f64();
    println!(
        "exon sweep on chr1: {queries} block queries in {el:.3} s ({:.1} ns/query), {found} hits",
        el * 1e9 / queries as f64
    );
    println!("peak RSS {:.1} MiB", mib(status("VmHWM:") * 1024));
}
