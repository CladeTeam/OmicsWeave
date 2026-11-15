//! Throughput on a synthetic single-cell library: writes a coordinate-sorted
//! BAM of paired reads with `CB` / `UB` / `NH` tags (and a GTF of genes over
//! the same references), then runs coverage, single-cell gene counts and
//! the ATAC products at several thread counts.
//!
//! ```text
//! cargo run --release -p omicsweave-quant --example quant_bench -- DIR [PAIRS] [THREADS,...]
//! ```

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use omicsweave_bgzf::Writer;
use omicsweave_core::exec::Pool;
use omicsweave_gtf::{Annotation, ReadOptions};
use omicsweave_quant::{Atac, Budget, Counts, Coverage, Product, Task, human, plan, run};

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }
}

const CHROMS: usize = 24;
const CHROM_LEN: u32 = 50_000_000;
const READ_LEN: u32 = 90;

/// `hts_reg2bin(beg, end, 14, 5)`.
fn reg2bin(beg: u32, end: u32) -> u16 {
    let end = end - 1;
    for (shift, offset) in [(14u32, 4681u32), (17, 585), (20, 73), (23, 9), (26, 1)] {
        if beg >> shift == end >> shift {
            return (offset + (beg >> shift)) as u16;
        }
    }
    0
}

struct Rec {
    pos: u32,
    bytes: Vec<u8>,
}

#[allow(clippy::too_many_arguments)]
fn record(
    tid: i32,
    pos: u32,
    flag: u16,
    cigar: &[(u32, u8)],
    mpos: u32,
    tlen: i32,
    name: &[u8],
    tags: &[u8],
    rng: &mut Rng,
) -> Rec {
    let span: u32 = cigar
        .iter()
        .filter(|(_, o)| matches!(o, 0 | 2 | 3))
        .map(|(l, _)| l)
        .sum();
    let l_seq = READ_LEN as usize;
    let mut b = Vec::with_capacity(64 + name.len() + 4 * cigar.len() + l_seq * 2 + tags.len());
    b.extend_from_slice(&[0; 4]);
    b.extend_from_slice(&tid.to_le_bytes());
    b.extend_from_slice(&(pos as i32).to_le_bytes());
    b.push(name.len() as u8 + 1);
    b.push(60);
    b.extend_from_slice(&reg2bin(pos, pos + span.max(1)).to_le_bytes());
    b.extend_from_slice(&(cigar.len() as u16).to_le_bytes());
    b.extend_from_slice(&flag.to_le_bytes());
    b.extend_from_slice(&(l_seq as i32).to_le_bytes());
    b.extend_from_slice(&tid.to_le_bytes());
    b.extend_from_slice(&(mpos as i32).to_le_bytes());
    b.extend_from_slice(&tlen.to_le_bytes());
    b.extend_from_slice(name);
    b.push(0);
    for &(l, o) in cigar {
        b.extend_from_slice(&(l << 4 | o as u32).to_le_bytes());
    }
    for _ in 0..l_seq.div_ceil(2) {
        let x = rng.next();
        b.push([0x12, 0x48, 0x84, 0x21][(x & 3) as usize]);
    }
    b.extend(std::iter::repeat_n(30u8, l_seq));
    b.extend_from_slice(tags);
    let size = (b.len() - 4) as i32;
    b[..4].copy_from_slice(&size.to_le_bytes());
    Rec { pos, bytes: b }
}

fn cigar(rng: &mut Rng) -> Vec<(u32, u8)> {
    if rng.below(100) < 20 {
        let a = 20 + rng.below(50) as u32;
        vec![(a, 0), (100 + rng.below(5000) as u32, 3), (READ_LEN - a, 0)]
    } else {
        vec![(READ_LEN, 0)]
    }
}

fn barcode(i: u64) -> Vec<u8> {
    let mut v: Vec<u8> = (0..16)
        .map(|k| b"ACGT"[(i >> (2 * k) & 3) as usize])
        .collect();
    v.extend_from_slice(b"-1");
    v
}

/// Writes `pairs` read pairs as `dir/bench.bam`, on `threads` compressors.
fn write_bam(dir: &Path, pairs: u64, threads: usize) -> PathBuf {
    let path = dir.join("bench.bam");
    let pool = Pool::new(threads);
    let file = BufWriter::with_capacity(4 << 20, File::create(&path).unwrap());
    let mut w = Writer::new(file, 1, Box::new(pool.ordered(4 * threads))).unwrap();
    let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for c in 0..CHROMS {
        text.push_str(&format!("@SQ\tSN:chr{}\tLN:{CHROM_LEN}\n", c + 1));
    }
    let mut head = b"BAM\x01".to_vec();
    head.extend_from_slice(&(text.len() as i32).to_le_bytes());
    head.extend_from_slice(text.as_bytes());
    head.extend_from_slice(&(CHROMS as i32).to_le_bytes());
    for c in 0..CHROMS {
        let n = format!("chr{}", c + 1);
        head.extend_from_slice(&(n.len() as i32 + 1).to_le_bytes());
        head.extend_from_slice(n.as_bytes());
        head.push(0);
        head.extend_from_slice(&CHROM_LEN.to_le_bytes());
    }
    w.write(&head).unwrap();
    w.flush_block().unwrap();
    let mut rng = Rng(42);
    let per_chrom = pairs / CHROMS as u64;
    let mut serial = 0u64;
    for c in 0..CHROMS {
        let mut recs = Vec::with_capacity(2 * per_chrom as usize);
        // Reads pile up in "genes": a third of the genome holds most of them.
        for _ in 0..per_chrom {
            serial += 1;
            let hot = rng.below(100) < 80;
            let s = if hot {
                let locus = rng.below(CHROM_LEN as u64 / 30_000) * 30_000;
                (locus + rng.below(10_000)) as u32
            } else {
                rng.below(CHROM_LEN as u64 - 20_000) as u32
            };
            let frag = 150 + rng.below(450) as u32;
            let (c1, c2) = (cigar(&mut rng), cigar(&mut rng));
            let span2: u32 = c2.iter().filter(|(_, o)| *o != 1).map(|(l, _)| l).sum();
            let p2 = s + frag.max(span2) - span2;
            let end = (s + READ_LEN).max(p2 + span2);
            let tlen = (end - s) as i32;
            let mut tags = Vec::new();
            tags.extend_from_slice(b"CBZ");
            tags.extend_from_slice(&barcode(rng.below(5000)));
            tags.push(0);
            tags.extend_from_slice(b"UBZ");
            for _ in 0..12 {
                tags.push(b"ACGT"[rng.below(4) as usize]);
            }
            tags.push(0);
            tags.extend_from_slice(b"NHC");
            tags.push(if rng.below(100) < 5 { 2 } else { 1 });
            let name = format!("r{serial}");
            let r1_left = rng.below(2) == 0;
            let (f1, f2) = if r1_left { (0x40, 0x80) } else { (0x80, 0x40) };
            recs.push(record(
                c as i32,
                s,
                0x1 | 0x2 | 0x20 | f1,
                &c1,
                p2,
                tlen,
                name.as_bytes(),
                &tags,
                &mut rng,
            ));
            recs.push(record(
                c as i32,
                p2,
                0x1 | 0x2 | 0x10 | f2,
                &c2,
                s,
                -tlen,
                name.as_bytes(),
                &tags,
                &mut rng,
            ));
        }
        recs.sort_by_key(|r| r.pos);
        for r in recs {
            w.write_record(&r.bytes).unwrap();
        }
    }
    let done = w.finish().unwrap();
    done.out.into_inner().unwrap().flush().unwrap();
    path
}

/// About 20k genes, one per 60 kb, two transcripts of up to eight exons.
fn write_gtf(dir: &Path) -> PathBuf {
    let path = dir.join("bench.gtf");
    let mut w = BufWriter::new(File::create(&path).unwrap());
    let mut rng = Rng(7);
    let mut g = 0;
    for c in 0..CHROMS {
        let mut at = 5_000u32;
        while at + 60_000 < CHROM_LEN {
            g += 1;
            let strand = if rng.below(2) == 0 { '+' } else { '-' };
            let len = 5_000 + rng.below(40_000) as u32;
            let a = format!("gene_id \"G{g}\"; gene_name \"N{g}\"; gene_type \"protein_coding\";");
            writeln!(
                w,
                "chr{}\tb\tgene\t{}\t{}\t.\t{strand}\t.\t{a}",
                c + 1,
                at + 1,
                at + len
            )
            .unwrap();
            for t in 0..2 {
                let k = 2 + rng.below(7) as u32;
                let piece = len / k;
                writeln!(
                    w,
                    "chr{}\tb\ttranscript\t{}\t{}\t.\t{strand}\t.\t{a} transcript_id \"T{g}.{t}\";",
                    c + 1,
                    at + 1,
                    at + len
                )
                .unwrap();
                for e in 0..k {
                    let s = at + e * piece + rng.below(piece as u64 / 4) as u32;
                    let x = s + 100 + rng.below(piece as u64 / 2) as u32;
                    writeln!(
                        w,
                        "chr{}\tb\texon\t{}\t{}\t.\t{strand}\t.\t{a} transcript_id \"T{g}.{t}\";",
                        c + 1,
                        s + 1,
                        x.min(at + len)
                    )
                    .unwrap();
                }
            }
            at += 60_000;
        }
    }
    path
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let dir = PathBuf::from(
        args.get(1)
            .expect("usage: quant_bench DIR [PAIRS] [THREADS,...]"),
    );
    let pairs: u64 = args.get(2).map_or(10_000_000, |s| s.parse().unwrap());
    let threads: Vec<usize> = args.get(3).map_or(vec![1, 8, 16], |s| {
        s.split(',').map(|t| t.parse().unwrap()).collect()
    });
    std::fs::create_dir_all(&dir).unwrap();
    let bam = dir.join("bench.bam");
    if !bam.exists() {
        let t = Instant::now();
        write_bam(&dir, pairs, 16);
        eprintln!(
            "wrote {} pairs in {:.1} s",
            pairs,
            t.elapsed().as_secs_f64()
        );
    }
    let gtf = write_gtf(&dir);
    let ann = Arc::new(Annotation::read(&gtf, &ReadOptions::default()).unwrap());
    let size = std::fs::metadata(&bam).unwrap().len();
    println!(
        "input: {} ({}), {} genes",
        bam.display(),
        human(size),
        ann.genes().len()
    );
    let tasks: Vec<(&str, Vec<Product>)> = vec![
        ("coverage", vec![Product::Coverage(Coverage::default())]),
        ("sc-genes", vec![Product::Counts(Counts::rna_single_cell())]),
        (
            "atac",
            vec![Product::Atac(Atac {
                filter: omicsweave_quant::Filter {
                    min_mapq: 0,
                    ..Atac::default().filter
                },
                ..Atac::default()
            })],
        ),
        (
            "all",
            vec![
                Product::Coverage(Coverage::default()),
                Product::Counts(Counts::rna_single_cell()),
                Product::Atac(Atac::default()),
            ],
        ),
    ];
    for (name, products) in tasks {
        for &t in &threads {
            let out = dir.join(format!("out-{name}-{t}"));
            let _ = std::fs::remove_dir_all(&out);
            let mut task = Task::new(vec![bam.clone()], &out);
            task.annotation = Some(Arc::clone(&ann));
            task.products = products.clone();
            let p = plan(task, Budget::new(t, 16 << 30)).unwrap();
            let start = Instant::now();
            let r = run(&p, None).unwrap();
            let s = start.elapsed().as_secs_f64();
            println!(
                "{name:>9} threads {t:>2}: {:6.2} s (scan {:5.2} s), {:5.2} M records/s, {} compressed/s",
                s,
                r.scan.as_secs_f64(),
                r.records as f64 / s / 1e6,
                human((size as f64 / s) as u64)
            );
            let _ = std::fs::remove_dir_all(&out);
        }
    }
}
