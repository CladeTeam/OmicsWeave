//! Allele counts against direct computation on read pairs drawn from a
//! generated genome with a phased sample's alleles at SNVs (molecules of
//! either haplotype, PCR copies, sequencing errors, `N`s, low qualities,
//! introns and deletions over sites, filtered records): the variant ×
//! barcode matrices in UMIs and the variant × file matrices in read pairs,
//! the per-variant table, gene-level haplotype counts; the same files at
//! any thread count and with spilled buffers.

mod common;

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use common::{
    Gene, Read, Rng, bcf_of_vcf, files, gene_hits, genes, genome, gtf, gunzip, read_10x_layer,
    reverse_for, write_bam,
};
use omicsweave_gtf::{Annotation, ReadOptions};
use omicsweave_quant::{
    AlleleGenes, Alleles, Budget, MatrixFormat, Product, Strandedness, Task, plan, quantify, run,
};

const REFS: [(&str, u32); 2] = [("chr1", 40_000), ("chr2", 25_000)];

/// A site of the panel: position, REF, ALT, and the genotype of the
/// first sample (`"0|1"`, `"0/1"`, ...).
#[derive(Debug, Clone)]
struct Snv {
    tid: usize,
    pos: u32,
    reference: u8,
    alt: u8,
    gt: &'static str,
}

impl Snv {
    /// The allele (0 or 1) of haplotype `h`.
    fn allele(&self, h: usize) -> usize {
        match self.gt.as_bytes()[2 * h] {
            b'1' => 1,
            _ => 0,
        }
    }

    /// Haplotype 1's allele at a phased heterozygous genotype.
    fn hap1(&self) -> Option<usize> {
        match self.gt {
            "0|1" => Some(0),
            "1|0" => Some(1),
            _ => None,
        }
    }
}

struct Fixture {
    dir: tempfile::TempDir,
    refs: Vec<(String, u32)>,
    reads: Vec<Read>,
    snvs: Vec<Snv>,
    genes: Vec<Gene>,
    ann: Arc<Annotation>,
    bam: PathBuf,
    bcf: PathBuf,
}

fn cigar(rng: &mut Rng, len: u32) -> Vec<(u32, char)> {
    let mut out = Vec::new();
    let mut left = len;
    if rng.chance(10) {
        let s = rng.range(1, 6) as u32;
        out.push((s, 'S'));
        left -= s;
    }
    while left > 0 {
        let m = (rng.range(15, 70) as u32).min(left);
        out.push((m, 'M'));
        left -= m;
        if left == 0 {
            break;
        }
        match rng.range(0, 8) {
            0 => out.push((rng.range(1, 5) as u32, 'D')),
            1 => out.push((rng.range(50, 400) as u32, 'N')),
            2 if left > 4 => {
                out.push((2, 'I'));
                left -= 2;
            }
            _ => {}
        }
    }
    out
}

fn fixture(seed: u64) -> Fixture {
    let dir = tempfile::tempdir().unwrap();
    let refs: Vec<(String, u32)> = REFS.iter().map(|(n, l)| (n.to_string(), *l)).collect();
    let mut rng = Rng::new(seed);
    let seqs = genome(&mut rng, &refs);
    let gs = genes(&mut rng, &refs, 6);
    let mut at = BTreeSet::new();
    for _ in 0..500 {
        let t = rng.range(0, 2) as usize;
        at.insert((t, rng.range(5, refs[t].1 as u64 - 5) as u32));
    }
    const GTS: [&str; 8] = ["0|1", "1|0", "0|1", "1|0", "0/1", "0|0", "1|1", "./."];
    let snvs: Vec<Snv> = at
        .iter()
        .map(|&(tid, pos)| {
            let reference = seqs[tid][pos as usize];
            let alt = *rng.pick(
                &b"ACGT"
                    .iter()
                    .copied()
                    .filter(|&b| b != reference)
                    .collect::<Vec<_>>(),
            );
            Snv {
                tid,
                pos,
                reference,
                alt,
                gt: GTS[rng.range(0, 8) as usize],
            }
        })
        .collect();
    let by_pos: HashMap<(usize, u32), &Snv> = snvs.iter().map(|s| ((s.tid, s.pos), s)).collect();

    // The panel: the SNVs, and records that are not biallelic SNVs.
    let mut vcf = String::from("##fileformat=VCFv4.2\n");
    for (n, l) in &refs {
        vcf.push_str(&format!("##contig=<ID={n},length={l}>\n"));
    }
    vcf.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    vcf.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n");
    for s in &snvs {
        vcf.push_str(&format!(
            "{}\t{}\t.\t{}\t{}\t.\t.\t.\tGT\t{}\t0|0\n",
            refs[s.tid].0,
            s.pos + 1,
            s.reference as char,
            s.alt as char,
            s.gt
        ));
        if rng.chance(5) {
            // An indel or a multiallelic site right after it.
            let r = seqs[s.tid][s.pos as usize + 1];
            let others: Vec<char> = b"ACGT"
                .iter()
                .filter(|&&b| b != r)
                .map(|&b| b as char)
                .collect();
            let r = r as char;
            let line = if rng.chance(50) {
                format!("{r}\t{r}A")
            } else {
                format!("{r}\t{},{}", others[0], others[1])
            };
            vcf.push_str(&format!(
                "{}\t{}\t.\t{line}\t.\t.\t.\tGT\t0|1\t0|0\n",
                refs[s.tid].0,
                s.pos + 2
            ));
        }
    }
    let bcf = bcf_of_vcf(dir.path(), "panel", &vcf);

    let mut reads = Vec::new();
    let mut n = 0;
    for m in 0..1500 {
        let tid = rng.range(0, 2) as usize;
        let start = rng.range(0, (refs[tid].1 - 3000) as u64) as u32;
        let hap = rng.range(0, 2) as usize;
        let cb = rng.chance(95).then(|| format!("C{}", rng.range(0, 8)));
        let ub: Option<String> = rng
            .chance(97)
            .then(|| (0..6).map(|_| *rng.pick(&['A', 'C', 'G', 'T'])).collect());
        let reverse = rng.chance(50);
        for _ in 0..rng.range(1, 4) {
            n += 1;
            let h = if rng.chance(4) { 1 - hap } else { hap };
            let s = start + rng.range(0, 40) as u32;
            let mut mates = Vec::new();
            for k in 0..2 {
                let len = rng.range(50, 110) as u32;
                let c = cigar(&mut rng, len);
                let pos = s + k * rng.range(20, 200) as u32;
                let mut bases = Vec::new();
                let mut at = pos;
                for &(l, o) in &c {
                    match o {
                        'M' => {
                            for i in 0..l {
                                let p = at + i;
                                let mut b = match by_pos.get(&(tid, p)) {
                                    Some(v) if v.allele(h) == 1 => v.alt,
                                    _ => seqs[tid][p as usize],
                                };
                                if rng.chance(2) {
                                    b = *rng.pick(b"ACGT");
                                } else if rng.chance(1) {
                                    b = b'N';
                                }
                                bases.push(b);
                            }
                            at += l;
                        }
                        'I' | 'S' => bases.extend((0..l).map(|_| *rng.pick(b"ACGT"))),
                        _ => at += l,
                    }
                }
                let qual = (!rng.chance(10)).then(|| {
                    (0..bases.len())
                        .map(|_| b'!' + rng.range(5, 41) as u8)
                        .collect()
                });
                mates.push((pos, c, bases, qual));
            }
            let dup = rng.chance(3);
            let ends: Vec<u32> = mates.iter().map(|m| m.0).collect();
            for (k, (pos, c, bases, qual)) in mates.into_iter().enumerate() {
                let mut flag = 0x1 | 0x2 | if k == 0 { 0x40 } else { 0x80 };
                if (k == 1) != reverse {
                    flag |= 0x10;
                } else {
                    flag |= 0x20;
                }
                if dup {
                    flag |= 0x400;
                }
                match rng.range(0, 100) {
                    0 => flag |= 0x100,
                    1 => flag |= 0x200,
                    _ => {}
                }
                let mut tags = Vec::new();
                if let Some(cb) = &cb {
                    tags.push(format!("CB:Z:{cb}"));
                }
                if let Some(ub) = &ub {
                    tags.push(format!("UB:Z:{ub}"));
                }
                reads.push(Read {
                    name: format!("m{m}r{n}"),
                    flag,
                    tid,
                    pos,
                    mapq: if rng.chance(5) { 10 } else { 60 },
                    cigar: c,
                    mtid: tid as i32,
                    mpos: ends[1 - k],
                    tlen: 0,
                    tags,
                    seq: Some(bases),
                    qual,
                });
            }
        }
    }
    reads.sort_by_key(|r| (r.tid, r.pos));
    let bam = write_bam(dir.path(), "rna", &refs, &reads, seed);
    let g = dir.path().join("genes.gtf");
    std::fs::write(&g, gtf(&refs, &gs)).unwrap();
    let ann = Arc::new(Annotation::read(&g, &ReadOptions::default()).unwrap());
    Fixture {
        dir,
        refs,
        reads,
        snvs,
        genes: gs,
        ann,
        bam,
        bcf,
    }
}

fn site_id(f: &Fixture, s: &Snv) -> String {
    format!(
        "{}:{}:{}>{}",
        f.refs[s.tid].0,
        s.pos + 1,
        s.reference as char,
        s.alt as char
    )
}

/// The aligned query position over reference position `p`.
fn query_at(r: &Read, p: u32) -> Option<usize> {
    common::aligned_pairs(r)
        .into_iter()
        .find(|x| x.0 == p)
        .map(|x| x.1)
}

type Layers = [BTreeMap<(String, String), u32>; 3];

/// The expected `ref`, `alt`, `other` layers and gene `hap1`, `hap2`,
/// `both` layers, units by `unit` (UMI or read name), rows by `row`.
fn expected(
    f: &Fixture,
    exclude: u16,
    row: &dyn Fn(&Read) -> Option<String>,
    unit: &dyn Fn(&Read) -> Option<String>,
    genes: bool,
) -> (Layers, Layers) {
    let mut units: BTreeSet<(String, String, usize, String)> = BTreeSet::new();
    let mut haps: BTreeMap<(String, String, String), u32> = BTreeMap::new();
    for r in &f.reads {
        if r.flag & exclude != 0 || r.mapq < 20 {
            continue;
        }
        let (Some(key), Some(sub)) = (row(r), unit(r)) else {
            continue;
        };
        let mut bits = 0;
        for s in f.snvs.iter().filter(|s| s.tid == r.tid) {
            if s.pos < r.pos || s.pos >= r.end() {
                continue;
            }
            let Some(q) = query_at(r, s.pos) else {
                continue;
            };
            if r.qual.as_ref().is_some_and(|v| v[q] - 33 < 20) {
                continue;
            }
            let b = r.seq.as_ref().unwrap()[q];
            let allele = if b == s.reference {
                0
            } else if b == s.alt {
                1
            } else if b"ACGT".contains(&b) {
                2
            } else {
                continue;
            };
            units.insert((key.clone(), site_id(f, s), allele, sub.clone()));
            if let (Some(h), true) = (s.hap1(), allele < 2) {
                bits |= if allele == h { 1 } else { 2 };
            }
        }
        if genes && bits != 0 {
            let hits = gene_hits(&f.genes, r.tid, &r.blocks(false), Some(reverse_for(r, 1)));
            if hits.len() == 1 {
                *haps
                    .entry((key.clone(), f.genes[hits[0]].id.clone(), sub.clone()))
                    .or_insert(0) |= bits;
            }
        }
    }
    let mut v: Layers = Default::default();
    for (k, s, a, _) in units {
        *v[a].entry((k, s)).or_insert(0) += 1;
    }
    let mut g: Layers = Default::default();
    for ((k, gene, _), bits) in haps {
        *g[bits as usize - 1].entry((k, gene)).or_insert(0) += 1;
    }
    (v, g)
}

fn task(f: &Fixture, out: &Path, products: Vec<Product>) -> Task {
    let mut t = Task::new(vec![f.bam.clone()], out);
    t.annotation = Some(Arc::clone(&f.ann));
    t.matrix = MatrixFormat::Mtx;
    t.products = products;
    t
}

fn products(f: &Fixture) -> Vec<Product> {
    let mut sc = Alleles::single_cell(f.bcf.clone());
    sc.genes = Some(AlleleGenes {
        sample: None,
        strand: Strandedness::Forward,
    });
    let mut bulk = Alleles::bulk(f.bcf.clone());
    bulk.name = "bulk".into();
    vec![Product::Alleles(sc), Product::Alleles(bulk)]
}

#[test]
fn allele_matrices_and_table_match_the_reads() {
    let f = fixture(91);
    let out = f.dir.path().join("out");
    quantify(
        task(&f, &out, products(&f)),
        Budget::new(3, 256 << 20),
        None,
    )
    .unwrap();

    let cb = |r: &Read| r.tag("CB").map(str::to_owned);
    let ub = |r: &Read| r.tag("UB").map(str::to_owned);
    let (want, want_genes) = expected(&f, 0x304, &cb, &ub, true);
    assert!(want.iter().all(|l| !l.is_empty()));
    let dir = out.join("alleles");
    for (i, layer) in ["ref", "alt", "other"].iter().enumerate() {
        assert_eq!(read_10x_layer(&dir, layer), want[i], "{layer}");
    }
    assert!(want_genes[0].len() > 5 && want_genes[1].len() > 5);
    for (i, layer) in ["hap1", "hap2", "both"].iter().enumerate() {
        assert_eq!(
            read_10x_layer(&dir.join("genes"), layer),
            want_genes[i],
            "{layer}"
        );
    }
    // Every panel SNV is a column, in file order, the others skipped.
    let ids: Vec<String> = gunzip(&dir.join("features.tsv.gz"))
        .lines()
        .map(|l| l.split('\t').next().unwrap().to_owned())
        .collect();
    assert_eq!(
        ids,
        f.snvs.iter().map(|s| site_id(&f, s)).collect::<Vec<_>>()
    );
    // The table sums each layer over the rows.
    let table = gunzip(&dir.join("variants.tsv.gz"));
    let mut lines = table.lines();
    assert_eq!(
        lines.next().unwrap(),
        "chrom\tpos\tref\talt\tref_count\talt_count\tother_count\talt_fraction"
    );
    for (line, s) in lines.zip(&f.snvs) {
        let id = site_id(&f, s);
        let n = |l: usize| -> u64 {
            want[l]
                .iter()
                .filter(|((_, x), _)| *x == id)
                .map(|(_, v)| *v as u64)
                .sum()
        };
        let (a, b, c) = (n(0), n(1), n(2));
        let af = if a + b > 0 {
            format!("{:.6}", b as f64 / (a + b) as f64)
        } else {
            String::new()
        };
        assert_eq!(
            line,
            format!(
                "{}\t{}\t{}\t{}\t{a}\t{b}\t{c}\t{af}",
                f.refs[s.tid].0,
                s.pos + 1,
                s.reference as char,
                s.alt as char
            )
        );
    }

    // Bulk: read pairs by name, per file, duplicates excluded.
    let bam = f.bam.display().to_string();
    let file = |_: &Read| Some(bam.clone());
    let name = |r: &Read| Some(r.name.clone());
    let (want, _) = expected(&f, 0x704, &file, &name, false);
    for (i, layer) in ["ref", "alt", "other"].iter().enumerate() {
        assert_eq!(read_10x_layer(&out.join("bulk"), layer), want[i], "{layer}");
    }
    let summary = std::fs::read_to_string(out.join("summary.json")).unwrap();
    assert!(summary.contains(&format!("\"sites\": {},", f.snvs.len())));
    assert!(summary.contains("\"sites_skipped_not_biallelic_snv\""));
}

#[test]
fn allele_outputs_do_not_depend_on_threads_or_spilling() {
    let f = fixture(92);
    let a = f.dir.path().join("a");
    quantify(task(&f, &a, products(&f)), Budget::new(1, 1 << 30), None).unwrap();
    let b = f.dir.path().join("b");
    let mut p = plan(task(&f, &b, products(&f)), Budget::new(6, 256 << 20)).unwrap();
    for s in &mut p.stages {
        s.1 = 8 << 10;
    }
    run(&p, None).unwrap();
    let summary = std::fs::read_to_string(b.join("summary.json")).unwrap();
    assert!(
        summary
            .lines()
            .any(|l| l.contains("spilled_bytes") && !l.contains(": 0,")),
        "the small buffers spill"
    );
    let strip = |v: Vec<(String, Vec<u8>)>| -> Vec<(String, Vec<u8>)> {
        v.into_iter().filter(|(n, _)| n != "summary.json").collect()
    };
    assert_eq!(strip(files(&a)), strip(files(&b)));
    // A sample that is not in the panel is an error before any output.
    let mut sc = Alleles::single_cell(f.bcf.clone());
    sc.genes = Some(AlleleGenes {
        sample: Some("S9".into()),
        strand: Strandedness::Forward,
    });
    let err = plan(
        task(&f, &f.dir.path().join("c"), vec![Product::Alleles(sc)]),
        Budget::new(1, 256 << 20),
    )
    .unwrap_err();
    assert!(err.to_string().contains("S9"), "{err}");
}

#[test]
fn placed_unmapped_spans_and_absent_sequences_do_not_add_observations() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("reads.sam");
    std::fs::write(
        &sam,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
u\t4\tchr1\t101\t60\t4M\t*\t0\t0\tACAA\t*\n\
m\t0\tchr1\t201\t60\t4M\t*\t0\t0\t*\t*\n\
n\t0\tchr1\t301\t60\t4M\t*\t0\t0\tACAA\t*\n",
    )
    .unwrap();
    let bam = dir.path().join("reads.bam");
    common::bam_of(&sam, &bam);
    let mut vcf = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    for p in [100, 101, 200, 201, 300, 301] {
        let (r, a) = if p % 100 == 0 { ('A', 'C') } else { ('C', 'A') };
        vcf.push_str(&format!("chr1\t{}\t.\t{r}\t{a}\t.\t.\t.\n", p + 1));
    }
    let panel = bcf_of_vcf(dir.path(), "sites", &vcf);
    let mut spec = Alleles::bulk(panel);
    spec.filter.exclude = 0;
    let out = dir.path().join("out");
    let mut task = Task::new(vec![bam], &out);
    task.matrix = MatrixFormat::Mtx;
    task.products = vec![Product::Alleles(spec)];
    quantify(task, Budget::new(1, 128 << 20), None).unwrap();
    let layer = read_10x_layer(&out.join("alleles"), "ref");
    let features: Vec<_> = layer
        .into_iter()
        .map(|((_, feature), n)| (feature, n))
        .collect();
    assert_eq!(
        features,
        [
            ("chr1:101:A>C".into(), 1),
            ("chr1:301:A>C".into(), 1),
            ("chr1:302:C>A".into(), 1)
        ]
    );
}
