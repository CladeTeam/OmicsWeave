//! Metabolic labeling against direct computation on reads drawn from a
//! generated genome with T>C conversions, sequencing errors and SNPs: the
//! per-read histogram, the mismatch table by cycle, labeled layers of
//! single-cell and bulk gene matrices; the same files from a bgzip FASTA,
//! a 2bit genome and `MD` tags, at any thread count.

mod common;

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use common::{
    Gene, Read, Rng, aligned_pairs, complement, gene_hits, genes, genome, gtf, gunzip, md_tag,
    read_10x_layer, reverse_for, write_bam, write_bcf, write_fasta, write_twobit,
};
use omicsweave_gtf::{Annotation, ReadOptions};
use omicsweave_quant::{
    Budget, Conversion, Counts, Labeled, Labeling, MatrixFormat, Product, Task, quantify,
};

const REFS: [(&str, u32); 2] = [("chr1", 40_000), ("chr2", 25_000)];

struct Fixture {
    dir: tempfile::TempDir,
    reads: Vec<Read>,
    genes: Vec<Gene>,
    seqs: Vec<Vec<u8>>,
    snps: BTreeSet<(usize, u32)>,
    bam: PathBuf,
    md_bam: PathBuf,
    fasta: PathBuf,
    twobit: PathBuf,
    bcf: PathBuf,
    ann: Arc<Annotation>,
}

/// A CIGAR of `len` query bases: matches with the odd insertion, deletion,
/// intron and soft clip.
fn cigar(rng: &mut Rng, len: u32) -> Vec<(u32, char)> {
    let mut out = Vec::new();
    let mut left = len;
    if rng.chance(20) {
        let s = rng.range(1, 6) as u32;
        out.push((s, 'S'));
        left -= s;
    }
    while left > 0 {
        let m = (rng.range(10, 70) as u32).min(left);
        out.push((m, 'M'));
        left -= m;
        if left == 0 {
            break;
        }
        match rng.range(0, 10) {
            0 => out.push((rng.range(1, 6) as u32, 'D')),
            1 => out.push((rng.range(50, 800) as u32, 'N')),
            2 if left > 4 => {
                let i = rng.range(1, 3) as u32;
                out.push((i, 'I'));
                left -= i;
            }
            3 if left < 8 => {
                out.push((left, 'S'));
                left = 0;
            }
            _ => {}
        }
    }
    if rng.chance(5) {
        out.insert(0, (rng.range(1, 10) as u32, 'H'));
    }
    out
}

fn ref_span(c: &[(u32, char)]) -> u32 {
    c.iter()
        .filter(|(_, o)| matches!(o, 'M' | 'D' | 'N'))
        .map(|(l, _)| l)
        .sum()
}

/// Bases of a read over `seq` at `pos` with `cigar`, on a transcript of
/// strand `reverse`: the sample's SNP alleles, T>C conversions of labeled
/// molecules (transcript sense), random errors and Ns.
#[allow(clippy::too_many_arguments)]
fn bases(
    rng: &mut Rng,
    seq: &[u8],
    tid: usize,
    pos: u32,
    cigar: &[(u32, char)],
    reverse: bool,
    labeled: bool,
    snps: &HashMap<(usize, u32), u8>,
) -> Vec<u8> {
    let mut out = Vec::new();
    let mut at = pos;
    for &(l, o) in cigar {
        match o {
            'M' => {
                for i in 0..l {
                    let p = at + i;
                    let mut b = *snps.get(&(tid, p)).unwrap_or(&seq[p as usize]);
                    let sense = if reverse { complement(b) } else { b };
                    if sense == b'T' && labeled && rng.chance(15) {
                        b = if reverse { b'G' } else { b'C' };
                    } else if rng.chance(1) {
                        b = *rng.pick(b"ACGT");
                    } else if rng.chance(1) {
                        b = b'N';
                    }
                    out.push(b);
                }
                at += l;
            }
            'I' | 'S' => out.extend((0..l).map(|_| *rng.pick(b"ACGT"))),
            'D' | 'N' => at += l,
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
    let genes = genes(&mut rng, &refs, 8);
    let mut snp_alt = HashMap::new();
    for _ in 0..300 {
        let t = rng.range(0, refs.len() as u64) as usize;
        let p = rng.range(0, refs[t].1 as u64) as u32;
        let r = seqs[t][p as usize];
        let alt = *rng.pick(
            &b"ACGT"
                .iter()
                .copied()
                .filter(|&b| b != r)
                .collect::<Vec<_>>(),
        );
        snp_alt.insert((t, p), alt);
    }
    let mut reads = Vec::new();
    for i in 0..3000 {
        let tid = rng.range(0, refs.len() as u64) as usize;
        let len = rng.range(40, 120) as u32;
        let c1 = cigar(&mut rng, len);
        let len2 = rng.range(40, 120) as u32;
        let c2 = cigar(&mut rng, len2);
        let room = ref_span(&c1).max(300 + ref_span(&c2)) + 1;
        let pos = rng.range(0, (refs[tid].1 - room) as u64) as u32;
        let transcript_reverse = rng.chance(50);
        let labeled = rng.chance(40);
        let mut tags = vec![
            format!("CB:Z:CB{}-1", rng.range(0, 8)),
            format!(
                "UB:Z:{}",
                ["AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT"][rng.range(0, 4) as usize]
            ),
            format!("NH:i:{}", if rng.chance(5) { 2 } else { 1 }),
        ];
        if rng.chance(5) {
            tags.retain(|t| !t.starts_with("UB"));
        }
        let qual = |rng: &mut Rng, n: usize| -> Option<Vec<u8>> {
            if rng.chance(5) {
                None
            } else {
                Some((0..n).map(|_| b'!' + rng.range(2, 41) as u8).collect())
            }
        };
        let dup = if rng.chance(3) { 0x400 } else { 0 };
        if rng.chance(70) {
            // Single-end, forward-stranded: the read lies on the transcript.
            let seq = bases(
                &mut rng,
                &seqs[tid],
                tid,
                pos,
                &c1,
                transcript_reverse,
                labeled,
                &snp_alt,
            );
            let q = qual(&mut rng, seq.len());
            reads.push(Read {
                name: format!("r{i}"),
                flag: dup | if transcript_reverse { 0x10 } else { 0 },
                tid,
                pos,
                mapq: 60,
                cigar: c1,
                mtid: -1,
                mpos: 0,
                tlen: 0,
                tags,
                seq: Some(seq),
                qual: q.map(|q| q.into_iter().take(len as usize).collect()),
            });
        } else {
            // A pair: read 1 on the transcript's strand, read 2 opposite.
            let p2 = pos + rng.range(0, 300) as u32;
            let s1 = bases(
                &mut rng,
                &seqs[tid],
                tid,
                pos,
                &c1,
                transcript_reverse,
                labeled,
                &snp_alt,
            );
            let s2 = bases(
                &mut rng,
                &seqs[tid],
                tid,
                p2,
                &c2,
                transcript_reverse,
                labeled,
                &snp_alt,
            );
            let (f1, f2) = if transcript_reverse {
                (0x10, 0x20)
            } else {
                (0x20, 0x10)
            };
            let end = (pos + ref_span(&c1)).max(p2 + ref_span(&c2));
            let tl = (end - pos) as i64;
            let q1 = qual(&mut rng, s1.len());
            let q2 = qual(&mut rng, s2.len());
            reads.push(Read {
                name: format!("r{i}"),
                flag: dup | 0x1 | 0x2 | 0x40 | f1,
                tid,
                pos,
                mapq: 60,
                cigar: c1,
                mtid: tid as i32,
                mpos: p2,
                tlen: tl,
                tags: tags.clone(),
                seq: Some(s1),
                qual: q1,
            });
            reads.push(Read {
                name: format!("r{i}"),
                flag: dup | 0x1 | 0x2 | 0x80 | f2,
                tid,
                pos: p2,
                mapq: 60,
                cigar: c2,
                mtid: tid as i32,
                mpos: pos,
                tlen: -tl,
                tags,
                seq: Some(s2),
                qual: q2,
            });
        }
    }
    reads.sort_by_key(|r| (r.tid, r.pos));
    let bam = write_bam(dir.path(), "in", &refs, &reads, seed);
    let with_md: Vec<Read> = reads
        .iter()
        .map(|r| {
            let mut r = r.clone();
            let md = md_tag(&r, &seqs[r.tid]);
            r.tags.push(format!("MD:Z:{md}"));
            r
        })
        .collect();
    let md_bam = write_bam(dir.path(), "md", &refs, &with_md, seed);
    let fasta = write_fasta(dir.path(), "genome", &refs, &seqs);
    let twobit = write_twobit(dir.path(), "genome", &refs, &seqs);
    let sites: Vec<(usize, u32, u8, u8)> = snp_alt
        .iter()
        .map(|(&(t, p), &a)| (t, p, seqs[t][p as usize], a))
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let bcf = write_bcf(dir.path(), "snps", &refs, &sites);
    let g = dir.path().join("genes.gtf");
    std::fs::write(&g, gtf(&refs, &genes)).unwrap();
    let ann = Arc::new(Annotation::read(&g, &ReadOptions::default()).unwrap());
    Fixture {
        dir,
        reads,
        genes,
        seqs,
        snps: snp_alt.keys().copied().collect(),
        bam,
        md_bam,
        fasta,
        twobit,
        bcf,
        ann,
    }
}

const MIN_QUALITY: u8 = 20;

/// Per transcript strand (forward, reverse): coverage and conversions of
/// a read, as the module defines them.
fn conversions(f: &Fixture, r: &Read, masked: bool) -> [(u32, u32); 2] {
    let seq = r.seq.as_ref().unwrap();
    let mut out = [(0, 0); 2];
    for (p, q) in aligned_pairs(r) {
        let g = f.seqs[r.tid][p as usize];
        let b = seq[q];
        let ok_qual = r.qual.as_ref().is_none_or(|qs| qs[q] - 33 >= MIN_QUALITY);
        if !ok_qual || b == b'N' || (masked && f.snps.contains(&(r.tid, p))) {
            continue;
        }
        for (s, (from, to)) in [(b'T', b'C'), (b'A', b'G')].into_iter().enumerate() {
            if g == from {
                out[s].0 += 1;
                out[s].1 += (b == to) as u32;
            }
        }
    }
    out
}

fn labeling_filter(r: &Read) -> bool {
    r.flag & 0xf04 == 0 && r.nh() == 1
}

/// `conversions.tsv.gz` and `mismatches.tsv` computed from the reads.
fn expected_tables(f: &Fixture, masked: bool) -> (String, String) {
    let mut hist: BTreeMap<(u32, u32), u64> = BTreeMap::new();
    // (mate, cycle, ref) → read base counts.
    let mut cycles: BTreeMap<(u8, u32, u8), [u64; 4]> = BTreeMap::new();
    for r in &f.reads {
        if !labeling_filter(r) {
            continue;
        }
        let rev = reverse_for(r, 1);
        let c = conversions(f, r, masked)[rev as usize];
        *hist.entry(c).or_insert(0) += 1;
        let seq = r.seq.as_ref().unwrap();
        let lead = match r.cigar.first() {
            Some(&(l, 'H')) => l,
            _ => 0,
        };
        let tail = match r.cigar.last() {
            Some(&(l, 'H')) if r.cigar.len() > 1 => l,
            _ => 0,
        };
        for (p, q) in aligned_pairs(r) {
            let g = f.seqs[r.tid][p as usize];
            let b = seq[q];
            let ok_qual = r.qual.as_ref().is_none_or(|qs| qs[q] - 33 >= MIN_QUALITY);
            if !ok_qual || b == b'N' || (masked && f.snps.contains(&(r.tid, p))) {
                continue;
            }
            let (g, b) = if rev {
                (complement(g), complement(b))
            } else {
                (g, b)
            };
            let cycle = if r.reverse() {
                tail + (seq.len() - 1 - q) as u32
            } else {
                lead + q as u32
            };
            let mate = if r.flag & 0x81 == 0x81 { 2 } else { 1 };
            let row = cycles.entry((mate, cycle, g)).or_insert([0; 4]);
            row[b"ACGT".iter().position(|&x| x == b).unwrap()] += 1;
        }
    }
    let mut h = "T_coverage\tT>C_conversions\treads\n".to_string();
    for ((a, b), n) in hist {
        h.push_str(&format!("{a}\t{b}\t{n}\n"));
    }
    let mut m = "mate\tcycle\treference\tA\tC\tG\tT\n".to_string();
    // Rows go by mate, cycle, then reference base in ACGT order.
    for ((mate, cycle, g), v) in cycles {
        m.push_str(&format!(
            "{mate}\t{cycle}\t{}\t{}\t{}\t{}\t{}\n",
            g as char, v[0], v[1], v[2], v[3]
        ));
    }
    (h, m)
}

fn task(f: &Fixture, bam: &Path, out: &Path, reference: Option<&Path>, mask: bool) -> Task {
    let mut sc = Counts::rna_single_cell();
    sc.labeled = Some(Labeled {
        conversion: Conversion::default(),
        min_conversions: 2,
    });
    let mut bulk = Counts::bulk_rna();
    bulk.name = "bulk".into();
    bulk.labeled = Some(Labeled::default());
    let mut t = Task::new(vec![bam.to_owned()], out);
    t.annotation = Some(Arc::clone(&f.ann));
    t.reference = reference.map(Path::to_owned);
    t.variants = mask.then(|| f.bcf.clone());
    t.matrix = MatrixFormat::Mtx;
    t.products = vec![
        Product::Labeling(Labeling::default()),
        Product::Counts(sc),
        Product::Counts(bulk),
    ];
    t
}

fn product_files(dir: &Path) -> Vec<(String, Vec<u8>)> {
    let mut out = Vec::new();
    for sub in ["labeling", "genes", "bulk"] {
        let mut names: Vec<_> = std::fs::read_dir(dir.join(sub))
            .unwrap()
            .map(|e| e.unwrap().path())
            .collect();
        names.sort();
        for p in names {
            out.push((
                p.strip_prefix(dir).unwrap().display().to_string(),
                std::fs::read(&p).unwrap(),
            ));
        }
    }
    out
}

/// The labeled layers computed from the reads: single-cell molecules
/// (barcode, gene, UMI) with a read of two or more conversions, and bulk
/// reads with one or more on the strand of their gene.
type Labels = (
    BTreeMap<(String, String), u32>,
    BTreeMap<(String, String), u32>,
);

fn expected_labels(f: &Fixture, bam: &Path, masked: bool) -> Labels {
    let mut molecules = BTreeSet::new();
    let mut bulk = BTreeMap::new();
    for r in &f.reads {
        if r.nh() > 1 {
            continue;
        }
        let c = conversions(f, r, masked);
        if r.flag & 0x704 == 0
            && let (Some(cb), Some(ub)) = (r.tag("CB"), r.tag("UB"))
        {
            let rev = reverse_for(r, 1);
            let hits = gene_hits(&f.genes, r.tid, &r.blocks(false), Some(rev));
            if hits.len() == 1 && c[rev as usize].1 >= 2 {
                molecules.insert((cb.to_owned(), f.genes[hits[0]].id.clone(), ub.to_owned()));
            }
        }
        if r.flag & 0x4 == 0 {
            let hits = gene_hits(&f.genes, r.tid, &r.blocks(false), None);
            if hits.len() == 1 {
                let g = &f.genes[hits[0]];
                if c[(g.strand == '-') as usize].1 >= 1 {
                    *bulk
                        .entry((bam.display().to_string(), g.id.clone()))
                        .or_insert(0) += 1;
                }
            }
        }
    }
    let mut sc = BTreeMap::new();
    for (cb, g, _) in molecules {
        *sc.entry((cb, g)).or_insert(0) += 1;
    }
    (sc, bulk)
}

#[test]
fn conversions_and_mismatches_match_the_reads() {
    let f = fixture(51);
    for masked in [true, false] {
        let out = f.dir.path().join(format!("fa{masked}"));
        quantify(
            task(&f, &f.bam, &out, Some(&f.fasta), masked),
            Budget::new(3, 256 << 20),
            None,
        )
        .unwrap();
        let (hist, mism) = expected_tables(&f, masked);
        assert_eq!(gunzip(&out.join("labeling/conversions.tsv.gz")), hist);
        assert_eq!(
            std::fs::read_to_string(out.join("labeling/mismatches.tsv")).unwrap(),
            mism
        );
        let (sc, bulk) = expected_labels(&f, &f.bam, masked);
        assert!(!sc.is_empty() && !bulk.is_empty());
        assert_eq!(read_10x_layer(&out.join("genes"), "labeled"), sc);
        assert_eq!(read_10x_layer(&out.join("bulk"), "labeled"), bulk);
    }
}

#[test]
fn twobit_md_tags_and_threads_give_the_same_files() {
    let f = fixture(52);
    let reference = f.dir.path().join("ref");
    quantify(
        task(&f, &f.bam, &reference, Some(&f.fasta), true),
        Budget::new(1, 256 << 20),
        None,
    )
    .unwrap();
    let want = product_files(&reference);
    let twobit = f.dir.path().join("twobit");
    quantify(
        task(&f, &f.bam, &twobit, Some(&f.twobit), true),
        Budget::new(6, 256 << 20),
        None,
    )
    .unwrap();
    assert_eq!(product_files(&twobit), want);
    assert_eq!(
        std::fs::read(twobit.join("summary.json")).unwrap(),
        std::fs::read(reference.join("summary.json")).unwrap()
    );
    let md = f.dir.path().join("md");
    quantify(
        task(&f, &f.md_bam, &md, None, true),
        Budget::new(4, 256 << 20),
        None,
    )
    .unwrap();
    // The bulk rows name the input file, which differs here.
    let rename = |v: Vec<(String, Vec<u8>)>| -> Vec<(String, Vec<u8>)> {
        v.into_iter()
            .filter(|(n, _)| n != "bulk/barcodes.tsv.gz" && n != "bulk/counts.tsv")
            .collect()
    };
    let got = rename(product_files(&md));
    let want_md = rename(want.clone());
    let differ: Vec<&String> = got
        .iter()
        .zip(&want_md)
        .filter(|(a, b)| a != b)
        .map(|(a, _)| &a.0)
        .collect();
    assert!(
        differ.is_empty() && got.len() == want_md.len(),
        "{differ:?}"
    );
    let (_, bulk) = expected_labels(&f, &f.md_bam, true);
    assert_eq!(read_10x_layer(&md.join("bulk"), "labeled"), bulk);
    // Without a reference or MD tags, nothing is labeled.
    let none = f.dir.path().join("none");
    quantify(
        task(&f, &f.bam, &none, None, false),
        Budget::new(2, 256 << 20),
        None,
    )
    .unwrap();
    assert!(read_10x_layer(&none.join("genes"), "labeled").is_empty());
    let summary = std::fs::read_to_string(none.join("summary.json")).unwrap();
    assert!(summary.contains("\"unassigned_no_reference\""));
}
