//! Splicing products against direct computation on generated reads: the
//! junction table (STAR `SJ.out.tab` columns, motifs from a bgzip FASTA
//! and a 2bit genome), the junction × barcode matrix, and the velocity
//! layers of gene matrices; outputs independent of threads and spilling.

mod common;

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use common::{
    Gene, Read, Rng, Spec, genes, genome, gtf, gunzip, read_10x, read_10x_layer, reverse_for,
    synth, write_bam, write_fasta, write_twobit,
};
use omicsweave_gtf::{Annotation, ReadOptions};
use omicsweave_quant::{
    Budget, Counts, JunctionMatrix, Junctions, MatrixFormat, Product, Task, plan, quantify, run,
};

struct Fixture {
    dir: tempfile::TempDir,
    spec: Spec,
    reads: Vec<Read>,
    genes: Vec<Gene>,
    seqs: Vec<Vec<u8>>,
    bam: PathBuf,
    fasta: PathBuf,
    twobit: PathBuf,
    ann: Arc<Annotation>,
}

/// The merged exon blocks of each transcript.
fn transcripts(g: &Gene) -> Vec<Vec<(u32, u32)>> {
    g.transcripts
        .iter()
        .map(|ex| {
            let mut ex = ex.clone();
            ex.sort();
            let mut out: Vec<(u32, u32)> = Vec::new();
            for (s, e) in ex {
                match out.last_mut() {
                    Some(l) if s <= l.1 => l.1 = l.1.max(e),
                    _ => out.push((s, e)),
                }
            }
            out
        })
        .collect()
}

/// Reads spliced exactly at annotated introns: single reads, and pairs
/// whose mates both cross the same intron.
fn on_introns(rng: &mut Rng, genes: &[Gene], n: usize) -> Vec<Read> {
    let mut out = Vec::new();
    let mut serial = 0;
    while out.len() < n {
        let g = rng.pick(genes);
        let t = rng.pick(&transcripts(g)).clone();
        if t.len() < 2 {
            continue;
        }
        let i = rng.range(0, t.len() as u64 - 1) as usize;
        let (donor, acceptor) = (t[i].1, t[i + 1].0);
        let left = (t[i].1 - t[i].0).min(rng.range(8, 50) as u32);
        let right = (t[i + 1].1 - t[i + 1].0).min(rng.range(8, 50) as u32);
        if left < 3 || right < 3 {
            continue;
        }
        serial += 1;
        let tags = vec![
            format!("CB:Z:AAAACCCC-{}", rng.range(0, 3)),
            format!(
                "UB:Z:{}",
                ["ACGTACGTAA", "TTGGCCAATT"][rng.range(0, 2) as usize]
            ),
            format!("NH:i:{}", if rng.chance(10) { 2 } else { 1 }),
        ];
        let intron = acceptor - donor;
        let mk = |name: String, flag: u16, pos: u32, l: u32, r: u32| Read {
            name,
            flag,
            tid: g.tid,
            pos,
            mapq: 60,
            cigar: vec![(l, 'M'), (intron, 'N'), (r, 'M')],
            mtid: -1,
            mpos: 0,
            tlen: 0,
            tags: tags.clone(),
            seq: None,
            qual: None,
        };
        if rng.chance(50) {
            let flag = if rng.chance(50) { 16 } else { 0 };
            out.push(mk(format!("j{serial}"), flag, donor - left, left, right));
        } else {
            let l2 = rng.range(2, left as u64) as u32;
            let r2 = right.min(rng.range(3, 40) as u32);
            let (p1, p2) = (donor - left, donor - l2);
            let end = (acceptor + right).max(acceptor + r2);
            let mut a = mk(
                format!("j{serial}"),
                0x1 | 0x2 | 0x20 | 0x40,
                p1,
                left,
                right,
            );
            let mut b = mk(format!("j{serial}"), 0x1 | 0x2 | 0x10 | 0x80, p2, l2, r2);
            a.mtid = g.tid as i32;
            a.mpos = p2;
            a.tlen = (end - p1) as i64;
            b.mtid = g.tid as i32;
            b.mpos = p1;
            b.tlen = -((end - p1) as i64);
            out.push(a);
            out.push(b);
        }
    }
    out
}

/// The junctions of a read: `(start, end, overhang)`.
fn junctions_of(r: &Read) -> Vec<(u32, u32, u32)> {
    let mut segs = Vec::new();
    let mut found = Vec::new();
    let (mut run, mut at) = (0u32, r.pos);
    for &(l, o) in &r.cigar {
        match o {
            'M' | '=' | 'X' => {
                run += l;
                at += l;
            }
            'N' => {
                segs.push(run);
                found.push((at, at + l, segs.len() - 1));
                run = 0;
                at += l;
            }
            'D' => {
                segs.push(run);
                run = 0;
                at += l;
            }
            'I' | 'S' => {
                segs.push(run);
                run = 0;
            }
            _ => {}
        }
    }
    segs.push(run);
    found
        .into_iter()
        .map(|(s, e, i)| (s, e, segs[i].min(segs[i + 1])))
        .collect()
}

fn fixture(seed: u64) -> Fixture {
    let dir = tempfile::tempdir().unwrap();
    let spec = Spec {
        fragments: 4000,
        ..Spec::default()
    };
    let mut rng = Rng::new(seed);
    let mut reads = synth(&mut rng, &spec);
    let genes = genes(&mut rng, &spec.refs, 10);
    reads.extend(on_introns(&mut rng, &genes, 600));
    reads.sort_by_key(|r| (r.tid, r.pos));
    let mut seqs = genome(&mut rng, &spec.refs);
    // Canonical motifs at most junctions, so that every motif is seen.
    let motifs: [(&[u8; 2], &[u8; 2]); 6] = [
        (b"GT", b"AG"),
        (b"CT", b"AC"),
        (b"GC", b"AG"),
        (b"CT", b"GC"),
        (b"AT", b"AC"),
        (b"GT", b"AT"),
    ];
    let mut all = BTreeSet::new();
    for r in &reads {
        for (s, e, _) in junctions_of(r) {
            all.insert((r.tid, s, e));
        }
    }
    for (t, s, e) in all {
        if e - s >= 4 && rng.chance(70) {
            let (d, a) = rng.pick(&motifs);
            seqs[t][s as usize..s as usize + 2].copy_from_slice(*d);
            seqs[t][e as usize - 2..e as usize].copy_from_slice(*a);
        }
    }
    let bam = write_bam(dir.path(), "in", &spec.refs, &reads, seed);
    let fasta = write_fasta(dir.path(), "genome", &spec.refs, &seqs);
    let twobit = write_twobit(dir.path(), "genome", &spec.refs, &seqs);
    let g = dir.path().join("genes.gtf");
    std::fs::write(&g, gtf(&spec.refs, &genes)).unwrap();
    let ann = Arc::new(Annotation::read(&g, &ReadOptions::default()).unwrap());
    Fixture {
        dir,
        spec,
        reads,
        genes,
        seqs,
        bam,
        fasta,
        twobit,
        ann,
    }
}

fn task(f: &Fixture, out: &Path, products: Vec<Product>, reference: Option<&Path>) -> Task {
    let mut t = Task::new(vec![f.bam.clone()], out);
    t.annotation = Some(Arc::clone(&f.ann));
    t.reference = reference.map(Path::to_owned);
    t.products = products;
    t
}

/// A mate pair: its name with both positions.
type PairKey = (String, (usize, u32), (i32, u32));

/// The mate pair a record belongs to.
fn pair_key(r: &Read) -> PairKey {
    let me = (r.tid, r.pos);
    let mate = (r.mtid, r.mpos);
    if (me.0 as i32, me.1) <= mate {
        (r.name.clone(), me, mate)
    } else {
        (
            r.name.clone(),
            (r.mtid as usize, r.mpos),
            (r.tid as i32, r.pos),
        )
    }
}

fn expects_mate(r: &Read) -> bool {
    r.flag & 1 != 0 && r.flag & 0xc == 0 && r.mtid >= 0
}

/// A counted unit: its records, in record order.
fn units(reads: &[Read], keep: impl Fn(&Read) -> bool) -> Vec<Vec<&Read>> {
    let mut out: Vec<Vec<&Read>> = Vec::new();
    let mut open: HashMap<PairKey, usize> = HashMap::new();
    for r in reads {
        if !keep(r) {
            continue;
        }
        if expects_mate(r) {
            let k = pair_key(r);
            match open.remove(&k) {
                Some(i) => out[i].push(r),
                None => {
                    open.insert(k, out.len());
                    out.push(vec![r]);
                }
            }
        } else {
            out.push(vec![r]);
        }
    }
    out
}

/// The union of a unit's junctions, largest overhang kept.
fn unit_junctions(u: &[&Read]) -> BTreeMap<(usize, u32, u32), u32> {
    let mut m = BTreeMap::new();
    for r in u {
        for (s, e, o) in junctions_of(r) {
            let x = m.entry((r.tid, s, e)).or_insert(0);
            *x = (*x).max(o);
        }
    }
    m
}

fn table_filter(r: &Read) -> bool {
    r.flag & (0x4 | 0x200 | 0x800) == 0 && !junctions_of(r).is_empty()
}

/// `SJ.out.tab` lines computed from the reads, the genes and the genome.
fn expected_table(f: &Fixture, with_motifs: bool) -> Vec<String> {
    let mut rows: BTreeMap<(usize, u32, u32), (u32, u32, u32)> = BTreeMap::new();
    for u in units(&f.reads, table_filter) {
        let multi = u.iter().any(|r| r.nh() > 1);
        for (k, o) in unit_junctions(&u) {
            let x = rows.entry(k).or_insert((0, 0, 0));
            if multi {
                x.1 += 1;
            } else {
                x.0 += 1;
            }
            x.2 = x.2.max(o);
        }
    }
    let mut introns: HashMap<(usize, u32, u32), char> = HashMap::new();
    for g in &f.genes {
        for t in transcripts(g) {
            for w in t.windows(2) {
                let e = introns.entry((g.tid, w[0].1, w[1].0)).or_insert(g.strand);
                if g.strand == '+' {
                    *e = '+';
                }
            }
        }
    }
    rows.into_iter()
        .map(|((t, s, e), (uq, mu, o))| {
            let ann = introns.get(&(t, s, e));
            let motif = if with_motifs && e - s >= 4 {
                let seq = &f.seqs[t];
                let d = &seq[s as usize..s as usize + 2];
                let a = &seq[e as usize - 2..e as usize];
                match (d, a) {
                    (b"GT", b"AG") => 1,
                    (b"CT", b"AC") => 2,
                    (b"GC", b"AG") => 3,
                    (b"CT", b"GC") => 4,
                    (b"AT", b"AC") => 5,
                    (b"GT", b"AT") => 6,
                    _ => 0,
                }
            } else {
                0
            };
            let strand = match (motif, ann) {
                (0, Some('+')) => 1,
                (0, Some(_)) => 2,
                (0, None) => 0,
                (m, _) if m % 2 == 1 => 1,
                _ => 2,
            };
            format!(
                "{}\t{}\t{}\t{strand}\t{motif}\t{}\t{uq}\t{mu}\t{o}",
                f.spec.refs[t].0,
                s + 1,
                e,
                ann.is_some() as u8
            )
        })
        .collect()
}

fn lines(p: &Path) -> Vec<String> {
    std::fs::read_to_string(p)
        .unwrap()
        .lines()
        .map(str::to_owned)
        .collect()
}

#[test]
fn junction_table_matches_star_columns() {
    let f = fixture(41);
    let want = expected_table(&f, true);
    assert!(want.iter().any(|l| l.split('\t').nth(5) == Some("1")));
    for m in 1..=6 {
        assert!(
            want.iter()
                .any(|l| l.split('\t').nth(4) == Some(&m.to_string())),
            "motif {m} seen"
        );
    }
    let mut outputs = Vec::new();
    for (i, reference) in [&f.fasta, &f.twobit].into_iter().enumerate() {
        let out = f.dir.path().join(format!("j{i}"));
        let t = task(
            &f,
            &out,
            vec![Product::Junctions(Junctions::default())],
            Some(reference),
        );
        quantify(t, Budget::new(1 + 4 * i, 256 << 20), None).unwrap();
        let got = lines(&out.join("junctions/junctions.tab"));
        assert_eq!(got, want, "reference {}", reference.display());
        outputs.push(std::fs::read(out.join("junctions/junctions.tab")).unwrap());
    }
    assert_eq!(outputs[0], outputs[1]);
    // Without a reference: no motifs, strands from the annotation.
    let out = f.dir.path().join("nomotif");
    let t = task(
        &f,
        &out,
        vec![Product::Junctions(Junctions::default())],
        None,
    );
    quantify(t, Budget::new(3, 256 << 20), None).unwrap();
    assert_eq!(
        lines(&out.join("junctions/junctions.tab")),
        expected_table(&f, false)
    );
}

#[test]
fn junction_matrix_counts_umis_per_barcode() {
    let f = fixture(42);
    let out = f.dir.path().join("m");
    let j = Junctions {
        matrix: Some(JunctionMatrix::single_cell()),
        ..Junctions::default()
    };
    let mut t = task(&f, &out, vec![Product::Junctions(j)], Some(&f.fasta));
    t.matrix = MatrixFormat::Mtx;
    quantify(t, Budget::new(4, 256 << 20), None).unwrap();

    let matrix_filter = |r: &Read| r.flag & 0x704 == 0 && r.nh() == 1;
    // Barcodes in order of first appearance among reads counted in the
    // matrix; units take the row and UMI of their first such read.
    let mut rows: Vec<String> = Vec::new();
    for r in &f.reads {
        if table_filter(r)
            && matrix_filter(r)
            && let Some(cb) = r.tag("CB")
            && !rows.iter().any(|x| x == cb)
        {
            rows.push(cb.to_owned());
        }
    }
    let mut seen = BTreeSet::new();
    for u in units(&f.reads, table_filter) {
        let Some(first) = u
            .iter()
            .find(|r| matrix_filter(r) && r.tag("CB").is_some() && r.tag("UB").is_some())
        else {
            continue;
        };
        let (cb, ub) = (first.tag("CB").unwrap(), first.tag("UB").unwrap());
        for (t, s, e) in unit_junctions(&u).into_keys() {
            seen.insert((
                cb.to_owned(),
                format!("{}:{s}-{e}", f.spec.refs[t].0),
                ub.to_owned(),
            ));
        }
    }
    let mut want = BTreeMap::new();
    for (cb, j, _) in seen {
        *want.entry((cb, j)).or_insert(0u32) += 1;
    }
    let dir = out.join("junctions");
    assert_eq!(read_10x(&dir), want);
    let barcodes: Vec<String> = gunzip(&dir.join("barcodes.tsv.gz"))
        .lines()
        .map(str::to_owned)
        .collect();
    assert_eq!(barcodes, rows);
    // One column per junction of the table, in its order.
    let table = lines(&dir.join("junctions.tab"));
    let features = gunzip(&dir.join("features.tsv.gz"));
    assert_eq!(features.lines().count(), table.len());
    for (l, t) in features.lines().zip(&table) {
        let fv: Vec<&str> = l.split('\t').collect();
        let tv: Vec<&str> = t.split('\t').collect();
        assert_eq!(fv[2], "Splice Junction");
        assert_eq!(fv[3], tv[0]);
        assert_eq!(
            fv[4].parse::<u32>().unwrap() + 1,
            tv[1].parse::<u32>().unwrap()
        );
        assert_eq!(fv[5], tv[2]);
    }
}

/// velocyto's rules on the generated reads: per read, exonic or intronic
/// against each transcript of its one gene holding it; per molecule, the
/// OR of its reads.
fn classify(g: &Gene, blocks: &[(u32, u32)]) -> u32 {
    let mut bits = 0;
    for t in transcripts(g) {
        let (ts, te) = (t[0].0, t.last().unwrap().1);
        if blocks[0].0 < ts || blocks.last().unwrap().1 > te {
            continue;
        }
        let exonic = blocks
            .iter()
            .all(|&(s, e)| t.iter().any(|&(a, b)| a <= s && e <= b));
        bits |= if exonic { 1 } else { 2 };
    }
    bits
}

type Layers = [BTreeMap<(String, String), u32>; 3];

fn expected_velocity(f: &Fixture, single_cell: bool) -> Layers {
    let mut molecules: BTreeMap<(String, String, String), u32> = BTreeMap::new();
    for (i, r) in f.reads.iter().enumerate() {
        if r.nh() > 1 {
            continue;
        }
        let (row, sub) = if single_cell {
            if r.flag & 0x704 != 0 {
                continue;
            }
            let (Some(cb), Some(ub)) = (r.tag("CB"), r.tag("UB")) else {
                continue;
            };
            (cb.to_owned(), ub.to_owned())
        } else {
            if r.flag & 0x4 != 0 {
                continue;
            }
            (f.bam.display().to_string(), i.to_string())
        };
        let blocks = r.blocks(false);
        let strand = single_cell.then(|| reverse_for(r, 1));
        let hits: Vec<&Gene> = f
            .genes
            .iter()
            .filter(|g| g.tid == r.tid)
            .filter(|g| strand.is_none_or(|rev| (g.strand == '-') == rev))
            .filter(|g| blocks.iter().any(|&(s, e)| s < g.end && g.start < e))
            .collect();
        if hits.len() != 1 {
            continue;
        }
        let bits = classify(hits[0], &blocks);
        if bits > 0 {
            *molecules.entry((row, hits[0].id.clone(), sub)).or_insert(0) |= bits;
        }
    }
    let mut out: Layers = Default::default();
    for ((row, gene, _), bits) in molecules {
        let layer = match bits {
            1 => 0,
            2 => 1,
            _ => 2,
        };
        *out[layer].entry((row, gene)).or_insert(0) += 1;
    }
    out
}

#[test]
fn velocity_layers_follow_velocyto_rules() {
    let f = fixture(43);
    let out = f.dir.path().join("v");
    let mut sc = Counts::rna_single_cell();
    sc.velocity = true;
    let mut bulk = Counts::bulk_rna();
    bulk.name = "bulk".into();
    bulk.velocity = true;
    bulk.table = false;
    let t = task(
        &f,
        &out,
        vec![Product::Counts(sc), Product::Counts(bulk)],
        None,
    );
    quantify(t, Budget::new(5, 256 << 20), None).unwrap();
    for (name, single_cell) in [("genes", true), ("bulk", false)] {
        let want = expected_velocity(&f, single_cell);
        assert!(
            want.iter().all(|l| !l.is_empty()),
            "{name}: every class seen"
        );
        for (i, layer) in ["spliced", "unspliced", "ambiguous"].iter().enumerate() {
            assert_eq!(
                read_10x_layer(&out.join(name), layer),
                want[i],
                "{name}/{layer}"
            );
        }
    }
}

fn files(dir: &Path) -> Vec<(String, Vec<u8>)> {
    let mut out = Vec::new();
    let mut stack = vec![dir.to_owned()];
    while let Some(d) = stack.pop() {
        for e in std::fs::read_dir(&d).unwrap() {
            let p = e.unwrap().path();
            if p.is_dir() {
                stack.push(p);
            } else {
                out.push((
                    p.strip_prefix(dir).unwrap().display().to_string(),
                    std::fs::read(&p).unwrap(),
                ));
            }
        }
    }
    out.sort();
    out
}

#[test]
fn splicing_outputs_do_not_depend_on_threads_or_spilling() {
    let f = fixture(44);
    let products = || {
        let j = Junctions {
            matrix: Some(JunctionMatrix::single_cell()),
            ..Junctions::default()
        };
        let mut sc = Counts::rna_single_cell();
        sc.velocity = true;
        vec![Product::Junctions(j), Product::Counts(sc)]
    };
    let a = f.dir.path().join("a");
    quantify(
        task(&f, &a, products(), Some(&f.fasta)),
        Budget::new(1, 1 << 30),
        None,
    )
    .unwrap();
    let b = f.dir.path().join("b");
    let mut p = plan(
        task(&f, &b, products(), Some(&f.fasta)),
        Budget::new(7, 256 << 20),
    )
    .unwrap();
    for s in &mut p.stages {
        s.1 = 16 << 10;
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
}
