//! DNA methylation against direct computation on bisulfite read pairs drawn
//! from a generated genome (directional library, both original strands,
//! molecules methylated per cytosine at random, sequencing errors, `N`s,
//! indels and clips, overlapping mates, filtered records, known SNPs) and on
//! long reads with `MM` / `ML` tags: calls from the reference, from `MD`
//! tags and from Bismark `XM` tags; the per-cytosine table (with the tabix
//! index htslib builds for it), the beta and coverage tracks, M-bias,
//! region × barcode and × file matrices; the same files at any thread count
//! and with spilled buffers.

mod common;

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

use common::{
    Read, Rng, assert_tabix_as_htslib, bigwig_values, complement, files, genome, gunzip, md_tag,
    read_10x_layer, write_bam, write_bcf, write_fasta,
};
use omicsweave_quant::{
    Budget, CHG, CHH, CPG, CallSource, FeatureSource, MethylMatrix, Methylation, Product, RowKey,
    Task, plan, quantify, run,
};

const REFS: [(&str, u32); 2] = [("chr1", 30_000), ("chr2", 18_000)];
const CONTEXTS: [&str; 3] = ["CpG", "CHG", "CHH"];

#[test]
fn several_input_matrices_keep_file_rows_and_output_restrictions() {
    let dir = tempfile::tempdir().unwrap();
    let mut inputs = Vec::new();
    for name in ["one", "two"] {
        let sam = dir.path().join(format!("{name}.sam"));
        std::fs::write(&sam,"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nr\t0\tchr1\t101\t60\t3M\t*\t0\t0\tCGC\tIII\tXM:Z:Z.z\tXG:Z:CT\n").unwrap();
        let bam = dir.path().join(format!("{name}.bam"));
        common::bam_of(&sam, &bam);
        inputs.push(bam);
    }
    let out = dir.path().join("out");
    let mut task = Task::new(inputs, &out);
    task.matrix = omicsweave_quant::MatrixFormat::Mtx;
    let spec = Methylation {
        source: CallSource::Bismark,
        table: false,
        bigwig: false,
        overlap: false,
        matrices: vec![MethylMatrix::new("bins", FeatureSource::Bins(1000))],
        ..Methylation::default()
    };
    task.products = vec![Product::Methylation(spec.clone())];
    quantify(task.clone(), Budget::new(2, 128 << 20), None).unwrap();
    let path = out.join("methylation/bins");
    for (layer, want) in [("methylated", 1), ("total", 2)] {
        let matrix = read_10x_layer(&path, layer);
        assert_eq!(matrix.len(), 2);
        assert!(matrix.values().all(|&n| n == want));
        assert_eq!(
            matrix.keys().map(|(r, _)| r).collect::<BTreeSet<_>>().len(),
            2
        );
    }
    let mut forbidden = spec;
    forbidden.bigwig = true;
    task.products = vec![Product::Methylation(forbidden)];
    assert!(
        plan(task, Budget::new(1, 128 << 20))
            .unwrap_err()
            .to_string()
            .contains("one BAM")
    );
}

struct Fixture {
    dir: tempfile::TempDir,
    refs: Vec<(String, u32)>,
    seqs: Vec<Vec<u8>>,
    reads: Vec<Read>,
    snps: BTreeSet<(usize, u32)>,
    bam: PathBuf,
    fasta: PathBuf,
    bcf: PathBuf,
}

/// Matches with the odd insertion, deletion and soft clip, and now and
/// then a hard clip in front; `len` query bases.
fn cigar(rng: &mut Rng, len: u32) -> Vec<(u32, char)> {
    let mut out = Vec::new();
    if rng.chance(5) {
        out.push((rng.range(1, 8) as u32, 'H'));
    }
    let mut left = len;
    if rng.chance(15) {
        let s = rng.range(1, 6) as u32;
        out.push((s, 'S'));
        left -= s;
    }
    while left > 0 {
        let m = (rng.range(10, 60) as u32).min(left);
        out.push((m, 'M'));
        left -= m;
        if left == 0 {
            break;
        }
        match rng.range(0, 8) {
            0 => out.push((rng.range(1, 5) as u32, 'D')),
            1 if left > 4 => {
                let i = rng.range(1, 3) as u32;
                out.push((i, 'I'));
                left -= i;
            }
            2 if left < 8 => {
                out.push((left, 'S'));
                left = 0;
            }
            _ => {}
        }
    }
    out
}

fn ref_span(c: &[(u32, char)]) -> u32 {
    c.iter()
        .filter(|(_, o)| matches!(o, 'M' | 'D' | 'N'))
        .map(|(l, _)| l)
        .sum()
}

/// `(ref_pos, query_pos)` of every `M` base (hard clips ignored).
fn aligned(r: &Read) -> Vec<(u32, usize)> {
    common::aligned_pairs(r)
}

/// The context index of a cytosine followed by `n1`, `n2` on its strand.
fn context(n1: u8, n2: u8) -> Option<u8> {
    let h = |b: u8| b == b'A' || b == b'C' || b == b'T';
    if n1 == b'G' {
        Some(0)
    } else if h(n1) && n2 == b'G' {
        Some(1)
    } else if h(n1) && h(n2) {
        Some(2)
    } else {
        None
    }
}

fn base_at(seq: &[u8], p: i64) -> u8 {
    if p < 0 || p as usize >= seq.len() {
        b'N'
    } else {
        seq[p as usize]
    }
}

/// The genome context of the cytosine at `p` on the top or bottom strand.
fn genome_context(seq: &[u8], p: u32, top: bool) -> Option<u8> {
    let p = p as i64;
    if top {
        context(base_at(seq, p + 1), base_at(seq, p + 2))
    } else {
        context(
            complement(base_at(seq, p - 1)),
            complement(base_at(seq, p - 2)),
        )
    }
}

/// Bases of a read of the top (or bottom) strand at `pos` with `cigar`:
/// unmethylated cytosines of the molecule converted, random errors, `N`s.
#[allow(clippy::too_many_arguments)]
fn bases(
    rng: &mut Rng,
    seq: &[u8],
    pos: u32,
    cigar: &[(u32, char)],
    top: bool,
    state: &mut HashMap<u32, bool>,
    tid: usize,
    rng_state: &mut Rng,
) -> Vec<u8> {
    let mut out = Vec::new();
    let mut at = pos;
    for &(l, o) in cigar {
        match o {
            'M' => {
                for i in 0..l {
                    let p = at + i;
                    let g = seq[p as usize];
                    let mut b = g;
                    let target = if top { b'C' } else { b'G' };
                    if g == target {
                        let meth = *state
                            .entry(p)
                            .or_insert_with(|| rng_state.chance(level(tid, p, top)));
                        if !meth {
                            b = if top { b'T' } else { b'A' };
                        }
                    }
                    if rng.chance(1) {
                        b = *rng.pick(b"ACGT");
                    } else if rng.next().is_multiple_of(200) {
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

/// The methylation level (percent) of a cytosine.
fn level(tid: usize, p: u32, top: bool) -> u64 {
    let h =
        (tid as u64 * 1_000_003 + p as u64 * 2 + top as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    [0, 15, 85, 100][(h >> 60) as usize % 4]
}

/// Bismark's `XM` for a read: each aligned cytosine of its strand by
/// genome context (`U` unknown), uppercase when read unconverted.
fn xm_tag(r: &Read, seq: &[u8], top: bool) -> String {
    let bases = r.seq.as_ref().unwrap();
    let mut xm = vec![b'.'; bases.len()];
    for (p, q) in aligned(r) {
        let (target, conv) = if top { (b'C', b'T') } else { (b'G', b'A') };
        if seq[p as usize] != target {
            continue;
        }
        let letter = match genome_context(seq, p, top) {
            Some(c) => b"zxh"[c as usize],
            None => b'u',
        };
        if bases[q] == target {
            xm[q] = letter.to_ascii_uppercase();
        } else if bases[q] == conv {
            xm[q] = letter;
        }
    }
    String::from_utf8(xm).unwrap()
}

fn fixture(seed: u64) -> Fixture {
    let dir = tempfile::tempdir().unwrap();
    let refs: Vec<(String, u32)> = REFS.iter().map(|(n, l)| (n.to_string(), *l)).collect();
    let mut rng = Rng::new(seed);
    let seqs = genome(&mut rng, &refs);
    let mut states = Rng::new(seed ^ 0xabc);
    let mut reads = Vec::new();
    for i in 0..2600 {
        let tid = rng.range(0, 2) as usize;
        let len = refs[tid].1;
        let flen = rng.range(120, 420) as u32;
        let start = rng.range(0, (len - flen - 1) as u64) as u32;
        let top = rng.chance(50);
        let mut state = HashMap::new();
        let mut mk = |rng: &mut Rng, left: bool| -> (u32, Vec<(u32, char)>, Vec<u8>) {
            let n = rng.range(50, 130) as u32;
            let mut c = cigar(rng, n);
            while ref_span(&c) > flen {
                let n = rng.range(40, 100) as u32;
                c = cigar(rng, n);
            }
            let pos = if left {
                start
            } else {
                start + flen - ref_span(&c)
            };
            let b = bases(rng, &seqs[tid], pos, &c, top, &mut state, tid, &mut states);
            (pos, c, b)
        };
        let (lp, lc, lb) = mk(&mut rng, true);
        let (rp, rc, rb) = mk(&mut rng, false);
        let cb = if rng.chance(95) {
            Some(format!("C{}", rng.range(0, 6)))
        } else {
            None
        };
        let dup = rng.chance(3);
        let yd = rng.chance(50);
        let name = format!("f{i}");
        // OT: read 1 forward on the left; OB: read 1 reverse on the right.
        for (left, (pos, c, b)) in [(true, (lp, lc.clone(), lb)), (false, (rp, rc.clone(), rb))] {
            let read1 = left == top;
            let mut flag = 0x1 | 0x2 | if read1 { 0x40 } else { 0x80 };
            flag |= if left { 0x20 } else { 0x10 };
            if dup {
                flag |= 0x400;
            }
            if rng.chance(2) {
                flag |= 0x200;
            }
            let qual = if rng.chance(10) {
                None
            } else {
                Some(
                    (0..b.len())
                        .map(|_| b'!' + rng.range(2, 41) as u8)
                        .collect(),
                )
            };
            let mut r = Read {
                name: name.clone(),
                flag,
                tid,
                pos,
                mapq: if rng.chance(5) { 3 } else { 60 },
                cigar: c,
                mtid: tid as i32,
                mpos: if left { rp } else { lp },
                tlen: if left { flen as i64 } else { -(flen as i64) },
                tags: Vec::new(),
                seq: Some(b),
                qual,
            };
            r.tags
                .push(format!("XG:Z:{}", if top { "CT" } else { "GA" }));
            if yd {
                r.tags.push(format!("YD:Z:{}", if top { "f" } else { "r" }));
            }
            r.tags.push(format!("XM:Z:{}", xm_tag(&r, &seqs[tid], top)));
            r.tags.push(format!("MD:Z:{}", md_tag(&r, &seqs[tid])));
            if let Some(cb) = &cb {
                r.tags.push(format!("CB:Z:{cb}"));
            }
            reads.push(r);
        }
    }
    reads.sort_by_key(|r| (r.tid, r.pos));
    let bam = write_bam(dir.path(), "bs", &refs, &reads, seed);
    let fasta = write_fasta(dir.path(), "genome", &refs, &seqs);
    let mut sites = BTreeSet::new();
    for _ in 0..150 {
        let t = rng.range(0, 2) as usize;
        sites.insert((t, rng.range(0, refs[t].1 as u64) as u32));
    }
    let bcf_sites: Vec<(usize, u32, u8, u8)> = sites
        .iter()
        .map(|&(t, p)| {
            let r = seqs[t][p as usize];
            (t, p, r, if r == b'A' { b'C' } else { b'A' })
        })
        .collect();
    let bcf = write_bcf(dir.path(), "snps", &refs, &bcf_sites);
    Fixture {
        dir,
        refs,
        seqs,
        reads,
        snps: sites,
        bam,
        fasta,
        bcf,
    }
}

/// How the expected calls are made.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Mode {
    /// Reference bases from the genome.
    Genome,
    /// Reference bases from `MD`, so only aligned neighbours are known.
    Md,
    /// The `XM` tag.
    Xm,
}

/// One call: position, code (context | reverse << 2), methylated, query
/// position.
type Call = (u32, u8, bool, usize);

fn passes(r: &Read) -> bool {
    r.flag & 0xf04 == 0 && r.mapq >= 10
}

fn qual_ok(r: &Read, q: usize) -> bool {
    r.qual.as_ref().is_none_or(|v| v[q] - 33 >= 5)
}

fn top_of(r: &Read) -> bool {
    match r.tag("YD") {
        Some("f") => return true,
        Some("r") => return false,
        _ => {}
    }
    r.tag("XG") == Some("CT")
}

fn calls_of(f: &Fixture, r: &Read, mode: Mode, masked: bool) -> Vec<Call> {
    let seq = &f.seqs[r.tid];
    let bases = r.seq.as_ref().unwrap();
    let pairs = aligned(r);
    let covered: BTreeSet<u32> = pairs.iter().map(|x| x.0).collect();
    let mask = |p: u32| masked && f.snps.contains(&(r.tid, p));
    let mut out = Vec::new();
    if mode == Mode::Xm {
        let xm = r.tag("XM").unwrap().as_bytes();
        let rev = r.tag("XG") == Some("GA");
        for (p, q) in pairs {
            let ctx = match xm[q].to_ascii_lowercase() {
                b'z' => 0,
                b'x' => 1,
                b'h' => 2,
                _ => continue,
            };
            if qual_ok(r, q) && !mask(p) {
                out.push((p, ctx | (rev as u8) << 2, xm[q].is_ascii_uppercase(), q));
            }
        }
        return out;
    }
    let top = top_of(r);
    let (target, conv) = if top { (b'C', b'T') } else { (b'G', b'A') };
    let near = |p: i64| -> u8 {
        match mode {
            Mode::Md if p < 0 || !covered.contains(&(p as u32)) => b'N',
            _ => base_at(seq, p),
        }
    };
    for (p, q) in pairs {
        if seq[p as usize] != target || mask(p) || !qual_ok(r, q) {
            continue;
        }
        let meth = match bases[q] {
            b if b == target => true,
            b if b == conv => false,
            _ => continue,
        };
        let pi = p as i64;
        let ctx = if top {
            context(near(pi + 1), near(pi + 2))
        } else {
            context(complement(near(pi - 1)), complement(near(pi - 2)))
        };
        if let Some(c) = ctx {
            out.push((p, c | (!top as u8) << 2, meth, q));
        }
    }
    out
}

/// Each counted read's calls after overlapping mates are resolved: read 2's
/// calls inside read 1's span dropped.
fn counted(f: &Fixture, mode: Mode, masked: bool) -> Vec<(&Read, Vec<Call>)> {
    let kept: Vec<&Read> = f.reads.iter().filter(|r| passes(r)).collect();
    let mut by_name: HashMap<&str, Vec<&Read>> = HashMap::new();
    for r in &kept {
        by_name.entry(&r.name).or_default().push(r);
    }
    let mut out = Vec::new();
    for r in kept {
        let mut calls = calls_of(f, r, mode, masked);
        let mates = &by_name[r.name.as_str()];
        if r.flag & 0x80 != 0
            && let Some(one) = mates.iter().find(|m| m.flag & 0x40 != 0)
        {
            let (s, e) = (one.pos, one.end());
            calls.retain(|c| c.0 < s || c.0 >= e);
        }
        out.push((r, calls));
    }
    out
}

/// `(tid, pos, code) → (methylated, unmethylated)` over `contexts`.
fn cytosines(
    f: &Fixture,
    mode: Mode,
    masked: bool,
    contexts: u8,
) -> BTreeMap<(usize, u32, u8), (u32, u32)> {
    let mut out: BTreeMap<(usize, u32, u8), (u32, u32)> = BTreeMap::new();
    for (r, calls) in counted(f, mode, masked) {
        for (p, code, meth, _) in calls {
            if contexts & (1 << (code & 3)) == 0 {
                continue;
            }
            let e = out.entry((r.tid, p, code)).or_default();
            if meth {
                e.0 += 1;
            } else {
                e.1 += 1;
            }
        }
    }
    out
}

/// Table lines `(tid, start, end, context, strand, meth, unmeth)`, the
/// CpGs of both strands merged with `merge`.
type Line = (usize, u32, u32, usize, char, u32, u32);

fn lines(c: &BTreeMap<(usize, u32, u8), (u32, u32)>, merge: bool) -> Vec<Line> {
    // Lines at one position come forward strand first, then by context.
    type Key = (usize, u32, u32, bool, usize, char);
    let mut m: BTreeMap<Key, (u32, u32)> = BTreeMap::new();
    for (&(t, p, code), &(a, b)) in c {
        let ctx = (code & 3) as usize;
        let rev = code & 4 != 0;
        let key = if merge && ctx == 0 {
            let s = if rev { p - 1 } else { p };
            (t, s, s + 2, false, 0, '.')
        } else {
            (t, p, p + 1, rev, ctx, if rev { '-' } else { '+' })
        };
        let e = m.entry(key).or_default();
        e.0 += a;
        e.1 += b;
    }
    m.into_iter()
        .map(|((t, s, e, _, c, st), (a, b))| (t, s, e, c, st, a, b))
        .collect()
}

fn text(f: &Fixture, lines: &[Line]) -> String {
    let mut s = String::new();
    for &(t, a, b, c, st, m, u) in lines {
        s.push_str(&format!(
            "{}\t{a}\t{b}\t{}\t{st}\t{m}\t{u}\n",
            f.refs[t].0, CONTEXTS[c]
        ));
    }
    s
}

/// Equal texts, or a panic naming the first line that differs.
fn same_text(got: &str, want: &str) {
    if got != want {
        let (g, w): (Vec<&str>, Vec<&str>) = (got.lines().collect(), want.lines().collect());
        let i = g
            .iter()
            .zip(&w)
            .position(|(a, b)| a != b)
            .unwrap_or(g.len().min(w.len()));
        panic!(
            "{} vs {} lines; first difference at {i}: got {:?}, want {:?}",
            g.len(),
            w.len(),
            g.get(i.saturating_sub(1)..(i + 2).min(g.len())),
            w.get(i.saturating_sub(1)..(i + 2).min(w.len()))
        );
    }
}

fn product(contexts: u8, merge: bool, source: CallSource) -> Methylation {
    Methylation {
        contexts,
        merge_cpg: merge,
        source,
        matrices: vec![
            MethylMatrix {
                rows: RowKey::Tag(*b"CB"),
                ..MethylMatrix::new("bins", FeatureSource::Bins(2000))
            },
            MethylMatrix {
                contexts: CPG | CHH,
                ..MethylMatrix::new("files", FeatureSource::Bins(5000))
            },
        ],
        ..Methylation::default()
    }
}

fn task(f: &Fixture, out: &Path, m: Methylation, reference: bool, masked: bool) -> Task {
    let mut t = Task::new(vec![f.bam.clone()], out);
    t.reference = reference.then(|| f.fasta.clone());
    t.variants = masked.then(|| f.bcf.clone());
    t.matrix = omicsweave_quant::MatrixFormat::Mtx;
    t.products = vec![Product::Methylation(m)];
    t
}

/// A matrix layer: `(row, feature id) → count`.
type Layer = BTreeMap<(String, String), u32>;

/// The expected matrix layers over `bin`-base bins: `(row, bin id) →
/// (methylated, total)`.
fn matrix(
    f: &Fixture,
    mode: Mode,
    masked: bool,
    bin: u32,
    contexts: u8,
    row: impl Fn(&Read) -> Option<String>,
) -> (Layer, Layer) {
    let (mut meth, mut total) = (BTreeMap::new(), BTreeMap::new());
    for (r, calls) in counted(f, mode, masked) {
        let Some(key) = row(r) else { continue };
        for (p, code, m, _) in calls {
            if contexts & (1 << (code & 3)) == 0 {
                continue;
            }
            let s = p / bin * bin;
            let e = (s + bin).min(f.refs[r.tid].1);
            let id = format!("{}:{s}-{e}", f.refs[r.tid].0);
            *total.entry((key.clone(), id.clone())).or_insert(0) += 1;
            if m {
                *meth.entry((key.clone(), id)).or_insert(0) += 1;
            }
        }
    }
    (meth, total)
}

fn mbias(f: &Fixture, mode: Mode, masked: bool) -> String {
    let mut m: BTreeMap<(usize, usize, usize), (u64, u64)> = BTreeMap::new();
    for r in f.reads.iter().filter(|r| passes(r)) {
        let l_seq = r.seq.as_ref().unwrap().len();
        let lead = match r.cigar.first() {
            Some(&(l, 'H')) => l as usize,
            _ => 0,
        };
        for (_, code, meth, q) in calls_of(f, r, mode, masked) {
            let cycle = if r.reverse() { l_seq - 1 - q } else { lead + q };
            let mate = (r.flag & 0x80 != 0) as usize;
            let e = m.entry(((code & 3) as usize, mate, cycle)).or_default();
            if meth {
                e.0 += 1;
            } else {
                e.1 += 1;
            }
        }
    }
    let mut s = String::from("context\tmate\tposition\tmethylated\tunmethylated\n");
    // Every position up to the longest one seen is listed.
    let mut longest: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for &(c, mate, p) in m.keys() {
        let l = longest.entry((c, mate)).or_insert(0);
        *l = (*l).max(p + 1);
    }
    for (c, name) in CONTEXTS.iter().enumerate() {
        for mate in 0..2 {
            for p in 0..longest.get(&(c, mate)).copied().unwrap_or(0) {
                let (a, b) = m.get(&(c, mate, p)).copied().unwrap_or((0, 0));
                s.push_str(&format!("{name}\t{}\t{}\t{a}\t{b}\n", mate + 1, p + 1));
            }
        }
    }
    s
}

#[test]
fn reference_calls_match_the_reads() {
    let f = fixture(71);
    let out = f.dir.path().join("ref");
    quantify(
        task(
            &f,
            &out,
            product(CPG | CHG | CHH, false, CallSource::Reference),
            true,
            true,
        ),
        Budget::new(3, 256 << 20),
        None,
    )
    .unwrap();
    let cyt = cytosines(&f, Mode::Genome, true, CPG | CHG | CHH);
    assert!(cyt.len() > 5000);
    let want = lines(&cyt, false);
    let table = out.join("methylation/cytosines.tsv.gz");
    assert_eq!(gunzip(&table), text(&f, &want));
    assert_tabix_as_htslib(&table, f.dir.path());
    // The tracks hold each line's beta and coverage.
    let beta = bigwig_values(&out.join("methylation/beta.bw"), &f.refs);
    let cov = bigwig_values(&out.join("methylation/coverage.bw"), &f.refs);
    let mut want_beta: Vec<Vec<f32>> = f.refs.iter().map(|r| vec![0.0; r.1 as usize]).collect();
    let mut want_cov = want_beta.clone();
    for &(t, s, _, _, _, m, u) in &want {
        want_cov[t][s as usize] = (m + u) as f32;
        want_beta[t][s as usize] = (m as f64 / (m + u) as f64) as f32;
    }
    assert_eq!(cov, want_cov);
    assert_eq!(beta, want_beta);
    assert_eq!(
        std::fs::read_to_string(out.join("methylation/mbias.tsv")).unwrap(),
        mbias(&f, Mode::Genome, true)
    );
    let (meth, total) = matrix(&f, Mode::Genome, true, 2000, CPG, |r| {
        r.tag("CB").map(str::to_owned)
    });
    assert!(!meth.is_empty());
    assert_eq!(
        read_10x_layer(&out.join("methylation/bins"), "methylated"),
        meth
    );
    assert_eq!(
        read_10x_layer(&out.join("methylation/bins"), "total"),
        total
    );
    let bam = f.bam.display().to_string();
    let (meth, total) = matrix(&f, Mode::Genome, true, 5000, CPG | CHH, |_| {
        Some(bam.clone())
    });
    assert_eq!(
        read_10x_layer(&out.join("methylation/files"), "methylated"),
        meth
    );
    assert_eq!(
        read_10x_layer(&out.join("methylation/files"), "total"),
        total
    );

    let summary = std::fs::read_to_string(out.join("summary.json")).unwrap();
    let json: Vec<&str> = summary.lines().collect();
    let methylated: u64 = cyt.values().map(|v| v.0 as u64).sum();
    assert!(summary.contains("\"mbias\""));
    assert!(
        json.iter()
            .any(|l| l.trim() == format!("\"cytosines\": {},", cyt.len())),
        "{summary}"
    );
    assert!(methylated > 0);
    let dropped = json
        .iter()
        .find_map(|l| l.trim().strip_prefix("\"overlap_calls_dropped\": "))
        .unwrap();
    assert_ne!(dropped, "0,", "overlapping mates count once");

    // Without the mask, the SNP positions are called too.
    let open = f.dir.path().join("open");
    quantify(
        task(
            &f,
            &open,
            product(CPG, true, CallSource::Reference),
            true,
            false,
        ),
        Budget::new(2, 256 << 20),
        None,
    )
    .unwrap();
    let want = lines(&cytosines(&f, Mode::Genome, false, CPG), true);
    assert!(want.iter().any(|l| l.4 == '.'));
    assert_eq!(
        gunzip(&open.join("methylation/cytosines.tsv.gz")),
        text(&f, &want)
    );
}

#[test]
fn md_and_bismark_tags_give_their_calls() {
    let f = fixture(72);
    let md = f.dir.path().join("md");
    quantify(
        task(
            &f,
            &md,
            product(CPG | CHG | CHH, false, CallSource::Reference),
            false,
            false,
        ),
        Budget::new(4, 256 << 20),
        None,
    )
    .unwrap();
    let want = lines(&cytosines(&f, Mode::Md, false, CPG | CHG | CHH), false);
    assert_eq!(
        gunzip(&md.join("methylation/cytosines.tsv.gz")),
        text(&f, &want)
    );

    // Auto takes XM where it is, with or without a reference.
    for (name, reference) in [("xm", false), ("xm-ref", true)] {
        let out = f.dir.path().join(name);
        quantify(
            task(
                &f,
                &out,
                product(CPG | CHH, true, CallSource::Auto),
                reference,
                true,
            ),
            Budget::new(2, 256 << 20),
            None,
        )
        .unwrap();
        let want = lines(&cytosines(&f, Mode::Xm, true, CPG | CHH), true);
        assert_eq!(
            gunzip(&out.join("methylation/cytosines.tsv.gz")),
            text(&f, &want)
        );
        assert_eq!(
            std::fs::read_to_string(out.join("methylation/mbias.tsv")).unwrap(),
            mbias(&f, Mode::Xm, true)
        );
    }
}

#[test]
fn outputs_do_not_depend_on_threads_or_spilling() {
    let f = fixture(73);
    let m = || product(CPG | CHG, true, CallSource::Reference);
    let a = f.dir.path().join("a");
    quantify(task(&f, &a, m(), true, true), Budget::new(1, 1 << 30), None).unwrap();
    let b = f.dir.path().join("b");
    let mut p = plan(task(&f, &b, m(), true, true), Budget::new(6, 256 << 20)).unwrap();
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
}

/// A long read with `MM` / `ML`: (read, per C of the read as sequenced
/// its ML value or none, explicit mode).
struct Modified {
    read: Read,
    probs: Vec<Option<u8>>,
    explicit: bool,
}

#[test]
fn modification_tags_give_calls() {
    let dir = tempfile::tempdir().unwrap();
    let refs: Vec<(String, u32)> = REFS.iter().map(|(n, l)| (n.to_string(), *l)).collect();
    let mut rng = Rng::new(74);
    let seqs = genome(&mut rng, &refs);
    let mut mods = Vec::new();
    for i in 0..500 {
        let tid = rng.range(0, 2) as usize;
        let n = rng.range(200, 900) as u32;
        let mut c = cigar(&mut rng, n);
        c.retain(|x| x.1 != 'H');
        let span = ref_span(&c);
        let pos = rng.range(0, (refs[tid].1 - span - 1) as u64) as u32;
        let reverse = rng.chance(50);
        let mut seq = Vec::new();
        let mut at = pos;
        for &(l, o) in &c {
            match o {
                'M' => {
                    for k in 0..l {
                        let g = seqs[tid][(at + k) as usize];
                        seq.push(if rng.chance(2) { *rng.pick(b"ACGT") } else { g });
                    }
                    at += l;
                }
                'I' | 'S' => seq.extend((0..l).map(|_| *rng.pick(b"ACGT"))),
                'D' => at += l,
                _ => {}
            }
        }
        // Cs of the read as sequenced: its forward bases, or the
        // complement of the stored ones from the end.
        let n = seq.len();
        let original: Vec<u8> = if reverse {
            (0..n).map(|k| complement(seq[n - 1 - k])).collect()
        } else {
            seq.clone()
        };
        let n_c = original.iter().filter(|&&b| b == b'C').count();
        let explicit = rng.chance(40);
        let probs: Vec<Option<u8>> = (0..n_c)
            .map(|_| rng.chance(70).then(|| rng.range(0, 256) as u8))
            .collect();
        let mut skips = Vec::new();
        let mut ml = Vec::new();
        let mut gap = 0;
        for p in &probs {
            match p {
                Some(v) => {
                    skips.push(gap.to_string());
                    ml.push(v.to_string());
                    gap = 0;
                }
                None => gap += 1,
            }
        }
        if n_c == 0 {
            continue;
        }
        let mode = if explicit {
            "?"
        } else if rng.chance(50) {
            "."
        } else {
            ""
        };
        // A 5hmC list (the first C) comes before the 5mC one, its ML value
        // first.
        let tags = vec![
            format!(
                "MM:Z:C+h?,0;C+m{mode}{}{};",
                if skips.is_empty() { "" } else { "," },
                skips.join(",")
            ),
            format!(
                "ML:B:C,7{}{}",
                if ml.is_empty() { "" } else { "," },
                ml.join(",")
            ),
        ];
        let qual: Vec<u8> = (0..n).map(|_| b'!' + rng.range(0, 41) as u8).collect();
        mods.push(Modified {
            read: Read {
                name: format!("m{i}"),
                flag: if reverse { 0x10 } else { 0 },
                tid,
                pos,
                mapq: 60,
                cigar: c,
                mtid: -1,
                mpos: 0,
                tlen: 0,
                tags,
                seq: Some(seq),
                qual: Some(qual),
            },
            probs,
            explicit,
        });
    }
    mods.sort_by_key(|m| (m.read.tid, m.read.pos));
    let reads: Vec<Read> = mods.iter().map(|m| m.read.clone()).collect();
    let bam = write_bam(dir.path(), "mm", &refs, &reads, 74);
    let fasta = write_fasta(dir.path(), "genome", &refs, &seqs);

    let expect = |reference: bool, p: f64| -> BTreeMap<(usize, u32, u8), (u32, u32)> {
        let mut out: BTreeMap<(usize, u32, u8), (u32, u32)> = BTreeMap::new();
        for m in &mods {
            let r = &m.read;
            let bases = r.seq.as_ref().unwrap();
            let n = bases.len();
            let rev = r.reverse();
            let to_ref: HashMap<usize, u32> = aligned(r).into_iter().map(|(p, q)| (q, p)).collect();
            let mut k_c = 0;
            for k in 0..n {
                let b = if rev {
                    complement(bases[n - 1 - k])
                } else {
                    bases[k]
                };
                if b != b'C' {
                    continue;
                }
                let prob = m.probs[k_c];
                k_c += 1;
                let meth = match prob {
                    Some(v) => {
                        let x = (v as f64 + 0.5) / 256.0;
                        if x >= p {
                            true
                        } else if x <= 1.0 - p {
                            false
                        } else {
                            continue;
                        }
                    }
                    None if m.explicit => continue,
                    None => false,
                };
                let q = if rev { n - 1 - k } else { k };
                let Some(&pos) = to_ref.get(&q) else { continue };
                if r.qual.as_ref().unwrap()[q] - 33 < 5 {
                    continue;
                }
                let seq = &seqs[r.tid];
                let ctx = if reference {
                    let want = if rev { b'G' } else { b'C' };
                    if seq[pos as usize] != want {
                        continue;
                    }
                    genome_context(seq, pos, !rev)
                } else {
                    let read = |d: i64| base_at(bases, q as i64 + d);
                    if rev {
                        context(complement(read(-1)), complement(read(-2)))
                    } else {
                        context(read(1), read(2))
                    }
                };
                let Some(c) = ctx else { continue };
                let e = out.entry((r.tid, pos, c | (rev as u8) << 2)).or_default();
                if meth {
                    e.0 += 1;
                } else {
                    e.1 += 1;
                }
            }
        }
        out
    };
    let fx = Fixture {
        dir,
        refs: refs.clone(),
        seqs: seqs.clone(),
        reads: Vec::new(),
        snps: BTreeSet::new(),
        bam: bam.clone(),
        fasta: fasta.clone(),
        bcf: PathBuf::new(),
    };
    for (name, reference, p) in [("mm-ref", true, 0.5), ("mm-read", false, 0.8)] {
        let out = fx.dir.path().join(name);
        let m = Methylation {
            contexts: CPG | CHG | CHH,
            min_probability: p,
            matrices: Vec::new(),
            ..Methylation::default()
        };
        quantify(
            task(&fx, &out, m, reference, false),
            Budget::new(3, 256 << 20),
            None,
        )
        .unwrap();
        let want = lines(&expect(reference, p), false);
        assert!(want.len() > 1000);
        same_text(
            &gunzip(&out.join("methylation/cytosines.tsv.gz")),
            &text(&fx, &want),
        );
    }
}
