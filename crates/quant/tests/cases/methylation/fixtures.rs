use super::*;

pub(super) const REFS: [(&str, u32); 2] = [("chr1", 30_000), ("chr2", 18_000)];
pub(super) const CONTEXTS: [&str; 3] = ["CpG", "CHG", "CHH"];

pub(super) struct Fixture {
    pub(super) dir: tempfile::TempDir,
    pub(super) refs: Vec<(String, u32)>,
    pub(super) seqs: Vec<Vec<u8>>,
    pub(super) reads: Vec<Read>,
    pub(super) snps: BTreeSet<(usize, u32)>,
    pub(super) bam: PathBuf,
    pub(super) fasta: PathBuf,
    pub(super) bcf: PathBuf,
}

/// Matches with the odd insertion, deletion and soft clip, and now and
/// then a hard clip in front; `len` query bases.
pub(super) fn cigar(rng: &mut Rng, len: u32) -> Vec<(u32, char)> {
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

pub(super) fn ref_span(c: &[(u32, char)]) -> u32 {
    c.iter()
        .filter(|(_, o)| matches!(o, 'M' | 'D' | 'N'))
        .map(|(l, _)| l)
        .sum()
}

/// `(ref_pos, query_pos)` of every `M` base (hard clips ignored).
pub(super) fn aligned(r: &Read) -> Vec<(u32, usize)> {
    common::aligned_pairs(r)
}

/// The context index of a cytosine followed by `n1`, `n2` on its strand.
pub(super) fn context(n1: u8, n2: u8) -> Option<u8> {
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

pub(super) fn base_at(seq: &[u8], p: i64) -> u8 {
    if p < 0 || p as usize >= seq.len() {
        b'N'
    } else {
        seq[p as usize]
    }
}

/// The genome context of the cytosine at `p` on the top or bottom strand.
pub(super) fn genome_context(seq: &[u8], p: u32, top: bool) -> Option<u8> {
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
pub(super) fn bases(
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
pub(super) fn level(tid: usize, p: u32, top: bool) -> u64 {
    let h =
        (tid as u64 * 1_000_003 + p as u64 * 2 + top as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    [0, 15, 85, 100][(h >> 60) as usize % 4]
}

/// Bismark's `XM` for a read: each aligned cytosine of its strand by
/// genome context (`U` unknown), uppercase when read unconverted.
pub(super) fn xm_tag(r: &Read, seq: &[u8], top: bool) -> String {
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

pub(super) fn fixture(seed: u64) -> Fixture {
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

pub(super) fn product(contexts: u8, merge: bool, source: CallSource) -> Methylation {
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

pub(super) fn task(f: &Fixture, out: &Path, m: Methylation, reference: bool, masked: bool) -> Task {
    let mut t = Task::new(vec![f.bam.clone()], out);
    t.reference = reference.then(|| f.fasta.clone());
    t.variants = masked.then(|| f.bcf.clone());
    t.matrix = omicsweave_quant::MatrixFormat::Mtx;
    t.products = vec![Product::Methylation(m)];
    t
}
