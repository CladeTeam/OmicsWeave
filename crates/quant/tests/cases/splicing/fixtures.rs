use super::*;

pub(super) struct Fixture {
    pub(super) dir: tempfile::TempDir,
    pub(super) spec: Spec,
    pub(super) reads: Vec<Read>,
    pub(super) genes: Vec<Gene>,
    pub(super) seqs: Vec<Vec<u8>>,
    pub(super) bam: PathBuf,
    pub(super) fasta: PathBuf,
    pub(super) twobit: PathBuf,
    pub(super) ann: Arc<Annotation>,
}

/// The merged exon blocks of each transcript.
pub(super) fn transcripts(g: &Gene) -> Vec<Vec<(u32, u32)>> {
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
pub(super) fn on_introns(rng: &mut Rng, genes: &[Gene], n: usize) -> Vec<Read> {
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

pub(super) fn fixture(seed: u64) -> Fixture {
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

pub(super) fn task(
    f: &Fixture,
    out: &Path,
    products: Vec<Product>,
    reference: Option<&Path>,
) -> Task {
    let mut t = Task::new(vec![f.bam.clone()], out);
    t.annotation = Some(Arc::clone(&f.ann));
    t.reference = reference.map(Path::to_owned);
    t.products = products;
    t
}
