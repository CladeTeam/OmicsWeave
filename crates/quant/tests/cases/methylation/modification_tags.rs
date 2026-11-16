use super::*;

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
