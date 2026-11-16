use super::*;

/// The junctions of a read: `(start, end, overhang)`.
pub(super) fn junctions_of(r: &Read) -> Vec<(u32, u32, u32)> {
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

/// A mate pair: its name with both positions.
pub(super) type PairKey = (String, (usize, u32), (i32, u32));

/// The mate pair a record belongs to.
pub(super) fn pair_key(r: &Read) -> PairKey {
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

pub(super) fn expects_mate(r: &Read) -> bool {
    r.flag & 1 != 0 && r.flag & 0xc == 0 && r.mtid >= 0
}

/// A counted unit: its records, in record order.
pub(super) fn units(reads: &[Read], keep: impl Fn(&Read) -> bool) -> Vec<Vec<&Read>> {
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
pub(super) fn unit_junctions(u: &[&Read]) -> BTreeMap<(usize, u32, u32), u32> {
    let mut m = BTreeMap::new();
    for r in u {
        for (s, e, o) in junctions_of(r) {
            let x = m.entry((r.tid, s, e)).or_insert(0);
            *x = (*x).max(o);
        }
    }
    m
}

pub(super) fn table_filter(r: &Read) -> bool {
    r.flag & (0x4 | 0x200 | 0x800) == 0 && !junctions_of(r).is_empty()
}

/// `SJ.out.tab` lines computed from the reads, the genes and the genome.
pub(super) fn expected_table(f: &Fixture, with_motifs: bool) -> Vec<String> {
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

/// velocyto's rules on the generated reads: per read, exonic or intronic
/// against each transcript of its one gene holding it; per molecule, the
/// OR of its reads.
pub(super) fn classify(g: &Gene, blocks: &[(u32, u32)]) -> u32 {
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

pub(super) type Layers = [BTreeMap<(String, String), u32>; 3];

pub(super) fn expected_velocity(f: &Fixture, single_cell: bool) -> Layers {
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
