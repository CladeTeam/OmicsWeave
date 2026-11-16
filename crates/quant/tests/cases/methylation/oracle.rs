use super::*;

/// How the expected calls are made.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum Mode {
    /// Reference bases from the genome.
    Genome,
    /// Reference bases from `MD`, so only aligned neighbours are known.
    Md,
    /// The `XM` tag.
    Xm,
}

/// One call: position, code (context | reverse << 2), methylated, query
/// position.
pub(super) type Call = (u32, u8, bool, usize);

pub(super) fn passes(r: &Read) -> bool {
    r.flag & 0xf04 == 0 && r.mapq >= 10
}

pub(super) fn qual_ok(r: &Read, q: usize) -> bool {
    r.qual.as_ref().is_none_or(|v| v[q] - 33 >= 5)
}

pub(super) fn top_of(r: &Read) -> bool {
    match r.tag("YD") {
        Some("f") => return true,
        Some("r") => return false,
        _ => {}
    }
    r.tag("XG") == Some("CT")
}

pub(super) fn calls_of(f: &Fixture, r: &Read, mode: Mode, masked: bool) -> Vec<Call> {
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
pub(super) fn counted(f: &Fixture, mode: Mode, masked: bool) -> Vec<(&Read, Vec<Call>)> {
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
pub(super) fn cytosines(
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

/// A matrix layer: `(row, feature id) → count`.
pub(super) type Layer = BTreeMap<(String, String), u32>;

/// The expected matrix layers over `bin`-base bins: `(row, bin id) →
/// (methylated, total)`.
pub(super) fn matrix(
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

pub(super) fn mbias(f: &Fixture, mode: Mode, masked: bool) -> String {
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
