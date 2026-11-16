//! Step ⑦ and execution: needs are linked into kernel modes, sequences are
//! packed into tiles, each tile is encoded once and scanned by every motif.

use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc;
use std::sync::{Condvar, Mutex};

use omicsweave_core::{Error, Result, Strings};
use omicsweave_ffm::INVALID;

use crate::lower::{Block, Lowered, Packing, Reader};
use crate::model::{StrandModel, Threshold};
use crate::program::Program;
use crate::reduce::{Demand, Match, Need, Reducer, ScanInfo, SequenceChunk, Window, Windows};

/// Scan execution settings; results do not depend on them.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ScanOptions {
    /// Worker threads; `None` uses the available parallelism.
    pub threads: Option<usize>,
    /// Symbols per tile (each tile is encoded once for all motifs).
    pub tile_len: usize,
}

impl Default for ScanOptions {
    fn default() -> ScanOptions {
        ScanOptions {
            threads: None,
            tile_len: 1 << 16,
        }
    }
}

/// What the kernel computes for one motif on one strand.
#[derive(Clone, Copy, Debug)]
struct Mode {
    all: bool,
    /// Emit windows scoring at least this.
    tau: Option<i32>,
    /// Emit each segment's best window.
    best: bool,
    /// Blocks summed for every window before pruning.
    split: usize,
}

/// Links the reducers' demands into a mode per (motif, strand); `None`
/// skips.
fn link(program: &Program, demands: &[Demand]) -> Vec<Vec<Option<Mode>>> {
    (0..program.len())
        .map(|i| {
            let part = program.part_index(i);
            let context = program.part_of(i).background().order();
            let needs: Vec<Need> = demands
                .iter()
                .filter(|d| d.part.is_none_or(|p| p == part))
                .map(|d| d.need)
                .collect();
            let model = program.motif(i);
            model
                .strands()
                .iter()
                .zip(program.lowered(i))
                .enumerate()
                .map(|(si, (sm, low))| {
                    let all = needs.iter().any(|n| matches!(n, Need::All));
                    let best = needs.iter().any(|n| matches!(n, Need::Best));
                    let above: Vec<Threshold> = needs
                        .iter()
                        .filter_map(|n| match n {
                            Need::Above(t) | Need::Ordered(t) => Some(*t),
                            _ => None,
                        })
                        .collect();
                    if !all && !best && above.is_empty() {
                        return None;
                    }
                    let tau = above.iter().filter_map(|t| sm.full().cutoff(*t)).min();
                    let ambiguity = program.part_of(i).scoring().ambiguity.rule
                        != crate::scoring::AmbiguityRule::Skip;
                    if !all && !best && tau.is_none() && context == 0 && !ambiguity {
                        return None;
                    }
                    let nb = low.blocks.len();
                    let split = if all {
                        nb
                    } else {
                        let target = tau.or_else(|| sm.full().cutoff(Threshold::PValue(1e-3)));
                        choose_split(program, i, low, target, (i * 2 + si) as u64)
                    };
                    Some(Mode {
                        all,
                        tau,
                        best,
                        split,
                    })
                })
                .collect()
        })
        .collect()
}

/// Picks how many blocks to sum densely before pruning, from the survival
/// rate of background windows (a speed choice only).
fn choose_split(
    program: &Program,
    motif: usize,
    low: &Lowered,
    target: Option<i32>,
    seed: u64,
) -> usize {
    let nb = low.blocks.len();
    if nb <= 1 {
        return nb;
    }
    let Some(tau) = target else { return 1 };
    let packing = program.packing();
    let letters = program.part_of(motif).background().letters();
    let mut cdf = Vec::with_capacity(letters.len());
    let mut acc = 0.0;
    for &p in letters {
        acc += p;
        cdf.push(acc);
    }
    let footprint = (low.width as i32 - low.lo) as usize + packing.k;
    let samples = 512;
    let mut survive = vec![0usize; nb + 1];
    let mut rng = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    let mut codes = vec![0u8; footprint];
    for _ in 0..samples {
        for c in codes.iter_mut() {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            let u = (rng >> 11) as f64 / (1u64 << 53) as f64 * acc;
            *c = cdf.iter().position(|&x| u < x).unwrap_or(cdf.len() - 1) as u8;
        }
        let mut s = 0i32;
        for (j, b) in low.blocks.iter().enumerate() {
            if s + low.rest_max[j] >= tau {
                survive[j] += 1;
            }
            s += lookup_codes(b, &codes, (-low.lo) as usize, packing);
        }
        if s >= tau {
            survive[nb] += 1;
        }
    }
    let mut best = (f64::INFINITY, 1usize);
    for (j, &alive) in survive.iter().enumerate().skip(1) {
        let rate = alive as f64 / samples as f64;
        let cost = j as f64 + rate * (nb - j) as f64 * 3.0;
        if cost < best.0 {
            best = (cost, j);
        }
    }
    best.1
}

fn lookup_codes(b: &Block, codes: &[u8], at: usize, packing: &Packing) -> i32 {
    let base = (at as isize + b.offset as isize) as usize;
    match &b.reader {
        Reader::Kmer { shift } => {
            let len = packing.k - (*shift / packing.bits) as usize;
            let idx = codes[base..base + len]
                .iter()
                .fold(0usize, |acc, &c| (acc << packing.bits) | c as usize);
            b.table[idx]
        }
        Reader::Sparse { positions } => {
            let idx = positions.iter().fold(0usize, |acc, &p| {
                acc * packing.n + codes[(base as isize + p as isize) as usize] as usize
            });
            b.table[idx]
        }
    }
}

/// Owned windows of one sequence stretch inside a tile.
#[derive(Clone, Copy, Debug)]
struct Segment {
    sequence: u32,
    /// Sequence position of the first owned window.
    first: u32,
    /// Owned windows (window starts `first .. first + count`).
    count: u32,
    /// Tile index of the first owned window.
    at: usize,
    /// Sequence bytes copied and where they go in the tile.
    bytes: (usize, usize),
    copy_at: usize,
}

#[derive(Clone, Debug, Default)]
struct TilePlan {
    segments: Vec<Segment>,
    len: usize,
}

/// Packs sequences into tiles: each tile starts with `lead` unknown
/// symbols, each segment carries up to `context` symbols before its first
/// window and the bases its last window reads, segments are separated by an
/// unknown symbol, and the tile ends with `tail` unknown symbols.
fn plan_tiles(
    sequences: &Strings<'_>,
    tile_len: usize,
    context: usize,
    max_width: usize,
    tail: usize,
) -> Vec<TilePlan> {
    let lead = context + 1;
    let cap = tile_len.max(1024);
    let min_chunk = 256usize;
    let mut tiles = Vec::new();
    let mut cur = TilePlan {
        segments: Vec::new(),
        len: lead,
    };
    for i in 0..sequences.len() {
        let len = sequences.get(i).len();
        let mut a = 0usize;
        while a < len {
            let pre = a.min(context);
            let fixed = pre + max_width.saturating_sub(1) + 1;
            let avail = cap.saturating_sub(cur.len + fixed + tail);
            if !cur.segments.is_empty() && avail < min_chunk.min(len - a) {
                tiles.push(std::mem::replace(
                    &mut cur,
                    TilePlan {
                        segments: Vec::new(),
                        len: lead,
                    },
                ));
                continue;
            }
            let take = (len - a).min(avail.max(1));
            let b = a + take;
            let lo = a - pre;
            let hi = len.min(b + max_width.saturating_sub(1));
            let copy_at = cur.len;
            cur.segments.push(Segment {
                sequence: i as u32,
                first: a as u32,
                count: take as u32,
                at: copy_at + pre,
                bytes: (lo, hi),
                copy_at,
            });
            cur.len = copy_at + (hi - lo) + 1;
            a = b;
        }
    }
    if !cur.segments.is_empty() {
        tiles.push(cur);
    }
    for t in &mut tiles {
        t.len += tail;
    }
    tiles
}

/// A tile encoded once for every motif.
#[derive(Default)]
struct Encoded {
    /// Global codes, `INVALID` for unknown symbols.
    codes: Vec<u8>,
    /// The same with unknown symbols as 0, for lookups whose window is
    /// discarded anyway.
    zeroed: Vec<u8>,
    /// Packed k-mers starting at each index (unknown symbols packed as 0).
    kmer: Vec<u32>,
    /// Known symbols in a row starting at each index.
    run: Vec<u32>,
    /// The segment owning the window at each index, or `NO_OWNER`.
    owner: Vec<u32>,
    /// Owned indices with less background context than the order.
    edges: Vec<usize>,
    /// Ambiguity code index at each index (`0xFF` for none).
    ambiguity: Vec<u8>,
    /// Alphabet or ambiguity symbols in a row starting at each index.
    soft_run: Vec<u32>,
    /// Indices holding an ambiguity code, ascending.
    ambiguous_at: Vec<usize>,
    /// Known symbols in a row ending just before each index, at most the
    /// background order.
    before: Vec<u8>,
}

fn encode(
    plan: &TilePlan,
    sequences: &Strings<'_>,
    lut: &[u8; 256],
    ambiguity_lut: &[u8; 256],
    packing: &Packing,
    context: usize,
    out: &mut Encoded,
) {
    let len = plan.len;
    out.codes.clear();
    out.codes.resize(len, INVALID);
    out.ambiguity.clear();
    out.ambiguity.resize(len, 0xFF);
    out.ambiguous_at.clear();
    for s in &plan.segments {
        let src = &sequences.get(s.sequence as usize)[s.bytes.0..s.bytes.1];
        for (i, &b) in src.iter().enumerate() {
            let at = s.copy_at + i;
            out.codes[at] = lut[b as usize];
            let a = ambiguity_lut[b as usize];
            if a != 0xFF {
                out.ambiguity[at] = a;
                out.ambiguous_at.push(at);
            }
        }
    }
    out.soft_run.clear();
    out.soft_run.resize(len + 1, 0);
    if !out.ambiguous_at.is_empty() {
        for i in (0..len).rev() {
            out.soft_run[i] = if out.codes[i] == INVALID && out.ambiguity[i] == 0xFF {
                0
            } else {
                out.soft_run[i + 1] + 1
            };
        }
    }
    out.owner.clear();
    out.owner.resize(len, NO_OWNER);
    for (si, seg) in plan.segments.iter().enumerate() {
        out.owner[seg.at..seg.at + seg.count as usize].fill(si as u32);
    }
    out.zeroed.clear();
    out.zeroed
        .extend(out.codes.iter().map(|&c| if c == INVALID { 0 } else { c }));
    let bits = packing.bits;
    let mask = packing.mask();
    let k = packing.k;
    out.kmer.clear();
    out.kmer.resize(len, 0);
    let mut v = 0u32;
    for (i, &c) in out.codes.iter().enumerate() {
        let c = if c == INVALID { 0 } else { c as u32 };
        v = ((v << bits) | c) & mask;
        if i + 1 >= k {
            out.kmer[i + 1 - k] = v;
        }
    }
    out.run.clear();
    out.run.resize(len + 1, 0);
    for i in (0..len).rev() {
        out.run[i] = if out.codes[i] == INVALID {
            0
        } else {
            out.run[i + 1] + 1
        };
    }
    out.edges.clear();
    if context > 0 {
        out.before.clear();
        out.before.resize(len, 0);
        for i in 1..len {
            out.before[i] = if out.codes[i - 1] == INVALID {
                0
            } else {
                (out.before[i - 1] as usize + 1).min(context) as u8
            };
            if (out.before[i] as usize) < context
                && out.owner[i] != NO_OWNER
                && out.codes[i] != INVALID
            {
                out.edges.push(i);
            }
        }
    }
}

/// No segment owns this tile index.
const NO_OWNER: u32 = u32::MAX;

/// One tile's ordered-stage matches, tagged with their sequence.
type TileMatches = (usize, Vec<(u32, Match)>);

/// Scans a tile's owned windows for one motif on one strand.
#[allow(clippy::too_many_arguments)]
fn kernel(
    plan: &TilePlan,
    enc: &Encoded,
    low: &Lowered,
    model: &StrandModel,
    mode: &Mode,
    packing: &Packing,
    context: usize,
    ambiguous: bool,
    acc: &mut Vec<i32>,
    out: &mut Vec<Window>,
) {
    let (Some(first), Some(last)) = (plan.segments.first(), plan.segments.last()) else {
        return;
    };
    let lo = first.at;
    let count = last.at + last.count as usize - lo;
    let width = low.width;
    acc.clear();
    acc.resize(count, 0);
    for b in &low.blocks[..mode.split] {
        let base = (lo as isize + b.offset as isize) as usize;
        match &b.reader {
            Reader::Kmer { shift } => {
                let kmers = &enc.kmer[base..base + count];
                for (a, &km) in acc.iter_mut().zip(kmers) {
                    // SAFETY: lower allocates 2^(bits*len) cells for this
                    // block; encode masks each k-mer to bits*k bits, and
                    // shift removes bits*(k-len), including unused codes.
                    *a += unsafe { *b.table.get_unchecked((km >> shift) as usize) };
                }
            }
            Reader::Sparse { .. } => {
                for (i, a) in acc.iter_mut().enumerate() {
                    *a += lookup(b, enc, lo + i, packing);
                }
            }
        }
    }
    let window = |o: u32, t: usize, context: usize, score: i32| {
        let seg = &plan.segments[o as usize];
        Window {
            sequence: seg.sequence,
            start: seg.first + (t - seg.at) as u32,
            context: context as u8,
            score,
        }
    };
    // Windows with truncated background context use their own models.
    for &t in &enc.edges {
        if context > 0 && (enc.before[t] as usize) < context && enc.run[t] >= width {
            let c = enc.before[t] as usize;
            if let Some(s) = model.variant(c).int_model().score(&enc.codes, t) {
                out.push(window(enc.owner[t], t, c, s));
            }
        }
    }
    // Windows holding ambiguity codes (and nothing unknown) score through
    // the codes' cells; they are never ordinary windows.
    if ambiguous && !enc.ambiguous_at.is_empty() {
        let w = width as usize;
        let int = model.full().int_model();
        let mut next = 0usize;
        for &a in &enc.ambiguous_at {
            for t in a.saturating_sub(w - 1).max(next)..=a {
                if enc.owner[t] == NO_OWNER || enc.run[t] >= width || enc.soft_run[t] < width {
                    continue;
                }
                if let Some(s) = int.score_ambiguous(&enc.codes, &enc.ambiguity, t) {
                    out.push(window(enc.owner[t], t, context, s));
                }
            }
            next = a + 1;
        }
    }
    let ordinary = |t: usize| {
        enc.owner[t] != NO_OWNER
            && enc.run[t] >= width
            && (context == 0 || enc.before[t] as usize >= context)
    };
    let rest = &low.blocks[mode.split..];
    // Completes a window's score; `None` when it cannot reach `floor`.
    let complete = |t: usize, partial: i32, floor: i32| -> Option<i32> {
        let mut s = partial;
        for (j, b) in rest.iter().enumerate() {
            s += lookup(b, enc, t, packing);
            if s + low.rest_max[mode.split + j + 1] < floor {
                return None;
            }
        }
        (s >= 0).then_some(s)
    };
    let need = mode.tau.unwrap_or(i32::MAX);
    if mode.all {
        for (i, &partial) in acc.iter().enumerate() {
            let t = lo + i;
            if partial >= 0 && ordinary(t) {
                out.push(window(enc.owner[t], t, context, partial));
            }
        }
        return;
    }
    if !mode.best {
        // Few windows survive the first blocks: filter them in a tight loop.
        let entry = need.saturating_sub(low.rest_max[mode.split]);
        for (i, &partial) in acc.iter().enumerate() {
            if partial < entry {
                continue;
            }
            let t = lo + i;
            if ordinary(t)
                && let Some(s) = complete(t, partial, need)
            {
                out.push(window(enc.owner[t], t, context, s));
            }
        }
        return;
    }
    // The current segment's best window: (score, tile index, emitted).
    let mut best: Option<(i32, usize, bool)> = None;
    let mut best_owner = NO_OWNER;
    for (i, &partial) in acc.iter().enumerate() {
        let t = lo + i;
        if !ordinary(t) {
            continue;
        }
        let o = enc.owner[t];
        if o != best_owner {
            if let Some((s, bt, false)) = best {
                out.push(window(best_owner, bt, context, s));
            }
            best = None;
            best_owner = o;
        }
        let floor = match best {
            Some((b, _, _)) => need.min(b.saturating_add(1)),
            None => 0,
        };
        if partial + low.rest_max[mode.split] < floor {
            continue;
        }
        let Some(s) = complete(t, partial, floor) else {
            continue;
        };
        let emit = s >= need;
        if emit {
            out.push(window(o, t, context, s));
        }
        if best.is_none_or(|(b, _, _)| s > b) {
            best = Some((s, t, emit));
        }
    }
    if let Some((s, bt, false)) = best {
        out.push(window(best_owner, bt, context, s));
    }
}

#[inline]
fn lookup(b: &Block, enc: &Encoded, t: usize, packing: &Packing) -> i32 {
    let base = (t as isize + b.offset as isize) as usize;
    match &b.reader {
        Reader::Kmer { shift } => {
            // SAFETY: kernel visits tile indices between its first and last
            // owned window. plan_tiles adds context before the first and
            // width+k after the last, covering every lowered block offset.
            // encode masks k-mers to bits*k, and lower's table has exactly
            // 2^(bits*k-shift) cells, so the shifted index is in bounds.
            unsafe {
                let km = *enc.kmer.get_unchecked(base);
                *b.table.get_unchecked((km >> shift) as usize)
            }
        }
        Reader::Sparse { positions } => {
            let idx = positions.iter().fold(0usize, |acc, &p| {
                acc * packing.n + enc.zeroed[(base as isize + p as isize) as usize] as usize
            });
            // SAFETY: zeroed codes are <n, and lower clones the checked
            // n^positions.len() factor table. Every prefix index fits usize
            // because that table was successfully allocated by compilation.
            unsafe { *b.table.get_unchecked(idx) }
        }
    }
}

pub(crate) fn check_sequences(sequences: &Strings<'_>) -> Result<()> {
    let bad = |detail: String| Err(Error::invalid("motif scan sequences", detail));
    if sequences.offsets.is_empty() {
        return bad("offsets are empty".into());
    }
    if sequences.len() > u32::MAX as usize {
        return bad("more than 2^32 − 1 sequences".into());
    }
    let mut prev = sequences.offsets[0];
    if prev < 0 {
        return bad("negative offset".into());
    }
    for &o in &sequences.offsets[1..] {
        if o < prev {
            return bad("offsets decrease".into());
        }
        if o - prev > u32::MAX as i64 {
            return bad("a sequence is longer than 2^32 − 1".into());
        }
        prev = o;
    }
    if prev as u64 > sequences.data.len() as u64 {
        return bad("offsets run past the data".into());
    }
    Ok(())
}

pub(crate) fn run<R: Reducer>(
    program: &Program,
    sequences: Strings<'_>,
    mut reducer: R,
    options: &ScanOptions,
) -> Result<R::Output> {
    check_sequences(&sequences)?;
    let info = ScanInfo { sequences };
    reducer.prepare(program, &info)?;
    let mut demands = Vec::new();
    reducer.needs(&mut demands);
    let modes = link(program, &demands);
    let packing = *program.packing();
    let context = program.context();
    let max_width = (0..program.len())
        .map(|i| program.motif(i).width() as usize)
        .max()
        .unwrap_or(1);
    let tiles = plan_tiles(
        &sequences,
        options.tile_len,
        context,
        max_width,
        max_width + packing.k,
    );
    // Work units are (tile, group of motifs), enough to keep every thread
    // busy even when the input fits in a few tiles.
    let threads = options
        .threads
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()))
        .max(1);
    let motif_count = program.len().max(1);
    let groups = threads
        .saturating_mul(4)
        .div_ceil(tiles.len().max(1))
        .clamp(1, motif_count);
    let per_group = motif_count.div_ceil(groups);
    let units = tiles.len() * groups;
    let threads = threads.min(units.max(1));
    let ordered: Vec<(Option<usize>, Threshold)> = demands
        .iter()
        .filter_map(|d| match d.need {
            Need::Ordered(t) => Some((d.part, t)),
            _ => None,
        })
        .collect();
    // Tiles the ordered stage has consumed; workers stay within a window
    // ahead of it so buffered results stay bounded.
    let gate = (Mutex::new(0usize), Condvar::new());
    let ahead = threads * 2 + 2;
    let next = AtomicUsize::new(0);
    let lut = program.alphabet().lut();
    let ambiguity_lut = program.ambiguity_lut();
    let work = |mut local: R, tx: Option<mpsc::Sender<TileMatches>>| -> R {
        let mut enc = Encoded::default();
        let mut encoded = usize::MAX;
        let mut acc = Vec::new();
        let mut buf = Vec::new();
        loop {
            let unit = next.fetch_add(1, Ordering::Relaxed);
            if unit >= units {
                break;
            }
            let (ti, group) = (unit / groups, unit % groups);
            if tx.is_some() {
                let mut done = gate.0.lock().expect("scan gate");
                while ti >= *done + ahead {
                    done = gate.1.wait(done).expect("scan gate");
                }
            }
            let plan = &tiles[ti];
            if encoded != ti {
                encode(
                    plan,
                    &sequences,
                    lut,
                    ambiguity_lut,
                    &packing,
                    context,
                    &mut enc,
                );
                encoded = ti;
            }
            let mut matches = Vec::new();
            let motifs = group * per_group..((group + 1) * per_group).min(program.len());
            for mi in motifs {
                let model = program.motif(mi);
                let part = program.part_index(mi);
                let order = program.part_of(mi).background().order();
                let ambiguous = order == 0
                    && program.part_of(mi).scoring().ambiguity.rule
                        != crate::scoring::AmbiguityRule::Skip;
                for ((sm, low), mode) in model
                    .strands()
                    .iter()
                    .zip(program.lowered(mi))
                    .zip(&modes[mi])
                {
                    let Some(mode) = mode else { continue };
                    buf.clear();
                    kernel(
                        plan, &enc, low, sm, mode, &packing, order, ambiguous, &mut acc, &mut buf,
                    );
                    if buf.is_empty() {
                        continue;
                    }
                    let batch = Windows {
                        motif: mi as u32,
                        part: part as u32,
                        strand: sm.strand(),
                        windows: &buf,
                        model: sm,
                        width: model.width(),
                    };
                    local.windows(&batch);
                    if tx.is_some() {
                        for w in &buf {
                            if ordered
                                .iter()
                                .any(|(p, t)| p.is_none_or(|p| p == part) && batch.meets(w, *t))
                            {
                                matches.push((
                                    w.sequence,
                                    Match {
                                        score: batch.score(w),
                                        pvalue: batch.pvalue(w),
                                        motif: mi as u32,
                                        start: w.start,
                                        end: w.start + model.width(),
                                        strand: sm.strand(),
                                        context: w.context,
                                        int_score: w.score,
                                    },
                                ));
                            }
                        }
                    }
                }
            }
            if let Some(tx) = &tx
                && tx.send((ti, matches)).is_err()
            {
                break;
            }
        }
        local
    };
    let forks: Vec<R> = (0..threads).map(|_| reducer.fork()).collect();
    let done: Vec<R> = if ordered.is_empty() {
        if threads == 1 {
            forks.into_iter().map(|f| work(f, None)).collect()
        } else {
            std::thread::scope(|scope| {
                let handles: Vec<_> = forks
                    .into_iter()
                    .map(|f| scope.spawn(|| work(f, None)))
                    .collect();
                handles
                    .into_iter()
                    .map(|h| h.join().expect("scan worker panicked"))
                    .collect()
            })
        }
    } else {
        let (tx, rx) = mpsc::channel();
        std::thread::scope(|scope| {
            let handles: Vec<_> = forks
                .into_iter()
                .map(|f| {
                    let tx = tx.clone();
                    let work = &work;
                    scope.spawn(move || work(f, Some(tx)))
                })
                .collect();
            drop(tx);
            // The ordered stage: tiles in order, each once all its groups
            // have arrived.
            let mut pending: HashMap<usize, (usize, Vec<(u32, Match)>)> = HashMap::new();
            let mut next_tile = 0usize;
            let mut chunk = Vec::new();
            for (ti, mut matches) in rx {
                let entry = pending.entry(ti).or_default();
                entry.0 += 1;
                entry.1.append(&mut matches);
                while pending.get(&next_tile).is_some_and(|e| e.0 == groups) {
                    let (_, mut all) = pending.remove(&next_tile).expect("complete tile");
                    all.sort_by_key(|(s, m)| (*s, m.start, m.motif, m.strand));
                    let mut at = 0usize;
                    for seg in &tiles[next_tile].segments {
                        let from = at;
                        while at < all.len() && all[at].0 == seg.sequence {
                            at += 1;
                        }
                        chunk.clear();
                        chunk.extend(all[from..at].iter().map(|(_, m)| *m));
                        let len = sequences.get(seg.sequence as usize).len() as u32;
                        reducer.sequence(
                            program,
                            &SequenceChunk {
                                sequence: seg.sequence,
                                windows: seg.first..seg.first + seg.count,
                                first: seg.first == 0,
                                last: seg.first + seg.count == len,
                                matches: &chunk,
                            },
                        );
                    }
                    next_tile += 1;
                    *gate.0.lock().expect("scan gate") = next_tile;
                    gate.1.notify_all();
                }
            }
            handles
                .into_iter()
                .map(|h| h.join().expect("scan worker panicked"))
                .collect()
        })
    };
    for d in done {
        reducer.merge(d);
    }
    Ok(reducer.finish(program, &info))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tiles_own_every_window_once() {
        let data: Vec<u8> = (0..5000).map(|i| b"ACGT"[i % 4]).collect();
        let offsets = [0i64, 10, 4000, 4000, 5000];
        let s = Strings {
            offsets: &offsets,
            data: &data,
        };
        let tiles = plan_tiles(&s, 1024, 2, 12, 17);
        let mut owned: Vec<Vec<bool>> = (0..4).map(|i| vec![false; s.get(i).len()]).collect();
        for t in &tiles {
            assert!(t.len <= 1024 + 17 || t.segments.len() == 1);
            for seg in &t.segments {
                for p in seg.first..seg.first + seg.count {
                    let o = &mut owned[seg.sequence as usize][p as usize];
                    assert!(!*o);
                    *o = true;
                }
            }
        }
        assert!(owned.iter().all(|o| o.iter().all(|&x| x)));
    }
}
