use super::scan_fixtures::{Rng, reader, sequences};
use omicsweave_core::{Strand, Strings};
use omicsweave_ffm::{Alphabet, Motif, MotifRef};
use omicsweave_motif::{
    Background, Demand, Dense, Hits, Need, Program, Reducer, ScanInfo, ScanOptions, Scoring,
    SequenceBest, SequenceChunk, SequenceCounts, Strands, Threshold, Windows,
};
/// Every reported window as (motif, strand, sequence, start, context, score).
#[derive(Default)]
pub(crate) struct Collect(pub(crate) Vec<(u32, Strand, u32, u32, u8, i32)>);

impl Reducer for Collect {
    type Output = Vec<(u32, Strand, u32, u32, u8, i32)>;
    fn needs(&self, out: &mut Vec<Demand>) {
        out.push(Need::All.into());
    }
    fn fork(&self) -> Collect {
        Collect::default()
    }
    fn windows(&mut self, batch: &Windows<'_>) {
        for w in batch.windows {
            self.0.push((
                batch.motif,
                batch.strand,
                w.sequence,
                w.start,
                w.context,
                w.score,
            ));
        }
    }
    fn merge(&mut self, mut other: Collect) {
        self.0.append(&mut other.0);
    }
    fn finish(mut self, _: &Program, _: &ScanInfo<'_>) -> Self::Output {
        self.0.sort();
        self.0
    }
}

/// The ordered stage's chunks: (sequence, window range, first, last) and
/// every match as (sequence, motif, start, strand, int score).
struct Ordered {
    threshold: Threshold,
    chunks: Vec<(u32, u32, u32, bool, bool)>,
    matches: Vec<(u32, u32, u32, Strand, i32)>,
}

impl Reducer for Ordered {
    type Output = Ordered;
    fn needs(&self, out: &mut Vec<Demand>) {
        out.push(Need::Ordered(self.threshold).into());
    }
    fn fork(&self) -> Ordered {
        Ordered {
            threshold: self.threshold,
            chunks: Vec::new(),
            matches: Vec::new(),
        }
    }
    fn sequence(&mut self, _: &Program, chunk: &SequenceChunk<'_>) {
        self.chunks.push((
            chunk.sequence,
            chunk.windows.start,
            chunk.windows.end,
            chunk.first,
            chunk.last,
        ));
        let mut prev = None;
        for m in chunk.matches {
            let key = (m.start, m.motif, m.strand);
            assert!(prev.is_none_or(|p| p <= key), "matches out of order");
            prev = Some(key);
            assert!(chunk.windows.contains(&m.start));
            self.matches
                .push((chunk.sequence, m.motif, m.start, m.strand, m.int_score));
        }
    }
    fn merge(&mut self, _: Ordered) {}
    fn finish(self, _: &Program, _: &ScanInfo<'_>) -> Ordered {
        self
    }
}

/// Window-by-window evaluation of the compiled integer models.
pub(crate) fn naive(
    program: &Program,
    seqs: &Strings<'_>,
) -> Vec<(u32, Strand, u32, u32, u8, i32)> {
    let lut = program.alphabet().lut();
    let mut out = Vec::new();
    for si in 0..seqs.len() {
        let bytes = seqs.get(si);
        for mi in 0..program.len() {
            let part = program.part_of(mi);
            let m = part.background().order();
            let mut padded = vec![0u8; m];
            padded.extend(bytes.iter().map(|&b| lut[b as usize]));
            let known = |i: usize| padded[i + m] != 0xFF;
            let model = program.motif(mi);
            let w = model.width() as usize;
            let amb = &part.scoring().ambiguity;
            let code_of = |b: u8| {
                amb.codes.iter().position(|(c, _)| *c == b).or_else(|| {
                    (lut[b as usize] == 0xFF)
                        .then_some(amb.unknown)
                        .flatten()
                        .and_then(|u| amb.codes.iter().position(|(c, _)| *c == u))
                })
            };
            let ambiguity: Vec<u8> = bytes
                .iter()
                .map(|&b| code_of(b).map_or(0xFF, |k| k as u8))
                .collect();
            for t in 0..bytes.len() {
                if t + w > bytes.len() {
                    continue;
                }
                if !(t..t + w).all(known) {
                    // Ambiguity codes score only under a 0-order background.
                    let scored = m == 0
                        && amb.rule != omicsweave_motif::AmbiguityRule::Skip
                        && (t..t + w).all(|i| known(i) || ambiguity[i] != 0xFF);
                    if scored {
                        for sm in model.strands() {
                            if let Some(s) = sm
                                .full()
                                .int_model()
                                .score_ambiguous(&padded, &ambiguity, t)
                            {
                                out.push((mi as u32, sm.strand(), si as u32, t as u32, 0, s));
                            }
                        }
                    }
                    continue;
                }
                let c = (1..=m.min(t)).take_while(|&d| known(t - d)).count();
                for sm in model.strands() {
                    if let Some(s) = sm.variant(c).int_model().score(&padded, t + m) {
                        out.push((mi as u32, sm.strand(), si as u32, t as u32, c as u8, s));
                    }
                }
            }
        }
    }
    out.sort();
    out
}

pub(crate) fn check(motifs: &[Motif], background: &Background, scoring: &Scoring, seed: u64) {
    let alphabet = Alphabet::dna();
    let r = reader(motifs);
    let refs: Vec<MotifRef> = (0..r.len()).map(|i| r.motif(i).unwrap()).collect();
    let views: Vec<_> = refs.iter().map(|m| m.view()).collect();
    let program =
        Program::compile(&views, &alphabet, background, scoring, Strands::Both, None).unwrap();
    let mut rng = Rng(seed);
    let (data, offsets) = sequences(&mut rng, 12);
    let seqs = Strings {
        offsets: &offsets,
        data: &data,
    };
    let expect = naive(&program, &seqs);
    let options = [
        ScanOptions {
            threads: Some(1),
            tile_len: 1024,
        },
        ScanOptions {
            threads: Some(4),
            tile_len: 3000,
        },
        ScanOptions::default(),
    ];
    let p = 2e-3;
    for opt in &options {
        let got = program.scan(seqs, Collect::default(), opt).unwrap();
        assert_eq!(
            got, expect,
            "seed={seed}, background={background:?}, scoring={scoring:?}, all windows, {opt:?}"
        );

        let (hits, best, counts, dense) = program
            .scan(
                seqs,
                (
                    Hits::pvalue(p),
                    SequenceBest::new(),
                    SequenceCounts::new(Threshold::PValue(p)),
                    Dense::new(),
                ),
                opt,
            )
            .unwrap();
        let variant = |mi: u32, strand: Strand, c: u8| {
            program
                .motif(mi as usize)
                .strand(strand)
                .unwrap()
                .variant(c as usize)
        };
        let mut want_hits: Vec<_> = expect
            .iter()
            .filter(|w| variant(w.0, w.1, w.4).meets(w.5, Threshold::PValue(p)))
            .map(|w| (w.0, w.2, w.3, w.1))
            .collect();
        want_hits.sort();
        let got_hits: Vec<_> = hits
            .iter()
            .map(|h| (h.motif, h.sequence, h.start, h.strand))
            .collect();
        assert_eq!(
            got_hits, want_hits,
            "seed={seed}, background={background:?}, scoring={scoring:?}, hits, {opt:?}"
        );

        // The ordered stage sees the same windows, by sequence and position,
        // with every non-empty sequence covered exactly once, in order.
        let ordered = program
            .scan(
                seqs,
                Ordered {
                    threshold: Threshold::PValue(p),
                    chunks: Vec::new(),
                    matches: Vec::new(),
                },
                opt,
            )
            .unwrap();
        let mut want_ordered: Vec<_> = expect
            .iter()
            .filter(|w| variant(w.0, w.1, w.4).meets(w.5, Threshold::PValue(p)))
            .map(|w| (w.2, w.0, w.3, w.1, w.5))
            .collect();
        want_ordered.sort_by_key(|w| (w.0, w.2, w.1, w.3));
        assert_eq!(
            ordered.matches, want_ordered,
            "seed={seed}, background={background:?}, scoring={scoring:?}, ordered, {opt:?}"
        );
        let mut covered = 0u32;
        let mut last_seq = None;
        for &(s, a, b, first, last) in &ordered.chunks {
            if last_seq != Some(s) {
                assert!(first && a == 0, "chunk starts a sequence");
                covered = 0;
                last_seq = Some(s);
            }
            assert_eq!(a, covered);
            covered = b;
            assert_eq!(last, b as usize == seqs.get(s as usize).len());
        }
        let nonempty = (0..seqs.len()).filter(|&i| !seqs.get(i).is_empty()).count();
        assert_eq!(
            ordered.chunks.iter().filter(|c| c.4).count(),
            nonempty,
            "every sequence ends once"
        );
        for h in &hits {
            assert_eq!(h.end - h.start, program.motif(h.motif as usize).width());
        }

        // Best window per (sequence, motif).
        let mut want_best: std::collections::BTreeMap<(u32, u32), (f64, u32, Strand)> =
            Default::default();
        for w in &expect {
            let raw = variant(w.0, w.1, w.4).int_model().raw(w.5);
            let e = want_best.entry((w.2, w.0)).or_insert((raw, w.3, w.1));
            if raw > e.0 || (raw == e.0 && (w.3, w.1) < (e.1, e.2)) {
                *e = (raw, w.3, w.1);
            }
        }
        let got_best: Vec<_> = best
            .iter()
            .map(|b| ((b.sequence, b.motif), (b.start, b.strand)))
            .collect();
        let want_best: Vec<_> = want_best
            .into_iter()
            .map(|(k, v)| (k, (v.1, v.2)))
            .collect();
        assert_eq!(
            got_best, want_best,
            "seed={seed}, background={background:?}, scoring={scoring:?}, best, {opt:?}"
        );

        // Counts.
        let mut want_counts = std::collections::BTreeMap::<(u32, u32), u32>::new();
        for h in &want_hits {
            *want_counts.entry((h.1, h.0)).or_default() += 1;
        }
        let mut got_counts = Vec::new();
        for s in 0..seqs.len() {
            for j in counts.indptr[s] as usize..counts.indptr[s + 1] as usize {
                got_counts.push(((s as u32, counts.motifs[j]), counts.counts[j]));
            }
        }
        assert_eq!(
            got_counts,
            want_counts.into_iter().collect::<Vec<_>>(),
            "seed={seed}, counts"
        );

        // Dense scores at every reported window, NaN elsewhere.
        let mut filled = 0usize;
        for w in &expect {
            let len = seqs.get(w.2 as usize).len() as u64;
            let strand = if w.1 == Strand::Reverse { 1 } else { 0 };
            let at = dense.offsets[w.2 as usize] + (w.0 as u64 * 2 + strand) * len + w.3 as u64;
            let v = dense.values[at as usize];
            assert_eq!(v, variant(w.0, w.1, w.4).score(w.5) as f32);
            filled += 1;
        }
        assert_eq!(dense.values.iter().filter(|v| !v.is_nan()).count(), filled);
    }
}
