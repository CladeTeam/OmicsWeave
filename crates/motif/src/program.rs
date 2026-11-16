//! A compiled motif collection: every motif's models and lowered blocks,
//! ready to scan sequences.

use std::ops::Range;
use std::sync::atomic::{AtomicUsize, Ordering};

use omicsweave_core::{Error, Result, Strand, Strings};
use omicsweave_ffm::{Alphabet, MotifView};

use crate::background::Background;
use crate::lower::{Lowered, Packing, lower};
use crate::model::{MotifModel, check_strands};
use crate::reduce::Reducer;
use crate::scan::{ScanOptions, run};
use crate::scoring::{Scoring, Strands};

type Compiled = Result<(MotifModel, Vec<Lowered>)>;

/// Motifs compiled under one background, scoring policy and strand choice.
#[derive(Clone, Debug)]
pub struct Part {
    background: Background,
    scoring: Scoring,
    strands: Strands,
    motifs: Range<usize>,
}

impl Part {
    /// The background after the scoring policy's preparation.
    pub fn background(&self) -> &Background {
        &self.background
    }

    pub fn scoring(&self) -> &Scoring {
        &self.scoring
    }

    pub fn strands(&self) -> Strands {
        self.strands
    }

    /// The program's motif indices that belong to this part.
    pub fn motifs(&self) -> Range<usize> {
        self.motifs.clone()
    }
}

/// Motifs compiled for scanning: one or more [`Part`]s over one alphabet,
/// scanned together so each stretch of sequence is encoded once for all.
#[derive(Clone, Debug)]
pub struct Program {
    alphabet: Alphabet,
    packing: Packing,
    parts: Vec<Part>,
    /// Part of each motif.
    part_of: Vec<u32>,
    motifs: Vec<MotifModel>,
    /// Per motif, per scanned strand.
    lowered: Vec<Vec<Lowered>>,
    /// Ambiguity code index of each byte (`0xFF` for none), shared by every
    /// part that scores ambiguity codes.
    ambiguity_lut: [u8; 256],
}

/// The byte table of the parts' ambiguity codes; parts that score them must
/// agree on the codes.
fn ambiguity_lut(alphabet: &Alphabet, parts: &[Part]) -> Result<[u8; 256]> {
    let mut shared: Option<&crate::scoring::Ambiguity> = None;
    for p in parts {
        let a = &p.scoring.ambiguity;
        if a.rule == crate::scoring::AmbiguityRule::Skip {
            continue;
        }
        match shared {
            Some(s) if s.codes != a.codes || s.unknown != a.unknown => {
                return Err(Error::invalid(
                    "motif program",
                    "parts that score ambiguity codes must share them",
                ));
            }
            _ => shared = Some(a),
        }
    }
    let mut lut = [0xFFu8; 256];
    let Some(a) = shared else { return Ok(lut) };
    for (k, (b, _)) in a.codes.iter().enumerate() {
        lut[*b as usize] = k as u8;
    }
    if let Some(u) = a.unknown {
        let k = a
            .codes
            .iter()
            .position(|(b, _)| *b == u)
            .expect("checked code") as u8;
        for (byte, slot) in lut.iter_mut().enumerate() {
            if *slot == 0xFF && alphabet.lut()[byte] == omicsweave_ffm::INVALID {
                *slot = k;
            }
        }
    }
    Ok(lut)
}

impl Program {
    /// Compiles `motifs` (steps ①–⑥). The background is prepared by
    /// `scoring`; `threads = None` uses the available parallelism.
    pub fn compile(
        motifs: &[MotifView<'_>],
        alphabet: &Alphabet,
        background: &Background,
        scoring: &Scoring,
        strands: Strands,
        threads: Option<usize>,
    ) -> Result<Program> {
        scoring.check()?;
        scoring.ambiguity.check(alphabet)?;
        check_strands(alphabet, strands)?;
        if alphabet.len() > 255 {
            return Err(Error::invalid("motif scan", "alphabet larger than 255"));
        }
        let background = scoring.prepare_background(background, alphabet, strands)?;
        if background.order() > 32 {
            return Err(Error::invalid("motif background", "order above 32"));
        }
        let packing = Packing::new(alphabet.len());
        let compile_one = |i: usize| -> Compiled {
            let model = MotifModel::build(&motifs[i], alphabet, &background, scoring, strands)
                .map_err(|e| in_motif(i, e))?;
            let lowered = model
                .strands()
                .iter()
                .map(|s| lower(s.full().int_model(), &packing))
                .collect::<Result<Vec<_>>>()
                .map_err(|e| in_motif(i, e))?;
            Ok((model, lowered))
        };
        let threads = threads
            .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()))
            .clamp(1, motifs.len().max(1));
        let compiled: Vec<Compiled> = if threads == 1 {
            (0..motifs.len()).map(compile_one).collect()
        } else {
            let next = AtomicUsize::new(0);
            let mut parts: Vec<Vec<(usize, Compiled)>> = std::thread::scope(|scope| {
                let handles: Vec<_> = (0..threads)
                    .map(|_| {
                        scope.spawn(|| {
                            let mut out = Vec::new();
                            loop {
                                let i = next.fetch_add(1, Ordering::Relaxed);
                                if i >= motifs.len() {
                                    break;
                                }
                                out.push((i, compile_one(i)));
                            }
                            out
                        })
                    })
                    .collect();
                handles
                    .into_iter()
                    .map(|h| h.join().expect("compile worker panicked"))
                    .collect()
            });
            let mut all: Vec<(usize, Compiled)> = parts.drain(..).flatten().collect();
            all.sort_by_key(|(i, _)| *i);
            all.into_iter().map(|(_, r)| r).collect()
        };
        let mut models = Vec::with_capacity(motifs.len());
        let mut lowered = Vec::with_capacity(motifs.len());
        for r in compiled {
            let (m, l) = r?;
            models.push(m);
            lowered.push(l);
        }
        let parts = vec![Part {
            background,
            scoring: scoring.clone(),
            strands,
            motifs: 0..models.len(),
        }];
        Ok(Program {
            alphabet: alphabet.clone(),
            packing,
            ambiguity_lut: ambiguity_lut(alphabet, &parts)?,
            parts,
            part_of: vec![0; models.len()],
            motifs: models,
            lowered,
        })
    }

    /// One program scanning every motif of `programs`, which must share an
    /// alphabet. Part `i` of the result is the `i`-th part met; motif
    /// indices follow the same order.
    pub fn combine(programs: impl IntoIterator<Item = Program>) -> Result<Program> {
        let mut iter = programs.into_iter();
        let mut out = iter
            .next()
            .ok_or_else(|| Error::invalid("motif program", "nothing to combine"))?;
        for p in iter {
            if p.alphabet != out.alphabet {
                return Err(Error::invalid(
                    "motif program",
                    "combined programs must share an alphabet",
                ));
            }
            let base = out.motifs.len();
            let part_base = out.parts.len() as u32;
            for mut part in p.parts {
                part.motifs = part.motifs.start + base..part.motifs.end + base;
                out.parts.push(part);
            }
            out.part_of.extend(p.part_of.iter().map(|&q| q + part_base));
            out.motifs.extend(p.motifs);
            out.lowered.extend(p.lowered);
        }
        out.ambiguity_lut = ambiguity_lut(&out.alphabet, &out.parts)?;
        Ok(out)
    }

    /// Scans `sequences` (ASCII bytes, one row per sequence) once for every
    /// motif and returns what `reducer` keeps.
    pub fn scan<R: Reducer>(
        &self,
        sequences: Strings<'_>,
        reducer: R,
        options: &ScanOptions,
    ) -> Result<R::Output> {
        run(self, sequences, reducer, options)
    }

    /// Scans one sequence.
    pub fn scan_sequence<R: Reducer>(
        &self,
        sequence: &[u8],
        reducer: R,
        options: &ScanOptions,
    ) -> Result<R::Output> {
        let offsets = [0i64, sequence.len() as i64];
        self.scan(
            Strings {
                offsets: &offsets,
                data: sequence,
            },
            reducer,
            options,
        )
    }

    pub fn len(&self) -> usize {
        self.motifs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.motifs.is_empty()
    }

    pub fn motif(&self, i: usize) -> &MotifModel {
        &self.motifs[i]
    }

    pub fn alphabet(&self) -> &Alphabet {
        &self.alphabet
    }

    pub fn parts(&self) -> &[Part] {
        &self.parts
    }

    /// The index of the part motif `motif` belongs to.
    pub fn part_index(&self, motif: usize) -> usize {
        self.part_of[motif] as usize
    }

    /// The part motif `motif` belongs to.
    pub fn part_of(&self, motif: usize) -> &Part {
        &self.parts[self.part_of[motif] as usize]
    }

    /// The largest background order of any part.
    pub(crate) fn context(&self) -> usize {
        self.parts
            .iter()
            .map(|p| p.background.order())
            .max()
            .unwrap_or(0)
    }

    /// Strands any part scans, forward first.
    pub fn strand_list(&self) -> Vec<Strand> {
        let mut out = Vec::with_capacity(2);
        if self.parts.iter().any(|p| p.strands.forward()) {
            out.push(Strand::Forward);
        }
        if self.parts.iter().any(|p| p.strands.reverse()) {
            out.push(Strand::Reverse);
        }
        out
    }

    pub(crate) fn ambiguity_lut(&self) -> &[u8; 256] {
        &self.ambiguity_lut
    }

    pub(crate) fn packing(&self) -> &Packing {
        &self.packing
    }

    pub(crate) fn lowered(&self, i: usize) -> &[Lowered] {
        &self.lowered[i]
    }
}

fn in_motif(i: usize, e: Error) -> Error {
    match e {
        Error::Invalid { what, detail } => Error::Invalid {
            what,
            detail: format!("motif {i}: {detail}"),
        },
        other => other,
    }
}
