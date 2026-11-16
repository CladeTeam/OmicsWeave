use omicsweave_ffm::{
    Alphabet, Factor, FactorKind, FfmReader, FfmWriter, Motif, Schema, Support, WriterOptions,
};

pub(crate) use crate::rng::Rng;

impl Rng {
    pub(crate) fn table(&mut self, len: usize, zeros: bool) -> Vec<f64> {
        (0..len)
            .map(|_| {
                if zeros && self.below(6) == 0 {
                    0.0
                } else {
                    (self.unit() * 100.0).round() + 1.0
                }
            })
            .collect()
    }
}

pub(crate) fn pwm(rng: &mut Rng, width: usize) -> Motif {
    let rows: Vec<Vec<f64>> = (0..width).map(|_| rng.table(4, true)).collect();
    Motif::pfm(&rows, Some(10.0 + rng.below(40) as f64))
}

pub(crate) fn dinucleotide(rng: &mut Rng, width: usize) -> Motif {
    Motif {
        width: width as u32,
        support: Support::All,
        sites: None,
        factors: (0..width - 1)
            .map(|i| Factor {
                kind: FactorKind::Joint,
                positions: vec![i as u32, i as u32 + 1],
                table: rng.table(16, true),
            })
            .collect(),
        gaps: None,
    }
}

pub(crate) fn markov(rng: &mut Rng, width: usize) -> Motif {
    Motif {
        width: width as u32,
        support: Support::All,
        sites: Some(30.0),
        factors: (0..width)
            .map(|i| {
                let positions: Vec<u32> = (i.saturating_sub(2)..=i).map(|p| p as u32).collect();
                let len = 4usize.pow(positions.len() as u32);
                Factor {
                    kind: FactorKind::Conditional,
                    positions,
                    table: rng.table(len, false),
                }
            })
            .collect(),
        gaps: None,
    }
}

pub(crate) fn gapped(rng: &mut Rng) -> Motif {
    let mut m = pwm(rng, 10);
    m.factors.push(Factor {
        kind: FactorKind::Potential,
        positions: vec![0, 9],
        table: (0..16).map(|_| 0.5 + rng.unit()).collect(),
    });
    m
}

pub(crate) fn subset(rng: &mut Rng) -> Motif {
    // Models A and G individually; C and T form the rest.
    Motif {
        width: 6,
        support: Support::Codes(vec![0, 2]),
        sites: Some(15.0),
        factors: (0..6)
            .map(|i| Factor {
                kind: FactorKind::Joint,
                positions: vec![i],
                table: rng.table(3, true),
            })
            .collect(),
        gaps: None,
    }
}

pub(crate) fn kmer_potential(rng: &mut Rng) -> Motif {
    Motif {
        width: 4,
        support: Support::All,
        sites: None,
        factors: vec![Factor {
            kind: FactorKind::Potential,
            positions: vec![0, 1, 2, 3],
            table: (0..256).map(|_| rng.unit() * 3.0).collect(),
        }],
        gaps: None,
    }
}

pub(crate) fn reader(motifs: &[Motif]) -> FfmReader<Vec<u8>> {
    let mut w = FfmWriter::new(
        Vec::new(),
        Alphabet::dna(),
        Schema::default(),
        WriterOptions::default(),
    )
    .unwrap();
    for (i, m) in motifs.iter().enumerate() {
        w.push(&format!("m{i}"), m, &[]).unwrap();
    }
    FfmReader::open(w.finish().unwrap()).unwrap()
}

pub(crate) fn sequences(rng: &mut Rng, count: usize) -> (Vec<u8>, Vec<i64>) {
    let mut data = Vec::new();
    let mut offsets = vec![0i64];
    for i in 0..count {
        let len = match i {
            0 => 0,
            1 => 3,
            2 => 20_000,
            _ => rng.below(1500),
        };
        for _ in 0..len {
            let b = match rng.below(100) {
                0 | 1 => b'N',
                2 => b'R',
                3 => b'-',
                4..=12 => b"acgt"[rng.below(4)],
                _ => b"ACGT"[rng.below(4)],
            };
            data.push(b);
        }
        offsets.push(data.len() as i64);
    }
    (data, offsets)
}

pub(crate) fn motif_zoo(rng: &mut Rng) -> Vec<Motif> {
    vec![
        pwm(rng, 1),
        pwm(rng, 4),
        pwm(rng, 12),
        pwm(rng, 21),
        dinucleotide(rng, 7),
        markov(rng, 9),
        gapped(rng),
        subset(rng),
        kmer_potential(rng),
    ]
}
