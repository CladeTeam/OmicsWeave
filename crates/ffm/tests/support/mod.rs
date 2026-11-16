use omicsweave_ffm::{
    Alphabet, Factor, FactorKind, FfmWriter, Motif, Schema, Support, WriterOptions,
};

pub mod wire;

pub fn pwm(width: usize, seed: u64) -> Motif {
    let rows = (0..width)
        .map(|p| {
            (0..4)
                .map(|b| ((seed + p as u64 * 7 + b) % 13) as f64 + 0.5)
                .collect()
        })
        .collect::<Vec<Vec<f64>>>();
    Motif::pfm(&rows, Some(20.0))
}

pub fn write(
    alphabet: Alphabet,
    schema: Schema,
    opts: WriterOptions,
    motifs: &[(String, Motif, Vec<String>)],
) -> Vec<u8> {
    let mut w = FfmWriter::new(Vec::new(), alphabet, schema, opts).unwrap();
    for (id, m, ann) in motifs {
        let ann: Vec<&str> = ann.iter().map(String::as_str).collect();
        w.push(id, m, &ann).unwrap();
    }
    w.finish().unwrap()
}

pub fn schema() -> Schema {
    Schema::new(vec![
        omicsweave_ffm::Key::new("name", omicsweave_ffm::KeyType::String),
        omicsweave_ffm::Key::new("n", omicsweave_ffm::KeyType::Integer),
        omicsweave_ffm::Key::new("p", omicsweave_ffm::KeyType::Float),
    ])
    .unwrap()
}

pub fn pfm() -> Motif {
    Motif::pfm(
        &[vec![1.0, 2.0, 3.0, 4.0], vec![4.0, 3.0, 2.0, 1.0]],
        Some(10.0),
    )
}

pub fn pair() -> Motif {
    Motif {
        width: 3,
        support: Support::Codes(vec![0]),
        sites: None,
        factors: vec![Factor {
            kind: FactorKind::Conditional,
            positions: vec![0, 2],
            table: vec![1.0, 2.0, 3.0, 4.0],
        }],
        gaps: None,
    }
}
