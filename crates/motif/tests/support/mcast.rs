pub(crate) use crate::rng::Rng;
use omicsweave_ffm::{Alphabet, FfmReader, FfmWriter, Motif, MotifRef, Schema, WriterOptions};
use omicsweave_motif::mcast;
use omicsweave_motif::{Background, Program, Strands};

pub(crate) fn program(motifs: &[Motif]) -> Program {
    let alphabet = Alphabet::dna();
    let mut w = FfmWriter::new(
        Vec::new(),
        alphabet.clone(),
        Schema::default(),
        WriterOptions::default(),
    )
    .unwrap();
    for (i, m) in motifs.iter().enumerate() {
        w.push(&format!("M{i}"), m, &[]).unwrap();
    }
    let reader = FfmReader::open(w.finish().unwrap()).unwrap();
    let refs: Vec<MotifRef> = (0..reader.len())
        .map(|i| reader.motif(i).unwrap())
        .collect();
    let views: Vec<_> = refs.iter().map(|r| r.view()).collect();
    let bg = Background::nrdb(&alphabet).unwrap();
    let scoring = mcast::scoring(&bg, &alphabet).unwrap();
    Program::compile(&views, &alphabet, &bg, &scoring, Strands::Both, Some(1)).unwrap()
}
