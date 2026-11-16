use omicsweave_ffm::{Alphabet, FfmReader, FfmWriter, Motif, Schema, WriterOptions};

pub(crate) fn resident(motif: &Motif) -> FfmReader<Vec<u8>> {
    let mut writer = FfmWriter::new(
        Vec::new(),
        Alphabet::dna(),
        Schema::default(),
        WriterOptions::default(),
    )
    .unwrap();
    writer.push("m", motif, &[]).unwrap();
    FfmReader::open(writer.finish().unwrap()).unwrap()
}
