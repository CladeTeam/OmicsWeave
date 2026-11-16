use omicsweave_convert::{Conversion, Format, Options};

#[test]
fn sidecars_are_only_options_of_bbi_and_hic_outputs() {
    let options = Options {
        sidecar: true,
        ..Options::default()
    };
    for (from, to) in [
        (Format::Sam, Format::Bam),
        (Format::Vcf, Format::Bcf),
        (Format::Fasta, Format::BgzipFasta),
        (Format::Fasta, Format::TwoBit),
    ] {
        let mut conversion = Conversion::between(from, to, "input", "output").unwrap();
        let error = conversion.configure(&options).unwrap_err();
        assert!(error.to_string().contains("sidecar"), "{error}");
    }
    for (from, to) in [
        (Format::Bed, Format::BigBed),
        (Format::BedGraph, Format::BigWig),
        (Format::Wig, Format::BigWig),
        (Format::BedMethyl, Format::BigBed),
        (Format::BedMethyl, Format::BigWig),
        (Format::Cool, Format::Hic),
        (Format::Pairs, Format::Hic),
    ] {
        let mut conversion = Conversion::between(from, to, "input", "output").unwrap();
        conversion.configure(&options).unwrap();
    }
}
