//! Sparse text conversions preserve coordinates beyond the TBI limit.

use omicsweave_bbi::BbiReader;
use omicsweave_convert::{BedGraphToBigWig, BedToBigBed, Budget, ChromSizes, Conversion, convert};
use omicsweave_core::{
    FileSource, IntervalRead, IntervalsOut, Locus, Region, StringColumn, TrackRead,
};

#[test]
fn sparse_bigwig_and_bigbed_preserve_long_chromosome_ends() {
    let dir = tempfile::tempdir().unwrap();
    for len in [830_000_000u32, u32::MAX] {
        let sizes = dir.path().join("chrom.sizes");
        std::fs::write(&sizes, format!("long\t{len}\n")).unwrap();
        let start = len - 10;
        let region = Region::new("long", start, len);
        let bg = dir.path().join("in.bedgraph");
        std::fs::write(&bg, format!("long\t{start}\t{len}\t3.5\n")).unwrap();
        let bw = dir.path().join("out.bw");
        convert(
            Conversion::BedGraphToBigWig(BedGraphToBigWig::new(
                &bg,
                &bw,
                ChromSizes::File(sizes.clone()),
            )),
            Budget::new(2, 64 << 20),
            None,
        )
        .unwrap();
        let mut reader = BbiReader::open(FileSource::open(&bw).unwrap()).unwrap();
        assert_eq!(reader.chrom_len("long").unwrap(), len);
        let mut values = [0.0; 10];
        reader.values(&Locus::from(region), &mut values).unwrap();
        assert_eq!(values, [3.5; 10]);

        let bed = dir.path().join("in.bed");
        std::fs::write(&bed, format!("long\t{start}\t{len}\tfeature\n")).unwrap();
        let bb = dir.path().join("out.bb");
        convert(
            Conversion::BedToBigBed(BedToBigBed::new(&bed, &bb, ChromSizes::File(sizes))),
            Budget::new(2, 64 << 20),
            None,
        )
        .unwrap();
        let mut reader = BbiReader::open(FileSource::open(&bb).unwrap()).unwrap();
        let size = reader.query(region, false).unwrap();
        assert_eq!(size.entries, 1);
        let mut entries = vec![Default::default(); size.entries];
        let mut offsets = vec![0; size.entries + 1];
        let mut data = vec![0; size.annotation_bytes];
        reader
            .write(IntervalsOut {
                entries: &mut entries,
                annotations: StringColumn {
                    offsets: &mut offsets,
                    data: &mut data,
                },
            })
            .unwrap();
        assert_eq!((entries[0].start, entries[0].end), (start, len));
        assert_eq!(data, b"feature");
    }
}
