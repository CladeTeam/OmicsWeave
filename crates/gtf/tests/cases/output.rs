//! Lossless and canonical GTF output, compression and BED12 geometry.
use crate::support::{
    fixtures::*,
    generator::{Rng, random_gtf},
    gtf_text, read,
};
use omicsweave_gtf::{Annotation, AttributePolicy, Dialect, GtfLayout, ReadOptions, Row};
use std::io::Write;

#[test]
fn gtf_round_trips_byte_for_byte() {
    let odd = "c\ts\tgene\t1\t9\t5\t+\t1\tgene_id \"g\";\n#mid\nc\ts\ttranscript\t1\t9\t.\t?\t2\tgene_id \"g\"; transcript_id \"t\";\n";
    for text in [GENCODE_GTF, ENSEMBL_GTF, REFSEQ_GTF, UCSC_GTF, odd] {
        let ann = read(text);
        assert!(ann.has_original_gtf());
        assert_eq!(gtf_text(&ann, GtfLayout::Original), text);
        assert_eq!(gtf_text(&ann, GtfLayout::Auto), text);
    }
}

#[test]
fn model_layout_reads_back_to_the_same_model() {
    for text in [
        GENCODE_GTF,
        ENSEMBL_GTF,
        REFSEQ_GTF,
        UCSC_GTF,
        ENSEMBL_GFF3,
        REFSEQ_GFF3,
    ] {
        let a = read(text);
        let model = gtf_text(&a, GtfLayout::Model);
        let b = read(&model);
        assert_eq!(b.dialect(), Dialect::Gtf);
        assert_eq!(a.genes().len(), b.genes().len());
        assert_eq!(a.transcripts().len(), b.transcripts().len());
        assert_eq!(a.features().len(), b.features().len());
        for i in 0..a.genes().len() {
            let (ga, gb) = (a.genes(), b.genes());
            assert_eq!(a.str(ga.id[i]), b.str(gb.id[i]));
            assert_eq!(a.str(ga.name[i]), b.str(gb.name[i]));
            assert_eq!(
                (ga.start[i], ga.end[i], ga.strand[i]),
                (gb.start[i], gb.end[i], gb.strand[i])
            );
        }
        assert_eq!(a.transcripts().start, b.transcripts().start);
        assert_eq!(a.transcripts().end, b.transcripts().end);
        assert_eq!(a.transcripts().strand, b.transcripts().strand);
        assert_eq!(a.transcripts().gene, b.transcripts().gene);
        assert_eq!(a.transcripts().tag_bits, b.transcripts().tag_bits);
        for i in 0..a.transcripts().len() {
            assert_eq!(
                a.str(a.transcripts().id[i]),
                b.str(b.transcripts().id[i]),
                "transcript {i}"
            );
            assert_eq!(
                a.str(a.transcripts().name[i]),
                b.str(b.transcripts().name[i]),
                "transcript {i}"
            );
        }
        assert_eq!(a.features().start, b.features().start);
        assert_eq!(a.features().end, b.features().end);
        assert_eq!(a.features().strand, b.features().strand);
        assert_eq!(a.features().phase, b.features().phase);
        assert_eq!(a.features().transcript, b.features().transcript);
        assert_eq!(a.features().kind, b.features().kind);
        // And writing the re-read model again is a fixed point.
        assert_eq!(gtf_text(&b, GtfLayout::Original), model);
    }
}

#[test]
fn compression_is_told_by_content() {
    let text = GENCODE_GTF.as_bytes();
    let plain = read(GENCODE_GTF);
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gz.write_all(text).unwrap();
    let gz = gz.finish().unwrap();
    // BGZF: concatenated gzip members with the BC extra field.
    let mut bgzf = Vec::new();
    for chunk in text.chunks(700) {
        let mut out = Vec::new();
        omicsweave_bgzf::deflate_block_at(6, chunk, &mut out).unwrap();
        bgzf.extend_from_slice(&out);
    }
    bgzf.extend_from_slice(&omicsweave_bgzf::EOF_BLOCK);
    let dir = tempfile::tempdir().unwrap();
    for (name, bytes) in [
        ("a.gtf", text.to_vec()),
        ("a.gtf.gz", gz),
        ("b.gtf.gz", bgzf),
    ] {
        let path = dir.path().join(name);
        std::fs::write(&path, &bytes).unwrap();
        for threads in [1, 3] {
            let opts = ReadOptions {
                threads,
                batch_bytes: 512,
                ..ReadOptions::default()
            };
            let a = Annotation::read(&path, &opts).unwrap();
            assert_eq!(gtf_text(&a, GtfLayout::Original), GENCODE_GTF, "{name}");
            let b = Annotation::parse(&bytes, &opts).unwrap();
            assert_eq!(b.genes().start, plain.genes().start);
        }
    }
}

#[test]
fn output_does_not_depend_on_threads_or_batches() {
    let mut rng = Rng(7);
    let text = random_gtf(&mut rng, 4, 300, 20_000);
    let base = Annotation::parse(
        text.as_bytes(),
        &ReadOptions {
            threads: 1,
            batch_bytes: 1 << 30,
            ..ReadOptions::default()
        },
    )
    .unwrap();
    let want = (
        gtf_text(&base, GtfLayout::Original),
        gtf_text(&base, GtfLayout::Model),
    );
    assert_eq!(want.0, text);
    for (threads, batch) in [(1, 300), (4, 97), (8, 4096)] {
        let a = Annotation::parse(
            text.as_bytes(),
            &ReadOptions {
                threads,
                batch_bytes: batch,
                ..ReadOptions::default()
            },
        )
        .unwrap();
        assert_eq!(
            (
                gtf_text(&a, GtfLayout::Original),
                gtf_text(&a, GtfLayout::Model)
            ),
            want
        );
        assert_eq!(a.genes().id.len(), base.genes().id.len());
        // String ids are assigned in file order too.
        assert_eq!(a.genes().id, base.genes().id);
        assert_eq!(a.features().record, base.features().record);
    }
}

#[test]
fn attribute_policies() {
    let only = ReadOptions {
        attributes: AttributePolicy::Only(vec!["exon_id".into()]),
        ..ReadOptions::default()
    };
    let a = Annotation::parse(GENCODE_GTF.as_bytes(), &only).unwrap();
    // The hierarchy and names do not depend on kept attributes.
    assert_eq!(a.str(a.genes().name[2]), "OR4F5");
    let f = a.transcripts().features(0).start;
    let kept: Vec<(&str, &str)> = a.attributes(Row::Feature(f)).collect();
    assert_eq!(kept, [("exon_id", "ENSE00002234944.1")]);
    assert!(!a.has_original_gtf());
    assert!(a.write_gtf(&mut Vec::new(), GtfLayout::Original).is_err());
    let none = ReadOptions {
        attributes: AttributePolicy::None,
        ..ReadOptions::default()
    };
    let n = Annotation::parse(GENCODE_GTF.as_bytes(), &none).unwrap();
    assert!(n.attribute_table().is_empty());
    assert!(n.heap_bytes() < read(GENCODE_GTF).heap_bytes());
    // The model layout falls back to the model's names and tags.
    let model = gtf_text(&n, GtfLayout::Model);
    assert!(model.contains(
        "transcript_name \"OR4F5-201\"; transcript_type \"protein_coding\"; tag \"basic\";"
    ));
}

#[test]
fn bed12_coding_span_stays_within_exon_bounds() {
    for (cds_start, cds_end, expected) in [(21, 30, (10, 10)), (1, 5, (4, 5)), (1, 3, (10, 10))] {
        let text = format!(
            "c\t.\texon\t5\t10\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\nc\t.\tCDS\t{cds_start}\t{cds_end}\t.\t+\t0\tgene_id \"g\"; transcript_id \"t\";\n"
        );
        let ann = read(&text);
        let mut bed = Vec::new();
        ann.write_bed12(&mut bed, None).unwrap();
        let text = String::from_utf8(bed).unwrap();
        let columns: Vec<_> = text.trim_end().split('\t').collect();
        assert_eq!(
            (
                columns[6].parse::<u32>().unwrap(),
                columns[7].parse::<u32>().unwrap()
            ),
            expected
        );
    }
}

#[test]
fn bed12_of_transcripts() {
    let ann = read(GENCODE_GTF);
    let mut out = Vec::new();
    let tx = ann.transcript("ENST00000641515.2").unwrap();
    ann.write_bed12(
        &mut out,
        Some(&[tx, ann.transcript("ENST00000488147.1").unwrap()]),
    )
    .unwrap();
    let text = String::from_utf8(out).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(
        lines[0],
        "chr1\t65418\t71585\tENST00000641515.2\t0\t+\t65564\t70008\t0\t3\t15,54,2549,\t0,101,3618,"
    );
    assert_eq!(
        lines[1],
        "chr1\t14403\t29570\tENST00000488147.1\t0\t-\t29570\t29570\t0\t3\t98,154,37,\t0,10334,15130,"
    );
}

#[test]
fn exports_propagate_sink_errors_and_allow_retry_with_a_fresh_sink() {
    struct FailingSink;
    impl Write for FailingSink {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("annotation output unavailable"))
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let ann = read(GENCODE_GTF);
    for layout in [GtfLayout::Original, GtfLayout::Model, GtfLayout::Auto] {
        let err = ann.write_gtf(&mut FailingSink, layout).unwrap_err();
        assert!(
            err.to_string().contains("annotation output unavailable"),
            "layout={layout:?}: {err}"
        );
        let text = gtf_text(&ann, layout);
        let parsed = read(&text);
        assert_eq!(parsed.genes().len(), 4);
        assert_eq!(parsed.transcripts().len(), 4);
    }
    let err = ann.write_bed12(&mut FailingSink, None).unwrap_err();
    assert!(
        err.to_string().contains("annotation output unavailable"),
        "{err}"
    );
    let mut bed = Vec::new();
    ann.write_bed12(&mut bed, Some(&[])).unwrap();
    assert!(bed.is_empty());
    ann.write_bed12(&mut bed, None).unwrap();
    assert_eq!(String::from_utf8(bed).unwrap().lines().count(), 4);
}
