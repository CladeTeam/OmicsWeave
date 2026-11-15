//! Reading each source's dialect, round trips, compression, determinism
//! across threads, validation and aliases.

mod common;

use std::io::Write;

use common::*;
use omicsweave_core::Region;
use omicsweave_gtf::{
    Aliases, Annotation, AttributePolicy, Dialect, FeatureKind, GtfLayout, Level, ProblemKind,
    ReadOptions, Row, Strand, Tags,
};

fn gtf_text(ann: &Annotation, layout: GtfLayout) -> String {
    let mut out = Vec::new();
    ann.write_gtf(&mut out, layout).unwrap();
    String::from_utf8(out).unwrap()
}

#[test]
fn gff3_fasta_tail_is_not_parsed_as_annotation_across_batches() {
    let mut text = format!(
        "##gff-version 3\n#{}\nchr1\t.\tgene\t1\t20\t.\t+\t.\tID=g\n##FASTA\n>chr1\nACGTACGTACGTACGTACGT\ninvalid annotation text\n",
        "long comment".repeat(32_768)
    ).into_bytes();
    text.extend_from_slice(b"\0\xff\xfe\n");
    for threads in [1, 2, 4] {
        for batch_bytes in [1, 64, 1024, 1 << 20] {
            let ann = Annotation::parse(
                &text,
                &ReadOptions {
                    threads,
                    batch_bytes,
                    ..ReadOptions::default()
                },
            )
            .unwrap();
            assert_eq!(ann.genes().len(), 1);
            assert_eq!(ann.str(ann.genes().id[0]), "g");
        }
    }
    let invalid = b"##gff-version 3\ninvalid annotation text\n##FASTA\n>chr1\nACGT\n";
    for threads in [1, 4] {
        assert!(
            Annotation::parse(
                invalid,
                &ReadOptions {
                    threads,
                    batch_bytes: 1,
                    ..ReadOptions::default()
                }
            )
            .is_err()
        );
    }
}

fn gene_by_id(ann: &Annotation, id: &str) -> u32 {
    ann.gene(id).unwrap_or_else(|| panic!("no gene {id}"))
}

#[test]
fn gencode_hierarchy_names_and_tags() {
    let ann = read(GENCODE_GTF);
    assert_eq!(ann.dialect(), Dialect::Gtf);
    assert!(ann.report().is_clean(), "{:?}", ann.report().problems);
    let (g, t, f) = (ann.genes(), ann.transcripts(), ann.features());
    assert_eq!((g.len(), t.len()), (4, 4));
    assert_eq!(ann.report().synthesized_genes, 0);
    // Genes sorted by seqid then start.
    let ids: Vec<&str> = g.id.iter().map(|&i| ann.str(i)).collect();
    assert_eq!(
        ids,
        [
            "ENSG00000290825.1",
            "ENSG00000227232.5",
            "ENSG00000186092.7",
            "ENSG00000198888.2"
        ]
    );
    let or4f5 = gene_by_id(&ann, "ENSG00000186092.7");
    assert_eq!(ann.str(g.name[or4f5 as usize]), "OR4F5");
    assert_eq!(ann.str(g.biotype[or4f5 as usize]), "protein_coding");
    assert_eq!(
        (g.start[or4f5 as usize], g.end[or4f5 as usize]),
        (65418, 71585)
    );
    let tx = g.transcripts(or4f5).start;
    assert_eq!(ann.str(t.name[tx as usize]), "OR4F5-201");
    let bits = t.tag_bits[tx as usize];
    assert!(bits.contains(Tags::BASIC | Tags::MANE_SELECT | Tags::CCDS | Tags::APPRIS_PRINCIPAL));
    assert_eq!(t.tags(tx).len(), 5);
    // Features sorted by (start, end, kind).
    let kinds: Vec<FeatureKind> = t.features(tx).map(|i| f.kind[i as usize]).collect();
    assert_eq!(kinds[0], FeatureKind::Exon);
    assert_eq!(kinds[1], FeatureKind::Utr);
    let starts: Vec<u32> = t.features(tx).map(|i| f.start[i as usize]).collect();
    assert!(starts.windows(2).all(|w| w[0] <= w[1]));
    assert_eq!(ann.transcript_lengths()[tx as usize], 15 + 54 + 2549);
    assert_eq!(ann.cds_lengths()[tx as usize], 9 + 969);
    // Minus-strand exons are position-sorted, not by exon_number.
    let wash = gene_by_id(&ann, "ENSG00000227232.5");
    let wt = g.transcripts(wash).start;
    let ws: Vec<u32> = t.features(wt).map(|i| f.start[i as usize]).collect();
    assert_eq!(ws, [14403, 24737, 29533]);
    assert_eq!(t.strand[wt as usize], Strand::Reverse);
    // Attributes keep order and quoting.
    let row = Row::Gene(or4f5);
    assert_eq!(ann.attribute(row, "level"), Some("2"));
    assert_eq!(ann.attributes(row).count(), 4);
}

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
        assert_eq!(a.features().start, b.features().start);
        assert_eq!(a.features().kind, b.features().kind);
        // And writing the re-read model again is a fixed point.
        assert_eq!(gtf_text(&b, GtfLayout::Original), model);
    }
}

#[test]
fn ensembl_and_refseq_gtf_fields() {
    let ens = read(ENSEMBL_GTF);
    assert!(ens.report().is_clean());
    let g = ens.gene("ENSG00000186092").unwrap() as usize;
    assert_eq!(ens.str(ens.genes().biotype[g]), "protein_coding");
    let t = ens.transcript("ENST00000641515").unwrap() as usize;
    assert_eq!(ens.str(ens.transcripts().biotype[t]), "protein_coding");
    assert!(ens.transcripts().tag_bits[t].contains(Tags::ENSEMBL_CANONICAL));
    let kinds: Vec<FeatureKind> = ens
        .transcripts()
        .features(t as u32)
        .map(|i| ens.features().kind[i as usize])
        .collect();
    assert!(
        kinds.contains(&FeatureKind::FivePrimeUtr) && kinds.contains(&FeatureKind::ThreePrimeUtr)
    );

    let rs = read(REFSEQ_GTF);
    assert!(rs.report().is_clean(), "{:?}", rs.report().problems);
    // `transcript_id ""` on gene lines does not make a transcript.
    assert_eq!(rs.transcripts().len(), 2);
    let g = rs.gene("DDX11L1").unwrap() as usize;
    assert_eq!(rs.str(rs.genes().name[g]), "DDX11L1");
    assert_eq!(rs.str(rs.genes().biotype[g]), "transcribed_pseudogene");
    let t = rs.transcript("NM_001005484.2").unwrap() as usize;
    assert!(rs.transcripts().tag_bits[t].contains(Tags::MANE_SELECT | Tags::REFSEQ_SELECT));
}

#[test]
fn ucsc_genes_and_transcripts_are_synthesized() {
    let ann = read(UCSC_GTF);
    let r = ann.report();
    assert!(r.is_clean(), "{:?}", r.problems);
    assert_eq!((r.synthesized_genes, r.synthesized_transcripts), (3, 3));
    let (g, t) = (ann.genes(), ann.transcripts());
    // One id on chrX and chrY: two genes, two transcripts.
    assert_eq!(ann.genes_by_id("PARG").len(), 2);
    assert_eq!(ann.transcripts_by_id("PART").len(), 2);
    let x = ann.genes_by_id("PARG")[0] as usize;
    assert_eq!(ann.seqid_name(g.seqid[x]), "chrX");
    assert_eq!(
        (g.start[x], g.end[x], g.strand[x]),
        (400, 1100, Strand::Reverse)
    );
    let d = ann.gene("ENST00000456328.2").unwrap() as usize;
    assert_eq!(ann.str(g.name[d]), "DDX11L2");
    assert_eq!((g.start[d], g.end[d]), (11868, 14409));
    let tx = g.transcripts(d as u32).start as usize;
    assert!(t.is_synthesized(tx as u32) && g.is_synthesized(d as u32));
    assert_eq!(ann.str(t.feature_type[tx]), "transcript");
    // Scores are kept as written.
    let f = t.features(tx as u32).start as usize;
    assert_eq!(ann.str(ann.features().score[f]), "0.000000");
}

#[test]
fn ensembl_gff3_hierarchy() {
    let ann = read(ENSEMBL_GFF3);
    assert_eq!(ann.dialect(), Dialect::Gff3);
    assert!(ann.report().is_clean(), "{:?}", ann.report().problems);
    let (g, t) = (ann.genes(), ann.transcripts());
    assert_eq!((g.len(), t.len()), (2, 2));
    // gene_id is preferred over ID.
    let lnc = ann.gene("ENSG00000290825").unwrap() as usize;
    assert_eq!(ann.str(g.feature_type[lnc]), "ncRNA_gene");
    assert_eq!(ann.str(g.name[lnc]), "DDX11L2");
    assert_eq!(ann.str(g.biotype[lnc]), "lncRNA");
    let tx = ann.transcript("ENST00000641515").unwrap();
    assert_eq!(ann.str(t.name[tx as usize]), "OR4F5-201");
    assert!(t.tag_bits[tx as usize].contains(Tags::MANE_SELECT));
    // Two CDS lines sharing an ID are two rows.
    let cds = t
        .features(tx)
        .filter(|&i| ann.features().kind[i as usize] == FeatureKind::Cds)
        .count();
    assert_eq!(cds, 2);
    // chromosome and biological_region are outside the hierarchy.
    assert_eq!(ann.others().len(), 2);
    let o = ann.others();
    assert_eq!(ann.str(o.score[1]), "1.3e+03");
    assert_eq!(ann.str(o.id[0]), "chromosome:1");
}

#[test]
fn refseq_gff3_hierarchy() {
    let ann = read(REFSEQ_GFF3);
    let r = ann.report();
    assert!(r.is_clean(), "{:?}", r.problems);
    let (g, t) = (ann.genes(), ann.transcripts());
    assert_eq!(g.len(), 4);
    let ddx = ann.gene("gene-DDX11L1").unwrap() as usize;
    assert_eq!(ann.str(g.feature_type[ddx]), "pseudogene");
    assert_eq!(ann.str(g.biotype[ddx]), "transcribed_pseudogene");
    assert_eq!(ann.str(g.name[ddx]), "DDX11L1");
    assert_eq!(
        ann.transcript("NR_046018.2").map(|x| t.gene[x as usize]),
        Some(ddx as u32)
    );
    // A V segment between gene and exons is the transcript.
    let v = ann.gene("gene-IGHV1-2").unwrap();
    let vt = g.transcripts(v).start as usize;
    assert_eq!(ann.str(t.feature_type[vt]), "V_gene_segment");
    assert_eq!(t.features(vt as u32).len(), 2);
    // A CDS directly under a gene gets a synthesized transcript.
    let nd1 = ann.gene("gene-ND1").unwrap();
    let nt = g.transcripts(nd1);
    assert_eq!(nt.len(), 1);
    assert!(t.is_synthesized(nt.start));
    assert_eq!(r.synthesized_transcripts, 1);
    // ##FASTA ends the records.
    assert_eq!(r.records, 19);
    assert_eq!(ann.others().len(), 1);
    let tx = ann.transcript("NM_001005484.2").unwrap();
    assert!(t.tag_bits[tx as usize].contains(Tags::MANE_SELECT));
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
fn problems_are_reported_or_raised() {
    let text = "\
c\ts\tgene\t1\t100\t.\t+\t.\tgene_id \"g1\";
c\ts\ttranscript\t1\t50\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";
c\ts\texon\t1\t60\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";
c\ts\texon\t70\t80\t.\t-\t.\tgene_id \"g1\"; transcript_id \"t1\";
c\ts\texon\t70\t80\t.\t+\t.\ttranscript_id \"t9\";
c\ts\tgene\t1\t100\t.\t+\t.\tgene_id \"g1\";
c\ts\texon\t1\t10\t.\t+\t.\tgene_id \"g2\"; transcript_id \"t1\";
";
    let ann = read(text);
    let r = ann.report();
    for kind in [
        ProblemKind::OutsideParent,
        ProblemKind::StrandConflict,
        ProblemKind::MissingGeneId,
        ProblemKind::DuplicateId,
        ProblemKind::GeneConflict,
    ] {
        assert!(r.count(kind) > 0, "{kind:?}: {:?}", r.problems);
    }
    // The orphan and the duplicate gene line are kept as others.
    assert_eq!(ann.others().len(), 2);
    assert_eq!(gtf_text(&ann, GtfLayout::Original), text);
    let strict = ReadOptions {
        strict: true,
        ..ReadOptions::default()
    };
    let err = Annotation::parse(text.as_bytes(), &strict).unwrap_err();
    assert!(
        matches!(err, omicsweave_core::Error::Format { line: 3, .. }),
        "{err}"
    );

    let gff = "##gff-version 3\nc\ts\tgene\t1\t100\t.\t+\t.\tID=g\nc\ts\texon\t1\t10\t.\t+\t.\tParent=nope\nc\ts\texon\t1\t10\t.\t+\t.\tID=e\nc\ts\tmRNA\t1\t100\t.\t+\t.\tID=g\n";
    let ann = read(gff);
    let r = ann.report();
    assert_eq!(r.count(ProblemKind::UnknownParent), 1);
    assert_eq!(r.count(ProblemKind::MissingParent), 1);
    assert_eq!(r.count(ProblemKind::DuplicateId), 1);

    for bad in [
        "c\ts\texon\t0\t10\t.\t+\t.\tgene_id \"g\";\n",
        "c\ts\texon\t10\t9\t.\t+\t.\tgene_id \"g\";\n",
        "c\ts\texon\t1\t9\t.\tx\t.\tgene_id \"g\";\n",
        "c\ts\texon\t1\t9\t.\t+\t3\tgene_id \"g\";\n",
        "c\ts\texon\t1\t9\tabc\t+\t.\tgene_id \"g\";\n",
        "c\ts\texon\t1\t9\n",
        "c\ts\texon\t1\t9\t.\t+\t.\tgene_id \"g;\n",
    ] {
        assert!(
            Annotation::parse(bad.as_bytes(), &ReadOptions::default()).is_err(),
            "{bad}"
        );
    }
}

#[test]
fn aliases_rename_and_match() {
    let opts = ReadOptions {
        aliases: Some(Aliases::from_pairs([("1", "chr1"), ("MT", "chrM")])),
        ..ReadOptions::default()
    };
    let ann = Annotation::parse(ENSEMBL_GTF.as_bytes(), &opts).unwrap();
    let names: Vec<&str> = ann.seqid_names().collect();
    assert_eq!(names, ["chr1", "chrM"]);
    let mut hits = Vec::new();
    ann.overlaps(Level::Gene, Region::new("chr1", 65000, 66000), &mut hits);
    assert_eq!(hits.len(), 1);

    let plain = read(ENSEMBL_GTF);
    let bam = ["chr1", "chr2", "chrM"];
    assert_eq!(plain.match_seqids(&bam, None), [Some(0), None, Some(1)]);
    let lengths = plain.seqid_lengths(&[("chr1", 248956422), ("chrM", 16569)], None);
    assert_eq!(lengths, [Some(248956422), Some(16569)]);
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
