//! Dialect interpretation, hierarchy, identifiers and sequence names.
use crate::support::{fixtures::*, read};
use omicsweave_core::Region;
use omicsweave_gtf::{
    Aliases, Annotation, Dialect, FeatureKind, Level, ReadOptions, Row, Strand, Tags,
};

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
fn closed_single_base_coordinates_and_duplicate_ids_stay_seqid_local() {
    let text = "chr1\ts\texon\t1\t1\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\"; gene_name \"shared\";\n\
                chr2\ts\texon\t4294967295\t4294967295\t.\t-\t.\tgene_id \"g\"; transcript_id \"t\"; gene_name \"shared\";\n";
    let ann = read(text);
    assert!(ann.report().is_clean(), "{:?}", ann.report());
    assert_eq!(ann.genes_by_id("g"), [0, 1]);
    assert_eq!(ann.genes_by_name("shared"), [0, 1]);
    assert_eq!(ann.transcripts_by_id("t"), [0, 1]);
    assert_eq!(ann.gene("g"), Some(0));
    assert_eq!(ann.transcript("t"), Some(0));
    assert!(ann.genes_by_id("absent").is_empty());
    assert_eq!(ann.transcript("absent"), None);
    assert_eq!(ann.features().start, [0, u32::MAX - 1]);
    assert_eq!(ann.features().end, [1, u32::MAX]);
    assert_eq!(ann.transcript_lengths(), [1, 1]);
    for (seqid, start, end, expected) in [
        ("chr1", 0, 1, vec![0]),
        ("chr1", 1, 2, vec![]),
        ("chr2", u32::MAX - 1, u32::MAX, vec![1]),
    ] {
        let mut hits = Vec::new();
        ann.overlaps(Level::Exon, Region::new(seqid, start, end), &mut hits);
        assert_eq!(hits, expected, "{seqid}:{start}-{end}");
    }
    assert_eq!(
        crate::support::gtf_text(&ann, omicsweave_gtf::GtfLayout::Original),
        text
    );
}

#[test]
fn explicit_aliases_match_both_directions_and_exact_names_take_priority() {
    let aliases = Aliases::parse("# aliases\n\nNC_1\tchr1\nMT chrM\n").unwrap();
    assert_eq!(aliases.rename("NC_1"), "chr1");
    assert_eq!(aliases.rename("unknown"), "unknown");
    assert_eq!(
        omicsweave_gtf::match_names(&["NC_1", "MT"], &["chr1", "chrM", "absent"], Some(&aliases)),
        [Some(0), Some(1), None]
    );
    assert_eq!(
        omicsweave_gtf::match_names(&["chr1", "chrM"], &["NC_1", "MT"], Some(&aliases)),
        [Some(0), Some(1)]
    );
    assert_eq!(
        omicsweave_gtf::match_names(
            &["NC_1", "chr1", "1"],
            &["chr1", "NC_1", "1"],
            Some(&aliases)
        ),
        [Some(1), Some(0), Some(2)]
    );
    let err = Aliases::parse("# aliases\n\nmissing_target\n").unwrap_err();
    assert!(
        matches!(err, omicsweave_core::Error::Format { line: 3, .. }),
        "{err}"
    );
}
