//! Combined filters and the scope of derived views.
use crate::support::{base_sets::rows, fixtures::GENCODE_GTF, random, read};
use omicsweave_gtf::{Flank, Level, Select, Selection, Tags};

#[test]
fn selections_restrict_views() {
    let a = random(5);
    let sel = a.select(&Select::new().transcript_biotypes(["tt0"]));
    let tx = a.transcripts();
    assert!(!sel.transcripts.is_empty());
    for &t in &sel.transcripts {
        assert_eq!(a.str(tx.biotype[t as usize]), "tt0");
        assert!(sel.genes.binary_search(&tx.gene[t as usize]).is_ok());
    }
    let want = (0..tx.len() as u32)
        .filter(|&t| a.str(tx.biotype[t as usize]) == "tt0")
        .count();
    assert_eq!(sel.transcripts.len(), want);
    let introns = a.derive().selection(&sel).introns();
    assert!(
        introns
            .owner
            .iter()
            .all(|t| sel.transcripts.binary_search(t).is_ok())
    );
    let all = a
        .introns()
        .filter_owners(|t| sel.transcripts.binary_search(&t).is_ok());
    assert_eq!(rows(&introns), rows(&all));

    let by_gene = a.select(&Select::new().genes(["G0_1", "G2_3"]).seqids(["c0"]));
    assert_eq!(by_gene.genes, [a.gene("G0_1").unwrap()]);
    let bt = a.select(&Select::new().gene_biotypes(["bt1"]));
    assert!(
        bt.genes
            .iter()
            .all(|&g| a.str(a.genes().biotype[g as usize]) == "bt1")
    );

    let gencode = read(GENCODE_GTF);
    let mane = gencode.select(&Select::new().tags(Tags::MANE_SELECT));
    assert_eq!(
        mane.transcripts,
        [gencode.transcript("ENST00000641515.2").unwrap()]
    );
    let basic = gencode.select(&Select::new().tag_values(["basic"]));
    assert_eq!(basic.transcripts.len(), 4);
}

#[test]
fn sparse_gene_candidates_match_ids_names_and_combined_filters() {
    let ann = read(GENCODE_GTF);
    let wanted = ["OR4F5", "ENSG00000290825.1", "OR4F5", "absent"];
    for seqids in [None, Some(vec!["chr1"]), Some(vec!["chrM"])] {
        let mut filter = Select::new().genes(wanted);
        if let Some(seqids) = &seqids {
            filter = filter.seqids(seqids.iter().copied());
        }
        let got = ann.select(&filter);
        let g = ann.genes();
        let genes: Vec<u32> = (0..g.len() as u32)
            .filter(|&r| {
                let i = r as usize;
                (wanted.contains(&ann.str(g.id[i])) || wanted.contains(&ann.str(g.name[i])))
                    && seqids
                        .as_ref()
                        .is_none_or(|s| s.contains(&ann.seqid_name(g.seqid[i])))
            })
            .collect();
        let transcripts: Vec<u32> = genes
            .iter()
            .flat_map(|&g| ann.genes().transcripts(g))
            .collect();
        assert_eq!(got.genes, genes);
        assert_eq!(got.transcripts, transcripts);
    }
    let empty = ann.select(&Select::new().genes(Vec::<String>::new()));
    assert!(empty.genes.is_empty());
    assert!(empty.transcripts.is_empty());
}

#[test]
fn empty_filter_values_and_empty_selection_produce_empty_views() {
    let ann = read(GENCODE_GTF);
    let all = ann.select(&Select::new());
    assert_eq!(all.genes.len(), ann.genes().len());
    assert_eq!(all.transcripts.len(), ann.transcripts().len());
    for filter in [
        Select::new().genes(Vec::<String>::new()),
        Select::new().seqids(Vec::<String>::new()),
        Select::new().gene_biotypes(Vec::<String>::new()),
        Select::new().transcript_biotypes(Vec::<String>::new()),
        Select::new().tag_values(Vec::<String>::new()),
        Select::new().genes(["absent"]),
    ] {
        let empty = ann.select(&filter);
        assert_eq!(empty, Selection::default(), "filter={filter:?}");
        empty.check(&ann).unwrap();
        let scope = ann.derive().selection(&empty);
        for set in [
            scope.introns(),
            scope.exonic_regions(),
            scope.tss(Level::Gene),
            scope.tes(Level::Transcript),
            scope.promoters(Level::Gene, Flank::new(10, 5)),
            scope.tes_windows(Level::Transcript, Flank::new(5, 10)),
            scope.gene_bodies(Level::Gene, Flank::new(10, 0)),
        ] {
            assert!(set.is_empty(), "filter={filter:?}, set={set:?}");
            assert_eq!(set.seqid_count(), ann.seqid_count());
        }
        assert!(scope.junctions().is_empty());
        for aggregate in [false, true] {
            let bins = scope.exon_bins(aggregate);
            assert!(bins.set.is_empty());
            assert!(bins.transcripts.is_empty());
            assert!(bins.genes.is_empty());
        }
    }
    assert_eq!(ann.select(&Select::new()), all);
    assert!(!ann.introns().is_empty());
}
