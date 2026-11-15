//! Subsets of an annotation as row lists, without copying its data.

use crate::annotation::Annotation;
use crate::model::Tags;
use crate::strings::NONE;

/// Selected gene and transcript rows, ascending. Every selected
/// transcript's gene is selected.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Selection {
    pub genes: Vec<u32>,
    pub transcripts: Vec<u32>,
}

/// Filters for [`Annotation::select`]; filters combine with "and", values
/// within one filter with "or".
#[derive(Debug, Clone, Default)]
pub struct Select {
    gene_biotypes: Option<Vec<String>>,
    transcript_biotypes: Option<Vec<String>>,
    all_tags: Tags,
    any_tags: Option<Vec<String>>,
    genes: Option<Vec<String>>,
    seqids: Option<Vec<String>>,
}

impl Select {
    pub fn new() -> Select {
        Select::default()
    }

    /// Genes of these biotypes.
    pub fn gene_biotypes<S: Into<String>>(mut self, v: impl IntoIterator<Item = S>) -> Self {
        self.gene_biotypes = Some(v.into_iter().map(Into::into).collect());
        self
    }

    /// Transcripts of these biotypes (and the genes having any).
    pub fn transcript_biotypes<S: Into<String>>(mut self, v: impl IntoIterator<Item = S>) -> Self {
        self.transcript_biotypes = Some(v.into_iter().map(Into::into).collect());
        self
    }

    /// Transcripts carrying all of these known tags.
    pub fn tags(mut self, tags: Tags) -> Self {
        self.all_tags |= tags;
        self
    }

    /// Transcripts carrying any of these tag values (exact text).
    pub fn tag_values<S: Into<String>>(mut self, v: impl IntoIterator<Item = S>) -> Self {
        self.any_tags = Some(v.into_iter().map(Into::into).collect());
        self
    }

    /// Genes with these ids or names.
    pub fn genes<S: Into<String>>(mut self, v: impl IntoIterator<Item = S>) -> Self {
        self.genes = Some(v.into_iter().map(Into::into).collect());
        self
    }

    /// Genes on these chromosomes.
    pub fn seqids<S: Into<String>>(mut self, v: impl IntoIterator<Item = S>) -> Self {
        self.seqids = Some(v.into_iter().map(Into::into).collect());
        self
    }
}

impl Annotation {
    /// The genes and transcripts passing `filter`. A gene is selected when
    /// it passes the gene filters and, if there are transcript filters,
    /// has a transcript passing them.
    pub fn select(&self, filter: &Select) -> Selection {
        let ids = |v: &Option<Vec<String>>| -> Option<Vec<u32>> {
            v.as_ref().map(|v| {
                let mut ids: Vec<u32> = v
                    .iter()
                    .filter_map(|s| self.pool.get(s.as_bytes()))
                    .collect();
                ids.sort_unstable();
                ids
            })
        };
        let has = |set: &Option<Vec<u32>>, id: u32| {
            set.as_ref()
                .is_none_or(|s| id != NONE && s.binary_search(&id).is_ok())
        };
        let gene_types = ids(&filter.gene_biotypes);
        let tx_types = ids(&filter.transcript_biotypes);
        let any_tags = ids(&filter.any_tags);
        let seqids: Option<Vec<u32>> = filter.seqids.as_ref().map(|v| {
            let mut s: Vec<u32> = v.iter().filter_map(|n| self.seqid(n)).collect();
            s.sort_unstable();
            s
        });
        let gene_rows: Option<Vec<u32>> = filter.genes.as_ref().map(|v| {
            let mut rows: Vec<u32> = v
                .iter()
                .flat_map(|s| {
                    self.genes_by_id(s)
                        .iter()
                        .chain(self.genes_by_name(s))
                        .copied()
                })
                .collect();
            rows.sort_unstable();
            rows.dedup();
            rows
        });
        let tx_filtered =
            tx_types.is_some() || any_tags.is_some() || filter.all_tags != Tags::default();
        let g = &self.genes;
        let t = &self.transcripts;
        let mut out = Selection::default();
        let candidates = gene_rows.as_ref().map_or(g.len(), Vec::len);
        for candidate in 0..candidates {
            let gene = gene_rows
                .as_ref()
                .map_or(candidate as u32, |rows| rows[candidate]);
            let i = gene as usize;
            if !has(&gene_types, g.biotype[i])
                || seqids
                    .as_ref()
                    .is_some_and(|s| s.binary_search(&g.seqid[i]).is_err())
            {
                continue;
            }
            let before = out.transcripts.len();
            for tx in g.transcripts(gene) {
                let j = tx as usize;
                if !has(&tx_types, t.biotype[j]) || !t.tag_bits[j].contains(filter.all_tags) {
                    continue;
                }
                if let Some(any) = &any_tags
                    && !t.tags(tx).iter().any(|v| any.binary_search(v).is_ok())
                {
                    continue;
                }
                out.transcripts.push(tx);
            }
            if !tx_filtered || out.transcripts.len() > before {
                out.genes.push(gene);
            }
        }
        out
    }
}
