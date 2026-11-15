//! The columnar tables of an annotation and their small value types.

use std::ops::Range;

use crate::strings::NONE;

/// Strand of a record, normalized from the GTF / GFF3 strand column.
pub use omicsweave_core::Strand;

/// What a sub-transcript feature is, from its type column.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum FeatureKind {
    Exon = 0,
    Cds = 1,
    FivePrimeUtr = 2,
    ThreePrimeUtr = 3,
    /// `UTR` without a side (GENCODE GTF).
    Utr = 4,
    StartCodon = 5,
    StopCodon = 6,
    /// `Selenocysteine` / `stop_codon_redefined_as_selenocysteine`.
    Selenocysteine = 7,
    /// Any other type found below a transcript (`intron`, `polyA_site`, ...).
    Other = 8,
}

impl FeatureKind {
    /// The kind of a feature type name; `None` for names that are not
    /// sub-transcript features in GFF3 (`gene`, `mRNA`, `region`, ...).
    pub fn of_type(t: &[u8]) -> Option<FeatureKind> {
        Some(match t {
            b"exon" | b"noncoding_exon" | b"pseudogenic_exon" => FeatureKind::Exon,
            b"CDS" => FeatureKind::Cds,
            b"five_prime_UTR" | b"5UTR" | b"five_prime_utr" | b"5'UTR" => FeatureKind::FivePrimeUtr,
            b"three_prime_UTR" | b"3UTR" | b"three_prime_utr" | b"3'UTR" => {
                FeatureKind::ThreePrimeUtr
            }
            b"UTR" => FeatureKind::Utr,
            b"start_codon" => FeatureKind::StartCodon,
            b"stop_codon" => FeatureKind::StopCodon,
            b"Selenocysteine" | b"stop_codon_redefined_as_selenocysteine" => {
                FeatureKind::Selenocysteine
            }
            b"intron"
            | b"polyA_site"
            | b"polyA_signal_sequence"
            | b"transcription_start_site"
            | b"TSS"
            | b"splice_acceptor_site"
            | b"splice_donor_site"
            | b"intron_retention" => FeatureKind::Other,
            _ => return None,
        })
    }

    /// Coding parts (CDS and codons), which make up the thick part of BED12.
    pub fn is_coding(self) -> bool {
        matches!(
            self,
            FeatureKind::Cds | FeatureKind::StartCodon | FeatureKind::StopCodon
        )
    }

    /// Parts that lie on the mature transcript: exon pieces and the coding
    /// and UTR pieces exons consist of.
    pub fn is_exonic(self) -> bool {
        !matches!(self, FeatureKind::Other)
    }
}

/// Well-known transcript tags as bits; every tag is also kept by name
/// ([`Transcripts::tags`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Tags(pub u32);

impl Tags {
    pub const BASIC: Tags = Tags(1);
    pub const ENSEMBL_CANONICAL: Tags = Tags(1 << 1);
    pub const MANE_SELECT: Tags = Tags(1 << 2);
    pub const MANE_PLUS_CLINICAL: Tags = Tags(1 << 3);
    pub const CCDS: Tags = Tags(1 << 4);
    pub const GENCODE_PRIMARY: Tags = Tags(1 << 5);
    pub const APPRIS_PRINCIPAL: Tags = Tags(1 << 6);
    pub const APPRIS_ALTERNATIVE: Tags = Tags(1 << 7);
    pub const REFSEQ_SELECT: Tags = Tags(1 << 8);
    pub const MRNA_START_NF: Tags = Tags(1 << 9);
    pub const MRNA_END_NF: Tags = Tags(1 << 10);
    pub const CDS_START_NF: Tags = Tags(1 << 11);
    pub const CDS_END_NF: Tags = Tags(1 << 12);
    pub const READTHROUGH: Tags = Tags(1 << 13);

    /// (bit, canonical name) of every known tag.
    pub const KNOWN: [(Tags, &'static str); 14] = [
        (Tags::BASIC, "basic"),
        (Tags::ENSEMBL_CANONICAL, "Ensembl_canonical"),
        (Tags::MANE_SELECT, "MANE_Select"),
        (Tags::MANE_PLUS_CLINICAL, "MANE_Plus_Clinical"),
        (Tags::CCDS, "CCDS"),
        (Tags::GENCODE_PRIMARY, "GENCODE_Primary"),
        (Tags::APPRIS_PRINCIPAL, "appris_principal"),
        (Tags::APPRIS_ALTERNATIVE, "appris_alternative"),
        (Tags::REFSEQ_SELECT, "RefSeq_Select"),
        (Tags::MRNA_START_NF, "mRNA_start_NF"),
        (Tags::MRNA_END_NF, "mRNA_end_NF"),
        (Tags::CDS_START_NF, "cds_start_NF"),
        (Tags::CDS_END_NF, "cds_end_NF"),
        (Tags::READTHROUGH, "readthrough_transcript"),
    ];

    /// The bit of a tag value. Case and `' '` versus `'_'` do not matter
    /// (RefSeq writes `MANE Select`); `appris_principal_1` ... `_5` and
    /// `appris_alternative_*` map to their family.
    pub fn of_name(name: &[u8]) -> Tags {
        let mut buf = [0u8; 32];
        if name.len() > buf.len() {
            return Tags(0);
        }
        for (d, &c) in buf.iter_mut().zip(name) {
            *d = if c == b' ' {
                b'_'
            } else {
                c.to_ascii_lowercase()
            };
        }
        let n = &buf[..name.len()];
        if n.starts_with(b"appris_principal") {
            return Tags::APPRIS_PRINCIPAL;
        }
        if n.starts_with(b"appris_alternative") {
            return Tags::APPRIS_ALTERNATIVE;
        }
        for (bit, known) in Tags::KNOWN {
            if known.len() == n.len() && known.as_bytes().eq_ignore_ascii_case(n) {
                return bit;
            }
        }
        Tags(0)
    }

    /// The bit of a canonical or alias tag name.
    pub fn parse(name: &str) -> Option<Tags> {
        let t = Tags::of_name(name.as_bytes());
        (t.0 != 0).then_some(t)
    }

    pub fn contains(self, other: Tags) -> bool {
        self.0 & other.0 == other.0
    }

    pub fn intersects(self, other: Tags) -> bool {
        self.0 & other.0 != 0
    }
}

impl std::ops::BitOr for Tags {
    type Output = Tags;
    fn bitor(self, o: Tags) -> Tags {
        Tags(self.0 | o.0)
    }
}

impl std::ops::BitOrAssign for Tags {
    fn bitor_assign(&mut self, o: Tags) {
        self.0 |= o.0;
    }
}

/// The text dialect of an annotation file.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Dialect {
    /// GTF 2.2 (`key "value";` attributes, hierarchy by `gene_id` /
    /// `transcript_id`).
    Gtf,
    /// GFF3 (`key=value` attributes, hierarchy by `ID` / `Parent`).
    Gff3,
}

/// Which attributes (column 9) are kept. The hierarchy, names, biotypes
/// and tags are read whatever the policy.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum AttributePolicy {
    /// Every attribute in file order: lossless, and needed to write the
    /// original GTF back.
    #[default]
    All,
    /// Only these keys, in file order.
    Only(Vec<String>),
    /// No attribute table.
    None,
}

/// Gene rows, sorted by (seqid, start, end). Gene `g`'s transcripts are
/// rows `transcripts(g)` of [`Transcripts`].
///
/// String columns hold ids into the annotation's
/// [`StringPool`](crate::StringPool), [`NONE`] when absent; `record` is the
/// row's line in file order ([`NONE`] for a record synthesized from its
/// children).
#[derive(Debug, Clone, Default)]
pub struct Genes {
    pub seqid: Vec<u32>,
    pub start: Vec<u32>,
    pub end: Vec<u32>,
    pub strand: Vec<Strand>,
    /// `gene_id` (GTF), or `gene_id` / `ID` (GFF3).
    pub id: Vec<u32>,
    /// `gene_name`, `gene`, or GFF3 `Name`.
    pub name: Vec<u32>,
    /// `gene_type`, `gene_biotype`, or GFF3 `biotype`.
    pub biotype: Vec<u32>,
    pub source: Vec<u32>,
    /// The type column (`gene`, `ncRNA_gene`, ...).
    pub feature_type: Vec<u32>,
    /// The score column as written, [`NONE`] for `.`.
    pub score: Vec<u32>,
    /// The phase column (`0`..`2`, or [`PHASE_NONE`]).
    pub phase: Vec<u8>,
    pub record: Vec<u32>,
    /// `transcript_offsets[g]..transcript_offsets[g + 1]`: rows in
    /// [`Transcripts`]; one more element than genes.
    pub transcript_offsets: Vec<u32>,
}

impl Genes {
    pub fn len(&self) -> usize {
        self.start.len()
    }

    pub fn is_empty(&self) -> bool {
        self.start.is_empty()
    }

    /// The transcripts of gene `g`.
    #[inline]
    pub fn transcripts(&self, g: u32) -> Range<u32> {
        self.transcript_offsets[g as usize]..self.transcript_offsets[g as usize + 1]
    }

    /// Whether gene `g` had no record of its own.
    pub fn is_synthesized(&self, g: u32) -> bool {
        self.record[g as usize] == NONE
    }
}

/// Transcript rows: those of each gene are contiguous and sorted by
/// (start, end). Transcript `t`'s features are rows `features(t)` of
/// [`Features`].
#[derive(Debug, Clone, Default)]
pub struct Transcripts {
    pub seqid: Vec<u32>,
    pub start: Vec<u32>,
    pub end: Vec<u32>,
    pub strand: Vec<Strand>,
    pub gene: Vec<u32>,
    /// `transcript_id` (GTF), or `transcript_id` / `ID` (GFF3).
    pub id: Vec<u32>,
    /// `transcript_name` or GFF3 `Name`.
    pub name: Vec<u32>,
    /// `transcript_type`, `transcript_biotype`, or GFF3 `biotype`.
    pub biotype: Vec<u32>,
    pub source: Vec<u32>,
    pub feature_type: Vec<u32>,
    pub score: Vec<u32>,
    pub phase: Vec<u8>,
    pub record: Vec<u32>,
    /// Known tags as bits.
    pub tag_bits: Vec<Tags>,
    /// Every tag value (string ids), `tag_offsets[t]..tag_offsets[t + 1]`.
    pub tag_values: Vec<u32>,
    pub tag_offsets: Vec<u32>,
    /// `feature_offsets[t]..feature_offsets[t + 1]`: rows in [`Features`].
    pub feature_offsets: Vec<u32>,
}

impl Transcripts {
    pub fn len(&self) -> usize {
        self.start.len()
    }

    pub fn is_empty(&self) -> bool {
        self.start.is_empty()
    }

    #[inline]
    pub fn features(&self, t: u32) -> Range<u32> {
        self.feature_offsets[t as usize]..self.feature_offsets[t as usize + 1]
    }

    /// Tag values of transcript `t` as string ids.
    pub fn tags(&self, t: u32) -> &[u32] {
        &self.tag_values
            [self.tag_offsets[t as usize] as usize..self.tag_offsets[t as usize + 1] as usize]
    }

    pub fn is_synthesized(&self, t: u32) -> bool {
        self.record[t as usize] == NONE
    }
}

/// Sub-transcript features (exons, CDS, UTRs, codons, ...): those of each
/// transcript are contiguous and sorted by (start, end, kind). A GFF3
/// record with several transcript parents has a row under each.
#[derive(Debug, Clone, Default)]
pub struct Features {
    pub transcript: Vec<u32>,
    pub kind: Vec<FeatureKind>,
    pub start: Vec<u32>,
    pub end: Vec<u32>,
    pub strand: Vec<Strand>,
    /// `0`, `1`, `2`, or [`PHASE_NONE`] for `.`.
    pub phase: Vec<u8>,
    pub feature_type: Vec<u32>,
    pub source: Vec<u32>,
    pub score: Vec<u32>,
    pub record: Vec<u32>,
}

/// The phase of a feature written as `.`.
pub const PHASE_NONE: u8 = u8::MAX;

impl Features {
    pub fn len(&self) -> usize {
        self.start.len()
    }

    pub fn is_empty(&self) -> bool {
        self.start.is_empty()
    }
}

/// Records outside the gene → transcript → feature hierarchy (GFF3
/// `region`, `chromosome`, records with problems), in file order.
#[derive(Debug, Clone, Default)]
pub struct Others {
    pub seqid: Vec<u32>,
    pub start: Vec<u32>,
    pub end: Vec<u32>,
    pub strand: Vec<Strand>,
    pub phase: Vec<u8>,
    pub feature_type: Vec<u32>,
    pub source: Vec<u32>,
    pub score: Vec<u32>,
    /// `ID` (GFF3), or `gene_id` (GTF).
    pub id: Vec<u32>,
    pub record: Vec<u32>,
}

impl Others {
    pub fn len(&self) -> usize {
        self.start.len()
    }

    pub fn is_empty(&self) -> bool {
        self.start.is_empty()
    }
}

/// Which table a file record went to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Row {
    Gene(u32),
    Transcript(u32),
    Feature(u32),
    Other(u32),
}

/// A [`Row`] packed in 32 bits: table in the top two bits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct RowRef(pub u32);

pub(crate) const ROW_LIMIT: u32 = 1 << 30;

impl RowRef {
    pub fn new(row: Row) -> RowRef {
        RowRef(match row {
            Row::Gene(i) => i,
            Row::Transcript(i) => (1 << 30) | i,
            Row::Feature(i) => (2 << 30) | i,
            Row::Other(i) => (3 << 30) | i,
        })
    }

    pub fn get(self) -> Row {
        let i = self.0 & (ROW_LIMIT - 1);
        match self.0 >> 30 {
            0 => Row::Gene(i),
            1 => Row::Transcript(i),
            2 => Row::Feature(i),
            _ => Row::Other(i),
        }
    }
}

/// Column 9: each distinct (key, value) pair once, and per file record
/// the pairs in file order: record `r`'s pairs are
/// `ids[offsets[r]..offsets[r + 1]]`, indexes into `pairs`.
#[derive(Debug, Clone, Default)]
pub struct Attributes {
    /// (key, value) string ids; the key's top bit marks a GTF value written
    /// without quotes.
    pub(crate) pairs: Vec<(u32, u32)>,
    pub(crate) ids: Vec<u32>,
    pub(crate) offsets: Vec<u32>,
}

pub(crate) const UNQUOTED: u32 = 1 << 31;

impl Attributes {
    /// The (key, value) string ids of file record `record`.
    pub fn of(&self, record: u32) -> impl ExactSizeIterator<Item = (u32, u32)> + '_ {
        self.ids[self.raw(record)].iter().map(|&p| {
            let (k, v) = self.pairs[p as usize];
            (k & !UNQUOTED, v)
        })
    }

    /// Positions in `ids` of record `record`'s pairs.
    pub(crate) fn raw(&self, record: u32) -> Range<usize> {
        if record == NONE || record as usize + 1 >= self.offsets.len() {
            return 0..0;
        }
        self.offsets[record as usize] as usize..self.offsets[record as usize + 1] as usize
    }

    /// The pair at position `i` of `ids`: (key with quoting bit, value).
    #[inline]
    pub(crate) fn at(&self, i: usize) -> (u32, u32) {
        self.pairs[self.ids[i] as usize]
    }

    /// Attribute occurrences held.
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }

    /// Distinct (key, value) pairs.
    pub fn distinct(&self) -> usize {
        self.pairs.len()
    }

    pub fn heap_bytes(&self) -> usize {
        8 * self.pairs.capacity() + 4 * (self.ids.capacity() + self.offsets.capacity())
    }
}
