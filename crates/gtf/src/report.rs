//! What reading found wrong with an annotation, without failing it.

use crate::strings::NONE;

/// A kind of problem in an annotation file.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum ProblemKind {
    /// A GTF record without `gene_id`; kept among [`Others`](crate::Others).
    MissingGeneId,
    /// A GTF `transcript` record without `transcript_id`.
    MissingTranscriptId,
    /// A GFF3 sub-transcript feature (exon, CDS, ...) without `Parent`.
    MissingParent,
    /// A GFF3 `Parent` naming no record's `ID`.
    UnknownParent,
    /// A second gene or transcript record with the same identity (GTF:
    /// same seqid and id; GFF3: same `ID`); kept among the others.
    DuplicateId,
    /// Records of one transcript naming different genes.
    GeneConflict,
    /// A feature on another strand than its transcript, or a transcript on
    /// another strand than its gene.
    StrandConflict,
    /// A GFF3 child on another seqid than its parent; kept among the others.
    SeqidConflict,
    /// A feature outside its transcript's span, or a transcript outside
    /// its gene's.
    OutsideParent,
    /// A GFF3 hierarchy that is not gene → transcript → feature (an exon
    /// directly under a gene that also has transcripts, a transcript under
    /// a transcript, or cyclic Parent links). Records that cannot be
    /// assigned to the hierarchy are kept among the others.
    BadNesting,
    /// A GFF3 transcript with several parents; the first is used.
    MultipleParents,
}

impl ProblemKind {
    pub const ALL: [ProblemKind; 11] = [
        ProblemKind::MissingGeneId,
        ProblemKind::MissingTranscriptId,
        ProblemKind::MissingParent,
        ProblemKind::UnknownParent,
        ProblemKind::DuplicateId,
        ProblemKind::GeneConflict,
        ProblemKind::StrandConflict,
        ProblemKind::SeqidConflict,
        ProblemKind::OutsideParent,
        ProblemKind::BadNesting,
        ProblemKind::MultipleParents,
    ];

    /// A snake_case name.
    pub fn name(self) -> &'static str {
        match self {
            ProblemKind::MissingGeneId => "missing_gene_id",
            ProblemKind::MissingTranscriptId => "missing_transcript_id",
            ProblemKind::MissingParent => "missing_parent",
            ProblemKind::UnknownParent => "unknown_parent",
            ProblemKind::DuplicateId => "duplicate_id",
            ProblemKind::GeneConflict => "gene_conflict",
            ProblemKind::StrandConflict => "strand_conflict",
            ProblemKind::SeqidConflict => "seqid_conflict",
            ProblemKind::OutsideParent => "outside_parent",
            ProblemKind::BadNesting => "bad_nesting",
            ProblemKind::MultipleParents => "multiple_parents",
        }
    }

    pub fn describe(self) -> &'static str {
        match self {
            ProblemKind::MissingGeneId => "record without gene_id",
            ProblemKind::MissingTranscriptId => "transcript record without transcript_id",
            ProblemKind::MissingParent => "transcript feature without Parent",
            ProblemKind::UnknownParent => "Parent names no ID",
            ProblemKind::DuplicateId => "duplicate gene or transcript identity",
            ProblemKind::GeneConflict => "transcript records name different genes",
            ProblemKind::StrandConflict => "strand differs from parent",
            ProblemKind::SeqidConflict => "seqid differs from parent",
            ProblemKind::OutsideParent => "extends outside parent",
            ProblemKind::BadNesting => "not a gene > transcript > feature hierarchy",
            ProblemKind::MultipleParents => "transcript with several parents",
        }
    }
}

/// One problem: its kind, the 1-based line of the record, and the id
/// involved (a string id; [`NONE`] when there is none).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Problem {
    pub kind: ProblemKind,
    pub line: u64,
    pub id: u32,
}

/// Counts of what was read, and every problem found, in line order.
#[derive(Debug, Clone, Default)]
pub struct Report {
    /// Feature lines read (comments and blank lines excluded).
    pub records: usize,
    pub comments: usize,
    pub genes: usize,
    pub transcripts: usize,
    /// Feature rows (a GFF3 record with several parents counts once per
    /// parent).
    pub features: usize,
    pub others: usize,
    /// Genes and transcripts built from their children, without a record.
    pub synthesized_genes: usize,
    pub synthesized_transcripts: usize,
    pub problems: Vec<Problem>,
}

impl Report {
    pub fn count(&self, kind: ProblemKind) -> usize {
        self.problems.iter().filter(|p| p.kind == kind).count()
    }

    pub fn is_clean(&self) -> bool {
        self.problems.is_empty()
    }

    pub(crate) fn push(&mut self, kind: ProblemKind, line: u64, id: u32) {
        self.problems.push(Problem { kind, line, id });
    }

    pub(crate) fn finish(&mut self) {
        self.problems.sort_by_key(|p| (p.line, p.kind));
        self.problems.dedup();
    }
}

impl Problem {
    pub fn has_id(&self) -> bool {
        self.id != NONE
    }
}
