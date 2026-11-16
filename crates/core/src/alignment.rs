//! The alignment family's output layout: fixed fields of each record as a
//! `#[repr(C)]` row, and the variable-length fields a query asks for as
//! Arrow `LargeString`-style columns, all in caller-allocated memory.

use std::ops::BitOr;

use crate::error::{Error, Result, check_output};
use crate::layout::StringColumn;
use crate::region::Region;

/// The fixed fields of one alignment record.
///
/// numpy: `np.dtype([("start", "<i8"), ("end", "<i8"), ("mate_start",
/// "<i4"), ("mate_tid", "<i4"), ("tlen", "<i4"), ("flag", "<u2"), ("mapq",
/// "u1"), ("_pad", "u1")])`.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct AlignmentEntry {
    /// 0-based leftmost position (`POS - 1`); -1 when the record has none.
    pub start: i64,
    /// End of the reference span, as htslib's `bam_endpos`: `start` plus
    /// the CIGAR's reference length, or plus 1 when the record is unmapped
    /// or its CIGAR consumes no reference.
    pub end: i64,
    /// 0-based position of the mate (`PNEXT - 1`); -1 when none.
    pub mate_start: i32,
    /// Reference index of the mate; -1 when none.
    pub mate_tid: i32,
    /// Observed template length (`TLEN`).
    pub tlen: i32,
    pub flag: u16,
    pub mapq: u8,
    pub _pad: u8,
}

const _: () = assert!(
    std::mem::size_of::<AlignmentEntry>() == 32 && std::mem::align_of::<AlignmentEntry>() == 8
);

/// One base modification of a record, as its `MM` / `ML` tags give it
/// (SAMtags §1.7), at a base of the record's SEQ.
///
/// numpy: `np.dtype([("ref_pos", "<i8"), ("qpos", "<i4"), ("code", "<i4"),
/// ("canonical", "u1"), ("probability", "u1"), ("flags", "u1"), ("_pad",
/// "u1", 5)])`.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct ModificationEntry {
    /// 0-based reference position of the base; -1 when it is not aligned
    /// (inserted, soft-clipped, or the record is unmapped).
    pub ref_pos: i64,
    /// 0-based index of the base in SEQ as stored, including soft clips.
    /// FLAG REVERSE determines its relation to the sequenced orientation.
    pub qpos: i32,
    /// The modification: a single-letter code as its ASCII value (`'m'`
    /// is 109), or a ChEBI identifier negated (`-17596` for `N+17596`),
    /// as htslib's `hts_base_mod.modified_base` holds it.
    pub code: i32,
    /// The unmodified base of the `MM` entry as ASCII (`A`, `C`, `G`, `T`
    /// or `N`; `U` reads as `T`), on the strand the read was sequenced.
    pub canonical: u8,
    /// The `ML` value, the probability scaled to 0–255; 0 when the record
    /// has no `ML` tag ([`ModificationEntry::UNSCORED`]).
    pub probability: u8,
    /// [`ModificationEntry::MINUS`], [`ModificationEntry::EXPLICIT`] and
    /// [`ModificationEntry::UNSCORED`] bits.
    pub flags: u8,
    pub _pad: [u8; 5],
}

impl ModificationEntry {
    /// The modification is on the strand opposite the sequenced one (`-`).
    pub const MINUS: u8 = 1;
    /// Its `MM` entry is explicit (`?`): bases it does not list were not
    /// examined, rather than unmodified (`.` or no mark).
    pub const EXPLICIT: u8 = 2;
    /// The record has no `ML` tag, so the probability is unknown.
    pub const UNSCORED: u8 = 4;
}

const _: () = assert!(
    std::mem::size_of::<ModificationEntry>() == 24
        && std::mem::align_of::<ModificationEntry>() == 8
);

/// Which variable-length fields a query also returns; the fixed fields of
/// [`AlignmentEntry`] always come back.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct AlignmentFields(u8);

impl AlignmentFields {
    pub const NONE: AlignmentFields = AlignmentFields(0);
    /// Read names, without their NUL.
    pub const NAME: AlignmentFields = AlignmentFields(1);
    /// Bases as ASCII (`=ACMGRSVTWYHKDBN`), one per query position.
    pub const SEQ: AlignmentFields = AlignmentFields(2);
    /// Phred base qualities as stored, one per query position (`0xff`
    /// throughout when the record has none).
    pub const QUAL: AlignmentFields = AlignmentFields(4);
    /// CIGAR operations as BAM stores them: `len << 4 | op`.
    pub const CIGAR: AlignmentFields = AlignmentFields(8);
    /// Optional fields in their binary BAM encoding.
    pub const AUX: AlignmentFields = AlignmentFields(16);
    /// Base modifications decoded from the `MM` / `ML` tags, as
    /// [`ModificationEntry`] rows.
    pub const MODIFICATIONS: AlignmentFields = AlignmentFields(32);
    pub const ALL: AlignmentFields = AlignmentFields(63);

    pub fn from_bits(bits: u8) -> Result<AlignmentFields> {
        if bits & !Self::ALL.0 != 0 {
            return Err(Error::invalid(
                "alignment fields",
                format!("unknown field bits {bits:#x}"),
            ));
        }
        Ok(AlignmentFields(bits))
    }

    pub fn bits(self) -> u8 {
        self.0
    }

    pub fn contains(self, other: AlignmentFields) -> bool {
        self.0 & other.0 == other.0
    }

    pub fn is_empty(self) -> bool {
        self.0 == 0
    }
}

impl BitOr for AlignmentFields {
    type Output = AlignmentFields;

    fn bitor(self, rhs: AlignmentFields) -> AlignmentFields {
        AlignmentFields(self.0 | rhs.0)
    }
}

/// Sizes of an alignment query's output, to allocate before filling. Sizes
/// of fields the query did not ask for are zero.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AlignmentSize {
    pub records: usize,
    pub name_bytes: usize,
    /// Query positions over all records: the bytes of the SEQ and of the
    /// QUAL column.
    pub seq_bytes: usize,
    pub cigar_ops: usize,
    pub aux_bytes: usize,
    pub modifications: usize,
}

/// CIGAR operations of every record: record `i`'s are
/// `ops[offsets[i]..offsets[i + 1]]`.
pub struct CigarColumn<'a> {
    pub offsets: &'a mut [i64],
    pub ops: &'a mut [u32],
}

/// Base modifications of every record: record `i`'s are
/// `entries[offsets[i]..offsets[i + 1]]`, in SEQ order.
pub struct ModificationColumn<'a> {
    pub offsets: &'a mut [i64],
    pub entries: &'a mut [ModificationEntry],
}

/// Caller-allocated output of an alignment query, sized by
/// [`AlignmentSize`]: one buffer for every field the query asked for, and
/// `None` for the others.
pub struct AlignmentsOut<'a> {
    pub entries: &'a mut [AlignmentEntry],
    pub names: Option<StringColumn<'a>>,
    pub seqs: Option<StringColumn<'a>>,
    pub quals: Option<StringColumn<'a>>,
    pub cigars: Option<CigarColumn<'a>>,
    pub aux: Option<StringColumn<'a>>,
    pub modifications: Option<ModificationColumn<'a>>,
}

fn check_strings(
    what: &'static str,
    column: &Option<StringColumn<'_>>,
    wanted: bool,
    rows: usize,
    bytes: usize,
) -> Result<()> {
    match (column, wanted) {
        (None, false) => Ok(()),
        (Some(c), true) => {
            check_output(what, c.offsets.len(), rows + 1)?;
            check_output(what, c.data.len(), bytes)
        }
        (Some(_), false) => Err(Error::invalid(what, "given but not queried")),
        (None, true) => Err(Error::invalid(what, "queried but not given")),
    }
}

impl AlignmentsOut<'_> {
    /// Fails unless the buffers match the query's `fields` and `size`
    /// exactly.
    pub fn check(&self, fields: AlignmentFields, size: AlignmentSize) -> Result<()> {
        let rows = size.records;
        check_output("alignment entries", self.entries.len(), rows)?;
        check_strings(
            "alignment names",
            &self.names,
            fields.contains(AlignmentFields::NAME),
            rows,
            size.name_bytes,
        )?;
        check_strings(
            "alignment sequences",
            &self.seqs,
            fields.contains(AlignmentFields::SEQ),
            rows,
            size.seq_bytes,
        )?;
        check_strings(
            "alignment qualities",
            &self.quals,
            fields.contains(AlignmentFields::QUAL),
            rows,
            size.seq_bytes,
        )?;
        check_strings(
            "alignment aux",
            &self.aux,
            fields.contains(AlignmentFields::AUX),
            rows,
            size.aux_bytes,
        )?;
        match (&self.cigars, fields.contains(AlignmentFields::CIGAR)) {
            (None, false) => {}
            (Some(c), true) => {
                check_output("alignment cigars", c.offsets.len(), rows + 1)?;
                check_output("alignment cigars", c.ops.len(), size.cigar_ops)?;
            }
            (Some(_), false) => {
                return Err(Error::invalid("alignment cigars", "given but not queried"));
            }
            (None, true) => {
                return Err(Error::invalid("alignment cigars", "queried but not given"));
            }
        }
        match (
            &self.modifications,
            fields.contains(AlignmentFields::MODIFICATIONS),
        ) {
            (None, false) => Ok(()),
            (Some(m), true) => {
                check_output("alignment modifications", m.offsets.len(), rows + 1)?;
                check_output(
                    "alignment modifications",
                    m.entries.len(),
                    size.modifications,
                )
            }
            (Some(_), false) => Err(Error::invalid(
                "alignment modifications",
                "given but not queried",
            )),
            (None, true) => Err(Error::invalid(
                "alignment modifications",
                "queried but not given",
            )),
        }
    }
}

/// Alignments (BAM or CRAM). A query returns the records whose reference span
/// overlaps the region, in file order, exactly those htslib's region
/// iterator returns: unmapped records placed on the reference included.
/// Unknown chromosomes are errors; regions may run past a chromosome's
/// end; an empty region yields none.
pub trait AlignmentRead {
    /// Length of `chrom`; [`Error::UnknownChrom`] if absent.
    fn chrom_len(&self, chrom: &str) -> Result<u32>;

    /// Decodes the records of `region` once and holds what `fields` asks
    /// for until the next query. Returns the output sizes, for the caller
    /// to allocate. Base modifications are kept only for the codes in
    /// `modification_codes` ([`ModificationEntry::code`]), or all of them
    /// when it is empty; other records' fields are unaffected.
    fn query_filtered(
        &mut self,
        region: Region<'_>,
        fields: AlignmentFields,
        modification_codes: &[i32],
    ) -> Result<AlignmentSize>;

    /// [`AlignmentRead::query_filtered`] keeping every modification.
    fn query(&mut self, region: Region<'_>, fields: AlignmentFields) -> Result<AlignmentSize> {
        self.query_filtered(region, fields, &[])
    }

    /// Copies the records of the last [`AlignmentRead::query`] into `out`,
    /// sized exactly by it. Only copies (and expands packed bases): nothing
    /// is read again.
    fn write(&self, out: AlignmentsOut<'_>) -> Result<()>;
}
