//! [`IntervalSet`]: a columnar set of stranded intervals with an owner and
//! a name per row, its own [`IntervalIndex`] and [`SweepCursor`]. Derived
//! views of an annotation (introns, exonic regions, promoters, ...) are
//! interval sets; so can be peaks, fixed bins or any other features a
//! counting layer indexes.

use std::sync::{Arc, OnceLock};

use omicsweave_core::{Error, Result};

use crate::index::{IntervalIndex, SweepCursor};
use crate::model::Strand;
use crate::strings::{IdMap, NONE, StringPool, id_map};

/// What the rows of an interval set belong to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Level {
    Gene,
    Transcript,
    /// Rows of [`Features`](crate::Features) (exons).
    Exon,
}

impl Level {
    pub fn name(self) -> &'static str {
        match self {
            Level::Gene => "gene",
            Level::Transcript => "transcript",
            Level::Exon => "exon",
        }
    }

    pub fn parse(s: &str) -> Option<Level> {
        Some(match s {
            "gene" => Level::Gene,
            "transcript" => Level::Transcript,
            "exon" => Level::Exon,
            _ => return None,
        })
    }
}

/// Per-row lists in CSR layout: row `i`'s values are
/// `values[offsets[i]..offsets[i + 1]]`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Csr {
    pub offsets: Vec<u32>,
    pub values: Vec<u32>,
}

impl Csr {
    pub fn new() -> Csr {
        Csr {
            offsets: vec![0],
            values: Vec::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.offsets.len().saturating_sub(1)
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    #[inline]
    pub fn get(&self, i: u32) -> &[u32] {
        &self.values[self.offsets[i as usize] as usize..self.offsets[i as usize + 1] as usize]
    }

    pub(crate) fn push_row(&mut self, values: impl IntoIterator<Item = u32>) {
        self.values.extend(values);
        self.offsets.push(
            u32::try_from(self.values.len()).expect("more than 2^32 derived interval memberships"),
        );
    }

    pub fn heap_bytes(&self) -> usize {
        4 * (self.offsets.capacity() + self.values.capacity())
    }
}

/// Stranded intervals as columns: row `i` is `seqid[i]`,
/// `[start[i], end[i])` (0-based, half-open), `strand[i]`, belonging to
/// `owner[i]` (a row of the [`Level`] the set names, or [`NONE`]) and
/// called `name[i]` (a string id, or [`NONE`]).
///
/// Seqids number the set's chromosome names ([`IntervalSet::seqid_name`]);
/// sets derived from an annotation share its numbering and string pool.
/// The index is built on first use.
#[derive(Debug, Clone)]
pub struct IntervalSet {
    pub(crate) pool: Arc<StringPool>,
    pub(crate) seqid_names: Arc<[u32]>,
    pub(crate) owner_level: Option<Level>,
    pub seqid: Vec<u32>,
    pub start: Vec<u32>,
    pub end: Vec<u32>,
    pub strand: Vec<Strand>,
    pub owner: Vec<u32>,
    pub name: Vec<u32>,
    index: OnceLock<IntervalIndex>,
}

impl IntervalSet {
    pub(crate) fn empty(
        pool: Arc<StringPool>,
        seqid_names: Arc<[u32]>,
        owner_level: Option<Level>,
    ) -> IntervalSet {
        IntervalSet {
            pool,
            seqid_names,
            owner_level,
            seqid: Vec::new(),
            start: Vec::new(),
            end: Vec::new(),
            strand: Vec::new(),
            owner: Vec::new(),
            name: Vec::new(),
            index: OnceLock::new(),
        }
    }

    /// A set from columns: `seqid[i]` numbers `chroms`; `owner` and
    /// `names` may be absent (all [`NONE`]).
    #[allow(clippy::too_many_arguments)]
    pub fn from_columns(
        chroms: &[&str],
        seqid: Vec<u32>,
        start: Vec<u32>,
        end: Vec<u32>,
        strand: Vec<Strand>,
        owner: Option<Vec<u32>>,
        names: Option<&[&str]>,
        owner_level: Option<Level>,
    ) -> Result<IntervalSet> {
        let n = start.len();
        let mut pool = StringPool::new();
        let mut seqid_names = Vec::with_capacity(chroms.len());
        for c in chroms {
            let id = pool.intern(c.as_bytes())?;
            // Only chromosome names have been interned in this fresh pool;
            // a new name receives the next dense ID in constant time.
            if id as usize != seqid_names.len() {
                return Err(Error::invalid(
                    "interval set",
                    format!("chromosome {c:?} twice"),
                ));
            }
            seqid_names.push(id);
        }
        let owner = owner.unwrap_or_else(|| vec![NONE; n]);
        if seqid.len() != n || end.len() != n || strand.len() != n || owner.len() != n {
            return Err(Error::invalid("interval set", "columns differ in length"));
        }
        if let Some(i) = (0..n).find(|&i| start[i] > end[i]) {
            return Err(Error::InvalidRegion {
                start: start[i],
                end: end[i],
            });
        }
        if seqid.iter().any(|&s| s as usize >= chroms.len()) {
            return Err(Error::invalid(
                "interval set",
                "seqid outside the chromosome list",
            ));
        }
        let name = match names {
            Some(v) if v.len() != n => {
                return Err(Error::invalid("interval set", "names differ in length"));
            }
            Some(v) => v
                .iter()
                .map(|s| pool.intern(s.as_bytes()))
                .collect::<Result<Vec<u32>>>()?,
            None => vec![NONE; n],
        };
        Ok(IntervalSet {
            pool: Arc::new(pool),
            seqid_names: Arc::from(seqid_names),
            owner_level,
            seqid,
            start,
            end,
            strand,
            owner,
            name,
            index: OnceLock::new(),
        })
    }

    #[inline]
    pub(crate) fn push(
        &mut self,
        seqid: u32,
        start: u32,
        end: u32,
        strand: Strand,
        owner: u32,
        name: u32,
    ) {
        self.seqid.push(seqid);
        self.start.push(start);
        self.end.push(end);
        self.strand.push(strand);
        self.owner.push(owner);
        self.name.push(name);
    }

    pub fn len(&self) -> usize {
        self.start.len()
    }

    pub fn is_empty(&self) -> bool {
        self.start.is_empty()
    }

    /// What `owner` refers to.
    pub fn owner_level(&self) -> Option<Level> {
        self.owner_level
    }

    /// Chromosome names known to the set.
    pub fn seqid_count(&self) -> usize {
        self.seqid_names.len()
    }

    pub fn seqid_name(&self, seqid: u32) -> &str {
        self.pool.str(self.seqid_names[seqid as usize])
    }

    /// The seqid called `name`.
    pub fn seqid_of(&self, name: &str) -> Option<u32> {
        let id = self.pool.get(name.as_bytes())?;
        self.seqid_names
            .iter()
            .position(|&n| n == id)
            .map(|i| i as u32)
    }

    /// The name of row `i`.
    pub fn name(&self, i: usize) -> Option<&str> {
        self.pool.get_str(self.name[i])
    }

    /// The strings behind the name column.
    pub fn pool(&self) -> &StringPool {
        &self.pool
    }

    /// The overlap index over the rows (ids are row numbers).
    /// After editing public columns, call [`IntervalSet::clear_index`]
    /// before querying again; edits are not detected automatically.
    pub fn index(&self) -> &IntervalIndex {
        self.index.get_or_init(|| {
            IntervalIndex::build(self.seqid_names.len(), &self.seqid, &self.start, &self.end)
        })
    }

    /// Releases the cached index. The next query builds it from the
    /// current columns; call this after editing coordinates or seqids.
    /// Existing cursors borrow the set and must be dropped first.
    pub fn clear_index(&mut self) {
        self.index.take();
    }

    /// Appends the rows overlapping `[start, end)` of `seqid`.
    pub fn overlaps(&self, seqid: u32, start: u32, end: u32, out: &mut Vec<u32>) {
        self.index().overlaps(seqid, start, end, out);
    }

    /// A sweep cursor over `seqid`.
    pub fn cursor(&self, seqid: u32) -> SweepCursor<'_> {
        self.index().cursor(seqid)
    }

    /// The rows `rows` (in that order) as a new set.
    pub fn select(&self, rows: &[u32]) -> IntervalSet {
        let mut s = IntervalSet::empty(
            Arc::clone(&self.pool),
            Arc::clone(&self.seqid_names),
            self.owner_level,
        );
        for &r in rows {
            let r = r as usize;
            s.push(
                self.seqid[r],
                self.start[r],
                self.end[r],
                self.strand[r],
                self.owner[r],
                self.name[r],
            );
        }
        s
    }

    /// The rows whose owner passes `keep`.
    pub fn filter_owners(&self, mut keep: impl FnMut(u32) -> bool) -> IntervalSet {
        let rows: Vec<u32> = (0..self.len() as u32)
            .filter(|&r| keep(self.owner[r as usize]))
            .collect();
        self.select(&rows)
    }

    pub fn heap_bytes(&self) -> usize {
        4 * (self.seqid.capacity()
            + self.start.capacity()
            + self.end.capacity()
            + self.owner.capacity()
            + self.name.capacity())
            + self.strand.capacity()
            + self.index.get().map_or(0, IntervalIndex::heap_bytes)
    }
}

/// Builds an [`IntervalSet`] row by row from names (peaks, bins, features
/// of any source).
#[derive(Debug, Default)]
pub struct IntervalSetBuilder {
    pool: StringPool,
    seqids: Vec<u32>,
    seqid_of: IdMap<u32, u32>,
    rows: Option<IntervalSet>,
    owner_level: Option<Level>,
}

impl IntervalSetBuilder {
    pub fn new() -> IntervalSetBuilder {
        IntervalSetBuilder {
            seqid_of: id_map(64),
            ..IntervalSetBuilder::default()
        }
    }

    /// Declares chromosomes up front, fixing their seqid numbering.
    pub fn chroms<'a>(mut self, names: impl IntoIterator<Item = &'a str>) -> Result<Self> {
        for n in names {
            self.seqid(n)?;
        }
        Ok(self)
    }

    /// What owners refer to.
    pub fn owner_level(mut self, level: Option<Level>) -> Self {
        self.owner_level = level;
        self
    }

    /// The seqid of `chrom`, numbering it if new.
    pub fn seqid(&mut self, chrom: &str) -> Result<u32> {
        let id = self.pool.intern(chrom.as_bytes())?;
        let n = self.seqids.len() as u32;
        let s = *self.seqid_of.entry(id).or_insert(n);
        if s == n {
            self.seqids.push(id);
        }
        Ok(s)
    }

    /// Adds `chrom:[start, end)`; `owner` may be [`NONE`].
    pub fn push(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
        strand: Strand,
        owner: u32,
        name: Option<&str>,
    ) -> Result<()> {
        if start > end {
            return Err(Error::InvalidRegion { start, end });
        }
        let seqid = self.seqid(chrom)?;
        let name = match name {
            Some(n) => self.pool.intern(n.as_bytes())?,
            None => NONE,
        };
        let rows = self.rows.get_or_insert_with(|| {
            IntervalSet::empty(Arc::new(StringPool::new()), Arc::from(Vec::new()), None)
        });
        rows.push(seqid, start, end, strand, owner, name);
        Ok(())
    }

    pub fn build(self) -> IntervalSet {
        let mut set = self.rows.unwrap_or_else(|| {
            IntervalSet::empty(Arc::new(StringPool::new()), Arc::from(Vec::new()), None)
        });
        set.pool = Arc::new(self.pool);
        set.seqid_names = Arc::from(self.seqids);
        set.owner_level = self.owner_level;
        set
    }
}
