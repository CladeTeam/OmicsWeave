//! Factor frequency models: a motif window scored by local factors, each a
//! small table over the symbols of a few positions.

use crate::error::{Result, invalid};

/// Highest factor order (positions per factor).
pub const MAX_ORDER: usize = 16;

/// Most cells one factor table may hold.
pub const MAX_TABLE: usize = 1 << 24;

/// Motifs addressable by the four-column u32 index stream. The IDs stream
/// imposes an additional limit on the combined offsets and ID bytes.
pub const MAX_MOTIFS: usize = (u32::MAX / 4) as usize;

/// How a factor's table is read.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum FactorKind {
    /// Joint frequency of the symbol combination at the factor's positions;
    /// normalized over the whole table.
    Joint,
    /// Frequency of the last position's symbol given the earlier positions;
    /// normalized over each run of `arity` cells sharing a context.
    Conditional,
    /// Unnormalized weight; the score contribution is `ln(value)`, and a zero
    /// forbids the combination.
    Potential,
}

impl FactorKind {
    pub(crate) fn id(self) -> u32 {
        match self {
            FactorKind::Joint => 0,
            FactorKind::Conditional => 1,
            FactorKind::Potential => 2,
        }
    }

    pub(crate) fn from_id(id: u32) -> Option<FactorKind> {
        match id {
            0 => Some(FactorKind::Joint),
            1 => Some(FactorKind::Conditional),
            2 => Some(FactorKind::Potential),
            _ => None,
        }
    }
}

/// The symbols a motif models individually.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Support {
    /// Every symbol of the alphabet; there is no rest category.
    All,
    /// A strict subset of global codes, ascending. The other symbols of the
    /// alphabet form the rest category, the last local code.
    Codes(Vec<u8>),
}

impl Support {
    pub fn codes(&self) -> Option<&[u8]> {
        match self {
            Support::All => None,
            Support::Codes(c) => Some(c),
        }
    }
}

/// Local codes per position: the support, plus one for the rest category
/// when the support does not cover the alphabet.
pub fn arity(support: Option<&[u8]>, alphabet_len: usize) -> usize {
    match support {
        None => alphabet_len,
        Some(codes) => codes.len() + 1,
    }
}

/// Cells of a table over `order` positions of `arity` local codes.
pub fn table_len(arity: usize, order: usize) -> Result<usize> {
    if !(1..=crate::MAX_SYMBOLS).contains(&arity) || !(1..=MAX_ORDER).contains(&order) {
        return Err(invalid(
            "FFM factor",
            "arity or order outside format limits",
        ));
    }
    let mut n = 1usize;
    for _ in 0..order {
        n = n.saturating_mul(arity);
    }
    if n > MAX_TABLE {
        return Err(invalid(
            "FFM factor",
            format!("table of {arity}^{order} cells exceeds {MAX_TABLE}"),
        ));
    }
    Ok(n)
}

#[derive(Clone, Debug, PartialEq)]
pub struct Factor {
    pub kind: FactorKind,
    /// Strictly increasing positions within the motif window.
    pub positions: Vec<u32>,
    /// `arity^order` cells; the last position varies fastest.
    pub table: Vec<f64>,
}

/// Gap counts of a gapped profile: an alignment of `sites` sequences whose
/// `width` columns are the motif positions, with residues inserted between
/// columns and columns deleted in some sequences.
#[derive(Clone, Debug, PartialEq)]
pub struct Gaps {
    /// Per column: sequences in which the column is deleted (`width` values).
    pub deletions: Vec<f64>,
    /// Per junction between column `i` and `i + 1`: residues inserted there
    /// over all sequences (`width − 1` values).
    pub insertions: Vec<f64>,
}

/// One motif's model. Its ID and annotation are stored beside it, not in it.
#[derive(Clone, Debug, PartialEq)]
pub struct Motif {
    pub width: u32,
    pub support: Support,
    /// Observations behind the frequencies, when known.
    pub sites: Option<f64>,
    pub factors: Vec<Factor>,
    /// Present for a gapped profile, whose factors are then one
    /// single-position joint factor per column, in column order, holding
    /// the residue counts of the column.
    pub gaps: Option<Gaps>,
}

impl Motif {
    /// A position frequency matrix: one joint factor per row of `rows`, each
    /// row holding a frequency or count per symbol of the alphabet.
    pub fn pfm(rows: &[Vec<f64>], sites: Option<f64>) -> Motif {
        Motif {
            width: rows.len() as u32,
            support: Support::All,
            sites,
            factors: rows
                .iter()
                .enumerate()
                .map(|(p, row)| Factor {
                    kind: FactorKind::Joint,
                    positions: vec![p as u32],
                    table: row.clone(),
                })
                .collect(),
            gaps: None,
        }
    }

    /// Local codes per position under an alphabet of `alphabet_len` symbols.
    pub fn arity(&self, alphabet_len: usize) -> usize {
        arity(self.support.codes(), alphabet_len)
    }

    /// Checks every constraint the file format places on a model.
    pub fn validate(&self, alphabet_len: usize) -> Result<()> {
        let support = self.support.codes();
        check_header(self.width, support, self.sites, alphabet_len)?;
        if self.factors.is_empty() {
            return Err(invalid("FFM motif", "no factors"));
        }
        let r = arity(support, alphabet_len);
        for (i, f) in self.factors.iter().enumerate() {
            let deleted = self
                .gaps
                .as_ref()
                .is_some_and(|g| g.deletions.get(i).is_some_and(|&d| d > 0.0));
            check_factor(f.kind, &f.positions, &f.table, self.width, r, deleted)?;
        }
        if let Some(g) = &self.gaps {
            check_gaps(
                self.width,
                &g.deletions,
                &g.insertions,
                self.factors
                    .iter()
                    .map(|f| (f.kind, f.positions.as_slice())),
            )?;
        }
        Ok(())
    }
}

/// Checks a gapped profile's counts and that its factors are its columns.
pub(crate) fn check_gaps<'a>(
    width: u32,
    deletions: &[f64],
    insertions: &[f64],
    factors: impl ExactSizeIterator<Item = (FactorKind, &'a [u32])>,
) -> Result<()> {
    let w = width as usize;
    if deletions.len() != w || insertions.len() != w - 1 {
        return Err(invalid(
            "FFM gaps",
            format!(
                "{} deletion and {} insertion counts for width {w}, need {w} and {}",
                deletions.len(),
                insertions.len(),
                w - 1
            ),
        ));
    }
    if let Some(v) = deletions
        .iter()
        .chain(insertions)
        .find(|v| !(v.is_finite() && **v >= 0.0))
    {
        return Err(invalid(
            "FFM gaps",
            format!("count {v} is not finite and >= 0"),
        ));
    }
    let columns = factors.len() == w;
    if !columns
        || !factors
            .enumerate()
            .all(|(i, (kind, positions))| kind == FactorKind::Joint && positions == [i as u32])
    {
        return Err(invalid(
            "FFM gaps",
            "a gapped profile has one single-position joint factor per column, in order",
        ));
    }
    Ok(())
}

pub(crate) fn check_header(
    width: u32,
    support: Option<&[u8]>,
    sites: Option<f64>,
    alphabet_len: usize,
) -> Result<()> {
    if !(1..=crate::MAX_SYMBOLS).contains(&alphabet_len) {
        return Err(invalid("FFM motif", "alphabet size must be in 1..=255"));
    }
    if width == 0 {
        return Err(invalid("FFM motif", "width is 0"));
    }
    if let Some(codes) = support {
        check_support(codes, alphabet_len)?;
    }
    if let Some(s) = sites
        && !(s.is_finite() && s > 0.0)
    {
        return Err(invalid(
            "FFM motif",
            format!("sites {s} is not finite and > 0"),
        ));
    }
    Ok(())
}

pub(crate) fn check_support<T: Copy + Into<u32>>(codes: &[T], alphabet_len: usize) -> Result<()> {
    if codes.is_empty() || codes.len() >= alphabet_len {
        return Err(invalid(
            "FFM motif",
            format!(
                "support of {} codes; a subset needs 1..{alphabet_len} (use all for the whole alphabet)",
                codes.len()
            ),
        ));
    }
    if codes.windows(2).any(|w| w[0].into() >= w[1].into())
        || codes[codes.len() - 1].into() as usize >= alphabet_len
    {
        return Err(invalid(
            "FFM motif",
            "support codes must be ascending and inside the alphabet",
        ));
    }
    Ok(())
}

/// Checks one factor. `deleted`: the factor is a column of a gapped profile
/// that some sequences delete, which may then have no residues at all.
pub(crate) fn check_factor(
    kind: FactorKind,
    positions: &[u32],
    table: &[f64],
    width: u32,
    arity: usize,
    deleted: bool,
) -> Result<()> {
    let order = positions.len();
    if order == 0 || order > MAX_ORDER {
        return Err(invalid(
            "FFM factor",
            format!("order {order}, need 1..={MAX_ORDER}"),
        ));
    }
    if positions.windows(2).any(|w| w[0] >= w[1]) || positions[order - 1] >= width {
        return Err(invalid(
            "FFM factor",
            format!("positions {positions:?} are not ascending inside width {width}"),
        ));
    }
    let need = table_len(arity, order)?;
    if table.len() != need {
        return Err(invalid(
            "FFM factor",
            format!(
                "table has {} cells, need {arity}^{order} = {need}",
                table.len()
            ),
        ));
    }
    if let Some(v) = table.iter().find(|v| !(v.is_finite() && **v >= 0.0)) {
        return Err(invalid(
            "FFM factor",
            format!("value {v} is not finite and >= 0"),
        ));
    }
    if deleted && table.iter().all(|&v| v == 0.0) {
        return Ok(());
    }
    let groups = match kind {
        FactorKind::Joint => table.len(),
        FactorKind::Conditional => arity,
        FactorKind::Potential => return Ok(()),
    };
    // Non-negative finite cells have positive total mass iff one is positive;
    // testing that directly also accepts counts whose sum overflows f64.
    if table.chunks(groups).any(|g| !g.iter().any(|&v| v > 0.0)) {
        return Err(invalid(
            "FFM factor",
            format!("a {kind:?} normalization group sums to 0"),
        ));
    }
    Ok(())
}
