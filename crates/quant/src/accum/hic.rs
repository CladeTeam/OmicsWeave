//! Hi-C / Micro-C contacts from a coordinate-sorted BAM whose records
//! carry their mate's CIGAR (`MC`) and MAPQ (`MQ`), as `samtools fixmate
//! -m` (or bwa) writes them.
//!
//! Definitions (4DN pairs, pairtools):
//! - A **contact** is a read pair passing the filter (default: paired,
//!   both mates mapped, primary, not QC-failed, supplementary or
//!   duplicate), taken from its leftmost mate in coordinate order (at one
//!   position, read 1) with the mate's position, strand, `MC` and `MQ`:
//!   no mate is buffered. Its two ends are the mates' 5' positions
//!   (1-based; a reverse mate's from its `MC` CIGAR), ordered by reference
//!   index then position (`chr1`, `pos1` first).
//! - **MAPQ**: both mates at least `min_mapq` (the mate's from `MQ`; a
//!   record without `MQ` passes only when `min_mapq` is 0).
//! - Within a reference, contacts closer than `min_distance` are dropped;
//!   with `self_ligation`, so are inward (`+ -`) and outward (`- +`)
//!   contacts closer than it — dangling ends and self-circles; same-strand
//!   contacts cannot be either and are kept.
//! - **Duplicates**: contacts with the same two ends and strands (and
//!   barcode) are one; the first in file order is kept.
//! - Per barcode (`barcode_tag`): pairs, duplicates, contacts, cis /
//!   trans, cis at 1 kb, 10 kb and 20 kb or more (pairtools stats).
//!
//! Products: `<name>/contacts.pairs.gz` (BGZF; header with `#chromsize`
//! lines in reference order, sorted `chr1-chr2-pos1-pos2`), the Juicer
//! `<name>/contacts.hic` (v9, raw counts) binned at `resolutions` as the
//! Juicer tools' `pre` bins the pairs, and `<name>/cell_metrics.csv`.

use std::io::Read;
use std::path::PathBuf;

use omicsweave_bam::flag;
use omicsweave_convert::stage::{Sorter, Spill, read_or_end};
use omicsweave_core::{Error, Result};

use super::counts::{GlobalRows, LocalRows, Rows};
use super::coverage::Tally;
use super::{Accumulator, BatchInfo, Delivery, MemoryNeed, Needs, PlanInfo, RunCtx};
use crate::model::{BarcodeMetrics, BarcodeTable, Contact, ContactTable};
use crate::prelude::{Filter, ReadView, coord_key};
use crate::sink::json::Json;

/// A contact product; see the module documentation.
#[derive(Debug, Clone)]
pub struct Contacts {
    /// The output directory (`hic`).
    pub name: String,
    pub filter: Filter,
    /// Contacts within a reference closer than this are dropped (0: none).
    pub min_distance: u32,
    /// Inward and outward contacts closer than this are dropped (0: none).
    pub self_ligation: u32,
    pub dedup: bool,
    /// Write `contacts.pairs.gz`.
    pub pairs: bool,
    /// `.hic` base-pair resolutions (none: no `.hic`).
    pub resolutions: Vec<u32>,
    /// The genome assembly named in both files.
    pub genome: String,
    /// Per-barcode statistics (and a `barcode` column) by this tag.
    pub barcode_tag: Option<[u8; 2]>,
    /// Only these barcodes, rows in this order.
    pub whitelist: Option<PathBuf>,
    /// Also write the `.hic` file's `.sidecar`.
    pub sidecar: bool,
}

impl Default for Contacts {
    /// MAPQ 30 on both mates, duplicates removed, no distance filter, pairs
    /// and a `.hic` at the Juicer tools' resolutions.
    fn default() -> Contacts {
        Contacts {
            name: "hic".into(),
            filter: Filter {
                require: flag::PAIRED,
                exclude: flag::UNMAPPED
                    | flag::MATE_UNMAPPED
                    | flag::SECONDARY
                    | flag::QC_FAIL
                    | flag::SUPPLEMENTARY
                    | flag::DUPLICATE,
                min_mapq: 30,
                unique: false,
            },
            min_distance: 0,
            self_ligation: 0,
            dedup: true,
            pairs: true,
            resolutions: omicsweave_convert::JUICER_RESOLUTIONS.to_vec(),
            genome: "unknown".into(),
            barcode_tag: None,
            whitelist: None,
            sidecar: false,
        }
    }
}

impl Spill for Contact {
    type Key = (u32, u32, u32, u32, bool, bool, u32, u32, u64);
    fn key(&self) -> Self::Key {
        (
            self.c1,
            self.c2,
            self.p1,
            self.p2,
            self.s1,
            self.s2,
            self.row,
            self.input,
            self.ordinal,
        )
    }
    fn mem_size(&self) -> usize {
        self.name.capacity()
    }
    fn encode(&self, out: &mut Vec<u8>) {
        for v in [self.c1, self.p1, self.c2, self.p2, self.row, self.input] {
            out.extend_from_slice(&v.to_le_bytes());
        }
        out.extend_from_slice(&self.ordinal.to_le_bytes());
        out.extend_from_slice(&[self.s1 as u8 | (self.s2 as u8) << 1, self.mapq1, self.mapq2]);
        out.extend_from_slice(&(self.name.len() as u16).to_le_bytes());
        out.extend_from_slice(&self.name);
    }
    fn decode(input: &mut impl Read) -> Result<Option<Contact>> {
        let mut b = [0u8; 37];
        if !read_or_end(input, &mut b)? {
            return Ok(None);
        }
        let u = |i: usize| u32::from_le_bytes(b[4 * i..4 * i + 4].try_into().unwrap());
        let n = u16::from_le_bytes([b[35], b[36]]) as usize;
        let mut name = vec![0u8; n];
        input.read_exact(&mut name)?;
        Ok(Some(Contact {
            c1: u(0),
            p1: u(1),
            c2: u(2),
            p2: u(3),
            row: u(4),
            input: u(5),
            ordinal: u64::from_le_bytes(b[24..32].try_into().unwrap()),
            s1: b[32] & 1 != 0,
            s2: b[32] & 2 != 0,
            mapq1: b[33],
            mapq2: b[34],
            name,
        }))
    }
}

/// Why pairs were or were not taken.
#[derive(Debug, Clone, Copy, Default)]
pub struct ContactStatus {
    pub pairs: u64,
    /// Mates whose pair is taken from the other mate.
    pub not_leader: u64,
    /// Mates on no reference.
    pub mate_unplaced: u64,
    pub no_mate_mapq: u64,
    pub mate_mapq: u64,
    /// Reverse mates without `MC`, whose 5' end is unknown.
    pub no_mate_cigar: u64,
    pub too_close: u64,
    pub self_ligation: u64,
    pub no_barcode: u64,
}

impl ContactStatus {
    fn add(&mut self, o: &ContactStatus) {
        self.pairs += o.pairs;
        self.not_leader += o.not_leader;
        self.mate_unplaced += o.mate_unplaced;
        self.no_mate_mapq += o.no_mate_mapq;
        self.mate_mapq += o.mate_mapq;
        self.no_mate_cigar += o.no_mate_cigar;
        self.too_close += o.too_close;
        self.self_ligation += o.self_ligation;
        self.no_barcode += o.no_barcode;
    }
}

/// A [`Contacts`] product resolved against a task.
pub struct ContactStage {
    pub(crate) spec: Contacts,
    pub(crate) rows: Option<Rows>,
    pub(crate) table_bytes: u64,
}

pub struct ContactBatch {
    contacts: Vec<Contact>,
    local: LocalRows,
    tally: Tally,
    status: ContactStatus,
}

pub struct ContactTotal<'p> {
    sorter: Sorter<'p, Contact>,
    rows: GlobalRows,
    tally: Tally,
    status: ContactStatus,
}

/// The reference length of a CIGAR string (`M`, `D`, `N`, `=`, `X`).
fn cigar_ref_len(text: &[u8]) -> Option<u32> {
    let mut n = 0u32;
    let mut len = 0u32;
    let mut digits = false;
    for &c in text {
        if c.is_ascii_digit() {
            n = n.checked_mul(10)?.checked_add((c - b'0') as u32)?;
            digits = true;
            continue;
        }
        if !digits {
            return None;
        }
        match c {
            b'M' | b'D' | b'N' | b'=' | b'X' => len = len.checked_add(n)?,
            b'I' | b'S' | b'H' | b'P' => {}
            _ => return None,
        }
        n = 0;
        digits = false;
    }
    (!digits).then_some(len)
}

/// Per-barcode statistics, in [`CELL_METRICS`] order.
const CELL_METRICS: [&str; 8] = [
    "pairs",
    "duplicates",
    "contacts",
    "cis",
    "trans",
    "cis_1kb+",
    "cis_10kb+",
    "cis_20kb+",
];

impl Accumulator for ContactStage {
    type Batch = ContactBatch;
    type Total<'p> = ContactTotal<'p>;
    type Product = (ContactTable, Tally, ContactStatus, u64);

    fn name(&self) -> &str {
        &self.spec.name
    }

    fn needs(&self) -> Needs {
        Needs::CIGAR | Needs::TAGS
    }

    fn estimate(&self, info: &PlanInfo) -> MemoryNeed {
        MemoryNeed {
            fixed: self.table_bytes,
            // A contact of about 80 bytes for two records.
            per_batch: info.batch_bytes / 4,
            min_buffer: 4 << 20,
            weight: 4,
        }
    }

    fn batch(&self, _info: BatchInfo) -> ContactBatch {
        ContactBatch {
            contacts: Vec::new(),
            local: LocalRows::default(),
            tally: Tally::default(),
            status: ContactStatus::default(),
        }
    }

    fn observe(&self, b: &mut ContactBatch, r: &ReadView<'_>) -> Result<()> {
        b.tally.seen += 1;
        if let Some(why) = self.spec.filter.reject(r) {
            b.tally.reject(why);
            return Ok(());
        }
        let f = r.flag();
        if f & flag::PAIRED == 0 || f & (flag::UNMAPPED | flag::MATE_UNMAPPED) != 0 {
            b.status.mate_unplaced += 1;
            return Ok(());
        }
        let rec = &r.rec;
        if rec.mate_tid < 0 || rec.mate_pos < 0 {
            b.status.mate_unplaced += 1;
            return Ok(());
        }
        let me = coord_key(rec.tid, rec.pos);
        let mate = coord_key(rec.mate_tid, rec.mate_pos);
        if me > mate || (me == mate && f & flag::READ1 == 0) {
            b.status.not_leader += 1;
            return Ok(());
        }
        let mq = r.tag(*b"MQ").and_then(|a| a.int());
        let mapq2 = match mq {
            Some(q) => {
                if q < self.spec.filter.min_mapq as i64 {
                    b.status.mate_mapq += 1;
                    return Ok(());
                }
                q.clamp(0, 255) as u8
            }
            None if self.spec.filter.min_mapq > 0 => {
                b.status.no_mate_mapq += 1;
                return Ok(());
            }
            None => 255,
        };
        let mate_reverse = f & flag::MATE_REVERSE != 0;
        let mate_five = if mate_reverse {
            match r.tag_bytes(*b"MC").and_then(cigar_ref_len) {
                Some(l) => rec.mate_pos as u32 + l.max(1),
                None => {
                    b.status.no_mate_cigar += 1;
                    return Ok(());
                }
            }
        } else {
            rec.mate_pos as u32 + 1
        };
        let own = (r.tid(), r.five_prime() + 1, r.is_reverse(), r.mapq());
        let other = (rec.mate_tid as u32, mate_five, mate_reverse, mapq2);
        let (a, z) = if (own.0, own.1) <= (other.0, other.1) {
            (own, other)
        } else {
            (other, own)
        };
        if a.0 == z.0 {
            let d = z.1 - a.1;
            if d < self.spec.min_distance {
                b.status.too_close += 1;
                return Ok(());
            }
            if a.2 != z.2 && d < self.spec.self_ligation {
                b.status.self_ligation += 1;
                return Ok(());
            }
        }
        let row = match &self.rows {
            Some(rows) => match rows.of(r, &mut b.local) {
                Some(x) => x,
                None => {
                    b.status.no_barcode += 1;
                    return Ok(());
                }
            },
            None => u32::MAX,
        };
        b.status.pairs += 1;
        b.contacts.push(Contact {
            c1: a.0,
            p1: a.1,
            c2: z.0,
            p2: z.1,
            s1: a.2,
            s2: z.2,
            mapq1: a.3,
            mapq2: z.3,
            row,
            input: r.input,
            ordinal: r.ordinal,
            name: if self.spec.pairs {
                r.name().to_vec()
            } else {
                Vec::new()
            },
        });
        Ok(())
    }

    fn total<'p>(&'p self, rc: &RunCtx<'p>) -> Result<ContactTotal<'p>> {
        Ok(ContactTotal {
            sorter: Sorter::new(rc.pool, rc.buffer, &rc.tmp),
            rows: GlobalRows::default(),
            tally: Tally::default(),
            status: ContactStatus::default(),
        })
    }

    fn merge<'p>(
        &'p self,
        t: &mut ContactTotal<'p>,
        b: ContactBatch,
        _frontier: u64,
        rc: &RunCtx<'p>,
    ) -> Result<()> {
        t.tally.add(&b.tally);
        t.status.add(&b.status);
        let ContactBatch {
            mut contacts,
            local,
            ..
        } = b;
        if self.rows.as_ref().is_some_and(Rows::local) {
            let map = t.rows.remap(local);
            for c in &mut contacts {
                c.row = map[c.row as usize];
            }
            if t.rows.bytes() > rc.buffer.max(64 << 20) {
                return Err(Error::invalid(
                    "budget",
                    format!(
                        "{}: {} barcodes outgrow the memory budget; give a whitelist or raise \
                         the memory",
                        self.spec.name,
                        t.rows.names.len()
                    ),
                ));
            }
        }
        t.sorter.extend(contacts)
    }

    fn finish<'p>(
        &'p self,
        t: ContactTotal<'p>,
        rc: &RunCtx<'p>,
    ) -> Result<(ContactTable, Tally, ContactStatus, u64)> {
        let spilled = t.sorter.spilled();
        let names = match self.rows.as_ref().and_then(Rows::fixed_names) {
            Some(n) => n.to_vec(),
            None => t.rows.names,
        };
        Ok((
            ContactTable {
                chroms: rc.inputs[0].refs.clone(),
                barcodes: BarcodeTable { names },
                contacts: t.sorter.finish()?,
            },
            t.tally,
            t.status,
            spilled,
        ))
    }

    fn deliver(
        &self,
        (table, tally, status, spilled): (ContactTable, Tally, ContactStatus, u64),
        d: &mut Delivery<'_>,
    ) -> Result<Json> {
        let ContactTable {
            chroms,
            barcodes,
            mut contacts,
        } = table;
        let name = &self.spec.name;
        let mut files = Vec::new();
        let barcode_column = self.rows.is_some();
        let mut pairs = match self.spec.pairs {
            true => {
                let rel = format!("{name}/contacts.pairs.gz");
                let path = d.path(&rel)?;
                files.push(Json::from(rel));
                Some((
                    crate::sink::pairs::PairsWriter::create(
                        d,
                        &path,
                        &chroms,
                        &self.spec.genome,
                        barcode_column,
                    )?,
                    path,
                ))
            }
            false => None,
        };
        let mut hic = match self.spec.resolutions.is_empty() {
            false => Some(crate::sink::hic::HicPixels::new(
                &chroms,
                &self.spec.genome,
                &self.spec.resolutions,
                d,
            )?),
            true => None,
        };
        let mut cells: Vec<[u64; 8]> = vec![[0; 8]; barcodes.len()];
        let mut totals = [0u64; 8];
        let mut last: Option<(u32, u32, u32, u32, bool, bool, u32)> = None;
        while let Some(c) = contacts.next_record()? {
            let mut m = [0u64; 8];
            m[0] = 1;
            let ends = c.ends();
            if self.spec.dedup && last == Some(ends) {
                m[1] = 1;
            } else {
                last = Some(ends);
                m[2] = 1;
                if c.c1 == c.c2 {
                    let dist = c.p2 - c.p1;
                    m[3] = 1;
                    m[5] = (dist >= 1_000) as u64;
                    m[6] = (dist >= 10_000) as u64;
                    m[7] = (dist >= 20_000) as u64;
                } else {
                    m[4] = 1;
                }
                if let Some((w, _)) = &mut pairs {
                    let bc = (c.row != u32::MAX).then(|| barcodes.names[c.row as usize].as_slice());
                    w.push(&c, bc)?;
                }
                if let Some(h) = &mut hic {
                    h.push(c.c1, c.p1, c.c2, c.p2)?;
                }
            }
            for (t, v) in totals.iter_mut().zip(m) {
                *t += v;
            }
            if c.row != u32::MAX {
                for (t, v) in cells[c.row as usize].iter_mut().zip(m) {
                    *t += v;
                }
            }
        }
        if let Some((w, _)) = pairs {
            w.finish()?;
        }
        let mut pixels = Json::Null;
        if let Some(h) = hic {
            let rel = format!("{name}/contacts.hic");
            let path = d.path(&rel)?;
            let n = h.finish(&path, self.spec.sidecar, d)?;
            pixels = Json::from(n);
            files.push(Json::from(rel.clone()));
            if self.spec.sidecar {
                files.push(Json::from(format!("{rel}.sidecar")));
            }
        }
        let mut fields = vec![
            ("product", Json::from("hic")),
            ("pairs", Json::from(totals[0])),
            ("duplicates", Json::from(totals[1])),
            ("contacts", Json::from(totals[2])),
            ("cis", Json::from(totals[3])),
            ("trans", Json::from(totals[4])),
            ("cis_1kb+", Json::from(totals[5])),
            ("cis_10kb+", Json::from(totals[6])),
            ("cis_20kb+", Json::from(totals[7])),
            ("hic_pixels", pixels),
            (
                "mates_counted_with_their_pair",
                Json::from(status.not_leader),
            ),
            ("unassigned_mate_unplaced", Json::from(status.mate_unplaced)),
            ("unassigned_mate_mapq", Json::from(status.mate_mapq)),
            ("unassigned_no_mate_mapq", Json::from(status.no_mate_mapq)),
            ("unassigned_no_mate_cigar", Json::from(status.no_mate_cigar)),
            ("unassigned_too_close", Json::from(status.too_close)),
            ("unassigned_self_ligation", Json::from(status.self_ligation)),
            ("unassigned_no_barcode", Json::from(status.no_barcode)),
            ("spilled_bytes", Json::from(spilled)),
        ];
        if barcode_column {
            let rel = format!("{name}/cell_metrics.csv");
            let path = d.path(&rel)?;
            let metrics = BarcodeMetrics {
                barcodes: barcodes.clone(),
                columns: CELL_METRICS
                    .iter()
                    .enumerate()
                    .map(|(i, n)| (n.to_string(), cells.iter().map(|c| c[i]).collect()))
                    .collect(),
            };
            crate::sink::table::write_metrics(&metrics, &path, d)?;
            files.push(Json::from(rel));
            let per_cell = barcodes
                .names
                .iter()
                .zip(&cells)
                .map(|(b, c)| {
                    let mut f = vec![(
                        "barcode",
                        Json::from(String::from_utf8_lossy(b).into_owned()),
                    )];
                    f.extend(
                        CELL_METRICS
                            .iter()
                            .zip(c)
                            .map(|(n, &v)| (*n, Json::from(v))),
                    );
                    Json::obj(f)
                })
                .collect();
            fields.push(("cells", Json::Arr(per_cell)));
        }
        fields.extend(tally.json());
        fields.push(("files", Json::Arr(files)));
        Ok(Json::obj(fields))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mate_cigars_give_reference_lengths() {
        assert_eq!(cigar_ref_len(b"10S30M2I5D20M"), Some(55));
        assert_eq!(cigar_ref_len(b"50M100N50M"), Some(200));
        assert_eq!(cigar_ref_len(b"5H40M"), Some(40));
        assert_eq!(cigar_ref_len(b"M"), None);
        assert_eq!(cigar_ref_len(b"40"), None);
        assert_eq!(cigar_ref_len(b"*"), None);
    }
}
