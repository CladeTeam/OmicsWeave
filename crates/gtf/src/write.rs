//! Writing GTF and BED12.

use std::io::Write;

use omicsweave_core::{Error, Result};

use crate::annotation::Annotation;
use crate::model::{AttributePolicy, Dialect, PHASE_NONE, Row, Strand, UNQUOTED};
use crate::strings::NONE;

/// The record order and attribute layout of a written GTF.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GtfLayout {
    /// `Original` when possible, else `Model`.
    #[default]
    Auto,
    /// The file as read: every record and `#` line in file order, every
    /// column as written, attributes in order with their quoting. Needs a
    /// GTF read with [`AttributePolicy::All`]; records synthesized at read
    /// are not written, seqids are written as renamed by read aliases.
    /// Byte-identical to the input when its attributes use the canonical
    /// layout (`key "value";` separated by one space, as GENCODE, Ensembl,
    /// RefSeq and UCSC write); blank lines and `\r` line ends are not kept.
    Original,
    /// The model: per gene, a `gene` line, then per transcript a
    /// `transcript` line and its features. Each line starts with the
    /// model's fields (`gene_id`, `transcript_id`, and on gene and
    /// transcript lines `gene_name`, `gene_type`, `transcript_name`,
    /// `transcript_type`, `tag`), followed by the other kept attributes
    /// (GFF3 `ID` / `Parent` dropped). Records outside the hierarchy are
    /// not written.
    Model,
}

const FLUSH: usize = 1 << 20;

/// Keys the model layout writes from the model, not from the attributes.
const MODEL_KEYS: [&[u8]; 10] = [
    b"gene_id",
    b"transcript_id",
    b"gene_name",
    b"gene_type",
    b"gene_biotype",
    b"transcript_name",
    b"transcript_type",
    b"transcript_biotype",
    b"ID",
    b"Parent",
];

#[inline]
fn put_u32(buf: &mut Vec<u8>, mut v: u32) {
    let mut tmp = [0u8; 10];
    let mut i = tmp.len();
    loop {
        i -= 1;
        tmp[i] = b'0' + (v % 10) as u8;
        v /= 10;
        if v == 0 {
            break;
        }
    }
    buf.extend_from_slice(&tmp[i..]);
}

fn io(e: std::io::Error) -> Error {
    e.into()
}

struct Line<'a> {
    seqid: u32,
    source: u32,
    feature_type: u32,
    start: u32,
    end: u32,
    score: u32,
    strand: Strand,
    phase: u8,
    ann: &'a Annotation,
}

impl Line<'_> {
    fn columns(&self, buf: &mut Vec<u8>) {
        let a = self.ann;
        buf.extend_from_slice(a.seqid_name(self.seqid).as_bytes());
        buf.push(b'\t');
        buf.extend_from_slice(dot(a, self.source));
        buf.push(b'\t');
        buf.extend_from_slice(dot(a, self.feature_type));
        buf.push(b'\t');
        put_u32(buf, self.start + 1);
        buf.push(b'\t');
        put_u32(buf, self.end);
        buf.push(b'\t');
        buf.extend_from_slice(dot(a, self.score));
        buf.push(b'\t');
        buf.push(self.strand.as_char() as u8);
        buf.push(b'\t');
        if self.phase == PHASE_NONE {
            buf.push(b'.');
        } else {
            buf.push(b'0' + self.phase);
        }
        buf.push(b'\t');
    }
}

fn dot(a: &Annotation, id: u32) -> &[u8] {
    if id == NONE { b"." } else { a.pool.bytes(id) }
}

fn put_attr(buf: &mut Vec<u8>, first: &mut bool, key: &[u8], value: &[u8], quoted: bool) {
    if !*first {
        buf.push(b' ');
    }
    *first = false;
    buf.extend_from_slice(key);
    buf.push(b' ');
    if quoted {
        buf.push(b'"');
        buf.extend_from_slice(value);
        buf.push(b'"');
    } else {
        buf.extend_from_slice(value);
    }
    buf.push(b';');
}

impl Annotation {
    /// Whether [`GtfLayout::Original`] can be written.
    pub fn has_original_gtf(&self) -> bool {
        self.dialect == Dialect::Gtf && self.policy == AttributePolicy::All
    }

    fn line(&self, row: Row) -> Line<'_> {
        match row {
            Row::Gene(i) => {
                let (g, i) = (&self.genes, i as usize);
                Line {
                    seqid: g.seqid[i],
                    source: g.source[i],
                    feature_type: g.feature_type[i],
                    start: g.start[i],
                    end: g.end[i],
                    score: g.score[i],
                    strand: g.strand[i],
                    phase: g.phase[i],
                    ann: self,
                }
            }
            Row::Transcript(i) => {
                let (t, i) = (&self.transcripts, i as usize);
                Line {
                    seqid: t.seqid[i],
                    source: t.source[i],
                    feature_type: t.feature_type[i],
                    start: t.start[i],
                    end: t.end[i],
                    score: t.score[i],
                    strand: t.strand[i],
                    phase: t.phase[i],
                    ann: self,
                }
            }
            Row::Feature(f) => {
                let (fs, i) = (&self.features, f as usize);
                Line {
                    seqid: self.feature_seqid(f),
                    source: fs.source[i],
                    feature_type: fs.feature_type[i],
                    start: fs.start[i],
                    end: fs.end[i],
                    score: fs.score[i],
                    strand: fs.strand[i],
                    phase: fs.phase[i],
                    ann: self,
                }
            }
            Row::Other(i) => {
                let (o, i) = (&self.others, i as usize);
                Line {
                    seqid: o.seqid[i],
                    source: o.source[i],
                    feature_type: o.feature_type[i],
                    start: o.start[i],
                    end: o.end[i],
                    score: o.score[i],
                    strand: o.strand[i],
                    phase: o.phase[i],
                    ann: self,
                }
            }
        }
    }

    /// Writes the annotation as GTF.
    pub fn write_gtf(&self, w: &mut impl Write, layout: GtfLayout) -> Result<()> {
        let original = match layout {
            GtfLayout::Auto => self.has_original_gtf(),
            GtfLayout::Original if !self.has_original_gtf() => {
                return Err(Error::Unsupported(
                    "the original GTF layout needs a GTF read with all attributes",
                ));
            }
            GtfLayout::Original => true,
            GtfLayout::Model => false,
        };
        let mut buf = Vec::with_capacity(FLUSH + 4096);
        if original {
            self.write_original(w, &mut buf)?;
        } else {
            self.write_model(w, &mut buf)?;
        }
        w.write_all(&buf).map_err(io)
    }

    fn write_original(&self, w: &mut impl Write, buf: &mut Vec<u8>) -> Result<()> {
        let mut c = 0;
        for (r, row) in self.rows().enumerate() {
            while c < self.comments.len() && self.comments[c].0 as usize <= r {
                buf.extend_from_slice(self.pool.bytes(self.comments[c].1));
                buf.push(b'\n');
                c += 1;
            }
            self.line(row).columns(buf);
            let mut first = true;
            let attrs = self.attributes.raw(r as u32);
            for i in attrs {
                let (k, v) = self.attributes.at(i);
                put_attr(
                    buf,
                    &mut first,
                    self.pool.bytes(k & !UNQUOTED),
                    self.pool.bytes(v),
                    k & UNQUOTED == 0,
                );
            }
            buf.push(b'\n');
            if buf.len() >= FLUSH {
                w.write_all(buf).map_err(io)?;
                buf.clear();
            }
        }
        for &(_, text) in &self.comments[c..] {
            buf.extend_from_slice(self.pool.bytes(text));
            buf.push(b'\n');
        }
        Ok(())
    }

    /// Writes `fields` (skipping absent ones), then the kept attributes of
    /// `record` other than those the model writes itself.
    fn model_attrs(&self, buf: &mut Vec<u8>, record: u32, fields: &[(&[u8], u32)], tags: bool) {
        let mut first = true;
        for &(k, v) in fields {
            if v != NONE {
                put_attr(buf, &mut first, k, self.pool.bytes(v), true);
            }
        }
        for i in self.attributes.raw(record) {
            let (k, v) = self.attributes.at(i);
            let key = self.pool.bytes(k & !UNQUOTED);
            if MODEL_KEYS.contains(&key) || (tags && key == b"tag") {
                continue;
            }
            put_attr(
                buf,
                &mut first,
                key,
                self.pool.bytes(v),
                k & UNQUOTED == 0 || self.dialect == Dialect::Gff3,
            );
        }
        buf.push(b'\n');
    }

    fn write_model(&self, w: &mut impl Write, buf: &mut Vec<u8>) -> Result<()> {
        let g = &self.genes;
        let t = &self.transcripts;
        let gene_t = self.pool.get(b"gene");
        let tx_t = self.pool.get(b"transcript");
        let mut fields: Vec<(&[u8], u32)> = Vec::new();
        for gene in 0..g.len() as u32 {
            let gi = gene as usize;
            let gid = g.id[gi];
            let mut line = self.line(Row::Gene(gene));
            line.feature_type = gene_t.unwrap_or(line.feature_type);
            line.columns(buf);
            self.model_attrs(
                buf,
                g.record[gi],
                &[
                    (b"gene_id", gid),
                    (b"gene_name", g.name[gi]),
                    (b"gene_type", g.biotype[gi]),
                ],
                false,
            );
            for tx in g.transcripts(gene) {
                let ti = tx as usize;
                let ids: [(&[u8], u32); 2] = [(b"gene_id", gid), (b"transcript_id", t.id[ti])];
                let mut line = self.line(Row::Transcript(tx));
                line.feature_type = tx_t.unwrap_or(line.feature_type);
                line.columns(buf);
                fields.clear();
                fields.extend_from_slice(&ids);
                fields.extend_from_slice(&[
                    (b"gene_name".as_slice(), g.name[gi]),
                    (b"gene_type", g.biotype[gi]),
                    (b"transcript_name", t.name[ti]),
                    (b"transcript_type", t.biotype[ti]),
                ]);
                fields.extend(t.tags(tx).iter().map(|&v| (b"tag".as_slice(), v)));
                self.model_attrs(buf, t.record[ti], &fields, true);
                for f in t.features(tx) {
                    self.line(Row::Feature(f)).columns(buf);
                    self.model_attrs(buf, self.features.record[f as usize], &ids, false);
                }
                if buf.len() >= FLUSH {
                    w.write_all(buf).map_err(io)?;
                    buf.clear();
                }
            }
        }
        Ok(())
    }

    /// Writes transcripts as BED12 (all, or `transcripts`): name is the
    /// transcript id, score 0; blocks are the exon blocks
    /// ([`Annotation::exon_blocks`]), the thick part spans CDS and codon
    /// rows (UCSC convention `thickStart = thickEnd = chromEnd` without
    /// them). Transcripts without exon blocks are written as one block.
    pub fn write_bed12(&self, w: &mut impl Write, transcripts: Option<&[u32]>) -> Result<()> {
        let t = &self.transcripts;
        let f = &self.features;
        let mut buf = Vec::with_capacity(FLUSH + 4096);
        let mut blocks = Vec::new();
        let all: Vec<u32>;
        let list = match transcripts {
            Some(l) => l,
            None => {
                all = (0..t.len() as u32).collect();
                &all
            }
        };
        for &tx in list {
            let ti = tx as usize;
            self.exon_blocks(tx, &mut blocks);
            if blocks.is_empty() {
                blocks.push((t.start[ti], t.end[ti]));
            }
            let (cs, ce) = (blocks[0].0, blocks[blocks.len() - 1].1);
            let mut thick: Option<(u32, u32)> = None;
            for i in t.features(tx) {
                let i = i as usize;
                if f.kind[i].is_coding() {
                    thick = Some(match thick {
                        Some((s, e)) => (s.min(f.start[i]), e.max(f.end[i])),
                        None => (f.start[i], f.end[i]),
                    });
                }
            }
            let (ts, te) = thick.map_or((ce, ce), |(s, e)| (s.max(cs), e.min(ce)));
            buf.extend_from_slice(self.seqid_name(t.seqid[ti]).as_bytes());
            buf.push(b'\t');
            put_u32(&mut buf, cs);
            buf.push(b'\t');
            put_u32(&mut buf, ce);
            buf.push(b'\t');
            buf.extend_from_slice(dot(self, t.id[ti]));
            buf.extend_from_slice(b"\t0\t");
            buf.push(match t.strand[ti] {
                Strand::Forward => b'+',
                Strand::Reverse => b'-',
                _ => b'.',
            });
            buf.push(b'\t');
            put_u32(&mut buf, ts);
            buf.push(b'\t');
            put_u32(&mut buf, te);
            buf.extend_from_slice(b"\t0\t");
            put_u32(&mut buf, blocks.len() as u32);
            buf.push(b'\t');
            for &(s, e) in &blocks {
                put_u32(&mut buf, e - s);
                buf.push(b',');
            }
            buf.push(b'\t');
            for &(s, _) in &blocks {
                put_u32(&mut buf, s - cs);
                buf.push(b',');
            }
            buf.push(b'\n');
            if buf.len() >= FLUSH {
                w.write_all(&buf).map_err(io)?;
                buf.clear();
            }
        }
        w.write_all(&buf).map_err(io)
    }
}
