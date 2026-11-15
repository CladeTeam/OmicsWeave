//! A whole CRAM file container by container: slices located without
//! decoding them, decoded independently into BAM records (conversion), and
//! the `.crai` htslib's `cram_index_build` writes.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use omicsweave_core::{Error, ReadAt, Result, SequenceRead};

use crate::comp::{CompressionHeader, ds};
use crate::container::{Block, ContainerHeader, content};
use crate::header::{CramHeader, read_container_body, read_container_header};
use crate::index::CraiEntry;
use crate::reader::CramOptions;
use crate::record::{BamFields, Naming};
use crate::reference::ByName;
use crate::slice::{Decode, Decoded, References, SliceInput};

/// One container: its header, compression header and the byte ranges of
/// its slices within its body.
pub struct Container {
    /// File offset of the container.
    pub offset: u64,
    pub header: ContainerHeader,
    comp: Arc<CompressionHeader>,
    body: Arc<[u8]>,
    slices: Vec<std::ops::Range<usize>>,
    methods: Vec<u8>,
}

/// The containers of a CRAM file in order, from its first data container.
pub struct Scanner<S> {
    source: S,
    offset: u64,
    head: Vec<u8>,
    verify: bool,
}

impl<S: ReadAt> Scanner<S> {
    pub fn new(source: S, header: &CramHeader, verify_checksums: bool) -> Scanner<S> {
        Scanner {
            source,
            offset: header.first_container(),
            head: Vec::new(),
            verify: verify_checksums,
        }
    }

    /// Bytes of the source consumed so far.
    pub fn position(&self) -> u64 {
        self.offset
    }

    /// The next container; `None` at the end of the file. Each slice's
    /// header and block headers are read to locate it, as
    /// `cram_index_container` reads them, and must match the container's
    /// landmarks and length.
    pub fn next_container(&mut self) -> Result<Option<Container>> {
        let Some(c) = read_container_header(&self.source, self.offset, &mut self.head)? else {
            return Ok(None);
        };
        let mut body = Vec::new();
        read_container_body(&self.source, &c, &mut body)?;
        let mut at = 0;
        let comp_block = Block::parse(&body, &mut at)?;
        if comp_block.content_type != content::COMPRESSION_HEADER {
            return Err(Error::invalid(
                "CRAM container",
                "no compression header block",
            ));
        }
        let mut scratch = Vec::new();
        let comp = CompressionHeader::parse(comp_block.uncompress(&mut scratch, self.verify)?)?;
        let mut slices = Vec::with_capacity(c.header.landmarks.len());
        let mut methods = vec![comp_block.method];
        for &landmark in &c.header.landmarks {
            if at as i64 != landmark as i64 {
                return Err(Error::invalid(
                    "CRAM container",
                    format!("slice at {at} does not match landmark {landmark}"),
                ));
            }
            let start = at;
            let slice = SliceInput::read(&body, &mut at, &mut scratch)?;
            methods.extend(slice.blocks.iter().map(|b| b.method));
            slices.push(start..at);
        }
        if at != body.len() {
            return Err(Error::invalid(
                "CRAM container",
                "container length does not match its blocks",
            ));
        }
        let offset = self.offset;
        self.offset += c.header.header_len as u64 + body.len() as u64;
        Ok(Some(Container {
            offset,
            header: c.header,
            comp: Arc::new(comp),
            body: body.into(),
            slices,
            methods,
        }))
    }
}

/// Where a BAM record goes in a BAM index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BamEntry {
    pub tid: i32,
    pub beg: i64,
    pub end: i64,
    pub mapped: bool,
}

/// The BAM records of a slice: their bytes back to back, each with its
/// `block_size`, and the length and index entry of each.
#[derive(Debug, Default)]
pub struct BamRecords {
    pub bytes: Vec<u8>,
    pub records: Vec<(u32, BamEntry)>,
}

/// The references slices are decoded against, shared by decoding
/// threads.
pub struct SharedReference<R> {
    reader: Mutex<R>,
}

impl<R: SequenceRead> SharedReference<R> {
    pub fn new(reader: R) -> SharedReference<R> {
        SharedReference {
            reader: Mutex::new(reader),
        }
    }
}

struct Locked<'a, R> {
    shared: &'a SharedReference<R>,
    header: &'a CramHeader,
}

impl<R: SequenceRead> References for Locked<'_, R> {
    fn ref_len(&mut self, tid: usize) -> Result<i64> {
        let mut r = self.shared.reader.lock().unwrap_or_else(|e| e.into_inner());
        ByName {
            reader: &mut *r,
            refs: self.header.references(),
        }
        .ref_len(tid)
    }

    fn bases(&mut self, tid: usize, start: i64, end: i64, out: &mut Vec<u8>) -> Result<()> {
        let mut r = self.shared.reader.lock().unwrap_or_else(|e| e.into_inner());
        ByName {
            reader: &mut *r,
            refs: self.header.references(),
        }
        .bases(tid, start, end, out)
    }
}

impl Container {
    /// Number of slices.
    pub fn slices(&self) -> usize {
        self.slices.len()
    }

    /// The compression method of each block of the container's
    /// compression header and slices (0 raw, 1 gzip, 2 bzip2, 3 lzma,
    /// 4 rANS 4x8, 5 rANS Nx16, 6 arithmetic, 7 fqzcomp, 8 name
    /// tokeniser).
    pub fn block_methods(&self) -> &[u8] {
        &self.methods
    }

    /// Records of the container.
    pub fn records(&self) -> i32 {
        self.header.num_records
    }

    /// Decodes slice `i` with every field into the BAM records htslib
    /// reads from it, appended to `out`.
    pub fn decode_bam<R: SequenceRead>(
        &self,
        i: usize,
        header: &CramHeader,
        reference: Option<&SharedReference<R>>,
        options: &CramOptions,
        out: &mut BamRecords,
    ) -> Result<()> {
        let mut scratch = Vec::new();
        let mut at = self.slices[i].start;
        let slice = SliceInput::read(&self.body, &mut at, &mut scratch)?;
        let decode = Decode {
            comp: &self.comp,
            ref_lens: &header.ref_lens,
            unknown_rg: header.unknown_rg,
            wanted: ds::ALL,
            bases: true,
            decode_md: options.decode_md,
            mates: true,
            holds_refs: true,
            range_ref: -2,
            verify: options.verify_checksums,
        };
        let mut locked;
        let refs: Option<&mut dyn References> = match reference {
            Some(shared) => {
                locked = Locked { shared, header };
                Some(&mut locked)
            }
            None => None,
        };
        let mut dec = Decoded::default();
        decode.run(&slice, refs, &mut dec)?;
        let naming = Naming {
            prefix: &options.name_prefix,
            read_groups: &header.read_groups,
        };
        let mut f = BamFields::default();
        for r in 0..dec.recs.len() {
            dec.bam_fields(r, &naming, true, true, &mut f)?;
            let at = out.bytes.len();
            dec.push_bam(r, &f, &mut out.bytes)?;
            out.records.push((
                (out.bytes.len() - at) as u32,
                BamEntry {
                    tid: f.tid,
                    beg: f.entry.start,
                    end: f.entry.end,
                    mapped: f.entry.flag & 0x4 == 0,
                },
            ));
        }
        Ok(())
    }

    /// The `.crai` lines of slice `i` (`cram_index_slice`).
    fn index_slice(
        &self,
        i: usize,
        header: &CramHeader,
        verify: bool,
        out: &mut Vec<CraiEntry>,
    ) -> Result<()> {
        let range = self.slices[i].clone();
        let mut scratch = Vec::new();
        let mut at = range.start;
        let slice = SliceInput::read(&self.body, &mut at, &mut scratch)?;
        let sz = (range.end - range.start) as i32;
        let landmark = self.header.landmarks[i];
        let h = &slice.header;
        if h.ref_seq_id != -2 {
            out.push(CraiEntry {
                ref_id: h.ref_seq_id,
                start: h.ref_seq_start,
                span: h.ref_seq_span,
                container: self.offset as i64,
                landmark,
                slice_len: sz,
            });
            return Ok(());
        }
        let decode = Decode {
            comp: &self.comp,
            ref_lens: &header.ref_lens,
            unknown_rg: header.unknown_rg,
            wanted: ds::RI | ds::BF | ds::AP | ds::CIGAR,
            bases: false,
            decode_md: false,
            mates: false,
            holds_refs: false,
            range_ref: -2,
            verify,
        };
        let mut dec = Decoded::default();
        decode.run(&slice, None, &mut dec)?;
        let entry = |r: i32, start: i64, end: i64| CraiEntry {
            ref_id: r,
            start,
            span: end - start + 1,
            container: self.offset as i64,
            landmark,
            slice_len: sz,
        };
        let (mut cur, mut start, mut end) = (-2i32, 0i64, i32::MIN as i64);
        let (mut last_ref, mut last_pos) = (-9i32, -9i64);
        for cr in &dec.recs {
            if cr.ref_id == last_ref && cr.apos < last_pos {
                return Err(Error::invalid(
                    "CRAM index",
                    "the file is not sorted by reference and position",
                ));
            }
            last_ref = cr.ref_id;
            last_pos = cr.apos as i32 as i64;
            if cr.ref_id == cur {
                end = end.max(cr.aend);
                continue;
            }
            if cur != -2 {
                out.push(entry(cur, start, end));
            }
            cur = cr.ref_id;
            start = cr.apos;
            end = cr.aend;
        }
        if cur != -2 {
            out.push(entry(cur, start, end));
        }
        Ok(())
    }
}

/// The `.crai` entries `samtools index` writes for the CRAM `source`
/// (`cram_index_build`). Only multi-reference slices are decoded, to find
/// where each reference's records lie.
pub fn index_entries<S: ReadAt>(source: S) -> Result<Vec<CraiEntry>> {
    let header = CramHeader::read(&source)?;
    let mut scan = Scanner::new(source, &header, true);
    let mut out = Vec::new();
    let (mut last_ref, mut last_start) = (-9i64, -9i64);
    while let Some(c) = scan.next_container()? {
        let h = &c.header;
        if h.ref_seq_id as i64 == last_ref && h.ref_seq_start < last_start {
            return Err(Error::invalid(
                "CRAM index",
                "the file is not sorted by reference and position",
            ));
        }
        last_ref = h.ref_seq_id as i64;
        last_start = h.ref_seq_start;
        for i in 0..c.slices() {
            c.index_slice(i, &header, true, &mut out)?;
        }
    }
    Ok(out)
}

/// The `.crai` of the CRAM `source` as `samtools index` writes it: its
/// lines, gzip-compressed.
pub fn build_index<S: ReadAt>(source: S) -> Result<Vec<u8>> {
    let entries = index_entries(source)?;
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    for e in &entries {
        gz.write_all(e.line().as_bytes())?;
    }
    Ok(gz.finish()?)
}

/// [`build_index`] of the CRAM file at `cram`, saved beside it as
/// `<cram>.crai` (replacing one there). Returns its path.
pub fn build_index_file(cram: &Path) -> Result<PathBuf> {
    let source = omicsweave_core::FileSource::open(cram)?;
    let bytes = build_index(source)?;
    let mut path = cram.as_os_str().to_owned();
    path.push(".crai");
    let path = PathBuf::from(path);
    omicsweave_core::sidecar::write_file(&path, &bytes)?;
    Ok(path)
}
