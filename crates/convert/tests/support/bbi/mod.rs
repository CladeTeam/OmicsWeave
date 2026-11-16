//! Test-side composition of the crate's steps into whole-file queries, with
//! plain `pread` I/O and owned buffers. This is what a caller of the crate
//! writes; it is kept out of the public API on purpose.

#![allow(dead_code)]

pub mod libbigwig;

use std::fs::File;
use std::os::unix::fs::FileExt;
use std::path::Path;

use omicsweave_bbi::{
    BedBlock, BlockRef, Chroms, HEADER_LEN, Header, Index, Inflater, Layout, Level, NodeInfo,
    Summary, TotalSummary, WigBlock, ZoomBlock, ZoomLevel, parse_auto_sql,
};

fn read(file: &File, range: std::ops::Range<u64>) -> Vec<u8> {
    let mut buf = vec![0u8; (range.end - range.start) as usize];
    file.read_exact_at(&mut buf, range.start).expect("read");
    buf
}

pub struct Opened {
    file: File,
    pub header: Header,
    pub zooms: Vec<ZoomLevel>,
    pub layout: Layout,
    pub chroms: Chroms,
    pub summary: Option<TotalSummary>,
    pub auto_sql: Option<String>,
    /// (file offset, bytes) of the base index and of each zoom index.
    indexes: Vec<(u64, Vec<u8>)>,
    inflater: Inflater,
    raw: Vec<u8>,
    block: Vec<u8>,
}

impl Opened {
    pub fn open(path: &Path) -> Opened {
        let file = File::open(path).expect("open");
        let file_len = file.metadata().expect("stat").len();
        let header = Header::parse(&read(&file, 0..HEADER_LEN as u64)).expect("header");
        let zooms =
            ZoomLevel::parse_all(&header, &read(&file, header.zoom_headers())).expect("zooms");
        let layout = Layout::new(&header, &zooms, file_len).expect("layout");
        let range = layout.chrom_tree();
        let chroms = Chroms::parse(&read(&file, range.clone()), range.start).expect("chroms");
        let summary = layout
            .total_summary()
            .map(|r| TotalSummary::parse(&read(&file, r)).expect("summary"));
        let auto_sql = layout
            .auto_sql()
            .map(|r| parse_auto_sql(&read(&file, r)).expect("autoSql").to_owned());
        let indexes = std::iter::once(Level::Base)
            .chain((0..zooms.len()).map(Level::Zoom))
            .map(|level| {
                let range = layout.index(level).expect("index range");
                (range.start, read(&file, range))
            })
            .collect();
        let block = vec![0u8; header.block_buffer_len()];
        Opened {
            file,
            header,
            zooms,
            layout,
            chroms,
            summary,
            auto_sql,
            indexes,
            inflater: Inflater::new(),
            raw: Vec::new(),
            block,
        }
    }

    pub fn index(&self, level: Level) -> Index<'_> {
        let i = match level {
            Level::Base => 0,
            Level::Zoom(z) => z + 1,
        };
        let (base, bytes) = &self.indexes[i];
        Index::new(bytes, *base).expect("index")
    }

    pub fn blocks(&self, level: Level, chrom: u32, start: u32, end: u32) -> Vec<BlockRef> {
        let index = self.index(level);
        let raw: Vec<_> = index
            .query(chrom, start, end)
            .collect::<Result<_, _>>()
            .expect("query");
        let mut nodes: Vec<NodeInfo> = Vec::new();
        let ready = index.prepare(&mut nodes).expect("prepare index");
        let prepared: Vec<_> = ready
            .query(chrom, start, end)
            .collect::<Result<_, _>>()
            .expect("prepared query");
        assert_eq!(raw, prepared);
        prepared
    }

    /// Reads one block and returns its decompressed bytes.
    fn load(&mut self, block: BlockRef) -> &[u8] {
        self.raw.resize(block.size as usize, 0);
        self.file
            .read_exact_at(&mut self.raw, block.offset)
            .expect("read block");
        if self.header.compressed() {
            let n = self
                .inflater
                .inflate(&self.raw, &mut self.block)
                .expect("inflate");
            &self.block[..n]
        } else {
            &self.raw
        }
    }

    /// `values` as pyBigWig returns them: NaN where uncovered.
    pub fn values(&mut self, chrom: &str, start: u32, end: u32) -> Vec<f32> {
        let mut out = vec![f32::NAN; (end - start) as usize];
        let Some(id) = self.chroms.id(chrom) else {
            return out;
        };
        for block in self.blocks(Level::Base, id, start, end) {
            WigBlock::parse(self.load(block))
                .expect("wig block")
                .prepare()
                .fill(id, start, &mut out);
        }
        out
    }

    pub fn intervals(&mut self, chrom: &str, start: u32, end: u32) -> Vec<(u32, u32, f32)> {
        let id = self.chroms.id(chrom).expect("chrom");
        let mut out = Vec::new();
        for block in self.blocks(Level::Base, id, start, end) {
            out.extend(
                WigBlock::parse(self.load(block))
                    .expect("wig block")
                    .prepare()
                    .overlapping(id, start, end),
            );
        }
        out
    }

    pub fn entries(&mut self, chrom: &str, start: u32, end: u32) -> Vec<(u32, u32, Vec<u8>)> {
        let id = self.chroms.id(chrom).expect("chrom");
        let mut out = Vec::new();
        for block in self.blocks(Level::Base, id, start, end) {
            let raw = BedBlock::parse(self.load(block)).expect("bed block");
            let mut offsets = Vec::new();
            for entry in raw
                .prepare(&mut offsets)
                .expect("prepare bed")
                .overlapping(id, start, end)
            {
                out.push((entry.start, entry.end, entry.rest.to_vec()));
            }
        }
        out
    }

    pub fn summaries(&mut self, zoom: usize, chrom: u32, start: u32, end: u32) -> Vec<Summary> {
        let mut out = Vec::new();
        for block in self.blocks(Level::Zoom(zoom), chrom, start, end) {
            out.extend(
                ZoomBlock::parse(self.load(block))
                    .expect("zoom block")
                    .prepare()
                    .overlapping(chrom, start, end),
            );
        }
        out
    }
}

/// Values equal, or both uncovered.
pub fn same(a: f32, b: f32) -> bool {
    a == b || (a.is_nan() && b.is_nan())
}
