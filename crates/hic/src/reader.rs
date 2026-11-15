//! [`HicReader`]: standard [`ContactRead`] access to a `.hic` file over a
//! caller-provided byte source.

use std::collections::HashMap;
use std::ops::Range;
use std::sync::Arc;

use omicsweave_core::cache::{BlockCache, BlockKey, CacheLease, SourceToken};
use omicsweave_core::{
    ContactEntry, ContactQuery, ContactRead, ContactUnit, Error, ReadAt, Result, check_output,
};

use crate::block::{Inflater, decode_block};
use crate::cache::{HicBlockRequest, HicCacheLimits, PreparedBlock};
use crate::cursor::Cursor;
use crate::matrix::{BlockEntry, ZoomData};
use crate::meta::{Depth, Footer, Header};

type ZoomKey = (u32, u32, ContactUnit, u32);
const CACHE_BYPASS: &str = "hic cache admission bypass";

fn cache_limit_error(error: &Error) -> bool {
    matches!(
        error,
        Error::Unsupported(CACHE_BYPASS) | Error::Inflate("block output limit exceeded")
    ) || matches!(error, Error::Invalid { what: "hic cache", detail } if detail == "block record limit exceeded" || detail == "cached block exceeds current limits")
}

/// A `.hic` file (v8 or v9). Opening reads the header and the footer's
/// master index; each chromosome pair's block index and each fragment site
/// list is read the first time a query needs it and kept. A query reads and
/// decodes each block it overlaps once. Values are the stored contact
/// counts; stored normalization and expected vectors are not read.
pub struct HicReader<S> {
    source: S,
    header: Arc<Header>,
    footer: Arc<Footer>,
    depth: Depth,
    chrom_ids: Arc<HashMap<String, u32>>,
    zooms: HashMap<ZoomKey, Option<Arc<ZoomData>>>,
    sites: HashMap<u32, Arc<[u32]>>,
    inflater: Option<Inflater>,
    raw: Vec<u8>,
    block: Vec<u8>,
    blocks: Vec<BlockEntry>,
    entries: Vec<ContactEntry>,
    entries_ready: bool,
    cache: Option<Arc<BlockCache>>,
    cache_namespace: SourceToken,
    cache_limits: HicCacheLimits,
    strict_cache_limits: bool,
    selected: Vec<usize>,
    cache_requests: Vec<HicBlockRequest>,
    cached_records_visited: usize,
    cache_preparations: u64,
    block_inflations: u64,
}

/// A query resolved against the file: which stored matrix to read, which
/// bins of its x (chr1) and y (chr2) axes, and how cells map onto the tile.
struct Plan {
    chr1: u32,
    chr2: u32,
    x: Range<u32>,
    y: Range<u32>,
    /// The tile's rows run along y (rows lie on the higher-index chromosome).
    transposed: bool,
}

fn missing(what: &str, detail: String) -> Error {
    Error::invalid("hic query", format!("no {what} {detail}"))
}

impl<S: ReadAt> HicReader<S> {
    pub fn open(source: S) -> Result<HicReader<S>> {
        let mut cursor = Cursor::new(&source, 0);
        let header = Header::parse(&mut cursor)?;
        let footer = Footer::parse(&mut cursor, &header, source.len())?;
        let mut chrom_ids = HashMap::new();
        for (i, chrom) in header.chroms.iter().enumerate() {
            if chrom_ids.insert(chrom.name.clone(), i as u32).is_some() {
                return Err(Error::invalid("hic header", "duplicate chromosome name"));
            }
        }
        Ok(HicReader {
            source,
            depth: header.depth(),
            header: Arc::new(header),
            footer: Arc::new(footer),
            chrom_ids: Arc::new(chrom_ids),
            zooms: HashMap::new(),
            sites: HashMap::new(),
            inflater: Some(Inflater::new()),
            raw: Vec::new(),
            block: Vec::new(),
            blocks: Vec::new(),
            entries: Vec::new(),
            entries_ready: false,
            cache: None,
            cache_namespace: SourceToken::new(),
            cache_limits: HicCacheLimits::default(),
            strict_cache_limits: false,
            selected: Vec::new(),
            cache_requests: Vec::new(),
            cached_records_visited: 0,
            cache_preparations: 0,
            block_inflations: 0,
        })
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn footer(&self) -> &Footer {
        &self.footer
    }

    /// Attaches externally managed block storage. The source and its indexes
    /// must remain unchanged until this source session is closed.
    pub fn with_cache(mut self, cache: Arc<BlockCache>) -> Self {
        self.cache = Some(cache);
        self
    }

    pub fn set_cache(&mut self, cache: Option<Arc<BlockCache>>) {
        self.cache = cache;
    }

    pub fn set_cache_limits(&mut self, limits: HicCacheLimits) -> Result<()> {
        limits.check()?;
        self.cache_limits = limits;
        self.strict_cache_limits = true;
        Ok(())
    }

    /// Sets preparation/admission thresholds. Larger valid blocks retain the
    /// ordinary streaming decoder instead of being rejected as query inputs.
    pub fn set_cache_admission_limits(&mut self, limits: HicCacheLimits) -> Result<()> {
        limits.check()?;
        self.cache_limits = limits;
        self.strict_cache_limits = false;
        Ok(())
    }

    pub fn cache_namespace(&self) -> SourceToken {
        self.cache_namespace
    }

    fn context_with_source(&self, source: S) -> Self {
        Self {
            source,
            header: Arc::clone(&self.header),
            footer: Arc::clone(&self.footer),
            depth: self.depth,
            chrom_ids: Arc::clone(&self.chrom_ids),
            zooms: self.zooms.clone(),
            sites: self.sites.clone(),
            inflater: Some(Inflater::new()),
            raw: Vec::new(),
            block: Vec::new(),
            blocks: Vec::new(),
            entries: Vec::new(),
            entries_ready: false,
            cache: self.cache.clone(),
            cache_namespace: self.cache_namespace,
            cache_limits: self.cache_limits,
            strict_cache_limits: self.strict_cache_limits,
            selected: Vec::new(),
            cache_requests: Vec::new(),
            cached_records_visited: 0,
            cache_preparations: 0,
            block_inflations: 0,
        }
    }

    /// Shares the immutable source session/indexes/cache while starting an
    /// independent decoder and query scratch. The source clone must refer to
    /// the same unchanged content; `Arc<S>` makes that sharing explicit.
    pub fn clone_context(&self) -> Self
    where
        S: Clone,
    {
        self.context_with_source(self.source.clone())
    }

    /// Records selected from prepared blocks in the most recent cached query.
    /// This excludes unrelated records skipped by the row/column index.
    pub fn cached_records_visited(&self) -> usize {
        self.cached_records_visited
    }

    /// Cumulative block-index preparations performed by this context.
    pub fn cache_preparations(&self) -> u64 {
        self.cache_preparations
    }

    /// Cumulative block inflation requests, including uncached queries.
    pub fn block_inflations(&self) -> u64 {
        self.block_inflations
    }

    fn cache_key(&self, query: &ContactQuery<'_>, plan: &Plan, block: &BlockEntry) -> BlockKey {
        BlockKey::new(
            self.cache_namespace,
            3,
            [
                (plan.chr1 as u64) << 32 | plan.chr2 as u64,
                (query.resolution as u64) << 1 | u64::from(query.unit == ContactUnit::Frag),
            ],
            block.number as u64,
            (1u64 << 32) | self.header.version as u64,
        )
    }

    /// Describes independently reusable blocks, without loading payload data.
    pub fn cache_plan(
        &mut self,
        query: ContactQuery<'_>,
        requests: &mut Vec<HicBlockRequest>,
    ) -> Result<()> {
        requests.clear();
        let plan = self.prepare(&query)?;
        let mut blocks = std::mem::take(&mut self.blocks);
        if let Some(zoom) = &self.zooms[&(plan.chr1, plan.chr2, query.unit, query.resolution)] {
            zoom.blocks_for_into(
                self.header.version,
                self.depth,
                plan.x.clone(),
                plan.y.clone(),
                &mut blocks,
            );
        } else {
            blocks.clear();
        }
        requests.extend(blocks.iter().map(|entry| HicBlockRequest {
            key: self.cache_key(&query, &plan, entry),
            entry: *entry,
        }));
        self.blocks = blocks;
        Ok(())
    }

    pub fn describe_cache(
        &mut self,
        query: ContactQuery<'_>,
        keys: &mut Vec<BlockKey>,
    ) -> Result<()> {
        let mut requests = Vec::new();
        self.cache_plan(query, &mut requests)?;
        keys.clear();
        keys.extend(requests.iter().map(HicBlockRequest::key));
        Ok(())
    }

    /// Loads one planned dependency into the external cache. Callers may use
    /// this to control synchronous prefetch order without materializing tiles.
    pub fn load_cached_block(
        &mut self,
        request: &HicBlockRequest,
    ) -> Result<Option<CacheLease<PreparedBlock>>> {
        match self.load_cached_range(request, request.entry.range(), &mut None, &mut None) {
            Ok(block) => Ok(Some(block)),
            Err(Error::Unsupported(CACHE_BYPASS)) => Ok(None),
            Err(error) => Err(error),
        }
    }

    fn load_cached_range(
        &mut self,
        request: &HicBlockRequest,
        read_range: Range<u64>,
        loaded: &mut Option<Range<u64>>,
        inflated: &mut Option<BlockKey>,
    ) -> Result<CacheLease<PreparedBlock>> {
        if request.key.source != self.cache_namespace {
            return Err(Error::invalid(
                "hic cache",
                "request belongs to another source session",
            ));
        }
        let cache = self
            .cache
            .clone()
            .ok_or_else(|| Error::invalid("hic cache", "no external cache attached"))?;
        let limits = self.cache_limits;
        let (max_bytes, max_entries) = cache.limits();
        if max_bytes == 0 || max_entries == 0 {
            return Err(Error::Unsupported(CACHE_BYPASS));
        }
        let mut inflater = self.inflater.take().unwrap_or_default();
        let mut raw = std::mem::take(&mut self.raw);
        let mut block = std::mem::take(&mut self.block);
        let result = cache.get_or_load(request.key, || {
            let entry = request.entry;
            if entry
                .position
                .checked_add(entry.size as u64)
                .is_none_or(|end| end > self.source.len())
            {
                return Err(Error::invalid("hic block", "file range outside source"));
            }
            if !loaded.as_ref().is_some_and(|range| {
                range.start <= entry.position && range.end >= entry.range().end
            }) {
                if read_range.start > entry.position
                    || read_range.end < entry.range().end
                    || read_range.end > self.source.len()
                {
                    return Err(Error::invalid("hic cache", "invalid block read batch"));
                }
                let size = usize::try_from(read_range.end - read_range.start)
                    .map_err(|_| Error::invalid("hic cache", "read batch exceeds usize"))?;
                raw.resize(size, 0);
                self.source.read_exact_at(&mut raw, read_range.start)?;
                *loaded = Some(read_range.clone());
            }
            let at =
                (entry.position - loaded.as_ref().expect("read batch populated").start) as usize;
            *inflated = None;
            self.block_inflations = self.block_inflations.saturating_add(1);
            if let Err(error) = inflater.inflate_limited(
                &raw[at..at + entry.size as usize],
                &mut block,
                limits.max_inflated_bytes,
            ) {
                if !self.strict_cache_limits
                    && error == Error::Inflate("block output limit exceeded")
                {
                    return Err(Error::Unsupported(CACHE_BYPASS));
                }
                return Err(error);
            }
            *inflated = Some(request.key);
            if PreparedBlock::estimated_heap_bytes(self.header.version, &block)
                .and_then(|heap| heap.checked_add(512))
                .is_none_or(|charge| charge > max_bytes)
            {
                return Err(Error::Unsupported(CACHE_BYPASS));
            }
            self.cache_preparations = self.cache_preparations.saturating_add(1);
            let Some(prepared) =
                PreparedBlock::decode(self.header.version, &block, limits.max_records)?
            else {
                return Err(if self.strict_cache_limits {
                    Error::invalid("hic cache", "block record limit exceeded")
                } else {
                    Error::Unsupported(CACHE_BYPASS)
                });
            };
            let size = prepared.allocated_bytes();
            Ok((prepared, size))
        });
        self.inflater = Some(inflater);
        self.raw = raw;
        self.block = block;
        let prepared = result?;
        if self.strict_cache_limits && !prepared.within(limits) {
            return Err(Error::invalid(
                "hic cache",
                "cached block exceeds current limits",
            ));
        }
        Ok(prepared)
    }

    fn visit_cached(
        &mut self,
        query: &ContactQuery<'_>,
        plan: &Plan,
        mut emit: impl FnMut(u32, u32, f32),
    ) -> Result<()> {
        let mut requests = std::mem::take(&mut self.cache_requests);
        self.cache_plan(*query, &mut requests)?;
        let mut selected = std::mem::take(&mut self.selected);
        let mut loaded: Option<Range<u64>> = None;
        let mut inflated = None;
        let mut batch = 0..0;
        let mut batch_end = 0;
        self.cached_records_visited = 0;
        let result = (|| {
            for (i, request) in requests.iter().enumerate() {
                if i == batch_end {
                    let start = request.entry.position;
                    let mut end = start
                        .checked_add(request.entry.size as u64)
                        .filter(|&end| end <= self.source.len())
                        .ok_or_else(|| Error::invalid("hic block", "file range outside source"))?;
                    batch_end = i + 1;
                    while batch_end < requests.len() {
                        let next = requests[batch_end].entry;
                        if next.position != end {
                            break;
                        }
                        let next_end = end
                            .checked_add(next.size as u64)
                            .filter(|&end| end <= self.source.len())
                            .ok_or_else(|| {
                                Error::invalid("hic block", "file range outside source")
                            })?;
                        if next_end - start > 1 << 20 {
                            break;
                        }
                        end = next_end;
                        batch_end += 1;
                    }
                    batch = start..end;
                }
                let block = match self.load_cached_range(
                    request,
                    batch.clone(),
                    &mut loaded,
                    &mut inflated,
                ) {
                    Ok(block) => block,
                    Err(error) if cache_limit_error(&error) => {
                        self.visit_uncached_dependency(
                            request,
                            &mut loaded,
                            &mut inflated,
                            plan,
                            &mut emit,
                        )?;
                        continue;
                    }
                    Err(error) => return Err(error),
                };
                if self.strict_cache_limits && !block.within(self.cache_limits) {
                    return Err(Error::invalid(
                        "hic cache",
                        "cached block exceeds current limits",
                    ));
                }
                self.cached_records_visited =
                    self.cached_records_visited.saturating_add(block.visit(
                        &plan.x,
                        &plan.y,
                        plan.chr1 == plan.chr2,
                        &mut selected,
                        |x, y, value| {
                            if plan.transposed {
                                emit(y, x, value);
                            } else {
                                emit(x, y, value);
                            }
                        },
                    ));
            }
            Ok(())
        })();
        self.selected = selected;
        self.cache_requests = requests;
        result
    }

    fn visit_uncached_dependency(
        &mut self,
        request: &HicBlockRequest,
        loaded: &mut Option<Range<u64>>,
        inflated: &mut Option<BlockKey>,
        plan: &Plan,
        mut emit: impl FnMut(u32, u32, f32),
    ) -> Result<()> {
        let entry = request.entry;
        let end = entry
            .position
            .checked_add(entry.size as u64)
            .filter(|&end| end <= self.source.len())
            .ok_or_else(|| Error::invalid("hic block", "file range outside source"))?;
        if !loaded
            .as_ref()
            .is_some_and(|range| range.start <= entry.position && range.end >= end)
        {
            self.raw.resize(entry.size as usize, 0);
            self.source.read_exact_at(&mut self.raw, entry.position)?;
            *loaded = Some(entry.position..end);
        }
        if *inflated != Some(request.key) {
            let at =
                (entry.position - loaded.as_ref().expect("read batch populated").start) as usize;
            let inflater = self
                .inflater
                .as_mut()
                .expect("decoder restored after cache load");
            self.block_inflations = self.block_inflations.saturating_add(1);
            if self.strict_cache_limits {
                inflater.inflate_limited(
                    &self.raw[at..at + entry.size as usize],
                    &mut self.block,
                    self.cache_limits.max_inflated_bytes,
                )?;
            } else {
                inflater.inflate(&self.raw[at..at + entry.size as usize], &mut self.block)?;
            }
            *inflated = Some(request.key);
        } else if self.strict_cache_limits
            && self.block.len() > self.cache_limits.max_inflated_bytes
        {
            return Err(Error::Inflate("block output limit exceeded"));
        }
        let mut count = 0;
        let mut overflow = false;
        decode_block(self.header.version, &self.block, |x, y, value| {
            let (Ok(x), Ok(y)) = (u32::try_from(x), u32::try_from(y)) else {
                return;
            };
            if !value.is_finite() || value == 0.0 {
                return;
            }
            count += 1;
            if self.strict_cache_limits && count > self.cache_limits.max_records {
                overflow = true;
                return;
            }
            let mut inside = |x, y| {
                if plan.x.contains(&x) && plan.y.contains(&y) {
                    if plan.transposed {
                        emit(y, x, value);
                    } else {
                        emit(x, y, value);
                    }
                }
            };
            inside(x, y);
            if plan.chr1 == plan.chr2 && x != y {
                inside(y, x);
            }
        })?;
        if overflow {
            return Err(Error::invalid("hic cache", "block record limit exceeded"));
        }
        self.cached_records_visited += count;
        Ok(())
    }

    /// Index of a chromosome by name.
    pub fn chrom_index(&self, name: &str) -> Result<u32> {
        self.chrom_ids
            .get(name)
            .copied()
            .ok_or_else(|| Error::UnknownChrom(name.to_owned()))
    }

    /// A chromosome's fragment count (its number of restriction sites).
    pub fn fragment_count(&self, chrom: &str) -> Result<u32> {
        let c = self.chrom_index(chrom)?;
        self.header
            .sites
            .get(c as usize)
            .map(|s| s.count)
            .ok_or_else(|| missing("fragment sites", format!("for {chrom}")))
    }

    /// A chromosome's restriction site positions, read once.
    pub fn sites(&mut self, chrom: u32) -> Result<&[u32]> {
        if !self.sites.contains_key(&chrom) {
            let list = *self
                .header
                .sites
                .get(chrom as usize)
                .ok_or_else(|| missing("fragment sites", format!("for chromosome {chrom}")))?;
            let buf = self.read_range(list.position, list.count as u64 * 4, "fragment sites")?;
            let sites = buf
                .chunks_exact(4)
                .map(|b| {
                    u32::try_from(i32::from_le_bytes(b.try_into().unwrap()))
                        .map_err(|_| Error::invalid("fragment sites", "negative position"))
                })
                .collect::<Result<Vec<_>>>()?;
            if sites.windows(2).any(|w| w[0] > w[1]) {
                return Err(Error::invalid("fragment sites", "unsorted positions"));
            }
            self.sites.insert(chrom, Arc::from(sites));
        }
        Ok(&self.sites[&chrom])
    }

    /// The block index of a chromosome pair (`chr1 <= chr2`) at `bin_size`
    /// in `unit`, read once; `None` when the file has no such map.
    pub fn zoom(
        &mut self,
        chr1: u32,
        chr2: u32,
        unit: ContactUnit,
        bin_size: u32,
    ) -> Result<Option<&ZoomData>> {
        let key = (chr1, chr2, unit, bin_size);
        if !self.zooms.contains_key(&key) {
            let zoom = match self.footer.matrices.get(&(chr1, chr2)) {
                Some(&position) => ZoomData::parse(
                    &mut Cursor::new(&self.source, position),
                    position,
                    unit,
                    bin_size,
                )?,
                None => None,
            };
            if zoom
                .as_ref()
                .is_some_and(|z| (z.chr1, z.chr2) != (chr1, chr2))
            {
                return Err(Error::invalid(
                    "hic matrix",
                    "chromosome pair disagrees with master index",
                ));
            }
            self.zooms.insert(key, zoom.map(Arc::new));
        }
        Ok(self.zooms[&key].as_deref())
    }

    fn read_range(&self, position: u64, nbytes: u64, what: &'static str) -> Result<Vec<u8>> {
        if position
            .checked_add(nbytes)
            .is_none_or(|end| end > self.source.len())
        {
            return Err(Error::invalid(what, "range outside source"));
        }
        let n = usize::try_from(nbytes).map_err(|_| Error::invalid(what, "range exceeds usize"))?;
        let mut buf = vec![0u8; n];
        self.source.read_exact_at(&mut buf, position)?;
        Ok(buf)
    }

    /// Resolves the query and reads its block index.
    fn prepare(&mut self, query: &ContactQuery<'_>) -> Result<Plan> {
        query.check()?;
        let (unit, res) = (query.unit, query.resolution);
        let (r, c) = (
            self.chrom_index(query.rows.chrom)?,
            self.chrom_index(query.cols.chrom)?,
        );
        let plan = if r <= c {
            Plan {
                chr1: r,
                chr2: c,
                x: query.row_bins(),
                y: query.col_bins(),
                transposed: false,
            }
        } else {
            Plan {
                chr1: c,
                chr2: r,
                x: query.col_bins(),
                y: query.row_bins(),
                transposed: true,
            }
        };
        let stored = self.zoom(plan.chr1, plan.chr2, unit, res)?.is_some();
        if !stored && !self.header.resolutions(unit).contains(&res) {
            return Err(missing("resolution", format!("{res} {}", unit.as_str())));
        }
        Ok(plan)
    }

    /// Calls `f(row_bin, col_bin, count)` for every stored cell of the tile
    /// with a finite, nonzero count, decoding each overlapping block once.
    fn visit(
        &mut self,
        query: &ContactQuery<'_>,
        plan: &Plan,
        mut f: impl FnMut(u32, u32, f32),
    ) -> Result<()> {
        self.cached_records_visited = 0;
        if plan.x.is_empty() || plan.y.is_empty() {
            return Ok(());
        }
        if self.cache.as_ref().is_some_and(|cache| {
            let (bytes, entries) = cache.limits();
            bytes != 0 && entries != 0
        }) {
            return self.visit_cached(query, plan, f);
        }
        let strict_limits =
            (self.cache.is_some() && self.strict_cache_limits).then_some(self.cache_limits);
        let (unit, res) = (query.unit, query.resolution);
        let version = self.header.version;
        let mut blocks = std::mem::take(&mut self.blocks);
        if let Some(zoom) = &self.zooms[&(plan.chr1, plan.chr2, unit, res)] {
            zoom.blocks_for_into(
                version,
                self.depth,
                plan.x.clone(),
                plan.y.clone(),
                &mut blocks,
            );
        } else {
            blocks.clear();
        }
        let mut inflater = self.inflater.take().unwrap_or_default();
        let mut raw = std::mem::take(&mut self.raw);
        let mut block = std::mem::take(&mut self.block);
        let result = (|| {
            let intra = plan.chr1 == plan.chr2;
            let mut emit = |x: i32, y: i32, v: f32| {
                let (Ok(x), Ok(y)) = (u32::try_from(x), u32::try_from(y)) else {
                    return;
                };
                if !plan.x.contains(&x) || !plan.y.contains(&y) {
                    return;
                }
                if !v.is_finite() || v == 0.0 {
                    return;
                }
                if plan.transposed {
                    f(y, x, v)
                } else {
                    f(x, y, v)
                }
            };
            let mut first = 0;
            while first < blocks.len() {
                let start = blocks[first].position;
                let end_of = |entry: &BlockEntry| {
                    entry
                        .position
                        .checked_add(entry.size as u64)
                        .filter(|&end| end <= self.source.len())
                        .ok_or_else(|| Error::invalid("hic block", "file range outside source"))
                };
                let mut end = end_of(&blocks[first])?;
                let mut last = first + 1;
                while last < blocks.len() && blocks[last].position == end {
                    let next = end_of(&blocks[last])?;
                    if next - start > 1 << 20 {
                        break;
                    }
                    end = next;
                    last += 1;
                }
                let len = usize::try_from(end - start)
                    .map_err(|_| Error::invalid("hic block", "read exceeds usize"))?;
                raw.resize(len, 0);
                self.source.read_exact_at(&mut raw, start)?;
                for entry in &blocks[first..last] {
                    let lo = (entry.position - start) as usize;
                    self.block_inflations = self.block_inflations.saturating_add(1);
                    if let Some(limits) = strict_limits {
                        inflater.inflate_limited(
                            &raw[lo..lo + entry.size as usize],
                            &mut block,
                            limits.max_inflated_bytes,
                        )?;
                        let mut records = 0usize;
                        decode_block(version, &block, |x, y, value| {
                            if x >= 0 && y >= 0 && value.is_finite() && value != 0.0 {
                                records = records.saturating_add(1);
                                if records > limits.max_records {
                                    return;
                                }
                            }
                            emit(x, y, value);
                            if intra && x != y {
                                emit(y, x, value);
                            }
                        })?;
                        if records > limits.max_records {
                            return Err(Error::invalid("hic cache", "block record limit exceeded"));
                        }
                        continue;
                    }
                    inflater.inflate(&raw[lo..lo + entry.size as usize], &mut block)?;
                    decode_block(version, &block, |x, y, v| {
                        emit(x, y, v);
                        if intra && x != y {
                            emit(y, x, v);
                        }
                    })?;
                }
                first = last;
            }
            Ok(())
        })();
        self.blocks = blocks;
        self.inflater = Some(inflater);
        self.raw = raw;
        self.block = block;
        result
    }
}

impl HicReader<omicsweave_core::FileSource> {
    /// Duplicates the open file handle, sharing session identity and metadata.
    pub fn try_clone_context(&self) -> Result<Self> {
        Ok(self.context_with_source(self.source.try_clone()?))
    }
}

impl<S: ReadAt> ContactRead for HicReader<S> {
    fn chrom_len(&self, chrom: &str) -> Result<u32> {
        let len = self.header.chroms[self.chrom_index(chrom)? as usize].len;
        u32::try_from(len)
            .map_err(|_| Error::invalid("chromosome", format!("{chrom} longer than 2^32 - 1")))
    }

    fn resolutions(&self, unit: ContactUnit) -> &[u32] {
        self.header.resolutions(unit)
    }

    fn matrix(&mut self, query: ContactQuery<'_>, out: &mut [f32]) -> Result<()> {
        query.check()?;
        let (rows, cols) = query.shape();
        let cells = rows
            .checked_mul(cols)
            .ok_or_else(|| Error::invalid("contact matrix", "shape overflows usize"))?;
        check_output("contact matrix", out.len(), cells)?;
        let plan = self.prepare(&query)?;
        if cells == 0 {
            return Ok(());
        }
        let (r0, c0) = (query.row_bins().start, query.col_bins().start);
        out.fill(0.0);
        self.visit(&query, &plan, |r, c, v| {
            let at = (r - r0) as usize * cols + (c - c0) as usize;
            // SAFETY: visit clips emitted bins to the query. The checked
            // product rows*cols equals out.len(), so at is in range.
            unsafe {
                *out.get_unchecked_mut(at) = v;
            }
        })
    }

    fn query_entries(&mut self, query: ContactQuery<'_>) -> Result<usize> {
        self.entries_ready = false;
        let plan = self.prepare(&query)?;
        let mut staged = std::mem::take(&mut self.entries);
        staged.clear();
        let res = query.resolution;
        let result = self.visit(&query, &plan, |r, c, value| {
            staged.push(ContactEntry {
                row: r.saturating_mul(res),
                col: c.saturating_mul(res),
                value,
            })
        });
        self.entries = staged;
        result?;
        let staged = &mut self.entries;
        staged.sort_unstable_by_key(|e| (e.row, e.col));
        let n = staged.len();
        self.entries_ready = true;
        Ok(n)
    }

    fn write_entries(&self, out: &mut [ContactEntry]) -> Result<()> {
        if !self.entries_ready {
            return Err(Error::NotQueried);
        }
        let staged = &self.entries;
        check_output("contact entries", out.len(), staged.len())?;
        out.copy_from_slice(staged);
        Ok(())
    }
}
