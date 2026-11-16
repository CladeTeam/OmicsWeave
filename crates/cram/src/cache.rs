//! Immutable container metadata and decoded slices in the caller's cache.

use std::ops::{Deref, Range};
use std::sync::Arc;

use omicsweave_core::cache::{BlockCache, BlockKey, CacheLease, SourceToken};
use omicsweave_core::{AlignmentFields, Error, ReadAt, Result};

use crate::comp::{CompressionHeader, SliceHeader};
use crate::container::{Block, ContainerHeader, content};
use crate::header::{ReadContainer, read_container_body, read_container_header};
use crate::reader::{CramOptions, CramStats};
use crate::slice::{Decoded, SliceInput};

pub(crate) enum Loaded<T: 'static> {
    Owned(T),
    Cached(CacheLease<T>),
}

impl<T: Send + Sync> Deref for Loaded<T> {
    type Target = T;
    fn deref(&self) -> &T {
        match self {
            Self::Owned(v) => v,
            Self::Cached(v) => v,
        }
    }
}

pub(crate) fn load<T: Send + Sync, F: FnOnce() -> Result<(T, usize)>>(
    cache: &Option<Arc<BlockCache>>,
    key: BlockKey,
    f: F,
) -> Result<Loaded<T>> {
    match cache {
        Some(cache) => cache.get_or_load(key, f).map(Loaded::Cached),
        None => f().map(|(v, _)| Loaded::Owned(v)),
    }
}

pub(crate) fn meta_key(source: SourceToken, offset: u64) -> BlockKey {
    BlockKey::new(source, 6, [0, 0], offset, 0)
}

pub(crate) fn body_key(source: SourceToken, offset: u64) -> BlockKey {
    BlockKey::new(source, 6, [0, 1], offset, 0)
}

pub(crate) fn slice_key(
    source: SourceToken,
    offset: u64,
    index: usize,
    tid: i32,
    fields: AlignmentFields,
    options: &CramOptions,
) -> BlockKey {
    let variant = ((tid as u32 as u64) << 32)
        | fields.bits() as u64
        | ((options.verify_checksums as u64) << 8)
        | ((options.decode_md as u64) << 9);
    BlockKey::new(source, 6, [offset, 2], index as u64, variant)
}

pub(crate) struct Meta {
    pub header: Option<ContainerHeader>,
}

impl Meta {
    pub fn read(source: &impl ReadAt, offset: u64, head: &mut Vec<u8>) -> Result<(Self, usize)> {
        let header = read_container_header(source, offset, head)?.map(|v| v.header);
        if header.as_ref().is_some_and(|h| h.length < 0) {
            return Err(Error::invalid("CRAM container", "negative length"));
        }
        let charge = header.as_ref().map_or(0, |h| h.landmarks.capacity() * 4);
        Ok((Self { header }, charge))
    }
}

pub(crate) struct Part {
    pub header: SliceHeader,
    /// Starts after the header block, so cache hits need not inflate it.
    pub blocks: Range<usize>,
}

pub(crate) struct PreparedContainer {
    pub body: Vec<u8>,
    pub comp: CompressionHeader,
    pub parts: Vec<Part>,
}

impl PreparedContainer {
    pub fn read(
        source: &impl ReadAt,
        offset: u64,
        header: &ContainerHeader,
        scratch: &mut Vec<u8>,
        verify: bool,
        stats: &mut CramStats,
        mut body: Vec<u8>,
    ) -> Result<(Self, usize)> {
        read_container_body(
            source,
            &ReadContainer {
                offset,
                header: header.clone(),
            },
            &mut body,
        )?;
        stats.containers += 1;
        let mut at = 0;
        let comp = Block::parse(&body, &mut at)?;
        if comp.content_type != content::COMPRESSION_HEADER {
            return Err(Error::invalid("CRAM container", "no compression header"));
        }
        let comp = CompressionHeader::parse(comp.uncompress(scratch, verify)?)?;
        let mut parts = Vec::with_capacity(header.landmarks.len());
        for &landmark in &header.landmarks {
            if at as i64 != landmark as i64 {
                return Err(Error::invalid(
                    "CRAM container",
                    "slice does not match landmark",
                ));
            }
            let input = Block::parse(&body, &mut at)?;
            let decoded = input.uncompress(scratch, verify)?;
            let sh = SliceHeader::parse(input.content_type, decoded)?;
            let start = at;
            for _ in 0..sh.num_blocks {
                Block::parse(&body, &mut at)?;
            }
            parts.push(Part {
                header: sh,
                blocks: start..at,
            });
        }
        if at != body.len() {
            return Err(Error::invalid(
                "CRAM container",
                "container length does not match blocks",
            ));
        }
        let charge = body.capacity() + comp.heap_bytes() + parts.capacity() * size_of::<Part>();
        Ok((Self { body, comp, parts }, charge))
    }

    pub fn slice(&self, part: &Part) -> Result<SliceInput<'_>> {
        let mut at = part.blocks.start;
        let mut blocks = Vec::with_capacity(part.header.num_blocks as usize);
        for _ in 0..part.header.num_blocks {
            blocks.push(Block::parse(&self.body, &mut at)?);
        }
        Ok(SliceInput {
            header: part.header.clone(),
            blocks,
        })
    }
}

/// Uncached containers return only empty allocation capacity to the reader.
pub(crate) struct ContainerUse<'a> {
    pub loaded: Loaded<PreparedContainer>,
    pub spare: &'a mut Vec<u8>,
}

impl Deref for ContainerUse<'_> {
    type Target = PreparedContainer;
    fn deref(&self) -> &PreparedContainer {
        &self.loaded
    }
}

impl Drop for ContainerUse<'_> {
    fn drop(&mut self) {
        if let Loaded::Owned(container) = &mut self.loaded {
            container.body.clear();
            *self.spare = std::mem::take(&mut container.body);
        }
    }
}

/// Decoded results belong to the cache; otherwise only cleared scratch survives.
pub(crate) struct SliceUse<'a> {
    pub loaded: Loaded<PreparedSlice>,
    pub spare: &'a mut Decoded,
}

impl Deref for SliceUse<'_> {
    type Target = PreparedSlice;
    fn deref(&self) -> &PreparedSlice {
        &self.loaded
    }
}

impl Drop for SliceUse<'_> {
    fn drop(&mut self) {
        if let Loaded::Owned(slice) = &mut self.loaded {
            slice.decoded.recs.clear();
            slice.decoded.names.clear();
            slice.decoded.cigars.clear();
            slice.decoded.seqs.clear();
            slice.decoded.quals.clear();
            slice.decoded.aux.clear();
            *self.spare = std::mem::take(&mut slice.decoded);
        }
    }
}

pub(crate) struct PreparedSlice {
    pub decoded: Decoded,
    prefix_end: Vec<i64>,
}

impl PreparedSlice {
    pub fn new(decoded: Decoded, cache: &Option<Arc<BlockCache>>) -> (Self, usize) {
        let mut charge = decoded.recs.capacity() * size_of::<crate::slice::Rec>()
            + decoded.names.capacity()
            + 4 * decoded.cigars.capacity()
            + decoded.seqs.capacity()
            + decoded.quals.capacity()
            + decoded.aux.capacity();
        let eligible = cache.as_ref().is_some_and(|c| {
            let (bytes, entries) = c.limits();
            entries > 0
                && charge
                    .checked_add(decoded.recs.len() * 8)
                    .is_some_and(|n| n < bytes)
        });
        let indexed = eligible
            && decoded
                .recs
                .windows(2)
                .all(|w| w[0].ref_id == w[1].ref_id && w[0].apos <= w[1].apos);
        let mut prefix_end = Vec::new();
        if indexed {
            prefix_end.reserve_exact(decoded.recs.len());
            let mut end = i64::MIN;
            for rec in &decoded.recs {
                end = end.max(rec.aend);
                prefix_end.push(end);
            }
        }
        charge += prefix_end.capacity() * 8;
        (
            Self {
                decoded,
                prefix_end,
            },
            charge,
        )
    }

    pub fn candidates(&self, tid: i32, start: i64, end: i64) -> Range<usize> {
        if self.prefix_end.is_empty() || self.decoded.recs[0].ref_id != tid || tid == -1 {
            return 0..self.decoded.recs.len();
        }
        let lo = self.prefix_end.partition_point(|&p| p < start);
        let hi = self.decoded.recs.partition_point(|r| r.apos <= end);
        lo.min(hi)..hi
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::slice::Rec;

    #[test]
    fn evicted_slice_lease_keeps_its_owned_bytes_charged_until_release() {
        let cache = Arc::new(BlockCache::new(1 << 20, 8));
        let key = meta_key(SourceToken::new(), 100);
        let lease = cache
            .get_or_load(key, || {
                let decoded = Decoded {
                    recs: vec![Rec::default(); 100],
                    names: vec![b'A'; 4096],
                    ..Decoded::default()
                };
                Ok(PreparedSlice::new(decoded, &Some(cache.clone())))
            })
            .unwrap();
        assert!(lease.is_cached());
        let live = cache.stats().live_bytes;
        assert!(live >= 4096 + 100 * size_of::<Rec>());
        cache.clear();
        assert_eq!(cache.stats().resident_bytes, 0);
        assert_eq!(cache.stats().live_bytes, live);
        assert_eq!(lease.decoded.names, vec![b'A'; 4096]);
        drop(lease);
        assert_eq!(cache.stats().live_bytes, 0);
        assert_eq!(cache.stats().live_entries, 0);
    }

    #[test]
    fn cached_candidates_match_full_scan_and_skip_expired_prefixes() {
        for nested in [false, true] {
            let mut decoded = Decoded::default();
            for i in 0..10_000 {
                decoded.recs.push(Rec {
                    ref_id: 0,
                    apos: i * 10,
                    aend: i * 10 + if nested && i % 31 == 0 { 10_000 } else { 2 },
                    ..Rec::default()
                });
            }
            let cache = Some(Arc::new(BlockCache::new(16 << 20, 128)));
            let (prepared, _) = PreparedSlice::new(decoded, &cache);
            for start in [0, 8, 999, 50_000, 99_000, 100_000] {
                let end = start + 2;
                let want: Vec<_> = prepared
                    .decoded
                    .recs
                    .iter()
                    .enumerate()
                    .filter_map(|(i, r)| (r.apos <= end && r.aend >= start).then_some(i))
                    .collect();
                let positions = prepared.candidates(0, start, end);
                let got: Vec<_> = positions
                    .clone()
                    .filter(|&i| prepared.decoded.recs[i].aend >= start)
                    .collect();
                assert_eq!(got, want);
                if !nested {
                    assert!(positions.len() <= 2);
                }
            }
        }
    }
}
