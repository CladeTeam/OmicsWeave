use crate::allocator::counted;
use omicsweave_bbi::{
    BBI_SIDECAR_VERSION, BIGWIG_MAGIC, BbiSidecar, BigWigWriter, BlockInfo, Chrom, Error, Kind,
    Level, LevelBlocks, WriteOptions,
};
use omicsweave_core::sidecar::{self, SidecarKind, SourceIdentity};
use omicsweave_core::{ReadAt, Result};
use std::cell::Cell;
use std::io::Cursor;

fn fixture(n: u32, items: u32, zooms: Vec<u32>) -> (Vec<u8>, BbiSidecar) {
    let mut writer = BigWigWriter::create(
        Cursor::new(Vec::new()),
        vec![Chrom {
            name: "c".into(),
            len: n * 8 + 8,
        }],
        WriteOptions {
            items_per_slot: Some(items),
            block_size: 4,
            zoom_reductions: Some(zooms),
            ..WriteOptions::default()
        },
    )
    .unwrap();
    for i in 0..n {
        writer.add(0, i * 8, i * 8 + 4, (i % 97) as f32).unwrap();
    }
    let (out, sidecar) = writer.finish().unwrap();
    (out.into_inner(), sidecar)
}

fn measured_build(bytes: &[u8]) -> (omicsweave_bbi::Result<BbiSidecar>, usize) {
    let (result, largest, _) = counted(|| BbiSidecar::build(bytes));
    (result, largest)
}
#[test]
fn sidecar_small_blocks_ignore_a_large_allocation_hint_and_keep_zoom_items() {
    let (mut bytes, mut expected) = fixture(64, 4, vec![16, 128]);
    assert_eq!(BbiSidecar::build(&bytes).unwrap(), expected);
    for ceiling in [u32::MAX, 256 * 1024, 256 * 1024 + 1] {
        bytes[52..56].copy_from_slice(&ceiling.to_le_bytes());
        expected.uncompress_buf_size = ceiling;
        let (got, largest) = measured_build(&bytes);
        assert_eq!(got.unwrap(), expected);
        assert!(largest <= 256 * 1024, "requested {largest} bytes");
    }
    assert!(expected.levels.iter().all(|level| !level.blocks.is_empty()));
}

#[test]
fn sidecar_grows_to_actual_output_and_keeps_the_declared_ceiling() {
    let (bytes, expected) = fixture(12_000, 12_000, vec![]);
    assert!(expected.levels[0].blocks[0].inflated > 128 * 1024);
    let (got, largest) = measured_build(&bytes);
    assert_eq!(got.unwrap(), expected);
    assert!(largest <= 512 * 1024, "requested {largest} bytes");

    let mut too_small = bytes.clone();
    too_small[52..56].copy_from_slice(&32_768u32.to_le_bytes());
    assert!(BbiSidecar::build(&too_small).is_err());
    let block = &expected.levels[0].blocks[0];
    let mut corrupt = bytes.clone();
    corrupt[block.offset as usize] ^= 0xff;
    assert!(BbiSidecar::build(&corrupt).is_err());
    assert_eq!(BbiSidecar::build(&bytes).unwrap(), expected);
}

#[test]
fn empty_sidecar_does_not_allocate_an_inflate_buffer() {
    let (mut bytes, mut expected) = fixture(0, 4, vec![]);
    bytes[52..56].copy_from_slice(&u32::MAX.to_le_bytes());
    expected.uncompress_buf_size = u32::MAX;
    let (got, largest) = measured_build(&bytes);
    assert_eq!(got.unwrap(), expected);
    assert!(largest < 64 * 1024, "requested {largest} bytes");
    assert!(expected.level(Level::Zoom(usize::MAX)).is_none());
    assert!(expected.level(Level::Zoom(0)).is_none());
}

struct CountReads<'a> {
    bytes: &'a [u8],
    calls: Cell<usize>,
}
impl ReadAt for CountReads<'_> {
    fn len(&self) -> u64 {
        self.bytes.len() as u64
    }
    fn read_exact_at(&self, out: &mut [u8], at: u64) -> Result<()> {
        self.calls.set(self.calls.get() + 1);
        self.bytes.read_exact_at(out, at)
    }
}

#[test]
fn zoom_table_must_fit_before_allocating_or_reading_it() {
    let mut bytes = [0u8; 64];
    bytes[..4].copy_from_slice(&BIGWIG_MAGIC.to_le_bytes());
    bytes[6..8].copy_from_slice(&u16::MAX.to_le_bytes());
    let source = CountReads {
        bytes: &bytes,
        calls: Cell::new(0),
    };
    let (result, largest, _) = counted(|| BbiSidecar::build(&source));
    assert!(matches!(result, Err(Error::OutOfRange { .. })));
    assert_eq!(source.calls.get(), 1);
    assert!(largest <= 128, "largest allocation: {largest}");
}

#[test]
fn sidecar_headers_and_block_ranges_are_validated_against_available_bytes() {
    let source = SourceIdentity {
        len: 100,
        fingerprint: 1,
    };
    let mut body = vec![0, 0];
    body.extend(u16::MAX.to_le_bytes());
    body.extend(0u32.to_le_bytes());
    let sealed = sidecar::seal(SidecarKind::Bbi, BBI_SIDECAR_VERSION, source, &body);
    let (result, largest, _) = counted(|| BbiSidecar::parse(&sealed, source));
    assert!(result.is_err());
    assert!(largest < 4096, "largest allocation: {largest}");
    for (offset, size, valid) in [
        (90, 10, true),
        (100, 0, true),
        (99, 2, false),
        (101, 0, false),
        (u64::MAX, 1, false),
    ] {
        let sidecar = BbiSidecar {
            kind: Kind::BigWig,
            uncompress_buf_size: 0,
            levels: vec![LevelBlocks {
                reduction: 0,
                blocks: vec![BlockInfo {
                    start: (0, 0),
                    end: (0, 1),
                    offset,
                    size,
                    inflated: size,
                    items: 1,
                }],
            }],
        };
        let result = BbiSidecar::parse(&sidecar.to_bytes(source), source);
        if valid {
            assert_eq!(result.unwrap(), sidecar);
            assert_eq!(
                sidecar.levels[0].blocks[0].range(),
                offset..offset + size as u64
            );
        } else {
            assert!(matches!(result, Err(Error::OutOfRange { .. })));
        }
    }
}
