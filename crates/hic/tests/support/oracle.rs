//! Bins every resolution directly from the original input, independently of builder accumulators.
use omicsweave_core::ContactUnit;
use omicsweave_hic::{BlockLayout, HicHeader, HicSidecar, HicWriter, Record};
use std::collections::BTreeMap;
use std::io::Cursor;
pub type Contact = (u32, u32, u32, u32, f32);

// Bins each resolution directly from the input and hands the unmerged
// records to the public writer, in the header's unit/resolution/pair order.
pub fn binned_file(
    header: &HicHeader,
    contacts: &[Contact],
    blocks: BlockLayout,
) -> (Vec<u8>, HicSidecar) {
    let mut writer =
        HicWriter::create_with(Cursor::new(Vec::new()), header.clone(), 1, blocks).unwrap();
    for unit in [ContactUnit::Bp, ContactUnit::Frag] {
        for &res in header.resolutions(unit) {
            let mut pairs: BTreeMap<(u32, u32), Vec<Record>> = BTreeMap::new();
            for &(a, p, b, q, value) in contacts {
                let (a, p, b, q) = if a <= b { (a, p, b, q) } else { (b, q, a, p) };
                let position = |c: u32, p: u32| {
                    let p = p.min(header.chroms[c as usize].len as u32);
                    match unit {
                        ContactUnit::Bp => p,
                        ContactUnit::Frag => {
                            header.sites[c as usize].iter().filter(|&&s| s <= p).count() as u32
                        }
                    }
                };
                let (mut x, mut y) = (position(a, p) / res, position(b, q) / res);
                if a == b && x > y {
                    std::mem::swap(&mut x, &mut y);
                }
                pairs
                    .entry((a, b))
                    .or_default()
                    .push(Record { x, y, value });
            }
            for ((a, b), records) in pairs {
                writer.write_matrix(a, b, unit, res, &records).unwrap();
            }
        }
    }
    let (out, sidecar) = writer.finish_with_sidecar().unwrap();
    (out.into_inner(), sidecar)
}
