//! Shared Juicer-layout inputs for Python bindings and upstream straw.
use std::io::Cursor;
use std::path::PathBuf;

use omicsweave_core::{ContactUnit, sidecar::checksum};
use omicsweave_hic::{BlockLayout, Chrom, HicHeader, HicWriter, Record};

fn header(version: i32, fragments: bool) -> HicHeader {
    let mut h = HicHeader::new(
        "round8",
        ["All", "chr1", "chr2"]
            .into_iter()
            .map(|name| Chrom {
                name: name.into(),
                len: if name == "All" { 8_000_000 } else { 4_000_000 },
            })
            .collect(),
        if fragments { vec![] } else { vec![1000] },
    );
    h.version = version;
    if fragments {
        h.frag_resolutions = vec![1, 2];
        h.sites = vec![
            vec![],
            vec![0, 17, 100, 700, 4_000_000],
            vec![0, 17, 100, 700, 4_000_000],
        ];
    }
    h
}

fn fixture(version: i32, fragments: bool) -> Vec<u8> {
    let mut writer = HicWriter::create_with(
        Cursor::new(Vec::new()),
        header(version, fragments),
        1,
        BlockLayout::juicer(),
    )
    .unwrap();
    if fragments {
        for resolution in [1, 2] {
            for (chr2, cells) in [
                (1, vec![(1, 1, 1.0), (2, 3, 3.0), (4, 4, 4.0), (5, 5, 8.0)]),
                (2, vec![(1, 2, 3.0), (3, 4, 5.5), (5, 5, 7.0)]),
            ] {
                let records: Vec<_> = cells
                    .into_iter()
                    .map(|(x, y, value)| Record {
                        x: x / resolution,
                        y: y / resolution,
                        value,
                    })
                    .collect();
                writer
                    .write_matrix(1, chr2, ContactUnit::Frag, resolution, &records)
                    .unwrap();
            }
            writer
                .write_matrix(2, 2, ContactUnit::Frag, resolution, &[])
                .unwrap();
        }
    } else {
        let value = |x, y| ((x * 3 + y * 7) % 97 + 1) as f32 / 8.0;
        let mut cis: Vec<_> = (0..4000)
            .flat_map(|x| {
                (x..(x + 48).min(4000)).map(move |y| Record {
                    x,
                    y,
                    value: value(x, y),
                })
            })
            .collect();
        cis.push(Record {
            x: 4000,
            y: 4000,
            value: 0.5,
        });
        writer
            .write_matrix(1, 1, ContactUnit::Bp, 1000, &cis)
            .unwrap();
        let mut trans: Vec<_> = (0..4000)
            .flat_map(|x| {
                (0..12).map(move |offset| {
                    let y = (x * 17 + offset * 31) % 4000;
                    Record {
                        x,
                        y,
                        value: value(x, y),
                    }
                })
            })
            .collect();
        trans.push(Record {
            x: 4000,
            y: 4000,
            value: 1.5,
        });
        writer
            .write_matrix(1, 2, ContactUnit::Bp, 1000, &trans)
            .unwrap();
        writer
            .write_matrix(2, 2, ContactUnit::Bp, 1000, &[])
            .unwrap();
    }
    writer.finish().unwrap().into_inner()
}

fn main() {
    let dir = PathBuf::from(std::env::args_os().nth(1).expect("output directory"));
    std::fs::create_dir_all(&dir).unwrap();
    for version in [8, 9] {
        for fragments in [false, true] {
            let unit = if fragments { "frag" } else { "bp" };
            let path = dir.join(format!("round8-{unit}-v{version}.hic"));
            let bytes = fixture(version, fragments);
            std::fs::write(&path, &bytes).unwrap();
            println!(
                "{}\t{}\t{:016x}",
                path.display(),
                bytes.len(),
                checksum(&bytes)
            );
        }
    }
}
