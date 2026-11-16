use std::io::Cursor;

use omicsweave_bbi::{BedSchema, BigBedWriter, Chrom, WriteOptions};
use omicsweave_core::{Region, Strand};
use omicsweave_msa::MsaReader;

pub const SQL: &str = "table bigMaf\n\"alignment\"\n(\nstring chrom;\nuint chromStart;\nuint chromEnd;\nlstring mafBlock;\n)\n";

pub fn bigmaf(records: &[(u32, u32, String)], sql: Option<&str>, fields: u16) -> Vec<u8> {
    let mut writer = BigBedWriter::create(
        Cursor::new(Vec::new()),
        vec![Chrom {
            name: "c".into(),
            len: 100,
        }],
        BedSchema {
            field_count: fields,
            defined_field_count: 3,
            auto_sql: sql.map(str::to_owned),
        },
        WriteOptions {
            items_per_slot: Some(2),
            zoom_levels: 0,
            ..Default::default()
        },
    )
    .unwrap();
    for (start, end, block) in records {
        writer.add(0, *start, *end, block.as_bytes()).unwrap();
    }
    writer.finish().unwrap().0.into_inner()
}

pub fn reader(
    records: &[(u32, u32, String)],
    species: Option<Vec<String>>,
    fill: u8,
) -> MsaReader<Vec<u8>> {
    MsaReader::open(bigmaf(records, Some(SQL), 4), species, fill).unwrap()
}

pub fn fetch<S: omicsweave_core::ReadAt>(
    r: &mut MsaReader<S>,
    start: u32,
    end: u32,
    strand: Strand,
) -> Vec<u8> {
    let size = r.query(Region::new("c", start, end)).unwrap();
    let mut out = vec![0; size.elements().unwrap()];
    r.write(strand, &mut out).unwrap();
    out
}
