use super::{Rec, Rng};
use omicsweave_bam::{BamReader, BinIndex};
use omicsweave_core::{
    AlignmentFields, AlignmentRead, AlignmentsOut, CigarColumn, FileSource, ModificationColumn,
    Region, StringColumn,
};
use std::path::Path;

pub fn open(bam: &Path, idx: &Path) -> BamReader<FileSource> {
    BamReader::open(
        FileSource::open(bam).unwrap(),
        BinIndex::read_file(idx).unwrap(),
    )
    .unwrap()
}

/// Queries `region` with every field and rebuilds htslib-style records.
pub fn query_all(r: &mut BamReader<FileSource>, region: Region<'_>) -> Vec<Rec> {
    let size = r.query(region, AlignmentFields::ALL).unwrap();
    let mut entries = vec![Default::default(); size.records];
    let n = size.records + 1;
    let (mut no, mut nd) = (vec![0i64; n], vec![0u8; size.name_bytes]);
    let (mut so, mut sd) = (vec![0i64; n], vec![0u8; size.seq_bytes]);
    let (mut qo, mut qd) = (vec![0i64; n], vec![0u8; size.seq_bytes]);
    let (mut co, mut cd) = (vec![0i64; n], vec![0u32; size.cigar_ops]);
    let (mut ao, mut ad) = (vec![0i64; n], vec![0u8; size.aux_bytes]);
    let (mut mo, mut md) = (vec![0i64; n], vec![Default::default(); size.modifications]);
    r.write(AlignmentsOut {
        entries: &mut entries,
        names: Some(StringColumn {
            offsets: &mut no,
            data: &mut nd,
        }),
        seqs: Some(StringColumn {
            offsets: &mut so,
            data: &mut sd,
        }),
        quals: Some(StringColumn {
            offsets: &mut qo,
            data: &mut qd,
        }),
        cigars: Some(CigarColumn {
            offsets: &mut co,
            ops: &mut cd,
        }),
        aux: Some(StringColumn {
            offsets: &mut ao,
            data: &mut ad,
        }),
        modifications: Some(ModificationColumn {
            offsets: &mut mo,
            entries: &mut md,
        }),
    })
    .unwrap();
    let tid = r.header().tid(region.chrom).unwrap() as i32;
    let span = |o: &[i64], i: usize| o[i] as usize..o[i + 1] as usize;
    entries
        .iter()
        .enumerate()
        .map(|(i, e)| Rec {
            tid,
            pos: e.start,
            end: e.end,
            flag: e.flag,
            mapq: e.mapq,
            mtid: e.mate_tid,
            mpos: e.mate_start as i64,
            isize: e.tlen as i64,
            name: nd[span(&no, i)].to_vec(),
            cigar: cd[span(&co, i)].to_vec(),
            seq: sd[span(&so, i)].to_vec(),
            qual: qd[span(&qo, i)].to_vec(),
            aux: ad[span(&ao, i)].to_vec(),
        })
        .collect()
}

pub fn regions(
    rng: &mut Rng,
    refs: &[(&'static str, u64)],
    n: usize,
) -> Vec<(&'static str, u32, u32)> {
    let mut out = Vec::new();
    for &(name, len) in refs {
        let len = len as u32;
        out.push((name, 0, len));
        out.push((name, 0, 1));
        out.push((name, len - 1, len));
        out.push((name, len, len + 1000));
        out.push((name, len - 50, len + 50));
        out.push((name, 16_383, 16_385));
    }
    for _ in 0..n {
        let &(name, len) = &refs[rng.range(0, refs.len() as u64) as usize];
        let start = rng.range(0, len) as u32;
        let width = match rng.range(0, 4) {
            0 => rng.range(1, 10),
            1 => rng.range(10, 2_000),
            2 => rng.range(2_000, 100_000),
            _ => rng.range(100_000, 1_500_000),
        } as u32;
        out.push((name, start, start.saturating_add(width)));
    }
    out
}
