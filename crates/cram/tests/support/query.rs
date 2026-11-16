//! Allocation and record comparison for region-query assertions.

use super::{Fixture, Hts, Rec, Rng, fasta_reader};
use omicsweave_core::{
    AlignmentFields, AlignmentRead, AlignmentsOut, CigarColumn, FileSource, ModificationColumn,
    Region, StringColumn,
};
use omicsweave_cram::{CramIndex, CramOptions, CramReader};
use omicsweave_fasta::FastaReader;
use std::path::Path;

pub type Reader = CramReader<FileSource, FastaReader<FileSource>>;

pub fn open(cram: &Path, crai: &Path, fasta: Option<&Path>) -> Reader {
    let header = omicsweave_cram::CramHeader::read(&FileSource::open(cram).unwrap()).unwrap();
    let index = CramIndex::read_file(crai, header.references().len()).unwrap();
    let options = CramOptions {
        name_prefix: cram
            .file_name()
            .unwrap()
            .to_str()
            .unwrap()
            .as_bytes()
            .to_vec(),
        ..CramOptions::default()
    };
    CramReader::open(
        FileSource::open(cram).unwrap(),
        index,
        fasta.map(fasta_reader),
        options,
    )
    .unwrap()
}

/// Queries `region` with `fields` and rebuilds htslib-style records
/// (fields not asked for stay empty).
pub fn query<R: AlignmentRead>(
    r: &mut R,
    tid: i32,
    region: Region<'_>,
    fields: AlignmentFields,
) -> Vec<Rec> {
    let size = r.query(region, fields).unwrap_or_else(|error| {
        panic!("query tid={tid}, region={region:?}, fields={fields:?}: {error}")
    });
    let mut entries = vec![Default::default(); size.records];
    let n = size.records + 1;
    let want = |f| fields.contains(f);
    let (mut no, mut nd) = (vec![0i64; n], vec![0u8; size.name_bytes]);
    let (mut so, mut sd) = (vec![0i64; n], vec![0u8; size.seq_bytes]);
    let (mut qo, mut qd) = (vec![0i64; n], vec![0u8; size.seq_bytes]);
    let (mut co, mut cd) = (vec![0i64; n], vec![0u32; size.cigar_ops]);
    let (mut ao, mut ad) = (vec![0i64; n], vec![0u8; size.aux_bytes]);
    let (mut mo, mut md) = (vec![0i64; n], vec![Default::default(); size.modifications]);
    r.write(AlignmentsOut {
        entries: &mut entries,
        names: want(AlignmentFields::NAME).then(|| StringColumn {
            offsets: &mut no,
            data: &mut nd,
        }),
        seqs: want(AlignmentFields::SEQ).then(|| StringColumn {
            offsets: &mut so,
            data: &mut sd,
        }),
        quals: want(AlignmentFields::QUAL).then(|| StringColumn {
            offsets: &mut qo,
            data: &mut qd,
        }),
        cigars: want(AlignmentFields::CIGAR).then(|| CigarColumn {
            offsets: &mut co,
            ops: &mut cd,
        }),
        aux: want(AlignmentFields::AUX).then(|| StringColumn {
            offsets: &mut ao,
            data: &mut ad,
        }),
        modifications: want(AlignmentFields::MODIFICATIONS).then(|| ModificationColumn {
            offsets: &mut mo,
            entries: &mut md,
        }),
    })
    .unwrap();
    for (field, offsets, data_len) in [
        (AlignmentFields::NAME, &no, nd.len()),
        (AlignmentFields::SEQ, &so, sd.len()),
        (AlignmentFields::QUAL, &qo, qd.len()),
        (AlignmentFields::CIGAR, &co, cd.len()),
        (AlignmentFields::AUX, &ao, ad.len()),
        (AlignmentFields::MODIFICATIONS, &mo, md.len()),
    ] {
        if want(field) {
            assert_eq!(offsets[0], 0, "{region:?} {field:?}: first offset");
            assert!(
                offsets.windows(2).all(|pair| pair[0] <= pair[1]),
                "{region:?} {field:?}: offsets {offsets:?}"
            );
            assert_eq!(
                *offsets.last().unwrap(),
                data_len as i64,
                "{region:?} {field:?}: terminal offset"
            );
        }
    }
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
            name: if want(AlignmentFields::NAME) {
                nd[span(&no, i)].to_vec()
            } else {
                Vec::new()
            },
            cigar: if want(AlignmentFields::CIGAR) {
                cd[span(&co, i)].to_vec()
            } else {
                Vec::new()
            },
            seq: if want(AlignmentFields::SEQ) {
                sd[span(&so, i)].to_vec()
            } else {
                Vec::new()
            },
            qual: if want(AlignmentFields::QUAL) {
                qd[span(&qo, i)].to_vec()
            } else {
                Vec::new()
            },
            aux: if want(AlignmentFields::AUX) {
                ad[span(&ao, i)].to_vec()
            } else {
                Vec::new()
            },
        })
        .collect()
}

/// `want` with the fields not in `fields` emptied.
pub fn only(want: &[Rec], fields: AlignmentFields) -> Vec<Rec> {
    want.iter()
        .map(|r| {
            let keep = |f| fields.contains(f);
            Rec {
                tid: r.tid,
                pos: r.pos,
                end: r.end,
                flag: r.flag,
                mapq: r.mapq,
                mtid: r.mtid,
                mpos: r.mpos,
                isize: r.isize,
                name: if keep(AlignmentFields::NAME) {
                    r.name.clone()
                } else {
                    Vec::new()
                },
                cigar: if keep(AlignmentFields::CIGAR) {
                    r.cigar.clone()
                } else {
                    Vec::new()
                },
                seq: if keep(AlignmentFields::SEQ) {
                    r.seq.clone()
                } else {
                    Vec::new()
                },
                qual: if keep(AlignmentFields::QUAL) {
                    r.qual.clone()
                } else {
                    Vec::new()
                },
                aux: if keep(AlignmentFields::AUX) {
                    r.aux.clone()
                } else {
                    Vec::new()
                },
            }
        })
        .collect()
}

fn regions(rng: &mut Rng, refs: &[(String, u64)], n: usize) -> Vec<(String, u32, u32)> {
    let mut out = Vec::new();
    for (name, len) in refs {
        let len = *len as u32;
        out.push((name.clone(), 0, len));
        out.push((name.clone(), 0, 1));
        out.push((name.clone(), len - 1, len));
        out.push((name.clone(), len, len + 1000));
        out.push((name.clone(), len.saturating_sub(5_000), len + 50));
    }
    for _ in 0..n {
        let (name, len) = &refs[rng.range(0, refs.len() as u64) as usize];
        let start = rng.range(0, *len) as u32;
        let width = match rng.range(0, 4) {
            0 => rng.range(1, 10),
            1 => rng.range(10, 2_000),
            2 => rng.range(2_000, 50_000),
            _ => rng.range(50_000, 400_000),
        } as u32;
        out.push((name.clone(), start, start.saturating_add(width)));
    }
    out
}

/// Every region of `regions` through htslib and through our reader, with
/// every field and with subsets of them.
pub fn check_queries(f: &Fixture, cram: &Path, crai: &Path, with_ref: bool, seed: u64) -> usize {
    let reference = with_ref.then_some(f.fasta.as_path());
    let hts = Hts::open(cram, crai, Some(f.fasta.as_path()).filter(|_| with_ref));
    let mut ours = open(cram, crai, reference);
    let mut rng = Rng::new(seed);
    let mut total = 0;
    let subsets = [
        AlignmentFields::NONE,
        AlignmentFields::CIGAR,
        AlignmentFields::NAME | AlignmentFields::QUAL,
        AlignmentFields::SEQ,
        AlignmentFields::AUX,
    ];
    for (k, (chrom, start, end)) in regions(&mut rng, &f.refs, 60).iter().enumerate() {
        let tid = hts.tid(chrom);
        let region = Region::new(chrom, *start, *end);
        let want = hts.query(tid, *start as i64, *end as i64);
        let got = query(&mut ours, tid, region, AlignmentFields::ALL);
        assert_eq!(got.len(), want.len(), "{cram:?} {chrom}:{start}-{end}");
        for (record, (g, w)) in got.iter().zip(&want).enumerate() {
            assert_eq!(
                g, w,
                "{cram:?} seed={seed} region#{k} {chrom}:{start}-{end} record#{record}"
            );
        }
        let fields = subsets[k % subsets.len()];
        assert_eq!(
            query(&mut ours, tid, region, fields),
            only(&want, fields),
            "{cram:?} seed={seed} region#{k} {chrom}:{start}-{end} {fields:?}"
        );
        total += got.len();
    }
    total
}
