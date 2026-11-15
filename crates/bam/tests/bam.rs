//! BAM against htslib: region queries return the records `sam_itr_queryi`
//! returns, field for field, through BAI and CSI indexes; indexes built
//! here hold what `samtools index` builds.

mod common;

use std::path::{Path, PathBuf};
use std::sync::Arc;

use common::{Hts, Rec, Rng, SamSpec, blocks, parse_index, random_sam, read_all};
use omicsweave_bam::{BamReader, BinIndex, IndexChoice, IndexKind, Record, build_index};
use omicsweave_core::cache::BlockCache;
use omicsweave_core::exec::{Inline, Pool};
use omicsweave_core::{
    AlignmentFields, AlignmentRead, AlignmentSize, AlignmentsOut, CigarColumn, Error, FileSource,
    Region, StringColumn,
};

struct Fixture {
    _dir: tempfile::TempDir,
    bam: PathBuf,
    bai: PathBuf,
    csi: PathBuf,
    csi12: PathBuf,
    refs: Vec<(&'static str, u64)>,
}

fn fixture(seed: u64, per_ref: usize, gap: u64, level: u32) -> Fixture {
    let dir = tempfile::tempdir().unwrap();
    let refs = vec![
        ("chr1", 3_000_000),
        ("chr2", 1_200_000),
        ("chrM", 16_569),
        ("empty", 50_000),
    ];
    let mut rng = Rng::new(seed);
    let mut spec_refs = refs.clone();
    spec_refs.pop();
    let text = random_sam(
        &mut rng,
        &SamSpec {
            refs: spec_refs,
            per_ref,
            gap,
            long_cigars: true,
        },
    );
    // The empty reference stays in the header without records.
    let mut header_fix = String::from_utf8(text).unwrap();
    header_fix = header_fix.replacen("@RG\t", "@SQ\tSN:empty\tLN:50000\n@RG\t", 1);
    let sam = dir.path().join("in.sam");
    std::fs::write(&sam, header_fix).unwrap();
    let bam = dir.path().join("in.bam");
    common::bam_of(&sam, &bam, level);
    let bai = dir.path().join("in.bam.bai");
    common::index(&bam, &bai, 0);
    let csi = dir.path().join("in.bam.csi");
    common::index(&bam, &csi, 14);
    let csi12 = dir.path().join("in12.csi");
    common::index(&bam, &csi12, 12);
    Fixture {
        _dir: dir,
        bam,
        bai,
        csi,
        csi12,
        refs,
    }
}

fn open(bam: &Path, idx: &Path) -> BamReader<FileSource> {
    BamReader::open(
        FileSource::open(bam).unwrap(),
        BinIndex::read_file(idx).unwrap(),
    )
    .unwrap()
}

/// Queries `region` with every field and rebuilds htslib-style records.
fn query_all(r: &mut BamReader<FileSource>, region: Region<'_>) -> Vec<Rec> {
    let size = r.query(region, AlignmentFields::ALL).unwrap();
    let mut entries = vec![Default::default(); size.records];
    let n = size.records + 1;
    let (mut no, mut nd) = (vec![0i64; n], vec![0u8; size.name_bytes]);
    let (mut so, mut sd) = (vec![0i64; n], vec![0u8; size.seq_bytes]);
    let (mut qo, mut qd) = (vec![0i64; n], vec![0u8; size.seq_bytes]);
    let (mut co, mut cd) = (vec![0i64; n], vec![0u32; size.cigar_ops]);
    let (mut ao, mut ad) = (vec![0i64; n], vec![0u8; size.aux_bytes]);
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

fn regions(rng: &mut Rng, refs: &[(&'static str, u64)], n: usize) -> Vec<(&'static str, u32, u32)> {
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

#[test]
fn queries_match_htslib_through_bai_and_csi() {
    for (seed, per_ref, gap, level) in [(1u64, 6_000, 30, 6), (2, 2_500, 400, 1), (3, 12_000, 5, 0)]
    {
        let f = fixture(seed, per_ref, gap, level);
        let mut rng = Rng::new(seed + 100);
        let regions = regions(&mut rng, &f.refs, 120);
        let blocks = common::blocks(&std::fs::read(&f.bam).unwrap()).len();
        assert!(blocks > 50, "{blocks} blocks");
        for idx in [&f.bai, &f.csi, &f.csi12] {
            let hts = Hts::open(&f.bam, idx);
            let mut ours = open(&f.bam, idx);
            let (mut total, mut nonempty) = (0, 0);
            for &(chrom, start, end) in &regions {
                let want = hts.query(hts.tid(chrom), start as i64, end as i64);
                let got = query_all(&mut ours, Region::new(chrom, start, end));
                assert_eq!(got.len(), want.len(), "{idx:?} {chrom}:{start}-{end}");
                assert_eq!(got, want, "{idx:?} {chrom}:{start}-{end}");
                total += got.len();
                nonempty += !got.is_empty() as usize;
            }
            assert!(
                total > 20_000 && nonempty > 60,
                "{total} records, {nonempty} regions"
            );
        }
    }
}

#[test]
fn fields_are_optional_and_sizes_exact() {
    let f = fixture(7, 3_000, 50, 6);
    let mut r = open(&f.bam, &f.bai);
    let region = Region::new("chr1", 0, 3_000_000);
    let all = r.query(region, AlignmentFields::ALL).unwrap();
    let only = r.query(region, AlignmentFields::NONE).unwrap();
    assert_eq!(
        only,
        AlignmentSize {
            records: all.records,
            ..Default::default()
        }
    );
    let mut entries = vec![Default::default(); only.records];
    r.write(AlignmentsOut {
        entries: &mut entries,
        names: None,
        seqs: None,
        quals: None,
        cigars: None,
        aux: None,
    })
    .unwrap();
    let names = r.query(region, AlignmentFields::NAME).unwrap();
    let mut offsets = vec![0i64; names.records + 1];
    let mut data = vec![0u8; names.name_bytes];
    let mut entries = vec![Default::default(); names.records];
    // A buffer for a field the query did not ask for is refused.
    let mut so = vec![0i64; names.records + 1];
    let err = r
        .write(AlignmentsOut {
            entries: &mut entries,
            names: Some(StringColumn {
                offsets: &mut offsets,
                data: &mut data,
            }),
            seqs: Some(StringColumn {
                offsets: &mut so,
                data: &mut [],
            }),
            quals: None,
            cigars: None,
            aux: None,
        })
        .unwrap_err();
    assert!(matches!(err, Error::Invalid { .. }), "{err}");
    let mut short = vec![Default::default(); names.records.saturating_sub(1)];
    let err = r
        .write(AlignmentsOut {
            entries: &mut short,
            names: Some(StringColumn {
                offsets: &mut offsets,
                data: &mut data,
            }),
            seqs: None,
            quals: None,
            cigars: None,
            aux: None,
        })
        .unwrap_err();
    assert!(matches!(err, Error::OutputSize { .. }), "{err}");
}

#[test]
fn unknown_and_empty_regions() {
    let f = fixture(8, 500, 100, 6);
    let mut r = open(&f.bam, &f.bai);
    assert!(matches!(
        r.query(Region::new("nope", 0, 10), AlignmentFields::NONE),
        Err(Error::UnknownChrom(_))
    ));
    assert_eq!(
        r.query(Region::new("chr1", 500, 500), AlignmentFields::NONE)
            .unwrap()
            .records,
        0
    );
    assert!(matches!(
        r.query(Region::new("chr1", 5, 1), AlignmentFields::NONE),
        Err(Error::InvalidRegion { .. })
    ));
    // A failed query leaves nothing to write.
    let mut entries = Vec::new();
    assert!(matches!(
        r.write(AlignmentsOut {
            entries: &mut entries,
            names: None,
            seqs: None,
            quals: None,
            cigars: None,
            aux: None,
        }),
        Err(Error::NotQueried)
    ));
    assert_eq!(r.chrom_len("chrM").unwrap(), 16_569);
}

#[test]
fn built_indexes_match_samtools_index() {
    for (seed, level) in [(11u64, 6), (12, 0), (13, 9)] {
        let f = fixture(seed, 8_000, 20, level);
        let theirs_bai = parse_index(&std::fs::read(&f.bai).unwrap());
        let theirs_csi = parse_index(&std::fs::read(&f.csi).unwrap());
        let theirs_csi12 = parse_index(&std::fs::read(&f.csi12).unwrap());
        let raw = std::fs::read(&f.bam).unwrap();
        let pool = Pool::new(4);
        for (choice, kind, want) in [
            (IndexChoice::Bai, IndexKind::Bai, &theirs_bai),
            (IndexChoice::Auto, IndexKind::Bai, &theirs_bai),
            (
                IndexChoice::Csi { min_shift: 14 },
                IndexKind::Csi,
                &theirs_csi,
            ),
            (
                IndexChoice::Csi { min_shift: 12 },
                IndexKind::Csi,
                &theirs_csi12,
            ),
        ] {
            let (k, bytes) = build_index(&raw[..], Box::new(Inline::new()), choice).unwrap();
            assert_eq!(k, kind);
            assert_eq!(&parse_index(&bytes), want, "{choice:?}");
            let (_, threaded) = build_index(&raw[..], Box::new(pool.ordered(8)), choice).unwrap();
            assert_eq!(threaded, bytes);
            // Our index reads back as htslib's does.
            BinIndex::parse(&bytes).unwrap();
        }
    }
}

/// A copy of `bam` with an empty BGZF block after every `every`th block.
fn with_empty_blocks(bam: &Path, out: &Path, every: usize) {
    let raw = std::fs::read(bam).unwrap();
    let mut spliced = Vec::new();
    for (i, b) in blocks(&raw).into_iter().enumerate() {
        spliced.extend_from_slice(b);
        if i % every == 0 {
            spliced.extend_from_slice(&omicsweave_bgzf::EOF_BLOCK);
        }
    }
    std::fs::write(out, spliced).unwrap();
}

#[test]
fn empty_blocks_inside_the_file() {
    let f = fixture(21, 6_000, 25, 6);
    let spliced = f.bam.with_file_name("spliced.bam");
    with_empty_blocks(&f.bam, &spliced, 3);
    let bai = f.bam.with_file_name("spliced.bam.bai");
    common::index(&spliced, &bai, 0);
    let raw = std::fs::read(&spliced).unwrap();
    let (_, ours) = build_index(&raw[..], Box::new(Inline::new()), IndexChoice::Bai).unwrap();
    assert_eq!(
        parse_index(&ours),
        parse_index(&std::fs::read(&bai).unwrap())
    );
    let hts = Hts::open(&spliced, &bai);
    let mut r = open(&spliced, &bai);
    let mut rng = Rng::new(5);
    for (chrom, start, end) in regions(&mut rng, &f.refs, 60) {
        let want = hts.query(hts.tid(chrom), start as i64, end as i64);
        assert_eq!(query_all(&mut r, Region::new(chrom, start, end)), want);
    }
}

#[test]
fn records_with_long_cigars_use_their_cg_tag() {
    let f = fixture(41, 4_000, 30, 6);
    let all = read_all(&f.bam);
    assert!(
        all.iter().any(|r| r.cigar.len() > 65_535),
        "fixture has a long CIGAR"
    );
    let raw = std::fs::read(&f.bam).unwrap();
    let mut scan = omicsweave_bam::Scanner::new(&raw[..], Box::new(Inline::new()));
    let header = scan.header().unwrap();
    let n_refs = header.references().len();
    let mut i = 0;
    while let Some(range) = scan.next_record().unwrap() {
        let rec = Record::parse(scan.slice(range), n_refs).unwrap();
        let want = &all[i];
        assert_eq!(rec.n_cigar(), want.cigar.len());
        assert_eq!(rec.end(), want.end);
        let (a, b) = rec.aux();
        assert_eq!([a, b].concat(), want.aux);
        i += 1;
    }
    assert_eq!(i, all.len());
}

#[test]
fn external_blocks_prefetch_clone_and_budget_fallback_preserve_records() {
    let f = fixture(73, 4_000, 5, 6);
    let region = Region::new("chr1", 0, 100_000);
    let mut reader = open(&f.bam, &f.bai);
    let expected = query_all(&mut reader, region);
    assert!(!expected.is_empty());
    reader.reset_stats();
    assert_eq!(query_all(&mut reader, region), expected);
    assert!(
        reader.stats().inflate_calls > 0,
        "default queries do not retain decoded blocks"
    );

    let cache = Arc::new(BlockCache::new(64 << 20, 8192));
    reader.set_cache(Some(Arc::clone(&cache)));
    let mut keys = Vec::new();
    reader.describe_cache(region, &mut keys).unwrap();
    assert!(!keys.is_empty());
    assert!(
        keys.iter()
            .all(|key| key.source == reader.cache_namespace() && key.kind == 5)
    );
    assert!(keys.windows(2).all(|pair| pair[0] != pair[1]));
    reader.prefetch_cache(region).unwrap();
    assert!(cache.stats().loads > 0);
    reader.reset_stats();
    assert_eq!(query_all(&mut reader, region), expected);
    assert_eq!(reader.stats().read_calls, 0);
    assert_eq!(reader.stats().inflate_calls, 0);
    let mut clone = reader.try_clone().unwrap();
    assert!(std::ptr::eq(reader.header(), clone.header()));
    assert!(std::ptr::eq(reader.index(), clone.index()));
    assert_eq!(reader.cache_namespace(), clone.cache_namespace());
    drop(reader);
    assert_eq!(query_all(&mut clone, region), expected);
    assert_eq!(clone.stats().read_calls, 0);
    assert_eq!(clone.stats().inflate_calls, 0);
    cache.clear_source(clone.cache_namespace());
    clone.reset_stats();
    assert_eq!(query_all(&mut clone, region), expected);
    assert!(clone.stats().inflate_calls > 0);
    for (bytes, entries) in [(0, 0), (1, 1)] {
        let limited = Arc::new(BlockCache::new(bytes, entries));
        clone.set_cache(Some(Arc::clone(&limited)));
        clone.reset_stats();
        assert_eq!(query_all(&mut clone, region), expected);
        assert!(clone.stats().inflate_calls > 0);
        assert_eq!(limited.stats().resident_entries, 0);
    }
}

#[test]
fn wrong_cached_payload_fails_safely_and_can_be_evicted() {
    let f = fixture(77, 500, 30, 6);
    let mut reader = open(&f.bam, &f.bai);
    let region = Region::new("chr1", 0, 100_000);
    let expected = query_all(&mut reader, region);
    let cache = Arc::new(BlockCache::new(64 << 20, 8192));
    reader.set_cache(Some(Arc::clone(&cache)));
    let mut keys = Vec::new();
    reader.describe_cache(region, &mut keys).unwrap();
    assert!(!keys.is_empty());
    cache.get_or_load(keys[0], || Ok((7u64, 0))).unwrap();
    assert!(reader.query(region, AlignmentFields::ALL).is_err());
    let mut entries = Vec::new();
    assert!(matches!(
        reader.write(AlignmentsOut {
            entries: &mut entries,
            names: None,
            seqs: None,
            quals: None,
            cigars: None,
            aux: None,
        }),
        Err(Error::NotQueried)
    ));
    assert!(cache.evict(keys[0]));
    assert_eq!(query_all(&mut reader, region), expected);
}

#[test]
fn selected_columns_match_htslib_with_guarded_outputs() {
    let f = fixture(101, 24, 50, 0);
    let hts = Hts::open(&f.bam, &f.bai);
    let mut reader = open(&f.bam, &f.bai);
    for region in [
        Region::new("chr1", 0, 3_000_000),
        Region::new("empty", 0, 50_000),
    ] {
        let expected = hts.query(
            hts.tid(region.chrom),
            region.start as i64,
            region.end as i64,
        );
        for bits in 0..=AlignmentFields::ALL.bits() {
            let fields = AlignmentFields::from_bits(bits).unwrap();
            let size = reader.query(region, fields).unwrap();
            assert_eq!(size.records, expected.len());
            let n = size.records + 1;
            let mut entries = vec![omicsweave_core::AlignmentEntry::default(); size.records + 2];
            entries[0].start = -777;
            entries[size.records + 1].start = -777;
            let (mut no, mut nd) = (vec![-777; n + 2], vec![0x5a; size.name_bytes + 2]);
            let (mut so, mut sd) = (vec![-777; n + 2], vec![0x5a; size.seq_bytes + 2]);
            let (mut qo, mut qd) = (vec![-777; n + 2], vec![0x5a; size.seq_bytes + 2]);
            let (mut ao, mut ad) = (vec![-777; n + 2], vec![0x5a; size.aux_bytes + 2]);
            let (mut co, mut cd) = (vec![-777; n + 2], vec![0xdeadbeef; size.cigar_ops + 2]);
            reader
                .write(AlignmentsOut {
                    entries: &mut entries[1..size.records + 1],
                    names: fields
                        .contains(AlignmentFields::NAME)
                        .then_some(StringColumn {
                            offsets: &mut no[1..n + 1],
                            data: &mut nd[1..size.name_bytes + 1],
                        }),
                    seqs: fields
                        .contains(AlignmentFields::SEQ)
                        .then_some(StringColumn {
                            offsets: &mut so[1..n + 1],
                            data: &mut sd[1..size.seq_bytes + 1],
                        }),
                    quals: fields
                        .contains(AlignmentFields::QUAL)
                        .then_some(StringColumn {
                            offsets: &mut qo[1..n + 1],
                            data: &mut qd[1..size.seq_bytes + 1],
                        }),
                    aux: fields
                        .contains(AlignmentFields::AUX)
                        .then_some(StringColumn {
                            offsets: &mut ao[1..n + 1],
                            data: &mut ad[1..size.aux_bytes + 1],
                        }),
                    cigars: fields
                        .contains(AlignmentFields::CIGAR)
                        .then_some(CigarColumn {
                            offsets: &mut co[1..n + 1],
                            ops: &mut cd[1..size.cigar_ops + 1],
                        }),
                })
                .unwrap();
            for (flag, offsets, data, part) in [
                (
                    AlignmentFields::NAME,
                    &no,
                    &nd,
                    (|r: &Rec| &r.name) as fn(&Rec) -> &Vec<u8>,
                ),
                (
                    AlignmentFields::SEQ,
                    &so,
                    &sd,
                    (|r: &Rec| &r.seq) as fn(&Rec) -> &Vec<u8>,
                ),
                (
                    AlignmentFields::QUAL,
                    &qo,
                    &qd,
                    (|r: &Rec| &r.qual) as fn(&Rec) -> &Vec<u8>,
                ),
                (
                    AlignmentFields::AUX,
                    &ao,
                    &ad,
                    (|r: &Rec| &r.aux) as fn(&Rec) -> &Vec<u8>,
                ),
            ] {
                assert_eq!((offsets[0], offsets[n + 1]), (-777, -777));
                assert_eq!((data[0], data[data.len() - 1]), (0x5a, 0x5a));
                if fields.contains(flag) {
                    let mut at = 0;
                    assert_eq!(offsets[1], 0);
                    for (i, record) in expected.iter().enumerate() {
                        let bytes = part(record);
                        assert_eq!(&data[1 + at..1 + at + bytes.len()], bytes);
                        at += bytes.len();
                        assert_eq!(offsets[i + 2], at as i64);
                    }
                    assert_eq!(at, data.len() - 2);
                }
            }
            assert_eq!((co[0], co[n + 1]), (-777, -777));
            assert_eq!((cd[0], cd[cd.len() - 1]), (0xdeadbeef, 0xdeadbeef));
            if fields.contains(AlignmentFields::CIGAR) {
                let want: Vec<_> = expected
                    .iter()
                    .flat_map(|r| r.cigar.iter().copied())
                    .collect();
                assert_eq!(&cd[1..cd.len() - 1], want);
                let mut at = 0;
                assert_eq!(co[1], 0);
                for (i, rec) in expected.iter().enumerate() {
                    at += rec.cigar.len();
                    assert_eq!(co[i + 2], at as i64);
                }
            }
            assert_eq!(
                (entries[0].start, entries[size.records + 1].start),
                (-777, -777)
            );
            for (entry, record) in entries[1..size.records + 1].iter().zip(&expected) {
                assert_eq!((entry.start, entry.end), (record.pos, record.end));
            }
        }
    }
}

#[test]
fn long_cigar_and_empty_sequence_outputs_match_htslib() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("long.sam");
    let bam = dir.path().join("long.bam");
    let bai = dir.path().join("long.bam.bai");
    let cigar = "1M1I".repeat(33_000) + "1M";
    let sequence = "AC".repeat(33_000) + "G";
    std::fs::write(
        &sam,
        format!(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\nlong\t0\tchr1\t1\t60\t{cigar}\t*\t0\t0\t{sequence}\t*\tNM:i:33000\nmissing\t0\tchr1\t2\t60\t1M\t*\t0\t0\t*\t*\n"
        ),
    )
    .unwrap();
    common::bam_of(&sam, &bam, 0);
    common::index(&bam, &bai, 0);
    let hts = Hts::open(&bam, &bai);
    let expected = hts.query(0, 0, 100_000);
    assert_eq!(expected[0].cigar.len(), 66_001);
    assert!(expected[1].seq.is_empty());
    let mut reader = open(&bam, &bai);
    assert_eq!(
        query_all(&mut reader, Region::new("chr1", 0, 100_000)),
        expected
    );

    let raw = std::fs::read(&bam).unwrap();
    let mut scan = omicsweave_bam::Scanner::new(&raw[..], Box::new(Inline::new()));
    scan.header().unwrap();
    let range = scan.next_record().unwrap().unwrap();
    let data = scan.slice(range);
    for cut in [0, 31, 32, 33, data.len() / 2, data.len() - 1] {
        assert!(Record::parse(&data[..cut], 1).is_err(), "cut {cut}");
    }
}
