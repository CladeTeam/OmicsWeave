//! CRAM against htslib: region queries return the records
//! `sam_itr_queryi` returns, field for field, across CRAM versions,
//! codecs, embedded and absent references and multi-reference slices;
//! indexes built here hold what `samtools index` builds.

mod common;

use std::path::{Path, PathBuf};

use common::{
    Hts, Rec, Rng, SamSpec, crai_of, cram_of, gunzip, random_sam, reference, write_fasta,
};
use omicsweave_core::{
    AlignmentFields, AlignmentRead, AlignmentsOut, CigarColumn, Error, FileSource, Region,
    StringColumn,
};
use omicsweave_cram::{CramIndex, CramOptions, CramReader, NoReference, Scanner, index_entries};
use omicsweave_fasta::FastaReader;

struct Fixture {
    dir: tempfile::TempDir,
    refs: Vec<(String, u64)>,
    fasta: PathBuf,
    sam: PathBuf,
}

fn fixture(seed: u64, refs: Vec<(String, u64)>, spec: SamSpec) -> Fixture {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(seed);
    let seqs = reference(&mut rng, &refs);
    let fasta = dir.path().join("ref.fa");
    write_fasta(&fasta, &refs, &seqs);
    let text = random_sam(&mut rng, &refs, &seqs, &spec);
    let sam = dir.path().join("in.sam");
    std::fs::write(&sam, text).unwrap();
    Fixture {
        dir,
        refs,
        fasta,
        sam,
    }
}

fn standard_refs() -> Vec<(String, u64)> {
    vec![
        ("chr1".into(), 600_000),
        ("chr2".into(), 250_000),
        ("chrM".into(), 16_569),
        ("empty".into(), 30_000),
    ]
}

impl Fixture {
    /// The CRAM written with `options` and its htslib index.
    fn cram(&self, name: &str, options: &[&str], with_ref: bool) -> (PathBuf, PathBuf) {
        let cram = self.dir.path().join(format!("{name}.cram"));
        let reference = with_ref.then_some(self.fasta.as_path());
        cram_of(&self.sam, &cram, reference, options);
        let crai = self.dir.path().join(format!("{name}.cram.crai"));
        crai_of(&cram, &crai);
        (cram, crai)
    }
}

type Reader = CramReader<FileSource, FastaReader<FileSource>>;

fn open(cram: &Path, crai: &Path, fasta: Option<&Path>) -> Reader {
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
        fasta.map(common::fasta_reader),
        options,
    )
    .unwrap()
}

/// Queries `region` with `fields` and rebuilds htslib-style records
/// (fields not asked for stay empty).
fn query<R: AlignmentRead>(
    r: &mut R,
    tid: i32,
    region: Region<'_>,
    fields: AlignmentFields,
) -> Vec<Rec> {
    let size = r.query(region, fields).unwrap();
    let mut entries = vec![Default::default(); size.records];
    let n = size.records + 1;
    let want = |f| fields.contains(f);
    let (mut no, mut nd) = (vec![0i64; n], vec![0u8; size.name_bytes]);
    let (mut so, mut sd) = (vec![0i64; n], vec![0u8; size.seq_bytes]);
    let (mut qo, mut qd) = (vec![0i64; n], vec![0u8; size.seq_bytes]);
    let (mut co, mut cd) = (vec![0i64; n], vec![0u32; size.cigar_ops]);
    let (mut ao, mut ad) = (vec![0i64; n], vec![0u8; size.aux_bytes]);
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
    })
    .unwrap();
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
fn only(want: &[Rec], fields: AlignmentFields) -> Vec<Rec> {
    want.iter()
        .map(|r| {
            let mut r = r.clone();
            let keep = |f| fields.contains(f);
            if !keep(AlignmentFields::NAME) {
                r.name.clear();
            }
            if !keep(AlignmentFields::CIGAR) {
                r.cigar.clear();
            }
            if !keep(AlignmentFields::SEQ) {
                r.seq.clear();
            }
            if !keep(AlignmentFields::QUAL) {
                r.qual.clear();
            }
            if !keep(AlignmentFields::AUX) {
                r.aux.clear();
            }
            r
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
fn check_queries(f: &Fixture, cram: &Path, crai: &Path, with_ref: bool, seed: u64) -> usize {
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
        for (g, w) in got.iter().zip(&want) {
            assert_eq!(g, w, "{cram:?} {chrom}:{start}-{end}");
        }
        let fields = subsets[k % subsets.len()];
        assert_eq!(
            query(&mut ours, tid, region, fields),
            only(&want, fields),
            "{cram:?} {chrom}:{start}-{end} {fields:?}"
        );
        total += got.len();
    }
    total
}

#[test]
fn queries_match_htslib_across_versions_and_codecs() {
    let f = fixture(
        1,
        standard_refs(),
        SamSpec {
            per_ref: 3_000,
            gap: 50,
            long_cigars: true,
            unplaced: 20,
        },
    );
    let configs: [(&str, &[&str]); 9] = [
        ("v30", &["version=3.0"]),
        ("v31", &["version=3.1"]),
        ("gzip", &["version=3.0", "use_rans=0", "level=1"]),
        (
            "bz2lzma",
            &[
                "version=3.0",
                "use_rans=0",
                "use_bzip2=1",
                "use_lzma=1",
                "level=9",
            ],
        ),
        (
            "arith",
            &[
                "version=3.1",
                "use_arith=1",
                "use_fqz=1",
                "use_tok=1",
                "level=9",
            ],
        ),
        ("archive", &["version=3.1", "archive", "level=9"]),
        (
            "small",
            &[
                "version=3.1",
                "seqs_per_slice=300",
                "slices_per_container=3",
            ],
        ),
        ("fast", &["version=3.1", "fast", "use_tok=0"]),
        (
            "fqz",
            &["version=3.1", "use_rans=0", "use_fqz=1", "level=1"],
        ),
    ];
    let mut methods = std::collections::BTreeSet::new();
    for (i, (name, opts)) in configs.iter().enumerate() {
        let (cram, crai) = f.cram(name, opts, true);
        methods.extend(block_methods(&cram));
        let n = check_queries(&f, &cram, &crai, true, 100 + i as u64);
        assert!(n > 5_000, "{name}: {n} records");
    }
    // Every block compression method was decoded at least once.
    assert_eq!(
        methods.into_iter().collect::<Vec<_>>(),
        (0..=8).collect::<Vec<u8>>(),
        "methods covered"
    );
}

/// Compression methods of every block of a file.
fn block_methods(cram: &Path) -> Vec<u8> {
    let source = FileSource::open(cram).unwrap();
    let header = omicsweave_cram::CramHeader::read(&source).unwrap();
    let mut scan = Scanner::new(source, &header, true);
    let mut out = Vec::new();
    while let Some(c) = scan.next_container().unwrap() {
        out.extend_from_slice(c.block_methods());
    }
    out
}

#[test]
fn embedded_and_absent_references() {
    let f = fixture(
        2,
        standard_refs(),
        SamSpec {
            per_ref: 1_500,
            gap: 200,
            long_cigars: false,
            unplaced: 5,
        },
    );
    let (cram, crai) = f.cram("embed", &["version=3.1", "embed_ref=1"], true);
    assert!(check_queries(&f, &cram, &crai, false, 7) > 1_000);
    let (cram, crai) = f.cram("noref", &["version=3.0", "no_ref=1"], false);
    assert!(check_queries(&f, &cram, &crai, false, 8) > 1_000);
}

#[test]
fn multi_reference_slices() {
    let refs: Vec<(String, u64)> = (0..40)
        .map(|i| (format!("ctg{i}"), 6_000 + 97 * i))
        .collect();
    let f = fixture(
        3,
        refs,
        SamSpec {
            per_ref: 12,
            gap: 100,
            long_cigars: false,
            unplaced: 10,
        },
    );
    for (name, opts) in [
        (
            "multi",
            &["version=3.1", "multi_seq_per_slice=1", "seqs_per_slice=150"][..],
        ),
        (
            "multi30",
            &[
                "version=3.0",
                "multi_seq_per_slice=1",
                "seqs_per_slice=1000",
                "embed_ref=1",
            ],
        ),
    ] {
        let (cram, crai) = f.cram(name, opts, true);
        assert!(check_queries(&f, &cram, &crai, true, 9) > 300, "{name}");
        let want = String::from_utf8(gunzip(&std::fs::read(&crai).unwrap())).unwrap();
        let got: String = index_entries(FileSource::open(&cram).unwrap())
            .unwrap()
            .iter()
            .map(|e| e.line())
            .collect();
        assert_eq!(got, want, "{name}");
        assert!(
            want.lines().any(|l| l.starts_with("-1\t")),
            "{name}: unplaced entries"
        );
    }
}

#[test]
fn indexes_match_samtools_index() {
    let f = fixture(
        4,
        standard_refs(),
        SamSpec {
            per_ref: 2_000,
            gap: 30,
            long_cigars: true,
            unplaced: 30,
        },
    );
    for (name, opts) in [
        ("a", &["version=3.0"][..]),
        (
            "b",
            &[
                "version=3.1",
                "seqs_per_slice=200",
                "slices_per_container=4",
            ],
        ),
        ("c", &["version=3.1", "no_ref=1"]),
    ] {
        let (cram, crai) = f.cram(name, opts, name != "c");
        let want = gunzip(&std::fs::read(&crai).unwrap());
        let built = omicsweave_cram::build_index(FileSource::open(&cram).unwrap()).unwrap();
        assert_eq!(gunzip(&built), want, "{name}");
        let header = omicsweave_cram::CramHeader::read(&FileSource::open(&cram).unwrap()).unwrap();
        let parsed = CramIndex::parse(&built, header.references().len()).unwrap();
        assert_eq!(
            parsed.entries().len(),
            want.split(|&c| c == b'\n')
                .filter(|l| !l.is_empty())
                .count()
        );
    }
}

#[test]
fn fixed_fields_need_no_reference() {
    let f = fixture(
        5,
        standard_refs(),
        SamSpec {
            per_ref: 1_000,
            gap: 80,
            long_cigars: true,
            unplaced: 0,
        },
    );
    let (cram, crai) = f.cram("plain", &["version=3.1"], true);
    let hts = Hts::open(&cram, &crai, Some(&f.fasta));
    let mut bare: CramReader<FileSource, NoReference> = {
        let header = omicsweave_cram::CramHeader::read(&FileSource::open(&cram).unwrap()).unwrap();
        let index = CramIndex::read_file(&crai, header.references().len()).unwrap();
        CramReader::open(
            FileSource::open(&cram).unwrap(),
            index,
            None,
            CramOptions::default(),
        )
        .unwrap()
    };
    let fields = AlignmentFields::CIGAR | AlignmentFields::QUAL;
    for (chrom, len) in &f.refs {
        let tid = hts.tid(chrom);
        let want = hts.query(tid, 0, *len as i64);
        let got = query(&mut bare, tid, Region::new(chrom, 0, *len as u32), fields);
        assert_eq!(got, only(&want, fields), "{chrom}");
    }
    let err = bare
        .query(Region::new("chr1", 0, 600_000), AlignmentFields::SEQ)
        .unwrap_err();
    assert!(matches!(err, Error::Invalid { .. }), "{err}");
}

#[test]
fn reference_checks_and_errors() {
    let f = fixture(
        6,
        standard_refs(),
        SamSpec {
            per_ref: 300,
            gap: 80,
            long_cigars: false,
            unplaced: 0,
        },
    );
    let (cram, crai) = f.cram("x", &["version=3.0"], true);
    // A different reference fails the slice MD5 check, unless told not to.
    let mut rng = Rng::new(99);
    let other = reference(&mut rng, &f.refs);
    let wrong = f.dir.path().join("wrong.fa");
    write_fasta(&wrong, &f.refs, &other);
    let mut r = open(&cram, &crai, Some(&wrong));
    let err = r
        .query(Region::new("chr1", 0, 600_000), AlignmentFields::SEQ)
        .unwrap_err();
    assert!(err.to_string().contains("MD5"), "{err}");
    assert!(
        r.query(Region::new("chr1", 0, 600_000), AlignmentFields::NONE)
            .is_ok()
    );
    assert!(matches!(
        r.query(Region::new("nope", 0, 10), AlignmentFields::NONE),
        Err(Error::UnknownChrom(_))
    ));
    assert_eq!(
        r.query(Region::new("chr1", 500, 500), AlignmentFields::ALL)
            .unwrap()
            .records,
        0
    );
    assert!(
        r.write(AlignmentsOut {
            entries: &mut [],
            names: None,
            seqs: None,
            quals: None,
            cigars: Some(CigarColumn {
                offsets: &mut [0],
                ops: &mut []
            }),
            aux: None,
        })
        .is_err()
    );

    // Versions other than 3.0 / 3.1 are unsupported.
    let mut bytes = std::fs::read(&cram).unwrap();
    for (major, minor) in [(2u8, 1u8), (4, 0)] {
        bytes[4] = major;
        bytes[5] = minor;
        let err = omicsweave_cram::CramHeader::read(&bytes[..]).unwrap_err();
        assert!(matches!(err, Error::Unsupported(_)), "{err}");
    }
}
