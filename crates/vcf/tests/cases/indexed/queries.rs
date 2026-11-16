//! BCF region queries against htslib: the records `bcf_itr_queryi` returns
//! through a CSI index, with their genotypes, and the genotype matrix built
//! from them.

use crate::indexed_fixture::{Fixture, fixture};
use std::ffi::{CStr, CString, c_int, c_void};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use omicsweave_bgzf::BinIndex;
use omicsweave_core::cache::BlockCache;
use omicsweave_core::{Error, FileSource, ReadAt, Region, Result};
use omicsweave_vcf::{BcfReader, GT_ABSENT, GT_MISSING, Record, gt_allele, gt_is_phased};

const CONTIGS: [(&str, u64); 3] = [("chr1", 2_000_000), ("chr2", 600_000), ("chr3", 50_000)];
const SAMPLES: usize = 6;
const MISSING: i32 = i32::MIN;
const VECTOR_END: i32 = i32::MIN + 1;

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

fn c(s: &str) -> CString {
    CString::new(s).unwrap()
}

fn bases(rng: &mut Rng, n: usize) -> String {
    (0..n)
        .map(|_| b"ACGT"[rng.below(4) as usize] as char)
        .collect()
}

/// A sorted VCF of SNPs, indels, multiallelic sites, symbolic deletions
/// with END, mixed ploidy, missing and phased genotypes, and records
/// without genotypes.
fn random_vcf(seed: u64, per_contig: usize, gap: u64) -> String {
    let mut rng = Rng(seed | 1);
    let mut text = String::from("##fileformat=VCFv4.2\n");
    for (name, len) in CONTIGS {
        text.push_str(&format!("##contig=<ID={name},length={len}>\n"));
    }
    text.push_str("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n");
    text.push_str("##ALT=<ID=DEL,Description=\"Deletion\">\n");
    text.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    text.push_str("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n");
    text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for s in 0..SAMPLES {
        text.push_str(&format!("\ts{s}"));
    }
    text.push('\n');
    for (name, len) in CONTIGS {
        let mut pos = 1 + rng.below(gap);
        for i in 0..per_contig {
            if pos >= len - 2000 {
                break;
            }
            let kind = rng.below(10);
            let (reference, alt, info, n_alt) = match kind {
                0 => {
                    let n = 1 + rng.below(30) as usize;
                    let r = bases(&mut rng, n);
                    let a = bases(&mut rng, 1);
                    (r, a, ".".to_owned(), 1)
                }
                1 => {
                    let r = bases(&mut rng, 1);
                    let end = pos + 1 + rng.below(3000);
                    (r, "<DEL>".to_owned(), format!("END={end}"), 1)
                }
                2 => {
                    let r = bases(&mut rng, 1);
                    (r, "C,G,TT".to_owned(), ".".to_owned(), 3)
                }
                3 => (bases(&mut rng, 1), ".".to_owned(), ".".to_owned(), 0),
                _ => {
                    let r = bases(&mut rng, 1);
                    let a = if r == "A" { "G" } else { "A" };
                    (r, a.to_owned(), ".".to_owned(), 1)
                }
            };
            let no_gt = rng.below(15) == 0;
            text.push_str(&format!(
                "{name}\t{pos}\tv{i}\t{reference}\t{alt}\t50\tPASS\t{info}"
            ));
            if no_gt {
                text.push_str("\tDP");
                for _ in 0..SAMPLES {
                    text.push_str(&format!("\t{}", rng.below(40)));
                }
            } else {
                text.push_str("\tGT:DP");
                let haploid = rng.below(8) == 0;
                for _ in 0..SAMPLES {
                    let allele = |rng: &mut Rng| {
                        if rng.below(12) == 0 {
                            ".".to_owned()
                        } else {
                            rng.below(n_alt + 1).to_string()
                        }
                    };
                    let gt = if haploid || rng.below(20) == 0 {
                        allele(&mut rng)
                    } else {
                        let sep = if rng.below(2) == 0 { '|' } else { '/' };
                        format!("{}{sep}{}", allele(&mut rng), allele(&mut rng))
                    };
                    text.push_str(&format!("\t{gt}:{}", rng.below(40)));
                }
            }
            text.push('\n');
            pos += rng.below(gap);
        }
    }
    text
}

/// What a test compares of a record: position, decoded span, alleles and
/// genotype values in htslib's 32-bit encoding.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Rec {
    pos: i64,
    rlen: i64,
    alleles: Vec<Vec<u8>>,
    gt: Option<(usize, Vec<i32>)>,
}

fn htslib_query(f: &Fixture, chrom: &str, beg: i64, end: i64) -> Vec<Rec> {
    let mut out = Vec::new();
    // SAFETY: handles come from htslib and are checked before use; GT
    // values are read within the count htslib returns.
    unsafe {
        let path = c(f.bcf.to_str().unwrap());
        let fp = htslib_sys::hts_open(path.as_ptr(), c("r").as_ptr());
        assert!(!fp.is_null());
        let hdr = htslib_sys::bcf_hdr_read(fp);
        let idx = htslib_sys::ow_bcf_index_load(path.as_ptr(), c(f.csi.to_str().unwrap()).as_ptr());
        assert!(!idx.is_null());
        let tid = htslib_sys::ow_bcf_hdr_name2id(hdr, c(chrom).as_ptr());
        assert!(tid >= 0);
        let itr = htslib_sys::ow_bcf_itr_queryi(idx, tid, beg, end);
        assert!(!itr.is_null());
        let rec = htslib_sys::bcf_init();
        let mut dst: *mut c_void = std::ptr::null_mut();
        let mut ndst: c_int = 0;
        loop {
            let r = htslib_sys::ow_bcf_itr_next(fp, itr, rec);
            if r == -1 {
                break;
            }
            assert!(r >= 0, "htslib iterator error {r}");
            htslib_sys::bcf_unpack(rec, htslib_sys::BCF_UN_ALL);
            let alleles = (0..htslib_sys::ow_bcf_n_allele(rec))
                .map(|i| {
                    CStr::from_ptr(htslib_sys::ow_bcf_allele(rec, i))
                        .to_bytes()
                        .to_vec()
                })
                .collect();
            let n = htslib_sys::bcf_get_format_values(
                hdr,
                rec,
                c("GT").as_ptr(),
                &mut dst,
                &mut ndst,
                htslib_sys::BCF_HT_INT,
            );
            let gt = (n > 0).then(|| {
                let values = std::slice::from_raw_parts(dst as *const i32, n as usize).to_vec();
                (n as usize / SAMPLES, values)
            });
            out.push(Rec {
                pos: htslib_sys::ow_bcf_pos(rec),
                rlen: htslib_sys::ow_bcf_rlen(rec),
                alleles,
                gt,
            });
        }
        htslib_sys::free(dst);
        htslib_sys::bcf_destroy(rec);
        htslib_sys::ow_hts_itr_destroy(itr);
        htslib_sys::ow_hts_idx_destroy(idx);
        htslib_sys::bcf_hdr_destroy(hdr);
        htslib_sys::hts_close(fp);
    }
    out
}

fn open(f: &Fixture) -> BcfReader<FileSource> {
    let index = BinIndex::read_file(&f.csi).unwrap();
    BcfReader::open(FileSource::open(&f.bcf).unwrap(), index).unwrap()
}

fn ours_query(r: &mut BcfReader<FileSource>, chrom: &str, beg: u32, end: u32) -> Vec<Rec> {
    let mut out = Vec::new();
    let mut rec = Record::default();
    r.visit(Region::new(chrom, beg, end), &mut rec, |rec| {
        out.push(Rec {
            pos: rec.pos,
            rlen: rec.rlen,
            alleles: rec.alleles().map(<[u8]>::to_vec).collect(),
            gt: rec.genotypes().map(|g| (g.width, g.values.to_vec())),
        });
        Ok(())
    })
    .unwrap();
    out
}

fn regions(rng: &mut Rng, n: usize) -> Vec<(&'static str, u32, u32)> {
    (0..n)
        .map(|i| {
            let (name, len) = CONTIGS[rng.below(CONTIGS.len() as u64) as usize];
            let start = rng.below(len) as u32;
            let span = match i % 4 {
                0 => 1 + rng.below(100),
                1 => 1 + rng.below(20_000),
                2 => 1 + rng.below(400_000),
                _ => len,
            };
            let end = (start as u64 + span).min(len) as u32;
            (name, start, end)
        })
        .collect()
}

#[test]
fn queries_match_htslib_through_csi() {
    for (seed, per_contig, gap, min_shift) in
        [(7, 4000, 600, 14), (11, 3000, 400, 10), (13, 200, 9000, 14)]
    {
        let f = fixture(&random_vcf(seed, per_contig, gap), min_shift);
        let mut reader = open(&f);
        let mut rng = Rng(seed * 31);
        let mut seen = 0;
        for (chrom, beg, end) in regions(&mut rng, 300) {
            let want = htslib_query(&f, chrom, beg as i64, end as i64);
            let got = ours_query(&mut reader, chrom, beg, end);
            assert_eq!(got, want, "{chrom}:{beg}-{end} seed {seed}");
            seen += want.len();
        }
        assert!(seen > 0);
    }
}

#[test]
fn genotype_matrix_follows_the_records() {
    let f = fixture(&random_vcf(5, 2000, 300), 14);
    let mut reader = open(&f);
    let samples = [3, 0, 5];
    let mut rng = Rng(99);
    for (chrom, beg, end) in regions(&mut rng, 60) {
        let recs = ours_query(&mut reader, chrom, beg, end);
        let m = reader
            .genotypes(Region::new(chrom, beg, end), &samples)
            .unwrap();
        assert_eq!(m.len(), recs.len());
        assert_eq!(m.samples, samples.len());
        let ploidy = recs
            .iter()
            .filter_map(|r| r.gt.as_ref().map(|g| g.0))
            .max()
            .unwrap_or(0);
        assert_eq!(m.ploidy, ploidy);
        assert_eq!(m.gt.len(), recs.len() * samples.len() * ploidy);
        for (i, r) in recs.iter().enumerate() {
            assert_eq!(m.start[i], r.pos);
            assert_eq!(m.end[i], r.pos + r.rlen);
            assert_eq!(m.reference.get(i), &r.alleles[0][..]);
            assert_eq!(m.alt.get(i), &r.alleles[1..].join(&b','));
            for (k, &s) in samples.iter().enumerate() {
                let row = &m.gt[(i * samples.len() + k) * ploidy..][..ploidy];
                let Some((width, values)) = &r.gt else {
                    assert!(row.iter().all(|&a| a == GT_ABSENT));
                    assert!(!m.phased[i * samples.len() + k]);
                    continue;
                };
                let values = &values[s * width..(s + 1) * width];
                let mut want = vec![GT_ABSENT; ploidy];
                let mut phases = Vec::new();
                for (j, &v) in values.iter().enumerate() {
                    want[j] = match v {
                        VECTOR_END => continue,
                        MISSING => GT_MISSING,
                        v => gt_allele(v) as i8,
                    };
                    phases.push(gt_is_phased(v));
                }
                let phased = match phases.len() {
                    0 => false,
                    1 => phases[0],
                    _ => phases[1..].iter().all(|&p| p),
                };
                assert_eq!(row, &want[..]);
                assert_eq!(m.phased[i * samples.len() + k], phased);
            }
        }
    }
}

#[test]
fn unknown_contigs_samples_and_indexes() {
    let f = fixture(&random_vcf(3, 50, 1000), 14);
    let mut reader = open(&f);
    let err = reader
        .genotypes(Region::new("chrZ", 0, 10), &[0])
        .unwrap_err();
    assert!(matches!(err, Error::UnknownChrom(_)), "{err:?}");
    let err = reader
        .genotypes(Region::new("chr1", 0, 10), &[SAMPLES])
        .unwrap_err();
    assert!(matches!(err, Error::Invalid { .. }), "{err:?}");
    let empty = reader
        .genotypes(Region::new("chr1", 10, 10), &[0, 1])
        .unwrap();
    assert!(empty.is_empty() && empty.ploidy == 0);
    let bai = BinIndex::parse(&bai_bytes()).unwrap();
    assert!(BcfReader::open(FileSource::open(&f.bcf).unwrap(), bai).is_err());
}

/// A BAI of one empty reference.
fn bai_bytes() -> Vec<u8> {
    let mut b = b"BAI\x01".to_vec();
    b.extend_from_slice(&1i32.to_le_bytes());
    b.extend_from_slice(&0i32.to_le_bytes());
    b.extend_from_slice(&0i32.to_le_bytes());
    b
}

struct Counted {
    source: FileSource,
    reads: AtomicUsize,
}

impl ReadAt for Counted {
    fn len(&self) -> u64 {
        self.source.len()
    }
    fn read_exact_at(&self, out: &mut [u8], offset: u64) -> Result<()> {
        self.reads.fetch_add(1, Ordering::Relaxed);
        self.source.read_exact_at(out, offset)
    }
}

fn collect<S: ReadAt>(reader: &mut BcfReader<S>, region: Region<'_>) -> Vec<Rec> {
    let mut out = Vec::new();
    reader
        .visit(region, &mut Record::default(), |record| {
            out.push(Rec {
                pos: record.pos,
                rlen: record.rlen,
                alleles: record.alleles().map(<[u8]>::to_vec).collect(),
                gt: record.genotypes().map(|gt| (gt.width, gt.values.to_vec())),
            });
            Ok(())
        })
        .unwrap();
    out
}

#[test]
fn external_block_cache_prefetch_and_clone_preserve_indexed_record_queries() {
    let f = fixture(&random_vcf(181, 4000, 300), 14);
    let source = Arc::new(Counted {
        source: FileSource::open(&f.bcf).unwrap(),
        reads: AtomicUsize::new(0),
    });
    let index = BinIndex::read_file(&f.csi).unwrap();
    let mut reader = BcfReader::open(Arc::clone(&source), index).unwrap();
    let region = Region::new("chr1", 5_000, 240_000);
    let wanted = htslib_query(&f, "chr1", 5_000, 240_000);
    assert_eq!(collect(&mut reader, region), wanted);
    let previous = source.reads.load(Ordering::Relaxed);
    let previous_inflations = reader.stream_stats().inflate_calls;
    assert_eq!(collect(&mut reader, region), wanted);
    assert!(
        source.reads.load(Ordering::Relaxed) > previous,
        "uncached query reused decoded blocks"
    );
    assert!(reader.stream_stats().inflate_calls > previous_inflations);
    let cache = Arc::new(BlockCache::new(32 << 20, 4096));
    reader.set_cache(Some(Arc::clone(&cache)));
    let mut keys = Vec::new();
    reader.describe_cache(region, &mut keys).unwrap();
    assert!(!keys.is_empty());
    assert!(
        keys.iter()
            .all(|key| key.source == reader.cache_namespace() && key.kind == 5)
    );
    reader.prefetch_cache(region).unwrap();
    let loaded = cache.stats().loads;
    let reads = source.reads.load(Ordering::Relaxed);
    reader.reset_stats();
    assert_eq!(collect(&mut reader, region), wanted);
    assert_eq!(source.reads.load(Ordering::Relaxed), reads);
    assert_eq!(cache.stats().loads, loaded);
    assert_eq!(reader.stream_stats().inflate_calls, 0);
    assert_eq!(reader.stream_stats().read_calls, 0);
    let mut clone = reader.clone_context();
    assert_eq!(clone.cache_namespace(), reader.cache_namespace());
    drop(reader);
    assert_eq!(collect(&mut clone, region), wanted);
    assert_eq!(source.reads.load(Ordering::Relaxed), reads);
    let matrix = clone.genotypes(region, &[5, 0, 2]).unwrap();
    assert_eq!(matrix.len(), wanted.len());
    assert_eq!(source.reads.load(Ordering::Relaxed), reads);
    cache.clear_source(clone.cache_namespace());
    assert_eq!(collect(&mut clone, region), wanted);
    assert!(source.reads.load(Ordering::Relaxed) > reads);
    clone.set_cache(None);
    let before = source.reads.load(Ordering::Relaxed);
    assert_eq!(collect(&mut clone, region), wanted);
    assert!(source.reads.load(Ordering::Relaxed) > before);
}

#[test]
fn file_context_clone_zero_budget_and_empty_dependencies_are_safe() {
    let f = fixture(&random_vcf(71, 200, 600), 10);
    let mut reader = open(&f);
    let cache = Arc::new(BlockCache::new(0, 0));
    reader.set_cache(Some(Arc::clone(&cache)));
    let mut clone = reader.try_clone().unwrap();
    let region = Region::new("chr1", 500, 90_000);
    let wanted = htslib_query(&f, "chr1", 500, 90_000);
    assert_eq!(ours_query(&mut clone, "chr1", 500, 90_000), wanted);
    assert_eq!(cache.stats().resident_entries, 0);
    cache.resize(32 << 20, 4096);
    reader.prefetch_cache(region).unwrap();
    let mut keys = cache.keys();
    clone
        .describe_cache(Region::new("chr1", 1, 1), &mut keys)
        .unwrap();
    assert!(keys.is_empty());
    let independent = open(&f);
    assert_ne!(independent.cache_namespace(), clone.cache_namespace());
    drop(reader);
    assert_eq!(ours_query(&mut clone, "chr1", 500, 90_000), wanted);
    let strings = omicsweave_vcf::Strings::default();
    assert!(strings.is_empty());
    assert_eq!(strings.offsets, &[0]);
}

fn typed_gt_fixture(ty: u8) -> Fixture {
    typed_gt_fixture_layout(ty, "4.4", 3)
}

fn typed_gt_fixture_layout(ty: u8, version: &str, ploidy: usize) -> Fixture {
    use omicsweave_bgzf::Writer;
    use omicsweave_core::exec::Inline;
    use omicsweave_vcf::{Header, bcf_header_bytes};

    let mut text = "##fileformat=VCFvVERSION\n##contig=<ID=chr1,length=2000000>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String>\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        .replace("VERSION", version);
    for i in 0..SAMPLES {
        text.push_str(&format!("\ts{i}"));
    }
    text.push('\n');
    let header = Header::from_vcf_text(text.as_bytes()).unwrap();
    let key = header.key_id("GT").unwrap();
    assert!(key < 128);
    let mut shared = Vec::new();
    for value in [0i32, 100, 1, 30.0f32.to_bits() as i32] {
        shared.extend_from_slice(&value.to_le_bytes());
    }
    shared.extend_from_slice(&0u16.to_le_bytes());
    shared.extend_from_slice(&3u16.to_le_bytes());
    shared.extend_from_slice(&((1u32 << 24) | SAMPLES as u32).to_le_bytes());
    shared.extend_from_slice(&[0x17, b'.', 0x17, b'A', 0x17, b'C', 0x17, b'G', 0]);
    let mut indiv = vec![0x11, key as u8, ((ploidy as u8) << 4) | ty];
    for row in [
        [2, 5, VECTOR_END],
        [MISSING, 2, 7],
        [VECTOR_END, 5, 2],
        [3, VECTOR_END, 7],
        [0, 1, VECTOR_END],
        [7, 4, 3],
    ] {
        for value in row.into_iter().take(ploidy) {
            match ty {
                1 => indiv.push(match value {
                    MISSING => i8::MIN as u8,
                    VECTOR_END => (i8::MIN + 1) as u8,
                    v => v as u8,
                }),
                2 => indiv.extend_from_slice(
                    &match value {
                        MISSING => i16::MIN,
                        VECTOR_END => i16::MIN + 1,
                        v => v as i16,
                    }
                    .to_le_bytes(),
                ),
                3 => indiv.extend_from_slice(&value.to_le_bytes()),
                _ => unreachable!(),
            }
        }
    }
    let mut record = Vec::new();
    record.extend_from_slice(&(shared.len() as u32).to_le_bytes());
    record.extend_from_slice(&(indiv.len() as u32).to_le_bytes());
    record.extend_from_slice(&shared);
    record.extend_from_slice(&indiv);
    let dir = tempfile::tempdir().unwrap();
    let bcf = dir.path().join("typed.bcf");
    let mut writer = Writer::new(Vec::new(), 1, Box::new(Inline::new())).unwrap();
    writer.write(&bcf_header_bytes(&header).unwrap()).unwrap();
    writer.write(&record).unwrap();
    std::fs::write(&bcf, writer.finish().unwrap().out).unwrap();
    // SAFETY: the NUL-terminated path stays live during the call; htslib owns
    // and releases the index builder's internal handles.
    assert_eq!(
        unsafe {
            htslib_sys::bcf_index_build3(c(bcf.to_str().unwrap()).as_ptr(), std::ptr::null(), 14, 1)
        },
        0
    );
    let csi = dir.path().join("typed.bcf.csi");
    Fixture {
        _dir: dir,
        bcf,
        csi,
    }
}

#[test]
fn selected_genotypes_match_all_typed_widths_sentinels_and_duplicate_samples() {
    for (version, ty, ploidy) in ["4.2", "4.4"].into_iter().flat_map(|v| {
        [1, 2, 3]
            .into_iter()
            .flat_map(move |ty| [1, 2, 3].map(|p| (v, ty, p)))
    }) {
        let fixture = typed_gt_fixture_layout(ty, version, ploidy);
        let mut reader = open(&fixture);
        let wanted = htslib_query(&fixture, "chr1", 0, 1000);
        assert_eq!(wanted.len(), 1);
        assert_eq!(ours_query(&mut reader, "chr1", 0, 1000), wanted);
        let mut record = Record::default();
        reader
            .visit(Region::new("chr1", 0, 1000), &mut record, |rec| {
                assert_eq!(rec.genotype_type(), Some(ty));
                Ok(())
            })
            .unwrap();
        let samples = [5, 0, 2, 0];
        let matrix = reader
            .genotypes(Region::new("chr1", 0, 1000), &samples)
            .unwrap();
        assert_eq!(matrix.ploidy, ploidy);
        let (width, gt) = wanted[0].gt.as_ref().unwrap();
        for (i, &sample) in samples.iter().enumerate() {
            let row = &gt[sample * width..(sample + 1) * width];
            let expected: Vec<_> = row
                .iter()
                .map(|&v| match v {
                    MISSING => GT_MISSING,
                    VECTOR_END => GT_ABSENT,
                    v => gt_allele(v) as i8,
                })
                .collect();
            assert_eq!(
                &matrix.gt[i * width..(i + 1) * width],
                expected,
                "type {ty}"
            );
            let phases: Vec<_> = row
                .iter()
                .filter(|&&v| v != VECTOR_END)
                .map(|&v| gt_is_phased(v))
                .collect();
            let phased = match phases.len() {
                0 => false,
                1 => phases[0],
                _ => phases[1..].iter().all(|&v| v),
            };
            assert_eq!(matrix.phased[i], phased, "type {ty} sample {sample}");
        }
        let empty = reader.genotypes(Region::new("chr1", 0, 1000), &[]).unwrap();
        assert_eq!(empty.len(), 1);
        assert!(empty.gt.is_empty() && empty.phased.is_empty());
    }
}

#[test]
fn cached_record_spanning_bgzf_members_preserves_full_fields_and_genotypes() {
    let id = "r".repeat(180_000);
    let text = random_vcf(27, 1, 30).replacen("\tv0\t", &format!("\t{id}\t"), 1);
    let f = fixture(&text, 14);
    let source = Arc::new(Counted {
        source: FileSource::open(&f.bcf).unwrap(),
        reads: AtomicUsize::new(0),
    });
    let mut reader =
        BcfReader::open(Arc::clone(&source), BinIndex::read_file(&f.csi).unwrap()).unwrap();
    let region = Region::new("chr1", 0, 100_000);
    let wanted = htslib_query(&f, "chr1", 0, 100_000);
    assert_eq!(collect(&mut reader, region), wanted);
    let cache = Arc::new(BlockCache::new(4 << 20, 128));
    reader.set_cache(Some(Arc::clone(&cache)));
    reader.prefetch_cache(region).unwrap();
    let reads = source.reads.load(Ordering::Relaxed);
    let matrix = reader.genotypes(region, &[0]).unwrap();
    assert_eq!(matrix.id.get(0), id.as_bytes());
    assert_eq!(collect(&mut reader, region), wanted);
    assert_eq!(source.reads.load(Ordering::Relaxed), reads);
    assert!(cache.stats().resident_entries >= 3);
}

#[test]
fn invalid_csi_and_oversized_record_lengths_fail_without_trusting_allocations() {
    use omicsweave_bgzf::Writer;
    use omicsweave_core::exec::Inline;
    use std::io::Read;

    for bytes in [
        b"CSI\x01".as_slice(),
        b"CSI\x01\x0e\0\0\0\x05\0\0\0\xff\xff\xff\xff".as_slice(),
    ] {
        assert!(BinIndex::parse(bytes).is_err());
    }
    let f = typed_gt_fixture(1);
    let mut bytes = Vec::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(&f.bcf).unwrap())
        .read_to_end(&mut bytes)
        .unwrap();
    let header = 9 + u32::from_le_bytes(bytes[5..9].try_into().unwrap()) as usize;
    bytes[header..header + 8].copy_from_slice(&[0xff; 8]);
    let mut writer = Writer::new(Vec::new(), 1, Box::new(Inline::new())).unwrap();
    writer.write(&bytes).unwrap();
    std::fs::write(&f.bcf, writer.finish().unwrap().out).unwrap();
    let mut reader = open(&f);
    let error = reader
        .genotypes(Region::new("chr1", 0, 1000), &[0])
        .unwrap_err();
    assert!(
        matches!(error, Error::Invalid { .. } | Error::Truncated { .. }),
        "{error:?}"
    );
}

#[test]
fn visit_clears_record_on_empty_queries_and_failed_callbacks() {
    let fixture = typed_gt_fixture(1);
    let mut reader = open(&fixture);
    let mut rec = Record::default();
    let region = Region::new("chr1", 0, 1000);
    reader.visit(region, &mut rec, |_| Ok(())).unwrap();
    assert!(rec.n_alleles() > 0);
    reader
        .visit(Region::new("chr1", 1, 1), &mut rec, |_| {
            panic!("empty query")
        })
        .unwrap();
    assert_eq!(rec, Record::default());
    reader.visit(region, &mut rec, |_| Ok(())).unwrap();
    assert!(
        reader
            .visit(region, &mut rec, |_| Err(Error::invalid(
                "callback", "stop"
            )))
            .is_err()
    );
    assert_eq!(rec, Record::default());
    assert!(
        reader
            .visit(Region::new("missing", 0, 1), &mut rec, |_| Ok(()))
            .is_err()
    );
    assert_eq!(rec, Record::default());
    assert_eq!(
        ours_query(&mut reader, "chr1", 0, 1000),
        htslib_query(&fixture, "chr1", 0, 1000)
    );
}
