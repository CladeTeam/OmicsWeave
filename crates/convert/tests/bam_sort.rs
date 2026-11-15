//! BAM → sorted BAM: the file `samtools sort --no-PG` writes, byte for byte
//! at level 0 and block for block at other levels, for every budget, with
//! indexes equal to the ones htslib builds.

mod common;

use std::ffi::{CString, c_char, c_int};
use std::path::Path;
use std::sync::Mutex;

use common::{Rng, hts, inflate};
use omicsweave_convert::{
    BamIndex, BamSort, Budget, Conversion, Format, Index, Options, SamToBam, convert, plan,
};

/// Runs `samtools sort --no-PG -o dst args.. src`; returns its status.
fn samtools_sort(src: &Path, dst: &Path, args: &[&str]) -> i32 {
    // samtools keeps process-global state (getopt, the sort order).
    static LOCK: Mutex<()> = Mutex::new(());
    let _guard = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let prefix = dst.with_extension("samtools-tmp");
    let mut words = vec!["sort", "--no-PG", "-o", dst.to_str().unwrap()];
    words.extend(["-T", prefix.to_str().unwrap()]);
    words.extend(args);
    words.push(src.to_str().unwrap());
    let owned: Vec<CString> = words.iter().map(|w| CString::new(*w).unwrap()).collect();
    let mut argv: Vec<*mut c_char> = owned.iter().map(|c| c.as_ptr() as *mut c_char).collect();
    unsafe { samtools_sys::ow_samtools(argv.len() as c_int, argv.as_mut_ptr()) }
}

const REFS: [(&str, u64); 4] = [
    ("chr1", 60_000),
    ("chr2", 30_000),
    ("chrM", 16_569),
    ("unused", 1_000),
];

fn bases(rng: &mut Rng, n: usize) -> (String, String) {
    let seq = (0..n).map(|_| *rng.pick(b"ACGT") as char).collect();
    let qual = (0..n)
        .map(|_| (b'!' + rng.range(0, 40) as u8) as char)
        .collect();
    (seq, qual)
}

/// SAM records in no order over `refs`: piles of reads at the same
/// position on both strands, unmapped reads placed at their mate's
/// position and unplaced ones.
fn random_records(rng: &mut Rng, n: usize, refs: &[(&str, u64)], out: &mut String) {
    let placed = refs.len().saturating_sub(1).max(1);
    for i in 0..n {
        let l = rng.range(20, 120) as usize;
        let (seq, qual) = bases(rng, l);
        let kind = rng.range(0, 100);
        if kind < 8 || refs.is_empty() {
            let flag = *rng.pick(&[4, 4 | 16, 77, 141]);
            out.push_str(&format!(
                "u{i}\t{flag}\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}\n"
            ));
            continue;
        }
        let (name, len) = refs[rng.range(0, placed as u64) as usize];
        let pos = if rng.chance(30) {
            rng.range(1, 40)
        } else {
            rng.range(1, len - 30_000.min(len / 2))
        };
        if kind < 20 {
            let flag = 4 | 1 | *rng.pick(&[0, 16, 32, 48]) | *rng.pick(&[64, 128]);
            out.push_str(&format!(
                "p{i}\t{flag}\t{name}\t{pos}\t0\t*\t=\t{pos}\t0\t{seq}\t{qual}\n"
            ));
            continue;
        }
        let flag = *rng.pick(&[0, 16, 99, 147, 83, 163, 256, 272, 1024, 2048, 2064, 528]);
        let a = rng.range(1, l as u64 - 1);
        let cigar = match rng.range(0, 4) {
            0 => format!("{l}M"),
            1 => format!("{a}M{}D{}M", rng.range(1, 50), l as u64 - a),
            2 => format!("{a}S{}M", l as u64 - a),
            _ => format!("{a}M{}N{}M", rng.range(100, 25_000), l as u64 - a),
        };
        out.push_str(&format!(
            "m{i}\t{flag}\t{name}\t{pos}\t{}\t{cigar}\t*\t0\t0\t{seq}\t{qual}\tNM:i:{}\n",
            rng.range(0, 61),
            rng.range(0, 9)
        ));
    }
}

fn sam_header(refs: &[(&str, u64)]) -> String {
    let mut s = String::from("@HD\tVN:1.6\tSO:unsorted\n");
    for (name, len) in refs {
        s.push_str(&format!("@SQ\tSN:{name}\tLN:{len}\n"));
    }
    s.push_str("@RG\tID:rg1\tSM:x\n@PG\tID:aligner\tPN:aligner\n@CO\tunsorted reads\n");
    s
}

/// An unsorted BAM of `sam` text at `dst`.
fn bam_of(dir: &Path, sam: &str, dst: &Path) {
    let src = dir.join("unsorted.sam");
    std::fs::write(&src, sam).unwrap();
    let mut c = SamToBam::new(&src, dst);
    c.level = 1;
    convert(Conversion::SamToBam(c), Budget::new(2, 1 << 30), None).unwrap();
}

/// A BAM taken apart: header text, references, records with their
/// `block_size`.
#[derive(Debug, PartialEq)]
struct Bam {
    text: Vec<u8>,
    refs: Vec<(Vec<u8>, u32)>,
    records: Vec<Vec<u8>>,
}

fn u32_at(b: &[u8], at: usize) -> u32 {
    u32::from_le_bytes(b[at..at + 4].try_into().unwrap())
}

fn read_bam(path: &Path) -> Bam {
    let raw = inflate(path);
    assert_eq!(&raw[..4], b"BAM\x01");
    let l_text = u32_at(&raw, 4) as usize;
    let text = raw[8..8 + l_text].to_vec();
    let mut at = 8 + l_text;
    let n_ref = u32_at(&raw, at);
    at += 4;
    let mut refs = Vec::new();
    for _ in 0..n_ref {
        let l = u32_at(&raw, at) as usize;
        let name = raw[at + 4..at + 4 + l - 1].to_vec();
        refs.push((name, u32_at(&raw, at + 4 + l)));
        at += 8 + l;
    }
    let mut records = Vec::new();
    while at < raw.len() {
        let n = 4 + u32_at(&raw, at) as usize;
        records.push(raw[at..at + n].to_vec());
        at += n;
    }
    Bam {
        text,
        refs,
        records,
    }
}

fn write_bam(path: &Path, bam: &Bam) {
    let mut head = b"BAM\x01".to_vec();
    head.extend((bam.text.len() as u32).to_le_bytes());
    head.extend(&bam.text);
    head.extend((bam.refs.len() as u32).to_le_bytes());
    for (name, len) in &bam.refs {
        head.extend((name.len() as u32 + 1).to_le_bytes());
        head.extend(name);
        head.push(0);
        head.extend(len.to_le_bytes());
    }
    let mut w = omicsweave_bgzf::Writer::new(
        Vec::new(),
        1,
        Box::new(omicsweave_core::exec::Inline::new()),
    )
    .unwrap();
    w.write(&head).unwrap();
    w.flush_block().unwrap();
    for r in &bam.records {
        w.write_record(r).unwrap();
    }
    std::fs::write(path, w.finish().unwrap().out).unwrap();
}

fn name_of(record: &[u8]) -> String {
    let l = (u32_at(record, 12) & 0xff) as usize;
    String::from_utf8_lossy(&record[36..36 + l - 1]).into_owned()
}

/// Uncompressed sizes of the BGZF blocks of a file.
fn block_sizes(path: &Path) -> Vec<u32> {
    let raw = std::fs::read(path).unwrap();
    let mut at = 0;
    let mut out = Vec::new();
    while at < raw.len() {
        let size = u16::from_le_bytes([raw[at + 16], raw[at + 17]]) as usize + 1;
        out.push(u32_at(&raw, at + size - 4));
        at += size;
    }
    out
}

/// The same header, records and blocks; and, with `bytes`, the same file.
fn assert_same(ours: &Path, theirs: &Path, bytes: bool) {
    let (a, b) = (read_bam(ours), read_bam(theirs));
    assert_eq!(
        String::from_utf8_lossy(&a.text),
        String::from_utf8_lossy(&b.text)
    );
    assert_eq!(a.refs, b.refs);
    assert_eq!(a.records.len(), b.records.len());
    if let Some(i) = (0..a.records.len()).find(|&i| a.records[i] != b.records[i]) {
        panic!(
            "record {i} differs: ours {} theirs {}",
            name_of(&a.records[i]),
            name_of(&b.records[i])
        );
    }
    assert_eq!(block_sizes(ours), block_sizes(theirs));
    if bytes {
        assert!(std::fs::read(ours).unwrap() == std::fs::read(theirs).unwrap());
    }
}

fn sort(src: &Path, dst: &Path, level: i32, index: BamIndex, budget: Budget) -> u64 {
    let mut c = BamSort::new(src, dst);
    c.level = level;
    c.index = index;
    convert(Conversion::BamSort(c), budget, None)
        .unwrap()
        .spilled
}

/// Roomy and tight budgets on one and several threads; the tight ones
/// spill the reads of [`random_records`] many times over.
fn budgets() -> Vec<(Budget, bool)> {
    vec![
        (Budget::new(1, 1 << 30), false),
        (Budget::new(5, 1 << 30), false),
        (Budget::new(1, 20 << 20), true),
        (Budget::new(3, 24 << 20), true),
    ]
}

#[test]
fn sorts_as_samtools_byte_for_byte() {
    let dir = tempfile::tempdir().unwrap();
    for seed in 0..2u64 {
        let mut rng = Rng::new(seed);
        let mut sam = sam_header(&REFS);
        random_records(&mut rng, 80_000, &REFS, &mut sam);
        let src = dir.path().join("in.bam");
        bam_of(dir.path(), &sam, &src);
        let theirs = dir.path().join("theirs.bam");
        assert_eq!(samtools_sort(&src, &theirs, &["-l", "0", "-m", "1M"]), 0);
        let mut first: Option<Vec<u8>> = None;
        for (budget, spills) in budgets() {
            let ours = dir.path().join("ours.bam");
            let spilled = sort(&src, &ours, 0, BamIndex::None, budget.clone());
            assert_eq!(spilled > 0, spills, "{budget:?}");
            assert_same(&ours, &theirs, true);
            let bytes = std::fs::read(&ours).unwrap();
            match &first {
                None => first = Some(bytes),
                Some(f) => assert!(*f == bytes, "{budget:?}"),
            }
        }
    }
    assert!(
        std::fs::read_dir(dir.path()).unwrap().all(|e| !e
            .unwrap()
            .file_name()
            .to_string_lossy()
            .contains("spill")),
        "spill directories are removed"
    );
}

#[test]
fn compressed_output_has_samtools_blocks_and_records() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(7);
    let mut sam = sam_header(&REFS);
    random_records(&mut rng, 50_000, &REFS, &mut sam);
    let src = dir.path().join("in.bam");
    bam_of(dir.path(), &sam, &src);
    let theirs = dir.path().join("theirs.bam");
    // Several threads and in-memory blocks merged with spilled files.
    assert_eq!(
        samtools_sort(&src, &theirs, &["-l", "6", "-m", "1M", "-@", "3"]),
        0
    );
    for (budget, _) in budgets() {
        let ours = dir.path().join("ours.bam");
        sort(&src, &ours, 6, BamIndex::None, budget);
        assert_same(&ours, &theirs, false);
    }
}

#[test]
fn indexes_match_htslib() {
    let dir = tempfile::tempdir().unwrap();
    let long = [("chr1", 60_000), ("big", 900_000_000), ("chrM", 16_569)];
    for (refs, index, min_shift, suffix) in [
        (&REFS[..], BamIndex::Auto, 0, "bai"),
        (&REFS[..], BamIndex::Csi { min_shift: 14 }, 14, "csi"),
        (&long[..], BamIndex::Auto, 14, "csi"),
    ] {
        let mut rng = Rng::new(3);
        let mut sam = sam_header(refs);
        random_records(&mut rng, 30_000, refs, &mut sam);
        let src = dir.path().join("in.bam");
        bam_of(dir.path(), &sam, &src);
        let theirs = dir.path().join("theirs.bam");
        assert_eq!(samtools_sort(&src, &theirs, &["-l", "0"]), 0);
        for (budget, _) in budgets() {
            let ours = dir.path().join("ours.bam");
            sort(&src, &ours, 0, index, budget);
            assert_same(&ours, &theirs, true);
            let built = dir.path().join("built.idx");
            hts::bam_index(&ours, &built, min_shift);
            assert_eq!(
                hts::parse_index(&dir.path().join(format!("ours.bam.{suffix}"))),
                hts::parse_index(&built)
            );
        }
    }
}

/// Sorts `bam` with both and compares the files byte for byte.
fn check_against_samtools(dir: &Path, bam: &Bam) {
    let src = dir.join("case.bam");
    write_bam(&src, bam);
    let theirs = dir.join("case-theirs.bam");
    assert_eq!(samtools_sort(&src, &theirs, &["-l", "0"]), 0);
    for (budget, _) in &budgets()[..2] {
        let ours = dir.join("case-ours.bam");
        sort(&src, &ours, 0, BamIndex::Auto, budget.clone());
        assert_same(&ours, &theirs, true);
    }
}

fn small_bam(dir: &Path, seed: u64, refs: &[(&str, u64)], n: usize) -> Bam {
    let mut rng = Rng::new(seed);
    let mut sam = sam_header(refs);
    random_records(&mut rng, n, refs, &mut sam);
    let src = dir.join("small.bam");
    bam_of(dir, &sam, &src);
    read_bam(&src)
}

#[test]
fn headers_are_rewritten_as_samtools_rewrites_them() {
    let dir = tempfile::tempdir().unwrap();
    let base = small_bam(dir.path(), 11, &REFS, 500);
    let texts: [&[u8]; 10] = [
        // No @HD.
        b"@SQ\tSN:chr1\tLN:60000\n@CO\tno HD\n",
        // @HD after other lines, with a group order and a sub-sort.
        b"@CO\tfirst\n@HD\tVN:1.4\tGO:query\tSO:queryname\tSS:queryname:natural\n@RG\tID:a\n",
        // Two @HD lines.
        b"@HD\tVN:1.6\n@SQ\tSN:chr2\tLN:30000\n@HD\tVN:1.0\tSO:x\n",
        // No text: every reference gets an @SQ line.
        b"",
        // NUL padding, and a last line without a newline.
        b"@HD\tVN:1.6\tSO:unsorted\n@PG\tID:p\tPN:p\0\0\0\0",
        b"@HD\tVN:1.6\n@CO\tlast line\twithout newline",
        // A NUL inside the text, which htslib reads as a newline.
        b"@HD\tVN:1.6\0@SQ\tSN:chr1\tLN:60000\n",
        // A length differing from the reference list, and alternative names.
        b"@SQ\tSN:chrM\tLN:5\tAN:MT,M\n@SQ\tSN:chr1\tLN:60000\tAN:1\n",
        // An @SQ line for a reference the list lacks.
        b"@SQ\tSN:extra\tLN:777\n@SQ\tSN:chr2\tLN:30000\n",
        // Carriage returns stay in the fields.
        b"@HD\tVN:1.6\r\n@CO\tdos\r\n",
    ];
    for text in texts {
        let bam = Bam {
            text: text.to_vec(),
            refs: base.refs.clone(),
            records: base.records.clone(),
        };
        check_against_samtools(dir.path(), &bam);
    }
}

#[test]
fn bins_are_recomputed_as_htslib_reads_them() {
    let dir = tempfile::tempdir().unwrap();
    let mut bam = small_bam(dir.path(), 12, &REFS, 2_000);
    for r in &mut bam.records {
        r[14..16].copy_from_slice(&0u16.to_le_bytes());
    }
    check_against_samtools(dir.path(), &bam);
}

#[test]
fn empty_and_unmapped_files() {
    let dir = tempfile::tempdir().unwrap();
    // No records at all.
    let empty = small_bam(dir.path(), 1, &REFS, 0);
    assert!(empty.records.is_empty());
    check_against_samtools(dir.path(), &empty);
    // Only unplaced reads, no references.
    let unplaced = small_bam(dir.path(), 2, &[], 3_000);
    assert!(unplaced.refs.is_empty());
    check_against_samtools(dir.path(), &unplaced);
    // Only unmapped reads, some placed with their mates.
    let mut placed = small_bam(dir.path(), 3, &REFS, 6_000);
    placed
        .records
        .retain(|r| (u32_at(r, 16) >> 16) as u16 & omicsweave_bam::flag::UNMAPPED != 0);
    assert!(placed.records.iter().any(|r| u32_at(r, 4) != u32::MAX));
    check_against_samtools(dir.path(), &placed);
}

#[test]
fn a_bam_converts_to_a_bam_only_by_sorting() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("in.bam");
    bam_of(dir.path(), &sam_header(&REFS), &src);
    let dst = dir.path().join("out.bam");
    assert!(Conversion::routes().any(|r| r == (Format::Bam, Format::Bam)));
    let mut c = Conversion::infer(&src, &dst).unwrap();
    let err = c.configure(&Options::default()).unwrap_err();
    assert!(err.to_string().contains("pass sort"), "{err}");
    let err = c
        .configure(&Options {
            index: Some(Index::Auto),
            ..Options::default()
        })
        .unwrap_err();
    assert!(
        err.to_string().contains("index needs a coordinate-sorted"),
        "{err}"
    );
    let sorting = Options {
        sort: true,
        ..Options::default()
    };
    let mut sam = Conversion::infer(dir.path().join("in.sam"), &dst).unwrap();
    let err = sam.configure(&sorting).unwrap_err();
    assert!(
        err.to_string().contains("does not apply to SAM to BAM"),
        "{err}"
    );
    c.configure(&Options {
        index: Some(Index::Bai),
        level: Some(0),
        ..sorting
    })
    .unwrap();
    let err = plan(c.clone(), Budget::new(4, 8 << 20)).unwrap_err();
    assert!(
        err.to_string().contains("sorting a BAM needs at least"),
        "{err}"
    );
    let report = convert(c, Budget::new(2, 64 << 20), None).unwrap();
    assert_eq!(
        report.files,
        vec![
            (Format::Bam, dst.clone()),
            (Format::Bai, dir.path().join("out.bam.bai"))
        ]
    );
    assert_eq!(report.records, 0);
}

#[test]
fn malformed_headers_fail_as_in_samtools() {
    let dir = tempfile::tempdir().unwrap();
    let base = small_bam(dir.path(), 5, &REFS, 10);
    for text in [
        &b"@SQ\tSN:chr1\n"[..],
        b"@HD\tVN:1.6\t\n",
        b"@RG\tSM:x\n",
        b"@HD\tVN:1.6\n\n",
        b"@HD\tVN:1.6\nab\n",
    ] {
        let src = dir.path().join("bad.bam");
        write_bam(
            &src,
            &Bam {
                text: text.to_vec(),
                refs: base.refs.clone(),
                records: base.records.clone(),
            },
        );
        assert_ne!(samtools_sort(&src, &dir.path().join("t.bam"), &[]), 0);
        let ours = dir.path().join("o.bam");
        let mut c = BamSort::new(&src, &ours);
        c.level = 0;
        assert!(convert(Conversion::BamSort(c), Budget::new(1, 1 << 30), None).is_err());
        assert!(!ours.exists());
    }
}
