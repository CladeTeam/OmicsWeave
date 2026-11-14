//! SAM → BAM: the uncompressed stream htslib writes, the same file for
//! every budget, and BAI / CSI indexes equal to the ones htslib builds.

mod common;

use common::{Rng, STORAGES, budgets, hts, inflate, write_input};
use omicsweave_convert::{BamIndex, Conversion, SamToBam, convert};

fn random_sam(rng: &mut Rng, n: usize, long_ref: bool) -> Vec<u8> {
    let refs: Vec<(&str, u64)> = if long_ref {
        vec![("chr1", 248_956_422), ("big", 900_000_000)]
    } else {
        vec![
            ("chr1", 248_956_422),
            ("chr2", 242_193_529),
            ("chrM", 16_569),
        ]
    };
    let mut s = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for (name, len) in &refs {
        s.push_str(&format!("@SQ\tSN:{name}\tLN:{len}\n"));
    }
    s.push_str("@RG\tID:rg1\tSM:x\n");
    let per = n / refs.len();
    let bases = b"ACGT";
    for (name, len) in &refs {
        let mut pos = 1u64;
        for i in 0..per {
            // Clustered reads, as real data, so bins fill and merge.
            pos += if rng.chance(95) {
                rng.range(0, 40)
            } else {
                rng.range(1, 2_000_000)
            };
            if pos + 200 >= *len {
                break;
            }
            let l = rng.range(30, 160) as usize;
            let seq: String = (0..l).map(|_| *rng.pick(bases) as char).collect();
            let qual: String = (0..l)
                .map(|_| (b'!' + rng.range(0, 40) as u8) as char)
                .collect();
            let (cigar, flag) = if rng.chance(3) {
                ("*".to_string(), 4)
            } else if rng.chance(20) {
                (
                    format!("{}M{}D{}M", l / 2, rng.range(1, 500), l - l / 2),
                    99,
                )
            } else {
                (format!("{l}M"), 0)
            };
            s.push_str(&format!(
                "r{i}\t{flag}\t{name}\t{pos}\t{}\t{cigar}\t=\t{pos}\t0\t{seq}\t{qual}\tRG:Z:rg1\tNM:i:{}\n",
                rng.range(0, 61),
                rng.range(0, 300)
            ));
        }
    }
    for i in 0..rng.range(0, 50) {
        s.push_str(&format!("u{i}\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n"));
    }
    s.into_bytes()
}

#[test]
fn bam_matches_htslib_with_its_indexes() {
    let dir = tempfile::tempdir().unwrap();
    for (seed, long_ref) in [(0u64, false), (1, false), (2, true)] {
        let mut rng = Rng::new(seed);
        let text = random_sam(&mut rng, 60_000, long_ref);
        let plain = dir.path().join("ref.sam");
        std::fs::write(&plain, &text).unwrap();
        let reference = dir.path().join("ref.bam");
        hts::bam_of(&plain, &reference);
        let want = inflate(&reference);
        let mut first: Option<Vec<u8>> = None;
        for storage in STORAGES {
            let src = dir.path().join("in.sam");
            write_input(&src, &text, storage);
            for budget in budgets() {
                let dst = dir.path().join("out.bam");
                let mut conversion = SamToBam::new(&src, &dst);
                conversion.index = BamIndex::Auto;
                convert(Conversion::SamToBam(conversion), budget.clone(), None).unwrap();
                assert_eq!(inflate(&dst), want, "{storage:?} {budget:?}");
                let (ours, min_shift) = if long_ref {
                    (dir.path().join("out.bam.csi"), 14)
                } else {
                    (dir.path().join("out.bam.bai"), 0)
                };
                let theirs = dir.path().join("theirs.idx");
                hts::bam_index(&dst, &theirs, min_shift);
                assert_eq!(
                    hts::parse_index(&ours),
                    hts::parse_index(&theirs),
                    "{storage:?} {budget:?}"
                );
                let bytes = std::fs::read(&dst).unwrap();
                match &first {
                    None => first = Some(bytes),
                    Some(f) => assert_eq!(f, &bytes, "{storage:?} {budget:?}"),
                }
            }
        }
    }
}

#[test]
fn an_index_needs_sorted_input() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("a.sam");
    std::fs::write(
        &src,
        b"@SQ\tSN:c\tLN:1000\nr1\t0\tc\t500\t9\t4M\t*\t0\t0\tACGT\t*\nr2\t0\tc\t10\t9\t4M\t*\t0\t0\tACGT\t*\n",
    )
    .unwrap();
    let dst = dir.path().join("a.bam");
    let mut conversion = SamToBam::new(&src, &dst);
    conversion.index = BamIndex::Bai;
    let err = convert(Conversion::SamToBam(conversion), budgets()[0].clone(), None).unwrap_err();
    assert!(err.to_string().contains("unsorted"), "{err}");
    assert!(!dst.exists());
    // Without an index the order is kept.
    convert(
        Conversion::SamToBam(SamToBam::new(&src, &dst)),
        budgets()[0].clone(),
        None,
    )
    .unwrap();
    assert!(dst.exists());
}
