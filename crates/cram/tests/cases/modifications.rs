//! Base modifications of CRAM records: a region query returns what the
//! BAM decoder makes of the records htslib reads from the same CRAM, with
//! bases rebuilt from the reference, for CRAM 3.0 and 3.1.

use std::path::Path;

use crate::support::{Rng, crai_of, cram_of, fasta_reader, read_all, reference, write_fasta};
use omicsweave_bam::mods::{ModDecoder, ModSeq, decode_record};
use omicsweave_core::{
    AlignmentFields, AlignmentRead, AlignmentSize, AlignmentsOut, FileSource, ModificationColumn,
    ModificationEntry, Region,
};
use omicsweave_cram::{CramIndex, CramOptions, CramReader};

const LEN: u64 = 200_000;

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        x => x,
    }
}

/// Reads drawn from `genome` with mismatches, an insertion now and then,
/// soft clips and well-formed `MM` / `ML` tags (one or two entries, both
/// strands, explicit and implicit, ChEBI codes).
fn sam(rng: &mut Rng, genome: &[u8]) -> String {
    let mut text = format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{LEN}\n");
    let mut pos = 1u64;
    for i in 0..3000 {
        pos += rng.range(0, 60);
        let len = rng.range(30, 160) as usize;
        if pos as usize + len >= genome.len() {
            break;
        }
        let mut seq: Vec<u8> = genome[pos as usize - 1..pos as usize - 1 + len]
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();
        for b in seq.iter_mut() {
            if rng.chance(3) {
                *b = b"ACGT"[rng.range(0, 4) as usize];
            }
        }
        let clip = if rng.chance(30) {
            rng.range(1, 10) as usize
        } else {
            0
        };
        let cigar = if clip > 0 {
            format!("{clip}S{}M", len - clip)
        } else if rng.chance(20) {
            format!("10M2I{}M", len - 12)
        } else {
            format!("{len}M")
        };
        let reverse = rng.chance(50);
        let original: Vec<u8> = if reverse {
            seq.iter().rev().map(|&b| complement(b)).collect()
        } else {
            seq.clone()
        };
        let mut mm = String::new();
        let mut n_ml = 0;
        for (base, codes, n_codes) in [
            (b'C', ["m", "mh", "m?", "h."][rng.range(0, 4) as usize], 0),
            (b'A', ["a", "17596"][rng.range(0, 2) as usize], 1),
            (b'G', ["m"][0], 1),
        ] {
            if base != b'C' && rng.chance(50) {
                continue;
            }
            let n_codes = if n_codes == 0 {
                codes.chars().filter(char::is_ascii_alphabetic).count()
            } else {
                n_codes
            };
            let strand = if base == b'G' { '-' } else { '+' };
            mm.push_str(&format!("{}{strand}{codes}", base as char));
            let mut left = original.iter().filter(|&&b| b == base).count();
            while left > 0 && rng.chance(85) {
                let d = rng.range(0, (left as u64).min(4)) as usize;
                mm.push_str(&format!(",{d}"));
                n_ml += n_codes;
                left -= d + 1;
            }
            mm.push(';');
        }
        let mut tags = format!("\tMM:Z:{mm}");
        if n_ml > 0 && rng.chance(90) {
            let v: Vec<String> = (0..n_ml).map(|_| rng.range(0, 256).to_string()).collect();
            tags.push_str(&format!("\tML:B:C,{}", v.join(",")));
        }
        let flag = if reverse { 16 } else { 0 };
        text.push_str(&format!(
            "r{i}\t{flag}\tchr1\t{pos}\t60\t{cigar}\t*\t0\t0\t{}\t{}{tags}\n",
            String::from_utf8(seq).unwrap(),
            "I".repeat(len)
        ));
    }
    text
}

fn modifications(
    r: &mut impl AlignmentRead,
    region: Region<'_>,
) -> (Vec<i64>, Vec<ModificationEntry>) {
    let size = r.query(region, AlignmentFields::MODIFICATIONS).unwrap();
    write_modifications(r, size)
}

fn write_modifications(
    r: &impl AlignmentRead,
    size: AlignmentSize,
) -> (Vec<i64>, Vec<ModificationEntry>) {
    let mut entries = vec![Default::default(); size.records];
    let mut offsets = vec![0i64; size.records + 1];
    let mut rows = vec![ModificationEntry::default(); size.modifications];
    r.write(AlignmentsOut {
        entries: &mut entries,
        names: None,
        seqs: None,
        quals: None,
        cigars: None,
        aux: None,
        modifications: Some(ModificationColumn {
            offsets: &mut offsets,
            entries: &mut rows,
        }),
    })
    .unwrap();
    (offsets, rows)
}

#[test]
fn cram_modifications_match_the_records_htslib_reads() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(5);
    let refs = vec![("chr1".to_string(), LEN)];
    let genome = reference(&mut rng, &refs);
    let fasta = dir.path().join("ref.fa");
    write_fasta(&fasta, &refs, &genome);
    let sam_path = dir.path().join("in.sam");
    std::fs::write(&sam_path, sam(&mut rng, &genome[0])).unwrap();
    for version in ["version=3.0", "version=3.1"] {
        let cram = dir.path().join(format!("{version}.cram"));
        cram_of(&sam_path, &cram, Some(&fasta), &[version]);
        let crai = dir.path().join(format!("{version}.cram.crai"));
        crai_of(&cram, &crai);
        let records = read_all(&cram, Some(&fasta));
        let mut decoder = ModDecoder::new();
        let mut want = Vec::new();
        let mut want_offsets = vec![0i64];
        for r in &records {
            decode_record(
                &mut decoder,
                &[&r.aux],
                ModSeq::Ascii(&r.seq),
                r.flag,
                r.pos,
                r.cigar.iter().map(|&c| ((c & 0xf) as u8, c >> 4)),
                &r.name,
                &[],
                &mut want,
            )
            .unwrap();
            want_offsets.push(want.len() as i64);
        }
        assert!(want.len() > 20_000, "{}", want.len());
        assert!(want.iter().any(|m| m.code == -17596));
        assert!(want.iter().any(|m| m.ref_pos < 0));

        let header = omicsweave_cram::CramHeader::read(&FileSource::open(&cram).unwrap()).unwrap();
        let index = CramIndex::read_file(&crai, header.references().len()).unwrap();
        let mut reader = CramReader::open(
            FileSource::open(&cram).unwrap(),
            index,
            Some(fasta_reader(&fasta)),
            CramOptions::default(),
        )
        .unwrap();
        let (offsets, rows) = modifications(&mut reader, Region::new("chr1", 0, LEN as u32));
        assert_eq!(offsets, want_offsets, "{version}");
        assert_eq!(rows, want, "{version}");

        for codes in [
            &[b'h' as i32][..],
            &[b'h' as i32, -17596][..],
            &[i32::MAX][..],
        ] {
            let size = reader
                .query_filtered(
                    Region::new("chr1", 0, LEN as u32),
                    AlignmentFields::MODIFICATIONS,
                    codes,
                )
                .unwrap();
            let mut expected_rows = Vec::new();
            let mut expected_offsets = vec![0];
            for record in want_offsets.windows(2) {
                expected_rows.extend(
                    want[record[0] as usize..record[1] as usize]
                        .iter()
                        .copied()
                        .filter(|row| codes.contains(&row.code)),
                );
                expected_offsets.push(expected_rows.len() as i64);
            }
            let (offsets, rows) = write_modifications(&reader, size);
            assert_eq!(size.records, records.len(), "{version}, codes={codes:?}");
            assert_eq!(
                size.modifications,
                expected_rows.len(),
                "{version}, codes={codes:?}"
            );
            assert_eq!(offsets, expected_offsets, "{version}, codes={codes:?}");
            assert_eq!(rows, expected_rows, "{version}, codes={codes:?}");
        }
    }
}

#[test]
fn modifications_need_the_reference_bases() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(9);
    let refs = vec![("chr1".to_string(), LEN)];
    let genome = reference(&mut rng, &refs);
    let fasta = dir.path().join("ref.fa");
    write_fasta(&fasta, &refs, &genome);
    let sam_path = dir.path().join("in.sam");
    std::fs::write(&sam_path, sam(&mut rng, &genome[0])).unwrap();
    let cram = dir.path().join("in.cram");
    cram_of(&sam_path, &cram, Some(&fasta), &[]);
    let crai = dir.path().join("in.cram.crai");
    crai_of(&cram, &crai);
    let header = omicsweave_cram::CramHeader::read(&FileSource::open(&cram).unwrap()).unwrap();
    let open = |path: &Path| {
        CramReader::<FileSource, omicsweave_fasta::FastaReader<FileSource>>::open(
            FileSource::open(path).unwrap(),
            CramIndex::read_file(&crai, header.references().len()).unwrap(),
            None,
            CramOptions::default(),
        )
        .unwrap()
    };
    let mut reader = open(&cram);
    let region = Region::new("chr1", 0, 10_000);
    assert!(
        reader
            .query(region, AlignmentFields::MODIFICATIONS)
            .is_err()
    );
    assert!(reader.query(region, AlignmentFields::CIGAR).is_ok());
}
