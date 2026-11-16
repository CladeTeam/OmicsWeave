use crate::support::modifications::{
    HtsMod, as_hts, htslib_mods, random_sam, spec_lines, write_bam,
};
use omicsweave_bam::{BamReader, BinIndex};
use omicsweave_core::{
    AlignmentFields, AlignmentRead, AlignmentsOut, FileSource, ModificationColumn,
    ModificationEntry, Region,
};

/// Reference positions of every query position, by walking `cigar`.
fn ref_of_query(cigar: &[u32], pos: i64, l_seq: usize) -> Vec<i64> {
    let mut out = vec![-1; l_seq];
    let (mut q, mut r) = (0usize, pos);
    for &c in cigar {
        let (op, len) = (c & 0xf, (c >> 4) as usize);
        match op {
            0 | 7 | 8 => {
                // A record without SEQ still has a CIGAR.
                for k in 0..len {
                    if let Some(o) = out.get_mut(q + k) {
                        *o = r + k as i64;
                    }
                }
                q += len;
                r += len as i64;
            }
            1 | 4 => q += len,
            2 | 3 => r += len as i64,
            _ => {}
        }
    }
    out
}

#[test]
fn region_queries_place_modifications_on_the_reference() {
    let dir = tempfile::tempdir().unwrap();
    let mut lines = spec_lines();
    lines.extend(random_sam(47, 4000));
    lines.sort_by_key(|l| l.split('\t').nth(3).unwrap().parse::<u64>().unwrap());
    let probe = write_bam(dir.path(), "probe", &lines);
    // Keep the records htslib decodes: one bad record fails a query.
    let kept: Vec<String> = htslib_mods(&probe)
        .iter()
        .zip(&lines)
        .filter(|(m, _)| m.is_some())
        .map(|(_, l)| l.clone())
        .collect();
    let bam = write_bam(dir.path(), "kept", &kept);
    let bai = dir.path().join("kept.bam.bai");
    crate::support::index(&bam, &bai, 0);
    let expected = htslib_mods(&bam);
    let records = crate::support::read_all(&bam);
    let mut reader = BamReader::open(
        FileSource::open(&bam).unwrap(),
        BinIndex::read_file(&bai).unwrap(),
    )
    .unwrap();
    for codes in [vec![], vec![b'm' as i32], vec![b'h' as i32, -76792]] {
        let region = Region::new("chr1", 0, 1_000_000);
        let size = reader
            .query_filtered(region, AlignmentFields::MODIFICATIONS, &codes)
            .unwrap();
        assert_eq!(size.records, records.len());
        let mut entries = vec![Default::default(); size.records];
        let mut offsets = vec![-777i64; size.records + 3];
        let sentinel = ModificationEntry {
            qpos: -777,
            ..Default::default()
        };
        let mut rows = vec![sentinel; size.modifications + 2];
        reader
            .write(AlignmentsOut {
                entries: &mut entries,
                names: None,
                seqs: None,
                quals: None,
                cigars: None,
                aux: None,
                modifications: Some(ModificationColumn {
                    offsets: &mut offsets[1..size.records + 2],
                    entries: &mut rows[1..size.modifications + 1],
                }),
            })
            .unwrap();
        assert_eq!((offsets[0], offsets[size.records + 2]), (-777, -777));
        assert_eq!(
            (rows[0], rows[size.modifications + 1]),
            (sentinel, sentinel)
        );
        let offsets = &offsets[1..size.records + 2];
        let rows = &rows[1..size.modifications + 1];
        assert_eq!(
            (offsets[0], offsets[size.records]),
            (0, size.modifications as i64)
        );
        let mut total = 0;
        for (i, (rec, theirs)) in records.iter().zip(&expected).enumerate() {
            let (mods, _) = theirs.as_ref().unwrap();
            let got = &rows[offsets[i] as usize..offsets[i + 1] as usize];
            let want: Vec<HtsMod> = mods
                .iter()
                .copied()
                .filter(|m| codes.is_empty() || codes.contains(&m.1))
                .collect();
            assert_eq!(as_hts(got), want, "record {i}");
            let unmapped = rec.flag & 0x4 != 0;
            let refs = ref_of_query(&rec.cigar, rec.pos, rec.seq.len());
            for m in got {
                let want = if unmapped { -1 } else { refs[m.qpos as usize] };
                assert_eq!(m.ref_pos, want, "record {i} at {}", m.qpos);
            }
            total += got.len();
        }
        assert_eq!(total, size.modifications);
        if codes.is_empty() {
            assert!(total > 10_000, "{total}");
        }
    }
}
