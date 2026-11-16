use crate::support::modifications::{
    as_hts, htslib_mods, our_mods, random_sam, spec_lines, types_as_hts, write_bam,
};

#[test]
fn records_decode_as_htslib_decodes_them() {
    let dir = tempfile::tempdir().unwrap();
    for (seed, lines) in [
        (0, spec_lines()),
        (1, random_sam(11, 3000)),
        (2, random_sam(29, 3000)),
    ] {
        let bam = write_bam(dir.path(), &format!("s{seed}"), &lines);
        let ours = our_mods(&bam);
        let theirs = htslib_mods(&bam);
        assert_eq!(ours.len(), theirs.len());
        assert_eq!(ours.len(), lines.len());
        let (mut good, mut bad) = (0, 0);
        for (i, (o, t)) in ours.iter().zip(&theirs).enumerate() {
            let o = o
                .as_ref()
                .map(|(rows, types)| (as_hts(rows), types_as_hts(types)));
            if t.is_some() {
                good += 1;
            } else {
                bad += 1;
            }
            assert_eq!(&o, t, "record {i}: {}", lines[i]);
        }
        if seed > 0 {
            assert!(good > 1000 && bad > 50, "{good} decoded, {bad} rejected");
        }
    }
}
