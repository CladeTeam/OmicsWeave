//! VCF → BCF: the uncompressed stream htslib writes, the same file for
//! every budget, and a CSI index equal to the one htslib builds for it.

mod common;

use common::{Rng, STORAGES, budgets, hts, inflate, write_input};
use omicsweave_convert::{Conversion, VcfToBcf, convert};

fn random_vcf(rng: &mut Rng, n: usize) -> Vec<u8> {
    let contigs = [
        ("chr1", 248_956_422u64),
        ("chr2", 242_193_529),
        ("chrM", 16_569),
    ];
    let samples = ["s1", "s2", "s3"];
    let mut v = String::from("##fileformat=VCFv4.2\n");
    for (c, l) in contigs {
        v.push_str(&format!("##contig=<ID={c},length={l}>\n"));
    }
    v.push_str("##FILTER=<ID=q10,Description=\"low\">\n");
    v.push_str("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n");
    v.push_str("##INFO=<ID=AF,Number=A,Type=Float,Description=\"a\">\n");
    v.push_str("##INFO=<ID=DB,Number=0,Type=Flag,Description=\"b\">\n");
    v.push_str("##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n");
    v.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n");
    v.push_str("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"ad\">\n");
    v.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for s in samples {
        v.push('\t');
        v.push_str(s);
    }
    v.push('\n');
    let per = n / contigs.len();
    for (ci, (c, len)) in contigs.iter().enumerate() {
        let mut pos = 1u64;
        for i in 0..per {
            pos += rng.range(1, (len / per as u64).max(2) * 2);
            if pos >= *len {
                break;
            }
            let alts = rng.range(1, 3) as usize;
            let alt: Vec<&str> = (0..alts)
                .map(|k| ["A", "CT", "<DEL>"][(k + ci) % 3])
                .collect();
            let mut info = format!("DP={}", rng.range(0, 70_000));
            info.push_str(&format!(";AF={}", vec!["0.5"; alts].join(",")));
            if rng.chance(30) {
                info.push_str(";DB");
            }
            if alt.contains(&"<DEL>") {
                info.push_str(&format!(";END={}", (pos + rng.range(1, 3000)).min(*len)));
            }
            let qual = if rng.chance(10) {
                ".".into()
            } else {
                format!("{}", rng.range(0, 1000))
            };
            let filter = *rng.pick(&["PASS", "q10", "."]);
            v.push_str(&format!(
                "{c}\t{pos}\trs{i}\tG\t{}\t{qual}\t{filter}\t{info}\tGT:AD",
                alt.join(",")
            ));
            for _ in samples {
                let gt = *rng.pick(&["0/1", "1|1", "./.", "0", "1/2"]);
                let ad: Vec<String> = (0..=alts).map(|_| rng.range(0, 300).to_string()).collect();
                v.push_str(&format!("\t{gt}:{}", ad.join(",")));
            }
            v.push('\n');
        }
    }
    v.into_bytes()
}

#[test]
fn bcf_matches_htslib_with_its_index() {
    let dir = tempfile::tempdir().unwrap();
    for seed in 0..3u64 {
        let mut rng = Rng::new(seed);
        let text = random_vcf(&mut rng, 30_000);
        let plain = dir.path().join("ref.vcf");
        std::fs::write(&plain, &text).unwrap();
        let reference = dir.path().join("ref.bcf");
        hts::bcf_of(&plain, &reference);
        let want = inflate(&reference);
        let mut first: Option<(Vec<u8>, Vec<u8>)> = None;
        for storage in STORAGES {
            let src = dir.path().join("in.vcf");
            write_input(&src, &text, storage);
            for budget in budgets() {
                let dst = dir.path().join("out.bcf");
                let mut conversion = VcfToBcf::new(&src, &dst);
                conversion.index = true;
                convert(Conversion::VcfToBcf(conversion), budget.clone(), None).unwrap();
                assert_eq!(inflate(&dst), want, "{storage:?} {budget:?}");
                let ours = dir.path().join("out.bcf.csi");
                let theirs = dir.path().join("theirs.csi");
                hts::bcf_index(&dst, &theirs);
                assert_eq!(hts::parse_index(&ours), hts::parse_index(&theirs));
                let files = (std::fs::read(&dst).unwrap(), std::fs::read(&ours).unwrap());
                match &first {
                    None => first = Some(files),
                    Some(f) => assert_eq!(f, &files, "{storage:?} {budget:?}"),
                }
            }
        }
    }
}

#[test]
fn missing_contigs_are_declared_on_request() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(42);
    let full = random_vcf(&mut rng, 3_000);
    let text = String::from_utf8(full).unwrap();
    // The same file without its ##contig lines, and the header our first
    // pass makes of it (contigs in order of use, no lengths).
    let bare: String = text
        .lines()
        .filter(|l| !l.starts_with("##contig"))
        .map(|l| format!("{l}\n"))
        .collect();
    let declared: String = bare
        .lines()
        .map(|l| {
            if l.starts_with("#CHROM") {
                format!("##contig=<ID=chr1>\n##contig=<ID=chr2>\n##contig=<ID=chrM>\n{l}\n")
            } else {
                format!("{l}\n")
            }
        })
        .collect();
    let src = dir.path().join("bare.vcf");
    std::fs::write(&src, &bare).unwrap();
    let reference = dir.path().join("declared.vcf");
    std::fs::write(&reference, &declared).unwrap();
    let want_path = dir.path().join("want.bcf");
    hts::bcf_of(&reference, &want_path);
    let dst = dir.path().join("out.bcf");
    assert!(
        convert(
            Conversion::VcfToBcf(VcfToBcf::new(&src, &dst)),
            budgets()[1].clone(),
            None
        )
        .is_err()
    );
    let mut conversion = VcfToBcf::new(&src, &dst);
    conversion.add_missing_contigs = true;
    conversion.index = true;
    convert(Conversion::VcfToBcf(conversion), budgets()[1].clone(), None).unwrap();
    assert_eq!(inflate(&dst), inflate(&want_path));
}

#[test]
fn undeclared_contigs_fail_without_output() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("a.vcf");
    std::fs::write(
        &src,
        b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr9\t5\t.\tA\tC\t.\t.\t.\n",
    )
    .unwrap();
    let dst = dir.path().join("a.bcf");
    let err = convert(
        Conversion::VcfToBcf(VcfToBcf::new(&src, &dst)),
        budgets()[1].clone(),
        None,
    )
    .unwrap_err();
    assert!(err.to_string().contains("line 3"), "{err}");
    assert!(!dst.exists());
}
