use super::{Overlap, VcfStore};
use std::fmt::Write;

fn fixture(monotone: bool) -> VcfStore {
    let mut text = String::from(
        "##fileformat=VCFv4.3\n##contig=<ID=c,length=100000>\n\
         ##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    // An unsorted stream, duplicate positions, a long first interval and
    // many short intervals exercise both search paths and stable ties.
    for i in (0..1024u64).rev() {
        let pos = (i * 311) % 1024;
        let span = if monotone {
            7
        } else if pos == 0 {
            100_000
        } else {
            3 + ((i * 7919) % 199)
        };
        let alt = if i % 3 == 0 { "AAC" } else { "C" };
        writeln!(
            text,
            "c\t{}\t.\tAAA\t{alt}\t.\t.\tEND={}",
            pos / 2 + 1,
            pos / 2 + span
        )
        .unwrap();
    }
    VcfStore::read(text.as_bytes(), None).unwrap()
}

#[test]
fn interval_index_matches_full_scan_with_duplicates_and_nested_spans() {
    for monotone in [false, true] {
        let store = fixture(monotone);
        let c = &store.contigs[0];
        assert_eq!(c.subtree_last.is_empty(), monotone);
        for start in (-2..1200).step_by(7) {
            for len in [0, 1, 2, 8, 100, 400] {
                let end = start + len - 1;
                for overlap in [Overlap::Start, Overlap::Record, Overlap::Variant] {
                    let mut got = vec![usize::MAX];
                    store.overlapping(0, start, end, overlap, &mut got);
                    let mut want = vec![usize::MAX];
                    for (i, site) in c.sites.iter().enumerate() {
                        let mut begin = site.pos;
                        if overlap == Overlap::Variant {
                            let reference = store.allele(site, 0);
                            let mut shared = site.rlen;
                            for alt in 1..site.n_allele as usize {
                                let allele = store.allele(site, alt);
                                let mut equal = 0;
                                while equal < reference.len()
                                    && equal < allele.len()
                                    && reference[equal] == allele[equal]
                                {
                                    equal += 1;
                                }
                                shared = shared.min(equal as i64);
                            }
                            begin += shared;
                        }
                        let takes = end >= start
                            && site.pos <= end
                            && match overlap {
                                Overlap::Start => site.pos >= start,
                                Overlap::Record => site.pos + site.rlen > start,
                                Overlap::Variant => site.pos + site.rlen > start && begin <= end,
                            };
                        if takes {
                            want.push(i);
                        }
                    }
                    assert_eq!(got, want, "monotone={monotone} {start}..={end} {overlap:?}");
                }
            }
        }
    }
}

#[test]
fn interval_index_preserves_file_order_at_equal_positions() {
    let store = VcfStore::read(
        &b"##fileformat=VCFv4.3\n##contig=<ID=c>\n\
          #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
          c\t2\t.\tAA\tT\t.\t.\t.\n\
          c\t1\t.\tAAA\tC\t.\t.\t.\n\
          c\t2\t.\tA\tG\t.\t.\t.\n"[..],
        None,
    )
    .unwrap();
    let mut indices = Vec::new();
    store.overlapping(0, 1, 1, Overlap::Record, &mut indices);
    let alleles: Vec<_> = indices
        .iter()
        .map(|&i| store.allele(store.site(0, i), 1))
        .collect();
    assert_eq!(alleles, [b"C", b"T", b"G"]);
}

#[test]
fn metadata_requires_representable_exclusive_reference_end() {
    let text = b"##fileformat=VCFv4.3\n##contig=<ID=c>\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
        c\t1\t.\tA\tC\t.\t.\t.\n";
    let mut reader = omicsweave_vcf::VariantReader::new(&text[..]).unwrap();
    let mut record = omicsweave_vcf::Record::default();
    assert!(reader.read(&mut record).unwrap());
    let mut store = VcfStore::read(&text[..], None).unwrap();
    record.pos = 2;
    record.rlen = i64::MAX;
    assert!(store.push(&record).is_err());
    record.rlen = i64::MAX - 2;
    assert_eq!(store.push(&record).unwrap().last(), i64::MAX - 1);
    record.pos = -1;
    record.rlen = 0;
    assert_eq!(store.push(&record).unwrap().last(), -2);
    record.rlen = -1;
    assert!(store.push(&record).is_err());
}

#[test]
fn variant_prefix_cache_matches_alleles_and_reference_span() {
    let header = "##fileformat=VCFv4.3\n##contig=<ID=c>\n\
        ##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let cases = [
        ("AAAAC", "AAAAT,AAAGC", 3),
        ("AAAAC", "AAAACG,AAAAC", 5),
        ("AAAAC", "A,AAAAT", 1),
        ("AAAAC", "C,AAAAC", 0),
        ("AAAAC", "<NON_REF>", 0),
        ("AAAAC", ".", 17),
    ];
    let mut text = header.to_owned();
    for (i, (reference, alt, _)) in cases.iter().enumerate().rev() {
        writeln!(
            text,
            "c\t{}\t.\t{reference}\t{alt}\t.\t.\tEND={}",
            i * 20 + 1,
            i * 20 + 17
        )
        .unwrap();
    }
    let store = VcfStore::read(text.as_bytes(), None).unwrap();
    for (i, (reference, alt, prefix)) in cases.iter().enumerate() {
        let site = store.site(0, i);
        assert_eq!(store.allele(site, 0), reference.as_bytes());
        assert_eq!(store.variant_begin(site), site.pos + prefix);
        if *alt != "." {
            for (j, expected) in alt.split(',').enumerate() {
                assert_eq!(store.allele(site, j + 1), expected.as_bytes());
            }
        }
    }

    // Public record metadata can give a shorter span than its allele
    // spelling; overlap mode 2 limits its shared prefix to that span.
    let text = format!("{header}c\t1\t.\tAAAA\tAAAAT\t.\t.\t.\n");
    let mut reader = omicsweave_vcf::VariantReader::new(text.as_bytes()).unwrap();
    let mut record = omicsweave_vcf::Record::default();
    assert!(reader.read(&mut record).unwrap());
    let mut store = VcfStore::read(header.as_bytes(), None).unwrap();
    for span in [0, 1, 2, 4, 5, (1i64 << 32) + 1, i64::MAX] {
        record.rlen = span;
        let site = store.push(&record).unwrap();
        assert_eq!(store.variant_begin(&site), span.min(4));
    }
}
