//! Public indexed tables: dense padding, sample order, filters and reuse.

use crate::indexed_fixture::{Fixture, fixture};
use omicsweave_bgzf::BinIndex;
use omicsweave_core::{Error, FileSource, Region};
use omicsweave_vcf::{
    BcfReader, FLOAT_MISSING, FLOAT_VECTOR_END, GT_ABSENT, INT32_MISSING, INT32_VECTOR_END,
    RecordFilter, RecordQuery, Values, VariantTable,
};

const VCF: &str = "##fileformat=VCFv4.4\n\
    ##contig=<ID=chr1,length=2000>\n\
    ##contig=<ID=chr2,length=2000>\n\
    ##FILTER=<ID=q10,Description=\"Low quality\">\n\
    ##FILTER=<ID=Low,Description=\"Low depth\">\n\
    ##INFO=<ID=I,Number=.,Type=Integer,IDX=20>\n\
    ##INFO=<ID=F,Number=.,Type=Float>\n\
    ##INFO=<ID=S,Number=1,Type=String>\n\
    ##INFO=<ID=B,Number=0,Type=Flag>\n\
    ##FORMAT=<ID=GT,Number=1,Type=String>\n\
    ##FORMAT=<ID=IV,Number=.,Type=Integer>\n\
    ##FORMAT=<ID=FV,Number=.,Type=Float>\n\
    ##FORMAT=<ID=SV,Number=1,Type=String>\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
    chr1\t10\tfirst\tA\tC\t12.5\tPASS\tI=1,.;F=0.5;S=alpha;B\tGT:IV:FV:SV\t0|1:1,.:0.25:aa\t1:9:.,1.5:bbb\n\
    chr1\t20\tmiddle\tAT\tA\t.\tq10;Low\t.\tGT:SV\t.:x\t0/0:.\n\
    chr1\t30\tlast\tA\tG\t99\t.\tI=-200;F=.,2;S=omega\tGT:IV\t1/1:4\t0:5,6,7\n";

fn open(fixture: &Fixture) -> BcfReader<FileSource> {
    BcfReader::open(
        FileSource::open(&fixture.bcf).unwrap(),
        BinIndex::read_file(&fixture.csi).unwrap(),
    )
    .unwrap()
}

fn query() -> RecordQuery {
    RecordQuery {
        info: ["S", "I", "B", "F", "I"].map(str::to_owned).to_vec(),
        format: ["SV", "IV", "FV", "IV"].map(str::to_owned).to_vec(),
    }
}

fn strings(values: &omicsweave_vcf::Strings) -> Vec<&[u8]> {
    (0..values.len()).map(|i| values.get(i)).collect()
}

fn float_bits(values: &Values) -> (usize, Vec<u32>) {
    let Values::Float { width, values } = values else {
        panic!("expected float column, got {values:?}");
    };
    (*width, values.iter().map(|value| value.to_bits()).collect())
}

#[test]
fn records_preserve_selected_sample_order_duplicate_columns_and_dense_sentinels() {
    let fixture = fixture(VCF, 14);
    let mut reader = open(&fixture);
    let samples = [1, 0, 1];
    let region = Region::new("chr1", 0, 100);
    let table = reader.records(region, &samples, &query(), None).unwrap();
    assert_eq!(table.genotypes.start, [9, 19, 29]);
    assert_eq!(table.genotypes.end, [10, 21, 30]);
    assert_eq!((table.genotypes.samples, table.genotypes.ploidy), (3, 2));
    assert_eq!(
        &table.genotypes.gt[..6],
        &[1, GT_ABSENT, 0, 1, 1, GT_ABSENT]
    );
    assert_eq!(table.genotypes, reader.genotypes(region, &samples).unwrap());
    assert_eq!(
        strings(&table.genotypes.id),
        [b"first".as_slice(), b"middle", b"last"]
    );
    assert_eq!(
        table.qual.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        [12.5f32.to_bits(), FLOAT_MISSING, 99f32.to_bits()]
    );
    assert_eq!(
        strings(&table.filters),
        [b"PASS".as_slice(), b"q10;Low", b""]
    );
    assert_eq!(
        table
            .info
            .iter()
            .map(|c| c.name.as_str())
            .collect::<Vec<_>>(),
        ["S", "I", "B", "F", "I"]
    );
    assert_eq!(
        table
            .format
            .iter()
            .map(|c| c.name.as_str())
            .collect::<Vec<_>>(),
        ["SV", "IV", "FV", "IV"]
    );
    let Values::Str(ref info_strings) = table.info[0].values else {
        panic!("INFO/S")
    };
    assert_eq!(strings(info_strings), [b"alpha".as_slice(), b"", b"omega"]);
    let ve = INT32_VECTOR_END;
    let missing = INT32_MISSING;
    assert_eq!(
        table.info[1].values,
        Values::Int {
            width: 2,
            values: vec![1, missing, ve, ve, -200, ve]
        }
    );
    assert_eq!(table.info[1], table.info[4]);
    assert_eq!(table.info[2].values, Values::Flag(vec![true, false, false]));
    assert_eq!(
        float_bits(&table.info[3].values),
        (
            2,
            vec![
                0.5f32.to_bits(),
                FLOAT_VECTOR_END,
                FLOAT_VECTOR_END,
                FLOAT_VECTOR_END,
                FLOAT_MISSING,
                2f32.to_bits()
            ]
        )
    );
    let Values::Str(ref format_strings) = table.format[0].values else {
        panic!("FORMAT/SV")
    };
    assert_eq!(
        strings(format_strings),
        [
            b"bbb".as_slice(),
            b"aa",
            b"bbb",
            b".",
            b"x",
            b".",
            b"",
            b"",
            b""
        ]
    );
    assert_eq!(
        table.format[1].values,
        Values::Int {
            width: 3,
            values: vec![
                9, ve, ve, 1, missing, ve, 9, ve, ve, ve, ve, ve, ve, ve, ve, ve, ve, ve, 5, 6, 7,
                4, ve, ve, 5, 6, 7
            ]
        }
    );
    assert_eq!(table.format[1], table.format[3]);
    let mut expected = vec![
        FLOAT_MISSING,
        1.5f32.to_bits(),
        0.25f32.to_bits(),
        FLOAT_VECTOR_END,
        FLOAT_MISSING,
        1.5f32.to_bits(),
    ];
    expected.extend([FLOAT_VECTOR_END; 12]);
    assert_eq!(float_bits(&table.format[2].values), (2, expected));
}

#[test]
fn include_exclude_and_empty_selection_keep_all_columns_aligned() {
    let fixture = fixture(VCF, 14);
    let mut reader = open(&fixture);
    for include in [true, false] {
        let mut filter = if include {
            RecordFilter::include(reader.header(), "INFO/B=1").unwrap()
        } else {
            RecordFilter::exclude(reader.header(), "INFO/B=1").unwrap()
        };
        for _ in 0..2 {
            let table = reader
                .records(
                    Region::new("chr1", 0, 100),
                    &[],
                    &query(),
                    Some(&mut filter),
                )
                .unwrap();
            let expected = if include { vec![9] } else { vec![19, 29] };
            assert_eq!(table.genotypes.start, expected);
            assert_eq!(table.qual.len(), expected.len());
            assert_eq!(table.filters.len(), expected.len());
            assert_eq!(table.genotypes.samples, 0);
            assert!(table.genotypes.gt.is_empty());
            for column in &table.format {
                match &column.values {
                    Values::Int { values, .. } => assert!(values.is_empty()),
                    Values::Float { values, .. } => assert!(values.is_empty()),
                    Values::Str(values) => assert!(values.is_empty()),
                    Values::Flag(_) => panic!("FORMAT flag"),
                }
            }
        }
    }
    for region in [Region::new("chr1", 50, 50), Region::new("chr2", 0, 100)] {
        let empty = reader.records(region, &[1, 0], &query(), None).unwrap();
        assert!(empty.genotypes.is_empty());
        assert_eq!(empty.genotypes.samples, 2);
        assert!(empty.qual.is_empty() && empty.filters.is_empty());
        assert_eq!(empty.info.len(), 5);
        assert_eq!(empty.format.len(), 4);
    }
}

#[test]
fn invalid_queries_fail_and_reader_can_return_owned_results_afterwards() {
    let fixture = fixture(VCF, 14);
    let mut reader = open(&fixture);
    let region = Region::new("chr1", 0, 100);
    for bad in [
        RecordQuery {
            info: vec!["unknown".into()],
            format: vec![],
        },
        RecordQuery {
            info: vec![],
            format: vec!["unknown".into()],
        },
        RecordQuery {
            info: vec![],
            format: vec!["GT".into()],
        },
    ] {
        assert!(matches!(
            reader.records(region, &[0], &bad, None),
            Err(Error::Invalid { .. })
        ));
    }
    assert!(matches!(
        reader.records(region, &[2], &query(), None),
        Err(Error::Invalid { .. })
    ));
    assert!(matches!(
        reader.records(Region::new("missing", 0, 10), &[0], &query(), None),
        Err(Error::UnknownChrom(_))
    ));
    let first: VariantTable = reader.records(region, &[0], &query(), None).unwrap();
    let mut clone = reader.try_clone().unwrap();
    let later = clone
        .records(Region::new("chr1", 29, 30), &[1], &query(), None)
        .unwrap();
    assert_eq!(later.genotypes.start, [29]);
    assert_eq!(first.genotypes.start, [9, 19, 29]);
    assert_eq!(
        strings(&first.genotypes.id),
        [b"first".as_slice(), b"middle", b"last"]
    );
}
