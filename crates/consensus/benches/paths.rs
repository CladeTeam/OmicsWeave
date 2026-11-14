//! Public-API, memory-source workloads. Fixtures and exact output checks
//! are outside the timed groups; the identical file runs on the baseline.

use std::fmt::Write;
use std::hint::black_box;
use std::sync::Arc;
use std::time::{Duration, Instant};

use omicsweave_consensus::{
    Alleles, Consensus, ConsensusOptions, Mask, MaskWith, Overlap, VcfStore,
};
use omicsweave_core::{Region, Result, SequenceRead};

// These adapters select the unfiltered constructor across the two public
// signatures. Store creation is outside every timed group.
struct WithoutFilter;
struct WithFilter<T>(std::marker::PhantomData<T>);
trait ReadUnfiltered<'a, Signature> {
    fn read_unfiltered(self, bytes: &'a [u8]) -> Result<VcfStore>;
}
impl<'a, F> ReadUnfiltered<'a, WithoutFilter> for F
where
    F: FnOnce(&'a [u8]) -> Result<VcfStore>,
{
    fn read_unfiltered(self, bytes: &'a [u8]) -> Result<VcfStore> {
        self(bytes)
    }
}
impl<'a, F, T> ReadUnfiltered<'a, WithFilter<T>> for F
where
    F: FnOnce(&'a [u8], Option<T>) -> Result<VcfStore>,
{
    fn read_unfiltered(self, bytes: &'a [u8]) -> Result<VcfStore> {
        self(bytes, None)
    }
}
fn read_store(bytes: &[u8]) -> VcfStore {
    VcfStore::read.read_unfiltered(bytes).unwrap()
}

struct Reference(Vec<u8>);

impl SequenceRead for Reference {
    fn chrom_len(&self, _: &str) -> Result<u32> {
        Ok(self.0.len() as u32)
    }

    fn sequence(&mut self, region: Region<'_>, out: &mut [u8]) -> Result<()> {
        out.copy_from_slice(&self.0[region.start as usize..region.end as usize]);
        Ok(())
    }
}

fn measure(name: &str, mut run: impl FnMut()) {
    if std::env::var("CONSENSUS_BENCH_FILTER").is_ok_and(|f| !name.contains(&f)) {
        return;
    }
    run();
    let mut samples = [0.0; 5];
    for sample in &mut samples {
        let start = Instant::now();
        let mut n = 0u64;
        while start.elapsed() < Duration::from_millis(100) {
            run();
            n += 1;
        }
        *sample = start.elapsed().as_nanos() as f64 / n as f64;
    }
    samples.sort_by(f64::total_cmp);
    println!("{name}\t{:.2}", samples[2]);
}

fn fixture(n: u32, long_first: bool, missing: bool) -> String {
    let mut text = String::from(
        "##fileformat=VCFv4.3\n##contig=<ID=c,length=1000000>\n\
         ##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
    );
    if long_first {
        writeln!(text, "c\t1\t.\tA\t<NON_REF>\t.\t.\tEND={}\tGT\t0/0", n * 4).unwrap();
    }
    let gt = if missing { "./." } else { "1/1" };
    for i in 0..n {
        writeln!(text, "c\t{}\t.\tA\tC\t.\t.\t.\tGT\t{gt}", i * 2 + 1).unwrap();
    }
    text
}

fn output(consensus: &mut Consensus) -> Vec<u8> {
    let n = consensus
        .query(Alleles::Haplotype {
            sample: "S",
            haplotype: "1",
        })
        .unwrap();
    let mut out = vec![0; n];
    consensus.write(&mut out).unwrap();
    out
}

fn loads(long_first: bool, overlap: Overlap, width: u32) {
    const N: u32 = 131_072;
    let text = fixture(N, long_first, false);
    let store = Arc::new(read_store(text.as_bytes()));
    let mut consensus = Consensus::new(
        store,
        ConsensusOptions {
            regions_overlap: overlap,
            ..ConsensusOptions::default()
        },
    );
    let mut reference = Reference(vec![b'A'; (N * 4) as usize]);
    let region = Region::new("c", N * 2 - width * 2, N * 2 - width);
    consensus.load(&mut reference, region).unwrap();
    let expected: Vec<_> = (region.start..region.end)
        .map(|p| if p < N * 2 && p % 2 == 0 { b'C' } else { b'A' })
        .collect();
    assert_eq!(output(&mut consensus), expected);
    let shape = if long_first { "nested" } else { "monotone" };
    measure(&format!("{shape}.{overlap:?}.load.{width}"), || {
        consensus
            .load(black_box(&mut reference), black_box(region))
            .unwrap();
        black_box(&consensus);
    });
}

fn queries(n: u32, missing: bool) {
    let text = fixture(n, false, missing);
    let store = Arc::new(read_store(text.as_bytes()));
    let options = ConsensusOptions {
        absent: missing.then_some(b'.'),
        ..ConsensusOptions::default()
    };
    let mut consensus = Consensus::new(store, options);
    let mut reference = Reference(vec![b'A'; (n * 4) as usize]);
    consensus
        .load(&mut reference, Region::new("c", 0, n * 4))
        .unwrap();
    let expected: Vec<_> = (0..n * 4)
        .map(|p| {
            if missing {
                b'.'
            } else if p < n * 2 && p % 2 == 0 {
                b'C'
            } else {
                b'A'
            }
        })
        .collect();
    assert_eq!(output(&mut consensus), expected);
    let shape = if missing { "absent-skipped" } else { "snv" };
    measure(&format!("{shape}.query.{n}"), || {
        black_box(
            consensus
                .query(black_box(Alleles::Haplotype {
                    sample: "S",
                    haplotype: "1",
                }))
                .unwrap(),
        );
    });
}

fn shared_prefix_loads(prefix: usize, records: usize) {
    let mut text = String::from(
        "##fileformat=VCFv4.3\n##contig=<ID=c>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    let reference = "A".repeat(prefix + 1);
    let alt = format!("{}C", "A".repeat(prefix));
    for _ in 0..records {
        writeln!(text, "c\t1\t.\t{reference}\t{alt}\t.\t.\t.").unwrap();
    }
    let store = Arc::new(read_store(text.as_bytes()));
    let mut consensus = Consensus::new(
        store,
        ConsensusOptions {
            regions_overlap: Overlap::Variant,
            ..ConsensusOptions::default()
        },
    );
    let mut reference = Reference(vec![b'A'; prefix + 1]);
    let region = Region::new("c", (prefix - 8) as u32, (prefix - 4) as u32);
    consensus.load(&mut reference, region).unwrap();
    assert_eq!(consensus.query(Alleles::AllAlt).unwrap(), 4);
    let mut out = [0; 4];
    consensus.write(&mut out).unwrap();
    assert_eq!(&out, b"AAAA");
    measure(
        &format!("shared-prefix.Variant.load.{prefix}.{records}"),
        || {
            consensus
                .load(black_box(&mut reference), black_box(region))
                .unwrap();
            black_box(&consensus);
        },
    );
    measure(&format!("shared-prefix.store.{prefix}.{records}"), || {
        black_box(read_store(black_box(text.as_bytes())));
    });
}

fn long_reference_queries(width: usize) {
    let reference = "A".repeat(width);
    let alt = "C".repeat(width);
    let text = format!(
        "##fileformat=VCFv4.3\n##contig=<ID=c>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         c\t1\t.\t{reference}\t{alt}\t.\t.\t.\n"
    );
    let store = Arc::new(read_store(text.as_bytes()));
    let mut consensus = Consensus::new(store, ConsensusOptions::default());
    let mut reference = Reference(vec![b'A'; width]);
    consensus
        .load(&mut reference, Region::new("c", 0, width as u32))
        .unwrap();
    assert_eq!(consensus.query(Alleles::AllAlt).unwrap(), width);
    let mut out = vec![0; width];
    consensus.write(&mut out).unwrap();
    assert_eq!(out, alt.as_bytes());
    measure(&format!("long-reference.query.{width}"), || {
        black_box(consensus.query(black_box(Alleles::AllAlt)).unwrap());
    });
}

fn iupac_queries(width: usize, samples: bool, short_alt: bool) {
    const RECORDS: usize = 32;
    let mut text = String::from(
        "##fileformat=VCFv4.3\n##contig=<ID=c>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
    );
    let reference = "A".repeat(width);
    let alt_width = if short_alt { 1 } else { width };
    let alt = "C".repeat(alt_width);
    for i in 0..RECORDS {
        writeln!(
            text,
            "c\t{}\t.\t{reference}\t{alt}\t.\t.\t.\tGT\t0/1",
            i * width + 1
        )
        .unwrap();
    }
    let store = Arc::new(read_store(text.as_bytes()));
    let options = ConsensusOptions {
        iupac_codes: true,
        ..ConsensusOptions::default()
    };
    let mut consensus = Consensus::new(store.clone(), options.clone());
    let mut reference = Reference(vec![b'A'; width * RECORDS]);
    let region = Region::new("c", 0, (width * RECORDS) as u32);
    let mode = if samples {
        Alleles::Samples(&["S"])
    } else {
        Alleles::AllAlt
    };
    consensus.load(&mut reference, region).unwrap();
    let n = consensus.query(mode).unwrap();
    let mut out = vec![0; n];
    consensus.write(&mut out).unwrap();
    assert_eq!(out, vec![b'M'; alt_width * RECORDS]);
    let mode_name = if samples { "samples" } else { "all-alt" };
    let shape = if short_alt { "short-alt" } else { "equal-alt" };
    measure(&format!("iupac.{mode_name}.{shape}.warm.{width}"), || {
        black_box(consensus.query(black_box(mode)).unwrap());
    });
    // Includes builder allocation, loading, the first query, and dropping
    // its buffers; store parsing and the reference remain outside.
    measure(&format!("iupac.{mode_name}.{shape}.first.{width}"), || {
        let mut fresh = Consensus::new(store.clone(), options.clone());
        fresh.load(&mut reference, region).unwrap();
        black_box(fresh.query(black_box(mode)).unwrap());
    });
}

fn checksum(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0u64, |sum, &b| {
        sum.wrapping_mul(31).wrapping_add(u64::from(b))
    })
}

fn masked_loads(records: u32, stride: u32, mask_stride: u32) {
    let width = records * stride;
    let mut text = String::from(
        "##fileformat=VCFv4.3\n##contig=<ID=c>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    for i in 0..records {
        writeln!(text, "c\t{}\t.\tA\tC\t.\t.\t.", i * stride + 1).unwrap();
    }
    let mut bed = String::new();
    for begin in (0..width).step_by(mask_stride as usize) {
        writeln!(bed, "c\t{begin}\t{}", (begin + 2).min(width)).unwrap();
    }
    let store = Arc::new(read_store(text.as_bytes()));
    let mut consensus = Consensus::new(
        store,
        ConsensusOptions {
            mask: Some(Arc::new(
                Mask::parse(bed.as_bytes(), MaskWith::Char(b'N')).unwrap(),
            )),
            ..ConsensusOptions::default()
        },
    );
    let mut reference = Reference(vec![b'A'; width as usize]);
    let region = Region::new("c", 0, width);
    consensus.load(&mut reference, region).unwrap();
    let n = consensus.query(Alleles::AllAlt).unwrap();
    let mut out = vec![0; n];
    consensus.write(&mut out).unwrap();
    let expected: Vec<_> = (0..width)
        .map(|p| {
            if p % mask_stride < 2 {
                b'N'
            } else if p % stride == 0 {
                b'C'
            } else {
                b'A'
            }
        })
        .collect();
    assert_eq!(out, expected);
    let name = format!("round7.mask.load.{records}.{stride}.{mask_stride}");
    println!("{name}.checksum\t{}", checksum(&out));
    measure(&name, || {
        consensus
            .load(black_box(&mut reference), black_box(region))
            .unwrap();
        black_box(&consensus);
    });
}

fn gapped_queries(records: u32, chain: bool) {
    let mut text = String::from(
        "##fileformat=VCFv4.3\n##contig=<ID=c>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
    );
    for i in 0..records {
        writeln!(text, "c\t{}\t.\tA\tAC\t.\t.\t.\tGT\t0/1", i * 4 + 1).unwrap();
    }
    let store = Arc::new(read_store(text.as_bytes()));
    let mut consensus = Consensus::new(
        store,
        ConsensusOptions {
            chain,
            ..ConsensusOptions::default()
        },
    );
    let mut reference = Reference(vec![b'A'; (records * 4) as usize]);
    consensus
        .load(&mut reference, Region::new("c", 0, records * 4))
        .unwrap();
    for haplotype in ["1", "2"] {
        let mode = Alleles::Haplotype {
            sample: "S",
            haplotype,
        };
        let n = consensus.query(mode).unwrap();
        let mut out = vec![0; n];
        consensus.write(&mut out).unwrap();
        let unit: &[u8] = if haplotype == "1" { b"AAAA" } else { b"ACAAA" };
        assert_eq!(out, unit.repeat(records as usize));
        let name = format!("round7.indel.query.{records}.chain-{chain}.hap-{haplotype}");
        println!(
            "{name}.checksum\t{}",
            checksum(&out).wrapping_add(consensus.chain().map_or(0, |c| checksum(c.as_bytes())))
        );
        measure(&name, || {
            black_box(consensus.query(black_box(mode)).unwrap());
            black_box(consensus.chain());
        });
    }
}

fn main() {
    println!("case\tns/op (median of 5 x 100ms; memory source, warm buffers)");
    masked_loads(16_384, 4, 12);
    masked_loads(128, 4096, 8);
    for chain in [false, true] {
        for records in [8, 4096] {
            gapped_queries(records, chain);
        }
    }
    for long_first in [false, true] {
        for overlap in [Overlap::Start, Overlap::Record, Overlap::Variant] {
            for width in [8, 1024] {
                loads(long_first, overlap, width);
            }
        }
    }
    for n in [4096, 16384] {
        queries(n, false);
        queries(n, true);
    }
    for prefix in [64, 8192] {
        shared_prefix_loads(prefix, 128);
    }
    for width in [1, 64, 8192] {
        long_reference_queries(width);
        for samples in [false, true] {
            iupac_queries(width, samples, false);
            if width > 1 {
                iupac_queries(width, samples, true);
            }
        }
    }
    let text = fixture(4096, false, false);
    measure("snv.store.4096", || {
        black_box(read_store(black_box(text.as_bytes())));
    });
}
