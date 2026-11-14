//! Integer workspace and shared-dictionary fixtures for real htslib oracles.

pub fn fixture(
    samples: usize,
    width: usize,
    keys: usize,
    records: usize,
    values: &[&str],
    version: &str,
) -> String {
    let mut text = format!("##fileformat=VCFv{version}\n##contig=<ID=chr1>\n");
    for key in 0..keys {
        text.push_str(&format!(
            "##INFO=<ID=I{key},Number=1,Type=Integer,Description=\"integer\">\n\
             ##FORMAT=<ID=I{key},Number=.,Type=Integer,Description=\"integers\">\n"
        ));
    }
    text.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n");
    text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if samples > 0 {
        text.push_str("\tFORMAT");
        for sample in 0..samples {
            text.push_str(&format!("\tS{sample}"));
        }
    }
    text.push('\n');
    for record in 0..records {
        text.push_str(&format!("chr1\t{}\t.\tA\tC\t.\tPASS\t", record + 1));
        for key in 0..keys {
            if key > 0 {
                text.push(';');
            }
            text.push_str(&format!("I{key}={}", values[(record + key) % values.len()]));
        }
        if samples > 0 {
            text.push_str("\tGT");
            for key in 0..keys {
                text.push_str(&format!(":I{key}"));
            }
            for sample in 0..samples {
                text.push('\t');
                text.push_str(["0|1", "0/0", ".", "1", "./."][(sample + record) % 5]);
                if sample % 7 == 6 {
                    continue;
                }
                for key in 0..keys {
                    text.push(':');
                    for i in 0..1 + (sample + record) % width {
                        if i > 0 {
                            text.push(',');
                        }
                        text.push_str(values[(sample + record + key + i) % values.len()]);
                    }
                }
            }
        }
        text.push('\n');
    }
    text
}
