mod support;

#[path = "support/htslib.rs"]
mod htslib;

#[path = "support/cram/mod.rs"]
mod cram_support;

#[path = "cases/sam_to_bam.rs"]
mod sam_to_bam;

#[path = "cases/bam_sort.rs"]
mod bam_sort;

#[path = "cases/vcf_to_bcf.rs"]
mod vcf_to_bcf;

#[path = "cases/cram_to_bam.rs"]
mod cram_to_bam;
