//! VCF and BCF variant files, read whole and in order or (BCF) by region,
//! decoded as htslib decodes them: contig, position, alleles, reference span (`rlen`, from
//! REF, INFO/END, INFO/SVLEN and FORMAT/LEN), genotypes in htslib's
//! encoding, and ID, QUAL, FILTER and every INFO and FORMAT value in the
//! typed form htslib holds them in (see [`InfoField`], [`FormatField`] and
//! the `Record::info_*`/`format_*` getters).
//!
//! [`VariantReader`] detects the format (VCF text or BCF) and compression
//! (none, gzip or BGZF) from the first bytes and fills a reusable
//! [`Record`] per call, so reading allocates only while records grow.
//!
//! [`BcfReader`] reads a BGZF BCF by region through its CSI index, as
//! htslib's `bcf_itr_queryi` does, and gives the [`GenotypeMatrix`] of a
//! region.
//!
//! [`read_vcf_header`], [`bcf_header_bytes`] and [`BcfEncoder`] convert VCF
//! text to the uncompressed BCF bytes htslib writes from it, a line at a
//! time and from any number of threads.

mod bcf;
mod bcfenc;
mod cnum;
mod fields;
mod filter;
mod header;
mod hrec;
mod indexed;
mod reader;
mod record;
mod text;

pub use bcfenc::{BcfEncoder, BcfScratch, BcfSpan, bcf_header_bytes, read_vcf_header};
pub use fields::{
    BT_CHAR, BT_FLOAT, BT_INT8, BT_INT16, BT_INT32, BT_INT64, BT_NULL, FLOAT_MISSING,
    FLOAT_VECTOR_END, FormatField, GetError, INT32_MISSING, INT32_VECTOR_END, InfoField, Typed,
    Value, type_width,
};
pub use filter::Filter;
pub use header::{Contig, Header, Key, LineKind, Number, ValueType};
pub use indexed::{BcfReader, GT_ABSENT, GT_MISSING, GenotypeMatrix, Strings};
pub use reader::{Format, VariantReader};
pub use record::{
    GT_MISSING_VALUE, GT_VECTOR_END, Genotypes, Record, VCF_BND, VCF_INDEL, VCF_MNP, VCF_OTHER,
    VCF_OVERLAP, VCF_REF, VCF_SNP, gt_allele, gt_is_missing, gt_is_phased, variant, variant_type,
    variant_types,
};
