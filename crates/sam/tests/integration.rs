//! SAM parsing and BAM encoding by public behavior, with an independent htslib oracle.

mod support;

mod cases {
    mod aux;
    mod cigar;
    mod encoder;
    mod header;
    mod record;
}
