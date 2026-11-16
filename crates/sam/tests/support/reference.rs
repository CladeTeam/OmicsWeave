use omicsweave_sam::BamSpan;

/// One record as an encoder produced it.
pub type Record = Option<(Vec<u8>, BamSpan)>;

/// A whole SAM text encoded: the BAM header (`None` if the header is
/// rejected) and, per data line, the record (`None` if the line is).
#[derive(Debug)]
pub struct Encoded {
    pub header: Option<Vec<u8>>,
    pub records: Vec<Record>,
}
