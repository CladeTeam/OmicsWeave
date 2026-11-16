//! The `.csi` of a BGZF BCF, built from the file as `bcftools index` builds
//! it (htslib's `bcf_index_build3`): one pass over the records, each pushed
//! with its stored reference span.

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use omicsweave_bgzf::{BinningIndex, IndexFormat, Scanner, inflate_queue, try_csi_levels};
use omicsweave_core::exec::Ordered;
use omicsweave_core::{Error, Result, sidecar};

use crate::bcf::{Cursor, Fixed, checked_key, parse_record_fields, span_int, width};
use crate::fields::{BT_CHAR, BT_INT8, BT_INT16, BT_INT32, BT_INT64, BT_NULL};
use crate::header::{Header, ValueType};
use crate::record::Record;

/// htslib's default `min_shift` for BCF indexes: bins of 16 KiB.
pub const CSI_MIN_SHIFT: u32 = 14;

fn bad(detail: impl Into<String>) -> Error {
    Error::invalid("bcf record", detail)
}

/// A section cannot consume the next section's bytes. Ordinary payloads
/// are discarded a BGZF window at a time; span repair retains actual bytes
/// only, without allocating from untrusted lengths in advance.
struct Section<'a, R: Read> {
    scan: &'a mut Scanner<R>,
    left: usize,
    capture: Option<Vec<u8>>,
}

impl<'a, R: Read> Section<'a, R> {
    fn new(scan: &'a mut Scanner<R>, len: usize, capture: bool) -> Self {
        Section {
            scan,
            left: len,
            capture: capture.then(Vec::new),
        }
    }

    fn take(&mut self, n: usize) -> Result<&[u8]> {
        if n > self.left {
            return Err(bad("section malformed or too short"));
        }
        let bytes = self.scan.take(n, "bcf record")?;
        self.left -= n;
        if let Some(out) = &mut self.capture {
            out.try_reserve(bytes.len())
                .map_err(|e| bad(e.to_string()))?;
            out.extend_from_slice(bytes);
        }
        Ok(bytes)
    }

    fn skip(&mut self, n: usize) -> Result<()> {
        if n > self.left {
            return Err(bad("section malformed or too short"));
        }
        if self.capture.is_none() {
            self.scan.skip(n, "bcf record")?;
            self.left -= n;
        } else {
            let mut left = n;
            while left > 0 {
                let available = self.scan.fill(1)?;
                if available == 0 {
                    return Err(Error::Truncated {
                        what: "bcf record",
                        need: n as u64,
                        have: (n - left) as u64,
                    });
                }
                let chunk = left.min(available);
                self.take(chunk)?;
                left -= chunk;
            }
        }
        Ok(())
    }

    fn descriptor(&mut self) -> Result<(u8, usize)> {
        // A descriptor holds at most an extra descriptor and an int64
        // length. Cursor supplies the same validation as record decoding.
        let mut bytes = [0; 10];
        bytes[0] = self.take(1)?[0];
        let mut len = 1;
        if bytes[0] >> 4 == 15 {
            bytes[1] = self.take(1)?[0];
            let w = width(bytes[1] & 0xf)?;
            bytes[2..2 + w].copy_from_slice(self.take(w)?);
            len = 2 + w;
        }
        let (ty, n) = Cursor {
            buf: &bytes[..len],
            at: 0,
        }
        .descriptor()?;
        if ty == BT_NULL && n != 0 {
            return Err(bad("nonempty null vector"));
        }
        width(ty)?;
        Ok((ty, n))
    }

    fn int(&mut self) -> Result<i64> {
        let (ty, n) = self.descriptor()?;
        if n != 1 || !matches!(ty, BT_INT8 | BT_INT16 | BT_INT32 | BT_INT64) {
            return Err(bad("expected a typed integer"));
        }
        let mut bytes = [0; 9];
        bytes[0] = 0x10 | ty;
        let w = width(ty)?;
        bytes[1..1 + w].copy_from_slice(self.take(w)?);
        Cursor {
            buf: &bytes[..1 + w],
            at: 0,
        }
        .int()
    }

    fn vector(&mut self) -> Result<(u8, usize, usize)> {
        let (ty, n) = self.descriptor()?;
        let size = n
            .checked_mul(width(ty)?)
            .ok_or_else(|| bad("vector size"))?;
        if size > self.left {
            return Err(bad("section malformed or too short"));
        }
        Ok((ty, n, size))
    }
}

fn record_span<R: Read>(
    scan: &mut Scanner<R>,
    header: &Header,
    fixed: &Fixed,
    shared_len: usize,
    indiv_len: usize,
    gt_id: Option<u32>,
) -> Result<i64> {
    let repair = fixed.rlen < 0;
    let mut shared = Section::new(scan, shared_len, repair);
    let (ty, _, bytes) = shared.vector()?;
    if ty != BT_CHAR {
        return Err(bad("ID is not a string"));
    }
    shared.skip(bytes)?;
    if fixed.n_allele == 0 {
        return Err(bad("missing REF allele"));
    }
    for _ in 0..fixed.n_allele {
        let (ty, _, bytes) = shared.vector()?;
        if ty != BT_CHAR {
            return Err(bad("allele is not a string"));
        }
        shared.skip(bytes)?;
    }
    let (ty, count, _) = shared.vector()?;
    if count != 0 {
        if !matches!(ty, BT_INT8 | BT_INT16 | BT_INT32 | BT_INT64) {
            return Err(bad("FILTER is not an integer vector"));
        }
        let w = width(ty)?;
        for _ in 0..count {
            let key = span_int(ty, shared.take(w)?).ok_or_else(|| bad("missing FILTER id"))?;
            checked_key(header, key)?;
        }
    }
    for _ in 0..fixed.n_info {
        checked_key(header, shared.int()?)?;
        let (_, _, bytes) = shared.vector()?;
        shared.skip(bytes)?;
    }
    if shared.left != 0 {
        return Err(bad("trailing bytes in shared block"));
    }
    let shared = shared.capture;

    // bcf_read1_core clears stale n_fmt on empty individual/sample data.
    let n_fmt = if indiv_len == 0 || fixed.n_sample == 0 {
        0
    } else {
        fixed.n_fmt
    };
    let samples = header.samples.len();
    if n_fmt > 0 && fixed.n_sample as usize != samples {
        return Err(bad(format!(
            "{} samples, the header has {samples}",
            fixed.n_sample
        )));
    }
    let mut indiv = Section::new(scan, indiv_len, repair);
    for _ in 0..n_fmt {
        let key = checked_key(header, indiv.int()?)?;
        let (ty, n) = indiv.descriptor()?;
        let bytes = n
            .checked_mul(width(ty)?)
            .and_then(|n| n.checked_mul(samples))
            .filter(|&n| n <= i32::MAX as usize)
            .ok_or_else(|| bad("FORMAT size"))?;
        if Some(key) == gt_id && !matches!(ty, BT_INT8 | BT_INT16 | BT_INT32) {
            return Err(bad("GT is not an integer vector"));
        }
        indiv.skip(bytes)?;
    }
    if n_fmt > 0 && indiv.left != 0 {
        return Err(bad("trailing bytes in FORMAT block"));
    }
    indiv.skip(indiv.left)?;
    if let Some(shared) = shared {
        let mut rec = Record::default();
        parse_record_fields::<false>(header, fixed, &shared, &indiv.capture.unwrap(), &mut rec)?;
        Ok(rec.rlen)
    } else {
        Ok(fixed.rlen as i64)
    }
}

/// The header of a BCF stream, read past its text.
fn read_header<R: Read>(scan: &mut Scanner<R>) -> Result<Header> {
    let magic = scan.take(5, "bcf header")?;
    if magic[..4] != *b"BCF\x02" || !matches!(magic[4], 1 | 2) {
        return Err(Error::invalid(
            "bcf header",
            "not a BCF 2 file (an indexed BCF is BGZF-compressed)",
        ));
    }
    let len = u32::from_le_bytes(scan.take(4, "bcf header")?.try_into().unwrap()) as usize;
    let text = scan.take(len, "bcf header")?;
    if text.last() != Some(&0) {
        return Err(Error::invalid("bcf header", "missing terminating NUL"));
    }
    let text = std::str::from_utf8(text)
        .map_err(|_| Error::invalid("bcf header", "header is not UTF-8"))?;
    Header::parse(text)
}

/// Builds the CSI of a BGZF BCF read from `input` (blocks inflated on
/// `exec`), with bins of `2^min_shift` bases and as many levels as the
/// longest contig needs (`2^31 - 1` bases when no contig has a length), as
/// `bcftools index` does; returns the index file's bytes. Records must be
/// sorted by contig and position.
pub fn build_index<R: Read>(
    input: R,
    exec: Box<dyn Ordered<Vec<u8>>>,
    min_shift: u32,
) -> Result<Vec<u8>> {
    let mut scan = Scanner::new(input, exec);
    let header = read_header(&mut scan)?;
    let gt_id = header.key_id("GT").filter(|&id| {
        header
            .key(id)
            .is_some_and(|key| key.format == Some(ValueType::String))
    });
    let max_len = header
        .contigs
        .iter()
        .filter_map(|c| c.len)
        .max()
        .filter(|&len| len != 0)
        .unwrap_or((1 << 31) - 1);
    let (min_shift, levels) = try_csi_levels(max_len, min_shift)?;
    let format = IndexFormat::Csi { min_shift, levels };
    let mut index = BinningIndex::try_new(format, header.contigs.len(), scan.tell())?;
    while scan.fill(1)? > 0 {
        let lens = scan.take(8, "bcf record")?;
        // SAFETY: take returned exactly eight initialized bytes. These
        // unaligned loads cover the first and second four-byte fields;
        // both pointers remain within the borrowed buffer.
        let (l_shared, l_indiv) = unsafe {
            let p = lens.as_ptr();
            (
                u32::from_le(p.cast::<u32>().read_unaligned()) as usize,
                u32::from_le(p.add(4).cast::<u32>().read_unaligned()) as usize,
            )
        };
        if l_shared < 24 {
            return Err(Error::invalid("bcf record", "shared block too short"));
        }
        (l_shared - 24)
            .checked_add(l_indiv)
            .ok_or_else(|| Error::invalid("bcf record", "record length overflows usize"))?;
        let fixed = Fixed::parse(scan.take(24, "bcf record")?.try_into().unwrap());
        if fixed.rid < 0
            || header
                .contigs
                .get(fixed.rid as usize)
                .is_none_or(|c| c.name.is_empty())
        {
            return Err(Error::invalid("bcf record", "invalid contig id"));
        }
        let span = record_span(&mut scan, &header, &fixed, l_shared - 24, l_indiv, gt_id)?;
        let beg = if fixed.pos == u32::MAX {
            -1
        } else {
            fixed.pos as i64
        };
        let end = beg
            .checked_add(span)
            .ok_or_else(|| bad("record end overflows i64"))?;
        index.push(fixed.rid, beg, end, scan.tell(), true)?;
    }
    let end = scan.tell();
    index.finish(end, scan.blocks())
}

/// [`build_index`] of the BCF file at `bcf`, inflated on `threads` workers
/// and saved beside it as `<bcf>.csi` (replacing one there); returns its
/// path.
pub fn build_index_file(bcf: &Path, min_shift: u32, threads: usize) -> Result<PathBuf> {
    let input = BufReader::with_capacity(1 << 20, File::open(bcf)?);
    let bytes = build_index(input, inflate_queue(threads), min_shift)?;
    let mut path = bcf.as_os_str().to_owned();
    path.push(".csi");
    let path = PathBuf::from(path);
    sidecar::write_file(&path, &bytes)?;
    Ok(path)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ffi::CString;

    use omicsweave_bgzf::{BinIndex, Writer};
    use omicsweave_core::exec::Inline;

    fn bcf(length: u64, spans: &[(u32, i32)], trim: usize) -> Vec<u8> {
        let text = format!(
            "##fileformat=VCFv4.2\n##contig=<ID=c,length={length}>\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        );
        let header = Header::from_vcf_text(text.as_bytes()).unwrap();
        let mut bytes = crate::bcfenc::bcf_header_bytes(&header).unwrap();
        for &(pos, len) in spans {
            let mut fixed = [0; 24];
            fixed[4..8].copy_from_slice(&pos.to_le_bytes());
            fixed[8..12].copy_from_slice(&len.to_le_bytes());
            fixed[18..20].copy_from_slice(&1u16.to_le_bytes());
            bytes.extend_from_slice(&29u32.to_le_bytes());
            bytes.extend_from_slice(&0u32.to_le_bytes());
            bytes.extend_from_slice(&fixed);
            bytes.extend_from_slice(&[0x17, b'.', 0x17, b'A', 0]);
        }
        bytes.truncate(bytes.len() - trim);
        let mut writer = Writer::new(Vec::new(), 1, Box::new(Inline::new())).unwrap();
        writer.write(&bytes).unwrap();
        writer.finish().unwrap().out
    }

    #[test]
    fn high_positions_and_telomeres_keep_their_unsigned_bcf_coordinates() {
        let high = (1u32 << 31) + 99;
        let raw = bcf(3_000_000_000, &[(u32::MAX, 1), (high, 1)], 0);
        let bytes = build_index(&raw[..], inflate_queue(1), CSI_MIN_SHIFT).unwrap();
        let index = BinIndex::parse(&bytes).unwrap();
        assert_eq!(index.counts(0), Some((2, 0)));
        assert!(!index.chunks(0, 0, 1).is_empty());
        assert!(!index.chunks(0, high as i64, high as i64 + 1).is_empty());
    }

    #[test]
    fn zero_contig_lengths_use_the_unknown_length_default() {
        let raw = bcf(0, &[(100_000, 1)], 0);
        let bytes = build_index(&raw[..], inflate_queue(1), CSI_MIN_SHIFT).unwrap();
        let index = BinIndex::parse(&bytes).unwrap();
        let (shift, levels) = try_csi_levels((1 << 31) - 1, CSI_MIN_SHIFT).unwrap();
        assert_eq!((index.min_shift(), index.levels()), (shift, levels));
    }

    #[test]
    fn negative_spans_are_repaired_and_truncated_record_bodies_are_rejected() {
        let raw = bcf(1_000, &[(10, -1)], 0);
        let bytes = build_index(&raw[..], inflate_queue(1), CSI_MIN_SHIFT).unwrap();
        let index = BinIndex::parse(&bytes).unwrap();
        assert!(!index.chunks(0, 10, 11).is_empty());
        assert_eq!(oracle(&raw), (0, 1));
        let raw = bcf(1_000, &[(10, 1)], 2);
        assert!(matches!(
            build_index(&raw[..], inflate_queue(1), CSI_MIN_SHIFT),
            Err(Error::Truncated {
                what: "bcf record",
                ..
            })
        ));
    }

    fn pack(header: &Header, fixed: &[u8; 24], shared: &[u8], indiv: &[u8]) -> Vec<u8> {
        let mut bytes = crate::bcfenc::bcf_header_bytes(header).unwrap();
        bytes.extend_from_slice(&((24 + shared.len()) as u32).to_le_bytes());
        bytes.extend_from_slice(&(indiv.len() as u32).to_le_bytes());
        bytes.extend_from_slice(fixed);
        bytes.extend_from_slice(shared);
        bytes.extend_from_slice(indiv);
        let mut writer = Writer::new(Vec::new(), 1, Box::new(Inline::new())).unwrap();
        writer.write(&bytes).unwrap();
        writer.finish().unwrap().out
    }

    fn oracle(raw: &[u8]) -> (i32, i64) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("record.bcf");
        std::fs::write(&path, raw).unwrap();
        let path = CString::new(path.to_str().unwrap()).unwrap();
        // SAFETY: the path and mode are live NUL-terminated strings. Every
        // opaque handle is checked before use and destroyed exactly once.
        unsafe {
            let file = htslib_sys::hts_open(path.as_ptr(), c"rb".as_ptr());
            assert!(!file.is_null());
            let header = htslib_sys::bcf_hdr_read(file);
            assert!(!header.is_null());
            let rec = htslib_sys::bcf_init();
            assert!(!rec.is_null());
            let status = htslib_sys::bcf_read(file, header, rec);
            let span = if status == 0 {
                htslib_sys::ow_bcf_rlen(rec)
            } else {
                0
            };
            htslib_sys::bcf_destroy(rec);
            htslib_sys::bcf_hdr_destroy(header);
            htslib_sys::hts_close(file);
            (status, span)
        }
    }

    fn fixture() -> (Header, [u8; 24], Vec<u8>) {
        let header = Header::from_vcf_text(
            b"##fileformat=VCFv4.4\n##contig=<ID=c,length=100000>\n\
              ##INFO=<ID=DP,Number=1,Type=Integer>\n\
              ##FORMAT=<ID=GT,Number=1,Type=String>\n\
              ##FORMAT=<ID=AD,Number=R,Type=Integer>\n\
              #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n",
        )
        .unwrap();
        let mut fixed = [0; 24];
        fixed[4..8].copy_from_slice(&9u32.to_le_bytes());
        fixed[8..12].copy_from_slice(&1i32.to_le_bytes());
        fixed[18..20].copy_from_slice(&1u16.to_le_bytes());
        (header, fixed, vec![0x17, b'.', 0x17, b'A', 0])
    }

    #[test]
    fn malformed_shared_fields_cannot_receive_a_successful_index() {
        let (header, fixed, shared) = fixture();
        let mut cases = Vec::new();
        // ID and REF type, FILTER id/type, and the descriptor length.
        for body in [
            vec![0x11, 0, 0x17, b'A', 0],
            vec![0x17, b'.', 0x11, 1, 0],
            vec![0x17, b'.', 0x17, b'A', 0x11, 99],
            vec![0x17, b'.', 0x17, b'A', 0x15, 0, 0, 0, 0],
            vec![0xf7, 0x11, 0xff],
            vec![0xf7, 0x13, 0xff, 0xff, 0xff, 0x7f],
        ] {
            cases.push((fixed, body));
        }
        let dp = header.key_id("DP").unwrap() as u8;
        for extra in [vec![0x11, 99, 0x11, 1], vec![0x11, dp, 0x10]] {
            let mut fields = shared.clone();
            fields.extend_from_slice(&extra);
            let mut fixed = fixed;
            fixed[16..18].copy_from_slice(&1u16.to_le_bytes());
            cases.push((fixed, fields));
        }
        for (fixed, shared) in cases {
            let raw = pack(&header, &fixed, &shared, &[]);
            assert!(build_index(&raw[..], inflate_queue(1), CSI_MIN_SHIFT).is_err());
            assert!(oracle(&raw).0 < -1);
        }
    }

    #[test]
    fn format_lengths_samples_and_dictionary_ids_are_checked() {
        let (header, mut fixed, shared) = fixture();
        fixed[20..24].copy_from_slice(&((1u32 << 24) | 2).to_le_bytes());
        let ad = header.key_id("AD").unwrap() as u8;
        let bodies = [
            vec![0x11, 99, 0x11, 1, 2],
            vec![0x11, ad, 0x10],
            vec![0x11, ad, 0xf1, 0x13, 0xff, 0xff, 0xff, 0x7f],
            vec![0x11, ad, 0x21, 1, 2],
        ];
        for indiv in bodies {
            let raw = pack(&header, &fixed, &shared, &indiv);
            assert!(build_index(&raw[..], inflate_queue(1), CSI_MIN_SHIFT).is_err());
            assert!(oracle(&raw).0 < -1);
        }
        fixed[20..24].copy_from_slice(&((1u32 << 24) | 1).to_le_bytes());
        let raw = pack(&header, &fixed, &shared, &[0x11, ad, 0x11, 1]);
        assert!(build_index(&raw[..], inflate_queue(1), CSI_MIN_SHIFT).is_err());
    }

    #[test]
    fn stale_format_counts_are_ignored_without_individual_or_sample_data() {
        let (header, mut fixed, shared) = fixture();
        for (samples, indiv) in [(2u32, &[][..]), (0, &[0xff][..])] {
            fixed[20..24].copy_from_slice(&((1u32 << 24) | samples).to_le_bytes());
            let raw = pack(&header, &fixed, &shared, indiv);
            assert!(build_index(&raw[..], inflate_queue(1), CSI_MIN_SHIFT).is_ok());
            assert_eq!(oracle(&raw), (0, 1));
        }
    }

    #[test]
    fn negative_spans_use_existing_end_svlen_and_format_len_repair() {
        let header = Header::from_vcf_text(
            b"##fileformat=VCFv4.5\n##contig=<ID=c,length=100000>\n\
              ##INFO=<ID=END,Number=1,Type=Integer>\n\
              ##INFO=<ID=SVLEN,Number=A,Type=Integer>\n\
              ##FORMAT=<ID=GT,Number=1,Type=String>\n\
              ##FORMAT=<ID=LEN,Number=1,Type=Integer>\n\
              #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n",
        )
        .unwrap();
        let encoder = crate::BcfEncoder::new(&header);
        let mut scratch = crate::BcfScratch::default();
        for (line, span) in [
            (
                b"c\t10\tv\tA\tC\t.\tPASS\tEND=25\tGT:LEN\t0|1:3\t1:7".as_slice(),
                16,
            ),
            (
                b"c\t10\tv\tA\t<DEL>\t.\tPASS\tEND=11;SVLEN=100\tGT:LEN\t0|1:3\t1:7",
                101,
            ),
            (
                b"c\t10\tv\tA\t<NON_REF>\t.\tPASS\tEND=11\tGT:LEN\t0|1:3\t1:70",
                70,
            ),
        ] {
            let mut bytes = Vec::new();
            encoder
                .encode_line(&mut scratch, line, 1, &mut bytes)
                .unwrap();
            let shared_len = u32::from_le_bytes(bytes[..4].try_into().unwrap()) as usize;
            let mut fixed: [u8; 24] = bytes[8..32].try_into().unwrap();
            fixed[8..12].copy_from_slice(&(-1i32).to_le_bytes());
            let shared = &bytes[32..8 + shared_len];
            let indiv = &bytes[8 + shared_len..];
            let raw = pack(&header, &fixed, shared, indiv);
            assert_eq!(oracle(&raw), (0, span));
            let mut scan = Scanner::new(&raw[..], inflate_queue(1));
            read_header(&mut scan).unwrap();
            scan.skip(32, "bcf record").unwrap();
            assert_eq!(
                record_span(
                    &mut scan,
                    &header,
                    &Fixed::parse(&fixed),
                    shared.len(),
                    indiv.len(),
                    header.key_id("GT")
                )
                .unwrap(),
                span
            );
            let bytes = build_index(&raw[..], inflate_queue(1), CSI_MIN_SHIFT).unwrap();
            assert_eq!(BinIndex::parse(&bytes).unwrap().counts(0), Some((1, 0)));
        }
    }

    #[test]
    fn repaired_spans_cannot_overflow_the_index_end() {
        let header = Header::from_vcf_text(
            b"##fileformat=VCFv4.5\n##contig=<ID=c,length=100000>\n\
              ##FORMAT=<ID=LEN,Number=1,Type=Integer>\n\
              #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n",
        )
        .unwrap();
        let mut fixed = [0; 24];
        fixed[4..8].copy_from_slice(&9u32.to_le_bytes());
        fixed[8..12].copy_from_slice(&(-1i32).to_le_bytes());
        fixed[18..20].copy_from_slice(&2u16.to_le_bytes());
        fixed[20..24].copy_from_slice(&((1u32 << 24) | 1).to_le_bytes());
        let mut shared = vec![0x17, b'.', 0x17, b'A', 0x97];
        shared.extend_from_slice(b"<NON_REF>");
        shared.push(0);
        let mut indiv = vec![0x11, header.key_id("LEN").unwrap() as u8, 0x14];
        indiv.extend_from_slice(&i64::MAX.to_le_bytes());
        let raw = pack(&header, &fixed, &shared, &indiv);
        let err = build_index(&raw[..], inflate_queue(1), CSI_MIN_SHIFT).unwrap_err();
        assert!(
            err.to_string().contains("record end overflows i64"),
            "{err}"
        );
    }
}
