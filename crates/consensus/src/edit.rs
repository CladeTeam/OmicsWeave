//! Applies records to a region's reference the way consensus.c does
//! (`apply_variant`, `apply_absent`, `freeze_ref`), state for state.
//!
//! consensus.c edits one buffer in place and shifts its tail at every
//! indel. Records arrive in order and never edit before an earlier edit's
//! end, so here the buffer is a modified prefix (`out`) followed by the
//! untouched rest of the reference (`src[cur..]`): an edit appends to the
//! prefix, and a region costs time linear in its length.

use omicsweave_core::{Error, Result};
use omicsweave_vcf::{VCF_INDEL, VCF_OTHER, variant};

use crate::chain::Chain;
use crate::consensus::Mark;
use crate::records::Site;
use crate::select::Choice;

/// What every record of a query is applied with.
pub(crate) struct Marks {
    pub missing: Option<u8>,
    pub absent: Option<u8>,
    pub mark_del: Option<u8>,
    pub mark_ins: Option<Mark>,
    pub mark_snv: Option<Mark>,
}

/// consensus.c's `fa_*` and `prev_*` state over the region's buffer.
pub(crate) struct Edit<'a> {
    src: &'a [u8],
    cur: usize,
    out: Vec<u8>,
    /// `fa_ori_pos`: the region's first base.
    ori_pos: i64,
    /// `fa_mod_off`: length change so far.
    mod_off: i64,
    /// `fa_frz_pos`: last reference base a record has claimed (-1 none).
    frz_pos: i64,
    /// `fa_frz_mod`: where `-a` filling may resume in the buffer.
    frz_mod: i64,
    prev_base: u8,
    prev_base_pos: i64,
    prev_is_insert: bool,
    /// Case of the first byte retained in bcftools' backing buffer after
    /// deleting its entire logical sequence; a later insertion can use it.
    empty_case: u8,
}

fn mismatch(pos: i64, reference: &[u8]) -> Error {
    Error::invalid(
        "consensus",
        format!(
            "the reference does not match REF {:?} at position {}",
            String::from_utf8_lossy(reference),
            pos + 1
        ),
    )
}

fn mark_bytes(bytes: &mut [u8], mark: Mark) {
    for b in bytes {
        *b = match mark {
            Mark::Upper => b.to_ascii_uppercase(),
            Mark::Lower => b.to_ascii_lowercase(),
            Mark::Char(c) => c,
        };
    }
}

/// `mark_snv`: marks the ALT bases that differ from REF where both have
/// one.
fn mark_snv(reference: &[u8], alt: &mut [u8], mark: Mark) {
    for (a, r) in alt.iter_mut().zip(reference) {
        if !a.eq_ignore_ascii_case(r) {
            mark_bytes(std::slice::from_mut(a), mark);
        }
    }
}

/// `mark_del`: `alt` (or, for `<DEL>`, REF) padded with `mark` to `rlen`.
fn mark_del(prefix: &[u8], rlen: i64, mark: u8, out: &mut Vec<u8>) {
    out.clear();
    out.extend_from_slice(&prefix[..prefix.len().min(rlen.max(0) as usize)]);
    out.resize(rlen.max(0) as usize, mark);
}

impl<'a> Edit<'a> {
    /// Starts on `src`, the region's reference beginning at `ori_pos`,
    /// building into `out` (cleared; its capacity is reused).
    pub fn new(src: &'a [u8], ori_pos: i64, mut out: Vec<u8>) -> Edit<'a> {
        out.clear();
        out.reserve(src.len());
        Edit {
            src,
            cur: 0,
            out,
            ori_pos,
            mod_off: 0,
            frz_pos: -1,
            frz_mod: -1,
            prev_base: 0,
            prev_base_pos: -1,
            prev_is_insert: false,
            empty_case: 0,
        }
    }

    /// The finished sequence.
    pub fn finish(mut self) -> Vec<u8> {
        self.out.extend_from_slice(&self.src[self.cur..]);
        self.out
    }

    /// Releases the buffer after a failure.
    pub fn abandon(self) -> Vec<u8> {
        self.out
    }

    fn len(&self) -> i64 {
        (self.out.len() + self.src.len() - self.cur) as i64
    }

    fn at(&self, i: i64) -> u8 {
        let i = i as usize;
        match i.checked_sub(self.out.len()) {
            None => self.out[i],
            Some(k) => self.src[self.cur + k],
        }
    }

    /// Whether `bytes` matches the buffer at `at`, ignoring case.
    fn matches(&self, at: i64, bytes: &[u8]) -> bool {
        let Ok(at) = usize::try_from(at) else {
            return false;
        };
        let (remaining, offset) = if at < self.out.len() {
            let n = bytes.len().min(self.out.len() - at);
            if !self.out[at..at + n].eq_ignore_ascii_case(&bytes[..n]) {
                return false;
            }
            (&bytes[n..], 0)
        } else {
            (bytes, at - self.out.len())
        };
        self.src[self.cur..]
            .get(offset..)
            .and_then(|suffix| suffix.get(..remaining.len()))
            .is_some_and(|suffix| suffix.eq_ignore_ascii_case(remaining))
    }

    /// Moves the boundary so that the prefix holds `[0, end)`.
    fn extend_to(&mut self, end: usize) {
        if let Some(n) = end.checked_sub(self.out.len()) {
            self.out
                .extend_from_slice(&self.src[self.cur..self.cur + n]);
            self.cur += n;
        }
    }

    /// Replaces `remove` bytes at `at` with `insert`.
    fn splice(&mut self, at: usize, remove: usize, insert: &[u8]) -> Result<()> {
        let len = self.len() as usize;
        if at > len || remove > len - at {
            return Err(Error::Unsupported(
                "consensus edit extends beyond the current sequence",
            ));
        }
        if at == 0 && insert.is_empty() && remove > 0 && remove == self.len() as usize {
            self.empty_case = self.at(0);
        }
        if at >= self.out.len() {
            self.extend_to(at);
            self.out.extend_from_slice(insert);
            self.cur += remove;
        } else {
            let in_out = remove.min(self.out.len() - at);
            self.out.splice(at..at + in_out, insert.iter().copied());
            self.cur += remove - in_out;
        }
        Ok(())
    }

    /// `apply_absent`: fills with `absent` from where the last record's
    /// claim ends up to `pos`.
    pub fn absent(&mut self, pos: i64, absent: u8) {
        let len = self.len();
        if len == 0 || pos <= self.frz_pos + 1 || pos <= self.ori_pos {
            return;
        }
        let at = pos
            .saturating_sub(self.ori_pos)
            .saturating_add(self.mod_off);
        let end = if pos != 0 && at < len { at } else { len };
        let begin = self.frz_mod.max(0);
        if end <= begin {
            return;
        }
        self.extend_to(end as usize);
        self.out[begin as usize..end as usize].fill(absent);
        // These bytes already have the query's absent character. A skipped
        // record need not fill them again; an applied edit or reference
        // claim establishes its own resume boundary.
        self.frz_mod = end;
    }

    /// `freeze_ref`: the reference over `[pos, pos + rlen)` stands.
    pub fn freeze(&mut self, pos: i64, rlen: i64) {
        if rlen <= 0 || self.frz_pos >= pos + rlen - 1 {
            return;
        }
        self.frz_pos = pos + rlen - 1;
        // A reference block may extend far beyond this loaded region;
        // its resume boundary only needs to remain beyond the output.
        self.frz_mod = (pos + rlen - self.ori_pos).saturating_add(self.mod_off);
    }

    /// Applies one record whose allele is chosen: `reference` is its REF,
    /// `alt` the bytes of the chosen allele (`Choice::Allele`), rewritten in
    /// place as consensus.c rewrites them; `scratch` is spare room.
    #[allow(clippy::too_many_arguments)]
    pub fn apply(
        &mut self,
        site: &Site,
        choice: Choice,
        reference: &[u8],
        alt: &mut Vec<u8>,
        scratch: &mut Vec<u8>,
        marks: &Marks,
        chain: Option<&mut Chain>,
    ) -> Result<()> {
        match choice {
            Choice::Skip => Ok(()),
            Choice::Ref => {
                if marks.absent.is_some() {
                    self.freeze(site.pos, site.rlen);
                }
                Ok(())
            }
            Choice::Missing => match marks.missing {
                Some(c) => self.missing(site, reference, c, marks.mark_snv),
                None => Ok(()),
            },
            Choice::Allele(ialt) => self.allele(site, ialt, reference, alt, scratch, marks, chain),
        }
    }

    /// The general path of `apply_variant` (lines 791-1058).
    #[allow(clippy::too_many_arguments)]
    fn allele(
        &mut self,
        site: &Site,
        ialt: usize,
        reference: &[u8],
        alt: &mut Vec<u8>,
        scratch: &mut Vec<u8>,
        marks: &Marks,
        chain: Option<&mut Chain>,
    ) -> Result<()> {
        if self.src.is_empty() {
            return Ok(());
        }
        let mut pos = site.pos;
        let mut rlen = site.rlen;
        let symbolic = alt.first() == Some(&b'<');
        let is_del = alt.eq_ignore_ascii_case(b"<DEL>");
        let is_ins = alt.len() >= 4 && alt[..4].eq_ignore_ascii_case(b"<INS");
        let is_gvcf = alt.eq_ignore_ascii_case(b"<*>") || alt.eq_ignore_ascii_case(b"<NON_REF>");

        // POS and REF include the base before an indel's event unless REF's
        // first base itself changes.
        let (var_type, var_n) = variant(reference, alt);
        let mut trim_beg = 0;
        let mut var_len = var_n.unwrap_or(0);
        if var_type & VCF_INDEL != 0 {
            trim_beg = (reference.first() == alt.first()) as i64;
        } else if var_type & VCF_OTHER != 0 && is_del {
            trim_beg = 1;
            var_len = 1 - rlen;
        } else if var_type & VCF_OTHER != 0 && is_ins {
            trim_beg = 1;
        }

        // Offsets into REF and ALT as consensus.c advances its pointers, and
        // ALT's end as it truncates.
        let mut ref_off: i64 = 0;
        let mut alt_off: i64 = 0;
        let mut alt_end = alt.len() as i64;

        // REF chosen over the end of an earlier deletion: what remains of it
        // after the claimed bases.
        if ialt == 0 && pos <= self.frz_pos && pos + rlen - 1 > self.frz_pos {
            let ntrim = self.frz_pos - pos + 1;
            let nref = reference.len() as i64;
            if ntrim >= nref {
                return Err(Error::invalid(
                    "consensus",
                    format!(
                        "failed to trim the overlapping variant at position {}",
                        pos + 1
                    ),
                ));
            }
            pos += ntrim;
            rlen -= ntrim;
            ref_off += ntrim;
            alt_off += ntrim;
        }
        // Starting on claimed bases is allowed only for an insertion right
        // after the claim that does not follow another insertion.
        if pos <= self.frz_pos
            && (pos < self.frz_pos || trim_beg == 0 || var_len == 0 || self.prev_is_insert)
        {
            return Ok(());
        }

        // Where the record starts in the buffer; before it when the record
        // starts before the region.
        let mut idx = pos - self.ori_pos + self.mod_off;
        if idx < 0 {
            let nref = reference.len() as i64 - ref_off;
            if symbolic {
                pos -= idx + 1;
                rlen += idx + 1;
                idx = -1;
            } else if nref < -idx {
                return Ok(());
            } else if alt_end - alt_off > -idx {
                pos -= idx;
                rlen += idx;
                ref_off -= idx;
                alt_off -= idx;
                idx = 0;
            } else {
                // Only REF reaches the region: keep one base before it.
                pos -= idx + 1;
                rlen += idx + 1;
                ref_off -= idx + 1;
                alt_off = alt_end - 1;
                idx = -1;
            }
        }
        let len = self.len();
        if rlen > len - idx {
            rlen = len - idx;
            if !symbolic && alt_end - alt_off > rlen {
                alt_end = alt_off + rlen;
            }
        }
        if idx > 0 && idx >= len {
            return Ok(());
        }

        let ref_at = |off: i64| &reference[(off.max(0) as usize).min(reference.len())..];
        let mut len_diff;
        if symbolic {
            if !is_del && !is_gvcf {
                return Err(Error::invalid(
                    "consensus",
                    format!(
                        "symbolic allele {:?} at position {}: only <DEL>, <*> and <NON_REF> apply",
                        String::from_utf8_lossy(alt),
                        pos + 1
                    ),
                ));
            }
            if is_gvcf {
                // A gVCF block: evidence for the reference throughout.
                self.freeze(pos, rlen);
                return Ok(());
            }
            if let Some(mark) = marks.mark_del {
                mark_del(ref_at(ref_off), rlen, mark, scratch);
                len_diff = 0;
            } else {
                // The base before the event stands.
                scratch.clear();
                scratch.extend_from_slice(ref_at(ref_off));
                scratch.truncate(1);
                len_diff = 1 - rlen;
            }
            std::mem::swap(alt, scratch);
            alt_off = 0;
            alt_end = alt.len() as i64;
        } else {
            let r = ref_at(ref_off);
            let n = rlen.max(0) as usize;
            if idx >= 0 && !(r.len() >= n && self.matches(idx, &r[..n])) {
                // REF's first base may be one an earlier deletion already
                // rewrote.
                let ok = self.prev_base_pos == pos
                    && r.first()
                        .is_some_and(|b| b.eq_ignore_ascii_case(&self.prev_base))
                    && (rlen == 1 || (r.len() >= n && self.matches(idx + 1, &r[1..n])));
                if !ok {
                    return Err(mismatch(pos, r));
                }
            }
            len_diff = (alt_end - alt_off) - rlen;
            if let (Some(mark), true) = (marks.mark_del, len_diff < 0) {
                mark_del(
                    &alt[alt_off as usize..alt_end as usize],
                    rlen,
                    mark,
                    scratch,
                );
                std::mem::swap(alt, scratch);
                alt_off = 0;
                alt_end = alt.len() as i64;
                len_diff = 0;
            }
        }
        let alen = alt_end - alt_off;
        let a = &mut alt[alt_off as usize..alt_end as usize];

        // The ALT takes the case of the reference base it starts on; any
        // byte that is not a lowercase letter counts as uppercase.
        let case = if len == 0 {
            self.empty_case
        } else {
            self.at(idx.max(0))
        };
        if case.is_ascii_lowercase() {
            a.make_ascii_lowercase();
        } else {
            a.make_ascii_uppercase();
        }
        let r = ref_at(ref_off);
        if let (Some(mark), true) = (marks.mark_ins, len_diff > 0) {
            let from = r.len().min(a.len());
            mark_bytes(&mut a[from..], mark);
        }
        if let Some(mark) = marks.mark_snv {
            mark_snv(r, a, mark);
        }

        if len_diff <= 0 {
            if idx + rlen >= 1 {
                self.prev_base = self.at(idx + rlen - 1);
            }
            self.prev_base_pos = pos + rlen - 1;
            self.prev_is_insert = false;
            self.frz_mod = idx + alen;
            let start = (idx + trim_beg).max(0);
            let from = ((start - idx) as usize).min(a.len());
            self.splice(
                start as usize,
                (idx + rlen - start).max(0) as usize,
                &a[from..],
            )?;
        } else {
            self.prev_is_insert = true;
            self.prev_base_pos = pos;
            // Leading bases the insertion keeps do not overwrite an earlier
            // record's base at the same position.
            let mut ibeg = 0;
            while ibeg < alen
                && r.get(ibeg as usize) == Some(&a[ibeg as usize])
                && pos + ibeg <= self.prev_base_pos
            {
                ibeg += 1;
            }
            let span_end = idx + rlen;
            if idx + ibeg <= span_end {
                let start = (idx + ibeg).max(0);
                let from = (start - idx) as usize;
                self.splice(start as usize, (span_end - start) as usize, &a[from..])?;
            } else {
                // The insertion opens a gap before what it writes.
                let mut fill = vec![0; (idx + ibeg - span_end) as usize];
                fill.extend_from_slice(&a[ibeg as usize..]);
                self.splice(span_end as usize, 0, &fill)?;
            }
            self.frz_mod = idx + alen - ibeg + 1;
        }

        if let (Some(chain), true) = (chain, len_diff != 0) {
            let same_first =
                matches!((r.first(), a.first()), (Some(x), Some(y)) if x.eq_ignore_ascii_case(y));
            if same_first {
                chain.push_gap(pos + 1, rlen - 1, pos + 1 + self.mod_off, alen - 1);
            } else {
                chain.push_gap(pos, rlen, pos + self.mod_off, alen);
            }
        }
        self.mod_off += len_diff;
        self.frz_pos = pos + rlen - 1;
        Ok(())
    }

    /// A missing genotype with `-M` (lines 766-788): consensus.c rewrites
    /// the alleles to REF and the missing character (over the whole span of
    /// a gVCF block) and applies that.
    fn missing(&mut self, site: &Site, reference: &[u8], c: u8, mark: Option<Mark>) -> Result<()> {
        if self.len() == 0 {
            return Ok(());
        }
        let mut pos = site.pos;
        let mut rlen = if site.ref_block { site.rlen } else { 1 };
        let mut idx = pos - self.ori_pos + self.mod_off;
        if pos <= self.frz_pos {
            return Ok(());
        }
        if idx < 0 {
            if rlen <= -idx {
                return Ok(());
            }
            pos -= idx;
            rlen += idx;
            idx = 0;
        }
        let len = self.len();
        if idx > 0 && idx >= len {
            return Ok(());
        }
        rlen = rlen.min(len - idx);
        if rlen <= 0 {
            return Ok(());
        }
        let mut written: smallvec::SmallVec<[u8; 16]> = smallvec::SmallVec::new();
        let mut original: smallvec::SmallVec<[u8; 16]> = smallvec::SmallVec::new();
        if site.ref_block {
            original.extend((idx..idx + rlen).map(|i| self.at(i)));
        } else {
            let base = reference.first().copied().unwrap_or(0);
            // With the alleles rewritten, htslib takes the span from
            // INFO/END again: beyond one base the one-base REF mismatches.
            let span = match site.info_end {
                Some(end) if end > site.pos => end - site.pos,
                _ => 1,
            };
            if span.min(len - idx) > 1 {
                return Err(mismatch(pos, &[base]));
            }
            let ok = self.at(idx).eq_ignore_ascii_case(&base)
                || (self.prev_base_pos == pos && base.eq_ignore_ascii_case(&self.prev_base));
            if !ok {
                return Err(mismatch(pos, &[base]));
            }
            original.push(base);
        }
        let last = self.at(idx + rlen - 1);
        written.resize(rlen as usize, c);
        if self.at(idx).is_ascii_lowercase() {
            written.make_ascii_lowercase();
        } else {
            written.make_ascii_uppercase();
        }
        if let Some(mark) = mark {
            mark_snv(&original, &mut written, mark);
        }
        self.splice(idx as usize, rlen as usize, &written)?;
        self.prev_base = last;
        self.prev_base_pos = pos + rlen - 1;
        self.prev_is_insert = false;
        self.frz_mod = idx + rlen;
        self.frz_pos = pos + rlen - 1;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::Edit;

    fn check_matches(edit: &Edit<'_>) {
        let mut logical = edit.out.clone();
        logical.extend_from_slice(&edit.src[edit.cur..]);
        for at in -2..logical.len() as i64 + 3 {
            for len in 0..logical.len() + 3 {
                let mut bytes = vec![b'A'; len];
                if let Ok(start) = usize::try_from(at)
                    && let Some(span) = logical.get(start..).and_then(|s| s.get(..len))
                {
                    bytes.copy_from_slice(span);
                    bytes.make_ascii_lowercase();
                }
                for mismatch in [false, true] {
                    if mismatch && !bytes.is_empty() {
                        bytes[len / 2] ^= 0x80;
                    }
                    let want = usize::try_from(at)
                        .ok()
                        .and_then(|at| logical.get(at..))
                        .and_then(|s| s.get(..len))
                        .is_some_and(|span| span.eq_ignore_ascii_case(&bytes));
                    assert_eq!(edit.matches(at, &bytes), want, "at={at} bytes={bytes:?}");
                }
            }
        }
        assert!(!edit.matches(i64::MAX, b"A"));
        assert!(!edit.matches(i64::MIN, b""));
    }

    #[test]
    fn reference_matching_handles_both_initialized_spans() {
        let mut edit = Edit::new(b"ACgtN\xffTtaC", 0, Vec::new());
        check_matches(&edit);
        edit.splice(2, 3, b"agCT").unwrap();
        check_matches(&edit);
        edit.splice(1, 2, b"T").unwrap();
        check_matches(&edit);
        edit.splice(0, edit.len() as usize, b"").unwrap();
        check_matches(&edit);
    }
}
