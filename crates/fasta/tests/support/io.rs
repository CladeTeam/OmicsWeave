//! Caller-side compositions of staged reads, independent from the reference C libraries.

use omicsweave_fasta::{
    Directory, FaiEntry, Gzi, Inflater, MAX_BLOCK_LEN, Parsed, Record, TWOBIT_HEADER_LEN,
    TwoBitHeader,
};

/// Bases `[start, end)` of a plain FASTA held in `text`.
pub fn fetch_plain(entry: &FaiEntry, text: &[u8], start: u32, end: u32) -> Vec<u8> {
    let range = entry.text_range(start, end).unwrap();
    let mut out = vec![0u8; (end - start) as usize];
    entry
        .decode(
            start,
            &text[range.start as usize..range.end as usize],
            &mut out,
        )
        .unwrap();
    out
}

/// A bgzip FASTA held in memory with its `.gzi`.
pub struct Bgzipped<'a> {
    pub bytes: &'a [u8],
    pub gzi: Gzi,
    inflater: Inflater,
}

impl<'a> Bgzipped<'a> {
    pub fn new(bytes: &'a [u8], gzi: Gzi) -> Bgzipped<'a> {
        Bgzipped {
            bytes,
            gzi,
            inflater: Inflater::new(),
        }
    }

    fn read(&self, range: std::ops::Range<u64>) -> &'a [u8] {
        &self.bytes[range.start as usize..(range.end as usize).min(self.bytes.len())]
    }

    /// One contiguous read through [`Gzi::span`].
    pub fn fetch_span(&mut self, entry: &FaiEntry, start: u32, end: u32) -> Vec<u8> {
        let text = entry.text_range(start, end).unwrap();
        let span = self
            .gzi
            .span(text.clone(), self.bytes.len() as u64)
            .expect("non-empty");
        let mut inflated = vec![0u8; span.max_text_len];
        let n = self
            .inflater
            .inflate_blocks(
                self.read(span.compressed.clone()),
                span.blocks,
                &mut inflated,
            )
            .unwrap();
        let skip = (text.start - span.text_start) as usize;
        let mut out = vec![0u8; (end - start) as usize];
        entry.decode(start, &inflated[skip..n], &mut out).unwrap();
        out
    }

    /// Block by block through [`Gzi::blocks`], as a caching loader would.
    pub fn fetch_blocks(&mut self, entry: &FaiEntry, start: u32, end: u32) -> Vec<u8> {
        let text = entry.text_range(start, end).unwrap();
        let blocks: Vec<_> = self
            .gzi
            .blocks(text.clone(), self.bytes.len() as u64)
            .collect();
        let mut inflated = Vec::new();
        let mut block = vec![0u8; MAX_BLOCK_LEN];
        for b in &blocks {
            assert_eq!(
                b.text_start,
                blocks[0].text_start + inflated.len() as u64,
                "blocks are contiguous"
            );
            let (_, n) = self
                .inflater
                .inflate_block(self.read(b.compressed.clone()), &mut block)
                .unwrap();
            if let Some(len) = b.text_len {
                assert_eq!(len, n as u64, "indexed block length");
            }
            inflated.extend_from_slice(&block[..n]);
        }
        let skip = (text.start - blocks[0].text_start) as usize;
        let mut out = vec![0u8; (end - start) as usize];
        entry.decode(start, &inflated[skip..], &mut out).unwrap();
        out
    }
}

/// A 2bit file opened through the crate's staged parsing.
pub struct OpenTwoBit<'a> {
    pub bytes: &'a [u8],
    pub directory: Directory,
    pub records: Vec<Record>,
    /// Read rounds each record header took.
    pub rounds: Vec<usize>,
}

impl<'a> OpenTwoBit<'a> {
    pub fn open(bytes: &'a [u8]) -> OpenTwoBit<'a> {
        let header = TwoBitHeader::parse(&bytes[..TWOBIT_HEADER_LEN]).unwrap();
        let bound = header.directory_bound(bytes.len() as u64);
        let directory =
            match Directory::parse(&header, &bytes[bound.start as usize..bound.end as usize])
                .unwrap()
            {
                Parsed::Done(d) => d,
                Parsed::Need(n) => panic!("directory bound too small: need {n}"),
            };
        let (mut records, mut rounds) = (Vec::new(), Vec::new());
        for (_, offset) in directory.entries() {
            let at = *offset as usize;
            let mut want = Record::FIRST_READ;
            let mut round = 1;
            let record = loop {
                match Record::parse(&bytes[at..at + want], *offset).unwrap() {
                    Parsed::Done(r) => break r,
                    Parsed::Need(n) => {
                        assert!(n > want, "a round must ask for more bytes");
                        want = n;
                        round += 1;
                    }
                }
            };
            records.push(record);
            rounds.push(round);
        }
        OpenTwoBit {
            bytes,
            directory,
            records,
            rounds,
        }
    }

    pub fn fetch(&self, i: usize, start: u32, end: u32, mask: bool) -> Vec<u8> {
        let record = &self.records[i];
        let range = record.packed_range(start, end).unwrap();
        let mut out = vec![0u8; (end - start) as usize];
        record
            .decode(
                start,
                &self.bytes[range.start as usize..range.end as usize],
                &mut out,
                mask,
            )
            .unwrap();
        out
    }
}

/// Observable positional source for read budgets and recoverable failure checks.
pub struct Source {
    pub bytes: std::cell::RefCell<Vec<u8>>,
    pub calls: std::cell::Cell<usize>,
    pub total: std::cell::Cell<usize>,
    pub largest: std::cell::Cell<usize>,
    pub requests: std::cell::RefCell<Vec<(u64, usize)>>,
    pub fail: std::cell::Cell<Option<usize>>,
}
impl Source {
    pub fn new(bytes: Vec<u8>) -> Self {
        Self {
            bytes: bytes.into(),
            calls: 0.into(),
            total: 0.into(),
            largest: 0.into(),
            requests: Vec::new().into(),
            fail: None.into(),
        }
    }
    pub fn reset(&self) {
        self.calls.set(0);
        self.total.set(0);
        self.largest.set(0);
        self.requests.borrow_mut().clear();
    }
}
impl omicsweave_core::ReadAt for Source {
    fn len(&self) -> u64 {
        self.bytes.borrow().len() as u64
    }
    fn read_exact_at(&self, out: &mut [u8], offset: u64) -> omicsweave_core::Result<()> {
        let call = self.calls.get() + 1;
        self.calls.set(call);
        self.total.set(self.total.get() + out.len());
        self.largest.set(self.largest.get().max(out.len()));
        self.requests.borrow_mut().push((offset, out.len()));
        if self.fail.get() == Some(call) {
            out.fill(0);
            return Err(omicsweave_core::Error::invalid(
                "source",
                "injected I/O failure",
            ));
        }
        self.bytes.borrow().as_slice().read_exact_at(out, offset)
    }
}
