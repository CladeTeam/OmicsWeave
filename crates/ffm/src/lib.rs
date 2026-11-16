//! FFM (Factor Frequency Model): a binary format for motif collections.
//!
//! A motif is a window of `width` positions scored by local factors. Each
//! factor covers a few strictly increasing positions and holds a table of
//! non-negative frequencies over the motif's local codes: the symbols it
//! models individually (its support) plus, when the support is a strict
//! subset of the alphabet, one rest category for every other symbol. A
//! position frequency matrix is the case of one single-position joint factor
//! per position; dinucleotide, inhomogeneous Markov, energy and k-mer models
//! use higher-order factors. Bytes outside the alphabet are [`INVALID`] and
//! never index a table.
//!
//! The file stores models faithfully (counts or frequencies); turning them
//! into scores with backgrounds and pseudocounts belongs to the matching
//! layer. Layout (see `docs/formats/ffm.md`):
//!
//! - a fixed 128-byte header;
//! - metadata back to back: alphabet, annotation schema, IDs, block
//!   directory and a per-motif index of block number and stream offsets;
//! - blocks of consecutive records, each a structure stream (`u32`), a value
//!   stream (`f64`) and an annotation stream (bytes), compressed separately
//!   with pcodec or zstd.
//!
//! [`FfmReader::open`] reads the header and the metadata in two reads and
//! builds an in-memory ID table. A motif is located in O(1) (expected O(1)
//! by ID); its block is read whole, and the annotation stream is decoded
//! only when asked for. Random access reuses decoded blocks through an
//! explicitly attached caller-owned [`omicsweave_core::cache::BlockCache`];
//! a sequential [`FfmReader::motifs`] traversal holds only its current block.
//! [`FfmWriter`] encodes blocks as motifs are pushed and writes the file on
//! [`FfmWriter::finish`].
//!
//! Small collections can be loaded with [`MotifSet::load`]. It validates all
//! records and packs tables, positions and annotations into shared arenas.
//! Lookups then return borrowed views without I/O, locks or allocation.
//!
//! ```
//! use omicsweave_ffm::{Alphabet, MotifSet, MotifSetBuilder, Schema, WriterOptions};
//!
//! let mut builder = MotifSetBuilder::new(Alphabet::dna(), Schema::default());
//! builder.push_pfm("TF1", &[1.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0, 1.0], None, &[])?;
//! let motifs = builder.finish();
//! let bytes = motifs.write_to(Vec::new(), WriterOptions::default())?;
//! let motifs = MotifSet::load(bytes)?;
//! let entry = motifs.get_by_id("TF1").unwrap();
//! let model = entry.model();
//! assert_eq!(model.pfm().unwrap().row(1), Some(&[4.0, 3.0, 2.0, 1.0][..]));
//! assert_eq!(model.factor(0).unwrap().positions, &[0]);
//! # Ok::<(), omicsweave_ffm::Error>(())
//! ```

mod alphabet;
mod block;
mod collection;
mod error;
mod layout;
mod meta;
mod model;
mod options;
mod reader;
mod schema;
mod stream;
mod view;
mod writer;

pub use alphabet::{Alphabet, INVALID, MAX_SYMBOLS};
pub use collection::{AnnotationView, MotifEntry, MotifSet, MotifSetBuilder};
pub use error::{Error, Result};
pub use layout::{MAGIC, MAJOR, MINOR};
pub use model::{
    Factor, FactorKind, Gaps, MAX_MOTIFS, MAX_ORDER, MAX_TABLE, Motif, Support, arity, table_len,
};
pub use options::ReaderOptions;
pub use reader::{AnnotationRef, FfmReader, MotifRef};
pub use schema::{AnnotationValue, Key, KeyType, Schema};
pub use stream::{Codec, Levels};
pub use view::{FactorView, Factors, GapsView, MotifEncoding, MotifView, PfmView};
pub use writer::{FfmWriter, WriterOptions};
