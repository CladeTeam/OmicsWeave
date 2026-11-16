//! Public BGZF behavior, independent wire fixtures and fault injection.
mod support;

#[path = "cases/blocks.rs"]
mod blocks;
#[path = "cases/index.rs"]
mod index;
#[path = "cases/index_compression.rs"]
mod index_compression;
#[path = "cases/reader.rs"]
mod reader;
#[path = "cases/scanner.rs"]
mod scanner;
#[path = "cases/seek_cache.rs"]
mod seek_cache;
#[path = "cases/writer.rs"]
mod writer;
