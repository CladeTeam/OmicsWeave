//! Reading: line batches parsed on a worker pool, merged in order.

use std::io::BufRead;
use std::path::Path;

use omicsweave_core::Result;
use omicsweave_core::exec::{Exec, Ordered, Pool};

use crate::alias::Aliases;
use crate::annotation::Annotation;
use crate::build::{BuildOptions, build};
use crate::input::{LineBatch, LineBatches, open_bytes, open_file};
use crate::model::{AttributePolicy, Dialect};
use crate::parse::{Keep, Raw, parse_batch, sniff_dialect};

/// How to read an annotation.
#[derive(Debug, Clone)]
pub struct ReadOptions {
    /// The dialect; `None` tells it from the text (a `##gff-version 3`
    /// directive, else the attribute syntax of the first record).
    pub dialect: Option<Dialect>,
    pub attributes: AttributePolicy,
    /// Seqids renamed at read (`1` → `chr1`, ...).
    pub aliases: Option<Aliases>,
    /// Fail on the first problem instead of reporting it.
    pub strict: bool,
    /// Parsing threads; `0` uses the available parallelism (at most 8).
    pub threads: usize,
    /// Bytes of text per parse batch.
    pub batch_bytes: usize,
}

impl Default for ReadOptions {
    fn default() -> ReadOptions {
        ReadOptions {
            dialect: None,
            attributes: AttributePolicy::All,
            aliases: None,
            strict: false,
            threads: 0,
            batch_bytes: 4 << 20,
        }
    }
}

impl ReadOptions {
    fn threads(&self) -> usize {
        if self.threads > 0 {
            return self.threads;
        }
        std::thread::available_parallelism()
            .map_or(1, |n| n.get())
            .min(8)
    }
}

impl Annotation {
    /// Reads a GTF or GFF3 file, plain, gzip or BGZF (told by content).
    pub fn read(path: impl AsRef<Path>, opts: &ReadOptions) -> Result<Annotation> {
        let threads = opts.threads();
        let pool = (threads > 1).then(|| Pool::new(threads));
        let exec = match &pool {
            Some(p) => Exec::Pool {
                pool: p,
                limit: 2 * threads,
            },
            None => Exec::Inline,
        };
        let text = open_file(path.as_ref(), exec)?;
        read_text(text, opts, exec)
    }

    /// Parses annotation text held in memory (plain, gzip or BGZF).
    pub fn parse(bytes: &[u8], opts: &ReadOptions) -> Result<Annotation> {
        let threads = opts.threads();
        let pool = (threads > 1).then(|| Pool::new(threads));
        let exec = match &pool {
            Some(p) => Exec::Pool {
                pool: p,
                limit: 2 * threads,
            },
            None => Exec::Inline,
        };
        read_text(open_bytes(bytes), opts, exec)
    }
}

fn read_text(text: impl BufRead, opts: &ReadOptions, exec: Exec<'_>) -> Result<Annotation> {
    let mut batches = LineBatches::new(text, opts.batch_bytes);
    let first = batches.next_batch()?;
    let dialect = opts.dialect.unwrap_or_else(|| {
        first
            .as_ref()
            .map_or(Dialect::Gtf, |b| sniff_dialect(&b.text))
    });
    let keep = Keep::new(&opts.attributes);
    let raw = match exec {
        Exec::Inline => {
            let mut raw = Raw::new();
            dispatch(&mut batches, first, dialect, &keep, exec, &mut |part| {
                raw.append(part)
            })?;
            raw
        }
        // Merging moves each batch's strings into the shared pool, which
        // only one thread can do; it runs beside reading and dispatching.
        Exec::Pool { .. } => std::thread::scope(|scope| -> Result<Raw> {
            let (tx, rx) = std::sync::mpsc::sync_channel::<Raw>(4);
            let merger = scope.spawn(move || -> Result<Raw> {
                let mut raw = Raw::new();
                for part in rx {
                    raw.append(part)?;
                }
                Ok(raw)
            });
            let sent = dispatch(&mut batches, first, dialect, &keep, exec, &mut |part| {
                tx.send(part)
                    .map_err(|_| omicsweave_core::Error::invalid("annotation merge", "stopped"))
            });
            drop(tx);
            let merged = merger.join().unwrap_or_else(|_| {
                Err(omicsweave_core::Error::invalid(
                    "annotation merge",
                    "panicked",
                ))
            });
            // A merge error is the cause of a failed send.
            let raw = merged?;
            sent?;
            Ok(raw)
        })?,
    };
    build(
        raw,
        BuildOptions {
            dialect,
            policy: opts.attributes.clone(),
            aliases: opts.aliases.as_ref(),
            strict: opts.strict,
        },
    )
}

/// Parses every batch on `exec` and hands the results to `sink` in input
/// order, stopping after a batch that ends at a `##FASTA` directive.
fn dispatch(
    batches: &mut LineBatches<impl BufRead>,
    mut pending: Option<LineBatch>,
    dialect: Dialect,
    keep: &Keep,
    exec: Exec<'_>,
    sink: &mut dyn FnMut(Raw) -> Result<()>,
) -> Result<()> {
    let mut queue: Box<dyn Ordered<Raw>> = exec.queue();
    let mut stop = false;
    let mut take = |part: Raw, stop: &mut bool| -> Result<()> {
        if *stop {
            return Ok(());
        }
        *stop = part.fasta;
        sink(part)
    };
    while !stop {
        let batch = match pending.take() {
            Some(b) => b,
            None => match batches.next_batch()? {
                Some(b) => b,
                None => break,
            },
        };
        while queue.full() {
            if let Some(done) = queue.next(true)? {
                take(done, &mut stop)?;
                if stop {
                    return Ok(());
                }
            }
        }
        let keep = keep.clone();
        queue.submit(Box::new(move || {
            parse_batch(&batch.text, batch.first_line, dialect, &keep)
        }))?;
        while let Some(done) = queue.next(false)? {
            take(done, &mut stop)?;
            if stop {
                return Ok(());
            }
        }
    }
    while let Some(done) = queue.next(true)? {
        take(done, &mut stop)?;
        if stop {
            return Ok(());
        }
    }
    Ok(())
}
