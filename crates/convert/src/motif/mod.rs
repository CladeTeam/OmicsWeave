//! Motif collections into FFM: every source format is parsed into factor
//! models over one alphabet and written as a single `.ffm`.
//!
//! Sources are a list of paths (as a shell expands `pwms/*.txt`); together,
//! in the order given and each file once, they form one collection.
//! Values are stored as the source states them: counts stay counts and
//! probabilities stay probabilities. Each motif keeps its ID and one
//! annotation, `name` (the alternate name, TF or display name the source
//! gives, empty when none). GLAM2 alignments become gapped profiles.

mod bamm;
mod cisbp;
mod glam2;
mod jaspar;
mod meme;
mod probound;
mod rows;
mod transfac;

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;
use omicsweave_core::{Error, FileSource, Result};
use omicsweave_ffm::{
    Alphabet, Factor, FactorKind, FfmWriter, Key, KeyType, Motif, MotifSet, ReaderOptions, Schema,
    Support, WriterOptions,
};

use crate::budget::{self, Budget, MemoryPlan};
use crate::conversion::{Context, Report};
use crate::format::Format;
use crate::pipeline::input::with_path;
use crate::pipeline::output::Outputs;

/// Motif source formats, in [`Format`] order.
pub(crate) const SOURCES: [Format; 14] = [
    Format::Meme,
    Format::Jaspar,
    Format::JasparPfm,
    Format::Transfac,
    Format::Homer,
    Format::CisBp,
    Format::Pcm,
    Format::Pwm,
    Format::Dpcm,
    Format::Dpwm,
    Format::Bamm,
    Format::ProBound,
    Format::Glam2,
    Format::Ffm,
];

/// One motif read from a source.
pub(crate) struct Entry {
    pub id: String,
    pub name: String,
    pub motif: Motif,
}

/// What one source file holds.
pub(crate) struct Parsed {
    pub alphabet: Alphabet,
    pub entries: Vec<Entry>,
    /// Records passed over without failing (empty matrices).
    pub skipped: Vec<String>,
}

impl Parsed {
    fn dna(entries: Vec<Entry>) -> Parsed {
        Parsed {
            alphabet: Alphabet::dna(),
            entries,
            skipped: Vec::new(),
        }
    }
}

enum Source {
    Text(Box<Parsed>),
    Resident(MotifSet),
}

impl Source {
    fn alphabet(&self) -> &Alphabet {
        match self {
            Source::Text(parsed) => &parsed.alphabet,
            Source::Resident(set) => set.alphabet(),
        }
    }
}

/// An error at `line` (1-based; 0 for the whole file) of `path`.
pub(crate) fn bad(path: &Path, line: usize, detail: impl std::fmt::Display) -> Error {
    let at = if line == 0 {
        path.display().to_string()
    } else {
        format!("{}:{line}", path.display())
    };
    Error::invalid("motif", format!("{at}: {detail}"))
}

/// Whitespace-separated numbers of a line.
pub(crate) fn numbers(line: &str) -> Option<Vec<f64>> {
    line.split_whitespace()
        .map(|w| w.parse::<f64>().ok())
        .collect()
}

/// `>id name...` split into the ID and the rest.
pub(crate) fn header(line: &str) -> (String, String) {
    let rest = line.trim_start_matches('>').trim();
    match rest.split_once(char::is_whitespace) {
        Some((id, name)) => (id.to_string(), name.trim().to_string()),
        None => (rest.to_string(), String::new()),
    }
}

/// The file name without its directory, `.gz` and last extension.
pub(crate) fn stem(path: &Path) -> String {
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    match name.rsplit_once('.') {
        Some((s, _)) if !s.is_empty() => s.to_string(),
        _ => name.to_string(),
    }
}

/// A motif of one single-position joint factor per row. `rows[p]` holds
/// the values of position `p` in alphabet order. When the rows are counts
/// (some position sums past 1) `sites` is the largest position sum.
pub(crate) fn per_position(rows: Vec<Vec<f64>>, sites: Option<f64>) -> Motif {
    let sites = sites.or_else(|| count_sites(&rows));
    Motif {
        width: rows.len() as u32,
        support: Support::All,
        sites,
        factors: rows
            .into_iter()
            .enumerate()
            .map(|(p, table)| Factor {
                kind: FactorKind::Joint,
                positions: vec![p as u32],
                table,
            })
            .collect(),
        gaps: None,
    }
}

/// The largest table sum, when the tables hold counts rather than
/// frequencies (some sum is clearly above 1).
pub(crate) fn count_sites(tables: &[Vec<f64>]) -> Option<f64> {
    let max = tables
        .iter()
        .map(|t| t.iter().sum::<f64>())
        .fold(0.0f64, f64::max);
    (max > 1.0 + 1e-3).then_some(max)
}

/// Potentials of log weights: a weight `w` scores as `ln(e^w) = w`.
pub(crate) fn potentials(weights: Vec<f64>) -> Vec<f64> {
    weights.into_iter().map(f64::exp).collect()
}

/// The text of a source, inflated when gzip compressed. Bytes that are not
/// UTF-8 (free text in some databases) are replaced.
pub(crate) fn read_text(path: &Path, budget: &Budget) -> Result<(String, u64)> {
    let file = File::open(path).map_err(|e| with_path(e, path))?;
    let len = file.metadata()?.len();
    let mut input = BufReader::with_capacity(64 << 10, file);
    let gzip = input.fill_buf()?.starts_with(&[0x1f, 0x8b]);
    let input: Box<dyn Read> = if gzip {
        Box::new(MultiGzDecoder::new(input))
    } else {
        Box::new(input)
    };
    let mut bytes = Vec::new();
    input
        .take(budget.memory.saturating_add(1))
        .read_to_end(&mut bytes)
        .map_err(|e| with_path(e, path))?;
    budget::require(budget, bytes.len() as u64, "decoded motif text")?;
    let text = String::from_utf8(bytes)
        .unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned());
    budget::require(budget, text.len() as u64, "decoded motif text")?;
    Ok((text, len))
}

/// The source files, in the order given, each once. Every one must exist.
pub(crate) fn expand(paths: &[PathBuf]) -> Result<Vec<PathBuf>> {
    if paths.is_empty() {
        return Err(Error::invalid(
            "motif",
            "at least one source file is required",
        ));
    }
    let mut seen = HashSet::new();
    let mut files = Vec::with_capacity(paths.len());
    for path in paths {
        if !path.is_file() {
            return Err(Error::Io {
                kind: std::io::ErrorKind::NotFound,
                detail: format!("{} does not exist", path.display()),
            });
        }
        if seen.insert(path.clone()) {
            files.push(path.clone());
        }
    }
    Ok(files)
}

/// The format of every source file, which must be one. Files are told apart
/// in parallel: collections of thousands of small files are bound by
/// file-system latency.
pub fn detect_all(paths: &[PathBuf]) -> Result<Format> {
    let files = expand(paths)?;
    let threads = std::thread::available_parallelism()
        .map_or(1, |n| n.get())
        .clamp(1, 16)
        .min(files.len());
    let chunk = files.len().div_ceil(threads);
    let formats: Vec<Format> = std::thread::scope(|scope| {
        let handles: Vec<_> = files
            .chunks(chunk)
            .map(|part| {
                scope.spawn(move || {
                    part.iter()
                        .map(|f| Format::detect(f))
                        .collect::<Result<Vec<_>>>()
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("detection thread"))
            .collect::<Result<Vec<Vec<_>>>>()
    })?
    .into_iter()
    .flatten()
    .collect();
    let mut found: Option<(Format, &Path)> = None;
    for (f, &format) in files.iter().zip(&formats) {
        match found {
            None => found = Some((format, f)),
            Some((g, first)) if g != format => {
                return Err(Error::invalid(
                    "motif",
                    format!(
                        "{} holds {g} but {} holds {format}; name the format",
                        first.display(),
                        f.display()
                    ),
                ));
            }
            Some(_) => {}
        }
    }
    Ok(found.expect("at least one file").0)
}

/// Converts motif collections (one or more files) into one FFM.
#[derive(Debug, Clone, PartialEq)]
pub struct MotifToFfm {
    /// The source files.
    pub sources: Vec<PathBuf>,
    pub dst: PathBuf,
    /// The source format every matched file holds.
    pub from: Format,
}

impl MotifToFfm {
    pub fn new(src: impl Into<PathBuf>, dst: impl Into<PathBuf>, from: Format) -> MotifToFfm {
        MotifToFfm {
            sources: vec![src.into()],
            dst: dst.into(),
            from,
        }
    }

    pub(crate) fn plan(&self, _budget: &Budget) -> Result<MemoryPlan> {
        let files = expand(&self.sources)?;
        if files.iter().any(|f| f == &self.dst) {
            return Err(Error::invalid(
                "convert",
                "source and destination are the same file",
            ));
        }
        // Whole-file parsers and compressed output remain resident. Source
        // sizes estimate this payload; gzip expansion is bounded while read.
        let fixed = files
            .iter()
            .try_fold(0u64, |size, f| {
                let len = std::fs::metadata(f).map_err(|e| with_path(e, f))?.len();
                size.checked_add(len)
                    .ok_or_else(|| Error::invalid("budget", "motif source sizes overflow"))
            })?
            .saturating_mul(4);
        Ok(MemoryPlan {
            fixed,
            pipeline: 0,
            buffer: 0,
        })
    }

    pub(crate) fn run(&self, ctx: &Context<'_>) -> Result<Report> {
        let files = expand(&self.sources)?;
        let total = files.iter().try_fold(0u64, |size, f| {
            let len = std::fs::metadata(f).map_err(|e| with_path(e, f))?.len();
            size.checked_add(len)
                .ok_or_else(|| Error::invalid("budget", "motif source sizes overflow"))
        })?;
        let schema = Schema::new(vec![Key::new("name", KeyType::String)])?;
        let mut outputs = Outputs::new();
        let mut writer: Option<(FfmWriter<File>, PathBuf)> = None;
        let mut read = 0u64;
        let mut skipped = Vec::new();
        for path in &files {
            ctx.progress("read", read, total);
            let (parsed, len) = self.parse(path, ctx.budget)?;
            read += len;
            let (w, first) = match &mut writer {
                Some(w) => w,
                None => writer.insert((
                    FfmWriter::new(
                        outputs.create(&self.dst, Format::Ffm)?,
                        parsed.alphabet().clone(),
                        schema.clone(),
                        WriterOptions::default(),
                    )?,
                    path.clone(),
                )),
            };
            if w.alphabet() != parsed.alphabet() {
                return Err(Error::invalid(
                    "motif",
                    format!(
                        "{} uses alphabet {:?} but {} uses {:?}; convert them separately",
                        path.display(),
                        parsed.alphabet().name(),
                        first.display(),
                        w.alphabet().name()
                    ),
                ));
            }
            match parsed {
                Source::Text(parsed) => {
                    skipped.extend(parsed.skipped);
                    for e in parsed.entries {
                        w.push(&e.id, &e.motif, &[&e.name])
                            .map_err(|err| bad(path, 0, format!("motif {:?}: {err}", e.id)))?;
                    }
                }
                Source::Resident(set) => {
                    let name_key = set.schema().index("name");
                    for entry in set.iter() {
                        let annotation = entry.annotation();
                        let name = name_key.and_then(|k| annotation.get(k)).unwrap_or("");
                        w.push_view(entry.id(), entry.model(), &[name])
                            .map_err(|err| {
                                bad(path, 0, format!("motif {:?}: {err}", entry.id()))
                            })?;
                    }
                }
            }
        }
        // Thousands of empty shells (CIS-BP) are one note, not thousands.
        if skipped.len() > 5 {
            ctx.warn(format!(
                "{} motifs were skipped, such as: {}",
                skipped.len(),
                skipped[..3].join("; ")
            ));
        } else {
            skipped.into_iter().for_each(|s| ctx.warn(s));
        }
        let (writer, _) = writer.expect("at least one file");
        let records = writer.len() as u64;
        ctx.progress("write", read, total);
        drop(writer.finish()?);
        Ok(Report {
            files: outputs.commit()?,
            records,
            bytes_read: read,
            spilled: 0,
            elapsed: Default::default(),
            warnings: Vec::new(),
        })
    }

    fn parse(&self, path: &Path, budget: &Budget) -> Result<(Source, u64)> {
        if self.from == Format::Ffm {
            return read_ffm(path, budget).map(|(set, len)| (Source::Resident(set), len));
        }
        let (text, len) = read_text(path, budget)?;
        let parsed = match self.from {
            Format::Meme => meme::parse(path, &text)?,
            Format::Jaspar => jaspar::parse_brackets(path, &text)?,
            Format::JasparPfm => jaspar::parse_raw(path, &text)?,
            Format::Transfac => transfac::parse(path, &text)?,
            Format::Homer => rows::parse_homer(path, &text)?,
            Format::CisBp => cisbp::parse(path, &text)?,
            Format::Pcm => rows::parse(path, &text, rows::Kind::Counts)?,
            Format::Pwm => rows::parse(path, &text, rows::Kind::Weights)?,
            Format::Dpcm => rows::parse(path, &text, rows::Kind::DiCounts)?,
            Format::Dpwm => rows::parse(path, &text, rows::Kind::DiWeights)?,
            Format::Bamm => bamm::parse(path, &text)?,
            Format::ProBound => probound::parse(path, &text)?,
            Format::Glam2 => glam2::parse(path, &text)?,
            other => {
                return Err(Error::invalid(
                    "motif",
                    format!("{other} is not a motif format"),
                ));
            }
        };
        Ok((Source::Text(Box::new(parsed)), len))
    }
}

/// Reads all models and annotations once into shared resident arenas.
fn read_ffm(path: &Path, budget: &Budget) -> Result<(MotifSet, u64)> {
    let source = FileSource::open(path).map_err(|e| bad(path, 0, e))?;
    let len = omicsweave_core::ReadAt::len(&source);
    let limit = usize::try_from(budget.memory).unwrap_or(usize::MAX);
    let defaults = ReaderOptions::default();
    let reader = MotifSet::load_with_options(
        source,
        ReaderOptions {
            max_metadata_bytes: defaults.max_metadata_bytes.min(limit),
            max_block_bytes: defaults.max_block_bytes.min(limit),
            max_stream_bytes: defaults.max_stream_bytes.min(limit),
            max_resident_bytes: limit,
        },
    )
    .map_err(|e| bad(path, 0, e))?;
    Ok((reader, len))
}

/// The motif format a text shows: `None` when it is not motif text. `head`
/// is the start of the file (inflated), up to 64 KiB. Motif text whose
/// layout cannot be told (square 4 x 4 matrices) is an error asking for the
/// format.
pub(crate) fn sniff(head: &[u8], complete: bool) -> Result<Option<Format>> {
    let text = String::from_utf8_lossy(head);
    let Some(first) = text.lines().map(str::trim).find(|l| !l.is_empty()) else {
        return Ok(None);
    };
    if first.starts_with("MEME version") {
        return Ok(Some(Format::Meme));
    }
    if glam2::is_glam2(first) {
        return Ok(Some(Format::Glam2));
    }
    if first.starts_with('{') && text.contains("\"coefficients\"") {
        return Ok(Some(Format::ProBound));
    }
    if first.split_whitespace().next() == Some("Pos")
        && first.split_whitespace().skip(1).all(|w| w.len() == 1)
    {
        return Ok(Some(Format::CisBp));
    }
    if transfac::is_transfac(&text) {
        return Ok(Some(Format::Transfac));
    }
    if first.starts_with('>') {
        return match rows::sniff_headed(&text, complete) {
            rows::Headed::Text => Ok(None),
            rows::Headed::Motif(f) => Ok(Some(f)),
            rows::Headed::Ambiguous => Err(Error::invalid(
                "format",
                "numeric matrices of 4 x 4 cannot be told apart; name the format \
                 (pcm, pwm, jaspar-pfm, homer)",
            )),
        };
    }
    if bamm::is_bamm(&text, complete) {
        return Ok(Some(Format::Bamm));
    }
    Ok(None)
}
