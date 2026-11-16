pub(crate) use crate::rng::Rng;
use omicsweave_core::Strings;
use omicsweave_ffm::{
    Alphabet, Factor, FactorKind, FfmReader, FfmWriter, Gaps, Motif, MotifRef, Schema, Support,
    WriterOptions,
};
use omicsweave_motif::profile::{
    Glam2Text, ProfileProgram, ProfileScoring, TopAlignments, glam2_background, write_glam2_text,
};
use omicsweave_motif::{ScanOptions, Strands};
use std::fmt::Write as _;

/// A GLAM2 motif as `glam2scan` reads it: the first line whose first word
/// starts and ends with `*` gives the columns (`*`) and insertion slots;
/// each following line with at least three words gives an aligned row as
/// its third word.
pub(crate) fn parse_glam2_motif(text: &str, alphabet: &Alphabet) -> Motif {
    let mut lines = text.lines();
    let key = loop {
        let line = lines.next().expect("a column line");
        if let Some(w) = line.split_whitespace().next()
            && w.starts_with('*')
            && w.ends_with('*')
        {
            break w.as_bytes().to_vec();
        }
    };
    let rows: Vec<Vec<u8>> = lines
        .map_while(|l| l.split_whitespace().nth(2).map(|w| w.as_bytes().to_vec()))
        .collect();
    let width = key.iter().filter(|&&c| c == b'*').count();
    let n = alphabet.len();
    let mut residues = vec![vec![0.0; n]; width];
    let mut deletions = vec![0.0; width];
    let mut insertions = vec![0.0; width - 1];
    for row in &rows {
        assert_eq!(row.len(), key.len());
        let mut k = 0;
        for (&slot, &c) in key.iter().zip(row) {
            if slot == b'*' {
                if c == b'.' {
                    deletions[k] += 1.0;
                } else {
                    residues[k][alphabet.code(c).expect("residue") as usize] += 1.0;
                }
                k += 1;
            } else if c != b'.' {
                insertions[k - 1] += 1.0;
            }
        }
    }
    Motif {
        width: width as u32,
        support: Support::All,
        sites: Some(rows.len() as f64),
        factors: residues
            .into_iter()
            .enumerate()
            .map(|(p, table)| Factor {
                kind: FactorKind::Joint,
                positions: vec![p as u32],
                table,
            })
            .collect(),
        gaps: Some(Gaps {
            deletions,
            insertions,
        }),
    }
}

/// Sequences parsed as `glam2scan` reads FASTA: records start at `>`, the
/// name is the title's first word, whitespace inside sequences is dropped.
pub(crate) struct Fasta {
    names: Vec<u8>,
    name_offsets: Vec<i64>,
    pub(crate) data: Vec<u8>,
    offsets: Vec<i64>,
}

impl Fasta {
    pub(crate) fn parse(text: &str) -> Fasta {
        let mut f = Fasta {
            names: Vec::new(),
            name_offsets: vec![0],
            data: Vec::new(),
            offsets: vec![0],
        };
        for record in text.split('>').skip(1) {
            let (title, body) = record.split_once('\n').unwrap_or((record, ""));
            f.names
                .extend_from_slice(title.split_whitespace().next().unwrap_or("").as_bytes());
            f.name_offsets.push(f.names.len() as i64);
            f.data
                .extend(body.bytes().filter(|b| !b.is_ascii_whitespace()));
            f.offsets.push(f.data.len() as i64);
        }
        f
    }
    pub(crate) fn names(&self) -> Strings<'_> {
        Strings {
            offsets: &self.name_offsets,
            data: &self.names,
        }
    }
    pub(crate) fn sequences(&self) -> Strings<'_> {
        Strings {
            offsets: &self.offsets,
            data: &self.data,
        }
    }
}

pub(crate) fn program(
    motifs: &[Motif],
    alphabet: &Alphabet,
    scoring: &ProfileScoring,
    strands: Strands,
) -> omicsweave_core::Result<ProfileProgram> {
    let mut w = FfmWriter::new(
        Vec::new(),
        alphabet.clone(),
        Schema::default(),
        WriterOptions::default(),
    )?;
    for (i, m) in motifs.iter().enumerate() {
        w.push(&format!("P{i}"), m, &[])?;
    }
    let reader = FfmReader::open(w.finish()?)?;
    let refs: Vec<MotifRef> = (0..reader.len())
        .map(|i| reader.motif(i))
        .collect::<omicsweave_core::Result<_>>()?;
    let views: Vec<_> = refs.iter().map(|r| r.view()).collect();
    ProfileProgram::compile(
        &views,
        alphabet,
        &glam2_background(alphabet)?,
        scoring,
        strands,
    )
}

/// Our hits for motif 0 in `glam2scan -t`'s layout (without the header).
pub(crate) fn hits_text(
    program: &ProfileProgram,
    fasta: &Fasta,
    n: usize,
    threads: Option<usize>,
) -> String {
    let top = program
        .scan(
            fasta.sequences(),
            TopAlignments::glam2(n, fasta.names()),
            &ScanOptions {
                threads,
                ..ScanOptions::default()
            },
        )
        .unwrap();
    let mut out = Vec::new();
    let text = Glam2Text::new(program.alphabet(), fasta.names(), fasta.sequences());
    write_glam2_text(&mut out, &top[0], &text).unwrap();
    String::from_utf8(out).unwrap()
}

pub(crate) const HAND_MOTIF: &str = "some header line
Score: 1  Columns: 6  Sequences: 8

                ***..***
r1            1 acg..tca 6 + 1
r2            1 acgg.tca 7 + 1
r3            1 acgggtca 8 + 1
r4            1 a.g..tca 5 + 1
r5            1 acg..t.a 5 + 1
r6            1 tcg..tca 6 + 1
r7            1 acgt.tca 7 + 1
r8            1 acgN.tca 7 + 1

                acg  tca
";

pub(crate) const HAND_FASTA: &str = "junk before
>q1 description here
ACGTCA
>q2
ttTGACGTnnACGGTCAgg
>q3
acguca-acgtca
>q1
tgacgt
>empty
>q5
aNNNNa
";

/// `glam2scan -t -2 -n 7 n` on the two texts above.
pub(crate) const HAND_HITS: &str = "                           ******
q1                       1 acgtca 6 + 9.57

                           ******
q1                       6 acgtca 1 - 9.57

                           ******
q2                       8 acgtca 3 - 9.57

                           ******
q3                       1 acgtca 6 + 9.57

                           ******
q3                       8 acgtca 13 + 9.57

                            ***.***
q2                       11 acggtca 17 + 8.09

                            ***..***
q3                       11 acgtNtga 4 - 2.40

";

pub(crate) fn random_motif(
    rng: &mut Rng,
    residues: &[u8],
    max_width: usize,
    max_rows: usize,
) -> String {
    let width = 1 + rng.below(max_width);
    let rows = 1 + rng.below(max_rows);
    let mut key = String::new();
    for i in 0..width {
        key.push('*');
        if i + 1 < width {
            for _ in 0..[0, 0, 1, 2][rng.below(4)] {
                key.push('.');
            }
        }
    }
    let mut column: Vec<u8> = residues.to_vec();
    column.push(b'.');
    let mut insert: Vec<u8> = residues.to_vec();
    insert.extend_from_slice(b"N...");
    let mut table: Vec<Vec<u8>> = (0..rows)
        .map(|_| {
            key.bytes()
                .map(|c| rng.pick(if c == b'*' { &column } else { &insert }))
                .collect()
        })
        .collect();
    // FFM stores a column's residue counts as a joint factor, which needs
    // some mass: every column keeps at least one residue.
    for (k, c) in key.bytes().enumerate() {
        if c == b'*' && table.iter().all(|row| row[k] == b'.') {
            table[0][k] = residues[0];
        }
    }
    let mut text = format!("header\n\n    {key}\n");
    for (r, row) in table.iter().enumerate() {
        writeln!(text, "n{r} 1 {} 9 + 1", std::str::from_utf8(row).unwrap()).unwrap();
    }
    text.push('\n');
    text
}

pub(crate) fn random_fasta(rng: &mut Rng, residues: &[u8], count: usize, max_len: usize) -> String {
    let mut text = String::new();
    for k in 0..count {
        let len = if k == 2 { 0 } else { rng.below(max_len + 1) };
        if k % 7 == 3 {
            writeln!(text, ">a_sequence_name_longer_than_24_columns_{k} x").unwrap();
        } else {
            writeln!(text, ">s{} some description", k % 9).unwrap();
        }
        let seq: Vec<u8> = (0..len).map(|_| rng.pick(residues)).collect();
        for chunk in seq.chunks(60) {
            text.push_str(std::str::from_utf8(chunk).unwrap());
            text.push('\n');
        }
    }
    text
}
