use omicsweave_core::{Strand, Strings};
use omicsweave_ffm::Alphabet;
use omicsweave_motif::HitEntry;
use omicsweave_motif::fimo::{FimoText, write_text};

#[test]
fn text_preserves_forward_bytes_and_complements_reverse_coordinates() {
    let alphabet = Alphabet::dna();
    let text = FimoText {
        alphabet: &alphabet,
        motif_ids: &["m"],
        motif_alt_ids: &[],
        sequence_names: Strings {
            offsets: &[0, 3],
            data: b"seq",
        },
        sequences: Strings {
            offsets: &[0, 5],
            data: b"aCgTN",
        },
        origins: Some(&[100]),
    };
    let forward = HitEntry {
        score: 2.5,
        pvalue: 0.001,
        sequence: 0,
        motif: 0,
        start: 0,
        end: 4,
        strand: Strand::Forward,
    };
    let reverse = HitEntry {
        start: 1,
        end: 5,
        strand: Strand::Reverse,
        ..forward
    };
    let mut out = Vec::new();
    write_text(&mut out, &[forward, reverse], &text).unwrap();
    assert_eq!(
        String::from_utf8(out).unwrap(),
        concat!(
            "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n",
            "m\t\tseq\t101\t104\t+\t2.5\t0.001\t\taCgT\n",
            "m\t\tseq\t102\t105\t-\t2.5\t0.001\t\tNACG\n",
        )
    );
    let mut empty = Vec::new();
    write_text(&mut empty, &[], &text).unwrap();
    assert!(empty.is_empty(), "no hits must not emit a header");
}

#[test]
fn text_propagates_output_failure() {
    use std::io::{self, Write};
    struct Full;
    impl Write for Full {
        fn write(&mut self, _: &[u8]) -> io::Result<usize> {
            Err(io::Error::from(io::ErrorKind::StorageFull))
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }
    let alphabet = Alphabet::dna();
    let text = FimoText {
        alphabet: &alphabet,
        motif_ids: &["m"],
        motif_alt_ids: &[],
        sequence_names: Strings {
            offsets: &[0, 1],
            data: b"s",
        },
        sequences: Strings {
            offsets: &[0, 1],
            data: b"A",
        },
        origins: None,
    };
    let hit = HitEntry {
        score: 1.0,
        pvalue: 0.25,
        sequence: 0,
        motif: 0,
        start: 0,
        end: 1,
        strand: Strand::Forward,
    };
    assert_eq!(
        write_text(Full, &[hit], &text).unwrap_err().kind(),
        io::ErrorKind::StorageFull
    );
}
