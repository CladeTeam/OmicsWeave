use omicsweave_ffm::{
    Alphabet, Codec, Error, FfmReader, FfmWriter, Key, KeyType, MAX_MOTIFS, Motif, MotifSet,
    ReaderOptions, Schema, WriterOptions,
};
use xxhash_rust::xxh3::xxh3_64;

use crate::support::{pfm, pwm, write};

#[test]
fn independent_file_rejects_every_truncated_prefix_and_trailing_bytes() {
    let bytes = crate::support::wire::independent_file();
    for end in 0..bytes.len() {
        let result = FfmReader::open(&bytes[..end]).and_then(|reader| reader.into_set());
        assert!(result.is_err(), "accepted truncated file ending at {end}");
    }
    let mut trailing = bytes;
    trailing.push(0);
    assert!(FfmReader::open(trailing).is_err());
}

#[test]
fn resealed_header_rejects_unsupported_version_flags_and_reserved_words() {
    let original = crate::support::wire::independent_file();
    for (at, value) in [(8, 2), (12, 1), (112, 1)] {
        let mut bytes = original.clone();
        bytes[at] = value;
        let checksum = xxh3_64(&bytes[..120]);
        bytes[120..128].copy_from_slice(&checksum.to_le_bytes());
        assert!(
            FfmReader::open(bytes).is_err(),
            "header byte {at} = {value}"
        );
    }
}

fn u64_at(bytes: &[u8], at: usize) -> u64 {
    u64::from_le_bytes(bytes[at..at + 8].try_into().unwrap())
}

fn fixture(annotation: bool) -> Vec<u8> {
    let schema = if annotation {
        Schema::new(vec![Key::new("name", KeyType::String)]).unwrap()
    } else {
        Schema::default()
    };
    let opts = WriterOptions {
        annotation_codec: Codec::Raw,
        ..WriterOptions::default()
    };
    let mut writer = FfmWriter::new(Vec::new(), Alphabet::dna(), schema, opts).unwrap();
    writer
        .push(
            "a",
            &Motif::pfm(&[vec![1.0; 4]], None),
            if annotation { &["abc"] } else { &[] },
        )
        .unwrap();
    writer.finish().unwrap()
}

fn reseal_stream(bytes: &mut [u8], at: usize) {
    let len = u64_at(bytes, at + 16) as usize;
    let mut checked = bytes[at..at + 24].to_vec();
    checked.extend_from_slice(&bytes[at + 32..at + 32 + len]);
    bytes[at + 24..at + 32].copy_from_slice(&xxh3_64(&checked).to_le_bytes());
}

fn block_streams(bytes: &[u8]) -> [usize; 3] {
    let directory = u64_at(bytes, 80) as usize;
    let structure = u64_at(bytes, directory) as usize + 8;
    let values = structure + 32 + u64_at(bytes, structure + 16) as usize;
    let annotation = values + 32 + u64_at(bytes, values + 16) as usize;
    [structure, values, annotation]
}

#[test]
fn index_count_is_checked_before_allocating_the_claimed_output() {
    let mut bytes = fixture(false);
    let index = u64_at(&bytes, 96) as usize;
    bytes[index + 8..index + 16].copy_from_slice(&(u32::MAX as u64).to_le_bytes());
    reseal_stream(&mut bytes, index);
    let error = match FfmReader::open(bytes) {
        Ok(_) => panic!("accepted invalid index"),
        Err(e) => e,
    };
    assert!(
        matches!(error, Error::Invalid { what: "FFM index", ref detail } if detail.contains("length"))
    );
}

#[test]
fn huge_numeric_length_is_bounded_even_with_valid_checksum() {
    let mut bytes = fixture(false);
    let values = block_streams(&bytes)[1];
    bytes[values + 8..values + 16].copy_from_slice(&(u32::MAX as u64).to_le_bytes());
    reseal_stream(&mut bytes, values);
    let reader = FfmReader::open(bytes).unwrap();
    let error = match reader.motif(0) {
        Ok(_) => panic!("accepted invalid values"),
        Err(e) => e,
    };
    assert!(
        matches!(error, Error::Invalid { ref detail, .. } if detail.contains("allocation limit"))
    );
}

#[test]
fn motif_count_must_fit_the_index_stream() {
    let mut bytes = fixture(false);
    bytes[16..24].copy_from_slice(&(MAX_MOTIFS as u64 + 1).to_le_bytes());
    let checksum = xxh3_64(&bytes[..120]);
    bytes[120..128].copy_from_slice(&checksum.to_le_bytes());
    assert!(FfmReader::open(bytes).is_err());
}

#[test]
fn resident_load_validates_annotation_even_when_file_reader_is_lazy() {
    let mut bytes = fixture(true);
    let annotation = block_streams(&bytes)[2];
    bytes[annotation + 32] = 255;
    reseal_stream(&mut bytes, annotation);
    let reader = FfmReader::open(&bytes).unwrap();
    assert!(reader.motif(0).is_ok());
    assert!(reader.annotation(0).is_err());
    assert!(reader.into_set().is_err());
    assert!(MotifSet::load(&bytes).is_err());
    let options = ReaderOptions {
        max_stream_bytes: 2,
        ..ReaderOptions::default()
    };
    assert!(FfmReader::open_with_options(&bytes, options).is_err());
}

#[test]
fn checked_block_shapes_fail_before_constructing_factor_views() {
    let opts = WriterOptions {
        structure_codec: Codec::Raw,
        values_codec: Codec::Raw,
        ..WriterOptions::default()
    };
    let mut writer = FfmWriter::new(Vec::new(), Alphabet::dna(), Schema::default(), opts).unwrap();
    writer
        .push(
            "model",
            &Motif::pfm(&[vec![1.0; 4], vec![2.0; 4]], None),
            &[],
        )
        .unwrap();
    let original = writer.finish().unwrap();
    let structure = block_streams(&original)[0];
    // width, support count, flags, factor count, then kind/order/position.
    for (word, value) in [
        (0, 0u32),
        (1, 0),
        (2, 8),
        (3, 0),
        (3, u32::MAX),
        (4, 3),
        (5, 0),
        (5, 17),
        (6, 2),
    ] {
        let mut bytes = original.clone();
        let at = structure + 32 + word * 4;
        bytes[at..at + 4].copy_from_slice(&value.to_le_bytes());
        reseal_stream(&mut bytes, structure);
        let reader = FfmReader::open(&bytes).unwrap();
        assert!(reader.motif(0).is_err(), "word {word}, value {value}");
        assert!(reader.motifs().next().unwrap().is_err());
        assert!(reader.into_set().is_err());
    }
}

#[test]
fn every_flipped_byte_is_detected() {
    let motifs: Vec<_> = (0..4)
        .map(|i| (format!("id{i}"), pwm(5, i), vec![]))
        .collect();
    let bytes = write(
        Alphabet::dna(),
        Schema::default(),
        WriterOptions::default(),
        &motifs,
    );
    for at in 0..bytes.len() {
        let mut bad = bytes.clone();
        bad[at] ^= 0x40;
        let read_all = FfmReader::open(bad)
            .and_then(|r| (0..r.len()).try_for_each(|i| r.motif(i).map(|_| ())));
        assert!(read_all.is_err(), "corruption at byte {at} went unnoticed");
    }
}

#[test]
fn reader_limits_apply_before_decode_and_resident_growth() {
    let mut w = FfmWriter::new(
        Vec::new(),
        Alphabet::dna(),
        Schema::default(),
        WriterOptions::default(),
    )
    .unwrap();
    w.push("a", &pfm(), &[]).unwrap();
    let bytes = w.finish().unwrap();
    let metadata = ReaderOptions {
        max_metadata_bytes: 1,
        ..ReaderOptions::default()
    };
    assert!(FfmReader::open_with_options(&bytes, metadata).is_err());
    let blocks = ReaderOptions {
        max_block_bytes: 1,
        ..ReaderOptions::default()
    };
    let reader = FfmReader::open_with_options(&bytes, blocks).unwrap();
    assert!(reader.motif(0).is_err());
    let streams = ReaderOptions {
        max_stream_bytes: 1,
        ..ReaderOptions::default()
    };
    assert!(FfmReader::open_with_options(&bytes, streams).is_err());
    let resident = ReaderOptions {
        max_resident_bytes: 1,
        ..ReaderOptions::default()
    };
    assert!(MotifSet::load_with_options(&bytes, resident).is_err());
}
