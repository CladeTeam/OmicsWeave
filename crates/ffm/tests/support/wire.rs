use xxhash_rust::xxh3::xxh3_64;

fn sealed(mut payload: Vec<u8>) -> Vec<u8> {
    payload.extend_from_slice(&xxh3_64(&payload).to_le_bytes());
    payload
}

fn raw(payload: &[u8], elements: usize) -> Vec<u8> {
    let mut descriptor = vec![0; 8];
    descriptor.extend_from_slice(&(elements as u64).to_le_bytes());
    descriptor.extend_from_slice(&(payload.len() as u64).to_le_bytes());
    let mut checked = descriptor.clone();
    checked.extend_from_slice(payload);
    descriptor.extend_from_slice(&xxh3_64(&checked).to_le_bytes());
    descriptor.extend_from_slice(payload);
    descriptor
}

/// A complete one-record file, assembled from the published byte layout.
/// No FFM writer, alphabet serializer or model helper participates in this oracle.
pub fn independent_file() -> Vec<u8> {
    let mut alphabet = vec![2, 1, 1, 0, b'x', b'A', b'B', 1, 0];
    let mut lut = [255; 256];
    lut[b'A' as usize] = 0;
    lut[b'B' as usize] = 1;
    alphabet.extend_from_slice(&lut);
    let alphabet = sealed(alphabet);
    let schema = sealed(vec![
        2, 0, 0, 5, 0, b'l', b'a', b'b', b'e', b'l', 1, 1, 0, b'n',
    ]);
    let mut id_bytes = vec![0, 0, 0, 0, 3, 0, 0, 0];
    id_bytes.extend_from_slice(b"ref");
    let ids = raw(&id_bytes, id_bytes.len());
    let index = raw(&[0; 16], 4);
    let structure: Vec<u8> = [2u32, 2, 1, 2, 0, 1, 0, 1, 2, 0, 1]
        .into_iter()
        .flat_map(u32::to_le_bytes)
        .collect();
    let values: Vec<u8> = [5.0f64, 1.0, 3.0, 2.0, 6.0, 4.0, 4.0]
        .into_iter()
        .flat_map(f64::to_le_bytes)
        .collect();
    let text = "line\nβ\t-7\0".as_bytes();
    let mut block = vec![1, 0, 0, 0, 0, 0, 0, 0];
    block.extend(raw(&structure, 11));
    block.extend(raw(&values, 7));
    block.extend(raw(text, text.len()));

    let directory_at = 128 + alphabet.len() + schema.len() + ids.len();
    let block_at = directory_at + 24 + index.len();
    let mut directory = (block_at as u64).to_le_bytes().to_vec();
    directory.extend_from_slice(&(block.len() as u32).to_le_bytes());
    directory.extend_from_slice(&1u32.to_le_bytes());
    let directory = sealed(directory);
    let sections = [alphabet, schema, ids, directory, index];
    let mut header = vec![0; 128];
    header[..8].copy_from_slice(b"\x89FFM\r\n\x1a\n");
    header[8..10].copy_from_slice(&1u16.to_le_bytes());
    header[16..24].copy_from_slice(&1u64.to_le_bytes());
    header[24..32].copy_from_slice(&1u64.to_le_bytes());
    let mut offset = 128u64;
    for (i, section) in sections.iter().enumerate() {
        let at = 32 + i * 16;
        header[at..at + 8].copy_from_slice(&offset.to_le_bytes());
        header[at + 8..at + 16].copy_from_slice(&(section.len() as u64).to_le_bytes());
        offset += section.len() as u64;
    }
    let checksum = xxh3_64(&header[..120]);
    header[120..].copy_from_slice(&checksum.to_le_bytes());
    for section in sections {
        header.extend(section);
    }
    header.extend(block);
    header
}
