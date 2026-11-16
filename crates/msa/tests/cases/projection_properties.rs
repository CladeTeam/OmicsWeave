use omicsweave_core::Strand;
use proptest::prelude::*;

use crate::support::*;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(80))]
    #[test]
    fn generated_reference_positions_match_projected_truth(
        blocks in prop::collection::vec(
            (0u32..85, prop::collection::vec((0usize..8, prop::option::of(0usize..8), 0usize..3), 1..12), any::<bool>()),
            1..20,
        ),
        start in 0u32..100,
        length in 0u32..30,
    ) {
        let alphabet = b"ACGTacgt";
        let mut blocks = blocks;
        blocks.sort_by_key(|block| block.0);
        let mut records = Vec::new();
        let mut truth = vec![vec![b'.'; 100]; 2];
        for (at, positions, include_dog) in blocks {
            let size = positions.len() as u32;
            let mut reference = String::new();
            let mut dog = String::new();
            // Truth is expressed in reference positions before alignment
            // insertions are encoded, independent of the reader's column map.
            for (offset, (a, b, insertions)) in positions.into_iter().enumerate() {
                let ref_base = alphabet[a];
                let dog_base = b.map_or(b'-', |index| alphabet[index]);
                truth[0][at as usize + offset] = ref_base;
                if include_dog {
                    truth[1][at as usize + offset] = dog_base;
                }
                reference.extend(std::iter::repeat_n('-', insertions));
                dog.extend(std::iter::repeat_n('N', insertions));
                reference.push(ref_base as char);
                dog.push(dog_base as char);
            }
            reference.push('-');
            dog.push('N');
            let mut block = format!("s ref.c {at} {size} + 100 {reference};");
            if include_dog {
                let bases = dog.bytes().filter(|&base| base != b'-').count();
                block.push_str(&format!("s dog.c 0 {bases} - 100 {dog};"));
            }
            records.push((at, at + size, block));
        }
        let mut expected = Vec::new();
        for row in &truth {
            expected.extend((start..start + length).map(|position| row.get(position as usize).copied().unwrap_or(b'.')));
        }
        let mut r = reader(&records, Some(vec!["ref".into(), "dog".into()]), b'.');
        prop_assert_eq!(fetch(&mut r, start, start + length, Strand::Forward), expected.clone());
        if length > 0 {
            for row in expected.chunks_mut(length as usize) {
                row.reverse();
                for base in row {
                    *base = match *base {
                        b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C',
                        b'a' => b't', b't' => b'a', b'c' => b'g', b'g' => b'c',
                        base => base,
                    };
                }
            }
        }
        prop_assert_eq!(fetch(&mut r, start, start + length, Strand::Reverse), expected);
    }
}
