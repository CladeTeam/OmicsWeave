use crate::support::{
    CHROMS, CountType, MapSpec, Storage, check_all, expected, normalize_pixels, write_and_open,
};
use omicsweave_cool::{PixelBatch, PixelScratch};
use proptest::prelude::*;

fn storage_strategy() -> impl Strategy<Value = Storage> {
    prop_oneof![
        (
            1usize..200,
            any::<bool>(),
            proptest::option::of(1u8..9),
            any::<bool>()
        )
            .prop_map(|(chunk, shuffle, gzip, fletcher)| Storage::Chunked {
                chunk,
                shuffle,
                gzip,
                fletcher,
            }),
        Just(Storage::Contiguous),
        (1usize..100).prop_map(|chunk| Storage::ScaleOffset { chunk }),
    ]
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 48,
        failure_persistence: Some(Box::new(proptest::test_runner::FileFailurePersistence::Direct(
            concat!(env!("CARGO_MANIFEST_DIR"), "/tests/cases/properties.proptest-regressions")
        ))),
        ..ProptestConfig::default()
    })]

    #[test]
    fn random_pixels_round_trip(
        raw in proptest::collection::vec((0u64..175, 0u64..175, 1u32..100_000), 0..1500),
        storage in [storage_strategy(), storage_strategy(), storage_strategy()],
        float in any::<bool>(),
        target in 1u64..600,
        ranges in proptest::collection::vec((0.0f64..1.0, 0.0f64..1.0), 4),
    ) {
        let dir = tempfile::tempdir().unwrap();
        let mut s = MapSpec::new(CHROMS, 100);
        let raw: Vec<(u64, u64, f64)> = raw
            .iter()
            .map(|&(a, b, c)| (a, b, if float { f64::from(c) / 8.0 } else { f64::from(c) }))
            .collect();
        s.pixels = normalize_pixels(&raw, true);
        s.count_type = if float { CountType::F64 } else { CountType::I32 };
        // Scale-offset is an integer filter.
        if float && matches!(storage[2], Storage::ScaleOffset { .. }) {
            s.storage = [storage[0], storage[1], Storage::Contiguous];
        } else {
            s.storage = storage;
        }
        let cool = write_and_open(dir.path(), &s);
        check_all(&cool, &s, target);
        let n = s.pixels.len() as u64;
        let mut scratch = PixelScratch::new();
        let mut out = PixelBatch::new();
        for (a, b) in ranges {
            let (a, b) = ((a * n as f64) as u64, (b * n as f64) as u64);
            let rows = a.min(b)..a.max(b);
            cool.read_pixels(rows.clone(), &mut scratch, &mut out).unwrap();
            prop_assert_eq!(&out, &expected(&s, rows));
        }
    }
}
