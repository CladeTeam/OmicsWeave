use omicsweave_ffm::{
    Alphabet, Factor, FactorKind, FfmReader, Motif, Schema, Support, WriterOptions, table_len,
};
use proptest::prelude::*;

use crate::support::write;

proptest! {
    #[test]
    fn random_factor_models_round_trip(
        specs in prop::collection::vec(
            (1u32..12, prop::option::of(1usize..4), prop::collection::vec((0u8..3, 1usize..4, 0u32..1000), 1..5)),
            1..40,
        ),
        block_values in 1usize..400,
    ) {
        let alphabet = Alphabet::new("six", b"ABCDEF").unwrap();
        let mut motifs = Vec::new();
        for (i, (width, subset, factors)) in specs.into_iter().enumerate() {
            let support = match subset {
                None => Support::All,
                Some(k) => Support::Codes((0..k as u8).map(|c| c * 2).collect()),
            };
            let r = match &support { Support::All => 6, Support::Codes(c) => c.len() + 1 };
            let factors = factors
                .into_iter()
                .map(|(kind, order, salt)| {
                    let order = order.min(width as usize);
                    let positions: Vec<u32> = (0..order as u32).map(|p| (p + salt) % width).collect::<std::collections::BTreeSet<_>>().into_iter().collect();
                    let cells = table_len(r, positions.len()).unwrap();
                    Factor {
                        kind: [FactorKind::Joint, FactorKind::Conditional, FactorKind::Potential][kind as usize],
                        table: (0..cells).map(|c| ((c as u32 + salt) % 5) as f64 + 0.25).collect(),
                        positions,
                    }
                })
                .collect();
            let m = Motif { width, support, sites: (i % 2 == 0).then_some(3.0 + i as f64), factors, gaps: None };
            motifs.push((format!("motif-{i}"), m, vec![]));
        }
        let opts = WriterOptions { block_values, ..WriterOptions::default() };
        let bytes = write(alphabet, Schema::default(), opts, &motifs);
        let reader = FfmReader::open(bytes).unwrap();
        for (i, (id, m, _)) in motifs.iter().enumerate() {
            prop_assert_eq!(reader.find(id), Some(i));
            prop_assert_eq!(&reader.motif(i).unwrap().view().to_motif(), m);
        }
        let resident = reader.into_set().unwrap();
        for (record, (id, m, _)) in resident.iter().zip(&motifs) {
            prop_assert_eq!(record.id(), id);
            prop_assert_eq!(&record.model().to_motif(), m);
            for (i, factor) in m.factors.iter().enumerate() {
                let view = record.model().factor(i).unwrap();
                prop_assert_eq!(view.positions, &factor.positions);
                prop_assert_eq!(view.table, &factor.table);
                prop_assert_eq!(view.kind, factor.kind);
            }
        }
    }
}
