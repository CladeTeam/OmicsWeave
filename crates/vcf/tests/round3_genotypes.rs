use omicsweave_vcf::Genotypes;

#[test]
fn genotype_sample_rejects_overflow_and_preserves_zero_width() {
    let values = [2, 4, 6, 8];
    let gt = Genotypes {
        values: &values,
        width: 2,
    };
    assert_eq!(gt.sample(0), &[2, 4]);
    assert_eq!(gt.sample(1), &[6, 8]);
    for i in [2, 1usize << (usize::BITS - 1), usize::MAX] {
        assert!(std::panic::catch_unwind(|| gt.sample(i)).is_err());
    }
    let addition_overflow = Genotypes {
        values: &values,
        width: 1,
    };
    assert!(std::panic::catch_unwind(|| addition_overflow.sample(usize::MAX)).is_err());

    let zero = Genotypes {
        values: &[],
        width: 0,
    };
    for i in [0, 1usize << (usize::BITS - 1), usize::MAX] {
        assert!(zero.sample(i).is_empty());
    }
}
