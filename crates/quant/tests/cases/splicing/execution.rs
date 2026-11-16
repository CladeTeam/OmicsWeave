use super::*;

#[test]
fn splicing_outputs_do_not_depend_on_threads_or_spilling() {
    let f = fixture(44);
    let products = || {
        let j = Junctions {
            matrix: Some(JunctionMatrix::single_cell()),
            ..Junctions::default()
        };
        let mut sc = Counts::rna_single_cell();
        sc.velocity = true;
        vec![Product::Junctions(j), Product::Counts(sc)]
    };
    support::assert_execution_invariant(
        f.dir.path(),
        |out| task(&f, out, products(), Some(&f.fasta)),
        16 << 10,
        1,
    );
}
