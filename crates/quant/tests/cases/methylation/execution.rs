use super::*;

#[test]
fn outputs_do_not_depend_on_threads_or_spilling() {
    let f = fixture(73);
    let m = || product(CPG | CHG, true, CallSource::Reference);
    support::assert_execution_invariant(
        f.dir.path(),
        |out| task(&f, out, m(), true, true),
        8 << 10,
        1,
    );
}
