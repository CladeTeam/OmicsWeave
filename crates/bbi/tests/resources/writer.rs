use crate::allocator::counted;
use crate::writer_fixture::{Discard, chroms, options};
use omicsweave_bbi::{BedSchema, BigBedWriter, BigWigWriter};
#[test]
fn disabled_zooms_do_not_retain_the_sampling_runs() {
    let mut wig = BigWigWriter::create(Discard::default(), chroms(), options(false, 0)).unwrap();
    let (_, largest, _) = counted(|| {
        for i in 0..20_000 {
            wig.add(0, i * 8, i * 8 + 4, i as f32).unwrap();
        }
    });
    assert!(largest < 128 * 1024, "largest allocation: {largest}");
    let mut bed = BigBedWriter::create(
        Discard::default(),
        chroms(),
        BedSchema::bed(4),
        options(false, 0),
    )
    .unwrap();
    let (_, largest, _) = counted(|| {
        for i in 0..20_000 {
            bed.add(0, i * 3, i * 3 + 50, b"x").unwrap();
        }
    });
    assert!(largest < 128 * 1024, "largest allocation: {largest}");
}
