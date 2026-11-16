use std::sync::Arc;

use omicsweave_core::cache::BlockCache;
use omicsweave_core::{Error, Region, Strand};
use omicsweave_msa::MsaReader;

use crate::support::*;

#[test]
fn cached_projection_and_clones_keep_staging_independent() {
    let records = [(2, 5, "s ref.c 2 3 + 100 AcG;s dog.c 1 2 + 100 C-N;".into())];
    let file = tempfile::NamedTempFile::new().unwrap();
    std::fs::write(file.path(), bigmaf(&records, Some(SQL), 4)).unwrap();
    let mut original = MsaReader::open(
        omicsweave_core::FileSource::open(file.path()).unwrap(),
        None,
        b'.',
    )
    .unwrap();
    let cache = Arc::new(BlockCache::new(1 << 20, 64));
    original.set_cache(Some(cache.clone())).unwrap();
    let shape = original.query(Region::new("c", 0, 7)).unwrap();
    let mut cloned = original.try_clone().unwrap();
    assert!(cloned.species().is_empty());
    cloned.query(Region::new("c", 90, 100)).unwrap();
    let mut want = vec![0; shape.elements().unwrap()];
    original.write(Strand::Forward, &mut want).unwrap();
    assert_eq!(&want, b"..AcG....C-N..");
    let loads = cache.stats().loads;
    assert!(loads > 0);
    for _ in 0..2 {
        let shape = original.query(Region::new("c", 0, 7)).unwrap();
        let mut got = vec![0; shape.elements().unwrap()];
        original.write(Strand::Forward, &mut got).unwrap();
        assert_eq!(got, want);
    }
    assert_eq!(cache.stats().loads, loads);
    assert!(cache.stats().hits > 0);
    original.set_cache(None).unwrap();
    assert_eq!(fetch(&mut original, 2, 5, Strand::Forward), b"AcGC-N");
    assert_eq!(cache.stats().loads, loads);
    let shape = cloned.query(Region::new("c", 0, 7)).unwrap();
    let mut got = vec![0; shape.elements().unwrap()];
    cloned.write(Strand::Forward, &mut got).unwrap();
    assert_eq!(got, want);
    assert_eq!(cache.stats().loads, loads);
    drop(original);
    cache.clear();
    assert_eq!(fetch(&mut cloned, 2, 5, Strand::Forward), b"AcGC-N");
    assert!(cache.stats().loads > loads);
}

#[test]
fn bounds_empty_queries_and_staging_errors() {
    let records = [(98, 100, "s ref.c 98 2 + 100 AC".into())];
    let mut r = reader(&records, Some(vec!["ref".into()]), b'.');
    assert!(matches!(
        r.write(Strand::Forward, &mut []),
        Err(Error::NotQueried)
    ));
    assert_eq!(fetch(&mut r, 98, 104, Strand::Forward), b"AC....");
    assert_eq!(fetch(&mut r, 100, 104, Strand::Forward), b"....");
    assert!(fetch(&mut r, 10, 10, Strand::Reverse).is_empty());
    r.query(Region::new("c", 0, 2)).unwrap();
    assert!(matches!(
        r.write(Strand::Forward, &mut []),
        Err(Error::OutputSize { .. })
    ));
    assert_eq!(r.query(Region::new("unknown", 0, 2)).unwrap().len, 2);
    let mut unknown = [0; 2];
    r.write(Strand::Forward, &mut unknown).unwrap();
    assert_eq!(&unknown, b"..");
    assert!(matches!(
        r.query(Region::new("c", 5, 2)),
        Err(Error::InvalidRegion { .. })
    ));
    assert!(matches!(
        r.write(Strand::Forward, &mut []),
        Err(Error::NotQueried)
    ));
    assert_eq!(fetch(&mut r, 98, 100, Strand::Forward), b"AC");
    let mut dynamic = reader(&records, None, b'-');
    assert!(fetch(&mut dynamic, 0, 1, Strand::Forward).is_empty());
    assert!(dynamic.species().is_empty());
    assert_eq!(fetch(&mut dynamic, 98, 100, Strand::Forward), b"AC");
    assert!(fetch(&mut dynamic, 0, 1, Strand::Forward).is_empty());
    assert!(dynamic.species().is_empty());
}

#[test]
fn exact_output_size_is_checked_before_writing_and_staging_is_reusable() {
    let mut r = reader(
        &[(2, 5, "s ref.c 2 3 + 100 AcG;".into())],
        Some(vec!["ref".into()]),
        b'.',
    );
    r.query(Region::new("c", 2, 5)).unwrap();
    for len in [0, 2, 4] {
        let mut output = vec![0xa5; len];
        assert!(
            matches!(r.write(Strand::Reverse, &mut output), Err(Error::OutputSize { have, need: 3, .. }) if have == len)
        );
        assert_eq!(
            output,
            vec![0xa5; len],
            "wrong-sized destination must be untouched"
        );
    }
    let mut guarded = [0xa5; 5];
    r.write(Strand::Reverse, &mut guarded[1..4]).unwrap();
    assert_eq!(guarded, [0xa5, b'C', b'g', b'T', 0xa5]);
    r.write(Strand::Forward, &mut guarded[1..4]).unwrap();
    assert_eq!(guarded, [0xa5, b'A', b'c', b'G', 0xa5]);
}

#[test]
fn failed_decode_invalidates_partial_output_and_a_later_query_recovers() {
    let mut r = reader(
        &[
            (0, 3, "s ref.c 0 3 + 100 ACG;s dog.x 0 3 + 10 TTT;".into()),
            (10, 13, "s ref.c 10 3 + 100 AAC;s bad.x 0 3 + 10 A-;".into()),
            (20, 23, "s ref.c 20 3 + 100 GAT;".into()),
        ],
        None,
        b'.',
    );
    assert_eq!(fetch(&mut r, 0, 3, Strand::Forward), b"ACGTTT");
    assert!(r.query(Region::new("c", 0, 13)).is_err());
    let mut out = [0xa5; 26];
    assert_eq!(r.write(Strand::Forward, &mut out), Err(Error::NotQueried));
    assert_eq!(out, [0xa5; 26]);
    assert_eq!(fetch(&mut r, 20, 23, Strand::Forward), b"GAT");
    assert_eq!(r.species(), ["ref"]);
}

#[test]
fn source_read_failure_leaves_no_output_and_can_be_retried() {
    use omicsweave_core::{ReadAt, Result};
    use std::sync::atomic::{AtomicBool, Ordering};

    struct FallibleSource {
        bytes: Vec<u8>,
        fail: Arc<AtomicBool>,
    }
    impl ReadAt for FallibleSource {
        fn len(&self) -> u64 {
            self.bytes.len() as u64
        }
        fn read_exact_at(&self, out: &mut [u8], at: u64) -> Result<()> {
            if self.fail.load(Ordering::Relaxed) {
                return Err(std::io::Error::other("injected source failure").into());
            }
            self.bytes.read_exact_at(out, at)
        }
    }
    let fail = Arc::new(AtomicBool::new(false));
    let source = FallibleSource {
        bytes: bigmaf(&[(2, 5, "s ref.c 2 3 + 100 ACG;".into())], Some(SQL), 4),
        fail: fail.clone(),
    };
    let mut r = MsaReader::open(source, None, b'.').unwrap();
    fail.store(true, Ordering::Relaxed);
    assert!(matches!(
        r.query(Region::new("c", 0, 7)),
        Err(Error::Io { .. })
    ));
    let mut out = [0xa5; 7];
    assert_eq!(r.write(Strand::Forward, &mut out), Err(Error::NotQueried));
    assert_eq!(out, [0xa5; 7]);
    fail.store(false, Ordering::Relaxed);
    assert_eq!(fetch(&mut r, 0, 7, Strand::Forward), b"..ACG..");
}

#[test]
fn fixed_clone_retains_configuration_but_needs_its_own_query() {
    let file = tempfile::NamedTempFile::new().unwrap();
    std::fs::write(
        file.path(),
        bigmaf(&[(2, 5, "s ref.c 2 3 + 100 AcG;".into())], Some(SQL), 4),
    )
    .unwrap();
    let mut original = MsaReader::open(
        omicsweave_core::FileSource::open(file.path()).unwrap(),
        Some(vec!["absent".into(), "ref".into()]),
        b'A',
    )
    .unwrap();
    assert_eq!(fetch(&mut original, 2, 5, Strand::Forward), b"AAAAcG");
    let mut cloned = original.try_clone().unwrap();
    assert_eq!(cloned.species(), ["absent", "ref"]);
    let mut out = [0xa5; 6];
    assert_eq!(
        cloned.write(Strand::Forward, &mut out),
        Err(Error::NotQueried)
    );
    assert_eq!(out, [0xa5; 6]);
    drop(original);
    assert_eq!(fetch(&mut cloned, 1, 6, Strand::Reverse), b"AAAAAACgTA");
}
