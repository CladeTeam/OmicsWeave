use omicsweave_bbi::BbiSidecar;
use omicsweave_core::sidecar::SourceIdentity;
use omicsweave_core::{Error, FileSource};
use std::path::Path;
pub(crate) fn check_sidecar(path: &Path, written: &BbiSidecar) {
    let source = FileSource::open(path).unwrap();
    let built = BbiSidecar::build(&source).unwrap();
    assert_eq!(
        &built, written,
        "writer sidecar differs from one built from the file"
    );
    let identity = SourceIdentity::of(&source).unwrap();
    let bytes = written.to_bytes(identity);
    assert_eq!(&BbiSidecar::parse(&bytes, identity).unwrap(), written);
    let other = SourceIdentity {
        fingerprint: identity.fingerprint ^ 1,
        ..identity
    };
    assert_eq!(BbiSidecar::parse(&bytes, other), Err(Error::StaleSidecar));
    let mut damaged = bytes.clone();
    *damaged.last_mut().unwrap() ^= 1;
    assert!(matches!(
        BbiSidecar::parse(&damaged, identity),
        Err(Error::Invalid { .. })
    ));
}
