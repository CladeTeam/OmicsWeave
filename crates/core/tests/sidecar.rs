use omicsweave_core::Error;
use omicsweave_core::sidecar::{
    SIDECAR_HEADER_LEN, SidecarKind, SourceIdentity, checksum, seal, sidecar_path, unseal,
};

#[test]
fn sealed_bodies_open_only_for_their_source_and_kind() {
    let source: Vec<u8> = (0..300_000u32).map(|i| (i * 7 % 251) as u8).collect();
    let id = SourceIdentity::of(source.as_slice()).unwrap();
    assert_eq!(id, SourceIdentity::of(source.as_slice()).unwrap());
    let body = b"0123456789abcdef".to_vec();
    let sealed = seal(SidecarKind::Bbi, 3, id, &body);
    assert_eq!(sealed.len(), SIDECAR_HEADER_LEN + body.len());
    assert_eq!(
        unseal(&sealed, SidecarKind::Bbi, id).unwrap(),
        (3, body.as_slice())
    );
    assert!(unseal(&sealed, SidecarKind::Hic, id).is_err());

    // A changed sampled byte, or a changed length, is another source; the
    // interior pages sampled start at multiples of len / 15.
    let interior = source.len() / 15 * 7 + 10;
    for edit in [0, interior, source.len() - 1] {
        let mut other = source.clone();
        other[edit] ^= 1;
        let other = SourceIdentity::of(other.as_slice()).unwrap();
        assert_eq!(
            unseal(&sealed, SidecarKind::Bbi, other),
            Err(Error::StaleSidecar)
        );
    }
    // Bytes between the samples are not read.
    let mut unsampled = source.clone();
    unsampled[source.len() / 2] ^= 1;
    assert_eq!(SourceIdentity::of(unsampled.as_slice()).unwrap(), id);
    let longer = SourceIdentity::of(&[source.as_slice(), &[0]].concat()[..]).unwrap();
    assert_eq!(
        unseal(&sealed, SidecarKind::Bbi, longer),
        Err(Error::StaleSidecar)
    );

    let mut damaged = sealed.clone();
    damaged[SIDECAR_HEADER_LEN] ^= 1;
    assert!(matches!(
        unseal(&damaged, SidecarKind::Bbi, id),
        Err(Error::Invalid { .. })
    ));
    assert!(unseal(&sealed[..sealed.len() - 1], SidecarKind::Bbi, id).is_err());
    assert!(matches!(
        unseal(
            b"not a sidecar at all, not at all..............",
            SidecarKind::Bbi,
            id
        ),
        Err(Error::BadMagic { .. })
    ));
    assert_ne!(checksum(b"a"), checksum(b"b"));
}

#[test]
fn small_sources_and_paths() {
    for len in [0usize, 1, 4095, 70_000] {
        let source = vec![1u8; len];
        let id = SourceIdentity::of(source.as_slice()).unwrap();
        assert_eq!(id.len, len as u64);
    }
    assert_eq!(sidecar_path("/a/b.hic").to_str(), Some("/a/b.hic.sidecar"));
}
