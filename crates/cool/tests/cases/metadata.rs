use crate::support::{
    CHROMS, MapSpec, check_all, normalize_pixels, sample_pixels, spec, write_and_open, write_cool,
    write_mcool,
};
use omicsweave_cool::{Cool, StorageMode, is_mcool, mcool_resolutions};
use omicsweave_core::Error;

#[test]
fn metadata() {
    let dir = tempfile::tempdir().unwrap();
    let s = spec(500, 1);
    let cool = write_and_open(dir.path(), &s);
    let chroms: Vec<(&str, u32)> = cool
        .chroms()
        .iter()
        .map(|c| (c.name.as_str(), c.length))
        .collect();
    assert_eq!(chroms, CHROMS);
    assert_eq!(cool.bin_size(), 100);
    assert_eq!(cool.n_bins(), 100 + 74 + 1);
    assert_eq!(cool.chrom_offsets(), &[0, 100, 174, 175]);
    assert_eq!(cool.n_pixels(), s.pixels.len() as u64);
    assert_eq!(cool.storage_mode(), StorageMode::SymmetricUpper);
    assert_eq!(cool.assembly(), Some("test"));
    assert_eq!(cool.locate(0), (0, 0));
    assert_eq!(cool.locate(99), (0, 99));
    assert_eq!(cool.locate(100), (1, 0));
    assert_eq!(cool.locate(173), (1, 73));
    assert_eq!(cool.locate(174), (2, 0));
}

#[test]
fn locate_skips_empty_chromosomes() {
    let dir = tempfile::tempdir().unwrap();
    let mut s = MapSpec::new(&[("a", 250), ("empty", 0), ("b", 100)], 100);
    s.pixels = sample_pixels(s.n_bins(), 20, 3, true);
    let cool = write_and_open(dir.path(), &s);
    assert_eq!(cool.chrom_offsets(), &[0, 3, 3, 4]);
    assert_eq!(cool.locate(2), (0, 2));
    assert_eq!(cool.locate(3), (2, 0));
}

#[test]
#[should_panic(expected = "out of range")]
fn locate_rejects_the_bin_sentinel() {
    let dir = tempfile::tempdir().unwrap();
    let cool = write_and_open(dir.path(), &spec(10, 81));
    cool.locate(cool.n_bins());
}

#[test]
fn plain_integer_bins_and_variable_length_names() {
    let dir = tempfile::tempdir().unwrap();
    let mut s = spec(500, 9);
    s.chrom_enum = false;
    s.fixed_names = false;
    s.format_version = 2;
    let cool = write_and_open(dir.path(), &s);
    assert_eq!(cool.chroms()[1].name, "chr2");
    check_all(&cool, &s, 100);
}

#[test]
fn square_storage_mode() {
    let dir = tempfile::tempdir().unwrap();
    let mut s = MapSpec::new(CHROMS, 100);
    let mut raw = sample_pixels(s.n_bins(), 800, 10, true);
    // Lower-triangle pixels exist in square maps.
    raw.extend(raw.clone().iter().map(|&(a, b, c)| (b, a, c + 1.0)));
    s.pixels = normalize_pixels(&raw, false);
    s.storage_mode = Some("square");
    let cool = write_and_open(dir.path(), &s);
    assert_eq!(cool.storage_mode(), StorageMode::Square);
    assert!(s.pixels.iter().any(|p| p.0 > p.1));
    check_all(&cool, &s, 100);
}

#[test]
fn mcool_resolution_groups() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("m.mcool");
    let mut fine = MapSpec::new(CHROMS, 100);
    fine.pixels = sample_pixels(fine.n_bins(), 1500, 11, true);
    let mut coarse = MapSpec::new(CHROMS, 1000);
    coarse.pixels = sample_pixels(coarse.n_bins(), 100, 12, true);
    // Written coarse first: resolutions come back sorted regardless.
    write_mcool(&path, &[coarse.clone(), fine.clone()]);

    assert!(is_mcool(&path).unwrap());
    assert_eq!(mcool_resolutions(&path).unwrap(), vec![100, 1000]);
    for s in [&fine, &coarse] {
        let cool = Cool::open(&path, &format!("/resolutions/{}", s.bin_size)).unwrap();
        assert_eq!(cool.bin_size(), s.bin_size);
        assert_eq!(cool.n_bins(), s.n_bins());
        check_all(&cool, s, 128);
    }
    // The root of an .mcool is not a map.
    assert!(matches!(Cool::open(&path, "/"), Err(Error::Invalid { .. })));
    assert!(matches!(
        Cool::open(&path, "/resolutions/5"),
        Err(Error::Invalid { .. })
    ));

    let single = dir.path().join("s.cool");
    write_cool(&single, &fine);
    assert!(!is_mcool(&single).unwrap());
    assert!(mcool_resolutions(&single).is_err());
}
