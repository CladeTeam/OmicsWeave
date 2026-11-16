use crate::support::{MapSpec, spec, write_cool};
use omicsweave_cool::Cool;
use omicsweave_core::Error;

#[test]
fn datasets_require_one_dimensional_columns() {
    for name in [
        "chroms/name",
        "chroms/length",
        "bins/chrom",
        "bins/start",
        "bins/end",
        "indexes/chrom_offset",
        "pixels/bin1_id",
        "pixels/bin2_id",
        "pixels/count",
    ] {
        for scalar in [false, true] {
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("invalid-rank.cool");
            write_cool(&path, &spec(10, 71));
            {
                let file = hdf5_metno::File::open_rw(&path).unwrap();
                file.unlink(name).unwrap();
                let builder = file.new_dataset::<i64>();
                if scalar {
                    builder.shape(()).create(name).unwrap();
                } else {
                    builder.shape((1, 1)).create(name).unwrap();
                }
            }
            let err = Cool::open(&path, "/").err().expect("invalid rank accepted");
            assert!(
                matches!(&err, Error::Invalid { what, .. } if *what == name),
                "dataset {name}, scalar {scalar}: {err}"
            );
        }
    }
}

#[test]
fn metadata_and_pixel_columns_reject_incompatible_types() {
    for name in [
        "chroms/length",
        "bins/chrom",
        "bins/start",
        "bins/end",
        "indexes/chrom_offset",
    ] {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("float-metadata.cool");
        write_cool(&path, &spec(10, 72));
        {
            let file = hdf5_metno::File::open_rw(&path).unwrap();
            let n = file.dataset(name).unwrap().size();
            file.unlink(name).unwrap();
            file.new_dataset::<f64>().shape(n).create(name).unwrap();
        }
        let err = Cool::open(&path, "/")
            .err()
            .expect("float metadata accepted");
        assert!(
            matches!(&err, Error::Invalid { what, .. } if *what == name),
            "{name}: {err}"
        );
    }
    for name in ["pixels/bin1_id", "pixels/bin2_id", "pixels/count"] {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("string-pixels.cool");
        write_cool(&path, &spec(10, 73));
        {
            let file = hdf5_metno::File::open_rw(&path).unwrap();
            let n = file.dataset(name).unwrap().size();
            file.unlink(name).unwrap();
            file.new_dataset::<hdf5_metno::types::FixedAscii<8>>()
                .shape(n)
                .create(name)
                .unwrap();
        }
        let err = Cool::open(&path, "/")
            .err()
            .expect("string pixels accepted");
        assert!(
            matches!(&err, Error::Invalid { what, .. } if *what == name),
            "{name}: {err}"
        );
    }
}

#[test]
fn related_columns_must_have_matching_lengths() {
    for (name, what) in [
        ("chroms/length", "cool chroms"),
        ("bins/start", "cool bins"),
        ("bins/end", "cool bins"),
        ("pixels/bin2_id", "cool pixels"),
        ("pixels/count", "cool pixels"),
        ("indexes/chrom_offset", "indexes/chrom_offset"),
    ] {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("unequal-lengths.cool");
        write_cool(&path, &spec(10, 74));
        {
            let file = hdf5_metno::File::open_rw(&path).unwrap();
            let n = file.dataset(name).unwrap().size();
            file.unlink(name).unwrap();
            file.new_dataset::<i64>().shape(n + 1).create(name).unwrap();
        }
        let err = Cool::open(&path, "/")
            .err()
            .expect("unequal lengths accepted");
        assert!(
            matches!(&err, Error::Invalid { what: actual, .. } if *actual == what),
            "{name}: {err}"
        );
    }
}

#[test]
fn variable_bins_are_unsupported() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("v.cool");
    write_cool(&path, &spec(10, 15));
    {
        let file = hdf5_metno::File::open_rw(&path).unwrap();
        file.attr("bin-type")
            .unwrap()
            .write_scalar(
                &"variable"
                    .parse::<hdf5_metno::types::VarLenUnicode>()
                    .unwrap(),
            )
            .unwrap();
    }
    assert!(matches!(Cool::open(&path, "/"), Err(Error::Unsupported(_))));
}

#[test]
fn bins_off_the_fixed_grid_are_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("b.cool");
    write_cool(&path, &spec(10, 16));
    {
        let file = hdf5_metno::File::open_rw(&path).unwrap();
        let start = file.dataset("bins/start").unwrap();
        start.write_slice(&[150i32], 1..2).unwrap();
    }
    let err = Cool::open(&path, "/").err().unwrap();
    assert!(
        matches!(
            err,
            Error::Invalid {
                what: "cool bins",
                ..
            }
        ),
        "{err}"
    );
}

#[test]
fn chromosome_offsets_are_checked_before_bin_table_access() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("invalid-offset.cool");
    let s = MapSpec::new(&[("first", 1), ("empty", 0)], 1);
    write_cool(&path, &s);
    {
        let file = hdf5_metno::File::open_rw(&path).unwrap();
        // The declared first chromosome spans two bins, but the table
        // holds one row and the index returns to that final row count.
        file.dataset("chroms/length")
            .unwrap()
            .write_raw(&[2i32, 0])
            .unwrap();
        file.dataset("indexes/chrom_offset")
            .unwrap()
            .write_raw(&[0i64, 2, 1])
            .unwrap();
    }
    let err = Cool::open(&path, "/").err().unwrap();
    assert!(
        matches!(
            err,
            Error::Invalid {
                what: "indexes/chrom_offset",
                ..
            }
        ),
        "{err}"
    );
}

#[test]
fn bin_validation_crosses_batches_and_empty_chromosomes() {
    for chrom_enum in [false, true] {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("batched.cool");
        let mut s = MapSpec::new(
            &[
                ("empty-first", 0),
                ("first", 65_535),
                ("empty-middle", 0),
                ("second", 65_539),
                ("empty-last", 0),
            ],
            1,
        );
        s.chrom_enum = chrom_enum;
        write_cool(&path, &s);
        let cool = Cool::open(&path, "/").unwrap();
        assert_eq!(cool.n_bins(), 131_074);
        assert_eq!(cool.locate(65_534), (1, 65_534));
        assert_eq!(cool.locate(65_535), (3, 0));
        assert_eq!(cool.locate(131_073), (3, 65_538));
        drop(cool);

        // Exercise both the first row of a new batch and the final short
        // batch; a valid prefix must not hide invalid later metadata.
        for row in [65_536usize, 131_073] {
            for (name, expected) in [("start", row - 65_535), ("end", row - 65_534)] {
                {
                    let file = hdf5_metno::File::open_rw(&path).unwrap();
                    file.dataset(&format!("bins/{name}"))
                        .unwrap()
                        .write_slice(&[-1i32], row..row + 1)
                        .unwrap();
                }
                assert!(matches!(
                    Cool::open(&path, "/"),
                    Err(Error::Invalid {
                        what: "cool bins",
                        ..
                    })
                ));
                let file = hdf5_metno::File::open_rw(&path).unwrap();
                file.dataset(&format!("bins/{name}"))
                    .unwrap()
                    .write_slice(&[expected as i32], row..row + 1)
                    .unwrap();
            }
        }
        assert!(Cool::open(&path, "/").is_ok());
    }
}

#[test]
fn float_bin_ids_are_rejected_even_with_custom_hdf5_layout() {
    use hdf5_metno as hdf5;
    use hdf5_metno_sys as sys;

    for custom_bias in [false, true] {
        for name in ["pixels/bin1_id", "pixels/bin2_id"] {
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("float-ids.cool");
            let mut s = MapSpec::new(&[("chr1", 10)], 1);
            s.pixels = vec![(0, 1, 2.0)];
            write_cool(&path, &s);
            {
                let file = hdf5::File::open_rw(&path).unwrap();
                file.unlink(name).unwrap();
                let dtype = {
                    let _lock = sys::LOCK.lock();
                    // SAFETY: an open file initialized HDF5; the global
                    // lock protects the predefined type and its new copy.
                    let id = unsafe { sys::h5t::H5Tcopy(*sys::h5t::H5T_IEEE_F32LE) };
                    assert!(id >= 0);
                    // SAFETY: transfer the newly owned valid type once to
                    // the RAII wrapper while holding the HDF5 lock.
                    let dtype: hdf5::Datatype = unsafe { hdf5::from_id(id) }.unwrap();
                    if custom_bias {
                        // SAFETY: this transient copy is exclusively owned
                        // and not yet attached to a dataset; lock is held.
                        assert!(unsafe { sys::h5t::H5Tset_ebias(dtype.id(), 126) } >= 0);
                    }
                    dtype
                };
                file.new_dataset_builder()
                    .empty_as(&dtype)
                    .shape(1)
                    .create(name)
                    .unwrap()
                    .write_raw(&[0.75f32])
                    .unwrap();
            }
            assert!(matches!(
                Cool::open(&path, "/"),
                Err(Error::Invalid { what, .. }) if what == name
            ));
        }
    }
}
