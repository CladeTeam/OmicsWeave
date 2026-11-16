use std::ffi::{c_char, c_int, c_uint, c_void};
use std::ptr::NonNull;

use htscodecs_sys as h;

#[repr(C)]
struct FqzSlice {
    num_records: c_int,
    len: *mut u32,
    flags: *mut u32,
}

unsafe extern "C" {
    fn free(p: *mut c_void);
    fn rans_compress(input: *mut u8, len: c_uint, out: *mut c_uint, order: c_int) -> *mut u8;
    fn rans_compress_4x16(input: *mut u8, len: c_uint, out: *mut c_uint, order: c_int) -> *mut u8;
    fn arith_compress(input: *mut u8, len: c_uint, out: *mut c_uint, order: c_int) -> *mut u8;
    fn fqz_compress(
        version: c_int,
        records: *mut FqzSlice,
        input: *mut c_char,
        len: usize,
        out: *mut usize,
        strategy: c_int,
        parameters: *mut c_void,
    ) -> *mut c_char;
    fn tok3_encode_names(
        input: *mut c_char,
        len: c_int,
        level: c_int,
        arithmetic: c_int,
        out: *mut c_int,
        consumed: *mut c_int,
    ) -> *mut u8;
}

struct CBuffer {
    ptr: NonNull<u8>,
    len: usize,
}

impl CBuffer {
    unsafe fn owned(ptr: *mut u8, len: usize) -> Option<Self> {
        NonNull::new(ptr).map(|ptr| Self { ptr, len })
    }

    fn bytes(&self) -> &[u8] {
        // SAFETY: private construction takes a successful codec allocation
        // and its returned initialized byte count; this owner has not freed it.
        unsafe { std::slice::from_raw_parts(self.ptr.as_ptr(), self.len) }
    }
}

impl Drop for CBuffer {
    fn drop(&mut self) {
        // SAFETY: the private owner holds one unique malloc allocation.
        unsafe { free(self.ptr.as_ptr().cast()) }
    }
}

#[derive(Clone, Copy, Debug)]
enum Codec {
    Rans,
    Nx16,
    Arithmetic,
    Fqz,
    Names,
}

const CODECS: [Codec; 5] = [
    Codec::Rans,
    Codec::Nx16,
    Codec::Arithmetic,
    Codec::Fqz,
    Codec::Names,
];

fn encode(codec: Codec, input: &[u8], order: c_int) -> CBuffer {
    let mut input = input.to_vec();
    let len = c_uint::try_from(input.len()).unwrap();
    let mut count = 0;
    // SAFETY: each encoder receives its actual mutable input extent and
    // distinct writable length fields. FQZ record lengths exactly partition
    // the qualities; name input consists of complete NUL-terminated names.
    let (ptr, count) = unsafe {
        match codec {
            Codec::Rans => (
                rans_compress(input.as_mut_ptr(), len, &mut count, order),
                count as usize,
            ),
            Codec::Nx16 => (
                rans_compress_4x16(input.as_mut_ptr(), len, &mut count, order),
                count as usize,
            ),
            Codec::Arithmetic => (
                arith_compress(input.as_mut_ptr(), len, &mut count, order),
                count as usize,
            ),
            Codec::Fqz => {
                let mut lengths = [input.len() as u32];
                let mut flags = [0];
                let mut records = FqzSlice {
                    num_records: 1,
                    len: lengths.as_mut_ptr(),
                    flags: flags.as_mut_ptr(),
                };
                let mut count = 0;
                (
                    fqz_compress(
                        3,
                        &mut records,
                        input.as_mut_ptr().cast(),
                        input.len(),
                        &mut count,
                        0,
                        std::ptr::null_mut(),
                    )
                    .cast(),
                    count,
                )
            }
            Codec::Names => {
                let (mut count, mut consumed) = (0, 0);
                let ptr = tok3_encode_names(
                    input.as_mut_ptr().cast(),
                    len as c_int,
                    1,
                    order,
                    &mut count,
                    &mut consumed,
                );
                assert_eq!(consumed as usize, input.len());
                assert!(count >= 0);
                (ptr, count as usize)
            }
        }
    };
    // SAFETY: on success these C encoders return malloced initialized buffers.
    unsafe { CBuffer::owned(ptr, count) }
        .unwrap_or_else(|| panic!("fixture must encode: {codec:?}, order={order}"))
}

fn decode(codec: Codec, input: &[u8], size: usize) -> Option<CBuffer> {
    let before = input.to_vec();
    let mut input = input.to_vec();
    let len = c_uint::try_from(input.len()).unwrap();
    if input.is_empty() {
        input.push(0);
    }
    let mut count = 0;
    // SAFETY: even an empty input has an allocated byte; actual input extents and
    // writable distinct length fields are passed. NULL output requests malloc.
    let (ptr, count) = unsafe {
        match codec {
            Codec::Rans => (
                h::rans_uncompress(input.as_mut_ptr(), len, &mut count),
                count as usize,
            ),
            Codec::Nx16 => (
                h::rans_uncompress_to_4x16(
                    input.as_mut_ptr(),
                    len,
                    std::ptr::null_mut(),
                    &mut count,
                ),
                count as usize,
            ),
            Codec::Arithmetic => (
                h::arith_uncompress_to(input.as_mut_ptr(), len, std::ptr::null_mut(), &mut count),
                count as usize,
            ),
            Codec::Fqz => {
                let mut count = size;
                (
                    h::fqz_decompress(
                        input.as_mut_ptr().cast(),
                        len as usize,
                        &mut count,
                        std::ptr::null_mut(),
                        0,
                    )
                    .cast(),
                    count,
                )
            }
            Codec::Names => (
                h::tok3_decode_names(input.as_mut_ptr(), len, &mut count),
                count as usize,
            ),
        }
    };
    assert_eq!(
        &input[..len as usize],
        before,
        "decoder modified its input: {codec:?}, encoded_len={len}"
    );
    // SAFETY: NULL means failure; otherwise the decoder supplies initialized len.
    unsafe { CBuffer::owned(ptr, count) }
}

fn fixture(codec: Codec) -> Vec<u8> {
    match codec {
        Codec::Names => (0..32)
            .flat_map(|i| format!("instrument:1:lane:{i}:ACGT\0").into_bytes())
            .collect(),
        Codec::Fqz => (0..512)
            .map(|i| 20 + ((i * 13 + i / 7) % 22) as u8)
            .collect(),
        _ => (0..2048).map(|i| ((i * 13 + i / 7) % 17) as u8).collect(),
    }
}

fn fixture_cases(codec: Codec) -> Vec<(&'static str, Vec<u8>)> {
    let (repeated, varied) = match codec {
        Codec::Names => (
            b"same:read:0001\0".repeat(32),
            (0..32)
                .flat_map(|i| format!("run-{i:04}:lane_{}:ACGT/{}\0", i % 3, i % 2).into_bytes())
                .collect(),
        ),
        Codec::Fqz => (
            vec![32; 512],
            (0..512).map(|i| ((i * 29 + i / 3) % 42) as u8).collect(),
        ),
        _ => (
            vec![42; 2048],
            (0..2048).map(|i| ((i * 29 + i / 3) % 256) as u8).collect(),
        ),
    };
    vec![
        ("structured", fixture(codec)),
        ("repeated", repeated),
        ("varied", varied),
    ]
}

#[test]
fn allocated_decoders_preserve_bytes_across_codecs_orders_and_data_shapes() {
    for codec in CODECS {
        for (shape, input) in fixture_cases(codec) {
            for order in [0, 1] {
                let encoded = encode(codec, &input, order);
                let decoded = decode(codec, encoded.bytes(), input.len()).unwrap_or_else(|| {
                    panic!("decode failed: {codec:?}, order={order}, shape={shape}")
                });
                assert_eq!(
                    decoded.bytes(),
                    input,
                    "{codec:?}, order={order}, shape={shape}"
                );
            }
        }
    }
}

#[test]
fn truncated_streams_preserve_input_and_do_not_poison_later_decodes() {
    for codec in CODECS {
        let input = fixture(codec);
        for order in [0, 1] {
            let encoded = encode(codec, &input, order);
            let mut cuts: Vec<_> = (0..encoded.len.min(9)).collect();
            cuts.extend([encoded.len / 2, encoded.len - 1]);
            cuts.sort_unstable();
            cuts.dedup();
            for cut in cuts {
                if let Some(decoded) = decode(codec, &encoded.bytes()[..cut], input.len()) {
                    // Some formats permit a prefix that describes a shorter stream.
                    assert!(
                        decoded.len <= input.len(),
                        "{codec:?}, order={order}, prefix={cut}, decoded_len={}",
                        decoded.len
                    );
                }
                let decoded = decode(codec, encoded.bytes(), input.len()).unwrap_or_else(|| {
                    panic!(
                        "full decode failed after prefix: {codec:?}, order={order}, prefix={cut}"
                    )
                });
                assert_eq!(
                    decoded.bytes(),
                    input,
                    "recovery: {codec:?}, order={order}, prefix={cut}"
                );
            }
        }
    }
}

#[test]
fn caller_buffers_respect_capacity_and_keep_allocation_ownership_with_the_caller() {
    for codec in [Codec::Nx16, Codec::Arithmetic] {
        let input = fixture(codec);
        for order in [0, 1, 4, 5, 0x40, 0x80, 0x20, 0x100 | 8] {
            let encoded = encode(codec, &input, order);
            let mut packed = encoded.bytes().to_vec();
            for capacity in [0, 1, input.len() - 1, input.len(), input.len() + 16] {
                let mut out = vec![0x9b; capacity + 2];
                let mut len = capacity as c_uint;
                // SAFETY: the interior output slice has exactly the stated
                // writable capacity, and input/length storage do not alias it.
                let ptr = unsafe {
                    let out = out.as_mut_ptr().add(1);
                    match codec {
                        Codec::Nx16 => h::rans_uncompress_to_4x16(
                            packed.as_mut_ptr(),
                            packed.len() as c_uint,
                            out,
                            &mut len,
                        ),
                        Codec::Arithmetic => h::arith_uncompress_to(
                            packed.as_mut_ptr(),
                            packed.len() as c_uint,
                            out,
                            &mut len,
                        ),
                        _ => unreachable!(),
                    }
                };
                assert_eq!(
                    (out[0], out[capacity + 1]),
                    (0x9b, 0x9b),
                    "output guards: {codec:?}, order={order}, capacity={capacity}"
                );
                assert_eq!(
                    packed,
                    encoded.bytes(),
                    "input changed: {codec:?}, order={order}, capacity={capacity}"
                );
                // Striped streams require the exact decoded size on entry.
                let striped = encoded.bytes()[0] & 0x08 != 0;
                let accepted = capacity == input.len() || (capacity > input.len() && !striped);
                if accepted {
                    assert_eq!(
                        ptr,
                        out[1..].as_mut_ptr(),
                        "output ownership: {codec:?}, order={order}, capacity={capacity}"
                    );
                    assert_eq!(
                        len as usize,
                        input.len(),
                        "decoded size: {codec:?}, order={order}, capacity={capacity}"
                    );
                    assert_eq!(
                        &out[1..input.len() + 1],
                        input,
                        "decoded bytes: {codec:?}, order={order}, capacity={capacity}"
                    );
                } else {
                    // Failure does not establish output bytes or an output length.
                    assert!(
                        ptr.is_null(),
                        "invalid output size accepted: {codec:?}, order={order}, capacity={capacity}"
                    );
                }
            }
        }
    }
}

#[test]
fn fqz_record_lengths_are_written_only_inside_the_requested_extent() {
    let input = fixture(Codec::Fqz);
    let encoded = encode(Codec::Fqz, &input, 0);
    for count in [0, 1] {
        let mut packed = encoded.bytes().to_vec();
        let mut lengths = [-97; 3];
        let mut size = input.len();
        // SAFETY: both buffers and size are disjoint; the interior lengths slot
        // is writable for count entries. The two exterior slots are guards.
        let ptr = unsafe {
            h::fqz_decompress(
                packed.as_mut_ptr().cast(),
                packed.len(),
                &mut size,
                lengths.as_mut_ptr().add(1),
                count,
            )
        };
        // SAFETY: successful fqz output is one malloc allocation of size bytes.
        let decoded = unsafe { CBuffer::owned(ptr.cast(), size) }
            .unwrap_or_else(|| panic!("FQZ decode failed with lengths capacity={count}"));
        assert_eq!(decoded.bytes(), input, "lengths capacity={count}");
        assert_eq!(packed, encoded.bytes(), "lengths capacity={count}");
        assert_eq!(lengths[0], -97, "leading guard, lengths capacity={count}");
        assert_eq!(lengths[2], -97, "trailing guard, lengths capacity={count}");
        assert_eq!(
            lengths[1],
            if count == 0 {
                -97
            } else {
                input.len() as c_int
            }
        );
    }
}

#[test]
fn concurrent_decoders_and_thread_local_scratch_preserve_independent_streams() {
    let start = std::sync::Barrier::new(8);
    std::thread::scope(|scope| {
        for worker in 0..8 {
            let start = &start;
            scope.spawn(move || {
                let cases: Vec<_> = CODECS
                    .into_iter()
                    .map(|codec| (codec, fixture_cases(codec).swap_remove(worker % 3).1))
                    .collect();
                start.wait();
                for (codec, input) in cases {
                    let orders: &[c_int] = match codec {
                        Codec::Nx16 => &[0, 1, 4, 5],
                        _ => &[0, 1],
                    };
                    for &order in orders {
                        let packed = encode(codec, &input, order);
                        let decoded =
                            decode(codec, packed.bytes(), input.len()).unwrap_or_else(|| {
                                panic!("worker={worker}, codec={codec:?}, order={order}")
                            });
                        assert_eq!(
                            decoded.bytes(),
                            input,
                            "worker={worker}, codec={codec:?}, order={order}"
                        );
                    }
                }
            });
        }
    });
}
