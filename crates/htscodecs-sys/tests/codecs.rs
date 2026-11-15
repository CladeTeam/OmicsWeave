use std::ffi::{CStr, c_char, c_int, c_uint, c_void};
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
    unsafe { CBuffer::owned(ptr, count) }.expect("valid fixture must encode")
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
    assert!(
        input[..len as usize] == before,
        "decoder modified its input"
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

#[test]
fn every_codec_roundtrips_and_short_inputs_leave_source_unchanged() {
    for codec in [
        Codec::Rans,
        Codec::Nx16,
        Codec::Arithmetic,
        Codec::Fqz,
        Codec::Names,
    ] {
        let input = fixture(codec);
        for order in [0, 1] {
            let encoded = encode(codec, &input, order);
            let decoded = decode(codec, encoded.bytes(), input.len()).unwrap();
            assert_eq!(decoded.bytes(), input, "{codec:?}, order {order}");
            let mut cuts: Vec<_> = (0..encoded.len.min(9)).collect();
            cuts.extend([encoded.len / 2, encoded.len - 1]);
            cuts.sort_unstable();
            cuts.dedup();
            for cut in cuts {
                if let Some(decoded) = decode(codec, &encoded.bytes()[..cut], input.len()) {
                    // A prefix can itself describe a valid shorter stream;
                    // the format layer must check its expected decoded size.
                    assert!(decoded.len <= input.len(), "prefix {cut}, {codec:?}");
                }
            }
        }
    }
}

#[test]
fn caller_buffers_preserve_guards_and_reject_insufficient_capacity() {
    for codec in [Codec::Nx16, Codec::Arithmetic] {
        let input = fixture(codec);
        for order in [0, 1, 4, 5, 0x40, 0x80, 0x20, 0x100 | 8] {
            let encoded = encode(codec, &input, order);
            let mut packed = encoded.bytes().to_vec();
            for capacity in [input.len(), input.len() - 1] {
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
                assert_eq!((out[0], out[capacity + 1]), (0x9b, 0x9b));
                if capacity == input.len() {
                    assert_eq!(ptr, out[1..].as_mut_ptr());
                    assert_eq!(len as usize, input.len());
                    assert_eq!(&out[1..capacity + 1], input);
                } else {
                    assert!(ptr.is_null(), "{codec:?}, order {order}");
                }
            }
        }
    }
}

#[test]
fn decoder_initialization_and_thread_local_scratch_are_independent() {
    std::thread::scope(|scope| {
        for _ in 0..8 {
            scope.spawn(|| {
                let input = fixture(Codec::Nx16);
                for order in [4, 5] {
                    let packed = encode(Codec::Nx16, &input, order);
                    assert_eq!(
                        decode(Codec::Nx16, packed.bytes(), input.len())
                            .unwrap()
                            .bytes(),
                        input
                    );
                }
            });
        }
    });
    // SAFETY: this no-argument function returns the vendor's static C string.
    let version = unsafe { CStr::from_ptr(h::htscodecs_version()) };
    assert_eq!(version.to_bytes(), b"1.6.7");
}
