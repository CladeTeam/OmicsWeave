//! Public record scanning: logical positions, EOF and error state.
use crate::support::wire::stored;
use omicsweave_bgzf::{EOF_BLOCK, Logical, Scanner, inflate_queue};
use omicsweave_core::Error;

#[test]
fn catalog_and_streaming_scanners_share_positions_across_empty_members() {
    let first = stored(b"ab");
    let second = stored(b"cde");
    let raw = [
        EOF_BLOCK.to_vec(),
        first.clone(),
        EOF_BLOCK.to_vec(),
        second.clone(),
        EOF_BLOCK.to_vec(),
    ]
    .concat();
    for threads in [0, 1, 3] {
        for streaming in [false, true] {
            let queue = inflate_queue(threads);
            let mut scan = if streaming {
                Scanner::new_streaming(&raw[..], queue)
            } else {
                Scanner::new(&raw[..], queue)
            };
            assert_eq!(scan.tell(), Logical::new(0, 0));
            assert!(scan.take(0, "empty").unwrap().is_empty());
            assert_eq!(scan.take(1, "first base").unwrap(), b"a");
            assert_eq!(scan.tell(), Logical::new(1, 1));
            assert_eq!(scan.take(1, "member end").unwrap(), b"b");
            assert_eq!(scan.tell(), Logical::new(2, 0));
            let range = scan.take_range(3, "second member").unwrap();
            assert_eq!(scan.slice(range), b"cde");
            assert_eq!(scan.tell(), Logical::new(4, 0));
            assert_eq!(scan.fill(1).unwrap(), 0);
            assert_eq!(scan.compressed_read(), raw.len() as u64);
            assert_eq!(scan.try_blocks().is_some(), !streaming);
            if let Some(blocks) = scan.try_blocks() {
                assert_eq!(blocks.starts.len(), 5);
                assert_eq!(blocks.text_len, 5);
                assert_eq!(
                    blocks.resolve(Logical::new(1, 1)),
                    (EOF_BLOCK.len() as u64) << 16 | 1
                );
                assert_eq!(
                    blocks.resolve(scan.tell()),
                    ((raw.len() - EOF_BLOCK.len()) as u64) << 16
                );
            }
        }
    }
}

#[test]
fn scanner_truncation_distinguishes_atomic_take_from_consuming_skip() {
    let raw = stored(b"abc");
    for skip in [false, true] {
        let mut scan = Scanner::new(&raw[..], inflate_queue(1));
        let error = if skip {
            scan.skip(5, "frame").unwrap_err()
        } else {
            scan.take(5, "frame").unwrap_err()
        };
        assert!(matches!(
            error,
            Error::Truncated {
                need: 5,
                have: 3,
                ..
            }
        ));
        if skip {
            assert_eq!(scan.tell(), Logical::new(1, 0));
            assert_eq!(scan.fill(1).unwrap(), 0);
        } else {
            assert_eq!(scan.tell(), Logical::new(0, 0));
            assert_eq!(scan.take(3, "remaining frame").unwrap(), b"abc");
        }
    }
}

#[test]
fn scanner_decode_failure_remains_terminal_after_a_valid_frame() {
    let prefix = stored(b"ok");
    let mut broken = stored(b"bad");
    let crc = broken.len() - 8;
    broken[crc] ^= 1;
    let raw = [prefix, broken, stored(b"unreachable")].concat();
    for streaming in [false, true] {
        let queue = inflate_queue(1);
        let mut scan = if streaming {
            Scanner::new_streaming(&raw[..], queue)
        } else {
            Scanner::new(&raw[..], queue)
        };
        assert_eq!(scan.take(2, "valid frame").unwrap(), b"ok");
        let at = scan.tell();
        assert!(scan.fill(1).is_err());
        assert!(scan.take(1, "failed frame").is_err());
        assert!(scan.skip(0, "failed frame").is_err());
        assert_eq!(scan.tell(), at);
    }
}
