//! Fuzz target: no `(seq, k, num_hashes, pos, mask)` may panic or hang.
//!
//! All three hashers are driven to completion over arbitrary bytes and
//! arbitrary mask strings. Construction is allowed to fail (it returns a
//! `Result`); what must never happen is a panic, an out-of-bounds slice, or a
//! non-terminating roll. This is the campaign-grade counterpart to the
//! `no_panic_*` property tests, and reproduces the two bugs fixed in the same
//! change set (SeedNtHash stopping at `N`, BlindNtHash short-input panic).

#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;

#[derive(Arbitrary, Debug)]
struct Input {
    seq: Vec<u8>,
    k: u16,
    num_hashes: u8,
    pos: u32,
    mask: String,
}

fuzz_target!(|input: Input| {
    let Input {
        seq,
        k,
        num_hashes,
        pos,
        mask,
    } = input;
    let pos = pos as usize;
    // A hard backstop on iteration count: roll() advances position monotonically,
    // so this can only ever be hit by a genuine non-termination bug.
    let budget = seq.len() + 8;

    // Contiguous NtHash.
    if let Ok(mut h) = nthash_rs::NtHash::new(&seq, k, num_hashes, pos) {
        let mut n = 0;
        while h.roll() && n < budget {
            n += 1;
            std::hint::black_box(h.hashes());
        }
    }

    // Streaming BlindNtHash: feed the bytes past the initial window.
    if let Ok(mut h) = nthash_rs::BlindNtHash::new(&seq, k, num_hashes, pos as isize) {
        for &c in seq.iter().skip(pos.saturating_add(k as usize)) {
            h.roll(c);
            std::hint::black_box(h.hashes());
        }
    }

    // Spaced-seed SeedNtHash with an arbitrary mask string.
    if let Ok(mut h) = nthash_rs::SeedNtHash::new(
        &seq,
        std::slice::from_ref(&mask),
        num_hashes as usize,
        k,
        pos,
    ) {
        let mut n = 0;
        while h.roll() && n < budget {
            n += 1;
            std::hint::black_box(h.hashes());
        }
    }
});
