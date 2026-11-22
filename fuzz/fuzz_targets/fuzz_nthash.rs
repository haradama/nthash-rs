//! Fuzz target: no `(seq, k, num_hashes, pos, mask)` may panic or hang.
//!
//! All three hashers are driven to completion over arbitrary bytes and
//! arbitrary mask strings. Construction is allowed to fail (it returns a
//! `Result`); what must never happen is a panic, an out-of-bounds slice, or a
//! non-terminating roll. This is the campaign-grade counterpart to the
//! `no_panic_*` property tests. It covers the backward (`roll_back`) and
//! `peek*` entry points and unbounded minimizer windows as well as the forward
//! pass.

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
    w: usize,
}

fuzz_target!(|input: Input| {
    let Input {
        seq,
        k,
        num_hashes,
        pos,
        mask,
        w,
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
        // Walk back to the start, peeking both ways at every stop. roll_back()
        // must move strictly backward, so it is bounded by the same budget.
        let mut n = 0;
        let mut prev = h.pos() + 1;
        loop {
            h.peek();
            h.peek_back();
            if !h.roll_back() {
                break;
            }
            assert!(h.pos() < prev, "roll_back did not move backward");
            prev = h.pos();
            n += 1;
            assert!(n < budget, "roll_back did not terminate");
        }
    }
    // The backward / peek entry points on a fresh hasher.
    if let Ok(mut h) = nthash_rs::NtHash::new(&seq, k, num_hashes, pos) {
        let mut n = 0;
        while h.roll_back() && n < budget {
            n += 1;
        }
    }
    if let Ok(mut h) = nthash_rs::NtHash::new(&seq, k, num_hashes, pos) {
        h.peek();
        h.peek_back();
    }

    // Streaming BlindNtHash: feed the bytes past the initial window.
    if let Ok(mut h) = nthash_rs::BlindNtHash::new(&seq, k, num_hashes, pos) {
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

    // Canonical (w, k)-minimizers over the same input.
    if let Ok(mins) = nthash_rs::MinimizerBuilder::new(&seq)
        .k(k)
        .w(w)
        .finish()
    {
        let mut n = 0;
        for m in mins {
            std::hint::black_box(m);
            n += 1;
            if n >= budget {
                break;
            }
        }
    }

    // Batch read hasher over reads of mixed lengths. A read shorter than k must
    // reject the whole batch before anything is emitted.
    if let Ok(hasher) = nthash_rs::ReadHasher::new(k) {
        let reads: [&[u8]; 5] = [
            &seq,
            &seq[seq.len() / 2..],
            &seq[..seq.len() / 3],
            &[],
            &seq,
        ];
        let mut calls = 0usize;
        let res = hasher.hash_reads(&reads, |ri, pos, hash| {
            assert!(pos + k as usize <= reads[ri].len(), "k-mer out of range");
            calls += 1;
            std::hint::black_box(hash);
        });
        assert!(res.is_ok() || calls == 0, "emitted hashes for a rejected batch");
        // Without the short reads the same input must be accepted.
        let long: Vec<&[u8]> = reads
            .iter()
            .copied()
            .filter(|r| r.len() >= k as usize)
            .collect();
        assert!(hasher.hash_reads(&long, |_, _, _| {}).is_ok());
    }
});

