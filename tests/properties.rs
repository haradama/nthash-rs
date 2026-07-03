//! Property-based tests encoding the mathematical invariants of ntHash.
//!
//! Unlike the golden-vector regression tests (which anchor the lookup tables to
//! the C++ reference on a handful of fixed inputs), these tests assert
//! *structural* properties that must hold for **every** input. They are an
//! oracle-free correctness net: they caught nothing to write home about on
//! their own, but they make whole classes of future regressions impossible to
//! merge — including the two bugs fixed alongside them
//! (`SeedNtHash` stopping at `N`, `BlindNtHash` panicking on short input).
//!
//! The five invariants:
//! 1. **Purity** — a k-mer's hash depends only on its bases, not its position
//!    or neighbours (rolling hash == isolated recomputation).
//! 2. **Reverse-complement invariance** — the canonical hash is identical for a
//!    k-mer and its reverse complement (the defining property of ntHash).
//! 3. **Strand symmetry** — `forward(kmer) == reverse(revcomp(kmer))`.
//! 4. **Roll invertibility** — `roll` then `roll_back` restores the state.
//! 5. **Cross-variant agreement** — `NtHash`, `BlindNtHash`, and an all-ones
//!    `SeedNtHash` compute the same values.
//!
//! A sixth test asserts **robustness**: no `(seq, k, pos, mask)` — however
//! malformed — may panic.

use nthash_rs::{
    canonical,
    kmer::{base_forward_hash, base_reverse_hash},
    BlindNtHash, MinimizerBuilder, NtHash, NtHashBuilder, SeedNtHash,
};
use proptest::prelude::*;
use std::collections::HashSet;

/// Reverse complement of an A/C/G/T sequence.
fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => other,
        })
        .collect()
}

/// Strategy: a non-empty A/C/G/T sequence paired with a valid `k` in `1..=len`.
fn seq_and_k() -> impl Strategy<Value = (Vec<u8>, u16)> {
    prop::collection::vec(prop::sample::select(&b"ACGT"[..]), 1..256).prop_flat_map(|seq| {
        let len = seq.len();
        (Just(seq), 1u16..=(len.min(u16::MAX as usize) as u16))
    })
}

/// Brute-force reference for canonical (w, k)-minimizers: the window minimum
/// (leftmost on ties) over each gap-free run of valid k-mers, deduplicated by
/// position, with the dedup state reset at each run boundary.
fn brute_minimizers(seq: &[u8], k: u16, w: usize) -> Vec<(usize, u64)> {
    let kmers: Vec<(usize, u64)> = NtHashBuilder::new(seq)
        .k(k)
        .num_hashes(1)
        .pos(0)
        .finish()
        .unwrap()
        .map(|(p, h)| (p, h[0]))
        .collect();

    let mut out = Vec::new();
    let mut i = 0;
    while i < kmers.len() {
        let mut j = i;
        while j + 1 < kmers.len() && kmers[j + 1].0 == kmers[j].0 + 1 {
            j += 1;
        }
        let run = &kmers[i..=j];
        let mut last: Option<usize> = None;
        for start in 0..=run.len().saturating_sub(w) {
            if run.len() < w {
                break;
            }
            let m = run[start..start + w]
                .iter()
                .copied()
                .min_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)))
                .unwrap();
            if last != Some(m.0) {
                last = Some(m.0);
                out.push(m);
            }
        }
        i = j + 1;
    }
    out
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(400))]

    /// (1) The rolling forward/reverse hash of each window equals the hash of
    /// that window computed in isolation — the hash is a pure function of the
    /// k-mer, independent of position and surrounding bases.
    #[test]
    fn purity_rolling_equals_isolated((seq, k) in seq_and_k()) {
        let mut h = NtHash::new(&seq, k, 1, 0).unwrap();
        while h.roll() {
            let w = &seq[h.pos()..h.pos() + k as usize];
            prop_assert_eq!(h.forward_hash(), base_forward_hash(w, k));
            prop_assert_eq!(h.reverse_hash(), base_reverse_hash(w, k));
        }
    }

    /// (2) The canonical hash is invariant under reverse complement, and (3)
    /// forward/reverse hashes swap under reverse complement (strand symmetry).
    #[test]
    fn reverse_complement_symmetry((seq, k) in seq_and_k()) {
        let mut h = NtHash::new(&seq, k, 1, 0).unwrap();
        while h.roll() {
            let w = &seq[h.pos()..h.pos() + k as usize];
            let rc = revcomp(w);
            let (fwd_rc, rev_rc) = (base_forward_hash(&rc, k), base_reverse_hash(&rc, k));

            // (2) canonical is strand-independent
            prop_assert_eq!(
                canonical(h.forward_hash(), h.reverse_hash()),
                canonical(fwd_rc, rev_rc)
            );
            // (3) forward(kmer) == reverse(revcomp(kmer)) and vice versa
            prop_assert_eq!(h.forward_hash(), rev_rc);
            prop_assert_eq!(h.reverse_hash(), fwd_rc);
        }
    }

    /// (4) Rolling forward then back returns to exactly the prior hashes.
    #[test]
    fn roll_back_inverts_roll((seq, k) in seq_and_k()) {
        // Need at least two windows to exercise a forward+back cycle.
        prop_assume!(seq.len() > k as usize);

        let mut h = NtHash::new(&seq, k, 1, 0).unwrap();
        prop_assert!(h.roll()); // first window
        let (f0, r0) = (h.forward_hash(), h.reverse_hash());
        prop_assert!(h.roll()); // advance one
        prop_assert!(h.roll_back()); // and back
        prop_assert_eq!(h.forward_hash(), f0);
        prop_assert_eq!(h.reverse_hash(), r0);
    }

    /// (5a) `NtHash` and `BlindNtHash` agree on every window of an N-free input.
    #[test]
    fn nthash_agrees_with_blind((seq, k) in seq_and_k()) {
        let mut nt = NtHash::new(&seq, k, 1, 0).unwrap();
        let mut bl = BlindNtHash::new(&seq, k, 1, 0).unwrap();

        prop_assert!(nt.roll()); // NtHash initialises on first roll
        prop_assert_eq!(nt.forward_hash(), bl.forward_hash());
        prop_assert_eq!(nt.reverse_hash(), bl.reverse_hash());

        for &c in &seq[k as usize..] {
            bl.roll(c);
            prop_assert!(nt.roll());
            prop_assert_eq!(nt.forward_hash(), bl.forward_hash());
            prop_assert_eq!(nt.reverse_hash(), bl.reverse_hash());
        }
    }

    /// (5b) An all-ones spaced seed reproduces `NtHash`'s canonical hash.
    #[test]
    fn seed_all_ones_equals_nthash((seq, k) in seq_and_k()) {
        let mask = "1".repeat(k as usize);
        let mut sd = SeedNtHash::new(&seq, &[mask], 1, k, 0).unwrap();
        let mut nt = NtHash::new(&seq, k, 1, 0).unwrap();
        while nt.roll() && sd.roll() {
            prop_assert_eq!(
                sd.hashes()[0],
                canonical(nt.forward_hash(), nt.reverse_hash())
            );
        }
    }

    /// Minimizers match an independent brute-force window-minimum reference,
    /// for sequences that may contain `N` (exercising the gap reset).
    #[test]
    fn minimizers_match_brute_force(
        seq in prop::collection::vec(prop::sample::select(&b"ACGTN"[..]), 1..256),
        k in 1u16..24,
        w in 1usize..8,
    ) {
        prop_assume!((k as usize) <= seq.len());
        let got: Vec<(usize, u64)> = MinimizerBuilder::new(&seq)
            .k(k)
            .w(w)
            .finish()
            .unwrap()
            .map(|m| (m.pos, m.hash))
            .collect();
        prop_assert_eq!(got, brute_minimizers(&seq, k, w));
    }

    /// The set of selected minimizer hashes is invariant under reverse
    /// complement — canonical minimizers.
    #[test]
    fn minimizer_hashes_invariant_under_revcomp((seq, k) in seq_and_k(), w in 1usize..8) {
        let fwd: HashSet<u64> = MinimizerBuilder::new(&seq)
            .k(k).w(w).finish().unwrap().map(|m| m.hash).collect();
        let rc = revcomp(&seq);
        let rev: HashSet<u64> = MinimizerBuilder::new(&rc)
            .k(k).w(w).finish().unwrap().map(|m| m.hash).collect();
        prop_assert_eq!(fwd, rev);
    }

    /// Robustness: arbitrary bytes / k / pos / num_hashes must never panic —
    /// construction either errors cleanly or yields a hasher that rolls to
    /// completion. This is a proptest stand-in for fuzzing that runs in normal
    /// CI; it fails on the `BlindNtHash` short-input panic fixed today.
    #[test]
    fn no_panic_on_arbitrary_input(
        seq in prop::collection::vec(any::<u8>(), 0..300),
        k in 0u16..300,
        num_hashes in 0u8..6,
        pos in 0usize..400,
        w in 0usize..12,
    ) {
        if let Ok(mut h) = NtHash::new(&seq, k, num_hashes, pos) {
            while h.roll() { std::hint::black_box(h.hashes()); }
        }
        if let Ok(mut h) = BlindNtHash::new(&seq, k, num_hashes, pos as isize) {
            // Feed the remaining bytes; the window slides without bounds panics.
            let start = pos + k as usize;
            for &c in seq.iter().skip(start) {
                h.roll(c);
                std::hint::black_box(h.hashes());
            }
        }
        if let Ok(mins) = MinimizerBuilder::new(&seq).k(k).w(w).finish() {
            for m in mins { std::hint::black_box(m); }
        }
    }

    /// Robustness for spaced seeds: arbitrary bytes and mask strings must not
    /// panic (invalid masks must surface as errors).
    #[test]
    fn no_panic_on_arbitrary_seed(
        seq in prop::collection::vec(any::<u8>(), 0..200),
        k in 1u16..64,
        mask in "[01]{0,80}",
    ) {
        if let Ok(mut h) = SeedNtHash::new(&seq, std::slice::from_ref(&mask), 1, k, 0) {
            while h.roll() { std::hint::black_box(h.hashes()); }
        }
    }
}
