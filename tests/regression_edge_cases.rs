//! Regression tests for edge cases that once panicked, hung, or returned the
//! wrong k-mer. Each test names the behaviour it pins down.

use nthash_rs::{
    BlindNtHash, BlindNtHashBuilder, MinimizerBuilder, NtHash, NtHashBuilder, NtHashError,
    ReadHasher, SeedNtHash, SeedNtHashBuilder, SpacedSeeds,
};

/// Canonical hash of the single k-mer `kmer`.
fn hash_of(kmer: &[u8]) -> u64 {
    let mut h = NtHash::new(kmer, kmer.len() as u16, 1, 0).unwrap();
    assert!(h.roll());
    h.hashes()[0]
}

/// `roll()` running off the end through an `N`-skip used to leave `pos` past
/// the last window, so the next `roll_back()` indexed out of bounds.
#[test]
fn failed_roll_keeps_hasher_on_last_kmer() {
    let mut h = NtHash::new(b"ACGTN", 4, 1, 0).unwrap();
    assert!(h.roll());
    let first = h.hashes()[0];
    assert!(!h.roll());
    assert_eq!(h.pos(), 0);
    assert_eq!(h.hashes()[0], first);
    assert!(!h.roll_back()); // already at the first k-mer; must not panic
    assert!(!h.roll()); // and stays exhausted
    assert_eq!(h.pos(), 0);
}

/// `roll_back()` across an `N` used to re-scan forward and land on the k-mer it
/// started from, so `while h.roll_back()` never terminated.
#[test]
fn roll_back_skips_backward_over_n() {
    let seq = b"ACGTTNACGTA";
    let mut h = NtHash::new(seq, 4, 1, 0).unwrap();
    let mut forward = Vec::new();
    while h.roll() {
        forward.push((h.pos(), h.hashes()[0]));
    }
    assert_eq!(
        forward.iter().map(|&(p, _)| p).collect::<Vec<_>>(),
        [0, 1, 6, 7]
    );

    let mut backward = vec![(h.pos(), h.hashes()[0])];
    while h.roll_back() {
        backward.push((h.pos(), h.hashes()[0]));
        assert!(
            backward.len() <= forward.len(),
            "roll_back did not terminate"
        );
    }
    backward.reverse();
    assert_eq!(backward, forward);
    for &(p, hash) in &forward {
        assert_eq!(hash, hash_of(&seq[p..p + 4]));
    }
}

/// With nothing valid before an `N`, `roll_back()` reports the start of the
/// sequence and leaves the hasher where it was.
#[test]
fn roll_back_with_no_valid_kmer_before_n() {
    let seq = b"ACNACGT";
    let mut h = NtHash::new(seq, 4, 1, 0).unwrap();
    assert!(h.roll());
    assert_eq!(h.pos(), 3);
    assert!(!h.roll_back());
    assert_eq!(h.pos(), 3);
    assert_eq!(h.hashes()[0], hash_of(b"ACGT"));
}

/// On a fresh hasher every entry point initializes on the k-mer at the start
/// position — `roll_back` and `peek` used to step one further.
#[test]
fn first_call_initializes_regardless_of_direction() {
    let seq = b"ACGTTGCA";
    let at2 = hash_of(&seq[2..6]);

    let mut h = NtHash::new(seq, 4, 1, 2).unwrap();
    assert!(h.roll_back());
    assert_eq!((h.pos(), h.hashes()[0]), (2, at2));

    let mut h = NtHash::new(seq, 4, 1, 2).unwrap();
    assert!(h.peek());
    assert_eq!((h.pos(), h.hashes()[0]), (2, at2));

    let mut h = NtHash::new(seq, 4, 1, 2).unwrap();
    assert!(h.peek_back());
    assert_eq!((h.pos(), h.hashes()[0]), (2, at2));

    // At position 0 there is a valid k-mer, so initialization succeeds.
    let mut h = NtHash::new(seq, 4, 1, 0).unwrap();
    assert!(h.roll_back());
    assert_eq!((h.pos(), h.hashes()[0]), (0, hash_of(&seq[..4])));
    assert!(!h.roll_back());
}

/// The explicit-byte peeks follow the same rule: on a fresh hasher they only
/// initialize, whatever byte they are handed.
#[test]
fn peek_char_on_fresh_hasher_initializes() {
    let seq = b"ACGTTGCA";
    let at2 = hash_of(&seq[2..6]);

    let mut h = NtHash::new(seq, 4, 1, 2).unwrap();
    assert!(h.peek_char(b'N'));
    assert_eq!((h.pos(), h.hashes()[0]), (2, at2));

    let mut h = NtHash::new(seq, 4, 1, 2).unwrap();
    assert!(h.peek_back_char(b'N'));
    assert_eq!((h.pos(), h.hashes()[0]), (2, at2));
}

/// Once initialized, `peek` shows the neighbouring k-mer without moving.
#[test]
fn peek_shows_neighbour_without_moving() {
    let seq = b"ACGTTGCA";
    let mut h = NtHash::new(seq, 4, 1, 2).unwrap();
    assert!(h.roll());
    assert!(h.peek());
    assert_eq!((h.pos(), h.hashes()[0]), (2, hash_of(&seq[3..7])));
    assert!(h.peek_back());
    assert_eq!((h.pos(), h.hashes()[0]), (2, hash_of(&seq[1..5])));
}

/// `for_each_canonical` agrees with the iterator across `N`-skips, including a
/// trailing `N` (the path that used to overshoot `pos`), a sequence that ends
/// on a valid k-mer, and one with no valid k-mer at all.
#[test]
fn for_each_canonical_handles_trailing_n() {
    for seq in [
        &b"ACGTN"[..],
        b"ACGTNACGTAN",
        b"NNACGTANN",
        b"ACGTACGTNNN",
        b"ACGTACGT",
        b"NNNNN",
    ] {
        let expected: Vec<(usize, u64)> = NtHashBuilder::new(seq)
            .k(4)
            .finish()
            .unwrap()
            .map(|(p, h)| (p, h[0]))
            .collect();
        let mut got = Vec::new();
        NtHashBuilder::new(seq)
            .k(4)
            .finish()
            .unwrap()
            .for_each_canonical(|p, h| got.push((p, h)));
        assert_eq!(got, expected);
    }
}

/// A window larger than the sequence (up to `usize::MAX`) used to overflow or
/// attempt a `w`-sized allocation.
#[test]
fn minimizer_window_larger_than_sequence() {
    for w in [10, 1 << 45, usize::MAX] {
        let mins = MinimizerBuilder::new(b"ACGTACGTACGT")
            .k(4)
            .w(w)
            .finish()
            .unwrap();
        assert_eq!(mins.count(), 0, "w={w}");
    }
    // Exactly as many k-mers as the window: one minimizer.
    let mins = MinimizerBuilder::new(b"ACGTACGTACGT")
        .k(4)
        .w(9)
        .finish()
        .unwrap();
    assert_eq!(mins.count(), 1);
}

/// Repeated care indices used to desynchronise the ambiguous-base counter, so
/// every window after the first `N` was dropped.
#[test]
fn duplicate_care_indices_are_ignored() {
    let seq = b"ACGTACGTNACGTACGTACGTACGT";
    let collect = |care: Vec<usize>| {
        let mut h = SeedNtHash::from_care_indices(seq, vec![care], 1, 4, 0).unwrap();
        let mut out = Vec::new();
        while h.roll() {
            out.push((h.pos(), h.hashes()[0]));
        }
        out
    };
    let clean = collect(vec![0, 1, 2, 3]);
    assert_eq!(clean.len(), 18);
    assert_eq!(collect(vec![0, 1, 2, 3, 3, 3, 3, 3]), clean);
    assert_eq!(collect(vec![3, 1, 0, 2, 1]), clean);
}

/// An exhausted iterator stays exhausted, whichever way it is driven next.
#[test]
fn exhausted_iterators_stay_exhausted() {
    let seq = b"ACGTNACGTA";

    let drained = || {
        let mut it = NtHashBuilder::new(seq).k(4).finish().unwrap();
        assert_eq!(it.by_ref().count(), 3);
        assert!(it.next().is_none());
        it
    };
    let mut calls = 0;
    drained().for_each_hash(|_, _| calls += 1);
    drained().for_each_canonical(|_, _| calls += 1);
    assert_eq!(calls, 0);
    assert_eq!(drained().canonical().count(), 0);

    let drained = || {
        let mut it = SeedNtHashBuilder::new(seq)
            .k(4)
            .masks(["1011"])
            .finish()
            .unwrap();
        assert!(it.by_ref().count() > 0);
        assert!(it.next().is_none());
        it
    };
    drained().for_each_hash(|_, _| calls += 1);
    assert_eq!(calls, 0);
}

#[test]
fn seed_errors_name_the_actual_problem() {
    let seq = b"ACGTACGT";
    assert_eq!(
        SeedNtHash::new(seq, &["1011"], 1, 0, 0).err(),
        Some(NtHashError::InvalidK)
    );
    assert_eq!(
        SeedNtHash::new(seq, &["1011"], 1, 4, 5).err(),
        Some(NtHashError::PositionOutOfRange { pos: 5, seq_len: 8 })
    );
    assert_eq!(
        SpacedSeeds::from_care_indices(vec![vec![0]], 0).err(),
        Some(NtHashError::InvalidK)
    );
    // The builder surfaces mask errors instead of panicking.
    assert_eq!(
        SeedNtHashBuilder::new(seq)
            .k(4)
            .masks(["10x1"])
            .finish()
            .err(),
        Some(NtHashError::InvalidMask)
    );
    // Sequence too short, through each way of supplying the seeds.
    let too_short = Some(NtHashError::SequenceTooShort { seq_len: 3, k: 4 });
    assert_eq!(
        SeedNtHash::from_care_indices(b"ACG", vec![vec![0, 1]], 1, 4, 0).err(),
        too_short
    );
    let seeds = SpacedSeeds::from_masks(&["1011"], 4).unwrap();
    assert_eq!(
        SeedNtHashBuilder::new(b"ACG").seeds(&seeds).finish().err(),
        too_short
    );
    assert!(SpacedSeeds::from_masks(&[""; 0], 4).unwrap().is_empty());
    assert!(!SpacedSeeds::from_masks(&["1011"], 4).unwrap().is_empty());
    assert_eq!(
        SeedNtHash::new(seq, &["101"], 1, 4, 0).err(),
        Some(NtHashError::InvalidMaskLength { mask_len: 3, k: 4 })
    );
    assert_eq!(
        SeedNtHash::new(seq, &["1x01"], 1, 4, 0).err(),
        Some(NtHashError::InvalidMask)
    );
    assert_eq!(
        SeedNtHash::from_care_indices(seq, vec![vec![0, 4]], 1, 4, 0).err(),
        Some(NtHashError::InvalidCareIndex { index: 4, k: 4 })
    );
    assert_eq!(
        SpacedSeeds::from_masks(&["1"; 0], 0).err(),
        Some(NtHashError::InvalidK)
    );
}

/// Precompiled seeds give the same hashes as per-hasher masks, through both
/// the constructor and the builder.
#[test]
fn shared_seeds_match_owned_masks() {
    let seq = b"ATCGTACGATGCATGNATGCTGACG";
    let masks = ["110011", "010101"];
    let expected: Vec<(usize, Vec<u64>)> = SeedNtHashBuilder::new(seq)
        .k(6)
        .masks(masks)
        .num_hashes(2)
        .finish()
        .unwrap()
        .collect();
    assert!(!expected.is_empty());

    let seeds = SpacedSeeds::from_masks(&masks, 6).unwrap();
    assert_eq!((seeds.k(), seeds.len()), (6, 2));
    let via_builder: Vec<(usize, Vec<u64>)> = SeedNtHashBuilder::new(seq)
        .seeds(&seeds)
        .num_hashes(2)
        .finish()
        .unwrap()
        .collect();
    assert_eq!(via_builder, expected);

    let mut h = SeedNtHash::with_seeds(seq, &seeds, 2, 0).unwrap();
    let mut direct = Vec::new();
    while h.roll() {
        direct.push((h.pos(), h.hashes().to_vec()));
    }
    assert_eq!(direct, expected);

    // A conflicting explicit k is rejected rather than silently ignored.
    assert_eq!(
        SeedNtHashBuilder::new(seq)
            .k(5)
            .seeds(&seeds)
            .finish()
            .err(),
        Some(NtHashError::InvalidK)
    );
    // Sequence shorter than the seeds' k.
    assert_eq!(
        SeedNtHash::with_seeds(b"ACG", &seeds, 1, 0).err(),
        Some(NtHashError::SequenceTooShort { seq_len: 3, k: 6 })
    );
}

/// The blind builder reports an invalid configuration instead of panicking,
/// and its callback traversals continue from a partially consumed iterator.
#[test]
fn blind_builder_errors_and_resumes() {
    assert_eq!(
        BlindNtHashBuilder::new(b"ACG").k(4).finish().err(),
        Some(NtHashError::SequenceTooShort { seq_len: 3, k: 4 })
    );

    let seq = b"ACGTTGCATG";
    let iter = || BlindNtHashBuilder::new(seq).k(4).finish().unwrap();
    let full: Vec<(usize, u64)> = iter().map(|(p, h)| (p, h[0])).collect();

    // Each traversal is driven from a fresh iterator and from one that has
    // already yielded `taken` items.
    for taken in [0, 3] {
        let mut it = iter();
        let mut got: Vec<(usize, u64)> = it.by_ref().take(taken).map(|(p, h)| (p, h[0])).collect();
        it.for_each_hash(|p, h| got.push((p, h[0])));
        assert_eq!(got, full, "for_each_hash after {taken}");

        let mut it = iter();
        let mut got: Vec<(usize, u64)> = it.by_ref().take(taken).map(|(p, h)| (p, h[0])).collect();
        it.for_each_canonical(|p, h| got.push((p, h)));
        assert_eq!(got, full, "for_each_canonical after {taken}");
    }
}

/// `k == 0` is rejected by every constructor.
#[test]
fn zero_k_is_rejected_everywhere() {
    let seq = b"ACGTACGT";
    let invalid_k = Some(NtHashError::InvalidK);
    assert_eq!(NtHash::new(seq, 0, 1, 0).err(), invalid_k);
    assert_eq!(NtHashBuilder::new(seq).finish().err(), invalid_k);
    assert_eq!(BlindNtHash::new(seq, 0, 1, 0).err(), invalid_k);
    assert_eq!(BlindNtHashBuilder::new(seq).finish().err(), invalid_k);
    assert_eq!(ReadHasher::new(0).err(), invalid_k);
    assert_eq!(MinimizerBuilder::new(seq).w(3).finish().err(), invalid_k);
    assert_eq!(
        SeedNtHashBuilder::new(seq).masks(["1"]).finish().err(),
        invalid_k
    );
}

/// A start position past the last k-mer is rejected by every constructor.
#[test]
fn start_position_out_of_range_is_rejected() {
    let seq = b"ACGTACGT";
    let out_of_range = Some(NtHashError::PositionOutOfRange { pos: 5, seq_len: 8 });
    assert_eq!(NtHash::new(seq, 4, 1, 5).err(), out_of_range);
    assert_eq!(
        NtHashBuilder::new(seq).k(4).pos(5).finish().err(),
        out_of_range
    );
    assert_eq!(
        BlindNtHashBuilder::new(seq).k(4).pos(5).finish().err(),
        out_of_range
    );
    assert!(NtHash::new(seq, 4, 1, 4).is_ok());
}

#[test]
fn blind_position_out_of_range_is_reported_as_given() {
    assert_eq!(
        BlindNtHash::new(b"ACGTACGT", 4, 1, 5).err(),
        Some(NtHashError::PositionOutOfRange { pos: 5, seq_len: 8 })
    );
}

/// `seq` with every `from` byte replaced by `to`.
fn swap(seq: &[u8], from: u8, to: u8) -> Vec<u8> {
    seq.iter()
        .map(|&b| if b == from { to } else { b })
        .collect()
}

/// RNA `U` hashes exactly like `T` in every hasher, as in the C++ reference.
/// The seed table used to lack `U`, so `NtHash` / `SeedNtHash` skipped every
/// k-mer containing it and `BlindNtHash` rolled to values that disagreed with
/// a from-scratch hash of the same window.
#[test]
fn rna_u_hashes_like_t() {
    let dna = b"ACGTACGTTGCATTNGATTACAT";
    for rna in [
        swap(dna, b'T', b'U'),
        swap(&dna.to_ascii_lowercase(), b't', b'u'),
    ] {
        let rna = rna.as_slice();
        for k in [1u16, 3, 4, 5, 8] {
            let nt = |s: &[u8]| -> Vec<(usize, Vec<u64>)> {
                NtHashBuilder::new(s)
                    .k(k)
                    .num_hashes(2)
                    .finish()
                    .unwrap()
                    .collect()
            };
            assert!(!nt(dna).is_empty());
            assert_eq!(nt(rna), nt(dna), "NtHash k={k}");

            let mask = "1".repeat(k as usize);
            let spaced = |s: &[u8], mask: &str| -> Vec<(usize, Vec<u64>)> {
                SeedNtHashBuilder::new(s)
                    .k(k)
                    .masks([mask])
                    .finish()
                    .unwrap()
                    .collect()
            };
            assert_eq!(spaced(rna, &mask), spaced(dna, &mask), "SeedNtHash k={k}");
            if k >= 3 {
                let gapped = format!("1{}1", "0".repeat(k as usize - 2));
                assert_eq!(spaced(rna, &gapped), spaced(dna, &gapped), "gapped k={k}");
            }

            let mins = |s: &[u8]| -> Vec<_> {
                MinimizerBuilder::new(s)
                    .k(k)
                    .w(3)
                    .finish()
                    .unwrap()
                    .collect()
            };
            assert_eq!(mins(rna), mins(dna), "Minimizers k={k}");
        }
    }

    // The blind and batch hashers take N-free input.
    let dna = b"ACGTACGTTGCATTGATTACAT";
    let rna = swap(dna, b'T', b'U');
    for k in [1u16, 4, 7] {
        let blind = |s: &[u8]| -> Vec<(usize, u64)> {
            BlindNtHashBuilder::new(s)
                .k(k)
                .finish()
                .unwrap()
                .canonical()
                .collect()
        };
        assert_eq!(blind(&rna), blind(dna), "BlindNtHash k={k}");
        // Rolling agrees with hashing each window from scratch.
        let recomputed: Vec<(usize, u64)> = (0..=rna.len() - k as usize)
            .map(|i| (i, blind(&rna[i..i + k as usize])[0].1))
            .collect();
        assert_eq!(
            blind(&rna),
            recomputed,
            "BlindNtHash roll vs recompute k={k}"
        );

        let batch = |s: &[u8]| {
            let mut out = Vec::new();
            ReadHasher::new(k)
                .unwrap()
                .hash_reads(&[s], |_, pos, h| out.push((pos, h)))
                .unwrap();
            out
        };
        assert_eq!(batch(&rna), batch(dna), "ReadHasher k={k}");
        assert_eq!(batch(&rna), blind(dna), "ReadHasher vs BlindNtHash k={k}");
    }
}
