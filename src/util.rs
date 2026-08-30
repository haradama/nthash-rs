//! Miscellaneous helpers shared across all ntHash variants.
//!
//! This module provides two small yet critical building blocks used by
//! every ntHash implementation (`NtHash`, `SeedNtHash`, `BlindNtHash`, etc.):
//!
//! - **`canonical`** - combine forward and reverse‐complement hashes into a
//!   strand‐independent value by wrapping addition.
//! - **`extend_hashes`** - generate a sequence of "extra" hash values from
//!   one canonical base hash, matching the C++ reference's multiplicative
//!   mixing and shift scheme.
//!
//! These functions are marked `#[inline]` for zero‐overhead calls in hot paths,
//! and the code is dependency‐free (only `core`/`std`), so it can be used
//! in no‐std contexts if needed.

use crate::constants::{MULTISEED, MULTISHIFT};

/// Combine forward and reverse‐complement strand hashes into a single
/// *canonical* k‑mer hash (strand‐independent).
///
/// The original ntHash definition simply **adds** the two 64‑bit words with
/// wrapping arithmetic to remain well‐defined on overflow.
///
/// # Examples
///
/// ```
/// # use nthash_rs::util::canonical;
/// let fwd = 0xFFFF_FFFF_FFFF_FFFF;
/// let rev = 1;
/// // wraps around to 0
/// assert_eq!(canonical(fwd, rev), 0);
/// ```
#[inline(always)]
pub const fn canonical(fwd: u64, rev: u64) -> u64 {
    fwd.wrapping_add(rev)
}

/// Expand a single canonical hash into a user‐provided slice of additional
/// hash values.
///
/// This implements the same scheme as the C++ ntHash reference:
/// each extra hash `h_i` (for `i ≥ 1`) is computed as:
///
/// ```text
///   mix   = (i as u64) ^ (k as u64 * MULTISEED)
///   h_i   = base.wrapping_mul(mix)
///   h_i  ^= h_i >> MULTISHIFT
/// ```
///
/// - `fwd`, `rev`  - forward and reverse‐complement strand hashes.
/// - `k`           - k‑mer span or seed weight, used in the mixing step.
/// - `hashes`      - output slice; the length determines how many values
///   (including the canonical hash at index 0) are generated.
///
/// If `hashes` is empty, this function returns immediately, avoiding any
/// unnecessary branching in callers.
///
/// # Examples
///
/// ```
/// # use nthash_rs::util::extend_hashes;
/// let mut out = [0u64; 4];
/// extend_hashes(0x1234, 0x5678, 5, &mut out);
/// assert_eq!(out[0], 0x1234u64.wrapping_add(0x5678));
/// // subsequent elements are nonzero mixes
/// assert!(out[1] != out[0]);
/// ```
#[inline]
pub fn extend_hashes(fwd: u64, rev: u64, k: u32, hashes: &mut [u64]) {
    extend_from_canonical(canonical(fwd, rev), k, hashes)
}

/// [`extend_hashes`] starting from an already-combined canonical hash.
///
/// Every extra hash is a pure function of the canonical value, so a driver that
/// already has it — every `for_each_*` loop in this crate — should not
/// re-derive it from the strand pair.
#[inline]
pub fn extend_from_canonical(base: u64, k: u32, hashes: &mut [u64]) {
    match hashes.len() {
        0 => return,
        1 => {
            hashes[0] = base;
            return;
        }
        _ => {}
    }
    hashes[0] = base;

    let seed = (k as u64).wrapping_mul(MULTISEED);

    // Compute extra hashes for i = 1 .. len−1
    for (i, slot) in hashes.iter_mut().enumerate().skip(1) {
        let mut h = base.wrapping_mul((i as u64) ^ seed);
        h ^= h >> MULTISHIFT;
        *slot = h;
    }
}

/// Const-sized [`extend_from_canonical`].
///
/// With `N` known at compile time the mixing loop unrolls and the buffer can
/// stay in registers, so a caller that knows its hash count statically pays
/// neither the trip through a heap slice nor a loop whose trip count the
/// optimizer cannot see. Used by
/// [`NtHashIter::for_each_hash_n`](crate::kmer::NtHashIter::for_each_hash_n).
#[inline(always)]
pub fn extend_from_canonical_n<const N: usize>(base: u64, k: u32, hashes: &mut [u64; N]) {
    if N == 0 {
        return;
    }
    hashes[0] = base;
    if N == 1 {
        return;
    }
    let seed = (k as u64).wrapping_mul(MULTISEED);
    for (i, slot) in hashes.iter_mut().enumerate().skip(1) {
        let mut h = base.wrapping_mul((i as u64) ^ seed);
        h ^= h >> MULTISHIFT;
        *slot = h;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_wraps_on_overflow() {
        let max = u64::MAX;
        assert_eq!(canonical(max, 1), 0);
    }

    #[test]
    fn extend_zero_length_slice() {
        let mut out: [u64; 0] = [];
        extend_hashes(123, 456, 7, &mut out);
        // no panic, no change
    }

    /// The const-sized form must agree with the slice form element for element.
    #[test]
    fn const_sized_matches_slice_form() {
        const F: u64 = 0xFEED_FACE_CAFE_BEEF;
        const R: u64 = 0x0BAD_C0DE_1337_D00D;
        let base = canonical(F, R);
        for k in [1u32, 21, 31] {
            let mut want = [0u64; 8];
            extend_hashes(F, R, k, &mut want);
            let mut got = [0u64; 8];
            extend_from_canonical_n::<8>(base, k, &mut got);
            assert_eq!(got, want, "k={k}");

            let mut one = [0u64; 1];
            extend_from_canonical_n::<1>(base, k, &mut one);
            assert_eq!(one[0], want[0]);
            extend_from_canonical_n::<0>(base, k, &mut []);
        }
    }

    #[test]
    fn extend_matches_cpp_reference() {
        const F: u64 = 0x1234_5678_9ABC_DEF0;
        const R: u64 = 0x0FED_CBA9_8765_4321;
        const K: u32 = 21;
        let mut v = [0u64; 8];
        extend_hashes(F, R, K, &mut v);
        let base = F.wrapping_add(R);
        for (i, &got) in v.iter().enumerate() {
            let expected = if i == 0 {
                base
            } else {
                let mut t = base.wrapping_mul((i as u64) ^ (K as u64).wrapping_mul(MULTISEED));
                t ^= t >> MULTISHIFT;
                t
            };
            assert_eq!(got, expected);
        }
    }
}
