//! Low‑level **split‑rotate** primitives for ntHash.
//!
//! This module implements the core bit‑twiddling operations used by all
//! ntHash variants.  The **split‑rotate** (`srol` / `sror`) operations rotate
//! a 64‑bit word in two independent halves (33 bits + 31 bits) to preserve
//! strand‑symmetry properties.  We also provide a lookup‑based variant
//! (`srol_table`) that applies a split‑rotate to a pre‑seeded constant
//! (A/C/G/T/N) and supports arbitrary rotation distances without branches.
//!
//! All functions are marked `#[inline(always)]` for maximum performance — each
//! compiles down to a handful of shifts, masks, and table lookups.

use crate::constants::{MS_TAB_31L, MS_TAB_33R};

/// One‑bit split‑rotate left (33 + 31 bit halves).
///
/// Conceptually, the 64‑bit word is split into:
/// - a 33‑bit "upper" half (bits 32–63)
/// - a 31‑bit "lower" half (bits 0–30)
///
/// Each half is rotated left by one bit, exchanging carry bits:
/// - bit 63 ➔ bit 33 (upper ➔ lower)
/// - bit 32 ➔ bit 0  (lower ➔ upper)
///
/// This preserves the strand‑symmetry invariants important to ntHash.
#[inline(always)]
pub const fn srol(x: u64) -> u64 {
    // extract the wrap bits from each half
    let m = ((x & 0x8000_0000_0000_0000) >> 30)   // bit 63 ➔ bit 33
          | ((x & 0x0000_0001_0000_0000) >> 32);  // bit 32 ➔ bit 0
    // shift left and re‑insert those bits
    ((x << 1) & 0xFFFF_FFFD_FFFF_FFFF) | m
}

/// Arbitrary‑distance split‑rotate left (0 ≤ d < 64).
///
/// This implements `d` repeated one‑bit split‑rotates efficiently:
/// 1. Perform a full 64‑bit rotate left by `d`.
/// 2. "Unscramble" any bits that crossed the 33/31 boundary to match
///    the effect of split‑rotating each half independently.
#[inline(always)]
pub const fn srol_n(x: u64, d: u32) -> u64 {
    if d == 0 {
        return x;
    }
    // full rotate
    let v = x.rotate_left(d);
    // detect bits that straddle the 33/31 boundary
    let y = (v ^ (v >> 33)) & (!0u64 >> (64 - d));
    // correct their placement
    v ^ (y | (y << 33))
}

/// One‑bit split‑rotate right (33 + 31 bit halves).
///
/// Inverse of [`srol`].  Rotates each half right by one bit, exchanging:
/// - bit 33 ➔ bit 63
/// - bit 0  ➔ bit 32
#[inline(always)]
pub const fn sror(x: u64) -> u64 {
    // extract wrap bits for right rotation
    let m = ((x & 0x0000_0002_0000_0000) << 30)   // bit 33 ➔ bit 63
          | ((x & 0x0000_0000_0000_0001) << 32);  // bit 0  ➔ bit 32
    ((x >> 1) & 0xFFFF_FFFE_FFFF_FFFF) | m
}

/// Lookup‑based split‑rotate left.
///
/// Applies a split‑rotate of distance `d` to the 64‑bit seed constant for
/// nucleotide `c` (A,C,G,T,N).  Internally indexes two pre‑computed tables:
/// - `MS_TAB_31L[c][d % 31]` for the 31‑bit lower half
/// - `MS_TAB_33R[c][d % 33]` for the 33‑bit upper half
///
/// This avoids any runtime loops or branching in the hot path.
#[inline(always)]
pub fn srol_table(c: u8, d: u32) -> u64 {
    let idx31 = (d % 31) as usize;
    let idx33 = (d % 33) as usize;
    MS_TAB_31L[c as usize][idx31] | MS_TAB_33R[c as usize][idx33]
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `srol` followed by `sror` repeatedly should restore the original value.
    #[test]
    fn srol_and_sror_inverse() {
        let mut x = 0xDEADBEEF_DEADBEEF_u64;
        for _ in 0..128 {
            x = srol(x);
            x = sror(x);
        }
        assert_eq!(x, 0xDEADBEEF_DEADBEEF);
    }

    /// `srol_n` should match `d` repeated calls to `srol`.
    #[test]
    fn srol_vs_srol_n() {
        let x = 0x1234_5678_9ABC_DEF0_u64;
        let mut y = x;
        for _ in 0..17 {
            y = srol(y);
        }
        assert_eq!(y, srol_n(x, 17));
    }

    /// `srol_table` must agree with `srol_n` when applied to the base seed.
    #[test]
    fn table_matches_srol_n() {
        let bases = [b'A', b'C', b'G', b'T', b'N'];
        for &b in &bases {
            let seed = srol_table(b, 0);
            for d in 0..64u32 {
                assert_eq!(srol_table(b, d), srol_n(seed, d));
            }
        }
    }
}
