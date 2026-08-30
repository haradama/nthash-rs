//! Low‑level **split‑rotate** primitives for ntHash.
//!
//! This module implements the core bit‑twiddling operations used by all
//! ntHash variants.  The **split‑rotate** (`srol` / `sror`) operations rotate
//! a 64‑bit word in two independent halves (33 bits + 31 bits) to preserve
//! strand‑symmetry properties.  We also provide a lookup‑based variant
//! (`srol_table`) that applies a split‑rotate to a pre‑seeded constant
//! (A/C/G/T/N) and supports arbitrary rotation distances without branches.
//!
//! All functions are marked `#[inline(always)]` for maximum performance - each
//! compiles down to a handful of shifts, masks, and table lookups.

use crate::constants::{CP_OFF, MS_TAB_31L, MS_TAB_33R, SEED_TAB};

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
          | ((x & 0x0000_0001_0000_0000) >> 32); // bit 32 ➔ bit 0
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

/// [`srol_n`] specialised to a compile-time distance `D` in `1..31`.
///
/// The generic form rotates the whole word and then repairs the bits that
/// crossed the 33/31 boundary — a seven-op chain six cycles deep. Rotating the
/// two fields directly is four cycles deep, and this sits on a *serial* chain:
/// [`base_forward_hash`](crate::kmer::base_forward_hash) is `k/4` of these back
/// to back, which is ~40% of batch hashing for 100 bp reads.
///
/// `D` must be in `1..31`; outside that range the 31-bit field would wrap more
/// than once and the masks below no longer describe the rotation. A unit test
/// pins every valid `D` against [`srol_n`].
#[inline(always)]
pub const fn srol_const<const D: u32>(x: u64) -> u64 {
    let m = (1u64 << D) - 1;
    // Everything that stays inside its own field, then each field's wrap-around.
    ((x << D) & !(m | (m << 33)))       // bits D..32 and 33+D..63
        | ((x >> (33 - D)) & m)         // lower field: bits 33-D..32 -> 0..D-1
        | ((x >> (31 - D)) & (m << 33)) // upper field: bits 64-D..63 -> 33..33+D-1
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
          | ((x & 0x0000_0000_0000_0001) << 32); // bit 0  ➔ bit 32
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

/// The only byte values whose split-rotate tables are non-zero: A/C/G/T (upper
/// and lower case), RNA `U`, and the low indices produced by `byte & CP_OFF`.
/// Every other byte selects the all-zero `N` tables, so `srol_table(c, k) == 0`.
const SROL_NONZERO_BYTES: [u8; 15] = [
    b'A', b'a', b'C', b'c', b'G', b'g', b'T', b't', b'U', b'u', 1, 3, 4, 5, 7,
];

/// Precompute `srol_table(c, k)` for every byte `c` into a 256-entry table.
///
/// Used by [`BlindNtHash`](crate::blind), which — unlike [`NtHash`](crate::kmer)
/// — does not skip ambiguous bases: it hashes whatever it is fed, and the
/// reference defines an `N` (or any non-ACGT byte) as a zero contribution. A
/// full byte-indexed table preserves that, where the 4-constant [`RollConsts`]
/// (which folds `N` into `G`) cannot. Only the ~15 meaningful bytes are
/// computed; the rest stay `0`, which is exactly `srol_table(c, k)` for them. A
/// debug assertion guards that the shortcut stays in sync with `MS_TAB_*`.
pub fn build_srol_k(k: u16) -> Box<[u64; 256]> {
    let d = k as u32;
    let mut t = Box::new([0u64; 256]);
    for &c in &SROL_NONZERO_BYTES {
        t[c as usize] = srol_table(c, d);
    }
    debug_assert!(
        (0..256).all(|c| t[c] == srol_table(c as u8, d)),
        "SROL_NONZERO_BYTES is out of sync with MS_TAB_*"
    );
    t
}

/// Branchless 2-bit base code: folds A/C/G/T (and lowercase / RNA `U`) to
/// `0..3`. Non-ACGT bytes are out of contract (callers gate on `N` first), so
/// their code is unspecified — in particular `N` folds onto `G`, so this must
/// not be used on inputs that may contain ambiguous bases.
///
/// The grouping matches the lookup tables exactly — `a`/`A` share a code and a
/// seed, `t`/`T`/`u`/`U` share a code — so selecting a constant by code is
/// identical to indexing `SEED_TAB`/`srol_table` by the raw byte for every
/// in-contract base.
#[inline(always)]
pub fn base_code(b: u8) -> usize {
    ((b >> 1) & 3) as usize
}

/// Per-base rolling constants for a fixed `k`, indexed by [`base_code`].
///
/// Supersedes the 256-entry `srol_k` table: only four distinct base
/// contributions exist, so four constants per table (128 bytes total, held in
/// registers) are selected by the 2-bit code. This is measurably faster in the
/// roll hot-path than a variable-index load into a 2 KiB table, and — unlike
/// the table — needs no heap allocation, which matters when a fresh hasher is
/// built per short read.
///
/// The four separate tables are what the *backward* roll needs (it pairs them
/// the other way round). The *forward* roll — the hot direction — instead uses
/// the pre-folded [`fwd_q`](Self::fwd_q) / [`rev_q`](Self::rev_q) pair; see
/// [`roll_index`].
#[derive(Clone)]
pub struct RollConsts {
    /// `SEED_TAB[base]` — forward incoming contribution.
    pub seed: [u64; 4],
    /// `srol_table(base, k)` — forward outgoing contribution.
    pub srolk: [u64; 4],
    /// `srol_table(base & CP_OFF, k)` — reverse incoming contribution.
    pub srolk_rc: [u64; 4],
    /// `SEED_TAB[base & CP_OFF]` — reverse outgoing contribution.
    pub seed_rc: [u64; 4],
    /// `srolk[out] ^ seed[in]`, indexed by [`roll_index`]`(out, in)`.
    pub fwd_q: [u64; 16],
    /// `srolk_rc[in] ^ seed_rc[out]`, indexed by the same [`roll_index`].
    pub rev_q: [u64; 16],
    /// `[fwd_q[q], sror(rev_q[q])]` interleaved, for the 128-bit dual-chain
    /// roll (`simd::roll_x2`). The reverse delta is pre-rotated because that
    /// form computes `sror(rev) ^ delta` rather than `sror(rev ^ delta)` —
    /// equal by linearity, and it takes the XOR off the dependency chain.
    pub fr_q: [u64; 32],
}

/// Index into [`RollConsts::fwd_q`] / [`RollConsts::rev_q`] for one roll step:
/// the outgoing base's 2-bit code in bits 2-3, the incoming base's in bits 0-1.
///
/// `base_code(b)` is `(b >> 1) & 3`, so `base_code(out) * 4` is just
/// `(out & 6) << 1` — no shift-then-shift round trip. Non-ACGT bytes are out of
/// contract exactly as for [`base_code`].
#[inline(always)]
pub const fn roll_index(char_out: u8, char_in: u8) -> usize {
    ((char_out as usize & 6) << 1) | ((char_in as usize >> 1) & 3)
}

impl RollConsts {
    /// Precompute the per-base and per-base-pair constants for k-mer length `k`.
    pub fn new(k: u16) -> Self {
        let mut c = RollConsts {
            seed: [0; 4],
            srolk: [0; 4],
            srolk_rc: [0; 4],
            seed_rc: [0; 4],
            fwd_q: [0; 16],
            rev_q: [0; 16],
            fr_q: [0; 32],
        };
        for &b in b"ACGT" {
            let i = base_code(b);
            c.seed[i] = SEED_TAB[b as usize];
            c.srolk[i] = srol_table(b, k as u32);
            c.srolk_rc[i] = srol_table(b & CP_OFF, k as u32);
            c.seed_rc[i] = SEED_TAB[(b & CP_OFF) as usize];
        }
        // Fold each (outgoing, incoming) pair into a single delta. The roll is
        // latency- not throughput-bound, so removing one XOR from each of the
        // forward and reverse dependency chains is the whole scalar win.
        for out in 0..4 {
            for inc in 0..4 {
                c.fwd_q[out * 4 + inc] = c.srolk[out] ^ c.seed[inc];
                c.rev_q[out * 4 + inc] = c.srolk_rc[inc] ^ c.seed_rc[out];
            }
        }
        for q in 0..16 {
            c.fr_q[2 * q] = c.fwd_q[q];
            c.fr_q[2 * q + 1] = sror(c.rev_q[q]);
        }
        c
    }
}

/// One forward roll step: `(srol(fwd) ^ fwd_q[q], sror(rev ^ rev_q[q]))`.
///
/// Both chains are one XOR shorter than the equivalent pair of
/// [`RollConsts::seed`]/[`RollConsts::srolk`] lookups.
#[inline(always)]
pub fn roll_forward(c: &RollConsts, fwd: u64, rev: u64, q: usize) -> (u64, u64) {
    debug_assert!(q < 16);
    // SAFETY-equivalent: `roll_index` can only produce 0..16, and the debug
    // assertion above guards hand-built indices.
    let fd = c.fwd_q[q & 15];
    let rd = c.rev_q[q & 15];
    (srol(fwd) ^ fd, sror(rev ^ rd))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn srol_boundaries() {
        // Case 1: all zeros -> zero
        assert_eq!(srol(0x0000_0000_0000_0000), 0x0000_0000_0000_0000);
        // Case 2: LSB only -> shifts to bit 1
        assert_eq!(srol(0x0000_0000_0000_0001), 0x0000_0000_0000_0002);
        // Case 3: bit 1 only -> shifts to bit 2
        assert_eq!(srol(0x0000_0000_0000_0002), 0x0000_0000_0000_0004);
        // Case 4: lower 32 bits all ones -> verify mask clears bit 32 and shifts rest
        assert_eq!(srol(0x0000_0000_FFFF_FFFF), 0x0000_0001_FFFF_FFFE);
        // Case 5: bit 32 only -> should wrap into LSB
        assert_eq!(srol(0x0000_0001_0000_0000), 0x0000_0000_0000_0001);
        // Case 6: bit 33 only -> normal left shift to bit 34
        assert_eq!(srol(0x0000_0002_0000_0000), 0x0000_0004_0000_0000);
        // Case 7: bits 0–62 all ones, bit 63 zero -> shifts and wrap LSB
        assert_eq!(srol(0x7FFF_FFFF_FFFF_FFFF), 0xFFFF_FFFD_FFFF_FFFF);
        // Case 8: bit 63 only -> should wrap into bit 33
        assert_eq!(srol(0x8000_0000_0000_0000), 0x0000_0002_0000_0000);
        // Case 9: bits 63, 32, 0 -> test combined wraps and shift
        assert_eq!(srol(0x8000_0001_0000_0001), 0x0000_0002_0000_0003);
        // Case 10: all ones -> rotates to all ones
        assert_eq!(srol(0xFFFF_FFFF_FFFF_FFFF), 0xFFFF_FFFF_FFFF_FFFF);
        // Case 11: random intermediate value
        assert_eq!(srol(0x0123456789ABCDEF), 0x0246_8ACD_1357_9BDF);
    }

    #[test]
    fn srol_n_boundaries() {
        // Representative x values for boundary testing:
        // Zero case:               0x0000_0000_0000_0000  // all bits zero
        // LSB only:                0x0000_0000_0000_0001  // only bit 0 set
        // Lower 32 bits all ones:  0x0000_0000_FFFF_FFFF  // bits 0–31 all set
        // Upper 32 bits all ones:  0xFFFF_FFFF_0000_0000  // bits 32–63 all set
        // Bit 32 only:             0x0000_0001_0000_0000  // 33‑bit split lower edge
        // Bit 33 only:             0x0000_0002_0000_0000  // 33‑bit split upper edge
        // MSB only:                0x8000_0000_0000_0000  // only most significant bit set
        // Random pattern:          0x0123_4567_89AB_CDEF  // generic mixed‐bit pattern
        //
        // Representative rotation distances (d):
        // d = 0:  no rotation (identity)
        // d = 1:  minimal 1‑bit rotate
        // d = 32: half‑word boundary rotate
        // d = 33: boundary + 1 rotate
        // d = 63: maximal rotate (equivalent to 1‑bit right rotate)

        // PICT-generated (x, d) → expected
        assert_eq!(srol_n(0x0000_0000_FFFF_FFFF, 1), 0x0000_0001_FFFF_FFFE);
        assert_eq!(srol_n(0x0000_0000_0000_0000, 32), 0x0000_0000_0000_0000);
        assert_eq!(srol_n(0xFFFF_FFFF_0000_0000, 32), 0xFFFF_FFFE_0000_0000);
        assert_eq!(srol_n(0x0000_0000_0000_0001, 0), 0x0000_0000_0000_0001);
        assert_eq!(srol_n(0x0000_0002_0000_0000, 33), 0x0000_0008_0000_0000);
        assert_eq!(srol_n(0x0000_0001_0000_0000, 63), 0x0000_0000_0000_0000);
        assert_eq!(srol_n(0x8000_0000_0000_0000, 63), 0x0000_0000_2000_0000);
        assert_eq!(srol_n(0x0000_0000_FFFF_FFFF, 33), 0x0000_0002_7FFF_FFFF);
        assert_eq!(srol_n(0x0123_4567_89AB_CDEF, 0), 0x0123_4567_89AB_CDEF);
        assert_eq!(srol_n(0x0000_0000_0000_0001, 1), 0x0000_0000_0000_0002);
        assert_eq!(srol_n(0x0000_0002_0000_0000, 0), 0x0000_0002_0000_0000);
        assert_eq!(srol_n(0xFFFF_FFFF_0000_0000, 33), 0xFFFF_FFFC_0000_0000);
        assert_eq!(srol_n(0xFFFF_FFFF_0000_0000, 0), 0xFFFF_FFFF_0000_0000);
        assert_eq!(srol_n(0x0000_0000_0000_0000, 0), 0x0000_0000_0000_0000);
        assert_eq!(srol_n(0x0000_0000_FFFF_FFFF, 63), 0xFFFF_FFFE_4000_0000);
        assert_eq!(srol_n(0x8000_0000_0000_0000, 32), 0x0000_0000_0000_0000);
        assert_eq!(srol_n(0x0123_4567_89AB_CDEF, 63), 0x892A_4D4C_4048_D159);
        assert_eq!(srol_n(0x0000_0000_0000_0000, 63), 0x0000_0000_0000_0000);
        assert_eq!(srol_n(0xFFFF_FFFF_0000_0000, 1), 0xFFFF_FFFE_0000_0001);
        assert_eq!(srol_n(0x0000_0000_0000_0001, 63), 0x0000_0000_4000_0000);
        assert_eq!(srol_n(0x8000_0000_0000_0000, 1), 0x0000_0002_0000_0000);
        assert_eq!(srol_n(0x0000_0002_0000_0000, 63), 0x0000_0000_0000_0000);
        assert_eq!(srol_n(0x0000_0000_0000_0001, 33), 0x0000_0000_0000_0001);
        assert_eq!(srol_n(0x0000_0000_0000_0001, 32), 0x0000_0001_0000_0000);
        assert_eq!(srol_n(0x0000_0000_0000_0000, 33), 0x0000_0000_0000_0000);
        assert_eq!(srol_n(0x0000_0001_0000_0000, 0), 0x0000_0001_0000_0000);
        assert_eq!(srol_n(0x0000_0002_0000_0000, 1), 0x0000_0004_0000_0000);
        assert_eq!(srol_n(0x0000_0001_0000_0000, 1), 0x0000_0000_0000_0001);
        assert_eq!(srol_n(0x8000_0000_0000_0000, 33), 0x0000_0000_0000_0000);
        assert_eq!(srol_n(0x0000_0001_0000_0000, 33), 0x0000_0004_0000_0000);
        assert_eq!(srol_n(0x0000_0000_FFFF_FFFF, 0), 0x0000_0000_FFFF_FFFF);
        assert_eq!(srol_n(0x0123_4567_89AB_CDEF, 1), 0x0246_8ACD_1357_9BDF);
        assert_eq!(srol_n(0x0123_4567_89AB_CDEF, 33), 0x048D_159E_09AB_CDEF);
        assert_eq!(srol_n(0x0000_0001_0000_0000, 32), 0x0000_0002_0000_0000);
        assert_eq!(srol_n(0x0123_4567_89AB_CDEF, 32), 0x0246_8ACF_44D5_E6F7);
        assert_eq!(srol_n(0x0000_0000_0000_0000, 1), 0x0000_0000_0000_0000);
        assert_eq!(srol_n(0xFFFF_FFFF_0000_0000, 63), 0x0000_0000_3FFF_FFFF);
        assert_eq!(srol_n(0x0000_0002_0000_0000, 32), 0x0000_0004_0000_0000);
        assert_eq!(srol_n(0x8000_0000_0000_0000, 0), 0x8000_0000_0000_0000);
        assert_eq!(srol_n(0x0000_0000_FFFF_FFFF, 32), 0x0000_0001_7FFF_FFFF);
    }

    #[test]
    fn sror_boundaries() {
        // Case 1: all zeros → zero
        assert_eq!(sror(0x0000_0000_0000_0000), 0x0000_0000_0000_0000);
        // Cas_eq!e 2: bit 0 only → moves into bit 32
        assert_eq!(sror(0x0000_0000_0000_0001), 0x0000_0001_0000_0000);
        // Cas_eq!e 3: bit 1 only → shifts into bit 0
        assert_eq!(sror(0x0000_0000_0000_0002), 0x0000_0000_0000_0001);
        // Cas_eq!e 4: lower 32 bits all ones
        assert_eq!(sror(0x0000_0000_FFFF_FFFF), 0x0000_0001_7FFF_FFFF);
        // Cas_eq!e 5: bit 32 only → shifts into bit 31
        assert_eq!(sror(0x0000_0001_0000_0000), 0x0000_0000_8000_0000);
        // Cas_eq!e 6: bit 33 only → wraps into bit 63
        assert_eq!(sror(0x0000_0002_0000_0000), 0x8000_0000_0000_0000);
        // Cas_eq!e 7: bits 0–62 all ones, bit 63 = 0
        assert_eq!(sror(0x7FFF_FFFF_FFFF_FFFF), 0xBFFF_FFFF_FFFF_FFFF);
        // Cas_eq!e 8: bit 63 only → shifts into bit 62
        assert_eq!(sror(0x8000_0000_0000_0000), 0x4000_0000_0000_0000);
        // Cas_eq!e 9: bits 63, 32, and 0 all set
        assert_eq!(sror(0x8000_0001_0000_0001), 0x4000_0001_8000_0000);
        // Cas_eq!e 10: all bits one → remains all ones
        assert_eq!(sror(0xFFFF_FFFF_FFFF_FFFF), 0xFFFF_FFFF_FFFF_FFFF);
        // Cas_eq!e 11: random intermediate
        // x =_eq! 0x0123_4567_89AB_CDEF → expect 0x8091_A2B3_C4D5_E6F7
        assert_eq!(sror(0x0123_4567_89AB_CDEF), 0x8091_A2B3_C4D5_E6F7);
    }

    #[test]
    fn srol_table_boundaries() {
        // Representative parameters for srol_table boundary testing:
        // Parameter c (table index):
        //   0  (N default)  – seed table for ambiguous base 'N'
        //   1  (T)          – seed table for base 'T'
        //   3  (G)          – seed table for base 'G'
        //   4  (A)          – seed table for base 'A' (4th entry in ASCII mapping)
        //   7  (C)          – seed table for base 'C'
        //
        // Parameter d (rotation count):
        //   0   (no‑op)     – no rotation, direct table lookup
        //   1   (minimal)   – single‑bit rotate
        //   30  (pre‑wrap)  – last valid index in 31‑length table (d < 31)
        //   31  (wrap‑31L)  – wraps A31L (31 % 31 == 0)
        //   32  (wrap+in‑33)– A31L wraps (32 % 31 == 1), A33R still in bounds (32 < 33)
        //   33  (both wrap) – wraps both tables (33 % 31 == 2, 33 % 33 == 0)
        //   64  (large)     – large rotation to force multi‑wrap behavior

        // PICT-generated (c, d) → expected
        assert_eq!(srol_table(0, 0), 0x0000_0000_0000_0000);
        assert_eq!(srol_table(3, 32), 0x4064_7DA0_412B_9192);
        assert_eq!(srol_table(4, 0), 0x3C8B_FBB3_95C6_0474);
        assert_eq!(srol_table(1, 0), 0x2955_49F5_4BE2_4456);
        assert_eq!(srol_table(7, 1), 0x6327_8308_C540_5699);
        assert_eq!(srol_table(1, 33), 0xA555_27D1_4BE2_4456);
        assert_eq!(srol_table(4, 33), 0xF22F_EEC9_95C6_0474);
        assert_eq!(srol_table(4, 30), 0x9E45_FDD9_32B8_C08E);
        assert_eq!(srol_table(0, 1), 0x0000_0000_0000_0000);
        assert_eq!(srol_table(0, 31), 0x0000_0000_0000_0000);
        assert_eq!(srol_table(7, 33), 0xC64F_0611_62A0_2B4C);
        assert_eq!(srol_table(1, 64), 0xA555_27D1_52F8_9115);
        assert_eq!(srol_table(1, 31), 0x2955_49F5_52F8_9115);
        assert_eq!(srol_table(0, 32), 0x0000_0000_0000_0000);
        assert_eq!(srol_table(3, 64), 0x80C8_FB40_2095_C8C9);
        assert_eq!(srol_table(7, 30), 0x18C9_E0C3_2C54_0569);
        assert_eq!(srol_table(7, 0), 0x3193_C185_62A0_2B4C);
        assert_eq!(srol_table(1, 1), 0x52AA_93E8_97C4_88AD);
        assert_eq!(srol_table(3, 0), 0x2032_3ED0_8257_2324);
        assert_eq!(srol_table(4, 64), 0xF22F_EEC8_6571_811D);
        assert_eq!(srol_table(4, 31), 0x3C8B_FBB2_6571_811D);
        assert_eq!(srol_table(3, 33), 0x80C8_FB40_8257_2324);
        assert_eq!(srol_table(0, 30), 0x0000_0000_0000_0000);
        assert_eq!(srol_table(0, 64), 0x0000_0000_0000_0000);
        assert_eq!(srol_table(7, 31), 0x3193_C184_58A8_0AD3);
        assert_eq!(srol_table(4, 1), 0x7917_F765_2B8C_08E9);
        assert_eq!(srol_table(3, 31), 0x2032_3ED0_2095_C8C9);
        assert_eq!(srol_table(1, 30), 0x14AA_A4FB_A97C_488A);
        assert_eq!(srol_table(1, 32), 0x52AA_93E8_A5F1_222B);
        assert_eq!(srol_table(7, 32), 0x6327_8308_B150_15A6);
        assert_eq!(srol_table(3, 1), 0x4064_7DA1_04AE_4648);
        assert_eq!(srol_table(7, 64), 0xC64F_0610_58A8_0AD3);
        assert_eq!(srol_table(4, 32), 0x7917_F764_CAE3_023A);
        assert_eq!(srol_table(3, 30), 0x1019_1F69_104A_E464);
        assert_eq!(srol_table(0, 33), 0x0000_0000_0000_0000);
    }

    /// `srol_const::<D>` must equal `srol_n(_, D)` for every valid `D`.
    #[test]
    fn srol_const_matches_srol_n() {
        const XS: [u64; 8] = [
            0,
            1,
            0xFFFF_FFFF_FFFF_FFFF,
            0x0000_0001_0000_0000,
            0x0000_0002_0000_0000,
            0x8000_0000_0000_0000,
            0x0123_4567_89AB_CDEF,
            0xDEAD_BEEF_CAFE_F00D,
        ];
        seq_macro_like(&XS);

        fn seq_macro_like(xs: &[u64; 8]) {
            macro_rules! check {
                ($($d:literal),*) => { $(
                    for &x in xs {
                        assert_eq!(srol_const::<$d>(x), srol_n(x, $d), "d={} x={x:#x}", $d);
                    }
                )* };
            }
            check!(
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                24, 25, 26, 27, 28, 29, 30
            );
        }
    }

    /// The folded `fwd_q`/`rev_q` tables must agree with the four separate
    /// per-base tables for every (outgoing, incoming) A/C/G/T pair.
    #[test]
    fn fused_tables_match_split_tables() {
        for k in [1u16, 2, 7, 31, 32, 33, 64, 1023, 1024] {
            let c = RollConsts::new(k);
            for &o in b"ACGTacgt" {
                for &i in b"ACGTacgtUu" {
                    let q = roll_index(o, i);
                    let (oc, ic) = (base_code(o), base_code(i));
                    assert_eq!(q, oc * 4 + ic, "k={k} o={o} i={i}");
                    assert_eq!(c.fwd_q[q], c.srolk[oc] ^ c.seed[ic], "fwd k={k}");
                    assert_eq!(c.rev_q[q], c.srolk_rc[ic] ^ c.seed_rc[oc], "rev k={k}");
                }
            }
        }
    }

    /// `roll_forward` must equal the un-fused split-table roll.
    #[test]
    fn roll_forward_matches_split_roll() {
        for k in [3u16, 21, 31] {
            let c = RollConsts::new(k);
            let (mut f, mut r) = (0xDEAD_BEEF_1234_5678u64, 0x0FED_CBA9_8765_4321u64);
            for (o, i) in b"ACGTACGTGGTTAACC".iter().zip(b"TTGGCCAAACGTACGT".iter()) {
                let want_f = srol(f) ^ c.seed[base_code(*i)] ^ c.srolk[base_code(*o)];
                let want_r = sror(r ^ c.srolk_rc[base_code(*i)] ^ c.seed_rc[base_code(*o)]);
                let (gf, gr) = roll_forward(&c, f, r, roll_index(*o, *i));
                assert_eq!((gf, gr), (want_f, want_r), "k={k}");
                f = gf;
                r = gr;
            }
        }
    }

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
}
