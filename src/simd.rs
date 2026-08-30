//! Shared **x86-64 vector primitives** for the lane-parallel rolling kernels.
//!
//! A single ntHash chain is serial, so every SIMD kernel in this crate gets its
//! parallelism the same way: `L` *independent* hash chains are rolled in
//! lockstep, one per lane. What differs is where the chains come from —
//! [`batch`](crate::batch) puts one read in each lane, [`wide`](crate::wide)
//! splits one long sequence into `L` disjoint ranges — so the per-step
//! arithmetic is identical and lives here.
//!
//! Two things make the AVX-512 kernel much faster than a blend-based AVX2 one:
//!
//! - **`vpermi2q` resolves the whole 16-entry roll table in one instruction.**
//!   Folding the four `RollConsts` tables into
//!   [`fwd_q`](crate::tables::RollConsts::fwd_q) / `rev_q` gives exactly 16
//!   `u64`s — two `zmm` registers — so a per-lane 4-bit index is a single
//!   permute. Emulating that with AVX2 blends needs 15 of them and is *slower
//!   than scalar*; AVX2 therefore uses a gather instead.
//! - **`vprolq` + `vpternlogq` cut `srol`/`sror` to five instructions**, versus
//!   eight for the shift/mask formulation (`srol` is a plain rotate with bits 0
//!   and 33 swapped, and the three-input XOR that repairs them is one op).

#![cfg(target_arch = "x86_64")]

use std::arch::x86_64::*;

// --------------------------------------------------------------------------
// AVX-512: eight lanes
// --------------------------------------------------------------------------

/// Lanes processed per step by the AVX-512 kernels.
pub const LANES_X8: usize = 8;

/// Eight-lane [`srol`](crate::tables::srol).
///
/// `srol(x)` is `x.rotate_left(1)` with bits 0 and 33 exchanged, so one
/// `vprolq` plus a two-`vpternlogq` bit swap suffices.
///
/// # Safety
/// Requires AVX-512F.
#[inline]
#[target_feature(enable = "avx512f")]
pub unsafe fn srol_x8(x: __m512i) -> __m512i {
    let v = _mm512_rol_epi64::<1>(x);
    // y = (v ^ (v >> 33)) & 1 — the XOR of the two bits that must swap.
    let y = _mm512_ternarylogic_epi64::<0x28>(v, _mm512_srli_epi64::<33>(v), _mm512_set1_epi64(1));
    _mm512_ternarylogic_epi64::<0x96>(v, y, _mm512_slli_epi64::<33>(y))
}

/// Eight-lane [`sror`](crate::tables::sror): `rotate_right(1)` with bits 32 and
/// 63 exchanged.
///
/// # Safety
/// Requires AVX-512F.
#[inline]
#[target_feature(enable = "avx512f")]
pub unsafe fn sror_x8(x: __m512i) -> __m512i {
    let v = _mm512_ror_epi64::<1>(x);
    let y = _mm512_ternarylogic_epi64::<0x28>(
        v,
        _mm512_srli_epi64::<31>(v),
        _mm512_set1_epi64(1i64 << 32),
    );
    _mm512_ternarylogic_epi64::<0x96>(v, y, _mm512_slli_epi64::<31>(y))
}

/// Per-lane [`roll_index`](crate::tables::roll_index) from the two byte streams.
///
/// `ow` / `iw` each hold, per lane, the next eight bases of that lane's
/// outgoing / incoming stream packed little-endian into a `u64`; only the low
/// byte is consulted, so a caller advances by shifting both right by 8.
///
/// # Safety
/// Requires AVX-512F.
#[inline]
#[target_feature(enable = "avx512f")]
pub unsafe fn roll_index_x8(ow: __m512i, iw: __m512i) -> __m512i {
    // ((out & 6) << 1) | ((in >> 1) & 3)
    let out4 = _mm512_slli_epi64::<1>(_mm512_and_si512(ow, _mm512_set1_epi64(6)));
    let inc = _mm512_srli_epi64::<1>(iw);
    _mm512_ternarylogic_epi64::<0xF8>(out4, inc, _mm512_set1_epi64(3))
}

/// One eight-lane forward roll against the folded roll tables.
///
/// `fq` / `rq` are the two halves of `RollConsts::fwd_q` and `rev_q` loaded as
/// `[lo, hi]`; `q` is the per-lane index from [`roll_index_x8`].
///
/// # Safety
/// Requires AVX-512F.
#[inline]
#[target_feature(enable = "avx512f")]
pub unsafe fn roll_x8(
    fwd: __m512i,
    rev: __m512i,
    q: __m512i,
    fq: [__m512i; 2],
    rq: [__m512i; 2],
) -> (__m512i, __m512i) {
    let f = _mm512_xor_si512(srol_x8(fwd), _mm512_permutex2var_epi64(fq[0], q, fq[1]));
    let r = sror_x8(_mm512_xor_si512(
        rev,
        _mm512_permutex2var_epi64(rq[0], q, rq[1]),
    ));
    (f, r)
}

/// Load a 16-`u64` roll table as the `[lo, hi]` pair `vpermi2q` selects from.
///
/// # Safety
/// Requires AVX-512F; `t` must be 16 elements.
#[inline]
#[target_feature(enable = "avx512f")]
pub unsafe fn load_table_x8(t: &[u64; 16]) -> [__m512i; 2] {
    [
        _mm512_loadu_si512(t.as_ptr() as *const __m512i),
        _mm512_loadu_si512(t.as_ptr().add(8) as *const __m512i),
    ]
}

/// Eight-lane split-rotate left by a **runtime** distance, applied field-wise.
///
/// `srol` rotates the 33-bit and 31-bit fields independently, so `srol^d` is
/// the low field rotated by `d % 33` and the high field by `d % 31`. Used once
/// per group to finish the reverse seed hash, so op count matters less than
/// handling any `d`.
///
/// # Safety
/// Requires AVX-512F.
#[inline]
#[target_feature(enable = "avx512f")]
pub unsafe fn srol_field_x8(x: __m512i, d: usize) -> __m512i {
    #[inline(always)]
    unsafe fn rot(v: __m512i, w: u32, by: u32, mask: __m512i) -> __m512i {
        if by == 0 {
            return v;
        }
        let l = _mm512_sll_epi64(v, _mm_cvtsi32_si128(by as i32));
        let r = _mm512_srl_epi64(v, _mm_cvtsi32_si128((w - by) as i32));
        _mm512_and_si512(_mm512_or_si512(l, r), mask)
    }
    let lo_mask = _mm512_set1_epi64(0x1_FFFF_FFFF);
    let hi_mask = _mm512_set1_epi64(0x7FFF_FFFF);
    let lo = rot(_mm512_and_si512(x, lo_mask), 33, (d % 33) as u32, lo_mask);
    let hi = rot(_mm512_srli_epi64::<33>(x), 31, (d % 31) as u32, hi_mask);
    _mm512_or_si512(lo, _mm512_slli_epi64::<33>(hi))
}

/// Eight-lane first-k-mer seed hashes: `(forward, reverse)` per lane.
///
/// The scalar [`base_forward_hash`](crate::kmer::base_forward_hash) /
/// [`base_reverse_hash`](crate::kmer::base_reverse_hash) pair is a serial
/// `O(k)` chain per sequence, and a lane group needs sixteen of them — ~28% of
/// batch hashing for 150 bp reads, ~39% for 100 bp. Here the eight lanes share
/// one chain.
///
/// The reverse hash is `XOR_j srol^j(seed_rc[c_j])`, which reads the bases
/// back to front. Rather than walk backwards, accumulate
/// `g = XOR_j sror^(k-1-j)(seed_rc[c_j])` in the same forward pass and finish
/// with a single `srol^(k-1)` — `srol^(k-1) ∘ sror^(k-1-j) = srol^j`.
///
/// `seed16` / `seed_rc16` are `RollConsts::seed` / `seed_rc` padded to sixteen
/// entries (only `0..4` are ever selected).
///
/// # Safety
/// Requires AVX-512F. `words(off)` must return, per lane, the eight bases at
/// byte offset `off`; every lane needs at least `max(k, 8)` bases, and `off + 8`
/// never exceeds `len`.
#[inline]
#[target_feature(enable = "avx512f")]
pub unsafe fn seed_x8<W: Fn(usize) -> __m512i>(
    k: usize,
    len: usize,
    words: W,
    seed16: &[u64; 16],
    seed_rc16: &[u64; 16],
) -> (__m512i, __m512i) {
    debug_assert!(len >= 8 && len >= k && k >= 1);
    let sd = load_table_x8(seed16);
    let rc = load_table_x8(seed_rc16);
    let three = _mm512_set1_epi64(3);

    let mut f = _mm512_setzero_si512();
    let mut g = _mm512_setzero_si512();
    let mut j = 0usize;
    while j < k {
        // Words are eight bases wide; near the end of a short lane, back the
        // load off to `len - 8` and shift the already-consumed bases out.
        let off = if j + 8 <= len { j } else { len - 8 };
        let mut w = words(off);
        let skip = j - off;
        if skip != 0 {
            w = _mm512_srl_epi64(w, _mm_cvtsi32_si128((skip * 8) as i32));
        }
        let n = core::cmp::min(8 - skip, k - j);
        for _ in 0..n {
            let c = _mm512_and_si512(_mm512_srli_epi64::<1>(w), three);
            f = _mm512_xor_si512(srol_x8(f), _mm512_permutex2var_epi64(sd[0], c, sd[1]));
            g = _mm512_xor_si512(sror_x8(g), _mm512_permutex2var_epi64(rc[0], c, rc[1]));
            w = _mm512_srli_epi64::<8>(w);
        }
        j += n;
    }
    (f, srol_field_x8(g, k - 1))
}

// --------------------------------------------------------------------------
// 128-bit: both chains of one stream
// --------------------------------------------------------------------------

/// Rotation amounts and masks for [`roll_x2`], built once per traversal.
pub struct Roll2 {
    rot1: __m128i,
    m1: __m128i,
    rot2: __m128i,
    m2: __m128i,
    m3: __m128i,
}

/// Build the [`Roll2`] constants. SSE2 only, so no feature gate is needed.
#[inline]
pub fn roll2() -> Roll2 {
    // SAFETY: SSE2 is part of the x86-64 baseline.
    unsafe {
        Roll2 {
            // lane 0 rotates left by 1 (`srol`), lane 1 right by 1 (`sror`).
            rot1: _mm_set_epi64x(63, 1),
            m1: _mm_set_epi64x(
                0x7FFF_FFFE_FFFF_FFFFu64 as i64,
                0xFFFF_FFFD_FFFF_FFFEu64 as i64,
            ),
            // The field that wraps: bit 63 -> 33 forward, bit 33 -> 63 reverse.
            rot2: _mm_set_epi64x(30, 34),
            m2: _mm_set_epi64x(1i64 << 63, 1i64 << 33),
            // The other wrap: bit 32 -> 0 forward, bit 0 -> 32 reverse.
            m3: _mm_set_epi64x(1i64 << 32, 1),
        }
    }
}

/// Roll **both** chains of a single stream, held as lane 0 = forward,
/// lane 1 = reverse, against the interleaved delta pair `d`
/// (`RollConsts::fr_q[2q..2q+2]`).
///
/// `srol` and `sror` run in opposite directions, which is why one register
/// cannot hold them under a plain shift — but each is three *masked rotations*
/// of the input, and `vprolvq` takes a per-lane rotation amount, so lane 0 can
/// rotate left while lane 1 rotates right in the same instruction. The three
/// masked terms occupy disjoint bits, so they combine with XOR, and the delta
/// folds into the same tree: two `vpternlogq`s finish it.
///
/// The result is seven ops and a three-deep chain for both chains, against
/// sixteen ops and five deep for the scalar pair — measured 1.4x on the
/// canonical scan, bit-identical.
///
/// # Safety
/// Requires AVX-512F and AVX-512VL.
#[inline]
#[target_feature(enable = "avx512f,avx512vl")]
pub unsafe fn roll_x2(x: __m128i, d: __m128i, c: &Roll2) -> __m128i {
    let a = _mm_and_si128(_mm_rolv_epi64(x, c.rot1), c.m1);
    let b = _mm_and_si128(_mm_rolv_epi64(x, c.rot2), c.m2);
    // (rot32(x) & m3) ^ d
    let e = _mm_ternarylogic_epi64::<0x6A>(_mm_rol_epi64::<32>(x), c.m3, d);
    // a ^ b ^ e
    _mm_ternarylogic_epi64::<0x96>(a, b, e)
}

/// Pack a `(forward, reverse)` pair into the [`roll_x2`] lane layout.
#[inline]
pub fn pack_x2(fwd: u64, rev: u64) -> __m128i {
    // SAFETY: SSE2 baseline.
    unsafe { _mm_set_epi64x(rev as i64, fwd as i64) }
}

/// The canonical hash of a [`roll_x2`] state: lane 0 plus lane 1, wrapping.
#[inline]
pub fn canonical_x2(x: __m128i) -> u64 {
    // SAFETY: SSE2 baseline.
    unsafe { _mm_cvtsi128_si64(_mm_add_epi64(x, _mm_unpackhi_epi64(x, x))) as u64 }
}

/// Unpack a [`roll_x2`] state back to `(forward, reverse)`.
#[inline]
pub fn unpack_x2(x: __m128i) -> (u64, u64) {
    // SAFETY: SSE2 baseline.
    unsafe {
        (
            _mm_cvtsi128_si64(x) as u64,
            _mm_cvtsi128_si64(_mm_unpackhi_epi64(x, x)) as u64,
        )
    }
}

// --------------------------------------------------------------------------
// AVX2: four lanes
// --------------------------------------------------------------------------

/// Lanes per 256-bit register. The AVX2 kernels drive two of these at once, so
/// they cover the same eight chains as the AVX-512 kernels.
pub const LANES_X4: usize = 4;

/// Four-lane [`srol`](crate::tables::srol).
///
/// AVX2 has no variable rotate, so this is the shift/mask formulation; the
/// masks are hoisted by the caller's loop, not rebuilt per step.
///
/// # Safety
/// Requires AVX2.
#[inline]
#[target_feature(enable = "avx2")]
pub unsafe fn srol_x4(x: __m256i) -> __m256i {
    let shifted = _mm256_and_si256(
        _mm256_slli_epi64::<1>(x),
        _mm256_set1_epi64x(0xFFFF_FFFD_FFFF_FFFFu64 as i64),
    );
    let hi = _mm256_srli_epi64::<30>(_mm256_and_si256(
        x,
        _mm256_set1_epi64x(0x8000_0000_0000_0000u64 as i64),
    ));
    let lo = _mm256_srli_epi64::<32>(_mm256_and_si256(x, _mm256_set1_epi64x(1i64 << 32)));
    _mm256_or_si256(shifted, _mm256_or_si256(hi, lo))
}

/// Four-lane [`sror`](crate::tables::sror).
///
/// # Safety
/// Requires AVX2.
#[inline]
#[target_feature(enable = "avx2")]
pub unsafe fn sror_x4(x: __m256i) -> __m256i {
    let shifted = _mm256_and_si256(
        _mm256_srli_epi64::<1>(x),
        _mm256_set1_epi64x(0xFFFF_FFFE_FFFF_FFFFu64 as i64),
    );
    let hi = _mm256_slli_epi64::<30>(_mm256_and_si256(x, _mm256_set1_epi64x(2i64 << 32)));
    let lo = _mm256_slli_epi64::<32>(_mm256_and_si256(x, _mm256_set1_epi64x(1)));
    _mm256_or_si256(shifted, _mm256_or_si256(hi, lo))
}

/// Four-lane [`roll_index_x8`].
///
/// # Safety
/// Requires AVX2.
#[inline]
#[target_feature(enable = "avx2")]
pub unsafe fn roll_index_x4(ow: __m256i, iw: __m256i) -> __m256i {
    let out4 = _mm256_slli_epi64::<1>(_mm256_and_si256(ow, _mm256_set1_epi64x(6)));
    let inc = _mm256_and_si256(_mm256_srli_epi64::<1>(iw), _mm256_set1_epi64x(3));
    _mm256_or_si256(out4, inc)
}

/// One four-lane forward roll. The 16-entry tables are read with a gather:
/// they are 128 bytes and stay L1-resident, and a blend tree over 16 entries
/// costs 15 `vpblendvb` — measurably slower than the scalar roll.
///
/// # Safety
/// Requires AVX2; `fq` and `rq` must each point at 16 `u64`s.
#[inline]
#[target_feature(enable = "avx2")]
pub unsafe fn roll_x4(
    fwd: __m256i,
    rev: __m256i,
    q: __m256i,
    fq: *const i64,
    rq: *const i64,
) -> (__m256i, __m256i) {
    let f = _mm256_xor_si256(srol_x4(fwd), _mm256_i64gather_epi64::<8>(fq, q));
    let r = sror_x4(_mm256_xor_si256(rev, _mm256_i64gather_epi64::<8>(rq, q)));
    (f, r)
}
