//! **Batch canonical hashing of many equal-length reads**, with an optional
//! SIMD backend.
//!
//! A single rolling-hash chain is inherently sequential — each hash depends on
//! the previous via a split-rotate — so SIMD cannot accelerate *one* stream.
//! The parallelism here is **across sequences**: `LANES` reads are rolled in
//! lockstep, one per SIMD lane. This is the natural shape of read processing in
//! sequencing, where inputs are millions of short, fixed-length reads.
//!
//! Because DNA has only four bases, the per-step table lookups
//! (`SEED_TAB[..]`, `srol_table(.., k)`) are resolved by **selecting among four
//! precomputed constants** rather than by a gather — cheap and fully
//! vectorisable. Reads are hashed in parallel via SIMD, detected at run time:
//! **AVX2** on `x86_64` (all four lanes in one 256-bit vector) and **NEON** on
//! `aarch64` (the four reads as two 128-bit lane-pairs). Everywhere else, and
//! for the ragged tail of reads that don't fill a lane group, a scalar path is
//! used. All paths are validated to produce **bit-identical** results.
//!
//! Like [`BlindNtHash`](crate::blind), this fast path assumes **N-free A/C/G/T**
//! input: it does not skip ambiguous bases. Use [`NtHash`](crate::kmer) when the
//! input may contain `N`.
//!
//! ```
//! use nthash_rs::batch::ReadHasher;
//! # fn main() -> nthash_rs::Result<()> {
//! let reads: &[&[u8]] = &[b"ACGTACGT", b"TTGGCCAA", b"GATTACAG", b"CCGGTTAA"];
//! let hasher = ReadHasher::new(4)?;
//! let mut count = 0;
//! hasher.hash_reads(reads, |_read_idx, _pos, _hash| count += 1)?;
//! assert_eq!(count, reads.len() * (8 - 4 + 1));
//! # Ok(()) }
//! ```

use crate::{
    kmer::{base_forward_hash, base_reverse_hash},
    tables::{base_code, srol, sror, RollConsts},
    util::canonical,
    NtHashError, Result,
};

/// Number of reads processed per SIMD lane group (AVX2 = 4×u64).
const LANES: usize = 4;

/// Batch hasher for many equal-length, N-free reads.
///
/// Construct once for a given `k`, then call [`hash_reads`](Self::hash_reads)
/// with each batch. Emits the canonical (strand-independent) hash of every
/// k-mer via a callback, allocating nothing per k-mer.
pub struct ReadHasher {
    k: u16,
    consts: RollConsts,
}

impl ReadHasher {
    /// Create a batch hasher for k-mer length `k` (`> 0`).
    pub fn new(k: u16) -> Result<Self> {
        if k == 0 {
            return Err(NtHashError::InvalidK);
        }
        Ok(Self {
            k,
            consts: RollConsts::new(k),
        })
    }

    /// Hash every k-mer of every read, invoking `f(read_idx, pos, hash)` where
    /// `hash` is the canonical hash of the k-mer at `pos` within `reads[read_idx]`.
    ///
    /// All reads must have the same length `≥ k`. Input is assumed N-free
    /// A/C/G/T (ambiguous bases are not skipped).
    ///
    /// # Errors
    /// [`NtHashError::SequenceTooShort`] if the common length is `< k`, or
    /// [`NtHashError::UnequalReadLengths`] if reads differ in length.
    pub fn hash_reads<F: FnMut(usize, usize, u64)>(&self, reads: &[&[u8]], mut f: F) -> Result<()> {
        if reads.is_empty() {
            return Ok(());
        }
        let len = reads[0].len();
        let k = self.k as usize;
        if len < k {
            return Err(NtHashError::SequenceTooShort {
                seq_len: len,
                k: self.k,
            });
        }
        for r in reads {
            if r.len() != len {
                return Err(NtHashError::UnequalReadLengths {
                    expected: len,
                    found: r.len(),
                });
            }
        }

        // Process full lane groups with the fastest available kernel, then the
        // ragged tail scalar.
        let mut base = 0;
        while base + LANES <= reads.len() {
            let group: [&[u8]; LANES] = [
                reads[base],
                reads[base + 1],
                reads[base + 2],
                reads[base + 3],
            ];
            self.hash_group(base, &group, &mut f);
            base += LANES;
        }
        for (offset, r) in reads[base..].iter().enumerate() {
            self.hash_one(base + offset, r, &mut f);
        }
        Ok(())
    }

    /// Scalar single-read roll (also the reference for the SIMD kernel).
    fn hash_one<F: FnMut(usize, usize, u64)>(&self, ri: usize, read: &[u8], f: &mut F) {
        let k = self.k as usize;
        let c = &self.consts;
        let mut fwd = base_forward_hash(read, self.k);
        let mut rev = base_reverse_hash(read, self.k);
        f(ri, 0, canonical(fwd, rev));
        for i in k..read.len() {
            let out = base_code(read[i - k]);
            let inc = base_code(read[i]);
            fwd = srol(fwd) ^ c.seed[inc] ^ c.srolk[out];
            rev = sror(rev ^ c.srolk_rc[inc] ^ c.seed_rc[out]);
            f(ri, i - k + 1, canonical(fwd, rev));
        }
    }

    /// Hash a full group of `LANES` reads, dispatching to SIMD when available.
    fn hash_group<F: FnMut(usize, usize, u64)>(
        &self,
        base: usize,
        group: &[&[u8]; LANES],
        f: &mut F,
    ) {
        #[cfg(target_arch = "x86_64")]
        {
            if std::arch::is_x86_feature_detected!("avx2") {
                // SAFETY: guarded by run-time AVX2 detection.
                unsafe { self.hash_group_avx2(base, group, f) };
                return;
            }
        }
        #[cfg(target_arch = "aarch64")]
        {
            if std::arch::is_aarch64_feature_detected!("neon") {
                // SAFETY: guarded by run-time NEON detection (always true on the
                // aarch64 baseline, checked for form).
                unsafe { self.hash_group_neon(base, group, f) };
                return;
            }
        }
        for (lane, read) in group.iter().enumerate() {
            self.hash_one(base + lane, read, f);
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    unsafe fn hash_group_avx2<F: FnMut(usize, usize, u64)>(
        &self,
        base: usize,
        group: &[&[u8]; LANES],
        f: &mut F,
    ) {
        use std::arch::x86_64::*;

        let k = self.k as usize;
        let len = group[0].len();
        let c = &self.consts;

        // Splat each 4-constant table into per-code AVX2 vectors.
        let splat = |t: &[u64; 4]| {
            [
                _mm256_set1_epi64x(t[0] as i64),
                _mm256_set1_epi64x(t[1] as i64),
                _mm256_set1_epi64x(t[2] as i64),
                _mm256_set1_epi64x(t[3] as i64),
            ]
        };
        let seed = splat(&c.seed);
        let srolk = splat(&c.srolk);
        let srolk_rc = splat(&c.srolk_rc);
        let seed_rc = splat(&c.seed_rc);

        // Select `t[code_lane]` per lane via two blends (2-bit code).
        #[inline(always)]
        unsafe fn select4(t: &[__m256i; 4], code: __m256i) -> __m256i {
            let one = _mm256_set1_epi64x(1);
            let m0 = _mm256_cmpeq_epi64(_mm256_and_si256(code, one), one);
            let lo = _mm256_blendv_epi8(t[0], t[1], m0);
            let hi = _mm256_blendv_epi8(t[2], t[3], m0);
            let m1 = _mm256_cmpeq_epi64(_mm256_and_si256(_mm256_srli_epi64(code, 1), one), one);
            _mm256_blendv_epi8(lo, hi, m1)
        }

        #[inline(always)]
        unsafe fn srol_avx2(x: __m256i) -> __m256i {
            let shifted = _mm256_and_si256(
                _mm256_slli_epi64(x, 1),
                _mm256_set1_epi64x(0xFFFF_FFFD_FFFF_FFFFu64 as i64),
            );
            let hi = _mm256_srli_epi64(
                _mm256_and_si256(x, _mm256_set1_epi64x(0x8000_0000_0000_0000u64 as i64)),
                30,
            );
            let lo = _mm256_srli_epi64(
                _mm256_and_si256(x, _mm256_set1_epi64x(0x0000_0001_0000_0000u64 as i64)),
                32,
            );
            _mm256_or_si256(shifted, _mm256_or_si256(hi, lo))
        }

        #[inline(always)]
        unsafe fn sror_avx2(x: __m256i) -> __m256i {
            let shifted = _mm256_and_si256(
                _mm256_srli_epi64(x, 1),
                _mm256_set1_epi64x(0xFFFF_FFFE_FFFF_FFFFu64 as i64),
            );
            let hi = _mm256_slli_epi64(
                _mm256_and_si256(x, _mm256_set1_epi64x(0x0000_0002_0000_0000u64 as i64)),
                30,
            );
            let lo = _mm256_slli_epi64(
                _mm256_and_si256(x, _mm256_set1_epi64x(0x0000_0000_0000_0001u64 as i64)),
                32,
            );
            _mm256_or_si256(shifted, _mm256_or_si256(hi, lo))
        }

        // Per-lane 2-bit codes for one column, as an AVX2 vector.
        #[inline(always)]
        unsafe fn codes_at(group: &[&[u8]; LANES], i: usize) -> __m256i {
            _mm256_set_epi64x(
                base_code(group[3][i]) as i64,
                base_code(group[2][i]) as i64,
                base_code(group[1][i]) as i64,
                base_code(group[0][i]) as i64,
            )
        }

        // Seed each lane with the first k-mer's forward/reverse hash.
        let mut fwd_arr = [0u64; LANES];
        let mut rev_arr = [0u64; LANES];
        for lane in 0..LANES {
            fwd_arr[lane] = base_forward_hash(group[lane], self.k);
            rev_arr[lane] = base_reverse_hash(group[lane], self.k);
            f(base + lane, 0, canonical(fwd_arr[lane], rev_arr[lane]));
        }
        let mut fwd = _mm256_loadu_si256(fwd_arr.as_ptr() as *const __m256i);
        let mut rev = _mm256_loadu_si256(rev_arr.as_ptr() as *const __m256i);

        for i in k..len {
            let out_code = codes_at(group, i - k);
            let in_code = codes_at(group, i);

            // fwd = srol(fwd) ^ seed[in] ^ srolk[out]
            fwd = _mm256_xor_si256(
                _mm256_xor_si256(srol_avx2(fwd), select4(&seed, in_code)),
                select4(&srolk, out_code),
            );
            // rev = sror(rev ^ srolk_rc[in] ^ seed_rc[out])
            let inner = _mm256_xor_si256(
                _mm256_xor_si256(rev, select4(&srolk_rc, in_code)),
                select4(&seed_rc, out_code),
            );
            rev = sror_avx2(inner);

            // canonical = fwd + rev (wrapping), per lane.
            let canon = _mm256_add_epi64(fwd, rev);
            _mm256_storeu_si256(fwd_arr.as_mut_ptr() as *mut __m256i, canon);
            let pos = i - k + 1;
            f(base, pos, fwd_arr[0]);
            f(base + 1, pos, fwd_arr[1]);
            f(base + 2, pos, fwd_arr[2]);
            f(base + 3, pos, fwd_arr[3]);
        }
    }

    #[cfg(target_arch = "aarch64")]
    #[target_feature(enable = "neon")]
    unsafe fn hash_group_neon<F: FnMut(usize, usize, u64)>(
        &self,
        base: usize,
        group: &[&[u8]; LANES],
        f: &mut F,
    ) {
        use core::arch::aarch64::*;

        let k = self.k as usize;
        let len = group[0].len();
        let c = &self.consts;

        // Splat each 4-constant table across both lanes of a uint64x2_t. NEON is
        // 128-bit (2×u64), so the 4-read group is processed as two lane-pairs.
        let splat = |t: &[u64; 4]| {
            [
                vdupq_n_u64(t[0]),
                vdupq_n_u64(t[1]),
                vdupq_n_u64(t[2]),
                vdupq_n_u64(t[3]),
            ]
        };
        let seed = splat(&c.seed);
        let srolk = splat(&c.srolk);
        let srolk_rc = splat(&c.srolk_rc);
        let seed_rc = splat(&c.seed_rc);

        // Select `t[code_lane]` per lane via two bit-selects (2-bit code).
        #[inline(always)]
        unsafe fn select4(t: &[uint64x2_t; 4], code: uint64x2_t) -> uint64x2_t {
            let one = vdupq_n_u64(1);
            let m0 = vceqq_u64(vandq_u64(code, one), one);
            let lo = vbslq_u64(m0, t[1], t[0]);
            let hi = vbslq_u64(m0, t[3], t[2]);
            let m1 = vceqq_u64(vandq_u64(vshrq_n_u64::<1>(code), one), one);
            vbslq_u64(m1, hi, lo)
        }

        #[inline(always)]
        unsafe fn srol_neon(x: uint64x2_t) -> uint64x2_t {
            let shifted = vandq_u64(vshlq_n_u64::<1>(x), vdupq_n_u64(0xFFFF_FFFD_FFFF_FFFF));
            let hi = vshrq_n_u64::<30>(vandq_u64(x, vdupq_n_u64(0x8000_0000_0000_0000)));
            let lo = vshrq_n_u64::<32>(vandq_u64(x, vdupq_n_u64(0x0000_0001_0000_0000)));
            vorrq_u64(shifted, vorrq_u64(hi, lo))
        }

        #[inline(always)]
        unsafe fn sror_neon(x: uint64x2_t) -> uint64x2_t {
            let shifted = vandq_u64(vshrq_n_u64::<1>(x), vdupq_n_u64(0xFFFF_FFFE_FFFF_FFFF));
            let hi = vshlq_n_u64::<30>(vandq_u64(x, vdupq_n_u64(0x0000_0002_0000_0000)));
            let lo = vshlq_n_u64::<32>(vandq_u64(x, vdupq_n_u64(0x0000_0000_0000_0001)));
            vorrq_u64(shifted, vorrq_u64(hi, lo))
        }

        // Roll one 2-lane vector pair forward by one base.
        #[inline(always)]
        #[allow(clippy::too_many_arguments)]
        unsafe fn roll_pair(
            fwd: uint64x2_t,
            rev: uint64x2_t,
            in_code: uint64x2_t,
            out_code: uint64x2_t,
            seed: &[uint64x2_t; 4],
            srolk: &[uint64x2_t; 4],
            srolk_rc: &[uint64x2_t; 4],
            seed_rc: &[uint64x2_t; 4],
        ) -> (uint64x2_t, uint64x2_t) {
            let fwd = veorq_u64(
                veorq_u64(srol_neon(fwd), select4(seed, in_code)),
                select4(srolk, out_code),
            );
            let inner = veorq_u64(
                veorq_u64(rev, select4(srolk_rc, in_code)),
                select4(seed_rc, out_code),
            );
            (fwd, sror_neon(inner))
        }

        // Per-lane 2-bit codes for two adjacent reads at column `i`.
        #[inline(always)]
        unsafe fn codes_pair(group: &[&[u8]; LANES], lane0: usize, i: usize) -> uint64x2_t {
            let arr = [
                base_code(group[lane0][i]) as u64,
                base_code(group[lane0 + 1][i]) as u64,
            ];
            vld1q_u64(arr.as_ptr())
        }

        // Seed each lane with the first k-mer's forward/reverse hash.
        let mut fwd_arr = [0u64; LANES];
        let mut rev_arr = [0u64; LANES];
        for lane in 0..LANES {
            fwd_arr[lane] = base_forward_hash(group[lane], self.k);
            rev_arr[lane] = base_reverse_hash(group[lane], self.k);
            f(base + lane, 0, canonical(fwd_arr[lane], rev_arr[lane]));
        }
        let mut fwd0 = vld1q_u64(fwd_arr.as_ptr());
        let mut fwd1 = vld1q_u64(fwd_arr.as_ptr().add(2));
        let mut rev0 = vld1q_u64(rev_arr.as_ptr());
        let mut rev1 = vld1q_u64(rev_arr.as_ptr().add(2));

        for i in k..len {
            let out0 = codes_pair(group, 0, i - k);
            let in0 = codes_pair(group, 0, i);
            let out1 = codes_pair(group, 2, i - k);
            let in1 = codes_pair(group, 2, i);

            let (nf0, nr0) = roll_pair(fwd0, rev0, in0, out0, &seed, &srolk, &srolk_rc, &seed_rc);
            let (nf1, nr1) = roll_pair(fwd1, rev1, in1, out1, &seed, &srolk, &srolk_rc, &seed_rc);
            fwd0 = nf0;
            rev0 = nr0;
            fwd1 = nf1;
            rev1 = nr1;

            // canonical = fwd + rev (wrapping), per lane.
            vst1q_u64(fwd_arr.as_mut_ptr(), vaddq_u64(fwd0, rev0));
            vst1q_u64(fwd_arr.as_mut_ptr().add(2), vaddq_u64(fwd1, rev1));
            let pos = i - k + 1;
            f(base, pos, fwd_arr[0]);
            f(base + 1, pos, fwd_arr[1]);
            f(base + 2, pos, fwd_arr[2]);
            f(base + 3, pos, fwd_arr[3]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BlindNtHashBuilder;
    use std::collections::BTreeMap;

    /// Reference: per-read canonical hashes via `BlindNtHash`.
    fn reference(reads: &[&[u8]], k: u16) -> BTreeMap<(usize, usize), u64> {
        let mut m = BTreeMap::new();
        for (ri, r) in reads.iter().enumerate() {
            BlindNtHashBuilder::new(r)
                .k(k)
                .num_hashes(1)
                .finish()
                .unwrap()
                .for_each_hash(|pos, h| {
                    m.insert((ri, pos), h[0]);
                });
        }
        m
    }

    fn batch(reads: &[&[u8]], k: u16) -> BTreeMap<(usize, usize), u64> {
        let mut m = BTreeMap::new();
        ReadHasher::new(k)
            .unwrap()
            .hash_reads(reads, |ri, pos, h| {
                m.insert((ri, pos), h);
            })
            .unwrap();
        m
    }

    #[test]
    fn matches_blind_across_lane_group_and_tail() {
        // 6 reads = one full AVX2 group of 4 + a 2-read scalar tail.
        let reads: &[&[u8]] = &[
            b"ACGTACGTACGT",
            b"TTGGCCAATTGG",
            b"GATTACAGATTA",
            b"CCGGTTAACCGG",
            b"AAAACCCCGGGG",
            b"TACGTACGTACG",
        ];
        for k in [1u16, 2, 4, 7, 12] {
            assert_eq!(batch(reads, k), reference(reads, k), "k={k}");
        }
    }

    #[test]
    fn unequal_lengths_error() {
        let reads: &[&[u8]] = &[b"ACGT", b"ACG"];
        assert_eq!(
            ReadHasher::new(2).unwrap().hash_reads(reads, |_, _, _| {}),
            Err(NtHashError::UnequalReadLengths {
                expected: 4,
                found: 3
            })
        );
    }

    #[test]
    fn too_short_errors() {
        let reads: &[&[u8]] = &[b"AC", b"GT"];
        assert!(matches!(
            ReadHasher::new(4).unwrap().hash_reads(reads, |_, _, _| {}),
            Err(NtHashError::SequenceTooShort { seq_len: 2, k: 4 })
        ));
    }

    #[test]
    fn empty_batch_ok() {
        assert_eq!(
            ReadHasher::new(4).unwrap().hash_reads(&[], |_, _, _| {}),
            Ok(())
        );
    }
}
