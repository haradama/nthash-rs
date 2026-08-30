//! **Batch canonical hashing of many equal-length reads**, with an optional
//! SIMD backend.
//!
//! A single rolling-hash chain is inherently sequential — each hash depends on
//! the previous via a split-rotate — so SIMD cannot accelerate *one* stream.
//! The parallelism here is **across sequences**: `LANES` reads are rolled in
//! lockstep, one per SIMD lane. This is the natural shape of read processing in
//! sequencing, where inputs are millions of short, fixed-length reads.
//!
//! Because DNA has only four bases, one roll step is fully described by the
//! `(outgoing, incoming)` pair — 16 possibilities — so the per-step table
//! lookup is a 16-entry index into [`RollConsts::fwd_q`] / `rev_q` rather than
//! a general gather. Backends, detected at run time:
//!
//! | backend | lanes | roll-table lookup |
//! |---|---|---|
//! | AVX-512F | 8 reads | `vpermi2q`, one instruction |
//! | AVX2 | 8 reads (2×256-bit) | 128-byte L1 gather |
//! | NEON (`aarch64`) | 4 reads (2×128-bit) | 2-bit blend tree |
//! | scalar | 1 read | direct index |
//!
//! Each lane's next eight bases are fetched as one `u64` and consumed by
//! shifting, so byte handling costs one load per lane per eight steps. All
//! paths are validated to produce **bit-identical** results.
//!
//! The callback order is **unspecified** — vector backends emit a whole lane
//! group per step, the scalar tail emits a whole read at a time.
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
    tables::{roll_forward, roll_index, RollConsts},
    util::canonical,
    NtHashError, Result,
};

#[cfg(target_arch = "aarch64")]
use crate::tables::base_code;

/// Number of reads per NEON lane group (two 128-bit lane-pairs).
#[cfg(target_arch = "aarch64")]
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

        // Dispatch once for the whole batch, not per lane group.
        #[cfg(target_arch = "x86_64")]
        {
            if std::arch::is_x86_feature_detected!("avx512f") {
                // SAFETY: guarded by run-time AVX-512F detection.
                unsafe { self.hash_reads_avx512(reads, &mut f) };
                return Ok(());
            }
            if std::arch::is_x86_feature_detected!("avx2") {
                // SAFETY: guarded by run-time AVX2 detection.
                unsafe { self.hash_reads_avx2(reads, &mut f) };
                return Ok(());
            }
        }
        #[cfg(target_arch = "aarch64")]
        {
            if std::arch::is_aarch64_feature_detected!("neon") {
                // SAFETY: guarded by run-time NEON detection (always true on the
                // aarch64 baseline, checked for form).
                unsafe { self.hash_reads_neon(reads, &mut f) };
                return Ok(());
            }
        }
        #[allow(unreachable_code)]
        {
            for (ri, r) in reads.iter().enumerate() {
                self.hash_one(ri, r, &mut f);
            }
            Ok(())
        }
    }

    /// Scalar single-read roll (also the reference for the SIMD kernels).
    fn hash_one<F: FnMut(usize, usize, u64)>(&self, ri: usize, read: &[u8], f: &mut F) {
        let k = self.k as usize;
        let fwd = base_forward_hash(read, self.k);
        let rev = base_reverse_hash(read, self.k);
        f(ri, 0, canonical(fwd, rev));
        self.roll_tail(ri, read, k, fwd, rev, f);
    }

    /// Continue a scalar roll from `(fwd, rev)`, whose k-mer ends just before
    /// `from`, to the end of the read. Shared by [`hash_one`](Self::hash_one)
    /// and every SIMD kernel's ragged tail.
    #[inline]
    fn roll_tail<F: FnMut(usize, usize, u64)>(
        &self,
        ri: usize,
        read: &[u8],
        from: usize,
        mut fwd: u64,
        mut rev: u64,
        f: &mut F,
    ) {
        let k = self.k as usize;
        let c = &self.consts;
        for i in from..read.len() {
            let (nf, nr) = roll_forward(c, fwd, rev, roll_index(read[i - k], read[i]));
            fwd = nf;
            rev = nr;
            f(ri, i - k + 1, canonical(fwd, rev));
        }
    }

    /// AVX-512F: roll eight reads per group.
    ///
    /// Each step resolves both 16-entry roll tables with one `vpermi2q` and
    /// rolls all eight lanes; the eight hashes of a step are buffered and
    /// handed to `f` a block at a time, so the callback never sits between the
    /// vector store and the scalar reads of the same value.
    ///
    /// # Safety
    /// Requires AVX-512F.
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx512f")]
    unsafe fn hash_reads_avx512<F: FnMut(usize, usize, u64)>(&self, reads: &[&[u8]], f: &mut F) {
        use crate::simd::{load_table_x8, roll_index_x8, roll_x8, seed_x8, LANES_X8};
        use std::arch::x86_64::*;

        /// Eight bases of each lane at `off`, one `u64` per lane.
        #[inline]
        #[target_feature(enable = "avx512f")]
        unsafe fn lane_words(group: &[&[u8]], off: usize) -> __m512i {
            let w = |lane: usize| -> i64 {
                (group.get_unchecked(lane).as_ptr().add(off) as *const u64)
                    .read_unaligned()
                    .to_le() as i64
            };
            _mm512_set_epi64(w(7), w(6), w(5), w(4), w(3), w(2), w(1), w(0))
        }

        let k = self.k as usize;
        let len = reads[0].len();
        let fq = load_table_x8(&self.consts.fwd_q);
        let rq = load_table_x8(&self.consts.rev_q);
        // `vpermi2q` selects from sixteen; only codes 0..4 are ever used.
        let mut seed16 = [0u64; 16];
        let mut seed_rc16 = [0u64; 16];
        seed16[..4].copy_from_slice(&self.consts.seed);
        seed_rc16[..4].copy_from_slice(&self.consts.seed_rc);

        let mut seed_f = [0u64; LANES_X8];
        let mut seed_r = [0u64; LANES_X8];
        let mut block = [0u64; LANES_X8 * 8];

        let mut base = 0;
        while base + LANES_X8 <= reads.len() {
            let group = &reads[base..base + LANES_X8];
            // Seed all eight lanes on one chain; `seed_x8` needs a whole word
            // to work with, so very short reads keep the scalar seeding.
            let (mut fwd, mut rev) = if len >= 8 {
                seed_x8(k, len, |off| lane_words(group, off), &seed16, &seed_rc16)
            } else {
                for (lane, r) in group.iter().enumerate() {
                    seed_f[lane] = base_forward_hash(r, self.k);
                    seed_r[lane] = base_reverse_hash(r, self.k);
                }
                (
                    _mm512_loadu_si512(seed_f.as_ptr() as *const __m512i),
                    _mm512_loadu_si512(seed_r.as_ptr() as *const __m512i),
                )
            };
            _mm512_storeu_si512(seed_f.as_mut_ptr() as *mut __m512i, fwd);
            _mm512_storeu_si512(seed_r.as_mut_ptr() as *mut __m512i, rev);
            for lane in 0..LANES_X8 {
                f(base + lane, 0, canonical(seed_f[lane], seed_r[lane]));
            }

            let mut i = k;
            while i + 8 <= len {
                // One load per lane per stream covers the next eight steps.
                let mut ow = lane_words(group, i - k);
                let mut iw = lane_words(group, i);
                for e in 0..8 {
                    let (nf, nr) = roll_x8(fwd, rev, roll_index_x8(ow, iw), fq, rq);
                    fwd = nf;
                    rev = nr;
                    _mm512_storeu_si512(
                        block.as_mut_ptr().add(e * LANES_X8) as *mut __m512i,
                        _mm512_add_epi64(fwd, rev),
                    );
                    ow = _mm512_srli_epi64::<8>(ow);
                    iw = _mm512_srli_epi64::<8>(iw);
                }
                for e in 0..8 {
                    let pos = i - k + 1 + e;
                    for lane in 0..LANES_X8 {
                        f(base + lane, pos, *block.get_unchecked(e * LANES_X8 + lane));
                    }
                }
                i += 8;
            }
            if i < len {
                _mm512_storeu_si512(seed_f.as_mut_ptr() as *mut __m512i, fwd);
                _mm512_storeu_si512(seed_r.as_mut_ptr() as *mut __m512i, rev);
                for (lane, r) in group.iter().enumerate() {
                    self.roll_tail(base + lane, r, i, seed_f[lane], seed_r[lane], f);
                }
            }
            base += LANES_X8;
        }
        for (offset, r) in reads[base..].iter().enumerate() {
            self.hash_one(base + offset, r, f);
        }
    }

    /// AVX2: the same eight-read group as the AVX-512 kernel, driven as two
    /// 256-bit halves. The roll table is read with a gather — 128 bytes, always
    /// L1-resident — because a 16-way blend tree costs 15 `vpblendvb` and is
    /// slower than rolling the reads one at a time.
    ///
    /// # Safety
    /// Requires AVX2.
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    unsafe fn hash_reads_avx2<F: FnMut(usize, usize, u64)>(&self, reads: &[&[u8]], f: &mut F) {
        use crate::simd::{roll_index_x4, roll_x4, LANES_X4};
        use std::arch::x86_64::*;

        const LANES: usize = 2 * LANES_X4;

        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn lane_words(group: &[&[u8]], half: usize, off: usize) -> __m256i {
            let w = |lane: usize| -> i64 {
                (group
                    .get_unchecked(half * LANES_X4 + lane)
                    .as_ptr()
                    .add(off) as *const u64)
                    .read_unaligned()
                    .to_le() as i64
            };
            _mm256_set_epi64x(w(3), w(2), w(1), w(0))
        }

        let k = self.k as usize;
        let len = reads[0].len();
        let fq = self.consts.fwd_q.as_ptr() as *const i64;
        let rq = self.consts.rev_q.as_ptr() as *const i64;

        let mut seed_f = [0u64; LANES];
        let mut seed_r = [0u64; LANES];
        let mut block = [0u64; LANES * 8];

        let mut base = 0;
        while base + LANES <= reads.len() {
            let group = &reads[base..base + LANES];
            for (lane, r) in group.iter().enumerate() {
                seed_f[lane] = base_forward_hash(r, self.k);
                seed_r[lane] = base_reverse_hash(r, self.k);
                f(base + lane, 0, canonical(seed_f[lane], seed_r[lane]));
            }
            let mut fwd = [
                _mm256_loadu_si256(seed_f.as_ptr() as *const __m256i),
                _mm256_loadu_si256(seed_f.as_ptr().add(LANES_X4) as *const __m256i),
            ];
            let mut rev = [
                _mm256_loadu_si256(seed_r.as_ptr() as *const __m256i),
                _mm256_loadu_si256(seed_r.as_ptr().add(LANES_X4) as *const __m256i),
            ];

            let mut i = k;
            while i + 8 <= len {
                let mut ow = [lane_words(group, 0, i - k), lane_words(group, 1, i - k)];
                let mut iw = [lane_words(group, 0, i), lane_words(group, 1, i)];
                for e in 0..8 {
                    for h in 0..2 {
                        let (nf, nr) = roll_x4(fwd[h], rev[h], roll_index_x4(ow[h], iw[h]), fq, rq);
                        fwd[h] = nf;
                        rev[h] = nr;
                        _mm256_storeu_si256(
                            block.as_mut_ptr().add(e * LANES + h * LANES_X4) as *mut __m256i,
                            _mm256_add_epi64(fwd[h], rev[h]),
                        );
                        ow[h] = _mm256_srli_epi64::<8>(ow[h]);
                        iw[h] = _mm256_srli_epi64::<8>(iw[h]);
                    }
                }
                for e in 0..8 {
                    let pos = i - k + 1 + e;
                    for lane in 0..LANES {
                        f(base + lane, pos, *block.get_unchecked(e * LANES + lane));
                    }
                }
                i += 8;
            }
            if i < len {
                for h in 0..2 {
                    _mm256_storeu_si256(
                        seed_f.as_mut_ptr().add(h * LANES_X4) as *mut __m256i,
                        fwd[h],
                    );
                    _mm256_storeu_si256(
                        seed_r.as_mut_ptr().add(h * LANES_X4) as *mut __m256i,
                        rev[h],
                    );
                }
                for (lane, r) in group.iter().enumerate() {
                    self.roll_tail(base + lane, r, i, seed_f[lane], seed_r[lane], f);
                }
            }
            base += LANES;
        }
        for (offset, r) in reads[base..].iter().enumerate() {
            self.hash_one(base + offset, r, f);
        }
    }

    /// NEON: roll four reads per group as two 128-bit lane-pairs.
    ///
    /// # Safety
    /// Requires NEON.
    #[cfg(target_arch = "aarch64")]
    #[target_feature(enable = "neon")]
    unsafe fn hash_reads_neon<F: FnMut(usize, usize, u64)>(&self, reads: &[&[u8]], f: &mut F) {
        let mut base = 0;
        while base + LANES <= reads.len() {
            let group: [&[u8]; LANES] = [
                reads[base],
                reads[base + 1],
                reads[base + 2],
                reads[base + 3],
            ];
            self.hash_group_neon(base, &group, f);
            base += LANES;
        }
        for (offset, r) in reads[base..].iter().enumerate() {
            self.hash_one(base + offset, r, f);
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

    /// Deterministic pseudo-random N-free reads.
    fn make_reads(n: usize, len: usize) -> Vec<Vec<u8>> {
        (0..n)
            .map(|r| {
                (0..len)
                    .map(|i| b"ACGT"[(r * 7 + i * 5 + i / 3) % 4])
                    .collect()
            })
            .collect()
    }

    fn collect_scalar(h: &ReadHasher, reads: &[&[u8]]) -> BTreeMap<(usize, usize), u64> {
        let mut m = BTreeMap::new();
        for (ri, r) in reads.iter().enumerate() {
            h.hash_one(ri, r, &mut |ri, pos, x| {
                m.insert((ri, pos), x);
            });
        }
        m
    }

    /// Every SIMD backend must agree with the scalar roll bit-for-bit, over
    /// read counts that exercise a full lane group and a ragged tail, and read
    /// lengths that exercise whole eight-step vector blocks plus a tail.
    #[test]
    fn simd_backends_match_scalar() {
        // (8, 5) and (16, 9) straddle the eight-base word the SIMD seeding
        // needs; (7, 33) and (3, 12) leave a ragged tail after the vector loop.
        for &(n_reads, len) in &[
            (20usize, 100usize),
            (8, 40),
            (7, 33),
            (3, 12),
            (16, 9),
            (8, 5),
            (9, 8),
        ] {
            let owned = make_reads(n_reads, len);
            let reads: Vec<&[u8]> = owned.iter().map(|r| r.as_slice()).collect();
            for k in [1u16, 2, 4, 7, 9, 31, 40] {
                if k as usize > len {
                    continue;
                }
                let h = ReadHasher::new(k).unwrap();
                let want = collect_scalar(&h, &reads);
                assert_eq!(
                    want,
                    reference(&reads, k),
                    "scalar vs blind k={k} len={len}"
                );

                let mut got = BTreeMap::new();
                h.hash_reads(&reads, |ri, pos, x| {
                    got.insert((ri, pos), x);
                })
                .unwrap();
                assert_eq!(want, got, "dispatched backend k={k} len={len} n={n_reads}");

                #[cfg(target_arch = "x86_64")]
                {
                    if std::arch::is_x86_feature_detected!("avx2") {
                        let mut v = BTreeMap::new();
                        unsafe {
                            h.hash_reads_avx2(&reads, &mut |ri, pos, x| {
                                v.insert((ri, pos), x);
                            })
                        };
                        assert_eq!(want, v, "avx2 k={k} len={len} n={n_reads}");
                    }
                    if std::arch::is_x86_feature_detected!("avx512f") {
                        let mut v = BTreeMap::new();
                        unsafe {
                            h.hash_reads_avx512(&reads, &mut |ri, pos, x| {
                                v.insert((ri, pos), x);
                            })
                        };
                        assert_eq!(want, v, "avx512 k={k} len={len} n={n_reads}");
                    }
                }
                #[cfg(target_arch = "aarch64")]
                {
                    if std::arch::is_aarch64_feature_detected!("neon") {
                        let mut v = BTreeMap::new();
                        unsafe {
                            h.hash_reads_neon(&reads, &mut |ri, pos, x| {
                                v.insert((ri, pos), x);
                            })
                        };
                        assert_eq!(want, v, "neon k={k} len={len} n={n_reads}");
                    }
                }
            }
        }
    }

    /// Every k-mer of every read must be emitted exactly once.
    #[test]
    fn emits_each_kmer_once() {
        let owned = make_reads(20, 100);
        let reads: Vec<&[u8]> = owned.iter().map(|r| r.as_slice()).collect();
        for k in [1u16, 7, 31] {
            let mut count = 0usize;
            let mut seen = BTreeMap::new();
            ReadHasher::new(k)
                .unwrap()
                .hash_reads(&reads, |ri, pos, _| {
                    count += 1;
                    assert!(seen.insert((ri, pos), ()).is_none(), "duplicate {ri}/{pos}");
                })
                .unwrap();
            assert_eq!(count, reads.len() * (100 - k as usize + 1), "k={k}");
        }
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
