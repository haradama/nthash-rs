//! **Lane-parallel canonical hashing of one long sequence.**
//!
//! A single rolling-hash chain is serial: `srol` feeds the next `srol`, so one
//! stream is bounded by that dependency (~2.5 ns/base here) no matter how much
//! issue width is left over. But ntHash's value is a pure function of the
//! k-mer — position-independent — so **disjoint ranges of the same sequence are
//! independent chains**. [`parallel`](crate::parallel) exploits that across
//! threads; this module exploits it across SIMD lanes, on one core, by cutting
//! the k-mer range into `L` equal parts and rolling all `L` in lockstep.
//!
//! Backends are detected at run time and share their per-step arithmetic with
//! [`batch`](crate::batch): AVX-512F (8 lanes, one `vpermi2q` per roll table),
//! AVX2 (8 lanes as two 256-bit halves, L1 gather), else scalar. All produce
//! **bit-identical** results.
//!
//! Splitting reorders emission, so the two entry points differ in what they
//! promise:
//!
//! - [`for_each_canonical_unordered`] hands every hash to a callback in an
//!   **unspecified order** — the fastest path, for consumers that fold (a
//!   sketch, a Bloom filter, a counter).
//! - [`hash_canonical_into`] / [`hash_canonical`] materialise **in position
//!   order**, `out[i]` being the hash of the k-mer at `i`.
//!
//! Like [`BlindNtHash`](crate::blind), this assumes **N-free A/C/G/T** input;
//! use [`NtHash`](crate::kmer) for sequences that may contain `N`.
//!
//! ```
//! use nthash_rs::wide::hash_canonical;
//! # fn main() -> nthash_rs::Result<()> {
//! let hashes = hash_canonical(b"ACGTACGTACGTACGT", 4)?;
//! assert_eq!(hashes.len(), 16 - 4 + 1);
//! # Ok(()) }
//! ```

use crate::{
    kmer::{base_forward_hash, base_reverse_hash},
    tables::{roll_forward, roll_index, RollConsts},
    util::canonical,
    NtHashError, Result,
};

/// Lanes the vector backends drive at once.
const LANES: usize = 8;

/// Minimum k-mers per lane before splitting is worth it.
///
/// The vector roll advances eight k-mers at a time, so fewer than nine per lane
/// leaves nothing for it to do and the whole group falls to the scalar tail.
const MIN_PER_LANE: usize = 9;

/// Number of k-mers in `seq`, or the reason there are none.
fn kmer_count(seq: &[u8], k: u16) -> Result<usize> {
    if k == 0 {
        return Err(NtHashError::InvalidK);
    }
    let ku = k as usize;
    if seq.len() < ku {
        return Err(NtHashError::SequenceTooShort {
            seq_len: seq.len(),
            k,
        });
    }
    Ok(seq.len() - ku + 1)
}

/// How many k-mers each lane covers: `0` means "roll the whole sequence
/// serially".
///
/// Beyond [`MIN_PER_LANE`], splitting also has to earn back its seed hash.
/// Serial rolling seeds *once*; the eight-lane seed costs about three times
/// that, because it walks one base per step where the scalar tetramer table
/// walks four. Measured break-even is near `k / 6` k-mers per lane — below it,
/// serial rolling wins (`k = 127` only turns over at ~24 k-mers per lane).
fn lane_span(n_kmers: usize, k: u16) -> usize {
    let per = n_kmers / LANES;
    if per < MIN_PER_LANE || 6 * per < k as usize {
        0
    } else {
        per
    }
}

/// Canonical hash of every k-mer, in position order (`out[i]` is the k-mer at
/// `i`). `out` is resized to the exact k-mer count — reuse one `Vec` across
/// calls to avoid reallocating.
///
/// # Errors
/// [`NtHashError::InvalidK`] if `k == 0`, [`NtHashError::SequenceTooShort`] if
/// `seq.len() < k`.
pub fn hash_canonical_into(seq: &[u8], k: u16, out: &mut Vec<u64>) -> Result<()> {
    let n = kmer_count(seq, k)?;
    // Grow only when short and never re-zero: the kernel writes every slot
    // exactly once (asserted by `matches_serial_across_lengths_and_k`), so
    // pre-filling a reused buffer would just be a second pass over the output —
    // at these rates, as expensive as the hashing itself.
    if out.len() < n {
        out.resize(n, 0);
    }
    out.truncate(n);
    fill(seq, k, out)
}

/// Allocating form of [`hash_canonical_into`].
///
/// # Errors
/// Same as [`hash_canonical_into`].
pub fn hash_canonical(seq: &[u8], k: u16) -> Result<Vec<u64>> {
    let n = kmer_count(seq, k)?;
    // `vec![0; n]` takes the `alloc_zeroed` path — fresh pages arrive zeroed, so
    // unlike `Vec::resize` this costs no pre-pass over the buffer.
    let mut out = vec![0u64; n];
    fill(seq, k, &mut out)?;
    Ok(out)
}

/// Write each k-mer's canonical hash to its own slot of `out`.
///
/// `out.len()` must be the k-mer count of `seq`; [`parallel`](crate::parallel)
/// calls this per thread-chunk so the two levels of splitting compose.
pub(crate) fn fill(seq: &[u8], k: u16, out: &mut [u64]) -> Result<()> {
    for_each_canonical_unordered(seq, k, |pos, h| {
        if let Some(slot) = out.get_mut(pos) {
            *slot = h;
        }
    })
}

/// Invoke `f(pos, canonical_hash)` once for every k-mer, in an **unspecified
/// order**. The fastest path: nothing is materialised, so a folding consumer
/// runs at the kernel's full rate.
///
/// # Errors
/// Same as [`hash_canonical_into`].
pub fn for_each_canonical_unordered<F: FnMut(usize, u64)>(
    seq: &[u8],
    k: u16,
    mut f: F,
) -> Result<()> {
    let n = kmer_count(seq, k)?;
    let consts = RollConsts::new(k);
    let per = lane_span(n, k);

    if per == 0 {
        roll_range(seq, k, &consts, 0, n, &mut f);
        return Ok(());
    }

    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx512f") {
            // SAFETY: guarded by run-time AVX-512F detection.
            unsafe { lanes_avx512(seq, k, &consts, per, &mut f) };
            roll_range(seq, k, &consts, LANES * per, n, &mut f);
            return Ok(());
        }
        if std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: guarded by run-time AVX2 detection.
            unsafe { lanes_avx2(seq, k, &consts, per, &mut f) };
            roll_range(seq, k, &consts, LANES * per, n, &mut f);
            return Ok(());
        }
    }

    #[allow(unreachable_code)]
    {
        roll_range(seq, k, &consts, 0, n, &mut f);
        Ok(())
    }
}

/// Scalar roll over the k-mer positions `start..end`, seeded from scratch.
/// Also the reference the SIMD kernels are checked against.
#[inline]
fn roll_range<F: FnMut(usize, u64)>(
    seq: &[u8],
    k: u16,
    c: &RollConsts,
    start: usize,
    end: usize,
    f: &mut F,
) {
    if start >= end {
        return;
    }
    let ku = k as usize;
    let mut fwd = base_forward_hash(&seq[start..], k);
    let mut rev = base_reverse_hash(&seq[start..], k);
    f(start, canonical(fwd, rev));
    for pos in (start + 1)..end {
        let (nf, nr) = roll_forward(c, fwd, rev, roll_index(seq[pos - 1], seq[pos - 1 + ku]));
        fwd = nf;
        rev = nr;
        f(pos, canonical(fwd, rev));
    }
}

/// Continue a lane's scalar roll from `state` over the k-mer positions `range`.
/// Finishes the ragged end of each SIMD lane.
#[cfg(target_arch = "x86_64")]
#[inline]
fn roll_lane_tail<F: FnMut(usize, u64)>(
    seq: &[u8],
    k: u16,
    c: &RollConsts,
    range: std::ops::Range<usize>,
    state: (u64, u64),
    f: &mut F,
) {
    let ku = k as usize;
    let (mut fwd, mut rev) = state;
    for pos in range {
        let (nf, nr) = roll_forward(c, fwd, rev, roll_index(seq[pos - 1], seq[pos - 1 + ku]));
        fwd = nf;
        rev = nr;
        f(pos, canonical(fwd, rev));
    }
}

/// AVX-512F: eight lanes, one `vpermi2q` per roll table per step.
///
/// # Safety
/// Requires AVX-512F. `per` must be the value [`lane_span`] returned (`> 0`),
/// so every byte the kernel touches is in bounds.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn lanes_avx512<F: FnMut(usize, u64)>(
    seq: &[u8],
    k: u16,
    c: &RollConsts,
    per: usize,
    f: &mut F,
) {
    use crate::simd::{load_table_x8, roll_index_x8, roll_x8, seed_x8};
    use std::arch::x86_64::*;

    let ku = k as usize;
    let base = seq.as_ptr();
    let fq = load_table_x8(&c.fwd_q);
    let rq = load_table_x8(&c.rev_q);

    // Eight bases of each lane's stream at `off` bytes past that lane's start.
    let words = |off: usize| -> __m512i {
        let w = |lane: usize| -> i64 {
            (base.add(lane * per + off) as *const u64)
                .read_unaligned()
                .to_le() as i64
        };
        _mm512_set_epi64(w(7), w(6), w(5), w(4), w(3), w(2), w(1), w(0))
    };

    // Seed all eight lanes on one chain. `per >= MIN_PER_LANE` guarantees every
    // lane has a whole word to read, and passing `k + 8` as the span tells
    // `seed_x8` it never has to back a load off the end.
    let mut seed16 = [0u64; 16];
    let mut seed_rc16 = [0u64; 16];
    seed16[..4].copy_from_slice(&c.seed);
    seed_rc16[..4].copy_from_slice(&c.seed_rc);
    let (mut fwd, mut rev) = seed_x8(ku, ku + 8, words, &seed16, &seed_rc16);

    let mut state_f = [0u64; LANES];
    let mut state_r = [0u64; LANES];
    _mm512_storeu_si512(state_f.as_mut_ptr() as *mut __m512i, fwd);
    _mm512_storeu_si512(state_r.as_mut_ptr() as *mut __m512i, rev);
    for lane in 0..LANES {
        f(lane * per, canonical(state_f[lane], state_r[lane]));
    }

    let mut block = [0u64; LANES * 8];
    let mut step = 0usize;
    while step + 8 < per {
        let mut ow = words(step);
        let mut iw = words(step + ku);
        for e in 0..8 {
            let (nf, nr) = roll_x8(fwd, rev, roll_index_x8(ow, iw), fq, rq);
            fwd = nf;
            rev = nr;
            _mm512_storeu_si512(
                block.as_mut_ptr().add(e * LANES) as *mut __m512i,
                _mm512_add_epi64(fwd, rev),
            );
            ow = _mm512_srli_epi64::<8>(ow);
            iw = _mm512_srli_epi64::<8>(iw);
        }
        // Emit lane-major: the multiply leaves the inner loop and each lane's
        // eight hashes land on consecutive positions, which is what an ordered
        // consumer wants to store.
        for lane in 0..LANES {
            let start = lane * per + step + 1;
            for e in 0..8 {
                f(start + e, *block.get_unchecked(e * LANES + lane));
            }
        }
        step += 8;
    }

    _mm512_storeu_si512(state_f.as_mut_ptr() as *mut __m512i, fwd);
    _mm512_storeu_si512(state_r.as_mut_ptr() as *mut __m512i, rev);
    for lane in 0..LANES {
        let s = lane * per;
        roll_lane_tail(
            seq,
            k,
            c,
            (s + step + 1)..(s + per),
            (state_f[lane], state_r[lane]),
            f,
        );
    }
}

/// AVX2: the same eight lanes as two 256-bit halves.
///
/// # Safety
/// Requires AVX2; `per` as for [`lanes_avx512`].
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn lanes_avx2<F: FnMut(usize, u64)>(
    seq: &[u8],
    k: u16,
    c: &RollConsts,
    per: usize,
    f: &mut F,
) {
    use crate::simd::{roll_index_x4, roll_x4, LANES_X4};
    use std::arch::x86_64::*;

    let ku = k as usize;
    let base = seq.as_ptr();
    let fq = c.fwd_q.as_ptr() as *const i64;
    let rq = c.rev_q.as_ptr() as *const i64;

    let words = |half: usize, off: usize| -> __m256i {
        let w = |lane: usize| -> i64 {
            (base.add((half * LANES_X4 + lane) * per + off) as *const u64)
                .read_unaligned()
                .to_le() as i64
        };
        _mm256_set_epi64x(w(3), w(2), w(1), w(0))
    };

    let mut state_f = [0u64; LANES];
    let mut state_r = [0u64; LANES];
    for lane in 0..LANES {
        let s = lane * per;
        state_f[lane] = base_forward_hash(&seq[s..], k);
        state_r[lane] = base_reverse_hash(&seq[s..], k);
        f(s, canonical(state_f[lane], state_r[lane]));
    }
    let mut fwd = [
        _mm256_loadu_si256(state_f.as_ptr() as *const __m256i),
        _mm256_loadu_si256(state_f.as_ptr().add(LANES_X4) as *const __m256i),
    ];
    let mut rev = [
        _mm256_loadu_si256(state_r.as_ptr() as *const __m256i),
        _mm256_loadu_si256(state_r.as_ptr().add(LANES_X4) as *const __m256i),
    ];

    let mut block = [0u64; LANES * 8];
    let mut step = 0usize;
    while step + 8 < per {
        let mut ow = [words(0, step), words(1, step)];
        let mut iw = [words(0, step + ku), words(1, step + ku)];
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
        // Emit lane-major: the multiply leaves the inner loop and each lane's
        // eight hashes land on consecutive positions, which is what an ordered
        // consumer wants to store.
        for lane in 0..LANES {
            let start = lane * per + step + 1;
            for e in 0..8 {
                f(start + e, *block.get_unchecked(e * LANES + lane));
            }
        }
        step += 8;
    }

    for h in 0..2 {
        _mm256_storeu_si256(
            state_f.as_mut_ptr().add(h * LANES_X4) as *mut __m256i,
            fwd[h],
        );
        _mm256_storeu_si256(
            state_r.as_mut_ptr().add(h * LANES_X4) as *mut __m256i,
            rev[h],
        );
    }
    for lane in 0..LANES {
        let s = lane * per;
        roll_lane_tail(
            seq,
            k,
            c,
            (s + step + 1)..(s + per),
            (state_f[lane], state_r[lane]),
            f,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BlindNtHashBuilder;

    /// Position-ordered canonical hashes via the single-stream reference.
    fn reference(seq: &[u8], k: u16) -> Vec<u64> {
        let mut v = Vec::new();
        BlindNtHashBuilder::new(seq)
            .k(k)
            .num_hashes(1)
            .finish()
            .unwrap()
            .for_each_canonical(|_pos, h| v.push(h));
        v
    }

    fn dna(len: usize) -> Vec<u8> {
        (0..len).map(|i| b"ACGT"[(i * 7 + i / 13) % 4]).collect()
    }

    /// Every backend, at lengths that straddle the split threshold, must match
    /// the serial rolling hasher exactly.
    #[test]
    fn matches_serial_across_lengths_and_k() {
        // 102/140/200 sit just above the split threshold for k = 31 (nine
        // k-mers per lane), 100 just below it; 513+ exercise long runs.
        for len in [16usize, 100, 102, 140, 200, 513, 1024, 4096, 5001] {
            let seq = dna(len);
            for k in [1u16, 2, 5, 16, 31, 63] {
                if k as usize > len {
                    continue;
                }
                let want = reference(&seq, k);
                assert_eq!(hash_canonical(&seq, k).unwrap(), want, "k={k} len={len}");

                // Unordered: same multiset, every position exactly once.
                let mut got = vec![None; want.len()];
                for_each_canonical_unordered(&seq, k, |pos, h| {
                    assert!(got[pos].is_none(), "duplicate pos {pos}");
                    got[pos] = Some(h);
                })
                .unwrap();
                let got: Vec<u64> = got.into_iter().map(|x| x.expect("missing pos")).collect();
                assert_eq!(got, want, "unordered k={k} len={len}");
            }
        }
    }

    /// The split must be exercised, not silently skipped, on a long sequence.
    #[test]
    fn long_sequence_actually_splits() {
        let seq = dna(100_000);
        let n = seq.len() - 31 + 1;
        assert!(lane_span(n, 31) > 0, "test would not cover the SIMD path");
        assert_eq!(hash_canonical(&seq, 31).unwrap(), reference(&seq, 31));
    }

    #[test]
    fn rejects_bad_arguments() {
        assert_eq!(hash_canonical(b"ACGT", 0), Err(NtHashError::InvalidK));
        assert_eq!(
            hash_canonical(b"ACG", 4),
            Err(NtHashError::SequenceTooShort { seq_len: 3, k: 4 })
        );
    }

    #[test]
    fn into_reuses_the_buffer() {
        let seq = dna(2048);
        let mut buf = vec![0u64; 9999];
        hash_canonical_into(&seq, 31, &mut buf).unwrap();
        assert_eq!(buf.len(), seq.len() - 31 + 1);
        assert_eq!(buf, reference(&seq, 31));
    }
}
