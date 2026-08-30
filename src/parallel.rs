//! **Multi-core canonical hashing** of a single long sequence
//! (requires the `parallel` feature).
//!
//! A single rolling-hash chain is serial — each hash depends on the previous —
//! so one thread is bounded by the `srol` dependency latency. But ntHash's hash
//! is a **pure function of the k-mer** (position-independent), so disjoint
//! ranges of the same sequence can be hashed completely independently.
//!
//! This module splits the k-mer range across all CPU cores with rayon; each
//! core writes to a disjoint region of the output (no locks, no atomics),
//! giving near-linear scaling with core count. The two levels of splitting
//! compose: each thread hands its chunk to [`wide`](crate::wide), which splits
//! it again across SIMD lanes.
//!
//! Like [`BlindNtHash`](crate::blind), this assumes **N-free A/C/G/T** input.
//!
//! ```
//! # #[cfg(feature = "parallel")] {
//! use nthash_rs::parallel::par_hash_canonical;
//! let seq = b"ACGTACGTACGTACGT";
//! let hashes = par_hash_canonical(seq, 4).unwrap();
//! assert_eq!(hashes.len(), seq.len() - 4 + 1);
//! # }
//! ```

use rayon::prelude::*;

use crate::{wide, NtHashError, Result};

/// Compute the canonical hash of **every k-mer** of `seq`, in parallel across
/// all CPU cores, returning them in position order (`out[i]` is the hash of the
/// k-mer starting at `i`).
///
/// Input is assumed N-free A/C/G/T (ambiguous bases are not skipped); use
/// [`NtHash`](crate::kmer) for sequences that may contain `N`.
///
/// # Errors
/// [`NtHashError::InvalidK`] if `k == 0`, or [`NtHashError::SequenceTooShort`]
/// if `seq.len() < k`.
pub fn par_hash_canonical(seq: &[u8], k: u16) -> Result<Vec<u64>> {
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
    let n_kmers = seq.len() - ku + 1;
    // `alloc_zeroed` pages, not a written-out pre-pass; every slot is
    // overwritten below.
    let mut out = vec![0u64; n_kmers];

    // A few chunks per thread balances load without excessive overhead. Each
    // chunk owns a disjoint output slice, so no synchronization is needed.
    let threads = rayon::current_num_threads().max(1);
    let chunk = n_kmers.div_ceil(threads * 4).max(1);

    out.par_chunks_mut(chunk)
        .enumerate()
        .for_each(|(ci, oslice)| {
            // A chunk of `m` k-mers starting at `start` spans `m + k - 1` bases.
            let start = ci * chunk;
            let span = &seq[start..start + oslice.len() + ku - 1];
            wide::fill(span, k, oslice).expect("chunk span is always >= k");
        });

    Ok(out)
}

/// XOR-fold every k-mer's canonical hash **in parallel without materializing**
/// them. This is the compute-bound counterpart of [`par_hash_canonical`] (no
/// per-k-mer output write), representative of the scaling an on-the-fly consumer
/// — a sketch, Bloom filter, or counter — can expect. Returns the XOR of all
/// canonical hashes. N-free A/C/G/T input.
///
/// # Errors
/// Same as [`par_hash_canonical`].
pub fn par_hash_canonical_xor(seq: &[u8], k: u16) -> Result<u64> {
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
    let n_kmers = seq.len() - ku + 1;
    let threads = rayon::current_num_threads().max(1);
    let chunk = n_kmers.div_ceil(threads * 4).max(1);
    let n_chunks = n_kmers.div_ceil(chunk);

    let acc = (0..n_chunks)
        .into_par_iter()
        .map(|ci| {
            let start = ci * chunk;
            let end = (start + chunk).min(n_kmers);
            let span = &seq[start..end + ku - 1];
            let mut x = 0u64;
            // XOR is order-independent, so the unordered (fastest) kernel is
            // exactly what this reducer wants.
            wide::for_each_canonical_unordered(span, k, |_pos, h| x ^= h)
                .expect("chunk span is always >= k");
            x
        })
        .reduce(|| 0u64, |a, b| a ^ b);
    Ok(acc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BlindNtHashBuilder;

    #[test]
    fn matches_serial_canonical() {
        let seq = b"ACGTACGTTTGGCCAAGATTACAGATTACACGCGCGTATATACGGGCCC";
        for k in [1u16, 2, 5, 16, 31] {
            if k as usize > seq.len() {
                continue;
            }
            let par = par_hash_canonical(seq, k).unwrap();
            let mut serial = Vec::new();
            BlindNtHashBuilder::new(seq)
                .k(k)
                .num_hashes(1)
                .finish()
                .unwrap()
                .for_each_canonical(|_pos, h| serial.push(h));
            assert_eq!(par, serial, "k={k} (materialized)");

            // The compute-bound reducer must equal the XOR of the same hashes.
            let xor = par_hash_canonical_xor(seq, k).unwrap();
            assert_eq!(xor, serial.iter().fold(0u64, |a, &h| a ^ h), "k={k} (xor)");
        }
    }
}
