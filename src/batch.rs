//! **Batch canonical hashing of many reads.**
//!
//! [`ReadHasher`] is built once for a given `k` and then hashes any number of
//! reads, emitting the canonical hash of every k-mer through a callback. All
//! per-`k` setup (the rolling constants) is paid once, and nothing is
//! allocated per read or per k-mer — the natural shape of read processing in
//! sequencing, where inputs are millions of short reads. Reads may differ in
//! length, but each must be at least `k` long.
//!
//! Each read is rolled independently, with the same forward roll
//! [`NtHash`](crate::kmer) uses.
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
    tables::RollConsts,
    util::canonical,
    NtHashError, Result,
};

/// Batch hasher for many N-free reads.
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
    /// Reads are processed one after another, each from its first k-mer to its
    /// last, so calls arrive ordered by `read_idx`, then `pos`.
    ///
    /// Reads may differ in length, but each must be at least `k` long. Input is
    /// assumed N-free A/C/G/T (ambiguous bases are not skipped).
    ///
    /// # Errors
    /// [`NtHashError::SequenceTooShort`] if any read is shorter than `k` — the
    /// same error [`NtHash::new`](crate::NtHash::new) reports for such a
    /// sequence. Lengths are checked up front, so `f` is never called for a
    /// rejected batch.
    pub fn hash_reads<R, F>(&self, reads: &[R], mut f: F) -> Result<()>
    where
        R: AsRef<[u8]>,
        F: FnMut(usize, usize, u64),
    {
        let k = self.k as usize;
        if let Some(r) = reads.iter().find(|r| r.as_ref().len() < k) {
            return Err(NtHashError::SequenceTooShort {
                seq_len: r.as_ref().len(),
                k: self.k,
            });
        }
        for (read_idx, read) in reads.iter().enumerate() {
            self.hash_one(read_idx, read.as_ref(), &mut f);
        }
        Ok(())
    }

    /// Roll one read (of length `>= k`), emitting every k-mer's canonical hash.
    #[inline]
    fn hash_one<F: FnMut(usize, usize, u64)>(&self, ri: usize, read: &[u8], f: &mut F) {
        let k = self.k as usize;
        let c = &self.consts;
        let mut fwd = base_forward_hash(read, self.k);
        let mut rev = base_reverse_hash(read, self.k);
        f(ri, 0, canonical(fwd, rev));
        // Zipping the outgoing and incoming bases avoids per-base bounds checks.
        for (pos, (&out, &inc)) in read.iter().zip(&read[k..]).enumerate() {
            (fwd, rev) = c.roll(fwd, rev, out, inc);
            f(ri, pos + 1, canonical(fwd, rev));
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
    fn matches_blind_for_every_read() {
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

    /// Reads of different lengths are each hashed in full, exactly as if they
    /// had been hashed alone.
    #[test]
    fn mixed_lengths_match_blind() {
        let reads: &[&[u8]] = &[b"ACGT", b"ACGTACGTACGTTGCA", b"TTGGC", b"GATTACAGATTACA"];
        for k in [1u16, 3, 4] {
            let got = batch(reads, k);
            assert_eq!(got, reference(reads, k), "k={k}");
            for (ri, r) in reads.iter().enumerate() {
                let n = got.keys().filter(|&&(i, _)| i == ri).count();
                assert_eq!(n, r.len() - k as usize + 1, "read {ri}, k={k}");
            }
        }
    }

    /// A read shorter than `k` is rejected with the same error `NtHash::new`
    /// gives, before anything is emitted — even for reads ahead of it.
    #[test]
    fn read_shorter_than_k_errors_before_emitting() {
        let reads: &[&[u8]] = &[b"ACGTAC", b"TTGGCC", b"ACG", b""];
        let mut calls = 0;
        let res = ReadHasher::new(5)
            .unwrap()
            .hash_reads(reads, |_, _, _| calls += 1);
        assert_eq!(res, Err(NtHashError::SequenceTooShort { seq_len: 3, k: 5 }));
        assert_eq!(calls, 0);
        assert_eq!(
            crate::NtHash::new(b"ACG", 5, 1, 0).err(),
            res.err(),
            "batch and single-sequence hashers must agree"
        );
    }

    #[test]
    fn empty_batch_ok() {
        let mut calls = 0;
        let res = ReadHasher::new(4)
            .unwrap()
            .hash_reads::<&[u8], _>(&[], |_, _, _| calls += 1);
        assert_eq!((res, calls), (Ok(()), 0));
    }

    /// Calls arrive ordered by read, then by position within the read.
    #[test]
    fn emits_in_read_then_position_order() {
        let reads = ["ACGTACGT", "TTGGCCAA", "GATTACAG", "CCGGTTAA", "AAAACCCC"];
        let mut seen = Vec::new();
        ReadHasher::new(5)
            .unwrap()
            .hash_reads(&reads, |ri, pos, _| seen.push((ri, pos)))
            .unwrap();
        let expected: Vec<(usize, usize)> = (0..reads.len())
            .flat_map(|ri| (0..4).map(move |pos| (ri, pos)))
            .collect();
        assert_eq!(seen, expected);
    }
}
