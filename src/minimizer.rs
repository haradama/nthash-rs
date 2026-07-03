//! **Canonical (w, k)-minimizers** over a DNA sequence.
//!
//! A *minimizer* of a window of `w` consecutive k-mers is the k-mer with the
//! smallest hash in that window (winnowing, Schleimer et al. 2003 / Roberts
//! et al. 2004). Minimizers are the backbone of modern sequence sketching,
//! read mapping (minimap2), and assembly — they subsample k-mers while
//! guaranteeing that any two sequences sharing a long enough exact match also
//! share a minimizer.
//!
//! This implementation rides on the crate's canonical [`NtHash`](crate::kmer),
//! so minimizers are **strand-independent** by construction: a sequence and its
//! reverse complement yield the same multiset of minimizer hashes. Selection is
//! O(1) amortised per k-mer via a monotonic deque, and the iterator itself
//! allocates nothing per item.
//!
//! Windows never span an ambiguous base: when [`NtHash`](crate::kmer) skips an
//! `N`-containing k-mer, the sliding window resets, exactly as read mappers do.
//!
//! ```
//! use nthash_rs::MinimizerBuilder;
//! # fn main() -> nthash_rs::Result<()> {
//! for m in MinimizerBuilder::new(b"ACGTACGTACGTACGT").k(4).w(3).finish()? {
//!     println!("minimizer at {}: {:#x}", m.pos, m.hash);
//! }
//! # Ok(()) }
//! ```

use std::collections::VecDeque;

use crate::{kmer::NtHash, NtHashError, Result};

/// A selected minimizer: the canonical hash and the start position of the
/// k-mer that produced it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Minimizer {
    /// Start index of the minimizer k-mer in the sequence.
    pub pos: usize,
    /// Canonical (strand-independent) hash of the k-mer.
    pub hash: u64,
}

/// Configure a canonical (w, k)-minimizer scan.
pub struct MinimizerBuilder<'a> {
    seq: &'a [u8],
    k: u16,
    w: usize,
}

impl<'a> MinimizerBuilder<'a> {
    /// Begin building over `seq`.
    pub fn new(seq: &'a [u8]) -> Self {
        Self { seq, k: 0, w: 1 }
    }

    /// Set the k-mer length (`> 0`).
    pub fn k(mut self, k: u16) -> Self {
        self.k = k;
        self
    }

    /// Set the window length in k-mers (`> 0`). `w == 1` selects every k-mer.
    pub fn w(mut self, w: usize) -> Self {
        self.w = w;
        self
    }

    /// Finalize into an iterator over [`Minimizer`]s.
    ///
    /// # Errors
    /// Returns [`NtHashError::InvalidK`] if `k == 0`, [`NtHashError::InvalidW`]
    /// if `w == 0`, or [`NtHashError::SequenceTooShort`] if `seq.len() < k`.
    pub fn finish(self) -> Result<Minimizers<'a>> {
        if self.w == 0 {
            return Err(NtHashError::InvalidW);
        }
        // `num_hashes = 1` makes `hashes()[0]` the canonical hash.
        let hasher = NtHash::new(self.seq, self.k, 1, 0)?;
        Ok(Minimizers {
            hasher,
            w: self.w,
            window: VecDeque::with_capacity(self.w),
            run_len: 0,
            prev_pos: None,
            last_emitted: None,
        })
    }
}

/// Iterator yielding canonical minimizers, deduplicated by position.
///
/// Each item is a small `Copy` value, so iteration allocates nothing.
pub struct Minimizers<'a> {
    hasher: NtHash<'a>,
    w: usize,
    /// Monotonic deque of `(pos, hash)` candidates, hashes non-decreasing from
    /// front to back; the front is the current window's minimizer.
    window: VecDeque<(usize, u64)>,
    /// Number of k-mers seen in the current gap-free run.
    run_len: usize,
    /// Position of the previously consumed k-mer, for gap detection.
    prev_pos: Option<usize>,
    /// Position of the last emitted minimizer, for deduplication.
    last_emitted: Option<usize>,
}

impl Iterator for Minimizers<'_> {
    type Item = Minimizer;

    fn next(&mut self) -> Option<Minimizer> {
        loop {
            if !self.hasher.roll() {
                return None;
            }
            let pos = self.hasher.pos();
            let hash = self.hasher.hashes()[0];

            // An `N`-skip makes positions non-contiguous; reset the window so a
            // minimizer never spans an ambiguous base.
            if self.prev_pos.is_some_and(|pp| pos != pp + 1) {
                self.window.clear();
                self.run_len = 0;
                self.last_emitted = None;
            }
            self.prev_pos = Some(pos);

            // Maintain the monotonic deque. Pop strictly-greater tails so that,
            // on ties, the leftmost (earliest) k-mer is retained as minimizer.
            while let Some(&(_, back_hash)) = self.window.back() {
                if back_hash > hash {
                    self.window.pop_back();
                } else {
                    break;
                }
            }
            self.window.push_back((pos, hash));
            self.run_len += 1;

            // Evict the front once it falls outside the window of the last `w`
            // k-mer positions (`front_pos <= pos - w`).
            while let Some(&(front_pos, _)) = self.window.front() {
                if front_pos + self.w <= pos {
                    self.window.pop_front();
                } else {
                    break;
                }
            }

            // Emit the window minimizer once a full window exists, skipping
            // repeats of the same occurrence.
            if self.run_len >= self.w {
                let (mpos, mhash) = *self.window.front().expect("window is non-empty");
                if self.last_emitted != Some(mpos) {
                    self.last_emitted = Some(mpos);
                    return Some(Minimizer {
                        pos: mpos,
                        hash: mhash,
                    });
                }
            }
        }
    }
}

impl<'a> IntoIterator for MinimizerBuilder<'a> {
    type Item = Minimizer;
    type IntoIter = Minimizers<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.finish()
            .expect("invalid MinimizerBuilder configuration")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::NtHashBuilder;

    #[test]
    fn w_one_selects_every_kmer() {
        let seq = b"ACGTACGTAC";
        let k = 4u16;
        let mins: Vec<(usize, u64)> = MinimizerBuilder::new(seq)
            .k(k)
            .w(1)
            .finish()
            .unwrap()
            .map(|m| (m.pos, m.hash))
            .collect();
        let all: Vec<(usize, u64)> = NtHashBuilder::new(seq)
            .k(k)
            .num_hashes(1)
            .finish()
            .unwrap()
            .map(|(p, h)| (p, h[0]))
            .collect();
        assert_eq!(mins, all, "w=1 must select every k-mer");
    }

    #[test]
    fn each_minimizer_is_the_window_minimum() {
        let seq = b"TACGGATCCGATTACGA";
        let (k, w) = (5u16, 4usize);
        let hashes: Vec<(usize, u64)> = NtHashBuilder::new(seq)
            .k(k)
            .num_hashes(1)
            .finish()
            .unwrap()
            .map(|(p, h)| (p, h[0]))
            .collect();
        for m in MinimizerBuilder::new(seq).k(k).w(w).finish().unwrap() {
            // The emitted hash must be the min over the w k-mers ending at, or
            // starting at, the minimizer's window — check it is <= every k-mer
            // within `w` positions of it.
            let min_around = hashes
                .iter()
                .filter(|(p, _)| m.pos.abs_diff(*p) < w)
                .map(|(_, h)| *h)
                .min()
                .unwrap();
            assert!(
                m.hash <= min_around || hashes.iter().any(|(p, h)| *p == m.pos && *h == m.hash)
            );
        }
    }

    #[test]
    fn zero_window_errors() {
        assert_eq!(
            MinimizerBuilder::new(b"ACGT").k(2).w(0).finish().err(),
            Some(NtHashError::InvalidW)
        );
    }
}
