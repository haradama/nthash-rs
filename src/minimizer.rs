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

use crate::{
    kmer::{NtHashBuilder, NtHashCanonicalIter},
    NtHashError, Result,
};

/// Monotonic ring-buffer deque of `(pos, hash)` candidates, hashes non-decreasing
/// from front to back. The window holds at most `w` positions, plus one
/// transient slot before the front is evicted, so a fixed `cap + 1` capacity
/// never overflows — no reallocation, no capacity checks, and index wrap-around
/// by conditional subtraction (never a `%`) in the hot path.
struct MonoWindow {
    buf: Box<[(usize, u64)]>,
    head: usize,
    len: usize,
}

impl MonoWindow {
    /// `cap` is the largest number of candidates the window can hold at once.
    fn new(cap: usize) -> Self {
        Self {
            buf: vec![(0usize, 0u64); cap + 1].into_boxed_slice(),
            head: 0,
            len: 0,
        }
    }

    #[inline(always)]
    fn is_empty(&self) -> bool {
        self.len == 0
    }

    #[inline(always)]
    fn clear(&mut self) {
        self.head = 0;
        self.len = 0;
    }

    #[inline(always)]
    fn front(&self) -> (usize, u64) {
        self.buf[self.head]
    }

    #[inline(always)]
    fn back(&self) -> (usize, u64) {
        let cap = self.buf.len();
        let mut i = self.head + self.len - 1;
        if i >= cap {
            i -= cap;
        }
        self.buf[i]
    }

    #[inline(always)]
    fn push_back(&mut self, x: (usize, u64)) {
        let cap = self.buf.len();
        let mut i = self.head + self.len;
        if i >= cap {
            i -= cap;
        }
        self.buf[i] = x;
        self.len += 1;
    }

    #[inline(always)]
    fn pop_back(&mut self) {
        self.len -= 1;
    }

    #[inline(always)]
    fn pop_front(&mut self) {
        self.head += 1;
        if self.head == self.buf.len() {
            self.head = 0;
        }
        self.len -= 1;
    }
}

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
        // One canonical hash per k-mer, with `N`-containing k-mers skipped.
        let hashes = NtHashBuilder::new(self.seq).k(self.k).finish()?.canonical();
        // The deque holds distinct k-mer positions, so it never needs more slots
        // than the sequence has k-mers — which keeps an oversized `w` from
        // turning into an oversized (or overflowing) allocation.
        let num_kmers = self.seq.len() - self.k as usize + 1;
        Ok(Minimizers {
            hashes,
            w: self.w,
            window: MonoWindow::new(self.w.min(num_kmers)),
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
    /// Canonical hash of every valid k-mer, in order.
    hashes: NtHashCanonicalIter<'a>,
    w: usize,
    /// Monotonic window; the front is the current window's minimizer.
    window: MonoWindow,
    /// Number of k-mers seen in the current gap-free run.
    run_len: usize,
    /// Position of the previously consumed k-mer, for gap detection.
    prev_pos: Option<usize>,
    /// Position of the last emitted minimizer, for deduplication.
    last_emitted: Option<usize>,
}

impl Minimizers<'_> {
    /// Drive the scan to completion, invoking `f` for each selected minimizer in
    /// a single tight loop — no per-item `Iterator::next` boundary. Equivalent
    /// to `for m in self { f(m) }`.
    #[inline]
    pub fn for_each<F: FnMut(Minimizer)>(mut self, mut f: F) {
        while let Some(m) = self.step() {
            f(m);
        }
    }

    /// Advance until the next minimizer is emitted (or the sequence ends).
    /// Shared by [`for_each`](Self::for_each) and the [`Iterator`] impl.
    #[inline(always)]
    fn step(&mut self) -> Option<Minimizer> {
        loop {
            let (pos, hash) = self.hashes.next()?;

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
            while !self.window.is_empty() && self.window.back().1 > hash {
                self.window.pop_back();
            }
            self.window.push_back((pos, hash));
            self.run_len += 1;

            // Evict the front once it falls outside the window of the last `w`
            // k-mer positions (`pos - front_pos >= w`; the front never exceeds
            // `pos`, and this form cannot overflow for a huge `w`).
            while pos - self.window.front().0 >= self.w {
                self.window.pop_front();
            }

            // Emit the window minimizer once a full window exists, skipping
            // repeats of the same occurrence.
            if self.run_len >= self.w {
                let (mpos, mhash) = self.window.front();
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

impl Iterator for Minimizers<'_> {
    type Item = Minimizer;

    #[inline]
    fn next(&mut self) -> Option<Minimizer> {
        self.step()
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
        let hashes: Vec<u64> = NtHashBuilder::new(seq)
            .k(k)
            .num_hashes(1)
            .finish()
            .unwrap()
            .map(|(_, h)| h[0])
            .collect();
        // Leftmost minimum of every window of `w` consecutive k-mers, deduplicated.
        let mut expected: Vec<(usize, u64)> = Vec::new();
        for (start, win) in hashes.windows(w).enumerate() {
            let min = *win.iter().min().unwrap();
            let off = win.iter().position(|&h| h == min).unwrap();
            if expected.last() != Some(&(start + off, min)) {
                expected.push((start + off, min));
            }
        }
        let got: Vec<(usize, u64)> = MinimizerBuilder::new(seq)
            .k(k)
            .w(w)
            .finish()
            .unwrap()
            .map(|m| (m.pos, m.hash))
            .collect();
        assert_eq!(got, expected);
    }

    /// A window far larger than the sequence must neither overflow nor try to
    /// allocate `w` slots; with fewer than `w` k-mers nothing is selected.
    #[test]
    fn oversized_window_is_harmless() {
        for w in [usize::MAX, usize::MAX - 1, 1 << 45] {
            let n = MinimizerBuilder::new(b"ACGTACGTACGT")
                .k(4)
                .w(w)
                .finish()
                .unwrap()
                .count();
            assert_eq!(n, 0, "w={w}");
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
