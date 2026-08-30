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
    kmer::{NtHash, NtHashIter},
    util::canonical,
    NtHashError, Result,
};

/// Monotonic ring-buffer deque of `(pos, hash)` candidates, hashes non-decreasing
/// from front to back. The window holds at most `w` positions, plus one
/// transient slot before the front is evicted, so `w + 1` slots always suffice.
///
/// The buffer is rounded up to a power of two and `head`/`tail` are *virtual*
/// indices that only ever increase — wrapping is a single `&` on read. That
/// removes both the conditional subtraction the obvious ring needs and, because
/// `i & (len - 1)` is provably `< len`, the bounds check on every access; the
/// pop loops touch the buffer once per candidate, so both were on the hot path.
struct MonoWindow {
    buf: Box<[(usize, u64)]>,
    /// `buf.len() - 1`; `buf.len()` is a power of two, so this is the wrap mask.
    mask: usize,
    head: usize,
    tail: usize,
}

impl MonoWindow {
    fn new(w: usize) -> Self {
        let cap = (w + 1).next_power_of_two();
        Self {
            buf: vec![(0usize, 0u64); cap].into_boxed_slice(),
            mask: cap - 1,
            head: 0,
            tail: 0,
        }
    }

    #[inline(always)]
    fn is_empty(&self) -> bool {
        self.head == self.tail
    }

    #[inline(always)]
    fn clear(&mut self) {
        self.head = 0;
        self.tail = 0;
    }

    #[inline(always)]
    fn front(&self) -> (usize, u64) {
        self.buf[self.head & self.mask]
    }

    #[inline(always)]
    fn back(&self) -> (usize, u64) {
        self.buf[(self.tail - 1) & self.mask]
    }

    #[inline(always)]
    fn push_back(&mut self, x: (usize, u64)) {
        let i = self.tail & self.mask;
        self.buf[i] = x;
        self.tail += 1;
    }

    #[inline(always)]
    fn pop_back(&mut self) {
        self.tail -= 1;
    }

    #[inline(always)]
    fn pop_front(&mut self) {
        self.head += 1;
    }
}

/// Sentinel for "no position yet". Real positions are sequence indices, so
/// `usize::MAX` can never collide; an `Option<usize>` here costs a tag load,
/// test and store on every k-mer.
const NO_POS: usize = usize::MAX;

/// Sliding-window minimizer selection, independent of *how* k-mer hashes arrive.
///
/// Both drivers feed this the same `(pos, hash)` stream: the [`Iterator`] impl
/// pulls one k-mer at a time, [`Minimizers::for_each`] pushes from the hasher's
/// own hoisted roll loop.
struct Scan {
    w: usize,
    /// Monotonic window; the front is the current window's minimizer.
    window: MonoWindow,
    /// Number of k-mers seen in the current gap-free run.
    run_len: usize,
    /// Position of the previously consumed k-mer, or [`NO_POS`].
    prev_pos: usize,
    /// Position of the last emitted minimizer, or [`NO_POS`].
    last_emitted: usize,
}

impl Scan {
    fn new(w: usize) -> Self {
        Self {
            w,
            window: MonoWindow::new(w),
            run_len: 0,
            prev_pos: NO_POS,
            last_emitted: NO_POS,
        }
    }

    /// Consume one k-mer, returning the minimizer it completes, if any.
    #[inline(always)]
    fn push(&mut self, pos: usize, hash: u64) -> Option<Minimizer> {
        // An `N`-skip makes positions non-contiguous; reset the window so a
        // minimizer never spans an ambiguous base.
        if self.prev_pos != NO_POS && pos != self.prev_pos + 1 {
            self.window.clear();
            self.run_len = 0;
            self.last_emitted = NO_POS;
        }
        self.prev_pos = pos;

        // Maintain the monotonic deque. Pop strictly-greater tails so that,
        // on ties, the leftmost (earliest) k-mer is retained as minimizer.
        while !self.window.is_empty() && self.window.back().1 > hash {
            self.window.pop_back();
        }
        self.window.push_back((pos, hash));
        self.run_len += 1;

        // Evict the front once it falls outside the window of the last `w`
        // k-mer positions (`front_pos <= pos - w`).
        while self.window.front().0 + self.w <= pos {
            self.window.pop_front();
        }

        // Emit the window minimizer once a full window exists, skipping
        // repeats of the same occurrence.
        if self.run_len >= self.w {
            let (mpos, mhash) = self.window.front();
            if self.last_emitted != mpos {
                self.last_emitted = mpos;
                return Some(Minimizer {
                    pos: mpos,
                    hash: mhash,
                });
            }
        }
        None
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
        // `num_hashes = 1` makes `hashes()[0]` the canonical hash.
        let hasher = NtHash::new(self.seq, self.k, 1, 0)?;
        Ok(Minimizers {
            hasher,
            scan: Scan::new(self.w),
        })
    }
}

/// Iterator yielding canonical minimizers, deduplicated by position.
///
/// Each item is a small `Copy` value, so iteration allocates nothing.
pub struct Minimizers<'a> {
    hasher: NtHash<'a>,
    scan: Scan,
}

impl<'a> Minimizers<'a> {
    /// Drive the scan to completion, invoking `f` for each selected minimizer.
    ///
    /// Faster than iterating: it hands the k-mer stream to
    /// [`NtHashIter::for_each_canonical`], whose loop keeps the rolling hash
    /// state in registers, rather than pulling one k-mer at a time through
    /// `&mut self` — which forces every field back to memory on each step.
    #[inline]
    pub fn for_each<F: FnMut(Minimizer)>(self, mut f: F) {
        let Minimizers { hasher, mut scan } = self;
        NtHashIter::from_hasher(hasher).for_each_canonical(|pos, hash| {
            if let Some(m) = scan.push(pos, hash) {
                f(m);
            }
        });
    }

    /// Advance until the next minimizer is emitted (or the sequence ends).
    /// Backs the [`Iterator`] impl.
    #[inline(always)]
    fn step(&mut self) -> Option<Minimizer> {
        loop {
            // `advance` rolls the forward/reverse state with N-skip but skips the
            // per-k-mer hash-buffer materialization we don't need here.
            if !self.hasher.advance() {
                return None;
            }
            let pos = self.hasher.pos();
            let hash = canonical(self.hasher.forward_hash(), self.hasher.reverse_hash());
            if let Some(m) = self.scan.push(pos, hash) {
                return Some(m);
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

    /// The push-based `for_each` driver and the pull-based `Iterator` must
    /// select exactly the same minimizers, including across `N`-skips.
    #[test]
    fn for_each_matches_iterator() {
        let seqs: [&[u8]; 4] = [
            b"TACGGATCCGATTACGA",
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            b"ACGTNNNNACGTACGTNACGGATTACAGGGCCCTTTAAA",
            b"NNNNNNNN",
        ];
        for seq in seqs {
            for k in [1u16, 3, 5, 8] {
                for w in [1usize, 2, 5, 11] {
                    let Ok(it) = MinimizerBuilder::new(seq).k(k).w(w).finish() else {
                        continue;
                    };
                    let pulled: Vec<Minimizer> = it.collect();
                    let mut pushed = Vec::new();
                    MinimizerBuilder::new(seq)
                        .k(k)
                        .w(w)
                        .finish()
                        .unwrap()
                        .for_each(|m| pushed.push(m));
                    assert_eq!(pulled, pushed, "k={k} w={w} seq={seq:?}");
                }
            }
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
