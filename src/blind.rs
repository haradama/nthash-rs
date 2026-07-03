//! **Streaming (“blind”) ntHash** for *contiguous* k‑mers.
//!
//! Unlike [`kmer::NtHash`](crate::kmer), which owns an entire DNA string and
//! skips windows containing  'N', **`BlindNtHash` works on a pre‑cleaned input
//! where every window of length *k* is known to be valid**.  
//!
//! The hasher maintains an *exact* k‑base sliding window in a small ring
//! buffer and lets the caller **feed the next / previous character** to move
//! the window forward or backward.
//!
//! Heavy bit‑twiddling is delegated to the `tables` (split‑rotate) and
//! `constants` (lookup tables) modules, plus `util::extend_hashes` for
//! generating extra hash values per window.
//!
//! A Rust‑idiomatic **builder + iterator** facade
//! (`BlindNtHashBuilder` / `BlindNtHashIter`) is included for ergonomic
//! streaming over an already‑sanitised sequence.

use std::collections::VecDeque;

use crate::{
    constants::{CP_OFF, SEED_TAB},
    kmer::{base_forward_hash, base_reverse_hash},
    tables::{build_srol_k, srol, sror},
    util::extend_hashes,
    NtHashError, Result,
};

/// Rolling hash over a *fixed‑width* window that the caller rolls manually.
///
/// The window is stored in a `VecDeque<u8>`:
/// - `roll()` removes the **front** base and pushes a new base at the **back**.
/// - `roll_back()` does the opposite.
/// - `peek()` / `peek_back()` compute hashes for the next / previous window
///   **without** mutating internal state.
pub struct BlindNtHash {
    window: VecDeque<u8>,
    k: u16,
    pos: isize,
    fwd_hash: u64,
    rev_hash: u64,
    hashes: Vec<u64>,
    /// Precomputed `srol_table(c, k)` for every byte value `c`. Unlike
    /// [`NtHash`](crate::kmer), the blind variant does not gate out `N`, so it
    /// needs the full byte-indexed table (where `N` maps to `0`) rather than the
    /// 4-constant `RollConsts` (which would fold `N` onto `G`).
    srol_k: Box<[u64; 256]>,
}

impl BlindNtHash {
    /// Create a new `BlindNtHash` whose initial window is `seq[pos..pos+k]`.
    ///
    /// The caller must guarantee that the slice contains **no ambiguous bases
    /// ('N')** – the blind variant will not skip over invalid windows.
    ///
    /// # Errors
    ///
    /// Returns if `k == 0`, `seq.len() < k`, or `pos` too large.
    pub fn new(seq: &[u8], k: u16, num_hashes: u8, pos: isize) -> Result<Self> {
        if k == 0 {
            return Err(NtHashError::InvalidK);
        }
        let len = seq.len();
        let k_usz = k as usize;

        // Guard before `len - k_usz` below, which would otherwise underflow.
        if len < k_usz {
            return Err(NtHashError::SequenceTooShort { seq_len: len, k });
        }

        if pos < 0 || (pos as usize) > len - k_usz {
            return Err(NtHashError::PositionOutOfRange {
                pos: pos as usize,
                seq_len: len,
            });
        }

        let slice = &seq[(pos as usize)..(pos as usize + k_usz)];
        let mut window = VecDeque::with_capacity(k_usz);
        window.extend(slice.iter().copied());

        let fwd_hash = base_forward_hash(slice, k);
        let rev_hash = base_reverse_hash(slice, k);

        let mut hashes = vec![0; num_hashes as usize];
        extend_hashes(fwd_hash, rev_hash, k as u32, &mut hashes);

        Ok(Self {
            window,
            k,
            pos,
            fwd_hash,
            rev_hash,
            hashes,
            srol_k: build_srol_k(k),
        })
    }

    /// Returns `true` if a new valid hash was produced.
    pub fn roll(&mut self, char_in: u8) -> bool {
        let char_out = self
            .window
            .pop_front()
            .expect("window length is always k > 0");
        self.window.push_back(char_in);

        self.fwd_hash = next_forward_hash(&self.srol_k, self.fwd_hash, char_out, char_in);
        self.rev_hash = next_reverse_hash(&self.srol_k, self.rev_hash, char_out, char_in);
        extend_hashes(
            self.fwd_hash,
            self.rev_hash,
            self.k as u32,
            &mut self.hashes,
        );
        self.pos += 1;
        true
    }

    pub fn roll_back(&mut self, char_in: u8) -> bool {
        debug_assert_eq!(self.window.len(), self.k as usize);
        let char_out = self
            .window
            .pop_back()
            .expect("window length is always k > 0");
        self.window.push_front(char_in);

        self.fwd_hash = prev_forward_hash(&self.srol_k, self.fwd_hash, char_out, char_in);
        self.rev_hash = prev_reverse_hash(&self.srol_k, self.rev_hash, char_out, char_in);
        extend_hashes(
            self.fwd_hash,
            self.rev_hash,
            self.k as u32,
            &mut self.hashes,
        );
        self.pos -= 1;
        true
    }

    /// Compute hashes for the **next** window without mutating `self`.
    pub fn peek(&mut self, char_in: u8) {
        let char_out = *self.window.front().unwrap();
        let fwd = next_forward_hash(&self.srol_k, self.fwd_hash, char_out, char_in);
        let rev = next_reverse_hash(&self.srol_k, self.rev_hash, char_out, char_in);
        extend_hashes(fwd, rev, self.k as u32, &mut self.hashes);
    }

    pub fn peek_back(&mut self, char_in: u8) {
        let char_out = *self.window.back().unwrap();
        let fwd = prev_forward_hash(&self.srol_k, self.fwd_hash, char_out, char_in);
        let rev = prev_reverse_hash(&self.srol_k, self.rev_hash, char_out, char_in);
        extend_hashes(fwd, rev, self.k as u32, &mut self.hashes);
    }

    #[inline(always)]
    pub fn hashes(&self) -> &[u64] {
        &self.hashes
    }

    #[inline(always)]
    pub fn pos(&self) -> isize {
        self.pos
    }

    #[inline(always)]
    pub fn forward_hash(&self) -> u64 {
        self.fwd_hash
    }

    #[inline(always)]
    pub fn reverse_hash(&self) -> u64 {
        self.rev_hash
    }
}

#[inline(always)]
fn next_forward_hash(srol_k: &[u64; 256], prev: u64, char_out: u8, char_in: u8) -> u64 {
    srol(prev) ^ SEED_TAB[char_in as usize] ^ srol_k[char_out as usize]
}

#[inline(always)]
fn prev_forward_hash(srol_k: &[u64; 256], prev: u64, char_out: u8, char_in: u8) -> u64 {
    sror(prev ^ srol_k[char_in as usize] ^ SEED_TAB[char_out as usize])
}

#[inline(always)]
fn next_reverse_hash(srol_k: &[u64; 256], prev: u64, char_out: u8, char_in: u8) -> u64 {
    sror(prev ^ srol_k[(char_in & CP_OFF) as usize] ^ SEED_TAB[(char_out & CP_OFF) as usize])
}

#[inline(always)]
fn prev_reverse_hash(srol_k: &[u64; 256], prev: u64, char_out: u8, char_in: u8) -> u64 {
    srol(prev) ^ SEED_TAB[(char_in & CP_OFF) as usize] ^ srol_k[(char_out & CP_OFF) as usize]
}

pub struct BlindNtHashBuilder<'a> {
    seq: &'a [u8],
    k: u16,
    num_hashes: u8,
    start_pos: usize,
}

impl<'a> BlindNtHashBuilder<'a> {
    pub fn new(seq: &'a [u8]) -> Self {
        Self {
            seq,
            k: 0,
            num_hashes: 1,
            start_pos: 0,
        }
    }

    pub fn k(mut self, k: u16) -> Self {
        self.k = k;
        self
    }

    pub fn num_hashes(mut self, m: u8) -> Self {
        self.num_hashes = m;
        self
    }

    pub fn pos(mut self, pos: usize) -> Self {
        self.start_pos = pos;
        self
    }

    pub fn finish(self) -> Result<BlindNtHashIter<'a>> {
        let hasher = BlindNtHash::new(self.seq, self.k, self.num_hashes, self.start_pos as isize)?;
        let end = self.seq.len() - self.k as usize;
        Ok(BlindNtHashIter {
            seq: self.seq,
            end,
            hasher,
            first: true,
        })
    }
}

pub struct BlindNtHashIter<'a> {
    seq: &'a [u8],
    end: usize,
    hasher: BlindNtHash,
    first: bool,
}

impl<'a> BlindNtHashIter<'a> {
    /// Drive the hasher to completion, invoking `f(pos, &hashes)` for each
    /// window **without allocating**.
    ///
    /// The [`Iterator`] impl yields an owned `Vec<u64>` per window (one heap
    /// allocation each); this method lends the internal buffer to `f` instead,
    /// so the traversal allocates nothing. See
    /// [`NtHashIter::for_each_hash`](crate::kmer::NtHashIter::for_each_hash).
    #[inline]
    pub fn for_each_hash<F: FnMut(usize, &[u64])>(mut self, mut f: F) {
        if self.first {
            self.first = false;
            f(self.hasher.pos() as usize, self.hasher.hashes());
        }
        while (self.hasher.pos() as usize) < self.end {
            let cur = self.hasher.pos() as usize;
            let incoming = self.seq[cur + self.hasher.k as usize];
            self.hasher.roll(incoming);
            f(self.hasher.pos() as usize, self.hasher.hashes());
        }
    }
}

impl<'a> Iterator for BlindNtHashIter<'a> {
    type Item = (usize, Vec<u64>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.first {
            self.first = false;
            return Some((self.hasher.pos() as usize, self.hasher.hashes().to_vec()));
        }

        let cur = self.hasher.pos() as usize;
        if cur >= self.end {
            return None;
        }

        let incoming = self.seq[cur + self.hasher.k as usize];
        self.hasher.roll(incoming);

        Some((self.hasher.pos() as usize, self.hasher.hashes().to_vec()))
    }
}

impl<'a> IntoIterator for BlindNtHashBuilder<'a> {
    type Item = (usize, Vec<u64>);
    type IntoIter = BlindNtHashIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.finish()
            .expect("invalid BlindNtHashBuilder configuration")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn too_short_sequence_errors_without_panicking() {
        match BlindNtHash::new(b"ACG", 4, 1, 0) {
            Err(NtHashError::SequenceTooShort { seq_len: 3, k: 4 }) => {}
            other => panic!("expected SequenceTooShort, got {:?}", other.map(|_| ())),
        }
    }
}
