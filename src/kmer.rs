//! Canonical **ntHash** implementation for *contiguous* k‑mers.
//!
//! This module provides a panic‑free, idiomatic Rust port of the C++ `NtHash`
//! class.  It computes rolling hashes over DNA k‑mers (A/C/G/T/N) in **O(1)**
//! time per base (after an initial **O(k)** seed computation), skipping over
//! any windows that contain 'N' exactly as the reference does.
//!
//! All heavy bit‑twiddling is delegated to the `tables` (split‑rotate) and
//! `constants` (lookup tables) modules, plus `util::extend_hashes` for
//! generating extra hash values per k‑mer.
//!
//! Additionally, a Rust‑idiomatic **builder + iterator** facade
//! (`NtHashBuilder` / `NtHashIter`) is provided.

use crate::{
    constants::*,
    tables::{base_code, srol, srol_n, sror, RollConsts},
    util::{canonical, extend_hashes},
    NtHashError, // unified crate-level error
};

/// Convenient alias for fallible operations in this module.
pub type Result<T> = crate::Result<T>;

/// Rolling k‑mer hasher over a contiguous DNA sequence.
///
/// - Initialization is deferred until the first valid k‑mer (skips any
///   windows containing `N`).
/// - `roll()` / `roll_back()` advance by one base, handling skips transparently.
/// - Each valid k‑mer emits `num_hashes` values: the canonical hash plus
///   extra mixes.
pub struct NtHash<'a> {
    seq: &'a [u8],
    k: u16,
    pos: usize,
    initialized: bool,
    fwd_hash: u64,
    rev_hash: u64,
    hashes: Vec<u64>,
    /// Per-base rolling constants for this `k` (see [`RollConsts`]).
    consts: RollConsts,
}

impl<'a> NtHash<'a> {
    /// Create a new `NtHash` starting at `pos`.
    ///
    /// # Arguments
    ///
    /// * `seq` – full DNA sequence (`A,C,G,T,N` recognized; others treated as `N`)
    /// * `k` – k‑mer length (> 0)
    /// * `num_hashes` – how many hash values per k‑mer
    /// * `pos` – starting index
    ///
    /// # Errors
    ///
    /// Returns if `k == 0`, `seq.len() < k`, or `pos` too large.
    pub fn new(seq: &'a [u8], k: u16, num_hashes: u8, pos: usize) -> Result<Self> {
        if k == 0 {
            return Err(NtHashError::InvalidK);
        }
        let len = seq.len();
        let k_usz = k as usize;
        if len < k_usz {
            return Err(NtHashError::SequenceTooShort { seq_len: len, k });
        }
        if pos > len - k_usz {
            return Err(NtHashError::PositionOutOfRange { pos, seq_len: len });
        }
        Ok(Self {
            seq,
            k,
            pos,
            initialized: false,
            fwd_hash: 0,
            rev_hash: 0,
            hashes: vec![0; num_hashes as usize],
            consts: RollConsts::new(k),
        })
    }

    /// Advance forward by one base, skipping over k‑mers with `N`.
    /// Returns `true` if a new valid hash was produced.
    pub fn roll(&mut self) -> bool {
        if self.advance() {
            self.update_hashes();
            true
        } else {
            false
        }
    }

    /// Advance by one base (with `N`-skip), updating the forward/reverse hashes
    /// and position but **not** the materialized hash buffer.
    ///
    /// This is the state-only core of [`roll`](Self::roll) (which is `advance`
    /// plus [`update_hashes`](Self::update_hashes)). The canonical-only fast
    /// path uses it to skip the per-k-mer `extend_hashes` + heap buffer entirely.
    #[inline]
    pub(crate) fn advance(&mut self) -> bool {
        if !self.initialized {
            return self.init_advance();
        }
        let k_usz = self.k as usize;
        if self.pos >= self.seq.len() - k_usz {
            return false;
        }
        let incoming = self.seq[self.pos + k_usz];
        if SEED_TAB[incoming as usize] == SEED_N {
            self.pos += k_usz;
            return self.init_advance();
        }
        let outgoing = self.seq[self.pos];
        self.fwd_hash = next_forward_hash(&self.consts, self.fwd_hash, outgoing, incoming);
        self.rev_hash = next_reverse_hash(&self.consts, self.rev_hash, outgoing, incoming);
        self.pos += 1;
        true
    }

    /// Move backward by one base, skipping over k‑mers with `N`.
    pub fn roll_back(&mut self) -> bool {
        if !self.initialized && !self.init() {
            return false;
        }
        if self.pos == 0 {
            return false;
        }
        let incoming = self.seq[self.pos - 1];
        if SEED_TAB[incoming as usize] == SEED_N {
            if self.pos < self.k as usize {
                return false;
            }
            self.pos -= self.k as usize;
            return self.init();
        }
        let outgoing = self.seq[self.pos + self.k as usize - 1];
        self.fwd_hash = prev_forward_hash(&self.consts, self.fwd_hash, outgoing, incoming);
        self.rev_hash = prev_reverse_hash(&self.consts, self.rev_hash, outgoing, incoming);
        self.update_hashes();
        self.pos -= 1;
        true
    }

    /// Peek the next k‑mer without mutating self.
    pub fn peek(&mut self) -> bool {
        if self.pos >= self.seq.len() - self.k as usize {
            return false;
        }
        let incoming = self.seq[self.pos + self.k as usize];
        self.peek_char(incoming)
    }

    /// Peek with an explicit incoming byte.
    pub fn peek_char(&mut self, incoming: u8) -> bool {
        if !self.initialized && !self.init() {
            return false;
        }
        if SEED_TAB[incoming as usize] == SEED_N {
            return false;
        }
        let outgoing = self.seq[self.pos];
        let fwd = next_forward_hash(&self.consts, self.fwd_hash, outgoing, incoming);
        let rev = next_reverse_hash(&self.consts, self.rev_hash, outgoing, incoming);
        self.fill_hash_buffer(fwd, rev);
        true
    }

    /// Peek the previous k‑mer without mutating self.
    pub fn peek_back(&mut self) -> bool {
        if self.pos == 0 {
            return false;
        }
        let incoming = self.seq[self.pos - 1];
        self.peek_back_char(incoming)
    }

    /// Peek backward with explicit incoming byte.
    pub fn peek_back_char(&mut self, incoming: u8) -> bool {
        if !self.initialized && !self.init() {
            return false;
        }
        if SEED_TAB[incoming as usize] == SEED_N {
            return false;
        }
        let outgoing = self.seq[self.pos + self.k as usize - 1];
        let fwd = prev_forward_hash(&self.consts, self.fwd_hash, outgoing, incoming);
        let rev = prev_reverse_hash(&self.consts, self.rev_hash, outgoing, incoming);
        self.fill_hash_buffer(fwd, rev);
        true
    }

    /// Returns the most recent hash buffer.
    #[inline(always)]
    pub fn hashes(&self) -> &[u64] {
        &self.hashes
    }

    /// Returns the current k‑mer start index.
    #[inline(always)]
    pub fn pos(&self) -> usize {
        self.pos
    }

    /// Returns the forward‑strand hash.
    #[inline(always)]
    pub fn forward_hash(&self) -> u64 {
        self.fwd_hash
    }

    /// Returns the reverse‑complement hash.
    #[inline(always)]
    pub fn reverse_hash(&self) -> u64 {
        self.rev_hash
    }

    /// Initialize on the first valid k‑mer (materializing the hash buffer).
    fn init(&mut self) -> bool {
        if self.init_advance() {
            self.update_hashes();
            true
        } else {
            false
        }
    }

    /// Position on the first valid k‑mer, updating the forward/reverse hashes but
    /// not the materialized buffer.
    fn init_advance(&mut self) -> bool {
        let k_usz = self.k as usize;
        while self.pos <= self.seq.len() - k_usz {
            let mut skip = 0;
            if has_invalid_base(&self.seq[self.pos..], k_usz, &mut skip) {
                self.pos += skip + 1;
                continue;
            }
            self.fwd_hash = base_forward_hash(&self.seq[self.pos..], self.k);
            self.rev_hash = base_reverse_hash(&self.seq[self.pos..], self.k);
            self.initialized = true;
            return true;
        }
        false
    }

    #[inline(always)]
    fn update_hashes(&mut self) {
        extend_hashes(
            self.fwd_hash,
            self.rev_hash,
            self.k as u32,
            &mut self.hashes,
        );
    }

    #[inline(always)]
    fn fill_hash_buffer(&mut self, fwd: u64, rev: u64) {
        extend_hashes(fwd, rev, self.k as u32, &mut self.hashes);
    }
}

#[inline(always)]
pub fn has_invalid_base(seq: &[u8], k: usize, pos_n: &mut usize) -> bool {
    if let Some(idx) = seq[..k]
        .iter()
        .rposition(|&c| SEED_TAB[c as usize] == SEED_N)
    {
        *pos_n = idx;
        true
    } else {
        false
    }
}

#[inline]
pub fn base_forward_hash(seq: &[u8], k: u16) -> u64 {
    let k = k as usize;
    let mut h = 0_u64;

    for chunk in seq[..k - k % 4].chunks_exact(4) {
        h = srol_n(h, 4);

        // build 0‑255 index with 8‑bit wrapping
        let idx = (CONVERT_TAB[chunk[0] as usize] as usize) * 64
            + (CONVERT_TAB[chunk[1] as usize] as usize) * 16
            + (CONVERT_TAB[chunk[2] as usize] as usize) * 4
            + CONVERT_TAB[chunk[3] as usize] as usize;
        h ^= TETRAMER_TAB[idx & 0xFF];
    }

    h = srol_n(h, (k % 4) as u32);
    match k % 4 {
        3 => {
            let idx = (CONVERT_TAB[seq[k - 3] as usize] as usize) * 16
                + (CONVERT_TAB[seq[k - 2] as usize] as usize) * 4
                + CONVERT_TAB[seq[k - 1] as usize] as usize;
            h ^= TRIMER_TAB[idx & 0x3F];
        }
        2 => {
            let idx = (CONVERT_TAB[seq[k - 2] as usize] as usize) * 4
                + CONVERT_TAB[seq[k - 1] as usize] as usize;
            h ^= DIMER_TAB[idx & 0x0F];
        }
        1 => h ^= SEED_TAB[seq[k - 1] as usize],
        _ => {}
    }
    h
}

#[inline]
pub fn base_reverse_hash(seq: &[u8], k: u16) -> u64 {
    let k = k as usize;
    let mut h = 0_u64;

    // Handle the 'tail' (k % 4 = 1,2,3)
    match k % 4 {
        3 => {
            let idx = (RC_CONVERT_TAB[seq[k - 1] as usize] as usize) * 16
                + (RC_CONVERT_TAB[seq[k - 2] as usize] as usize) * 4
                + RC_CONVERT_TAB[seq[k - 3] as usize] as usize;
            h ^= TRIMER_TAB[idx & 0x3F];
        }
        2 => {
            let idx = (RC_CONVERT_TAB[seq[k - 1] as usize] as usize) * 4
                + RC_CONVERT_TAB[seq[k - 2] as usize] as usize;
            h ^= DIMER_TAB[idx & 0x0F];
        }
        1 => {
            let c = seq[k - 1] & CP_OFF;
            h ^= SEED_TAB[c as usize];
        }
        _ => {}
    }

    // Process full 4‑mer chunks in reverse order
    let mut i = k - k % 4;
    while i >= 4 {
        // split‑rotate the accumulator by 4
        h = srol_n(h, 4);

        // build 4‑mer index, mask to 8 bits
        let idx = (RC_CONVERT_TAB[seq[i - 1] as usize] as usize) * 64
            + (RC_CONVERT_TAB[seq[i - 2] as usize] as usize) * 16
            + (RC_CONVERT_TAB[seq[i - 3] as usize] as usize) * 4
            + RC_CONVERT_TAB[seq[i - 4] as usize] as usize;
        h ^= TETRAMER_TAB[idx & 0xFF];

        i -= 4;
    }
    h
}

#[inline(always)]
fn next_forward_hash(c: &RollConsts, prev: u64, char_out: u8, char_in: u8) -> u64 {
    let mut h = srol(prev);
    h ^= c.seed[base_code(char_in)];
    h ^= c.srolk[base_code(char_out)];
    h
}

#[inline(always)]
fn prev_forward_hash(c: &RollConsts, prev: u64, char_out: u8, char_in: u8) -> u64 {
    let mut h = prev ^ c.srolk[base_code(char_in)];
    h ^= c.seed[base_code(char_out)];
    sror(h)
}

#[inline(always)]
fn next_reverse_hash(c: &RollConsts, prev: u64, char_out: u8, char_in: u8) -> u64 {
    let mut h = prev ^ c.srolk_rc[base_code(char_in)];
    h ^= c.seed_rc[base_code(char_out)];
    sror(h)
}

#[inline(always)]
fn prev_reverse_hash(c: &RollConsts, prev: u64, char_out: u8, char_in: u8) -> u64 {
    let mut h = srol(prev);
    h ^= c.seed_rc[base_code(char_in)];
    h ^= c.srolk_rc[base_code(char_out)];
    h
}

// -------------------------------------------------------------------------
// Builder + Iterator facade
// -------------------------------------------------------------------------

/// Configure and consume a rolling‐hash computation as an iterator.
pub struct NtHashBuilder<'a> {
    seq: &'a [u8],
    k: u16,
    num_hashes: u8,
    pos: usize,
}

impl<'a> NtHashBuilder<'a> {
    /// Begin building over `seq`.
    pub fn new(seq: &'a [u8]) -> Self {
        NtHashBuilder {
            seq,
            k: 0,
            num_hashes: 1,
            pos: 0,
        }
    }

    /// Set the k‑mer length.
    pub fn k(mut self, k: u16) -> Self {
        self.k = k;
        self
    }

    /// Set how many hashes per k‑mer.
    pub fn num_hashes(mut self, m: u8) -> Self {
        self.num_hashes = m;
        self
    }

    /// Set the starting position.
    pub fn pos(mut self, pos: usize) -> Self {
        self.pos = pos;
        self
    }

    /// Finalize into an iterator.
    pub fn finish(self) -> Result<NtHashIter<'a>> {
        let hasher = NtHash::new(self.seq, self.k, self.num_hashes, self.pos)?;
        Ok(NtHashIter {
            hasher,
            done: false,
        })
    }
}

/// Iterator yielding `(pos, Vec<u64>)` for each valid k‑mer.
pub struct NtHashIter<'a> {
    hasher: NtHash<'a>,
    done: bool,
}

impl<'a> NtHashIter<'a> {
    /// Drive the hasher to completion, invoking `f(pos, &hashes)` for each
    /// valid k‑mer **without allocating**.
    ///
    /// The [`Iterator`] impl yields an owned `Vec<u64>` per k‑mer, which means
    /// one heap allocation for every window — the dominant cost when hashing
    /// long sequences. This method instead lends the hasher's internal buffer
    /// to `f`, so the whole traversal allocates nothing. Prefer it whenever the
    /// hashes are consumed immediately rather than collected.
    ///
    /// # Examples
    ///
    /// ```
    /// # use nthash_rs::NtHashBuilder;
    /// # fn main() -> nthash_rs::Result<()> {
    /// let mut count = 0usize;
    /// NtHashBuilder::new(b"ACGTACGTACGT").k(4).finish()?
    ///     .for_each_hash(|_pos, hashes| count += hashes.len());
    /// assert!(count > 0);
    /// # Ok(()) }
    /// ```
    #[inline]
    pub fn for_each_hash<F: FnMut(usize, &[u64])>(mut self, mut f: F) {
        if self.done {
            return;
        }
        while self.hasher.roll() {
            f(self.hasher.pos(), self.hasher.hashes());
        }
    }

    /// Drive the hasher to completion, invoking `f(pos, canonical_hash)` for each
    /// valid k‑mer — the **fastest** way to obtain the single canonical hash per
    /// k‑mer.
    ///
    /// Unlike [`for_each_hash`](Self::for_each_hash), this bypasses the per-k‑mer
    /// `extend_hashes` mixing and the heap hash buffer entirely, computing the
    /// canonical hash straight from the rolling forward/reverse state. It ignores
    /// `num_hashes` (always yielding exactly the canonical hash), so use it when
    /// one hash per k‑mer suffices — e.g. minimizers, sketching, membership.
    #[inline]
    pub fn for_each_canonical<F: FnMut(usize, u64)>(mut self, mut f: F) {
        if self.done {
            return;
        }
        let h = &mut self.hasher;
        let seq = h.seq;
        let k = h.k as usize;
        let len = seq.len();
        if len < k {
            return;
        }
        let last = len - k; // last valid k-mer start index

        // Prime once: position on the first valid k-mer. This hoists the
        // `initialized` check out of the hot loop below (it recurs only on an
        // actual N-skip), leaving a tight roll with each byte read exactly once.
        if !h.advance() {
            return;
        }
        loop {
            f(h.pos, canonical(h.fwd_hash, h.rev_hash));
            if h.pos >= last {
                return;
            }
            let incoming = seq[h.pos + k];
            if SEED_TAB[incoming as usize] == SEED_N {
                // N-skip: jump past the ambiguous base and re-init, exactly as
                // `roll` does; the loop top then emits the next run's first k-mer.
                h.pos += k;
                h.initialized = false;
                if !h.advance() {
                    return;
                }
                continue;
            }
            let outgoing = seq[h.pos];
            h.fwd_hash = next_forward_hash(&h.consts, h.fwd_hash, outgoing, incoming);
            h.rev_hash = next_reverse_hash(&h.consts, h.rev_hash, outgoing, incoming);
            h.pos += 1;
        }
    }
}

impl<'a> Iterator for NtHashIter<'a> {
    type Item = (usize, Vec<u64>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        if !self.hasher.roll() {
            self.done = true;
            return None;
        }
        let out = (self.hasher.pos(), self.hasher.hashes().to_owned());
        Some(out)
    }
}

impl<'a> IntoIterator for NtHashBuilder<'a> {
    type Item = (usize, Vec<u64>);
    type IntoIter = NtHashIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.finish().expect("invalid NtHashBuilder configuration")
    }
}
