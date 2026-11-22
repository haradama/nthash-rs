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
    NtHashError, Result,
};

/// Rolling k‑mer hasher over a contiguous DNA sequence.
///
/// - Initialization is deferred until the first valid k‑mer (skips any
///   windows containing `N`).
/// - `roll()` / `roll_back()` advance by one base, handling skips transparently.
/// - Each valid k‑mer emits `num_hashes` values: the canonical hash plus
///   extra mixes.
/// - On a fresh hasher, the first call to any of `roll`, `roll_back`, `peek*`
///   only *initializes*: it positions on the first valid k‑mer at or after the
///   starting `pos` and reports whether one exists.
/// - A call that returns `false` leaves the hasher on its current k‑mer, so
///   `pos()`, `hashes()` and further rolling in either direction stay valid.
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
            // Every window overlapping the ambiguous base is invalid; resume
            // scanning just past it.
            return match resume_after(self.seq, self.k, self.pos + k_usz) {
                Some((pos, fwd, rev)) => {
                    (self.pos, self.fwd_hash, self.rev_hash) = (pos, fwd, rev);
                    true
                }
                None => false,
            };
        }
        let outgoing = self.seq[self.pos];
        (self.fwd_hash, self.rev_hash) =
            self.consts
                .roll(self.fwd_hash, self.rev_hash, outgoing, incoming);
        self.pos += 1;
        true
    }

    /// Move backward by one base, skipping over k‑mers with `N`.
    /// Returns `true` if a new valid hash was produced.
    pub fn roll_back(&mut self) -> bool {
        if !self.initialized {
            return self.init();
        }
        if self.pos == 0 {
            return false;
        }
        let incoming = self.seq[self.pos - 1];
        if SEED_TAB[incoming as usize] == SEED_N {
            return self.init_back();
        }
        let outgoing = self.seq[self.pos + self.k as usize - 1];
        self.fwd_hash = prev_forward_hash(&self.consts, self.fwd_hash, outgoing, incoming);
        self.rev_hash = prev_reverse_hash(&self.consts, self.rev_hash, outgoing, incoming);
        self.update_hashes();
        self.pos -= 1;
        true
    }

    /// Peek the next k‑mer: write its hashes to the buffer returned by
    /// [`hashes`](Self::hashes) without moving the window. `pos()`,
    /// `forward_hash()` and `reverse_hash()` keep describing the current k‑mer.
    ///
    /// Returns `false` (buffer untouched) at the end of the sequence or when the
    /// incoming base is ambiguous.
    pub fn peek(&mut self) -> bool {
        if !self.initialized {
            return self.init();
        }
        if self.pos >= self.seq.len() - self.k as usize {
            return false;
        }
        let incoming = self.seq[self.pos + self.k as usize];
        self.peek_char(incoming)
    }

    /// Like [`peek`](Self::peek), with an explicit incoming byte.
    pub fn peek_char(&mut self, incoming: u8) -> bool {
        if !self.initialized {
            return self.init();
        }
        if SEED_TAB[incoming as usize] == SEED_N {
            return false;
        }
        let outgoing = self.seq[self.pos];
        let (fwd, rev) = self
            .consts
            .roll(self.fwd_hash, self.rev_hash, outgoing, incoming);
        self.fill_hash_buffer(fwd, rev);
        true
    }

    /// Peek the previous k‑mer; the backward counterpart of
    /// [`peek`](Self::peek).
    pub fn peek_back(&mut self) -> bool {
        if !self.initialized {
            return self.init();
        }
        if self.pos == 0 {
            return false;
        }
        let incoming = self.seq[self.pos - 1];
        self.peek_back_char(incoming)
    }

    /// Like [`peek_back`](Self::peek_back), with an explicit incoming byte.
    pub fn peek_back_char(&mut self, incoming: u8) -> bool {
        if !self.initialized {
            return self.init();
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

    /// Position on the first valid k‑mer at or after the current position,
    /// updating the forward/reverse hashes but not the materialized buffer.
    #[inline]
    fn init_advance(&mut self) -> bool {
        self.init_from(self.pos)
    }

    /// Position on the first valid k‑mer starting at or after `start`.
    ///
    /// State is committed only on success: when no valid k‑mer remains, `pos`
    /// and the hashes are left untouched, so the hasher still describes the
    /// k‑mer it was on (and can be rolled back from it).
    fn init_from(&mut self, start: usize) -> bool {
        match first_valid_window(self.seq, self.k, start) {
            Some((pos, fwd, rev)) => {
                (self.pos, self.fwd_hash, self.rev_hash) = (pos, fwd, rev);
                self.initialized = true;
                true
            }
            None => false,
        }
    }

    /// Position on the last valid k‑mer starting strictly before the current
    /// position (the backward counterpart of [`init_from`](Self::init_from)),
    /// materializing the hash buffer. Leaves the state untouched on failure.
    fn init_back(&mut self) -> bool {
        let k_usz = self.k as usize;
        let mut pos = self.pos;
        while pos > 0 {
            pos -= 1;
            let Some(idx) = first_invalid_base(&self.seq[pos..pos + k_usz]) else {
                self.set_window(pos);
                self.update_hashes();
                return true;
            };
            // Every window covering the ambiguous base at `pos + idx` is
            // invalid; the next candidate ends just before it.
            let bad = pos + idx;
            if bad < k_usz {
                return false;
            }
            pos = bad - k_usz + 1; // decremented at the loop top
        }
        false
    }

    /// Recompute the forward/reverse hashes from scratch for the (valid) window
    /// starting at `pos` and move there.
    fn set_window(&mut self, pos: usize) {
        self.fwd_hash = base_forward_hash(&self.seq[pos..], self.k);
        self.rev_hash = base_reverse_hash(&self.seq[pos..], self.k);
        self.pos = pos;
        self.initialized = true;
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

/// Find the first valid k‑mer starting at or after `start` and hash it from
/// scratch, returning `(pos, forward_hash, reverse_hash)`.
///
/// Requires `seq.len() >= k`.
fn first_valid_window(seq: &[u8], k: u16, start: usize) -> Option<(usize, u64, u64)> {
    let k_usz = k as usize;
    let last = seq.len() - k_usz;
    let mut pos = start;
    while pos <= last {
        if let Some(idx) = last_invalid_base(&seq[pos..pos + k_usz]) {
            pos += idx + 1;
            continue;
        }
        let window = &seq[pos..];
        return Some((
            pos,
            base_forward_hash(window, k),
            base_reverse_hash(window, k),
        ));
    }
    None
}

/// Out-of-line N-skip for the rolling loops: the first valid k‑mer after the
/// ambiguous base at index `n_pos` (every window overlapping it is invalid).
///
/// Deliberately never inlined, and taking and returning plain values rather
/// than `&mut NtHash`. Folding the O(k) re-initialization into a hot loop — or
/// letting a pointer to the rolling state escape into a call — makes the
/// optimizer keep that state on the stack and reload it on every base, which
/// puts a store-to-load round trip on the loop-carried hash dependency.
#[cold]
#[inline(never)]
fn resume_after(seq: &[u8], k: u16, n_pos: usize) -> Option<(usize, u64, u64)> {
    first_valid_window(seq, k, n_pos + 1)
}

/// Index of the last ambiguous (non‑ACGT) base in `window`, if any.
#[inline(always)]
fn last_invalid_base(window: &[u8]) -> Option<usize> {
    window.iter().rposition(|&c| SEED_TAB[c as usize] == SEED_N)
}

/// Index of the first ambiguous (non‑ACGT) base in `window`, if any.
#[inline(always)]
fn first_invalid_base(window: &[u8]) -> Option<usize> {
    window.iter().position(|&c| SEED_TAB[c as usize] == SEED_N)
}

/// Forward‑strand hash of the k‑mer `seq[..k]`, computed from scratch in O(k).
///
/// # Panics
/// Panics if `seq.len() < k`.
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

/// Reverse‑complement hash of the k‑mer `seq[..k]`, computed from scratch in
/// O(k).
///
/// # Panics
/// Panics if `seq.len() < k`.
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
fn prev_forward_hash(c: &RollConsts, prev: u64, char_out: u8, char_in: u8) -> u64 {
    let mut h = prev ^ c.srolk[base_code(char_in)];
    h ^= c.seed[base_code(char_out)];
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
///
/// Every item owns a freshly allocated `Vec`. When one hash per k‑mer is
/// enough, [`canonical`](Self::canonical) gives an allocation-free iterator;
/// to borrow all `num_hashes` values without allocating, use
/// [`for_each_hash`](Self::for_each_hash).
pub struct NtHashIter<'a> {
    hasher: NtHash<'a>,
    done: bool,
}

impl<'a> NtHashIter<'a> {
    /// Turn this into an iterator over `(pos, canonical_hash)` — one `u64` per
    /// valid k‑mer and **no allocation**, several times faster than iterating
    /// the `Vec`-yielding [`Iterator`] impl.
    ///
    /// Like [`for_each_canonical`](Self::for_each_canonical), it ignores
    /// `num_hashes` and yields exactly the canonical hash (what `hashes()[0]`
    /// would be), but as a regular iterator that composes with adapters.
    ///
    /// # Examples
    ///
    /// ```
    /// # use nthash_rs::NtHashBuilder;
    /// # fn main() -> nthash_rs::Result<()> {
    /// let min = NtHashBuilder::new(b"ACGTNACGTACG").k(4).finish()?
    ///     .canonical()
    ///     .map(|(_pos, hash)| hash)
    ///     .min();
    /// assert!(min.is_some());
    /// # Ok(()) }
    /// ```
    #[inline]
    pub fn canonical(self) -> NtHashCanonicalIter<'a> {
        let h = self.hasher;
        NtHashCanonicalIter {
            seq: h.seq,
            k: h.k,
            consts: h.consts,
            pos: h.pos,
            fwd: h.fwd_hash,
            rev: h.rev_hash,
            state: if self.done {
                CanonicalState::Done
            } else if h.initialized {
                CanonicalState::Rolling
            } else {
                CanonicalState::Fresh
            },
        }
    }

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

        // Prime once: position on the first valid k-mer (or, for a hasher that
        // was already rolled, step to the next one). This hoists the
        // `initialized` check out of the hot loop below, leaving a tight roll
        // with each byte read exactly once.
        if !h.advance() {
            return;
        }
        // Keep the rolling state in locals so it lives in registers for the whole
        // traversal; the hasher itself is not touched again.
        let consts = &h.consts;
        let (mut pos, mut fwd, mut rev) = (h.pos, h.fwd_hash, h.rev_hash);
        loop {
            f(pos, canonical(fwd, rev));
            if pos >= last {
                return;
            }
            let incoming = seq[pos + k];
            if SEED_TAB[incoming as usize] == SEED_N {
                // N-skip: jump past the ambiguous base and re-init, exactly as
                // `roll` does; the loop top then emits the next run's first k-mer.
                match resume_after(seq, h.k, pos + k) {
                    Some(next) => (pos, fwd, rev) = next,
                    None => return,
                }
                continue;
            }
            (fwd, rev) = consts.roll(fwd, rev, seq[pos], incoming);
            pos += 1;
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

/// Allocation-free iterator yielding `(pos, canonical_hash)` for each valid
/// k‑mer. Created by [`NtHashIter::canonical`].
///
/// It carries its own copy of the rolling state (rather than wrapping an
/// [`NtHash`]) so that no pointer to it ever leaves `next`, which lets the
/// optimizer keep the hashes in registers across a `for` loop.
pub struct NtHashCanonicalIter<'a> {
    seq: &'a [u8],
    k: u16,
    consts: RollConsts,
    /// Start of the current k‑mer (or the start position while `Fresh`).
    pos: usize,
    fwd: u64,
    rev: u64,
    state: CanonicalState,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum CanonicalState {
    /// Nothing emitted yet: the first k‑mer at or after `pos` is still to be found.
    Fresh,
    /// `pos` / `fwd` / `rev` describe the k‑mer emitted last.
    Rolling,
    /// The sequence is exhausted.
    Done,
}

impl Iterator for NtHashCanonicalIter<'_> {
    type Item = (usize, u64);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let k = self.k as usize;
        let next = match self.state {
            CanonicalState::Rolling => {
                if self.pos >= self.seq.len() - k {
                    None
                } else {
                    let incoming = self.seq[self.pos + k];
                    if SEED_TAB[incoming as usize] == SEED_N {
                        resume_after(self.seq, self.k, self.pos + k)
                    } else {
                        let (fwd, rev) =
                            self.consts
                                .roll(self.fwd, self.rev, self.seq[self.pos], incoming);
                        Some((self.pos + 1, fwd, rev))
                    }
                }
            }
            CanonicalState::Fresh => {
                self.state = CanonicalState::Rolling;
                first_valid_window(self.seq, self.k, self.pos)
            }
            CanonicalState::Done => None,
        };
        match next {
            Some((pos, fwd, rev)) => {
                (self.pos, self.fwd, self.rev) = (pos, fwd, rev);
                Some((pos, canonical(fwd, rev)))
            }
            None => {
                self.state = CanonicalState::Done;
                None
            }
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        // At most one k-mer per remaining start position; `N`s may remove any.
        let last = self.seq.len() - self.k as usize;
        let remaining = match self.state {
            CanonicalState::Fresh => last - self.pos + 1,
            CanonicalState::Rolling => last - self.pos,
            CanonicalState::Done => 0,
        };
        (0, Some(remaining))
    }
}

impl std::iter::FusedIterator for NtHashCanonicalIter<'_> {}
