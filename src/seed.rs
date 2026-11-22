//! **Spaced-seed ntHash** for *non-contiguous* k‑mers.
//!
//! **`SeedNtHash` computes hashes using spaced seeds**, where only selected
//! positions in the k‑mer are considered (“care sites”). Several seeds can be
//! hashed at once, each yielding `num_hashes` values per window.
//!
//! Windows are **rolled**, not recomputed: sliding by one base only changes the
//! care membership at the *run boundaries* of the mask, so the per-window cost
//! is O(#care-runs). Seeds whose mask is too scattered for that to pay off fall
//! back to an O(#care) recompute, chosen per seed.
//!
//! The per-seed lookup tables depend only on the masks and `k`, never on the
//! sequence. They live in [`SpacedSeeds`], which can be built once and shared
//! by any number of hashers via [`SeedNtHash::with_seeds`] — the right shape
//! when hashing many short reads with the same seeds.
//!
//! Bit-level operations are delegated to `tables`, `constants`, and
//! `util::extend_hashes` for efficient hash computation.
//!
//! A Rust‑idiomatic **builder + iterator** (`SeedNtHashBuilder` / `SeedNtHashIter`)
//! provides ergonomic traversal over valid k‑mers.

use crate::{
    constants::{CP_OFF, SEED_N, SEED_TAB},
    tables::{srol, srol_table, sror},
    util::extend_hashes,
    NtHashError, Result,
};

/// Parses a spaced-seed mask string composed of '0' and '1' characters
/// into a list of indices indicating which positions should be used ("care positions").
///
/// # Errors
/// Returns an error if the mask length does not match `k`, or contains characters other than '0' or '1'.
fn parse_seed_string(mask: &str, k: u16) -> Result<Vec<usize>> {
    if mask.len() != k as usize {
        return Err(NtHashError::InvalidMaskLength {
            mask_len: mask.len(),
            k,
        });
    }
    if !mask.bytes().all(|b| b == b'0' || b == b'1') {
        return Err(NtHashError::InvalidMask);
    }
    Ok(mask
        .bytes()
        .enumerate()
        .filter_map(|(i, b)| if b == b'1' { Some(i) } else { None })
        .collect())
}

/// A base-indexed `[fwd, rev]` contribution table for one care position, used at
/// a run boundary during rolling. Indexed by the raw base byte (256 entries).
struct BoundaryTab {
    off: usize,
    tab: Vec<[u64; 2]>,
}

/// A single spaced seed: the care positions plus everything needed to both
/// recompute a window from scratch and *roll* between windows.
///
/// - `rot` — flattened `care.len() × 256` table `[fwd, rev]` for recomputing a
///   window from scratch (the O(#care) path, used at init and for seeds where
///   rolling doesn't pay).
/// - `left` / `right` — the **run boundaries** of the care pattern: the only
///   positions whose care-membership changes when the window slides by one.
///   Rolling is `srol(fwd)` / `sror(rev)` plus a correction at each boundary, so
///   the per-window cost is O(#runs) instead of O(#care).
struct SeedTab {
    /// Care positions, strictly increasing.
    care: Vec<usize>,
    rot: Vec<[u64; 2]>,
    left: Vec<BoundaryTab>,
    right: Vec<BoundaryTab>,
    /// Whether rolling (O(#boundaries)) is cheaper than recompute (O(#care)) for
    /// this seed. Scattered seeds (many short runs) can have more boundaries than
    /// care positions, where recompute wins — so we pick per seed and never
    /// regress relative to plain recompute.
    rollable: bool,
}

/// Build a 256-entry `[srol_table(base, d_fwd), srol_table(base & CP, d_rev)]`
/// contribution table for one boundary position.
fn boundary_tab(off: usize, d_fwd: u32, d_rev: u32) -> BoundaryTab {
    let mut tab = vec![[0u64; 2]; 256];
    for c in 0u32..256 {
        tab[c as usize] = [
            srol_table(c as u8, d_fwd),
            srol_table((c as u8) & CP_OFF, d_rev),
        ];
    }
    BoundaryTab { off, tab }
}

/// Precompute the recompute and roll tables for one seed's care positions,
/// which must be strictly increasing and `< k`.
fn build_seed_tab(care: Vec<usize>, k: usize) -> SeedTab {
    let mut rot = vec![[0u64; 2]; care.len() * 256];
    for (i, &p) in care.iter().enumerate() {
        let d_fwd = (k - 1 - p) as u32;
        let d_rev = p as u32;
        for c in 0u32..256 {
            rot[i * 256 + c as usize] = [
                srol_table(c as u8, d_fwd),
                srol_table((c as u8) & CP_OFF, d_rev),
            ];
        }
    }

    // Run boundaries: a forward slide turns `srol(fwd)`/`sror(rev)` of the old
    // hash into the new hash, except at care positions whose neighbour on one
    // side is *not* a care position — there the contribution must be corrected.
    // `srol_table` has period lcm(31, 33) = 1023 in its rotation argument, so a
    // rotation of `p - 1` at `p == 0` is taken mod 1023 (→ 1022).
    let mut in_care = vec![false; k];
    for &p in &care {
        in_care[p] = true;
    }
    let mut left = Vec::new();
    let mut right = Vec::new();
    for &p in &care {
        if p == 0 || !in_care[p - 1] {
            // Base leaving the care coverage on a forward slide.
            let d_rev = (p as i64 - 1).rem_euclid(1023) as u32;
            left.push(boundary_tab(p, (k - p) as u32, d_rev));
        }
        if p + 1 == k || !in_care[p + 1] {
            // Base entering the care coverage on a forward slide.
            right.push(boundary_tab(p, (k - 1 - p) as u32, p as u32));
        }
    }

    let rollable = left.len() + right.len() < care.len();
    SeedTab {
        care,
        rot,
        left,
        right,
        rollable,
    }
}

impl SeedTab {
    /// Compute the state for the window at `pos` from scratch — the O(#care)
    /// path, used at init and for seeds where rolling doesn't pay.
    #[inline]
    fn recompute(&self, seq: &[u8], pos: usize) -> SeedState {
        let mut fwd = 0u64;
        let mut rev = 0u64;
        let mut n_invalid = 0usize;
        for (i, &p) in self.care.iter().enumerate() {
            let c = seq[pos + p] as usize;
            let e = self.rot[i * 256 + c];
            fwd ^= e[0];
            rev ^= e[1];
            if SEED_TAB[c] == SEED_N {
                n_invalid += 1;
            }
        }
        SeedState {
            fwd,
            rev,
            n_invalid,
        }
    }

    /// Roll the state `s` of the window at `pos` to the window at `pos + 1`,
    /// correcting only at the care pattern's run boundaries — the
    /// O(#boundaries) path.
    #[inline]
    fn roll(&self, seq: &[u8], pos: usize, s: SeedState) -> SeedState {
        let mut fwd = srol(s.fwd);
        let mut rev = sror(s.rev);
        let mut n_invalid = s.n_invalid;
        // Bases leaving the care coverage.
        for lb in &self.left {
            let c = seq[pos + lb.off] as usize;
            let e = lb.tab[c];
            fwd ^= e[0];
            rev ^= e[1];
            if SEED_TAB[c] == SEED_N {
                n_invalid -= 1;
            }
        }
        // Bases entering the care coverage.
        for rb in &self.right {
            let c = seq[pos + 1 + rb.off] as usize;
            let e = rb.tab[c];
            fwd ^= e[0];
            rev ^= e[1];
            if SEED_TAB[c] == SEED_N {
                n_invalid += 1;
            }
        }
        SeedState {
            fwd,
            rev,
            n_invalid,
        }
    }
}

/// A precompiled set of spaced seeds for a fixed `k`.
///
/// Building the per-seed lookup tables is by far the most expensive part of
/// setting up a [`SeedNtHash`], and depends only on the masks — not on the
/// sequence. Compile the seeds once and hand a reference to each hasher with
/// [`SeedNtHash::with_seeds`] (or [`SeedNtHashBuilder::seeds`]) to hash many
/// sequences without rebuilding them.
///
/// ```
/// use nthash_rs::seed::{SeedNtHash, SpacedSeeds};
/// # fn main() -> nthash_rs::Result<()> {
/// let seeds = SpacedSeeds::from_masks(&["110011", "101101"], 6)?;
/// for read in [&b"ACGTACGTAC"[..], &b"TTGACCATGA"[..]] {
///     let mut h = SeedNtHash::with_seeds(read, &seeds, 1, 0)?;
///     while h.roll() {
///         let _ = (h.pos(), h.hashes());
///     }
/// }
/// # Ok(()) }
/// ```
pub struct SpacedSeeds {
    k: u16,
    tabs: Vec<SeedTab>,
}

impl SpacedSeeds {
    /// Compile seeds from mask strings of length `k`, where `'1'` marks a care
    /// position and `'0'` a don't-care position.
    ///
    /// # Errors
    /// [`NtHashError::InvalidK`] if `k == 0`,
    /// [`NtHashError::InvalidMaskLength`] if a mask is not `k` characters long,
    /// or [`NtHashError::InvalidMask`] if it contains anything but `'0'`/`'1'`.
    pub fn from_masks<S: AsRef<str>>(masks: &[S], k: u16) -> Result<Self> {
        if k == 0 {
            return Err(NtHashError::InvalidK);
        }
        let tabs = masks
            .iter()
            .map(|m| parse_seed_string(m.as_ref(), k).map(|care| build_seed_tab(care, k as usize)))
            .collect::<Result<_>>()?;
        Ok(Self { k, tabs })
    }

    /// Compile seeds from pre-parsed care indices (one list per seed). The
    /// order of the indices does not matter and duplicates are ignored.
    ///
    /// # Errors
    /// [`NtHashError::InvalidK`] if `k == 0`, or
    /// [`NtHashError::InvalidCareIndex`] if an index is `>= k`.
    pub fn from_care_indices(seeds: Vec<Vec<usize>>, k: u16) -> Result<Self> {
        if k == 0 {
            return Err(NtHashError::InvalidK);
        }
        let tabs = seeds
            .into_iter()
            .map(|mut care| {
                if let Some(&index) = care.iter().find(|&&i| i >= k as usize) {
                    return Err(NtHashError::InvalidCareIndex { index, k });
                }
                // The tables assume a set of positions: a repeated index would
                // be XOR-ed in twice and counted twice as an ambiguous base.
                care.sort_unstable();
                care.dedup();
                Ok(build_seed_tab(care, k as usize))
            })
            .collect::<Result<_>>()?;
        Ok(Self { k, tabs })
    }

    /// The k‑mer span these seeds were compiled for.
    #[inline]
    pub fn k(&self) -> u16 {
        self.k
    }

    /// Number of seeds.
    #[inline]
    pub fn len(&self) -> usize {
        self.tabs.len()
    }

    /// Whether the set contains no seeds.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.tabs.is_empty()
    }
}

/// The seeds a hasher rolls with: compiled for this hasher alone, or borrowed
/// from the caller so the tables are shared across hashers.
enum SeedSet<'a> {
    Owned(SpacedSeeds),
    Shared(&'a SpacedSeeds),
}

impl SeedSet<'_> {
    #[inline(always)]
    fn tabs(&self) -> &[SeedTab] {
        match self {
            SeedSet::Owned(s) => &s.tabs,
            SeedSet::Shared(s) => &s.tabs,
        }
    }
}

/// Rolling state for one seed: its forward/reverse hashes for the current
/// window and how many of its care positions currently hold an ambiguous base.
#[derive(Clone, Copy)]
struct SeedState {
    fwd: u64,
    rev: u64,
    n_invalid: usize,
}

/// Rolling spaced-seed ntHash. Handles multiple seeds and multiple hashes per
/// k-mer, advancing between windows in O(#care-runs) rather than recomputing.
pub struct SeedNtHash<'a> {
    /// Input nucleotide sequence.
    seq: &'a [u8],
    /// k-mer span.
    k: usize,
    /// Number of hashes per seed.
    num_hashes: usize,
    /// Care indices + precomputed roll/recompute tables.
    seeds: SeedSet<'a>,
    /// Per-seed rolling hash state.
    state: Vec<SeedState>,
    /// Current position in the sequence.
    pos: usize,
    /// Hash results, flattened as `seeds × num_hashes`.
    hashes: Vec<u64>,
    /// Whether the starting window has been evaluated yet.
    initialized: bool,
}

impl<'a> SeedNtHash<'a> {
    /// Creates a new hasher from a sequence and spaced-seed masks.
    ///
    /// The masks are compiled for this hasher alone; to hash many sequences
    /// with the same seeds, compile them once with [`SpacedSeeds::from_masks`]
    /// and use [`with_seeds`](Self::with_seeds).
    ///
    /// `num_hashes_per_seed` is clamped to at least 1.
    ///
    /// # Errors
    /// Returns an error if `k` is zero, the sequence is too short, `start_pos`
    /// is out of range, or a mask is invalid.
    pub fn new<S: AsRef<str>>(
        seq: &'a [u8],
        seed_masks: &[S],
        num_hashes_per_seed: usize,
        k: u16,
        start_pos: usize,
    ) -> Result<Self> {
        Self::check_bounds(seq, k, start_pos)?;
        let seeds = SpacedSeeds::from_masks(seed_masks, k)?;
        Ok(Self::build(
            seq,
            SeedSet::Owned(seeds),
            num_hashes_per_seed,
            start_pos,
        ))
    }

    /// Alternative constructor using pre-parsed care indices (skips mask
    /// parsing). The order of the indices does not matter and duplicates are
    /// ignored.
    ///
    /// # Errors
    /// As [`new`](Self::new), with [`NtHashError::InvalidCareIndex`] for an
    /// index `>= k`.
    pub fn from_care_indices(
        seq: &'a [u8],
        seeds: Vec<Vec<usize>>,
        num_hashes_per_seed: usize,
        k: u16,
        start_pos: usize,
    ) -> Result<Self> {
        Self::check_bounds(seq, k, start_pos)?;
        let seeds = SpacedSeeds::from_care_indices(seeds, k)?;
        Ok(Self::build(
            seq,
            SeedSet::Owned(seeds),
            num_hashes_per_seed,
            start_pos,
        ))
    }

    /// Creates a hasher over `seq` that borrows precompiled `seeds`, skipping
    /// the table construction [`new`](Self::new) performs on every call. The
    /// k‑mer span is [`seeds.k()`](SpacedSeeds::k).
    ///
    /// # Errors
    /// Returns an error if the sequence is shorter than `k` or `start_pos` is
    /// out of range.
    pub fn with_seeds(
        seq: &'a [u8],
        seeds: &'a SpacedSeeds,
        num_hashes_per_seed: usize,
        start_pos: usize,
    ) -> Result<Self> {
        Self::check_bounds(seq, seeds.k, start_pos)?;
        Ok(Self::build(
            seq,
            SeedSet::Shared(seeds),
            num_hashes_per_seed,
            start_pos,
        ))
    }

    fn check_bounds(seq: &[u8], k: u16, start_pos: usize) -> Result<()> {
        if k == 0 {
            return Err(NtHashError::InvalidK);
        }
        let k_usz = k as usize;
        if seq.len() < k_usz {
            return Err(NtHashError::SequenceTooShort {
                seq_len: seq.len(),
                k,
            });
        }
        if start_pos > seq.len() - k_usz {
            return Err(NtHashError::PositionOutOfRange {
                pos: start_pos,
                seq_len: seq.len(),
            });
        }
        Ok(())
    }

    /// Assemble a hasher from validated parts.
    fn build(seq: &'a [u8], seeds: SeedSet<'a>, num_hashes: usize, pos: usize) -> Self {
        let num_hashes = num_hashes.max(1);
        let (num_seeds, k) = match &seeds {
            SeedSet::Owned(s) => (s.len(), s.k),
            SeedSet::Shared(s) => (s.len(), s.k),
        };
        Self {
            seq,
            k: k as usize,
            num_hashes,
            seeds,
            state: vec![
                SeedState {
                    fwd: 0,
                    rev: 0,
                    n_invalid: 0,
                };
                num_seeds
            ],
            pos,
            hashes: vec![0; num_seeds * num_hashes],
            initialized: false,
        }
    }

    /// Returns the current position in the sequence.
    #[inline(always)]
    pub fn pos(&self) -> usize {
        self.pos
    }

    /// Returns the current set of hash values, flattened seed by seed
    /// (`num_hashes` values each).
    #[inline(always)]
    pub fn hashes(&self) -> &[u64] {
        &self.hashes
    }

    /// Advances the iterator by one position.
    /// On first call, searches for the first valid k-mer (initialization).
    ///
    /// Windows containing an ambiguous base at a care position are **skipped**
    /// (the position advances to the next valid k-mer), matching the behaviour
    /// of the contiguous [`NtHash`](crate::kmer::NtHash). `false` is returned
    /// only when the end of the sequence is reached.
    pub fn roll(&mut self) -> bool {
        let last = self.seq.len() - self.k;
        if !self.initialized {
            self.initialized = true;
            // Recompute the state for the starting window, then emit it if valid.
            let tabs = self.seeds.tabs();
            for (s, tab) in self.state.iter_mut().zip(tabs) {
                *s = tab.recompute(self.seq, self.pos);
            }
            if self.all_valid() {
                self.emit();
                return true;
            }
            // Otherwise fall through and roll forward to the first valid window.
        }
        while self.pos < last {
            self.roll_step();
            if self.all_valid() {
                self.emit();
                return true;
            }
        }
        false // End of sequence
    }

    /// Whether every seed's current window is free of ambiguous care bases.
    #[inline(always)]
    fn all_valid(&self) -> bool {
        self.state.iter().all(|s| s.n_invalid == 0)
    }

    /// Advance every seed by one position, rolling or recomputing per seed.
    fn roll_step(&mut self) {
        let pos = self.pos;
        let tabs = self.seeds.tabs();
        for (s, tab) in self.state.iter_mut().zip(tabs) {
            *s = if tab.rollable {
                tab.roll(self.seq, pos, *s)
            } else {
                tab.recompute(self.seq, pos + 1)
            };
        }
        self.pos = pos + 1;
    }

    /// Materialize the hash outputs for the current (valid) window.
    fn emit(&mut self) {
        for (s, out) in self
            .state
            .iter()
            .zip(self.hashes.chunks_exact_mut(self.num_hashes))
        {
            extend_hashes(s.fwd, s.rev, self.k as u32, out);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Builder + Iterator façade for ergonomic traversal of spaced-seed hashes
// ─────────────────────────────────────────────────────────────────────────────

/// Builder for creating a `SeedNtHashIter`, providing ergonomic configuration.
///
/// Example:
/// ```rust
/// use nthash_rs::{SeedNtHashBuilder, Result};
///
/// # fn main() -> Result<()> {
/// let seq   = b"ATCGTACGATGCATGCATGCTGACG";
/// let masks = vec!["000111", "010101"];
///
/// for (pos, hashes) in SeedNtHashBuilder::new(seq)
///                        .k(6)
///                        .masks(masks)
///                        .num_hashes(2)
///                        .finish()? {
///     println!("{pos:2}  {:016x}", hashes[0]);
/// }
/// # Ok(()) }
/// ```
pub struct SeedNtHashBuilder<'a> {
    seq: &'a [u8],
    masks: Vec<String>,
    seeds: Option<&'a SpacedSeeds>,
    k: u16,
    num_hashes: usize,
    start_pos: usize,
}

impl<'a> SeedNtHashBuilder<'a> {
    /// Starts building a new ntHash configuration from the given sequence.
    pub fn new(seq: &'a [u8]) -> Self {
        Self {
            seq,
            masks: Vec::new(),
            seeds: None,
            k: 0,
            num_hashes: 1,
            start_pos: 0,
        }
    }

    /// Sets the k-mer size.
    pub fn k(mut self, k: u16) -> Self {
        self.k = k;
        self
    }

    /// Adds seed masks where '1' indicates positions to hash.
    pub fn masks<S: Into<String>, I: IntoIterator<Item = S>>(mut self, m: I) -> Self {
        self.masks = m.into_iter().map(Into::into).collect();
        self
    }

    /// Uses precompiled seeds instead of [`masks`](Self::masks), so the lookup
    /// tables are shared rather than rebuilt for this sequence. Takes precedence
    /// over `masks`; `k` may be omitted, but if set it must equal
    /// [`seeds.k()`](SpacedSeeds::k).
    pub fn seeds(mut self, seeds: &'a SpacedSeeds) -> Self {
        self.seeds = Some(seeds);
        self
    }

    /// Specifies number of hashes per spaced seed.
    pub fn num_hashes(mut self, n: usize) -> Self {
        self.num_hashes = n;
        self
    }

    /// Sets the start position in the sequence.
    pub fn pos(mut self, p: usize) -> Self {
        self.start_pos = p;
        self
    }

    /// Finalizes the builder and returns an iterator over the hashes.
    pub fn finish(self) -> Result<SeedNtHashIter<'a>> {
        let hasher = match self.seeds {
            Some(seeds) => {
                if self.k != 0 && self.k != seeds.k() {
                    return Err(NtHashError::InvalidK);
                }
                SeedNtHash::with_seeds(self.seq, seeds, self.num_hashes, self.start_pos)?
            }
            None => SeedNtHash::new(
                self.seq,
                &self.masks,
                self.num_hashes,
                self.k,
                self.start_pos,
            )?,
        };
        Ok(SeedNtHashIter {
            hasher,
            done: false,
        })
    }
}

/// Iterator for traversing valid k-mers and yielding spaced-seed hashes.
pub struct SeedNtHashIter<'a> {
    hasher: SeedNtHash<'a>,
    done: bool,
}

impl<'a> SeedNtHashIter<'a> {
    /// Drive the hasher to completion, invoking `f(pos, &hashes)` for each
    /// valid k‑mer **without allocating**.
    ///
    /// The [`Iterator`] impl yields an owned `Vec<u64>` per k‑mer (one heap
    /// allocation each); this method lends the internal buffer to `f` instead,
    /// so the traversal allocates nothing. Windows with an ambiguous base at a
    /// care position are skipped, mirroring the [`Iterator`].
    #[inline]
    pub fn for_each_hash<F: FnMut(usize, &[u64])>(mut self, mut f: F) {
        if self.done {
            return;
        }
        while self.hasher.roll() {
            f(self.hasher.pos(), self.hasher.hashes());
        }
    }
}

impl<'a> Iterator for SeedNtHashIter<'a> {
    type Item = (usize, Vec<u64>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        if !self.hasher.roll() {
            self.done = true;
            return None;
        }
        Some((self.hasher.pos(), self.hasher.hashes().to_vec()))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Basic Unit Test
// ─────────────────────────────────────────────────────────────────────────────
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_spaced_seed() {
        let seq = b"ATCGTACGATGCATGCATGCTGACG";
        let mut h = SeedNtHash::new(seq, &["000111", "010101"], 2, 6, 0).unwrap();
        assert!(h.roll()); // first valid
        let first = h.hashes()[0];
        assert!(h.roll()); // next valid
        assert_ne!(first, h.hashes()[0]); // hashes should differ
    }

    /// A mid-sequence `N` at a care position must be **skipped**, not treated as
    /// end-of-sequence. This mirrors the contiguous `NtHash`, which resumes
    /// after the ambiguous window.
    #[test]
    fn skips_ambiguous_windows_and_resumes() {
        // 'N' at index 8; valid 4-mers exist both before and after it.
        let seq = b"ACGTACGTNACGTACGTACGT";
        let k = 4u16;

        // Reference set of valid start positions = what the contiguous hasher emits.
        let expected: Vec<usize> = crate::NtHashBuilder::new(seq)
            .k(k)
            .finish()
            .unwrap()
            .map(|(p, _)| p)
            .collect();

        // All-care mask ("1111") makes each spaced window equal the contiguous one.
        let got: Vec<usize> = SeedNtHashBuilder::new(seq)
            .k(k)
            .masks(vec!["1111".to_string()])
            .finish()
            .unwrap()
            .map(|(p, _)| p)
            .collect();

        assert_eq!(
            got, expected,
            "SeedNtHash must skip N-windows and resume, like NtHash"
        );
    }
}
