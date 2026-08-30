//! **Streaming spaced-seed ntHash** for *non-contiguous* k‑mers.
//!
//! **`SeedNtHash` computes hashes using spaced seeds**, where only selected
//! positions in the k‑mer are considered (“care sites”).
//!
//! Sliding the window by one changes a care position's membership only at the
//! **run boundaries** of the mask, so a window costs `O(#runs)` rather than
//! `O(#care)` — unless the mask is so scattered that it has more boundaries
//! than care sites, in which case recomputing wins and is chosen per seed.
//!
//! Two things keep the constant small. Contribution tables are indexed by
//! *byte class* (see `ClassMap`) rather than by the raw byte, which is the
//! difference between 256 B and 4 KiB per position — a `1110`-repeating mask
//! over `k = 31` needs sixteen of them, so the byte-indexed form used to spill
//! out of L1. And a single seed — the common case — gets a driver that keeps
//! the rolling state in locals instead of walking three heap `Vec`s per k‑mer.
//!
//! Bit-level operations are delegated to `tables`, `constants`, and
//! `util::extend_hashes` for efficient hash computation.
//!
//! A Rust‑idiomatic **builder + iterator** (`SeedNtHashBuilder` / `SeedNtHashIter`)
//! provides ergonomic traversal over valid k‑mers.

use std::sync::OnceLock;

use crate::{
    constants::{CP_OFF, SEED_N, SEED_TAB},
    tables::{srol, srol_table, sror},
    util::extend_hashes,
    NtHashError, Result,
};

// ─────────────────────────────────────────────────────────────────────────────
// Byte classes
// ─────────────────────────────────────────────────────────────────────────────

/// Upper bound on distinct byte classes (the map builds 15; see [`ClassMap`]).
const NCLASS: usize = 16;

/// Partition of the 256 byte values by everything a spaced seed can observe
/// about a base.
///
/// `MS_TAB_31L` / `MS_TAB_33R` give every byte one of five rows (A/C/G/T/N), so
/// `srol_table(b, d)` depends on `b` only through its row — for *every* `d`. A
/// seed position needs the forward row, the row of the complement `b & CP_OFF`,
/// and whether `SEED_TAB[b]` is `SEED_N`; those three together partition the
/// byte space into 15 classes.
///
/// This shrinks every per-position contribution table from 256 entries (4 KiB)
/// to 16 (256 B). That is the difference between fitting in L1 and not: a
/// `1110`-repeating mask over `k = 31` has 16 run boundaries, which used to
/// mean 64 KiB of separately-allocated tables against a 48 KiB L1D.
struct ClassMap {
    /// Byte → class index.
    class: [u8; 256],
    /// Bit `i` set iff class `i` is ambiguous (`SEED_TAB[b] == SEED_N`).
    n_mask: u16,
    /// A byte representing each class, used to build contribution tables.
    rep: [u8; NCLASS],
}

impl ClassMap {
    fn build() -> Self {
        let mut m = ClassMap {
            class: [0; 256],
            n_mask: 0,
            rep: [0; NCLASS],
        };
        // `srol_table(_, 0)` is the row's seed, and the five seeds are
        // distinct, so the triple below identifies the class exactly.
        let mut keys: Vec<(u64, u64, bool)> = Vec::new();
        for (b, &seed) in SEED_TAB.iter().enumerate() {
            let c = b as u8;
            let key = (srol_table(c, 0), srol_table(c & CP_OFF, 0), seed == SEED_N);
            let idx = keys.iter().position(|k| *k == key).unwrap_or_else(|| {
                keys.push(key);
                keys.len() - 1
            });
            assert!(idx < NCLASS, "byte classes outgrew NCLASS");
            m.class[b] = idx as u8;
            m.rep[idx] = c;
            if key.2 {
                m.n_mask |= 1 << idx;
            }
        }
        m
    }
}

fn class_map() -> &'static ClassMap {
    static MAP: OnceLock<ClassMap> = OnceLock::new();
    MAP.get_or_init(ClassMap::build)
}

/// Parses a spaced-seed mask string composed of '0' and '1' characters
/// into a list of indices indicating which positions should be used ("care positions").
///
/// # Errors
/// Returns an error if the mask length does not match `k`, or contains characters other than '0' or '1'.
fn parse_seed_string(mask: &str, k: usize) -> Result<Vec<usize>> {
    if mask.len() != k {
        return Err(NtHashError::InvalidK);
    }
    if !mask.bytes().all(|b| b == b'0' || b == b'1') {
        return Err(NtHashError::InvalidSequence);
    }
    Ok(mask
        .bytes()
        .enumerate()
        .filter_map(|(i, b)| if b == b'1' { Some(i) } else { None })
        .collect())
}

/// A class-indexed `[fwd, rev]` contribution table for one care position, used
/// at a run boundary during rolling. Indexed by [`ClassMap::class`], so the
/// table is 16 entries and lives inline in its owner rather than behind its own
/// allocation.
struct BoundaryTab {
    off: usize,
    tab: [[u64; 2]; NCLASS],
}

/// A single spaced seed: the care positions plus everything needed to both
/// recompute a window from scratch and *roll* between windows.
///
/// - `rot` — flattened `care.len() × NCLASS` table `[fwd, rev]` for recomputing a
///   window (used at init and after every roll would be too slow — it is the
///   O(care) path).
/// - `left` / `right` — the **run boundaries** of the care pattern: the only
///   positions whose care-membership changes when the window slides by one.
///   Rolling is `srol(fwd)` / `sror(rev)` plus a correction at each boundary, so
///   the per-window cost is O(#runs) instead of O(#care).
struct SeedTab {
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

/// Build a class-indexed `[srol_table(base, d_fwd), srol_table(base & CP, d_rev)]`
/// contribution table for one boundary position.
fn boundary_tab(off: usize, d_fwd: u32, d_rev: u32) -> BoundaryTab {
    let cm = class_map();
    let mut tab = [[0u64; 2]; NCLASS];
    for (cls, slot) in tab.iter_mut().enumerate() {
        let r = cm.rep[cls];
        *slot = [srol_table(r, d_fwd), srol_table(r & CP_OFF, d_rev)];
    }
    BoundaryTab { off, tab }
}

/// Precompute the recompute and roll tables for one seed's care positions.
fn build_seed_tab(care: Vec<usize>, k: usize) -> SeedTab {
    let cm = class_map();
    let mut rot = vec![[0u64; 2]; care.len() * NCLASS];
    for (i, &p) in care.iter().enumerate() {
        let d_fwd = (k - 1 - p) as u32;
        let d_rev = p as u32;
        for cls in 0..NCLASS {
            let r = cm.rep[cls];
            rot[i * NCLASS + cls] = [srol_table(r, d_fwd), srol_table(r & CP_OFF, d_rev)];
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

/// Recompute one seed's `(fwd, rev, n_invalid)` for the window at `pos` — the
/// O(#care) path, shared by the generic driver and the hoisted single-seed loop.
#[inline]
fn recompute_one(seq: &[u8], st: &SeedTab, pos: usize) -> (u64, u64, usize) {
    let cm = class_map();
    let mut fwd = 0u64;
    let mut rev = 0u64;
    let mut n_invalid = 0usize;
    for (i, &p) in st.care.iter().enumerate() {
        let cls = cm.class[seq[pos + p] as usize] as usize;
        let e = st.rot[i * NCLASS + cls];
        fwd ^= e[0];
        rev ^= e[1];
        // Branchless: the ambiguity of a class is one bit of `n_mask`.
        n_invalid += ((cm.n_mask >> cls) & 1) as usize;
    }
    (fwd, rev, n_invalid)
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
    seq: &'a [u8],         // Input nucleotide sequence
    k: usize,              // k-mer size
    num_hashes: usize,     // Number of hashes per seed
    seeds: Vec<SeedTab>,   // Care indices + precomputed roll/recompute tables
    state: Vec<SeedState>, // Per-seed rolling hash state
    pos: usize,            // Current position in the sequence
    hashes: Vec<u64>,      // Hash results (flattened)
    initialised: bool,     // Whether the hasher has found the first valid k-mer
}

impl<'a> SeedNtHash<'a> {
    /// Creates a new hasher from a sequence and spaced-seed masks.
    ///
    /// # Errors
    /// Returns an error if `k` is zero, the sequence is too short, or a mask is invalid.
    pub fn new(
        seq: &'a [u8],
        seed_masks: &[String],
        num_hashes_per_seed: usize,
        k: u16,
        start_pos: usize,
    ) -> Result<Self> {
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

        let mut seeds = Vec::with_capacity(seed_masks.len());
        for m in seed_masks {
            let care = parse_seed_string(m, k_usz)?;
            seeds.push(build_seed_tab(care, k_usz));
        }

        let state = vec![
            SeedState {
                fwd: 0,
                rev: 0,
                n_invalid: 0,
            };
            seeds.len()
        ];
        Ok(Self {
            seq,
            k: k_usz,
            num_hashes: num_hashes_per_seed.max(1),
            seeds,
            state,
            pos: start_pos,
            hashes: vec![0; seed_masks.len() * num_hashes_per_seed.max(1)],
            initialised: false,
        })
    }

    /// Alternative constructor using pre-parsed care indices (skips mask parsing).
    pub fn from_care_indices(
        seq: &'a [u8],
        seeds: Vec<Vec<usize>>,
        num_hashes_per_seed: usize,
        k: u16,
        start_pos: usize,
    ) -> Result<Self> {
        let k_usz = k as usize;
        if seeds.iter().any(|v| v.iter().any(|&i| i >= k_usz)) {
            return Err(NtHashError::InvalidWindowOffsets);
        }
        Self::new(
            seq,
            &vec![String::from_utf8(vec![b'0'; k_usz]).unwrap(); seeds.len()], // dummy masks
            num_hashes_per_seed,
            k,
            start_pos,
        )
        .map(|mut s| {
            // Replace the dummy (empty-care) seeds with the caller's care indices,
            // rebuilding the precomputed rotation tables for each.
            s.seeds = seeds
                .into_iter()
                .map(|care| build_seed_tab(care, k_usz))
                .collect();
            s
        })
    }

    /// Returns the current position in the sequence.
    #[inline(always)]
    pub fn pos(&self) -> usize {
        self.pos
    }

    /// Returns the current set of hash values.
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
        if !self.initialised {
            self.initialised = true;
            // Recompute the state for the starting window, then emit it if valid.
            self.compute_state_at(self.pos);
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

    /// Recompute one seed's state from scratch for the window at `pos` — the
    /// O(#care) path, used at init and for seeds where rolling doesn't pay.
    fn recompute_seed(&mut self, si: usize, pos: usize) {
        let (fwd, rev, n_invalid) = recompute_one(self.seq, &self.seeds[si], pos);
        self.state[si] = SeedState {
            fwd,
            rev,
            n_invalid,
        };
    }

    /// Roll one seed's state from the window at `pos` to `pos + 1`, correcting
    /// only at the care pattern's run boundaries — the O(#boundaries) path.
    fn roll_seed(&mut self, si: usize, pos: usize) {
        let cm = class_map();
        let st = &self.seeds[si];
        let s = &mut self.state[si];
        let mut fwd = srol(s.fwd);
        let mut rev = sror(s.rev);
        let mut n_invalid = s.n_invalid;
        // Bases leaving the care coverage.
        for lb in &st.left {
            let cls = cm.class[self.seq[pos + lb.off] as usize] as usize;
            let e = lb.tab[cls];
            fwd ^= e[0];
            rev ^= e[1];
            n_invalid -= ((cm.n_mask >> cls) & 1) as usize;
        }
        // Bases entering the care coverage.
        for rb in &st.right {
            let cls = cm.class[self.seq[pos + 1 + rb.off] as usize] as usize;
            let e = rb.tab[cls];
            fwd ^= e[0];
            rev ^= e[1];
            n_invalid += ((cm.n_mask >> cls) & 1) as usize;
        }
        *s = SeedState {
            fwd,
            rev,
            n_invalid,
        };
    }

    /// Recompute every seed's state from scratch for the window at `pos`.
    fn compute_state_at(&mut self, pos: usize) {
        for si in 0..self.seeds.len() {
            self.recompute_seed(si, pos);
        }
    }

    /// Advance every seed by one position, rolling or recomputing per seed.
    fn roll_step(&mut self) {
        let pos = self.pos;
        for si in 0..self.seeds.len() {
            if self.seeds[si].rollable {
                self.roll_seed(si, pos);
            } else {
                self.recompute_seed(si, pos + 1);
            }
        }
        self.pos = pos + 1;
    }

    /// Materialize the hash outputs for the current (valid) window.
    fn emit(&mut self) {
        for si in 0..self.seeds.len() {
            let s = self.state[si];
            let slice = &mut self.hashes[si * self.num_hashes..(si + 1) * self.num_hashes];
            extend_hashes(s.fwd, s.rev, self.k as u32, slice);
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
        let hasher = SeedNtHash::new(
            self.seq,
            &self.masks,
            self.num_hashes,
            self.k,
            self.start_pos,
        )?;
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
    /// so the traversal allocates nothing. Traversal stops at the first window
    /// with an ambiguous base at a care position, mirroring the [`Iterator`].
    #[inline]
    pub fn for_each_hash<F: FnMut(usize, &[u64])>(mut self, mut f: F) {
        if self.done {
            return;
        }
        // One seed is the common case, and it is the one worth specialising:
        // the generic driver walks three separate heap `Vec`s (seeds, state,
        // hashes) for every k-mer, so the rolling state never reaches a
        // register. The hoisted loop keeps it in locals.
        if self.hasher.seeds.len() == 1 {
            self.for_each_hash_single(f);
            return;
        }
        while self.hasher.roll() {
            f(self.hasher.pos(), self.hasher.hashes());
        }
    }

    /// [`for_each_hash`](Self::for_each_hash) specialised to a single seed.
    fn for_each_hash_single<F: FnMut(usize, &[u64])>(mut self, mut f: F) {
        let SeedNtHash {
            seq,
            k,
            seeds,
            state,
            pos,
            hashes,
            initialised,
            ..
        } = &mut self.hasher;
        let seq: &[u8] = seq;
        let k = *k;
        let last = seq.len() - k;
        let st = &seeds[0];
        let cm = class_map();
        let class = &cm.class;
        let n_mask = cm.n_mask;

        let mut p = *pos;
        let (mut fwd, mut rev, mut n_inv);
        if *initialised {
            let s = state[0];
            fwd = s.fwd;
            rev = s.rev;
            n_inv = s.n_invalid;
        } else {
            *initialised = true;
            (fwd, rev, n_inv) = recompute_one(seq, st, p);
            if n_inv == 0 {
                extend_hashes(fwd, rev, k as u32, hashes);
                f(p, hashes);
            }
        }

        if st.rollable {
            let left = &st.left[..];
            let right = &st.right[..];
            while p < last {
                fwd = srol(fwd);
                rev = sror(rev);
                // Bases leaving the care coverage.
                for lb in left {
                    let cls = class[seq[p + lb.off] as usize] as usize;
                    let e = lb.tab[cls];
                    fwd ^= e[0];
                    rev ^= e[1];
                    n_inv -= ((n_mask >> cls) & 1) as usize;
                }
                // Bases entering the care coverage.
                for rb in right {
                    let cls = class[seq[p + 1 + rb.off] as usize] as usize;
                    let e = rb.tab[cls];
                    fwd ^= e[0];
                    rev ^= e[1];
                    n_inv += ((n_mask >> cls) & 1) as usize;
                }
                p += 1;
                if n_inv == 0 {
                    extend_hashes(fwd, rev, k as u32, hashes);
                    f(p, hashes);
                }
            }
        } else {
            while p < last {
                p += 1;
                (fwd, rev, n_inv) = recompute_one(seq, st, p);
                if n_inv == 0 {
                    extend_hashes(fwd, rev, k as u32, hashes);
                    f(p, hashes);
                }
            }
        }

        *pos = p;
        state[0] = SeedState {
            fwd,
            rev,
            n_invalid: n_inv,
        };
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

impl<'a> IntoIterator for SeedNtHashBuilder<'a> {
    type Item = (usize, Vec<u64>);
    type IntoIter = SeedNtHashIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.finish()
            .expect("invalid SeedNtHashBuilder configuration")
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Basic Unit Test
// ─────────────────────────────────────────────────────────────────────────────
#[cfg(test)]
mod tests {
    /// The hoisted single-seed driver must select exactly the same windows and
    /// hashes as the generic pull-based one, including across `N`-skips and for
    /// masks that take the recompute path rather than rolling.
    #[test]
    fn single_seed_driver_matches_iterator() {
        let seqs: [&[u8]; 4] = [
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            b"ACGTNNNNACGTACGTNACGGATTACAGGGCCCTTTAAACGT",
            b"NNNNNNNNNNNNNNNN",
            b"TACGGATCCGATTACGATTACAGGGCCCTTTAAACGTACGTA",
        ];
        for k in [4u16, 8, 12] {
            let ku = k as usize;
            let masks = [
                "1".repeat(ku),
                (0..ku)
                    .map(|i| if i % 4 == 3 { '0' } else { '1' })
                    .collect(),
                (0..ku)
                    .map(|i| if i % 3 == 0 { '1' } else { '0' })
                    .collect(),
                format!("1{}1", "0".repeat(ku - 2)),
            ];
            for mask in &masks {
                for m in [1usize, 3] {
                    for seq in seqs {
                        if seq.len() < ku {
                            continue;
                        }
                        let pulled: Vec<(usize, Vec<u64>)> = SeedNtHashBuilder::new(seq)
                            .k(k)
                            .masks(vec![mask.clone()])
                            .num_hashes(m)
                            .finish()
                            .unwrap()
                            .collect();
                        let mut pushed = Vec::new();
                        SeedNtHashBuilder::new(seq)
                            .k(k)
                            .masks(vec![mask.clone()])
                            .num_hashes(m)
                            .finish()
                            .unwrap()
                            .for_each_hash(|p, h| pushed.push((p, h.to_vec())));
                        assert_eq!(pulled, pushed, "k={k} mask={mask} m={m}");
                    }
                }
            }
        }
    }

    /// The class partition must be exact: two bytes may share a class only if
    /// they are indistinguishable to a seed position at *every* rotation.
    #[test]
    fn class_map_is_exact() {
        let cm = class_map();
        // `srol_table` has period lcm(31, 33) = 1023 in its rotation argument.
        for (b, &seed) in SEED_TAB.iter().enumerate() {
            let c = b as u8;
            let r = cm.rep[cm.class[b] as usize];
            for d in 0..1023u32 {
                assert_eq!(srol_table(c, d), srol_table(r, d), "fwd b={b} d={d}");
                assert_eq!(
                    srol_table(c & CP_OFF, d),
                    srol_table(r & CP_OFF, d),
                    "rev b={b} d={d}"
                );
            }
            let is_n = (cm.n_mask >> cm.class[b]) & 1 == 1;
            assert_eq!(is_n, seed == SEED_N, "n-flag b={b}");
        }
    }

    use super::*;

    #[test]
    fn basic_spaced_seed() {
        let seq = b"ATCGTACGATGCATGCATGCTGACG";
        let masks = vec!["000111".to_string(), "010101".to_string()];
        let mut h = SeedNtHash::new(seq, &masks, 2, 6, 0).unwrap();
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
