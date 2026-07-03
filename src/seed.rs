//! **Streaming spaced-seed ntHash** for *non-contiguous* k‑mers.
//!
//! **`SeedNtHash` computes hashes using spaced seeds**, where only selected
//! positions in the k‑mer are considered (“care sites”).
//!
//! Hashes are re‑computed per window rather than rolled, allowing support
//! for multiple seeds and arbitrary binary masks.
//!
//! Bit-level operations are delegated to `tables`, `constants`, and
//! `util::extend_hashes` for efficient hash computation.
//!
//! A Rust‑idiomatic **builder + iterator** (`SeedNtHashBuilder` / `SeedNtHashIter`)
//! provides ergonomic traversal over valid k‑mers.

use crate::{
    constants::{CP_OFF, SEED_N, SEED_TAB},
    tables::srol_table,
    util::extend_hashes,
    NtHashError, Result,
};

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

/// A single spaced seed with its per-care-position rotation tables precomputed.
///
/// For a fixed seed and `k`, the split-rotate distance at each care position is
/// constant (`k-1-p` forward, `p` reverse); only the base at that position
/// varies. We therefore precompute `srol_table` for every byte value once, so
/// the per-window inner loop is a single array load instead of two `srol_table`
/// calls (each a `% 31` / `% 33` reduction plus a `MS_TAB_*` pointer-chase).
/// `rot` is a flattened `care.len() × 256` table indexed as
/// `[care_idx * 256 + byte]`; each entry packs the forward rotation in `[0]` and
/// the reverse rotation in `[1]`, so both are fetched from one cache line.
struct SeedTab {
    care: Vec<usize>,
    rot: Vec<[u64; 2]>,
}

/// Precompute the rotation tables for one seed's care positions.
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
    SeedTab { care, rot }
}

/// Struct for computing spaced-seed ntHash values in a re-computational manner.
/// Can handle multiple seeds and generates multiple hashes per k-mer.
pub struct SeedNtHash<'a> {
    seq: &'a [u8],       // Input nucleotide sequence
    k: usize,            // k-mer size
    num_hashes: usize,   // Number of hashes per seed
    seeds: Vec<SeedTab>, // Care indices + precomputed rotation tables per seed
    pos: usize,          // Current position in the sequence
    hashes: Vec<u64>,    // Hash results (flattened)
    initialised: bool,   // Whether the hasher has found the first valid k-mer
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

        Ok(Self {
            seq,
            k: k_usz,
            num_hashes: num_hashes_per_seed.max(1),
            seeds,
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
        if !self.initialised {
            return self.init();
        }

        while self.pos < self.seq.len() - self.k {
            self.pos += 1;
            if self.compute_current() {
                return true;
            }
        }
        false // End of sequence
    }

    /// Computes hashes for the k-mer at the current position.
    /// Returns false if any ambiguous base is found.
    fn compute_current(&mut self) -> bool {
        let win = &self.seq[self.pos..self.pos + self.k];
        for st in &self.seeds {
            if st.care.iter().any(|&p| SEED_TAB[win[p] as usize] == SEED_N) {
                return false;
            }
        }

        for (i_seed, st) in self.seeds.iter().enumerate() {
            let mut fwd = 0u64;
            let mut rev = 0u64;
            for (i, &p) in st.care.iter().enumerate() {
                let e = st.rot[i * 256 + win[p] as usize];
                fwd ^= e[0];
                rev ^= e[1];
            }
            let slice = &mut self.hashes[i_seed * self.num_hashes..(i_seed + 1) * self.num_hashes];
            extend_hashes(fwd, rev, self.k as u32, slice);
        }
        true
    }

    /// Initializes by finding the first valid k-mer in the sequence.
    fn init(&mut self) -> bool {
        while self.pos <= self.seq.len() - self.k {
            if self.compute_current() {
                self.initialised = true;
                return true;
            }
            self.pos += 1;
        }
        false
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
