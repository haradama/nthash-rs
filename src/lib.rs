//! # ntHash‑rs
//!
//! An idiomatic, pure‑Rust port of the *ntHash* rolling‑hash suite for DNA
//! sequences, bit‑compatible with the C++ reference
//! ([bcgsc/ntHash](https://github.com/bcgsc/ntHash)): the test suite compares
//! [`NtHash`], [`BlindNtHash`] and [`SeedNtHash`] against vectors generated
//! from the C++ library. Inputs on which the reference has no well-defined
//! output — `k < 3`, or spaced seeds over a sequence containing `N` — are the
//! exception; see `tests/README.md` in the repository.
//!
//! | Type | Use it for |
//! |---|---|
//! | [`NtHash`] / [`NtHashBuilder`] | contiguous k‑mers; windows containing `N` (or any other non‑ACGTU byte) are skipped |
//! | [`BlindNtHash`] / [`BlindNtHashBuilder`] | streaming over pre‑cleaned input, the caller feeds each base |
//! | [`SeedNtHash`] / [`SeedNtHashBuilder`] / [`SpacedSeeds`] | spaced seeds (one or more `0`/`1` masks) |
//! | [`MinimizerBuilder`] / [`Minimizers`] | canonical (w, k)-minimizers |
//! | [`ReadHasher`] | many N‑free reads with one reusable hasher |
//!
//! All heavy bit‑twiddling is delegated to the low‑level `tables` and
//! `constants` modules, which mirror the reference implementation, with
//! [`util`] providing canonicalization and hash extension.
//!
//! ## Example
//!
//! ```rust
//! use nthash_rs::{NtHash, NtHashBuilder, Result};
//!
//! fn main() -> Result<()> {
//!     // Low-level: roll by hand. k=4, two hashes per k‑mer, start at pos=0.
//!     let mut hasher = NtHash::new(b"ACGTNACGT", 4, 2, 0)?;
//!     // Each roll() moves to the next valid k‑mer (the first call finds the
//!     // first one) and returns false at the end of the sequence.
//!     while hasher.roll() {
//!         let hashes = hasher.hashes();
//!         println!("{}: {:#x}, {:#x}", hasher.pos(), hashes[0], hashes[1]);
//!     }
//!
//!     // High-level: one canonical hash per k‑mer, no allocation.
//!     for (pos, hash) in NtHashBuilder::new(b"ACGTNACGT").k(4).finish()?.canonical() {
//!         println!("{pos}: {hash:#x}");
//!     }
//!     Ok(())
//! }
//! ```

/// Low‑level random seeds, split‑rotate tables, and numeric constants.
// Not re‑exported directly.
mod constants;
mod tables;

/// Batch canonical hashing of many reads with one reusable hasher.
pub mod batch;
/// Streaming hasher over pre‑cleaned input, rolled one caller‑fed base at a time.
pub mod blind;
/// High‑level contiguous k‑mer rolling hasher.
/// Skips over non‑ACGT bases exactly as the original reference.
pub mod kmer;
/// Canonical (w, k)-minimizers built on the contiguous rolling hasher.
pub mod minimizer;
/// Spaced‑seed hasher for non‑contiguous k‑mers.
pub mod seed;
pub mod util;

// README examples are compiled and run as doctests so they cannot rot.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeDoctests;

// --------------------------------------------------------------------------
// Re‑exports: public API surface
// --------------------------------------------------------------------------

/// One‑bit split‑rotate left (33 + 31 halves).
pub use tables::srol;
/// Arbitrary split‑rotate via lookup tables.
pub use tables::srol_table;
/// One‑bit split‑rotate right (33 + 31 halves).
pub use tables::sror;

/// Combine forward and reverse hashes into a strand‑independent value.
pub use util::canonical;
/// Derive multiple hash values from a single canonical hash.
pub use util::extend_hashes;

/// Primary rolling k‑mer hasher.
///
/// See [`kmer::NtHash`] for full documentation.
pub use kmer::NtHash;
pub use kmer::NtHashBuilder;
pub use kmer::NtHashCanonicalIter;
pub use kmer::NtHashIter;

pub use blind::BlindNtHash;
pub use blind::BlindNtHashBuilder;
pub use blind::BlindNtHashCanonicalIter;
pub use blind::BlindNtHashIter;

pub use seed::SeedNtHash;
pub use seed::SeedNtHashBuilder;
pub use seed::SpacedSeeds;

pub use minimizer::Minimizer;
pub use minimizer::MinimizerBuilder;
pub use minimizer::Minimizers;

pub use batch::ReadHasher;

// --------------------------------------------------------------------------
// Crate‑wide result and error types
// --------------------------------------------------------------------------

/// Shorthand `Result` alias for this crate's operations.
pub type Result<T, E = NtHashError> = std::result::Result<T, E>;

/// Errors common to all ntHash k‑mer hashers.
#[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum NtHashError {
    /// `k` was zero (or disagreed with the `k` of precompiled seeds).
    #[error("k must be > 0")]
    InvalidK,

    /// Minimizer window `w` was zero.
    #[error("w must be > 0")]
    InvalidW,

    /// Provided sequence length is shorter than `k`.
    #[error("sequence length ({seq_len}) < k ({k})")]
    SequenceTooShort { seq_len: usize, k: u16 },

    /// Starting `pos` is beyond the last valid window (`seq.len() - k`).
    #[error("position ({pos}) exceeds the last k-mer start of a sequence of length {seq_len}")]
    PositionOutOfRange { pos: usize, seq_len: usize },

    /// A spaced-seed mask is not exactly `k` characters long.
    #[error("spaced-seed mask length ({mask_len}) != k ({k})")]
    InvalidMaskLength { mask_len: usize, k: u16 },

    /// A spaced-seed mask contains a character other than `'0'` or `'1'`.
    #[error("spaced-seed mask may only contain '0' and '1'")]
    InvalidMask,

    /// A spaced-seed care index does not fall inside the k‑mer.
    #[error("care index ({index}) must be < k ({k})")]
    InvalidCareIndex { index: usize, k: u16 },
}

// --------------------------------------------------------------------------
// Basic smoke tests
// --------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sanity_kmer() {
        // Create hasher over "ACGTACGT", k=4, 1 hash per k‑mer, start at 0
        let mut h = NtHash::new("ACGTACGT".as_bytes(), 4, 1, 0).unwrap();
        // First valid k‑mer should be produced
        assert!(h.roll());
        assert_eq!(h.hashes().len(), 1);
    }
}
