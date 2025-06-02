//! # ntHash‑rs
//!
//! An idiomatic, pure‑Rust port of the classic *ntHash* rolling‑hash suite,
//! focused on contiguous k‑mer hashing for DNA sequences.
//!
//! This crate currently provides:
//! - [`kmer::NtHash`]: the canonical contiguous‑k‑mer hasher that skips over
//!   non‑ACGT bases (`N` or other characters).
//!
//! All heavy bit‑twiddling is delegated to low‑level modules (`tables` and
//! `constants`), which mirror the original C++ reference implementation, and
//! helper functionality in `util` for canonicalization and hash extension.
//!
//! ## Example
//!
//! ```rust
//! use nthash_rs::{NtHash, Result};
//!
//! fn main() -> Result<()> {
//!     // Create a new NtHash over "ACGTNACGT", k=4, emit 2 hashes per k‑mer, start at pos=0
//!     let mut hasher = NtHash::new(b"ACGTNACGT", 4, 2, 0)?;
//!
//!     // First call to roll() initializes and returns true if a valid k‑mer was found
//!     assert!(hasher.roll());
//!     // Retrieve the two hash values for the first valid 4‑mer
//!     let hashes = hasher.hashes();
//!     println!("First k‑mer hashes: {:#x}, {:#x}", hashes[0], hashes[1]);
//!
//!     // Advance through the sequence
//!     while hasher.roll() {
//!         let h = hasher.hashes()[0];
//!         println!("Next k‑mer forward hash: {:#x}", h);
//!     }
//!     Ok(())
//! }
//! ```

// Uncomment to build with `no_std` support
// #![cfg_attr(not(feature = "std"), no_std)]

/// Low‑level random seeds, split‑rotate tables, and numeric constants.
// Not re‑exported directly.
mod constants;
mod tables;

pub mod util;
/// High‑level contiguous k‑mer rolling hasher.
/// Skips over non‑ACGT bases exactly as the original reference.
pub mod kmer;
pub mod blind;
pub mod seed;

// ──────────────────────────────────────────────────────────────
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
pub use kmer::NtHashIter;

pub use blind::BlindNtHash;
pub use blind::BlindNtHashBuilder;

pub use seed::SeedNtHash;
pub use seed::SeedNtHashBuilder;

// ──────────────────────────────────────────────────────────────
// Crate‑wide result and error types
// --------------------------------------------------------------------------

/// Shorthand `Result` alias for this crate’s operations.
pub type Result<T, E = NtHashError> = std::result::Result<T, E>;

/// Errors common to all ntHash k‑mer hashers.
#[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
pub enum NtHashError {
    /// `k` was zero.
    #[error("k must be > 0")]
    InvalidK,

    /// Provided sequence length is shorter than `k`.
    #[error("sequence length ({seq_len}) < k ({k})")]
    SequenceTooShort { seq_len: usize, k: u16 },

    /// Starting `pos` is beyond the last valid window (`seq.len() - k`).
    #[error("position ({pos}) exceeds sequence length ({seq_len})")]
    PositionOutOfRange { pos: usize, seq_len: usize },

    #[error("invalid sequence")]
    InvalidSequence,

    #[error("invalid window offsets")]
    InvalidWindowOffsets,
}

// ──────────────────────────────────────────────────────────────
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
