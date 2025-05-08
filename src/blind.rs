//! Blind‐mode k‑mer hashing — **BlindNtHash**
//!
//! A "blind" hasher never skips k‑mers that contain non‑ACGT characters;
//! instead such characters behave like an `'N'` (`SEED_N == 0`).  
//!
//! Compared with the canonical [`NtHash`](crate::NtHash):
//! * the caller explicitly provides the next / previous base when rolling
//!   (`roll()` / `roll_back()`), matching the C++ API;
//! * `peek()` / `peek_back()` return the would‑be hash array without mutating
//!   any internal state.
//!
//! **初期計算のみ**は contiguous ntHash のまま呼び出して完全一致させ、
//! その後は差分更新を適用します。

use crate::{
    constants::*, kmer::{base_forward_hash, base_reverse_hash}, tables::{srol, srol_n, srol_table, sror}, util::extend_hashes
};

/// Construction / runtime errors
#[derive(Debug, thiserror::Error)]
pub enum BlindNtHashError {
    #[error("k must be > 0")]
    ZeroK,
    #[error("sequence length ({seq_len}) < pos({pos}) + k({k})")]
    SequenceTooShort {
        seq_len: usize,
        pos: usize,
        k: usize,
    },
}

/// Convenience alias
pub type Result<T> = core::result::Result<T, BlindNtHashError>;

/// Blind version of NtHash
///
/// After construction we store a ring‑buffer (`buf` + `head`) of length `k`.
/// `roll()` / `roll_back()` は差分更新だけ、初回だけ contigous ロジックを呼びます。
pub struct BlindNtHash {
    buf: Vec<u8>,
    head: usize,
    k: u16,
    num_hashes: u8,
    fwd_hash: u64,
    rev_hash: u64,
    hashes: Vec<u64>,
    pos: isize,
}

impl BlindNtHash {
    #[inline(always)]
    pub fn hashes(&self) -> &[u64] {
        &self.hashes
    }
    #[inline(always)]
    pub fn forward_hash(&self) -> u64 {
        self.fwd_hash
    }
    #[inline(always)]
    pub fn reverse_hash(&self) -> u64 {
        self.rev_hash
    }
    #[inline(always)]
    pub fn k(&self) -> u16 {
        self.k
    }
    #[inline(always)]
    pub fn pos(&self) -> isize {
        self.pos
    }

    /// Build a new *blind* NtHash centred on `seq[pos..pos+k]`.
    pub fn new(seq: &str, k: u16, num_hashes: u8, pos: usize) -> Result<Self> {
        if k == 0 {
            return Err(BlindNtHashError::ZeroK);
        }
        let ksz = k as usize;
        if seq.len() < pos + ksz {
            return Err(BlindNtHashError::SequenceTooShort {
                seq_len: seq.len(),
                pos,
                k: ksz,
            });
        }

        // --- リングバッファに k-mer を格納 ---
        let mut buf = Vec::with_capacity(ksz);
        buf.extend_from_slice(&seq.as_bytes()[pos..pos + ksz]);
        let head = 0;

        let fwd_hash = base_forward_hash_buf(&buf, head, k);
        let rev_hash = base_reverse_hash_buf(&buf, head, k);

        let mut hashes = vec![0u64; num_hashes as usize];
        extend_hashes(fwd_hash, rev_hash, k as u32, &mut hashes);

        Ok(Self {
            buf,
            head,
            k,
            num_hashes,
            fwd_hash,
            rev_hash,
            hashes,
            pos: pos as isize,
        })
    }

    /// Roll the window forward by one base.
    pub fn roll(&mut self, char_in: u8) {
        let idx = self.head;
        let char_out = self.buf[idx];
        self.buf[idx] = char_in;
        self.head = if idx + 1 == self.k as usize {
            0
        } else {
            idx + 1
        };

        self.fwd_hash = next_forward_hash(self.fwd_hash, self.k, char_out, char_in);
        self.rev_hash = next_reverse_hash(self.rev_hash, self.k, char_out, char_in);
        extend_hashes(
            self.fwd_hash,
            self.rev_hash,
            self.k as u32,
            &mut self.hashes,
        );

        self.pos += 1;
    }

    /// Roll the window backward by one base.
    pub fn roll_back(&mut self, char_in: u8) {
        let back = if self.head == 0 {
            self.k as usize - 1
        } else {
            self.head - 1
        };
        let char_out = self.buf[back];
        self.buf[back] = char_in;
        self.head = back;

        self.fwd_hash = prev_forward_hash(self.fwd_hash, self.k, char_out, char_in);
        self.rev_hash = prev_reverse_hash(self.rev_hash, self.k, char_out, char_in);
        extend_hashes(
            self.fwd_hash,
            self.rev_hash,
            self.k as u32,
            &mut self.hashes,
        );

        self.pos -= 1;
    }

    /// Preview forward roll
    #[inline]
    pub fn peek<'a>(&self, char_in: u8, out: &'a mut [u64]) -> &'a [u64] {
        let char_out = self.buf[self.head];
        let fwd = next_forward_hash(self.fwd_hash, self.k, char_out, char_in);
        let rev = next_reverse_hash(self.rev_hash, self.k, char_out, char_in);
        extend_hashes(fwd, rev, self.k as u32, out);
        &out[..]
    }

    /// Preview backward roll
    #[inline]
    pub fn peek_back<'a>(&self, char_in: u8, out: &'a mut [u64]) -> &'a [u64] {
        let back = if self.head == 0 {
            self.k as usize - 1
        } else {
            self.head - 1
        };
        let char_out = self.buf[back];
        let fwd = prev_forward_hash(self.fwd_hash, self.k, char_out, char_in);
        let rev = prev_reverse_hash(self.rev_hash, self.k, char_out, char_in);
        extend_hashes(fwd, rev, self.k as u32, out);
        &out[..]
    }
}

// ──────────────────────────────────────────────────────────────
// 差分更新のヘルパーはそのまま
// -------------------------------------------------------------------------
#[inline(always)]
fn next_forward_hash(prev: u64, k: u16, char_out: u8, char_in: u8) -> u64 {
    let mut h = srol(prev);
    h ^= SEED_TAB[char_in as usize];
    h ^= srol_table(char_out, k as u32);
    h
}
#[inline(always)]
fn prev_forward_hash(prev: u64, k: u16, char_out: u8, char_in: u8) -> u64 {
    let mut h = prev ^ srol_table(char_in, k as u32);
    h ^= SEED_TAB[char_out as usize];
    sror(h)
}
#[inline(always)]
pub fn next_reverse_hash(prev: u64, k: u16, char_out: u8, char_in: u8) -> u64 {
    let mut h = prev ^ srol_table(char_in & CP_OFF, k as u32);
    h ^= SEED_TAB[(char_out & CP_OFF) as usize];
    sror(h)
}
#[inline(always)]
pub fn prev_reverse_hash(prev: u64, k: u16, char_out: u8, char_in: u8) -> u64 {
    let mut h = srol(prev);
    h ^= SEED_TAB[(char_in & CP_OFF) as usize];
    h ^= srol_table(char_out & CP_OFF, k as u32);
    h
}

//――――――――――――――――――――――――――――――――――――――――――――――
// BlindNtHash 初期計算用: バッファ全体を 1bit split‑rotate + XOR(SEED_TAB)
#[inline(always)]
fn base_forward_hash_buf(buf: &[u8], head: usize, k: u16) -> u64 {
    let ksz = k as usize;
    let mut h = 0u64;
    // one‐bit split‑rotate + XOR(SEED_TAB) per base
    for i in 0..ksz {
        h = srol(h);
        // SEED_TAB maps A,C,G,T→seed_A/C/G/T, all other ASCII→0
        let b = buf[(head + i) % ksz];
        h ^= SEED_TAB[b as usize];
    }
    h
}

#[inline(always)]
fn base_reverse_hash_buf(buf: &[u8], head: usize, k: u16) -> u64 {
    let ksz = k as usize;
    let mut h = 0u64;
    // same idea, but read from the right and mask to the complement
    for i in 0..ksz {
        h = srol(h);
        let raw = buf[(head + ksz - 1 - i) % ksz] & CP_OFF;
        // SEED_TAB[raw] is SEED_T/A/C/G when raw is one of those,
        // or zero if the original byte was not A/C/G/T.
        h ^= SEED_TAB[raw as usize];
    }
    h
}

