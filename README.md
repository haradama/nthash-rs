# nthsash‑rs

[<img alt="github" src="https://img.shields.io/badge/github-haradama/nthash__rs-8da0cb?style=for-the-badge&labelColor=555555&logo=github" height="20">](https://github.com/haradama/nthash-rs)
[<img alt="crates.io" src="https://img.shields.io/crates/v/nthash-rs.svg?style=for-the-badge&color=fc8d62&logo=rust" height="20">](https://crates.io/crates/nthash-rs)
[<img alt="docs.rs" src="https://img.shields.io/badge/docs.rs-nthash__rs-66c2a5?style=for-the-badge&labelColor=555555&logo=docs.rs" height="20">](https://docs.rs/nthash-rs)
[<img alt="build status" src="https://img.shields.io/github/actions/workflow/status/haradama/nthash-rs/rust.yml?branch=master&style=for-the-badge" height="20">](https://github.com/haradama/nthash-rs/actions)

Pure‑Rust port of [ntHash](https://github.com/bcgsc/ntHash) rolling‑hash suite, focused on contiguous k‑mer hashing for DNA sequences.

## Installation

```shell
cargo add nthash-rs
```

## Quick Start

```rust
use nthash_rs::{NtHashBuilder, NtHashError};

fn main() -> Result<(), NtHashError> {
    let seq = b"ACGTCAGTNNNNACGTACGT";
    let k = 4u16;
    let m = 2u8; // number of hashes per k-mer

    // Build an iterator over all valid k-mers
    let iter = NtHashBuilder::new(seq)
        .k(k)
        .num_hashes(m)
        .pos(0)
        .finish()?;

    for (pos, hashes) in iter {
        // slice out the current k-mer
        let kmer = &seq[pos..pos + k as usize];
        println!("{} → {:x?}", kmer, hashes);
    }

    Ok(())
}
```

### Low‑Level API

If you prefer to manage the rolling yourself:

```rust
use nthash_rs::NtHash;

let seq = b"ACGTCAGTNNNNACGTACGT";
let k = 4u16;
let m = 2u8;

let mut h = NtHash::new(seq, k, m, 0)?;
if h.roll() {
    println!("first hash: {:#x}", h.hashes()[0]);
}
while h.roll() {
    println!("next hash: {:#x}", h.forward_hash());
}
```

## Fast paths

The rolling hash is serial by construction, so a single chain is bound by the
`srol` dependency latency. Where the workload allows independent chains, this
crate rolls several at once in SIMD lanes (AVX-512F / AVX2 on `x86_64`, NEON on
`aarch64`, detected at run time; all backends are bit-identical to the scalar
one).

```rust
use nthash_rs::{wide, NtHashError, ReadHasher};

fn main() -> Result<(), NtHashError> {
    let seq = b"ACGTACGTACGTACGTACGTACGT";

    // One long sequence, split across lanes. Unordered is the fastest form —
    // use it when the consumer folds (a sketch, a Bloom filter, a counter).
    let mut xor = 0u64;
    wide::for_each_canonical_unordered(seq, 4, |_pos, hash| xor ^= hash)?;

    // Or materialise in position order: `out[i]` is the k-mer starting at `i`.
    let out = wide::hash_canonical(seq, 4)?;
    assert_eq!(out.len(), seq.len() - 4 + 1);

    // Many equal-length reads: one read per lane.
    let reads: &[&[u8]] = &[b"ACGTACGT", b"TTGGCCAA", b"GATTACAG", b"CCGGTTAA"];
    ReadHasher::new(4)?.hash_reads(reads, |_read, _pos, _hash| {})?;

    Ok(())
}
```

Both assume **N-free A/C/G/T** input; use `NtHash` when the sequence may contain
ambiguous bases. With the `parallel` feature, `parallel::par_hash_canonical`
spreads the same kernels across cores.

## License

This project is MIT‑licensed (see [LICENSE](LICENSE)).
