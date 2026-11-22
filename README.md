# nthash‑rs

[<img alt="github" src="https://img.shields.io/badge/github-haradama/nthash__rs-8da0cb?style=for-the-badge&labelColor=555555&logo=github" height="20">](https://github.com/haradama/nthash-rs)
[<img alt="crates.io" src="https://img.shields.io/crates/v/nthash-rs.svg?style=for-the-badge&color=fc8d62&logo=rust" height="20">](https://crates.io/crates/nthash-rs)
[<img alt="docs.rs" src="https://img.shields.io/badge/docs.rs-nthash__rs-66c2a5?style=for-the-badge&labelColor=555555&logo=docs.rs" height="20">](https://docs.rs/nthash-rs)
[<img alt="build status" src="https://img.shields.io/github/actions/workflow/status/haradama/nthash-rs/rust.yml?branch=master&style=for-the-badge" height="20">](https://github.com/haradama/nthash-rs/actions)

Pure‑Rust port of the [ntHash](https://github.com/bcgsc/ntHash) rolling‑hash suite for DNA sequences, bit‑compatible with the C++ reference.

Compatibility is checked, not assumed: the test suite compares `NtHash`, `BlindNtHash` and `SeedNtHash` against vectors generated from the C++ library, and CI regenerates those vectors from a pinned commit of the reference on every run. The few inputs on which the reference has no well-defined output (for example `k < 3`, or spaced seeds over a sequence containing `N`) are listed in [tests/README.md](tests/README.md).

| Type | Use it for |
|---|---|
| `NtHash` / `NtHashBuilder` | contiguous k‑mers; windows containing `N` (or any other non‑ACGTU byte) are skipped |
| `BlindNtHash` / `BlindNtHashBuilder` | streaming over pre‑cleaned input, the caller feeds each base |
| `SeedNtHash` / `SeedNtHashBuilder` / `SpacedSeeds` | spaced seeds (one or more `0`/`1` masks) |
| `MinimizerBuilder` | canonical (w, k)-minimizers |
| `ReadHasher` | many N‑free reads with one reusable hasher |

Bases are case-insensitive, and RNA `U` hashes exactly like `T`.

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

    // Iterate over all valid k-mers (windows containing `N` are skipped).
    let iter = NtHashBuilder::new(seq).k(k).num_hashes(m).finish()?;

    for (pos, hashes) in iter {
        let kmer = String::from_utf8_lossy(&seq[pos..pos + k as usize]);
        println!("{kmer} → {hashes:x?}");
    }

    Ok(())
}
```

That iterator yields an owned `Vec<u64>` per k‑mer, which costs one allocation
each. Prefer the allocation-free forms when the hashes are consumed on the spot:

```rust
use nthash_rs::NtHashBuilder;

# fn main() -> nthash_rs::Result<()> {
let seq = b"ACGTCAGTNNNNACGTACGT";

// One canonical hash per k-mer, as a regular iterator.
for (pos, hash) in NtHashBuilder::new(seq).k(4).finish()?.canonical() {
    println!("{pos}: {hash:#x}");
}

// The same as a callback — the tightest loop.
NtHashBuilder::new(seq).k(4).finish()?
    .for_each_canonical(|pos, hash| println!("{pos}: {hash:#x}"));

// All `num_hashes` values per k-mer, borrowed from the hasher's buffer.
NtHashBuilder::new(seq).k(4).num_hashes(2).finish()?
    .for_each_hash(|pos, hashes| println!("{pos}: {hashes:x?}"));
# Ok(()) }
```

### Low‑Level API

If you prefer to manage the rolling yourself:

```rust
use nthash_rs::NtHash;

# fn main() -> nthash_rs::Result<()> {
let seq = b"ACGTCAGTNNNNACGTACGT";

let mut h = NtHash::new(seq, 4, 2, 0)?; // k = 4, 2 hashes per k-mer, start at 0
while h.roll() {
    println!("{}: canonical {:#x}, forward {:#x}", h.pos(), h.hashes()[0], h.forward_hash());
}
# Ok(()) }
```

### Spaced seeds

```rust
use nthash_rs::{SeedNtHash, SeedNtHashBuilder, SpacedSeeds};

# fn main() -> nthash_rs::Result<()> {
// One-off: masks are compiled for this sequence.
for (pos, hashes) in SeedNtHashBuilder::new(b"ATCGTACGATGCATGCATGCTGACG")
    .k(6)
    .masks(["110011", "101101"])
    .finish()?
{
    println!("{pos}: {hashes:x?}"); // one hash per mask
}

// Many sequences, same seeds: compile the lookup tables once and share them.
let seeds = SpacedSeeds::from_masks(&["110011", "101101"], 6)?;
for read in [&b"ACGTACGTAC"[..], &b"TTGACCATGA"[..]] {
    let mut h = SeedNtHash::with_seeds(read, &seeds, 1, 0)?;
    while h.roll() {
        println!("{}: {:x?}", h.pos(), h.hashes());
    }
}
# Ok(()) }
```

### Minimizers

```rust
use nthash_rs::MinimizerBuilder;

# fn main() -> nthash_rs::Result<()> {
for m in MinimizerBuilder::new(b"ACGTACGTTGCATGCA").k(5).w(4).finish()? {
    println!("minimizer at {}: {:#x}", m.pos, m.hash);
}
# Ok(()) }
```

### Batch hashing of reads

```rust
use nthash_rs::ReadHasher;

# fn main() -> nthash_rs::Result<()> {
// N-free reads of any length >= k; the per-k setup is paid once for the batch.
let reads = ["ACGTACGT", "TTGGCCAATT", "GATTAC", "CCGG"];
ReadHasher::new(4)?.hash_reads(&reads, |read_idx, pos, hash| {
    println!("read {read_idx} @ {pos}: {hash:#x}");
})?;
# Ok(()) }
```

## License

This project is MIT‑licensed (see [LICENSE](LICENSE)).
