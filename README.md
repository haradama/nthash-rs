# ntHash‑rs

Pure‑Rust port of [ntHash](https://github.com/bcgsc/ntHash) rolling‑hash suite, focused on contiguous k‑mer hashing for DNA sequences.

> **Note**: At the moment **only the canonical contiguous k‑mer hasher** (`kmer::NtHash`) is implemented.  
> Blind‐mode (`BlindNtHash`) and spaced‑seed variants (`SeedNtHash`, `BlindSeedNtHash`) are planned but not yet available.

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
nthash-rs = "0.1"
````

## Quick Start

```rust
use nthash_rs::{NtHashBuilder, NtHashError};

fn main() -> Result<(), NtHashError> {
    let seq = "ACGTCAGTNNNNACGTACGT";
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
use nthash_rs::kmer::NtHash;

let mut h = NtHash::new("ACGTCAGTACGT", 5, 1, 0)?;
if h.roll() {
    println!("first hash: {:#x}", h.hashes()[0]);
}
while h.roll() {
    println!("next hash: {:#x}", h.forward_hash());
}
```

## License

This project is MIT‑licensed (see [LICENSE](LICENSE)).
