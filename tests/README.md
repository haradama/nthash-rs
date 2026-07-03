# Correctness strategy

`nthash-rs` is validated on three independent layers. Together they pin the
output to the C++ reference *and* prove structural correctness on every input.

## 1. Golden vectors (absolute correctness vs. the C++ reference)

`regression_nthash.rs`, `regression_blindnthash.rs`, and `regression_seednthash.rs`
assert bit-exact equality against hash values produced by the upstream C++
[ntHash](https://github.com/bcgsc/ntHash) / btllib reference. These are the only
tests that can catch an error in the lookup **tables** themselves (a wrong seed
or split-rotate constant), since every Rust variant shares those tables.

### Regenerating / extending the vectors

The vectors are generated once from the reference and committed. To add more
(new `k`, longer sequences, more spaced seeds):

```bash
# 1. Build the reference
git clone https://github.com/bcgsc/btllib && cd btllib
meson setup build && ninja -C build

# 2. Emit (sequence, k, seed, hashes) tuples with a small driver that calls
#    btllib::NtHash / SeedNtHash, then paste them into the regression tests.
```

A future CI job can build btllib and diff a large random panel automatically;
until then the committed panel is the anchor. Keep additions small and readable
— the point is to cover *new shapes* (k around the 4-mer chunk boundary, k=1,
k>sequence edge cases, seeds with N at care sites), not volume.

## 2. Property tests (structural correctness on every input)

`properties.rs` uses `proptest` to assert the invariants that hold for *all*
inputs, so they cannot be satisfied by an accidental constant:

| Property | Why it matters |
|---|---|
| Purity: rolling == isolated recomputation | the hash is a pure function of the k-mer |
| Canonical invariant under reverse-complement | the defining property of ntHash |
| Strand symmetry: `fwd(kmer) == rev(revcomp)` | forward/reverse tables are consistent |
| `roll` ∘ `roll_back` = identity | the incremental update is exactly invertible |
| `NtHash` == `BlindNtHash`; all-ones `SeedNtHash` == `NtHash` | the three code paths agree |

These are an *oracle-free* net: no reference implementation is needed to run
them, so they gate every PR in ordinary CI.

## 3. Fuzzing (robustness — no panics, no hangs)

`../fuzz/` holds a `cargo-fuzz` target driving all three hashers over arbitrary
bytes, `k`, `pos`, and mask strings. Construction may return `Err`, but a panic,
an out-of-bounds slice, or a non-terminating roll is a bug.

```bash
cargo +nightly fuzz run fuzz_nthash
```

CI build-checks the target on every push (`fuzz-build` job) so it can't rot; run
campaigns out-of-band before releases. The `no_panic_*` cases in `properties.rs`
are a lightweight stand-in that runs in normal CI.

> Both bugs fixed in the change that introduced these tests — `SeedNtHash`
> silently stopping at the first `N`, and `BlindNtHash` panicking on a sequence
> shorter than `k` — are now covered by layers 2 and 3.
