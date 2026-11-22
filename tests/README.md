# Correctness strategy

`nthash-rs` is validated on independent layers. Together they pin the
output to the C++ reference *and* prove structural correctness on every input.

## 1. Reference vectors (bit-compatibility with the C++ implementation)

`reference_vectors.rs` checks every hasher that has a C++ counterpart —
`NtHash`, `BlindNtHash` and `SeedNtHash` — against output of the upstream C++
[ntHash](https://github.com/bcgsc/ntHash) library: positions, forward and
reverse strand hashes and every output hash, for every k-mer the reference
yields and no others. These are the only tests that can catch an error in the
lookup **tables** themselves (a wrong seed or split-rotate constant), since
every Rust variant shares those tables.

The expected values live in `reference/vectors.txt`. They are not written by
hand: `reference/generate.sh` fetches the reference at a pinned commit,
compiles `reference/gen_vectors.cpp` against it and captures the output.

```bash
tests/reference/generate.sh            # rewrite vectors.txt
tests/reference/generate.sh --check    # fail if vectors.txt is out of date
```

CI runs the `--check` form (the `reference` job), so the committed vectors are
guaranteed to be exactly what the reference produces. To cover a new shape, add
a case to `gen_vectors.cpp` and regenerate; the Rust test picks it up.

`regression_nthash.rs`, `regression_blindnthash.rs` and
`regression_seednthash.rs` are older, hand-transcribed vectors from the same
reference, kept as a small readable example of each hasher.

### Where the reference cannot be followed

The panel leaves out inputs for which the reference has no well-defined output.
On these the Rust implementation is deliberately *not* compatible:

| Input | Reference | This crate |
|---|---|---|
| `k < 3` (contiguous hashers) | reads out of bounds (`k - 3` underflows) | supported |
| Spaced seeds over a sequence containing `N` | hashes `N` as a zero contribution, or skips `k` positions, depending on where the `N` sits | skips exactly the windows with `N` at a care position |
| `BlindNtHash` constructed at `pos > 0` | hashes the first `k` bases regardless of `pos` | hashes the window at `pos` |
| `roll_back` across an `N` | re-scans forward and returns to where it started | steps to the previous valid k-mer |
| A read shorter than `k`, `k == 0`, ... | throws | returns `Err` |

`regression_edge_cases.rs` pins down edge cases that once panicked, hung, or
returned the wrong k-mer (rolling back across or after an `N`, the first call on
a fresh hasher, oversized minimizer windows, duplicate care indices).

## 2. Property tests (structural correctness on every input)

`properties.rs` uses `proptest` to assert the invariants that hold for *all*
inputs, so they cannot be satisfied by an accidental constant:

| Property | Why it matters |
|---|---|
| Purity: rolling == isolated recomputation | the hash is a pure function of the k-mer |
| Canonical invariant under reverse-complement | the defining property of ntHash |
| Strand symmetry: `fwd(kmer) == rev(revcomp)` | forward/reverse tables are consistent |
| `roll` ∘ `roll_back` = identity | the incremental update is exactly invertible |
| `roll_back` from the end retraces the forward pass (with `N`) | backward `N`-skipping visits the same k-mers and terminates |
| `peek` / `peek_back` == the neighbouring k-mer's hashes | peeking is consistent with rolling and does not move the window |
| care indices behave as a set; shared `SpacedSeeds` == per-hasher masks | seed precompilation is input-order independent and reusable |
| `NtHash` == `BlindNtHash`; all-ones `SeedNtHash` == `NtHash` | the independent code paths agree |

These are an *oracle-free* net: no reference implementation is needed to run
them, so they gate every PR in ordinary CI.

## 3. Fuzzing (robustness — no panics, no hangs)

`../fuzz/` holds a `cargo-fuzz` target driving every hasher over arbitrary
bytes, `k`, `pos`, and mask strings. Construction may return `Err`, but a panic,
an out-of-bounds slice, or a non-terminating roll is a bug.

```bash
cargo +nightly fuzz run fuzz_nthash
```

CI build-checks the target on every push (`fuzz-build` job) so it can't rot; run
campaigns out-of-band before releases. The `no_panic_*` cases in `properties.rs`
are a lightweight stand-in that runs in normal CI.

## Coverage

CI requires the unit, property and regression tests together to execute every
line of the library (the `coverage` job). To check locally:

```bash
cargo install cargo-llvm-cov        # needs the llvm-tools-preview component
CARGO_PROFILE_DEV_OPT_LEVEL=0 CARGO_PROFILE_TEST_OPT_LEVEL=0 \
  cargo llvm-cov --fail-under-lines 100 --show-missing-lines
```

Things to know when the gate fails unexpectedly:

- Measure at opt-level 0, as above. With optimizations on, a function inlined
  into several test binaries is instrumented differently in each; the report
  then warns that "functions have mismatched data" and drops their counts, so
  the result varies from run to run.
- Do not rely on a property test to reach a line by chance (a random `k`
  happening to be `0`, say). Give every error path a deterministic test.
- Use one target directory per toolchain; objects left by another compiler
  version are picked up by the report and show as unexecuted.

- A generic function is measured per instantiation, and the summary takes the
  best single one. A method taking a closure therefore needs *one* call site
  that reaches all of its lines; covering the halves from different call sites
  does not count.
- Code that cannot be reached is a finding in itself: remove it rather than
  excluding it from the report.
