use std::hash::BuildHasher;
use std::hash::Hasher;

use ahash::RandomState;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use fnv::FnvHasher;
use fxhash::FxHasher64;
use xxhash_rust::xxh3::xxh3_64;

use nthash_rs::SeedNtHashBuilder;
use nthash_rs::{kmer::NtHashBuilder, BlindNtHashBuilder};

/// Generate a DNA sequence of length `len` by cycling through "ACGT".
///
/// The sequence is deliberately **N‑free**: the three ntHash variants handle
/// ambiguous bases differently (`NtHash` skips N‑windows, `BlindNtHash` assumes
/// none exist, `SeedNtHash` stops at the first N on a care site), so only an
/// N‑free input yields an apples‑to‑apples throughput comparison in which every
/// hasher processes all `len − k + 1` windows. N‑handling correctness is
/// covered by the regression tests, not here.
fn generate_dna(len: usize) -> String {
    const PATTERN: &[u8] = b"ACGT";
    (0..len)
        .map(|i| PATTERN[i % PATTERN.len()] as char)
        .collect()
}

fn bench_nthash(c: &mut Criterion) {
    let seq = generate_dna(1_000_000);
    let k: u16 = 31;
    let m: u8 = 1;

    let mut group = c.benchmark_group("nthash_vs_others");
    group.throughput(Throughput::Bytes(seq.len() as u64));

    group.bench_with_input(BenchmarkId::new("NtHash", seq.len()), &seq, |b, seq| {
        b.iter(|| {
            // build a new rolling iterator each iteration and consume it via the
            // zero-allocation path (no Vec<u64> per k-mer)
            NtHashBuilder::new(seq.as_bytes())
                .k(k)
                .num_hashes(m)
                .pos(0)
                .finish()
                .unwrap()
                .for_each_hash(|_pos, hashes| {
                    black_box(hashes);
                });
        })
    });

    group.finish();
}

fn bench_blindnthash(c: &mut Criterion) {
    let seq = generate_dna(1_000_000);
    let k: u16 = 31;
    let m: u8 = 1;

    let mut group = c.benchmark_group("nthash_vs_others");
    group.throughput(Throughput::Bytes(seq.len() as u64));

    group.bench_with_input(
        BenchmarkId::new("BlindNtHash", seq.len()),
        &seq,
        |b, seq| {
            b.iter(|| {
                BlindNtHashBuilder::new(seq.as_bytes())
                    .k(k)
                    .num_hashes(m)
                    .pos(0)
                    .finish()
                    .unwrap()
                    .for_each_hash(|_pos, hashes| {
                        black_box(hashes);
                    });
            })
        },
    );

    group.finish();
}

fn bench_seednthash(c: &mut Criterion) {
    let seq = generate_dna(1_000_000);
    let k: u16 = 31;
    let m: usize = 1;

    let mut group = c.benchmark_group("nthash_vs_others");
    group.throughput(Throughput::Bytes(seq.len() as u64));

    // All-ones mask: every position is a care site (contiguous k-mer), so the
    // spaced-seed path does real work comparable to the other variants. The
    // previous all-zeros mask had no care sites and hashed nothing.
    let mask = "1".repeat(k as usize);
    group.bench_with_input(BenchmarkId::new("SeedNtHash", seq.len()), &seq, |b, seq| {
        b.iter(|| {
            SeedNtHashBuilder::new(seq.as_bytes())
                .k(k)
                .masks(vec![mask.clone()])
                .num_hashes(m)
                .pos(0)
                .finish()
                .unwrap()
                .for_each_hash(|_pos, hashes| {
                    black_box(hashes);
                });
        })
    });

    group.finish();
}

fn bench_xxh3(c: &mut Criterion) {
    let seq = generate_dna(1_000_000);
    let k: usize = 31;

    let mut group = c.benchmark_group("nthash_vs_others");
    group.throughput(Throughput::Bytes(seq.len() as u64));

    group.bench_with_input(BenchmarkId::new("xxh3_64", seq.len()), &seq, |b, seq| {
        b.iter(|| {
            let bytes = seq.as_bytes();
            // slide a k‑mer window and hash each one with xxh3_64
            for i in 0..=bytes.len().saturating_sub(k) {
                black_box(xxh3_64(&bytes[i..i + k]));
            }
        })
    });

    group.finish();
}

fn bench_ahash(c: &mut Criterion) {
    let seq = generate_dna(1_000_000);
    let k: usize = 31;

    let mut group = c.benchmark_group("nthash_vs_others");
    group.throughput(Throughput::Bytes(seq.len() as u64));

    group.bench_with_input(BenchmarkId::new("ahash", seq.len()), &seq, |b, seq| {
        let state = RandomState::new();
        b.iter(|| {
            let bytes = seq.as_bytes();
            for i in 0..=bytes.len().saturating_sub(k) {
                let mut hasher = state.build_hasher();
                hasher.write(&bytes[i..i + k]);
                black_box(hasher.finish());
            }
        })
    });

    group.finish();
}

fn bench_fnv(c: &mut Criterion) {
    let seq = generate_dna(1_000_000);
    let k: usize = 31;

    let mut group = c.benchmark_group("nthash_vs_others");
    group.throughput(Throughput::Bytes(seq.len() as u64));

    group.bench_with_input(BenchmarkId::new("fnv64", seq.len()), &seq, |b, seq| {
        b.iter(|| {
            let bytes = seq.as_bytes();
            for i in 0..=bytes.len().saturating_sub(k) {
                let mut hasher = FnvHasher::default();
                hasher.write(&bytes[i..i + k]);
                black_box(hasher.finish());
            }
        })
    });

    group.finish();
}

fn bench_fxhash(c: &mut Criterion) {
    let seq = generate_dna(1_000_000);
    let k: usize = 31;

    let mut group = c.benchmark_group("nthash_vs_others");
    group.throughput(Throughput::Bytes(seq.len() as u64));

    group.bench_with_input(BenchmarkId::new("fxhash64", seq.len()), &seq, |b, seq| {
        b.iter(|| {
            let bytes = seq.as_bytes();
            for i in 0..=bytes.len().saturating_sub(k) {
                let mut hasher = FxHasher64::default();
                hasher.write(&bytes[i..i + k]);
                black_box(hasher.finish());
            }
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_nthash,
    bench_blindnthash,
    bench_seednthash,
    bench_xxh3,
    bench_fnv,
    bench_fxhash,
    bench_ahash
);
criterion_main!(benches);
