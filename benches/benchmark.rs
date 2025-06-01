use std::hash::BuildHasher;
use std::hash::Hasher;

use ahash::RandomState;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use nthash_rs::{kmer::NtHashBuilder, BlindNtHashBuilder};
use xxhash_rust::xxh3::xxh3_64;

/// Generate a pseudo‐random DNA sequence of length `len` by
/// repeating "ACGT" and inserting occasional "N"s.
fn generate_dna(len: usize) -> String {
    const PATTERN: &str = "ACGTN";
    let mut s = String::with_capacity(len);
    let mut i = 0;
    while s.len() < len {
        s.push(PATTERN.as_bytes()[i % PATTERN.len()] as char);
        i += 1;
    }
    s.truncate(len);
    s
}

fn bench_nthash(c: &mut Criterion) {
    let seq = generate_dna(1_000_000);
    let k: u16 = 31;
    let m: u8 = 1;

    let mut group = c.benchmark_group("nthash_vs_others");
    group.throughput(Throughput::Bytes(seq.len() as u64));

    group.bench_with_input(
        BenchmarkId::new("NtHash", seq.len()),
        &seq,
        |b, seq| {
            b.iter(|| {
                // build a new rolling iterator each iteration
                let mut iter = NtHashBuilder::new(seq.as_bytes())
                    .k(k)
                    .num_hashes(m)
                    .pos(0)
                    .finish()
                    .unwrap();
                // consume it
                while let Some((_pos, _hashes)) = iter.next() {
                    // no-op
                }
            })
        },
    );

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
                let mut iter = BlindNtHashBuilder::new(seq.as_bytes())
                    .k(k)
                    .num_hashes(m)
                    .pos(0)
                    .finish()
                    .unwrap();
                // consume it
                while let Some((_pos, _hashes)) = iter.next() {
                    // no-op
                }
            })
        },
    );

    group.finish();
}

fn bench_xxh3(c: &mut Criterion) {
    let seq = generate_dna(1_000_000);
    let k: usize = 31;

    let mut group = c.benchmark_group("nthash_vs_others");
    group.throughput(Throughput::Bytes(seq.len() as u64));

    group.bench_with_input(
        BenchmarkId::new("xxh3_64", seq.len()),
        &seq,
        |b, seq| {
            b.iter(|| {
                let bytes = seq.as_bytes();
                // slide a k‑mer window and hash each one with xxh3_64
                for i in 0..=bytes.len().saturating_sub(k) {
                    let _h = xxh3_64(&bytes[i..i + k]);
                }
            })
        },
    );

    group.finish();
}

fn bench_ahash(c: &mut Criterion) {
    let seq = generate_dna(1_000_000);
    let k: usize = 31;

    let mut group = c.benchmark_group("nthash_vs_others");
    group.throughput(Throughput::Bytes(seq.len() as u64));

    group.bench_with_input(
        BenchmarkId::new("ahash", seq.len()),
        &seq,
        |b, seq| {
            let state = RandomState::new();
            b.iter(|| {
                let bytes = seq.as_bytes();
                for i in 0..=bytes.len().saturating_sub(k) {
                    let mut hasher = state.build_hasher();
                    hasher.write(&bytes[i..i + k]);
                    let _h = hasher.finish();
                }
            })
        },
    );

    group.finish();
}

criterion_group!(benches, bench_nthash, bench_blindnthash, bench_xxh3, bench_ahash);
criterion_main!(benches);
