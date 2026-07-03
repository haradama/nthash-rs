use nthash_rs::{NtHashError, ReadHasher};

fn main() -> Result<(), NtHashError> {
    println!("# Batch read hashing (AVX2 across reads where available)");

    // Many equal-length, N-free reads — the shape of short-read sequencing data.
    let reads: &[&[u8]] = &[
        b"ACGTACGTACGT",
        b"TTGGCCAATTGG",
        b"GATTACAGATTA",
        b"CCGGTTAACCGG",
        b"AAAACCCCGGGG",
    ];
    let k: u16 = 5;

    let hasher = ReadHasher::new(k)?;
    hasher.hash_reads(reads, |read_idx, pos, hash| {
        println!("read {read_idx}  pos {pos:>2}  {hash:#018x}");
    })?;

    Ok(())
}
