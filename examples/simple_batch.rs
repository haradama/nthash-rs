use nthash_rs::{NtHashError, ReadHasher};

fn main() -> Result<(), NtHashError> {
    println!("# Batch read hashing");

    // Many N-free reads — the shape of short-read sequencing data. Lengths may
    // differ, but a read shorter than k is an error.
    let reads: &[&[u8]] = &[
        b"ACGTACGTACGT",
        b"TTGGCCAATTGG",
        b"GATTACAGATTACAGA",
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
