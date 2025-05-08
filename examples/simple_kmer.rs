use nthash::{kmer::NtHashBuilder, NtHashError};

fn main() -> Result<(), NtHashError> {
    let seq = "ATCGTACGATGCATGCATGCTGACG";
    let kmer_size: u16 = 6;
    let num_hashes: u8 = 3;

    // build the iterator
    let iter = NtHashBuilder::new(seq)
        .k(kmer_size)
        .num_hashes(num_hashes)
        .pos(0)
        .finish()?;

    for (pos, hashes) in iter {
        // make a usize end index for slicing
        let end = pos + kmer_size as usize;
        // print the k‑mer & the hashes (in hex debug)
        println!("{} {:x?}", &seq[pos..end], hashes);
    }

    Ok(())
}
