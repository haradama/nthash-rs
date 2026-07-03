use nthash_rs::kmer::{NtHash, NtHashBuilder};
use nthash_rs::NtHashError;

fn main() -> Result<(), NtHashError> {
    println!("# NtHash");
    let seq = "ATCGTACGATGCATGCATGCTGACG";
    let kmer_size: u16 = 6;
    let num_hashes: u8 = 3;

    println!("## NtHash Low-Level API");
    let mut h = NtHash::new(seq.as_bytes(), kmer_size, num_hashes, 0)?;
    while h.roll() {
        let pos = h.pos();
        let end = pos + kmer_size as usize;
        let kmer = &seq[pos..end];
        let hashes = h.hashes();
        println!("{} {:x?}", kmer, hashes);
    }

    println!("## NtHashBuilder");
    let iter = NtHashBuilder::new(seq.as_bytes())
        .k(kmer_size)
        .num_hashes(num_hashes)
        .pos(0)
        .finish()?;

    for (pos, hashes) in iter {
        let end = pos + kmer_size as usize;
        let kmer = &seq[pos..end];
        println!("{} {:x?}", kmer, hashes);
    }

    Ok(())
}
