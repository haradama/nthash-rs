use nthash_rs::{Result, SeedNtHash, SeedNtHashBuilder};

fn main() -> Result<()> {
    println!("# SeedNtHash");
    let seq = "ATCGTACGATGCATGCATGCTGACG";
    let seed_masks = vec![
        "000111".to_string(),
        "010101".to_string(),
    ];
    let k   = 6u16;
    let m2  = 2usize;

    println!("## NtHash Low-Level API");
    let mut h = SeedNtHash::new(seq.as_bytes(), &seed_masks, m2, k, 0)?;
    while h.roll() {
        let pos   = h.pos() as usize;
        let end = pos + k as usize;
        let kmer  = &seq[pos..end];
        let hashes = h.hashes();
        println!("{} {:x?}", kmer, hashes);
    }

    println!("## SeedNtHashBuilder");
    let iter = SeedNtHashBuilder::new(seq.as_bytes())
        .k(k)
        .masks(seed_masks)
        .num_hashes(m2)
        .pos(0)
        .finish()?;

    for (pos, hashes) in iter {
        let end = pos + k as usize;
        let kmer  = &seq[pos..end];
        println!("{} {:x?}", kmer, hashes);
    }

    Ok(())
}
