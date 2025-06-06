use nthash_rs::NtHashError;
use nthash_rs::blind::{BlindNtHash, BlindNtHashBuilder};

fn main() -> Result<(), NtHashError> {
        println!("# BlindNtHash");
        let seq = "ATCGTACGNNNNNNNNATGCTGACG";
        let kmer_size: u16 = 6;
        let num_hashes: u8 = 3;

        println!("## BlindNtHash Low-Level API");
        let mut h = BlindNtHash::new(seq.as_bytes(), kmer_size, num_hashes, 0)?;
        for incoming in seq.as_bytes()[kmer_size as usize..].iter().copied() {
            h.roll(incoming);
    
            let pos   = h.pos() as usize;
            let end = pos + kmer_size as usize;
            let kmer  = &seq[pos..end];
            let hashes = h.hashes();
            println!("{} {:x?}", kmer, hashes);
        }

        println!("## BlindNtHashBuilder");
        let iter = BlindNtHashBuilder::new(seq.as_bytes())
            .k(kmer_size)
            .num_hashes(num_hashes)
            .pos(0)
            .finish()?;

        for (pos, hashes) in iter {
            let end = pos + kmer_size as usize;
            println!("{} {:x?}", &seq[pos..end], hashes);
        }

    Ok(())
}
