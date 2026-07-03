use nthash_rs::{MinimizerBuilder, NtHashError};

fn main() -> Result<(), NtHashError> {
    println!("# Canonical (w, k)-minimizers");
    let seq = "ACGTCAGTNNNNACGTACGTTGCA";
    let k: u16 = 5;
    let w: usize = 4; // window of 4 consecutive k-mers

    // Minimizers are strand-independent and skip windows spanning 'N'.
    for m in MinimizerBuilder::new(seq.as_bytes()).k(k).w(w).finish()? {
        let kmer = &seq[m.pos..m.pos + k as usize];
        println!("pos {:>2}  {}  {:#018x}", m.pos, kmer, m.hash);
    }

    // The same sequence and its reverse complement share the same minimizer
    // hash set — the property that makes canonical minimizers useful for
    // strand-agnostic mapping and sketching.
    Ok(())
}
