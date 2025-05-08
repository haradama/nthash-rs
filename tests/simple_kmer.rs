use nthash::kmer::NtHashBuilder;

const SEQ: &str = "ATCGTACGATGCATGCATGCTGACG";
const K: u16 = 6;
const M: u8 = 3;

#[test]
fn regression_simple_kmer() {
    // build our iterator
    let iter = NtHashBuilder::new(SEQ)
        .k(K)
        .num_hashes(M)
        .pos(0)
        .finish()
        .expect("builder should succeed");

    // expected windows
    let expected_kmers = [
        "ATCGTA", "TCGTAC", "CGTACG", "GTACGA", "TACGAT", "ACGATG", "CGATGC", "GATGCA", "ATGCAT",
        "TGCATG", "GCATGC", "CATGCA", "ATGCAT", "TGCATG", "GCATGC", "CATGCT", "ATGCTG", "TGCTGA",
        "GCTGAC", "CTGACG",
    ];

    // expected hashes for each window (hex literals)
    let expected_hashes: &[[u64; 3]] = &[
        [0x245f429174d6e9b1, 0x43def5f731c6a724, 0x683e389db9281069],
        [0xadc1da5b636f030c, 0xe58f7285f9e27675, 0x93514ce6c152b6fd],
        [0x959782747cd43a12, 0x7193a7d0cb009c55, 0x72b2a53932debf2],
        [0xadc1da5b636f030c, 0xe58f7285f9e27675, 0x93514ce6c152b6fd],
        [0x245f429174d6e9b1, 0x43def5f731c6a724, 0x683e389db9281069],
        [0x2879d1c755ac3c6d, 0x3cbb942540251d43, 0x653565e68baf68c3],
        [0x1f0b21c53d089b2c, 0x9fb0763a80534ae, 0x29062922f4ab50d2],
        [0x8f39107ffe942fd4, 0x9510d6508cf4ae41, 0x2449e6c6a4be4d9d],
        [0x276f4e13fd6666b8, 0x5af4471d5d9a7618, 0x8263953a4c0a5ed0],
        [0x94ace44e1d5e2ac2, 0x86f85a0ad91ffc3a, 0x1ba53e6b57d52738],
        [0x5f0d08530ea31f76, 0x736d77a48b344f82, 0xd27a7fe7ba7225ea],
        [0x94ace44e1d5e2ac2, 0x86f85a0ad91ffc3a, 0x1ba53e6b57d52738],
        [0x276f4e13fd6666b8, 0x5af4471d5d9a7618, 0x8263953a4c0a5ed0],
        [0x94ace44e1d5e2ac2, 0x86f85a0ad91ffc3a, 0x1ba53e6b57d52738],
        [0x5f0d08530ea31f76, 0x736d77a48b344f82, 0xd27a7fe7ba7225ea],
        [0x194dca4e9ac1ef13, 0x5410d41335d8f751, 0x6d5e9e65f957ae70],
        [0x1ec97a4a7b81aadf, 0x7a13c12e97c673b9, 0x98dd3b7f4ae76fe8],
        [0x855cf714947b7558, 0x88728b28e52873c5, 0xdcf824fc604c39f],
        [0x39f25bd0505b6512, 0x3b597bc698f7109, 0x3da7f38bdb669a11],
        [0xfc2267e8f5d65148, 0x8e6aaa7c9b150e82, 0x8a8d12471db4deb9],
    ];

    let k_usize = K as usize;
    let results: Vec<(usize, Vec<u64>)> = iter.collect();
    assert_eq!(results.len(), expected_kmers.len());

    for (i, (pos, hashes)) in results.iter().enumerate() {
        // check the sequence window
        let window = &SEQ[*pos..*pos + k_usize];
        assert_eq!(window, expected_kmers[i], "window at pos {}", pos);

        // check the three hash values
        assert_eq!(
            hashes.as_slice(),
            &expected_hashes[i][..],
            "hashes at pos {} (window {})",
            pos,
            window
        );
    }
}
