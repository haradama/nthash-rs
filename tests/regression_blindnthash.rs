use nthash_rs::BlindNtHashBuilder;

#[test]
fn regression_simple_nthash() {
    let seq: &str = "ATCGTACGNNNNNNNNATGCTGACG";
    let k: u16 = 6;
    let m: u8 = 3;

    // build our iterator
    let iter = BlindNtHashBuilder::new(seq.as_bytes())
        .k(k)
        .num_hashes(m)
        .pos(0)
        .finish()
        .expect("builder should succeed");

    // expected windows
    let expected_kmers = [
        "ATCGTA", "TCGTAC", "CGTACG", "GTACGN", "TACGNN", "ACGNNN", "CGNNNN", "GNNNNN", "NNNNNN",
        "NNNNNN", "NNNNNN", "NNNNNA", "NNNNAT", "NNNATG", "NNATGC", "NATGCT", "ATGCTG", "TGCTGA",
        "GCTGAC", "CTGACG",
    ];

    // expected hashes for each window (hex literals)
    let expected_hashes: &[[u64; 3]] = &[
        [0x245f429174d6e9b1, 0x43def5f731c6a724, 0x683e389db9281069],
        [0xadc1da5b636f030c, 0xe58f7285f9e27675, 0x93514ce6c152b6fd],
        [0x959782747cd43a12, 0x7193a7d0cb009c55, 0x072b2a53932debf2],
        [0xe6e4e36432fd8854, 0x3c4252d79fc71686, 0x232736302e6c1251],
        [0x50b812075ad47599, 0x1c09f4d6436c42ad, 0x6cc206d08567960d],
        [0x3950cad04647d29c, 0xbcac852b4fce2337, 0xf5fd50139f0c56ec],
        [0x74719b60b88ed18f, 0xec144f2823ffabc7, 0x6085ea9a4ab84dc9],
        [0x37db9b8dad848fd4, 0x1147067318718822, 0x4922a1f7fa41ea03],
        [0x0000000000000000, 0x0000000000000000, 0x0000000000000000],
        [0x0000000000000000, 0x0000000000000000, 0x0000000000000000],
        [0x0000000000000000, 0x0000000000000000, 0x0000000000000000],
        [0x67353a3f120e8f48, 0x25bca433de634fc5, 0x8cf1de67e23d96bf],
        [0x546ea79a675318b2, 0x314f1d89b5c0c659, 0x85bdc53ab4ce351c],
        [0x312507f99a02c6c9, 0x6adb14a798bd6bdd, 0x9c001cb7dde151c5],
        [0x8f6e7e75c91742b7, 0x52747fc29e2d4de1, 0xe1e2fe22a5e63061],
        [0x06f7af81ccfc7b6d, 0x322f4c4e6d7b751c, 0x3926fbced1500eb8],
        [0x1ec97a4a7b81aadf, 0x7a13c12e97c673b9, 0x98dd3b7f4ae76fe8],
        [0x855cf714947b7558, 0x88728b28e52873c5, 0x0dcf824fc604c39f],
        [0x39f25bd0505b6512, 0x03b597bc698f7109, 0x3da7f38bdb669a11],
        [0xfc2267e8f5d65148, 0x8e6aaa7c9b150e82, 0x8a8d12471db4deb9],
    ];

    let k_usize = k as usize;
    let results: Vec<(usize, Vec<u64>)> = iter.collect();
    assert_eq!(results.len(), expected_kmers.len());

    for (i, (pos, hashes)) in results.iter().enumerate() {
        // check the sequence window
        let window = &seq[*pos..*pos + k_usize];
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


