use nthash_rs::SeedNtHashBuilder;


#[test]
fn regression_simple_seednthash() {
    let seq = "ATCGTACGATGCATGCATGCTGACG";
    let seed_masks = vec![
        "000111".to_string(),
        "010101".to_string(),
    ];
    let k   = 6u16;
    let m2  = 2usize;

    let iter = SeedNtHashBuilder::new(seq.as_bytes())
        .k(k)
        .masks(seed_masks)
        .num_hashes(m2)
        .finish()
        .expect("builder should succeed");

    // expected windows
    let expected_kmers = [
        "ATCGTA", "TCGTAC", "CGTACG", "GTACGA", "TACGAT", "ACGATG", "CGATGC", "GATGCA", "ATGCAT",
        "TGCATG", "GCATGC", "CATGCA", "ATGCAT", "TGCATG", "GCATGC", "CATGCT", "ATGCTG", "TGCTGA",
        "GCTGAC", "CTGACG",
    ];

    // expected hashes for each window (hex literals)
    let expected_hashes: &[[u64; 4]] = &[
        [0x5d721caa40879845, 0x4eeedc1f3039a84c, 0x083865846584a5e7, 0x7e89a5c357dcdcfb],
        [0x651daa0fc1953543, 0x9dcb12ccfc9c403f, 0x3077784e47a59043, 0x28338c7283342427],
        [0x2d2be53a3e74ddd5, 0xa1ce7e5cc9bfaeff, 0xed343943be941e6d, 0x5d87d853f940b810],
        [0xeca5505260dc164b, 0x6da63c1524e034ad, 0xc1bfa252e6e0874b, 0x3c6078fe44975b86],
        [0x59402af97d1851c4, 0x7e427fd7ee96f930, 0xb06df31758e3ebdc, 0xcb3cf11867f830cc],
        [0x312507f99a02c6c9, 0x6adb14a798bd6bdd, 0x03d7caee0f0693a7, 0xa1b57910bbf6c4ba],
        [0xcfc6bbd2185e7043, 0x108fdbfe3f847552, 0x19248fe289ef1a09, 0x4ca1d5cbe41d248b],
        [0x847963a7a616c171, 0x8fe6d1d45ed0a139, 0x699a5d15bf8827c5, 0xb54ea82ee37d8b18],
        [0x9bc7d809ff4cb45e, 0xa38cb88768eb2d45, 0x6448484a013c4b60, 0xd4c2e85c8a6f3922],
        [0x312507f99a02c6c9, 0x6adb14a798bd6bdd, 0xed343943be941e6d, 0x5d87d853f940b810],
        [0xcfc6bbd2185e7043, 0x108fdbfe3f847552, 0x3077784e47a59043, 0x28338c7283342427],
        [0x847963a7a616c171, 0x8fe6d1d45ed0a139, 0x699a5d15bf8827c5, 0xb54ea82ee37d8b18],
        [0x9bc7d809ff4cb45e, 0xa38cb88768eb2d45, 0x6448484a013c4b60, 0xd4c2e85c8a6f3922],
        [0x312507f99a02c6c9, 0x6adb14a798bd6bdd, 0xed343943be941e6d, 0x5d87d853f940b810],
        [0xcfc6bbd2185e7043, 0x108fdbfe3f847552, 0x3077784e47a59043, 0x28338c7283342427],
        [0xe97cc92710b28516, 0xb31b6d7b9076f840, 0xb06df31758e3ebdc, 0xcb3cf11867f830cc],
        [0xb02e2852b9ef3eb3, 0x5b62f3dd45fa6b42, 0xbe27d8f1242443d9, 0x0ab639eef1c398a1],
        [0xb1859d23b75b711f, 0xb8c82cdd6236a58c, 0x48ce9177b6762755, 0x5577ac399ea0dd07],
        [0xe7da171022322abb, 0x9e3345eecb7493ce, 0x7324ce0d9d498bbb, 0x8a20451c8f3d48d6],
        [0x2d2be53a3e74ddd5, 0xa1ce7e5cc9bfaeff, 0x490ed7a78c06bb67, 0xe990dd1f2bdad4a8],
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
