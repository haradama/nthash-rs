//! Bit-compatibility with the C++ reference implementation.
//!
//! `tests/reference/vectors.txt` holds the output of the C++ ntHash library
//! (bcgsc/ntHash, at the commit recorded in the file's header) for a panel of
//! deterministic inputs; `tests/reference/generate.sh` rebuilds the reference
//! and regenerates it, and CI checks that the committed file is exactly what
//! the reference produces. Every hasher that has a C++ counterpart must
//! reproduce those values bit for bit — positions, strand hashes and every
//! output hash, for every k-mer the reference yields and no others.
//!
//! The file format is described at the top of `tests/reference/gen_vectors.cpp`.

use nthash_rs::{
    BlindNtHash, BlindNtHashBuilder, NtHash, NtHashBuilder, ReadHasher, SeedNtHash,
    SeedNtHashBuilder, SpacedSeeds,
};

const VECTORS: &str = include_str!("reference/vectors.txt");

/// One k-mer as reported by the reference: its position, the strand hashes
/// (absent for spaced seeds) and the output hashes.
#[derive(Debug, PartialEq, Eq)]
struct Row {
    pos: i64,
    strands: Option<(u64, u64)>,
    hashes: Vec<u64>,
}

struct Case<'a> {
    kind: &'a str,
    k: u16,
    num_hashes: usize,
    masks: Vec<&'a str>,
    seq: &'a [u8],
    rows: Vec<Row>,
    line: usize,
}

fn hex(s: &str) -> u64 {
    u64::from_str_radix(s, 16).expect("hash is hexadecimal")
}

fn parse() -> Vec<Case<'static>> {
    let mut cases: Vec<Case> = Vec::new();
    for (i, line) in VECTORS.lines().enumerate() {
        let mut fields = line.split(' ');
        match fields.next().unwrap() {
            tag if tag.starts_with('#') => {}
            "=" => {
                let case = cases.last_mut().expect("row before any case");
                let pos = fields.next().unwrap().parse().unwrap();
                let strands = (case.kind != "S")
                    .then(|| (hex(fields.next().unwrap()), hex(fields.next().unwrap())));
                cases.last_mut().unwrap().rows.push(Row {
                    pos,
                    strands,
                    hashes: fields.map(hex).collect(),
                });
            }
            kind => {
                let k = fields.next().unwrap().parse().unwrap();
                let num_hashes = fields.next().unwrap().parse().unwrap();
                let masks = if kind == "S" {
                    fields.next().unwrap().split(',').collect()
                } else {
                    Vec::new()
                };
                cases.push(Case {
                    kind,
                    k,
                    num_hashes,
                    masks,
                    seq: fields.next().unwrap().as_bytes(),
                    rows: Vec::new(),
                    line: i + 1,
                });
            }
        }
    }
    cases
}

fn cases(kind: &str) -> Vec<Case<'static>> {
    let selected: Vec<_> = parse().into_iter().filter(|c| c.kind == kind).collect();
    assert!(!selected.is_empty(), "no `{kind}` cases in vectors.txt");
    selected
}

fn nthash_row(h: &NtHash) -> Row {
    Row {
        pos: h.pos() as i64,
        strands: Some((h.forward_hash(), h.reverse_hash())),
        hashes: h.hashes().to_vec(),
    }
}

fn blind_row(h: &BlindNtHash) -> Row {
    Row {
        pos: h.pos() as i64,
        strands: Some((h.forward_hash(), h.reverse_hash())),
        hashes: h.hashes().to_vec(),
    }
}

/// `(pos, canonical hash)` of every reference row.
fn canonical_rows(case: &Case) -> Vec<(usize, u64)> {
    case.rows
        .iter()
        .map(|r| (r.pos as usize, r.hashes[0]))
        .collect()
}

#[test]
fn nthash_forward_matches_reference() {
    for case in cases("N") {
        let at = format!("vectors.txt:{} (k={})", case.line, case.k);
        let mut h = NtHash::new(case.seq, case.k, case.num_hashes as u8, 0).unwrap();
        let mut got = Vec::new();
        while h.roll() {
            got.push(nthash_row(&h));
        }
        assert_eq!(got, case.rows, "{at}");

        // The other ways of driving the same hasher agree with the reference.
        let builder = || {
            NtHashBuilder::new(case.seq)
                .k(case.k)
                .num_hashes(case.num_hashes as u8)
                .finish()
                .unwrap()
        };
        let canonical = canonical_rows(&case);
        assert_eq!(builder().canonical().collect::<Vec<_>>(), canonical, "{at}");
        let mut via_callback = Vec::new();
        builder().for_each_canonical(|pos, hash| via_callback.push((pos, hash)));
        assert_eq!(via_callback, canonical, "{at}");
        let via_iter: Vec<(usize, Vec<u64>)> = builder().collect();
        let expected: Vec<(usize, Vec<u64>)> = case
            .rows
            .iter()
            .map(|r| (r.pos as usize, r.hashes.clone()))
            .collect();
        assert_eq!(via_iter, expected, "{at}");

        // On input the reference hashes in full, the batch hasher agrees too.
        if case.rows.len() == case.seq.len() - case.k as usize + 1 {
            let mut batch = Vec::new();
            ReadHasher::new(case.k)
                .unwrap()
                .hash_reads(&[case.seq], |_, pos, hash| batch.push((pos, hash)))
                .unwrap();
            assert_eq!(batch, canonical, "{at}");
        }
    }
}

#[test]
fn nthash_backward_matches_reference() {
    for case in cases("B") {
        let mut h = NtHash::new(case.seq, case.k, case.num_hashes as u8, 0).unwrap();
        while h.roll() {}
        let mut got = Vec::new();
        while h.roll_back() {
            got.push(nthash_row(&h));
        }
        assert_eq!(got, case.rows, "vectors.txt:{} (k={})", case.line, case.k);
    }
}

#[test]
fn blind_forward_matches_reference() {
    for case in cases("L") {
        let at = format!("vectors.txt:{} (k={})", case.line, case.k);
        let mut h = BlindNtHash::new(case.seq, case.k, case.num_hashes as u8, 0).unwrap();
        let mut got = vec![blind_row(&h)];
        for &incoming in &case.seq[case.k as usize..] {
            h.roll(incoming);
            got.push(blind_row(&h));
        }
        assert_eq!(got, case.rows, "{at}");

        let canonical: Vec<(usize, u64)> = BlindNtHashBuilder::new(case.seq)
            .k(case.k)
            .finish()
            .unwrap()
            .canonical()
            .collect();
        assert_eq!(canonical, canonical_rows(&case), "{at}");
    }
}

#[test]
fn blind_backward_matches_reference() {
    for case in cases("R") {
        let k = case.k as usize;
        let mut h = BlindNtHash::new(case.seq, case.k, case.num_hashes as u8, 0).unwrap();
        for &incoming in &case.seq[k..] {
            h.roll(incoming);
        }
        let mut got = vec![blind_row(&h)];
        for &incoming in case.seq[..case.seq.len() - k].iter().rev() {
            h.roll_back(incoming);
            got.push(blind_row(&h));
        }
        assert_eq!(got, case.rows, "vectors.txt:{} (k={})", case.line, case.k);
    }
}

#[test]
fn seed_matches_reference() {
    for case in cases("S") {
        let at = format!("vectors.txt:{} (k={})", case.line, case.k);
        let mut h = SeedNtHash::new(case.seq, &case.masks, case.num_hashes, case.k, 0).unwrap();
        let mut got = Vec::new();
        while h.roll() {
            got.push(Row {
                pos: h.pos() as i64,
                strands: None,
                hashes: h.hashes().to_vec(),
            });
        }
        assert_eq!(got, case.rows, "{at}");

        // Precompiled seeds, through the builder, give the same values.
        let seeds = SpacedSeeds::from_masks(&case.masks, case.k).unwrap();
        let shared: Vec<(usize, Vec<u64>)> = SeedNtHashBuilder::new(case.seq)
            .seeds(&seeds)
            .num_hashes(case.num_hashes)
            .finish()
            .unwrap()
            .collect();
        let expected: Vec<(usize, Vec<u64>)> = case
            .rows
            .iter()
            .map(|r| (r.pos as usize, r.hashes.clone()))
            .collect();
        assert_eq!(shared, expected, "{at}");
    }
}
