// Emits reference vectors from the C++ ntHash library (bcgsc/ntHash).
//
// Built and run by generate.sh against a pinned commit of the reference; the
// output is committed as vectors.txt and checked by tests/reference_vectors.rs.
// Everything here is deterministic: the inputs come from a fixed-seed
// generator and are written into the output next to the expected values, so
// the Rust side never has to reproduce the generator.
//
// Output format — one case per block, fields separated by single spaces:
//
//   N <k> <num_hashes> <seq>            NtHash, rolled forward to the end
//   B <k> <num_hashes> <seq>            NtHash, rolled to the end, then back
//   L <k> <num_hashes> <seq>            BlindNtHash, fed seq[k..] forward
//   R <k> <num_hashes> <seq>            BlindNtHash, fed forward to the end, then backward
//   S <k> <num_hashes> <masks> <seq>    SeedNtHash (masks comma-separated)
//   = <pos> [<fwd> <rev>] <hash>...     one row per k-mer the reference yields
//
// Rows of N, B, L and R carry the forward and reverse strand hashes before the
// output hashes; SeedNtHash does not expose them. All hashes are 16 hex digits.

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "nthash/nthash.hpp"

namespace {

// xorshift64: small, fixed, and independent of the C++ standard library.
struct Rng
{
  uint64_t s;
  uint64_t next()
  {
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    return s;
  }
  unsigned below(unsigned n) { return (next() >> 33) % n; }
};

std::string
random_seq(Rng& rng, size_t len, const char* alphabet)
{
  const size_t n = std::string(alphabet).size();
  std::string s(len, 'A');
  for (auto& c : s) {
    c = alphabet[rng.below(n)];
  }
  return s;
}

void
print_hashes(const uint64_t* h, unsigned n)
{
  for (unsigned i = 0; i < n; i++) {
    std::printf(" %016llx", (unsigned long long)h[i]);
  }
  std::printf("\n");
}

template<class H>
void
row(const H& h, long pos, unsigned num_hashes)
{
  std::printf("= %ld %016llx %016llx",
              pos,
              (unsigned long long)h.get_forward_hash(),
              (unsigned long long)h.get_reverse_hash());
  print_hashes(h.hashes(), num_hashes);
}

void
nthash_forward(const std::string& seq, unsigned k, unsigned m)
{
  std::printf("N %u %u %s\n", k, m, seq.c_str());
  nthash::NtHash h(seq, m, k);
  while (h.roll()) {
    row(h, (long)h.get_pos(), m);
  }
}

// Only for sequences without ambiguous bases.
void
nthash_backward(const std::string& seq, unsigned k, unsigned m)
{
  std::printf("B %u %u %s\n", k, m, seq.c_str());
  nthash::NtHash h(seq, m, k);
  while (h.roll()) {
  }
  while (h.roll_back()) {
    row(h, (long)h.get_pos(), m);
  }
}

void
blind_forward(const std::string& seq, unsigned k, unsigned m)
{
  std::printf("L %u %u %s\n", k, m, seq.c_str());
  nthash::BlindNtHash h(seq.data(), m, k);
  row(h, (long)h.get_pos(), m);
  for (size_t i = k; i < seq.size(); i++) {
    h.roll(seq[i]);
    row(h, (long)h.get_pos(), m);
  }
}

// Starts from the first window and rolls to the last before turning round:
// the reference's constructor hashes seq[0..k] whatever `pos` is given, so a
// BlindNtHash constructed at a later position starts out inconsistent.
void
blind_backward(const std::string& seq, unsigned k, unsigned m)
{
  std::printf("R %u %u %s\n", k, m, seq.c_str());
  const long last = (long)(seq.size() - k);
  nthash::BlindNtHash h(seq.data(), m, k);
  for (size_t i = k; i < seq.size(); i++) {
    h.roll(seq[i]);
  }
  row(h, (long)h.get_pos(), m);
  for (long i = last - 1; i >= 0; i--) {
    h.roll_back(seq[i]);
    row(h, (long)h.get_pos(), m);
  }
}

void
seed_forward(const std::string& seq,
             const std::vector<std::string>& masks,
             unsigned k,
             unsigned m)
{
  std::printf("S %u %u ", k, m);
  for (size_t i = 0; i < masks.size(); i++) {
    std::printf("%s%s", i ? "," : "", masks[i].c_str());
  }
  std::printf(" %s\n", seq.c_str());
  nthash::SeedNtHash h(seq, masks, m, k);
  while (h.roll()) {
    std::printf("= %ld", (long)h.get_pos());
    print_hashes(h.hashes(), m * (unsigned)masks.size());
  }
}

// A sequence with ambiguous bases: one near the start, and one further in
// when the sequence is long enough to leave valid k-mers on both sides.
std::string
with_n(std::string seq, unsigned k)
{
  seq[1] = 'N';
  if (seq.size() >= 2 * (size_t)k + 4) {
    seq[k + 2] = 'N';
  }
  return seq;
}

std::string
random_mask(Rng& rng, unsigned k, unsigned zero_one_in)
{
  std::string m(k, '1');
  for (unsigned i = 1; i + 1 < k; i++) {
    if (rng.below(zero_one_in) == 0) {
      m[i] = '0';
    }
  }
  return m;
}

} // namespace

int
main()
{
  Rng rng{ 0x9E3779B97F4A7C15ULL };
  const unsigned extra = 14; // windows per sequence beyond the first

  // Around the 4-base chunk boundary of the initial hash, around the 31/33/64
  // bit rotation periods, and well past them.
  //
  // k = 1 and k = 2 are absent on purpose: the reference's base_forward_hash
  // loops while `i < k - 3` on an unsigned k, which underflows for k < 3 and
  // reads out of bounds. There is no reference output to be compatible with.
  const unsigned ks[] = { 3,  4,  5,  6,  7,  8,  9,   12,  16,   17,  30, 31, 32,
                          33, 34, 48, 63, 64, 65, 66,  96,  128,  200, 1023, 1024 };
  unsigned i = 0;
  for (unsigned k : ks) {
    const unsigned m = 1 + i++ % 3;
    nthash_forward(random_seq(rng, k + extra, "ACGT"), k, m);
    nthash_forward(with_n(random_seq(rng, 2 * k + extra, "ACGT"), k), k, m);
    nthash_forward(random_seq(rng, k + extra, "ACGTacgt"), k, m);
    nthash_forward(random_seq(rng, k + extra, "ACGUacgu"), k, m);
    if (k <= 8) {
      // Dense ambiguity, including runs with no valid k-mer in between.
      nthash_forward(random_seq(rng, 6 * k + extra, "ACGTACGTACGTN"), k, m);
    }
  }
  nthash_forward(std::string(40, 'N'), 5, 1);

  const unsigned back_ks[] = { 3, 4, 5, 16, 31, 33, 64, 65, 200 };
  for (unsigned k : back_ks) {
    const unsigned m = 1 + i++ % 3;
    nthash_backward(random_seq(rng, k + extra, "ACGTacgt"), k, m);
    blind_backward(random_seq(rng, k + extra, "ACGT"), k, m);
  }

  for (unsigned k : ks) {
    const unsigned m = 1 + i++ % 3;
    blind_forward(random_seq(rng, k + extra, "ACGT"), k, m);
    blind_forward(random_seq(rng, k + extra, "ACGUacgu"), k, m);
  }

  // Spaced seeds are compared on unambiguous input only. With an `N` the
  // reference is not a usable oracle: its window initialisation tests each
  // base against the numeric SEED_N (0) rather than looking it up, so an `N`
  // is hashed as a zero contribution there, while a roll skips k positions
  // only when the `N` arrives as the last base of the window.
  const unsigned seed_ks[] = { 2, 3, 4, 5, 8, 12, 16, 21, 31, 32, 33, 64, 65, 100 };
  for (unsigned k : seed_ks) {
    const unsigned m = 1 + i++ % 3;
    const std::string solid(k, '1');
    const std::string sparse = random_mask(rng, k, 2);
    const std::string blocky = random_mask(rng, k, 6);
    seed_forward(random_seq(rng, k + extra, "ACGT"), { solid }, k, m);
    seed_forward(random_seq(rng, k + extra, "ACGTacgu"), { sparse }, k, m);
    seed_forward(random_seq(rng, k + extra, "ACGT"), { blocky, sparse, solid }, k, m);
    seed_forward(random_seq(rng, 2 * k + extra, "ACGUacgt"), { sparse, blocky }, k, m);
  }
  return 0;
}
