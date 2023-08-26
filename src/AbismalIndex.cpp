/* Copyright (C) 2018-2023 Andrew D. Smith and Guilherme Sena
 *
 * Authors: Andrew D. Smith and Guilherme Sena
 *
 * This file is part of abismal.
 *
 * abismal is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * abismal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#include "AbismalIndex.hpp"

#include <algorithm>
#include <chrono>
#include <deque>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <utility>
#include <array>

#include <omp.h>

#include "dna_four_bit.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

using std::cerr;
using std::endl;
using std::ifstream;
using std::inclusive_scan;
using std::min;
using std::ofstream;
using std::runtime_error;
using std::sort;
using std::string;
using std::to_string;
using std::vector;

using std::chrono::duration;
using std::chrono::steady_clock;
using std::chrono::time_point;

template<typename T> using num_lim = std::numeric_limits<T>;

bool AbismalIndex::VERBOSE = false;
string AbismalIndex::internal_identifier = "AbismalIndex";

using genome_iterator = genome_four_bit_itr;

static string
delta_seconds(const time_point<steady_clock> &a,
              const time_point<steady_clock> &b) {
  std::ostringstream oss;
  oss.setf(std::ios::fixed);
  oss.precision(2);
  oss << "[time: " << duration<double>(b - a).count() << "s]";
  return oss.str();
}

void
AbismalIndex::create_index(const string &genome_file) {
  auto s_time = steady_clock::now();
  if (VERBOSE) cerr << "[loading genome]";
  vector<uint8_t> orig_genome;
  load_genome(genome_file, orig_genome, cl);
  if (VERBOSE) cerr << delta_seconds(s_time, steady_clock::now()) << endl;

  s_time = steady_clock::now();
  if (VERBOSE) cerr << "[encoding genome]";
  const auto orig_size = orig_genome.size();
  const auto compressed_size = orig_size / 16 + (orig_size % 16 != 0);
  genome.resize(compressed_size);
  encode_dna_four_bit(begin(orig_genome), end(orig_genome), begin(genome));
  vector<uint8_t>().swap(orig_genome);
  if (VERBOSE) cerr << delta_seconds(s_time, steady_clock::now()) << endl;

  // creat genome-wide mask of positions to keep
  keep.clear();
  keep.resize(cl.get_genome_size(), true);

  initialize_bucket_sizes<false>();
  select_two_letter_positions();
  compress_dp();
  initialize_bucket_sizes<true>();
  hash_genome();
  sort_buckets();
}

template<const bool use_mask> void
AbismalIndex::get_bucket_sizes_two() {
  constexpr auto genome_start = seed::padding_size;

  counter_size = (1ull << seed::key_weight);
  counter.clear();
  counter.resize(counter_size + 1, 0);  // one extra entry for convenience

  uint32_t hash_key = 0;
  auto gi = genome_iterator(cbegin(genome)) + genome_start;

  // start spooling the hash key
  const auto gi_lim = gi + (seed::key_weight - 1);
  while (gi != gi_lim) shift_hash_key(*gi++, hash_key);

  auto itl_itr = cbegin(is_two_let) + genome_start;
  auto keep_itr = cbegin(keep) + genome_start;

  // general loop to count positions for corresponding buckets
  const auto lim = cl.get_genome_size() - seed::key_weight - seed::padding_size;
  const auto itl_lim = cbegin(is_two_let) + lim;
  for (; itl_itr != itl_lim; ++itl_itr) {
    shift_hash_key(*gi++, hash_key);
    if (*keep_itr++) counter[hash_key] += (!use_mask || *itl_itr);
  }
}

template<const three_conv_type the_conv, const bool use_mask> void
AbismalIndex::get_bucket_sizes_three() {
  counter_size_three = seed::hash_mask_three;

  if (the_conv == c_to_t)
    counter_t.clear();
  else
    counter_a.clear();

  if (the_conv == c_to_t)
    counter_t.resize(counter_size_three + 1, 0);
  else
    counter_a.resize(counter_size_three + 1, 0);

  constexpr auto genome_start = seed::padding_size;
  const auto lim =
    cl.get_genome_size() - seed::key_weight_three - seed::padding_size;

  uint32_t hash_key = 0;
  auto gi = genome_iterator(cbegin(genome)) + genome_start;

  // start building up the hash key
  const auto gi_lim = gi + (seed::key_weight_three - 1);
  while (gi != gi_lim) shift_three_key<the_conv>(*gi++, hash_key);

  auto itl_itr = cbegin(is_two_let) + genome_start;
  auto keep_itr = cbegin(keep) + genome_start;

  // general loop to count positions for corresponding buckets
  const auto itl_lim = cbegin(is_two_let) + lim;
  for (; itl_itr != itl_lim; ++itl_itr) {
    shift_three_key<the_conv>(*gi++, hash_key);
    if (*keep_itr++) {
      if (the_conv == c_to_t)
        counter_t[hash_key] += !(use_mask && *itl_itr);
      else
        counter_a[hash_key] += !(use_mask && *itl_itr);
    }
  }
}

static inline uint32_t
two_letter_cost(const uint32_t count) {
  return count;
}

static inline uint64_t
three_letter_cost(const uint32_t count_t, const uint32_t count_a) {
  return (count_t + count_a) >> 1;
}

template<const bool use_mask> void
AbismalIndex::initialize_bucket_sizes() {
  auto s_time = steady_clock::now();
  if (VERBOSE) cerr << "[compute bucket sizes]";
#pragma omp parallel sections
  {
#pragma omp section
  get_bucket_sizes_two<use_mask>();
#pragma omp section
  get_bucket_sizes_three<c_to_t, use_mask>();
#pragma omp section
  get_bucket_sizes_three<g_to_a, use_mask>();
  }
  if (VERBOSE) cerr << delta_seconds(s_time, steady_clock::now()) << endl;
}

void
AbismalIndex::select_two_letter_positions() {
  constexpr size_t block_size = 1000000ul;
  auto s_time = steady_clock::now();
  if (VERBOSE) cerr << "[selecting 2-letter positions: choice]";

  // choose which have lower count under two-letters
  // assert(is_two_let.empty());
  is_two_let.resize(cl.get_genome_size(), false);
  const auto itl_beg = begin(is_two_let);
  const auto g_beg = genome_iterator(begin(genome));

  const auto lim = cl.get_genome_size() - seed::key_weight - seed::padding_size;

#pragma omp parallel for
  for (size_t i = seed::padding_size; i < lim; i += block_size) {

    uint32_t hash_two = 0;
    auto gi_two = g_beg + i;

    // spool the two letter hash
    const auto gi_lim_two = gi_two + (seed::key_weight - 1);
    while (gi_two != gi_lim_two) shift_hash_key(*gi_two++, hash_two);

    uint32_t hash_t = 0, hash_a = 0;
    auto gi_three = g_beg + i;

    // spool the three letter hash
    const auto gi_lim_three = gi_three + (seed::key_weight_three - 1);
    while (gi_three != gi_lim_three) {
      const auto x = *gi_three++;
      shift_three_key<c_to_t>(x, hash_t);
      shift_three_key<g_to_a>(x, hash_a);
    }

    // general loop to decide whether a position is represented in the
    // hash using the two or three letter encoding
    const auto itl_lim = itl_beg + min(i + block_size, lim);
    for (auto itl = itl_beg + i; itl != itl_lim; ++itl) {
      shift_hash_key(*gi_two++, hash_two);

      const auto gi3_val = *gi_three++;
      shift_three_key<c_to_t>(gi3_val, hash_t);
      shift_three_key<g_to_a>(gi3_val, hash_a);

      *itl = two_letter_cost(counter[hash_two]) <=
             three_letter_cost(counter_t[hash_t], counter_a[hash_a]);
    }
  }
  if (VERBOSE) cerr << delta_seconds(s_time, steady_clock::now()) << endl;
}

void
AbismalIndex::hash_genome() {
  // count k-mers under each encoding with masking
  auto s_time = steady_clock::now();
  if (VERBOSE) cerr << "[initializing index structures]";
#pragma omp parallel sections
  {
#pragma omp section
  inclusive_scan(begin(counter), end(counter), begin(counter));
#pragma omp section
  inclusive_scan(begin(counter_t), end(counter_t), begin(counter_t));
#pragma omp section
  inclusive_scan(begin(counter_a), end(counter_a), begin(counter_a));
  }

  index_size = counter[counter_size];
  index_size_three = counter_t[counter_size_three];

  index.resize(index_size, 0);
  index_t.resize(index_size_three, 0);
  index_a.resize(index_size_three, 0);
  if (VERBOSE)
    cerr << delta_seconds(s_time, steady_clock::now()) << endl
         << "[index sizes: " << index_size << " " << index_size_three << "]"
         << endl;
  s_time = steady_clock::now();

  if (VERBOSE) cerr << "[counting k-mers]";
  const auto lim = cl.get_genome_size() - seed::key_weight - seed::padding_size;

#pragma omp parallel sections
  {
#pragma omp section
    {
      auto gi_two = genome_iterator(cbegin(genome)) + seed::padding_size;
      const auto gi_lim = gi_two + (seed::key_weight - 1);
      uint32_t hash_two = 0;
      while (gi_two != gi_lim) shift_hash_key(*gi_two++, hash_two);
      for (size_t i = seed::padding_size; i < lim; ++i) {
        shift_hash_key(*gi_two++, hash_two);
        if (keep[i] && is_two_let[i]) index[--counter[hash_two]] = i;
      }
    }
#pragma omp section
    {
      auto gi_three = genome_iterator(cbegin(genome)) + seed::padding_size;
      const auto gi_lim = gi_three + (seed::key_weight_three - 1);
      uint32_t hash_t = 0;
      while (gi_three != gi_lim) shift_three_key<c_to_t>(*gi_three++, hash_t);
      for (size_t i = seed::padding_size; i < lim; ++i) {
        shift_three_key<c_to_t>(*gi_three++, hash_t);
        if (keep[i] && !is_two_let[i]) index_t[--counter_t[hash_t]] = i;
      }
    }
#pragma omp section
    {
      auto gi_three = genome_iterator(begin(genome)) + seed::padding_size;
      const auto gi_lim = gi_three + (seed::key_weight_three - 1);
      uint32_t hash_a = 0;
      while (gi_three != gi_lim) shift_three_key<g_to_a>(*gi_three++, hash_a);
      for (size_t i = seed::padding_size; i < lim; ++i) {
        shift_three_key<g_to_a>(*gi_three++, hash_a);
        if (keep[i] && !is_two_let[i]) index_a[--counter_a[hash_a]] = i;
      }
    }
  }
  if (VERBOSE) cerr << delta_seconds(s_time, steady_clock::now()) << endl;
}

struct dp_sol {
  static const uint64_t sentinel = num_lim<uint64_t>::max();
  uint64_t cost{sentinel};
  uint64_t prev{sentinel};
};

// ADS: Occupancy of the deque can be 0 through window_size, so
// (seed::window_size + 1). However, if we only use this size, then we can't
// tell the difference using `f` and `b` of whether the deque is empty
// or full. We need to ensure that the size of this deque ring buffer
// is at least (window_size + 2).
template<class T> struct fixed_ring_buffer {
  constexpr static size_t qsz = 32ul;  // must be >= (seed::window_size + 2)
  constexpr static size_t qsz_msk = 31ul;  // mask to keep positions in range

  std::array<T, qsz> h;

  auto size() const -> size_t {
    return (b == f)  ? 0 : (b > f) ? b - f : qsz - (f - b);
  }

  auto empty() const -> bool { return f == b; }

  auto back() -> T & { return h[(b - 1) & qsz_msk]; }

  auto front() -> T & { return h[f]; }

  auto pop_back() -> void {
    --b;
    b &= qsz_msk;
  }

  auto pop_front() -> void {
    ++f;
    f &= qsz_msk;
  }

  auto emplace_back(T &&x) -> void {
    h[b] = std::move(x);
    ++b;
    b &= qsz_msk;
  }

  uint32_t f{0};
  uint32_t b{0};
};

using deque = fixed_ring_buffer<dp_sol>;

static inline void
add_sol(deque &helper, const uint64_t prev, const uint64_t cost) {
  while (!helper.empty() && helper.back().cost > cost) helper.pop_back();

  helper.emplace_back({cost, prev});

  while (helper.front().prev + seed::window_size <= prev) helper.pop_front();

  assert(!helper.empty() && helper.size() <= seed::window_size);
}

static inline uint64_t
hybrid_cost(const bool is_two_let, const uint32_t count_two,
            const uint32_t count_t, const uint32_t count_a) {
  return is_two_let ? two_letter_cost(count_two)
                    : three_letter_cost(count_t, count_a);
}

void
AbismalIndex::compress_dp() {
  constexpr auto block_size = 1000000ul;
  auto s_time = steady_clock::now();
  if (VERBOSE) cerr << "[dynamic programming to optimize seed selection]";

  // by default no position is indexed
  std::fill(begin(keep), end(keep), false);

  const auto lim = cl.get_genome_size() - seed::padding_size - seed::key_weight;

#pragma omp parallel for schedule(static, 1)
  for (size_t block_start = seed::padding_size; block_start < lim;
       block_start += block_size) {

    const size_t sz = min(block_start + block_size, lim) - block_start;

    uint32_t hash_two = 0;

    // spool the hash for two letter encoding
    auto gi_two = genome_iterator(begin(genome)) + block_start;
    const auto gi_lim_two(gi_two + (seed::key_weight - 1));
    while (gi_two != gi_lim_two) shift_hash_key(*gi_two++, hash_two);

    uint32_t hash_t = 0, hash_a = 0;

    // spool the hashes for three letter encodings
    auto gi_three = genome_iterator(begin(genome)) + block_start;
    const auto gi_lim_three(gi_three + (seed::key_weight_three - 1));
    while (gi_three != gi_lim_three) {
      const auto gi3_val = *gi_three++;
      shift_three_key<c_to_t>(gi3_val, hash_t);
      shift_three_key<g_to_a>(gi3_val, hash_a);
    }

    vector<dp_sol> opt(block_size + 1);
    deque helper;

    auto opt_itr = begin(opt);
    const auto opt_lim = cbegin(opt) + sz;
    auto itl = cbegin(is_two_let) + block_start;

    // do the base cases for the dynamic programming: the solutions
    // for the first "window size" positions
    size_t i = 0;
    for (; i < seed::window_size; ++opt_itr) {
      shift_hash_key(*gi_two++, hash_two);
      const auto gi3_val = *gi_three++;
      shift_three_key<c_to_t>(gi3_val, hash_t);
      shift_three_key<g_to_a>(gi3_val, hash_a);

      const auto c = hybrid_cost(*itl++, counter[hash_two], counter_t[hash_t],
                                 counter_a[hash_a]);
      *opt_itr = {c, dp_sol::sentinel};
      add_sol(helper, i++, opt_itr->cost);
    }

    // general dynamic programming
    for (; opt_itr != opt_lim; ++opt_itr) {
      shift_hash_key(*gi_two++, hash_two);
      const auto gi3_val = *gi_three++;
      shift_three_key<c_to_t>(gi3_val, hash_t);
      shift_three_key<g_to_a>(gi3_val, hash_a);

      const auto c = hybrid_cost(*itl++, counter[hash_two], counter_t[hash_t],
                                 counter_a[hash_a]);
      *opt_itr = helper.front();
      opt_itr->cost += c;
      add_sol(helper, i++, opt_itr->cost);
    }

    // get the final solution
    uint64_t opt_ans = num_lim<uint64_t>::max();
    uint64_t last = num_lim<uint64_t>::max();
    for (size_t i = sz - 1; i >= sz - seed::window_size; --i) {
      const auto cand_cost = opt[i].cost;
      if (cand_cost < opt_ans) {
        opt_ans = cand_cost;
        last = i;
      }
    }

    *opt_itr = {opt_ans, last};
    const auto block_keep = begin(keep) + block_start;
    const auto opt_beg = begin(opt);
    while (opt_itr->prev != dp_sol::sentinel) {
      *(block_keep + opt_itr->prev) = true;
      opt_itr = opt_beg + opt_itr->prev;
    }
  }

  max_candidates = 100u;  // GS: this is a heuristic
  if (VERBOSE) cerr << delta_seconds(s_time, steady_clock::now()) << endl;
}

struct BucketLess {
  BucketLess(const Genome &g): g_start(begin(g)) {}

  bool operator()(const uint32_t a, const uint32_t b) const {
    auto idx1(g_start + a + seed::key_weight);
    const auto lim1(g_start + a + seed::n_sorting_positions);
    auto idx2(g_start + b + seed::key_weight);
    while (idx1 != lim1) {
      const uint8_t c1 = get_bit(*idx1++);
      const uint8_t c2 = get_bit(*idx2++);
      if (c1 != c2) return c1 < c2;
    }
    return false;
  }

  const genome_iterator g_start;
};

template<const three_conv_type the_conv>
static inline three_letter_t
get_three_letter_num_srt(const uint8_t nt) {
  // C=T=0, A=1, G=4
  // A=G=0, C=2, T=8
  return the_conv == c_to_t ? nt & 5 : nt & 10;
}

template<const three_conv_type the_type> struct BucketLessThree {
  BucketLessThree(const Genome &g): g_start(begin(g)) {}

  bool operator()(const uint32_t a, const uint32_t b) const {
    auto idx1(g_start + a + seed::key_weight_three);
    const auto lim1(g_start + a + seed::n_sorting_positions);
    auto idx2(g_start + b + seed::key_weight_three);
    while (idx1 != lim1) {
      const uint8_t c1 = get_three_letter_num_srt<the_type>(*idx1++);
      const uint8_t c2 = get_three_letter_num_srt<the_type>(*idx2++);
      if (c1 != c2) return c1 < c2;
    }
    return false;
  }

  const genome_iterator g_start;
};

void
AbismalIndex::sort_buckets() {
  {
    auto s_time = steady_clock::now();
    if (VERBOSE) cerr << "[sorting two-letter buckets]";
    const auto b = begin(index);
    const BucketLess bucket_less(genome);
#pragma omp parallel for
    for (size_t i = 0; i < counter_size; ++i)
      if (counter[i + 1] > counter[i] + 1)
        sort(b + counter[i], b + counter[i + 1], bucket_less);
    if (VERBOSE) cerr << delta_seconds(s_time, steady_clock::now()) << endl;
  }
  {
    auto s_time = steady_clock::now();
    if (VERBOSE) cerr << "[sorting three-letter T buckets]";
    const auto b = begin(index_t);
    const BucketLessThree<c_to_t> bucket_less(genome);
#pragma omp parallel for
    for (size_t i = 0; i < counter_size_three; ++i)
      if (counter_t[i + 1] > counter_t[i] + 1)
        sort(b + counter_t[i], b + counter_t[i + 1], bucket_less);
    if (VERBOSE) cerr << delta_seconds(s_time, steady_clock::now()) << endl;
  }
  {
    auto s_time = steady_clock::now();
    if (VERBOSE) cerr << "[sorting three-letter A buckets]";
    const auto b = begin(index_a);
    const BucketLessThree<g_to_a> bucket_less(genome);
#pragma omp parallel for
    for (size_t i = 0; i < counter_size_three; ++i)
      if (counter_a[i + 1] > counter_a[i] + 1)
        sort(b + counter_a[i], b + counter_a[i + 1], bucket_less);
    if (VERBOSE) cerr << delta_seconds(s_time, steady_clock::now()) << endl;
  }
}

static void
write_internal_identifier(FILE *out) {
  if (fwrite((char *)&AbismalIndex::internal_identifier[0], 1,
             AbismalIndex::internal_identifier.size(),
             out) != AbismalIndex::internal_identifier.size())
    throw runtime_error("failed writing index identifier");
}

void
seed::read(FILE *in) {
  static const std::string error_msg("failed to read seed data");
  uint32_t _key_weight;
  uint32_t _window_size;
  uint32_t _n_sorting_positions;

  // key_weight
  if (fread((char *)&_key_weight, sizeof(uint32_t), 1, in) != 1)
    throw runtime_error(error_msg);

  if (_key_weight != key_weight) {
    throw runtime_error(
      "inconsistent k-mer size. Expected: " + to_string(key_weight) +
      ", got: " + to_string(_key_weight));
  }

  // window_size
  if (fread((char *)&_window_size, sizeof(uint32_t), 1, in) != 1)
    throw runtime_error(error_msg);

  if (_window_size != window_size) {
    throw runtime_error(
      "inconsistent window size size. Expected: " + to_string(window_size) +
      ", got: " + to_string(_window_size));
  }

  // n_sorting_positions
  if (fread((char *)&_n_sorting_positions, sizeof(uint32_t), 1, in) != 1)
    throw runtime_error(error_msg);

  if (_n_sorting_positions != n_sorting_positions) {
    throw runtime_error("inconsistent sorting size size. Expected: " +
                        to_string(n_sorting_positions) +
                        ", got: " + to_string(_n_sorting_positions));
  }
}

void
seed::write(FILE *out) {
  static const std::string error_msg("failed to write seed data");
  if (fwrite((char *)&seed::key_weight, sizeof(uint32_t), 1, out) != 1 ||
      fwrite((char *)&seed::window_size, sizeof(uint32_t), 1, out) != 1 ||
      fwrite((char *)&seed::n_sorting_positions, sizeof(uint32_t), 1, out) != 1)
    throw runtime_error(error_msg);
}

void
AbismalIndex::write(const string &index_file) const {
  FILE *out = fopen(index_file.c_str(), "wb");
  if (!out) throw runtime_error("cannot open output file " + index_file);

  write_internal_identifier(out);
  seed::write(out);
  cl.write(out);

  if (fwrite((char *)&genome[0], sizeof(element_t), genome.size(), out) !=
        genome.size() ||
      fwrite((char *)&max_candidates, sizeof(uint32_t), 1, out) != 1 ||
      fwrite((char *)&counter_size, sizeof(size_t), 1, out) != 1 ||
      fwrite((char *)&counter_size_three, sizeof(size_t), 1, out) != 1 ||
      fwrite((char *)&index_size, sizeof(size_t), 1, out) != 1 ||
      fwrite((char *)&index_size_three, sizeof(size_t), 1, out) != 1 ||
      fwrite((char *)(&counter[0]), sizeof(uint32_t), counter_size + 1, out) !=
        (counter_size + 1) ||
      fwrite((char *)(&counter_t[0]), sizeof(uint32_t), counter_size_three + 1,
             out) != (counter_size_three + 1) ||
      fwrite((char *)(&counter_a[0]), sizeof(uint32_t), counter_size_three + 1,
             out) != (counter_size_three + 1) ||

      fwrite((char *)(&index[0]), sizeof(uint32_t), index_size, out) !=
        (index_size) ||
      fwrite((char *)(&index_t[0]), sizeof(uint32_t), index_size_three, out) !=
        (index_size_three) ||
      fwrite((char *)(&index_a[0]), sizeof(uint32_t), index_size_three, out) !=
        (index_size_three))
    throw runtime_error("failed writing index");

  if (fclose(out) != 0)
    throw runtime_error("problem closing file: " + index_file);
}

static bool
check_internal_identifier(FILE *in) {
  string id_found;
  while (id_found.size() < AbismalIndex::internal_identifier.size())
    id_found.push_back(getc(in));

  return (id_found == AbismalIndex::internal_identifier);
}

void
AbismalIndex::read(const string &index_file) {
  static const string error_msg("failed loading index file");

  FILE *in = fopen(index_file.c_str(), "rb");
  if (!in) throw runtime_error("cannot open input file " + index_file);

  if (!check_internal_identifier(in))
    throw runtime_error("index file format problem: " + index_file);

  seed::read(in);
  cl.read(in);

  const size_t genome_to_read = (cl.get_genome_size() + 15) / 16;
  // read the 4-bit encoded genome
  genome.resize(genome_to_read);
  if (fread((char *)&genome[0], sizeof(element_t), genome_to_read, in) !=
      genome_to_read)
    throw runtime_error(error_msg);

  // read the sizes of counter and index vectors
  if (fread((char *)&max_candidates, sizeof(uint32_t), 1, in) != 1 ||
      fread((char *)&counter_size, sizeof(size_t), 1, in) != 1 ||
      fread((char *)&counter_size_three, sizeof(size_t), 1, in) != 1 ||
      fread((char *)&index_size, sizeof(size_t), 1, in) != 1 ||
      fread((char *)&index_size_three, sizeof(size_t), 1, in) != 1)
    throw runtime_error(error_msg);

  // allocate then read the counter vector
  counter = vector<uint32_t>(counter_size + 1);
  if (fread((char *)(&counter[0]), sizeof(uint32_t), (counter_size + 1), in) !=
      (counter_size + 1))
    throw runtime_error(error_msg);

  counter_t = vector<uint32_t>(counter_size_three + 1);
  if (fread((char *)(&counter_t[0]), sizeof(uint32_t), (counter_size_three + 1),
            in) != (counter_size_three + 1))
    throw runtime_error(error_msg);

  counter_a = vector<uint32_t>(counter_size_three + 1);
  if (fread((char *)(&counter_a[0]), sizeof(uint32_t), (counter_size_three + 1),
            in) != (counter_size_three + 1))
    throw runtime_error(error_msg);

  // allocate the read the index vector
  index = vector<uint32_t>(index_size);
  if (fread((char *)(&index[0]), sizeof(uint32_t), index_size, in) !=
      index_size)
    throw runtime_error(error_msg);

  index_t = vector<uint32_t>(index_size_three);
  if (fread((char *)(&index_t[0]), sizeof(uint32_t), index_size_three, in) !=
      index_size_three)
    throw runtime_error(error_msg);

  index_a = vector<uint32_t>(index_size_three);
  if (fread((char *)(&index_a[0]), sizeof(uint32_t), index_size_three, in) !=
      index_size_three)
    throw runtime_error(error_msg);

  if (fclose(in) != 0)
    throw runtime_error("problem closing file: " + index_file);
}

std::ostream &
ChromLookup::write(std::ostream &out) const {
  const uint32_t n_chroms = names.size();
  out.write((char *)&n_chroms, sizeof(uint32_t));
  for (size_t i = 0; i < n_chroms; ++i) {
    const uint32_t name_size = names[i].length();
    out.write((char *)&name_size, sizeof(uint32_t));
    out.write(names[i].c_str(), name_size);
  }
  out.write((char *)(&starts[0]), sizeof(uint32_t) * (n_chroms + 1));
  return out;
}

FILE *
ChromLookup::write(FILE *out) const {
  const uint32_t n_chroms = names.size();
  fwrite((char *)&n_chroms, sizeof(uint32_t), 1, out);
  for (size_t i = 0; i < n_chroms; ++i) {
    const uint32_t name_size = names[i].length();
    fwrite((char *)&name_size, sizeof(uint32_t), 1, out);
    fwrite(names[i].c_str(), 1, name_size, out);
  }
  fwrite((char *)(&starts[0]), sizeof(uint32_t), n_chroms + 1, out);
  return out;
}

void
ChromLookup::write(const string &outfile) const {
  std::ofstream out(outfile, std::ios::binary);
  if (!out) throw runtime_error("cannot open output file " + outfile);
  write(out);
}

std::istream &
ChromLookup::read(std::istream &in) {
  // read the number of chroms
  uint32_t n_chroms = 0;
  in.read((char *)&n_chroms, sizeof(uint32_t));

  // allocate the number of chroms
  names.resize(n_chroms);

  // get each chrom name
  for (size_t i = 0; i < n_chroms; ++i) {
    uint32_t name_size = 0;
    // get the size of the chrom name
    in.read((char *)&name_size, sizeof(uint32_t));
    // allocate the chrom name
    names[i].resize(name_size);
    // read the chrom name
    in.read((char *)&names[i][0], name_size);
  }

  // allocate then read the starts vector
  starts = vector<uint32_t>(n_chroms + 1);
  in.read((char *)(&starts[0]), sizeof(uint32_t) * (n_chroms + 1));

  return in;
}

FILE *
ChromLookup::read(FILE *in) {
  constexpr auto error_msg = "failed loading chrom info from index";

  // read the number of chroms in the reference
  uint32_t n_chroms = 0;
  if (fread((char *)&n_chroms, sizeof(uint32_t), 1, in) != 1)
    throw runtime_error(error_msg);

  names.resize(n_chroms);  // allocate a name for each chrom

  // get each chrom name
  for (size_t i = 0; i < n_chroms; ++i) {
    uint32_t name_size = 0;
    // get size of chrom name
    if (fread((char *)&name_size, sizeof(uint32_t), 1, in) != 1)
      throw runtime_error(error_msg);
    // allocate the chrom name
    names[i].resize(name_size);
    // read the chrom name
    if (fread((char *)&names[i][0], 1, name_size, in) != name_size)
      throw runtime_error(error_msg);
  }

  // allocate then read the starts vector
  starts = vector<uint32_t>(n_chroms + 1);
  if (fread((char *)(&starts[0]), sizeof(uint32_t), n_chroms + 1, in) !=
      n_chroms + 1)
    throw runtime_error(error_msg);

  return in;
}

void
ChromLookup::read(const std::string &infile) {
  std::ifstream in(infile, std::ios::binary);
  if (!in) throw runtime_error("cannot open input file " + infile);
  read(in);
}

std::ostream &
operator<<(std::ostream &out, const ChromLookup &cl) {
  return out << cl.tostring();
}

string
ChromLookup::tostring() const {
  std::ostringstream iss;
  for (size_t i = 0; i < names.size(); ++i)
    iss << i << '\t' << names[i] << '\t' << starts[i + 1] << '\n';
  return iss.str();
}

void
ChromLookup::get_chrom_idx_and_offset(const uint32_t pos, uint32_t &chrom_idx,
                                      uint32_t &offset) const {
  auto idx = upper_bound(cbegin(starts), cend(starts), pos);

  assert(idx != cbegin(starts));

  --idx;

  chrom_idx = std::distance(cbegin(starts), idx);
  offset = pos - starts[chrom_idx];
}

uint32_t
ChromLookup::get_pos(const string &chrom, const uint32_t offset) const {
  const auto itr = find(cbegin(names), cend(names), chrom);
  return itr == cend(names)
           ? num_lim<uint32_t>::max()
           : starts[std::distance(cbegin(names), itr)] + offset;
}

bool
ChromLookup::get_chrom_idx_and_offset(const uint32_t pos,
                                      const uint32_t readlen,
                                      uint32_t &chrom_idx,
                                      uint32_t &offset) const {
  auto idx = upper_bound(cbegin(starts), cend(starts), pos);

  if (idx == cbegin(starts)) return false;  // read is before any chrom

  --idx;

  chrom_idx = std::distance(cbegin(starts), idx);
  offset = pos - starts[chrom_idx];
  return (pos + readlen <= starts[chrom_idx + 1]);
}
