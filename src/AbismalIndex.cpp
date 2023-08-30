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

#include <omp.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <deque>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <utility>

#include "dna_four_bit.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

using std::clog;
using std::cerr;
using std::endl;
using std::ifstream;
using std::inclusive_scan;
using std::min;
using std::ofstream;
using std::pair;
using std::runtime_error;
using std::sort;
using std::string;
using std::to_string;
using std::vector;
using std::begin;

using std::chrono::duration;
using std::chrono::steady_clock;
using std::chrono::time_point;

template<typename T> using num_lim = std::numeric_limits<T>;

bool AbismalIndex::VERBOSE = false;
const string AbismalIndex::internal_identifier = "AbismalIndex";

using genome_iterator = genome_four_bit_itr;

static string
delta_seconds(const time_point<steady_clock> &a) {
  const auto b = steady_clock::now();
  std::ostringstream oss;
  oss.setf(std::ios::fixed);
  oss.precision(2);
  oss << "[time: " << duration<double>(b - a).count() << "s]";
  return oss.str();
}

struct g_interval {
  string chrom{};
  size_t start_pos{};
  size_t end_pos{};

  g_interval(const string &c, size_t s, size_t e) :
    chrom{c}, start_pos{s}, end_pos{e} {}
  g_interval() {}
  bool operator<(const g_interval &rhs) const {
    const int x = chrom.compare(rhs.chrom);
    return (x < 0 || (x == 0 && (start_pos < rhs.start_pos ||
                                 (start_pos == rhs.start_pos &&
                                  (end_pos == rhs.end_pos)))));
  }
};

template <class T> T&
operator<<(T &out, const g_interval &g) {
  return out << g.chrom << '\t' << g.start_pos << '\t' << g.end_pos;
}

auto
load_target_regions(const string &target_filename) -> vector<g_interval> {
  constexpr auto msg = "failed parsing target region";
  ifstream in(target_filename);
  if (!in) throw runtime_error("failed reading target file");

  vector<g_interval> targets{};
  string line;
  string chrom;
  size_t start_pos = 0, end_pos = 0;
  while (getline(in, line)) {

    std::istringstream iss(line);
    if (!(iss >> chrom)) throw runtime_error(msg);
    if (!(iss >> start_pos)) throw runtime_error(msg);
    if (!(iss >> end_pos)) throw runtime_error(msg);

    targets.emplace_back(chrom, start_pos, end_pos);
  }

  return targets;
}

auto
mask_non_target(const vector<pair<uint32_t, uint32_t>> &targets,
                vector<uint8_t> &genome) {
  const auto target_end = cend(targets);
  auto target_itr = cbegin(targets);
  for (size_t i = 0; i < genome.size(); ++i) {
    if (target_itr == target_end || i < target_itr->first)
      genome[i] = 'N';
    else
      while (target_itr != target_end && target_itr->second <= i)
        ++target_itr;
  }
}

auto
contiguous_n(const vector<uint8_t> &genome) -> vector<pair<size_t, size_t>> {
  using std::distance;
  vector<pair<size_t, size_t>> r{};
  const auto g_beg = cbegin(genome);
  auto g_left = g_beg;
  bool inside_block = false;
  for (auto g = g_beg; g != cend(genome); ++g) {
    if (*g == 'N' && !inside_block) {
      g_left = g;
      inside_block = true;
    }
    else if (*g != 'N' && inside_block) {
      r.emplace_back(distance(g_beg, g_left), distance(g_beg, g));
      inside_block = false;
    }
  }
  if (inside_block) r.emplace_back(distance(g_beg, g_left), genome.size());
  return r;
}

template<typename G> static inline auto
get_exclude_itrs(const G &genome, const vector<pair<size_t, size_t>> &exclude)
  -> vector<pair<genome_iterator, genome_iterator>> {
  vector<pair<genome_iterator, genome_iterator>> exclude_itrs;
  const auto g_beg = genome_iterator(cbegin(genome));
  for (auto &i : exclude)
    exclude_itrs.emplace_back(g_beg + i.first, g_beg + i.second);
  return exclude_itrs;
}

template<typename G> static inline void
replace_included_n(const vector<pair<size_t, size_t>> &exclude, G &genome) {
  size_t j = 0;
  for (size_t i = 0; i < genome.size(); ++i) {
    if (genome[i] == 'N' && i < exclude[j].first) genome[i] = random_base();
    if (exclude[j].second <= i) ++j;
  }
}

static inline size_t
get_compressed_size(const size_t original_size) {
  const size_t nt_per_word = 2 * sizeof(original_size);
  return original_size / nt_per_word + (original_size % nt_per_word != 0);
}

static auto
sort_by_chrom(const vector<string> &names, const vector<g_interval> &u)
  -> vector<g_interval> {
  vector<g_interval> t(u.size());
  auto t_itr = begin(t);
  for (size_t i = 0; i < names.size(); ++i) {
    const auto u_itr = find_if(cbegin(u), cend(u), [&](const g_interval &x) {
      return x.chrom == names[i];
    });
    const auto v_itr = find_if(
      u_itr, cend(u), [&](const g_interval &x) { return x.chrom != names[i]; });

    if (!std::is_sorted(u_itr, v_itr))
      throw runtime_error("target regions not sorted");

    std::copy(u_itr, v_itr, t_itr);
    t_itr += std::distance(u_itr, v_itr);
  }
  if (t_itr != end(t))
    throw runtime_error("target regions unsorted or include extra chromosomes");
  return t;
}

void
AbismalIndex::create_index(const string &targets_file,
                           const string &genome_file) {

  auto s_time = steady_clock::now();
  if (VERBOSE) clog << "[loading genome]";
  vector<uint8_t> orig_genome;
  load_genome(genome_file, orig_genome, cl);
  if (VERBOSE) clog << delta_seconds(s_time) << endl;

  s_time = steady_clock::now();
  if (VERBOSE) clog << "[loading target regions]";
  vector<g_interval> orig_targets = load_target_regions(targets_file);
  orig_targets = sort_by_chrom(cl.names, orig_targets);
  if (VERBOSE) clog << delta_seconds(s_time) << endl;

  s_time = steady_clock::now();
  if (VERBOSE) clog << "[cleaning reference genome]";

  vector<pair<uint32_t, uint32_t>> targets;
  for (auto &i : orig_targets)
    targets.emplace_back(cl.get_pos(i.chrom, i.start_pos),
                         cl.get_pos(i.chrom, i.end_pos));

  mask_non_target(targets, orig_genome);

  exclude = contiguous_n(orig_genome);

  auto rm_itr = std::remove_if(begin(exclude), end(exclude),
                               [](const pair<size_t, size_t> &p) {
                                 return p.second - p.first <= max_n_count;
                               });
  exclude.erase(rm_itr, cend(exclude));

  for (size_t i = 0; i < exclude.size(); ++i)
    if (exclude[i].second - exclude[i].first <= max_n_count) {
      clog << "wtf" << endl;
      exit(1);
    }
  replace_included_n(exclude, orig_genome);
  if (VERBOSE) clog << delta_seconds(s_time) << endl;

  s_time = steady_clock::now();

  if (VERBOSE) clog << "[encoding genome]";
  genome.resize(get_compressed_size(orig_genome.size()));
  encode_dna_four_bit(begin(orig_genome), end(orig_genome), begin(genome));
  vector<uint8_t>().swap(orig_genome);
  if (VERBOSE) clog << delta_seconds(s_time) << endl;

  exclude_itr = get_exclude_itrs(genome, exclude);

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

void
AbismalIndex::create_index(const string &genome_file) {

  auto s_time = steady_clock::now();
  if (VERBOSE) clog << "[loading genome]";
  vector<uint8_t> orig_genome;
  load_genome(genome_file, orig_genome, cl);
  if (VERBOSE) clog << delta_seconds(s_time) << endl;

  s_time = steady_clock::now();
  if (VERBOSE) clog << "[cleaning reference genome]";
  exclude = contiguous_n(orig_genome);

  auto rm_itr = std::remove_if(begin(exclude), end(exclude),
                               [](const pair<size_t, size_t> &p) {
                                 return p.second - p.first <= max_n_count;
                               });
  exclude.erase(rm_itr, cend(exclude));

  replace_included_n(exclude, orig_genome);
  if (VERBOSE) clog << delta_seconds(s_time) << endl;

  s_time = steady_clock::now();

  if (VERBOSE) clog << "[encoding genome]";
  genome.resize(get_compressed_size(orig_genome.size()));
  encode_dna_four_bit(begin(orig_genome), end(orig_genome), begin(genome));
  vector<uint8_t>().swap(orig_genome);
  if (VERBOSE) clog << delta_seconds(s_time) << endl;

  exclude_itr = get_exclude_itrs(genome, exclude);

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
  constexpr auto genome_start = 0ul; // seed::padding_size;

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
  auto nidx_itr = cbegin(exclude);

  // general loop to count positions for corresponding buckets
  const auto lim = cl.get_genome_size() - seed::key_weight - genome_start;
  // const auto itl_lim = cbegin(is_two_let) + lim;
  for (size_t i = 0; i < lim; ++i) {
    shift_hash_key(*gi++, hash_key);
    assert(nidx_itr != cend(exclude));
    if (i < nidx_itr->first && *keep_itr)
      counter[hash_key] += (!use_mask || *itl_itr);
    if (nidx_itr->second <= i) ++nidx_itr;
    ++keep_itr;
    ++itl_itr;
  }
}

template<const three_conv_type the_conv, const bool use_mask> void
AbismalIndex::get_bucket_sizes_three() {
  constexpr auto genome_start = 0ul; // seed::padding_size;

  counter_size_three = seed::hash_mask_three;

  if (the_conv == c_to_t)
    counter_t.clear();
  else
    counter_a.clear();

  if (the_conv == c_to_t)
    counter_t.resize(counter_size_three + 1, 0);
  else
    counter_a.resize(counter_size_three + 1, 0);

  const auto lim = cl.get_genome_size() - seed::key_weight_three - genome_start;

  uint32_t hash_key = 0;
  auto gi = genome_iterator(cbegin(genome)) + genome_start;

  // start building up the hash key
  const auto gi_lim = gi + (seed::key_weight_three - 1);
  while (gi != gi_lim) shift_three_key<the_conv>(*gi++, hash_key);

  auto itl_itr = cbegin(is_two_let) + genome_start;
  auto keep_itr = cbegin(keep) + genome_start;
  auto nidx_itr = cbegin(exclude);

  // general loop to count positions for corresponding buckets
  for (size_t i = genome_start; i < lim; ++i) {
    shift_three_key<the_conv>(*gi++, hash_key);
    assert(nidx_itr != cend(exclude));
    if (i < nidx_itr->first && *keep_itr) {
      if (the_conv == c_to_t)
        counter_t[hash_key] += (!use_mask || !(*itl_itr));
      else
        counter_a[hash_key] += (!use_mask || !(*itl_itr));
    }
    if (nidx_itr->second <= i) ++nidx_itr;
    ++keep_itr;
    ++itl_itr;
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
  if (VERBOSE) clog << "[computing bucket sizes]";
#pragma omp parallel sections
  {
#pragma omp section
    get_bucket_sizes_two<use_mask>();
#pragma omp section
    get_bucket_sizes_three<c_to_t, use_mask>();
#pragma omp section
    get_bucket_sizes_three<g_to_a, use_mask>();
  }
  if (VERBOSE) clog << delta_seconds(s_time) << endl;
}

static inline auto
get_block_bounds(const size_t start_pos, const size_t step_size,
                 const size_t end_pos,
                 const vector<pair<size_t, size_t>> &exclude)
  -> vector<pair<size_t, size_t>> {
  vector<pair<size_t, size_t>> blocks;

  auto block_start = start_pos;
  auto i = cbegin(exclude);
  while (block_start < end_pos && i != cend(exclude)) {
    if (block_start < i->first) {
      const auto block_end = min({i->first, block_start + step_size, end_pos});
      blocks.emplace_back(block_start, block_end);
      block_start += step_size;

      // this might move block start backward
      if (block_start >= i->second) block_start = (*i++).second;
    }
    else
      block_start = (*i++).second;
  }
  while (block_start < end_pos) {
    const auto block_end = min({block_start + step_size, end_pos});
    blocks.emplace_back(block_start, block_end);
    block_start += step_size;
  }

  // for (auto &i : blocks) {
  //   if (i.second - i.first < AbismalIndex::max_n_count) {
  //     cerr << i.first << '\t' << i.second << endl;
  //     exit(1);
  //   }
  // }

  return blocks;
}

void
AbismalIndex::select_two_letter_positions() {
  constexpr auto genome_start = 0ul;
  constexpr size_t block_size = 1000000ul;
  auto s_time = steady_clock::now();
  if (VERBOSE) clog << "[selecting two-letter positions]";

  // choose which have lower count under two-letters
  is_two_let.resize(cl.get_genome_size(), false);
  const auto itl_beg = begin(is_two_let);
  const auto g_beg = genome_iterator(begin(genome));

  const auto lim = cl.get_genome_size() - seed::key_weight - genome_start;

  const auto blocks = get_block_bounds(genome_start, block_size, lim, exclude);

#pragma omp parallel for
  for (auto &block : blocks) {

    uint32_t hash_two = 0;
    auto gi_two = g_beg + block.first;

    // spool the two letter hash
    const auto gi_lim_two = gi_two + (seed::key_weight - 1);
    while (gi_two != gi_lim_two) shift_hash_key(*gi_two++, hash_two);

    uint32_t hash_t = 0, hash_a = 0;
    auto gi_three = g_beg + block.first;

    // spool the three letter hash
    const auto gi_lim_three = gi_three + (seed::key_weight_three - 1);
    while (gi_three != gi_lim_three) {
      const auto x = *gi_three++;
      shift_three_key<c_to_t>(x, hash_t);
      shift_three_key<g_to_a>(x, hash_a);
    }

    // general loop to decide whether a position is represented in the
    // hash using the two or three letter encoding
    const auto itl_lim = itl_beg + block.second;
    for (auto itl = itl_beg + block.first; itl != itl_lim; ++itl) {
      shift_hash_key(*gi_two++, hash_two);

      const auto gi3_val = *gi_three++;
      shift_three_key<c_to_t>(gi3_val, hash_t);
      shift_three_key<g_to_a>(gi3_val, hash_a);

      *itl = two_letter_cost(counter[hash_two]) <=
             three_letter_cost(counter_t[hash_t], counter_a[hash_a]);
    }
  }
  if (VERBOSE) clog << delta_seconds(s_time) << endl;
}

void
AbismalIndex::hash_genome() {
  constexpr auto genome_start = 0ul; // seed::padding_size;
  // count k-mers under each encoding with masking
  auto s_time = steady_clock::now();
  if (VERBOSE) clog << "[initializing index structures]";
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
    clog << delta_seconds(s_time) << endl
         << "[index sizes: two-letter=" << index_size << " "
         << "three-letter=" << index_size_three << "]" << endl;
  s_time = steady_clock::now();

  if (VERBOSE) clog << "[counting k-mers]";
  const auto lim = cl.get_genome_size() - seed::key_weight - genome_start;

#pragma omp parallel sections
  {
#pragma omp section
    {
      auto gi_two = genome_iterator(cbegin(genome)) + genome_start;
      const auto gi_lim = gi_two + (seed::key_weight - 1);
      uint32_t hash_two = 0;
      while (gi_two != gi_lim) shift_hash_key(*gi_two++, hash_two);
      auto nidx_itr = cbegin(exclude);
      for (size_t i = genome_start; i < lim; ++i) {
        shift_hash_key(*gi_two++, hash_two);
        assert(nidx_itr != cend(exclude));
        if (i < nidx_itr->first && keep[i] && is_two_let[i]) {
          assert(counter[hash_two] > 0);
          index[--counter[hash_two]] = i;
        }
        if (nidx_itr->second <= i) ++nidx_itr;
      }
    }
#pragma omp section
    {
      auto gi_three = genome_iterator(cbegin(genome)) + genome_start;
      const auto gi_lim = gi_three + (seed::key_weight_three - 1);
      uint32_t hash_t = 0;
      while (gi_three != gi_lim) shift_three_key<c_to_t>(*gi_three++, hash_t);
      auto nidx_itr = cbegin(exclude);
      for (size_t i = genome_start; i < lim; ++i) {
        shift_three_key<c_to_t>(*gi_three++, hash_t);
        assert(nidx_itr != cend(exclude));
        if (i < nidx_itr->first && keep[i] && !is_two_let[i]) {
          assert(counter_t[hash_t] > 0);
          index_t[--counter_t[hash_t]] = i;
        }
        if (nidx_itr->second <= i) ++nidx_itr;
      }
    }
#pragma omp section
    {
      auto gi_three = genome_iterator(begin(genome)) + genome_start;
      const auto gi_lim = gi_three + (seed::key_weight_three - 1);
      uint32_t hash_a = 0;
      while (gi_three != gi_lim) shift_three_key<g_to_a>(*gi_three++, hash_a);
      auto nidx_itr = cbegin(exclude);
      for (size_t i = genome_start; i < lim; ++i) {
        shift_three_key<g_to_a>(*gi_three++, hash_a);
        assert(nidx_itr != cend(exclude));
        if (i < nidx_itr->first && keep[i] && !is_two_let[i]) {
          assert(counter_a[hash_a] > 0);
          index_a[--counter_a[hash_a]] = i;
        }
        if (nidx_itr->second <= i) ++nidx_itr;
      }
    }
  }
  if (VERBOSE) clog << delta_seconds(s_time) << endl;
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

  auto push_back(T &&x) -> void {
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

  helper.push_back({cost, prev});

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
  constexpr auto genome_start = 0ul; // seed::padding_size;
  constexpr auto block_size = 1000000ul;
  auto s_time = steady_clock::now();
  if (VERBOSE) clog << "[dynamic programming to optimize seed selection]";

  // default is no position will be indexed
  std::fill(begin(keep), end(keep), false);

  const auto lim = cl.get_genome_size() - seed::key_weight - genome_start;

  const auto blocks = get_block_bounds(genome_start, block_size, lim, exclude);

#pragma omp parallel for
  for (auto &block : blocks) {

    const size_t block_start = block.first;
    const uint32_t current_block_size = block.second - block.first;

    if (current_block_size < seed::window_size) continue;

    uint32_t hash_two = 0;

    // spool the hash for two letter encoding
    auto gi_two = genome_iterator(begin(genome)) + block_start;
    const auto gi_lim_two =
      gi_two + min(current_block_size, seed::key_weight - 1);
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
    const auto opt_lim = cbegin(opt) + current_block_size;
    auto itl = cbegin(is_two_let) + block_start;

    // do the base cases for the dynamic programming: the solutions
    // for the first "window size" positions
    size_t i = 0;
    for (; i < seed::window_size; ++opt_itr) {
      shift_hash_key(*gi_two++, hash_two);
      const auto gi3_val = *gi_three++;
      shift_three_key<c_to_t>(gi3_val, hash_t);
      shift_three_key<g_to_a>(gi3_val, hash_a);

      assert(hash_two < counter.size());
      assert(hash_t < counter_t.size());
      assert(hash_a < counter_a.size());

      const auto c = hybrid_cost(*itl++, counter[hash_two], counter_t[hash_t],
                                 counter_a[hash_a]);
      assert(opt_itr != end(opt));
      *opt_itr = {c, dp_sol::sentinel};
      add_sol(helper, i++, opt_itr->cost);
    }

    // general dynamic programming
    for (; opt_itr != opt_lim; ++opt_itr) {
      shift_hash_key(*gi_two++, hash_two);
      const auto gi3_val = *gi_three++;
      shift_three_key<c_to_t>(gi3_val, hash_t);
      shift_three_key<g_to_a>(gi3_val, hash_a);

      assert(hash_two < counter.size());
      assert(hash_t < counter_t.size());
      assert(hash_a < counter_a.size());

      const auto c = hybrid_cost(*itl++, counter[hash_two], counter_t[hash_t],
                                 counter_a[hash_a]);
      assert(helper.size() > 0);
      assert(opt_itr != end(opt));
      *opt_itr = helper.front();
      opt_itr->cost += c;
      add_sol(helper, i++, opt_itr->cost);
    }

    // get solution for start of traceback
    uint64_t opt_ans = num_lim<uint64_t>::max();
    uint64_t last = num_lim<uint64_t>::max();
    for (size_t i = current_block_size - 1;
         i >= current_block_size - seed::window_size; --i) {
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
  if (VERBOSE) clog << delta_seconds(s_time) << endl;
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
three_letter_num_srt(const uint8_t nt) {
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
      const uint8_t c1 = three_letter_num_srt<the_type>(*idx1++);
      const uint8_t c2 = three_letter_num_srt<the_type>(*idx2++);
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
    if (VERBOSE) clog << "[sorting two-letter buckets]";
    const auto b = begin(index);
    const BucketLess bucket_less(genome);
#pragma omp parallel for
    for (size_t i = 0; i < counter_size; ++i)
      if (counter[i + 1] > counter[i] + 1)
        sort(b + counter[i], b + counter[i + 1], bucket_less);
    if (VERBOSE) clog << delta_seconds(s_time) << endl;
  }
  {
    auto s_time = steady_clock::now();
    if (VERBOSE) clog << "[sorting three-letter T buckets]";
    const auto b = begin(index_t);
    const BucketLessThree<c_to_t> bucket_less(genome);
#pragma omp parallel for
    for (size_t i = 0; i < counter_size_three; ++i)
      if (counter_t[i + 1] > counter_t[i] + 1)
        sort(b + counter_t[i], b + counter_t[i + 1], bucket_less);
    if (VERBOSE) clog << delta_seconds(s_time) << endl;
  }
  {
    auto s_time = steady_clock::now();
    if (VERBOSE) clog << "[sorting three-letter A buckets]";
    const auto b = begin(index_a);
    const BucketLessThree<g_to_a> bucket_less(genome);
#pragma omp parallel for
    for (size_t i = 0; i < counter_size_three; ++i)
      if (counter_a[i + 1] > counter_a[i] + 1)
        sort(b + counter_a[i], b + counter_a[i + 1], bucket_less);
    if (VERBOSE) clog << delta_seconds(s_time) << endl;
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
  uint32_t _key_weight = 0;
  uint32_t _window_size = 0;
  uint32_t _n_sorting_positions = 0;

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
  for (auto i = 0u; i < names.size(); ++i)
    iss << i << '\t'
        << names[i] << '\t'
        << starts[i] << '\t'
        << starts[i + 1] << '\n';
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

template<class G> void
load_genome_impl(const string &genome_file, G &genome, ChromLookup &cl) {
  std::ifstream in(genome_file);
  if (!in) throw std::runtime_error("bad genome file: " + genome_file);

  namespace fs = std::filesystem;
  const size_t file_size = fs::file_size(fs::path(genome_file));

  genome.clear();
  // reserve space for padding at both ends
  genome.reserve(file_size + 2*seed::padding_size);
  auto g_ins = std::back_inserter(genome);

  // pad the start of the concatenated sequence
  cl.names.push_back("pad_start");
  cl.starts.push_back(0);
  fill_n(g_ins, seed::padding_size, 'N');

  std::string line;
  while (getline(in, line))
    if (line[0] != '>')
      copy(std::cbegin(line), std::cend(line), g_ins);
    else {
      cl.names.push_back(line.substr(1, line.find_first_of(" \t") - 1));
      cl.starts.push_back(genome.size());
    }

  if (cl.names.size() < 2)
    throw std::runtime_error("no names found in genome file");

  // now pad the end of the concatenated sequence
  cl.names.push_back("pad_end");
  cl.starts.push_back(genome.size());
  std::fill_n(g_ins, seed::padding_size, 'N');

  // this one additional "start" is the end of all chroms
  cl.starts.push_back(genome.size());
}

void
load_genome(const string &genome_file, string &genome, ChromLookup &cl) {
  load_genome_impl(genome_file, genome, cl);
}

void
load_genome(const string &genome_file, vector<uint8_t> &genome,
            ChromLookup &cl) {
  load_genome_impl(genome_file, genome, cl);
}
