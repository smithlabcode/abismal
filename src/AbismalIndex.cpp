/* Copyright (C) 2018-2025 Andrew D. Smith and Guilherme Sena
 *
 * Authors: Andrew D. Smith and Guilherme Sena
 *
 * This file is part of abismal.
 *
 * abismal is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * abismal is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 */

#include "AbismalIndex.hpp"
#include "dna_four_bit.hpp"

#include "bamxx.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <utility>

// NOLINTBEGIN

using abismal_clock = std::chrono::steady_clock;
using std::chrono::time_point;

template <typename T> using num_lim = std::numeric_limits<T>;

bool AbismalIndex::VERBOSE = false;
std::size_t AbismalIndex::n_threads_global = 1;
const std::string AbismalIndex::internal_identifier = "AbismalIndex";

using genome_iterator = genome_four_bit_itr;

static std::string
delta_seconds(const time_point<abismal_clock> &a) {
  const auto b = abismal_clock::now();
  std::ostringstream oss;
  oss.setf(std::ios::fixed);
  oss.precision(2);
  oss << "[time: " << std::chrono::duration<double>(b - a).count() << "s]";
  return oss.str();
}

struct g_interval {
  std::string chrom{};
  std::size_t start_pos{};
  std::size_t end_pos{};

  g_interval(const std::string &c, std::size_t s, std::size_t e) :
    chrom{c}, start_pos{s}, end_pos{e} {}
  g_interval() {}
  bool
  operator<(const g_interval &rhs) const {
    const int x = chrom.compare(rhs.chrom);
    return (x < 0 || (x == 0 && (start_pos < rhs.start_pos ||
                                 (start_pos == rhs.start_pos &&
                                  (end_pos < rhs.end_pos)))));
  }
};

template <class T>
T &
operator<<(T &out, const g_interval &g) {
  return out << g.chrom << '\t' << g.start_pos << '\t' << g.end_pos;
}

auto
load_target_regions(const std::string &target_filename)
  -> std::vector<g_interval> {
  constexpr auto msg = "failed parsing target region";
  std::ifstream in(target_filename);
  if (!in)
    throw std::runtime_error("failed reading target file");

  std::vector<g_interval> targets{};
  std::string line;
  std::string chrom;
  std::size_t start_pos = 0, end_pos = 0;
  while (getline(in, line)) {
    std::istringstream iss(line);
    if (!(iss >> chrom))
      throw std::runtime_error(msg);
    if (!(iss >> start_pos))
      throw std::runtime_error(msg);
    if (!(iss >> end_pos))
      throw std::runtime_error(msg);
    targets.emplace_back(chrom, start_pos, end_pos);
  }

  return targets;
}

auto
mask_non_target(
  const std::vector<std::pair<std::uint32_t, std::uint32_t>> &targets,
  std::vector<std::uint8_t> &genome) {
  const auto target_end = std::cend(targets);
  auto target_itr = std::cbegin(targets);
  for (std::size_t i = 0; i < std::size(genome); ++i) {
    if (target_itr == target_end || i < target_itr->first)
      genome[i] = 'N';
    else
      while (target_itr != target_end && target_itr->second <= i)
        ++target_itr;
  }
}

auto
contiguous_n(const std::vector<std::uint8_t> &genome)
  -> std::vector<std::pair<std::size_t, std::size_t>> {
  std::vector<std::pair<std::size_t, std::size_t>> r{};
  const auto g_beg = std::cbegin(genome);
  auto g_left = g_beg;
  bool inside_block = false;
  for (auto g = g_beg; g != std::cend(genome); ++g) {
    if (*g == 'N' && !inside_block) {
      g_left = g;
      inside_block = true;
    }
    else if (*g != 'N' && inside_block) {
      r.emplace_back(std::distance(g_beg, g_left), std::distance(g_beg, g));
      inside_block = false;
    }
  }
  if (inside_block)
    r.emplace_back(std::distance(g_beg, g_left), genome.size());
  return r;
}

template <typename G>
static inline auto
get_exclude_itrs(
  const G &genome,
  const std::vector<std::pair<std::size_t, std::size_t>> &exclude)
  -> std::vector<std::pair<genome_iterator, genome_iterator>> {
  std::vector<std::pair<genome_iterator, genome_iterator>> exclude_itrs;
  const auto g_beg = genome_iterator(std::cbegin(genome));
  for (auto &i : exclude)
    exclude_itrs.emplace_back(g_beg + i.first, g_beg + i.second);
  return exclude_itrs;
}

template <typename G>
static inline void
replace_included_n(
  const std::vector<std::pair<std::size_t, std::size_t>> &exclude, G &genome) {
  std::size_t j = 0;
  for (std::size_t i = 0; i < std::size(genome); ++i) {
    if (genome[i] == 'N' && i < exclude[j].first)
      genome[i] = random_base();
    if (exclude[j].second <= i)
      ++j;
  }
}

static inline std::size_t
get_compressed_size(const std::size_t original_size) {
  const std::size_t nt_per_word = 2 * sizeof(original_size);
  return original_size / nt_per_word + (original_size % nt_per_word != 0);
}

static auto
sort_by_chrom(const std::vector<std::string> &names,
              const std::vector<g_interval> &u) -> std::vector<g_interval> {
  // This function will exclude any target regions that are not
  // corresponding to chromosomes in the given reference genome.
  std::vector<g_interval> t(u.size());
  auto t_itr = std::begin(t);
  for (std::size_t i = 0; i < names.size(); ++i) {
    const auto u_itr =
      find_if(std::cbegin(u), std::cend(u),
              [&](const g_interval &x) { return x.chrom == names[i]; });
    const auto v_itr = find_if(u_itr, std::cend(u), [&](const g_interval &x) {
      return x.chrom != names[i];
    });

    if (!std::is_sorted(u_itr, v_itr))
      throw std::runtime_error("target regions not sorted");

    std::copy(u_itr, v_itr, t_itr);
    t_itr += std::distance(u_itr, v_itr);
  }
  t.resize(std::distance(std::begin(t), t_itr));
  return t;
}

void
AbismalIndex::create_index(const std::string &targets_file,
                           const std::string &genome_file) {

  auto s_time = abismal_clock::now();
  if (VERBOSE)
    std::clog << "[loading genome]";
  std::vector<std::uint8_t> orig_genome;
  load_genome(genome_file, orig_genome, cl);
  if (VERBOSE)
    std::clog << delta_seconds(s_time) << '\n';

  s_time = abismal_clock::now();
  if (VERBOSE)
    std::clog << "[loading target regions]";
  std::vector<g_interval> orig_targets = load_target_regions(targets_file);
  orig_targets = sort_by_chrom(cl.names, orig_targets);
  if (VERBOSE)
    std::clog << delta_seconds(s_time) << '\n';

  s_time = abismal_clock::now();
  if (VERBOSE)
    std::clog << "[cleaning reference genome]";

  std::vector<std::pair<std::uint32_t, std::uint32_t>> targets;
  for (auto &i : orig_targets)
    targets.emplace_back(cl.get_pos(i.chrom, i.start_pos),
                         cl.get_pos(i.chrom, i.end_pos));

  mask_non_target(targets, orig_genome);

  exclude = contiguous_n(orig_genome);

  auto rm_itr =
    std::remove_if(std::begin(exclude), std::end(exclude),
                   [](const std::pair<std::size_t, std::size_t> &p) {
                     return p.second - p.first <= max_n_count;
                   });
  exclude.erase(rm_itr, std::cend(exclude));

  replace_included_n(exclude, orig_genome);
  if (VERBOSE)
    std::clog << delta_seconds(s_time) << '\n';

  s_time = abismal_clock::now();

  if (VERBOSE)
    std::clog << "[encoding genome]";
  genome.resize(get_compressed_size(orig_genome.size()));
  encode_dna_four_bit(std::cbegin(orig_genome), std::cend(orig_genome),
                      std::begin(genome));
  std::vector<std::uint8_t>().swap(orig_genome);
  if (VERBOSE)
    std::clog << delta_seconds(s_time) << '\n';

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
AbismalIndex::create_index(const std::string &genome_file) {

  auto s_time = abismal_clock::now();
  if (VERBOSE)
    std::clog << "[loading genome]";
  std::vector<std::uint8_t> orig_genome;
  load_genome(genome_file, orig_genome, cl);
  if (VERBOSE)
    std::clog << "[" << cl.names.size() << " targets]" << delta_seconds(s_time)
              << '\n';

  s_time = abismal_clock::now();
  if (VERBOSE)
    std::clog << "[cleaning reference genome]";
  exclude = contiguous_n(orig_genome);

  auto rm_itr =
    std::remove_if(std::begin(exclude), std::end(exclude),
                   [](const std::pair<std::size_t, std::size_t> &p) {
                     return p.second - p.first <= max_n_count;
                   });
  exclude.erase(rm_itr, std::cend(exclude));

  replace_included_n(exclude, orig_genome);
  if (VERBOSE)
    std::clog << delta_seconds(s_time) << '\n';

  s_time = abismal_clock::now();

  if (VERBOSE)
    std::clog << "[encoding genome]";
  genome.resize(get_compressed_size(orig_genome.size()));
  encode_dna_four_bit(std::cbegin(orig_genome), std::cend(orig_genome),
                      std::begin(genome));
  std::vector<std::uint8_t>().swap(orig_genome);
  if (VERBOSE)
    std::clog << delta_seconds(s_time) << '\n';

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

template <const bool use_mask>
void
AbismalIndex::get_bucket_sizes_two() {

  counter_size = (1ull << seed::key_weight);
  counter.clear();
  counter.resize(counter_size + 1, 0);  // one extra entry for convenience

  std::uint32_t hash_key = 0;
  auto gi = genome_iterator(std::cbegin(genome));

  // start spooling the hash key
  const auto gi_lim = gi + (seed::key_weight - 1);
  while (gi != gi_lim)
    shift_hash_key(*gi++, hash_key);

  auto itl_itr = std::cbegin(is_two_let);
  auto keep_itr = std::cbegin(keep);
  auto nidx_itr = std::cbegin(exclude);

  // general loop to count positions for corresponding buckets
  const auto lim = cl.get_genome_size() - seed::key_weight + 1;
  // const auto itl_lim = std::cbegin(is_two_let) + lim;
  for (std::size_t i = 0; i < lim; ++i) {
    shift_hash_key(*gi++, hash_key);
    assert(nidx_itr != std::cend(exclude));
    if (i < nidx_itr->first && *keep_itr)
      counter[hash_key] += (!use_mask || *itl_itr);
    if (nidx_itr->second <= i)
      ++nidx_itr;
    ++keep_itr;
    ++itl_itr;
  }
}

template <const three_conv_type the_conv, const bool use_mask>
void
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

  std::uint32_t hash_key = 0;
  auto gi = genome_iterator(std::cbegin(genome));

  // start building up the hash key
  const auto gi_lim = gi + (seed::key_weight_three - 1);
  while (gi != gi_lim)
    shift_three_key<the_conv>(*gi++, hash_key);

  auto itl_itr = std::cbegin(is_two_let);
  auto keep_itr = std::cbegin(keep);
  auto nidx_itr = std::cbegin(exclude);

  // general loop to count positions for corresponding buckets
  const auto lim = cl.get_genome_size() - seed::key_weight_three + 1;
  for (std::size_t i = 0; i < lim; ++i) {
    shift_three_key<the_conv>(*gi++, hash_key);
    assert(nidx_itr != std::cend(exclude));
    if (i < nidx_itr->first && *keep_itr) {
      if (the_conv == c_to_t)
        counter_t[hash_key] += (!use_mask || !(*itl_itr));
      else
        counter_a[hash_key] += (!use_mask || !(*itl_itr));
    }
    if (nidx_itr->second <= i)
      ++nidx_itr;
    ++keep_itr;
    ++itl_itr;
  }
}

static inline std::uint32_t
two_letter_cost(const std::uint32_t count) {
  return count;
}

static inline std::uint64_t
three_letter_cost(const std::uint32_t count_t, const std::uint32_t count_a) {
  return (count_t + count_a) >> 1;
}

template <const bool use_mask>
void
AbismalIndex::initialize_bucket_sizes() {
  const auto s_time = abismal_clock::now();
  if (VERBOSE)
    std::clog << "[computing bucket sizes]";
  std::thread bucket_two_ltr([&] { get_bucket_sizes_two<use_mask>(); });
  std::thread bucket_ct([&] { get_bucket_sizes_three<c_to_t, use_mask>(); });
  std::thread bucket_ga([&] { get_bucket_sizes_three<g_to_a, use_mask>(); });
  bucket_two_ltr.join();
  bucket_ct.join();
  bucket_ga.join();
  if (VERBOSE)
    std::clog << delta_seconds(s_time) << '\n';
}

static inline auto
get_block_bounds(
  const std::size_t start_pos, const std::size_t step_size,
  const std::size_t end_pos,
  const std::vector<std::pair<std::size_t, std::size_t>> &exclude)
  -> std::vector<std::pair<std::size_t, std::size_t>> {
  std::vector<std::pair<std::size_t, std::size_t>> blocks;

  auto block_start = start_pos;
  auto i = std::cbegin(exclude);
  while (block_start < end_pos && i != std::cend(exclude)) {
    if (block_start < i->first) {
      const auto block_end =
        std::min({i->first, block_start + step_size, end_pos});
      blocks.emplace_back(block_start, block_end);
      block_start += step_size;

      // this might move block start backward
      if (block_start >= i->second)
        block_start = (*i++).second;
    }
    else
      block_start = (*i++).second;
  }
  while (block_start < end_pos) {
    const auto block_end = std::min({block_start + step_size, end_pos});
    blocks.emplace_back(block_start, block_end);
    block_start += step_size;
  }

  return blocks;
}

void
AbismalIndex::select_two_letter_positions() {
  constexpr std::size_t block_size = 1000000ul;
  const auto s_time = abismal_clock::now();
  if (VERBOSE)
    std::clog << "[selecting two-letter positions]";

  // choose which have lower count under two-letters
  is_two_let.resize(cl.get_genome_size(), false);
  const auto itl_beg = std::begin(is_two_let);
  const auto g_beg = genome_iterator(std::cbegin(genome));

  // ADS: this works below because seed::key_weight is always at least
  // seed::key_weight_three
  const auto lim = cl.get_genome_size() - seed::key_weight + 1;

  const auto blocks = get_block_bounds(0, block_size, lim, exclude);
  const auto n_blocks = std::size(blocks);

  const auto n_threads = std::min(n_blocks, n_threads_global);

  const auto blocks_beg = std::cbegin(blocks);
  const std::uint32_t n_per = (n_blocks + n_threads - 1) / n_threads;

  std::vector<std::thread> threads;
  for (auto i = 0ul; i < n_threads; ++i) {
    const auto curr_beg = blocks_beg + i * n_per;
    const auto curr_end = blocks_beg + std::min(n_blocks, (i + 1) * n_per);
    threads.emplace_back([&, curr_beg, curr_end] {
      // parallelism over blocks; each thread does a set
      for (auto block_itr = curr_beg; block_itr != curr_end; ++block_itr) {
        auto block = *block_itr;

        std::uint32_t hash_two = 0;
        auto gi_two = g_beg + block.first;

        // spool the two letter hash
        const auto gi_lim_two = gi_two + (seed::key_weight - 1);
        while (gi_two != gi_lim_two)
          shift_hash_key(*gi_two++, hash_two);

        std::uint32_t hash_t = 0, hash_a = 0;
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
    });
  }
  for (auto &thread : threads)
    thread.join();
  if (VERBOSE)
    std::clog << delta_seconds(s_time) << '\n';
}

void
AbismalIndex::hash_genome() {
  // count k-mers under each encoding with masking
  auto s_time = abismal_clock::now();
  if (VERBOSE)
    std::clog << "[initializing index structures]";
  const auto inc_scan = [](auto &v) {
    std::inclusive_scan(std::cbegin(v), std::cend(v), std::begin(v));
  };
  std::thread scanner([&] { inc_scan(counter); });
  std::thread scanner_t([&] { inc_scan(counter_t); });
  std::thread scanner_a([&] { inc_scan(counter_a); });
  scanner.join();
  scanner_t.join();
  scanner_a.join();

  index_size = counter[counter_size];
  index_size_three = counter_t[counter_size_three];

  index.resize(index_size, 0);
  index_t.resize(index_size_three, 0);
  index_a.resize(index_size_three, 0);
  if (VERBOSE)
    std::clog << delta_seconds(s_time) << '\n'
              << "[index sizes: two-letter=" << index_size << " "
              << "three-letter=" << index_size_three << "]\n";
  s_time = abismal_clock::now();

  if (VERBOSE)
    std::clog << "[counting k-mers]";
  const auto lim = cl.get_genome_size() - seed::key_weight + 1;

  std::thread ltr_counter([&] {
    auto gi_two = genome_iterator(std::cbegin(genome));
    const auto gi_lim = gi_two + (seed::key_weight - 1);
    std::uint32_t hash_two = 0;
    while (gi_two != gi_lim)
      shift_hash_key(*gi_two++, hash_two);
    auto nidx_itr = std::cbegin(exclude);
    for (std::size_t i = 0; i < lim; ++i) {
      shift_hash_key(*gi_two++, hash_two);
      assert(nidx_itr != std::cend(exclude));
      if (i < nidx_itr->first && keep[i] && is_two_let[i]) {
        assert(counter[hash_two] > 0);
        index[--counter[hash_two]] = i;
      }
      if (nidx_itr->second <= i)
        ++nidx_itr;
    }
  });

  std::thread ltr_counter_t([&] {
    auto gi_three = genome_iterator(std::cbegin(genome));
    const auto gi_lim = gi_three + (seed::key_weight_three - 1);
    std::uint32_t hash_t = 0;
    while (gi_three != gi_lim)
      shift_three_key<c_to_t>(*gi_three++, hash_t);
    auto nidx_itr = std::cbegin(exclude);
    for (std::size_t i = 0; i < lim; ++i) {
      shift_three_key<c_to_t>(*gi_three++, hash_t);
      assert(nidx_itr != std::cend(exclude));
      if (i < nidx_itr->first && keep[i] && !is_two_let[i]) {
        assert(counter_t[hash_t] > 0);
        index_t[--counter_t[hash_t]] = i;
      }
      if (nidx_itr->second <= i)
        ++nidx_itr;
    }
  });

  std::thread ltr_counter_a([&] {
    auto gi_three = genome_iterator(std::cbegin(genome));
    const auto gi_lim = gi_three + (seed::key_weight_three - 1);
    std::uint32_t hash_a = 0;
    while (gi_three != gi_lim)
      shift_three_key<g_to_a>(*gi_three++, hash_a);
    auto nidx_itr = std::cbegin(exclude);
    for (std::size_t i = 0; i < lim; ++i) {
      shift_three_key<g_to_a>(*gi_three++, hash_a);
      assert(nidx_itr != std::cend(exclude));
      if (i < nidx_itr->first && keep[i] && !is_two_let[i]) {
        assert(counter_a[hash_a] > 0);
        index_a[--counter_a[hash_a]] = i;
      }
      if (nidx_itr->second <= i)
        ++nidx_itr;
    }
  });

  ltr_counter.join();
  ltr_counter_t.join();
  ltr_counter_a.join();

  if (VERBOSE)
    std::clog << delta_seconds(s_time) << '\n';
}

struct dp_sol {
  static const std::uint64_t sentinel = num_lim<std::uint64_t>::max();
  std::uint64_t cost{sentinel};
  std::uint64_t prev{sentinel};
};

// ADS: Occupancy of the deque can be 0 through window_size, so
// (seed::window_size + 1). However, if we only use this size, then we can't
// tell the difference using `f` and `b` of whether the deque is empty
// or full. We need to ensure that the size of this deque ring buffer
// is at least (window_size + 2).
template <class T> struct fixed_ring_buffer {
  constexpr static std::size_t qsz =
    32ul;  // must be >= (seed::window_size + 2)
  constexpr static std::size_t qsz_msk =
    31ul;  // mask to keep positions in range

  std::array<T, qsz> h;

  auto
  size() const -> std::size_t {
    return (b == f) ? 0 : (b > f) ? b - f : qsz - (f - b);
  }

  auto
  empty() const -> bool {
    return f == b;
  }

  auto
  back() -> T & {
    return h[(b - 1) & qsz_msk];
  }

  auto
  front() -> T & {
    return h[f];
  }

  auto
  pop_back() -> void {
    --b;
    b &= qsz_msk;
  }

  auto
  pop_front() -> void {
    ++f;
    f &= qsz_msk;
  }

  auto
  push_back(T &&x) -> void {
    h[b] = std::move(x);
    ++b;
    b &= qsz_msk;
  }

  std::uint32_t f{0};
  std::uint32_t b{0};
};

using deque = fixed_ring_buffer<dp_sol>;

static inline void
add_sol(deque &helper, const std::uint64_t prev, const std::uint64_t cost) {
  while (!helper.empty() && helper.back().cost > cost)
    helper.pop_back();

  helper.push_back({cost, prev});

  while (helper.front().prev + seed::window_size <= prev)
    helper.pop_front();

  assert(!helper.empty() && helper.size() <= seed::window_size);
}

static inline std::uint64_t
hybrid_cost(const bool is_two_let, const std::uint32_t count_two,
            const std::uint32_t count_t, const std::uint32_t count_a) {
  return is_two_let ? two_letter_cost(count_two)
                    : three_letter_cost(count_t, count_a);
}

void
AbismalIndex::compress_dp() {
  constexpr auto block_size = 1000000ul;
  const auto s_time = abismal_clock::now();
  if (VERBOSE)
    std::clog << "[dynamic programming to optimize seed selection]";

  // default is no position will be indexed
  std::fill(std::begin(keep), std::end(keep), false);

  const auto lim = cl.get_genome_size() - seed::key_weight + 1;

  const auto blocks = get_block_bounds(0, block_size, lim, exclude);
  const auto n_blocks = std::size(blocks);

  const auto n_threads = std::min(n_blocks, n_threads_global);
  std::vector<std::thread> threads;

  const auto blocks_beg = std::cbegin(blocks);
  const auto n_per = (n_blocks + n_threads - 1) / n_threads;

  for (auto i = 0u; i < n_threads; ++i) {
    const auto curr_beg = blocks_beg + i * n_per;
    const auto curr_end = blocks_beg + std::min(n_blocks, (i + 1) * n_per);
    threads.emplace_back([&, curr_beg, curr_end] {
      for (auto block_itr = curr_beg; block_itr != curr_end; ++block_itr) {
        auto block = *block_itr;
        const std::size_t block_start = block.first;
        const std::uint32_t current_block_size = block.second - block.first;

        if (current_block_size < seed::window_size)
          continue;

        std::uint32_t hash_two = 0;

        // spool the hash for two letter encoding
        auto gi_two = genome_iterator(std::cbegin(genome)) + block_start;
        const auto gi_lim_two =
          gi_two + std::min(current_block_size, seed::key_weight - 1);
        while (gi_two != gi_lim_two)
          shift_hash_key(*gi_two++, hash_two);

        std::uint32_t hash_t = 0, hash_a = 0;

        // spool the hashes for three letter encodings
        auto gi_three = genome_iterator(std::cbegin(genome)) + block_start;
        const auto gi_lim_three(gi_three + (seed::key_weight_three - 1));
        while (gi_three != gi_lim_three) {
          const auto gi3_val = *gi_three++;
          shift_three_key<c_to_t>(gi3_val, hash_t);
          shift_three_key<g_to_a>(gi3_val, hash_a);
        }

        std::vector<dp_sol> opt(block_size + 1);
        deque helper;

        auto opt_itr = std::begin(opt);
        const auto opt_lim = std::cbegin(opt) + current_block_size;
        auto itl = std::cbegin(is_two_let) + block_start;

        // do the base cases for the dynamic programming: the solutions
        // for the first "window size" positions
        std::size_t i = 0;
        for (; i < seed::window_size; ++opt_itr) {
          shift_hash_key(*gi_two++, hash_two);
          const auto gi3_val = *gi_three++;
          shift_three_key<c_to_t>(gi3_val, hash_t);
          shift_three_key<g_to_a>(gi3_val, hash_a);

          assert(hash_two < counter.size());
          assert(hash_t < counter_t.size());
          assert(hash_a < counter_a.size());

          const auto c = hybrid_cost(*itl++, counter[hash_two],
                                     counter_t[hash_t], counter_a[hash_a]);
          assert(opt_itr != std::end(opt));
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

          const auto c = hybrid_cost(*itl++, counter[hash_two],
                                     counter_t[hash_t], counter_a[hash_a]);
          assert(helper.size() > 0);
          assert(opt_itr != std::end(opt));
          *opt_itr = helper.front();
          opt_itr->cost += c;
          add_sol(helper, i++, opt_itr->cost);
        }

        // get solution for start of traceback
        std::uint64_t opt_ans = num_lim<std::uint64_t>::max();
        std::uint64_t last = num_lim<std::uint64_t>::max();
        for (std::size_t i = current_block_size - 1;
             i >= current_block_size - seed::window_size; --i) {
          const auto cand_cost = opt[i].cost;
          if (cand_cost < opt_ans) {
            opt_ans = cand_cost;
            last = i;
          }
        }

        *opt_itr = {opt_ans, last};
        const auto block_keep = std::begin(keep) + block_start;
        const auto opt_beg = std::begin(opt);
        while (opt_itr->prev != dp_sol::sentinel) {
          *(block_keep + opt_itr->prev) = true;
          opt_itr = opt_beg + opt_itr->prev;
        }
      }
    });
  }

  for (auto &thread : threads)
    thread.join();

  max_candidates = 100u;  // GS: this is a heuristic
  if (VERBOSE)
    std::clog << delta_seconds(s_time) << '\n';
}

struct BucketLess {
  BucketLess(const Genome &g) : g_start(std::begin(g)) {}

  bool
  operator()(const std::uint32_t a, const std::uint32_t b) const {
    auto idx1(g_start + a + seed::key_weight);
    const auto lim1(g_start + a + seed::n_sorting_positions);
    auto idx2(g_start + b + seed::key_weight);
    while (idx1 != lim1) {
      const std::uint8_t c1 = get_bit(*idx1++);
      const std::uint8_t c2 = get_bit(*idx2++);
      if (c1 != c2)
        return c1 < c2;
    }
    return false;
  }

  const genome_iterator g_start;
};

template <const three_conv_type the_conv>
static inline three_letter_t
three_letter_num_srt(const std::uint8_t nt) {
  // C=T=0, A=1, G=4
  // A=G=0, C=2, T=8
  return the_conv == c_to_t ? nt & 5 : nt & 10;
}

template <const three_conv_type the_type> struct BucketLessThree {
  BucketLessThree(const Genome &g) : g_start(std::begin(g)) {}

  bool
  operator()(const std::uint32_t a, const std::uint32_t b) const {
    auto idx1(g_start + a + seed::key_weight_three);
    const auto lim1(g_start + a + seed::n_sorting_positions);
    auto idx2(g_start + b + seed::key_weight_three);
    while (idx1 != lim1) {
      const std::uint8_t c1 = three_letter_num_srt<the_type>(*idx1++);
      const std::uint8_t c2 = three_letter_num_srt<the_type>(*idx2++);
      if (c1 != c2)
        return c1 < c2;
    }
    return false;
  }

  const genome_iterator g_start;
};

void
AbismalIndex::sort_buckets() {
  {
    const auto s_time = abismal_clock::now();
    if (VERBOSE)
      std::clog << "[sorting two-letter buckets]";
    const auto b = std::begin(index);
    const BucketLess bucket_less(genome);
    const auto n_threads = std::min(counter_size, n_threads_global);
    std::vector<std::thread> threads;
    const auto n_per = (counter_size + n_threads - 1) / n_threads;
    for (auto i = 0u; i < n_threads; ++i) {
      const auto start_idx = i * n_per;
      const auto stop_idx = std::min((i + 1) * n_per, counter_size);
      threads.push_back(std::thread([&, start_idx, stop_idx]() {
        for (std::size_t i = start_idx; i < stop_idx; ++i)
          if (counter[i + 1] > counter[i] + 1)
            std::sort(b + counter[i], b + counter[i + 1], bucket_less);
      }));
    }
    for (auto &thread : threads)
      thread.join();
    if (VERBOSE)
      std::clog << delta_seconds(s_time) << '\n';
  }
  {
    const auto s_time = abismal_clock::now();
    if (VERBOSE)
      std::clog << "[sorting three-letter T buckets]";
    const auto b = std::begin(index_t);
    const BucketLessThree<c_to_t> bucket_less(genome);
    const auto n_threads = std::min(counter_size_three, n_threads_global);
    std::vector<std::thread> threads;
    const auto n_per = (counter_size_three + n_threads - 1) / n_threads;
    for (auto i = 0u; i < n_threads; ++i) {
      const auto start_idx = i * n_per;
      const auto stop_idx = std::min((i + 1) * n_per, counter_size_three);
      threads.push_back(std::thread([&, start_idx, stop_idx]() {
        for (std::size_t i = start_idx; i < stop_idx; ++i)
          if (counter_t[i + 1] > counter_t[i] + 1)
            std::sort(b + counter_t[i], b + counter_t[i + 1], bucket_less);
      }));
    }
    for (auto &thread : threads)
      thread.join();
    if (VERBOSE)
      std::clog << delta_seconds(s_time) << '\n';
  }
  {
    const auto s_time = abismal_clock::now();
    if (VERBOSE)
      std::clog << "[sorting three-letter A buckets]";
    const auto b = std::begin(index_a);
    const BucketLessThree<g_to_a> bucket_less(genome);
    const auto n_threads = std::min(counter_size_three, n_threads_global);
    std::vector<std::thread> threads;
    const auto n_per = (counter_size_three + n_threads - 1) / n_threads;
    for (auto i = 0u; i < n_threads; ++i) {
      const auto start_idx = i * n_per;
      const auto stop_idx = std::min((i + 1) * n_per, counter_size_three);
      threads.push_back(std::thread([&, start_idx, stop_idx]() {
        for (std::size_t i = start_idx; i < stop_idx; ++i)
          if (counter_a[i + 1] > counter_a[i] + 1)
            std::sort(b + counter_a[i], b + counter_a[i + 1], bucket_less);
      }));
    }
    for (auto &thread : threads)
      thread.join();
    if (VERBOSE)
      std::clog << delta_seconds(s_time) << '\n';
  }
}

static void
write_internal_identifier(FILE *out) {
  if (fwrite((char *)&AbismalIndex::internal_identifier[0], 1,
             AbismalIndex::internal_identifier.size(),
             out) != AbismalIndex::internal_identifier.size())
    throw std::runtime_error("failed writing index identifier");
}

void
seed::read(FILE *in) {
  static constexpr auto error_msg{"failed to read seed data"};

  // key_weight
  std::uint32_t key_weight_from_file = 0;
  if (fread((char *)&key_weight_from_file, sizeof(std::uint32_t), 1, in) != 1)
    throw std::runtime_error(error_msg);

  if (key_weight_from_file != key_weight) {
    throw std::runtime_error(
      "inconsistent k-mer size. Expected: " + std::to_string(key_weight) +
      ", got: " + std::to_string(key_weight_from_file));
  }

  // window_size
  std::uint32_t window_size_from_file = 0;
  if (fread((char *)&window_size_from_file, sizeof(std::uint32_t), 1, in) != 1)
    throw std::runtime_error(error_msg);

  if (window_size_from_file != window_size) {
    throw std::runtime_error("inconsistent window size size. Expected: " +
                             std::to_string(window_size) +
                             ", got: " + std::to_string(window_size_from_file));
  }

  // n_sorting_positions
  std::uint32_t n_sorting_positions_from_file = 0;
  if (fread((char *)&n_sorting_positions_from_file, sizeof(std::uint32_t), 1,
            in) != 1)
    throw std::runtime_error(error_msg);

  if (n_sorting_positions_from_file != n_sorting_positions) {
    throw std::runtime_error("inconsistent sorting size size. Expected: " +
                             std::to_string(n_sorting_positions) + ", got: " +
                             std::to_string(n_sorting_positions_from_file));
  }
}

void
seed::write(FILE *out) {
  static const std::string error_msg("failed to write seed data");
  if (fwrite((char *)&seed::key_weight, sizeof(std::uint32_t), 1, out) != 1 ||
      fwrite((char *)&seed::window_size, sizeof(std::uint32_t), 1, out) != 1 ||
      fwrite((char *)&seed::n_sorting_positions, sizeof(std::uint32_t), 1,
             out) != 1)
    throw std::runtime_error(error_msg);
}

void
AbismalIndex::write(const std::string &index_file) const {
  FILE *out = fopen(index_file.c_str(), "wb");
  if (!out)
    throw std::runtime_error("cannot open output file " + index_file);

  write_internal_identifier(out);
  seed::write(out);
  cl.write(out);

  if (fwrite((char *)&genome[0], sizeof(element_t), genome.size(), out) !=
        genome.size() ||
      fwrite((char *)&max_candidates, sizeof(std::uint32_t), 1, out) != 1 ||
      fwrite((char *)&counter_size, sizeof(std::size_t), 1, out) != 1 ||
      fwrite((char *)&counter_size_three, sizeof(std::size_t), 1, out) != 1 ||
      fwrite((char *)&index_size, sizeof(std::size_t), 1, out) != 1 ||
      fwrite((char *)&index_size_three, sizeof(std::size_t), 1, out) != 1 ||
      fwrite((char *)(&counter[0]), sizeof(std::uint32_t), counter_size + 1,
             out) != (counter_size + 1) ||
      fwrite((char *)(&counter_t[0]), sizeof(std::uint32_t),
             counter_size_three + 1, out) != (counter_size_three + 1) ||
      fwrite((char *)(&counter_a[0]), sizeof(std::uint32_t),
             counter_size_three + 1, out) != (counter_size_three + 1) ||

      fwrite((char *)(&index[0]), sizeof(std::uint32_t), index_size, out) !=
        (index_size) ||
      fwrite((char *)(&index_t[0]), sizeof(std::uint32_t), index_size_three,
             out) != (index_size_three) ||
      fwrite((char *)(&index_a[0]), sizeof(std::uint32_t), index_size_three,
             out) != (index_size_three))
    throw std::runtime_error("failed writing index");

  if (fclose(out) != 0)
    throw std::runtime_error("problem closing file: " + index_file);
}

static bool
check_internal_identifier(FILE *in) {
  std::string id_found;
  while (id_found.size() < AbismalIndex::internal_identifier.size())
    id_found.push_back(getc(in));

  return (id_found == AbismalIndex::internal_identifier);
}

void
AbismalIndex::read(const std::string &index_file) {
  static const std::string error_msg("failed loading index file");

  FILE *in = fopen(index_file.data(), "rb");
  if (!in)
    throw std::runtime_error("cannot open input file " + index_file);

  if (!check_internal_identifier(in))
    throw std::runtime_error("index file format problem: " + index_file);

  seed::read(in);
  cl.read(in);

  const std::size_t genome_to_read = (cl.get_genome_size() + 15) / 16;
  // read the 4-bit encoded genome
  genome.resize(genome_to_read);
  if (fread((char *)&genome[0], sizeof(element_t), genome_to_read, in) !=
      genome_to_read)
    throw std::runtime_error(error_msg);

  // read the sizes of counter and index vectors
  if (fread((char *)&max_candidates, sizeof(std::uint32_t), 1, in) != 1 ||
      fread((char *)&counter_size, sizeof(std::size_t), 1, in) != 1 ||
      fread((char *)&counter_size_three, sizeof(std::size_t), 1, in) != 1 ||
      fread((char *)&index_size, sizeof(std::size_t), 1, in) != 1 ||
      fread((char *)&index_size_three, sizeof(std::size_t), 1, in) != 1)
    throw std::runtime_error(error_msg);

  // allocate then read the counter vector
  counter = std::vector<std::uint32_t>(counter_size + 1);
  if (fread((char *)(&counter[0]), sizeof(std::uint32_t), (counter_size + 1),
            in) != (counter_size + 1))
    throw std::runtime_error(error_msg);

  counter_t = std::vector<std::uint32_t>(counter_size_three + 1);
  if (fread((char *)(&counter_t[0]), sizeof(std::uint32_t),
            (counter_size_three + 1), in) != (counter_size_three + 1))
    throw std::runtime_error(error_msg);

  counter_a = std::vector<std::uint32_t>(counter_size_three + 1);
  if (fread((char *)(&counter_a[0]), sizeof(std::uint32_t),
            (counter_size_three + 1), in) != (counter_size_three + 1))
    throw std::runtime_error(error_msg);

  // allocate the read the index vector
  index = std::vector<std::uint32_t>(index_size);
  if (fread((char *)(&index[0]), sizeof(std::uint32_t), index_size, in) !=
      index_size)
    throw std::runtime_error(error_msg);

  index_t = std::vector<std::uint32_t>(index_size_three);
  if (fread((char *)(&index_t[0]), sizeof(std::uint32_t), index_size_three,
            in) != index_size_three)
    throw std::runtime_error(error_msg);

  index_a = std::vector<std::uint32_t>(index_size_three);
  if (fread((char *)(&index_a[0]), sizeof(std::uint32_t), index_size_three,
            in) != index_size_three)
    throw std::runtime_error(error_msg);

  if (fclose(in) != 0)
    throw std::runtime_error("problem closing file: " + index_file);
}

std::ostream &
ChromLookup::write(std::ostream &out) const {
  const std::uint32_t n_chroms = names.size();
  out.write((char *)&n_chroms, sizeof(std::uint32_t));
  for (std::size_t i = 0; i < n_chroms; ++i) {
    const std::uint32_t name_size = names[i].length();
    out.write((char *)&name_size, sizeof(std::uint32_t));
    out.write(names[i].c_str(), name_size);
  }
  out.write((char *)(&starts[0]), sizeof(std::uint32_t) * (n_chroms + 1));
  return out;
}

FILE *
ChromLookup::write(FILE *out) const {
  const std::uint32_t n_chroms = names.size();
  fwrite((char *)&n_chroms, sizeof(std::uint32_t), 1, out);
  for (std::size_t i = 0; i < n_chroms; ++i) {
    const std::uint32_t name_size = names[i].length();
    fwrite((char *)&name_size, sizeof(std::uint32_t), 1, out);
    fwrite(names[i].c_str(), 1, name_size, out);
  }
  fwrite((char *)(&starts[0]), sizeof(std::uint32_t), n_chroms + 1, out);
  return out;
}

void
ChromLookup::write(const std::string &outfile) const {
  std::ofstream out(outfile, std::ios::binary);
  if (!out)
    throw std::runtime_error("cannot open output file " + outfile);
  write(out);
}

std::istream &
ChromLookup::read(std::istream &in) {
  // read the number of chroms
  std::uint32_t n_chroms = 0;
  in.read((char *)&n_chroms, sizeof(std::uint32_t));

  // allocate the number of chroms
  names.resize(n_chroms);

  // get each chrom name
  for (std::size_t i = 0; i < n_chroms; ++i) {
    std::uint32_t name_size = 0;
    // get the size of the chrom name
    in.read((char *)&name_size, sizeof(std::uint32_t));
    // allocate the chrom name
    names[i].resize(name_size);
    // read the chrom name
    in.read((char *)&names[i][0], name_size);
  }

  // allocate then read the starts vector
  starts = std::vector<std::uint32_t>(n_chroms + 1);
  in.read((char *)(&starts[0]), sizeof(std::uint32_t) * (n_chroms + 1));

  return in;
}

FILE *
ChromLookup::read(FILE *in) {
  constexpr auto error_msg = "failed loading chrom info from index";

  // read the number of chroms in the reference
  std::uint32_t n_chroms = 0;
  if (fread((char *)&n_chroms, sizeof(std::uint32_t), 1, in) != 1)
    throw std::runtime_error(error_msg);

  names.resize(n_chroms);  // allocate a name for each chrom

  // get each chrom name
  for (std::size_t i = 0; i < n_chroms; ++i) {
    std::uint32_t name_size = 0;
    // get size of chrom name
    if (fread((char *)&name_size, sizeof(std::uint32_t), 1, in) != 1)
      throw std::runtime_error(error_msg);
    // allocate the chrom name
    names[i].resize(name_size);
    // read the chrom name
    if (fread((char *)&names[i][0], 1, name_size, in) != name_size)
      throw std::runtime_error(error_msg);
  }

  // allocate then read the starts vector
  starts = std::vector<std::uint32_t>(n_chroms + 1);
  if (fread((char *)(&starts[0]), sizeof(std::uint32_t), n_chroms + 1, in) !=
      n_chroms + 1)
    throw std::runtime_error(error_msg);

  return in;
}

void
ChromLookup::read(const std::string &infile) {
  std::ifstream in(infile, std::ios::binary);
  if (!in)
    throw std::runtime_error("cannot open input file " + infile);
  read(in);
}

std::ostream &
operator<<(std::ostream &out, const ChromLookup &cl) {
  return out << cl.tostring();
}

std::string
ChromLookup::tostring() const {
  std::ostringstream iss;
  for (auto i = 0u; i < names.size(); ++i)
    iss << i << '\t' << names[i] << '\t' << starts[i] << '\t' << starts[i + 1]
        << '\n';
  return iss.str();
}

void
ChromLookup::get_chrom_idx_and_offset(const std::uint32_t pos,
                                      std::uint32_t &chrom_idx,
                                      std::uint32_t &offset) const {
  auto idx = upper_bound(std::cbegin(starts), std::cend(starts), pos);

  assert(idx != std::cbegin(starts));

  --idx;

  chrom_idx = std::distance(std::cbegin(starts), idx);
  offset = pos - starts[chrom_idx];
}

std::uint32_t
ChromLookup::get_pos(const std::string &chrom,
                     const std::uint32_t offset) const {
  const auto itr = find(std::cbegin(names), std::cend(names), chrom);
  return itr == std::cend(names)
           ? num_lim<std::uint32_t>::max()
           : starts[std::distance(std::cbegin(names), itr)] + offset;
}

bool
ChromLookup::get_chrom_idx_and_offset(const std::uint32_t pos,
                                      const std::uint32_t readlen,
                                      std::uint32_t &chrom_idx,
                                      std::uint32_t &offset) const {
  auto idx = upper_bound(std::cbegin(starts), std::cend(starts), pos);

  if (idx == std::cbegin(starts))
    return false;  // read is before any chrom

  --idx;

  chrom_idx = std::distance(std::cbegin(starts), idx);
  offset = pos - starts[chrom_idx];
  return (pos + readlen <= starts[chrom_idx + 1]);
}

template <class G>
void
load_genome_impl(const std::string &genome_file, G &genome, ChromLookup &cl) {
  const auto file_size = std::filesystem::file_size(genome_file);

  bamxx::bgzf_file in(genome_file, "r");
  if (!in)
    throw std::runtime_error("failed to open genome file: " + genome_file);

  genome.clear();
  // reserve space for padding at both ends
  genome.reserve(file_size + 2 * seed::padding_size);
  auto g_ins = std::back_inserter(genome);

  // pad the start of the concatenated sequence
  cl.names.push_back("pad_start");
  cl.starts.push_back(0);
  fill_n(g_ins, seed::padding_size, 'N');

  std::string line;
  while (getline(in, line))
    if (line[0] != '>')
      std::copy(std::cbegin(line), std::cend(line), g_ins);
    else {
      cl.names.push_back(line.substr(1, line.find_first_of(" \t") - 1));
      cl.starts.push_back(genome.size());
    }

  if (cl.names.size() < 2)
    throw std::runtime_error("no names found in genome file");

  // now pad the end of the concatenated sequence
  cl.names.push_back("pad_end");
  cl.starts.push_back(genome.size());
  fill_n(g_ins, seed::padding_size, 'N');

  // this one additional "start" is the end of all chroms
  cl.starts.push_back(genome.size());
}

void
load_genome(const std::string &genome_file, std::string &genome,
            ChromLookup &cl) {
  load_genome_impl(genome_file, genome, cl);
}

void
load_genome(const std::string &genome_file, std::vector<std::uint8_t> &genome,
            ChromLookup &cl) {
  load_genome_impl(genome_file, genome, cl);
}

// NOLINTEND
