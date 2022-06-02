/* Copyright (C) 2018 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
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

#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"
#include "dna_four_bit.hpp"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <iterator>
#include <algorithm>
#include <deque>
#include <fstream>
#include <utility>

using std::pair;
using std::make_pair;
using std::ofstream;
using std::ifstream;
using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::runtime_error;
using std::sort;
using std::cout;
using std::min;
using std::deque;
using std::fill;
using std::to_string;
using std::max;

#define _DO_DP
#define _TWOLETTER
#define _THREELETTER

#if defined(_TWOLETTER) || defined(_THREELETTER)
  #define _HYBRID
#endif

bool AbismalIndex::VERBOSE = false;
string AbismalIndex::internal_identifier = "AbismalIndex";

using genome_iterator = genome_four_bit_itr;

void
AbismalIndex::create_index(const string &genome_file) {
  vector<uint8_t> inflated_genome;
  if (VERBOSE)
    cerr << "[loading genome]" << endl;
  load_genome(genome_file, inflated_genome, cl);
  encode_genome(inflated_genome);
  vector<uint8_t>().swap(inflated_genome);

  // creat genome-wide mask of positions to keep
  keep.resize(cl.get_genome_size());
  fill(begin(keep), end(keep), true);

  select_two_letter_positions();
  compress_dp();
  hash_genome();
  sort_buckets();
}


void
AbismalIndex::encode_genome(const vector<uint8_t> &input_genome) {
  if (VERBOSE)
    cerr << "[encoding genome]" << endl;

  genome.resize(input_genome.size()/32 + (input_genome.size()%32 != 0));
  encode_dna_four_bit(begin(input_genome), end(input_genome), begin(genome));
}

template<const bool use_mask>
void
AbismalIndex::get_bucket_sizes() {
  counter.clear();

  counter_size = (1ull << seed::key_weight);

  // the "counter" has an additional entry for convenience
  counter.resize(counter_size + 1, 0);

  const size_t genome_st = seed::padding_size;
  const size_t lim = cl.get_genome_size() - seed::key_weight - seed::padding_size;
  ProgressBar progress(lim, "counting " + to_string(seed::key_weight) +
                            "-bit words");

  if (VERBOSE)
    progress.report(cerr, 0);

  // start building up the hash key
  genome_iterator gi(begin(genome));
  gi = gi + genome_st;

  const auto gi_lim(gi + (seed::key_weight - 1));
  uint32_t hash_key = 0;
  while (gi != gi_lim)
    shift_hash_key(*gi++, hash_key);

  for (size_t i = genome_st; i < lim; ++i) {
    shift_hash_key(*gi++, hash_key);

    if (keep[i]) {
      if (VERBOSE && progress.time_to_report(i))
        progress.report(cerr, i);

      const bool count_base =  (!use_mask || is_two_letter[i]);
      counter[hash_key] += count_base;
    }
  }
  if (VERBOSE)
    progress.report(cerr, lim);
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

  const size_t genome_st = seed::padding_size;
  const size_t lim = cl.get_genome_size() - seed::key_weight_three - seed::padding_size;
  ProgressBar progress(lim, "counting " + to_string(seed::key_weight_three) +
                            "-3 letters");

  if (VERBOSE)
    progress.report(cerr, 0);

  // start building up the hash key
  genome_iterator gi(begin(genome));
  gi = gi + genome_st;

  const auto gi_lim(gi + (seed::key_weight_three - 1));
  uint32_t hash_key = 0;
  while (gi != gi_lim)
    shift_three_key<the_conv>(*gi++, hash_key);

  for (size_t i = genome_st; i < lim; ++i) {
    shift_three_key<the_conv>(*gi++, hash_key);
    if (keep[i]) {
      if (VERBOSE && progress.time_to_report(i))
        progress.report(cerr, i);

      const bool count_base = (!use_mask || !is_two_letter[i]);

      if (the_conv == c_to_t)
        counter_t[hash_key] += count_base;
      else
        counter_a[hash_key] += count_base;
    }
  }
  if (VERBOSE)
    progress.report(cerr, lim);
}

inline uint32_t
two_letter_cost(const uint32_t count) {
  return (count);
}

inline uint32_t
three_letter_cost(const uint32_t count_t, const uint32_t count_a) {
  return ((count_t + count_a) >> 1);
}

void
AbismalIndex::select_two_letter_positions() {
  // first get statistics on the full genome
#pragma omp parallel for
  for (size_t i = 0; i < 3; ++i) {
#if defined(_TWOLETTER)
    if (i == 0)
      get_bucket_sizes<false>();
#endif
#if defined(_THREELETTER)
    if (i == 1)
      get_bucket_sizes_three<c_to_t, false>();
    if (i == 2)
      get_bucket_sizes_three<g_to_a, false >();
#endif
  }

  // now choose which have lower count under two-letters
  is_two_letter.resize(cl.get_genome_size());
  fill(begin(is_two_letter), end(is_two_letter), false);
#if defined(_THREELETTER) && !defined(_TWOLETTER)
  return;
#endif

#if defined(_TWOLETTER) && !defined(_THREELETTER)
  fill(begin(is_two_letter), end(is_two_letter), true);
  return;
#endif

#if defined(_HYBRID)
  const size_t genome_st = seed::padding_size;
  const size_t lim = cl.get_genome_size() - seed::key_weight - seed::padding_size;
  ProgressBar progress(lim, "building hybrid index");

  if (VERBOSE)
    progress.report(cerr, 0);

  // start building up the hash key
  genome_iterator gi_two(begin(genome));
  genome_iterator gi_three(begin(genome));
  gi_two = gi_two + genome_st;
  gi_three = gi_three + genome_st;

  uint32_t hash_two = 0;
  uint32_t hash_t = 0;
  uint32_t hash_a = 0;

  const auto gi_lim_two(gi_two + (seed::key_weight - 1));
  const auto gi_lim_three(gi_three + (seed::key_weight_three - 1));
  while (gi_two != gi_lim_two) {
    shift_hash_key(*gi_two++, hash_two);
  }

  while (gi_three != gi_lim_three) {
    shift_three_key<c_to_t>(*gi_three, hash_t);
    shift_three_key<g_to_a>(*gi_three, hash_a);
    ++gi_three;
  }

  size_t tot_two = 0;
  size_t tot_three = 0;
  for (size_t i = genome_st; i < lim; ++i, ++gi_two, ++gi_three) {
    shift_hash_key(*gi_two, hash_two);
    shift_three_key<c_to_t>(*gi_three, hash_t);
    shift_three_key<g_to_a>(*gi_three, hash_a);

    if (keep[i]) {
      if (VERBOSE && progress.time_to_report(i))
        progress.report(cerr, i);

      is_two_letter[i] = (
        two_letter_cost(counter[hash_two]) <=
        three_letter_cost(counter_t[hash_t], counter_a[hash_a])
      );

      tot_two += is_two_letter[i];
      tot_three += !is_two_letter[i];
    }
  }
  cerr << "\ntotal two letter: " << tot_two << "\n";
  cerr << "total three letter: " << tot_three << "\n";
#endif
}

void
AbismalIndex::hash_genome() {

  // count k-mers under each encoding with masking
#pragma omp parallel for
  for (size_t i = 0; i < 3; ++i) {
    if (i == 0)
    get_bucket_sizes<true>();
    else if (i == 1)
      get_bucket_sizes_three<c_to_t, true>();
    else if (i == 2)
      get_bucket_sizes_three<g_to_a,  true>();
  }

  /*
  //============== DEBUG ===============
  ofstream of_two("counts-two.txt");
  ofstream of_three_t("counts-three-t.txt");
  ofstream of_three_a("counts-three-a.txt");
  cerr << "[saving 2-letter counts to counts-two.txt]\n";
  for (size_t i = 0; i < counter.size(); ++i)
    of_two << i << " " << counter[i] << "\n";

  cerr << "[saving 3-letter counts to counts-three-t.txt and counts-three-a.txt]\n";
  for (size_t i = 0; i < counter_t.size(); ++i) {
    of_three_t << i << " " << counter_t[i] << "\n";
    of_three_a << i << " " << counter_a[i] << "\n";
  }
  //============== END DEBUG ===============*/

  if (VERBOSE)
    cerr << "[allocating hash tables]" << endl;

  std::partial_sum(begin(counter), end(counter), begin(counter));
  std::partial_sum(begin(counter_t), end(counter_t), begin(counter_t));
  std::partial_sum(begin(counter_a), end(counter_a), begin(counter_a));

  index_size = counter[counter_size];
  index_size_three = counter_t[counter_size_three];

  index.resize(index_size, 0);
  index_t.resize(index_size_three, 0);
  index_a.resize(index_size_three, 0);
  cerr << "[index sizes: " << index_size << " " << index_size_three << "]\n";

  const size_t genome_st = seed::padding_size;
  const size_t lim = cl.get_genome_size() - seed::key_weight - seed::padding_size;
  ProgressBar progress(lim, "hashing genome");

  // start building up the hash key
  genome_iterator gi_two(begin(genome));
  genome_iterator gi_three(begin(genome));
  gi_two = gi_two + genome_st;
  gi_three = gi_three + genome_st;

  const auto gi_lim_two(gi_two + (seed::key_weight - 1));
  const auto gi_lim_three(gi_three + (seed::key_weight_three - 1));

  uint32_t hash_two = 0;
  uint32_t hash_t = 0;
  uint32_t hash_a = 0;

  while (gi_two != gi_lim_two) {
    shift_hash_key(*gi_two++, hash_two);
  }

  while (gi_three != gi_lim_three) {
    shift_three_key<c_to_t>(*gi_three, hash_t);
    shift_three_key<g_to_a>(*gi_three, hash_a);
    ++gi_three;
  }

  for (size_t i = genome_st; i < lim; ++i, ++gi_two, ++gi_three) {
    shift_hash_key(*gi_two, hash_two);
    shift_three_key<c_to_t>(*gi_three, hash_t);
    shift_three_key<g_to_a>(*gi_three, hash_a);

    if (keep[i]) {
      if (VERBOSE && progress.time_to_report(i))
        progress.report(cerr, i);

      if (is_two_letter[i])
        index[--counter[hash_two]] = i;
      else {
        index_t[--counter_t[hash_t]] = i;
        index_a[--counter_a[hash_a]] = i;
      }
    }
  }
  if (VERBOSE)
    progress.report(cerr, lim);
}

struct dp_sol {
  size_t cost;
  uint32_t prev;
  uint32_t num;
  dp_sol() { reset(); }
  void reset() { cost = INF; prev = NIL; num = 0; }
  dp_sol(const uint32_t c) : cost(c), prev(NIL), num(0) {}
  dp_sol(const uint32_t c, const uint32_t p) : cost(c), prev(p), num(0) {}

  static const size_t INF = std::numeric_limits<size_t>::max();
  static const uint32_t NIL = std::numeric_limits<uint32_t>::max();
};

typedef pair<dp_sol, uint32_t> helper_sol;

void
add_sol(deque<helper_sol> &helper, const uint32_t pos, const dp_sol &cand) {
#if defined(_DO_DP)
  while (!helper.empty() && helper.back().first.cost > cand.cost)
    helper.pop_back();

  helper.push_back(make_pair(cand, pos));

  while (helper.front().second + seed::window_size <= pos)
    helper.pop_front();
  //assert(helper.size() <= seed::window_size);
  //assert(!helper.empty());
  //assert(is_helper_sorted(helper));

#else
  if (helper.size() == seed::window_size)
    helper.pop_front();
  helper.push_back(make_pair(cand, pos));
#endif

}

inline uint32_t
get_hybrid_cost(const bool is_two_letter,
                const uint32_t count_two,
                const uint32_t count_t, const uint32_t count_a) {
#if defined(_TWOLETTER) && !defined(_THREELETTER)
  return two_letter_cost(count_two);
#endif

#if defined(_THREELETTER) && !defined(_TWOLETTER)
  return three_letter_cost(count_t, count_a);
#endif

#if defined(_HYBRID)
  return (is_two_letter ? two_letter_cost(count_two) :
                          three_letter_cost(count_t, count_a));
#endif
}

void
AbismalIndex::compress_dp() {
  // no position is indexed
  fill(begin(keep), end(keep), false);

  const size_t lim =
    cl.get_genome_size() - seed::padding_size - seed::key_weight;

  // the dp memory allocation
  static const size_t BLOCK_SIZE = 10000000;
  deque<helper_sol> helper;
  vector<dp_sol> opt(BLOCK_SIZE);


  ProgressBar progress(lim, "solving dynamic prog.");

  genome_iterator gi_two(begin(genome));
  genome_iterator gi_three(begin(genome));
  uint32_t hash_two = 0;
  uint32_t hash_t = 0;
  uint32_t hash_a = 0;
  size_t hit_rate = 0;
  size_t num_bases = 0;

  // fast forward padding positions
  gi_two = gi_two + seed::padding_size;
  gi_three = gi_three + seed::padding_size;

  // build the first hash key minus last base
  const auto gi_lim_two(gi_three + (seed::key_weight - 1));
  while (gi_two != gi_lim_two)
    shift_hash_key(*gi_two++, hash_two);

  const auto gi_lim_three(gi_three + (seed::key_weight_three - 1));
  while (gi_three != gi_lim_three) {
    shift_three_key<c_to_t>(*gi_three, hash_t);
    shift_three_key<g_to_a>(*gi_three++, hash_a);
  }

  size_t beg = seed::padding_size;
  while (beg < lim) {
    // the position up to which we will solve
    const size_t fin = min(beg + BLOCK_SIZE, lim);

    // the problem size
    const size_t sz = fin - beg;

    // resets without needing to reallocate
    for (size_t i = 0; i < sz; ++i)
      opt[i].reset();

    // get the first w solutions
    helper.clear();
    for (size_t i = 0; i < seed::window_size; ++i) {
      shift_hash_key(*gi_two++, hash_two);
      shift_three_key<c_to_t>(*gi_three, hash_t);
      shift_three_key<g_to_a>(*gi_three++, hash_a);

      opt[i].cost =
        get_hybrid_cost(is_two_letter[beg + i],
            counter[hash_two], counter_t[hash_t], counter_a[hash_a]);
      opt[i].prev = dp_sol::NIL;
      opt[i].num = 1;
      add_sol(helper, i, opt[i]);
    }

    // main dp loop
    for (size_t i = seed::window_size; i < sz; ++i) {
      if (VERBOSE && progress.time_to_report(beg + i))
        progress.report(cerr, beg + i);

      shift_hash_key(*gi_two++, hash_two);
      shift_three_key<c_to_t>(*gi_three, hash_t);
      shift_three_key<g_to_a>(*gi_three++, hash_a);

      opt[i].cost = helper.front().first.cost +
        get_hybrid_cost(is_two_letter[beg + i],
                         counter[hash_two], counter_t[hash_t], counter_a[hash_a]);

      opt[i].prev = helper.front().second;
      opt[i].num = helper.front().first.num + 1;
      add_sol(helper, i, opt[i]);
    }

    size_t opt_ans = std::numeric_limits<size_t>::max();
    size_t opt_num = std::numeric_limits<size_t>::max();
    uint32_t last = dp_sol::NIL;

    // get the final solution
    for (size_t i = sz - 1; i >= sz - seed::window_size; --i) {
      const size_t cand_cost = opt[i].cost;
      const size_t cand_num = opt[i].num;
      if (cand_cost < opt_ans || (cand_cost == opt_ans && cand_num < opt_num)) {
        opt_ans = cand_cost;
        opt_num = cand_num;
        last = i;
      }
    }

    hit_rate += opt_ans;

    // build traceback
    while (last != dp_sol::NIL) {
      keep[beg + last] = true;
      ++num_bases;
      last = opt[last].prev;
    }

    // go to next block
    beg = fin;
  }

  // GS: this is a heuristic
  max_candidates = 100u;
  if (VERBOSE)
    progress.report(cerr, lim);

}

struct BucketLess {
  BucketLess(const Genome &g) : g_start(begin(g)) {}
  bool operator()(const uint32_t a, const uint32_t b) const {
    auto idx1(g_start + a + seed::key_weight);
    const auto lim1(g_start + a + seed::n_sorting_positions);
    auto idx2(g_start + b + seed::key_weight);
    while (idx1 != lim1) {
      const uint8_t c1 = get_bit(*(idx1++));
      const uint8_t c2 = get_bit(*(idx2++));
      if (c1 != c2) return c1 < c2;
    }
    return false;
  }
  const genome_iterator g_start;
};

template<const three_conv_type the_type>
struct BucketLessThree {
  BucketLessThree(const Genome &g) : g_start(begin(g)) {}
  bool operator()(const uint32_t a, const uint32_t b) const {
    auto idx1(g_start + a + seed::key_weight_three);
    const auto lim1(g_start + a + seed::n_sorting_positions);
    auto idx2(g_start + b + seed::key_weight_three);
    while (idx1 != lim1) {
      const uint8_t c1 = get_three_letter_num<the_type>(*(idx1++));
      const uint8_t c2 = get_three_letter_num<the_type>(*(idx2++));
      if (c1 != c2) return c1 < c2;
    }
    return false;
  }
  const genome_iterator g_start;
};


void
AbismalIndex::sort_buckets() {
  if (VERBOSE)
    cerr << "[sorting two-letter buckets]" << endl;
  const vector<uint32_t>::iterator b(begin(index));
  const BucketLess bucket_less(genome);
#pragma omp parallel for
  for (size_t i = 0; i < counter_size; ++i)
    if (counter[i + 1] > counter[i] + 1) {
      sort(b + counter[i], b + counter[i + 1], bucket_less);
    }

  cerr << "[sorting three-letter buckets C-to-T]" << endl;
  const vector<uint32_t>::iterator b_t(begin(index_t));
  const BucketLessThree<c_to_t> bucket_less_t(genome);
#pragma omp parallel for
  for (size_t i = 0; i < counter_size_three; ++i)
    if (counter_t[i + 1] > counter_t[i] + 1) {
      sort(b_t + counter_t[i], b_t + counter_t[i + 1], bucket_less_t);
    }

  cerr << "[sorting three-letter buckets G-to-A]" << endl;
  const vector<uint32_t>::iterator b_a(begin(index_a));
  const BucketLessThree<g_to_a> bucket_less_a(genome);
#pragma omp parallel for
  for (size_t i = 0; i < counter_size_three; ++i)
    if (counter_a[i + 1] > counter_a[i] + 1) {
      sort(b_a + counter_a[i], b_a + counter_a[i + 1], bucket_less_a);
    }
}
static void
write_internal_identifier(FILE *out) {
  if (fwrite((char*)&AbismalIndex::internal_identifier[0], 1,
             AbismalIndex::internal_identifier.size(), out) !=
      AbismalIndex::internal_identifier.size())
    throw runtime_error("failed writing index identifier");
}

void seed::read(FILE* in) {
  static const std::string error_msg("failed to read seed data");
  uint32_t _key_weight;
  uint32_t _window_size;
  uint32_t _n_sorting_positions;

  // key_weight
  if (fread((char*)&_key_weight, sizeof(uint32_t), 1, in) != 1)
    throw runtime_error(error_msg);

  if(_key_weight != key_weight) {
    throw runtime_error("inconsistent k-mer size. Expected: " +
        to_string(key_weight) + ", got: " + to_string(_key_weight));
  }

  // window_size
  if (fread((char*)&_window_size, sizeof(uint32_t), 1, in) != 1)
    throw runtime_error(error_msg);

  if (_window_size != window_size) {
    throw runtime_error("inconsistent window size size. Expected: " +
        to_string(window_size) + ", got: " + to_string(_window_size));
  }

  // n_sorting_positions
  if (fread((char*)&_n_sorting_positions, sizeof(uint32_t), 1, in) != 1)
    throw runtime_error(error_msg);

  if (_n_sorting_positions != n_sorting_positions) {
    throw runtime_error("inconsistent sorting size size. Expected: " +
        to_string(n_sorting_positions) + ", got: " +
        to_string(_n_sorting_positions));
  }
}

void seed::write(FILE*out) {
  static const std::string error_msg("failed to write seed data");
  if (fwrite((char*)&seed::key_weight, sizeof(uint32_t), 1, out) != 1 ||
      fwrite((char*)&seed::window_size, sizeof(uint32_t), 1, out) != 1 ||
      fwrite((char*)&seed::n_sorting_positions, sizeof(uint32_t), 1, out) != 1)
    throw runtime_error(error_msg);
}

void
AbismalIndex::write(const string &index_file) const {

  FILE *out = fopen(index_file.c_str(), "wb");
  if (!out)
    throw runtime_error("cannot open output file " + index_file);

  write_internal_identifier(out);
  seed::write(out);
  cl.write(out);

  if (fwrite((char*)&genome[0], sizeof(element_t), genome.size(), out) != genome.size() ||
      fwrite((char*)&max_candidates, sizeof(uint32_t), 1, out) != 1 ||
      fwrite((char*)&counter_size, sizeof(size_t), 1, out) != 1 ||
      fwrite((char*)&counter_size_three, sizeof(size_t), 1, out) != 1 ||
      fwrite((char*)&index_size, sizeof(size_t), 1, out) != 1 ||
      fwrite((char*)&index_size_three, sizeof(size_t), 1, out) != 1 ||
      fwrite((char*)(&counter[0]), sizeof(uint32_t),
             counter_size + 1, out) != (counter_size + 1) ||
      fwrite((char*)(&counter_t[0]), sizeof(uint32_t),
             counter_size_three + 1, out) != (counter_size_three + 1) ||
      fwrite((char*)(&counter_a[0]), sizeof(uint32_t),
             counter_size_three + 1, out) != (counter_size_three + 1) ||

      fwrite((char*)(&index[0]), sizeof(uint32_t),
             index_size, out) != (index_size) ||
      fwrite((char*)(&index_t[0]), sizeof(uint32_t),
             index_size_three, out) != (index_size_three) ||
      fwrite((char*)(&index_a[0]), sizeof(uint32_t),
             index_size_three, out) != (index_size_three))
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
  if (!in)
    throw runtime_error("cannot open input file " + index_file);

  if (!check_internal_identifier(in))
    throw runtime_error("index file format problem: " + index_file);

  seed::read(in);
  cl.read(in);

  const size_t genome_to_read = (cl.get_genome_size() + 31)/32;
  // read the 4-bit encoded genome
  genome.resize(genome_to_read);
  if (fread((char*)&genome[0], sizeof(element_t), genome_to_read, in) != genome_to_read)
    throw runtime_error(error_msg);

  // read the sizes of counter and index vectors
  if (fread((char*)&max_candidates, sizeof(uint32_t), 1, in) != 1 ||
      fread((char*)&counter_size, sizeof(size_t), 1, in) != 1 ||
      fread((char*)&counter_size_three, sizeof(size_t), 1, in) != 1 ||
      fread((char*)&index_size, sizeof(size_t), 1, in) != 1 ||
      fread((char*)&index_size_three, sizeof(size_t), 1, in) != 1)
    throw runtime_error(error_msg);

  // allocate then read the counter vector
  counter = vector<uint32_t>(counter_size + 1);
  if (fread((char*)(&counter[0]), sizeof(uint32_t),
            (counter_size + 1), in) != (counter_size + 1))
    throw runtime_error(error_msg);

  counter_t = vector<uint32_t>(counter_size_three + 1);
  if (fread((char*)(&counter_t[0]), sizeof(uint32_t),
            (counter_size_three + 1), in) != (counter_size_three + 1))
    throw runtime_error(error_msg);

  counter_a = vector<uint32_t>(counter_size_three + 1);
  if (fread((char*)(&counter_a[0]), sizeof(uint32_t),
            (counter_size_three + 1), in) != (counter_size_three + 1))
    throw runtime_error(error_msg);



  // allocate the read the index vector
  index = vector<uint32_t>(index_size);
  if (fread((char*)(&index[0]), sizeof(uint32_t), index_size, in) != index_size)
    throw runtime_error(error_msg);

  index_t = vector<uint32_t>(index_size_three);
  if (fread((char*)(&index_t[0]), sizeof(uint32_t), index_size_three, in)
      != index_size_three)
    throw runtime_error(error_msg);

  index_a = vector<uint32_t>(index_size_three);
  if (fread((char*)(&index_a[0]), sizeof(uint32_t), index_size_three, in)
      != index_size_three)
    throw runtime_error(error_msg);

  if (fclose(in) != 0)
    throw runtime_error("problem closing file: " + index_file);
}

std::ostream &
ChromLookup::write(std::ostream &out) const {
  const uint32_t n_chroms = names.size();
  out.write((char*)&n_chroms, sizeof(uint32_t));
  for (size_t i = 0; i < n_chroms; ++i) {
    const uint32_t name_size = names[i].length();
    out.write((char*)&name_size, sizeof(uint32_t));
    out.write(names[i].c_str(), name_size);
  }
  out.write((char*)(&starts[0]), sizeof(uint32_t)*(n_chroms + 1));
  return out;
}

FILE *
ChromLookup::write(FILE *out) const {
  const uint32_t n_chroms = names.size();
  fwrite((char*)&n_chroms, sizeof(uint32_t), 1, out);
  for (size_t i = 0; i < n_chroms; ++i) {
    const uint32_t name_size = names[i].length();
    fwrite((char*)&name_size, sizeof(uint32_t), 1, out);
    fwrite(names[i].c_str(), 1, name_size, out);
  }
  fwrite((char*)(&starts[0]), sizeof(uint32_t), n_chroms + 1, out);
  return out;
}

void
ChromLookup::write(const string &outfile) const {
  std::ofstream out(outfile, std::ios::binary);
  if (!out)
    throw runtime_error("cannot open output file " + outfile);
  write(out);
}

std::istream &
ChromLookup::read(std::istream &in) {
  /* read the number of chroms */
  uint32_t n_chroms = 0;
  in.read((char*)&n_chroms, sizeof(uint32_t));

  /* allocate the number of chroms */
  names.resize(n_chroms);

  /* get each chrom name */
  for (size_t i = 0; i < n_chroms; ++i) {
    uint32_t name_size = 0;
    /* get the size of the chrom name */
    in.read((char*)&name_size, sizeof(uint32_t));
    /* allocate the chrom name */
    names[i].resize(name_size);
    /* read the chrom name */
    in.read((char*)&names[i][0], name_size);
  }

  /* allocate then read the starts vector */
  starts = vector<uint32_t>(n_chroms + 1);
  in.read((char*)(&starts[0]), sizeof(uint32_t)*(n_chroms + 1));

  return in;
}

FILE *
ChromLookup::read(FILE *in) {

  static const string error_msg("failed loading chrom info from index");

  /* read the number of chroms */
  uint32_t n_chroms = 0;
  if (fread((char*)&n_chroms, sizeof(uint32_t), 1, in) != 1)
    throw runtime_error(error_msg);

  /* allocate the number of chroms */
  names.resize(n_chroms);

  /* get each chrom name */
  for (size_t i = 0; i < n_chroms; ++i) {
    uint32_t name_size = 0;
    /* get size of chrom name */
    if (fread((char*)&name_size, sizeof(uint32_t), 1, in) != 1)
      throw runtime_error(error_msg);
    /* allocate the chrom name */
    names[i].resize(name_size);
    /* read the chrom name */
    if (fread((char*)&names[i][0], 1, name_size, in) != name_size)
      throw runtime_error(error_msg);
  }

  /* allocate then read the starts vector */
  starts = vector<uint32_t>(n_chroms + 1);
  if (fread((char*)(&starts[0]), sizeof(uint32_t), n_chroms + 1, in)
      != n_chroms + 1)
    throw runtime_error(error_msg);

  return in;
}

void
ChromLookup::read(const std::string &infile) {
  std::ifstream in(infile, std::ios::binary);
  if (!in)
    throw runtime_error("cannot open input file " + infile);
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
ChromLookup::get_chrom_idx_and_offset(const uint32_t pos,
                                      uint32_t &chrom_idx,
                                      uint32_t &offset) const {

  vector<uint32_t>::const_iterator idx =
    upper_bound(begin(starts), end(starts), pos);

  //assert(idx != begin(starts));

  --idx;

  chrom_idx = idx - begin(starts);
  offset = pos - starts[chrom_idx];
}

uint32_t
ChromLookup::get_pos(const string &chrom, const uint32_t offset) const {
  vector<string>::const_iterator itr(find(begin(names), end(names), chrom));
  return (itr == end(names)) ?
    std::numeric_limits<uint32_t>::max() : starts[itr - begin(names)] + offset;
}

bool
ChromLookup::get_chrom_idx_and_offset(const uint32_t pos,
                                      const uint32_t readlen,
                                      uint32_t &chrom_idx,
                                      uint32_t &offset) const {

  vector<uint32_t>::const_iterator idx =
    upper_bound(begin(starts), end(starts), pos);

  // read fell behind the beginning of the chromosome
  if (idx == begin(starts))
    return false;

  --idx;

  chrom_idx = idx - begin(starts);
  offset = pos - starts[chrom_idx];
  return (pos + readlen <= starts[chrom_idx + 1]);
}
