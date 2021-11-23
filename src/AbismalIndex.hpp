/*  Copyright (C) 2018-2019 Andrew D. Smith
 *
 *  Authors: Andrew D. Smith
 *
 *  This file is part of ABISMAL.
 *
 *  ABISMAL is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ABISMAL is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 */

#ifndef ABISMAL_INDEX_HPP
#define ABISMAL_INDEX_HPP

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <algorithm>
#include <deque>
#include <bitset>
#include <cassert>

#include "smithlab_utils.hpp"

typedef std::vector<size_t> Genome;
static inline char random_base() {return "ACGT"[rand() & 3];}

namespace seed {
  // number of positions in the hashed portion of the seed
  static const uint32_t key_weight = 25u;
  static const uint32_t n_dp_positions = 30u;

  // window in which we select the best k-mer. The longer it is,
  // the longer the minimum read length that guarantees an exact
  // match will be mapped
  static const uint32_t window_size = 10u;

  // number of positions to sort within buckets
  static const uint32_t n_sorting_positions = 200u;

  static const size_t hash_mask = (1ull << seed::key_weight) - 1;

  // the purpose of padding the left and right ends of the
  // concatenated genome is so that later we can avoid having to check
  // the (unlikely) case that a read maps partly off either end of the
  // genome.
  static const size_t padding_size = 1048576ull; // 2^20

  void read(FILE* in);
  void write(FILE* out);
};

struct ChromLookup {
  std::vector<std::string> names;
  std::vector<uint32_t> starts;

  void
  get_chrom_idx_and_offset(const uint32_t pos,
                           uint32_t &chrom_idx,
                           uint32_t &offset) const;
  bool
  get_chrom_idx_and_offset(const uint32_t pos,
                           const uint32_t readlen,
                           uint32_t &chrom_idx,
                           uint32_t &offset) const;

  uint32_t
  get_pos(const std::string &chrom, const uint32_t offset) const;
  uint32_t
  get_genome_size() const {return starts.back();}

  FILE * read(FILE *in);
  std::istream & read(std::istream &in);
  void read(const std::string &infile);

  FILE * write(FILE *out) const;
  std::ostream & write(std::ostream &out) const;
  void write(const std::string &outfile) const;

  std::string tostring() const;
};

template <class G>
void
load_genome(const std::string &genome_file, G &genome, ChromLookup &cl) {
  std::ifstream in(genome_file);
  if (!in)
    throw std::runtime_error("bad genome file: " + genome_file);

  const auto begin_pos = in.tellg();
  in.seekg(0, std::ios_base::end);
  const size_t file_size = in.tellg() - begin_pos;
  in.seekg(0, std::ios_base::beg);

  genome.clear();
  // pad on at start; the space for padding at the end will be
  // available because of the newlines and chromosome names
  genome.reserve(file_size + seed::padding_size);

  // pad the start of the concatenated sequence
  cl.names.push_back("pad_start");
  for (size_t i = 0; i < seed::padding_size; ++i)
    genome.push_back('Z');
  cl.starts.push_back(genome.size());

  std::string line;
  while (getline(in, line))
    if (line[0] != '>') {
      for (auto it(begin(line)); it != end(line); ++it) {
        if (base2int(*it) == 4) // non-acgts become random bases
          *it = random_base();
      }
      copy(std::begin(line), std::end(line), std::back_inserter(genome));
    }
    else {
      cl.names.push_back(line.substr(1, line.find_first_of(" \t") - 1));
      cl.starts.push_back(genome.size());
    }

  // now pad the end of the concatenated sequence
  cl.names.push_back("pad_end");
  cl.starts.push_back(genome.size());
  for (size_t i = 0; i < seed::padding_size; ++i)
    genome.push_back('Z');
  cl.starts.push_back(genome.size());
}

std::ostream &
operator<<(std::ostream &out, const ChromLookup &cl);

struct AbismalIndex {

  static bool VERBOSE;

  uint32_t max_candidates; // number of hits at which seeds are expanded

  // the default PE heap size estimated from the genome k-mer
  // frequencies
  uint32_t pe_heap_size; // number of candidates kept on PE reads

  size_t counter_size; // number of kmers indexed
  size_t index_size; // size of the index

  std::vector<uint32_t> index; // genome positions for each k-mer
  std::vector<uint32_t> counter; // offset of each k-mer in "index"
  Genome genome; // the genome
  ChromLookup cl;

  void create_index(const std::string &genome_file);

  /* count how many positions must be stored for each hash value */
  void get_bucket_sizes(std::vector<bool> &keep, const uint32_t word_size);

  /* get index statistics used in mapping */
  void calc_mapping_parameters(const bool sens);

  /* how much RAM is needed to map reads*/
  double estimate_ram();

  /* select genome positions by dp*/
  void compress_dp(std::vector<bool> &keep);

  /* put genome positions in the appropriate buckets */
  void hash_genome(std::vector<bool> &keep);

  /* Sort each bucket, if the seed length is more than 26, then use
   * binary search for the rest part of the seed */
  void sort_buckets();

  /* convert the genome to 4-bit encoding */
  void encode_genome(const std::vector<uint8_t> &input_genome);

  void write(const std::string &index_file) const;
  void read(const std::string &index_file);

  static std::string internal_identifier;
  AbismalIndex() {}
};

// A/T nucleotide to 1-bit value (0100 | 0001 = 5) is for A or G.
inline uint32_t
get_bit(const uint8_t nt) {return (nt & 5) == 0;}

inline void
shift_hash_key(const uint8_t c, uint32_t &hash_key) {
  hash_key = (((hash_key << 1) | get_bit(c)) & seed::hash_mask);
}

// get the hash value for a k-mer (specified as some iterator/pointer)
// and the encoding for the function above
template <class T>
inline void
get_1bit_hash(T r, uint32_t &k) {
  const auto lim = r + seed::key_weight;
  k = 0;
  while (r != lim) {
    k = ((k << 1) | get_bit(*r));
    ++r;
  }
}

#endif
