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

#include <unordered_set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "AbismalSeed.hpp"

typedef std::vector<char> Genome;

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

  const size_t begin_pos = in.tellg();
  in.seekg(0, std::ios_base::end);
  const size_t file_size = in.tellg() - begin_pos;
  in.seekg(0, std::ios_base::beg);

  genome.clear();
  genome.reserve(file_size);

  std::string line;
  while (getline(in, line))
    if (line[0] != '>')
      copy(std::begin(line), std::end(line), std::back_inserter(genome));
    else {
      cl.names.push_back(line.substr(1, line.find_first_of(" \t")));
      cl.starts.push_back(genome.size());
    }
  cl.starts.push_back(genome.size());
}

std::ostream &
operator<<(std::ostream &out, const ChromLookup &cl);

struct AbismalIndex {

  static bool VERBOSE;
  static uint32_t valid_bucket_limit;
  static uint32_t max_N_per_seed;

  uint32_t counter_size; // number of kmers indexed
  uint32_t index_size; // size of the index

  std::vector<uint32_t> index; // genome positions for each k-mer
  std::vector<uint32_t> counter; // offset of each k-mer in "index"
  std::vector<char> genome; // the genome
  ChromLookup cl;

  /* count how many positions must be stored for each hash value */
  void
  get_bucket_sizes(std::unordered_set<uint32_t> &big_buckets);

  /* put genome positions in the appropriate buckets */
  void
  hash_genome(const std::unordered_set<uint32_t>& big_buckets);

  /* Sort each bucket, if the seed length is more than 12, then use
   * binary search for the rest part of the seed */
  void sort_buckets();

  void write(const std::string &index_file) const;
  void read(const std::string &index_file);

  static std::string internal_identifier;
};

#endif
