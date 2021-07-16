/* Copyright (C) 2018 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * This file is part of ABISMAL.
 *
 * ABISMAL is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ABISMAL is distributed in the hope that it will be useful, but
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
#include <unordered_set>
#include <stdexcept>
#include <iterator>
#include <algorithm>
#include <deque>
#include <parallel/algorithm>

using std::string;
using std::vector;
using std::unordered_set;
using std::cerr;
using std::endl;
using std::runtime_error;
using std::sort;
using std::cout;
using std::min;
using std::deque;
using std::fill;
using std::to_string;

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

  vector<bool> pos_to_keep;
  pos_to_keep.resize(cl.get_genome_size());
  fill(begin(pos_to_keep), end(pos_to_keep), true);
  fill(begin(pos_to_keep), begin(pos_to_keep) + seed::padding_size, false);
  fill(end(pos_to_keep) - seed::padding_size, end(pos_to_keep), false);

  compress_minimizers(pos_to_keep);
  hash_genome(pos_to_keep);
  sort_buckets();
}


void
AbismalIndex::encode_genome(const vector<uint8_t> &input_genome) {
  if (VERBOSE)
    cerr << "[encoding genome]" << endl;

  genome.resize(input_genome.size()/16 + (input_genome.size()%16 != 0));
  encode_dna_four_bit(begin(input_genome), end(input_genome), begin(genome));
}

double
estimate_ram(const size_t counter_size,
             const size_t bases_indexed,
             const size_t genome_size) {
  size_t num_bytes =
    sizeof(size_t)*(genome_size) +
    sizeof(uint32_t)*(counter_size + bases_indexed);
           
  return num_bytes/(1024.0*1024.0*1024.0);
}

void
AbismalIndex::get_bucket_sizes(vector<bool> &keep) {
  counter.clear();

  // the "counter" has an additional entry for convenience
  counter_size = (1ull << seed::key_weight);
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
  size_t hash_key = 0;
  while (gi != gi_lim)
    shift_hash_key(*gi++, hash_key);

  for (size_t i = genome_st; i < lim; ++i) {
    if (VERBOSE && progress.time_to_report(i))
      progress.report(cerr, i);

    shift_hash_key(*gi++, hash_key);
    counter[hash_key] += keep[i];
  }
  if (VERBOSE)
    progress.report(cerr, lim);

  double mean = 0;
  size_t count = 0;
  if (VERBOSE)
    cerr << "[calculating index statistics]\n";

  // mean
  auto counter_end(end(counter));
  for (auto it(begin(counter)); it != counter_end; ++it) {
    mean += *it;
    count += (*it != 0);
  }

  mean = mean/static_cast<double>(count);
  vector<uint32_t> copy_of_counter;
  copy(begin(counter), end(counter), std::back_inserter(copy_of_counter));

  __gnu_parallel::sort(begin(copy_of_counter), end(copy_of_counter));
  max_candidates = copy_of_counter[
    static_cast<size_t>((1 - seed::overrep_kmer_quantile)*
    static_cast<double>(counter.size()))
  ];

  if (VERBOSE) {
    const size_t bases_indexed = std::accumulate(begin(counter), end(counter), 0ULL);
    cerr << std::fixed << std::setprecision(3)
         << "[bases indexed: " << bases_indexed/1000000 << "M"
         << ", mean: " << static_cast<size_t>(mean)
         << ", required RAM: "
         << estimate_ram(counter.size(), bases_indexed, genome.size()) << " GB"
         << ", max count: " << copy_of_counter.back()
         << ", cutoff: " << max_candidates << "]" << endl;
  }
}

void
AbismalIndex::hash_genome(vector<bool> &keep) {
  get_bucket_sizes(keep);

  if (VERBOSE)
    cerr << "[allocating hash table]" << endl;

  std::partial_sum(begin(counter), end(counter), begin(counter));
  index_size = counter[counter_size];
  index.resize(index_size, 0);

  const size_t genome_st = seed::padding_size;
  const size_t lim = cl.get_genome_size() - seed::key_weight - seed::padding_size;
  ProgressBar progress(lim, "hashing genome");

  // start building up the hash key
  genome_iterator gi(begin(genome));
  gi = gi + genome_st;

  const auto gi_lim(gi + (seed::key_weight - 1));
  size_t hash_key = 0;
  while (gi != gi_lim)
    shift_hash_key(*gi++, hash_key);

  for (size_t i = genome_st; i < lim; ++i) {
    if (VERBOSE && progress.time_to_report(i))
      progress.report(cerr, i);

    shift_hash_key(*gi++, hash_key);
    if (keep[i])
      index[--counter[hash_key]] = i;
  }
  if (VERBOSE)
    progress.report(cerr, lim);
}

void
AbismalIndex::compress_minimizers(vector<bool> &keep) {
  fill(begin(keep), end(keep), false);
  deque<kmer_loc> window_kmers;

  const size_t genome_st = seed::padding_size;
  const size_t lim = cl.get_genome_size() - seed::key_weight - seed::padding_size;

  ProgressBar progress(lim, "indexing (" +
    to_string(seed::window_size) + "," + to_string(seed::key_weight) +
    ") mins"
  );

  // start building up the hash key
  genome_iterator gi(begin(genome));
  gi = gi + genome_st;
  const auto gi_lim(gi + (seed::key_weight - 1));
  size_t hash_key = 0;
  while (gi != gi_lim)
    shift_hash_key(*gi++, hash_key);

  // get the first W minimizers
  for (size_t i = genome_st; i < genome_st + seed::window_size; ++i) {
    shift_hash_key(*gi++, hash_key);
    add_kmer(kmer_loc(hash_key, i), seed::window_size, window_kmers);
  }

  for (size_t i = genome_st + seed::window_size; i < lim; ++i) {
    if (VERBOSE && progress.time_to_report(i))
      progress.report(cerr, i);
    
    keep[window_kmers.front().loc] = true;

    // update minimizers using k-mer from current base
    shift_hash_key(*gi++, hash_key);
    add_kmer(kmer_loc(hash_key, i), seed::window_size, window_kmers);
  }

  for (auto it(begin(window_kmers)); it != end(window_kmers) && *it == window_kmers.front(); ++it)
    keep[it->loc] = true;

  if (VERBOSE)
    progress.report(cerr, lim);
}

struct BucketLess {
  BucketLess(const Genome &g) : g_start(begin(g)) {}
  bool operator()(const uint32_t a, const uint32_t b) const {
    auto idx1(g_start + a + seed::key_weight);
    auto lim1(g_start + a + seed::n_sorting_positions);
    auto idx2(g_start + b + seed::key_weight);
    while (idx1 != lim1) {
      const char c1 = get_bit(*(idx1++));
      const char c2 = get_bit(*(idx2++));
      if (c1 != c2) return c1 < c2;
    }
    return false;
  }
  const genome_iterator g_start;
};

void
AbismalIndex::sort_buckets() {
  if (VERBOSE)
    cerr << "[sorting buckets]" << endl;
  const vector<uint32_t>::iterator b(begin(index));
  const BucketLess bucket_less(genome);
#pragma omp parallel for
  for (size_t i = 0; i < counter_size; ++i)
    if (counter[i + 1] > counter[i] + 1) {
      sort(b + counter[i], b + counter[i + 1], bucket_less);
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

  if (fwrite((char*)&genome[0], sizeof(size_t), genome.size(), out) != genome.size() ||
      fwrite((char*)&max_candidates, sizeof(size_t), 1, out) != 1 ||
      fwrite((char*)&counter_size, sizeof(size_t), 1, out) != 1 ||
      fwrite((char*)&index_size, sizeof(size_t), 1, out) != 1 ||
      fwrite((char*)(&counter[0]), sizeof(uint32_t),
             counter_size + 1, out) != (counter_size + 1) ||
      fwrite((char*)(&index[0]), sizeof(uint32_t),
             index_size, out) != index_size)
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

  const size_t genome_to_read = (cl.get_genome_size() + 15)/16;
  // read the 4-bit encoded genome
  genome.resize(genome_to_read);
  if (fread((char*)&genome[0], sizeof(size_t), genome_to_read, in) != genome_to_read)
    throw runtime_error(error_msg);

  // read k-mer mean and standard deviation
  if (fread((char*)&max_candidates, sizeof(size_t), 1, in) != 1)
    throw runtime_error(error_msg);

  // read the sizes of counter and index vectors
  if (fread((char*)&counter_size, sizeof(size_t), 1, in) != 1 ||
      fread((char*)&index_size, sizeof(size_t), 1, in) != 1)
    throw runtime_error(error_msg);

  // allocate then read the counter vector
  counter = vector<uint32_t>(counter_size + 1);
  if (fread((char*)(&counter[0]), sizeof(uint32_t),
            (counter_size + 1), in) != (counter_size + 1))
    throw runtime_error(error_msg);

  // allocate the read the index vector
  index = vector<uint32_t>(index_size);
  if (fread((char*)(&index[0]), sizeof(uint32_t), index_size, in) != index_size)
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

  assert(idx != begin(starts));

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
