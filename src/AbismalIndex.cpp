/* Cmpyright (C) 2018 Andrew D. Smith
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

using std::string;
using std::vector;
using std::unordered_set;
using std::cerr;
using std::endl;
using std::runtime_error;
using std::sort;
using std::cout;
using std::min;

bool AbismalIndex::VERBOSE = false;
uint32_t AbismalIndex::max_invalid_per_seed = 0;
uint32_t AbismalIndex::deadzone_kmer_length = 100;
uint32_t AbismalIndex::valid_bucket_limit = 500000;

string AbismalIndex::internal_identifier = "AbismalIndex";

using genome_iterator = genome_four_bit_itr;

template <typename nuc_type>
bool
invalid_base(const nuc_type &x) {return x == 0;/*Z*/}

inline void
shift_hash_key_4bit(const uint8_t c, size_t &hash_key) {
  hash_key = ((hash_key << 1) | get_bit_4bit(c)) & seed::hash_mask;
}

void
AbismalIndex::encode_genome() {
  if (VERBOSE)
    cerr << "[encoding genome]" << endl;
  auto updated_end =
    encode_dna_four_bit(begin(genome), end(genome), begin(genome));
  genome.erase(updated_end, end(genome));
  vector<uint8_t>(genome).swap(genome); // shrink to fit
}

void
AbismalIndex::get_bucket_sizes(unordered_set<uint32_t> &big_buckets) {
  if (VERBOSE)
    cerr << "[allocating bucket counter]" << endl;
  counter_size = (1u << seed::key_weight);
  // the "counter" has an additional entry for convenience
  counter.resize(counter_size + 1, 0);

  const size_t lim = cl.get_genome_size() - seed::n_seed_positions;
  ProgressBar progress(lim, "counting bucket sizes");
  if (VERBOSE)
    progress.report(cerr, 0);

  // start counting the Ns in the seed region
  genome_iterator start_invalid_counter(begin(genome));
  auto end_invalid_counter(start_invalid_counter);
  uint32_t invalid_count = 0;
  while (end_invalid_counter !=
         start_invalid_counter + (seed::n_seed_positions - 1))
    invalid_count += invalid_base(*end_invalid_counter++);

  // start building up the hash key
  genome_iterator gi(begin(genome));
  const auto gi_lim(gi + (seed::key_weight - 1));
  size_t hash_key = 0;
  while (gi != gi_lim)
    shift_hash_key_4bit(*gi++, hash_key);

  for (size_t i = 0; i < lim; ++i) {
    if (VERBOSE && progress.time_to_report(i))
      progress.report(cerr, i);
    shift_hash_key_4bit(*gi++, hash_key);
    invalid_count += invalid_base(*end_invalid_counter++);
    if (invalid_count <= max_invalid_per_seed)
      counter[hash_key]++;
    invalid_count -= invalid_base(*start_invalid_counter++);
  }
  if (VERBOSE)
    progress.report(cerr, lim);

  if (VERBOSE)
    cerr << "[erasing big buckets]" << endl;
  for (uint32_t i = 0; i < counter_size; ++i)
    if (counter[i] > valid_bucket_limit) {
      counter[i] = 0;
      big_buckets.insert(i);
    }
  if (VERBOSE)
    cerr << "[erased " << big_buckets.size() << " buckets]" << endl;

  if (VERBOSE)
    cerr << "[computing bucket locations]" << endl;
  std::partial_sum(begin(counter), end(counter), begin(counter));
  index_size = counter[counter_size];
}

void
AbismalIndex::hash_genome(const unordered_set<uint32_t> &big_buckets) {
  if (VERBOSE)
    cerr << "[allocating hash table]" << endl;
  index.resize(index_size, 0);

  // ADS: make sure this works even if genome super small
  const size_t lim = cl.get_genome_size() - seed::n_seed_positions;
  ProgressBar progress(lim, "hashing genome");

  // start counting Ns in the seed windows
  genome_iterator start_invalid_counter(begin(genome));
  auto end_invalid_counter(start_invalid_counter);
  uint32_t invalid_count = 0;
  while (end_invalid_counter !=
         start_invalid_counter + (seed::n_seed_positions - 1))
    invalid_count += invalid_base(*end_invalid_counter++);

  // start building up the hash key
  genome_iterator gi(begin(genome));
  const auto gi_lim(gi + (seed::key_weight - 1));
  size_t hash_key = 0;
  while (gi != gi_lim)
    shift_hash_key_4bit(*gi++, hash_key);

  // begin(counter) + hash_key
  for (size_t i = 0; i < lim; ++i) {
    if (VERBOSE && progress.time_to_report(i))
      progress.report(cerr, i);
    invalid_count += invalid_base(*end_invalid_counter++);
    shift_hash_key_4bit(*gi++, hash_key);
    if (invalid_count <= max_invalid_per_seed &&
        big_buckets.find(hash_key) == end(big_buckets))
      index[--counter[hash_key]] = i;

    invalid_count -= invalid_base(*start_invalid_counter++);
  }
  if (VERBOSE)
    progress.report(cerr, lim);
}

struct BucketLessFour {
  BucketLessFour(const Genome &g, const uint32_t off) : 
    g_start(begin(g)), offset(off) {}
  bool operator()(const uint32_t a, const uint32_t b) const {
    auto idx1(g_start + a);
    auto lim1(g_start + a + offset);
    auto idx2(g_start + b);
    while (idx1 != lim1) {
      if (*idx1 != *idx2) return *idx1 < *idx2;
      ++idx1, ++idx2;
    }
    return false;
  }
  const genome_iterator g_start;
  const uint32_t offset;
};


struct BucketLessTwo {
  BucketLessTwo(const Genome &g) : g_start(begin(g)) {}
  bool operator()(const uint32_t a, const uint32_t b) const {
    auto idx1(g_start + a + seed::key_weight);
    auto lim1(g_start + a + seed::n_solid_positions);
    auto idx2(g_start + b + seed::key_weight);
    while (idx1 != lim1) {
      const char c1 = get_bit_4bit(*(idx1++));
      const char c2 = get_bit_4bit(*(idx2++));
      if (c1 != c2) return c1 < c2;
    }
    return false;
  }
  const genome_iterator g_start;
};

void
AbismalIndex::sort_buckets(const sort_type st) {
  const vector<uint32_t>::iterator b(begin(index));
  if (VERBOSE)
    cerr << "[sorting buckets by " << ((st == four_letter) ? "4":"2")
         << " letters]" << endl;

  if (st == four_letter) {
    const BucketLessFour bucket_less(genome, deadzone_kmer_length);
#pragma omp parallel for
    for (size_t i = 0; i < counter_size; ++i)
      if (counter[i + 1] > counter[i] + 1)
        sort(b + counter[i], b + counter[i + 1], bucket_less);
  }

  else {
    const BucketLessTwo bucket_less(genome);
#pragma omp parallel for
    for (size_t i = 0; i < counter_size; ++i)
      if (counter[i + 1] > counter[i] + 1)
        sort(b + counter[i], b + counter[i + 1], bucket_less);
  }
}

// GS TODO: add this to option parser and pass as parameter, or make it
// part of the AbismalIndex class
struct BucketEqual {
  BucketEqual(const Genome &g, const uint32_t off) :
    g_start(begin(g)), offset(off) {}
  bool operator()(const uint32_t a, const uint32_t b) const {
    auto idx1(g_start + a + offset);
    const auto lim1(g_start + a);
    auto idx2(g_start + b + offset);
    while (idx1 != lim1)
      if (*(--idx1) != *(--idx2))
        return false;
    return true;
  }
  const genome_iterator g_start;
  const uint32_t offset;
};

void
AbismalIndex::remove_big_buckets(const size_t max_candidates) {
  // first mark the positions to keep
  if (VERBOSE)
    cerr << "[finding big buckets]" << endl;
  const BucketEqual bucket_equal(genome, deadzone_kmer_length);
  const auto b(begin(index));
  vector<bool> keep(index_size, true);

#pragma omp parallel for
  for (size_t i = 0; i < counter_size; ++i) {
    uint32_t j = counter[i];
    while (j < counter[i+1]) {
      uint32_t k = j + 1;
      const uint32_t first_idx = *(b + j);
      while (k < counter[i+1] && bucket_equal(first_idx, *(b + k))) ++k;
      if (k - j > max_candidates)
        fill(begin(keep) + j, begin(keep) + k, false);
      j = k;
    }
  }

  // remove the indices that correspond to large buckets
  size_t num_removed = 0;
  if (VERBOSE)
    cerr << "[removing big buckets]" << endl;
  size_t j = 0;
  for (size_t i = 0; i < index_size; ++i) {
    if (keep[i])
      index[j++] = index[i];
    num_removed += !keep[i];
  }
  if (VERBOSE)
    cerr << "[removed " << num_removed << " positions from the index]" << endl;
  index.resize(j);
  index_size = j;

  // update the hash table to the positions of retained indices
  uint32_t total = 0;
  auto k_beg(begin(keep));
  uint32_t prev_counter = 0;
  for (size_t i = 1; i < counter.size(); ++i) {
    total += count(k_beg + prev_counter, k_beg + counter[i], false);
    prev_counter = counter[i];
    counter[i] -= total;
  }
}


/* ADS: original io functions using streams didn't work on mac os, so
   removed. Using the C file I/O rather than the streams...
*/

static void
write_internal_identifier(FILE *out) {
  if (fwrite((char*)&AbismalIndex::internal_identifier[0], 1,
             AbismalIndex::internal_identifier.size(), out) !=
      AbismalIndex::internal_identifier.size())
    throw runtime_error("failed writing index identifier");
}

void
AbismalIndex::write(const string &index_file) const {

  FILE *out = fopen(index_file.c_str(), "wb");
  if (!out)
    throw runtime_error("cannot open output file " + index_file);

  write_internal_identifier(out);

  const AbismalSeed the_seed(seed::n_seed_positions,
                             seed::key_weight,
                             seed::n_solid_positions,
                             seed::solid_positions);

  the_seed.write(out);

  cl.write(out);
  if (fwrite((char*)&genome[0], 1, genome.size(), out) != genome.size() ||
      fwrite((char*)&counter_size, sizeof(uint32_t), 1, out) != 1 ||
      fwrite((char*)&index_size, sizeof(uint32_t), 1, out) != 1 ||
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

  AbismalSeed the_seed;
  the_seed.read(in);

  const AbismalSeed expected_seed(seed::n_seed_positions,
                                  seed::key_weight,
                                  seed::n_solid_positions,
                                  seed::solid_positions);

  if (!(the_seed == expected_seed))
    throw runtime_error("inconsistent seeds (expected, indexed):\n" +
                        expected_seed.tostring() + "\n" +
                        the_seed.tostring());

  cl.read(in);

  const size_t genome_to_read = (cl.get_genome_size() + 1)/2;
  vector<uint8_t> enc_genome(genome_to_read);
  /* read the 4-bit encoded genome */
  if (fread((char*)&enc_genome[0], 1, genome_to_read, in) != genome_to_read)
    throw runtime_error(error_msg);

  /* expand the 4-bit encoded genome to 8-bit */
  genome.resize(cl.get_genome_size() + (cl.get_genome_size() % 2));
  auto d_itr(begin(genome)); //decoded iterator
  for (auto e_itr(begin(enc_genome)); e_itr != end(enc_genome); ++e_itr) {
    *d_itr++ = *e_itr & 15u; // low
    *d_itr++ = (*e_itr >> 4) & 15u; // high
  }
  vector<uint8_t>().swap(enc_genome); // release space

  /* read the sizes of counter and index vectors */
  if (fread((char*)&counter_size, sizeof(uint32_t), 1, in) != 1 ||
      fread((char*)&index_size, sizeof(uint32_t), 1, in) != 1)
    throw runtime_error(error_msg);

  /* allocate then read the counter vector */
  counter = vector<uint32_t>(counter_size + 1);
  if (fread((char*)(&counter[0]), sizeof(uint32_t),
            (counter_size + 1), in) != (counter_size + 1))
    throw runtime_error(error_msg);

  /* allocate then read the index vector */
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

  if (idx == begin(starts))
    throw std::runtime_error("bad chrom position: " + std::to_string(pos));

  --idx;

  chrom_idx = idx - begin(starts);
  offset = pos - starts[chrom_idx];
  return (pos + readlen <= starts[chrom_idx + 1]);
}

