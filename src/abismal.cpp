/* Copyright (C) 2018-2019 Andrew D. Smith
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

#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"
#include "OptionParser.hpp"
#include "zlib_wrapper.hpp"
#include "MappedRead.hpp"

#include "AbismalIndex.hpp"

#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>

#include <omp.h>

using std::vector;
using std::runtime_error;
using std::string;
using std::cerr;
using std::endl;
using std::cout;
using std::transform;
using std::numeric_limits;
using std::ostream;
using std::ofstream;
using std::count;

// type returned from atomic (e.g. nucleotide) comparisons
typedef bool cmp_t;

struct se_map_result {
  uint32_t pos;
  uint16_t diffs;
  char strand;
  bool ambig;

  se_map_result() : pos(0), diffs(max_diffs + 1), strand('+'), ambig(false) {}

  se_map_result(const uint32_t p, const uint16_t d, const char s) :
    pos(p), diffs(d), strand(s), ambig(false) {}

  bool operator<(const se_map_result &rhs) const {return diffs < rhs.diffs;}

  bool valid() const {return diffs <= max_diffs;}

  void reset() {diffs = max_diffs + 1; ambig = false;}

  void
  update(const uint16_t proposed_diffs, const uint32_t proposed_pos) {
    if (proposed_diffs < diffs) {
      diffs = proposed_diffs;
      pos = proposed_pos;
      strand = global_strand;
      ambig = false;
    }
    else if (proposed_diffs == diffs && pos != proposed_pos)
      ambig = true;
  }

  bool optimal(uint32_t seed_number = 0) const {
    return ambig && (diffs == 0 || (diffs == 1 && seed_number > 0));
  }

  uint32_t get_cutoff() const {return diffs;}

  static char global_strand;
  static uint16_t max_diffs;
};
uint16_t se_map_result::max_diffs = 6;
char se_map_result::global_strand = '+';

struct ReadLoader {
  ReadLoader(const string &fn,
             const size_t bs = numeric_limits<size_t>::max()) :
    filename(fn), batch_size(bs) {
    in = new igzfstream(fn);
    if (!in || !(*in))
      throw runtime_error("bad reads file: " + filename);
  }
  ~ReadLoader() {delete in;}
  bool good() const {return bool(*in);}
  size_t get_current_byte() const {return gztell(in->fileobj);}

  void
  load_reads(vector<string> &names, vector<string> &reads) {
    static const size_t reserve_size = 250;

    reads.clear();
    names.clear();

    size_t line_count = 0;
    string line;
    line.reserve(reserve_size);
    while (line_count < 4*batch_size && bool(getline(*in, line))) {
      if (line_count % 4 == 0) {
        names.push_back(line.substr(1, line.find_first_of(" \t")));
      }
      else if (line_count % 4 == 1) {
        if (line.length() < min_length ||
            count(begin(line), end(line), 'N') > se_map_result::max_diffs)
          line.clear();
        reads.push_back(line);
      }
      ++line_count;
    }
    in->peek(); // needed in case batch_size exactly divides the
                // number of reads in the file
  }
  string filename;
  size_t batch_size;
  igzfstream *in;

  static uint32_t min_length;
};
uint32_t ReadLoader::min_length = 32;

char
flip_strand(const char s) {return s=='+' ? '-' : '+';}

void
format_se_map_result(const se_map_result &res, const ChromLookup &cl,
                     const string &read, const string &read_name,
                     const bool rc, ofstream &out) {
  uint32_t offset = 0, chrom_idx = 0;
  if (cl.get_chrom_idx_and_offset(res.pos, read.length(),
                                  chrom_idx, offset)) {
    char s;
    string seq;
    if (rc) {
      s = flip_strand(res.strand);
      seq = revcomp(read);
    }
    else {
      s = res.strand;
      seq = read;
    }
    out << cl.names[chrom_idx] << '\t'
        << offset << '\t'
        << offset + read.length() << '\t'
        << read_name << '\t'
        << res.diffs << '\t'
        << s << '\t'
        << seq << '\n';
  }
}

struct pe_map_result {
  uint32_t pos;
  uint16_t diffs;
  pe_map_result(const uint32_t p, const uint16_t d) : pos(p), diffs(d) {}
  pe_map_result() : pos(0), diffs(max_diffs + 1) {}
  bool operator==(const pe_map_result &rhs) const {
    return diffs == rhs.diffs && pos == rhs.pos;
  }
  bool operator<(const pe_map_result &rhs) const {return diffs < rhs.diffs;}
  void reset() {diffs = max_diffs + 1;}
  bool not_set() const {return (diffs == max_diffs + 1);}

  static bool
  pos_less(const pe_map_result &a, const pe_map_result &b) {
    return a.pos < b.pos;
  }

  static uint16_t max_diffs;
};
uint16_t pe_map_result::max_diffs = 6;


struct pe_map_result_merged {
  uint32_t pos1;
  uint32_t pos2;
  uint32_t frag_len;
  uint16_t diffs;
  char strand;
  bool ambig;
  pe_map_result_merged(const pe_map_result &a,
                       const pe_map_result &b, const char s,
                       const uint32_t l) :
    pos1(a.pos), pos2(b.pos), frag_len(l),
    diffs(a.diffs + b.diffs), strand(s), ambig(false) {}
  pe_map_result_merged() :
    pos1(0), pos2(0), frag_len(0), diffs(max_diffs + 1),
    strand('+'), ambig(false) {}
  bool not_set() const {return (diffs == max_diffs + 1);}
  bool operator<(const pe_map_result_merged &rhs) const {
    return diffs < rhs.diffs ||
                   (diffs == rhs.diffs && frag_len < rhs.frag_len);
  }
  bool operator==(const pe_map_result_merged &rhs) const {
    return diffs == rhs.diffs && frag_len > 0 && frag_len == rhs.frag_len;
  }
  void reset() {
    diffs = max_diffs + 1;
  }

  static uint16_t max_diffs;
  static uint32_t min_dist;
  static uint32_t max_dist;
};
uint16_t pe_map_result_merged::max_diffs = 6;
uint32_t pe_map_result_merged::min_dist = 32;
uint32_t pe_map_result_merged::max_dist = 3000;

void
format_pe_map_result_merged(const pe_map_result_merged &res,
                            const ChromLookup &cl,
                            string &read1, string &read2,
                            const string &read_name1,
                            ofstream &out) {

  uint32_t offset = 0, chrom_idx = 0;
  if (cl.get_chrom_idx_and_offset(res.pos1, res.frag_len, chrom_idx, offset)) {
    string read;
    if (res.frag_len < read1.length()) {
      read.swap(read1);
      read.resize(res.frag_len);
    }
    else {
      revcomp_inplace(read2);
      if (res.frag_len < read2.length()) {
        read.swap(read2);
        read.resize(res.frag_len);
      }
      else {
        read.resize(res.frag_len, 'N');
        copy(begin(read1), end(read1), begin(read));
        copy(begin(read2), end(read2), end(read) - read2.length());
      }
    }
    out << cl.names[chrom_idx] << '\t'
        << offset << '\t'
        << offset + read.length() << '\t'
        << "FRAG:" << read_name1 << '\t'
        << res.diffs << '\t'
        << res.strand << '\t'
        << read << '\n';
  }
}


struct pe_candidates {

  pe_candidates() : v(vector<pe_map_result>(max_size)), sz(0) {}

  bool empty() const {return sz == 0;}
  bool full() const {return sz == max_size;}
  void reset() {
    for (auto i(begin(v)); i != begin(v) + sz; ++i)
      i->reset();
    sz = 0;
  }
  uint16_t get_cutoff() const {return v.front().diffs;}

  void
  update(const uint16_t proposed_diffs, const uint32_t proposed_pos) {
    if (full()) {
      if (proposed_diffs < v.front().diffs) {
        std::pop_heap(begin(v), end(v));
        v.back().diffs = proposed_diffs;
        v.back().pos = proposed_pos;
        std::push_heap(begin(v), end(v));
      }
    }
    else if (proposed_diffs <= pe_map_result::max_diffs) {
      v[sz].diffs = proposed_diffs;
      v[sz].pos = proposed_pos;
      ++sz;
      std::push_heap(begin(v), begin(v) + sz);
    }
  }

  bool optimal(uint32_t seed_number = 0) const {
    return full() &&
      (v.front().diffs == 0 ||
       (v.front().diffs == 1 && seed_number > 0));
  }

  void prepare_for_mating() {
    std::sort(begin(v), begin(v) + sz, pe_map_result::pos_less);
    sz = std::unique(begin(v), begin(v) + sz) - begin(v);
  }

  vector<pe_map_result> v;
  uint32_t sz;

  static uint32_t max_size;
  static char global_strand;
};
uint32_t pe_candidates::max_size = 20;
char pe_candidates::global_strand = '+';


static void
get_best_se_map_result(const char strand, const pe_candidates &pres,
                       se_map_result &res) {
  auto i(begin(pres.v));
  auto lim(begin(pres.v) + pres.sz);
  for (; i != lim; ++i) {
    if (i->diffs < res.diffs)
      res = se_map_result(i->pos, i->diffs, strand);
    else if (i->diffs == res.diffs && i->pos != res.pos)
      res.ambig = true;
  }
}


static void
get_best_pair(const char strand, pe_candidates &res1, pe_candidates &res2,
              pe_map_result_merged &best, const size_t r2_len) {

  res1.prepare_for_mating();
  res2.prepare_for_mating();

  auto j1 = begin(res1.v);
  auto j1_end = begin(res1.v) + res1.sz;
  auto j2_end = begin(res2.v) + res2.sz;

  for (auto j2(begin(res2.v)); j2 != j2_end; ++j2) {

    while (j1 != j1_end &&
           j1->pos + pe_map_result_merged::max_dist < j2->pos + r2_len)
      ++j1;

    while (j1 != j1_end &&
           j1->pos + pe_map_result_merged::min_dist <= j2->pos + r2_len) {
      const pe_map_result_merged p(*j1, *j2, strand,
                                   j2->pos + r2_len - j1->pos);
      if (p < best)
        best = p;
      else if (p == best)
        best.ambig = true;
      ++j1;
    }
  }
}

struct se_map_stats {
  se_map_stats() :
    total_reads(0), unique_reads(0), ambiguous_reads(0),
    unmapped_reads(0), skipped_reads(0) {}
  uint32_t total_reads;
  uint32_t unique_reads;
  uint32_t ambiguous_reads;
  uint32_t unmapped_reads;
  uint32_t skipped_reads;

  void update(const se_map_result &res) {
    ++total_reads;
    if (res.valid()) {
      if (!res.ambig) ++unique_reads;
      else ++ambiguous_reads;
    }
    else ++unmapped_reads;
  }

  string tostring(const size_t n_tabs = 0) const {
    static const string tab = "    ";
    string t;
    for (size_t i = 0; i < n_tabs; ++i) t += tab;
    std::ostringstream oss;
    oss << t     << "total_reads: " << total_reads << endl
        << t     << "mapped: " << endl
        << t+tab     << "percent_mapped: "
        << (100.0*(unique_reads + ambiguous_reads))/total_reads << endl
        << t+tab << "unique: " << unique_reads << endl
        << t+tab << "percent_unique: "
        << 100.0*unique_reads/total_reads << endl
        << t+tab << "ambiguous: " << ambiguous_reads << endl
        << t     << "unmapped: " << unmapped_reads << endl
        << t     << "skipped: " << skipped_reads << endl;
    return oss.str();
  }
};

struct pe_map_stats {
  pe_map_stats(const uint32_t min_d, const uint32_t max_d) :
    total_pairs(0), unique_pairs(0), ambiguous_pairs(0), unmapped_pairs(0),
    min_dist(min_d) {
    frag_len_hist.resize(max_d-min_dist+1, 0);
  }
  uint32_t total_pairs;
  uint32_t unique_pairs;
  uint32_t ambiguous_pairs;
  uint32_t unmapped_pairs;
  uint32_t min_dist;
  vector<uint32_t> frag_len_hist;
  se_map_stats end1_stats;
  se_map_stats end2_stats;

  void update_pair(const pe_map_result_merged &res, const bool allow_ambig) {
    ++total_pairs;
    if (res.ambig) ++ambiguous_pairs;
    else ++unique_pairs;
    if (allow_ambig || !res.ambig)
      ++frag_len_hist[res.frag_len - min_dist];
  }

  string tostring() const {
    std::ostringstream oss;
    static const string t = "    ";
    oss << "pairs:" << endl
        << t << "total_read_pairs: " << total_pairs << endl
        << t << "mapped:" << endl
        << t+t << "percent_mapped: "
        << (100.0*(unique_pairs + ambiguous_pairs))/total_pairs << endl
        << t+t << "unique: " << unique_pairs << endl
        << t+t << "percent_unique: " << 100.0*unique_pairs/total_pairs << endl
        << t+t << "ambiguous: " << ambiguous_pairs << endl
        << t << "unmapped: " << unmapped_pairs << endl
        << "mate1:" << endl
        << end1_stats.tostring(1)
        << "mate2:" << endl
        << end2_stats.tostring(1);

    oss << "frag_len_distribution:" << endl;
    double total = 0.0;
    for (size_t i = 0; i < frag_len_hist.size(); ++i) {
      oss << t << i + min_dist << ": " << frag_len_hist[i] << endl;
      total += ((i + min_dist)*frag_len_hist[i]);
    }
    oss << "frag_len_mean: " << total/unique_pairs << endl;
    return oss.str();
  }
};


static void
select_output(const ChromLookup &cl,
              const pe_map_result_merged &best,
              const se_map_result &pos1, const se_map_result &neg1,
              const se_map_result &pos2, const se_map_result &neg2,
              string &read1, const string &read_name1,
              string &read2, const string &read_name2,
              const bool allow_ambig, const bool rc,
              pe_map_stats &pe_stats, ofstream &out) {
  if (!best.not_set()) {
    if (allow_ambig || !best.ambig)
      format_pe_map_result_merged(best, cl, read1, read2, read_name1, out);
    pe_stats.update_pair(best, allow_ambig);
  }
  else {
    ++pe_stats.total_pairs;
    ++pe_stats.unmapped_pairs;
    const se_map_result *best_mate1 = pos1 < neg1 ? &pos1 : &neg1;
    const se_map_result *best_mate2 = pos2 < neg2 ? &pos2 : &neg2;
    pe_stats.end1_stats.update(*best_mate1);
    pe_stats.end2_stats.update(*best_mate2);
    if (read1.length() == 0) ++pe_stats.end1_stats.skipped_reads;
    if (read2.length() == 0) ++pe_stats.end2_stats.skipped_reads;
    if (best_mate1->valid() && (allow_ambig || !best_mate1->ambig))
      format_se_map_result(*best_mate1, cl, read1, read_name1, rc, out);
    if (best_mate2->valid() && (allow_ambig || !best_mate1->ambig))
      format_se_map_result(*best_mate2, cl, read2, read_name2, !rc, out);
  }
}


template <cmp_t (*compare)(const char, const char)>
uint16_t
full_compare(const uint16_t best_diffs,
             string::const_iterator read_itr,
             const string::const_iterator &read_end,
             Genome::const_iterator genome_itr) {
  uint16_t d = 0;
  while (d <= best_diffs && read_itr != read_end)
    d += compare(*read_itr++, *genome_itr++);
  return d;
}

inline cmp_t comp_ct(const char a, const char b) {
  return (a != b && !(b == 'C' && a == 'T')); // && b != 'N');
}

inline cmp_t comp_ga(const char a, const char b) {
  return (a != b && !(b == 'G' && a == 'A')); // && b != 'N');
}


template <cmp_t (*F)(const char, const char), class T>
void
check_hits(vector<uint32_t>::const_iterator start_idx,
           const vector<uint32_t>::const_iterator end_idx,
           const string::const_iterator read_start,
           const string::const_iterator read_end,
           const Genome::const_iterator genome_start, T &res) {

  for (; start_idx < end_idx && !res.optimal(0); ++start_idx) {
    const uint32_t pos = *start_idx;
    const uint16_t diffs =
      full_compare<F>(res.get_cutoff(), read_start, read_end,
                      genome_start + pos);
    res.update(diffs, pos);
  }
}


struct compare_bases {
  compare_bases(const Genome::const_iterator g_) : g(g_) {}
  bool operator()(const uint32_t mid, const uint32_t chr) const {
    return (get_bit(*(g + mid)) < chr);
  }
  const Genome::const_iterator g;
};


void
find_candidates(const string::const_iterator read_start,
                const Genome::const_iterator genome_start,
                const uint32_t read_lim, // not necessarily read len
                vector<uint32_t>::const_iterator &low,
                vector<uint32_t>::const_iterator &high) {
  size_t p = seed::key_weight;
  const size_t lim = std::min(read_lim, seed::n_solid_positions);
  while (p < lim) {
    auto first_1 = lower_bound(low, high, 1, compare_bases(genome_start + p));
    if (get_bit(*(read_start + p)) == 0) {
      if (first_1 == high) return; // need 0s; whole range is 0s
      high = first_1;
    }
    else {
      if (first_1 == low) return; // need 1s; whole range is 1s
      low = first_1;
    }
    ++p;
  }
}


template <cmp_t (*F)(const char, const char), class T>
void
process_seeds(vector<uint32_t> &hits,
              const uint32_t genome_size,
              const uint32_t max_candidates,
              const AbismalIndex &abismal_index,
              const Genome::const_iterator genome_start,
              const string &read, T &res) {

  const uint32_t readlen = read.length();
  const uint32_t genome_lim = genome_size - readlen;

  const string::const_iterator read_start(begin(read));
  const string::const_iterator read_end(end(read));
  const vector<uint32_t>::const_iterator index_st(begin(abismal_index.index));
  const vector<uint32_t>::const_iterator counter_st(begin(abismal_index.counter));

  const size_t valid_starts =
    readlen > seed::n_solid_positions ? readlen - seed::n_solid_positions : 0;
  const size_t shift_size = std::max(1ul, valid_starts/(seed::n_shifts - 1));

  for (uint32_t i = 0; i <= valid_starts && !res.optimal(i); i += shift_size) {

    hits.clear(); // assert(hits.capacity() == max_candidates);

    uint32_t k_low, k_high;
    get_1bit_hash_low_high(read_start + i, readlen - i, k_low, k_high);

    for (uint32_t k = k_low; k < k_high; ++k) {
      auto start_idx(index_st + *(counter_st + k));
      auto end_idx(index_st + *(counter_st + k+1));

      if (start_idx < end_idx) {
        find_candidates(read_start + i, genome_start,
                        readlen - i, start_idx, end_idx);

        if (end_idx - start_idx <= max_candidates) // optimization
          for (auto j(start_idx); j != end_idx &&
                 hits.size() < max_candidates /*optimization*/ ; ++j)
            if (*j >= i && *j <= genome_lim)
              hits.push_back(*j - i);
      }
    }
    if (!hits.empty())
      check_hits<F>(begin(hits), end(hits),
                    read_start, read_end, genome_start, res);
  }
}


template <cmp_t (*pos_cmp)(const char, const char),
          cmp_t (*neg_cmp)(const char, const char), class T>
void
map_single_ended_batch(const bool VERBOSE,
                       const vector<string> &reads,
                       const uint32_t max_candidates,
                       const AbismalIndex &abismal_index, vector<T> &res) {

  const uint32_t genome_size = abismal_index.genome.size();
  const Genome::const_iterator genome_start(begin(abismal_index.genome));

  T::global_strand = '+';
#pragma omp parallel
  {
    vector<uint32_t> hits;
    hits.reserve(max_candidates);
#pragma omp for
    for (size_t i = 0; i < reads.size(); ++i)
      if (reads[i].length() > 0)
        process_seeds<pos_cmp>(hits, genome_size, max_candidates,
                               abismal_index, genome_start, reads[i], res[i]);
  }

  T::global_strand = '-';
#pragma omp parallel
  {
    vector<uint32_t> hits;
    hits.reserve(max_candidates);
#pragma omp for
    for (size_t i = 0; i < reads.size(); ++i) {
      if (reads[i].length() > 0) {
        const string read_rc(revcomp(reads[i]));
        process_seeds<neg_cmp>(hits, genome_size, max_candidates,
                               abismal_index, genome_start, read_rc, res[i]);
      }
    }
  }
}


template <cmp_t (*pos_cmp)(const char, const char),
          cmp_t (*neg_cmp)(const char, const char)>
void
map_single_ended(const bool VERBOSE,
                 const string &reads_file,
                 const size_t batch_size,
                 const size_t max_candidates,
                 const AbismalIndex &abismal_index,
                 const bool allow_ambig,
                 const bool rc,
                 se_map_stats &se_stats,
                 ofstream &out) {

  const double start_time = omp_get_wtime();

  ReadLoader rl(reads_file, batch_size);

  ProgressBar progress(get_filesize(reads_file), "mapping reads");
  if (VERBOSE)
    progress.report(cerr, 0);

  while (rl.good()) {

    if (VERBOSE && progress.time_to_report(rl.get_current_byte()))
      progress.report(cerr, rl.get_current_byte());

    vector<string> names, reads;
    rl.load_reads(names, reads);

    vector<se_map_result> res(reads.size());
    map_single_ended_batch<pos_cmp, neg_cmp>(VERBOSE, reads,
                                             max_candidates, abismal_index, res);

    for (size_t i = 0 ; i < reads.size(); ++i) {
      se_stats.update(res[i]);
      if (reads[i].length() == 0) ++se_stats.skipped_reads;
      if (res[i].valid() && (allow_ambig || !res[i].ambig))
        format_se_map_result(res[i], abismal_index.cl,
                             reads[i], names[i], rc, out);
    }
  }
  if (VERBOSE)
    progress.report(cerr, get_filesize(reads_file));

  const double end_time = omp_get_wtime();
  if (VERBOSE)
    cerr << "[total mapping time: "
         << end_time - start_time << "s]" << endl;
}


template <cmp_t (*cmp)(const char, const char), class T>
void
map_paired_ended_batch(const bool VERBOSE,
                       const vector<string> &reads1,
                       const vector<string> &reads2,
                       const uint32_t max_candidates,
                       const AbismalIndex &abismal_index,
                       vector<T> &res1, vector<T> &res2) {

  const uint32_t genome_size = abismal_index.genome.size();
  const Genome::const_iterator genome_start(begin(abismal_index.genome));

#pragma omp parallel
  {
    vector<uint32_t> hits;
    hits.reserve(max_candidates);
#pragma omp for
    for (size_t i = 0; i < reads1.size(); ++i) {
      if (reads1[i].length() > 0)
        process_seeds<cmp>(hits, genome_size, max_candidates,
                           abismal_index, genome_start, reads1[i], res1[i]);
    }
  }

#pragma omp parallel
  {
    vector<uint32_t> hits;
    hits.reserve(max_candidates);
#pragma omp for
    for (size_t i = 0; i < reads2.size(); ++i) {
      if (reads2[i].length() > 0) {
        const string read_rc(revcomp(reads2[i]));
        process_seeds<cmp>(hits, genome_size, max_candidates,
                           abismal_index, genome_start, read_rc, res2[i]);
      }
    }
  }
}


template <cmp_t (*pos_cmp)(const char, const char),
          cmp_t (*neg_cmp)(const char, const char)>
void
map_paired_ended(const bool VERBOSE,
                 const string &reads_file1,
                 const string &reads_file2,
                 const size_t batch_size,
                 const size_t max_candidates,
                 const AbismalIndex &abismal_index,
                 const bool allow_ambig,
                 const bool rc,
                 pe_map_stats &pe_stats,
                 ofstream &out) {

  const double start_time = omp_get_wtime();

  ReadLoader rl1(reads_file1, batch_size);
  ReadLoader rl2(reads_file2, batch_size);

  vector<pe_candidates> res1(batch_size), res2(batch_size);
  vector<string> names1, reads1, names2, reads2;
  reads1.reserve(batch_size);
  names1.reserve(batch_size);
  reads2.reserve(batch_size);
  names2.reserve(batch_size);

  vector<pe_map_result_merged> bests(batch_size);
  vector<se_map_result> res_se_p1(batch_size), res_se_p2(batch_size);
  vector<se_map_result> res_se_n1(batch_size), res_se_n2(batch_size);

  ProgressBar progress(get_filesize(reads_file1), "mapping reads");
  if (VERBOSE)
    progress.report(cerr, 0);

  while (rl1.good() && rl2.good()) {

    if (VERBOSE && progress.time_to_report(rl1.get_current_byte()))
      progress.report(cerr, rl1.get_current_byte());

    rl1.load_reads(names1, reads1);
    rl2.load_reads(names2, reads2);

    const size_t n_reads = reads1.size();

#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i) {
      res1[i].reset();
      res2[i].reset();
      res_se_p1[i].reset();
      res_se_p2[i].reset();
      bests[i].reset();
    }
    map_paired_ended_batch<pos_cmp>(VERBOSE, reads1, reads2,
                                    max_candidates, abismal_index, res1, res2);

#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i) {
      get_best_pair('+', res1[i], res2[i], bests[i], reads2[i].length());
      if (bests[i].not_set()) {
        get_best_se_map_result('+', res1[i], res_se_p1[i]);
        get_best_se_map_result('-', res2[i], res_se_p2[i]);
      }
    }

#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i) {
      res1[i].reset();
      res2[i].reset();
      res_se_n1[i].reset();
      res_se_n2[i].reset();
    }
    map_paired_ended_batch<neg_cmp>(VERBOSE, reads2, reads1,
                                    max_candidates, abismal_index, res2, res1);

#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i) {
      get_best_pair('-', res2[i], res1[i], bests[i], reads1[i].length());
      if (bests[i].not_set()) {
        get_best_se_map_result('-', res1[i], res_se_n1[i]);
        get_best_se_map_result('+', res2[i], res_se_n2[i]);
      }
    }

    for (size_t i = 0 ; i < n_reads; ++i)
      select_output(abismal_index.cl, bests[i],
                    res_se_p1[i], res_se_n1[i], res_se_p2[i], res_se_n2[i],
                    reads1[i], names1[i], reads2[i], names2[i], allow_ambig,
                    rc, pe_stats, out);
  }
  if (VERBOSE)
    progress.report(cerr, get_filesize(reads_file1));
  const double end_time = omp_get_wtime();
  if (VERBOSE)
    cerr << "[total mapping time: "
         << end_time - start_time << "s)]" << endl;
}


int main(int argc, const char **argv) {

  try {

    string index_file;
    string outfile;
    bool VERBOSE = false;
    bool GA_conversion = false;
    bool allow_ambig = false;
    bool pbat_mode = false;
    uint32_t max_candidates = 3000;
    size_t max_diffs = 6;
    size_t batch_size = 1000000;
    size_t n_threads = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "map bisulfite converted reads",
                           "<reads-fq1> [<reads-fq2>]");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("index", 'i', "index file", true, index_file);
    opt_parse.add_opt("outfile", 'o', "output file", true, outfile);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
    opt_parse.add_opt("mismatches", 'm', "max allowed mismatches",
                      false, max_diffs);
    opt_parse.add_opt("shifts", 's', "number of seed shifts",
                      false, seed::n_shifts);
    opt_parse.add_opt("batch", 'b', "reads to load at once",
                      false, batch_size);
    opt_parse.add_opt("candidates", 'c', "max candidates for full comparison",
                      false, max_candidates);
    opt_parse.add_opt("max-mates", 'p', "max candidates as mates (pe mode)",
                      false, pe_candidates::max_size);
    opt_parse.add_opt("min-frag", 'l', "min fragment size (pe mode)",
                      false, pe_map_result_merged::min_dist);
    opt_parse.add_opt("max-frag", 'L', "max fragment size (pe mode)",
                      false, pe_map_result_merged::max_dist);
    opt_parse.add_opt("ambig", 'a', "report a posn for ambiguous mappers",
                      false, allow_ambig);
    opt_parse.add_opt("pbat", 'P', "input data follow the PBAT protocol",
                      false, pbat_mode);
    opt_parse.add_opt("a-rich", 'A', "indicates reads are a-rich (se mode)",
                      false, GA_conversion);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() < 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string reads_file = leftover_args.front();
    string reads_file2;
    bool paired_end = false;
    if (leftover_args.size() == 2) {
      paired_end = true;
      reads_file2 = leftover_args.back();
    }
    /****************** END COMMAND LINE OPTIONS *****************/

    se_map_result::max_diffs = max_diffs;
    pe_map_result::max_diffs = se_map_result::max_diffs;
    pe_map_result_merged::max_diffs = 2*se_map_result::max_diffs;

    omp_set_num_threads(n_threads);

    AbismalIndex::VERBOSE = VERBOSE;

    if (VERBOSE)
      cerr << "[loading abismal index]" << endl;
    AbismalIndex abismal_index;
    const double start_time = omp_get_wtime();
    abismal_index.read(index_file);
    const double end_time = omp_get_wtime();
    if (VERBOSE)
      cerr << "[loading time: " << (end_time - start_time) << "]" << endl;

    if (VERBOSE)
      cerr << "[mapping "
           << (paired_end ? "paired" : "single")
           << " end: " << reads_file << "]" << endl;

    if (VERBOSE)
      cerr << "[output file: " << outfile << "]" << endl;

    // avoiding opening the stats output file until mapping is done
    se_map_stats se_stats;
    pe_map_stats pe_stats(pe_map_result_merged::min_dist,
                          pe_map_result_merged::max_dist);

    std::ofstream out(outfile);
    if (!out)
      throw runtime_error("failed to open output file: " + outfile);

    if (reads_file2.empty()) {
      if (GA_conversion || pbat_mode)
        map_single_ended<comp_ga, comp_ct>(VERBOSE, reads_file, batch_size,
                                           max_candidates,
                                           abismal_index, allow_ambig,
                                           true, se_stats, out);
      else
        map_single_ended<comp_ct, comp_ga>(VERBOSE, reads_file, batch_size,
                                           max_candidates,
                                           abismal_index, allow_ambig,
                                           false, se_stats, out);
    }
    else {
      if (pbat_mode)
        map_paired_ended<comp_ga, comp_ct>(VERBOSE, reads_file, reads_file2,
                                           batch_size,
                                           max_candidates, abismal_index,
                                           allow_ambig, true, pe_stats, out);
      else
        map_paired_ended<comp_ct, comp_ga>(VERBOSE, reads_file, reads_file2,
                                           batch_size,
                                           max_candidates, abismal_index,
                                           allow_ambig, false, pe_stats, out);
    }

    std::ofstream stat_out(outfile + ".mapstats");
    stat_out << (reads_file2.empty() ?
                 se_stats.tostring() : pe_stats.tostring());
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
