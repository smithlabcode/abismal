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

using std::max;
using std::min;


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

  void load_reads(vector<string> &names, vector<string> &reads) {
    static const size_t reserve_size = 250;

    reads.clear();
    names.clear();

    size_t line_count = 0;
    const size_t num_lines_to_read = 4*batch_size;
    string line;
    line.reserve(reserve_size);
    while (line_count < num_lines_to_read && bool(getline(*in, line))) {
      if (!(line_count & 3))
        names.push_back(line.substr(1, line.find_first_of(" \t")));
      else if ((line_count & 3) == 1) {
        if (count_if(begin(line), end(line),
                     [](const char c) {return c != 'N';}) < min_length)
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

typedef int16_t scr_t;

// type returned from atomic (e.g. nucleotide) comparisons
typedef bool cmp_t;

struct se_result {
  uint32_t pos;
  uint16_t diffs;
  char strand;
  bool ambig;

  se_result() : pos(0), diffs(max_diffs + 1), strand('+'), ambig(false) {}

  se_result(const uint32_t p, const uint16_t d, const char s) :
    pos(p), diffs(d), strand(s), ambig(false) {}

  bool operator==(const se_result &rhs) const {
   return diffs == rhs.diffs && pos == rhs.pos;
  }

  bool rc() const {return strand == '-';}
  bool operator<(const se_result &rhs) const {return diffs < rhs.diffs;}
  void reset() {diffs = max_diffs + 1; ambig = false;}
  bool valid_hit() const {return diffs < invalid_hit_diffs;}
  bool valid() const {return diffs <= max_diffs;}

  void update(const uint32_t p, const scr_t d, const char s) {
    if (d < diffs) {
      pos = p;
      diffs = d;
      strand = s;
      ambig = false;
    }
    else if (d == diffs && pos != p)
      ambig = true;
  }
  bool sure_ambig(uint32_t seed_number = 0) const {
    return ambig && (diffs == 0 || (diffs == 1 && seed_number > 0));
  }
  uint32_t get_cutoff() const {return diffs;}

  static char global_strand;
  static uint16_t max_diffs;
  static uint16_t invalid_hit_diffs;
};

uint16_t se_result::invalid_hit_diffs = 50;
uint16_t se_result::max_diffs = 6;

inline char
flip_strand(const char s) {return s=='+' ? '-' : '+';}

void
format_se(const se_result &res, const ChromLookup &cl,
          string &read, const string &read_name,
          const bool rc, ofstream &out) {
  uint32_t offset = 0, chrom_idx = 0;
  if (cl.get_chrom_idx_and_offset(res.pos, read.length(),
                                  chrom_idx, offset)) {

    if (rc) // when doing SE, this is only for GA conversion
      revcomp_inplace(read);
    out << cl.names[chrom_idx] << '\t'
        << offset << '\t'
        << offset + read.length() << '\t'
        << read_name << '\t'
        << res.diffs << '\t'
        << (rc ? flip_strand(res.strand) : res.strand) << '\t'
        << read << '\n';
  }
}

struct pe_result { // assert(sizeof(pe_result) == 16);
  se_result r1;
  se_result r2;
  pe_result() {}
  pe_result(const uint32_t p1, const uint32_t p2,
            const scr_t d1, const scr_t d2, const char s) :
    r1(se_result(p1, d1, s)), r2(se_result(p2, d2, s)) {}
  pe_result(const se_result &a, const se_result &b) : r1(a), r2(b) {}
  bool rc() const {return r1.strand == '-';}
  scr_t diffs() const {return r1.diffs + r2.diffs;}
  bool valid() const {return r1.diffs + r2.diffs < 2*se_result::max_diffs;}
  bool valid_hit() const {
    return r1.diffs + r2.diffs < se_result::invalid_hit_diffs;}
  size_t frag_len(const size_t r2_len) const {return r2.pos + r2_len - r1.pos;}
  void reset() {r1.reset(); r2.reset();}
  bool ambig() const {return r1.ambig;}
  void set_ambig() {r1.ambig = true;}
  char strand() const {return r1.strand;}
  bool is_better_than(const pe_result &rhs, const size_t r2_len) const {
    return diffs() < rhs.diffs() ||
                     (diffs() == rhs.diffs() &&
                      frag_len(r2_len) < rhs.frag_len(r2_len));
  }
  bool is_equal_to(const pe_result &rhs, const size_t r2_len) const {
    const scr_t d = diffs();
    const size_t f = frag_len(r2_len);
    return d == rhs.diffs() && f > 0 && f == rhs.frag_len(r2_len);
  }

  static uint32_t min_dist;
  static uint32_t max_dist;
};
uint32_t pe_result::min_dist = 32;
uint32_t pe_result::max_dist = 3000;

void
format_pe(const pe_result &res, const ChromLookup &cl,
          string &read1, string &read2,
          const string &name1, const string &name2, ofstream &out) {

  const size_t frag_len = res.frag_len(read2.length());
  uint32_t offset = 0, chrom_idx = 0;
  if (!cl.get_chrom_idx_and_offset(res.r1.pos, frag_len, chrom_idx, offset))
    return;
  string read;
  if (frag_len < read1.length()) {
    read.swap(read1);
    read.resize(frag_len);
  }
  else {
    revcomp_inplace(read2);
    if (frag_len < read2.length()) {
      read.swap(read2);
      read.resize(frag_len);
    }
    else {
      read.resize(frag_len, 'N');
      copy(begin(read1), end(read1), begin(read));
      copy(begin(read2), end(read2), end(read) - read2.length());
    }
  }
  out << cl.names[chrom_idx] << '\t'
      << offset << '\t'
      << offset + read.length() << '\t'
      << "FRAG:" + name1 << '\t'
      << res.diffs() << '\t'
      << res.strand() << '\t'
      << read << '\n';
}

struct pe_candidates {
  pe_candidates() : v(vector<se_result>(max_size)), sz(0) {}
  bool full() const {return sz == max_size;}
  void reset() {v.front().reset(); sz = 0;}
  uint16_t get_cutoff() const {return v.front().diffs;}
  void update(const uint32_t p, const scr_t d, const char s) {
    if (full()) {
      if (d < v.front().diffs) {
        std::pop_heap(begin(v), end(v));
        v.back() = se_result(p, d, s);
        std::push_heap(begin(v), end(v));
      }
    }
    else if (d < se_result::max_diffs) {
      v[sz++] = se_result(p, d, s);
      std::push_heap(begin(v), begin(v) + sz);
    }
  }
  bool sure_ambig(uint32_t seed_number = 0) const {
    return full() && (v[0].diffs == 0 || (v[0].diffs == 1 && seed_number > 0));
  }
  void prepare_for_mating() {
    sort(begin(v), begin(v) + sz, // no sort_heap here as heapify used "diffs"
         [](const se_result &a, const se_result &b){return a.pos < b.pos;});
    sz = unique(begin(v), begin(v) + sz) - begin(v);
  }

  vector<se_result> v;
  uint32_t sz;
  static uint32_t max_size;
};
uint32_t pe_candidates::max_size = 20;

template <const char strand>
static void
best_single(const pe_candidates &pres, se_result &res) {
  auto lim(begin(pres.v) + pres.sz);
  for (auto i(begin(pres.v)); i != lim; ++i) {
    if (i->diffs < res.diffs)
      res = se_result(i->pos, i->diffs, strand);
    else if (i->diffs == res.diffs && i->pos != res.pos)
      res.ambig = true;
  }
}

// ADS: currently read1 not used, but likely will be soon
template <const char strand>
static void
best_pair(const pe_candidates &res1, const pe_candidates &res2,
          const string &read1, const string &read2, pe_result &best) {

  const size_t r2_len = read2.length();
  auto j1 = begin(res1.v);
  const auto j1_end = j1 + res1.sz;
  const auto j2_end = begin(res2.v) + res2.sz;

  for (auto j2(begin(res2.v)); j2 != j2_end; ++j2) {

    const uint32_t lim = j2->pos + read2.length();
    while (j1 != j1_end && j1->pos + pe_result::max_dist < lim) ++j1;

    while (j1 != j1_end && j1->pos + pe_result::min_dist <= lim) {
      const pe_result p(*j1, *j2);
      if (p.is_better_than(best, r2_len))
        best = p;
      else if (p.is_equal_to(best, r2_len))
        best.set_ambig();
      ++j1;
    }
  }
}

inline double pct(const double a, const double b) {return 100.0*a/b;}
struct se_map_stats {
  se_map_stats() :
    tot_rds(0), uniq_rds(0), ambig_rds(0), unmapped_rds(0), skipped_rds(0) {}
  uint32_t tot_rds;
  uint32_t uniq_rds;
  uint32_t ambig_rds;
  uint32_t unmapped_rds;
  uint32_t skipped_rds;

  void update(const string &read, const se_result &res) {
    ++tot_rds;
    if (res.valid()) {
      if (!res.ambig) ++uniq_rds;
      else ++ambig_rds;
    }
    else ++unmapped_rds;
    skipped_rds += (read.length() == 0);
  }

  string tostring(const size_t n_tabs = 0) const {
    static const string tab = "    ";
    string t;
    for (size_t i = 0; i < n_tabs; ++i) t += tab;
    std::ostringstream oss;
    oss << t     << "total_reads: " << tot_rds << endl
        << t     << "mapped: " << endl
        << t+tab << "percent_mapped: " << pct(uniq_rds+ambig_rds, tot_rds) << endl
        << t+tab << "unique: " << uniq_rds << endl
        << t+tab << "percent_unique: " << pct(uniq_rds, tot_rds) << endl
        << t+tab << "ambiguous: " << ambig_rds << endl
        << t     << "unmapped: " << unmapped_rds << endl
        << t     << "skipped: " << skipped_rds << endl;
    return oss.str();
  }
};

struct pe_map_stats {
  pe_map_stats(const uint32_t min_d, const uint32_t max_d) :
    tot_pairs(0), uniq_pairs(0), ambig_pairs(0), unmapped_pairs(0),
    min_dist(min_d) {}
  uint32_t tot_pairs;
  uint32_t uniq_pairs;
  uint32_t ambig_pairs;
  uint32_t unmapped_pairs;
  uint32_t min_dist;
  se_map_stats end1_stats;
  se_map_stats end2_stats;

  void update_pair(const pe_result &res) {
    ++tot_pairs;
    if (res.valid()) {
      ambig_pairs += res.ambig();
      uniq_pairs += !res.ambig();
    }
    else ++unmapped_pairs;
  }

  string tostring() const {
    std::ostringstream oss;
    static const string t = "    ";
    oss << "pairs:" << endl
        << t   << "total_read_pairs: " << tot_pairs << endl
        << t   << "mapped:" << endl
        << t+t << "percent_mapped: "<< pct(uniq_pairs + ambig_pairs, tot_pairs) << endl
        << t+t << "unique: " << uniq_pairs << endl
        << t+t << "percent_unique: " << pct(uniq_pairs, tot_pairs) << endl
        << t+t << "ambiguous: " << ambig_pairs << endl
        << t   << "unmapped: " << unmapped_pairs << endl
        << "mate1:" << endl << end1_stats.tostring(1)
        << "mate2:" << endl << end2_stats.tostring(1);
    return oss.str();
  }
};

static void
update_pe_stats(const pe_result &best,
                const se_result &pos1, const se_result &neg1,
                const se_result &pos2, const se_result &neg2,
                const string &read1, const string &read2,
                pe_map_stats &pe_stats) {
  pe_stats.update_pair(best);
  if (!best.valid()) {
    pe_stats.end1_stats.update(read1, pos1 < neg1 ? pos1 : neg1);
    pe_stats.end2_stats.update(read2, pos2 < neg2 ? pos2 : neg2);
  }
}


static void
select_output(const ChromLookup &cl,
              const pe_result &best,
              const se_result &pos1, const se_result &neg1,
              const se_result &pos2, const se_result &neg2,
              string &read1, const string &name1,
              string &read2, const string &name2,
              const bool rc, ofstream &out) {
  if (best.valid()) {
    if (!best.ambig())
      format_pe(best, cl, read1, read2, name1, name2, out);
  }
  else {
    const se_result *best_mate1 = pos1 < neg1 ? &pos1 : &neg1;
    if (best_mate1->valid() && !best_mate1->ambig)
      format_se(*best_mate1, cl, read1, name1, rc, out);
    const se_result *best_mate2 = pos2 < neg2 ? &pos2 : &neg2;
    if (best_mate2->valid() && !best_mate2->ambig)
      format_se(*best_mate2, cl, read2, name2, !rc, out);
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

template <cmp_t (*F)(const char, const char), const char strand, class T>
void
check_hits(vector<uint32_t>::const_iterator start_idx,
           const vector<uint32_t>::const_iterator end_idx,
           const string::const_iterator read_start,
           const string::const_iterator read_end,
           const Genome::const_iterator genome_st, T &res) {

  for (; start_idx < end_idx && !res.sure_ambig(0); ++start_idx) {
    const scr_t diffs = full_compare<F>(res.get_cutoff(), read_start,
                                        read_end, genome_st + *start_idx);
    res.update(*start_idx, diffs, strand);
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

template <cmp_t (*F)(const char, const char), const char strand, class T>
void
process_seeds(const uint32_t genome_size,
              const uint32_t max_candidates,
              const AbismalIndex &abismal_index,
              const Genome::const_iterator genome_st,
              const string &read,
              vector<uint32_t> &hits, T &res) {

  const uint32_t readlen = read.length();
  const uint32_t genome_lim = genome_size - readlen;

  const auto read_start(begin(read));
  const auto read_end(end(read));
  const auto index_st(begin(abismal_index.index));
  const auto counter_st(begin(abismal_index.counter));

  const size_t n_starts =
    readlen > seed::n_solid_positions ? readlen - seed::n_solid_positions : 0;
  const size_t shift = std::max(1ul, n_starts/(seed::n_shifts - 1));

  for (uint32_t i = 0; i <= n_starts && !res.sure_ambig(i); i += shift) {
    hits.clear(); // hits.capacity() == max_candidates;

    uint32_t k_low, k_high;
    get_1bit_hash_low_high(read_start + i, readlen - i, k_low, k_high);

    for (uint32_t k = k_low; k < k_high; ++k) {
      auto s_idx(index_st + *(counter_st + k));
      auto e_idx(index_st + *(counter_st + k+1));

      if (s_idx < e_idx) {
        find_candidates(read_start + i, genome_st, readlen - i, s_idx, e_idx);

        if (e_idx - s_idx <= max_candidates)
          for (; s_idx != e_idx && hits.size() < max_candidates; ++s_idx)
            if (*s_idx > i && *s_idx <= genome_lim)
              hits.push_back(*s_idx - i);
      }
    }
    check_hits<F, strand>(begin(hits), end(hits),
                          read_start, read_end, genome_st, res);
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
                 const bool rc,
                 se_map_stats &se_stats,
                 ofstream &out) {

  const uint32_t genome_size = abismal_index.genome.size();
  const Genome::const_iterator genome_st(begin(abismal_index.genome));

  vector<string> names, reads;
  reads.reserve(batch_size);
  names.reserve(batch_size);

  vector<se_result> res(batch_size);

  ReadLoader rl(reads_file, batch_size);

  ProgressBar progress(get_filesize(reads_file), "mapping reads");
  if (VERBOSE)
    progress.report(cerr, 0);

  double total_mapping_time = 0;
  while (rl.good()) {

    if (VERBOSE && progress.time_to_report(rl.get_current_byte()))
      progress.report(cerr, rl.get_current_byte());

    rl.load_reads(names, reads);
    const size_t n_reads = reads.size();

#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i) res[i].reset();

    const double start_time = omp_get_wtime();
#pragma omp parallel
    {
      vector<uint32_t> hits;
      hits.reserve(max_candidates);
#pragma omp for
      for (size_t i = 0; i < reads.size(); ++i)
        if (!reads[i].empty())
          process_seeds<pos_cmp, '+'>(genome_size, max_candidates,
                                      abismal_index, genome_st,
                                      reads[i], hits, res[i]);
#pragma omp for
      for (size_t i = 0; i < reads.size(); ++i)
        if (!reads[i].empty()) {
          const string read_rc(revcomp(reads[i]));
          process_seeds<neg_cmp, '-'>(genome_size, max_candidates,
                                      abismal_index, genome_st,
                                      read_rc, hits, res[i]);
        }
    }
    total_mapping_time += (omp_get_wtime() - start_time);

    for (size_t i = 0 ; i < n_reads; ++i) {
      se_stats.update(reads[i], res[i]);
      if (res[i].valid())
        format_se(res[i], abismal_index.cl, reads[i], names[i], rc, out);
    }
  }
  if (VERBOSE) {
    progress.report(cerr, get_filesize(reads_file));
    cerr << "[total mapping time: " << total_mapping_time << endl;
  }
}


template <cmp_t (*cmp)(const char, const char), const char strand, class T>
void
map_pe_batch(const vector<string> &reads1,
             const vector<string> &reads2,
             const uint32_t max_candidates,
             const AbismalIndex &abismal_index,
             vector<T> &res1, vector<T> &res2) {

  const uint32_t genome_size = abismal_index.genome.size();
  const Genome::const_iterator genome_st(begin(abismal_index.genome));

#pragma omp parallel
  {
    vector<uint32_t> hits;
    hits.reserve(max_candidates);
#pragma omp for
    for (size_t i = 0; i < reads1.size(); ++i)
      if (!reads1[i].empty())
        process_seeds<cmp, strand>(genome_size, max_candidates, abismal_index,
                                   genome_st, reads1[i], hits, res1[i]);
#pragma omp for
    for (size_t i = 0; i < reads2.size(); ++i)
      if (!reads2[i].empty()) {
        const string read_rc(revcomp(reads2[i]));
        process_seeds<cmp, strand>(genome_size, max_candidates, abismal_index,
                                   genome_st, read_rc, hits, res2[i]);
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
                 const bool rc,
                 pe_map_stats &pe_stats,
                 ofstream &out) {

  double total_mapping_time = 0;

  ReadLoader rl1(reads_file1, batch_size);
  ReadLoader rl2(reads_file2, batch_size);

  vector<string> names1, reads1, names2, reads2;
  reads1.reserve(batch_size);
  names1.reserve(batch_size);
  reads2.reserve(batch_size);
  names2.reserve(batch_size);

  vector<pe_candidates> res1(batch_size), res2(batch_size);
  vector<pe_result> bests(batch_size);
  vector<se_result>
    res_se_p1(batch_size), res_se_p2(batch_size),
    res_se_n1(batch_size), res_se_n2(batch_size);

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
    double start_time = omp_get_wtime();
    map_pe_batch<pos_cmp, '+'>(reads1, reads2,
                               max_candidates, abismal_index, res1, res2);
    total_mapping_time += (omp_get_wtime() - start_time);

#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i) {
      res1[i].prepare_for_mating();
      res2[i].prepare_for_mating();
      best_pair<'+'>(res1[i], res2[i], reads1[i], reads2[i], bests[i]);
      if (!bests[i].valid()) {
        best_single<'+'>(res1[i], res_se_p1[i]);
        best_single<'-'>(res2[i], res_se_p2[i]);
      }
    }

#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i) {
      res1[i].reset();
      res2[i].reset();
      res_se_n1[i].reset();
      res_se_n2[i].reset();
    }
    start_time = omp_get_wtime();
    map_pe_batch<neg_cmp, '-'>(reads2, reads1,
                               max_candidates, abismal_index, res2, res1);
    total_mapping_time += (omp_get_wtime() - start_time);

#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i) {
      res1[i].prepare_for_mating();
      res2[i].prepare_for_mating();
      best_pair<'-'>(res2[i], res1[i], reads2[i], reads1[i], bests[i]);
      if (!bests[i].valid()) {
        best_single<'-'>(res1[i], res_se_n1[i]);
        best_single<'+'>(res2[i], res_se_n2[i]);
      }
    }
    for (size_t i = 0 ; i < n_reads; ++i)
      update_pe_stats(bests[i], res_se_p1[i], res_se_n1[i], res_se_p2[i],
                      res_se_n2[i], reads1[i], reads2[i], pe_stats);

    for (size_t i = 0 ; i < n_reads; ++i)
      select_output(abismal_index.cl, bests[i],
                    res_se_p1[i], res_se_n1[i], res_se_p2[i], res_se_n2[i],
                    reads1[i], names1[i], reads2[i], names2[i], rc, out);
  }
  if (VERBOSE) {
    progress.report(cerr, get_filesize(reads_file1));
    cerr << "[total mapping time: " << total_mapping_time << endl;
  }
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
                      false, pe_result::min_dist);
    opt_parse.add_opt("max-frag", 'L', "max fragment size (pe mode)",
                      false, pe_result::max_dist);
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

    se_result::max_diffs = max_diffs;

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
    pe_map_stats pe_stats(pe_result::min_dist, pe_result::max_dist);

    std::ofstream out(outfile);
    if (!out)
      throw runtime_error("failed to open output file: " + outfile);

    if (reads_file2.empty()) {
      if (GA_conversion || pbat_mode)
        map_single_ended<comp_ga, comp_ct>(VERBOSE, reads_file, batch_size,
                                           max_candidates,
                                           abismal_index, true, se_stats, out);
      else
        map_single_ended<comp_ct, comp_ga>(VERBOSE, reads_file, batch_size,
                                           max_candidates,
                                           abismal_index, false, se_stats, out);
    }
    else {
      if (pbat_mode)
        map_paired_ended<comp_ga, comp_ct>(VERBOSE, reads_file, reads_file2,
                                           batch_size,
                                           max_candidates, abismal_index,
                                           true, pe_stats, out);
      else map_paired_ended<comp_ct, comp_ga>(VERBOSE, reads_file, reads_file2,
                                              batch_size,
                                              max_candidates, abismal_index,
                                              false, pe_stats, out);
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
