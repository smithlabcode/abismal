/* Copyright (C) 2018-2020 Andrew D. Smith
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

#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <deque>

#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"
#include "OptionParser.hpp"
#include "zlib_wrapper.hpp"
#include "sam_record.hpp"
#include "bisulfite_utils.hpp"

#include "dna_four_bit_bisulfite.hpp"
#include "AbismalIndex.hpp"
#include "AbismalAlign.hpp"

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
using std::max;
using std::min;
using std::to_string;
using std::begin;
using std::end;
using std::deque;

typedef uint16_t flags_t; // every bit is a flag
typedef int16_t score_t; // aln score, edit distance, hamming distance
typedef vector<uint8_t> Read; //4-bit encoding of reads

enum conversion_type { t_rich = false, a_rich = true };

constexpr conversion_type
flip_conv(const conversion_type conv) {
  return conv == t_rich ? a_rich : t_rich;
}

constexpr flags_t
get_strand_code(const char strand, const conversion_type conv) {
  return (((strand == '-')  ? samflags::read_rc : 0) |
          ((conv == a_rich) ? bsflags::read_is_a_rich: 0));
}

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
  size_t get_current_byte() const {return gzoffset(in->fileobj);}
  void load_reads(vector<string> &names, vector<string> &reads) {
    static const size_t reserve_size = 250;

    reads.clear();
    names.clear();

    size_t line_count = 0;
    const size_t num_lines_to_read = 4*batch_size;
    string line;
    line.reserve(reserve_size);
    while (line_count < num_lines_to_read && bool(getline(*in, line))) {
      if (line_count % 4 == 0) {
        names.push_back(line.substr(1, line.find_first_of(" \t") - 1));
      }
      else if (line_count % 4 == 1) {
        if (count_if(begin(line), end(line), 
                     [](const char c) {return c != 'N';}) < min_read_length)
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

  static uint32_t min_read_length;
};

uint32_t ReadLoader::min_read_length = seed::n_seed_positions;

inline void
update_max_read_length(size_t &max_length, const vector<string> &reads) {
  for (auto it (begin(reads)); it != end(reads); ++it)
    max_length = std::max(max_length, it->size());
}

struct se_element { //assert(sizeof(se_element) == 8)
  uint32_t pos;
  score_t diffs;
  flags_t flags;

  se_element() :
    pos(0), diffs(numeric_limits<score_t>::max() - 1),  flags(0) {}

  se_element(const uint32_t p, const score_t d,
             const flags_t f) :
    pos(p), diffs(d), flags(f) {}

  bool operator==(const se_element &rhs) const {
    return diffs == rhs.diffs && pos == rhs.pos;
  }
  bool operator!=(const se_element &rhs) const {
    return diffs != rhs.diffs || pos != rhs.pos;
  }

  // this is used to keep PE candidates sorted in the max heap
  bool operator<(const se_element &rhs) const {
    return diffs < rhs.diffs;
  }

  inline bool rc() const {return samflags::check(flags, samflags::read_rc);}
  inline bool elem_is_a_rich() const {
    return samflags::check(flags, bsflags::read_is_a_rich);
  }
  inline score_t invalid_hit_threshold(const uint32_t readlen) const {
    return static_cast<score_t>(readlen * invalid_hit_frac);
  }
  void reset(const uint32_t readlen) {
    diffs = invalid_hit_threshold(readlen);
    pos = 0;
  }
  static double valid_frac;
  static double invalid_hit_frac;
};

inline bool valid(const score_t dist, const uint32_t readlen) {
  return dist <= static_cast<score_t>(se_element::valid_frac*readlen);
}

inline bool valid_hit(const score_t dist, const uint32_t readlen) {
  return dist < static_cast<score_t>(se_element::invalid_hit_frac*readlen);
}

double se_element::valid_frac = 0.1;
double se_element::invalid_hit_frac = 0.4;

struct se_result { //assert(sizeof(se_result) == 16)
  se_result() : valid(false), ambig(false), sz(1), best(se_element()),
                v(vector<se_element>(max_size)) {}

  inline bool full() const { return sz == max_size; };

  void update_best(const uint32_t p, const score_t d, const flags_t s) {
    if (d < best.diffs) {
      best = se_element(p, d, s);
      ambig = false;
    }
    else if (d == best.diffs && (p != best.pos || s != best.flags))
      ambig = true;
  }

  void update_cand(const uint32_t p, const score_t d, const flags_t s) {
    if (d < v.front().diffs) {
      if (full()) {
        std::pop_heap(begin(v), end(v));
        v.back() = se_element(p, d, s);
        std::push_heap(begin(v), end(v));
      }
      else {
        v[sz++] = se_element(p, d, s);
        std::push_heap(begin(v), begin(v) + sz);
      }
    }
  }

  void update(const uint32_t p, const score_t d, const flags_t s) {
    update_best(p, d, s);
    update_cand(p, d, s);
  }

  inline bool sure_ambig() const {
    return ambig && best.diffs == 0;
  }
  void reset(const uint32_t readlen)  {
    valid = false;
    ambig = false;
    sz = 1;
    best.reset(readlen);
    v.front().reset(readlen);
  }

  void prepare_for_alignments() {
    sort_heap(begin(v), begin(v) + sz);
    sz = unique(begin(v), begin(v) + sz) - begin(v);
  }

  uint32_t get_cutoff() const { return v.front().diffs + 1; }

  bool valid;
  bool ambig;
  uint32_t sz;
  se_element best;
  vector<se_element> v;

  static uint32_t max_size;
};

uint32_t se_result::max_size = 100;

inline bool
chrom_and_posn(const ChromLookup &cl, const string &cig, const uint32_t p,
               uint32_t &r_p, uint32_t &r_e, uint32_t &r_chr) {
  const uint32_t ref_ops = cigar_rseq_ops(cig);
  if (!cl.get_chrom_idx_and_offset(p, ref_ops, r_chr, r_p)) return false;
  r_e = r_p + ref_ops;
  return true;
}

enum map_type { map_unmapped, map_unique, map_ambig };
static map_type
format_se(const bool allow_ambig, const se_result &res, const ChromLookup &cl,
          const string &read, const string &read_name, const string &cigar,
          ostream &out) {
  if (!allow_ambig && res.ambig)
    return map_ambig;

  const se_element s = res.best;

  uint32_t ref_s = 0, ref_e = 0, chrom_idx = 0;
  if (!res.valid || !chrom_and_posn(cl, cigar, s.pos, ref_s, ref_e, chrom_idx))
    return map_unmapped;

  sam_rec sr(read_name, 0, cl.names[chrom_idx], ref_s + 1,
             255, cigar, "*", 0, 0, read, "*");
  if (s.rc())
    set_flag(sr, samflags::read_rc);

  if (allow_ambig && res.ambig)
    set_flag(sr, samflags::secondary_aln);

  sr.add_tag("NM:i:" + to_string(s.diffs));
  sr.add_tag(s.elem_is_a_rich() ? "CV:A:A" : "CV:A:T");

  out << sr << "\n";
  return res.ambig ? map_ambig : map_unique;
}

struct pe_result {
  pe_result() : valid(false), ambig(false), aln_score(0), r1(se_element()), r2(se_element()) {}
  pe_result(const se_element s1, const se_element s2) :
    valid(false), ambig(false), aln_score(0), r1(s1), r2(s2) {}

  score_t diffs() const { return r1.diffs + r2.diffs; }
  void reset(const uint32_t readlen1, const uint32_t readlen2) {
    aln_score = 0;
    valid = false;
    ambig = false;
    r1.reset(readlen1);
    r2.reset(readlen2);
  }

  bool valid;
  bool ambig;
  score_t aln_score;
  se_element r1;
  se_element r2;

  static uint32_t min_dist;
  static uint32_t max_dist;
};
uint32_t pe_result::min_dist = 32;
uint32_t pe_result::max_dist = 3000;

/* The results passed into format_pe should be on opposite strands
 * already, as those are the only valid pairings. They also should
 * have opposite "richness" for the same reason.
 *
 * On output, each read sequence is exactly as it appears in the input
 * FASTQ files. The strand is opposite for each end, and the richness
 * as well. Positions are incremented, since the SAM format is
 * 1-based. CIGAR strings are written just as they were constructed:
 * always starting with the first position on the reference,
 * regardless of the strand indicated among the flags.
 *
 * Among optional tags, we include "CV" as conversion, and it is
 * Alphanumeric with value 'A' or 'T' to show whether the C->T
 * conversion was used or the G->A (for PBAT or 2nd end of PE reads).
 */
static map_type
format_pe(const bool allow_ambig,
          const pe_result &p, const ChromLookup &cl,
          const string &read1, const string &read2,
          const string &name1, const string &name2,
          const string &cig1,  const string &cig2,
          ostream &out) {
  if (!allow_ambig && p.ambig)
    return map_ambig;

  uint32_t r_s1 = 0, r_e1 = 0, chr1 = 0; // positions in chroms (0-based)
  uint32_t r_s2 = 0, r_e2 = 0, chr2 = 0;

  // PE chromosomes differ or couldn't be found, treat read as unmapped
  if (!p.valid || !chrom_and_posn(cl, cig1, p.r1.pos, r_s1, r_e1, chr1) ||
      !chrom_and_posn(cl, cig2, p.r2.pos, r_s2, r_e2, chr2) || chr1 != chr2)
    return map_unmapped;

  const bool rc = p.r1.rc();

  // ADS: will this always evaluate correctly with unsigned
  // intermediate vals?
  const int tlen = rc ?
    (static_cast<int>(r_s1) - static_cast<int>(r_e2)) :
    (static_cast<int>(r_e2) - static_cast<int>(r_s1));

  // ADS: +1 to POS & PNEXT; "=" for RNEXT; "*" for QUAL
  sam_rec sr1(name1, 0, cl.names[chr1], r_s1 + 1, 255,
              cig1, "=", r_s2 + 1, tlen, read1, "*");

  sr1.add_tag("NM:i:" + to_string(p.r1.diffs));
  sr1.add_tag(p.r1.elem_is_a_rich() ? "CV:A:A" : "CV:A:T");
  set_flag(sr1, samflags::read_paired);
  set_flag(sr1, samflags::read_pair_mapped);
  set_flag(sr1, samflags::template_first);

  sam_rec sr2(name2, 0, cl.names[chr2], r_s2 + 1, 255,
              cig2, "=", r_s1 + 1, -tlen, read2, "*");

  sr2.add_tag("NM:i:" + to_string(p.r2.diffs));
  sr2.add_tag(p.r2.elem_is_a_rich() ? "CV:A:A" : "CV:A:T");
  set_flag(sr2, samflags::read_paired);
  set_flag(sr2, samflags::read_pair_mapped);
  set_flag(sr2, samflags::template_last);

  // second mate is reverse strand and richness of 1st mate
  if (rc) {
    set_flag(sr1, samflags::read_rc);
    set_flag(sr2, samflags::mate_rc);
  }
  else {
    set_flag(sr1, samflags::mate_rc);
    set_flag(sr2, samflags::read_rc);
  }

  if (allow_ambig && p.ambig) {
    set_flag(sr1, samflags::secondary_aln);
    set_flag(sr2, samflags::secondary_aln);
  }

  out << sr1.tostring() << "\n" << sr2.tostring() << "\n";

  return p.ambig ? map_ambig : map_unique;
}

struct pe_candidates {
  pe_candidates() : v(vector<se_element>(max_size)), sz(1) {}
  bool full() const {return sz == max_size;}
  void reset(const uint32_t readlen) {
    v.front().reset(readlen); sz = 1;
  }
  score_t get_cutoff() const {return v.front().diffs + 1;}
  void update(const uint32_t p, const score_t d, const flags_t s) {
    if (d < v.front().diffs) {
      if (full()) {
        std::pop_heap(begin(v), end(v));
        v.back() = se_element(p, d, s);
        std::push_heap(begin(v), end(v));
      }
      else {
        v[sz++] = se_element(p, d, s);
        std::push_heap(begin(v), begin(v) + sz);
      }
    }
  }
  bool sure_ambig() const {
    return full() && (v.front().diffs == 0);
  }
  void prepare_for_mating() {
    sort(begin(v), begin(v) + sz, // no sort_heap here as heapify used "diffs"
         [](const se_element &a, const se_element &b){return a.pos < b.pos;});
    sz = unique(begin(v), begin(v) + sz) - begin(v);
  }
  vector<se_element> v;
  uint32_t sz;
  static uint32_t max_size;
};

uint32_t pe_candidates::max_size = 100;

inline double pct(const double a, const double b) {return 100.0*a/b;}

struct se_map_stats {
  se_map_stats() :
    tot_rds(0), uniq_rds(0), ambig_rds(0), unmapped_rds(0), skipped_rds(0) {}
  uint32_t tot_rds;
  uint32_t uniq_rds;
  uint32_t ambig_rds;
  uint32_t unmapped_rds;
  uint32_t skipped_rds;

  void update(const se_result &s, const bool skipped) {
    ++tot_rds;
    uniq_rds += (s.valid && !s.ambig);
    ambig_rds += (s.valid && s.ambig);
    unmapped_rds += !s.valid;
    skipped_rds += skipped;
  }

  string tostring(const size_t n_tabs = 0) const {
    static const string tab = "    ";
    string t;
    for (size_t i = 0; i < n_tabs; ++i) t += tab;
    std::ostringstream oss;

    oss << t     << "total_reads: " << tot_rds << endl
        << t     << "mapped: " << endl
        << t+tab << "num_mapped: " << uniq_rds+ambig_rds << endl
        << t+tab << "percent_mapped: "
        << pct(uniq_rds+ambig_rds, tot_rds == 0 ? 1 : tot_rds) << endl
        << t+tab << "num_unique: " << uniq_rds << endl
        << t+tab << "percent_unique: "
        << pct(uniq_rds, tot_rds == 0 ? 1 : tot_rds) << endl
        << t+tab << "ambiguous: " << ambig_rds << endl
        << t     << "unmapped: " << unmapped_rds << endl
        << t     << "skipped: " << skipped_rds << endl;
    return oss.str();
  }
};

struct pe_map_stats {
  pe_map_stats() :
    tot_pairs(0), uniq_pairs(0), ambig_pairs(0), unmapped_pairs(0) {}
  uint32_t tot_pairs;
  uint32_t uniq_pairs;
  uint32_t ambig_pairs;
  uint32_t unmapped_pairs;
  uint32_t min_dist;
  se_map_stats end1_stats;
  se_map_stats end2_stats;

  void update(const bool allow_ambig, const pe_result &p,
              const se_result &s1, const bool skipped_se1,
              const se_result &s2, const bool skipped_se2) {
    ++tot_pairs;
    ambig_pairs += (p.valid && p.ambig);
    uniq_pairs += (p.valid && !p.ambig);

    if (!p.valid || (!allow_ambig && p.ambig)) {
      end1_stats.update(s1, skipped_se1);
      end2_stats.update(s2, skipped_se2);
    }
    else ++unmapped_pairs;
  }

  string tostring() const {
    std::ostringstream oss;
    static const string t = "    ";
    oss << "pairs:" << endl
        << t   << "total_read_pairs: " << tot_pairs << endl
        << t   << "mapped:" << endl
        << t+t << "num_mapped: " << uniq_pairs + ambig_pairs << endl
        << t+t << "percent_mapped: "
        << pct(uniq_pairs + ambig_pairs, tot_pairs) << endl
        << t+t << "num_unique: " << uniq_pairs << endl
        << t+t << "percent_unique: " << pct(uniq_pairs, tot_pairs) << endl
        << t+t << "ambiguous: " << ambig_pairs << endl
        << t   << "unmapped: " << unmapped_pairs << endl
        << "mate1:" << endl << end1_stats.tostring(1)
        << "mate2:" << endl << end2_stats.tostring(1);
    return oss.str();
  }
};

static void
select_output(const bool allow_ambig, const ChromLookup &cl,
              pe_result &best, se_result &se1, se_result &se2,
              const string &read1, const string &name1,
              const string &read2, const string &name2,
              const string &cig1, const string &cig2,
              ostream &out) {

  const map_type pe_map_type = format_pe(allow_ambig, best, cl,
                                         read1, read2, name1,
                                         name2, cig1, cig2, out);
  if (pe_map_type == map_unmapped ||
      (!allow_ambig && pe_map_type == map_ambig)) {
    // GS: do not report in mapstats a read that was not reported
    if (pe_map_type == map_unmapped)
      best.reset(read1.size(), read2.size());

    if (format_se(allow_ambig, se1, cl, read1, name1, cig1, out) ==
        map_unmapped)
      se1.reset(read1.size());
    if (format_se(allow_ambig, se2, cl, read2, name2, cig2, out) ==
        map_unmapped)
      se2.reset(read2.size());
  }
}

inline bool
the_comp(const char a, const char b) {
  return (a & b) == 0;
}

score_t
full_compare(const score_t cutoff,
             Read::const_iterator read_itr,
             const Read::const_iterator read_end,
             Genome::const_iterator genome_itr) {
  score_t d = 0;
  while (d < cutoff && read_itr != read_end) {
    d += mismatch_lookup[*read_itr & *genome_itr];
    ++read_itr, ++genome_itr;
  }
  return d;
}

template <const uint16_t strand_code, class result_type>
inline void
check_hits(vector<uint32_t>::const_iterator start_idx,
           const vector<uint32_t>::const_iterator end_idx,
           const Read::const_iterator even_read_st,
           const Read::const_iterator even_read_end,
           const Read::const_iterator odd_read_st,
           const Read::const_iterator odd_read_end,
           const Genome::const_iterator genome_st,
           const uint32_t offset,
           result_type &res) {
  for (; start_idx != end_idx && !res.sure_ambig(); ++start_idx) {
    // GS: adds the next candidate to cache while current is compared
    __builtin_prefetch(
        &(*(genome_st + ((*(start_idx + 1) - offset) >> 1)))
    );

    const uint32_t the_pos = (*start_idx - offset);
    const score_t diffs = (the_pos & 1) ?
       full_compare(res.get_cutoff(), odd_read_st, odd_read_end,
                    genome_st + (the_pos >> 1)):
       full_compare(res.get_cutoff(), even_read_st, even_read_end,
                    genome_st + (the_pos >> 1));
    res.update(the_pos, diffs, strand_code);
  }
}

template<const uint32_t seed_lim>
static void
get_minimizers(const uint32_t readlen,
               Read::const_iterator read_start,
               vector<kmer_loc> &kmers) {
  kmers.clear();
  deque<kmer_loc> window_kmers;

  const uint32_t shift_lim = (readlen >= seed_lim) ?
                             (readlen - seed_lim) : 0;

  size_t kmer = 0;
  get_1bit_hash<seed::key_weight>(read_start, kmer);
  add_kmer(window_kmers, kmer_loc(kmer, 0));

  read_start += seed::key_weight;

  const uint32_t min_lim = min(shift_lim, seed::minimizer_window_size);
  for (size_t i = 1; i < min_lim; ++i) {
    shift_hash_key<seed::key_weight>(*read_start++, kmer);
    add_kmer(window_kmers, kmer_loc(kmer, i));
  }

  // GS: this number is guaranteed to be different from all kmers
  size_t last_minimizer = (1ull << (seed::key_weight + 1ull));
  for (size_t i = seed::minimizer_window_size; i <= shift_lim; ++i) {
    // GS: optimization to avoid getting the same minimizer several times
    const kmer_loc best_minimizer = window_kmers.front();
    if (best_minimizer.kmer != last_minimizer) {
      last_minimizer = best_minimizer.kmer;
      kmers.push_back(best_minimizer);
    }
    shift_hash_key<seed::key_weight>(*read_start++, kmer);
    add_kmer(window_kmers, kmer_loc(kmer, i));
  }

  const kmer_loc best_minimizer = window_kmers.front();
  if (best_minimizer.kmer != last_minimizer) {
    last_minimizer = best_minimizer.kmer;
    kmers.push_back(best_minimizer);
  }
}

// ADS: probably should be a lambda function for brevity
struct compare_bases {
  compare_bases(const genome_iterator g_) : g(g_) {}
  bool operator()(const uint32_t mid, const uint32_t chr) const {
    return get_bit(*(g + mid)) < chr;
  }
  const genome_iterator g;
};

template<const uint32_t seed_lim>
static void
find_candidates(const Read::const_iterator read_start,
                const genome_iterator gi,
                const uint32_t read_lim, // not necessarily read len
                vector<uint32_t>::const_iterator &low,
                vector<uint32_t>::const_iterator &high) {
  const uint32_t lim = min(read_lim, seed_lim);
  for (uint32_t p = seed::key_weight; p != lim && low < high; ++p) {
    auto first_1 = lower_bound(low, high, 1, compare_bases(gi + p));
    if (get_bit(*(read_start + p)) == 0) {
      if (first_1 == high) return; // need 0s; whole range is 0s
      high = first_1;
    }
    else {
      if (first_1 == low) return; // need 1s; whole range is 1s
      low = first_1;
    }
  }
}


template <const uint16_t strand_code, class result_type>
void
process_seeds(const uint32_t max_candidates,
              const vector<uint32_t>::const_iterator counter_st,
              const vector<uint32_t>::const_iterator index_st,
              const genome_iterator genome_st,
              const Read &read_seed, const Read &read_even,
              const Read &read_odd,
              vector<kmer_loc> &kmers,
              result_type &res) {
  const uint32_t readlen = read_seed.size();

  // used to get positions in the genome
  const auto read_start(begin(read_seed));

  // used to compare even positions in the genome
  const auto even_read_start(begin(read_even));
  const auto even_read_end(end(read_even));

  // used to compare odd positions in the genome
  const auto odd_read_start(begin(read_odd));
  const auto odd_read_end(end(read_odd));
  size_t k = 0;
  get_1bit_hash<seed::key_weight>(read_start, k);

  // specific step: look for exact matches
  for (uint32_t j = 0; j != seed::minimizer_window_size; ++j) {
    auto s_idx(index_st + *(counter_st + k));
    auto e_idx(index_st + *(counter_st + k + 1));
    if (s_idx < e_idx) {
      find_candidates<seed::n_sorting_positions>(
        read_start + j, genome_st, readlen - j, s_idx, e_idx
      );
      if ((e_idx - s_idx) <= max_candidates){
        check_hits<strand_code>(s_idx, e_idx, even_read_start, even_read_end,
          odd_read_start, odd_read_end, genome_st.itr, j, res
        );
      }
    }
    if (res.sure_ambig()) return;
    shift_hash_key<seed::key_weight>(*(read_start + seed::key_weight + j), k);
  }

  // sensitive step: get positions using minimizers
  get_minimizers<seed::n_seed_positions>(readlen, read_start, kmers);
  const size_t num_minimizers = kmers.size();
  for (size_t i = 0; i < num_minimizers; ++i) {
    auto s_idx(index_st + *(counter_st + kmers[i].kmer));
    auto e_idx(index_st + *(counter_st + kmers[i].kmer + 1));
    if (s_idx < e_idx) {
      find_candidates<seed::n_seed_positions>(
        read_start + kmers[i].loc, genome_st, readlen - kmers[i].loc,
        s_idx, e_idx
      );
      if ((e_idx - s_idx) <= max_candidates) {
        check_hits<strand_code>(
          s_idx, e_idx, even_read_start, even_read_end, odd_read_start,
          odd_read_end, genome_st.itr, kmers[i].loc, res
        );
      }
    }
  }
}

template <const bool convert_a_to_g>
static void
prep_read(const string &r, Read &pread) {
  pread.resize(r.size());
  for (size_t i = 0; i != r.size(); ++i)
    pread[i] = (convert_a_to_g ?
                (encode_base_a_rich[static_cast<unsigned char>(r[i])]) :
                (encode_base_t_rich[static_cast<unsigned char>(r[i])]));
}

/* GS: this function encodes an ASCII character string into
 * two byte arrays, where each byte contains two bases, and
 * each base is represented in four bits, with active bits
 * representing the possible base matches to the genome. One
 * of the arrays(pread_even) will be used to compare the read
 * to even positions in the genome, while the other (pread_odd)
 * is used to compare odd bases. When the number of bases is odd,
 * or to account for the fact that pread_odd has a single base
 * in the first byte, reads can have a 1111 as the first or last
 * base. This will match with any base in the genome and will therefore
 * be disregarded for mismatch counting except if it aligns with Ns
 * in the genome (encoded as 0000 to mismatch with everything). These
 * cases should be unlikely in practice
 * */
static void
prep_for_seeds(const Read &pread_seed, Read &pread_even,
               Read &pread_odd) {
  static const uint8_t first_base_match_any = 0x0F;
  static const uint8_t last_base_match_any = 0xF0;
  const size_t sz = pread_seed.size();
  pread_even.resize((sz + 1)/2);
  pread_odd.resize(sz/2 + 1);

  // even encoding
  const auto s_idx(begin(pread_seed));
  const auto e_idx(end(pread_seed));
  size_t i = 0;
  for (auto it(s_idx); i != pread_even.size(); ++it, ++i) {
    pread_even[i] = *it;
    if (++it != e_idx) pread_even[i] |= ((*it) << 4);
    else pread_even[i] |= last_base_match_any;
  }

  // odd encoding
  pread_odd[0] = first_base_match_any |  ((*s_idx) << 4);
  i = 1;
  for (auto it(s_idx + 1); i != pread_odd.size(); ++it, ++i) {
    pread_odd[i] = *it;
    if (++it != e_idx) pread_odd[i] |= ((*it) << 4);
    else pread_odd[i] |= last_base_match_any;
  }
}

using AbismalAlignSimple =
      AbismalAlign<simple_aln::mismatch_score, simple_aln::indel>;

static score_t
align_read(se_element &res, const Read &pread,
         uint32_t &len, string &cigar, AbismalAlignSimple &aln) {
  const score_t readlen = static_cast<score_t>(pread.size());
  if (res.diffs == readlen ||
      res.diffs < simple_aln::min_diffs_to_align) {
    simple_aln::make_default_cigar(readlen, cigar);
    len = readlen;

    // the score without indels
    return simple_aln::default_score(readlen, res.diffs);
  }
  return aln.align(pread, res.pos, len, cigar);
}

static score_t
align_read(se_element &res, const string &read,
           Read &pread, uint32_t &len, string &cigar,
           AbismalAlignSimple &aln) {
  // re-encodes the read based on best match
  if (res.rc()) {
    const string read_rc(revcomp(read));
    // rc reverses richness of read
    if (res.elem_is_a_rich()) prep_read<false>(read_rc, pread);
    else prep_read<true>(read_rc, pread);
  }
  else {
    if (res.elem_is_a_rich()) prep_read<true>(read, pread);
    else prep_read<false>(read, pread);
  }
  return align_read(res, pread, len, cigar, aln);
}

static void
align_se_candidates(se_result &res, string &cigar,
                    const string &read, Read &pread,
                    AbismalAlignSimple &aln) {
  const score_t readlen = static_cast<score_t>(read.size());
  // GS: no improvement expected from alignment
  if (res.best.diffs <= simple_aln::min_diffs_to_align) {
    res.valid = valid(res.best.diffs, readlen);
    simple_aln::make_default_cigar(readlen, cigar);
    return;
  }

  // variables to store the best alignment score
  bool best_is_ambig = false;
  bool best_is_valid = false;
  uint32_t len = 0;
  uint32_t best_len = 0;
  score_t best_score = 0;
  se_element the_best;
  string cand_cigar;

  res.prepare_for_alignments();
  const auto lim(begin(res.v) + res.sz);
  for (auto it(begin(res.v)); it != lim && valid_hit(it->diffs, readlen); ++it) {
    se_element s = *it;
    const score_t scr = align_read(s, read, pread, len, cand_cigar, aln);
    const score_t diffs = simple_aln::edit_distance(scr, len, cand_cigar);
    if (valid(diffs, len)) {
      if (scr > best_score) {
        the_best = s; best_score = scr; best_len = len;
        cigar = cand_cigar;
        best_is_valid = true;
        best_is_ambig = false;
      }
      else if (scr == best_score && !(s == the_best))
        best_is_ambig = true;
    }
  }

  // now replace res.best with the_best if there is an edit d
  res.best = the_best;
  res.best.diffs = simple_aln::edit_distance(best_score, best_len, cigar);
  res.ambig = best_is_ambig;
  res.valid = best_is_valid;
  res.valid = valid(res.best.diffs, read.size());
}

template <const  conversion_type conv>
inline void
map_single_ended(const bool VERBOSE, const bool allow_ambig,
                 const size_t batch_size, const size_t max_candidates,
                 const AbismalIndex &abismal_index, ReadLoader &rl,
                 omp_lock_t &read_lock, omp_lock_t &write_lock,
                 se_map_stats &se_stats, ostream &out,
                 ProgressBar &progress) {
  const auto counter_st(begin(abismal_index.counter));
  const auto index_st(begin(abismal_index.index));
  const genome_iterator genome_st(begin(abismal_index.genome));

  vector<string> names;
  vector<string> reads;
  vector<string> cigar;

  vector<se_result> res;

  names.reserve(batch_size);
  reads.reserve(batch_size);

  cigar.resize(batch_size);
  res.resize(batch_size);

  Read pread, pread_even, pread_odd;
  vector<kmer_loc> kmers;

  while (rl.good()) {
    omp_set_lock(&read_lock);
    rl.load_reads(names, reads);
    const size_t the_byte = rl.get_current_byte();
    omp_unset_lock(&read_lock);

    size_t max_batch_read_length = 0;
    update_max_read_length(max_batch_read_length, reads);

    AbismalAlignSimple aln(genome_st, max_batch_read_length);
    kmers.reserve(max_batch_read_length);

    const size_t n_reads = reads.size();
    for (size_t i = 0; i < n_reads; ++i) {
      res[i].reset(reads[i].size());
      if (!reads[i].empty()) {
        prep_read<conv>(reads[i], pread);
        prep_for_seeds(pread, pread_even, pread_odd);
        process_seeds<get_strand_code('+', conv)>(max_candidates,
          counter_st, index_st, genome_st, pread,
          pread_even, pread_odd, kmers, res[i]
        );

        const string read_rc(revcomp(reads[i]));
        prep_read<!conv>(read_rc, pread);
        prep_for_seeds(pread, pread_even, pread_odd);
        process_seeds<get_strand_code('-', conv)>(max_candidates,
          counter_st, index_st, genome_st, pread,
          pread_even, pread_odd, kmers, res[i]
        );
        align_se_candidates(res[i], cigar[i], reads[i], pread, aln);
      }
    }

    omp_set_lock(&write_lock);
    for (size_t i = 0; i < n_reads; ++i) {
      if (format_se(allow_ambig, res[i], abismal_index.cl, reads[i],
            names[i], cigar[i], out) == map_unmapped)
        res[i].reset(reads[i].size());
      se_stats.update(res[i], reads[i].length() == 0);
    }
    omp_unset_lock(&write_lock);
    if (VERBOSE && progress.time_to_report(the_byte))
      progress.report(cerr, the_byte);
  }
}

inline void
map_single_ended_rand(const bool VERBOSE, const bool allow_ambig,
                      const size_t batch_size, const size_t max_candidates,
                      const AbismalIndex &abismal_index, ReadLoader &rl,
                      omp_lock_t &read_lock, omp_lock_t &write_lock,
                      se_map_stats &se_stats, ostream &out,
                      ProgressBar &progress) {
  const auto counter_st(begin(abismal_index.counter));
  const auto index_st(begin(abismal_index.index));
  const genome_iterator genome_st(begin(abismal_index.genome));

  vector<string> names;
  vector<string> reads;
  vector<string> cigar;
  vector<se_result> res;

  names.reserve(batch_size);
  reads.reserve(batch_size);
  cigar.resize(batch_size);
  res.resize(batch_size);

  Read pread, pread_even, pread_odd;
  vector<kmer_loc> kmers;

  while (rl.good()) {
    omp_set_lock(&read_lock);
    rl.load_reads(names, reads);
    const size_t the_byte = rl.get_current_byte();
    omp_unset_lock(&read_lock);

    size_t max_batch_read_length = 0;
    update_max_read_length(max_batch_read_length, reads);

    AbismalAlignSimple aln(genome_st, max_batch_read_length);
    kmers.reserve(max_batch_read_length);

    const size_t n_reads = reads.size();
    for (size_t i = 0; i < n_reads; ++i) {
      res[i].reset(reads[i].size());
      if (!reads[i].empty()) {
        prep_read<t_rich>(reads[i], pread);
        prep_for_seeds(pread, pread_even, pread_odd);
        process_seeds<get_strand_code('+', t_rich)>(max_candidates,
          counter_st, index_st, genome_st, pread,
          pread_even, pread_odd, kmers, res[i]
        );

        prep_read<a_rich>(reads[i], pread);
        prep_for_seeds(pread, pread_even, pread_odd);
        process_seeds<get_strand_code('+', a_rich)>(max_candidates,
          counter_st, index_st, genome_st, pread,
          pread_even, pread_odd, kmers, res[i]
        );

        const string read_rc(revcomp(reads[i]));
        prep_read<t_rich>(read_rc, pread);
        prep_for_seeds(pread, pread_even, pread_odd);
        process_seeds<get_strand_code('-', a_rich)>(max_candidates,
          counter_st, index_st, genome_st, pread, pread_even,
          pread_odd, kmers, res[i]
        );

        prep_read<a_rich>(read_rc, pread);
        prep_for_seeds(pread, pread_even, pread_odd);
        process_seeds<get_strand_code('-', t_rich)>(max_candidates,
          counter_st, index_st, genome_st, pread,
          pread_even, pread_odd, kmers, res[i]
        );

        align_se_candidates(res[i], cigar[i], reads[i], pread, aln);
      }
    }

    omp_set_lock(&write_lock);
    for (size_t i = 0; i < n_reads; ++i) {
      if (format_se(allow_ambig, res[i], abismal_index.cl, reads[i],
            names[i], cigar[i], out) == map_unmapped)
        res[i].reset(reads[i].size());
      se_stats.update(res[i], reads[i].length() == 0);
    }
    omp_unset_lock(&write_lock);
    if (VERBOSE && progress.time_to_report(the_byte))
      progress.report(cerr, the_byte);
  }
}

template <const conversion_type conv, const bool random_pbat>
void
run_single_ended(const bool VERBOSE,
                 const bool allow_ambig,
                 const string &reads_file,
                 const size_t batch_size,
                 const size_t max_candidates,
                 const AbismalIndex &abismal_index,
                 se_map_stats &se_stats,
                 ostream &out) {
  ReadLoader rl(reads_file, batch_size);
  ProgressBar progress(get_filesize(reads_file), "mapping reads");

  omp_lock_t read_lock;
  omp_lock_t write_lock;

  omp_init_lock(&read_lock);
  omp_init_lock(&write_lock);

  if (VERBOSE)
    progress.report(cerr, 0);

  double start_time = omp_get_wtime();
  if (VERBOSE && progress.time_to_report(rl.get_current_byte()))
    progress.report(cerr, rl.get_current_byte());

#pragma omp parallel for
  for (int i = 0; i < omp_get_num_threads(); ++i) {
    if (random_pbat)
      map_single_ended_rand(VERBOSE, allow_ambig, batch_size, max_candidates,
        abismal_index, rl, read_lock, write_lock, se_stats, out, progress);
    else
      map_single_ended<conv>(VERBOSE, allow_ambig, batch_size, max_candidates,
          abismal_index, rl, read_lock, write_lock, se_stats, out, progress);
  }
  if (VERBOSE) {
    //progress.report(cerr, get_filesize(reads_file));
    cerr << "-> total mapping time: " << omp_get_wtime() - start_time << ""
         << endl;
  }
}

static void
best_single(const pe_candidates &pres, se_result &res) {
  const auto lim(begin(pres.v) + pres.sz);

  // get best and second best by mismatch
  for (auto i(begin(pres.v)); i != lim && !res.sure_ambig(); ++i) {
    res.update(i->pos, i->diffs, i->flags);
  }
}

template <const bool swap_ends>
static void
best_pair(const pe_candidates &res1, const pe_candidates &res2,
          const Read &pread1, const Read &pread2,
          string &cig1, string &cig2,
          AbismalAlignSimple &aln, pe_result &best) {

  auto j1 = begin(res1.v);
  const auto j1_end = j1 + res1.sz;
  const auto j2_end = begin(res2.v) + res2.sz;

  uint32_t len1 = 0, len2 = 0;
  se_element s1, s2;
  string cand_cig1, cand_cig2;

  for (auto j2(begin(res2.v)); j2 != j2_end; ++j2) {
    s2 = *j2;
    if (valid_hit(s2.diffs, pread2.size())){
      const uint32_t unaligned_lim = s2.pos + pread2.size();
      for (j1 = begin(res1.v); j1 != j1_end &&
           j1->pos + pe_result::max_dist < unaligned_lim; ++j1);

      bool aligned_s2 = false;
      score_t scr2 = 0;
      score_t d2 = 0;
      uint32_t aligned_lim = 0;
      while (j1 != j1_end && j1->pos + pe_result::min_dist <= unaligned_lim) {
        if (valid_hit(j1->diffs, pread1.size())) {
          s1 = *j1;
          const score_t scr1 = align_read(s1, pread1, len1, cand_cig1, aln);
          const score_t d1 = simple_aln::edit_distance(scr1, len1, cand_cig1);

          // GS: guarantees that j2 is aligned only once
          if (!aligned_s2) {
            scr2 = align_read(s2, pread2, len2, cand_cig2, aln);
            d2 = simple_aln::edit_distance(scr2, len2, cand_cig2);
            aligned_lim = s2.pos + cigar_rseq_ops(cand_cig2);
            aligned_s2 = true;
          }

          // GS: only accept if length post alignment is still within limits
          if (valid(d1, len1) && valid(d2, len2) &&
              (s1.pos + pe_result::max_dist >= aligned_lim) &&
              (s1.pos + pe_result::min_dist <= aligned_lim)){
            if (scr1 + scr2 > best.aln_score) {
              best.aln_score = scr1 + scr2;
              cig1 = cand_cig1; cig2 = cand_cig2;
              best.r1 = swap_ends ? s2 : s1;
              best.r2 = swap_ends ? s1 : s2;
              best.r1.diffs = swap_ends ? d2 : d1;
              best.r2.diffs = swap_ends ? d1 : d2;
              best.valid = true;
              best.ambig = false;
            }
            else if (scr1 + scr2 == best.aln_score) {
              best.ambig = true;
            }
          }
        }
        ++j1;
      }
    }
  }
}

template <const bool swap_ends>
inline void
select_maps(const Read &pread1, const Read &pread2,
            string &cig1, string &cig2,
            pe_candidates &res1, pe_candidates &res2,
            se_result &res_se1, se_result &res_se2,
            AbismalAlignSimple &aln, pe_result &best) {
  res1.prepare_for_mating();
  res2.prepare_for_mating();

  best_pair<swap_ends>(res1, res2, pread1, pread2, cig1, cig2, aln, best);

  // GS: gets best single if best pair is reported as it later it
  // can become ambiguous and we want the best possible SE in that case
  best_single(res1, res_se1);
  best_single(res2, res_se2);
}

template <const bool cmp, const bool swap_ends,
          const uint16_t strand_code1, const uint16_t strand_code2>
inline void
map_fragments(const string &read1, const string &read2,
             Read &pread1,  Read &pread2,
             Read &pread_even, Read &pread_odd,
             string &cigar1, string &cigar2,
             const uint32_t max_candidates,
             const vector<uint32_t>::const_iterator counter_st,
             const vector<uint32_t>::const_iterator index_st,
             const genome_iterator genome_st,
             AbismalAlignSimple &aln,
             vector<kmer_loc> &kmers,
             pe_candidates &res1, pe_candidates &res2,
             se_result &res_se1, se_result &res_se2,
             pe_result &bests) {
  res1.reset(read1.size());
  res2.reset(read2.size());

  if (!read1.empty()) {
    prep_read<cmp>(read1, pread1);
    prep_for_seeds(pread1, pread_even, pread_odd);
    process_seeds<strand_code1>(max_candidates, counter_st,
      index_st, genome_st, pread1, pread_even, pread_odd, kmers,
      res1
    );
  }

  if (!read2.empty()) {
    const string read_rc(revcomp(read2));
    prep_read<cmp>(read_rc, pread2);
    prep_for_seeds(pread2, pread_even, pread_odd);
    process_seeds<strand_code2>(max_candidates, counter_st,
      index_st, genome_st, pread2, pread_even, pread_odd, kmers,
      res2
    );
  }

  select_maps<swap_ends>(pread1, pread2, cigar1, cigar2,
                         res1, res2, res_se1, res_se2, aln, bests);
}

template <const conversion_type conv>
inline void
map_paired_ended(const bool VERBOSE,
                 const bool allow_ambig, const size_t batch_size,
                 const size_t max_candidates, const AbismalIndex &abismal_index,
                 ReadLoader &rl1, ReadLoader &rl2,
                 omp_lock_t &read_lock, omp_lock_t &write_lock,
                 pe_map_stats &pe_stats, ostream &out,
                 ProgressBar &progress) {
  const auto counter_st(begin(abismal_index.counter));
  const auto index_st(begin(abismal_index.index));
  const genome_iterator genome_st(begin(abismal_index.genome));

  vector<string> names1, reads1, cigar1;
  vector<string> names2, reads2, cigar2;

  vector<pe_result> bests;
  vector<pe_candidates> res1, res2;
  vector<se_result> res_se1, res_se2;

  names1.reserve(batch_size);
  reads1.reserve(batch_size);
  cigar1.resize(batch_size);

  names2.reserve(batch_size);
  reads2.reserve(batch_size);
  cigar2.resize(batch_size);

  bests.resize(batch_size);

  res1.resize(batch_size);
  res2.resize(batch_size);

  res_se1.resize(batch_size);
  res_se2.resize(batch_size);

  vector<kmer_loc> kmers;

  Read pread1, pread2, pread_even, pread_odd;

  while (rl1.good() && rl2.good()) {
    omp_set_lock(&read_lock);

    rl1.load_reads(names1, reads1);
    rl2.load_reads(names2, reads2);
    const size_t the_byte = rl1.get_current_byte();
    omp_unset_lock(&read_lock);

    size_t max_batch_read_length = 0;
    update_max_read_length(max_batch_read_length, reads1);
    update_max_read_length(max_batch_read_length, reads2);

    AbismalAlignSimple aln(genome_st, max_batch_read_length);
    kmers.reserve(max_batch_read_length);

    const size_t n_reads = reads1.size();
    for (size_t i = 0 ; i < n_reads; ++i) {
      res_se1[i].reset(reads1[i].size());
      res_se2[i].reset(reads2[i].size());
      bests[i].reset(reads1[i].size(), reads2[i].size());

      map_fragments<conv, false,
                   get_strand_code('+',conv),
                   get_strand_code('-', flip_conv(conv))>(
         reads1[i], reads2[i], pread1, pread2, pread_even, pread_odd,
         cigar1[i], cigar2[i], max_candidates, counter_st,
         index_st, genome_st, aln, kmers, res1[i], res2[i],
         res_se1[i], res_se2[i], bests[i]
      );

      map_fragments<!conv, true,
                   get_strand_code('+', flip_conv(conv)),
                   get_strand_code('-', conv)>(
         reads2[i], reads1[i], pread1, pread2, pread_even, pread_odd,
         cigar2[i], cigar1[i], max_candidates, counter_st,
         index_st, genome_st, aln, kmers, res2[i], res1[i],
         res_se2[i], res_se1[i], bests[i]
      );

      if (!bests[i].valid || (!allow_ambig && bests[i].ambig)) {
        align_se_candidates(res_se1[i], cigar1[i], reads1[i], pread1, aln);
        align_se_candidates(res_se2[i], cigar2[i], reads2[i], pread2, aln);
      }
    }

    omp_set_lock(&write_lock);
    for (size_t i = 0; i < n_reads; ++i) {
      select_output(allow_ambig, abismal_index.cl, bests[i],
        res_se1[i], res_se2[i], reads1[i], names1[i], reads2[i], names2[i],
          cigar1[i], cigar2[i], out
        );
      pe_stats.update(allow_ambig, bests[i],
                      res_se1[i], reads1[i].length() == 0,
                      res_se2[i], reads2[i].length() == 0);
    }
    omp_unset_lock(&write_lock);
    if (VERBOSE && progress.time_to_report(the_byte))
      progress.report(cerr, the_byte);
  }
}

inline void
map_paired_ended_rand(const bool VERBOSE, const bool allow_ambig,
                      const size_t batch_size, const size_t max_candidates,
                      const AbismalIndex &abismal_index,
                      ReadLoader &rl1, ReadLoader &rl2,
                      omp_lock_t &read_lock, omp_lock_t &write_lock,
                      pe_map_stats &pe_stats, ostream &out,
                      ProgressBar &progress) {
  const auto counter_st(begin(abismal_index.counter));
  const auto index_st(begin(abismal_index.index));
  const genome_iterator genome_st(begin(abismal_index.genome));

  vector<string> names1, reads1, cigar1;
  vector<string> names2, reads2, cigar2;

  vector<pe_result> bests;
  vector<pe_candidates> res1, res2;
  vector<se_result> res_se1, res_se2;

  names1.reserve(batch_size);
  reads1.reserve(batch_size);
  cigar1.resize(batch_size);

  names2.reserve(batch_size);
  reads2.reserve(batch_size);
  cigar2.resize(batch_size);

  bests.resize(batch_size);

  res1.resize(batch_size);
  res2.resize(batch_size);

  res_se1.resize(batch_size);
  res_se2.resize(batch_size);

  vector<kmer_loc> kmers;

  Read pread1, pread2, pread_even, pread_odd;

  while (rl1.good() && rl2.good()) {
    omp_set_lock(&read_lock);
    rl1.load_reads(names1, reads1);
    rl2.load_reads(names2, reads2);
    const size_t the_byte = rl1.get_current_byte();
    omp_unset_lock(&read_lock);

    size_t max_batch_read_length = 0;
    update_max_read_length(max_batch_read_length, reads1);
    update_max_read_length(max_batch_read_length, reads2);

    AbismalAlignSimple aln(genome_st, max_batch_read_length);
    kmers.reserve(max_batch_read_length);

    const size_t n_reads = reads1.size();
    for (size_t i = 0 ; i < n_reads; ++i) {
      res_se1[i].reset(reads1[i].size());
      res_se2[i].reset(reads2[i].size());
      bests[i].reset(reads1[i].size(), reads2[i].size());

      // GS: (1) T/A-rich +/- strand
      map_fragments<t_rich, false,
                   get_strand_code('+', t_rich),
                   get_strand_code('-', a_rich)>(
         reads1[i], reads2[i], pread1, pread2, pread_even, pread_odd,
         cigar1[i], cigar2[i], max_candidates, counter_st,
         index_st, genome_st, aln, kmers, res1[i], res2[i],
         res_se1[i], res_se2[i], bests[i]
      );

      // GS: (2) T/A-rich, -/+ strand
      map_fragments<a_rich, true,
                   get_strand_code('+', a_rich),
                   get_strand_code('-', t_rich)>(
         reads2[i], reads1[i], pread2, pread1, pread_even, pread_odd,
         cigar2[i], cigar1[i], max_candidates, counter_st,
         index_st, genome_st, aln, kmers, res2[i], res1[i],
         res_se2[i], res_se1[i], bests[i]
      );

      // GS: (3) A/T-rich +/- strand
      map_fragments<a_rich, false,
                   get_strand_code('+', a_rich),
                   get_strand_code('-', t_rich)>(
         reads1[i], reads2[i], pread1, pread2, pread_even, pread_odd,
         cigar1[i], cigar2[i], max_candidates, counter_st, 
         index_st, genome_st, aln, kmers, res1[i], res2[i],
         res_se1[i], res_se2[i], bests[i]
      );

      // GS: (4) A/T-rich, -/+ strand
      map_fragments<t_rich, true,
                   get_strand_code('+', t_rich),
                   get_strand_code('-', a_rich)>(
         reads2[i], reads1[i], pread2, pread1, pread_even, pread_odd,
         cigar2[i], cigar1[i], max_candidates, counter_st,
         index_st, genome_st, aln, kmers, res2[i], res1[i],
         res_se2[i], res_se1[i], bests[i]
      );

      // GS: align best SE candidates if no concordant pairs found
      if (!bests[i].valid || bests[i].ambig) {
        align_se_candidates(res_se1[i], cigar1[i], reads1[i], pread1, aln);
        align_se_candidates(res_se2[i], cigar2[i], reads2[i], pread2, aln);
      }
    }

    omp_set_lock(&write_lock);
    for (size_t i = 0; i < n_reads; ++i) {
      select_output(allow_ambig, abismal_index.cl, bests[i],
        res_se1[i], res_se2[i], reads1[i], names1[i], reads2[i], names2[i],
        cigar1[i], cigar2[i], out
      );

      pe_stats.update(allow_ambig, bests[i],
                      res_se1[i], reads1[i].length() == 0,
                      res_se2[i], reads2[i].length() == 0);
    }
    omp_unset_lock(&write_lock);
    if (VERBOSE && progress.time_to_report(the_byte))
      progress.report(cerr, the_byte);
  }
}

template <const conversion_type conv, const bool random_pbat>
void
run_paired_ended(const bool VERBOSE,
                 const bool allow_ambig,
                 const string &reads_file1,
                 const string &reads_file2,
                 const size_t batch_size,
                 const size_t max_candidates,
                 const AbismalIndex &abismal_index,
                 pe_map_stats &pe_stats,
                 ostream &out) {
  const auto counter_st(begin(abismal_index.counter));
  const auto index_st(begin(abismal_index.index));
  const genome_iterator genome_st(begin(abismal_index.genome));

  ReadLoader rl1(reads_file1, batch_size);
  ReadLoader rl2(reads_file2, batch_size);
  ProgressBar progress(get_filesize(reads_file1), "mapping reads");

  omp_lock_t read_lock;
  omp_lock_t write_lock;

  omp_init_lock(&read_lock);
  omp_init_lock(&write_lock);

  if (VERBOSE)
    progress.report(cerr, 0);

  double start_time = omp_get_wtime();

#pragma omp parallel for
  for (int i = 0; i < omp_get_num_threads(); ++i) {
    if (random_pbat)
      map_paired_ended_rand(VERBOSE, allow_ambig, batch_size, max_candidates,
          abismal_index, rl1, rl2, read_lock, write_lock, pe_stats, out, progress);

    else
      map_paired_ended<conv>(VERBOSE, allow_ambig, batch_size, max_candidates,
          abismal_index, rl1, rl2, read_lock, write_lock, pe_stats, out, progress);
  }

  if (VERBOSE) {
    cerr << "-> total mapping time: " << omp_get_wtime() - start_time
         << "" << endl;
  }
}

static void
select_max_candidates(const uint32_t genome_size,
                      uint32_t &max_candidates) {
  static const double genome_frac = 1e-6;
  static const uint32_t min_max_candidates = 100u;

  const uint32_t c = static_cast<uint32_t>(genome_size * genome_frac);
  max_candidates = max(c, min_max_candidates);
}

int main(int argc, const char **argv) {

  try {
    static const string ABISMAL_VERSION = "0.2.0";
    string index_file;
    string outfile;
    string stats_outfile = "";
    bool VERBOSE = false;
    bool GA_conversion = false;
    bool allow_ambig = false;
    bool pbat_mode = false;
    bool random_pbat = false;
    uint32_t max_candidates = 0;
    size_t batch_size = 2000;
    int n_threads = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "map bisulfite converted reads",
                           "<reads-fq1> [<reads-fq2>]");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("index", 'i', "index file", true, index_file);
    opt_parse.add_opt("outfile", 'o', "output file", false, outfile);
    opt_parse.add_opt("mapstats", 'm',
                      "mapstats output file [stderr]",
                      false, stats_outfile);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
    opt_parse.add_opt("batch", 'b', "reads to load at once",
                      false, batch_size);
    opt_parse.add_opt("candidates", 'c', "max candidates for full comparison",
                      false, max_candidates);
    opt_parse.add_opt("max-mates", 'p', "max candidates as mates (pe mode)",
                      false, pe_candidates::max_size);
    opt_parse.add_opt("max-alignments", 's', "max single-end alignments",
                      false, se_result::max_size);
    opt_parse.add_opt("min-frag", 'l', "min fragment size (pe mode)",
                      false, pe_result::min_dist);
    opt_parse.add_opt("max-frag", 'L', "max fragment size (pe mode)",
                      false, pe_result::max_dist);
    opt_parse.add_opt("max-distance", 'M', "max fractional edit distance",
                      false, se_element::valid_frac);
    opt_parse.add_opt("ambig", 'a', "report a posn for ambiguous mappers",
                      false, allow_ambig);
    opt_parse.add_opt("pbat", 'P', "input follows the PBAT protocol",
                      false, pbat_mode);
    opt_parse.add_opt("random-pbat", 'R', "input follows random PBAT protocol",
                      false, random_pbat);
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
    if (leftover_args.size() != 1 && leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (n_threads <= 0) {
      cerr << "please choose a positive number of threads" << endl;
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

    omp_set_num_threads(n_threads);
    AbismalIndex::VERBOSE = VERBOSE;

    if (VERBOSE) {
      if (paired_end)
        cerr << "[mapping paired end: " << reads_file << " "
             << reads_file2 << "]\n";
      else
        cerr << "[mapping single end: " << reads_file << "]\n";
      cerr << "[output file: " << outfile << "]" << endl;
    }

    if (VERBOSE)
      cerr << "[loading abismal index]" << endl;
    AbismalIndex abismal_index;
    const double start_time = omp_get_wtime();
    abismal_index.read(index_file);

    if (VERBOSE)
      cerr << "[loading time: " << (omp_get_wtime() - start_time) << "]" << endl;

    if (max_candidates == 0)
      select_max_candidates(abismal_index.cl.get_genome_size(), max_candidates);

    // avoiding opening the stats output file until mapping is done
    se_map_stats se_stats;
    pe_map_stats pe_stats;

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str(), std::ios::binary);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    if (!out)
      throw runtime_error("failed to open output file: " + outfile);

    write_sam_header(abismal_index.cl.names, abismal_index.cl.starts,
                     "ABISMAL", ABISMAL_VERSION, argc, argv, out);

    if (reads_file2.empty()) {
      if (GA_conversion || pbat_mode)
        run_single_ended<a_rich, false>(VERBOSE, allow_ambig, reads_file,
            batch_size, max_candidates, abismal_index, se_stats, out);
      else if (random_pbat)
        run_single_ended<t_rich, true>(VERBOSE, allow_ambig, reads_file,
            batch_size, max_candidates, abismal_index, se_stats, out);
      else
        run_single_ended<t_rich, false>(VERBOSE, allow_ambig, reads_file,
            batch_size, max_candidates, abismal_index, se_stats, out);
    }
    else {
      if (pbat_mode)
        run_paired_ended<a_rich, false>(VERBOSE, allow_ambig, reads_file,
            reads_file2, batch_size, max_candidates, abismal_index, pe_stats,
            out);
      else if (random_pbat)
        run_paired_ended<t_rich, true>(VERBOSE,  allow_ambig, reads_file,
            reads_file2, batch_size, max_candidates, abismal_index, pe_stats,
            out);
      else
        run_paired_ended<t_rich, false>(VERBOSE, allow_ambig, reads_file,
            reads_file2, batch_size, max_candidates, abismal_index, pe_stats,
            out);
    }
    std::ofstream stats_of;
    if (!stats_outfile.empty())
      stats_of.open(stats_outfile.c_str(), std::ios::binary);
    std::ostream stats_out(stats_outfile.empty() ? std::cerr.rdbuf() : stats_of.rdbuf());
    stats_out << (reads_file2.empty() ?
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
