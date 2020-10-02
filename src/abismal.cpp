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
#include <array>
#include <stdexcept>

#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"
#include "OptionParser.hpp"
#include "zlib_wrapper.hpp"
#include "GenomicRegion.hpp"
#include "sam_record.hpp"
#include "bisulfite_utils.hpp"

#include "dna_four_bit_bisulfite.hpp"
#include "AbismalIndex.hpp"
#include "AbismalAlign.hpp"

#include <omp.h>

#define PREFETCH_LOOP 10U

using std::vector;
using std::array;
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
using std::to_string;

typedef uint16_t flags_t; // every bit is a flag
typedef int16_t score_t; // alignment score
typedef vector<uint8_t> Read; //4-bit encoding of reads
typedef genome_four_bit_itr genome_iterator; // iterates over 4 bits per byte

enum genome_pos_parity { pos_even = false, pos_odd = true };
enum conversion_type { t_rich = false, a_rich = true };

constexpr conversion_type
flip_conv(const conversion_type conv) {
  return conv == t_rich ? a_rich : t_rich;
}

// ADS: why did we decide to have the flag set to "1" for t-rich?

// functions to simultaneously get/set rc and a/t-richness flags
constexpr flags_t
get_strand_code(const char strand, const conversion_type conv) {
  return (((strand == '-')  ? samflags::read_rc : 0) |
          ((conv == t_rich) ? bsflags::read_is_t_rich: 0));
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
uint32_t ReadLoader::min_read_length = 32;

static void
update_max_read_length(size_t &max_length, const vector<string> &reads) {
  for (auto it (begin(reads)); it != end(reads); ++it)
    max_length = std::max(max_length, it->size());
}

struct se_element {
  uint32_t pos;
  score_t diffs;
  flags_t flags;

  se_element() :
    pos(0), diffs(invalid_hit_diffs),  flags(0) {}

  se_element(const uint32_t p, const score_t d,
             const flags_t f) :
    pos(p), diffs(d), flags(f) {}

  bool operator==(const se_element &rhs) const {
    return diffs == rhs.diffs && pos == rhs.pos;
  }

  // this is used to keep PE candidates sorted in the max heap
  bool operator<(const se_element &rhs) const {
    return diffs < rhs.diffs;
  }

  bool is_better_than (const se_element &rhs) const {
    return diffs < rhs.diffs;
  }

  bool rc() const {return samflags::check(flags, samflags::read_rc);}
  bool elem_is_a_rich() const {
    return !samflags::check(flags, bsflags::read_is_t_rich);
  }
  bool valid_hit() const {return diffs < invalid_hit_diffs;}
  bool valid() const {return diffs <= max_diffs;}
  void reset() { diffs = invalid_hit_diffs; }

  bool is_equal_to(const se_element &rhs) const {
    return (diffs == rhs.diffs) && (pos != rhs.pos);
  }
  static score_t invalid_hit_diffs;
  static score_t max_diffs;
  static uint16_t min_aligned_length;
};

// GS: must have a valid rationale for these numbers
score_t se_element::invalid_hit_diffs = 40;
score_t se_element::max_diffs = 10;
uint16_t se_element::min_aligned_length = 32;

struct se_result {
  se_element best;
  se_element second_best;
  se_result() : best(se_element()), second_best(se_element()) {}

  void update(const uint32_t p, const score_t d, const flags_t s) {
    // avoid having two copies of the best hit
    if (p == best.pos && s == best.flags) return;

    const se_element cand(p, d, s);
    if (cand.is_better_than(second_best)) {
      second_best = cand;
      if (second_best.is_better_than(best)) std::swap(best, second_best);
    }
  }

  // sorts the first and second best candidates. This is used after
  // aligning the two candidates to account for the possibility that
  // second has a better alignment score than first
  bool sort_candidates() {
    if (second_best.is_better_than(best)) {
      std::swap(best, second_best);
      return true;
    }
    return false;
  }

  bool ambig() const {
    return best.diffs == second_best.diffs;
  }

  bool sure_ambig(uint32_t seed_number = 0) const {
    return ambig() &&
      (best.diffs == 0 || (best.diffs == 1 && seed_number > 0));
  }

  bool should_report() const {
    return best.valid() && !ambig();
  }

  void reset() {best.reset(); second_best.reset();}
  uint32_t get_cutoff() const {return second_best.diffs;}
};


inline bool
chrom_and_posn(const ChromLookup &cl, const string &cig, const uint32_t p,
               uint32_t &r_p, uint32_t &r_e, uint32_t &r_chr) {
  const uint32_t ref_ops = cigar_rseq_ops(cig);
  if (!cl.get_chrom_idx_and_offset(p, ref_ops, r_chr, r_p)) return false;
  r_e = r_p + ref_ops;
  return true;
}

static bool
format_se(const bool allow_ambig, se_result res, const ChromLookup &cl,
          string &read, const string &read_name, string &cigar,
          ostream &out) {

  const se_element s = res.best;

  // GS: not the same as "should_report" because we need to account
  // for the allow_ambig flag
  if (!s.valid() || (!allow_ambig && res.ambig()))
    return false;

  uint32_t ref_s = 0, ref_e = 0, chrom_idx = 0;
  if (chrom_and_posn(cl, cigar, s.pos, ref_s, ref_e, chrom_idx)) {
    sam_rec sr(read_name, 0, cl.names[chrom_idx], ref_s + 1,
               255, cigar, "*", 0, 0, read, "*");
    if (s.rc())
      set_flag(sr, samflags::read_rc);

    sr.add_tag("NM:i:" + to_string(s.diffs));
    sr.add_tag(s.elem_is_a_rich() ? "CV:A:A" : "CV:A:T");
    out << sr << '\n';
    return true;
  }
  return false;
}

struct pe_element {
  se_element r1;
  se_element r2;

  pe_element() : r1(se_element()), r2(se_element()) {}
  pe_element(const se_element &s1, const se_element &s2) : r1(s1), r2(s2) {}

  bool rc() const { return r1.rc(); }
  bool elem_is_a_rich() const {return r1.elem_is_a_rich();}
  score_t diffs() const { return r1.diffs + r2.diffs; }

  // essentially checks for a concordant pair because "update" in
  // pe_result is only called for concordant pairs
  bool valid_hit() const {
    return r1.valid_hit() && r2.valid_hit();
  }

  bool valid() const {
    return r1.diffs + r2.diffs <= 2*se_element::max_diffs;
  }

  bool is_equal_to(const pe_element &rhs) const {
    return diffs() == rhs.diffs() && !(r1.pos == rhs.r1.pos &&
                                       r2.pos == rhs.r2.pos);
  }

  bool is_better_than(const pe_element &rhs) const {
    return diffs() < rhs.diffs();
  }

  void reset() {r1.reset(); r2.reset();}
  static uint32_t min_dist;
  static uint32_t max_dist;
};

uint32_t pe_element::min_dist = 32;
uint32_t pe_element::max_dist = 3000;

struct pe_result { // assert(sizeof(pe_result) == 16);
  pe_element best;
  pe_element second_best;
  pe_result() {}
  pe_result(const pe_element &a, const pe_element &b)
    : best(a), second_best(b) {}

  // true if best was updated. We use the value returned by this
  // function to update the reported cigar if necessary
  bool update(const pe_element &cand) {
    if (cand.is_better_than(second_best)) {
      second_best = cand;

      if (second_best.is_better_than(best)) {
        std::swap(second_best, best);
        return true;
      }
      return false;
    }
    return false;
  }
  void reset() {
    best.reset();
    second_best.reset();
  }

  bool ambig() const {
    return best.diffs() == second_best.diffs();
  }

  bool should_report() const {
    return (best.valid() && !ambig());
  }
};

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
bool
format_pe(const bool allow_ambig,
          const pe_result &res, const ChromLookup &cl,
          const string &read1, const string &read2,
          const string &name1, const string &name2,
          const string &cig1,  const string &cig2,
          ostream &out) {

  uint32_t r_s1 = 0, r_e1 = 0, chr1 = 0; // positions in chroms (0-based)
  uint32_t r_s2 = 0, r_e2 = 0, chr2 = 0;
  const pe_element p = res.best;

  // score too low or ambiguous not required (we do not report a
  // random mate for paired-end reads
  if (!p.valid() || res.ambig())
    return false;

  // PE chromosomes differ or couldn't be found, treat read as unmapped
  if (!chrom_and_posn(cl, cig1, p.r1.pos, r_s1, r_e1, chr1) ||
      !chrom_and_posn(cl, cig2, p.r2.pos, r_s2, r_e2, chr2) || chr1 != chr2)
    return false;

  const bool rc = p.rc();

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

  out << sr1 << '\n'
      << sr2 << '\n';

  return true;
}

struct pe_candidates {
  pe_candidates() : v(vector<se_element>(max_size)), sz(1) {}
  bool full() const {return sz == max_size;}
  void reset() {v.front().reset(); sz = 1;}
  score_t get_cutoff() const {return v.front().diffs;}
  void update(const uint32_t p, const score_t d, const flags_t s) {
    if (full()) {
      if (d < v.front().diffs) {
        std::pop_heap(begin(v), end(v));
        v.back() = se_element(p, d, s);
        std::push_heap(begin(v), end(v));
      }
    }
    else if (d < se_element::invalid_hit_diffs) {
      v[sz++] = se_element(p, d, s);
      std::push_heap(begin(v), begin(v) + sz);
    }
  }
  bool sure_ambig(uint32_t seed_number = 0) const {
    return full() && (v[0].diffs == 0 || (v[0].diffs == 1 && seed_number != 0));
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
uint32_t pe_candidates::max_size = 20;


inline double pct(const double a, const double b) {return 100.0*a/b;}

struct se_map_stats {
  se_map_stats() :
    tot_rds(0), uniq_rds(0), ambig_rds(0), unmapped_rds(0), skipped_rds(0) {}
  uint32_t tot_rds;
  uint32_t uniq_rds;
  uint32_t ambig_rds;
  uint32_t unmapped_rds;
  uint32_t skipped_rds;

  void update(const string &read, const se_result &res,
              const bool se_reported) {
    ++tot_rds;
    if (res.best.valid()) {
      uniq_rds += se_reported;
      ambig_rds += res.ambig();
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

  void update_pair(const pe_result &res,
                   const bool pe_reported) {
    ++tot_pairs;
    if (res.best.valid()) {
      ambig_pairs += res.ambig();
      uniq_pairs += pe_reported;
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
update_pe_stats(const pe_result &best,
                const se_result &se1, const se_result &se2,
                const string &read1, const string &read2,
                const bool pe_reported,
                const bool se1_reported,
                const bool se2_reported,
                pe_map_stats &pe_stats) {
  pe_stats.update_pair(best, pe_reported);
  if (!pe_reported) {
    pe_stats.end1_stats.update(read1, se1, se1_reported);
    pe_stats.end2_stats.update(read2, se2, se2_reported);
  }
}

static void
select_output(const bool allow_ambig, const ChromLookup &cl,
              pe_result &best, se_result &se1, se_result &se2,
              string &read1, const string &name1,
              string &read2, const string &name2,
              string &cig1, string &cig2,
              bool &pe_reported,
              bool &se1_reported,
              bool &se2_reported,
              ostream &out) {

  pe_reported = format_pe(allow_ambig, best, cl, read1, read2, name1, name2,
                          cig1, cig2, out);
  if (!pe_reported) {
    se1_reported = format_se(allow_ambig, se1, cl, read1, name1,
                             cig1, out);
    se2_reported = format_se(allow_ambig, se2, cl, read2, name2,
                             cig2, out);
  }
}


inline bool
the_comp(const char a, const char b) {
  return (a & b) == 0;
}

template<const genome_pos_parity parity>
score_t
full_compare(const score_t cutoff,
             Read::const_iterator read_itr,
             const Read::const_iterator &read_end_low,
             const Read::const_iterator &read_end_high,
             Genome::const_iterator genome_itr) {
  score_t d = 0;
  auto g_i = genome_itr;
  while (d <= cutoff && read_itr != read_end_low) {
    d += the_comp(*read_itr, *g_i);
    ++read_itr, ++g_i;
  }

  g_i = genome_itr;
  if (parity == pos_odd) ++g_i; // odd case: first pos does only 1 comp
  while (d <= cutoff && read_itr != read_end_high) {
    d += the_comp(*read_itr, *g_i);
    ++read_itr, ++g_i;
  }
  return d;
}

template <const uint16_t strand_code, class result_type>
void
check_hits(vector<uint32_t>::const_iterator start_idx,
           const vector<uint32_t>::const_iterator end_idx,
           const Read::const_iterator even_read_st,
           const Read::const_iterator even_read_mid,
           const Read::const_iterator even_read_end,
           const Read::const_iterator odd_read_st,
           const Read::const_iterator odd_read_mid,
           const Read::const_iterator odd_read_end,
           const Genome::const_iterator genome_st,
           const uint32_t offset,
           result_type &res) {
  for (; start_idx != end_idx /*&& !res.sure_ambig(offset != 0)*/; ++start_idx) {

    /* GS: This prefetch adds the genome portion that will be compared
     * 10 iterations ahead to the cache. *(start_idx + PREFETCH_LOOP)
     * is this position. We then correct the read offset and divide by
     * 2 for the same reason as described below, where we adjust to
     * the fact that 2 bases are stored in each byte */
    __builtin_prefetch(
        &(*(genome_st + ((*(start_idx + PREFETCH_LOOP) - offset) >> 1))),
      0, 0);

    const uint32_t pos = (*start_idx) - offset;

    // ADS: (reminder) the adjustment below is because 2 bases are
    // stored in each byte
    const score_t diffs =
      ((pos & 1) ?
       full_compare<pos_odd>(res.get_cutoff(),
                             odd_read_st, odd_read_mid, odd_read_end,
                             genome_st + (pos >> 1)) :
       full_compare<pos_even>(res.get_cutoff(),
                              even_read_st, even_read_mid, even_read_end,
                              genome_st + (pos >> 1)));
    res.update(pos, diffs, strand_code);
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

void
find_candidates(const Read::const_iterator read_start,
                const genome_iterator gi,
                const uint32_t read_lim, // not necessarily read len
                vector<uint32_t>::const_iterator &low,
                vector<uint32_t>::const_iterator &high) {
  size_t p = seed::key_weight;
  const size_t lim = std::min(read_lim, seed::n_sorting_positions);
  while (p != lim) {
    auto first_1 = lower_bound(low, high, 1, compare_bases(gi + p));
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

template <const uint16_t strand_code, class Read, class result_type>
void
process_seeds(const uint32_t max_candidates,
              const AbismalIndex &abismal_index,
              const Genome::const_iterator genome_st,
              const genome_iterator gi,
              const Read &read_seed,
              const Read &read_even,
              const Read &read_odd,
              result_type &res) {

  const uint32_t readlen = read_seed.size();

  // seed iterators
  const auto read_start(begin(read_seed));
  const auto read_end(end(read_seed));

  // even comparison iterators
  const auto even_read_start(begin(read_even));
  const auto even_read_mid(begin(read_even) + ((readlen + 1)/2));
  const auto even_read_end(end(read_even));

  // odd comparison iterators
  const auto odd_read_start(begin(read_odd));
  const auto odd_read_mid(begin(read_odd) + ((readlen + 1)/2));
  const auto odd_read_end(end(read_odd));

  // ADS: there is something wrong with this code if so many
  // conditions are needed
  const uint32_t shift_lim =
    readlen > (seed::n_seed_positions + seed::index_interval - 1) ?
    (readlen - seed::n_seed_positions - seed::index_interval + 1) : 0;

  const uint32_t shift = (seed::n_shifts == 1u) ? (shift_lim + 1):
    std::max(1u, shift_lim/(seed::n_shifts - 1));

  bool found_good_seed = false;

  // uniformly spaced seeds
  for (uint32_t i = 0; i <= shift_lim && !res.sure_ambig(i); i += shift) {
    std::array<uint32_t, seed::index_interval> k = {0};
    get_1bit_hash_compressed(read_start + i, k);
    for (uint32_t j = 0; j < seed::index_interval; ++j) {
      const uint32_t offset = i + j;
      auto s_idx(index_st + *(counter_st + k[j]));
      auto e_idx(index_st + *(counter_st + k[j] + 1));

      if (s_idx < e_idx) {
        find_candidates(read_start + offset, genome_st, readlen - offset, s_idx, e_idx);
        if (e_idx - s_idx < max_candidates) {
          found_good_seed = true;
          check_hits<strand_code>(s_idx, e_idx,
                                  even_read_start, even_read_mid, even_read_end,
                                  odd_read_start, odd_read_mid, odd_read_end,
                                  genome_st.itr, offset, res);
        }
      }
    }
  }

  // if no good seeds found, use entire read as seed
  if (!found_good_seed) {
    std::array<uint32_t, seed::index_interval> k = {0};
    get_1bit_hash_compressed(read_start, k);
    for (uint32_t j = 0; j < seed::index_interval; ++j) {
      auto s_idx(index_st + *(counter_st + k[j]));
      auto e_idx(index_st + *(counter_st + k[j] + 1));
      if (s_idx < e_idx) {
        find_candidates(read_start, genome_st,
                        readlen - seed::index_interval + 1, s_idx, e_idx);

        // GS: seed::n_shifts * max candidates is the max acceptable
        // number of searches for a read under no good seed condition
        if (e_idx - s_idx < 80000)
          check_hits<strand_code>(s_idx, e_idx,
                                  even_read_start, even_read_mid, even_read_end,
                                  odd_read_start, odd_read_mid, odd_read_end,
                                  genome_st.itr, j, res);
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

/* Creates reads meant for comparison on a compressed genome.
 * This function is used to accelerate comparison between reads and
 * genome positions using a reference genome encoded in four bits
 * per
 * base, that is, two bases are present in each byte, with even
 * positions in lower bits and odd positions in higher bits.
 *
 * The function encodes the read in two ways. Every pread_* object
 * is a vector of characters, where each element of the vector has
 * only one active bit representing A, T, G or C.
 *
 * In pread_seed (input), we have encoded bases using the lower 4
 * bits
 * In pread_even, which is used to compare the read to even genome
 * positions, we put the bases in even positions as the first
 * elements, and the odd positions shifted by 4. This allows the &
 * operator to be used on four-bit comparisons.
 *
 * In pread_odd, used to compare to odd positions in the genome, we
 * first put the odd bases in the start, shifted by 4, then the
 * even
 * bases.
 *
 * In both comparison cases, the number of mismatches can be
 * obtained
 * by iterating from start to finish through the pread_even (resp
 * pread_odd) objects. The genome, however, has to be "rewinded"
 * once
 * we reach half of the read. This logic is implemented in the
 * "full_compare" function above. */

static void
prep_for_seeds(const Read &pread_seed, Read &pread_even,
               Read &pread_odd) {
  const size_t sz = pread_seed.size();
  pread_even.resize(sz);
  pread_odd.resize(sz);
  size_t i, j = 0;
  for (i = 0; i < sz; i += 2) pread_even[j++] = pread_seed[i];
  for (i = 1; i < sz; i += 2) pread_even[j++] = (pread_seed[i] << 4);

  j = 0;
  for (i = 0; i < sz; i += 2) pread_odd[j++] = (pread_seed[i] << 4);
  for (i = 1; i < sz; i += 2) pread_odd[j++] = pread_seed[i];
}


// this lookup improves speed when running alignment because we do not
// have to check if bases are equal or different
static inline score_t
mismatch_score(const char q_base, const uint8_t t_base) {
  return simple_local_alignment::score_lookup[the_comp(q_base, t_base)];
}

using AbismalAlignSimple = AbismalAlign<mismatch_score,
                                        simple_local_alignment::indel>;

// GS: does not consider soft-clipping into edit distance
// but requires a minimum alignment length
static score_t
edit_distance(const score_t aln_score, const size_t len) {
  if (len < se_element::min_aligned_length)
    return se_element::invalid_hit_diffs;

  // GS: this conversion only works for 1 -1 -1
  return (static_cast<score_t>(len) - aln_score) / 2;
}

// This function alings the read and converts "diff" from Hamming
// distance to edit distance
static void
align_read(se_element &res, string &cigar, const string &read,
           Read &pread, AbismalAlignSimple &aln) {
  // ends early if mismatches are good enough
  if (res.valid()) cigar = std::to_string(read.length()) + "M";
  else {
    uint32_t len = 0; // the region of the read the alignment spans
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
    const score_t aln_score = aln.align(pread, res.pos, len, cigar);
    res.diffs = edit_distance(aln_score, len);
  }
}

template <const  conversion_type conv>
void
map_single_ended(const bool VERBOSE,
                 const bool allow_ambig,
                 const string &reads_file,
                 const size_t batch_size,
                 const size_t max_candidates,
                 const AbismalIndex &abismal_index,
                 se_map_stats &se_stats,
                 ostream &out) {

  const uint32_t genome_size = abismal_index.cl.get_genome_size();
  const Genome::const_iterator genome_st(begin(abismal_index.genome));
  const genome_iterator gi(genome_st);

  size_t max_batch_read_length;
  vector<string> names, reads;
  reads.reserve(batch_size);
  names.reserve(batch_size);
  vector<se_result> res(batch_size);
  vector<string> cigar(batch_size);

  ReadLoader rl(reads_file, batch_size);

  ProgressBar progress(get_filesize(reads_file), "mapping reads");
  if (VERBOSE)
    progress.report(cerr, 0);

  double total_mapping_time = 0;
  while (rl.good()) {

    if (VERBOSE && progress.time_to_report(rl.get_current_byte()))
      progress.report(cerr, rl.get_current_byte());

    rl.load_reads(names, reads);

    max_batch_read_length = 0;
    update_max_read_length(max_batch_read_length, reads);
    const size_t n_reads = reads.size();

#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i) res[i].reset();

    const double start_time = omp_get_wtime();

#pragma omp parallel
    {
      Read pread_seed, pread_even, pread_odd;
#pragma omp for
      for (size_t i = 0; i < reads.size(); ++i) {
        if (!reads[i].empty()) {

          prep_read<conv>(reads[i], pread_seed);
          prep_for_seeds(pread_seed, pread_even, pread_odd);
          process_seeds<get_strand_code('+', conv)>(max_candidates,
                                                    abismal_index, genome_st,
                                                    gi, pread_seed, pread_even,
                                                    pread_odd, res[i]);

          const string read_rc(revcomp(reads[i]));
          prep_read<!conv>(read_rc, pread_seed);
          prep_for_seeds(pread_seed, pread_even, pread_odd);
          process_seeds<get_strand_code('-', conv)>(max_candidates,
                                                    abismal_index, genome_st,
                                                    gi, pread_seed, pread_even,
                                                    pread_odd, res[i]);
        }
      }
    }
    total_mapping_time += (omp_get_wtime() - start_time);

#pragma omp parallel
    {
      Read pread;
      string tmp_cigar;
      AbismalAlignSimple aln(gi, genome_size, max_batch_read_length);

#pragma omp for
      for (size_t i = 0; i < n_reads; ++i) {
        if (res[i].best.valid_hit())
          align_read(res[i].best, cigar[i], reads[i], pread, aln);
        if (res[i].second_best.valid_hit())
          align_read(res[i].second_best, tmp_cigar, reads[i], pread, aln);

        if (res[i].sort_candidates())
          cigar[i] = tmp_cigar;
      }
    }

    for (size_t i = 0 ; i < n_reads; ++i) {
      const bool reported = format_se(allow_ambig, res[i], abismal_index.cl,
                                      reads[i], names[i], cigar[i], out);
      se_stats.update(reads[i], res[i], reported);
    }
  }

  if (VERBOSE) {
    progress.report(cerr, get_filesize(reads_file));
    cerr << "[total mapping time: " << total_mapping_time << "]" << endl;
  }
}

void
map_single_ended_rand(const bool VERBOSE,
                      const bool allow_ambig,
                      const string &reads_file,
                      const size_t batch_size,
                      const size_t max_candidates,
                      const AbismalIndex &abismal_index,
                      se_map_stats &se_stats,
                      ostream &out) {
  const uint32_t genome_size = abismal_index.cl.get_genome_size();
  const Genome::const_iterator genome_st(begin(abismal_index.genome));
  const genome_iterator gi(genome_st);

  size_t max_batch_read_length;
  vector<string> names, reads;
  reads.reserve(batch_size);
  names.reserve(batch_size);
  vector<string> cigar(batch_size);
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
    max_batch_read_length = 0;
    update_max_read_length(max_batch_read_length, reads);

    const size_t n_reads = reads.size();

#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i) res[i].reset();

    const double start_time = omp_get_wtime();
#pragma omp parallel
    {
      Read pread_seed, pread_even, pread_odd;
#pragma omp for
      for (size_t i = 0; i < n_reads; ++i)
        if (!reads[i].empty()) {
          prep_read<t_rich>(reads[i], pread_seed);
          prep_for_seeds(pread_seed, pread_even, pread_odd);
          process_seeds<get_strand_code('+', t_rich)>(max_candidates,
                                                      abismal_index, genome_st,
                                                      gi, pread_seed, pread_even,
                                                      pread_odd, res[i]);
          prep_read<a_rich>(reads[i], pread_seed);
          prep_for_seeds(pread_seed, pread_even, pread_odd);
          process_seeds<get_strand_code('+', a_rich)>(max_candidates,
                                                      abismal_index, genome_st,
                                                      gi, pread_seed, pread_even,
                                                      pread_odd, res[i]);

          const string read_rc(revcomp(reads[i]));
          prep_read<t_rich>(read_rc, pread_seed);
          prep_for_seeds(pread_seed, pread_even, pread_odd);
          process_seeds<get_strand_code('-', a_rich)>(max_candidates,
                                                      abismal_index, genome_st,
                                                      gi, pread_seed, pread_even,
                                                      pread_odd, res[i]);

          prep_read<a_rich>(read_rc, pread_seed);
          prep_for_seeds(pread_seed, pread_even, pread_odd);
          process_seeds<get_strand_code('-', t_rich)>(max_candidates,
                                                      abismal_index, genome_st,
                                                      gi, pread_seed, pread_even,
                                                      pread_odd, res[i]);
        }
    }
    total_mapping_time += (omp_get_wtime() - start_time);
#pragma omp parallel
    {
      Read pread;
      string tmp_cigar;
      AbismalAlignSimple aln(gi, genome_size, max_batch_read_length);

#pragma omp for
      for (size_t i = 0; i < reads.size(); ++i) {
        if (res[i].best.valid_hit())
          align_read(res[i].best, cigar[i], reads[i], pread, aln);

        if (res[i].second_best.valid_hit())
          align_read(res[i].second_best, tmp_cigar, reads[i], pread, aln);

        if (res[i].sort_candidates())
          cigar[i] = tmp_cigar;
      }
    }

    for (size_t i = 0 ; i < n_reads; ++i) {
      const bool reported = format_se(allow_ambig, res[i], abismal_index.cl,
                                     reads[i], names[i], cigar[i], out);
      se_stats.update(reads[i], res[i], reported);
    }
  }
  if (VERBOSE) {
    progress.report(cerr, get_filesize(reads_file));
    cerr << "[total mapping time: " << total_mapping_time << "]" << endl;
  }
}

template <const bool cmp,
          const uint16_t strand_code1, const uint16_t strand_code2,
          class result_type>
void
map_pe_batch(const vector<string> &reads1, const vector<string> &reads2,
             const uint32_t max_candidates,
             const AbismalIndex &abismal_index,
             const Genome::const_iterator genome_st,
             const genome_iterator gi,
             vector<result_type> &res1, vector<result_type> &res2) {
#pragma omp parallel
  {
    Read pread_seed, pread_even, pread_odd;

#pragma omp for
    for (size_t i = 0 ; i < reads1.size(); ++i) {
      res1[i].reset();
      res2[i].reset();

      if (!reads1[i].empty()) {
        prep_read<cmp>(reads1[i], pread_seed);
        prep_for_seeds(pread_seed, pread_even, pread_odd);
        process_seeds<strand_code1>(max_candidates, abismal_index, genome_st,
                                    gi, pread_seed, pread_even, pread_odd,
                                    res1[i]);
      }
      if (!reads2[i].empty()) {
        const string read_rc(revcomp(reads2[i]));
        prep_read<cmp>(read_rc, pread_seed);
        prep_for_seeds(pread_seed, pread_even, pread_odd);
        process_seeds<strand_code2>(max_candidates, abismal_index, genome_st,
                                    gi, pread_seed, pread_even, pread_odd,
                                    res2[i]);
      }
    }
  }
}

static void
best_single(const pe_candidates &pres, se_result &res) {
  const auto lim(begin(pres.v) + pres.sz);

  // get best and second best by mismatch
  for (auto i(begin(pres.v)); i != lim && !res.sure_ambig(0); ++i)
    res.update(i->pos, i->diffs, i->flags);
}

template <const bool swap_ends>
static void
best_pair(const pe_candidates &res1, const pe_candidates &res2,
          const string &read1, const string &read2,
          string &cig1, string &cig2,
          AbismalAlignSimple &aln, pe_result &best) {

  auto j1 = begin(res1.v);
  const auto j1_end = j1 + res1.sz;
  const auto j2_end = begin(res2.v) + res2.sz;
  Read pread;
  se_element s1, s2;
  string cand_cig1, cand_cig2;
  for (auto j2(begin(res2.v)); j2 != j2_end; ++j2) {
    s2 = *j2;
    if (s2.valid_hit()){
      bool aligned_s2 = false;
      uint32_t lim_s2 = 0;

      const uint32_t lim = j2->pos + read2.length();
      while (j1 != j1_end && j1->pos + pe_element::max_dist < lim)
        ++j1;

      while (j1 != j1_end && j1->pos + pe_element::min_dist <= lim) {
        s1 = *j1;
        if (s1.valid_hit()) {
          align_read(s1, cand_cig1, read1, pread, aln);
          if (!aligned_s2) {
            align_read(s2, cand_cig2, read2, pread, aln);
            lim_s2 = s2.pos + cigar_rseq_ops(cand_cig2);
            aligned_s2 = true;
          }

          // get length after alignment, and only accept if it is
          // still within fragment limits
          if ((s1.pos + pe_element::max_dist >= lim_s2) &&
              (s1.pos + pe_element::min_dist <= lim_s2)) {
            const pe_element p(swap_ends ? s2 : s1, swap_ends ? s1 : s2);
            if (best.update(p)) {
              cig1 = cand_cig1;
              cig2 = cand_cig2;
            }
          }
        }
        ++j1;
      }
    }
  }
}

template <const bool swap_ends>
static void
select_maps(const string &read1, const string &read2,
            string &cig1, string &cig2,
            pe_candidates &res1, pe_candidates &res2,
            se_result &res_se1, se_result &res_se2,
            AbismalAlignSimple &aln, pe_result &best) {
  res1.prepare_for_mating();
  res2.prepare_for_mating();
  best_pair<swap_ends>(res1, res2, read1, read2, cig1, cig2, aln, best);

  // if PE should not be reported, try to find the best single
  best_single(res1, res_se1);
  best_single(res2, res_se2);
}

template <const conversion_type conv>
void
map_paired_ended(const bool VERBOSE,
                 const bool allow_ambig,
                 const string &reads_file1,
                 const string &reads_file2,
                 const size_t batch_size,
                 const size_t max_candidates,
                 const AbismalIndex &abismal_index,
                 pe_map_stats &pe_stats,
                 ostream &out) {
  double total_mapping_time = 0;

  ReadLoader rl1(reads_file1, batch_size);
  ReadLoader rl2(reads_file2, batch_size);

  size_t max_batch_read_length;
  vector<string> names1(batch_size), reads1(batch_size), cigar1(batch_size);
  vector<string> names2(batch_size), reads2(batch_size), cigar2(batch_size);
  vector<pe_candidates> res1(batch_size), res2(batch_size);
  vector<pe_result> bests(batch_size);
  vector<se_result> res_se1(batch_size), res_se2(batch_size);

  const uint32_t genome_size = abismal_index.cl.get_genome_size();
  const Genome::const_iterator genome_st(begin(abismal_index.genome));
  const genome_iterator gi(genome_st);

  bool pe_reported, se1_reported, se2_reported;
  ProgressBar progress(get_filesize(reads_file1), "mapping reads");
  if (VERBOSE)
    progress.report(cerr, 0);

  while (rl1.good() && rl2.good()) {
    if (VERBOSE && progress.time_to_report(rl1.get_current_byte()))
      progress.report(cerr, rl1.get_current_byte());

    rl1.load_reads(names1, reads1);
    rl2.load_reads(names2, reads2);

    // used to get AbismalAlign size
    max_batch_read_length = 0;
    update_max_read_length(max_batch_read_length, reads1);
    update_max_read_length(max_batch_read_length, reads2);

    const size_t n_reads = reads1.size();
    const double start_time = omp_get_wtime();

#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i) {
      res_se1[i].reset();
      res_se2[i].reset();
      bests[i].reset();
    }

    map_pe_batch<conv,
                 get_strand_code('+', conv),
                 get_strand_code('-', flip_conv(conv))>(reads1, reads2,
                                                        max_candidates,
                                                        abismal_index,
                                                        genome_st, gi,
                                                        res1, res2);

#pragma omp parallel
    {
      AbismalAlignSimple aln(gi, genome_size, max_batch_read_length);

#pragma omp for
      for (size_t i = 0 ; i < n_reads; ++i)
        select_maps<false>(reads1[i], reads2[i], cigar1[i], cigar2[i],
                           res1[i], res2[i], res_se1[i], res_se2[i],
                           aln, bests[i]);
    }

    map_pe_batch<!conv,
                 get_strand_code('+', flip_conv(conv)),
                 get_strand_code('-', conv)>(reads2, reads1, max_candidates,
                                             abismal_index, genome_st, gi,
                                             res2, res1);

#pragma omp parallel
    {
      AbismalAlignSimple aln(gi, genome_size, max_batch_read_length);

#pragma omp for
      for (size_t i = 0 ; i < n_reads; ++i)
        // ADS: note res2 and res1 are swapped below
        select_maps<true>(reads2[i], reads1[i], cigar2[i], cigar1[i],
                          res2[i], res1[i], res_se2[i], res_se1[i],
                          aln, bests[i]);
    }

#pragma omp parallel
    {
      Read pread;
      string tmp_cigar;
      AbismalAlignSimple aln(gi, genome_size, max_batch_read_length);

#pragma omp for
      for (size_t i = 0; i < n_reads; ++i) {
        if (!bests[i].should_report()) {
          if (res_se1[i].best.valid_hit())
            align_read(res_se1[i].best, cigar1[i], reads1[i], pread, aln);

          if (res_se1[i].second_best.valid_hit())
            align_read(res_se1[i].second_best, tmp_cigar, reads1[i], pread, aln);

          if (res_se1[i].sort_candidates())
            cigar1[i] = tmp_cigar;

          if (res_se2[i].best.valid_hit())
            align_read(res_se2[i].best, cigar2[i], reads2[i], pread, aln);

          if (res_se2[i].second_best.valid_hit())
            align_read(res_se2[i].second_best, tmp_cigar, reads2[i], pread, aln);

          if (res_se2[i].sort_candidates())
            cigar2[i] = tmp_cigar;
        }
      }
    }

    pe_reported = false;
    se1_reported = false;
    se2_reported = false;
    for (size_t i = 0 ; i < n_reads; ++i) {
      select_output(allow_ambig, abismal_index.cl,
                    bests[i], res_se1[i], res_se2[i],
                    reads1[i], names1[i],
                    reads2[i], names2[i],
                    cigar1[i], cigar2[i], pe_reported,
                    se1_reported, se2_reported, out);

      update_pe_stats(bests[i], res_se1[i], res_se2[i], reads1[i],
                      reads2[i], pe_reported, se1_reported,
                      se2_reported, pe_stats);
    }
    total_mapping_time += (omp_get_wtime() - start_time);
  }

  if (VERBOSE) {
    progress.report(cerr, get_filesize(reads_file1));
    cerr << "[total mapping time: " << total_mapping_time << "]" << endl;
  }
}

static void
map_paired_ended_rand(const bool VERBOSE,
                      const bool allow_ambig,
                      const string &reads_file1,
                      const string &reads_file2,
                      const size_t batch_size,
                      const size_t max_candidates,
                      const AbismalIndex &abismal_index,
                      pe_map_stats &pe_stats,
                      ostream &out) {
  double total_mapping_time = 0;

  ReadLoader rl1(reads_file1, batch_size);
  ReadLoader rl2(reads_file2, batch_size);

  size_t max_batch_read_length;
  vector<string> names1(batch_size), reads1(batch_size), cigar1(batch_size);
  vector<string> names2(batch_size), reads2(batch_size), cigar2(batch_size);

  vector<pe_candidates> res1(batch_size), res2(batch_size);
  vector<pe_result> bests(batch_size);
  vector<se_result> res_se1(batch_size), res_se2(batch_size);

  const uint32_t genome_size = abismal_index.cl.get_genome_size();
  const Genome::const_iterator genome_st(begin(abismal_index.genome));
  const genome_iterator gi(genome_st);

  bool pe_reported, se1_reported, se2_reported;
  ProgressBar progress(get_filesize(reads_file1), "mapping reads");
  if (VERBOSE)
    progress.report(cerr, 0);

  while (rl1.good() && rl2.good()) {

    if (VERBOSE && progress.time_to_report(rl1.get_current_byte()))
      progress.report(cerr, rl1.get_current_byte());

    rl1.load_reads(names1, reads1);
    rl2.load_reads(names2, reads2);

    max_batch_read_length = 0;
    update_max_read_length(max_batch_read_length, reads1);
    update_max_read_length(max_batch_read_length, reads2);
    const size_t n_reads = reads1.size();

    const double start_time = omp_get_wtime();
#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i) {
      res_se1[i].reset();
      res_se2[i].reset();
      bests[i].reset();
    }

    // t-rich end1, pos-strand end1
    map_pe_batch<t_rich,
                 get_strand_code('+', t_rich),
                 get_strand_code('-', a_rich)>(reads1, reads2, max_candidates,
                                               abismal_index, genome_st, gi,
                                               res1, res2);
#pragma omp parallel
    {
      AbismalAlignSimple aln(gi, genome_size, max_batch_read_length);

#pragma omp for
      for (size_t i = 0 ; i < n_reads; ++i)
        select_maps<false>(reads1[i], reads2[i],
                           cigar1[i], cigar2[i],
                           res1[i], res2[i],
                           res_se1[i], res_se2[i],
                           aln, bests[i]);
    }
    // t-rich end1, neg-strand end1
    map_pe_batch<a_rich,
                 get_strand_code('+', a_rich),
                 get_strand_code('-', t_rich)>(reads2, reads1, max_candidates,
                                               abismal_index, genome_st, gi,
                                               res2, res1);
#pragma omp parallel
    {
      AbismalAlignSimple aln(gi, genome_size, max_batch_read_length);

#pragma omp for
      for (size_t i = 0 ; i < n_reads; ++i)
        select_maps<true>(reads2[i], reads1[i],
                          cigar2[i], cigar1[i],
                          res2[i], res1[i],
                          res_se2[i], res_se1[i],
                          aln, bests[i]);
    }
    // a-rich end1, pos-strand end1
    map_pe_batch<a_rich,
                 get_strand_code('+', a_rich),
                 get_strand_code('-', t_rich)>(reads1, reads2, max_candidates,
                                               abismal_index, genome_st, gi,
                                               res1, res2);
#pragma omp parallel
    {
      AbismalAlignSimple aln(gi, genome_size, max_batch_read_length);

#pragma omp for
      for (size_t i = 0 ; i < n_reads; ++i)
        select_maps<false>(reads1[i], reads2[i],
                           cigar1[i], cigar2[i],
                           res1[i], res2[i],
                           res_se1[i], res_se2[i],
                           aln, bests[i]);
    }
    // a-rich end1, neg-strand end1
    map_pe_batch<t_rich,
                 get_strand_code('+', t_rich),
                 get_strand_code('-', a_rich)>(reads2, reads1, max_candidates,
                                               abismal_index, genome_st, gi,
                                               res2, res1);
#pragma omp parallel
    {
      AbismalAlignSimple aln(gi, genome_size, max_batch_read_length);

#pragma omp for
      for (size_t i = 0 ; i < n_reads; ++i)
        select_maps<true>(reads2[i], reads1[i],
                          cigar2[i], cigar1[i],
                          res2[i], res1[i],
                          res_se2[i], res_se1[i],
                          aln, bests[i]);
    }

    // only align singles if no concordant pair
#pragma omp parallel
    {
      Read pread;
      string tmp_cigar;
      AbismalAlignSimple aln(gi, genome_size, max_batch_read_length);

#pragma omp for
      for (size_t i = 0; i < n_reads; ++i) {
        if (!bests[i].should_report()) {
          // read 1
          if (res_se1[i].best.valid_hit())
            align_read(res_se1[i].best, cigar1[i], reads1[i], pread, aln);

          if (res_se1[i].second_best.valid_hit())
            align_read(res_se1[i].second_best, tmp_cigar, reads1[i], pread, aln);

          if (res_se1[i].sort_candidates())
            cigar1[i] = tmp_cigar;

          // read 2
          if (res_se2[i].best.valid_hit())
            align_read(res_se2[i].best, cigar2[i], reads2[i], pread, aln);

          if (res_se2[i].second_best.valid_hit())
            align_read(res_se2[i].second_best, tmp_cigar, reads2[i], pread, aln);

          if (res_se2[i].sort_candidates())
            cigar2[i] = tmp_cigar;
        }
      }
    }

    pe_reported = false;
    se1_reported = false;
    se2_reported = false;
    for (size_t i = 0 ; i < n_reads; ++i) {
      select_output(allow_ambig, abismal_index.cl,
                    bests[i], res_se1[i], res_se2[i],
                    reads1[i], names1[i],
                    reads2[i], names2[i],
                    cigar1[i], cigar2[i],
                    pe_reported,
                    se1_reported, se2_reported, out);

      update_pe_stats(bests[i], res_se1[i], res_se2[i], reads1[i],
                      reads2[i], pe_reported, se1_reported, se2_reported,
                      pe_stats);
    }

    total_mapping_time += (omp_get_wtime() - start_time);
  }

  if (VERBOSE) {
    progress.report(cerr, get_filesize(reads_file1));
    cerr << "[total mapping time: " << total_mapping_time << "]" << endl;
  }
}

static void
write_sam_header(const ChromLookup &cl,
                 const int argc, const char **argv,
                 ostream &out) {

  static const string ABISMAL_VERSION = "0.1";

  out <<"@HD" << '\t'
      << "VN:1.0" << endl; // sam version
  // ADS: should the logic here be somehow hidden in ChromLookup?
  // ADS: subtracting below because of padding sequences
  const size_t n_chroms = cl.names.size() - 1;
  for (size_t i = 1; i < n_chroms; ++i) {
    out << "@SQ" << '\t'
        << "SN:" << cl.names[i] << '\t'
        << "LN:" << cl.starts[i+1] - cl.starts[i] << '\n';
  }

  // write how the abismal program was run
  out << "@PG" << '\t'
      << "ID:"
      << "ABISMAL" << '\t'
      << "VN:" << ABISMAL_VERSION << '\t'
      << "CL:";

  out << "\"" << argv[0];
  for (int i = 1; i < argc; ++i)
    out << " " << argv[i];
  out << "\"" << endl;
}

int main(int argc, const char **argv) {

  try {

    string index_file;
    string outfile;
    string stats_outfile = "";
    bool VERBOSE = false;
    bool GA_conversion = false;
    bool allow_ambig = false;
    bool pbat_mode = false;
    bool random_pbat = false;
    uint32_t max_candidates = 2000;
    size_t batch_size = 100000;
    size_t n_threads = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "map bisulfite converted reads",
                           "<reads-fq1> [<reads-fq2>]");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("index", 'i', "index file", true, index_file);
    opt_parse.add_opt("outfile", 'o', "output file", false, outfile);
    opt_parse.add_opt("mapstats", 'm',
                      "mapstats output file. If not provided, it "
                      "will be generated as .mapstats suffix to the "
                      "output file name", false, stats_outfile);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
    opt_parse.add_opt("shifts", 's', "number of seed shifts",
                      false, seed::n_shifts);
    opt_parse.add_opt("seed-pos", 'S', "seed length",
                      false, seed::n_seed_positions);
    opt_parse.add_opt("batch", 'b', "reads to load at once",
                      false, batch_size);
    opt_parse.add_opt("candidates", 'c', "max candidates for full comparison",
                      false, max_candidates);
    opt_parse.add_opt("max-mates", 'p', "max candidates as mates (pe mode)",
                      false, pe_candidates::max_size);
    opt_parse.add_opt("min-frag", 'l', "min fragment size (pe mode)",
                      false, pe_element::min_dist);
    opt_parse.add_opt("max-frag", 'L', "max fragment size (pe mode)",
                      false, pe_element::max_dist);
    opt_parse.add_opt("ambig", 'a', "report a posn for ambiguous mappers",
                      false, allow_ambig);
    opt_parse.add_opt("pbat", 'P', "input data follow the PBAT protocol",
                      false, pbat_mode);
    opt_parse.add_opt("random-pbat", 'R', "input data follow random PBAT",
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
    if (leftover_args.size() < 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (outfile.empty() && stats_outfile.empty()) {
      cerr << "please provide a file name for either the .sam"
           << "output (-o flag) or the stats file (-m flag)" << endl;

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

    if (VERBOSE)
      cerr << "[loading abismal index]" << endl;
    AbismalIndex abismal_index;
    const double start_time = omp_get_wtime();
    abismal_index.read(index_file);
    const double end_time = omp_get_wtime();
    if (VERBOSE)
      cerr << "[loading time: " << (end_time - start_time) << "]" << endl;

    if (VERBOSE)
      cerr << "[using " << n_threads << " threads for mapping]\n";

    if (VERBOSE) {
      if (paired_end)
        cerr << "[mapping paired end: "
             << reads_file << " " << reads_file2 << "]\n";
      else
        cerr << "[mapping single end: " << reads_file << "]\n";
    }
    if (VERBOSE)
      cerr << "[output file: " << outfile << "]" << endl;

    // avoiding opening the stats output file until mapping is done
    se_map_stats se_stats;
    pe_map_stats pe_stats(pe_element::min_dist, pe_element::max_dist);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    if (!out)
      throw runtime_error("failed to open output file: " + outfile);

    write_sam_header(abismal_index.cl, argc, argv, out);
    if (reads_file2.empty()) {
      if (GA_conversion || pbat_mode)
        map_single_ended<a_rich>(VERBOSE, allow_ambig, reads_file, batch_size,
                                 max_candidates, abismal_index, se_stats, out);
      else if (random_pbat)
        map_single_ended_rand(VERBOSE, allow_ambig, reads_file, batch_size,
                              max_candidates, abismal_index, se_stats, out);
      else
        map_single_ended<t_rich>(VERBOSE, allow_ambig, reads_file, batch_size,
                                 max_candidates, abismal_index,
                                 se_stats, out);
    }
    else {
      if (pbat_mode)
        map_paired_ended<a_rich>(VERBOSE, allow_ambig, reads_file, reads_file2,
                                 batch_size, max_candidates,
                                 abismal_index, pe_stats, out);
      else if (random_pbat)
        map_paired_ended_rand(VERBOSE,  allow_ambig, reads_file, reads_file2,
                              batch_size, max_candidates,
                              abismal_index, pe_stats, out);
      else
        map_paired_ended<t_rich>(VERBOSE, allow_ambig, reads_file, reads_file2,
                                 batch_size, max_candidates,
                                 abismal_index, pe_stats, out);
    }

    std::ofstream stat_out(stats_outfile.empty() ?
                           (outfile + ".mapstats") : stats_outfile);
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
