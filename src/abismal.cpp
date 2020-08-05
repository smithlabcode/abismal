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

#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"
#include "OptionParser.hpp"
#include "zlib_wrapper.hpp"
#include "GenomicRegion.hpp"
#include "dna_four_bit.hpp"
#include "sam_record.hpp"
#include "bisulfite_utils.hpp"

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
      if (line_count % 4 == 0)
        names.push_back(line.substr(1, line.find_first_of(" \t")));
      else if (line_count % 4 == 1) {
        if (count_if(begin(line), end(line),
                     [](const char c) {return c != 'N';}) < min_read_length)
          line.clear();
        std::replace(begin(line), end(line), 'N', 'Z');
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
  score_t aln_score;
  flags_t flags;

  se_element() :
    pos(0), diffs(invalid_hit_diffs), aln_score(0), flags(0) {}

  se_element(const uint32_t p, const score_t d, const score_t scr,
             const flags_t f) :
    pos(p), diffs(d), aln_score(scr), flags(f) {}

  bool operator==(const se_element &rhs) const {
    return diffs == rhs.diffs && pos == rhs.pos;
  }

  // this is used to keep PE candidates sorted in the max heap
  bool operator<(const se_element &rhs) const {
    return diffs < rhs.diffs;
  }

  bool is_better_hit_than (const se_element &rhs) const {
    return diffs < rhs.diffs;
  }

  bool is_better_aln_than (const se_element &rhs) const {
    return aln_score > rhs.aln_score;
  }

  bool rc() const {return samflags::check(flags, samflags::read_rc);}
  bool elem_is_a_rich() const {
    return !samflags::check(flags, bsflags::read_is_t_rich);
  }
  bool valid_hit() const {return diffs < invalid_hit_diffs;}
  void reset() {diffs = invalid_hit_diffs; aln_score = 0;}

  bool is_equal_to(const se_element &rhs) const {
    return (diffs == rhs.diffs) && (pos != rhs.pos);
  }
  static score_t invalid_hit_diffs;
  static uint16_t min_aligned_length;
};

score_t se_element::invalid_hit_diffs = 40;
uint16_t se_element::min_aligned_length = ReadLoader::min_read_length;

struct se_result {
  se_element best;
  se_element second_best;
  se_result() : best(se_element()), second_best(se_element()) {}

  void update_by_mismatch(const uint32_t p, const score_t d, const flags_t s) {
    // avoid having two copies of the best hit
    if (p == best.pos && s == best.flags) return;
    const se_element cand(p, d, 0, s); // 0 = no alignment performed
    if (cand.is_better_hit_than(second_best)) second_best = cand;
    if (second_best.is_better_hit_than(best)) std::swap(best, second_best);
  }

  bool sort_by_score() {
    if (second_best.is_better_aln_than(best)) {
      std::swap(best, second_best);
      return true;
    }
    return false;
  }

  uint8_t mapq() const {
    if (!second_best.valid_hit()) return unknown_mapq_score;
    return max_mapq_score*(best.aln_score - second_best.aln_score) /
      best.aln_score;
  }

  bool ambig() const {
    return best.aln_score == second_best.aln_score;
  }

  bool ambig_diffs() const {
    return best.diffs == second_best.diffs;
  }

  bool sure_ambig(uint32_t seed_number = 0) const {
    return ambig_diffs() &&
      (best.diffs == 0 || (best.diffs == 1 && seed_number > 0));
  }

  bool should_report(const bool allow_ambig) const {
    return (allow_ambig || !ambig())  // unique mapping
      &&  best.valid_hit(); // concordant pair
  }

  void reset() {best.reset(); second_best.reset();}
  uint32_t get_cutoff() const {return second_best.diffs;}

  static uint8_t max_mapq_score;
  static uint8_t unknown_mapq_score;
  static uint8_t min_mapq_score;
};

uint8_t se_result::min_mapq_score = 1;
uint8_t se_result::max_mapq_score = 250;
uint8_t se_result::unknown_mapq_score = 255;

inline bool
chrom_and_posn(const ChromLookup &cl, const string &cig, const uint32_t p,
               uint32_t &r_p, uint32_t &r_e, uint32_t &r_chr) {
  const uint32_t ref_ops = cigar_rseq_ops(cig);
  if (!cl.get_chrom_idx_and_offset(p, ref_ops, r_chr, r_p)) return false;
  r_e = r_p + ref_ops;
  return true;
}

static void
format_se(const bool allow_ambig, se_result res, const ChromLookup &cl,
          string &read, const string &read_name,
          string &cigar, ofstream &out) {

  uint32_t ref_s = 0, ref_e = 0, chrom_idx = 0;

  se_element s = res.best;
  if (res.should_report(allow_ambig) &&
      chrom_and_posn(cl, cigar, s.pos, ref_s, ref_e, chrom_idx)) {
    sam_rec sr(read_name, 0, cl.names[chrom_idx], ref_s + 1,
               res.mapq(), cigar, "*", 0, 0, read, "*");
    if (s.rc())
      set_flag(sr, samflags::read_rc);
    sr.add_tag("NM:i:" + to_string(s.diffs));
    sr.add_tag(s.elem_is_a_rich() ? "CV:A:A" : "CV:A:T");
    out << sr << '\n';
  }
}

struct pe_element {
  se_element r1;
  se_element r2;

  pe_element() : r1(se_element()), r2(se_element()) {}
  pe_element(const se_element &s1, const se_element &s2) : r1(s1), r2(s2) {}

  bool rc() const { return r1.rc(); }
  bool elem_is_a_rich() const {return r1.elem_is_a_rich();}
  score_t diffs() const { return r1.diffs + r2.diffs; }
  score_t score() const { return r1.aln_score + r2.aln_score; }

  bool valid_hit() const {
    return r1.diffs < se_element::invalid_hit_diffs &&
                      r2.diffs < se_element::invalid_hit_diffs;
  }

  bool is_equal_to(const pe_element &rhs) const {
    return diffs() == rhs.diffs() && !(r1.pos == rhs.r1.pos &&
                                       r2.pos == rhs.r2.pos);
  }

  bool is_better_hit_than(const pe_element &rhs) const {
    return diffs() < rhs.diffs();
  }

  bool is_better_aln_than(const pe_element &rhs) const {
    return (score() > rhs.score());
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

  void reset() {
    best.reset();
    second_best.reset();
  }

  bool update_by_score(const pe_element &p) {
    if (p.is_better_aln_than(second_best)) second_best = p;
    if (second_best.is_better_aln_than(best)) {
      std::swap(best, second_best);
      return true; //best has been updated
    }
    return false;
  }

  bool ambig() const {
    return (se_result::max_mapq_score*(best.score() - second_best.score()) <
            best.score()*se_result::min_mapq_score);
  }

  uint8_t mapq() const {
    if (!second_best.valid_hit()) return se_result::unknown_mapq_score;
    return se_result::max_mapq_score*(best.score() - second_best.score()) /
      best.score();
  }

  bool should_report(const bool allow_ambig) const {
    return (allow_ambig || !ambig())  // unique mapping
      &&  best.valid_hit(); // concordant pair
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
          const string &cig1, const string &cig2, ofstream &out) {

  uint32_t r_s1 = 0, r_e1 = 0, chr1 = 0; // positions in ref (0-based)
  uint32_t r_s2 = 0, r_e2 = 0, chr2 = 0;
  const pe_element p = res.best;

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

  // ADS: does it make sense to use this? And for both ends?
  const uint8_t mapq = res.mapq();

  // ADS: +1 to POS & PNEXT; "=" for RNEXT; "*" for QUAL
  sam_rec sr1(name1, 0, cl.names[chr1], r_s1 + 1, mapq,
              cig1, "=", r_s2 + 1, tlen, read1, "*");
  sr1.add_tag("NM:i:" + to_string(p.r1.diffs));
  sr1.add_tag(p.r1.elem_is_a_rich() ? "CV:A:A" : "CV:A:T");
  set_flag(sr1, samflags::read_paired);
  set_flag(sr1, samflags::read_pair_mapped);
  set_flag(sr1, samflags::template_first);

  sam_rec sr2(name2, 0, cl.names[chr2], r_s2 + 1, mapq,
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
  void update_by_mismatch(const uint32_t p, const score_t d, const flags_t s) {
    if (full()) {
      if (d < v.front().diffs) {
        std::pop_heap(begin(v), end(v));
        v.back() = se_element(p, d, 0, s);
        std::push_heap(begin(v), end(v));
      }
    }
    else if (d < se_element::invalid_hit_diffs) {
      v[sz++] = se_element(p, d, 0, s);
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

  void update(const string &read, const se_result &res) {
    ++tot_rds;
    if (res.best.valid_hit()) {
      if (!res.ambig()) ++uniq_rds;
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
        << t+tab << "percent_mapped: "
        << pct(uniq_rds+ambig_rds, tot_rds == 0 ? 1 : tot_rds) << endl
        << t+tab << "unique: " << uniq_rds << endl
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

  void update_pair(const pe_result &res) {
    ++tot_pairs;
    if (res.best.valid_hit()) {
      const bool ambig = res.ambig();
      ambig_pairs += ambig;
      uniq_pairs += !ambig;
    }
    else ++unmapped_pairs;
  }

  string tostring() const {
    std::ostringstream oss;
    static const string t = "    ";
    oss << "pairs:" << endl
        << t   << "total_read_pairs: " << tot_pairs << endl
        << t   << "mapped:" << endl
        << t+t << "percent_mapped: "
        << pct(uniq_pairs + ambig_pairs, tot_pairs) << endl
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
                const se_result &se1, const se_result &se2,
                const string &read1, const string &read2,
                pe_map_stats &pe_stats) {
  pe_stats.update_pair(best);
  if (!best.should_report(false)) {
    pe_stats.end1_stats.update(read1, se1);
    pe_stats.end2_stats.update(read2, se2);
  }
}


static void
select_output(const bool allow_ambig, const ChromLookup &cl,
              pe_result &best, se_result &se1, se_result &se2,
              string &read1, const string &name1,
              string &read2, const string &name2,
              string &cig1, string &cig2, ofstream &out) {

  if (best.should_report(false)) {
    if (!format_pe(allow_ambig, best, cl, read1, read2, name1, name2,
                   cig1, cig2, out)) {
      // if unable to fetch chromosome positions (e.g. due to read
      // mapping in between chromosomes or cigars breaking dovetail
      // reads), consider it unmapped
      best.reset();
      se1.reset();
      se2.reset();
    }
  }
  else {
    format_se(allow_ambig, se1, cl, read1, name1, cig1, out);
    format_se(allow_ambig, se2, cl, read2, name2, cig2, out);
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
  for (auto i(start_idx); i != end_idx && !res.sure_ambig(offset != 0); ++i) {
    const uint32_t pos = (*i) - offset;
    const score_t diffs =
      ((pos & 1) ?
       full_compare<pos_odd>(res.get_cutoff(),
                             odd_read_st, odd_read_mid, odd_read_end,
                             genome_st + (pos >> 1)) :
       full_compare<pos_even>(res.get_cutoff(),
                              even_read_st, even_read_mid, even_read_end,
                              genome_st + (pos >> 1)));
    res.update_by_mismatch(pos, diffs, strand_code);
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
                const uint32_t n_sorting_positions,
                vector<uint32_t>::const_iterator &low,
                vector<uint32_t>::const_iterator &high) {
  size_t p = seed::key_weight;
  const size_t lim = std::min(read_lim, n_sorting_positions);
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

  const auto index_st(begin(abismal_index.index));
  const auto counter_st(begin(abismal_index.counter));

  // ADS: there is something wrong with this code if so many
  // conditions are needed
  const size_t shift_lim =
    readlen > (seed::n_seed_positions + 1) ?
    (readlen - seed::n_seed_positions - 1) : 0;
  const size_t shift = (seed::n_shifts == 1) ? shift_lim :
    std::max(1ul, shift_lim/(seed::n_shifts - 1));

  bool found_good_seed = false;

  // uniformly spaced seeds
  for (uint32_t i = 0; i <= shift_lim && !res.sure_ambig(i); i += shift) {
    const uint32_t offset = i;
    uint32_t k = 0;
    get_1bit_hash(read_start + offset, k);
    auto s_idx(index_st + *(counter_st + k));
    auto e_idx(index_st + *(counter_st + k + 1));

    if (s_idx < e_idx) {
      find_candidates(read_start + offset, gi, readlen - offset,
                      seed::n_seed_positions, s_idx, e_idx);

      if (e_idx - s_idx < max_candidates) {
        found_good_seed = true;
        check_hits<strand_code>(s_idx, e_idx,
                                even_read_start, even_read_mid, even_read_end,
                                odd_read_start, odd_read_mid, odd_read_end,
                                genome_st, offset, res);
      }
    }
  }

  // if no good seeds found, use entire read as seed
  if (!found_good_seed) {
    uint32_t k = 0;
    get_1bit_hash(read_start, k);
    auto s_idx(index_st + *(counter_st + k));
    auto e_idx(index_st + *(counter_st + k + 1));
    if (s_idx < e_idx) {
      const uint32_t posns_to_compare = min(readlen, seed::n_sorting_positions);
      // ADS: is above "min" needed, given check inside find_candidates?
      find_candidates(read_start, gi, readlen, posns_to_compare, s_idx, e_idx);

      if (e_idx - s_idx < max_candidates)
        check_hits<strand_code>(s_idx, e_idx,
                                even_read_start, even_read_mid, even_read_end,
                                odd_read_start, odd_read_mid, odd_read_end,
                                genome_st, 0, res);
    }
  }
}


template <const bool convert_a_to_g>
static void
prep_read(const string &r, Read &pread) {
  pread.resize(r.size());
  for (size_t i = 0; i < r.size(); ++i)
    pread[i] = encode_dna_four_bit(convert_a_to_g ?
                                   (r[i] == 'A' ? 'R' : r[i]) :
                                   (r[i] == 'T' ? 'Y' : r[i]));
}

// creates reads meant for comparison on a compressed genome
static void
prep_for_seeds(const Read &pread_seed, Read &pread_even, Read &pread_odd) {
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

namespace local_aln {
  static const score_t match = 1;
  static const score_t mismatch = -1;
  static const score_t indel = -1;

  // this lookup improves speed when running alignment
  static const vector<score_t> score_lookup = {match, mismatch};
  static inline score_t mismatch_score(const char q_base,
                                       const uint8_t t_base) {
    return score_lookup[the_comp(q_base, t_base)];
  }
};

using AbismalAlignSimple = AbismalAlign<local_aln::mismatch_score,
                                        local_aln::indel>;

static void
align_read(se_element &res, string &cigar, const string &read,
           Read &pread, AbismalAlignSimple &aln) {
  // ends early if alignment is nearly diagonal
  if (res.diffs <= 3) {
    cigar = std::to_string(read.length()) + "M";
    // ADS: aln_score
    res.aln_score = local_aln::match*(read.length() - res.diffs) +
      local_aln::mismatch*res.diffs;
  }
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
    res.aln_score = aln.align(pread, res.pos, len, cigar);
    // will be used to print edit distance on sam. Only works with 1 -1 -1
    res.diffs = (len - res.aln_score) / 2;
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
                 ofstream &out) {

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

        if (res[i].sort_by_score())
          cigar[i] = tmp_cigar;
      }
    }

    for (size_t i = 0 ; i < n_reads; ++i) {
      se_stats.update(reads[i], res[i]);
      format_se(allow_ambig, res[i], abismal_index.cl,
                reads[i], names[i], cigar[i], out);
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
                      ofstream &out) {
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

        if (res[i].sort_by_score())
          cigar[i] = tmp_cigar;
      }
    }

    for (size_t i = 0 ; i < n_reads; ++i) {
      se_stats.update(reads[i], res[i]);
      format_se(allow_ambig, res[i], abismal_index.cl,
                reads[i], names[i], cigar[i], out);
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
    res.update_by_mismatch(i->pos, i->diffs, i->flags);
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
            if (best.update_by_score(p)) {
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
void
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
                 ofstream &out) {
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
        select_maps<false>(reads1[i], reads2[i],
                           cigar1[i], cigar2[i],
                           res1[i], res2[i],
                           res_se1[i], res_se2[i],
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
        select_maps<true>(reads2[i], reads1[i],
                          cigar2[i], cigar1[i],
                          res2[i], res1[i],
                          res_se2[i], res_se1[i],
                          aln, bests[i]);
    }

#pragma omp parallel
    {
      Read pread;
      string tmp_cigar;
      AbismalAlignSimple aln(gi, genome_size, max_batch_read_length);

#pragma omp for
      for (size_t i = 0; i < n_reads; ++i) {
        if (!bests[i].should_report(false)) {
          if (res_se1[i].best.valid_hit())
            align_read(res_se1[i].best, cigar1[i], reads1[i], pread, aln);

          if (res_se1[i].second_best.valid_hit())
            align_read(res_se1[i].second_best, tmp_cigar, reads1[i], pread, aln);

          if (res_se1[i].sort_by_score())
            cigar1[i] = tmp_cigar;

          if (res_se2[i].best.valid_hit())
            align_read(res_se2[i].best, cigar2[i], reads2[i], pread, aln);

          if (res_se2[i].second_best.valid_hit())
            align_read(res_se2[i].second_best, tmp_cigar, reads2[i], pread, aln);

          if (res_se2[i].sort_by_score())
            cigar2[i] = tmp_cigar;
        }
      }
    }

    for (size_t i = 0 ; i < n_reads; ++i)
      select_output(allow_ambig, abismal_index.cl,
                    bests[i], res_se1[i], res_se2[i],
                    reads1[i], names1[i],
                    reads2[i], names2[i],
                    cigar1[i], cigar2[i], out);

    for (size_t i = 0 ; i < n_reads; ++i)
      update_pe_stats(bests[i], res_se1[i], res_se2[i], reads1[i],
                      reads2[i], pe_stats);
    total_mapping_time += (omp_get_wtime() - start_time);
  }

  if (VERBOSE) {
    progress.report(cerr, get_filesize(reads_file1));
    cerr << "[total mapping time: " << total_mapping_time << "]" << endl;
  }
}

void
map_paired_ended_rand(const bool VERBOSE,
                      const bool allow_ambig,
                      const string &reads_file1,
                      const string &reads_file2,
                      const size_t batch_size,
                      const size_t max_candidates,
                      const AbismalIndex &abismal_index,
                      pe_map_stats &pe_stats,
                      ofstream &out) {
  double total_mapping_time = 0;

  ReadLoader rl1(reads_file1, batch_size);
  ReadLoader rl2(reads_file2, batch_size);

  size_t max_batch_read_length;
  vector<string> names1(batch_size), reads1(batch_size), cigar1(batch_size);
  vector<string> names2(batch_size), reads2(batch_size), cigar2(batch_size);

  vector<pe_candidates> res1(batch_size), res2(batch_size);
  vector<pe_result> bests(batch_size);
  vector<se_result> res_se1(batch_size), res_se2(batch_size);

  // alignment stuff
  const uint32_t genome_size = abismal_index.cl.get_genome_size();
  const Genome::const_iterator genome_st(begin(abismal_index.genome));
  const genome_iterator gi(genome_st);
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
        if (!bests[i].should_report(false)) {
          // read 1
          if (res_se1[i].best.valid_hit())
            align_read(res_se1[i].best, cigar1[i], reads1[i], pread, aln);

          if (res_se1[i].second_best.valid_hit())
            align_read(res_se1[i].second_best, tmp_cigar, reads1[i], pread, aln);

          if (res_se1[i].sort_by_score())
            cigar1[i] = tmp_cigar;

          // read 2
          if (res_se2[i].best.valid_hit())
            align_read(res_se2[i].best, cigar2[i], reads2[i], pread, aln);

          if (res_se2[i].second_best.valid_hit())
            align_read(res_se2[i].second_best, tmp_cigar, reads2[i], pread, aln);

          if (res_se2[i].sort_by_score())
            cigar2[i] = tmp_cigar;
        }
      }
    }

    for (size_t i = 0 ; i < n_reads; ++i)
      select_output(allow_ambig, abismal_index.cl,
                    bests[i], res_se1[i], res_se2[i],
                    reads1[i], names1[i],
                    reads2[i], names2[i],
                    cigar1[i], cigar2[i], out);

    for (size_t i = 0 ; i < n_reads; ++i)
      update_pe_stats(bests[i], res_se1[i], res_se2[i], reads1[i],
                      reads2[i], pe_stats);

    total_mapping_time += (omp_get_wtime() - start_time);
  }

  if (VERBOSE) {
    progress.report(cerr, get_filesize(reads_file1));
    cerr << "[total mapping time: " << total_mapping_time << endl;
  }
}

void
write_sam_header(const ChromLookup &cl,
                 const int argc,
                 const char **argv,
                 ofstream &out) {
  out <<"@HD\tVN:1.0\n"; // sam version
  const size_t sz = cl.names.size();
  size_t prev = cl.starts[0];

  // sequence lengths
  for (size_t i = 1; i < sz; ++i) {
    out << "@SQ\tSN:" << cl.names[i] << "\tLN:" << cl.starts[i] - prev << '\n';
    prev = cl.starts[i];
  }

  // function call
  out << "@PG\tID:ABISMAL\tVN:1.0.0\tCL:\"";
  out << string(*(argv));
  for (int i = 1; i < argc; ++i)
    out << " " << string(*(argv + i));
  out << "\"\n";

}

int main(int argc, const char **argv) {

  try {

    string index_file;
    string outfile;
    bool VERBOSE = false;
    bool GA_conversion = false;
    bool allow_ambig = false;
    bool pbat_mode = false;
    bool random_pbat = false;
    uint32_t max_candidates = 2000;
    size_t batch_size = 100000;
    size_t n_threads = omp_get_max_threads();

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "map bisulfite converted reads",
                           "<reads-fq1> [<reads-fq2>]");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("index", 'i', "index file", true, index_file);
    opt_parse.add_opt("outfile", 'o', "output file", true, outfile);
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
    uint32_t index_max_cand = 0;
    abismal_index.read(index_file, seed::n_sorting_positions, index_max_cand);
    if (VERBOSE)
      // ADS: this should be removed soon; not helpful to a user
      cerr << "[index: n_sorting_positions = "
           << seed::n_sorting_positions << "]" << endl
           << "[index: max_candidates = "
           << seed::n_sorting_positions << "]" << endl;

    const double end_time = omp_get_wtime();
    if (VERBOSE)
      cerr << "[loading time: " << (end_time - start_time) << "]" << endl;

    if (VERBOSE)
      cerr << "[using " << n_threads << " threads for mapping]\n";

    if (VERBOSE) {
      if (paired_end)
        cerr << "[mapping paired end: " << reads_file << " "
             << reads_file2 << "]\n";
      else
        cerr << "[mapping single end: " << reads_file << "]\n";
    }
    if (VERBOSE)
      cerr << "[output file: " << outfile << "]" << endl;

    // avoiding opening the stats output file until mapping is done
    se_map_stats se_stats;
    pe_map_stats pe_stats(pe_element::min_dist, pe_element::max_dist);

    std::ofstream out(outfile);
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
