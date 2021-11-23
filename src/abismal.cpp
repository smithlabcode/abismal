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
#include <chrono>
#include <numeric>

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

using std::__gcd;
using std::vector;
using std::runtime_error;
using std::string;
using std::cerr;
using std::endl;
using std::cout;
using std::numeric_limits;
using std::ostream;
using std::ofstream;
using std::max;
using std::min;
using std::to_string;
using std::begin;
using std::end;
using std::deque;
using std::push_heap;
using std::pop_heap;
using std::chrono::system_clock;
using std::ostringstream;

typedef uint16_t flags_t; // every bit is a flag
typedef int16_t score_t; // aln score, edit distance, hamming distance
typedef vector<uint8_t> Read; //4-bit encoding of reads
typedef vector<size_t> PackedRead; //4-bit encoding of reads

enum conversion_type { t_rich = false, a_rich = true };

static inline void
print_with_time(const string &s) {
  auto tmp = system_clock::to_time_t(system_clock::now());
  string time_fmt(std::ctime(&tmp));
  time_fmt.pop_back();
  cerr << "[" << time_fmt << "] " << s << endl;
}

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
  ReadLoader(const string &fn) :
    filename(fn) {
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
        if (line.empty())
          throw runtime_error("malformatted FASTQ file "
                              "contains an empty read name\n");
        names.push_back(line.substr(1, line.find_first_of(" \t") - 1));
      }
      else if (line_count % 4 == 1) {
        // read too long, may pass the end of the genome
        if (line.size() >= seed::padding_size) {
          throw runtime_error("found a read of size " + to_string(line.size()) +
              ", which is too long. Maximum allowed read size = " +
              to_string(seed::padding_size));
        }
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
  igzfstream *in;

  static const size_t batch_size;
  static const uint32_t min_read_length;
};

const size_t ReadLoader::batch_size = 2000;
const uint32_t ReadLoader::min_read_length =
  seed::key_weight + seed::window_size - 1;

inline void
update_max_read_length(size_t &max_length, const vector<string> &reads) {
  for (vector<string>::const_iterator it (begin(reads)); it != end(reads); ++it)
    max_length = max(max_length, it->size());
}

struct se_element { //assert(sizeof(se_element) == 8)
  uint32_t pos;
  score_t diffs;
  flags_t flags;

  se_element() :
    pos(0), diffs(numeric_limits<score_t>::max() - 1),  flags(0) {}

  se_element(const uint32_t p, const score_t d, const flags_t f) :
    pos(p), diffs(d), flags(f) {}

  bool operator==(const se_element &rhs) const {
    return pos == rhs.pos && flags == rhs.flags;
  }
  bool operator!=(const se_element &rhs) const {
    return pos != rhs.pos || flags != rhs.flags;
  }

  // this is used to keep PE candidates sorted in the max heap
  bool operator<(const se_element &rhs) const {
    return diffs < rhs.diffs;
  }

  inline bool rc() const {
    return samflags::check(flags, samflags::read_rc);
  }
  inline bool elem_is_a_rich() const {
    return samflags::check(flags, bsflags::read_is_a_rich);
  }
  inline bool ambig() const {
    return samflags::check(flags, samflags::secondary_aln);
  }
  inline void set_ambig() {
    samflags::set(flags, samflags::secondary_aln);
  }
  inline void unset_ambig() {
    samflags::unset(flags, samflags::secondary_aln);
  }
  inline bool empty() const { return pos == 0; }
  inline bool sure_ambig() const { return ambig() && diffs == 0; }

  inline void reset() {
    pos = 0;
    unset_ambig();
  }

  inline void reset(const uint32_t readlen) {
    reset();
    diffs = static_cast<score_t>(invalid_hit_frac*readlen);
  }
  static double valid_frac;
  static const double invalid_hit_frac;
};

double se_element::valid_frac = 0.1;
const double se_element::invalid_hit_frac = 0.4;

inline bool
valid(const se_element &s, const uint32_t readlen) {
  return s.diffs <= static_cast<score_t>(se_element::valid_frac*readlen);
}

inline bool
valid_pair(const se_element &s1, const se_element &s2,
           const uint32_t readlen1, const uint32_t readlen2) {
  const score_t d = s1.diffs + s2.diffs;
  const uint32_t len = readlen1 + readlen2;
  return d <= static_cast<score_t>(se_element::valid_frac*len);
}

inline bool
valid_hit(const se_element s, const uint32_t readlen) {
  return s.diffs < static_cast<score_t>(se_element::invalid_hit_frac*readlen);
}

struct se_candidates {
  se_candidates () :
    sz(1), max_size(min_size), best(se_element()),
    v(vector<se_element>(expand_size)) {}
  inline bool full() const { return sz == max_size; };
  void update_exact_match(const uint32_t p, const score_t d, const flags_t s) {
    const se_element cand(p, d, s);
    if (best.empty())
      best = cand; // cand has ambig flag set to false

    else if (cand != best)
      best.set_ambig();
  }

  void update_cand(const uint32_t p, const score_t d, const flags_t s) {
    if (full()) {
      pop_heap(begin(v), begin(v) + sz);
      v[sz - 1] = se_element(p, d, s);
    }
    else {
      v[sz++] = se_element(p, d, s);
    }
    push_heap(begin(v), begin(v) + sz);
  }

  void update(const uint32_t p, const score_t d, const flags_t s) {
    if (d == 0)
      update_exact_match(p, d, s);
    update_cand(p, d, s);
  }

  inline void expand() {
    max_size = expand_size;
  }

  inline bool sure_ambig() const {
    return best.sure_ambig();
  }

  void reset(const uint32_t readlen)  {
    v.front().reset(readlen);
    sz = 1;
    best.reset(readlen);
    max_size = min_size;
  }

  void prepare_for_alignments() {
    sort(begin(v), begin(v) + sz, // no sort_heap here as heapify used "diffs"
         [](const se_element &a, const se_element &b) {
           return (a.pos < b.pos) || (a.pos == b.pos && a.flags < b.flags);
         });
    sz = unique(begin(v), begin(v) + sz) - begin(v);
  }

  score_t get_cutoff() const { return v.front().diffs; }

  uint32_t sz;
  uint32_t max_size;
  se_element best;
  vector<se_element> v;

  static const uint32_t min_size = 10u;
  static const uint32_t expand_size = 100u;
};

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
format_se(const bool allow_ambig, const se_element &res, const ChromLookup &cl,
          const string &read, const string &read_name, const string &cigar,
          ostringstream &out) {
  const bool ambig = res.ambig();
  const bool valid = !res.empty();
  if (!allow_ambig && ambig)
    return map_ambig;

  uint32_t ref_s = 0, ref_e = 0, chrom_idx = 0;
  if (!valid || !chrom_and_posn(cl, cigar, res.pos, ref_s, ref_e, chrom_idx))
    return map_unmapped;

  sam_rec sr(read_name, 0, cl.names[chrom_idx], ref_s + 1,
             255, cigar, "*", 0, 0, read, "*");
  if (res.rc())
    set_flag(sr, samflags::read_rc);

  if (allow_ambig && ambig)
    set_flag(sr, samflags::secondary_aln);

  sr.add_tag("NM:i:" + to_string(res.diffs));
  sr.add_tag(res.elem_is_a_rich() ? "CV:A:A" : "CV:A:T");

  out << sr.tostring() << "\n";
  return ambig ? map_ambig : map_unique;
}

struct pe_element {
  pe_element() :
    aln_score(0), r1(se_element()), r2(se_element()) {}

  score_t diffs() const { return r1.diffs + r2.diffs; }
  void reset(const uint32_t readlen1, const uint32_t readlen2) {
    aln_score = 0;
    r1.reset(readlen1);
    r2.reset(readlen2);
    max_aln_score = simple_aln::best_pair_score(readlen1, readlen2);
  }
  void reset() {
    aln_score = 0;
    r1.reset();
    r2.reset();
  }

  bool update(const score_t scr,
      const se_element s1, const se_element s2) {

    if (scr > aln_score) {
      r1 = s1;
      r2 = s2;
      aln_score = scr;
      return true;
    }
    else if (scr == aln_score)
      r1.set_ambig();

    return false;
  }

  inline bool should_report(const bool allow_ambig) const {
    return !empty() && (allow_ambig || !ambig());
  }
  inline bool ambig() const { return r1.ambig(); }
  inline bool empty() const { return r1.empty(); }
  inline bool sure_ambig() const {
    return ambig() && (aln_score == max_aln_score);
  }

  score_t aln_score;
  score_t max_aln_score;
  se_element r1;
  se_element r2;

  static uint32_t min_dist;
  static uint32_t max_dist;
};

uint32_t pe_element::min_dist = 32;
uint32_t pe_element::max_dist = 3000;

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
          const pe_element &p, const ChromLookup &cl,
          const string &read1, const string &read2,
          const string &name1, const string &name2,
          const string &cig1,  const string &cig2,
          ostringstream &out) {
  const bool ambig = p.ambig();
  const bool valid = !p.empty();
  if (!allow_ambig && ambig)
    return map_ambig;

  uint32_t r_s1 = 0, r_e1 = 0, chr1 = 0; // positions in chroms (0-based)
  uint32_t r_s2 = 0, r_e2 = 0, chr2 = 0;

  // PE chromosomes differ or couldn't be found, treat read as unmapped
  if (!valid || !chrom_and_posn(cl, cig1, p.r1.pos, r_s1, r_e1, chr1) ||
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

  if (allow_ambig && ambig) {
    set_flag(sr1, samflags::secondary_aln);
    set_flag(sr2, samflags::secondary_aln);
  }

  out << sr1.tostring() << "\n" << sr2.tostring() << "\n";

  return ambig ? map_ambig : map_unique;
}

struct pe_candidates {
  pe_candidates(const uint32_t s) :
    sz(1), min_size (s),
    max_size(min_size),
    v(vector<se_element>(max(s, pe_element::max_dist - pe_element::min_dist + 1))) {}
  bool full() const { return sz == max_size; }

  inline void reset(const uint32_t readlen) {
    v.front().reset(readlen);
    sz = 1;
    max_size = min_size;
  }

  score_t get_cutoff() const { return v.front().diffs; }

  void update(const uint32_t p, const score_t d, const flags_t s) {
    if (full()) {
      pop_heap(begin(v), begin(v) + sz);
      v[sz - 1] = se_element(p, d, s);
    }
    else {
      v[sz++] = se_element(p, d, s);
    }
    push_heap(begin(v), begin(v) + sz);
  }

  inline bool sure_ambig() const {
    return full() && (v.front().diffs == 0);
  }

  inline void expand() {
    max_size = pe_element::max_dist - pe_element::min_dist + 1;
  }

  void prepare_for_mating() {
    sort(begin(v), begin(v) + sz, // no sort_heap here as heapify used "diffs"
         [](const se_element &a, const se_element &b) {return a.pos < b.pos;});
    sz = unique(begin(v), begin(v) + sz) - begin(v);
  }

  uint32_t sz;
  const uint32_t min_size;
  uint32_t max_size;
  vector<se_element> v;
};

inline double pct(const double a, const double b) {return ((b == 0) ? 0.0 : 100.0*a/b);}

struct se_map_stats {
  se_map_stats() :
    tot_rds(0), uniq_rds(0), ambig_rds(0), unmapped_rds(0), skipped_rds(0),
    edit_distance(0), total_bases(0) {}
  uint32_t tot_rds;
  uint32_t uniq_rds;
  uint32_t ambig_rds;
  uint32_t unmapped_rds;
  uint32_t skipped_rds;

  size_t edit_distance;
  size_t total_bases;

  void update(const bool allow_ambig,
              const string &read, const string &cigar,
              const se_element s) {
    ++tot_rds;
    const bool valid = !s.empty();
    const bool ambig = s.ambig();
    uniq_rds += (valid && !ambig);
    ambig_rds += (valid && ambig);
    unmapped_rds += !valid;
    skipped_rds += read.empty();

    if (valid && (allow_ambig || !ambig))
      update_error_rate(s.diffs, cigar);
  }

  void update_error_rate(const score_t diffs, const string &cigar) {
    edit_distance += diffs;
    total_bases += cigar_rseq_ops(cigar);
  }

  string tostring(const size_t n_tabs = 0) const {
    static const string tab = "    ";
    string t;
    for (size_t i = 0; i < n_tabs; ++i) t += tab;
    ostringstream oss;

    oss << t     << "total_reads: " << tot_rds << endl
        << t     << "mapped: " << endl
        << t+tab << "num_mapped: " << uniq_rds+ambig_rds << endl
        << t+tab << "num_unique: " << uniq_rds << endl
        << t+tab << "num_ambiguous: " << ambig_rds << endl
        << t+tab << "percent_mapped: "
        << pct(uniq_rds+ambig_rds, tot_rds == 0 ? 1 : tot_rds) << endl
        << t+tab << "percent_unique: "
        << pct(uniq_rds, tot_rds == 0 ? 1 : tot_rds) << endl
        << t+tab << "percent_ambiguous: " << pct(ambig_rds, tot_rds) << endl
        << t+tab << "unique_error:" << endl
        << t+tab+tab << "edits: " << edit_distance << endl
        << t+tab+tab << "total_bases: " << total_bases << endl
        << t+tab+tab << "error_rate: " << pct(edit_distance, total_bases) << endl
        << t     << "num_unmapped: " << unmapped_rds << endl
        << t     << "num_skipped: " << skipped_rds << endl
        << t     << "percent_unmapped: " << pct(unmapped_rds, tot_rds) << endl
        << t     << "percent_skipped: " << pct(skipped_rds, tot_rds) << endl;
    return oss.str();
  }
};

struct pe_map_stats {
  pe_map_stats() :
    tot_pairs(0), uniq_pairs(0), ambig_pairs(0),
    unmapped_pairs(0), skipped_pairs(0), edit_distance(0), total_bases(0) {}
  uint32_t tot_pairs;
  uint32_t uniq_pairs;
  uint32_t ambig_pairs;
  uint32_t unmapped_pairs;
  uint32_t skipped_pairs;

  size_t edit_distance;
  size_t total_bases;

  se_map_stats end1_stats;
  se_map_stats end2_stats;

  void update(const bool allow_ambig,
              const string &reads1, const string &reads2,
              const string &cig1, const string &cig2,
              const pe_element &p,
              const se_element s1, const se_element s2) {
    const bool valid = !p.empty();
    const bool ambig = p.ambig();
    ++tot_pairs;
    ambig_pairs += (valid && ambig);
    uniq_pairs += (valid && !ambig);
    unmapped_pairs += !valid;
    skipped_pairs += (reads1.empty() || reads2.empty());

    if (p.should_report(allow_ambig)) {
      update_error_rate(p.r1.diffs, p.r2.diffs, cig1, cig2);
    }
    else {
      end1_stats.update(false, reads1, cig1, s1);
      end2_stats.update(false, reads1, cig2, s2);
    }
  }

  void update_error_rate(const score_t d1, const score_t d2,
                         const string &cig1, const string &cig2) {
    edit_distance += d1 + d2;
    total_bases += cigar_rseq_ops(cig1) + cigar_rseq_ops(cig2);
  }

  string tostring(const bool allow_ambig) const {
    ostringstream oss;
    static const string t = "    ";
    oss << "pairs:" << endl
        << t   << "total_pairs: " << tot_pairs << endl
        << t   << "mapped:" << endl
        << t+t << "num_mapped: " << uniq_pairs + ambig_pairs << endl
        << t+t << "num_unique: " << uniq_pairs << endl
        << t+t << "num_ambiguous: " << ambig_pairs << endl
        << t+t << "percent_mapped: "
        << pct(uniq_pairs + ambig_pairs, tot_pairs) << endl
        << t+t << "percent_unique: " << pct(uniq_pairs, tot_pairs) << endl
        << t+t << "percent_ambiguous: " << pct(ambig_pairs, tot_pairs) << endl
        << t+t<< "unique_error:" << endl
        << t+t+t<< "edits: " << edit_distance << endl
        << t+t+t<< "total_bases: " << total_bases << endl
        << t+t+t<< "error_rate: " << pct(edit_distance, total_bases) << endl
        << t   << "num_unmapped: " << unmapped_pairs << endl
        << t   << "num_skipped: " << skipped_pairs << endl
        << t   << "percent_unmapped: " << pct(unmapped_pairs, tot_pairs) << endl
        << t   << "percent_skipped: " << pct(skipped_pairs, tot_pairs) << endl;

    if (!allow_ambig)
      oss << "mate1:" << endl << end1_stats.tostring(1)
          << "mate2:" << endl << end2_stats.tostring(1);
    return oss.str();
  }
};

static void
select_output(const bool allow_ambig,
              const score_t scr1, const score_t scr2,
              const ChromLookup &cl,
              const string &read1, const string &name1,
              const string &read2, const string &name2,
              const string &cig1, const string &cig2,
              pe_element &best, se_element &se1, se_element &se2,
              ostringstream &out) {
  const map_type pe_map_type = format_pe(allow_ambig, best, cl,
                                         read1, read2, name1,
                                         name2, cig1, cig2, out);
  if (pe_map_type == map_unmapped ||
      (!allow_ambig && pe_map_type == map_ambig)) {

    // GS: do not report in mapstats a read that was not reported
    if (pe_map_type == map_unmapped)
      best.reset();

    const bool report_se1 = (se2.ambig() || scr1 >= scr2);
    if (!report_se1 ||
        format_se(allow_ambig, se1, cl, read1, name1, cig1, out) ==
        map_unmapped)
      se1.reset();

    const bool report_se2 = (se1.ambig() || scr2 >= scr1);
    if (!report_se2 ||
        format_se(allow_ambig, se2, cl, read2, name2, cig2, out) ==
                        map_unmapped)
      se2.reset();
  }
}

inline bool
the_comp(const char a, const char b) {
  return (a & b) == 0;
}

inline score_t
popcount64(size_t x) {
  x = (x & 0x5555555555555555ULL) + ((x >> 1) & 0x5555555555555555ULL);
  x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
  x = (x & 0x0F0F0F0F0F0F0F0FULL) + ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL);
  return (x * 0x0101010101010101ULL) >> 56;
}

/* GS: this function counts mismatches between read and genome when
 * they are packed as 64-bit integers, with 16 characters per integer.
 * The number of ones in the AND operation is the number of matches,
 * and there are at most 16 matches since each genome base has only 1
 * out of 4 active bits. Subtracting 16 from the popcount gives the
 * number of mismatches for 16 read bases.
 * The variable "offset" is the remainder of the position modulo 16,
 * and is necessary to adjust the genome bases to align them with the
 * read. The reason why we pad with << (63 - offset) << 1 instead of
 * << 64 - offset is that, when offset = 0 (i.e., the read is aligned
 * with the genome), offsetting 64 positions leads to undefined
 * behavior in some hardware architectures, so the compiler ignores
 * the directive and the number remains unchanged. This is a
 * workaround and there is probably a better way to do it. */
score_t
full_compare(const score_t cutoff, const PackedRead::const_iterator read_end,
             const uint32_t offset, PackedRead::const_iterator read_itr,
             Genome::const_iterator genome_itr) {
  static const score_t max_matches = static_cast<score_t>(16);
  score_t d = 0;
  while (d < cutoff && read_itr != read_end) {
    d += max_matches - popcount64(
      (*read_itr) & /*16 bases from the read*/

      /*16 bases from the padded genome*/
      ((*genome_itr >> offset) | ((*++genome_itr << (63 - offset)) << 1))
    );
    ++read_itr;
  }
  return d;
}

template <const uint16_t strand_code, class result_type>
inline void
check_hits(const uint32_t offset,
           const PackedRead::const_iterator read_st,
           const PackedRead::const_iterator read_end,
           const Genome::const_iterator genome_st,
           const vector<uint32_t>::const_iterator end_idx,
           vector<uint32_t>::const_iterator start_idx,
           result_type &res) {
  for (; start_idx != end_idx && !res.sure_ambig(); ++start_idx) {
    // GS: adds the next candidate to cache while current is compared
    __builtin_prefetch(
      &(*(genome_st + ((*(start_idx + 10) - offset) >> 4)))
    );
    const uint32_t the_pos = *start_idx - offset;
    const score_t cutoff = res.get_cutoff();
    /* GS: the_pos & 15u tells if the position is a multiple of 16, in
     * which case it is aligned with the genome. Otherwise we need to
     * use the unaligned comparison function that offsets genome
     * position by the_pos (mod 16). Multiplied by 4 because each base
     * uses 4 bits */
    const score_t diffs = full_compare(
      cutoff, read_end, ((the_pos & 15u) << 2),
      read_st, genome_st + (the_pos >> 4)
    );
    if (diffs < cutoff)
      res.update(the_pos, diffs, strand_code);
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

static void
find_candidates(const Read::const_iterator read_start,
                const genome_iterator gi,
                const uint32_t read_lim,
                vector<uint32_t>::const_iterator &low,
                vector<uint32_t>::const_iterator &high) {
  uint32_t p = seed::key_weight;
  auto prev_low = low;
  auto prev_high = high;
  for (; p != read_lim && low < high; ++p) {
    // keep last interval with >0 candidates
    prev_low = low;
    prev_high = high;

    // pointer to first 1 in the range
    const vector<uint32_t>::const_iterator first_1 =
      lower_bound(low, high, 1, compare_bases(gi + p));

    const bool the_bit = get_bit(*(read_start + p));
    high = ((the_bit) ? (high) : (first_1));
    low = ((the_bit) ? (first_1) : (low));
  }

  // some bit narrows it down to 0 candidates, roll back to when we
  // had some
  if (low == high) {
    low = prev_low;
    high = prev_high;
  }
}

template <const uint16_t strand_code, class result_type>
void
process_seeds(const uint32_t max_candidates,
              const vector<uint32_t>::const_iterator counter_st,
              const vector<uint32_t>::const_iterator index_st,
              const genome_iterator genome_st,
              const Read &read_seed,
              const PackedRead &packed_read,
              result_type &res) {
  const uint32_t readlen = read_seed.size();
  const PackedRead::const_iterator packed_read_start(begin(packed_read));
  const PackedRead::const_iterator packed_read_end(end(packed_read));
  const Read::const_iterator read_start(begin(read_seed));

  uint32_t k = 0u;

  const uint32_t lim = readlen - seed::key_weight + 1;
  uint32_t read_lim = lim;

  get_1bit_hash(read_start, k);
  for (uint32_t i = 0; i < read_lim; ++i) {
    vector<uint32_t>::const_iterator s_idx(index_st + *(counter_st + k));
    vector<uint32_t>::const_iterator e_idx(index_st + *(counter_st + k + 1));

    // sensitive step with small seeds
    if (s_idx < e_idx && (e_idx - s_idx <= max_candidates))
      check_hits<strand_code>(
        i, packed_read_start, packed_read_end,
        genome_st.itr, e_idx, s_idx, res
      );

    // specific step with whole read as seed
    if (i <= seed::window_size) {
      find_candidates(
        read_start + i, genome_st, readlen - i, s_idx, e_idx
      );
      if (e_idx - s_idx > max_candidates) res.expand();
      if (s_idx < e_idx) {
        check_hits<strand_code>(
          i, packed_read_start, packed_read_end,
          genome_st.itr, e_idx, s_idx, res
        );
      }
    }

    if (res.sure_ambig()) return;
    shift_hash_key(*(read_start + seed::key_weight + i), k);
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

/* GS: this function simply converts the vector<uint8_t> pread
 * to a vector<uint64_t> by putting 16 bases in each element of
 * the packed read. If the read length does not divide 16, we add
 * 1111s to the remaining positions so it divides 16. The remaining
 * bases match all bases in the reference genome
 * */
static void
pack_read(const Read &pread, PackedRead &packed_pread) {
  static const size_t base_match_any = 0xFull;
  const size_t sz = pread.size();
  const size_t num_complete_pos = sz/16;

  // divide by 16 and add an extra position if remainder not 0
  packed_pread.resize((sz + 15)/16);
  PackedRead::iterator it(begin(packed_pread));

  // first add the complete positions (i.e. having all 16 bases)
  size_t pread_ind = 0;
  for (size_t i = 0; i < num_complete_pos; ++i) {
    *it = 0;
    for (size_t j = 0; j < 16; ++j)
      *it |= (static_cast<size_t>(pread[pread_ind++]) << (j << 2));
    ++it;
  }

  // do not fill the flanking position
  if (pread_ind == sz) return;

  // now put only the remaining bases in the last pos. The rest
  // should match any base in the reference
  *it = 0;
  size_t j = 0;
  while (pread_ind < sz)
    *it |= (static_cast<size_t>(pread[pread_ind++]) << ((j++) << 2));

  while (j < 16)
    *it |= base_match_any << ((j++) << 2);
}

using AbismalAlignSimple =
      AbismalAlign<simple_aln::mismatch_score, simple_aln::indel>;

static score_t
align_read(const Read &pread, se_element &res,
         uint32_t &len, string &cigar, AbismalAlignSimple &aln) {
  const score_t readlen = static_cast<score_t>(pread.size());
  const score_t ans = ((res.diffs < simple_aln::min_diffs_to_align || res.diffs == readlen) ?
    aln.trim_ends(pread, res.pos, len, cigar) :
    aln.align    (pread, res.pos, len, cigar));
  res.diffs = simple_aln::edit_distance(ans, len, cigar);
  return ans;
}

static score_t
align_se_candidates(const Read &pread_t, const Read &pread_t_rc,
                    const Read &pread_a, const Read &pread_a_rc,
                    se_candidates &res, se_element &best,
                    string &cigar,
                    AbismalAlignSimple &aln) {
  const score_t readlen = static_cast<score_t>(pread_t.size());
  score_t best_score = 0;
  uint32_t len = 0;
  uint32_t best_len = 0;
  se_element s;
  string cand_cigar;

  res.prepare_for_alignments();
  const vector<se_element>::const_iterator lim(begin(res.v) + res.sz);

  vector<se_element>::const_iterator it(begin(res.v));
  for (; it != lim && it->empty(); ++it);
  for (; it != lim; ++it) {
    s = *it;
    if (valid_hit(s, readlen)) {
      const score_t scr = align_read(
        s.rc() ? (s.elem_is_a_rich() ? pread_a_rc : pread_t_rc) :
                 (s.elem_is_a_rich() ? pread_a : pread_t)
        , s, len, cand_cigar, aln
      );
      if (scr > best_score) {
        best = s; // GS: ambig is unset on s
        best_len = len;
        best_score = scr;
        cigar = std::move(cand_cigar);
      }

      /* GS: need to check != here in case alignment caused the pos
       * to be the same, even if not excluded by unique(). This is
       * the case when there are indels, where the same position can
       * be one apart because of k-mers sampled before and after the
       * indel */
      else if (scr == best_score && s != best)
        best.set_ambig();
    }
  }

  if (!valid(best, best_len))
    best.reset();

  return best_score;
}

template <const  conversion_type conv>
inline void
map_single_ended(const bool VERBOSE, const bool allow_ambig,
                 const AbismalIndex &abismal_index, ReadLoader &rl,
                 se_map_stats &se_stats, ostream &out,
                 ProgressBar &progress) {
  const vector<uint32_t>::const_iterator counter_st(begin(abismal_index.counter));
  const vector<uint32_t>::const_iterator index_st(begin(abismal_index.index));
  const genome_iterator genome_st(begin(abismal_index.genome));
  const uint32_t max_candidates = abismal_index.max_candidates;

  // batch variables used in reporting the SAM entry
  vector<string> names;
  vector<string> reads;
  vector<string> cigar;
  vector<se_element> bests;

  names.reserve(ReadLoader::batch_size);
  reads.reserve(ReadLoader::batch_size);

  cigar.resize(ReadLoader::batch_size);
  bests.resize(ReadLoader::batch_size);

  // pre-allocated variabes used idependently in each read
  Read pread, pread_rc;
  PackedRead packed_pread;
  se_candidates res;

  size_t the_byte;

  string out_str;
  out_str.reserve(ReadLoader::batch_size*500);
  ostringstream out_stream(out_str);

  while (rl.good()) {
#pragma omp critical
    {
      rl.load_reads(names, reads);
      the_byte = rl.get_current_byte();
    }

    size_t max_batch_read_length = 0;
    update_max_read_length(max_batch_read_length, reads);

    AbismalAlignSimple aln(genome_st, max_batch_read_length);

    const size_t n_reads = reads.size();

    out_stream.clear();
    out_stream.str("");
    for (size_t i = 0; i < n_reads; ++i) {
      res.reset(reads[i].size());
      bests[i].reset();
      if (!reads[i].empty()) {
        prep_read<conv>(reads[i], pread);
        pack_read(pread, packed_pread);
        process_seeds<get_strand_code('+', conv)>(
          max_candidates, counter_st, index_st, genome_st, pread,
          packed_pread, res
        );

        const string read_rc(revcomp(reads[i]));
        prep_read<!conv>(read_rc, pread_rc);
        pack_read(pread_rc, packed_pread);
        process_seeds<get_strand_code('-', conv)>(
          max_candidates, counter_st, index_st, genome_st, pread_rc,
          packed_pread, res
        );
        align_se_candidates(
          pread, pread_rc, pread, pread_rc, res, bests[i], cigar[i], aln
        );
        if (format_se(allow_ambig, bests[i], abismal_index.cl, reads[i],
              names[i], cigar[i], out_stream) == map_unmapped)
          bests[i].reset();
      }
    }

#pragma omp critical
    {
      out.write(out_stream.str().c_str(), out_stream.str().size());
      for (size_t i = 0; i < n_reads; ++i)
        se_stats.update(allow_ambig, reads[i], cigar[i], bests[i]);

      if (VERBOSE && progress.time_to_report(the_byte))
        progress.report(cerr, the_byte);
    }
  }
}

inline void
map_single_ended_rand(const bool VERBOSE, const bool allow_ambig,
                      const AbismalIndex &abismal_index, ReadLoader &rl,
                      se_map_stats &se_stats, ostream &out,
                      ProgressBar &progress) {
  const vector<uint32_t>::const_iterator counter_st(begin(abismal_index.counter));
  const vector<uint32_t>::const_iterator index_st(begin(abismal_index.index));
  const uint32_t max_candidates = abismal_index.max_candidates;

  const genome_iterator genome_st(begin(abismal_index.genome));

  vector<string> names;
  vector<string> reads;
  vector<string> cigar;
  vector<se_element> bests;

  names.reserve(ReadLoader::batch_size);
  reads.reserve(ReadLoader::batch_size);
  cigar.resize(ReadLoader::batch_size);
  bests.resize(ReadLoader::batch_size);

  // GS: pre-allocated variables used once per read
  // and not used for reporting
  Read pread_t, pread_t_rc, pread_a, pread_a_rc;
  PackedRead packed_pread;
  se_candidates res;

  size_t the_byte;
  string out_str;
  out_str.reserve(ReadLoader::batch_size*500);
  ostringstream out_stream(out_str);

  while (rl.good()) {
#pragma omp critical
    {
      rl.load_reads(names, reads);
      the_byte = rl.get_current_byte();
    }

    size_t max_batch_read_length = 0;
    update_max_read_length(max_batch_read_length, reads);

    AbismalAlignSimple aln(genome_st, max_batch_read_length);

    const size_t n_reads = reads.size();
    out_stream.clear();
    out_stream.str("");
    for (size_t i = 0; i < n_reads; ++i) {
      res.reset(reads[i].size());
      bests[i].reset();
      if (!reads[i].empty()) {
        prep_read<t_rich>(reads[i], pread_t);
        pack_read(pread_t, packed_pread);
        process_seeds<get_strand_code('+', t_rich)>(
          max_candidates, counter_st, index_st, genome_st, pread_t,
          packed_pread, res
        );

        prep_read<a_rich>(reads[i], pread_a);
        pack_read(pread_a, packed_pread);
        process_seeds<get_strand_code('+', a_rich)>(
          max_candidates, counter_st, index_st, genome_st, pread_a,
          packed_pread, res
        );

        const string read_rc(revcomp(reads[i]));
        prep_read<t_rich>(read_rc, pread_t_rc);
        pack_read(pread_t_rc, packed_pread);
        process_seeds<get_strand_code('-', a_rich)>(
          max_candidates, counter_st, index_st, genome_st, pread_t_rc,
          packed_pread, res
        );

        prep_read<a_rich>(read_rc, pread_a_rc);
        pack_read(pread_a_rc, packed_pread);
        process_seeds<get_strand_code('-', t_rich)>(
          max_candidates, counter_st, index_st, genome_st, pread_a_rc,
          packed_pread, res
        );

        align_se_candidates(
          pread_t, pread_t_rc, pread_a, pread_a_rc, res, bests[i], cigar[i], aln
        );
        if (format_se(allow_ambig, bests[i], abismal_index.cl, reads[i],
              names[i], cigar[i], out_stream) == map_unmapped)
          bests[i].reset();
      }
    }
#pragma omp critical
    {
      out.write(out_stream.str().c_str(), out_stream.str().size());
      for (size_t i = 0; i < n_reads; ++i)
        se_stats.update(allow_ambig, reads[i], cigar[i], bests[i]);

      if (VERBOSE && progress.time_to_report(the_byte))
        progress.report(cerr, the_byte);
    }
  }
}

template <const conversion_type conv, const bool random_pbat>
void
run_single_ended(const bool VERBOSE,
                 const bool allow_ambig,
                 const string &reads_file,
                 const AbismalIndex &abismal_index,
                 se_map_stats &se_stats,
                 ostream &out) {
  ReadLoader rl(reads_file);
  ProgressBar progress(get_filesize(reads_file), "mapping reads");

  if (VERBOSE)
    progress.report(cerr, 0);

  double start_time = omp_get_wtime();
  if (VERBOSE && progress.time_to_report(rl.get_current_byte()))
    progress.report(cerr, rl.get_current_byte());

#pragma omp parallel for
  for (int i = 0; i < omp_get_num_threads(); ++i) {
    if (random_pbat)
      map_single_ended_rand(VERBOSE, allow_ambig,
        abismal_index, rl, se_stats, out, progress);
    else
      map_single_ended<conv>(VERBOSE, allow_ambig,
          abismal_index, rl, se_stats, out, progress);
  }
  if (VERBOSE) {
    cerr << "\n";
    print_with_time("total mapping time: " + to_string(omp_get_wtime() - start_time) + "s");
  }
}

static void
best_single(const pe_candidates &pres, se_candidates &res) {
  const vector<se_element>::const_iterator lim(begin(pres.v) + pres.sz);
  for (vector<se_element>::const_iterator i(begin(pres.v));
      i != lim && !res.sure_ambig(); ++i) {
    res.update(i->pos, i->diffs, i->flags);
  }
}

template <const bool swap_ends>
static void
best_pair(const pe_candidates &res1, const pe_candidates &res2,
          const Read &pread1, const Read &pread2,
          string &cig1, string &cig2,
          AbismalAlignSimple &aln, pe_element &best) {
  const vector<se_element>::const_iterator j1_beg(begin(res1.v));
  vector<se_element>::const_iterator j1(j1_beg);
  vector<se_element>::const_iterator j2(begin(res2.v));
  const vector<se_element>::const_iterator j1_end = j1 + res1.sz;
  const vector<se_element>::const_iterator j2_end = j2 + res2.sz;
  // skips the empty elements for faster nested loops

  const uint32_t readlen1 = pread1.size();
  const uint32_t readlen2 = pread2.size();

  uint32_t len1 = 0, len2 = 0;
  score_t scr, scr2;
  se_element s1, s2;
  string cand_cig1, cand_cig2;

  // GS: upper bound on cigar size
  cand_cig1.reserve(2*readlen1);
  cand_cig2.reserve(2*readlen2);

  for (; j2 != j2_end; ++j2) {
    s2 = *j2;
    if (valid_hit(s2, readlen2)) {
      uint32_t lim = s2.pos + readlen2;

      // rewind to first concordant position. Needed in case of
      // many-to-many equivalence between pairs.
      for (; j1 != j1_beg && j1->pos + pe_element::max_dist >= lim; --j1);
      for (; j1 != j1_end && j1->pos + pe_element::max_dist < lim; ++j1);

      scr2 = 0;
      for (; j1 != j1_end && j1->pos + pe_element::min_dist <= lim; ++j1) {
        s1 = *j1;
        if (valid_hit(s1, readlen1)) {
          // GS: guarantees that j2 is aligned only once
          if (scr2 == 0) {
            scr2 = align_read(pread2, s2, len2, cand_cig2, aln);
            lim = s2.pos + len2;
          }
          scr = scr2 + align_read(pread1, s1, len1, cand_cig1, aln);
          // GS: only accept if length post alignment is still within limits
          if (valid_pair(s1, s2, len1, len2) &&
              (s1.pos + pe_element::max_dist >= lim) &&
              (s1.pos + pe_element::min_dist <= lim)) {

            const se_element r1 = (swap_ends ? s2 : s1);
            const se_element r2 = (swap_ends ? s1 : s2);
            if (best.update(scr, r1, r2)) {
              cig1 = std::move(cand_cig1);
              if (!cand_cig2.empty())
                cig2 = std::move(cand_cig2);

              // found 2 exact matches, enough to stop mating
              if (best.sure_ambig())
                return;
            }
          }
        }
      }
    }
  }
}

template <const bool swap_ends>
void
select_maps(Read &pread1, Read &pread2,
            string &cig1, string &cig2,
            pe_candidates &res1, pe_candidates &res2,
            se_candidates &res_se1, se_candidates &res_se2,
            AbismalAlignSimple &aln,  pe_element &best) {
  res1.prepare_for_mating();
  res2.prepare_for_mating();
  best_pair<swap_ends>(res1, res2, pread1, pread2, cig1, cig2, aln, best);
  best_single(res1, res_se1);
  best_single(res2, res_se2);
}


template <const bool cmp, const bool swap_ends,
          const uint16_t strand_code1, const uint16_t strand_code2>
void
map_fragments(const uint32_t max_candidates,
              const string &read1, const string &read2,
              const vector<uint32_t>::const_iterator counter_st,
              const vector<uint32_t>::const_iterator index_st,
              const genome_iterator genome_st,
              Read &pread1, Read &pread2, PackedRead &packed_pread,
              string &cigar1, string &cigar2,
              AbismalAlignSimple &aln,
              pe_candidates &res1, pe_candidates &res2,
              se_candidates &res_se1, se_candidates &res_se2,
              pe_element &best) {
  res1.reset(read1.size());
  res2.reset(read2.size());

  if (!read1.empty() && !read2.empty()) {
    prep_read<cmp>(read1, pread1);
    pack_read(pread1, packed_pread);
    process_seeds<strand_code1>(max_candidates, counter_st,
      index_st, genome_st, pread1, packed_pread, res1
    );

    const string read_rc(revcomp(read2));
    prep_read<cmp>(read_rc, pread2);
    pack_read(pread2, packed_pread);
    process_seeds<strand_code2>(max_candidates, counter_st,
      index_st, genome_st, pread2, packed_pread, res2
    );

    select_maps<swap_ends>(pread1, pread2, cigar1, cigar2,
      res1, res2, res_se1, res_se2, aln, best);
  }
}

template <const conversion_type conv>
inline void
map_paired_ended(const bool VERBOSE,
                 const bool allow_ambig,
                 const AbismalIndex &abismal_index,
                 ReadLoader &rl1, ReadLoader &rl2,
                 pe_map_stats &pe_stats, ostream &out,
                 ProgressBar &progress) {
  const vector<uint32_t>::const_iterator counter_st(begin(abismal_index.counter));
  const vector<uint32_t>::const_iterator index_st(begin(abismal_index.index));
  const uint32_t max_candidates = abismal_index.max_candidates;

  const genome_iterator genome_st(begin(abismal_index.genome));

  vector<string> names1, reads1, cigar1;
  vector<string> names2, reads2, cigar2;

  vector<pe_element> bests;
  vector<se_element> bests_se1;
  vector<se_element> bests_se2;

  names1.reserve(ReadLoader::batch_size);
  reads1.reserve(ReadLoader::batch_size);
  cigar1.resize(ReadLoader::batch_size);

  names2.reserve(ReadLoader::batch_size);
  reads2.reserve(ReadLoader::batch_size);
  cigar2.resize(ReadLoader::batch_size);

  bests.resize(ReadLoader::batch_size);
  bests_se1.resize(ReadLoader::batch_size);
  bests_se2.resize(ReadLoader::batch_size);

  // GS: pre-allocated variables used once per read
  // and not used for reporting
  Read pread1, pread1_rc, pread2, pread2_rc;
  PackedRead packed_pread;
  pe_candidates res1(abismal_index.pe_heap_size);
  pe_candidates res2(abismal_index.pe_heap_size);
  se_candidates res_se1;
  se_candidates res_se2;
  score_t scr1;
  score_t scr2;

  size_t the_byte;
  string out_str;
  out_str.reserve(ReadLoader::batch_size*500);
  ostringstream out_stream(out_str);

  while (rl1.good() && rl2.good()) {
#pragma omp critical
    {
      rl1.load_reads(names1, reads1);
      rl2.load_reads(names2, reads2);
      the_byte = rl1.get_current_byte();
    }

    size_t max_batch_read_length = 0;
    update_max_read_length(max_batch_read_length, reads1);
    update_max_read_length(max_batch_read_length, reads2);

    AbismalAlignSimple aln(genome_st, max_batch_read_length);

    const size_t n_reads = reads1.size();
    out_stream.clear();
    out_stream.str("");
    for (size_t i = 0 ; i < n_reads; ++i) {
      const uint32_t readlen1 = reads1[i].size();
      const uint32_t readlen2 = reads2[i].size();

      res1.reset(readlen1);
      res2.reset(readlen2);
      res_se1.reset(readlen1);
      res_se2.reset(readlen2);

      bests[i].reset(readlen1, readlen2);
      bests_se1[i].reset(readlen1);
      bests_se2[i].reset(readlen2);

      scr1 = 0;
      scr2 = 0;

      map_fragments<conv, false,
                   get_strand_code('+',conv),
                   get_strand_code('-', flip_conv(conv))>(
        max_candidates, reads1[i], reads2[i], counter_st,
        index_st, genome_st, pread1, pread2_rc, packed_pread,
        cigar1[i], cigar2[i], aln, res1, res2,
        res_se1, res_se2, bests[i]
      );

      map_fragments<!conv, true,
                   get_strand_code('+', flip_conv(conv)),
                   get_strand_code('-', conv)>(
        max_candidates, reads2[i], reads1[i], counter_st,
        index_st, genome_st, pread2, pread1_rc, packed_pread,
        cigar2[i], cigar1[i], aln, res2, res1,
        res_se2, res_se1, bests[i]
      );
      if (!bests[i].should_report(allow_ambig)) {
        scr1 = align_se_candidates(pread1, pread1_rc, pread1, pread1_rc,
                            res_se1, bests_se1[i], cigar1[i], aln);

        scr2 = align_se_candidates(pread2, pread2_rc, pread2, pread2_rc,
                            res_se2, bests_se2[i], cigar2[i], aln);
      }

      select_output(allow_ambig, scr1, scr2, abismal_index.cl,
          reads1[i], names1[i],reads2[i], names2[i], cigar1[i], cigar2[i],
          bests[i], bests_se1[i], bests_se2[i], out_stream);
    }

#pragma omp critical
    {
      out.write(out_stream.str().c_str(), out_stream.str().size());
      for (size_t i = 0; i < n_reads; ++i) {
        pe_stats.update(allow_ambig, reads1[i], reads2[i],
            cigar1[i], cigar2[i], bests[i], bests_se1[i], bests_se2[i]);
      }
      if (VERBOSE && progress.time_to_report(the_byte))
        progress.report(cerr, the_byte);
    }
  }
}

inline void
map_paired_ended_rand(const bool VERBOSE, const bool allow_ambig,
                      const AbismalIndex &abismal_index,
                      ReadLoader &rl1, ReadLoader &rl2,
                      pe_map_stats &pe_stats, ostream &out,
                      ProgressBar &progress) {
  const vector<uint32_t>::const_iterator counter_st(begin(abismal_index.counter));
  const vector<uint32_t>::const_iterator index_st(begin(abismal_index.index));
  const uint32_t max_candidates = abismal_index.max_candidates;

  const genome_iterator genome_st(begin(abismal_index.genome));

  vector<string> names1, reads1, cigar1;
  vector<string> names2, reads2, cigar2;

  vector<pe_element> bests;
  vector<se_element> bests_se1;
  vector<se_element> bests_se2;

  names1.reserve(ReadLoader::batch_size);
  reads1.reserve(ReadLoader::batch_size);
  cigar1.resize(ReadLoader::batch_size);

  names2.reserve(ReadLoader::batch_size);
  reads2.reserve(ReadLoader::batch_size);
  cigar2.resize(ReadLoader::batch_size);

  bests.resize(ReadLoader::batch_size);
  bests_se1.resize(ReadLoader::batch_size);
  bests_se2.resize(ReadLoader::batch_size);

  Read pread1_t, pread1_t_rc, pread2_t, pread2_t_rc;
  Read pread1_a, pread1_a_rc, pread2_a, pread2_a_rc;
  PackedRead packed_pread;
  pe_candidates res1(abismal_index.pe_heap_size);
  pe_candidates res2(abismal_index.pe_heap_size);
  se_candidates res_se1;
  se_candidates res_se2;
  score_t scr1;
  score_t scr2;

  size_t the_byte;
  string out_str;
  out_str.reserve(ReadLoader::batch_size*500);
  ostringstream out_stream(out_str);

  while (rl1.good() && rl2.good()) {
#pragma omp critical
    {
      rl1.load_reads(names1, reads1);
      rl2.load_reads(names2, reads2);
      the_byte = rl1.get_current_byte();
    }

    size_t max_batch_read_length = 0;
    update_max_read_length(max_batch_read_length, reads1);
    update_max_read_length(max_batch_read_length, reads2);

    AbismalAlignSimple aln(genome_st, max_batch_read_length);

    const size_t n_reads = reads1.size();
    out_stream.clear();
    out_stream.str("");
    for (size_t i = 0 ; i < n_reads; ++i) {
      const uint32_t readlen1 = reads1[i].size();
      const uint32_t readlen2 = reads1[i].size();

      res1.reset(readlen1);
      res2.reset(readlen2);
      res_se1.reset(readlen1);
      res_se2.reset(readlen2);
      bests[i].reset(readlen1, readlen2);
      bests_se1[i].reset(readlen1);
      bests_se2[i].reset(readlen2);
      scr1 = 0;
      scr2 = 0;

      // GS: (1) T/A-rich +/- strand
      map_fragments<t_rich, false,
                   get_strand_code('+', t_rich),
                   get_strand_code('-', a_rich)>(
        max_candidates, reads1[i], reads2[i], counter_st,
        index_st, genome_st, pread1_t, pread2_a_rc, packed_pread,
        cigar1[i], cigar2[i], aln, res1, res2, res_se1, res_se2, bests[i]
      );

      // GS: (2) T/A-rich, -/+ strand
      map_fragments<a_rich, true,
                   get_strand_code('+', a_rich),
                   get_strand_code('-', t_rich)>(
        max_candidates, reads2[i], reads1[i], counter_st,
        index_st, genome_st, pread2_a, pread1_t_rc, packed_pread,
        cigar2[i], cigar1[i], aln, res2, res1, res_se2, res_se1, bests[i]
      );

      // GS: (3) A/T-rich +/- strand
      map_fragments<a_rich, false,
                   get_strand_code('+', a_rich),
                   get_strand_code('-', t_rich)>(
        max_candidates, reads1[i], reads2[i], counter_st,
        index_st, genome_st, pread1_a, pread2_t_rc, packed_pread,
        cigar1[i], cigar2[i], aln, res1, res2, res_se1, res_se2, bests[i]
      );

      // GS: (4) A/T-rich, -/+ strand
      map_fragments<t_rich, true,
                   get_strand_code('+', t_rich),
                   get_strand_code('-', a_rich)>(
        max_candidates, reads2[i], reads1[i], counter_st,
        index_st, genome_st, pread2_t, pread1_a_rc, packed_pread,
        cigar2[i], cigar1[i], aln, res2, res1, res_se2, res_se1, bests[i]
      );

      // GS: align best SE candidates if no concordant pairs found
      if (!bests[i].should_report(allow_ambig)) {
        scr1 = align_se_candidates(
          pread1_t, pread1_t_rc, pread1_a, pread1_a_rc,
          res_se1, bests_se1[i], cigar1[i], aln
        );
        scr2 = align_se_candidates(
          pread2_t, pread2_t_rc, pread2_a, pread2_a_rc,
          res_se2, bests_se2[i], cigar2[i], aln
        );
      }
      select_output(allow_ambig, scr1, scr2, abismal_index.cl,
        reads1[i], names1[i], reads2[i], names2[i], cigar1[i], cigar2[i],
        bests[i], bests_se1[i], bests_se2[i], out_stream);
    }

#pragma omp critical
    {
      out.write(out_stream.str().c_str(), out_stream.str().size());
      for (size_t i = 0; i < n_reads; ++i)
        pe_stats.update(allow_ambig, reads1[i], reads2[i],
            cigar1[i], cigar2[i], bests[i], bests_se1[i], bests_se2[i]);

      if (VERBOSE && progress.time_to_report(the_byte))
        progress.report(cerr, the_byte);
    }
  }
}

template <const conversion_type conv, const bool random_pbat>
void
run_paired_ended(const bool VERBOSE,
                 const bool allow_ambig,
                 const string &reads_file1,
                 const string &reads_file2,
                 const AbismalIndex &abismal_index,
                 pe_map_stats &pe_stats,
                 ostream &out) {
  ReadLoader rl1(reads_file1);
  ReadLoader rl2(reads_file2);
  ProgressBar progress(get_filesize(reads_file1), "mapping reads");

  if (VERBOSE)
    progress.report(cerr, 0);

  double start_time = omp_get_wtime();

#pragma omp parallel for
  for (int i = 0; i < omp_get_num_threads(); ++i) {
    if (random_pbat)
      map_paired_ended_rand(VERBOSE, allow_ambig,
          abismal_index, rl1, rl2, pe_stats, out, progress);

    else
      map_paired_ended<conv>(VERBOSE, allow_ambig,
          abismal_index, rl1, rl2, pe_stats, out, progress);
  }

  if (VERBOSE) {
    cerr << "\n";
    print_with_time("total mapping time: " + to_string(omp_get_wtime() - start_time) + "s");
  }
}

int main(int argc, const char **argv) {

  try {
    static const string ABISMAL_VERSION = "1.0.0";
    bool VERBOSE = false;
    bool GA_conversion = false;
    bool allow_ambig = false;
    bool pbat_mode = false;
    bool random_pbat = false;
    bool sensitive = false;
    int n_threads = 1;
    string index_file = "";
    string genome_file = "";
    string outfile = "";
    string stats_outfile = "";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "map bisulfite converted reads",
                           "<reads-fq1> [<reads-fq2>]");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("index", 'i', "index file", false, index_file);
    opt_parse.add_opt("genome", 'g', "genome file (FASTA)", false, genome_file);
    opt_parse.add_opt("outfile", 'o', "output file (SAM) [stdout]", false, outfile);
    opt_parse.add_opt("stats", 's', "map statistics file (YAML)",
                      false, stats_outfile);
    opt_parse.add_opt("sensitive", 'x', "run abismal in sensitive mode", false, sensitive);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
    opt_parse.add_opt("min-frag", 'l', "min fragment size (pe mode)",
                      false, pe_element::min_dist);
    opt_parse.add_opt("max-frag", 'L', "max fragment size (pe mode)",
                      false, pe_element::max_dist);
    opt_parse.add_opt("max-distance", 'm',
                      "max fractional edit distance",
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
      cerr << "help requested" << endl;
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << "about requested" << endl;
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << "missing required option" << endl;
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 1 && leftover_args.size() != 2) {
      cerr << "leftover size = " << leftover_args.size() << endl;
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (n_threads <= 0) {
      cerr << "please choose a positive number of threads" << endl;
      return EXIT_SUCCESS;
    }
    if (index_file.empty() == genome_file.empty()) {
      cerr << "please select either an index file (-i) or a genome file (-g)"
           << endl;
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
        print_with_time("input (PE): " + reads_file + ", " + reads_file2);
      else
        print_with_time("input (SE): " + reads_file);

      print_with_time("output (SAM): " + (outfile.empty() ? "[stdout]" : outfile));

      if (!stats_outfile.empty())
        print_with_time("map statistics (YAML): " + stats_outfile);
    }

    AbismalIndex abismal_index;

    const double start_time = omp_get_wtime();
    if (!index_file.empty()) {
      if (VERBOSE)
        print_with_time("loading index " + index_file);
      abismal_index.read(index_file);

      if (VERBOSE)
        print_with_time("loading time: " + to_string(omp_get_wtime() - start_time) + "s");
    }
    else {
      if (VERBOSE)
        print_with_time("indexing genome " + genome_file);
      abismal_index.create_index(genome_file);
      if (VERBOSE)
        print_with_time("indexing time: " + to_string(omp_get_wtime() - start_time) + "s");
    }

    abismal_index.calc_mapping_parameters(sensitive);

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
        run_single_ended<a_rich, false>(
          VERBOSE, allow_ambig, reads_file, abismal_index, se_stats, out
        );
      else if (random_pbat)
        run_single_ended<t_rich, true>(
          VERBOSE, allow_ambig, reads_file, abismal_index, se_stats, out
        );
      else
        run_single_ended<t_rich, false>(
          VERBOSE, allow_ambig, reads_file, abismal_index, se_stats, out
        );
    }
    else {
      if (pbat_mode)
        run_paired_ended<a_rich, false>(VERBOSE, allow_ambig, reads_file,
            reads_file2, abismal_index, pe_stats,
            out);
      else if (random_pbat)
        run_paired_ended<t_rich, true>(VERBOSE,  allow_ambig, reads_file,
            reads_file2, abismal_index, pe_stats,
            out);
      else
        run_paired_ended<t_rich, false>(VERBOSE, allow_ambig, reads_file,
            reads_file2, abismal_index, pe_stats,
            out);
    }

    if (!stats_outfile.empty()) {
      std::ofstream stats_of(stats_outfile.c_str(), std::ios::binary);
      stats_of << (reads_file2.empty() ?
                    se_stats.tostring() : pe_stats.tostring(allow_ambig));
    }
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
