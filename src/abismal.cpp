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

#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"
#include "OptionParser.hpp"
#include "zlib_wrapper.hpp"
#include "AbismalIndex.hpp"
#include "AbismalAlign.hpp"
#include "GenomicRegion.hpp"
#include "dna_four_bit.hpp"

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

// Aliases for different types used throughout abismal
typedef uint16_t flags_t; // every bit is a flag
typedef int16_t score_t; // alignment score
typedef bool cmp_t; // match/mismatch type
typedef vector<char> Read; //4-bit encoding of reads

enum conversion_type { t_rich = false, a_rich = true };

constexpr conversion_type
flip_conv(const conversion_type &conv) {
  return conv == t_rich ? a_rich : t_rich;
}

namespace SamFlags {
  static const flags_t read_paired = 0x1;
  static const flags_t read_pair_mapped = 0x2;
  static const flags_t read_unmapped = 0x4;
  static const flags_t mate_unmapped = 0x8;
  static const flags_t read_rc = 0x10;
  static const flags_t mate_rc = 0x20;
  static const flags_t template_first = 0x40;
  static const flags_t template_second = 0x40;
  static const flags_t secondary_aln = 0x100;
  static const flags_t below_quality = 0x200;
  static const flags_t pcr_duplicate = 0x400;
  static const flags_t supplementary_aln = 0x800;
};

namespace BSFlags {
  static const flags_t a_rich = 0x1000;
  static const flags_t ambig = 0x2000;
};

constexpr flags_t
get_strand_code(const char strand, const conversion_type conv) {
  return (((strand == '-')  ? SamFlags::read_rc : 0) |
          ((conv == a_rich) ? BSFlags::a_rich : 0));
}

constexpr flags_t
flip_strand_code(const flags_t sc) {
  return (sc ^ SamFlags::read_rc) ^ BSFlags::a_rich;
}

constexpr bool
is_a_rich(const flags_t flags) {
  return (flags & BSFlags::a_rich) != 0;
}

constexpr bool
is_rc(const flags_t flags) {
  return (flags & SamFlags::read_rc) != 0;
}

namespace align_scores {
  static const score_t match = 1,
    mismatch = -1,
    indel = -1;
};

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
                     [](const char c) {return c != 'N';}) < min_length)
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

  static uint32_t min_length;
};
uint32_t ReadLoader::min_length = 32;


struct se_element {
  uint32_t pos;
  score_t diffs;
  flags_t flags;

  se_element() : pos(0),
                 diffs(invalid_hit_diffs + 1),
                 flags(0) {}

  se_element(const uint32_t p, const score_t d, const flags_t f) :
    pos(p), diffs(d), flags(f) {}

  bool operator==(const se_element &rhs) const {
    return diffs == rhs.diffs && pos == rhs.pos;
  }

  // this is used to keep PE candidates sorted in
  // the max heap
  bool operator<(const se_element &rhs) const {
    return diffs < rhs.diffs;
  }

  // this is used for deciding on ambiguity, and can
  // be overriden. For now it is the same implementation
  // as above but will be changed eventually.
  bool is_better_than (const se_element &rhs) const {
    return diffs < rhs.diffs;
  }

  bool rc() const {return is_rc(flags);}
  bool a_rich() const {return is_a_rich(flags);}
  bool valid() const {return diffs <= max_diffs;}
  bool valid_hit() const {return diffs <= invalid_hit_diffs;}
  bool do_align() const {return !valid() && valid_hit();}
  char strand() const {return rc() ? '-' : '+';}
  void flip_strand() {flags = flip_strand_code(flags);}
  void reset() { diffs = invalid_hit_diffs + 1; }

  bool is_equal_to (const se_element &rhs) const {
    return (diffs == rhs.diffs) && (pos != rhs.pos);
  }
  static score_t max_diffs;
  static score_t invalid_hit_diffs;
  static uint16_t min_aligned_length;
};

score_t se_element::invalid_hit_diffs = 40;
score_t se_element::max_diffs = 6;
uint16_t se_element::min_aligned_length = ReadLoader::min_length;

struct se_result {
  se_element best, second_best;
  se_result() : best(se_element()), second_best(se_element()) {}

  se_result(const uint32_t p, const score_t d, const flags_t s) :
    best(se_element(p, d, s)), second_best(se_element()) {}

  bool operator==(const se_result &rhs) const {
    return best == rhs.best;
  }

  bool operator<(const se_result &rhs) const {
    return best < rhs.best;
  }
  bool valid() const {return best.valid();}
  void update(const uint32_t p, const score_t d, const flags_t s) {
    if (p == best.pos || p == second_best.pos) return;
    const se_element cand(p, d, s);
    if (cand.is_better_than(second_best)) second_best = cand;
    if (second_best.is_better_than(best)) std::swap(best, second_best);
  }

  bool ambig() const {return valid() && best.is_equal_to(second_best);}
  bool sure_ambig(uint32_t seed_number = 0) const {
    return ambig() &&
      (best.diffs == 0 || (best.diffs == 1 && seed_number > 0));
  }

  void reset() {best.reset(); second_best.reset();}
  uint32_t get_cutoff() const {return second_best.diffs;}
  flags_t flags() const {return best.flags;}
};

bool
format_se(se_result res, const ChromLookup &cl,
          string &read, const string &read_name,
          const string &cigar, ofstream &out) {
  uint32_t offset = 0, chrom_idx = 0;
  se_element s = res.best;
  if (s.valid() &&
      cl.get_chrom_idx_and_offset(s.pos, read.size(), chrom_idx, offset)) {
    if (s.a_rich()) { // since SE, this is only for GA conversion
      revcomp_inplace(read);
      s.flip_strand();
    }
    out << cl.names[chrom_idx] << '\t'
        << offset << '\t'
        << offset + read.length() << '\t'
        << read_name << '\t'
        << s.diffs << '\t'
        << s.strand() << '\t'
        << read << '\t'
        << cigar << '\n';
    //       << res.flags() << '\n';
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
  bool a_rich() const {return r1.a_rich();}
  char strand() const {return r1.strand();}
  score_t diffs() const { return r1.diffs + r2.diffs; }
  flags_t flags() const { return r1.flags; }

  bool valid() const {return r1.diffs + r2.diffs < 2*se_element::max_diffs;}
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
  pe_element best, second_best;
  pe_result() {}
  pe_result(const pe_element &a, const pe_element &b)
    : best(a), second_best(b) {}

  void reset() {
    best.reset();
    second_best.reset();
  }

  bool ambig() const {
    return best.is_equal_to(second_best);
  }
  bool valid() const{return best.valid();}
  flags_t flags() const { return best.flags(); }
  bool update(const pe_element &p) {
    if (p.is_better_than(second_best)) second_best = p;
    if (second_best.is_better_than(best)) {
      std::swap(best, second_best);
      return true; //best has been updated
    }
    return false;
  }
};


static int
get_spacer_rlen(const bool rc, const long int s1, const long int e1,
                const long int s2, const long int e2) {
  return rc ? s1 - e2 : s2 - e1;
}

static int
get_head_rlen(const bool rc, const long int s1, const long int e1,
              const long int s2, const long int e2) {
  return rc ? e1 - e2 : s2 - s1;
}

static int
get_overlap_rlen(const bool rc, const long int s1, const long int e1,
                 const long int s2, const long int e2) {
  return rc ? e1 - s2 : e2 - s1;
}

inline bool
chrom_and_posn(const ChromLookup &cl, const string &cig, const uint32_t p,
               uint32_t &r_p, uint32_t &r_e, uint32_t &r_chr) {
  const uint32_t ref_ops = cigar_rseq_ops(cig);
  if (!cl.get_chrom_idx_and_offset(p, ref_ops, r_chr, r_p)) return false;
  r_e = r_p + ref_ops;
  return true;
}

bool
get_pe_overlap(GenomicRegion &gr,
               const bool rc,
               uint32_t r_s1, uint32_t r_e1, uint32_t chr1,
               uint32_t r_s2, uint32_t r_e2, uint32_t chr2,
               string &read1, string &read2,
               string &cig1, string &cig2) {

  const int spacer = get_spacer_rlen(rc, r_s1, r_e1, r_s2, r_e2);
  if (spacer >= 0) {
    /* fragments longer than or equal to 2x read length: this size of
     * the spacer ("_") is determined based on the reference positions
     * of the two ends, and depends on whether the mapping is on the
     * negative strand of the genome.
     *
     * left                                                             right
     * r_s1                         r_e1   r_s2                         r_e2
     * [------------end1------------]______[------------end2------------]
     */
    gr.set_name("FRAG_L:" + gr.get_name()); // DEBUG
    read1 += string(spacer, 'N');
    read1 += read2;

    cig1 += std::to_string(spacer) + "N";
    cig1 += cig2;
  }
  else {
    const int head = get_head_rlen(rc, r_s1, r_e1, r_s2, r_e2);
    if (head >= 0) { //
      /* fragment longer than or equal to the read length, but shorter
       * than twice the read length: this is determined by obtaining the
       * size of the "head" in the diagram below: the portion of end1
       * that is not within [=]. If the read maps to the positive
       * strand, this depends on the reference start of end2 minus the
       * reference start of end1. For negative strand, this is reference
       * start of end1 minus reference start of end2.
       *
       * left                                                 right
       * r_s1                   r_s2   r_e1                   r_e2
       * [------------end1------[======]------end2------------]
       */
      gr.set_name("FRAG_M:" + gr.get_name()); // DEBUG
      truncate_cigar_r(cig1, head);
      read1.resize(cigar_qseq_ops(cig1));
      cig1 += cig2;
      merge_equal_neighbor_cigar_ops(cig1);
      read1 += read2;
    }
    else {
      /* dovetail fragments shorter than read length: this is
       * identified if the above conditions are not satisfied, but
       * there is still some overlap. The overlap will be at the 5'
       * ends of reads, which in theory shouldn't happen unless the
       * two ends are covering identical genomic intervals.
       *
       * left                                           right
       * r_s2             r_s1         r_e2             r_e1
       * [--end2----------[============]----------end1--]
       */
      const int overlap = get_overlap_rlen(rc, r_s1, r_e1, r_s2, r_e2);
      if (overlap > 0) {
        gr.set_name("FRAG_S:" + gr.get_name()); // DEBUG
        truncate_cigar_r(cig1, overlap);
        read1.resize(overlap);
      }

      else return false;
    }
  }
  return true;
}

bool
format_pe(const pe_result &res, const ChromLookup &cl,
          string &read1, string &read2,
          const string &name1, const string &name2,
          string &cig1, string &cig2,
          ofstream &out) {
  uint32_t r_s1 = 0, r_e1 = 0, chr1 = 0;
  uint32_t r_s2 = 0, r_e2 = 0, chr2 = 0;
  const pe_element p = res.best;
  if (res.valid() &&
      !res.ambig()) {

    // PE chromosomes differ or couldn't be found, treat read as unmapped
    if (!chrom_and_posn(cl, cig1, p.r1.pos, r_s1, r_e1, chr1) ||
        !chrom_and_posn(cl, cig2, p.r2.pos, r_s2, r_e2, chr2) ||
        chr1 != chr2) return false;

    revcomp_inplace(read2);
    reverse_cigar(cig2);

    // Select the end points based on orientation, which indicates which
    // end is to the left (first) in the genome. Set the strand and read
    // name based on the first end.
    auto gr = p.rc() ?
      GenomicRegion(cl.names[chr2], r_s2, r_e1, name2, p.diffs(),
                    p.strand()) :
      GenomicRegion(cl.names[chr1], r_s1, r_e2, name1, p.diffs(),
                    p.strand());

    // CIGAR makes dovetail reads no longer overlap, treat as unmapped
    if (!get_pe_overlap(gr, p.rc(), r_s1, r_e1, chr1, r_s2, r_e2, chr2,
                        read1, read2, cig1, cig2)) return false;

    if (p.rc()) { // final revcomp if the first end was a-rich
      //gr.set_strand(gr.get_strand() == '+' ? '-' : '+');
      revcomp_inplace(read1);
      reverse_cigar(cig1);
    }
    if (cig1.empty()) throw runtime_error ("empty cigar\n");
    out << gr << '\t' << read1 << '\t' << cig1 << '\n';
//      << res.flags() << '\n';
  }
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
    return full() && (v[0].diffs == 0 || (v[0].diffs == 1 && seed_number > 0));
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
    if (res.best.valid()) {
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
        << t+tab << "percent_mapped: " << pct(uniq_rds+ambig_rds,
                                              tot_rds == 0 ? 1 : tot_rds) << endl
        << t+tab << "unique: " << uniq_rds << endl
        << t+tab << "percent_unique: " << pct(uniq_rds,
                                              tot_rds == 0 ? 1 : tot_rds) << endl
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
                const se_result &se1, const se_result &se2,
                const string &read1, const string &read2,
                pe_map_stats &pe_stats) {
  pe_stats.update_pair(best);
  if (!best.best.valid()) {
    pe_stats.end1_stats.update(read1, se1);
    pe_stats.end2_stats.update(read2, se2);
  }
}


bool
select_output(const ChromLookup &cl,
              const pe_result &best,
              const se_result &se1, const se_result &se2,
              string &read1, const string &name1,
              string &read2, const string &name2,
              string &cig1, string &cig2,
              ofstream &out) {
  if (format_pe(best, cl, read1, read2, name1, name2, cig1, cig2, out))
    return true;

  // try both ends if no concordant pairs or concordant pair is defective
  const bool good_se1 = format_se(se1, cl, read1, name1, cig1, out),
    good_se2 = format_se(se2, cl, read2, name2, cig2, out);

  return good_se1 || good_se2;
}



inline cmp_t
the_comp(const char a, const char b) {
  return (a & b) == 0;
}

score_t
full_compare(const score_t cutoff,
             Read::const_iterator read_itr,
             const Read::const_iterator &read_end,
             Genome::const_iterator genome_itr) {

  score_t d = 0;
  while (d <= cutoff && read_itr != read_end) {
    d += the_comp(*read_itr, *genome_itr);
    ++read_itr, ++genome_itr;
  }
  return d;
}

template <const uint16_t strand_code,
          class result_type>
void
check_hits(vector<uint32_t>::const_iterator start_idx,
           const vector<uint32_t>::const_iterator end_idx,
           const Read::const_iterator read_start,
           const Read::const_iterator read_end,
           const Genome::const_iterator genome_st,
           result_type &res) {

  for (; start_idx != end_idx && !res.sure_ambig(0); ++start_idx) {
    const score_t diffs = full_compare(res.get_cutoff(),
                                       read_start,
                                       read_end,
                                       genome_st + *start_idx);
    res.update(*start_idx, diffs, strand_code);
  }
}


struct compare_bases {
  compare_bases(const Genome::const_iterator g_) : g(g_) {}
  bool operator()(const uint32_t mid, const uint32_t chr) const {
    return (get_bit_4bit(*(g + mid)) < chr);
  }
  const Genome::const_iterator g;
};


void
find_candidates(const Read::const_iterator read_start,
                const Genome::const_iterator genome_start,
                const uint32_t read_lim, // not necessarily read len
                const uint32_t n_solid_positions,
                vector<uint32_t>::const_iterator &low,
                vector<uint32_t>::const_iterator &high) {
  size_t p = seed::key_weight;
  const size_t lim = std::min(read_lim, n_solid_positions);
  while (p < lim) {
    auto first_1 = lower_bound(low, high, 1, compare_bases(genome_start + p));
    if (get_bit_4bit(*(read_start + p)) == 0) {
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

template <const uint16_t strand_code, class Read,
          class result_type>
void
process_seeds(const uint32_t genome_size,
              const uint32_t max_candidates,
              const AbismalIndex &abismal_index,
              const Genome::const_iterator genome_st,
              const Read &read,
              vector<uint32_t> &hits,
              result_type &res) {
  const uint32_t readlen = read.size();
  const auto read_start(begin(read));
  const auto read_end(end(read));
  const auto index_st(begin(abismal_index.index));
  const auto counter_st(begin(abismal_index.counter));
  const size_t shift_lim =
    readlen > seed::n_seed_positions ? readlen - seed::n_seed_positions : 0;
  const size_t shift = std::max(1ul, shift_lim/(seed::n_shifts - 1));
  uint32_t k;
  bool found_good_seed = false;
  for (uint32_t i = 0; i <= shift_lim && !res.sure_ambig(i); i += shift) {
    hits.clear(); //hits.capacity() == max_candidates
    k = 0;
    get_1bit_hash_4bit(read_start + i, k);

    auto s_idx(index_st + *(counter_st + k));
    auto e_idx(index_st + *(counter_st + k + 1));

    if (s_idx < e_idx) {
      find_candidates(read_start + i, genome_st, readlen - i,
                      seed::n_seed_positions, s_idx, e_idx);
      if (e_idx - s_idx < max_candidates) {
        found_good_seed = true;
        for (; s_idx != e_idx; ++s_idx)
          hits.push_back(*s_idx - i);
      }
    }
    check_hits<strand_code>(begin(hits), end(hits),
                            read_start, read_end, genome_st,  res);
  }

  // All seeds ambiguous, try to use the whole read as seed
  if (!found_good_seed) {
    k = 0;
    get_1bit_hash_4bit(read_start, k);

    auto s_idx(index_st + *(counter_st + k));
    auto e_idx(index_st + *(counter_st + k + 1));

    hits.clear();
    if (s_idx < e_idx) {
      find_candidates(read_start, genome_st,
                      readlen, seed::n_solid_positions, s_idx, e_idx);
      if (e_idx - s_idx >= max_candidates)
        e_idx = s_idx + max_candidates;
      for (; s_idx != e_idx; ++s_idx)
        hits.push_back(*s_idx);
    }
    check_hits<strand_code>(begin(hits), end(hits),
                            read_start, read_end, genome_st, res);
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

static inline score_t
mismatch_score(const char q_base, const uint8_t t_base) {
  return the_comp(q_base, t_base) ?
    align_scores::mismatch : align_scores::match;
}

static inline bool
accept_alignment(const uint32_t len,
                 const score_t cand_diffs,
                 const score_t cur_diffs) {
  return (len >= se_element::min_aligned_length &&
          (cand_diffs <= cur_diffs));
}

template <score_t (*scr_fun)(const char, const uint8_t),
          score_t indel_pen>
void
align_read(se_element &res, string &cigar, const string &read,
            Read &pread, AbismalAlign<scr_fun, indel_pen> &aln) {
  const bool rc = res.rc();
  const bool a_rich = res.a_rich();
  if (res.do_align()) {
    if (rc) {
      const string read_rc(revcomp(read));
      // rc reverses richness of read
      if (a_rich) prep_read<false>(read_rc, pread);
      else prep_read<true>(read_rc, pread);
    }
    else {
      if (a_rich) prep_read<true>(read, pread);
      else prep_read<false>(read, pread);
    }
    uint32_t len; // the region of the read the alignment spans
    string cand_cigar;
    const score_t the_score = aln.align(pread, res.pos, len, cand_cigar);

    // GS the /2 here is specific to the 1, -1, -1 scoring scheme
    const score_t cand_diffs = (static_cast<score_t>(len) - the_score) / 2;
    if (!accept_alignment(len, cand_diffs, res.diffs)) {
      //if (len >= se_element::min_aligned_length)
      //  throw runtime_error("alignment fall through");
      cigar = std::to_string(read.size()) + "M"; // match/mismatch cigar
    }
    else {
      res.diffs = cand_diffs;
      cigar = cand_cigar;
    }
  }

  // Default cigar if alignment was not performed
  else cigar = std::to_string(read.size()) + "M"; // match/mismatch cigar
}

template <const  conversion_type conv>
void
map_single_ended(const bool VERBOSE,
                 const string &reads_file,
                 const size_t batch_size,
                 const size_t max_candidates,
                 const AbismalIndex &abismal_index,
                 se_map_stats &se_stats,
                 ofstream &out) {

  const uint32_t genome_size = abismal_index.genome.size();
  const Genome::const_iterator genome_st(begin(abismal_index.genome));
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
    const size_t n_reads = reads.size();

#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i) res[i].reset();

    const double start_time = omp_get_wtime();
#pragma omp parallel
    {
      Read pread;
      vector<uint32_t> hits;
      hits.reserve(max_candidates);
#pragma omp for
      for (size_t i = 0; i < reads.size(); ++i) {
        if (!reads[i].empty()) {
          prep_read<conv>(reads[i], pread);
          process_seeds<get_strand_code('+', conv)>(genome_size, max_candidates,
                                                    abismal_index, genome_st,
                                                    pread, hits, res[i]);
        }

        if (!reads[i].empty()) {
          const string read_rc(revcomp(reads[i]));
          prep_read<!conv>(read_rc, pread);
          process_seeds<get_strand_code('-', conv)>(genome_size, max_candidates,
                                                    abismal_index, genome_st,
                                                    pread, hits, res[i]);

        }
      }
    }
    total_mapping_time += (omp_get_wtime() - start_time);

#pragma omp parallel
    {
      Read pread;
      string tmp_cigar;
      AbismalAlign<mismatch_score, align_scores::indel>
        aln(genome_st, genome_size);

#pragma omp for
      for (size_t i = 0; i < n_reads; ++i) {
        align_read(res[i].best, cigar[i], reads[i], pread, aln);
        align_read(res[i].second_best, tmp_cigar, reads[i], pread, aln);

      }
    }

    for (size_t i = 0 ; i < n_reads; ++i) {
      se_stats.update(reads[i], res[i]);
      format_se(res[i], abismal_index.cl, reads[i], names[i], cigar[i], out);
    }
  }

  if (VERBOSE) {
    progress.report(cerr, get_filesize(reads_file));
    cerr << "[total mapping time: " << total_mapping_time << endl;
  }
}

void
map_single_ended_rand(const bool VERBOSE,
                      const string &reads_file,
                      const size_t batch_size,
                      const size_t max_candidates,
                      const AbismalIndex &abismal_index,
                      se_map_stats &se_stats,
                      ofstream &out) {

  const uint32_t genome_size = abismal_index.genome.size();
  const Genome::const_iterator genome_st(begin(abismal_index.genome));

  vector<string> names, reads;
  reads.reserve(batch_size);
  names.reserve(batch_size);
  vector<string> cigar(batch_size);
  vector<se_result> res(batch_size);

  ReadLoader rl(reads_file, batch_size);

  ProgressBar progress(get_filesize(reads_file), "mapping reads");
  if (VERBOSE)
    progress.report(cerr, 0);

  AbismalAlign<mismatch_score, align_scores::indel> aln(genome_st, genome_size);

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
      Read pread;
      vector<uint32_t> hits;
      hits.reserve(max_candidates);
#pragma omp for
      for (size_t i = 0; i < reads.size(); ++i)
        if (!reads[i].empty()) {
          prep_read<t_rich>(reads[i], pread);
          process_seeds<get_strand_code('+', t_rich)>(genome_size, max_candidates,
                                                      abismal_index, genome_st,
                                                      pread, hits, res[i]);

          prep_read<a_rich>(reads[i], pread);
          process_seeds<get_strand_code('+', a_rich)>(genome_size, max_candidates,
                                                      abismal_index, genome_st,
                                                      pread, hits, res[i]);

          const string read_rc(revcomp(reads[i]));
          prep_read<t_rich>(read_rc, pread);
          process_seeds<get_strand_code('-', a_rich)>(genome_size, max_candidates,
                                                      abismal_index, genome_st,
                                                      pread, hits, res[i]);

          prep_read<a_rich>(read_rc, pread);
          process_seeds<get_strand_code('-', t_rich)>(genome_size, max_candidates,
                                                      abismal_index, genome_st,
                                                      pread, hits, res[i]);
        }
    }
    total_mapping_time += (omp_get_wtime() - start_time);
#pragma omp parallel
    {
      Read pread;
      string tmp_cigar;
      AbismalAlign<mismatch_score, align_scores::indel>
        aln(genome_st, genome_size);

#pragma omp for
      for (size_t i = 0; i < reads.size(); ++i) {
        align_read(res[i].best, cigar[i], reads[i], pread, aln);
        align_read(res[i].second_best, tmp_cigar, reads[i], pread, aln);
      }
    }


    for (size_t i = 0 ; i < n_reads; ++i) {
      se_stats.update(reads[i], res[i]);
      format_se(res[i], abismal_index.cl, reads[i], names[i], cigar[i], out);
    }
  }
  if (VERBOSE) {
    progress.report(cerr, get_filesize(reads_file));
    cerr << "[total mapping time: " << total_mapping_time << endl;
  }
}


template <const bool cmp,
          const uint16_t strand_code1, const uint16_t strand_code2,
          class result_type>
void
map_pe_batch(const vector<string> &reads1, const vector<string> &reads2,
             const uint32_t max_candidates,
             const AbismalIndex &abismal_index,
             vector<result_type> &res1, vector<result_type> &res2) {
  const uint32_t genome_size = abismal_index.genome.size();
  const Genome::const_iterator genome_st(begin(abismal_index.genome));

#pragma omp parallel
  {
    Read pread;
    vector<uint32_t> hits;
    hits.reserve(max_candidates);

#pragma omp for
    for (size_t i = 0 ; i < reads1.size(); ++i) {
      res1[i].reset();
      res2[i].reset();

      if (!reads1[i].empty()) {
        prep_read<cmp>(reads1[i], pread);
        process_seeds<strand_code1>(genome_size, max_candidates,
                                    abismal_index, genome_st, pread,
                                    hits, res1[i]);
      }

      if (!reads2[i].empty()) {
        const string read_rc(revcomp(reads2[i]));
        prep_read<cmp>(read_rc, pread);
        process_seeds<strand_code2>(genome_size, max_candidates,
                                    abismal_index, genome_st, pread,
                                    hits, res2[i]);
      }
    }
  }
}

static void
best_single(const pe_candidates &pres, se_result &res,
            const string &read, string &cigar,
            const Genome::const_iterator genome_st,
            const uint32_t genome_size) {
  auto lim(begin(pres.v) + pres.sz);
  for (auto i(begin(pres.v)); i != lim; ++i)
    res.update(i->pos, i->diffs, i->flags);
}

template <const bool swap_ends>
static void
best_pair(const pe_candidates &res1, const pe_candidates &res2,
          const string &read1, const string &read2,
          string &cig1, string &cig2,
          const Genome::const_iterator genome_st,
          const uint32_t genome_size,
          pe_result &best) {

  auto j1 = begin(res1.v);
  const auto j1_end = j1 + res1.sz;
  const auto j2_end = begin(res2.v) + res2.sz;
  se_element s1, s2;
  AbismalAlign<mismatch_score, align_scores::indel> aln(genome_st, genome_size);
  Read pread;
  string cand_cig1, cand_cig2;
  for (auto j2(begin(res2.v)); j2 != j2_end; ++j2) {
    s2 = *j2;
    bool aligned_s2 = false;
    const uint32_t lim = j2->pos + read2.length();
    while (j1 != j1_end && j1->pos + pe_element::max_dist < lim) ++j1;
    while (j1 != j1_end && j1->pos + pe_element::min_dist <= lim) {
      s1 = *j1;

      align_read(s1, cand_cig1, read1, pread, aln);
      if (!aligned_s2) {
        align_read(s2, cand_cig2, read2, pread, aln);
        aligned_s2 = true;
      }

      const pe_element p(swap_ends ? s2 : s1, swap_ends ? s1 : s2);
      if (best.update(p)) {
        cig1 = cand_cig1;
        cig2 = cand_cig2;
      }
      ++j1;
    }
  }
}

template <const bool swap_ends>
void
select_maps(const string &read1, const string &read2,
            string &cig1, string &cig2,
            pe_candidates &res1, pe_candidates &res2,
            se_result &res_se1, se_result &res_se2,
            const Genome::const_iterator &genome_st,
            const uint32_t genome_size,
            pe_result &best) {
  res1.prepare_for_mating();
  res2.prepare_for_mating();
  best_pair<swap_ends>(res1, res2, read1, read2, cig1, cig2,
                                  genome_st, genome_size, best);

  // GS: This condition is necessary to not override the cigar
  if (!best.valid()) {
    best_single(res1, res_se1, read1, cig1, genome_st, genome_size);
    best_single(res2, res_se2, read2, cig2, genome_st, genome_size);
  }
}


template <const conversion_type conv>
void
map_paired_ended(const bool VERBOSE,
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

  vector<string> names1(batch_size), reads1(batch_size), cigar1(batch_size),
    names2(batch_size), reads2(batch_size), cigar2(batch_size);
  vector<pe_candidates> res1(batch_size), res2(batch_size);
  vector<pe_result> bests(batch_size);
  vector<se_result> res_se1(batch_size), res_se2(batch_size);

  const uint32_t genome_size = abismal_index.genome.size();
  const Genome::const_iterator genome_st(begin(abismal_index.genome));

  ProgressBar progress(get_filesize(reads_file1), "mapping reads");
  if (VERBOSE)
    progress.report(cerr, 0);

  while (rl1.good() && rl2.good()) {

    if (VERBOSE && progress.time_to_report(rl1.get_current_byte()))
      progress.report(cerr, rl1.get_current_byte());

    rl1.load_reads(names1, reads1);
    rl2.load_reads(names2, reads2);

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
                                                        res1, res2);

#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i)
      select_maps<false>(reads1[i], reads2[i],
                         cigar1[i], cigar2[i],
                         res1[i], res2[i],
                         res_se1[i], res_se2[i],
                         genome_st, genome_size, bests[i]);

    map_pe_batch<!conv,
                 get_strand_code('+', flip_conv(conv)),
                 get_strand_code('-', conv)>(reads2, reads1, max_candidates,
                                             abismal_index, res2, res1);

#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i)
      select_maps<true>(reads2[i], reads1[i],
                        cigar2[i], cigar1[i],
                        res2[i], res1[i],
                        res_se2[i], res_se1[i],
                        genome_st, genome_size, bests[i]);

#pragma omp parallel
    {
      Read pread;
      string tmp_cigar;
      AbismalAlign<mismatch_score, align_scores::indel>
        aln(genome_st, genome_size);

#pragma omp for
      for (size_t i = 0; i < n_reads; ++i) {
        if (!bests[i].valid()) {
          align_read(res_se1[i].best, cigar1[i], reads1[i], pread, aln);
          align_read(res_se1[i].second_best, tmp_cigar, reads1[i], pread, aln);
          align_read(res_se2[i].best, cigar2[i], reads2[i], pread, aln);
          align_read(res_se2[i].second_best, tmp_cigar, reads2[i], pread, aln);
        }
      }
    }

    for (size_t i = 0 ; i < n_reads; ++i)
      if (!select_output(abismal_index.cl, bests[i], res_se1[i], res_se2[i],
                         reads1[i], names1[i], reads2[i], names2[i],
                         cigar1[i], cigar2[i], out)) {
        bests[i].reset();
        res_se1[i].reset();
        res_se2[i].reset();
      }

    for (size_t i = 0 ; i < n_reads; ++i)
      update_pe_stats(bests[i], res_se1[i], res_se2[i], reads1[i],
                      reads2[i], pe_stats);
    total_mapping_time += (omp_get_wtime() - start_time);
  }

  if (VERBOSE) {
    progress.report(cerr, get_filesize(reads_file1));
    cerr << "[total mapping time: " << total_mapping_time << "]" << "\n";
  }
}

void
map_paired_ended_rand(const bool VERBOSE,
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

  vector<string> names1(batch_size), reads1(batch_size), cigar1(batch_size),
    names2(batch_size), reads2(batch_size), cigar2(batch_size);

  vector<pe_candidates> res1(batch_size), res2(batch_size);
  vector<pe_result> bests(batch_size);
  vector<se_result> res_se1(batch_size), res_se2(batch_size);

  // alignment stuff
  const uint32_t genome_size = abismal_index.genome.size();
  const Genome::const_iterator genome_st(begin(abismal_index.genome));

  AbismalAlign<mismatch_score, align_scores::indel> aln(genome_st, genome_size);

  ProgressBar progress(get_filesize(reads_file1), "mapping reads");
  if (VERBOSE)
    progress.report(cerr, 0);

  while (rl1.good() && rl2.good()) {

    if (VERBOSE && progress.time_to_report(rl1.get_current_byte()))
      progress.report(cerr, rl1.get_current_byte());

    rl1.load_reads(names1, reads1);
    rl2.load_reads(names2, reads2);

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
                                               abismal_index, res1, res2);
#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i)
      select_maps<false>(reads1[i], reads2[i],
                         cigar1[i], cigar2[i],
                         res1[i], res2[i],
                         res_se1[i], res_se2[i],
                         genome_st, genome_size, bests[i]);

    // t-rich end1, neg-strand end1
    map_pe_batch<a_rich,
                 get_strand_code('+', a_rich),
                 get_strand_code('-', t_rich)>(reads2, reads1, max_candidates,
                                               abismal_index, res2, res1);
#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i)
      select_maps<true>(reads2[i], reads1[i],
                        cigar2[i], cigar1[i],
                        res2[i], res1[i],
                        res_se2[i], res_se1[i],
                        genome_st, genome_size, bests[i]);

    // a-rich end1, pos-strand end1
    map_pe_batch<a_rich,
                 get_strand_code('+', a_rich),
                 get_strand_code('-', t_rich)>(reads1, reads2, max_candidates,
                                               abismal_index, res1, res2);

#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i)
      select_maps<false>(reads1[i], reads2[i],
                         cigar1[i], cigar2[i],
                         res1[i], res2[i],
                         res_se1[i], res_se2[i],
                         genome_st, genome_size, bests[i]);

    // a-rich end1, neg-strand end1
    map_pe_batch<t_rich,
                 get_strand_code('+', t_rich),
                 get_strand_code('-', a_rich)>(reads2, reads1, max_candidates,
                                               abismal_index, res2, res1);


#pragma omp parallel for
    for (size_t i = 0 ; i < n_reads; ++i)
      select_maps<true>(reads2[i], reads1[i],
                        cigar2[i], cigar1[i],
                        res2[i], res1[i],
                        res_se2[i], res_se1[i],
                        genome_st, genome_size, bests[i]);

#pragma omp parallel
    {
      Read pread;
      string tmp_cigar;
      AbismalAlign<mismatch_score, align_scores::indel>
        aln(genome_st, genome_size);

#pragma omp for
      for (size_t i = 0; i < n_reads; ++i) {
        if (!bests[i].valid()) {
          align_read(res_se1[i].best, cigar1[i], reads1[i], pread, aln);
          align_read(res_se1[i].second_best, tmp_cigar, reads1[i], pread, aln);
          align_read(res_se2[i].best, cigar2[i], reads2[i], pread, aln);
          align_read(res_se2[i].second_best, tmp_cigar, reads2[i], pread, aln);
        }
      }
    }

    for (size_t i = 0 ; i < n_reads; ++i)
      if (!select_output(abismal_index.cl, bests[i], res_se1[i], res_se2[i],
                         reads1[i], names1[i], reads2[i], names2[i],
                         cigar1[i], cigar2[i], out)) {
        bests[i].reset();
        res_se1[i].reset();
        res_se2[i].reset();
      }

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


int main(int argc, const char **argv) {
  try {
    string index_file;
    string outfile;
    bool VERBOSE = false;
    bool GA_conversion = false;
    bool allow_ambig = false;
    bool pbat_mode = false;
    bool random_pbat = false;
    uint32_t max_candidates = 3000;
    size_t max_diffs = 6;
    size_t invalid_hit_diffs = 40;
    size_t batch_size = 100000;
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
    opt_parse.add_opt("hits", 'h', "max number of mismatches for a hit to be "
                                   "considered for subsequent alignment",
                      false, invalid_hit_diffs);
    opt_parse.add_opt("shifts", 's', "number of seed shifts",
                      false, seed::n_shifts);
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

    se_element::max_diffs = max_diffs;
    se_element::invalid_hit_diffs = max(max_diffs, invalid_hit_diffs);

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
    pe_map_stats pe_stats(pe_element::min_dist, pe_element::max_dist);

    std::ofstream out(outfile);
    if (!out)
      throw runtime_error("failed to open output file: " + outfile);

    if (reads_file2.empty()) {
      if (GA_conversion || pbat_mode)
        map_single_ended<a_rich>(VERBOSE, reads_file, batch_size,
                                 max_candidates, abismal_index,
                                 se_stats, out);
      else if (random_pbat)
        map_single_ended_rand(VERBOSE, reads_file, batch_size, max_candidates,
                              abismal_index, se_stats, out);
      else
        map_single_ended<t_rich>(VERBOSE, reads_file, batch_size,
                                 max_candidates, abismal_index,
                                 se_stats, out);
    }
    else {
      if (pbat_mode)
        map_paired_ended<a_rich>(VERBOSE, reads_file, reads_file2,
                                          batch_size, max_candidates,
                                          abismal_index, pe_stats, out);
      else if (random_pbat)
        map_paired_ended_rand(VERBOSE, reads_file, reads_file2,
                              batch_size, max_candidates,
                              abismal_index, pe_stats, out);
      else
        map_paired_ended<t_rich>(VERBOSE, reads_file, reads_file2,
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
