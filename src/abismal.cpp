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

typedef uint16_t flags_t; // every bit is a flag
typedef int16_t score_t; // aln score, edit distance, hamming distance
typedef vector<uint8_t> Read; //4-bit encoding of reads
typedef vector<size_t> PackedRead; //4-bit encoding of reads

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
  igzfstream *in;

  static const size_t batch_size;
  static const uint32_t min_read_length;
};

const size_t ReadLoader::batch_size = 2000;
const uint32_t ReadLoader::min_read_length = seed::n_seed_positions;

inline void
update_max_read_length(size_t &max_length, const vector<string> &reads) {
  for (auto it (begin(reads)); it != end(reads); ++it)
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

inline bool
valid(const se_element s, const uint32_t readlen) {
  return s.diffs <= static_cast<score_t>(se_element::valid_frac*readlen);
}

inline bool
valid_hit(const se_element s, const uint32_t readlen) {
  return s.diffs < static_cast<score_t>(se_element::invalid_hit_frac*readlen);
}

double se_element::valid_frac = 0.1;
const double se_element::invalid_hit_frac = 0.4;

struct se_candidates {
  se_candidates () : sz(1), best(se_element()), v(vector<se_element>(max_size)) {}
  inline bool full() const { return sz == max_size; };
  void update_exact_match(const uint32_t p, const score_t d, const flags_t s) {
    if (d == 0) {
      if (best.empty())
        best = se_element(p, d, s); // s has ambig flag set to false

      else if ((p != best.pos || s != best.flags))
        best.set_ambig();
    }
  }

  void update_cand(const uint32_t p, const score_t d, const flags_t s) {
    if (d < v.front().diffs) {
      if (full()) {
        pop_heap(begin(v), end(v));
        v.back() = se_element(p, d, s);
        push_heap(begin(v), end(v));
      }
      else {
        v[sz++] = se_element(p, d, s);
        push_heap(begin(v), begin(v) + sz);
      }
    }
  }

  void update(const uint32_t p, const score_t d, const flags_t s) {
    update_exact_match(p, d, s);
    update_cand(p, d, s);
  }

  inline bool sure_ambig() const {
    return best.sure_ambig();
  }
  void reset(const uint32_t readlen)  {
    v.front().reset(readlen);
    sz = 1;
    best.reset(readlen);
  }

  void prepare_for_alignments() {
   sort(begin(v), begin(v) + sz, // no sort_heap here as heapify used "diffs"
         [](const se_element &a, const se_element &b) {
           return (a.pos < b.pos) || (a.pos == b.pos && a.flags < b.flags);
         });
    sz = unique(begin(v), begin(v) + sz) - begin(v);
  }

  uint32_t get_cutoff() const { return v.front().diffs; }

  uint32_t sz;
  se_element best;
  vector<se_element> v;
  static const uint32_t max_size;
};

const uint32_t se_candidates::max_size = 100;

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
          ostream &out) {
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

  out << sr << "\n";
  return ambig ? map_ambig : map_unique;
}

struct pe_element {
  pe_element() : r1(se_element()), r2(se_element()) {}
  pe_element(const se_element s1, const se_element s2) : r1(s1), r2(s2) {}

  score_t diffs() const { return r1.diffs + r2.diffs; }
  void reset() {
    r1.reset();
    r2.reset();
  }
  inline void update(const se_element s1, const se_element s2) {
    r1 = s2;
    r2 = s2;
  }

  inline bool ambig() const { return r1.ambig(); }
  inline bool empty() const { return r1.empty(); }
  inline void set_ambig() { r1.set_ambig(); }

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
          ostream &out) {
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
  pe_candidates() : v(vector<se_element>(max_size)), sz(1) {}
  bool full() const {return sz == max_size;}
  void reset(const uint32_t readlen) {
    v.front().reset(readlen); sz = 1;
  }
  score_t get_cutoff() const { return v.front().diffs; }
  void update(const uint32_t p, const score_t d, const flags_t s) {
    if (d < v.front().diffs) {
      if (full()) {
        pop_heap(begin(v), end(v));
        v.back() = se_element(p, d, s);
        push_heap(begin(v), end(v));
      }
      else {
        v[sz++] = se_element(p, d, s);
        push_heap(begin(v), begin(v) + sz);
      }
    }
  }
  bool sure_ambig() const {
    return full() && (v.front().diffs == 0);
  }
  void prepare_for_mating() {
    sort(begin(v), begin(v) + sz, // no sort_heap here as heapify used "diffs"
         [](const se_element &a, const se_element &b) {return a.pos < b.pos;});
    sz = unique(begin(v), begin(v) + sz) - begin(v);
  }
  vector<se_element> v;
  uint32_t sz;
  static const uint32_t max_size;
};

const uint32_t pe_candidates::max_size = 1000;

inline double pct(const double a, const double b) {return 100.0*a/b;}

struct se_map_stats {
  se_map_stats() :
    tot_rds(0), uniq_rds(0), ambig_rds(0), unmapped_rds(0), skipped_rds(0) {}
  uint32_t tot_rds;
  uint32_t uniq_rds;
  uint32_t ambig_rds;
  uint32_t unmapped_rds;
  uint32_t skipped_rds;

  void update(const se_element &s, const bool skipped) {
    ++tot_rds;
    const bool valid = !s.empty();
    const bool ambig = s.ambig();
    uniq_rds += (valid && !ambig);
    ambig_rds += (valid && ambig);
    unmapped_rds += !valid;
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

  void update(const bool allow_ambig, const pe_element &p,
              const se_element &s1, const bool skipped_se1,
              const se_element &s2, const bool skipped_se2) {
    const bool valid = !p.empty();
    const bool ambig = p.ambig();
    ++tot_pairs;
    ambig_pairs += (valid && ambig);
    uniq_pairs += (valid && !ambig);

    if (p.empty() || (!allow_ambig && p.ambig())) {
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
              const string &read1, const string &name1,
              const string &read2, const string &name2,
              const string &cig1, const string &cig2,
              pe_element &best, se_element &se1, se_element &se2,
              ostream &out) {
  const map_type pe_map_type = format_pe(allow_ambig, best, cl,
                                         read1, read2, name1,
                                         name2, cig1, cig2, out);
  if (pe_map_type == map_unmapped ||
      (!allow_ambig && pe_map_type == map_ambig)) {
    // GS: do not report in mapstats a read that was not reported
    if (pe_map_type == map_unmapped)
      best.reset();

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
  score_t d = 0;
  while (d < cutoff && read_itr != read_end) {
    d += 16 - __builtin_popcountll((*read_itr) & 
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
    /* GS: the_pos & 15u tells if the position is a multiple of 16, in
     * which case it is aligned with the genome. Otherwise we need to
     * use the unaligned comparison function that offsets genome
     * position by the_pos (mod 16). Multiplied by 4 because each base
     * uses 4 bits */
    const score_t diffs = full_compare(
      res.get_cutoff(), read_end, ((the_pos & 15u) << 2),
      read_st, genome_st + (the_pos >> 4)
    );
    res.update(the_pos, diffs, strand_code);
  }
}

static void
get_minimizer_offsets(const uint32_t readlen,
                      Read::const_iterator read_start,
                      vector<kmer_loc> &offsets,
                      deque<kmer_loc> &window_kmers) {
  const uint32_t shift_lim = (readlen >= seed::n_seed_positions) ? 
                             (readlen - seed::n_seed_positions) : 0;
  size_t kmer = 0;

  offsets.clear();
  window_kmers.clear();

  get_1bit_hash(read_start, kmer);
  add_kmer<seed::w_map>(kmer_loc(kmer, 0), window_kmers);
  read_start += seed::key_weight;

  for (size_t i = 1; i < seed::w_map; ++i) {
    shift_hash_key(*read_start++, kmer);
    add_kmer<seed::w_map>(kmer_loc(kmer, i), window_kmers);
  }

  for (size_t i = seed::w_map; i <= shift_lim; ++i) {
    if (window_kmers.front() != offsets.back())
      offsets.push_back(window_kmers.front());
    shift_hash_key(*read_start++, kmer);
    add_kmer<seed::w_map>(kmer_loc(kmer, i), window_kmers);
  }

  if (window_kmers.front() != offsets.back())
    offsets.push_back(window_kmers.front());
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
              const Read &read_seed,
              const PackedRead &packed_read,
              vector<kmer_loc> &offsets,
              deque<kmer_loc> &window_kmers,
              result_type &res) {
  const uint32_t readlen = read_seed.size();
  const auto read_start(begin(read_seed));
  const auto packed_read_start(begin(packed_read));
  const auto packed_read_end(end(packed_read));

  size_t k = 0;
  get_1bit_hash(read_start, k);
  // specific step: look for exact matches, which can be at most
  // w_index bases apart
  for (uint32_t j = 0; j < seed::w_index; ++j) {
    auto s_idx(index_st + *(counter_st + k));
    auto e_idx(index_st + *(counter_st + k + 1));
    if (s_idx < e_idx) {
      find_candidates<seed::n_sorting_positions>(
        read_start + j, genome_st, readlen - j, s_idx, e_idx
      );
      if (e_idx - s_idx <= max_candidates) {
        check_hits<strand_code>(
          j, packed_read_start, packed_read_end,
          genome_st.itr, e_idx, s_idx, res
        );
      }
    }

    // GS: this needs to be more readable
    shift_hash_key(*(read_start + seed::key_weight + j), k);
  }

  // sensitive step: get positions using minimizers
  get_minimizer_offsets(readlen, read_start, offsets, window_kmers);

  const auto offset_end = end(offsets);
  size_t num_kmers = 0;
  for (auto it(begin(offsets)); it != offset_end; ++it) {
    auto s_idx(index_st + *(counter_st + it->kmer));
    auto e_idx(index_st + *(counter_st + it->kmer + 1));
    if (s_idx < e_idx) {
      find_candidates<seed::n_seed_positions>(
        read_start + it->loc, genome_st, readlen - it->loc, s_idx, e_idx
      );
      if ((e_idx - s_idx) <= max_candidates) {
        ++num_kmers;
        check_hits<strand_code>(
          it->loc, packed_read_start, packed_read_end,
          genome_st.itr, e_idx, s_idx, res
        );
      }
      if (num_kmers == seed::max_seeds) return;
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
  auto it(begin(packed_pread));

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

  const score_t ans = aln.align(pread, res.pos, len, cigar);
  res.diffs = simple_aln::edit_distance(ans, len, cigar);
  return ans;
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
align_se_candidates(const string &read,
                    se_candidates &res, se_element &best, 
                    string &cigar, Read &pread,
                    AbismalAlignSimple &aln) {
  /* GS: this is faster, but potentially prevents ends to be
   * soft-clipped, which is helpful on reads with bad ends that were
   * not trimmed prior to alignment */
  if (res.best.diffs <= simple_aln::min_diffs_to_align) {
    best = res.best;
    simple_aln::make_default_cigar(read.size(), cigar);
    return;
  }

  const score_t readlen = static_cast<score_t>(read.size());

  // variables to store the best alignment score
  score_t best_score = 0;
  uint32_t len = 0;
  se_element s;
  string cand_cigar;

  res.prepare_for_alignments();

  const auto lim(begin(res.v) + res.sz);
  for (auto it(begin(res.v)); it != lim; ++it) {
    s = *it;
    if (valid_hit(s, readlen)) {
      const score_t scr = align_read(s, read, pread, len, cand_cigar, aln);
      if (valid(s, len)) {
        if (scr > best_score) {
          best = s; // GS: ambig is unset on s
          best_score = scr;
          cigar = cand_cigar;
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
  }
}

template <const  conversion_type conv>
inline void
map_single_ended(const bool VERBOSE, const bool allow_ambig,
                 const size_t max_candidates,
                 const AbismalIndex &abismal_index, ReadLoader &rl,
                 se_map_stats &se_stats, ostream &out,
                 ProgressBar &progress) {
  const auto counter_st(begin(abismal_index.counter));
  const auto index_st(begin(abismal_index.index));
  const genome_iterator genome_st(begin(abismal_index.genome));

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
  Read pread;
  PackedRead packed_pread;
  vector<kmer_loc> offsets;
  deque<kmer_loc> window_kmers;
  se_candidates res;

  size_t the_byte;
  while (rl.good()) {
#pragma omp critical
    {
      rl.load_reads(names, reads);
      the_byte = rl.get_current_byte();
    }

    size_t max_batch_read_length = 0;
    update_max_read_length(max_batch_read_length, reads);

    AbismalAlignSimple aln(genome_st, max_batch_read_length);
    offsets.reserve(max_batch_read_length);

    const size_t n_reads = reads.size();
    for (size_t i = 0; i < n_reads; ++i) {
      res.reset(reads[i].size());
      bests[i].reset();
      if (!reads[i].empty()) {
        prep_read<conv>(reads[i], pread);
        pack_read(pread, packed_pread);
        process_seeds<get_strand_code('+', conv)>(max_candidates,
          counter_st, index_st, genome_st, pread,
          packed_pread, offsets, window_kmers, res
        );

        const string read_rc(revcomp(reads[i]));
        prep_read<!conv>(read_rc, pread);
        pack_read(pread, packed_pread);
        process_seeds<get_strand_code('-', conv)>(max_candidates,
          counter_st, index_st, genome_st, pread,
          packed_pread, offsets, window_kmers, res
        );
        align_se_candidates(reads[i], res, bests[i], cigar[i], pread, aln);
      }
    }

#pragma omp critical
    {
      for (size_t i = 0; i < n_reads; ++i) {
        if (format_se(allow_ambig, bests[i], abismal_index.cl, reads[i],
              names[i], cigar[i], out) == map_unmapped)
          bests[i].reset();
        se_stats.update(bests[i], reads[i].length() == 0);
      }
    }
#pragma omp critical
    {
      if (VERBOSE && progress.time_to_report(the_byte))
        progress.report(cerr, the_byte);
    }
  }
}

inline void
map_single_ended_rand(const bool VERBOSE, const bool allow_ambig,
                      const size_t max_candidates,
                      const AbismalIndex &abismal_index, ReadLoader &rl,
                      se_map_stats &se_stats, ostream &out,
                      ProgressBar &progress) {
  const auto counter_st(begin(abismal_index.counter));
  const auto index_st(begin(abismal_index.index));
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
  Read pread;
  PackedRead packed_pread;
  se_candidates res;
  vector<kmer_loc> offsets;
  deque<kmer_loc> window_kmers;

  size_t the_byte;
  while (rl.good()) {
#pragma omp critical
    {
      rl.load_reads(names, reads);
      the_byte = rl.get_current_byte();
    }

    size_t max_batch_read_length = 0;
    update_max_read_length(max_batch_read_length, reads);

    AbismalAlignSimple aln(genome_st, max_batch_read_length);
    offsets.reserve(max_batch_read_length);

    const size_t n_reads = reads.size();
    for (size_t i = 0; i < n_reads; ++i) {
      res.reset(reads[i].size());
      bests[i].reset();
      if (!reads[i].empty()) {
        prep_read<t_rich>(reads[i], pread);
        pack_read(pread, packed_pread);
        process_seeds<get_strand_code('+', t_rich)>(max_candidates,
          counter_st, index_st, genome_st, pread,
          packed_pread, offsets, window_kmers, res
        );

        prep_read<a_rich>(reads[i], pread);
        pack_read(pread, packed_pread);
        process_seeds<get_strand_code('+', a_rich)>(max_candidates,
          counter_st, index_st, genome_st, pread,
          packed_pread, offsets, window_kmers, res
        );

        const string read_rc(revcomp(reads[i]));
        prep_read<t_rich>(read_rc, pread);
        pack_read(pread, packed_pread);
        process_seeds<get_strand_code('-', a_rich)>(max_candidates,
          counter_st, index_st, genome_st, pread, packed_pread,
          offsets, window_kmers, res
        );

        prep_read<a_rich>(read_rc, pread);
        pack_read(pread, packed_pread);
        process_seeds<get_strand_code('-', t_rich)>(max_candidates,
          counter_st, index_st, genome_st, pread,
          packed_pread, offsets, window_kmers, res
        );

        align_se_candidates(reads[i], res, bests[i], cigar[i], pread, aln);
      }
    }
#pragma omp critical
    {
      for (size_t i = 0; i < n_reads; ++i) {
        if (format_se(allow_ambig, bests[i], abismal_index.cl, reads[i],
              names[i], cigar[i], out) == map_unmapped)
          bests[i].reset(reads[i].size());
        se_stats.update(bests[i], reads[i].length() == 0);
      }
    }
#pragma omp critical
    {
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
                 const size_t max_candidates,
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
      map_single_ended_rand(VERBOSE, allow_ambig, max_candidates,
        abismal_index, rl, se_stats, out, progress);
    else
      map_single_ended<conv>(VERBOSE, allow_ambig, max_candidates,
          abismal_index, rl, se_stats, out, progress);
  }
  if (VERBOSE) {
    cerr << "[total mapping time: " << omp_get_wtime() - start_time << "]"
         << endl;
  }
}
static void
best_single(const pe_candidates &pres, se_candidates &res) {
  const auto lim(begin(pres.v) + pres.sz);
  for (auto i(begin(pres.v)); i != lim && !res.sure_ambig(); ++i) {
    res.update(i->pos, i->diffs, i->flags);
  }
}

template <const bool swap_ends>
static void
best_pair(const pe_candidates &res1, const pe_candidates &res2,
          const Read &pread1, const Read &pread2,
          score_t &aln_score, string &cig1, string &cig2,
          AbismalAlignSimple &aln, pe_element &best) {
  auto j1 = begin(res1.v);
  const auto j1_end = j1 + res1.sz;
  const auto j2_end = begin(res2.v) + res2.sz;

  uint32_t len1 = 0, len2 = 0;
  se_element s1, s2;
  string cand_cig1, cand_cig2;

  for (auto j2(begin(res2.v)); j2 != j2_end; ++j2) {
    s2 = *j2;
    if (valid_hit(s2, pread2.size())) {
      const uint32_t unaligned_lim = s2.pos + pread2.size();
      for (j1 = begin(res1.v); j1 != j1_end &&
           j1->pos + pe_element::max_dist < unaligned_lim; ++j1);

      score_t scr2 = 0;
      uint32_t aligned_lim = 0;
      while (j1 != j1_end && j1->pos + pe_element::min_dist <= unaligned_lim) {
        s1 = *j1;
        if (valid_hit(s1, pread1.size())) {
          const score_t scr1 = align_read(s1, pread1, len1, cand_cig1, aln);

          // GS: guarantees that j2 is aligned only once
          if (scr2 == 0) {
            scr2 = align_read(s2, pread2, len2, cand_cig2, aln);
            aligned_lim = s2.pos + cigar_rseq_ops(cand_cig2);
          }

          // GS: only accept if length post alignment is still within limits
          if (valid(s1, len1) && valid(s2, len2) &&
              (s1.pos + pe_element::max_dist >= aligned_lim) &&
              (s1.pos + pe_element::min_dist <= aligned_lim)) {
            const score_t scr = scr1 + scr2;
            if (scr > aln_score) {
              aln_score = scr;
              cig1 = cand_cig1;
              cig2 = cand_cig2;

              // GS: unsets ambig
              best.update(swap_ends ? s2 : s1, swap_ends ? s1 : s2);
            }
            else if ((scr == aln_score) &&
                     (((swap_ends ? s2 : s1) != best.r1) || 
                      ((swap_ends ? s1 : s2) != best.r2)))
              best.set_ambig();
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
            score_t &aln_score, string &cig1, string &cig2,
            pe_candidates &res1, pe_candidates &res2,
            se_candidates &res_se1, se_candidates &res_se2,
            AbismalAlignSimple &aln,  pe_element &best) {
  res1.prepare_for_mating();
  res2.prepare_for_mating();

  best_pair<swap_ends>(
    res1, res2, pread1, pread2, aln_score, cig1, cig2, aln, best
  );

  best_single(res1, res_se1);
  best_single(res2, res_se2);
}

template <const bool cmp, const bool swap_ends,
          const uint16_t strand_code1, const uint16_t strand_code2>
inline void
map_fragments(const string &read1, const string &read2,
              const uint32_t max_candidates,
              const vector<uint32_t>::const_iterator counter_st,
              const vector<uint32_t>::const_iterator index_st,
              const genome_iterator genome_st,
              Read &pread1, Read &pread2, PackedRead &packed_pread,
              score_t &aln_score, string &cigar1, string &cigar2,
              AbismalAlignSimple &aln,
              vector<kmer_loc> &offsets,
              deque<kmer_loc> &window_kmers,
              pe_candidates &res1, pe_candidates &res2,
              se_candidates &res_se1, se_candidates &res_se2,
              pe_element &best) {
  res1.reset(read1.size());
  res2.reset(read2.size());

  if (!read1.empty()) {
    prep_read<cmp>(read1, pread1);
    pack_read(pread1, packed_pread);
    process_seeds<strand_code1>(max_candidates, counter_st,
      index_st, genome_st, pread1, packed_pread, offsets,
      window_kmers, res1
    );
  }

  if (!read2.empty()) {
    const string read_rc(revcomp(read2));
    prep_read<cmp>(read_rc, pread2);
    pack_read(pread2, packed_pread);
    process_seeds<strand_code2>(max_candidates, counter_st,
      index_st, genome_st, pread2, packed_pread, offsets,
      window_kmers, res2
    );
  }

  select_maps<swap_ends>(pread1, pread2, aln_score, cigar1, cigar2,
                         res1, res2, res_se1, res_se2, aln, best);
}

template <const conversion_type conv>
inline void
map_paired_ended(const bool VERBOSE,
                 const bool allow_ambig,
                 const size_t max_candidates, const AbismalIndex &abismal_index,
                 ReadLoader &rl1, ReadLoader &rl2,
                 pe_map_stats &pe_stats, ostream &out,
                 ProgressBar &progress) {
  const auto counter_st(begin(abismal_index.counter));
  const auto index_st(begin(abismal_index.index));
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
  score_t aln_score;
  Read pread1, pread2;
  PackedRead packed_pread;
  pe_candidates res1;
  pe_candidates res2;
  se_candidates res_se1;
  se_candidates res_se2;
  vector<kmer_loc> offsets;
  deque<kmer_loc> window_kmers;

  size_t the_byte;
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
    offsets.reserve(max_batch_read_length);

    const size_t n_reads = reads1.size();
    for (size_t i = 0 ; i < n_reads; ++i) {
      res1.reset(reads1[i].size());
      res2.reset(reads2[i].size());
      res_se1.reset(reads1[i].size());
      res_se2.reset(reads2[i].size());

      bests[i].reset();
      bests_se1[i].reset();
      bests_se2[i].reset();
      aln_score = 0;

      map_fragments<conv, false,
                   get_strand_code('+',conv),
                   get_strand_code('-', flip_conv(conv))>(
        reads1[i], reads2[i], max_candidates, counter_st,
        index_st, genome_st, pread1, pread2, packed_pread, aln_score,
        cigar1[i], cigar2[i], aln, offsets, window_kmers, res1, res2,
        res_se1, res_se2, bests[i]
      );

      map_fragments<!conv, true,
                   get_strand_code('+', flip_conv(conv)),
                   get_strand_code('-', conv)>(
        reads2[i], reads1[i], max_candidates, counter_st,
        index_st, genome_st, pread1, pread2, packed_pread, aln_score,
        cigar2[i], cigar1[i], aln, offsets, window_kmers, res2, res1,
        res_se2, res_se1, bests[i]
      );

      if (bests[i].empty() || (!allow_ambig && bests[i].ambig())) {
        align_se_candidates(reads1[i], res_se1, bests_se1[i], cigar1[i], pread1, aln);
        align_se_candidates(reads2[i], res_se2, bests_se2[i], cigar2[i], pread2, aln);
      }
    }

#pragma omp critical
    {
      for (size_t i = 0; i < n_reads; ++i) {
        select_output(allow_ambig, abismal_index.cl,
            reads1[i], names1[i], reads2[i], names2[i], cigar1[i], cigar2[i],
            bests[i], bests_se1[i], bests_se2[i], out
          );
        pe_stats.update(allow_ambig, bests[i],
                        bests_se1[i], reads1[i].length() == 0,
                        bests_se2[i], reads2[i].length() == 0);
      }
    }
#pragma omp critical
    {
      if (VERBOSE && progress.time_to_report(the_byte))
        progress.report(cerr, the_byte);
    }
  }
}

inline void
map_paired_ended_rand(const bool VERBOSE, const bool allow_ambig,
                      const size_t max_candidates,
                      const AbismalIndex &abismal_index,
                      ReadLoader &rl1, ReadLoader &rl2,
                      pe_map_stats &pe_stats, ostream &out,
                      ProgressBar &progress) {
  const auto counter_st(begin(abismal_index.counter));
  const auto index_st(begin(abismal_index.index));
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

  score_t aln_score;
  Read pread1, pread2;
  PackedRead packed_pread;
  pe_candidates res1;
  pe_candidates res2;
  se_candidates res_se1;
  se_candidates res_se2;
  vector<kmer_loc> offsets;
  deque<kmer_loc> window_kmers;


  size_t the_byte;
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
    offsets.reserve(max_batch_read_length);

    const size_t n_reads = reads1.size();
    for (size_t i = 0 ; i < n_reads; ++i) {
      res1.reset(reads1[i].size());
      res2.reset(reads2[i].size());
      res_se1.reset(reads1[i].size());
      res_se2.reset(reads2[i].size());

      bests[i].reset();
      bests_se1[i].reset();
      bests_se2[i].reset();
      aln_score = 0;

      // GS: (1) T/A-rich +/- strand
      map_fragments<t_rich, false,
                   get_strand_code('+', t_rich),
                   get_strand_code('-', a_rich)>(
        reads1[i], reads2[i], max_candidates, counter_st,
        index_st, genome_st, pread1, pread2, packed_pread, aln_score,
        cigar1[i], cigar2[i], aln, offsets, window_kmers, res1, res2,
        res_se1, res_se2, bests[i]
      );

      // GS: (2) T/A-rich, -/+ strand
      map_fragments<a_rich, true,
                   get_strand_code('+', a_rich),
                   get_strand_code('-', t_rich)>(
        reads2[i], reads1[i], max_candidates, counter_st,
        index_st, genome_st, pread2, pread1, packed_pread, aln_score,
        cigar2[i], cigar1[i], aln, offsets, window_kmers, res2, res1,
        res_se2, res_se1, bests[i]
      );

      // GS: (3) A/T-rich +/- strand
      map_fragments<a_rich, false,
                   get_strand_code('+', a_rich),
                   get_strand_code('-', t_rich)>(
        reads1[i], reads2[i], max_candidates, counter_st,
        index_st, genome_st, pread1, pread2, packed_pread, aln_score,
        cigar1[i], cigar2[i], aln, offsets, window_kmers, res1, res2,
        res_se1, res_se2, bests[i]
      );

      // GS: (4) A/T-rich, -/+ strand
      map_fragments<t_rich, true,
                   get_strand_code('+', t_rich),
                   get_strand_code('-', a_rich)>(
        reads2[i], reads1[i], max_candidates, counter_st,
        index_st, genome_st, pread2, pread1, packed_pread, aln_score,
        cigar2[i], cigar1[i], aln, offsets, window_kmers, res2, res1,
        res_se2, res_se1, bests[i]
      );

      // GS: align best SE candidates if no concordant pairs found
      if (bests[i].empty() || bests[i].ambig()) {
        align_se_candidates(reads1[i], res_se1, bests_se1[i], cigar1[i], pread1, aln);
        align_se_candidates(reads2[i], res_se2, bests_se2[i], cigar2[i], pread2, aln);
      }
    }

#pragma omp critical
    {
      for (size_t i = 0; i < n_reads; ++i) {
        select_output(allow_ambig, abismal_index.cl,
          reads1[i], names1[i], reads2[i], names2[i], cigar1[i], cigar2[i], 
          bests[i], bests_se1[i], bests_se2[i], out
        );

        pe_stats.update(allow_ambig, bests[i],
                        bests_se1[i], reads1[i].length() == 0,
                        bests_se2[i], reads2[i].length() == 0);
      }
    }
#pragma omp critical
    {
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
                 const size_t max_candidates,
                 const AbismalIndex &abismal_index,
                 pe_map_stats &pe_stats,
                 ostream &out) {
  const auto counter_st(begin(abismal_index.counter));
  const auto index_st(begin(abismal_index.index));
  const genome_iterator genome_st(begin(abismal_index.genome));

  ReadLoader rl1(reads_file1);
  ReadLoader rl2(reads_file2);
  ProgressBar progress(get_filesize(reads_file1), "mapping reads");

  if (VERBOSE)
    progress.report(cerr, 0);

  double start_time = omp_get_wtime();

#pragma omp parallel for
  for (int i = 0; i < omp_get_num_threads(); ++i) {
    if (random_pbat)
      map_paired_ended_rand(VERBOSE, allow_ambig, max_candidates,
          abismal_index, rl1, rl2, pe_stats, out, progress);

    else
      map_paired_ended<conv>(VERBOSE, allow_ambig, max_candidates,
          abismal_index, rl1, rl2, pe_stats, out, progress);
  }

  if (VERBOSE) {
    cerr << "[total mapping time: " << omp_get_wtime() - start_time
         << "]" << endl;
  }
}

static void
select_max_candidates(const uint32_t genome_size,
                      uint32_t &max_candidates) {
  static const double genome_frac = 1.5e-5;
  static const uint32_t min_max_candidates = 100u;

  const uint32_t c = static_cast<uint32_t>(genome_size*genome_frac);
  max_candidates = max(c, min_max_candidates);
}

int main(int argc, const char **argv) {

  try {
    static const string ABISMAL_VERSION = "0.2.0";
    bool VERBOSE = false;
    bool GA_conversion = false;
    bool allow_ambig = false;
    bool pbat_mode = false;
    bool random_pbat = false;
    uint32_t max_candidates = 0;
    int n_threads = 1;
    string index_file;
    string outfile;
    string stats_outfile = "";

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
    opt_parse.add_opt("candidates", 'c', "max candidates for full comparison",
                      false, max_candidates);
    opt_parse.add_opt("min-frag", 'l', "min fragment size (pe mode)",
                      false, pe_element::min_dist);
    opt_parse.add_opt("max-frag", 'L', "max fragment size (pe mode)",
                      false, pe_element::max_dist);
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
          max_candidates, abismal_index, se_stats, out);
      else if (random_pbat)
        run_single_ended<t_rich, true>(VERBOSE, allow_ambig, reads_file,
             max_candidates, abismal_index, se_stats, out);
      else
        run_single_ended<t_rich, false>(VERBOSE, allow_ambig, reads_file,
          max_candidates, abismal_index, se_stats, out);
    }
    else {
      if (pbat_mode)
        run_paired_ended<a_rich, false>(VERBOSE, allow_ambig, reads_file,
            reads_file2, max_candidates, abismal_index, pe_stats,
            out);
      else if (random_pbat)
        run_paired_ended<t_rich, true>(VERBOSE,  allow_ambig, reads_file,
            reads_file2, max_candidates, abismal_index, pe_stats,
            out);
      else
        run_paired_ended<t_rich, false>(VERBOSE, allow_ambig, reads_file,
            reads_file2, max_candidates, abismal_index, pe_stats,
            out);
    }
    std::ofstream stats_of;
    if (!stats_outfile.empty())
      stats_of.open(stats_outfile.c_str(), std::ios::binary);
    std::ostream stats_out(
      stats_outfile.empty() ? std::cerr.rdbuf() : stats_of.rdbuf()
    );
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
