/* Copyright (C) 2018-2025 Andrew D. Smith and Guilherme Sena
 *
 * Authors: Andrew D. Smith and Guilherme Sena
 *
 * This file is part of ABISMAL.
 *
 * ABISMAL is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * ABISMAL is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 */

#include "abismal.hpp"

#include <config.h>
#include <unistd.h>

#include <atomic>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "AbismalAlign.hpp"
#include "AbismalIndex.hpp"
#include "OptionParser.hpp"
#include "bamxx.hpp"
#include "bisulfite_utils.hpp"
#include "dna_four_bit_bisulfite.hpp"
#include "popcnt.hpp"
#include "sam_record.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

using abismal_clock = std::chrono::steady_clock;
using abismal_timepoint = std::chrono::time_point<abismal_clock>;

using bamxx::bam_rec;

using AbismalAlignSimple =
  AbismalAlign<simple_aln::mismatch_score, simple_aln::indel>;

typedef std::uint16_t flags_t;  // every bit is a flag
typedef std::int16_t score_t;   // aln score, edit distance, hamming distance
typedef std::vector<uint8_t> Read;          // 4-bit encoding of reads
typedef std::vector<element_t> PackedRead;  // 4-bit encoding of reads

namespace abismal_concurrency {
static constexpr std::uint32_t max_n_threads = 1024;
static std::uint32_t n_threads = 1;
static std::mutex read_mutex;
static std::mutex write_mutex;
static std::mutex report_mutex;
static bool
invalid_n_threads() {
  return n_threads == 0 || n_threads > max_n_threads;
}
}  // namespace abismal_concurrency

enum conversion_type { t_rich = false, a_rich = true };

static void
log_msg(const std::string &s) {
  using std::chrono::system_clock;
  auto tmp = system_clock::to_time_t(system_clock::now());
  std::string time_fmt(std::ctime(&tmp));
  time_fmt.pop_back();
  std::cerr << "[" << time_fmt << "] " << s << std::endl;
}

static constexpr conversion_type
flip_conv(const conversion_type conv) {
  return conv == t_rich ? a_rich : t_rich;
}

static constexpr flags_t
get_strand_code(const char strand, const conversion_type conv) {
  return (((strand == '-') ? samflags::read_rc : 0) |
          ((conv == a_rich) ? bsflags::read_is_a_rich : 0));
}

struct ReadLoader {
  ReadLoader(const std::string &fn) : cur_line{0}, filename{fn}, in{fn, "r"} {}

  bool
  good() const {
    return in;
  }

  operator bool() const { return in; }

  std::size_t
  get_current_read() const {
    return cur_line / 4;
  }

  std::size_t
  get_current_byte() const {
    return in.tellg();
  }

  void
  load_reads(std::vector<std::string> &names, std::vector<std::string> &reads) {
    reads.clear();
    names.clear();

    std::size_t line_count = 0;
    const std::size_t num_lines_to_read = 4 * batch_size;
    std::string line;
    while (line_count < num_lines_to_read && getline(in, line)) {
      if (line_count % 4 == 0) {
        if (line.empty())
          throw std::runtime_error("file " + filename + " contains an empty " +
                                   "read name at line " +
                                   std::to_string(cur_line));
        names.emplace_back(line.substr(1, line.find_first_of(" \t") - 1));
      }
      else if (line_count % 4 == 1) {
        // read too long, may pass the end of the genome
        if (std::size(line) >= seed::padding_size)
          throw std::runtime_error(
            "found a read of size " + std::to_string(line.size()) +
            ", which is too long. Maximum allowed read size = " +
            std::to_string(seed::padding_size));

        if (count_if(std::cbegin(line), std::cend(line),
                     [](const char c) { return c != 'N'; }) < min_read_length)
          line.clear();
        else {
          while (line.back() == 'N')
            line.pop_back();                               // remove Ns from 3'
          line = line.substr(line.find_first_of("ACGT"));  // removes Ns from 5'
        }
        reads.emplace_back(line);
      }
      ++line_count;
      ++cur_line;
    }
  }

  std::uint32_t cur_line{};
  std::string filename;
  bamxx::bgzf_file in;

  static const std::size_t batch_size;
  static const std::uint32_t min_read_length;
};

const std::size_t ReadLoader::batch_size = 1000;

// GS: minimum length which an exact match can be
// guaranteed to map
const std::uint32_t ReadLoader::min_read_length =
  seed::key_weight + seed::window_size - 1;

// GS: used to allocate the appropriate dimensions of the banded
// alignment matrix for a batch of reads
static inline void
update_max_read_length(std::size_t &max_length,
                       const std::vector<std::string> &reads) {
  for (auto &i : reads)
    max_length = std::max(max_length, i.size());
}

struct se_element {   // size = 8
  score_t diffs;      // 2 bytes
  flags_t flags;      // 2 bytes
  std::uint32_t pos;  // 4 bytes

  se_element() : diffs(MAX_DIFFS), flags(0), pos(0) {}

  se_element(const score_t d, const flags_t f, const std::uint32_t p) :
    diffs(d), flags(f), pos(p) {}

  bool
  operator==(const se_element &rhs) const {
    return pos == rhs.pos && flags == rhs.flags;
  }

  bool
  operator!=(const se_element &rhs) const {
    return pos != rhs.pos || flags != rhs.flags;
  }

  // this is used to keep PE candidates sorted in the max heap
  bool
  operator<(const se_element &rhs) const {
    return diffs < rhs.diffs;
  }

  inline bool
  rc() const {
    return samflags::check(flags, samflags::read_rc);
  }

  inline bool
  elem_is_a_rich() const {
    return samflags::check(flags, bsflags::read_is_a_rich);
  }

  inline bool
  ambig() const {
    return samflags::check(flags, samflags::secondary_aln);
  }

  inline void
  set_ambig() {
    samflags::set(flags, samflags::secondary_aln);
  }

  inline bool
  empty() const {
    return pos == 0;
  }

  inline bool
  sure_ambig() const {
    return ambig() && diffs == 0;
  }

  inline void
  reset() {
    pos = 0;
    diffs = MAX_DIFFS;
  }

  inline void
  reset(const std::uint32_t readlen) {
    reset();
    diffs = static_cast<score_t>(invalid_hit_frac * readlen);
  }

  static double valid_frac;
  static const double invalid_hit_frac;
  static const score_t MAX_DIFFS;
};

double se_element::valid_frac = 0.1;
const score_t se_element::MAX_DIFFS = std::numeric_limits<score_t>::max() - 1;

// a liberal number of mismatches accepted to
// align a read downstream
const double se_element::invalid_hit_frac = 0.4;

static inline score_t
valid_diffs_cutoff(const std::uint32_t readlen, const double cutoff) {
  return static_cast<score_t>(cutoff * readlen);
}

static inline bool
valid_len(const std::uint32_t aln_len, const std::uint32_t readlen) {
  static const double min_aln_frac = 1.0 - se_element::invalid_hit_frac;

  return aln_len >=
         std::max(ReadLoader::min_read_length,
                  static_cast<std::uint32_t>(min_aln_frac * readlen));
}

static inline bool
valid(const se_element &s, const std::uint32_t aln_len,
      const std::uint32_t readlen, const double cutoff) {
  return valid_len(aln_len, readlen) &&
         s.diffs <= valid_diffs_cutoff(readlen, cutoff);
}

static inline bool
valid_hit(const se_element s, const std::uint32_t readlen) {
  return s.diffs < static_cast<score_t>(se_element::invalid_hit_frac * readlen);
}

template <class T>
static inline T
max16(const T x, const T y) {
  return (x > y) ? x : y;
}

struct se_candidates {
  se_candidates() :
    sz(1), best(se_element()), v(std::vector<se_element>(max_size)) {}

  inline bool
  full() const {
    return sz == max_size;
  };

  inline bool
  has_exact_match() const {
    return !best.empty();
  };

  void
  update_exact_match(const flags_t s, const std::uint32_t p) {
    const se_element cand(0, s, p);
    if (best.empty())
      best = cand;  // cand has ambig flag set to false

    else if (cand != best)
      best.set_ambig();
  }

  bool
  enough_good_hits() const {
    return full() && good_diff(cutoff);
  }

  bool
  good_diff(const score_t d) const {
    return (d <= good_cutoff);
  }

  bool
  should_do_sensitive() const {
    return (!full() || !good_diff(cutoff));
  }

  inline void
  set_specific() {
    cutoff = good_cutoff;
  }

  inline void
  set_sensitive() {
    cutoff = v.front().diffs;
  }

  void
  update_cand(const score_t d, const flags_t s, const std::uint32_t p) {
    if (full()) {
      std::pop_heap(std::begin(v), std::begin(v) + sz);
      v[sz - 1] = se_element(d, s, p);
    }
    else {
      v[sz++] = se_element(d, s, p);
    }
    std::push_heap(std::begin(v), std::begin(v) + sz);
  }

  void
  update(const bool specific, const score_t d, const flags_t s,
         const std::uint32_t p) {
    if (d == 0)
      update_exact_match(s, p);
    else
      update_cand(d, s, p);

    sure_ambig = best.sure_ambig();
    cutoff = (specific ? min16(cutoff, v.front().diffs) : v.front().diffs);
  }

  void
  reset() {
    best.reset();
    v.front().reset();

    cutoff = v.front().diffs;

    sure_ambig = false;
    sz = 1;
  }

  void
  reset(const std::uint32_t readlen) {
    best.reset(readlen);
    v.front().reset(readlen);
    cutoff = v.front().diffs;
    good_cutoff = readlen / 10u;

    sure_ambig = false;
    sz = 1;
  }

  // in SE reads, we sort to exclude duplicates
  void
  prepare_for_alignments() {
    sort(std::begin(v),
         std::begin(v) + sz,  // no sort_heap here as heapify used "diffs"
         [](const se_element &a, const se_element &b) {
           return (a.pos < b.pos) || (a.pos == b.pos && a.flags < b.flags);
         });
    sz = unique(std::begin(v), std::begin(v) + sz) - std::begin(v);
  }

  bool sure_ambig;
  score_t good_cutoff;
  score_t cutoff;
  std::uint32_t sz;
  se_element best;
  std::vector<se_element> v;

  static const std::uint32_t max_size;
};

const std::uint32_t se_candidates::max_size = 50u;

static inline bool
cigar_eats_ref(const std::uint32_t c) {
  return bam_cigar_type(bam_cigar_op(c)) & 2;
}

static inline std::uint32_t
cigar_rseq_ops(const bam_cigar_t &cig) {
  return accumulate(std::begin(cig), std::end(cig), 0u,
                    [](const std::uint32_t total, const std::uint32_t x) {
                      return total +
                             (cigar_eats_ref(x) ? bam_cigar_oplen(x) : 0);
                    });
}

static inline bool
chrom_and_posn(const ChromLookup &cl, const bam_cigar_t &cig,
               const std::uint32_t p, std::uint32_t &r_p, std::uint32_t &r_e,
               std::uint32_t &r_chr) {
  const std::uint32_t ref_ops = cigar_rseq_ops(cig);
  if (!cl.get_chrom_idx_and_offset(p, ref_ops, r_chr, r_p))
    return false;
  r_e = r_p + ref_ops;
  return true;
}

enum map_type { map_unmapped, map_unique, map_ambig };

static map_type
format_se(const bool allow_ambig, const se_element &res, const ChromLookup &cl,
          const std::string &read, const std::string &read_name,
          const bam_cigar_t &cigar, bam_rec &sr) {
  const bool ambig = res.ambig();
  const bool valid = !res.empty();
  if (!allow_ambig && ambig)
    return map_ambig;

  std::uint32_t ref_s = 0, ref_e = 0, chrom_idx = 0;
  if (!valid || !chrom_and_posn(cl, cigar, res.pos, ref_s, ref_e, chrom_idx))
    return map_unmapped;

  // ADS: we might be doing format_se for a mate in paried reads
  std::uint16_t flag = 0;
  if (res.rc())
    flag |= BAM_FREVERSE;

  if (allow_ambig && ambig)
    flag |= BAM_FSECONDARY;

  // flag |= BAM_FREAD1;  // ADS: this might be wrong...

  sr.b = bam_init1();
  int ret = bam_set1(sr.b,
                     read_name.size(),  // size_t l_qname,
                     read_name.data(),  // const char *qname,
                     flag,              // uint16_t flag,
                     chrom_idx - 1,     // int32_t tid (-1 for padding)
                     ref_s,             // hts_pos_t pos,
                     255,               // uint8_t mapq,
                     cigar.size(),      // size_t n_cigar,
                     cigar.data(),      // const uint32_t *cigar,
                     -1,                // int32_t mtid,
                     -1,                //  hts_pos_t mpos,
                     0,                 // hts_pos_t isize,
                     read.size(),       // size_t l_seq,
                     read.data(),       // const char *seq,
                     nullptr,           // const char *qual,
                     16);               // size_t l_aux);
  if (ret < 0)
    throw std::runtime_error("failed to format bam");

  ret = bam_aux_update_int(sr.b, "NM", res.diffs);
  if (ret < 0)
    throw std::runtime_error("bam_aux_update_int");

  ret = bam_aux_append(sr.b, "CV", 'A', 1,
                       (uint8_t *)(res.elem_is_a_rich() ? "A" : "T"));
  if (ret < 0)
    throw std::runtime_error("bam_aux_append");

  return ambig ? map_ambig : map_unique;
}

struct pe_element {
  pe_element() : aln_score(0), r1(se_element()), r2(se_element()) {}

  score_t
  diffs() const {
    return r1.diffs + r2.diffs;
  }

  void
  reset(const std::uint32_t readlen1, const std::uint32_t readlen2) {
    aln_score = 0;
    r1.reset(readlen1);
    r2.reset(readlen2);
    max_aln_score = simple_aln::best_pair_score(readlen1, readlen2);
  }

  void
  reset() {
    aln_score = 0;
    r1.reset();
    r2.reset();
  }

  bool
  update(const score_t scr, const se_element &s1, const se_element &s2) {
    const auto rd = r1.diffs + r2.diffs;
    const auto sd = s1.diffs + s2.diffs;
    if (scr > aln_score || (scr == aln_score && sd < rd)) {
      r1 = s1;
      r2 = s2;
      aln_score = scr;
      return true;
    }
    else if (scr == aln_score && sd == rd) {
      r1.set_ambig();
      return false;
    }

    return false;
  }

  // GS: this is used to decide whether ends should be
  // mapped as SE independently
  inline bool
  should_report(const bool allow_ambig) const {
    return !empty() && (allow_ambig || !ambig());
  }

  inline bool
  ambig() const {
    return r1.ambig();
  }

  inline bool
  empty() const {
    return r1.empty();
  }

  inline bool
  sure_ambig() const {
    return ambig() && (aln_score == max_aln_score);
  }

  score_t aln_score;
  score_t max_aln_score;
  se_element r1;
  se_element r2;

  static std::uint32_t min_dist;
  static std::uint32_t max_dist;
};

std::uint32_t pe_element::min_dist = 32;
std::uint32_t pe_element::max_dist = 3000;

static inline bool
valid_pair(const pe_element &best, const std::uint32_t readlen1,
           const std::uint32_t readlen2, const std::uint32_t aln_len1,
           const std::uint32_t aln_len2) {
  return valid_len(aln_len1, readlen1) && valid_len(aln_len2, readlen2) &&
         best.diffs() <=
           static_cast<score_t>(se_element::valid_frac * (aln_len1 + aln_len2));
}

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
format_pe(const bool allow_ambig, const pe_element &p, const ChromLookup &cl,
          const std::string &read1, const std::string &read2,
          const std::string &name1, const std::string &name2,
          const bam_cigar_t &cig1, const bam_cigar_t &cig2, bam_rec &sr1,
          bam_rec &sr2) {
  static const uint8_t cv[2] = {'T', 'A'};

  if (p.empty())
    return map_unmapped;

  const bool ambig = p.ambig();
  if (!allow_ambig && ambig)
    return map_ambig;

  std::uint32_t r_s1 = 0, r_e1 = 0, chr1 = 0;  // positions in chroms (0-based)
  std::uint32_t r_s2 = 0, r_e2 = 0, chr2 = 0;
  // PE chromosomes differ or couldn't be found, treat read as unmapped
  if (!chrom_and_posn(cl, cig1, p.r1.pos, r_s1, r_e1, chr1) ||
      !chrom_and_posn(cl, cig2, p.r2.pos, r_s2, r_e2, chr2) || chr1 != chr2)
    return map_unmapped;

  const bool rc = p.r1.rc();
  const int isize = rc ? (static_cast<int>(r_s1) - static_cast<int>(r_e2))
                       : (static_cast<int>(r_e2) - static_cast<int>(r_s1));

  std::uint16_t flag1 = 0;
  std::uint16_t flag2 = 0;

  flag1 |= BAM_FPAIRED | BAM_FPROPER_PAIR;

  flag2 |= BAM_FPAIRED | BAM_FPROPER_PAIR;

  if (p.r1.rc()) {  // ADS: is p.r1.rc() always !p.r2.rc()?
    flag1 |= BAM_FREVERSE;
    flag2 |= BAM_FMREVERSE;
  }
  if (p.r2.rc()) {
    flag2 |= BAM_FREVERSE;
    flag1 |= BAM_FMREVERSE;
  }
  if (allow_ambig && ambig) {
    // ADS: mark ambig for both the same way?
    flag1 |= BAM_FSECONDARY;
    flag2 |= BAM_FSECONDARY;
  }
  flag1 |= BAM_FREAD1;
  flag2 |= BAM_FREAD2;

  sr1.b = bam_init1();
  int ret = bam_set1(sr1.b,
                     name1.size(),  // size_t l_qname,
                     name1.data(),  // const char *qname,
                     flag1,         // uint16_t flag,
                     chr1 - 1,      // (-1 for padding) int32_t tid
                     r_s1,          // hts_pos_t pos,
                     255,           // uint8_t mapq,
                     cig1.size(),   // size_t n_cigar,
                     cig1.data(),   // const uint32_t *cigar,
                     chr2 - 1,      // (-1 for padding) int32_t mtid,
                     r_s2,          //  hts_pos_t mpos,
                     isize,         // hts_pos_t isize,
                     read1.size(),  // size_t l_seq,
                     read1.data(),  // const char *seq,
                     nullptr,       // const char *qual,
                     16);           // size_t l_aux);
  if (ret < 0)
    throw std::runtime_error("error formatting bam");

  ret = bam_aux_update_int(sr1.b, "NM", p.r1.diffs);
  if (ret < 0)
    throw std::runtime_error("error adding aux field");

  ret = bam_aux_append(sr1.b, "CV", 'A', 1, cv + p.r1.elem_is_a_rich());
  if (ret < 0)
    throw std::runtime_error("error adding aux field");

  sr2.b = bam_init1();
  ret = bam_set1(sr2.b,
                 name2.size(),  // size_t l_qname,
                 name2.data(),  // const char *qname,
                 flag2,         // uint16_t flag,
                 chr2 - 1,      // (-1 for padding) int32_t tid
                 r_s2,          // hts_pos_t pos,
                 255,           // uint8_t mapq,
                 cig2.size(),   // size_t n_cigar,
                 cig2.data(),   // const uint32_t *cigar,
                 chr1 - 1,      // (-1 for padding) int32_t mtid,
                 r_s1,          //  hts_pos_t mpos,
                 -isize,        // hts_pos_t isize,
                 read2.size(),  // size_t l_seq,
                 read2.data(),  // const char *seq,
                 nullptr,       // const char *qual,
                 16);           // size_t l_aux);
  if (ret < 0)
    throw std::runtime_error("failed to format bam");

  ret = bam_aux_update_int(sr2.b, "NM", p.r2.diffs);
  if (ret < 0)
    throw std::runtime_error("error adding aux field");

  ret = bam_aux_append(sr2.b, "CV", 'A', 1, cv + p.r2.elem_is_a_rich());
  if (ret < 0)
    throw std::runtime_error("error adding aux field");

  return ambig ? map_ambig : map_unique;
}

struct pe_candidates {
  pe_candidates() : v(std::vector<se_element>(max_size_large)) {}

  inline void
  reset(const std::uint32_t readlen) {
    v.front().reset(readlen);
    sure_ambig = false;
    cutoff = v.front().diffs;
    good_cutoff = static_cast<score_t>(readlen / 10);
    sz = 1;
    capacity = max_size_small;
  }

  inline void
  set_specific() {
    cutoff = good_cutoff;
  }

  inline void
  set_sensitive() {
    cutoff = v.front().diffs;
  }

  inline bool
  should_align() {
    return (sz != max_size_large || cutoff != 0);
  }

  inline bool
  full() const {
    return sz == capacity;
  }

  inline bool
  good_diff(const score_t d) const {
    return (d <= good_cutoff);
  }

  inline bool
  enough_good_hits() const {
    return full() && good_diff(cutoff);
  }

  inline bool
  should_do_sensitive() const {
    return (capacity == max_size_small || !good_diff(cutoff));
  }

  void
  update(const bool specific, const score_t d, const flags_t s,
         const std::uint32_t p) {
    if (full()) {
      // doubles capacity if heap is filled with good matches
      if (specific && capacity != max_size_large && good_diff(d)) {
        ++capacity;
      }
      else
        std::pop_heap(std::begin(v), std::begin(v) + sz--);
    }

    v[sz++] = se_element(d, s, p);
    std::push_heap(std::begin(v), std::begin(v) + sz);

    // update cutoff and sure_ambig
    cutoff = (specific ? min16(cutoff, v.front().diffs) : v.front().diffs);
    sure_ambig = (full() && cutoff == 0);
  }

  void
  prepare_for_mating() {
    sort(
      std::begin(v),
      std::begin(v) + sz,  // no sort_heap here as heapify used "diffs"
      [](const se_element &a, const se_element &b) { return a.pos < b.pos; });
    sz = unique(std::begin(v), std::begin(v) + sz) - std::begin(v);
  }

  bool sure_ambig;
  score_t cutoff;
  score_t good_cutoff;
  std::uint32_t sz;
  std::uint32_t capacity;
  std::vector<se_element> v;

  static const std::uint32_t max_size_small = 32u;
  static const std::uint32_t max_size_large = (max_size_small) << 10u;
};

struct se_map_stats {
  // total_reads is the total number of reads from the fastq file for
  // mapping that are considered for mapping as single-end reads. This
  // is all reads if the input is single-end, or all read ends for
  // which concordant mapping is not found in the case of paired-end.
  std::atomic_uint32_t total_reads{};

  // reads_mapped_unique is the number of reads for which exactly one location
  // in the reference genome has the best mapping score and that score meets the
  // minimum critiera for a good match.
  std::atomic_uint32_t reads_mapped_unique{};

  // reads_mapped_unique_fraction is the ratio of reads_mapped_unique
  // over total_reads.
  double reads_mapped_unique_fraction{};

  // reads_mapped_ambiguous is the number of reads that have two equally good
  // mapping locations within the reference genome, both meeting the
  // minimum criteria for a match.
  std::atomic_uint32_t reads_mapped_ambiguous{};

  // reads_mapped_ambiguous_fraction is the ratio of
  // reads_mapped_ambiguous over total_reads.
  double reads_mapped_ambiguous_fraction{};

  // reads_mapped is the sum of reads_mapped_unique and
  // reads_mapped_ambiguous.
  std::uint32_t reads_mapped{};

  // reads_mapped_fraction is the ratio of reads_mapped over
  // total_reads.
  double reads_mapped_fraction{};

  // reads_unmapped is the number of reads that have no mapping location
  // in the reference genome with a score that meets the minimum
  // criteria.
  std::atomic_uint32_t reads_unmapped{};

  // percent_unmapped is the ratio of reads_unmapped over total_reads,
  // multiplied by 100.
  double percent_unmapped{};

  // skipped_reads is the number of reads that are skipped and for which
  // mapping is not attempted due to poor quality -- usually these are
  // reads that are too short after adaptor trimming.
  std::atomic_uint32_t skipped_reads{};

  // percent_skipped is the ratio of skipped_reads over total_reads,
  // multiplied by 100.
  double percent_skipped{};

  // edit_distance is the sum of edit distances over all reads that
  // should be repored. This includes uniquely mapped reads but also
  // ambiguously mapping reads if the user requests those among the
  // output. This is used to obtain an average value.
  std::atomic_uint64_t edit_distance{};

  // total_bases is the total number of reference genome bases covered
  // by all mapped reads. This is used to obtain an average value.
  std::atomic_uint64_t total_bases{};

  // edit_distance_mean is the ratio of edit_distance over
  // total_bases.
  double edit_distance_mean{};

  void
  update(const bool allow_ambig, const std::string &read,
         const bam_cigar_t &cigar, const se_element s) {
    ++total_reads;
    const bool valid = !s.empty();
    const bool ambig = s.ambig();
    reads_mapped_unique += (valid && !ambig);
    reads_mapped_ambiguous += (valid && ambig);
    reads_unmapped += !valid;
    skipped_reads += read.empty();

    if (valid && (allow_ambig || !ambig))
      update_error_rate(s.diffs, cigar);
  }

  void
  update_error_rate(const score_t diffs, const bam_cigar_t &cigar) {
    edit_distance += diffs;
    total_bases += cigar_rseq_ops(cigar);
  }

  void
  assign_values() {
    constexpr auto pct = [](const double a, const double b) {
      return ((b == 0) ? 0.0 : 100.0 * a / b);
    };

    reads_mapped = reads_mapped_unique + reads_mapped_ambiguous;
    const std::uint32_t total_reads_tmp = total_reads;
    const std::uint32_t denom = total_reads_tmp == 0 ? 1 : total_reads_tmp;
    reads_mapped_fraction = pct(reads_mapped, denom);
    reads_mapped_unique_fraction = pct(reads_mapped_unique, denom);
    reads_mapped_ambiguous_fraction = pct(reads_mapped_ambiguous, denom);
    edit_distance_mean = edit_distance / static_cast<double>(total_bases);
    percent_unmapped = pct(reads_unmapped, total_reads);
    percent_skipped = pct(skipped_reads, total_reads);
  }

  std::string
  tostring(const std::size_t n_tabs = 0) {
    static constexpr auto tab = "    ";

    assign_values();

    std::string t;
    for (std::size_t i = 0; i < n_tabs; ++i)
      t += tab;
    std::ostringstream oss;
    // clang-format off
    oss << t << "total_reads: " << total_reads << '\n'
        << t << "mapped: " << '\n'
        << t + tab << "num_mapped: " << reads_mapped << '\n'
        << t + tab << "num_unique: " << reads_mapped_unique << '\n'
        << t + tab << "num_ambiguous: " << reads_mapped_ambiguous << '\n'
        << t + tab << "percent_mapped: " << reads_mapped_fraction << '\n'
        << t + tab << "percent_unique: " << reads_mapped_unique_fraction << '\n'
        << t + tab << "percent_ambiguous: " << reads_mapped_ambiguous_fraction << '\n'
        << t + tab << "unique_error:" << '\n'
        << t + tab + tab << "edits: " << edit_distance << '\n'
        << t + tab + tab << "total_bases: " << total_bases << '\n'
        << t + tab + tab << "error_rate: " << edit_distance_mean << '\n'
        << t << "num_unmapped: " << reads_unmapped << '\n'
        << t << "num_skipped: " << skipped_reads << '\n'
        << t << "percent_unmapped: " << percent_unmapped << '\n'
        << t << "percent_skipped: " << percent_skipped << '\n';
    // clang-format on
    return oss.str();
  }
};

struct pe_map_stats {
  // total_read_pairs is the total number of read pairs in the pair of
  // input fastq files.
  std::atomic_uint32_t total_read_pairs{};

  // reads_mapped_unique is the number of read pairs for which exactly
  // one pair of concordant locations in the reference genome has the
  // best mapping score for the pair and that score meets the minimum
  // critiera for a good match.
  std::atomic_uint32_t read_pairs_mapped_unique{};

  // read_pairs_mapped_unique_fraction is the ratio of
  // read_pairs_mapped_unique over total_read_pairs. This value
  // should be between 0 and 1, but for historical reasons is scaled
  // for display as a percentage.
  double read_pairs_mapped_unique_fraction{};

  // read_pairs_mapped_ambiguous is the number of reads that have two
  // equally good concordant mapping location pairs in the reference
  // genome, both meeting the minimum criteria for a match.
  std::atomic_uint32_t read_pairs_mapped_ambiguous{};

  // read_pairs_mapped_ambiguous_fraction is the ratio of
  // read_pairs_mapped_ambiguous over total_read_pairs. This value
  // should be between 0 and 1, but for historical reasons is scaled
  // for display as a percentage.
  double read_pairs_mapped_ambiguous_fraction{};

  // read_pairs_mapped is the sum of read_pairs_mapped_unique and
  // read_pairs_mapped_ambiguous.
  std::uint32_t read_pairs_mapped{};

  // read_pairs_mapped_fraction is the ratio of read_pairs_mapped over
  // total_read_pairs.
  double read_pairs_mapped_fraction{};

  // read_pairs_unmapped is the number of read pairs for which no
  // concordant pair of locations satisfies the mapping criteria for
  // both ends.
  std::atomic_uint32_t read_pairs_unmapped{};

  // percent_unmapped_pairs is the ratio of reads_pairs_unmapped over
  // total_read_pairs, multiplied by 100.
  double percent_unmapped_pairs{};

  // read_pairs_skipped is the number of read pairs that are skipped
  // as a pair and for which concordant mapping is not attempted due
  // to poor quality -- usually this is because of the quality for one
  // member of the pair, as in the case of end2 being generally very
  // poor quality.
  std::atomic_uint32_t read_pairs_skipped{};

  // percent_skipped_pairs is the ratio of read_pairs_skipped over
  // total_read_pairs, multiplied by 100.
  double percent_skipped_pairs{};

  // edit_distance_pairs is the total number of reference genome bases
  // covered by all concordantly mapped read pairs. This is used to
  // obtain an average value.
  std::atomic_uint64_t edit_distance_pairs{};

  // total_bases_pairs is the total number of reference genome bases
  // covered by all concordantly mapped read pairs. This is used to
  // obtain an average value.
  std::atomic_uint64_t total_bases_pairs{};

  // edit_distance_pairs_mean is the ratio of edit_distance_pairs over
  // total_bases_pairs.
  double edit_distance_pairs_mean{};

  se_map_stats end1_stats{};
  se_map_stats end2_stats{};

  void
  update(const bool allow_ambig, const std::string &reads1,
         const std::string &reads2, const bam_cigar_t &cig1,
         const bam_cigar_t &cig2, const pe_element &p, const se_element s1,
         const se_element s2) {
    const bool valid = !p.empty();
    const bool ambig = p.ambig();
    total_read_pairs++;
    read_pairs_mapped_ambiguous += (valid && ambig);
    read_pairs_mapped_unique += (valid && !ambig);
    read_pairs_unmapped += !valid;
    read_pairs_skipped += (reads1.empty() || reads2.empty());

    if (p.should_report(allow_ambig)) {
      update_error_rate(p.r1.diffs, p.r2.diffs, cig1, cig2);
    }
    else {
      end1_stats.update(false, reads1, cig1, s1);
      end2_stats.update(false, reads1, cig2, s2);
    }
  }

  void
  update_error_rate(const score_t d1, const score_t d2, const bam_cigar_t &cig1,
                    const bam_cigar_t &cig2) {
    edit_distance_pairs += d1 + d2;
    total_bases_pairs += cigar_rseq_ops(cig1) + cigar_rseq_ops(cig2);
  }

  void
  assign_values() {
    constexpr auto pct = [](const double a, const double b) {
      return ((b == 0) ? 0.0 : 100.0 * a / b);
    };

    read_pairs_mapped = read_pairs_mapped_unique + read_pairs_mapped_ambiguous;
    const std::uint32_t total_read_pairs_tmp = total_read_pairs;
    const std::uint32_t denom =
      total_read_pairs_tmp == 0 ? 1 : total_read_pairs_tmp;
    read_pairs_mapped_fraction = pct(read_pairs_mapped, denom);
    read_pairs_mapped_unique_fraction = pct(read_pairs_mapped_unique, denom);
    read_pairs_mapped_ambiguous_fraction =
      pct(read_pairs_mapped_ambiguous, denom);
    edit_distance_pairs_mean =
      edit_distance_pairs / static_cast<double>(total_bases_pairs);
    percent_unmapped_pairs = pct(read_pairs_unmapped, total_read_pairs_tmp);
    percent_skipped_pairs = pct(read_pairs_skipped, total_read_pairs_tmp);
  }

  std::string
  tostring(const bool allow_ambig) {
    static std::string t = "    ";

    assign_values();

    std::ostringstream oss;
    oss << "pairs:" << '\n'
        << t << "total_pairs: " << total_read_pairs << '\n'
        << t << "mapped:" << '\n'
        << t + t << "num_mapped: " << read_pairs_mapped << '\n'
        << t + t << "num_unique: " << read_pairs_mapped_unique << '\n'
        << t + t << "num_ambiguous: " << read_pairs_mapped_ambiguous << '\n'
        << t + t << "percent_mapped: " << read_pairs_mapped_fraction << '\n'
        << t + t << "percent_unique: " << read_pairs_mapped_unique_fraction
        << '\n'
        << t + t
        << "percent_ambiguous: " << read_pairs_mapped_ambiguous_fraction << '\n'
        << t + t << "unique_error:" << '\n'
        << t + t + t << "edits: " << edit_distance_pairs << '\n'
        << t + t + t << "total_bases: " << total_bases_pairs << '\n'
        << t + t + t << "error_rate: " << edit_distance_pairs_mean << '\n'
        << t << "num_unmapped: " << read_pairs_unmapped << '\n'
        << t << "num_skipped: " << read_pairs_skipped << '\n'
        << t << "percent_unmapped: " << percent_unmapped_pairs << '\n'
        << t << "percent_skipped: " << percent_skipped_pairs << '\n';

    if (!allow_ambig)
      oss << "mate1:" << '\n'
          << end1_stats.tostring(1) << "mate2:" << '\n'
          << end2_stats.tostring(1);
    return oss.str();
  }
};

static void
select_output(const bool allow_ambig, const ChromLookup &cl,
              const std::string &read1, const std::string &name1,
              const std::string &read2, const std::string &name2,
              const bam_cigar_t &cig1, const bam_cigar_t &cig2,
              pe_element &best, se_element &se1, se_element &se2, bam_rec &sr1,
              bam_rec &sr2) {
  const map_type pe_map_type = format_pe(allow_ambig, best, cl, read1, read2,
                                         name1, name2, cig1, cig2, sr1, sr2);

  if (!best.should_report(allow_ambig) || pe_map_type == map_unmapped) {
    if (pe_map_type == map_unmapped)
      best.reset();
    if (format_se(allow_ambig, se1, cl, read1, name1, cig1, sr1) ==
        map_unmapped)
      se1.reset();

    if (format_se(allow_ambig, se2, cl, read2, name2, cig2, sr2) ==
        map_unmapped)
      se2.reset();
  }
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
static score_t
full_compare(const score_t cutoff, const PackedRead::const_iterator read_end,
             const std::uint32_t offset, PackedRead::const_iterator read_itr,
             Genome::const_iterator genome_itr) {
  // max number of matches per element
  static const score_t max_matches = static_cast<score_t>(16);
  score_t d = 0;
  for (; d <= cutoff && read_itr != read_end;
       d += max_matches - popcnt64((*read_itr) & /*16 bases from the read*/
                                   /*16 bases from the padded genome*/
                                   ((*genome_itr >> offset) |
                                    ((*++genome_itr << (63 - offset)) << 1))),
       ++read_itr)
    ;
  return d;
}

template <const std::uint16_t strand_code, const bool specific,
          class result_type>
static inline void
check_hits(const std::uint32_t offset, const PackedRead::const_iterator read_st,
           const PackedRead::const_iterator read_end,
           const Genome::const_iterator genome_st,
           const std::vector<std::uint32_t>::const_iterator &end_idx,
           std::vector<std::uint32_t>::const_iterator start_idx,
           result_type &res) {
  for (; start_idx != end_idx && !res.sure_ambig; ++start_idx) {
    // GS: adds the next candidate to L1d cache while current is compared
#ifdef __SSE__
    _mm_prefetch(&(*(genome_st + ((*(start_idx + 10) - offset) >> 4))),
                 _MM_HINT_T0);
#endif
    const std::uint32_t the_pos = *start_idx - offset;
    /* GS: the_pos & 15u tells if the position is a multiple of 16, in
     * which case it is aligned with the genome. Otherwise we need to
     * use the unaligned comparison function that offsets genome
     * position by the_pos (mod 32). Multiplied by 4 because each base
     * uses 4 bits */
    const score_t diffs =
      full_compare(res.cutoff, read_end, ((the_pos & 15u) << 2), read_st,
                   genome_st + (the_pos >> 4));

    if (diffs <= res.cutoff)
      res.update(specific, diffs, strand_code, the_pos);
  }
}

struct compare_bases {
  compare_bases(const genome_iterator g_) : g(g_) {}

  bool
  operator()(const std::uint32_t mid, const two_letter_t chr) const {
    return get_bit(*(g + mid)) < chr;
  }

  const genome_iterator g;
};

template <const std::uint32_t start_length>
static std::uint32_t
find_candidates(const std::uint32_t max_candidates,
                const Read::const_iterator read_start, const genome_iterator gi,
                const std::uint32_t read_lim,
                std::vector<std::uint32_t>::const_iterator &low,
                std::vector<std::uint32_t>::const_iterator &high) {
  std::uint32_t p = start_length;
  auto prev_low = low;
  auto prev_high = high;
  for (; p != read_lim && ((high - low) > max_candidates); ++p) {
    // keep last interval with >0 candidates
    prev_low = low;
    prev_high = high;

    // pointer to first 1 in the range
    const std::vector<std::uint32_t>::const_iterator first_1 =
      lower_bound(low, high, 1, compare_bases(gi + p));

    const two_letter_t the_bit = get_bit(*(read_start + p));
    high = the_bit ? high : first_1;
    low = the_bit ? first_1 : low;
  }

  // some bit narrows it down to 0 candidates, roll back to when we
  // had some candidates to work with.
  if (low == high) {
    --p;
    low = prev_low;
    high = prev_high;
  }
  return p;
}

template <const three_conv_type the_conv>
static inline three_letter_t
get_three_letter_num_fast(const uint8_t nt) {
  return (the_conv == c_to_t) ? nt & 5 :  // C=T=0, A=1, G=4
           nt & 10;                       // A=G=0, C=2, T=8
}

template <const three_conv_type the_conv> struct compare_bases_three {
  compare_bases_three(const genome_iterator g_) : g(g_) {}

  bool
  operator()(const std::uint32_t mid, const three_letter_t chr) const {
    return get_three_letter_num_fast<the_conv>(*(g + mid)) < chr;
  }

  const genome_iterator g;
};

template <const std::uint32_t start_length, const three_conv_type the_conv>
static std::uint32_t
find_candidates_three(const std::uint32_t max_candidates,
                      const Read::const_iterator read_start,
                      const genome_iterator gi, const std::uint32_t max_size,
                      std::vector<std::uint32_t>::const_iterator &low,
                      std::vector<std::uint32_t>::const_iterator &high) {
  std::uint32_t p = start_length;
  auto prev_low = low;
  auto prev_high = high;
  for (; p != max_size && ((high - low) > max_candidates); ++p) {
    // keep last interval with >0 candidates
    prev_low = low;
    prev_high = high;

    // pointer to first 1 in the range
    const auto first_1 = lower_bound(low, high, (the_conv == c_to_t ? 1 : 2),
                                     compare_bases_three<the_conv>(gi + p));

    const auto first_2 = lower_bound(low, high, (the_conv == c_to_t ? 4 : 8),
                                     compare_bases_three<the_conv>(gi + p));

    const three_letter_t the_num =
      get_three_letter_num_fast<the_conv>(*(read_start + p));

    if (the_conv == c_to_t) {
      high = (the_num == 0) ? first_1 : ((the_num == 1) ? first_2 : high);
      low = (the_num == 0) ? low : ((the_num == 1) ? first_1 : first_2);
    }
    else {
      high = (the_num == 0) ? first_1 : ((the_num == 2) ? first_2 : high);
      low = (the_num == 0) ? low : ((the_num == 2) ? first_1 : first_2);
    }
  }

  // some bit narrows it down to 0 candidates, roll back to when we
  // had some candidates to work with.
  if (low == high) {
    --p;
    low = prev_low;
    high = prev_high;
  }
  return p;
}

static constexpr three_conv_type
get_conv_type(const std::uint16_t strand_code) {
  return ((samflags::check(strand_code, bsflags::read_is_a_rich) ^
           samflags::check(strand_code, samflags::read_rc))
            ? (g_to_a)
            : (c_to_t));
}

template <const std::uint16_t strand_code, class result_type>
static void
process_seeds(const std::uint32_t max_candidates,
              const std::vector<std::uint32_t>::const_iterator counter_st,
              const std::vector<std::uint32_t>::const_iterator counter_three_st,
              const std::vector<std::uint32_t>::const_iterator index_st,
              const std::vector<std::uint32_t>::const_iterator index_three_st,
              const genome_iterator genome_st, const Read &read_seed,
              const PackedRead &packed_read, result_type &res) {
  static constexpr three_conv_type the_conv = get_conv_type(strand_code);

  const std::uint32_t readlen = read_seed.size();
  const PackedRead::const_iterator pack_s_idx(std::begin(packed_read));
  const PackedRead::const_iterator pack_e_idx(std::end(packed_read));

  std::uint32_t k = 0u;
  std::uint32_t k_three = 0u;
  std::uint32_t i = 0u;

  Read::const_iterator read_idx(std::begin(read_seed));
  std::vector<std::uint32_t>::const_iterator s_idx;
  std::vector<std::uint32_t>::const_iterator e_idx;
  std::vector<std::uint32_t>::const_iterator s_idx_three;
  std::vector<std::uint32_t>::const_iterator e_idx_three;

  std::uint32_t d_two = 0;
  std::uint32_t d_three = 0;
  std::uint32_t l_two = 0;
  std::uint32_t l_three = 0;

  get_1bit_hash(read_idx, k);
  get_base_3_hash<the_conv>(read_idx, k_three);

  const std::uint32_t specific_len =
    min16(readlen - seed::window_size, readlen >> 1u);
  const std::uint32_t specific_lim =
    max16(seed::window_size, static_cast<std::uint32_t>(readlen >> 1u));

  res.set_specific();
  for (i = 0; i < specific_lim && !res.sure_ambig; ++i, ++read_idx) {
    s_idx = index_st + *(counter_st + k);
    e_idx = index_st + *(counter_st + k + 1);
    l_two = find_candidates<seed::key_weight>(
      max_candidates, read_idx, genome_st, readlen - i, s_idx, e_idx);
    d_two = (e_idx - s_idx);

    s_idx_three = index_three_st + *(counter_three_st + k_three);
    e_idx_three = index_three_st + *(counter_three_st + k_three + 1);
    l_three = find_candidates_three<seed::key_weight_three, the_conv>(
      max_candidates, read_idx, genome_st, readlen - i, s_idx_three,
      e_idx_three);

    d_three = (e_idx_three - s_idx_three);

    // two-letter seeds
    if (d_two <= max_candidates || l_two >= specific_len)
      check_hits<strand_code, true>(i, pack_s_idx, pack_e_idx, genome_st.itr,
                                    e_idx, s_idx, res);

    // three-letter seeds
    if (d_three <= max_candidates || l_three >= specific_len)
      check_hits<strand_code, true>(i, pack_s_idx, pack_e_idx, genome_st.itr,
                                    e_idx_three, s_idx_three, res);

    shift_hash_key(*(read_idx + seed::key_weight), k);
    shift_three_key<the_conv>(*(read_idx + seed::key_weight_three), k_three);
  }

  if (!res.should_do_sensitive())
    return;

  read_idx = std::begin(read_seed);
  get_1bit_hash(read_idx, k);
  get_base_3_hash<the_conv>(read_idx, k_three);

  res.set_sensitive();

  const std::uint32_t lim_two = readlen - seed::key_weight + 1;

  // GS: this is to avoid chasing down uninformative two-letter
  // seeds when there is a sufficiently high number of three
  // letter seeds that is lower than the number of two-letter hits
  static const std::uint32_t MIN_FOLD_SIZE = 10;
  for (i = 0; i < lim_two && !res.sure_ambig; ++i, ++read_idx) {
    s_idx = index_st + *(counter_st + k);
    e_idx = index_st + *(counter_st + k + 1);
    d_two = (e_idx - s_idx);

    s_idx_three = index_three_st + *(counter_three_st + k_three);
    e_idx_three = index_three_st + *(counter_three_st + k_three + 1);
    d_three = (e_idx_three - s_idx_three);

    // two-letter seeds
    if (d_two != 0 && d_two <= max_candidates &&
        (d_three == 0 || d_two <= MIN_FOLD_SIZE * d_three))
      check_hits<strand_code, true>(i, pack_s_idx, pack_e_idx, genome_st.itr,
                                    e_idx, s_idx, res);

    // three-letter seeds
    if (d_three != 0 && d_three <= max_candidates)
      check_hits<strand_code, true>(i, pack_s_idx, pack_e_idx, genome_st.itr,
                                    e_idx_three, s_idx_three, res);

    shift_hash_key(*(read_idx + seed::key_weight), k);
    shift_three_key<the_conv>(*(read_idx + seed::key_weight_three), k_three);
  }
}

template <const bool convert_a_to_g>
static void
prep_read(const std::string &r, Read &pread) {
  pread.resize(r.size());
  for (std::size_t i = 0; i != r.size(); ++i)
    pread[i] =
      (convert_a_to_g ? (encode_base_a_rich[static_cast<unsigned char>(r[i])])
                      : (encode_base_t_rich[static_cast<unsigned char>(r[i])]));
}

/* GS: this function simply converts the vector<uint8_t> pread
 * to a vector<uint64_t> by putting 16 bases in each element of
 * the packed read. If the read length does not divide 16, we add
 * 1111s to the remaining positions so it divides 16. The remaining
 * bases match all bases in the reference genome
 * */
static void
pack_read(const Read &pread, PackedRead &packed_pread) {
  static const element_t base_match_any = static_cast<element_t>(0xF);
  static const std::size_t NUM_BASES_PER_ELEMENT = 16;
  const std::size_t sz = pread.size();
  const std::size_t num_complete_pos = sz / NUM_BASES_PER_ELEMENT;

  // divide by 16 and add an extra position if remainder not 0
  packed_pread.resize((sz + NUM_BASES_PER_ELEMENT - 1) / NUM_BASES_PER_ELEMENT);
  PackedRead::iterator it(std::begin(packed_pread));

  // first add the complete positions (i.e. having all 16 bases)
  std::size_t pread_ind = 0;
  for (std::size_t i = 0; i < num_complete_pos; ++i) {
    *it = 0;
    for (std::size_t j = 0; j < NUM_BASES_PER_ELEMENT; ++j)
      *it |= (static_cast<element_t>(pread[pread_ind++]) << (j << 2));
    ++it;
  }

  // do not fill the flanking position
  if (pread_ind == sz)
    return;

  // now put only the remaining bases in the last pos. The rest
  // should match any base in the reference
  *it = 0;
  std::size_t j = 0;
  while (pread_ind < sz)
    *it |= (static_cast<element_t>(pread[pread_ind++]) << ((j++) << 2));

  while (j < NUM_BASES_PER_ELEMENT)
    *it |= base_match_any << ((j++) << 2);
}

static inline bool
same_pos(const std::uint32_t pos1, const std::uint32_t pos2) {
  const std::uint32_t diff = (pos1 > pos2) ? (pos1 - pos2) : (pos2 - pos1);
  static const std::uint32_t MIN_DIFF_FOR_EQUAL = 3;
  return diff <= MIN_DIFF_FOR_EQUAL;
}

static void
align_se_candidates(const Read &pread_t, const Read &pread_t_rc,
                    const Read &pread_a, const Read &pread_a_rc,
                    const double cutoff, se_candidates &res, se_element &best,
                    bam_cigar_t &cigar, AbismalAlignSimple &aln) {
  const score_t readlen = static_cast<score_t>(pread_t.size());
  const score_t max_diffs = valid_diffs_cutoff(readlen, cutoff);
  const score_t max_scr = simple_aln::best_single_score(readlen);
  if (res.has_exact_match()) {  // exact match, no need to align
    best = res.best;            // ambig info also passed here
    make_default_cigar(readlen, cigar);
    return;
  }

  score_t best_scr = 0;
  std::uint32_t cand_pos = 0;
  std::uint32_t best_pos = 0;

  res.prepare_for_alignments();
  std::vector<se_element>::const_iterator it(std::begin(res.v));
  const std::vector<se_element>::const_iterator lim(it + res.sz);

  for (; it != lim && it->empty(); ++it)
    ;
  for (; it != lim; ++it) {
    if (valid_hit(*it, readlen)) {
      cand_pos = it->pos;
      const score_t cand_scr = aln.align<false>(
        it->diffs, max_diffs,
        ((it->rc()) ? ((it->elem_is_a_rich()) ? (pread_t_rc) : (pread_a_rc))
                    : ((it->elem_is_a_rich()) ? (pread_a) : (pread_t))),
        cand_pos);

      if (cand_scr > best_scr) {
        best = *it;  // ambig = false
        best_scr = cand_scr;
        best_pos = cand_pos;
      }
      else if (cand_scr == best_scr &&
               ((cand_scr == max_scr) ? (cand_pos != best_pos)
                                      : !same_pos(cand_pos, best_pos)))
        best.set_ambig();
    }
  }

  if (best.pos != 0) {
    // recovers traceback to build CIGAR
    aln.align<true>(best.diffs, max_diffs,
                    (best.rc())
                      ? ((best.elem_is_a_rich()) ? (pread_t_rc) : (pread_a_rc))
                      : ((best.elem_is_a_rich()) ? (pread_a) : (pread_t)),
                    best.pos);

    std::uint32_t len = 0;
    aln.build_cigar_len_and_pos(best.diffs, max_diffs, cigar, len, best.pos);
    best.diffs = simple_aln::edit_distance(best_scr, len, cigar);

    // do not report and count it as unmapped if not valid
    if (!valid(best, len, readlen, cutoff))
      best.reset();
  }
  else
    best.reset();
}

static inline bool
valid_bam_rec(const bam_rec &b) {
  return b.b;
}

static inline void
reset_bam_rec(bam_rec &b) {
  if (b.b)
    bam_destroy1(b.b);
  b.b = nullptr;
}

template <const conversion_type conv>
static void
map_single_ended(const bool show_progress, const bool allow_ambig,
                 const AbismalIndex &abismal_index, ReadLoader &rl,
                 se_map_stats &se_stats, bamxx::bam_header &hdr,
                 bamxx::bam_out &out, ProgressBar &progress) {
  const auto counter_st(std::begin(abismal_index.counter));
  const auto counter_t_st(std::begin(abismal_index.counter_t));
  const auto counter_a_st(std::begin(abismal_index.counter_a));

  const auto index_st(std::begin(abismal_index.index));
  const auto index_t_st(std::begin(abismal_index.index_t));
  const auto index_a_st(std::begin(abismal_index.index_a));

  const genome_iterator genome_st(std::begin(abismal_index.genome));
  const std::uint32_t max_candidates = abismal_index.max_candidates;

  // batch variables used in reporting the SAM entry
  std::vector<std::string> names;
  std::vector<std::string> reads;
  std::vector<bam_cigar_t> cigar;
  std::vector<se_element> bests;
  std::vector<bam_rec> mr;

  names.reserve(ReadLoader::batch_size);
  reads.reserve(ReadLoader::batch_size);

  cigar.resize(ReadLoader::batch_size);
  bests.resize(ReadLoader::batch_size);
  mr.resize(ReadLoader::batch_size);

  // pre-allocated variabes used idependently in each read
  Read pread, pread_rc;
  PackedRead packed_pread;
  se_candidates res;
  AbismalAlignSimple aln(genome_st);

  std::size_t the_byte = 0;

  while (rl) {
    {
      const std::lock_guard<std::mutex> lock(abismal_concurrency::read_mutex);
      rl.load_reads(names, reads);
      the_byte = rl.get_current_byte();
    }

    std::size_t max_batch_read_length = 0;
    update_max_read_length(max_batch_read_length, reads);

    aln.reset(max_batch_read_length);

    const std::size_t n_reads = reads.size();

    for (std::size_t i = 0; i < n_reads; ++i) {
      res.reset(reads[i].size());
      bests[i].reset();
      if (!reads[i].empty()) {
        prep_read<conv>(reads[i], pread);
        pack_read(pread, packed_pread);
        process_seeds<get_strand_code('+', conv)>(
          max_candidates, counter_st,
          ((conv == t_rich) ? (counter_t_st) : (counter_a_st)),

          index_st, ((conv == t_rich) ? (index_t_st) : (index_a_st)), genome_st,
          pread, packed_pread, res);

        const std::string read_rc(revcomp(reads[i]));
        prep_read<!conv>(read_rc, pread_rc);
        pack_read(pread_rc, packed_pread);

        process_seeds<get_strand_code('-', conv)>(
          max_candidates, counter_st,
          (conv == t_rich) ? counter_a_st : counter_t_st, index_st,
          (conv == t_rich) ? index_a_st : index_t_st, genome_st, pread_rc,
          packed_pread, res);

        align_se_candidates(pread, pread_rc, pread, pread_rc,
                            se_element::valid_frac, res, bests[i], cigar[i],
                            aln);
        if (format_se(allow_ambig, bests[i], abismal_index.cl, reads[i],
                      names[i], cigar[i], mr[i]) == map_unmapped)
          bests[i].reset();
      }
    }
    {
      const std::lock_guard<std::mutex> lock(abismal_concurrency::write_mutex);
      for (std::size_t i = 0; i < n_reads; ++i) {
        if (valid_bam_rec(mr[i]) && !out.write(hdr, mr[i]))
          throw std::runtime_error("failed to write bam");
      }
    }
    for (std::size_t i = 0; i < n_reads; ++i) {
      if (valid_bam_rec(mr[i]))
        reset_bam_rec(mr[i]);
      se_stats.update(allow_ambig, reads[i], cigar[i], bests[i]);
      cigar[i].clear();
    }
    if (show_progress) {
      const std::lock_guard<std::mutex> lock(abismal_concurrency::report_mutex);
      if (progress.time_to_report(the_byte))
        progress.report(std::cerr, the_byte);
    }
  }
}

static void
map_single_ended_rand(const bool show_progress, const bool allow_ambig,
                      const AbismalIndex &abismal_index, ReadLoader &rl,
                      se_map_stats &se_stats, bamxx::bam_header &hdr,
                      bamxx::bam_out &out, ProgressBar &progress) {
  const auto counter_st(std::cbegin(abismal_index.counter));
  const auto counter_t_st(std::cbegin(abismal_index.counter_t));
  const auto counter_a_st(std::cbegin(abismal_index.counter_a));

  const auto index_st(std::cbegin(abismal_index.index));
  const auto index_t_st(std::cbegin(abismal_index.index_t));
  const auto index_a_st(std::cbegin(abismal_index.index_a));

  const std::uint32_t max_candidates = abismal_index.max_candidates;

  const genome_iterator genome_st(std::cbegin(abismal_index.genome));

  std::vector<std::string> names;
  std::vector<std::string> reads;
  std::vector<bam_cigar_t> cigar;
  std::vector<se_element> bests;
  std::vector<bam_rec> mr;

  names.reserve(ReadLoader::batch_size);
  reads.reserve(ReadLoader::batch_size);
  cigar.resize(ReadLoader::batch_size);
  bests.resize(ReadLoader::batch_size);
  mr.resize(ReadLoader::batch_size);

  // GS: pre-allocated variables used once per read
  // and not used for reporting
  Read pread_t, pread_t_rc, pread_a, pread_a_rc;
  PackedRead packed_pread;
  se_candidates res;
  AbismalAlignSimple aln(genome_st);

  std::size_t the_byte = 0;

  while (rl) {
    {
      const std::lock_guard<std::mutex> lock(abismal_concurrency::read_mutex);
      rl.load_reads(names, reads);
      the_byte = rl.get_current_byte();
    }

    std::size_t max_batch_read_length = 0;
    update_max_read_length(max_batch_read_length, reads);

    aln.reset(max_batch_read_length);

    const std::size_t n_reads = reads.size();
    for (std::size_t i = 0; i < n_reads; ++i) {
      res.reset(reads[i].size());
      bests[i].reset();
      if (!reads[i].empty()) {
        // T-rich, + strand
        prep_read<t_rich>(reads[i], pread_t);
        pack_read(pread_t, packed_pread);
        process_seeds<get_strand_code('+', t_rich)>(
          max_candidates, counter_st, counter_t_st, index_st, index_t_st,
          genome_st, pread_t, packed_pread, res);

        // A-rich, + strand
        prep_read<a_rich>(reads[i], pread_a);
        pack_read(pread_a, packed_pread);
        process_seeds<get_strand_code('+', a_rich)>(
          max_candidates, counter_st, counter_a_st, index_st, index_a_st,
          genome_st, pread_a, packed_pread, res);

        // A-rich, - strand
        const std::string read_rc(revcomp(reads[i]));
        prep_read<t_rich>(read_rc, pread_t_rc);
        pack_read(pread_t_rc, packed_pread);
        process_seeds<get_strand_code('-', a_rich)>(
          max_candidates, counter_st, counter_t_st, index_st, index_t_st,
          genome_st, pread_t_rc, packed_pread, res);

        // T-rich, - strand
        prep_read<a_rich>(read_rc, pread_a_rc);
        pack_read(pread_a_rc, packed_pread);
        process_seeds<get_strand_code('-', t_rich)>(
          max_candidates, counter_st, counter_a_st, index_st, index_a_st,
          genome_st, pread_a_rc, packed_pread, res);

        align_se_candidates(pread_t, pread_t_rc, pread_a, pread_a_rc,
                            se_element::valid_frac, res, bests[i], cigar[i],
                            aln);
        if (format_se(allow_ambig, bests[i], abismal_index.cl, reads[i],
                      names[i], cigar[i], mr[i]) == map_unmapped)
          bests[i].reset();
      }
    }
    {
      const std::lock_guard<std::mutex> lock(abismal_concurrency::write_mutex);
      for (std::size_t i = 0; i < n_reads; ++i)
        if (valid_bam_rec(mr[i]) && !out.write(hdr, mr[i]))
          throw std::runtime_error("failed to write bam");
    }
    for (std::size_t i = 0; i < n_reads; ++i) {
      if (valid_bam_rec(mr[i]))
        reset_bam_rec(mr[i]);
      se_stats.update(allow_ambig, reads[i], cigar[i], bests[i]);
      cigar[i].clear();
    }
    if (show_progress) {
      const std::lock_guard<std::mutex> lock(abismal_concurrency::report_mutex);
      if (progress.time_to_report(the_byte))
        progress.report(std::cerr, the_byte);
    }
  }
}

static std::string
format_duration(const abismal_timepoint start_time,
                const abismal_timepoint stop_time) {
  const std::chrono::duration<double> elapsed_seconds{stop_time - start_time};
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(2) << elapsed_seconds.count() << "s";
  return oss.str();
}

template <const conversion_type conv, const bool random_pbat>
static void
run_single_ended(const bool show_progress, const bool allow_ambig,
                 const std::string &reads_file,
                 const AbismalIndex &abismal_index, se_map_stats &se_stats,
                 bamxx::bam_header &hdr, bamxx::bam_out &out) {
  ReadLoader rl(reads_file);
  ProgressBar progress(get_filesize(reads_file), "mapping reads");

  const auto start_time{abismal_clock::now()};

  std::vector<std::thread> threads;
  for (auto i = 0u; i < abismal_concurrency::n_threads; ++i)
    threads.emplace_back([&] {
      if (random_pbat)
        map_single_ended_rand(show_progress, allow_ambig, abismal_index, rl,
                              se_stats, hdr, out, progress);
      else
        map_single_ended<conv>(show_progress, allow_ambig, abismal_index, rl,
                               se_stats, hdr, out, progress);
    });
  for (auto &thread : threads)
    thread.join();

  if (show_progress) {
    std::cerr << "\n";
    const auto stop_time{abismal_clock::now()};
    log_msg("reads mapped: " + std::to_string(rl.get_current_read()));
    log_msg("total mapping time: " + format_duration(start_time, stop_time));
  }
}

static void
best_single(const pe_candidates &pres, se_candidates &res) {
  const auto lim(std::begin(pres.v) + pres.sz);
  for (auto i(std::begin(pres.v)); i != lim && !res.sure_ambig; ++i)
    res.update(false, i->diffs, i->flags, i->pos);
}

template <const bool swap_ends>
static void
best_pair(const pe_candidates &res1, const pe_candidates &res2,
          const Read &pread1, const Read &pread2, bam_cigar_t &cigar1,
          bam_cigar_t &cigar2, std::vector<score_t> &mem_scr1,
          AbismalAlignSimple &aln, pe_element &best) {
  std::vector<se_element>::const_iterator j1(std::begin(res1.v));
  std::vector<se_element>::const_iterator j2(std::begin(res2.v));

  const std::vector<se_element>::const_iterator j1_end = j1 + res1.sz;
  const std::vector<se_element>::const_iterator j2_end = j2 + res2.sz;
  const std::vector<se_element>::const_iterator j1_beg(j1);

  // remembers alignment info on end1 to avoid redoing work
  const auto a1_beg(std::begin(mem_scr1));
  const auto a1_end(a1_beg + res1.sz);
  auto a1 = a1_beg;
  std::fill(a1_beg, a1_end, 0);

  const std::uint32_t readlen1 = pread1.size();
  const std::uint32_t readlen2 = pread2.size();
  const score_t max_diffs1 =
    valid_diffs_cutoff(readlen1, se_element::valid_frac);
  const score_t max_diffs2 =
    valid_diffs_cutoff(readlen2, se_element::valid_frac);

  score_t scr1 = 0;
  score_t scr2 = 0;
  score_t best_scr1 = 0;
  score_t best_scr2 = 0;
  std::uint32_t best_pos1 = 0;
  std::uint32_t best_pos2 = 0;

  se_element s1;
  se_element s2;

  // GS: skips empty hits which are in the beginning
  // because empty hits, by definition, have pos = 0
  for (; j1 != j1_end && j1->empty(); ++j1, ++a1)
    ;
  for (; j2 != j2_end && j2->empty(); ++j2)
    ;

  for (; j2 != j2_end && !best.sure_ambig(); ++j2) {
    s2 = *j2;
    scr2 = 0;

    // rewind to first concordant position. Needed in case of
    // many-to-many concordance between candidates
    const std::uint32_t lim = s2.pos + readlen2;
    for (; (j1 == j1_end) ||
           (j1 != j1_beg && j1->pos + pe_element::max_dist >= lim);
         --j1, --a1)
      ;

    for (; j1 != j1_end && j1->pos + pe_element::max_dist < lim; ++j1, ++a1)
      ;
    for (; j1 != j1_end && j1->pos + pe_element::min_dist <= lim &&
           !best.sure_ambig();
         ++j1, ++a1) {
      s1 = *j1;

      if (scr2 == 0) {  // ensures elements in j2 are aligned only once
        scr2 = aln.align<false>(j2->diffs, max_diffs2, pread2, s2.pos);
      }

      if (*a1 == 0) {  // ensures elements in j1 are aligned only once
        scr1 = aln.align<false>(j1->diffs, max_diffs1, pread1, s1.pos);
        *a1 = scr1;
      }

      const score_t pair_scr = scr2 + *a1;
      if (swap_ends ? best.update(pair_scr, s2, s1)
                    : best.update(pair_scr, s1, s2)) {
        best_scr1 = scr1;
        best_scr2 = scr2;
        best_pos1 = j1->pos;
        best_pos2 = j2->pos;
      }

      // if (best.sure_ambig()) return;
    }
  }

  if (best_pos1 != 0) {  // a new better alignment was found

    s1 = (swap_ends) ? (best.r2) : (best.r1);
    s2 = (swap_ends) ? (best.r1) : (best.r2);

    // re-aligns pos 1 with traceback
    std::uint32_t len1 = 0;
    aln.align<true>(s1.diffs, max_diffs1, pread1, best_pos1);
    aln.build_cigar_len_and_pos(s1.diffs, max_diffs1, cigar1, len1, best_pos1);
    s1.pos = best_pos1;
    s1.diffs = simple_aln::edit_distance(best_scr1, len1, cigar1);

    // re-aligns pos 2 with traceback
    std::uint32_t len2 = 0;
    aln.align<true>(s2.diffs, max_diffs2, pread2, best_pos2);
    aln.build_cigar_len_and_pos(s2.diffs, max_diffs2, cigar2, len2, best_pos2);
    s2.pos = best_pos2;
    s2.diffs = simple_aln::edit_distance(best_scr2, len2, cigar2);

    // last check if, after alignment, mates are still concordant
    const std::uint32_t frag_end = best_pos2 + len2;
    if (frag_end >= best_pos1 + pe_element::min_dist &&
        frag_end <= best_pos1 + pe_element::max_dist) {
      best.r1 = (swap_ends) ? (s2) : (s1);
      best.r2 = (swap_ends) ? (s1) : (s2);
    }
    else
      best.reset();
  }
}

template <const bool swap_ends>
static bool
select_maps(const Read &pread1, const Read &pread2, bam_cigar_t &cig1,
            bam_cigar_t &cig2, pe_candidates &res1, pe_candidates &res2,
            std::vector<score_t> &mem_scr1, se_candidates &res_se1,
            se_candidates &res_se2, AbismalAlignSimple &aln, pe_element &best) {
  if (res1.should_align() && res2.should_align()) {
    res1.prepare_for_mating();
    res2.prepare_for_mating();
    best_pair<swap_ends>(res1, res2, pread1, pread2, cig1, cig2, mem_scr1, aln,
                         best);
  }
  best_single(res1, res_se1);
  best_single(res2, res_se2);
  return true;
}

template <const bool cmp, const bool swap_ends,
          const std::uint16_t strand_code1, const std::uint16_t strand_code2>
static inline bool
map_fragments(const std::uint32_t max_candidates, const std::string &read1,
              const std::string &read2,
              const std::vector<std::uint32_t>::const_iterator counter_st,
              const std::vector<std::uint32_t>::const_iterator counter_three_st,
              const std::vector<std::uint32_t>::const_iterator index_st,
              const std::vector<std::uint32_t>::const_iterator index_three_st,
              const genome_iterator genome_st, Read &pread1, Read &pread2,
              PackedRead &packed_pread, bam_cigar_t &cigar1,
              bam_cigar_t &cigar2, AbismalAlignSimple &aln, pe_candidates &res1,
              pe_candidates &res2, std::vector<score_t> &mem_scr1,
              se_candidates &res_se1, se_candidates &res_se2,
              pe_element &best) {
  res1.reset(read1.size());
  res2.reset(read2.size());

  if (read1.empty() && read2.empty())
    return false;

  if (!read1.empty()) {
    prep_read<cmp>(read1, pread1);
    pack_read(pread1, packed_pread);
    process_seeds<strand_code1>(max_candidates, counter_st, counter_three_st,
                                index_st, index_three_st, genome_st, pread1,
                                packed_pread, res1);
  }

  if (!read2.empty()) {
    const std::string read_rc(revcomp(read2));
    prep_read<cmp>(read_rc, pread2);
    pack_read(pread2, packed_pread);
    process_seeds<strand_code2>(max_candidates, counter_st, counter_three_st,
                                index_st, index_three_st, genome_st, pread2,
                                packed_pread, res2);
  }

  return select_maps<swap_ends>(pread1, pread2, cigar1, cigar2, res1, res2,
                                mem_scr1, res_se1, res_se2, aln, best);
}

template <const conversion_type conv>
static void
map_paired_ended(const bool show_progress, const bool allow_ambig,
                 const AbismalIndex &abismal_index, ReadLoader &rl1,
                 ReadLoader &rl2, pe_map_stats &pe_stats,
                 bamxx::bam_header &hdr, bamxx::bam_out &out,
                 ProgressBar &progress) {
  const auto counter_st(std::begin(abismal_index.counter));
  const auto counter_t_st(std::begin(abismal_index.counter_t));
  const auto counter_a_st(std::begin(abismal_index.counter_a));

  const auto index_st(std::begin(abismal_index.index));
  const auto index_t_st(std::begin(abismal_index.index_t));
  const auto index_a_st(std::begin(abismal_index.index_a));

  const std::uint32_t max_candidates = abismal_index.max_candidates;

  const genome_iterator genome_st(std::begin(abismal_index.genome));

  // GS: objects used to report reads, need as many copies as
  // the batch size
  std::vector<std::string> names1, reads1;
  std::vector<std::string> names2, reads2;

  std::vector<bam_cigar_t> cigar1;
  std::vector<bam_cigar_t> cigar2;

  std::vector<pe_element> bests;
  std::vector<se_element> bests_se1;
  std::vector<se_element> bests_se2;
  std::vector<bam_rec> mr1;
  std::vector<bam_rec> mr2;

  names1.reserve(ReadLoader::batch_size);
  reads1.reserve(ReadLoader::batch_size);
  cigar1.resize(ReadLoader::batch_size);

  names2.reserve(ReadLoader::batch_size);
  reads2.reserve(ReadLoader::batch_size);
  cigar2.resize(ReadLoader::batch_size);

  bests.resize(ReadLoader::batch_size);
  bests_se1.resize(ReadLoader::batch_size);
  bests_se2.resize(ReadLoader::batch_size);

  mr1.resize(ReadLoader::batch_size);
  mr2.resize(ReadLoader::batch_size);

  // GS: pre-allocated variables used once per read
  // and not used for reporting
  Read pread1, pread1_rc, pread2, pread2_rc;
  PackedRead packed_pread;

  AbismalAlignSimple aln(genome_st);

  pe_candidates res1;
  pe_candidates res2;
  std::vector<score_t> mem_scr1(res1.v.size());
  se_candidates res_se1;
  se_candidates res_se2;

  std::size_t the_byte = 0;

  while (rl1 && rl2) {
    {
      const std::lock_guard<std::mutex> lock(abismal_concurrency::read_mutex);
      rl1.load_reads(names1, reads1);
      rl2.load_reads(names2, reads2);
      the_byte = rl1.get_current_byte();
    }

    if (reads1.size() != reads2.size()) {
      throw std::runtime_error("paired-end batch sizes differ. Batch 1: " +
                               std::to_string(reads1.size()) +
                               ", batch 2: " + std::to_string(reads2.size()) +
                               ". Are you sure your paired-end inputs "
                               "have the same number of reads?");
    }

    std::size_t max_batch_read_length = 0;
    update_max_read_length(max_batch_read_length, reads1);
    update_max_read_length(max_batch_read_length, reads2);

    aln.reset(max_batch_read_length);

    const std::size_t n_reads = reads1.size();
    for (std::size_t i = 0; i < n_reads; ++i) {
      const std::uint32_t readlen1 = reads1[i].size();
      const std::uint32_t readlen2 = reads2[i].size();

      res1.reset(readlen1);
      res2.reset(readlen2);
      res_se1.reset(readlen1);
      res_se2.reset(readlen2);

      bests[i].reset(readlen1, readlen2);
      bests_se1[i].reset(readlen1);
      bests_se2[i].reset(readlen2);

      const bool strand_pm_success =
        map_fragments<conv, false, get_strand_code('+', conv),
                      get_strand_code('-', flip_conv(conv))>(
          max_candidates, reads1[i], reads2[i], counter_st,
          (conv == t_rich) ? counter_t_st : counter_a_st, index_st,
          (conv == t_rich) ? index_t_st : index_a_st, genome_st, pread1,
          pread2_rc, packed_pread, cigar1[i], cigar2[i], aln, res1, res2,
          mem_scr1, res_se1, res_se2, bests[i]);

      const bool strand_mp_success =
        map_fragments<!conv, true, get_strand_code('+', flip_conv(conv)),
                      get_strand_code('-', conv)>(
          max_candidates, reads2[i], reads1[i], counter_st,
          (conv == t_rich) ? counter_a_st : counter_t_st, index_st,
          (conv == t_rich) ? index_a_st : index_t_st, genome_st, pread2,
          pread1_rc, packed_pread, cigar2[i], cigar1[i], aln, res2, res1,
          mem_scr1, res_se2, res_se1, bests[i]);

      if (!strand_pm_success && !strand_mp_success) {
        bests[i].reset();
        res_se1.reset();
        res_se2.reset();
      }

      if (!valid_pair(bests[i], reads1[i].size(), reads2[i].size(),
                      cigar_rseq_ops(cigar1[i]), cigar_rseq_ops(cigar2[i])))
        bests[i].reset();

      if (!bests[i].should_report(allow_ambig)) {
        align_se_candidates(pread1, pread1_rc, pread1, pread1_rc,
                            se_element::valid_frac / 2.0, res_se1, bests_se1[i],
                            cigar1[i], aln);

        align_se_candidates(pread2, pread2_rc, pread2, pread2_rc,
                            se_element::valid_frac / 2.0, res_se2, bests_se2[i],
                            cigar2[i], aln);
      }

      select_output(allow_ambig, abismal_index.cl, reads1[i], names1[i],
                    reads2[i], names2[i], cigar1[i], cigar2[i], bests[i],
                    bests_se1[i], bests_se2[i], mr1[i], mr2[i]);
    }

    {
      const std::lock_guard<std::mutex> lock(abismal_concurrency::write_mutex);
      for (std::size_t i = 0; i < n_reads; ++i) {
        if (valid_bam_rec(mr1[i]) && !out.write(hdr, mr1[i]))
          throw std::runtime_error("failed to write bam");
        if (valid_bam_rec(mr2[i]) && !out.write(hdr, mr2[i]))
          throw std::runtime_error("failed to write bam");
      }
    }
    for (std::size_t i = 0; i < n_reads; ++i) {
      if (valid_bam_rec(mr1[i]))
        reset_bam_rec(mr1[i]);
      if (valid_bam_rec(mr2[i]))
        reset_bam_rec(mr2[i]);
      pe_stats.update(allow_ambig, reads1[i], reads2[i], cigar1[i], cigar2[i],
                      bests[i], bests_se1[i], bests_se2[i]);
      cigar1[i].clear();
      cigar2[i].clear();
    }
    if (show_progress) {
      const std::lock_guard<std::mutex> lock(abismal_concurrency::report_mutex);
      if (progress.time_to_report(the_byte))
        progress.report(std::cerr, the_byte);
    }
  }
}

static void
map_paired_ended_rand(const bool show_progress, const bool allow_ambig,
                      const AbismalIndex &abismal_index, ReadLoader &rl1,
                      ReadLoader &rl2, pe_map_stats &pe_stats,
                      bamxx::bam_header &hdr, bamxx::bam_out &out,
                      ProgressBar &progress) {
  const auto counter_st(std::begin(abismal_index.counter));
  const auto counter_t_st(std::begin(abismal_index.counter_t));
  const auto counter_a_st(std::begin(abismal_index.counter_a));

  const auto index_st(std::begin(abismal_index.index));
  const auto index_t_st(std::begin(abismal_index.index_t));
  const auto index_a_st(std::begin(abismal_index.index_a));

  const std::uint32_t max_candidates = abismal_index.max_candidates;

  const genome_iterator genome_st(std::begin(abismal_index.genome));

  std::vector<std::string> names1, reads1;
  std::vector<std::string> names2, reads2;

  std::vector<bam_cigar_t> cigar1;
  std::vector<bam_cigar_t> cigar2;

  std::vector<pe_element> bests;
  std::vector<se_element> bests_se1;
  std::vector<se_element> bests_se2;
  std::vector<bam_rec> mr1;
  std::vector<bam_rec> mr2;

  names1.reserve(ReadLoader::batch_size);
  reads1.reserve(ReadLoader::batch_size);
  cigar1.resize(ReadLoader::batch_size);

  names2.reserve(ReadLoader::batch_size);
  reads2.reserve(ReadLoader::batch_size);
  cigar2.resize(ReadLoader::batch_size);

  bests.resize(ReadLoader::batch_size);
  bests_se1.resize(ReadLoader::batch_size);
  bests_se2.resize(ReadLoader::batch_size);

  mr1.resize(ReadLoader::batch_size);
  mr2.resize(ReadLoader::batch_size);

  Read pread1_t, pread1_t_rc, pread2_t, pread2_t_rc;
  Read pread1_a, pread1_a_rc, pread2_a, pread2_a_rc;
  PackedRead packed_pread;

  pe_candidates res1;
  pe_candidates res2;
  AbismalAlignSimple aln(genome_st);
  std::vector<score_t> mem_scr1(res1.v.size());
  se_candidates res_se1;
  se_candidates res_se2;

  std::size_t the_byte = 0;

  while (rl1 && rl2) {
    {
      const std::lock_guard<std::mutex> lock(abismal_concurrency::read_mutex);
      rl1.load_reads(names1, reads1);
      rl2.load_reads(names2, reads2);
      the_byte = rl1.get_current_byte();
    }

    if (reads1.size() != reads2.size()) {
      throw std::runtime_error("paired-end batch sizes differ. Batch 1: " +
                               std::to_string(reads1.size()) +
                               ", batch 2: " + std::to_string(reads2.size()) +
                               ". Are you sure your paired-end inputs "
                               "have the same number of reads?");
    }

    std::size_t max_batch_read_length = 0;
    update_max_read_length(max_batch_read_length, reads1);
    update_max_read_length(max_batch_read_length, reads2);

    aln.reset(max_batch_read_length);

    const std::size_t n_reads = reads1.size();
    for (std::size_t i = 0; i < n_reads; ++i) {
      const std::uint32_t readlen1 = reads1[i].size();
      const std::uint32_t readlen2 = reads2[i].size();

      res1.reset(readlen1);
      res2.reset(readlen2);
      res_se1.reset(readlen1);
      res_se2.reset(readlen2);
      bests[i].reset(readlen1, readlen2);
      bests_se1[i].reset(readlen1);
      bests_se2[i].reset(readlen2);

      // GS: (1) T/A-rich +/- strand
      const bool richness_ta_strand_pm_success =
        map_fragments<t_rich, false, get_strand_code('+', t_rich),
                      get_strand_code('-', a_rich)>(
          max_candidates, reads1[i], reads2[i], counter_st, counter_t_st,
          index_st, index_t_st, genome_st, pread1_t, pread2_t_rc, packed_pread,
          cigar1[i], cigar2[i], aln, res1, res2, mem_scr1, res_se1, res_se2,
          bests[i]);
      // GS: (2) T/A-rich, -/+ strand
      const bool richness_ta_strand_mp_success =
        map_fragments<a_rich, true, get_strand_code('+', a_rich),
                      get_strand_code('-', t_rich)>(
          max_candidates, reads2[i], reads1[i], counter_st, counter_a_st,
          index_st, index_a_st, genome_st, pread2_a, pread1_a_rc, packed_pread,
          cigar2[i], cigar1[i], aln, res2, res1, mem_scr1, res_se2, res_se1,
          bests[i]);
      // GS: (3) A/T-rich +/- strand
      const bool richness_at_strand_pm_success =
        map_fragments<a_rich, false, get_strand_code('+', a_rich),
                      get_strand_code('-', t_rich)>(
          max_candidates, reads1[i], reads2[i], counter_st, counter_a_st,
          index_st, index_a_st, genome_st, pread1_a, pread2_a_rc, packed_pread,
          cigar1[i], cigar2[i], aln, res1, res2, mem_scr1, res_se1, res_se2,
          bests[i]);
      // GS: (4) A/T-rich, -/+ strand
      const bool richness_at_strand_mp_success =
        map_fragments<t_rich, true, get_strand_code('+', t_rich),
                      get_strand_code('-', a_rich)>(
          max_candidates, reads2[i], reads1[i], counter_st, counter_t_st,
          index_st, index_t_st, genome_st, pread2_t, pread1_t_rc, packed_pread,
          cigar2[i], cigar1[i], aln, res2, res1, mem_scr1, res_se2, res_se1,
          bests[i]);

      if (!richness_ta_strand_pm_success && !richness_ta_strand_mp_success &&
          !richness_at_strand_pm_success && !richness_at_strand_mp_success) {
        bests[i].reset();
        res_se1.reset();
        res_se2.reset();
      }

      if (!valid_pair(bests[i], reads1[i].size(), reads2[i].size(),
                      cigar_rseq_ops(cigar1[i]), cigar_rseq_ops(cigar2[i])))
        bests[i].reset();

      if (!bests[i].should_report(allow_ambig)) {
        align_se_candidates(pread1_t, pread1_t_rc, pread1_a, pread1_a_rc,
                            se_element::valid_frac / 2.0, res_se1, bests_se1[i],
                            cigar1[i], aln);
        align_se_candidates(pread2_t, pread2_t_rc, pread2_a, pread2_a_rc,
                            se_element::valid_frac / 2.0, res_se2, bests_se2[i],
                            cigar2[i], aln);
      }
      select_output(allow_ambig, abismal_index.cl, reads1[i], names1[i],
                    reads2[i], names2[i], cigar1[i], cigar2[i], bests[i],
                    bests_se1[i], bests_se2[i], mr1[i], mr2[i]);
    }

    {
      const std::lock_guard<std::mutex> lock(abismal_concurrency::write_mutex);
      for (std::size_t i = 0; i < n_reads; ++i) {
        if (valid_bam_rec(mr1[i]) && !out.write(hdr, mr1[i]))
          throw std::runtime_error("failed to write bam");
        if (valid_bam_rec(mr2[i]) && !out.write(hdr, mr2[i]))
          throw std::runtime_error("failed to write bam");
      }
    }
    for (std::size_t i = 0; i < n_reads; ++i) {
      if (valid_bam_rec(mr1[i]))
        reset_bam_rec(mr1[i]);
      if (valid_bam_rec(mr2[i]))
        reset_bam_rec(mr2[i]);
      pe_stats.update(allow_ambig, reads1[i], reads2[i], cigar1[i], cigar2[i],
                      bests[i], bests_se1[i], bests_se2[i]);
      cigar1[i].clear();
      cigar2[i].clear();
    }
    if (show_progress) {
      const std::lock_guard<std::mutex> lock(abismal_concurrency::report_mutex);
      if (progress.time_to_report(the_byte))
        progress.report(std::cerr, the_byte);
    }
  }
}

template <const conversion_type conv, const bool random_pbat>
static void
run_paired_ended(const bool show_progress, const bool allow_ambig,
                 const std::string &reads_file1, const std::string &reads_file2,
                 const AbismalIndex &abismal_index, pe_map_stats &pe_stats,
                 bamxx::bam_header &hdr, bamxx::bam_out &out) {
  ReadLoader rl1(reads_file1);
  ReadLoader rl2(reads_file2);
  ProgressBar progress(get_filesize(reads_file1), "mapping reads");

  const auto start_time{abismal_clock::now()};

  std::vector<std::thread> threads;
  for (auto i = 0u; i < abismal_concurrency::n_threads; ++i)
    threads.emplace_back([&] {
      if (random_pbat)
        map_paired_ended_rand(show_progress, allow_ambig, abismal_index, rl1,
                              rl2, pe_stats, hdr, out, progress);
      else
        map_paired_ended<conv>(show_progress, allow_ambig, abismal_index, rl1,
                               rl2, pe_stats, hdr, out, progress);
    });
  for (auto &thread : threads)
    thread.join();

  if (show_progress) {
    std::cerr << "\n";
    const auto stop_time{abismal_clock::now()};
    log_msg("reads mapped: " + std::to_string(rl1.get_current_read()));
    log_msg("total mapping time: " + format_duration(start_time, stop_time));
  }
}

// this is used to fail before reading the index if any input FASTQ
// file does not exist
static inline bool
file_exists(const std::string &filename) {
  return (access(filename.data(), F_OK) == 0);
}

static int
abismal_make_sam_header(const ChromLookup &cl, const int argc, char *argv[],
                        bamxx::bam_header &hdr) {
  assert(cl.names.size() > 2);  // two entries exist for the padding
  assert(cl.starts.size() == cl.names.size() + 1);
  const std::vector<std::string> names(std::begin(cl.names) + 1,
                                       std::end(cl.names) - 1);
  std::vector<std::size_t> sizes(names.size());
  for (std::size_t i = 0; i < names.size(); ++i)
    sizes[i] = cl.starts[i + 2] - cl.starts[i + 1];

  static const std::string SAM_VERSION = "1.0";

  std::ostringstream out;

  // sam version
  out << "@HD" << '\t' << "VN:" << SAM_VERSION << '\n';  // sam version

  // chromosome sizes
  const std::size_t n_chroms = names.size();
  for (std::size_t i = 0; i < n_chroms; ++i)
    out << "@SQ" << '\t' << "SN:" << names[i] << '\t' << "LN:" << sizes[i]
        << '\n';

  // program details
  out << "@PG" << '\t' << "ID:"
      << "ABISMAL" << '\t' << "VN:" << VERSION << '\t';

  // how the program was run
  std::ostringstream the_command;
  copy(argv, argv + argc,
       std::ostream_iterator<const char *>(the_command, " "));
  out << "CL:\"" << the_command.str() << "\"" << std::endl;

  hdr.h = sam_hdr_init();
  return sam_hdr_add_lines(hdr.h, out.str().data(), out.str().size());
}

int
abismal(int argc, char *argv[]) {
  try {

    const std::string version_str =
      std::string("(v") + VERSION + std::string(")");
    const std::string description =
      "map bisulfite converted reads " + version_str;

    bool verbose = false;
    bool GA_conversion = false;
    bool allow_ambig = false;
    bool pbat_mode = false;
    bool random_pbat = false;
    bool write_bam_fmt = false;

    std::uint32_t max_candidates = 0;
    std::string index_file = "";
    std::string genome_file = "";
    std::string outfile("-");
    std::string stats_outfile = "";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<reads-fq1> [<reads-fq2>]");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("index", 'i', "index file", false, index_file);
    opt_parse.add_opt("genome", 'g', "genome file (FASTA)", false, genome_file);
    opt_parse.add_opt("outfile", 'o', "output file", false, outfile);
    opt_parse.add_opt("bam", 'B', "output BAM format", false, write_bam_fmt);
    opt_parse.add_opt("stats", 's', "map statistics file (YAML)", false,
                      stats_outfile);
    opt_parse.add_opt("max-candidates", 'c',
                      "max candidates per seed "
                      "(0 = use index estimate)",
                      false, max_candidates);
    opt_parse.add_opt("min-frag", 'l', "min fragment size (pe mode)", false,
                      pe_element::min_dist);
    opt_parse.add_opt("max-frag", 'L', "max fragment size (pe mode)", false,
                      pe_element::max_dist);
    opt_parse.add_opt("max-distance", 'm', "max fractional edit distance",
                      false, se_element::valid_frac);
    opt_parse.add_opt("ambig", 'a', "report a posn for ambiguous mappers",
                      false, allow_ambig);
    opt_parse.add_opt("pbat", 'P', "input follows the PBAT protocol", false,
                      pbat_mode);
    opt_parse.add_opt("random-pbat", 'R', "input follows random PBAT protocol",
                      false, random_pbat);
    opt_parse.add_opt("a-rich", 'A', "indicates reads are a-rich (se mode)",
                      false, GA_conversion);
    opt_parse.add_opt("threads", 't', "number of threads", false,
                      abismal_concurrency::n_threads);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << std::endl;
      std::cerr << opt_parse.about_message() << std::endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      std::cerr << opt_parse.about_message() << std::endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << "Missing required argument" << std::endl;
      std::cerr << opt_parse.option_missing_message() << std::endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 1 && leftover_args.size() != 2) {
      std::cerr << opt_parse.help_message() << std::endl;
      std::cerr << opt_parse.about_message() << std::endl;
      return EXIT_SUCCESS;
    }
    if (abismal_concurrency::invalid_n_threads()) {
      std::cerr << "Please choose a valid number of threads" << std::endl;
      return EXIT_SUCCESS;
    }
    if (index_file.empty() == genome_file.empty()) {
      std::cerr << "Select one of index file (-i) or genome file (-g)\n";
      return EXIT_SUCCESS;
    }

    const std::string reads_file = leftover_args.front();
    std::string reads_file2;

    if (!file_exists(reads_file)) {
      std::cerr << "cannot open read 1 FASTQ file: " << reads_file << std::endl;
      return EXIT_FAILURE;
    }
    bool paired_end = false;
    if (leftover_args.size() == 2) {
      paired_end = true;
      reads_file2 = leftover_args.back();

      if (!file_exists(reads_file2)) {
        std::cerr << "cannot open read 2 FASTQ file: " << reads_file2
                  << std::endl;
        return EXIT_FAILURE;
      }
    }
    /****************** END COMMAND LINE OPTIONS *****************/

    /****************** BEGIN THREAD VALIDATION *****************/
    const auto n_cores = std::thread::hardware_concurrency();
    abismal_concurrency::n_threads =
      std::min(abismal_concurrency::n_threads, n_cores);
    if (verbose && abismal_concurrency::n_threads > n_cores)
      log_msg("[WARNING] requesting more threads than the "
              "maximum of " +
              std::to_string(n_cores) +
              " processors available in "
              "this device");

    if (verbose)
      log_msg("using " + std::to_string(abismal_concurrency::n_threads) +
              " threads to map reads");
    /****************** END THREAD VALIDATION *****************/

    const bool show_progress = verbose && isatty(fileno(stderr));

    AbismalIndex::VERBOSE = verbose;

    if (verbose) {
      if (paired_end)
        log_msg("input (PE): " + reads_file + ", " + reads_file2);
      else
        log_msg("input (SE): " + reads_file);

      std::string output_msg = "output ";
      output_msg += (write_bam_fmt ? "(BAM): " : "(SAM): ");
      output_msg += (outfile == "-" ? "[stdout]" : outfile);
      log_msg(output_msg);

      if (!stats_outfile.empty())
        log_msg("map statistics (YAML): " + stats_outfile);
    }

    AbismalIndex abismal_index;

    const auto start_time{abismal_clock::now()};
    if (!index_file.empty()) {
      if (verbose)
        log_msg("loading index " + index_file);
      abismal_index.read(index_file);
      const auto stop_time{abismal_clock::now()};
      if (verbose)
        log_msg("loading time: " + format_duration(start_time, stop_time));
    }
    else {
      if (verbose)
        log_msg("indexing genome " + genome_file);
      abismal_index.create_index(genome_file);
      const auto stop_time{abismal_clock::now()};
      if (verbose)
        log_msg("indexing time: " + format_duration(start_time, stop_time));
    }

    if (max_candidates != 0) {
      log_msg("manually setting max_candidates to " +
              std::to_string(max_candidates));
      abismal_index.max_candidates = max_candidates;
    }

    // avoiding opening the stats output file until mapping is done
    se_map_stats se_stats;
    pe_map_stats pe_stats;

    bamxx::bam_out out(outfile, write_bam_fmt);
    if (!out)
      throw std::runtime_error("failed to open output file: " + outfile);

    bamxx::bam_header hdr;
    int ret = abismal_make_sam_header(abismal_index.cl, argc, argv, hdr);

    if (ret < 0)
      throw std::runtime_error("error formatting header");

    if (!out.write(hdr))
      throw std::runtime_error("error writing header");

    if (reads_file2.empty()) {
      if (GA_conversion || pbat_mode)
        run_single_ended<a_rich, false>(show_progress, allow_ambig, reads_file,
                                        abismal_index, se_stats, hdr, out);
      else if (random_pbat)
        run_single_ended<t_rich, true>(show_progress, allow_ambig, reads_file,
                                       abismal_index, se_stats, hdr, out);
      else
        run_single_ended<t_rich, false>(show_progress, allow_ambig, reads_file,
                                        abismal_index, se_stats, hdr, out);
    }
    else {
      if (pbat_mode)
        run_paired_ended<a_rich, false>(show_progress, allow_ambig, reads_file,
                                        reads_file2, abismal_index, pe_stats,
                                        hdr, out);
      else if (random_pbat)
        run_paired_ended<t_rich, true>(show_progress, allow_ambig, reads_file,
                                       reads_file2, abismal_index, pe_stats,
                                       hdr, out);
      else
        run_paired_ended<t_rich, false>(show_progress, allow_ambig, reads_file,
                                        reads_file2, abismal_index, pe_stats,
                                        hdr, out);
    }

    if (!stats_outfile.empty()) {
      std::ofstream stats_of(stats_outfile);
      if (stats_of)
        stats_of << (reads_file2.empty() ? se_stats.tostring()
                                         : pe_stats.tostring(allow_ambig));
      else
        std::cerr << "failed to open stats output file: " << stats_outfile
                  << std::endl;
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
