/* Copyright (C) 2019-2023 Andrew D. Smith and Guil Sena
 *
 * Authors: Andrew D. Smith and Guil Sena
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

#ifndef ABISMAL_ALIGN_HPP
#define ABISMAL_ALIGN_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <iterator>
#include <string>
#include <vector>

#include "abismal_cigar_utils.hpp"
#include "dna_four_bit.hpp"

// AbismalAlign class has a templated function for the comparison
// operation between letters in a sequence. This function currently
// returns a boolean value, so only computes longest common
// subsequence.
typedef int16_t score_t;
typedef genome_four_bit_itr genome_iterator;
typedef std::vector<uint32_t> bam_cigar_t;

static inline score_t
count_deletions(const bam_cigar_t &cigar) {
  score_t ans = 0;
  for (const auto &x : cigar)
    if (abismal_bam_cigar_op(x) == ABISMAL_BAM_CDEL)
      ans += abismal_bam_cigar_oplen(x);
  return ans;
}

static inline score_t
count_insertions(const bam_cigar_t &cigar) {
  score_t ans = 0;
  for (const auto &x : cigar)
    if (abismal_bam_cigar_op(x) == ABISMAL_BAM_CINS)
      ans += abismal_bam_cigar_oplen(x);
  return ans;
}

// A specific namespace for simple match/mismatch scoring system and a
// 1 -1 -1 scoring scheme for edit distance.
namespace simple_aln {
  static const score_t match = 2;
  static const score_t mismatch = -3;
  static const score_t indel = -4;
  static const std::array<score_t, 2> score_lookup = {match, mismatch};

  inline score_t default_score(const uint32_t len, const score_t diffs) {
    return match * (len - diffs) + mismatch * diffs;
  }

  inline score_t mismatch_score(const uint8_t q_base, const uint8_t t_base) {
    return score_lookup[(q_base & t_base) == 0];
  }

  // edit distance as a function of aln_score and len
  inline score_t edit_distance(const score_t scr, const uint32_t len,
                               const bam_cigar_t &cigar) {
    if (scr == 0) return len;
    const score_t ins = count_insertions(cigar);
    const score_t del = count_deletions(cigar);

    // A = S - (indel_pen) = match*M + mismatch*m
    // B = len - ins = M + m
    // m = (match*(len - ins) - A)/(match - mismatch)
    const score_t A = scr - indel * (ins + del);
    const score_t mism = (match * (len - ins) - A) / (match - mismatch);

    return mism + ins + del;
  }

  inline score_t best_single_score(const uint32_t readlen) {
    return match * readlen;
  }

  inline score_t best_pair_score(const uint32_t readlen1,
                                 const uint32_t readlen2) {
    return best_single_score(readlen1) + best_single_score(readlen2);
  }
};  // namespace simple_aln

template<score_t (*scr_fun)(const uint8_t, const uint8_t),
         score_t indel_pen = -1>
struct AbismalAlign {
  explicit AbismalAlign(const genome_iterator &target_start)
      : target(target_start), bw(2 * max_off_diag + 1) {}

  template<const bool do_traceback>
  score_t align(const score_t diffs, const score_t max_diffs,
                const std::vector<uint8_t> &query, const uint32_t t_pos);

  void build_cigar_len_and_pos(const score_t diffs, const score_t max_diffs,
                               bam_cigar_t &cigar, uint32_t &len,
                               uint32_t &t_pos);

  void reset(const uint32_t max_read_length);

  std::vector<score_t> table;
  std::vector<int8_t> traceback;
  const genome_iterator target;
  const size_t bw;

  // these are kept because they are needed in both
  // align and build_cigar
  uint16_t q_sz_max;
  uint16_t q_sz;

  static const size_t max_off_diag = 30;
};

template<score_t (*scr_fun)(const uint8_t, const uint8_t), score_t indel_pen>
void
AbismalAlign<scr_fun, indel_pen>::reset(
  const uint32_t max_read_length) {  // uses cigar
  q_sz_max = max_read_length;

  // size of alignment matrix and traceback matrix is maximum query
  // length times the width of the band around the diagonal
  const size_t n_cells = (q_sz_max + bw) * bw;
  table.resize(n_cells);
  traceback.resize(n_cells, -1);  // ADS: -1 no meaning for traceback
}

// for making the CIGAR string
static const uint8_t left_symbol = ABISMAL_BAM_CINS;             // I
static const uint8_t above_symbol = ABISMAL_BAM_CDEL;            // D
static const uint8_t diag_symbol = ABISMAL_BAM_CMATCH;           // M
static const uint8_t soft_clip_symbol = ABISMAL_BAM_CSOFT_CLIP;  // S

static inline bool  // consumes reference
is_deletion(const uint8_t c) {
  return c == above_symbol;
}

static inline bool  // does not consume reference
is_insertion(const uint8_t c) {
  return c == left_symbol;
}

static inline void
get_traceback(const size_t n_col, const std::vector<score_t> &table,
              // const std::vector<uint8_t> &traceback,
              const std::vector<int8_t> &traceback,
              std::vector<uint32_t> &cigar, size_t &the_row, size_t &the_col) {
  int8_t prev_arrow = traceback[the_row * n_col + the_col];
  // int8_t prev_arrow = traceback[the_row * n_col + the_col];
  const bool is_del = is_deletion(prev_arrow);
  const bool is_ins = is_insertion(prev_arrow);
  the_row -= !is_ins;
  the_col -= is_ins;
  the_col += is_del;  // ADS: straight up IS diagonal!
  uint32_t n = 1;
  while (table[the_row * n_col + the_col] > 0) {
    const int8_t arrow = traceback[the_row * n_col + the_col];
    const bool is_del = is_deletion(arrow);
    const bool is_ins = is_insertion(arrow);
    the_row -= !is_ins;
    the_col -= is_ins;
    the_col += is_del;
    if (arrow != prev_arrow) {
      cigar.emplace_back(abismal_bam_cigar_gen(n, prev_arrow));
      n = 0;
    }
    ++n;
    prev_arrow = arrow;
  }
  cigar.emplace_back(abismal_bam_cigar_gen(n, prev_arrow));
}

template<class T> inline void
max16(T &a, const T b) {
  a = ((a > b) ? a : b);
}

template<class T> inline T
min16(const T a, const T b) {
  return ((a < b) ? a : b);
}

static score_t
get_best_score(const std::vector<score_t> &table, const size_t n_cells,
               const size_t n_col, const size_t t_shift, size_t &best_i,
               size_t &best_j) {
  auto best_cell_itr =
    std::max_element(std::begin(table), std::begin(table) + n_cells);
  const size_t best_cell = std::distance(std::begin(table), best_cell_itr);
  best_i = best_cell / n_col;
  best_j = best_cell % n_col;
  return *best_cell_itr;
}

// ADS: it seems like with g++-13, on macos ventura on intel hardware
// the dynamic vectorized optimization of -O3 might be too aggressive
// and makes this function have strange behavior. Placing this pragma
// here helps, and below we restore it to the `-O3` default.
#pragma GCC optimize("vect-cost-model=very-cheap")

template<score_t (*scr_fun)(const uint8_t, const uint8_t), class T,
         class QueryConstItr>
void
from_diag(T next_row, const T next_row_end, T cur_row, QueryConstItr query_seq,
          uint8_t ref_base) {
  while (next_row != next_row_end) {
    const score_t score = scr_fun(*query_seq++, ref_base) + *cur_row++;
    max16(*next_row++, score);
  }
}

template<score_t indel_pen, class T> void
from_above(T above_itr, const T above_end, T target) {
  while (above_itr != above_end) {
    const score_t score = *above_itr++ + indel_pen;
    max16(*target++, score);
  }
}

// ADS: from_left is the same function as from_above, but uses
// different order on arguments, so rewritten to be more intuitive.
template<score_t indel_pen, class T> void
from_left(T left_itr, T target, const T target_end) {
  while (target != target_end) {
    const score_t score = *left_itr++ + indel_pen;
    max16(*target++, score);
  }
}

/********* SAME FUNCTIONS AS ABOVE BUT WITH TRACEBACK ********/
template<score_t (*scr_fun)(const uint8_t, const uint8_t), class T,
         class QueryConstItr, class U>
void
from_diag(T next_row, const T next_row_end, T cur_row, QueryConstItr query_seq,
          uint8_t ref_base, U traceback) {
  while (next_row != next_row_end) {
    const score_t score = scr_fun(*query_seq, ref_base) + *cur_row;
    max16(*next_row, score);
    *traceback = (*next_row == score) ? (diag_symbol) : (*traceback);
    ++traceback;
    ++next_row;
    ++query_seq;
    ++cur_row;
  }
}

template<score_t indel_pen, class T, class U> void
from_above(T above_itr, const T above_end, T target, U traceback) {
  while (above_itr != above_end) {
    const score_t score = *above_itr + indel_pen;
    max16(*target, score);
    *traceback = (*target == score) ? (above_symbol) : (*traceback);
    ++traceback;
    ++above_itr;
    ++target;
  }
}

template<score_t indel_pen, class T, class U> void
from_left(T left_itr, T target, const T target_end, U traceback) {
  while (target != target_end) {
    const score_t score = *left_itr + indel_pen;
    max16(*target, score);
    *traceback = (*target == score) ? (left_symbol) : (*traceback);
    ++traceback;
    ++left_itr;
    ++target;
  }
}

#pragma GCC optimize("vect-cost-model=dynamic")

inline void
make_default_cigar(const uint32_t len, std::string &cigar) {
  cigar = std::to_string(len) + 'M';
}

inline void
make_default_cigar(const uint32_t len, bam_cigar_t &cigar) {
  // ADS: below is = { abismal_bam_cigar_gen(len, 0)};
  cigar = {(len << ABISMAL_BAM_CIGAR_SHIFT)};
}

template<score_t (*scr_fun)(const uint8_t, const uint8_t), score_t indel_pen>
template<const bool do_traceback> score_t
AbismalAlign<scr_fun, indel_pen>::align(const score_t diffs,
                                        const score_t max_diffs,
                                        const std::vector<uint8_t> &qseq,
                                        const uint32_t t_pos) {
  q_sz = qseq.size();
  // edge case: diffs = 0 so alignment is "trivial"
  if (diffs == 0) return simple_aln::best_single_score(q_sz);

  // if diffs is small bw can be reduced
  const size_t bandwidth =
    min16(bw, static_cast<size_t>(2 * min16(diffs, max_diffs) + 1));
  const size_t n_cells = (q_sz + bandwidth) * bandwidth;

  std::fill(std::begin(table), std::begin(table) + n_cells, 0);
  if (do_traceback)
    std::fill(std::begin(traceback), std::begin(traceback) + n_cells, -1);

  // GS: non-negative because of padding. The mapper
  // must ensure t_pos is large enough when calling the function
  const size_t t_beg = t_pos - ((bandwidth - 1) / 2);
  const size_t t_shift = q_sz + bandwidth;

  // points to relevant reference sequence positions
  genome_iterator t_itr = target + t_beg;
  const auto q_itr(std::begin(qseq));
  auto tb_cur(std::begin(traceback));

  // prev and cur point to rows in the alignment matrix
  auto prev(std::begin(table));
  auto cur(prev);

  for (size_t i = 1; i < t_shift; ++i) {
    const size_t left = (i < bandwidth ? bandwidth - i : 0);
    const size_t right = min16(bandwidth, t_shift - i);

    cur += bandwidth;  // next row in aln matrix
    if (do_traceback) {
      tb_cur += bandwidth;  // next row in traceback
      from_diag<scr_fun>(cur + left, cur + right, prev + left,
                         q_itr + (i > bandwidth ? i - bandwidth : 0), *t_itr,
                         tb_cur + left);
      from_above<indel_pen>(prev + left + 1, prev + right, cur + left,
                            tb_cur + left);
      from_left<indel_pen>(cur + left, cur + left + 1, cur + right,
                           tb_cur + left + 1);
    }
    else {
      from_diag<scr_fun>(cur + left, cur + right, prev + left,
                         q_itr + (i > bandwidth ? i - bandwidth : 0), *t_itr);
      from_above<indel_pen>(prev + left + 1, prev + right, cur + left);
      from_left<indel_pen>(cur + left, cur + left + 1, cur + right);
    }
    ++t_itr;
    prev += bandwidth;  // update previous row
  }

  // locate the end of the alignment as max score
  size_t the_row = 0, the_col = 0;
  return get_best_score(table, n_cells, bandwidth, q_sz + bandwidth, the_row,
                        the_col);
}

template<score_t (*scr_fun)(const uint8_t, const uint8_t), score_t indel_pen>
void
AbismalAlign<scr_fun, indel_pen>::build_cigar_len_and_pos(  // uses cigar
  const score_t diffs, const score_t max_diffs, bam_cigar_t &cigar,
  uint32_t &len, uint32_t &t_pos) {
  // locate the end of the alignment as max score
  const size_t bandwidth =
    min16(bw, static_cast<size_t>(2 * min16(diffs, max_diffs) + 1));
  const size_t n_cells = (q_sz + bandwidth) * bandwidth;
  size_t the_row = 0, the_col = 0;
  const score_t r = get_best_score(table, n_cells, bandwidth, q_sz + bandwidth,
                                   the_row, the_col);

  // GS: unlikely, but possible, case where the score = 0, which
  // degenerates CIGAR string below
  if (r == 0 || diffs == 0) {
    make_default_cigar(q_sz, cigar);
    len = q_sz;
    // t_pos does not change in this case
    return;
  }

  // soft clip "S" at the start of the (reverse) uncompressed cigar
  const size_t soft_clip_bottom =
    (q_sz + (bandwidth - 1)) - (the_row + the_col);

  // run traceback, the_row and the_col now point to start of tb
  cigar.clear();
  get_traceback(bandwidth, table, traceback, cigar, the_row, the_col);

  // soft clip "S" at the end of the (reverse) uncompressed cigar
  const size_t soft_clip_top = (the_row + the_col) - (bandwidth - 1);

  // if there is any soft clip at the top, now append it
  if (soft_clip_top > 0)
    cigar.emplace_back(abismal_bam_cigar_gen(soft_clip_top, soft_clip_symbol));

  // cigar was formed working backwards so reverse it
  std::reverse(std::begin(cigar), std::end(cigar));

  // if any soft clip at the bottm, append now after reversing
  if (soft_clip_bottom > 0)
    cigar.emplace_back(
      abismal_bam_cigar_gen(soft_clip_bottom, soft_clip_symbol));

  // length of alignment: query size minus soft clip on both ends
  len = q_sz - soft_clip_bottom - soft_clip_top;

  // ADS: should have documented this better the first time around...
  const size_t t_beg = t_pos - ((bandwidth - 1) / 2);
  t_pos = t_beg + the_row;
}

#endif
