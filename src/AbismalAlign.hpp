/* Copyright (C) 2019 Andrew D. Smith
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

#ifndef ABISMAL_ALIGN_HPP
#define ABISMAL_ALIGN_HPP

#include <string>
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <array>

#include "cigar_utils.hpp"
#include "dna_four_bit.hpp"

// AbismalAlign class has a templated function for the comparison
// operation between letters in a sequence. This function currently
// returns a boolean value, so only computes longest common
// subsequence.
typedef int16_t score_t;
typedef genome_four_bit_itr genome_iterator;

template <score_t (*scr_fun)(const uint8_t, const uint8_t),
          score_t indel_pen = -1>
struct AbismalAlign {

  AbismalAlign(const genome_iterator &target_start,
               const size_t max_read_length);
  score_t align(const std::vector<uint8_t> &query, uint32_t &t_pos,
                uint32_t &len, std::string &cigar);

  std::vector<score_t> table;
  std::vector<uint8_t> traceback;
  std::vector<uint8_t> cigar_scratch;
  const genome_iterator target;
  const size_t q_sz_max;
  const size_t bw;

  static const uint16_t max_off_diag = 4;
};

template <score_t (*scr_fun)(const uint8_t, const uint8_t),
          score_t indel_pen>
AbismalAlign<scr_fun, indel_pen>::AbismalAlign(
    const genome_iterator &target_start,
    const size_t max_read_length) :
            target(target_start),
            q_sz_max(max_read_length), bw(2*max_off_diag + 1) {
  // size of alignment matrix and traceback matrix is maximum query
  // length times the width of the band around the diagonal
  const size_t n_cells = (q_sz_max + bw)*bw;
  table.resize(n_cells);
  traceback.resize(n_cells, ' ');
  cigar_scratch.resize(2*q_sz_max);
}

// for making the CIGAR string
static const uint8_t left_symbol = 'I';
static const uint8_t above_symbol = 'D';
static const uint8_t diag_symbol = 'M';
static const uint8_t soft_clip_symbol = 'S';

static inline bool
is_deletion(const uint8_t c) {
  return c == above_symbol; // consumes reference
}
static inline bool
is_insertion(const uint8_t c) {
  return c == left_symbol; // does not consume reference
}

static inline void
get_traceback(const size_t n_col,
              const std::vector<score_t> &table,
              const std::vector<uint8_t> &traceback,
              std::vector<uint8_t>::iterator &c_itr,
              size_t &the_row, size_t &the_col) {
  score_t score = table[the_row*n_col + the_col];
  while (score != 0) {
    const uint8_t the_arrow = traceback[the_row*n_col + the_col];
    const bool is_del = is_deletion(the_arrow);
    const bool is_ins = is_insertion(the_arrow);
    the_row -= !is_ins;
    the_col -= is_ins;
    the_col += is_del;
    *c_itr = the_arrow;
    ++c_itr;
    score = table[the_row*n_col + the_col];
  }
}


static score_t
get_best_score(const std::vector<score_t> &table, const size_t n_col,
               const size_t t_shift,
               size_t &best_i, size_t &best_j) {
  const auto best_cell_itr = std::max_element(begin(table), end(table));
  const size_t best_cell = std::distance(std::begin(table), best_cell_itr);
  best_i = best_cell/n_col;
  best_j = best_cell % n_col;
  return *best_cell_itr;
}

template <score_t (*scr_fun)(const uint8_t, const uint8_t),
          class T, class QueryConstItr, class U>
void
from_diag(T next_row, const T next_row_end, T cur_row,
          QueryConstItr query_seq, uint8_t ref_base, U traceback) {
  while (next_row != next_row_end) {
    const score_t score = scr_fun(*query_seq, ref_base) + *cur_row;
    *next_row = std::max(score, *next_row);
    *traceback = (*next_row == score) ? (diag_symbol) : (*traceback);
    ++next_row; ++traceback; ++query_seq; ++cur_row;
  }
}

template <score_t indel_pen, class T, class U>
void
from_above(T above_itr, const T above_end, T target, U traceback) {
  while (above_itr != above_end) {
    const score_t score = *above_itr + indel_pen;
    *target = std::max(score, *target);
    *traceback = (*target == score) ? (above_symbol) : (*traceback);
    ++above_itr; ++target; ++traceback;
  }
}

// ADS: from_left is the same function as from_above, but uses
// different order on arguments, so rewritten to be more intuitive.
template <score_t indel_pen, class T, class U>
void
from_left(T left_itr, T target, const T target_end, U traceback) {
  while (target != target_end) {
    const score_t score = *left_itr + indel_pen;
    *target = std::max(score, *target);
    *traceback = (*target == score) ? (left_symbol) : (*traceback);
    ++left_itr; ++target; ++traceback;
  }
}

template <score_t (*scr_fun)(const uint8_t, const uint8_t), score_t indel_pen>
score_t
AbismalAlign<scr_fun, indel_pen>::align(const std::vector<uint8_t> &qseq,
                                        uint32_t &t_pos, uint32_t &len,
                                        std::string &cigar) {

  std::fill(std::begin(table), std::end(table), 0);
  std::fill(std::begin(traceback), std::end(traceback), ' ');

  const size_t q_sz = qseq.size();

  // GS: non-negative because of padding
  const size_t t_beg = t_pos - ((bw - 1)/2);
  const size_t t_shift = q_sz + bw;

  // points to relevant reference sequence positions
  genome_iterator t_itr = target + t_beg;
  const auto q_itr(std::begin(qseq));
  auto tb_cur(std::begin(traceback));

  // prev and cur point to rows in the alignment matrix
  auto prev(std::begin(table));
  auto cur(prev);
  for (size_t i = 1; i < t_shift; ++i) {
    const size_t left = (i < bw ? bw - i : 0);
    const size_t right = std::min(bw, t_shift - i);

    tb_cur += bw; // next row in traceback
    cur += bw; // next row in aln matrix
    from_diag<scr_fun>(cur + left, cur + right, prev + left,
                       q_itr + (i > bw ? i - bw : 0), *t_itr, tb_cur + left);

    from_above<indel_pen>(prev + left + 1, prev + right, cur + left, tb_cur + left);

    from_left<indel_pen>(cur + left, cur + left + 1, cur + right, tb_cur + left + 1);
    ++t_itr;
    prev += bw; // update previous row
  }

  // locate the end of the alignment as max score
  size_t the_row = 0, the_col = 0;
  const score_t r = get_best_score(table, bw, t_shift, the_row, the_col);

  // GS: unlikely, but possible, case where the score = 0
  if (r == 0) return r;

  // soft clip "S" at the start of the (reverse) uncompressed cigar
  auto c_itr(std::begin(cigar_scratch));
  const size_t soft_clip_bottom = (q_sz + (bw - 1)) - (the_row + the_col);
  std::fill_n(c_itr, soft_clip_bottom, soft_clip_symbol);
  c_itr += soft_clip_bottom;

  get_traceback(bw, table, traceback, c_itr, the_row, the_col);

  // soft clip "S" at the end of the (reverse) uncompressed cigar
  const size_t soft_clip_top = (the_row + the_col) - (bw - 1);
  std::fill_n(c_itr, soft_clip_top, soft_clip_symbol);
  len = qseq.size() - soft_clip_bottom - soft_clip_top;
  c_itr += soft_clip_top;

  // put the uncompressed cigar back in the forward orientation
  std::reverse(std::begin(cigar_scratch), c_itr);

  // GS: max cigar size = 1 operation per base, giving 2 characters each
  cigar.resize(2*q_sz_max);
  compress_cigar(begin(cigar_scratch), c_itr, cigar);
  t_pos = t_beg + the_row;
  return r;
}


// A specific namespace for simple match/mismatch scoring system and a
// 1 -1 -1 scoring scheme for edit distance.
namespace simple_aln {
  static const score_t match = 1;
  static const score_t mismatch = -1;
  static const score_t indel = -1;
  static const score_t min_diffs_to_align = 1;
  static const std::array<score_t, 2> score_lookup = {match, mismatch};

  inline score_t default_score(const uint32_t len, const score_t diffs) {
    return match*(len - diffs) + mismatch*diffs;
  }
  static score_t
  count_deletions(const std::string &cigar) {
    score_t ans = 0;
    for (auto it(begin(cigar)); it != end(cigar); ++it) {
      ans += (extract_op_count(it, end(cigar))) * (*it == 'D');
    }
    return ans;
  }

  static score_t
  count_insertions(const std::string &cigar) {
    score_t ans = 0;
    for (auto it(begin(cigar)); it != end(cigar); ++it) {
      ans += (extract_op_count(it, end(cigar))) * (*it == 'I');
    }
    return ans;
  }

  inline score_t
  mismatch_score(const uint8_t q_base, const uint8_t t_base) {
    return score_lookup[(q_base & t_base) == 0];
  }

  // edit distance as a function of aln_score and len
  inline score_t edit_distance(const score_t scr, const uint32_t len,
                               const std::string &cigar) {
    if (scr == 0) return len;
    const score_t ins = count_insertions(cigar);
    const score_t del = count_deletions(cigar);

    // A = S - (indel_pen) = match*M + mismatch*m
    // B = len - ins = M + m
    // m = (match*(len - ins) - A)/(match - mismatch)
    const score_t A = scr - indel*(ins + del);
    const score_t mism = (match*(len - ins) - A)/(match - mismatch);

    return mism + ins + del;
  }
  inline void make_default_cigar(const uint32_t len, std::string &cigar) {
    cigar = std::to_string(len) + 'M';
  }
};

#endif
