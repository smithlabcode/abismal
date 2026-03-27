/* Copyright (C) 2026 Andrew D. Smith
 *
 * This is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This software is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CIGAR_UTILS_HPP_
#define CIGAR_UTILS_HPP_

// ADS: best way to do error checking in the snprintf?

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <initializer_list>
#include <iterator>
#include <string>

inline constexpr bool
consumes_query(const char op) {
  return ((op == 'M') || (op == 'I') || (op == 'S') || (op == '=') ||
          (op == 'X'));
}

inline constexpr bool
consumes_reference(const char op) {
  return ((op == 'M') || (op == 'D') || (op == 'N') || (op == '=') ||
          (op == 'X'));
}

template <class InItr>
static std::size_t
extract_op_count(InItr &itr, const InItr last) {
  std::size_t x = 0;
  while (itr != last && std::isdigit(*itr))
    x = x * 10 + (*itr++ - '0');
  return x;
}

template <class InputItr>
std::size_t
cigar_total_ops(InputItr itr, const InputItr last) {
  std::size_t op_count = 0;
  while (itr != last) {
    op_count += extract_op_count(itr, last);
    if (itr != last)
      *itr++;  // eat the op code
  }
  return op_count;
}

template <class InputItr>
std::size_t
cigar_qseq_ops(InputItr itr, const InputItr last) {
  std::size_t op_count = 0;
  while (itr != last) {
    const std::size_t curr_val = extract_op_count(itr, last);
    op_count += curr_val * consumes_query(*itr++);
  }
  return op_count;
}

template <class T>
std::size_t
cigar_qseq_ops(const T &cigar) {
  return cigar_qseq_ops(std::begin(cigar), std::end(cigar));
}

template <class InputItr>
std::size_t
cigar_rseq_ops(InputItr itr, const InputItr last) {
  std::size_t op_count = 0;
  while (itr != last) {
    const std::size_t curr_val = extract_op_count(itr, last);
    op_count += curr_val * consumes_reference(*itr++);
  }
  return op_count;
}

template <class T>
constexpr std::size_t
cigar_rseq_ops(const T &cigar) {
  return cigar_rseq_ops(std::begin(cigar), std::end(cigar));
}

template <class InputItr>
void
reverse_cigar(InputItr left, const InputItr last) {
  std::reverse(left, last);
  InputItr right(left);
  while (left != last) {
    right = std::find_if(++right, last,  // ++right eats op in rev orientation
                         [](const char c) { return !std::isdigit(c); });
    std::reverse(left, right);
    left = right;
  }
}

template <class T>
void
reverse_cigar(T &cigar) {
  reverse_cigar(std::begin(cigar), std::end(cigar));
}

template <class InputItr>
InputItr
truncate_cigar(InputItr left, const InputItr last,
               const std::size_t target_ops) {
  /* truncate cigar after specified number of total operations */
  InputItr right(left);
  std::size_t prev_ops = 0, curr_ops = 0;
  char op = '\0';
  while (left != last && prev_ops + curr_ops <= target_ops) {
    left = right;
    prev_ops += curr_ops;
    curr_ops = extract_op_count(right, last);
    op = *right++;  // cigar always ends with an op
  }
  if (left != last && target_ops > prev_ops)
    left += std::snprintf(&(*left), 16, "%lu%c", target_ops - prev_ops, op);
  return left;
}

template <class InputItr>
InputItr
truncate_cigar_q(InputItr left, const InputItr last,
                 const std::size_t target_ops) {
  /* truncate cigar after specified number of query operations */
  InputItr right(left);
  std::size_t prev_ops{};
  std::size_t curr_ops{};
  char op{};
  while (left != last && prev_ops + curr_ops <= target_ops) {
    left = right;
    prev_ops += curr_ops;
    curr_ops = extract_op_count(right, last);
    op = *right++;  // cigar always ends with an op
    if (!consumes_query(op))
      curr_ops = 0;
  }
  if (left != last && target_ops > prev_ops)
    left += std::snprintf(&(*left), 16, "%lu%c", target_ops - prev_ops, op);
  return left;
}

template <class T>
void
truncate_cigar_q(T &cigar, const std::size_t target_ops) {
  /* truncate cigar after specified number of query operations */
  auto c_beg(std::begin(cigar));
  auto updated_end = truncate_cigar_q(c_beg, std::end(cigar), target_ops);
  cigar.resize(std::distance(c_beg, updated_end));
}

template <class InputItr>
InputItr
truncate_cigar_r(InputItr left, const InputItr last,
                 const std::size_t target_ops) {
  /* truncate cigar after specified number of reference operations */
  InputItr right(left);
  std::size_t prev_ops{};
  std::size_t curr_ops{};
  char op{};
  while (left != last && prev_ops + curr_ops <= target_ops) {
    left = right;
    prev_ops += curr_ops;
    curr_ops = extract_op_count(right, last);
    op = *right++;  // cigar always ends with an op
    if (!consumes_reference(op))
      curr_ops = 0;
  }
  if (left != last && target_ops > prev_ops)
    left += std::snprintf(&(*left), 16, "%lu%c", target_ops - prev_ops, op);
  return left;
}

template <class T>
void
truncate_cigar_r(T &cigar, const std::size_t target_ops) {
  /* truncate cigar after specified number of reference operations */
  auto c_beg(std::begin(cigar));
  auto updated_end = truncate_cigar_r(c_beg, std::end(cigar), target_ops);
  cigar.resize(std::distance(c_beg, updated_end));
}

template <class InputItr>
InputItr
merge_equal_neighbor_cigar_ops(InputItr rd_itr, const InputItr last) {
  /* renders a cigar valid by merging consecutive identical operations */
  InputItr wr_itr(rd_itr);
  std::size_t prev_op_count = extract_op_count(rd_itr, last);
  char prev_op = *rd_itr++;

  while (rd_itr != last) {
    const std::size_t curr_op_count = extract_op_count(rd_itr, last);
    const char op = *rd_itr++;
    if (op != prev_op) {
      wr_itr += std::snprintf(&(*wr_itr), 16, "%lu%c", prev_op_count, prev_op);
      prev_op_count = 0;
    }
    prev_op_count += curr_op_count;
    prev_op = op;
  }
  wr_itr += std::snprintf(&(*wr_itr), 16, "%lu%c", prev_op_count, prev_op);
  return wr_itr;
}

template <class T>
void
merge_equal_neighbor_cigar_ops(T &cigar) {
  /* renders a cigar valid by merging consecutive identical operations */
  auto c_beg(std::begin(cigar));
  auto updated_end = merge_equal_neighbor_cigar_ops(c_beg, std::end(cigar));
  cigar.resize(std::distance(c_beg, updated_end));
}

template <class InputItr>
void
internal_S_to_M(InputItr rd_itr, const InputItr last) {
  // converts internal soft-clip (S) symbols into match/mismatch (M) symbols
  extract_op_count(rd_itr, last);  // move past first op count
  rd_itr++;                        // move past first op
  while (rd_itr != last) {
    extract_op_count(rd_itr, last);  // past curr op count
    auto curr_op = rd_itr++;         // save curr op iterator
    if (rd_itr != last && *curr_op == 'S') {
      *curr_op = 'M';
    }
  }
}

template <class T>
void
internal_S_to_M(T &cigar) {
  // converts internal soft-clip (S) symbols into match/mismatch (M) symbols
  internal_S_to_M(std::begin(cigar), std::end(cigar));
}

template <class T>
void
terminal_S_to_M(T &cigar) {
  // converts a terminal soft-clip symbol to a match/mismatch (M) symbol
  if (cigar.back() == 'S')
    cigar.back() = 'M';
}

template <class InputItr>
void
initial_S_to_M(InputItr rd_itr, const InputItr last) {
  // converts an initial soft-clip symbol to a match/mismatch (M) symbol
  extract_op_count(rd_itr, last);
  if (*rd_itr == 'S')
    *rd_itr = 'M';
}

template <class T>
void
initial_S_to_M(T &cigar) {
  // converts an initial soft-clip symbol to a match/mismatch (M) symbol
  initial_S_to_M(std::begin(cigar), std::end(cigar));
}

template <class InputItr>
std::size_t
get_soft_clip_size(InputItr rd_itr, const InputItr last) {
  // determine the total size of the soft-clips at both ends
  std::size_t r = 0;
  std::size_t curr_op_count = extract_op_count(rd_itr, last);
  if (*rd_itr++ == 'S')
    r += curr_op_count;
  char curr_op = '\0';
  while (rd_itr != last) {
    curr_op_count = extract_op_count(rd_itr, last);
    curr_op = *rd_itr++;
  }
  return (curr_op == 'S') ? r + curr_op_count : r;
}

template <class T>
constexpr std::size_t
get_soft_clip_size(T &cigar) {
  // determine the total size of the soft-clips at both ends
  return get_soft_clip_size(std::begin(cigar), std::end(cigar));
}

template <class InputItr>
std::size_t
get_soft_clip_size_start(InputItr rd_itr, const InputItr last) {
  // determine the size of the soft-clip at the start of cigar
  const std::size_t curr_op_count = extract_op_count(rd_itr, last);
  return (*rd_itr == 'S') ? curr_op_count : 0;
}

template <class T>
constexpr std::size_t
get_soft_clip_size_start(T &cigar) {
  // determine the size of the soft-clip at the start of cigar
  return get_soft_clip_size_start(std::begin(cigar), std::end(cigar));
}

// turns cigar symbols into compressed format
template <class T1, class T2>
void
compress_cigar(T1 c_itr, const T1 c_end, T2 &cigar) {
  // convert a sequence of cigar symbols into a valid cigar string
  const std::size_t N = cigar.size();
  unsigned j = 0, n = 0;
  char op = *c_itr;
  for (; c_itr != c_end; ++c_itr, ++n) {
    if (*c_itr != op) {
      j += std::snprintf(&cigar[j], N - j, "%d%c", n, op);
      op = *c_itr;
      n = 0;
    }
  }
  j += std::snprintf(&cigar[j], N - j, "%d%c", n, op);
  cigar.resize(j);
}

// turns a cigar string into a string of cigar symbols
template <class T1, class T2>
void
uncompress_cigar(T1 c_itr, const T1 c_end, T2 &cigar) {
  // convert a valid cigar string into a sequence of cigar symbols
  cigar.clear();  // ADS: should this be assumed?
  while (c_itr != c_end) {
    const std::size_t op_count = extract_op_count(c_itr, c_end);
    std::fill_n(std::back_inserter(cigar), op_count, *c_itr++);
  }
}

#endif  // CIGAR_UTILS_HPP_
