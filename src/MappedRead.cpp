/*
 *    Part of SMITHLAB_CPP software
 *
 *    Copyright (C) 2010 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "MappedRead.hpp"
#include "smithlab_utils.hpp"
#include "cigar_utils.hpp"

#include <fstream>
#include <algorithm>
#include <sstream>
#include <string>
#include <iostream>

using std::string;
using std::runtime_error;

template <typename SeqType, typename CigarType>
static void
apply_cigar(SeqType &sequence, const CigarType &cigar) {
  auto cig_it = begin(cigar), cig_end = end(cigar);
  size_t total_ops = cigar_total_ops(begin(cigar), end(cigar));
  size_t ref_ops = cigar_rseq_ops(begin(cigar), end(cigar));

  if (total_ops == 0)
    return;
  size_t cnt, pos = 0;
  sequence.resize(total_ops);
  while (cig_it != cig_end) {
    cnt = extract_op_count(cig_it, cig_end);
    if (*cig_it == 'M' || *cig_it == 'N') pos += cnt;
    else if (*cig_it =='I' || *cig_it =='S') sequence.erase(pos, cnt);
    else if (*cig_it == 'D') {
      sequence.insert(pos, string(cnt, 'N'));
      pos += cnt;
    }
    ++cig_it; //skip op character
  }
  // erases all non-ref ops
  sequence.erase(sequence.begin() + pos, sequence.end());

  if (sequence.size() != ref_ops)
    throw runtime_error("sequence size = " + std::to_string(sequence.size()) +
                        ", rseq ops: " + std::to_string(ref_ops) + " "
                        + sequence + " " + cigar);

  // now that cigar has been applied, destroy it to avoid double application
}

inline bool
is_score(const string &tmp, const string &read) {
  if (tmp.size() != read.size()) return false;
  return
  find_if(tmp.begin(), tmp.end(),
      [](char c) {
        // characters not used in score ascii
        return c == 'M' || c == 'N' || c == 'S';
      }) == tmp.end();
}

MappedRead::MappedRead(const string &line) {
  std::istringstream is;
  is.rdbuf()->pubsetbuf(const_cast<char*>(line.c_str()), line.size());

  string chrom, name, tmp;
  size_t start = 0ul, end = 0ul;
  char strand = '\0';
  double score;
  if (is >> chrom >> start >> tmp) {
    if (find_if(tmp.begin(), tmp.end(),
                [](char c) {return !std::isdigit(c);}) == tmp.end()) {
      end = std::stol(tmp);
      if (!(is >> name >> score >> strand >> seq))
        throw runtime_error("bad GenomicRegion MappedRead file: " + line);

      if (is >> tmp) {
        if (is_score(tmp, seq)) scr = tmp;
        else cigar = tmp;
      }
    }
    else {
      name = tmp;
      if (!(is >> score >> strand >> seq >> tmp))
        throw runtime_error("bad line in MappedRead file: " + line);
      end = start + seq.length();
      if (is_score(tmp, seq)) scr = tmp;
      else cigar = tmp;

    }
    r = GenomicRegion(chrom, start, end, name, score, strand);
    
    // apply and check cigar if empty
    if (!cigar.empty()) {
      apply_cigar(seq, cigar);
      cigar.clear();
    }

    const size_t expected_size = end - start;
    if (seq.length() != expected_size) {
      std::cerr << "seq: " << seq << "\n";
      throw runtime_error("cigar inconsistent with coordinates on line "
                          + line + ". seq size = "
                          + std::to_string(seq.length()) + ". expected: "
                          + std::to_string(expected_size));
    }
  }
  else throw runtime_error("bad line in MappedRead file: " + line);
}

string
MappedRead::tostring() const {
  std::ostringstream oss;
  oss << r; // no chaining for the << of GenomicRegion
  oss << '\t' << seq;

  if (!cigar.empty())
    oss << '\t' << cigar;
  if (!scr.empty())
    oss << '\t' << scr;
  return oss.str();
}

