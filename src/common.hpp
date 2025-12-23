/* Copyright (C) 2025 Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 */

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

template <typename T>
inline auto
revcomp_inplace(T &s) {
  std::transform(std::cbegin(s), std::cend(s), std::begin(s), [](const auto c) {
    return "TNGNNNCNNNNNNNNNNNNANNNNNN"[c - 'A'];
    //      A_C___G____________T______
  });
  std::reverse(std::begin(s), std::end(s));
}

template <typename T>
inline auto
revcomp(T &s) -> T {
  T t(s);
  revcomp_inplace(t);
  return t;
}

struct progress_bar {
  progress_bar(const std::size_t total,
               const std::string message = "completion") :
    total{total}, mid_tag{message} {
    // the 3 below is for the default left_tag and right_tag printed width and
    // the 5 is for the width of the percent (up to 100) plus two pipes ('|')
    bar_width = max_bar_width - std::size(message) - 3 - 5;
    bar = std::string(bar_width, ' ');
  }

  bool
  time_to_report(const std::size_t i) const {
    return std::round((100.0 * std::min(i, total)) / total) > prev;
  }

  void
  report(std::ostream &out, const std::size_t i) {
    prev = std::round((100.0 * std::min(i, total)) / total);
    const std::size_t x =
      std::min(static_cast<std::size_t>(bar_width * (prev / 100.0)), bar_width);
    std::fill_n(std::begin(bar), x, '=');
    out << left_tag << mid_tag << "|" << bar << "|" << std::setw(3) << prev
        << right_tag;
    if (i >= total)
      out << '\n';
  }

  std::size_t total{};
  std::size_t prev{};
  std::size_t bar_width{};
  std::string left_tag = "\r[";
  std::string mid_tag;
  std::string bar;
  std::string right_tag = "%]";

  static const std::size_t max_bar_width = 72;
};

// from 30 April 2020 SAM documentation
// 1    0x1   template having multiple segments in sequencing
// 2    0x2   each segment properly aligned according to the aligner
// 4    0x4   segment unmapped
// 8    0x8   next segment in the template unmapped
// 16   0x10  SEQ being reverse complemented
// 32   0x20  SEQ of the next segment in the template being reverse complemented
// 64   0x40  the first segment in the template
// 128  0x80  the last segment in the template
// 256  0x100 secondary alignment
// 512  0x200 not passing filters, such as platform/vendor quality controls
// 1024 0x400 PCR or optical duplicate
// 2048 0x800 supplementary alignment

namespace samflags {
// ADS: names of flags adjusted to how we typically interpret
static const std::uint16_t read_paired = 0x1;
static const std::uint16_t read_pair_mapped = 0x2;
static const std::uint16_t read_unmapped = 0x4;
static const std::uint16_t mate_unmapped = 0x8;
static const std::uint16_t read_rc = 0x10;
static const std::uint16_t mate_rc = 0x20;
static const std::uint16_t template_first = 0x40;
static const std::uint16_t template_last = 0x80;
static const std::uint16_t secondary_aln = 0x100;
static const std::uint16_t below_quality = 0x200;
static const std::uint16_t pcr_duplicate = 0x400;
static const std::uint16_t supplementary_aln = 0x800;
constexpr auto
check(const std::uint16_t to_check, const std::uint16_t &f) -> bool {
  return to_check & f;
}
constexpr auto
set(std::uint16_t &to_set, const std::uint16_t f) {
  to_set |= f;
}
constexpr auto
unset(std::uint16_t &to_unset, const std::uint16_t f) {
  to_unset &= ~f;
}
}  // namespace samflags
