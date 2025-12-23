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
#include <cstdlib>
#include <iomanip>

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
