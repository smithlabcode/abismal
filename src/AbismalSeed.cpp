/* Copyright (C) 2018 Andrew D. Smith
 *
 * Author: Andrew D. Smith
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

#include "AbismalSeed.hpp"

#include <iterator>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>

using std::vector;
using std::string;
using std::count;
using std::to_string;
using std::ostream_iterator;
using std::endl;
using std::runtime_error;

uint32_t seed::n_shifts = 3;

AbismalSeed::AbismalSeed(const string &pattern_arg,
                   const uint32_t key_weight_,
                   const uint32_t n_solid_positions_) :
  key_weight(key_weight_),
  n_solid_positions(n_solid_positions_) {

  if (n_solid_positions < key_weight)
    throw runtime_error("key weight less than total solid bits");

  bool found_solid = false;
  vector<bool> period;
  for (size_t i = 0; i < pattern_arg.length(); ++i) {
    if (pattern_arg[i] == '0')
      period.push_back(false);
    else if (pattern_arg[i] == '1') {
      found_solid = true;
      period.push_back(true);
    }
    else throw runtime_error("bad pattern: " + pattern_arg);
  }
  if (!found_solid)
    throw runtime_error("bad pattern: " + pattern_arg);

  const size_t period_size = period.size();
  vector<bool> seed_pattern;
  size_t solid_pos_count = 0;
  for (size_t i = 0; solid_pos_count < n_solid_positions; ++i) {
    seed_pattern.push_back(period[i % period_size]);
    if (seed_pattern.back())
      ++solid_pos_count;
  }

  for (size_t i = 0; i < seed_pattern.size(); ++i)
    if (seed_pattern[i])
      solid_positions.push_back(i);

  n_seed_positions = solid_positions.back() + 1;
}


string
AbismalSeed::tostring() const {
  std::ostringstream oss;
  oss << "const uint32_t n_seed_positions = " << n_seed_positions << ";\n"
      << "const uint32_t key_weight = " << key_weight << ";\n"
      << "const uint32_t n_solid_positions = " << n_solid_positions << ";\n"
      << "const uint32_t solid_positions[] = {\n";
  oss << std::setw(4) << solid_positions[0];
  for (size_t i = 1; i < solid_positions.size(); ++i) {
    oss << ',';
    if (i % 10 == 0) oss << endl;
    oss << std::setw(4) << solid_positions[i];
  }
  oss << endl << "};";
  return oss.str();
}

std::ostream &
operator<<(std::ostream &out, const AbismalSeed &s) {
  return out << s.tostring();
}


void
AbismalSeed::write(const string &filename) const {
  std::ofstream out(filename, std::ios::binary);
  if (!out)
    throw runtime_error("cannot open output file " + filename);
  write(out);
}

void
AbismalSeed::write(std::ostream &out) const {
  if (!out.write((char*)&n_seed_positions, sizeof(uint32_t)) ||
      !out.write((char*)&key_weight, sizeof(uint32_t)) ||
      !out.write((char*)&n_solid_positions, sizeof(uint32_t)) ||
      !out.write((char*)&solid_positions[0],
                 sizeof(uint32_t)*n_solid_positions))
    throw runtime_error("failed writing seed");
}

void
AbismalSeed::write(FILE *out) const {
  if (fwrite((char*)&n_seed_positions, sizeof(uint32_t), 1, out) != 1 ||
      fwrite((char*)&key_weight, sizeof(uint32_t), 1, out) != 1 ||
      fwrite((char*)&n_solid_positions, sizeof(uint32_t), 1, out) != 1 ||
      fwrite((char*)&solid_positions[0], sizeof(uint32_t),
             n_solid_positions, out) != n_solid_positions)
    throw runtime_error("failed writing seed");
}

void
AbismalSeed::read(const string &filename) {
  std::ifstream in(filename, std::ios::binary);
  if (!in)
    throw runtime_error("cannot open input file " + filename);
  read(in);
}

void
AbismalSeed::read(std::istream &in) {
  static const string error_msg("failed loading seed");

  if (!in.read((char*)&n_seed_positions, sizeof(uint32_t)) ||
      !in.read((char*)&key_weight, sizeof(uint32_t)) ||
      !in.read((char*)&n_solid_positions, sizeof(uint32_t)))
    throw runtime_error(error_msg);

  solid_positions.resize(n_solid_positions);
  if (!in.read((char*)&solid_positions[0], sizeof(uint32_t)*n_solid_positions))
    throw runtime_error(error_msg);
}

void
AbismalSeed::read(FILE *in) {
  static const string error_msg("failed loading seed");

  if (fread((char*)&n_seed_positions, sizeof(uint32_t), 1, in) != 1 ||
      fread((char*)&key_weight, sizeof(uint32_t), 1, in) != 1 ||
      fread((char*)&n_solid_positions, sizeof(uint32_t), 1, in) != 1)
    throw runtime_error(error_msg);

  solid_positions.resize(n_solid_positions);
  if (fread((char*)&solid_positions[0], sizeof(uint32_t),
            n_solid_positions, in) != n_solid_positions)
    throw runtime_error(error_msg);
}
