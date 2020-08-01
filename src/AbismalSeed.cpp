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
uint32_t seed::n_seed_positions = 32;
uint32_t seed::n_sorting_positions = 1000;

string
AbismalSeed::tostring() const {
  std::ostringstream oss;
  oss << "const uint32_t n_solid_positions = " << n_solid_positions << ";\n"
      << "const uint32_t key_weight = " << key_weight << ";\n";
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
  if (!out.write((char*)&n_solid_positions, sizeof(uint32_t)) ||
      !out.write((char*)&key_weight, sizeof(uint32_t)))
    throw runtime_error("failed writing seed");
}

void
AbismalSeed::write(FILE *out) const {
  if (fwrite((char*)&n_solid_positions, sizeof(uint32_t), 1, out) != 1 ||
      fwrite((char*)&key_weight, sizeof(uint32_t), 1, out) != 1)
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

  if (!in.read((char*)&n_solid_positions, sizeof(uint32_t)) ||
      !in.read((char*)&key_weight, sizeof(uint32_t)))
    throw runtime_error(error_msg);
}

void
AbismalSeed::read(FILE *in) {
  static const string error_msg("failed loading seed");

  if (fread((char*)&n_solid_positions, sizeof(uint32_t), 1, in) != 1 ||
      fread((char*)&key_weight, sizeof(uint32_t), 1, in) != 1)
    throw runtime_error(error_msg);
}
