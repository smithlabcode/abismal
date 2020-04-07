/* Copyright (C) 2018 Andrew D. Smith
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

#ifndef ABISMAL_SEED_HPP
#define ABISMAL_SEED_HPP

#include <string>
#include <vector>
#include <ostream>
#include <istream>
#include <cstdio>

struct AbismalSeed {

  AbismalSeed(const uint32_t _n_solid_positions,
           const uint32_t _key_weight) :
    n_solid_positions(_n_solid_positions),
    key_weight(_key_weight) {
  }

  /*AbismalSeed(//const std::string &period,
           const uint32_t key_weight_,
           const uint32_t n_solid_positions_);*/

  AbismalSeed() {}

  uint32_t n_solid_positions;  // number of positions to sort the genome by
  uint32_t key_weight;  // number of hash positions

  std::string tostring() const;
  void write(const std::string &filename) const;
  void write(std::ostream &out) const;
  void write(FILE *out) const;

  void read(const std::string &filename);
  void read(std::istream &in);
  void read(FILE *in);

  bool operator==(const AbismalSeed &rhs) const {
    return key_weight == rhs.key_weight &&
      n_solid_positions == rhs.n_solid_positions;
  }
};

std::ostream &
operator<<(std::ostream &out, const AbismalSeed &s);

namespace seed {
  // the number of shifts is not related to the seed pattern, but is a
  // choice we can make for any seed pattern.
  extern uint32_t n_shifts;

  //number of positions covered by the seed
  extern uint32_t n_seed_positions;

  const uint32_t n_solid_positions = 200;

  const uint32_t key_weight = 26;
  const uint32_t max_read_length = 200;
  const size_t hash_mask = (1 << seed::key_weight) - 1;
};

// A/T nucleotide to 1-bit number
inline uint32_t
get_bit_4bit(const uint8_t nt) {return (nt & 5) == 0;}

template <class T>
inline void
get_1bit_hash_4bit(T r, uint32_t &k) {
  const auto lim = r + seed::key_weight;
  while (r != lim) {
    k <<= 1;
    k += get_bit_4bit(*r);
    ++r;
  }
}

#endif
