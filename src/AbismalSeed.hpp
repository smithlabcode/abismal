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
    key_weight(_key_weight) {}

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

  // number of positions in searching with the seed, and cannot be
  // longer than n_solid_positions below. ADS: this is allowed to
  // change for debug purposes, but not sure if should be adjustable
  // by the end user.
  extern uint32_t n_seed_positions;

  // number of positions in the hashed portion of the seed
  const uint32_t key_weight = 26;

  const size_t hash_mask = (1 << seed::key_weight) - 1;

  // number of positions used to sort positions based on
  // sequence. This is used in building the index. ADS: this is
  // allowed to change for debug purposes, but not sure if should be
  // adjustable by the end user.
  extern uint32_t n_sorting_positions;
};

// A/T nucleotide to 1-bit value (0100 | 0001 = 5) is for A or G.
inline uint32_t
get_bit(const uint8_t nt) {return (nt & 5) == 0;}

// get the hash value for a k-mer (specified as some iterator/pointer)
// and the encoding for the function above
template <class T>
inline void
get_1bit_hash(T r, uint32_t &k) {
  const auto lim = r + seed::key_weight;
  while (r != lim) {
    k <<= 1;
    k += get_bit(*r);
    ++r;
  }
}

#endif
