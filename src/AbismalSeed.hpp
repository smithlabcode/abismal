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

  AbismalSeed(const uint32_t _n_seed_positions,
           const uint32_t _key_weight,
           const uint32_t _n_solid_positions,
           const uint32_t *const _solid_positions) :
    n_seed_positions(_n_seed_positions),
    key_weight(_key_weight),
    n_solid_positions(_n_solid_positions) {
    solid_positions.resize(n_solid_positions);
    copy(_solid_positions, _solid_positions + n_solid_positions,
         std::begin(solid_positions));
  }

  AbismalSeed(const std::string &period,
           const uint32_t key_weight_,
           const uint32_t n_solid_positions_);

  AbismalSeed() {}

  uint32_t n_seed_positions;  // number of positions covered by the seed
  uint32_t key_weight;  // number of hash positions
  uint32_t n_solid_positions; // total solid positions
  std::vector<uint32_t> solid_positions;

  std::string tostring() const;
  void write(const std::string &filename) const;
  void write(std::ostream &out) const;
  void write(FILE *out) const;

  void read(const std::string &filename);
  void read(std::istream &in);
  void read(FILE *in);

  bool operator==(const AbismalSeed &rhs) const {
    return key_weight == rhs.key_weight &&
      solid_positions == rhs.solid_positions;
  }
};

std::ostream &
operator<<(std::ostream &out, const AbismalSeed &s);

namespace seed {
  // the number of shifts is not related to the seed pattern, but is a
  // choice we can make for any seed pattern.
  extern uint32_t n_shifts;

  const uint32_t n_seed_positions = 32;
  const uint32_t key_weight = 26;
  const size_t hash_mask = (1 << seed::key_weight) - 1;
  const uint32_t n_solid_positions = 100;
  const uint32_t solid_positions[] = {
     0,   1,   2,   3,   4,   5,   6,   7,   8,   9,
    10,  11,  12,  13,  14,  15,  16,  17,  18,  19,
    20,  21,  22,  23,  24,  25,  26,  27,  28,  29,
    30,  31,  32,  33,  34,  35,  36,  37,  38,  39,
    40,  41,  42,  43,  44,  45,  46,  47,  48,  49,
    50,  51,  52,  53,  54,  55,  56,  57,  58,  59,
    60,  61,  62,  63,  64,  65,  66,  67,  68,  69,
    70,  71,  72,  73,  74,  75,  76,  77,  78,  79,
    80,  81,  82,  83,  84,  85,  86,  87,  88,  89,
    90,  91,  92,  93,  94,  95,  96,  97,  98,  99
  };
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
