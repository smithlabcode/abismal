/* Copyright (C) 2023-2025 Andrew D. Smith and Guil Sena
 *
 * Authors: Andrew D. Smith and Guil Sena
 *
 * This file is part of ABISMAL.
 *
 * ABISMAL is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * ABISMAL is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 */

#ifndef ABISMAL_CIGAR_UTILS_HPP
#define ABISMAL_CIGAR_UTILS_HPP

#include <cstdint>

// This file has equivalent functions and definitions as in HTSlib but just
// what's needed so HTSlib shouldn't be needed for AbismalAlign compile

static constexpr std::int8_t ABISMAL_BAM_CMATCH = 0;
static constexpr std::int8_t ABISMAL_BAM_CINS = 1;
static constexpr std::int8_t ABISMAL_BAM_CDEL = 2;
static constexpr std::int8_t ABISMAL_BAM_CREF_SKIP = 3;
static constexpr std::int8_t ABISMAL_BAM_CSOFT_CLIP = 4;
static constexpr std::int8_t ABISMAL_BAM_CHARD_CLIP = 5;
static constexpr std::int8_t ABISMAL_BAM_CPAD = 6;
static constexpr std::int8_t ABISMAL_BAM_CEQUAL = 7;
static constexpr std::int8_t ABISMAL_BAM_CDIFF = 8;
static constexpr std::int8_t ABISMAL_BAM_CBACK = 9;
static constexpr std::uint8_t ABISMAL_BAM_CIGAR_MASK = 0xf;
static constexpr std::uint32_t ABISMAL_BAM_CIGAR_SHIFT = 4;

inline auto
abismal_bam_cigar_gen(const std::uint32_t l,
                      const std::int8_t o) -> std::uint32_t {
  // ADS: "o" can have -1 for invalid cigar ops
  return (l << ABISMAL_BAM_CIGAR_SHIFT | static_cast<std::uint32_t>(o));
}

inline auto
abismal_bam_cigar_op(const std::uint32_t c) -> std::uint8_t {
  return c & ABISMAL_BAM_CIGAR_MASK;
}

inline auto
abismal_bam_cigar_oplen(const std::uint32_t c) -> std::uint8_t {
  return c >> ABISMAL_BAM_CIGAR_SHIFT;
}

#endif  // ABISMAL_CIGAR_UTILS_HPP
