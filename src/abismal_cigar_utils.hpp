/* Copyright (C) 2023-2025 Andrew D. Smith and Guil Sena
 *
 * Authors: Andrew D. Smith and Guil Sena
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

#ifndef ABISMAL_CIGAR_UTILS_HPP
#define ABISMAL_CIGAR_UTILS_HPP

// This file has equivalent functions and definitions as in HTSlib but
// just what's needed so HTSlib shouldn't be needed for AbismalAlign
// compile

static const int8_t ABISMAL_BAM_CMATCH = 0;
static const int8_t ABISMAL_BAM_CINS = 1;
static const int8_t ABISMAL_BAM_CDEL = 2;
static const int8_t ABISMAL_BAM_CREF_SKIP = 3;
static const int8_t ABISMAL_BAM_CSOFT_CLIP = 4;
static const int8_t ABISMAL_BAM_CHARD_CLIP = 5;
static const int8_t ABISMAL_BAM_CPAD = 6;
static const int8_t ABISMAL_BAM_CEQUAL = 7;
static const int8_t ABISMAL_BAM_CDIFF = 8;
static const int8_t ABISMAL_BAM_CBACK = 9;
static const uint32_t ABISMAL_BAM_CIGAR_SHIFT = 4;

inline uint32_t
abismal_bam_cigar_gen(const uint32_t l, const int8_t o) {
  // ADS: "o" can have -1 for invalid cigar ops
  return (l << ABISMAL_BAM_CIGAR_SHIFT | static_cast<uint32_t>(o));
}

static const uint8_t ABISMAL_BAM_CIGAR_MASK = 0xf;

inline uint8_t
abismal_bam_cigar_op(const uint32_t c) {
  return c & ABISMAL_BAM_CIGAR_MASK;
}

inline uint8_t
abismal_bam_cigar_oplen(const uint32_t c) {
  return c >> ABISMAL_BAM_CIGAR_SHIFT;
}

#endif  // ABISMAL_CIGAR_UTILS_HPP
