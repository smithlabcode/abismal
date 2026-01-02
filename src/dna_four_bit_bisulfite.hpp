/* Copyright (C) 2018-2025 Andrew D. Smith
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

#ifndef DNA_FOUR_BIT_BISULFITE_HPP
#define DNA_FOUR_BIT_BISULFITE_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>  // IWYU pragma: keep
#include <utility>
#include <vector>

enum base_in_byte {
  left,
  right,
};

// clang-format off
constexpr auto dna_four_bit_decoding = std::array{
  'Z', // = 0000 =  0 = {}        = Zero bases
  'A', // = 0001 =  1 = {A}       = Adenine
  'C', // = 0010 =  2 = {C}       = Cytosine
  'M', // = 0011 =  3 = {C,A}     = aMino
  'G', // = 0100 =  4 = {G}       = Guanine
  'R', // = 0101 =  5 = {G,A}     = puRine
  'S', // = 0110 =  6 = {G,C}     = Strong
  'V', // = 0111 =  7 = {G,C,A}   = not T
  'T', // = 1000 =  8 = {T}       = Thymine
  'W', // = 1001 =  9 = {T,A}     = Weak
  'Y', // = 1010 = 10 = {T,C}     = pYramidine
  'H', // = 1011 = 11 = {T,C,A}   = not G
  'K', // = 1100 = 12 = {T,G}     = Keto
  'D', // = 1101 = 13 = {T,G,A}   = not C
  'B', // = 1110 = 14 = {T,G,C}   = not A
  'N'  // = 1111 = 15 = {T,G,C,A} = aNything
};

template <typename uint_type>
constexpr uint_type
get_nibble(const uint_type x, const std::size_t offset) {
  return (x >> (4 * offset)) & 15ul;
}

template <typename uint_type>
constexpr char
decode_dna_four_bit(const uint_type x, const std::size_t offset) {
  return dna_four_bit_decoding[get_nibble(x, offset)];
}

template <class InputItr, class OutputIt>
OutputIt
decode_dna_four_bit(InputItr first, InputItr last, OutputIt d_first) {
  // ADS: assume destination has enough space
  while (first != last) {
    for (std::size_t offset = 0; offset < 16; ++offset)
      *d_first++ = decode_dna_four_bit(*first, offset);
    ++first;
  }
  // if original sequence length is odd and encoding not padded at the front,
  // then the final element in dest will be 'Z'
  return d_first;
}

template <class InCtr, class OutCtr>
void
decode_dna_four_bit(const InCtr &source, OutCtr &dest) {
  // expand out the bytes as pairs (do this backwards in case source == dest)
  const std::size_t source_size = std::size(source);
  dest.resize(16 * source_size);
  std::size_t i = source_size;
  std::size_t j = std::size(dest);
  while (i > 0) {
    dest[--j] = source[--i];
    dest[--j] = source[i];
  }
  for (i = 0; i < std::size(dest); i += 16) {
    for (std::size_t offset = 0; offset < 16; ++offset)
      dest[i + offset] = decode_dna_four_bit(dest[i], offset);
  }
}

/* Sorted by letter
  A = 0001 =  1 = {A}       = Adenine
  B = 1110 = 14 = {T,G,C}   = not A
  C = 0010 =  2 = {C}       = Cytosine
  D = 1101 = 13 = {T,G,A}   = not C
  E        = 15 =
  F        = 15 =
  G = 0100 =  4 = {G}       = Guanine
  H = 1011 = 11 = {T,C,A}   = not G
  I        = 15 =
  J        = 15 =
  K = 1100 = 12 = {T,G}     = Keto
  L        = 15 =
  M = 0011 =  3 = {C,A}     = aMino
  N = 1111 = 15 = {T,G,C,A} = aNything
  O        = 15 =
  P        = 15 =
  Q        = 15 =
  R = 0101 =  5 = {G,A}     = puRine
  S = 0110 =  6 = {G,C}     = Strong
  T = 1000 =  8 = {T}       = Thymine
  U        = 15 =
  V = 0111 =  7 = {G,C,A}   = not T
  W = 1001 =  9 = {T,A}     = Weak
  X        = 15 =
  Y = 1010 = 10 = {T,C}     = pYramidine
  Z = 0000 =  0 = {}        = Zero
*/
constexpr auto dna_four_bit_encoding = std::array{
  /*first*/                                               /*last*/
  /*  0*/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  /* 15*/
  /* 16*/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  /* 31*/
  /* 32*/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  /* 47*/
  /* 48*/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  /* 63*/
  /* 64*/ 0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3, 0, 0,  /* 79*/
  /* 80*/ 0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0, 0,  /* 95*/
  /* 96*/ 0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3, 0, 0,  /*111*/
  /*112*/ 0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0, 0,  /*127*/
};
//      .  A  B  C  D  .  .  G  H  .  .  K  .  M  N  .
//      .  .  R  S  T  .  V  W  .  Y  Z
// clang-format on

template <typename uint_type>
constexpr std::size_t
encode_dna_four_bit(const uint_type x, const std::size_t offset) {
  return (static_cast<std::size_t>(
           dna_four_bit_encoding[static_cast<unsigned>(x)]))
         << (4 * offset);
}

template <class InputItr, class OutputIt>
OutputIt
encode_dna_four_bit(InputItr first, InputItr last, OutputIt d_first) {
  while (first != last) {
    *d_first = 0;
    for (std::size_t i = 0; i < 16 && first != last; ++i)
      *d_first |= encode_dna_four_bit(std::move(*first++), i);
    ++d_first;
  }
  return d_first;
}

// GS: intended to be used as pointer to 4-bit encoding of DNA within a vector
// of std::size_t values
struct genome_four_bit_itr {
  genome_four_bit_itr(const std::vector<std::size_t>::const_iterator itr,
                      const int offset = 0) : itr{itr}, offset{offset} {}

  std::size_t
  operator*() const {
    return (*itr >> (offset << 2)) & 15ul;
  }

  genome_four_bit_itr &
  operator++() {
    offset = (offset + 1) & 15ul;
    itr += (offset == 0);
    return *this;
  }

  genome_four_bit_itr
  operator++(int) {
    genome_four_bit_itr tmp(*this);
    offset = (offset + 1) & 15ul;
    itr += (offset == 0);
    return tmp;
  }

  genome_four_bit_itr &
  operator--() {
    itr -= (offset == 0);

    offset = (offset - 1) & 15ul;
    return *this;
  }

  genome_four_bit_itr
  operator--(int) {
    genome_four_bit_itr tmp(*this);
    itr -= (offset == 0);
    offset = (offset - 1) & 15ul;
    return tmp;
  }

  genome_four_bit_itr
  operator+(const std::size_t step) const {
    // whether the sum of offsets is >= 16
    const bool shift_one_pos =
      ((offset + (static_cast<int>(step) & 15)) & 16) >> 4;

    const int new_offset = (offset + step) & 15;
    return genome_four_bit_itr(itr + step / 16 + shift_one_pos, new_offset);
  }

  bool
  operator!=(const genome_four_bit_itr &rhs) const {
    return itr != rhs.itr || offset != rhs.offset;
  }

  bool
  operator<(const genome_four_bit_itr &rhs) const {
    return itr < rhs.itr || (itr == rhs.itr && offset < rhs.offset);
  }

  bool
  operator<=(const genome_four_bit_itr &rhs) const {
    return itr < rhs.itr || (itr == rhs.itr && offset <= rhs.offset);
  }

  std::vector<std::size_t>::const_iterator itr;
  int offset{};
};

// clang-format off
/* encoding of ASCII characters into T-rich bases, used
 * in encoding reads.
 * A: 0001 = 1
 * C: 0010 = 2
 * G: 0100 = 4
 * T: 1010 = 10 */
constexpr auto encode_base_t_rich = std::array{
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //0
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //17
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //33
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //49
  0, 1, 0, 2, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, //@,A-O
  0, 0, 0, 0, 10,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //P-Z
  0, 1, 0, 2, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, //`,a-o
  0, 0, 0, 0, 10,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //p-z
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

/* encoding of ASCII characters into A-rich bases
 * A: 0101 = 5
 * C: 0010 = 2
 * G: 0100 = 4
 * T: 1000 = 8 */
constexpr auto encode_base_a_rich = std::array{
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //0
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //17
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //33
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //49
  0, 5, 0, 2, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, //@,A-O
  0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //P-Z
  0, 5, 0, 2, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, //`,a-o
  0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //p-z
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

// clang-format on

#endif  // DNA_FOUR_BIT_BISULFITE_HPP
