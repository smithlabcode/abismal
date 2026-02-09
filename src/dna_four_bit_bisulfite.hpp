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

/* encoding of ASCII characters into T-rich bases, used
 * in encoding reads.
 * A: 0001 = 1
 * C: 0010 = 2
 * G: 0100 = 4
 * T: 1010 = 10 */
constexpr auto encode_base_t_rich = std::array{
  0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 0
  0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 17
  0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 33
  0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 49
  0, 1, 0, 2, 0,  0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0,  //@,A-O
  0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // P-Z
  0, 1, 0, 2, 0,  0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0,  //`,a-o
  0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // p-z
};

/* encoding of ASCII characters into A-rich bases
 * A: 0101 = 5
 * C: 0010 = 2
 * G: 0100 = 4
 * T: 1000 = 8 */
constexpr auto encode_base_a_rich = std::array{
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 0
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 17
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 33
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 49
  0, 5, 0, 2, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0,  //@,A-O
  0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // P-Z
  0, 5, 0, 2, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0,  //`,a-o
  0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // p-z
};

constexpr auto dna_four_bit_decoding = std::array{
  'Z',  // = 0000 =  0 = {}        = Zero bases
  'A',  // = 0001 =  1 = {A}       = Adenine
  'C',  // = 0010 =  2 = {C}       = Cytosine
  'M',  // = 0011 =  3 = {C,A}     = aMino
  'G',  // = 0100 =  4 = {G}       = Guanine
  'R',  // = 0101 =  5 = {G,A}     = puRine
  'S',  // = 0110 =  6 = {G,C}     = Strong
  'V',  // = 0111 =  7 = {G,C,A}   = not T
  'T',  // = 1000 =  8 = {T}       = Thymine
  'W',  // = 1001 =  9 = {T,A}     = Weak
  'Y',  // = 1010 = 10 = {T,C}     = pYramidine
  'H',  // = 1011 = 11 = {T,C,A}   = not G
  'K',  // = 1100 = 12 = {T,G}     = Keto
  'D',  // = 1101 = 13 = {T,G,A}   = not C
  'B',  // = 1110 = 14 = {T,G,C}   = not A
  'N'   // = 1111 = 15 = {T,G,C,A} = aNything
};

enum class base_in_byte : std::uint8_t {
  left,
  right,
};

static constexpr auto nibble_mask = 15u;
static constexpr auto nibble_per_word = 16u;

template <typename uint_type>
constexpr auto
get_nibble(const uint_type x, const std::size_t offset) -> uint_type {
  return (x >> (4 * offset)) & nibble_mask;
}

template <typename uint_type>
constexpr auto
decode_dna_four_bit(const uint_type x, const std::size_t offset) -> char {
  return dna_four_bit_decoding[get_nibble(x, offset)];
}

template <class InputItr, class OutputIt>
auto
decode_dna_four_bit(InputItr first, InputItr last, OutputIt d_first)
  -> OutputIt {
  // ADS: assume destination has enough space
  for (; first != last; ++first)
    for (std::size_t offset = 0; offset < nibble_per_word; ++offset)
      *d_first++ = decode_dna_four_bit(*first, offset);
  // if original sequence length is odd and encoding not padded at the front,
  // then the final element in dest will be 'Z'
  return d_first;
}

template <class InCtr, class OutCtr>
void
decode_dna_four_bit(const InCtr &source, OutCtr &dest) {
  // expand out the bytes as pairs (do this backwards in case source == dest)
  const std::size_t source_size = std::size(source);
  dest.resize(nibble_per_word * source_size);
  std::size_t i = source_size;
  std::size_t j = std::size(dest);
  while (i > 0) {
    dest[--j] = source[--i];
    dest[--j] = source[i];
  }
  for (i = 0; i < std::size(dest); i += nibble_per_word)
    for (std::size_t offset = 0; offset < nibble_per_word; ++offset)
      dest[i + offset] = decode_dna_four_bit(dest[i], offset);
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
  /*  0*/ 0, 0, 0,  0, 0,  0, 0, 0, 0,  0,  0, 0,  0, 0, 0, 0, /* 15*/
  /* 16*/ 0, 0, 0,  0, 0,  0, 0, 0, 0,  0,  0, 0,  0, 0, 0, 0, /* 31*/
  /* 32*/ 0, 0, 0,  0, 0,  0, 0, 0, 0,  0,  0, 0,  0, 0, 0, 0, /* 47*/
  /* 48*/ 0, 0, 0,  0, 0,  0, 0, 0, 0,  0,  0, 0,  0, 0, 0, 0, /* 63*/
  /* 64*/ 0, 1, 14, 2, 13, 0, 0, 4, 11, 0,  0, 12, 0, 3, 0, 0, /* 79*/
  /* 80*/ 0, 0, 5,  6, 8,  0, 7, 9, 0,  10, 0, 0,  0, 0, 0, 0, /* 95*/
  /* 96*/ 0, 1, 14, 2, 13, 0, 0, 4, 11, 0,  0, 12, 0, 3, 0, 0, /*111*/
  /*112*/ 0, 0, 5,  6, 8,  0, 7, 9, 0,  10, 0, 0,  0, 0, 0, 0, /*127*/
};
//        .  A  B   C  D   .  .  G  H   .   .  K   .  M  N  .
//        .  .  R   S  T   .  V  W  .   Y   Z

template <typename uint_type>
constexpr auto
encode_dna_four_bit(const uint_type x, const std::size_t offset)
  -> std::size_t {
  // NOLINTNEXTLINE(*-constant-array-index)
  return static_cast<std::size_t>(dna_four_bit_encoding[x]) << (4 * offset);
}

template <class InputItr, class OutputIt>
auto
encode_dna_four_bit(InputItr first, InputItr last, OutputIt d_first)
  -> OutputIt {
  for (; first != last; ++d_first) {
    *d_first = 0;
    for (std::size_t i = 0; i < nibble_per_word && first != last; ++i)
      *d_first |= encode_dna_four_bit(std::move(*first++), i);
  }
  return d_first;
}

// GS: intended to be used as pointer to 4-bit encoding of DNA within a vector
// of std::size_t values
struct genome_four_bit_itr {
  explicit genome_four_bit_itr(
    const std::vector<std::size_t>::const_iterator itr,
    const std::uint32_t offset = 0) :
    itr{itr},
    offset{offset} {}

  auto
  operator*() const -> std::size_t {
    return (*itr >> (offset << 2)) & nibble_mask;
  }

  auto
  operator++() -> genome_four_bit_itr & {
    offset = (offset + 1) & nibble_mask;
    itr += (offset == 0);
    return *this;
  }

  auto
  operator++(int) -> genome_four_bit_itr {
    genome_four_bit_itr tmp(*this);
    offset = (offset + 1) & nibble_mask;
    itr += (offset == 0);
    return tmp;
  }

  auto
  operator--() -> genome_four_bit_itr & {
    itr -= (offset == 0);
    offset = (offset - 1) & nibble_mask;
    return *this;
  }

  auto
  operator--(int) -> genome_four_bit_itr {
    genome_four_bit_itr tmp(*this);
    itr -= (offset == 0);
    offset = (offset - 1) & nibble_mask;
    return tmp;
  }

  // ADS: careful, this isn't going to work if a negative value is given to
  // the 'step' argument...
  auto
  operator+(const std::size_t step) const -> genome_four_bit_itr {
    // check if the sum of offsets is >= 16
    const bool high_nibble =
      ((offset + (step & nibble_mask)) & nibble_per_word) >> 4;
    const auto new_offset = (offset + step) & nibble_mask;
    const auto word_offset = static_cast<std::int64_t>(step / nibble_per_word);
    return genome_four_bit_itr(itr + word_offset + high_nibble, new_offset);
  }

  auto
  operator!=(const genome_four_bit_itr &rhs) const -> bool {
    return itr != rhs.itr || offset != rhs.offset;
  }

  auto
  operator<(const genome_four_bit_itr &rhs) const -> bool {
    return itr < rhs.itr || (itr == rhs.itr && offset < rhs.offset);
  }

  auto
  operator<=(const genome_four_bit_itr &rhs) const -> bool {
    return itr < rhs.itr || (itr == rhs.itr && offset <= rhs.offset);
  }

  std::vector<std::size_t>::const_iterator itr;
  std::uint32_t offset{};
};

#endif  // DNA_FOUR_BIT_BISULFITE_HPP
