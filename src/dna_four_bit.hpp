/*  Copyright (C) 2020 University of Southern California
 *                     and Andrew D. Smith
 *
 *  Authors: Andrew D. Smith
 *
 *  This is free software: you can redistribute it and/or modify it under the
 *  terms of the GNU General Public License as published by the Free Software
 *  Foundation, either version 3 of the License, or (at your option) any later
 *  version.
 *
 *  This software is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 *  more details.
 */

#ifndef DNA_FOUR_BIT_HPP
#define DNA_FOUR_BIT_HPP

#include <cstdint>
#include <iterator>
#include <vector>

enum base_in_byte { left, right };

extern char dna_four_bit_decoding[16];

template <typename uint_type> constexpr
uint_type
get_low_nibble(const uint_type x) {return x & 15u;}

template <typename uint_type> constexpr
uint_type
get_high_nibble(const uint_type x) {return (x >> 4) & 15u;}

template <typename uint_type> constexpr
char
decode_dna_four_bit_low(const uint_type x) {
  return dna_four_bit_decoding[get_low_nibble(x)];
}

template <typename uint_type> constexpr
char
decode_dna_four_bit_high(const uint_type x) {
  return dna_four_bit_decoding[get_high_nibble(x)];
}

template <typename uint_type> constexpr
char
decode_dna_four_bit(const uint_type x,
                    const base_in_byte b = base_in_byte::left) {
  return b == base_in_byte::left ?
    decode_dna_four_bit_low(x) :
    decode_dna_four_bit_high(x);
}

template<class InputItr, class OutputIt>
OutputIt
decode_dna_four_bit(InputItr first, InputItr last, OutputIt d_first) {
  // ADS: assume destination has enough space
  while (first != last) {
    *d_first++ = decode_dna_four_bit(*first, base_in_byte::left);
    *d_first++ = decode_dna_four_bit(*first, base_in_byte::right);
    ++first;
  }
  // if original sequence length is odd and encoding not padded at the front,
  // then the final element in dest will be 'Z'
  return d_first;
}

template<class InCtr, class OutCtr>
void
decode_dna_four_bit(const InCtr &source, OutCtr &dest) {
  // expand out the bytes as pairs (do this backwards in case source == dest)
  const size_t source_size = source.size();
  dest.resize(2*source_size);
  size_t i = source_size;
  size_t j = dest.size();
  while (i > 0) {
    dest[--j] = source[--i];
    dest[--j] = source[i];
  }
  for (i = 0; i < dest.size(); i += 2) {
    dest[i] = decode_dna_four_bit(dest[i], base_in_byte::left);
    dest[i+1] = decode_dna_four_bit(dest[i+1], base_in_byte::right);
  }
}

extern uint8_t dna_four_bit_encoding[128];

template <typename uint_type> constexpr
uint8_t
encode_dna_four_bit_low(const uint_type x) {
  return dna_four_bit_encoding[static_cast<unsigned>(x)];
}

template <typename uint_type> constexpr
uint8_t
encode_dna_four_bit_high(const uint_type x) {
  return dna_four_bit_encoding[static_cast<unsigned>(x)] << 4;
}

template <typename uint_type> constexpr
uint8_t
encode_dna_four_bit(const uint_type x,
                    const base_in_byte b = base_in_byte::left) {
  return b == base_in_byte::left ?
    encode_dna_four_bit_low(x) :
    encode_dna_four_bit_high(x);
}

template<class InputItr, class OutputIt>
OutputIt
encode_dna_four_bit(InputItr first, InputItr last, OutputIt d_first) {
  while (first != last) {
    *d_first  = encode_dna_four_bit(*first++, base_in_byte::left);
    *d_first |= (first == last ? 0 :
                 encode_dna_four_bit(*first++, base_in_byte::right));
    ++d_first;
  }
  return d_first;
}

// ADS: indented to be used as pointer to 4-bit encoding of DNA within a vector
// of uint8_t values
struct genome_four_bit_itr {
  genome_four_bit_itr(const std::vector<uint8_t>::const_iterator itr_,
                      const bool odd_ = false) : itr(itr_), itr_odd(odd_) {}

  uint8_t operator*() const {
    return (!itr_odd ? *itr : (*itr >> 4)) & 15;
  }
  genome_four_bit_itr& operator++() {
    itr += itr_odd;
    itr_odd ^= 1ul;
    return *this;
  }
  genome_four_bit_itr operator++(int) {
    genome_four_bit_itr tmp(*this);
    itr += itr_odd;
    itr_odd ^= 1ul;
    return tmp;
  }
  genome_four_bit_itr& operator--() {
    itr -= !itr_odd;
    itr_odd ^= 1ul;
    return *this;
  }
  genome_four_bit_itr operator--(int) {
    genome_four_bit_itr tmp(*this);
    itr -= !itr_odd;
    itr_odd ^= 1ul;
    return tmp;
  }
  genome_four_bit_itr operator+(const size_t offset) const {
    const size_t offset_odd = offset & 1ul;
    return genome_four_bit_itr(itr + offset/2 + (itr_odd & offset_odd),
                               itr_odd != offset_odd);
  }
  bool operator!=(const genome_four_bit_itr &rhs) const {
    return itr != rhs.itr || itr_odd != rhs.itr_odd;
  }
  std::vector<uint8_t>::const_iterator itr;
  size_t itr_odd;
};

#endif
