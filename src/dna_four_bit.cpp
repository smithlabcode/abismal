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

#include "dna_four_bit.hpp"

char dna_four_bit_decoding[] = {
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
uint8_t dna_four_bit_encoding[] = {
/*first*/                                                               /*last*/
/*  0*/ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,  /* 15*/
/* 16*/ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,  /* 31*/
/* 32*/ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,  /* 47*/
/* 48*/ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,  /* 63*/
/* 64*/ 15,  1, 14,  2, 13, 15, 15,  4, 11, 15, 15, 12, 15,  3, 15, 15,  /* 79*/
/* 80*/ 15, 15,  5,  6,  8, 15,  7,  9, 15, 10,  0, 15, 15, 15, 15, 15,  /* 95*/
/* 96*/ 15,  1, 14,  2, 13, 15, 15,  4, 11, 15, 15, 12, 15,  3, 15, 15,  /*111*/
/*112*/ 15, 15,  5,  6,  8, 15,  7,  9, 15, 10,  0, 15, 15, 15, 15, 15   /*127*/
};
//       .   A   B   C   D   .   .   G   H   .   .   K   .   M   N   .
//       .   .   R   S   T   .   V   W   .   Y   Z
