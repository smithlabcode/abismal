/*
 * Part of SMITHLAB_CPP software
 *
 * Copyright (C) 2019 Meng Zhou and Andrew Smith
 *
 * Authors: Meng Zhou and Andrew Smith
 *
 * This is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#ifndef HTSLIB_WRAPPER_HPP
#define HTSLIB_WRAPPER_HPP

#include "smithlab_utils.hpp"
#include "MappedRead.hpp"

#include <string>
#include <vector>
#include <fstream>

extern "C" {
#include <htslib/sam.h>
#include <htslib/hts.h>
}

extern "C" {char check_htslib_wrapper();}

struct SamFlags {
  static const uint16_t read_paired = 0x1;
  static const uint16_t read_pair_mapped = 0x2;
  static const uint16_t read_unmapped = 0x4;
  static const uint16_t mate_unmapped = 0x8;
  static const uint16_t read_rc = 0x10;
  static const uint16_t mate_rc = 0x20;
  static const uint16_t template_first = 0x40;
  static const uint16_t template_second = 0x40;
  static const uint16_t secondary_aln = 0x100;
  static const uint16_t below_quality = 0x200;
  static const uint16_t pcr_duplicate = 0x400;
  static const uint16_t supplementary_aln = 0x800;
};


struct SAMRecord {
  MappedRead mr;
  bool is_Trich;
  bool is_mapping_paired;
  bool is_primary;
  bool is_mapped;
  int seg_len;

  std::string get_name() const {return mr.r.get_name();}
  std::string get_chrom() const {return mr.r.get_chrom();}
};

class SAMReader {
public:
  SAMReader(const std::string &filename, const std::string &mapper);
  ~SAMReader();

  operator bool() const {return good;}

  friend SAMReader &
  operator>>(SAMReader& sam_stream, SAMRecord &samr);

private:
  // internal methods
  bool
  get_SAMRecord(const std::string&, SAMRecord&);
  bool
  get_SAMRecord_bsmap(const std::string&, SAMRecord&);
  bool
  get_SAMRecord_bismark(const std::string&, SAMRecord&);
  bool
  get_SAMRecord_bsseeker(const std::string&, SAMRecord&);
  bool
  get_SAMRecord_general(const std::string&, SAMRecord&);

  // data
  std::string filename;
  std::string mapper;
  bool good;

  htsFile* hts;
  bam_hdr_t *hdr;
  bam1_t *b;
};

SAMReader &
operator>>(SAMReader& sam_stream, SAMRecord &samr);

#endif
