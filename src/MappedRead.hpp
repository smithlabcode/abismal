/*
 *    Part of SMITHLAB_CPP software
 *
 *    Copyright (C) 2010 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MAPPED_READ_HPP
#define MAPPED_READ_HPP

#include "GenomicRegion.hpp"

struct MappedRead {
  MappedRead() {}
  explicit MappedRead(const std::string &line);
  bool has_cigar;
  GenomicRegion r;
  std::string seq;
  std::string scr;
  std::string cigar;
  std::string tostring() const;
};

template <class T> T&
operator>>(T &the_stream, MappedRead &mr) {
  std::string buffer;
  if (getline(the_stream, buffer)) {
    mr = MappedRead(buffer);
  }
  return the_stream;
}

template <class T> T&
operator<<(T &the_stream, const MappedRead &mr) {
  the_stream << mr.tostring();
  return the_stream;
}

#endif
