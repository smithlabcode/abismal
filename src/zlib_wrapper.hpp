/* Copyright (C) 2019 Meng Zhou and Andrew D. Smith
 *
 * Author: Meng Zhou and Andrew D. Smith
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#ifndef ZLIB_WRAPPER_HPP
#define ZLIB_WRAPPER_HPP

#include <vector>
#include <string>

// ADS: need some way to check if this is installed
#include <zlib.h>

struct igzfstream {
  igzfstream(const std::string filename) :
    fileobj(gzopen(filename.c_str(), "r")) {
      chunk_size = 16 * 1024;
      // gz internal buffer has default size of 8k
      gzbuffer(fileobj, chunk_size);
      buf.resize(chunk_size, ' ');
    }
  ~igzfstream() {gzclose_r(fileobj);}
  operator bool() const {return (fileobj != NULL) && !gzeof(fileobj);}
  char peek();

  size_t chunk_size;
  std::vector<char> buf;
  gzFile fileobj;
};

igzfstream&
getline(igzfstream &in, std::string &line);

igzfstream&
operator>>(igzfstream &in, std::string &line);


struct ogzfstream {
  ogzfstream(const std::string filename) :
    fileobj(gzopen(filename.c_str(), "w")) {}
  ~ogzfstream() {gzclose_w(fileobj);}

  /* Below: eof is not a good error indicator for output streams and
   * we often use this only when checking after opening a file. This
   * function needs improving */
  operator bool() const {return fileobj != NULL;}
  gzFile fileobj;
};

ogzfstream&
operator<<(ogzfstream &out, const std::string &line);

ogzfstream&
operator<<(ogzfstream &out, const char c);

bool has_gz_ext(const std::string &filename);

#endif
