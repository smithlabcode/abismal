# This file is part of abismal
#
# Copyright (C) 2018-2020: Andrew D. Smith
#
# Authors: Andrew D. Smith
#
# This is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

EXTRA_DIST = README.md LICENSE docs
ACLOCAL_AMFLAGS = -I m4

SUBDIRS := src/smithlab_cpp

install installdirs: SUBDIRS := $(filter-out src/smithlab_cpp, $(SUBDIRS))
AM_CPPFLAGS = -I $(top_srcdir)/src/smithlab_cpp

CXXFLAGS = -O3  # default has optimization on
AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

LDADD = src/smithlab_cpp/libsmithlab_cpp.a

bin_PROGRAMS = abismal abismalidx simreads

abismal_SOURCES = \
	src/abismal.cpp \
	src/AbismalIndex.cpp \
	src/AbismalAlign.hpp \
	src/AbismalIndex.hpp \
	src/dna_four_bit_bisulfite.hpp

abismalidx_SOURCES = \
	src/abismalidx.cpp \
	src/AbismalIndex.cpp \
	src/AbismalIndex.hpp

simreads_SOURCES = \
	src/simreads.cpp \
	src/AbismalIndex.cpp \
	src/AbismalIndex.hpp
