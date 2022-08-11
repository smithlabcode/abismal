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

if BINARIES_ONLY
SUBDIRS := src/smithlab_cpp
endif

install installdirs: SUBDIRS := $(filter-out src/smithlab_cpp, $(SUBDIRS))
AM_CPPFLAGS = -I $(top_srcdir)/src/smithlab_cpp

CXXFLAGS = -O3 $(OPENMP_CXXFLAGS) # default has optimization on

if LIBRARY_ONLY
lib_LIBRARIES = libabismal.a
else
noinst_LIBRARIES = libabismal.a
endif

libabismal_a_SOURCES = \
	src/abismal.cpp \
	src/abismalidx.cpp \
	src/AbismalIndex.cpp \
	src/simreads.cpp

libabismal_a_SOURCES += \
	src/abismal.hpp \
	src/abismalidx.hpp \
	src/simreads.hpp \
	src/AbismalAlign.hpp \
	src/AbismalIndex.hpp \
	src/dna_four_bit_bisulfite.hpp \
	src/popcnt.hpp

if BINARIES_ONLY
LDADD = libabismal.a src/smithlab_cpp/libsmithlab_cpp.a

bin_PROGRAMS = abismal abismalidx simreads
abismal_SOURCES = src/abismal_main.cpp
abismalidx_SOURCES = src/abismalidx_main.cpp
simreads_SOURCES = src/simreads_main.cpp
endif
