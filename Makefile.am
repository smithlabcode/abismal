# This file is part of abismal
#
# Copyright (C) 2018-2019: Andrew D. Smith
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

EXTRA_DIST = README.md
ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src/smithlab_cpp
AM_CPPFLAGS = -I $(top_srcdir)/src -I $(top_srcdir)/src/smithlab_cpp

CXXFLAGS = -O3 -Wall -funroll-loops # default has optimization on
AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

noinst_LIBRARIES = libsmithlab_cpp.a
libsmithlab_cpp_a_SOURCES = \
	src/smithlab_cpp/GenomicRegion.cpp \
	src/smithlab_cpp/MappedRead.cpp \
	src/smithlab_cpp/OptionParser.cpp \
	src/smithlab_cpp/QualityScore.cpp \
	src/smithlab_cpp/bisulfite_utils.cpp \
	src/smithlab_cpp/chromosome_utils.cpp \
	src/smithlab_cpp/sim_utils.cpp \
	src/smithlab_cpp/smithlab_os.cpp \
	src/smithlab_cpp/smithlab_utils.cpp \
	src/smithlab_cpp/dna_four_bit.cpp \
	src/smithlab_cpp/zlib_wrapper.cpp
libsmithlab_cpp_a_SOURCES += \
	src/smithlab_cpp/GenomicRegion.hpp \
	src/smithlab_cpp/MappedRead.hpp \
	src/smithlab_cpp/OptionParser.hpp \
	src/smithlab_cpp/QualityScore.hpp \
	src/smithlab_cpp/bisulfite_utils.hpp \
	src/smithlab_cpp/chromosome_utils.hpp \
	src/smithlab_cpp/sim_utils.hpp \
	src/smithlab_cpp/smithlab_os.hpp \
	src/smithlab_cpp/smithlab_utils.hpp \
	src/smithlab_cpp/dna_four_bit.hpp \
	src/smithlab_cpp/zlib_wrapper.hpp

LDADD = libsmithlab_cpp.a

bin_PROGRAMS = abismal abismalidx simreads seedtool

abismal_SOURCES = \
	src/abismal.cpp \
	src/AbismalIndex.cpp \
	src/AbismalSeed.cpp

abismalidx_SOURCES = \
	src/abismalidx.cpp \
	src/AbismalIndex.cpp \
	src/AbismalSeed.cpp

simreads_SOURCES = \
	src/simreads.cpp \
	src/AbismalIndex.cpp \
	src/AbismalSeed.cpp

seedtool_SOURCES = \
	src/seedtool.cpp \
	src/AbismalSeed.cpp
