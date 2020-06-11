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

ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I $(top_srcdir)/src

AM_CXXFLAGS = $(OPENMP_CXXFLAGS)
CXXFLAGS = -O3 -Wall -funroll-loops # default has optimization on

bin_PROGRAMS = abismal abismalidx simreads seedtool

abismal_SOURCES = src/abismal.cpp \
	src/AbismalSeed.hpp \
	src/AbismalSeed.cpp \
	src/AbismalIndex.hpp \
	src/AbismalIndex.cpp \
	src/smithlab_os.cpp \
	src/smithlab_os.hpp \
	src/smithlab_utils.cpp \
	src/smithlab_utils.hpp \
	src/OptionParser.cpp \
	src/OptionParser.hpp \
	src/GenomicRegion.cpp \
	src/GenomicRegion.hpp \
	src/dna_four_bit.cpp \
	src/dna_four_bit.hpp \
	src/zlib_wrapper.cpp \
	src/zlib_wrapper.hpp

abismalidx_SOURCES = src/abismalidx.cpp \
	src/AbismalSeed.hpp \
	src/AbismalSeed.cpp \
	src/AbismalIndex.hpp \
	src/AbismalIndex.cpp \
	src/smithlab_os.cpp \
	src/smithlab_os.hpp \
	src/smithlab_utils.cpp \
	src/smithlab_utils.hpp \
	src/OptionParser.cpp \
	src/OptionParser.hpp \
	src/GenomicRegion.cpp \
	src/GenomicRegion.hpp \
	src/dna_four_bit.cpp \
	src/dna_four_bit.hpp


simreads_SOURCES = src/simreads.cpp \
	src/AbismalSeed.hpp \
	src/AbismalSeed.cpp \
	src/AbismalIndex.hpp \
	src/AbismalIndex.cpp \
	src/smithlab_os.cpp \
	src/smithlab_os.hpp \
	src/smithlab_utils.cpp \
	src/smithlab_utils.hpp \
	src/OptionParser.cpp \
	src/OptionParser.hpp \
	src/GenomicRegion.cpp \
	src/GenomicRegion.hpp \
	src/dna_four_bit.cpp \
	src/dna_four_bit.hpp


seedtool_SOURCES = src/seedtool.cpp \
	src/AbismalSeed.hpp \
	src/AbismalSeed.cpp \
	src/AbismalIndex.hpp \
	src/AbismalIndex.cpp \
	src/smithlab_os.cpp \
	src/smithlab_os.hpp \
	src/smithlab_utils.cpp \
	src/smithlab_utils.hpp \
	src/OptionParser.cpp \
	src/OptionParser.hpp \
	src/GenomicRegion.cpp \
	src/GenomicRegion.hpp \
	src/dna_four_bit.cpp \
	src/dna_four_bit.hpp
