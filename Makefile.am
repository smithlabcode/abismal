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
CXXFLAGS = -O3 -Wall # default has optimization on

bin_PROGRAMS = abismal abismalidx simreads seedtool

abismal_SOURCES = src/abismal.cpp \
	src/AbismalSeed.hpp \
	src/AbismalSeed.cpp \
	src/AbismalIndex.hpp \
	src/AbismalIndex.cpp

abismalidx_SOURCES = src/abismalidx.cpp \
	src/AbismalSeed.hpp \
	src/AbismalSeed.cpp \
	src/AbismalIndex.hpp \
	src/AbismalIndex.cpp

simreads_SOURCES = src/simreads.cpp \
	src/AbismalSeed.hpp \
	src/AbismalSeed.cpp \
	src/AbismalIndex.hpp \
	src/AbismalIndex.cpp

seedtool_SOURCES = src/seedtool.cpp \
	src/AbismalSeed.hpp \
	src/AbismalSeed.cpp \
	src/AbismalIndex.hpp \
	src/AbismalIndex.cpp
