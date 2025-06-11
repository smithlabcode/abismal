# This file is part of abismal
#
# Copyright (C) 2018-2024: Andrew D. Smith and Guilherme de Sena Brandine
#
# Authors: Andrew D. Smith and Guilherme de Sena Brandine
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

EXTRA_DIST = \
	README.md \
	LICENSE \
	docs \
	data/tRex1.fa \
        data/md5sum.txt \
	test_scripts/test_abismalidx.test \
	test_scripts/test_simreads.test \
	test_scripts/test_abismal.test \
	test_scripts/test_simreads_pe.test \
	test_scripts/test_abismal_pe.test \
	test_scripts/test_simreads_pbat.test \
	test_scripts/test_abismal_pbat.test \
	test_scripts/test_simreads_rpbat.test \
	test_scripts/test_abismal_rpbat.test \
	test_scripts/test_abismal_threads.test

ACLOCAL_AMFLAGS = -I m4

SUBDIRS := src/smithlab_cpp
install installdirs: SUBDIRS := $(filter-out src/smithlab_cpp, $(SUBDIRS))
AM_CPPFLAGS = -I $(top_srcdir)/src/smithlab_cpp -I $(top_srcdir)/src/bamxx

AM_CXXFLAGS = -Wall -Wextra -Wpedantic -Wno-unknown-attributes
if ENABLE_SHORT
AM_CXXFLAGS += -DENABLE_SHORT
endif

CXXFLAGS = -O3 -DNDEBUG # default has optimization on

noinst_LIBRARIES = libabismal.a

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
	src/popcnt.hpp \
	src/abismal_cigar_utils.hpp \
	src/bamxx/bamxx.hpp

LDADD = libabismal.a src/smithlab_cpp/libsmithlab_cpp.a

bin_PROGRAMS = abismal abismalidx simreads
abismal_SOURCES = src/abismal_main.cpp
abismalidx_SOURCES = src/abismalidx_main.cpp
simreads_SOURCES = src/simreads_main.cpp

TESTS = test_scripts/test_abismalidx.test \
	test_scripts/test_simreads.test \
	test_scripts/test_abismal.test \
	test_scripts/test_simreads_pe.test \
	test_scripts/test_abismal_pe.test \
	test_scripts/test_simreads_pbat.test \
	test_scripts/test_abismal_pbat.test \
	test_scripts/test_simreads_rpbat.test \
	test_scripts/test_abismal_rpbat.test \
	test_scripts/test_abismal_threads.test

TEST_EXTENSIONS = .test

test_scripts/test_abismal.log: \
	test_scripts/test_abismalidx.log \
	test_scripts/test_simreads.log
test_scripts/test_abismal_pe.log: \
	test_scripts/test_abismalidx.log \
	test_scripts/test_simreads_pe.log
test_scripts/test_abismal_pbat.log: \
	test_scripts/test_abismalidx.log \
	test_scripts/test_simreads_pbat.log
test_scripts/test_abismal_rpbat.log: \
	test_scripts/test_abismalidx.log \
	test_scripts/test_simreads_rpbat.log
test_scripts/test_abismal_threads.log: \
	test_scripts/test_abismalidx.log \
	test_scripts/test_simreads.log \
	test_scripts/test_simreads_pe.log \
	test_scripts/test_simreads_pbat.log \
	test_scripts/test_simreads_rpbat.log

CLEANFILES = \
    tests/tRex1.idx \
    tests/reads_1.fq \
    tests/reads.mstats \
    tests/reads.sam \
    tests/reads_pe_1.fq \
    tests/reads_pe_2.fq \
    tests/reads_pe.mstats \
    tests/reads_pe.sam \
    tests/reads_pbat_pe_1.fq \
    tests/reads_pbat_pe_2.fq \
    tests/reads_pbat_pe.mstats \
    tests/reads_pbat_pe.sam \
    tests/reads_rpbat_pe_1.fq \
    tests/reads_rpbat_pe_2.fq \
    tests/reads_rpbat_pe.mstats \
    tests/reads_rpbat_pe.sam
