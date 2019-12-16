# Copyright (C) 2018-2019 Andrew D. Smith
#
# Authors: Andrew D. Smith
#
# This file is part of ABISMAL.
#
# ABISMAL is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ABISMAL is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.

ifndef SMITHLAB_CPP
$(error Must define SMITHLAB_CPP variable)
endif

PROGS = abismal abismalidx simreads seedtool

CXX = g++
CXXFLAGS = -Wall -std=c++11 -fopenmp
LDLIBS = -lz
CPPFLAGS = -I $(SMITHLAB_CPP)

ifdef OPT
CXXFLAGS += -O3
endif

ifdef DEBUG
CXXFLAGS += -g
endif

all: $(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(CPPFLAGS)

$(PROGS): AbismalIndex.o AbismalSeed.o \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o \
	smithlab_utils.o OptionParser.o zlib_wrapper.o \
	GenomicRegion.o MappedRead.o)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(CPPFLAGS) $(LDLIBS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~
.PHONY: clean
