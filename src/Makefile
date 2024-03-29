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

# check if recursive clone was done
SMITHLAB_CPP = $(abspath $(dir $(MAKEFILE_LIST)))/smithlab_cpp
STATIC_LIB = $(addprefix $(SRC_ROOT)/, libabismal.a)

BINARIES = abismal abismalidx simreads
OBJECTS = abismal.o abismalidx.o simreads.o AbismalIndex.o

ifeq (,$(wildcard $(SMITHLAB_CPP)/Makefile))
$(error src/smithlab_cpp does not have a Makefile. \
	Did you use --recursive when running git clone?)
endif

PROGS = $(BINARIES)
CXX = g++
CXXFLAGS = -Wall -std=c++11
LDLIBS = -lz
CPPFLAGS = -I $(SMITHLAB_CPP)
OPTFLAGS = -O3
DEBUGFLAGS = -g

COMPILER := $(shell $(CXX) --version)

# clang++ -> add -lomp to link openmp
ifneq '' '$(findstring clang, $(COMPILER))'
	LDLIBS += lomp

# GNU g++ -> add -fopenmp
else ifneq '' '$(findstring g++, $(COMPILER))'
	CXXFLAGS += -fopenmp
endif

LIBDEPS = $(STATIC_LIB) $(addprefix $(SMITHLAB_CPP)/, libsmithlab_cpp.a)

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
else
CXXFLAGS += $(OPTFLAGS)
endif

all: $(PROGS) $(LIBDEPS)
install: $(PROGS)
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(CPPFLAGS)

libabismal.a : $(STATIC_LIB)
$(STATIC_LIB) : $(OBJECTS)
	ar cr $(STATIC_LIB) $^

$(BINARIES): $(LIBDEPS)

%: %_main.cpp $(LIBDEPS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(CPPFLAGS) $(LDLIBS)

clean:
	@-rm -f $(PROGS) $(STATIC_LIB) *.o *.so *.a *~
.PHONY: clean
