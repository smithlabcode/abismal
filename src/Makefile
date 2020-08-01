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
SMITHLAB_CPP=$(abspath $(dir $(MAKEFILE_LIST)))/smithlab_cpp

ifeq (,$(wildcard $(SMITHLAB_CPP)/Makefile))
$(error src/smithlab_cpp does not have a Makefile. \
	      Did you use --recursive when running git clone?)
endif

PROGS = abismal abismalidx simreads

CXX = g++
CXXFLAGS = -Wall -std=c++11 -fopenmp
LDLIBS = -lz
CPPFLAGS = -I $(SMITHLAB_CPP)
OPTFLAGS = -O3
DEBUGFLAGS = -g

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
else
CXXFLAGS += $(OPTFLAGS)
endif

all: $(PROGS)
install: $(PROGS)
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(CPPFLAGS)

$(PROGS): AbismalIndex.o \
	$(addprefix $(SMITHLAB_CPP)/, libsmithlab_cpp.a)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(CPPFLAGS) $(LDLIBS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~
.PHONY: clean
