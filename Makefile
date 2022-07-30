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
#
SMITHLAB_CPP=$(abspath $(dir $(MAKEFILE_LIST)))/src/smithlab_cpp
SRC_ROOT=$(shell pwd)

	ifeq (,$(wildcard $(SMITHLAB_CPP)/Makefile))
	$(error src/smithlab_cpp does not have a Makefile. \
					Did you use --recursive when running git clone?)
	endif

libabismal.a :
	@$(MAKE) -C src SRC_ROOT=$(SRC_ROOT) libabismal.a

all:
	@$(MAKE) -C $(SMITHLAB_CPP) all
	@$(MAKE) -C src SRC_ROOT=$(SRC_ROOT) OPT=1 all

install:
	@$(MAKE) -C src SRC_ROOT=$(SRC_ROOT) OPT=1 install

clean:
	@$(MAKE) -C $(SMITHLAB_CPP) clean
	@$(MAKE) -C src clean
	@-rm -f *.o *.a *~
.PHONY: clean
