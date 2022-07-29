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

ifdef NO_MAIN
	ifeq (,$(wildcard $(SMITHLAB_CPP)/Makefile))
	$(error src/smithlab_cpp does not have a Makefile. \
					Did you use --recursive when running git clone?)
	endif
endif

ifndef NO_MAIN
all:
	@make -C $(SMITHLAB_CPP) all
	@make -C src OPT=1
else
all:
	@make -C src OPT=1 NO_MAIN=1
endif

install:
	@make -C src SRC_ROOT=$(SRC_ROOT) OPT=1 install

clean:
	@make -C $(SMITHLAB_CPP) clean
	@make -C src clean
.PHONY: clean
