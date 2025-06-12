dnl This file is part of abismal
dnl
dnl Copyright (C) 2019-2025: Andrew D. Smith and Guilherme de Sena Brandine
dnl
dnl Authors: Andrew D. Smith and Guilherme de Sena Brandine
dnl
dnl This is free software: you can redistribute it and/or modify it
dnl under the terms of the GNU General Public License as published by
dnl the Free Software Foundation, either version 3 of the License, or
dnl (at your option) any later version.
dnl
dnl This software is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl General Public License for more details.

AC_INIT([abismal], [3.2.4], [andrewds@usc.edu],
        [abismal], [https://github.com/smithlabcode/abismal])

dnl the config.h is not currently #included in the source, and only
dnl used to keep command lines short.
AC_CONFIG_HEADERS([config.h])
AC_PREREQ([2.63])  dnl 4-argument AC_CHECK_HEADER
AM_INIT_AUTOMAKE([subdir-objects foreign])

AC_CONFIG_MACRO_DIR([m4])
AC_LANG(C++)
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_17([noext], [mandatory])

AM_PROG_AR
AC_PROG_RANLIB

AC_CONFIG_SUBDIRS([src/smithlab_cpp])

dnl check for HTSLib if requested
hts_fail_msg="
Failed to locate HTSLib on your system. Please use the LDFLAGS and CPPFLAGS
variables to specify the directories where the HTSLib library and headers can
be found.
"

dnl check for required libraries
AC_SEARCH_LIBS([pthread_create], [pthread], [], [AC_MSG_FAILURE(["pthread library not found"])])
AC_SEARCH_LIBS([gzopen], [z], [], [AC_MSG_FAILURE(["Zlib library not found"])])
AC_SEARCH_LIBS([hts_version], [hts], [], [AC_MSG_FAILURE([$hts_fail_msg])], [-lz -lpthread])

dnl Optimize for shorter reads by compiling with a smaller window when making
dnl the index. When abismal loads the index, it compares the value with a
dnl constant so abismalidx and abismal must both be compiled this way to work
dnl together.
AC_ARG_ENABLE([short],
  [AS_HELP_STRING([--enable-short],
                  [optimize for very short reads @<:@no@:>@])],
  [enable_short=yes], [enable_short=no])
AM_CONDITIONAL([ENABLE_SHORT], [test "x$enable_short" = "xyes"])


AC_CONFIG_FILES([
Makefile
])
dnl make the test data files available in the build tree
AC_CONFIG_LINKS([
tests/md5sum.txt:data/md5sum.txt
tests/tRex1.fa:data/tRex1.fa
])
AC_OUTPUT

dnl ADS: use the flag below to get mostly static build
dnl LDFLAGS="-static-libgcc -static-libstdc++"
