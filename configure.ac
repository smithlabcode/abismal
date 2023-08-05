dnl This file is part of abismal
dnl
dnl Copyright (C) 2019-2023: Andrew D. Smith and Guilherme de Sena Brandine
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

AC_INIT([abismal], [3.2.0], [andrewds@usc.edu],
        [abismal], [https://github.com/smithlabcode/abismal])

dnl the config.h is not currently #included in the source, and only
dnl used to keep command lines short.
AC_CONFIG_HEADERS([config.h])
AC_PREREQ([2.63])  dnl 4-argument AC_CHECK_HEADER
AM_INIT_AUTOMAKE([subdir-objects foreign])

AC_CONFIG_MACRO_DIR([m4])
AC_LANG(C++)
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])

AM_PROG_AR
AC_PROG_RANLIB

AC_CONFIG_SUBDIRS([src/smithlab_cpp])

dnl check for HTSLib if requested
hts_fail_msg="
Failed to locate HTSLib on your system. Please use the LDFLAGS and
CPPFLAGS variables to specify the directories where the HTSLib library
and headers can be found.
"

dnl check for required libraries
AC_SEARCH_LIBS([hts_version], [hts], [], [AC_MSG_FAILURE([$hts_fail_msg])])


dnl OpenMP happens here
AC_OPENMP([C++])
AS_VAR_IF(OPENMP_CXXFLAGS, [], [
dnl check for the OpenMP library; can't be later
AC_SEARCH_LIBS([omp_get_num_threads], [gomp omp], [],
               [AC_MSG_FAILURE([OpenMP library not found])])
dnl now we get setup for the right OpenMP flags
ADS_OPENMP([], [AC_MSG_FAILURE([OpenMP must be installed to build abismal])])
])dnl end of OpenMP stuff

AC_CONFIG_FILES([
Makefile
])
dnl make the test data files available in the build tree
AC_CONFIG_LINKS([
tests/md5sum.txt:data/md5sum.txt
tests/tRex1.fa:data/tRex1.fa
])
AC_OUTPUT
