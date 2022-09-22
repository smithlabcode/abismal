dnl This file is part of abismal
dnl
dnl Copyright (C) 2019: Andrew D. Smith
dnl
dnl Authors: Andrew D. Smith
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

AC_INIT([abismal], [3.1.0], [andrewds@usc.edu],
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

AX_OPENMP([], [AC_MSG_FAILURE([OpenMP must be installed to build abismal])]) dnl make sure we have openmp for multi-core in abismal

dnl check for the Zlib library
AC_CHECK_LIB([z], [zlibVersion], [],
  [AC_MSG_FAILURE([Zlib must be installed to build abismal])], [])

dnl check for OpenMP
AC_CHECK_LIB([gomp], [omp_get_num_threads], [],
  [AC_MSG_FAILURE([OpenMP must be installed to build abismal])], [])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
