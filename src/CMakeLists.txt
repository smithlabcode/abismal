# Copyright (C) 2025 Andrew D Smith
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.

if(USE_LIBDEFLATE)
  find_package(LIBDEFLATE REQUIRED)
endif()

find_package(HTSLIB REQUIRED)
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

if(NOT TARGET smithlab_cpp)
  add_subdirectory(smithlab_cpp)
endif()

if(NOT TARGET bamxx)
  add_subdirectory(bamxx)
endif()

add_library(abismal_objs OBJECT
  abismal.cpp
  abismalidx.cpp
  AbismalIndex.cpp
  simreads.cpp
)
target_link_libraries(abismal_objs PUBLIC
  bamxx
  smithlab_cpp
  HTSLIB::HTSLIB
  ZLIB::ZLIB
  Threads::Threads
)
if(USE_LIBDEFLATE)
  target_link_libraries(abismal_objs PUBLIC
    LIBDEFLATE::LIBDEFLATE
  )
endif()
# ADS: CMAKE_BINARY_DIR for config.h
target_include_directories(abismal_objs PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}
)

if(NOT LIBRARY_ONLY)
  add_executable(abismal abismal_main.cpp)
  target_link_libraries(abismal PRIVATE abismal_objs)
endif()
