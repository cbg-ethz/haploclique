# Copyright 2012-2017 Tobias Marschall & Armin Töpfer
#
# This file is part of HaploClique.
#
# HaploClique is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HaploClique is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HaploClique. If not, see <http://www.gnu.org/licenses/>

cmake_policy(SET CMP0048 NEW)
project (HAPLOCLIQUE VERSION 1.3.0 LANGUAGES CXX C)
cmake_minimum_required (VERSION 3.2)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# Main project paths
set(HC_ROOTDIR     ${HAPLOCLIQUE_SOURCE_DIR})
set(HC_SRCDIR      ${HC_ROOTDIR}/src)
set(HC_EXTERNALDIR ${HC_ROOTDIR}/external)
set(HC_TESTDIR     ${HC_ROOTDIR}/test)

# Project configuration
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake ${CMAKE_MODULE_PATH})

# Fixed order, do not touch!
include(hc-ccache)
include(hc-compilerflags)
include(hc-dependencies)

add_subdirectory(${HC_EXTERNALDIR}/BamTools external/bamtools)
add_subdirectory(${HC_EXTERNALDIR}/docopt.cpp external/docopt.cpp)
add_subdirectory(src)
add_subdirectory(test)
