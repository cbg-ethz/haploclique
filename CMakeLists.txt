# Copyright 2012 Tobias Marschall & Armin Töpfer
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
# along with HaploClique.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required (VERSION 2.6) 
project (STRUCTVAR)
SET(CMAKE_C_COMPILER "/usr/bin/gcc-4.9")

list(APPEND CMAKE_MODULE_PATH .) 

IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin)

enable_testing()

IF(${STATIC_LINKING} MATCHES "1")
	message("Statically linked binaries requested.")
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    SET(Boost_USE_STATIC_LIBS ON)
ENDIF()

find_package(Boost 1.36 COMPONENTS program_options iostreams unit_test_framework REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -Wno-long-long -Wno-unknown-pragmas -g -ggdb -O0 -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -Wno-long-long -Wno-unknown-pragmas -O3 -DNDEBUG -std=c++11")
set(CMAKE_CXX_FLAGS_PROFILE "-Wall -Wextra -Wno-long-long -Wno-unknown-pragmas -O3 -g -ggdb -std=c++11")
include_directories(${Boost_INCLUDE_DIR} src external/BamTools/src)

if( "${CMAKE_SYSTEM_NAME}" MATCHES "SunOS" )
	add_definitions( -DSUN_OS )
	set(CMAKE_EXE_LINKER_FLAGS "-lsocket ${CMAKE_EXE_LINKER_FLAGS}")
endif()

add_subdirectory(external/BamTools/src/api)
add_subdirectory(src)

#install(PROGRAMS scripts/haploclique-assembly scripts/postprocess-predictions bwa/bwa samtools/samtools DESTINATION bin)
