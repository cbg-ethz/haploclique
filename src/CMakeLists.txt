# Copyright 2012 Tobias Marschall
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

# Get build flags
get_property(HC_COMPILE_FLAGS GLOBAL PROPERTY HC_COMPIPLE_FLAGS_GLOBAL)
get_property(HC_LINK_FLAGS GLOBAL PROPERTY HC_LINK_FLAGS_GLOBAL)
get_property(HC_STATIC GLOBAL PROPERTY HC_STATIC_GLOBAL)

# Common lib
file(GLOB_RECURSE HC_SRCS "*.cpp")
add_library(hc STATIC ${HC_SRCS})
target_include_directories(hc PUBLIC
	${HC_SRCDIR}
	${Boost_INCLUDE_DIR}
	${HC_EXTERNALDIR}/BamTools/src
)
target_link_libraries(hc
	${Boost_LIBRARIES}
)
if (HC_STATIC)
	target_link_libraries(hc BamTools-static)
else()
	target_link_libraries(hc BamTools)
endif()

set_target_properties(hc PROPERTIES COMPILE_FLAGS ${HC_COMPILE_FLAGS})
if (HC_LINK_FLAGS)
	set_target_properties(hc PROPERTIES LINK_FLAGS ${HC_LINK_FLAGS})
endif()

# Haploclique executable
add_executable (haploclique	haploclique.cpp)
set_target_properties(haploclique PROPERTIES COMPILE_FLAGS ${HC_COMPILE_FLAGS})
if (HC_LINK_FLAGS)
	set_target_properties(haploclique PROPERTIES LINK_FLAGS ${HC_LINK_FLAGS})
endif()

set_target_properties(haploclique PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

target_link_libraries(haploclique hc)

install(TARGETS haploclique RUNTIME DESTINATION bin)