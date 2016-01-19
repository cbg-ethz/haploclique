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

add_executable (haploclique
	haploclique.cpp
	CLEVER.cpp
        BronKerbosch.cpp
	Clique.cpp
	QuasispeciesEdgeCalculator.cpp
        NewEdgeCalculator.cpp
	GaussianEdgeCalculator.cpp
	AnyDistributionEdgeCalculator.cpp
	AlignmentRecord.cpp
	ShortDnaSequence.cpp
	HistogramBasedDistribution.cpp
	Distributions.cpp
	VariationIndex.cpp
	IndelLengthDistribution.cpp
	Variation.cpp
	NamedDnaSequence.cpp
    docopt/docopt.cpp
)

IF(${STATIC_LINKING} MATCHES "1")
	target_link_libraries(haploclique ${Boost_LIBRARIES} BamTools-static)
ELSE()
    target_link_libraries(haploclique ${Boost_LIBRARIES} BamTools)
ENDIF()

set(ALL_EXE_TARGETS haploclique)

install(TARGETS ${ALL_EXE_TARGETS} RUNTIME DESTINATION bin)
