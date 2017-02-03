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
	#NoMeEdgeCalculator.cpp
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

set(libhaplo_SRCS
	AlignmentRecord.h
	AnyDistributionEdgeCalculator.h
	BamHelper.h
	BamReader.h
	BronKerbosch.h
	CLEVER.h
	Clique.h
	CliqueCollector.h
	CliqueFinder.h
	CoverageMonitor.h
	Distributions.h
	EdgeCalculator.h
	EdgeWriter.h
	FastaReader.h
	GaussianEdgeCalculator.h
	Histogram.h
	HistogramBasedDistribution.h
	IndelLengthDistribution.h
	LogWriter.h
	NamedDnaSequence.h
	NewEdgeCalculator.h
	PositionSet.h
	QuasispeciesEdgeCalculator.h
	ReadGroupAwareEdgeCalculator.h
	ReadGroups.h
	ReadSetGenericTester.h
	ReadSetGroupWiseZTester.h
	ReadSetSignificanceTester.h
	ReadSetZTester.h
	ShortDnaSequence.h
	SingleTrackCoverageMonitor.h
	ThreadPool.h
	Types.h
	Variation.h
	VariationIndex.h
	AlignmentRecord.cpp
	AnyDistributionEdgeCalculator.cpp
	BamHelper.cpp
	BronKerbosch.cpp
	CLEVER.cpp
	Clique.cpp
	CoverageMonitor.cpp
	Distributions.cpp
	EdgeWriter.cpp
	GaussianEdgeCalculator.cpp
	Histogram.cpp
	HistogramBasedDistribution.cpp
	IndelLengthDistribution.cpp
	NamedDnaSequence.cpp
	NewEdgeCalculator.cpp
	PositionSet.cpp
	QuasispeciesEdgeCalculator.cpp
	ReadGroupAwareEdgeCalculator.cpp
	ReadGroups.cpp
	ReadSetGenericTester.cpp
	ReadSetGroupWiseZTester.cpp
	ReadSetZTester.cpp
	ShortDnaSequence.cpp
	SingleTrackCoverageMonitor.cpp
	Variation.cpp
	VariationIndex.cpp
	haploclique.cpp
    docopt/docopt.cpp
    docopt/docopt.cpp
    docopt/docopt.h
    docopt/docopt_private.h
    docopt/docopt_util.h
    docopt/docopt_value.h
)

# Declare the library
add_library(libhaplo STATIC
    ${libhaplo_SRCS}
)

# Specify here the include directories exported
# by this library
target_include_directories(libhaplo PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)