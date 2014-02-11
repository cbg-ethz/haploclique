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

add_executable ( bam-to-alignment-priors
	bam-to-alignment-priors.cpp
	HistogramBasedDistribution.cpp
	ASTagCalculator.cpp
	Histogram.cpp
	SortedBamReader.cpp
	GroupWiseBamReader.cpp
	BamHelper.cpp
	ShortDnaSequence.cpp
	VariationIndex.cpp
	PositionSet.cpp
	IndelLengthDistribution.cpp
	Variation.cpp
	NamedDnaSequence.cpp
)

add_executable (haploclique
	haploclique.cpp
	CliqueFinder.cpp
	Clique.cpp
	QuasispeciesEdgeCalculator.cpp
	GaussianEdgeCalculator.cpp
	CliqueWriter.cpp
	EdgeWriter.cpp
	AnyDistributionEdgeCalculator.cpp
	ReadSetZTester.cpp
	ReadSetGenericTester.cpp
	ReadGroupAwareEdgeCalculator.cpp
	ReadSetGroupWiseZTester.cpp
	ReadGroups.cpp
	AlignmentRecord.cpp
	ShortDnaSequence.cpp
	CoverageMonitor.cpp
	VariationCaller.cpp
	SingleTrackCoverageMonitor.cpp
	HistogramBasedDistribution.cpp
	Distributions.cpp
	BamHelper.cpp
	VariationIndex.cpp
	IndelLengthDistribution.cpp
	PositionSet.cpp
	Variation.cpp
	NamedDnaSequence.cpp
)

target_link_libraries(bam-to-alignment-priors ${Boost_LIBRARIES} BamTools ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(haploclique ${Boost_LIBRARIES} BamTools)

set(ALL_EXE_TARGETS bam-to-alignment-priors haploclique)

install(TARGETS ${ALL_EXE_TARGETS} RUNTIME DESTINATION bin)
