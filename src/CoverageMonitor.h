/* Copyright 2012 Tobias Marschall
 * 
 * This file is part of HaploClique.
 * 
 * HaploClique is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HaploClique is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HaploClique.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef COVERAGEMONITOR_H_
#define COVERAGEMONITOR_H_

#include <iostream>
#include <vector>

#include "AlignmentRecord.h"
#include "SingleTrackCoverageMonitor.h"
#include "ReadGroups.h"

/** Class that monitors the local coverage while a sorted set of alignments is being read. */
class CoverageMonitor {
private:
	SingleTrackCoverageMonitor global_monitor;
	std::vector<SingleTrackCoverageMonitor> readwise_monitors;
	const ReadGroups* read_groups;
public:
	CoverageMonitor(const ReadGroups* read_groups = 0);
	virtual ~CoverageMonitor();

	/** Adds the given alignment. Alignments must be added in sorted order. */
	void addAlignment(const AlignmentRecord& ap);

	/** Returns the maximal local coverage (i.e. "under" the given alignment) that would result
	 *  from adding this alignment pair. Coverage refers to the coverage by insert segments, not
	 *  by the reads themselves. */
	size_t probeAlignment(const AlignmentRecord& ap) const;

	/** Free all coverage information on positions left of the given one. */
	void pruneLeftOf(size_t pos);

	size_t getCoverage(size_t pos);

	/** Returns true if read-group-wise coverage is available. */
	bool hasReadGroups() const;

	/** Returns coverage of reads from the given read group. */
	size_t getReadGroupCoverage(size_t read_group, size_t pos);

	/** Returns vector of coverages for all read groups. */
	std::unique_ptr<std::vector<size_t> > getReadGroupCoverages(size_t pos);

	friend std::ostream& operator<<(std::ostream& os, const CoverageMonitor& cm);
};

#endif /* COVERAGEMONITOR_H_ */
