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


#ifndef SINGLETRACKCOVERAGEMONITOR_H
#define SINGLETRACKCOVERAGEMONITOR_H

#include <iostream>
#include <deque>
#include <set>

#include "AlignmentRecord.h"

/** Read-group-agnostic coverage monitor. */
class SingleTrackCoverageMonitor {
private:
	static const int MAX_LOOKAHEAD;
	size_t offset;
	// coverage[i] is the coverage at position offset+i
	std::deque<size_t> coverage;
	// end positions of insert segments that start in the range covered by
	std::multiset<size_t> overhanging;
	/** Extend coverage deque such that the given position is included. */
	void extendTo(size_t pos);
public:
	SingleTrackCoverageMonitor();
	virtual ~SingleTrackCoverageMonitor();

	/** Adds the given alignment. Alignments must be added in sorted order. */
	void addAlignment(const AlignmentRecord& ap);

	/** Returns the maximal local coverage (i.e. "under" the given alignment) that would result
	 *  from adding this alignment pair. Coverage refers to the coverage by insert segments, not
	 *  by the reads themselves. */
	size_t probeAlignment(const AlignmentRecord& ap) const;

	/** Free all coverage information on positions left of the given one. */
	void pruneLeftOf(size_t pos);

	size_t getCoverage(size_t pos);

	friend std::ostream& operator<<(std::ostream& os, const SingleTrackCoverageMonitor& cm);
};

#endif // SINGLETRACKCOVERAGEMONITOR_H
