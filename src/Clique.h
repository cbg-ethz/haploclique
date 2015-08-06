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

#ifndef CLIQUE_H_
#define CLIQUE_H_

#include <iostream>

#include "AlignmentRecord.h"
#include "Types.h"

class CliqueFinder;

class Clique {
private:
	size_t leftmost_segment_start;
	size_t rightmost_segment_end;
	size_t alignment_count;
	alignment_set_t* alignment_set;
	CliqueFinder& parent;
	/** Computes insert_start, insert_end, and rightmost_segment_end from all contained alignments. */
	void init();
	/** Debug function */
	void printSet(std::ostream& os);

public:
	Clique(CliqueFinder& parent, size_t index, size_t capacity);
	Clique(CliqueFinder& parent, std::unique_ptr<alignment_set_t>& alignments);
	virtual ~Clique();

	void add(size_t index);

	/** Returns the intersection between nodes in the clique and the given set of nodes. */
	std::unique_ptr<alignment_set_t> intersect(const alignment_set_t& set) const;

	/** Updates the current bit set representation by a new one, where
	 *  new_table[j] := old_table[translation_table[j]]. */
	void translate(size_t* translation_table, size_t table_size, size_t new_capacity);

	size_t leftmostSegmentStart() const { return leftmost_segment_start; }
	size_t rightmostSegmentEnd() const { return rightmost_segment_end; }

	bool contains(const Clique& c);

	size_t size() const { return alignment_count; }

	const alignment_set_t& getAlignmentSet() const { return *alignment_set; }

	std::unique_ptr<std::vector<const AlignmentRecord*> > getAllAlignments() const;

	void computeIntervalIntersection(unsigned int* insert_start, unsigned int* insert_end);

	/** Returns the total coverage of all alignments at the center of the represented intervals, not
	 *  only those that are in the clique. */
	size_t totalCenterCoverage();

	/** Returns coverage at putative breakpoint for all read groups separately. */
	std::unique_ptr<std::vector<size_t> > readGroupWiseCoverage();

	friend std::ostream& operator<<(std::ostream&, const Clique& clique);
};

#endif /* CLIQUE_H_ */
