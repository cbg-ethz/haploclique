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

#ifndef COVERAGEWRITER_H_
#define COVERAGEWRITER_H_

#include <fstream>

#include "CoverageMonitor.h"

/** Writes coverage to a given file (one line per position). */
class CoverageWriter {
	bool finished;
	std::ofstream ofs;
	CoverageMonitor coverage_monitor;
	size_t rightmost_pos;
	// next position to be written
	size_t pos;
public:
	CoverageWriter(const std::string& filename);
	virtual ~CoverageWriter();

	/** Adds the given alignment. Alignments must be added in sorted order. */
	void addAlignment(const AlignmentRecord& ap);

	/** Completes the started file. After that, no more alignments can be added. */
	void finish();
};

#endif /* COVERAGEWRITER_H_ */
