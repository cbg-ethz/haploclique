/* Copyright 2013 Tobias Marschall
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

#ifndef BAMREADER_H_
#define BAMREADER_H_

#include <iostream>
#include <vector>
#include <map>

/** Abstract base class for classes reading BAM files. */
class BamReader {
public:
	virtual ~BamReader() {}

	/** Returns true if there is another read pair left in the input file. */
	virtual bool hasNext() const = 0;

	/** Reads next read pair. Before calling this method, hasNext() should
	 *  be called. */
	virtual void advance() = 0;

	/** Returns the name of the currently processed read. */
	virtual const std::string& getReadName() const = 0;

	/** Returns a reference to alignments current read. May only be called when
	 *  paired = false was given at construction time.
	 *  The returned reference is valid until the next call to advance(). */
	virtual const std::vector<BamTools::BamAlignment*>& getAlignments() const = 0;

	/** Returns a reference to alignments of first read in current read pair.
	 *  The returned reference is valid until the next call to advance(). */
	virtual const std::vector<BamTools::BamAlignment*>& getAlignmentsFirst() const = 0;

	/** Returns a reference to alignments of second read in current read pair.
	 *  The returned reference is valid until the next call to advance(). */
	virtual const std::vector<BamTools::BamAlignment*>& getAlignmentsSecond() const = 0;

	/** Same as getAlignments, but transfers pointer ownerships. */
	virtual std::auto_ptr<std::vector<BamTools::BamAlignment*> >  releaseAlignments() = 0;

	/** Same as getAlignmentsFirst, but transfers pointer ownerships. */
	virtual std::auto_ptr<std::vector<BamTools::BamAlignment*> >  releaseAlignmentsFirst() = 0;
	
	/** Same as getAlignmentsSecond, but transfers pointer ownerships. */
	virtual std::auto_ptr<std::vector<BamTools::BamAlignment*> >  releaseAlignmentsSecond() = 0;

	virtual bool isUnmapped() const = 0;
	virtual bool isFirstUnmapped() const = 0;
	virtual bool isSecondUnmapped() const = 0;

	virtual bool hasMultipleMappingsButNoXA() const = 0;

	virtual void enableProgressMessages(std::ostream& os, int frequency) = 0;

	virtual long long getSkippedDuplicates() const = 0;
	virtual long long getNonPairedCount() const = 0;

	virtual const BamTools::RefVector& getReferenceData() const = 0;

	virtual BamTools::SamHeader getHeader() const = 0;
};

#endif /* BAMREADER_H_ */
