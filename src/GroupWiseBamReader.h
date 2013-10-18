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

#ifndef GROUPWISEBAMREADER_H_
#define GROUPWISEBAMREADER_H_

#include <iostream>
#include <vector>
#include <boost/unordered_set.hpp>

#include <bamtools/api/BamReader.h>

#include "BamHelper.h"
#include "BamReader.h"

/** Reads a file in BAM format and returns groups of alignments belonging to
 *  the same read. Expects alignments to be grouped by read name in input BAM
 *  and throws an exception otherwise. */
class GroupWiseBamReader : public BamReader {
private:
	BamTools::BamReader bam_reader;
	std::string read_name;
	// alignments of the first and second read in the current read pair
	std::vector<BamTools::BamAlignment*>* alignments1;
	std::vector<BamTools::BamAlignment*>* alignments2;
	// read names to verify that reads are properly grouped
	boost::unordered_set<std::string> read_names;
	// count all processed alignments
	long long counter;
	std::ostream* progress_messages_os;
	int progress_message_frequency;
	BamTools::BamAlignment* next_read_aln;
	bool finished;
	bool skip_duplicates;
	long long skipped_duplicates;
	bool paired;
	bool unmapped1;
	bool unmapped2;
	bool inconsistentXTags;
	std::string inconsistentXTagsMessage;
	bool multiplyMappedButNoXA;
	bool multiplyMappedButNoXA1;
	bool multiplyMappedButNoXA2;
	bool expand_xa;
	void initAlignmentVectors();
	void clearAlignmentVectors();
	void addRead(std::vector<BamTools::BamAlignment*>* alignments, boost::unordered_set<BamHelper::alignment_coordinate_t>& coordinates);
	/** Verifies whether the number of reads complies with the given X0/X1 tags (if given). */
	void checkX0X1Tags(const std::vector<BamTools::BamAlignment*>* alignments, int count);

public:
	/** If paired is true, then the input BAM is expected to be paired end data. In this case, the methods
	 *  getAlignmentsFirst and getAlignmentsSecond must be used (instead of getAlignments) to get the
	 *  alignments of first/second read end. */
	GroupWiseBamReader(const std::string& filename, bool paired, bool skip_duplicates = true, bool expand_xa = false);
	virtual ~GroupWiseBamReader();

	/** Returns true if there is another read pair left in the input file. */
	virtual bool hasNext() const;

	/** Reads next read pair. Before calling this method, hasNext() should
	 *  be called. */
	virtual void advance();

	/** Returns the name of the currently processed read. */
	virtual const std::string& getReadName() const { return read_name; }

	virtual const std::vector<BamTools::BamAlignment*>& getAlignments() const;
	virtual const std::vector<BamTools::BamAlignment*>& getAlignmentsFirst() const;
	virtual const std::vector<BamTools::BamAlignment*>& getAlignmentsSecond() const;
	virtual std::auto_ptr<std::vector<BamTools::BamAlignment*> >  releaseAlignments();
	virtual std::auto_ptr<std::vector<BamTools::BamAlignment*> >  releaseAlignmentsFirst();
	virtual std::auto_ptr<std::vector<BamTools::BamAlignment*> >  releaseAlignmentsSecond();

	virtual bool isUnmapped() const;
	virtual bool isFirstUnmapped() const;
	virtual bool isSecondUnmapped() const;

	virtual bool hasInconsistentXTags() const;
	virtual bool hasMultipleMappingsButNoXA() const;
	virtual bool hasMultipleMappingsButNoXARead1() const;
	virtual bool hasMultipleMappingsButNoXARead2() const;
	virtual const std::string& getInconsistentXTagsMessage() const;

	virtual void enableProgressMessages(std::ostream& os, int frequency);

	virtual long long getSkippedDuplicates() const;
	virtual long long getNonPairedCount() const { return 0; }

	virtual const BamTools::RefVector& getReferenceData() const { return bam_reader.GetReferenceData(); }

	virtual BamTools::SamHeader getHeader() const { return bam_reader.GetHeader(); }
};

#endif /* GROUPWISEBAMREADER_H_ */
