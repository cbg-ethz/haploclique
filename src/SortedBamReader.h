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

#ifndef SORTEDBAMREADER_H_
#define SORTEDBAMREADER_H_

#include <iostream>
#include <vector>
#include <map>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

#include <bamtools/api/BamReader.h>

#include "BamReader.h"
#include "BamHelper.h"

/** Reads a sorted BAM file and returns groups of alignments belonging to
 *  the same read. If present, interprets XA tags as written by BWA.
 */
class SortedBamReader : public BamReader {
private:
	BamTools::BamReader bam_reader;
	// names of processed reads to verify that they do not appear again
	boost::unordered_set<std::string> finished_read_names;
	// stored read records waiting to be paired 
	typedef boost::unordered_map<std::string,BamTools::BamAlignment*> aln_by_readnames_t;
	aln_by_readnames_t aln_by_readnames;
	// the same set of alignments but ordered by position
	typedef std::multimap<int32_t,BamTools::BamAlignment*> aln_by_pos_t;
	aln_by_pos_t aln_by_pos;
	// count all processed alignments
	long long counter;
	std::ostream* progress_messages_os;
	int progress_message_frequency;
	long long extra_records;
	
	typedef struct alignment_coordinate_t {
		int32_t ref_id;
		int32_t start;
		int32_t end;
		bool operator==(const alignment_coordinate_t& c) const { return (ref_id == c.ref_id) && (start == c.start) && (end == c.end); }
		alignment_coordinate_t() : ref_id(-1), start(0), end(0) {}
		alignment_coordinate_t(const BamTools::BamAlignment& aln) : ref_id(aln.RefID), start(aln.Position), end(aln.GetEndPosition()) {}
	} alignment_coordinates;
	typedef struct read_record_t {
		bool unmapped1;
		bool unmapped2;
		bool multiplyMappedButNoXA;
		std::vector<BamTools::BamAlignment*> alignments1;
		std::vector<BamTools::BamAlignment*> alignments2;
		boost::unordered_set<BamHelper::alignment_coordinate_t> coordinates1;
		boost::unordered_set<BamHelper::alignment_coordinate_t> coordinates2;
		read_record_t() : unmapped1(false), unmapped2(false), multiplyMappedButNoXA(false) {}
		~read_record_t() {
			for (size_t i=0; i<alignments1.size(); ++i) delete alignments1[i];
			for (size_t i=0; i<alignments2.size(); ++i) delete alignments2[i];
		}
	} read_record_t;
	
	std::auto_ptr<read_record_t> next_read;
	std::auto_ptr<read_record_t> current_read;
	
	long long skipped_duplicates;
	long long not_paired_count;
	bool paired;
	int previous_chromosome;
	int32_t previous_position;
	int max_span;
	bool expand_xa;

	long long hashReadStart(const BamTools::BamAlignment& aln);
	void addRead(std::vector<BamTools::BamAlignment>& alignments, boost::unordered_set<long long>& alignments_starts, const BamTools::BamAlignment& read_aln);
	/** Verifies whether the number of reads complies with the given X0/X1 tags (if given). */
// 	void checkX0X1Tags(const std::vector<BamTools::BamAlignment>& alignments, int count);
	void read_bam_record(BamTools::BamAlignment* read_aln, read_record_t* next_read);
	void processNextRead();
	void parseCigar(const std::string& cigar, std::vector<BamTools::CigarOp>* target);

public:
	/** If paired is true, then the input BAM is expected to be paired end data. In this case, the methods
	 *  getAlignmentsFirst and getAlignmentsSecond must be used (instead of getAlignments) to get the
	 *  alignments of first/second read end. */
	SortedBamReader(const std::string& filename, bool paired, int max_span, bool expand_xa);
	virtual ~SortedBamReader();

	/** Returns true if there is another read pair left in the input file. */
	virtual bool hasNext() const;

	/** Reads next read pair. Before calling this method, hasNext() should
	 *  be called. */
	virtual void advance();

	/** Returns the name of the currently processed read. */
	virtual const std::string& getReadName() const;

	/** Returns a reference to alignments current read. May only be called when
	 *  paired = false was given at construction time.
	 *  The returned reference is valid until the next call to advance(). */
	virtual const std::vector<BamTools::BamAlignment*>& getAlignments() const;

	/** Returns a reference to alignments of first read in current read pair.
	 *  The returned reference is valid until the next call to advance(). */
	virtual const std::vector<BamTools::BamAlignment*>& getAlignmentsFirst() const;

	/** Returns a reference to alignments of second read in current read pair.
	 *  The returned reference is valid until the next call to advance(). */
	virtual const std::vector<BamTools::BamAlignment*>& getAlignmentsSecond() const;

	/** Same as getAlignments, but transfers pointer ownerships. */
	virtual std::auto_ptr<std::vector<BamTools::BamAlignment*> >  releaseAlignments();
	/** Same as getAlignmentsFirst, but transfers pointer ownerships. */
	virtual std::auto_ptr<std::vector<BamTools::BamAlignment*> >  releaseAlignmentsFirst();
	/** Same as getAlignmentsSecond, but transfers pointer ownerships. */
	virtual std::auto_ptr<std::vector<BamTools::BamAlignment*> >  releaseAlignmentsSecond();

	virtual bool isUnmapped() const;
	virtual bool isFirstUnmapped() const;
	virtual bool isSecondUnmapped() const;

	virtual bool hasMultipleMappingsButNoXA() const;

	virtual void enableProgressMessages(std::ostream& os, int frequency);

	virtual long long getSkippedDuplicates() const;
	virtual long long getNonPairedCount() const;

	virtual const BamTools::RefVector& getReferenceData() const { return bam_reader.GetReferenceData(); }

	virtual BamTools::SamHeader getHeader() const { return bam_reader.GetHeader(); }
};

#endif /* SORTEDBAMREADER_H_ */
