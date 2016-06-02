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

#ifndef ALIGNMENTRECORD_H_
#define ALIGNMENTRECORD_H_

#include <string>
#include <deque>
#include <utility>
#include <set>
#include <tuple>
#include <unordered_map>

#include <api/BamAux.h>
#include <api/BamAlignment.h>

#include "Types.h"
#include "ShortDnaSequence.h"

class Clique;

/** Class that represents alignments of a read pair. */
class AlignmentRecord {
    public:
    /** Represents entry for a map containing the ref positions of an AlignmentRecord, the base,
     * the quality score and the position of the base in the read; */
    struct mapValue{
        int ref; //pos in ref
        char base;
        char qual; //phred score of base (QUALiy+33)
        double prob; //error probability for qual
        int pir; //position in read
        int read; //number of paired end read: 0 for first, 1 for second read
    };
private:
	std::string name;
	int phred_sum1;
	unsigned int start1;
	unsigned int end1;
	std::vector<BamTools::CigarOp> cigar1;
	std::vector<char> cigar1_unrolled;
	int length_incl_deletions1;
	int length_incl_longdeletions1;
	ShortDnaSequence sequence1;
    std::string al_sequence1;
	int phred_sum2;
	unsigned int start2;
	unsigned int end2;
	std::vector<BamTools::CigarOp> cigar2;
	std::vector<char> cigar2_unrolled;
    int length_incl_deletions2;
	int length_incl_longdeletions2;
	ShortDnaSequence sequence2;
    std::vector<mapValue> cov_pos;
    double probability;
	alignment_id_t id;
	bool single_end;
	std::set<int> readNames;
    std::vector<std::string>* readNameMap;

    /** merges sequences to superreads, i and j correspond to sequence/cigar 1 or 2*/
    void mergeAlignmentRecordsSingle(const AlignmentRecord& ar, int i, int j);
    void mergeAlignmentRecordsPaired(const AlignmentRecord& ar);
    void mergeAlignmentRecordsMixed(const AlignmentRecord& ar);
public:
    AlignmentRecord(const BamTools::BamAlignment& alignment, int id, std::vector<std::string>* readNameMap);
    AlignmentRecord(std::unique_ptr<std::vector<const AlignmentRecord*>>& alignments,unsigned int clique_id);
    /** merges overlapping paired end reads while reading in bam files*/
    void noOverlapMerge(std::string& dna, std::string& qualities, std::string& nucigar, int& c_pos, int& q_pos, int& ref_pos, int i) const;
    void noOverlapMerge(const BamTools::BamAlignment& alignment, std::string& dna, std::string& qualities, std::string& nucigar, std::vector<char>& cigar_temp_unrolled, int& c_pos, int& q_pos, int& ref_pos) const;
    void overlapMerge(const BamTools::BamAlignment& alignment, std::string& dna, std::string& qualities, std::string& nucigar, std::vector<char>& cigar_temp_unrolled, int& c_pos1, int& c_pos2, int& q_pos1, int& q_pos2, int& ref_pos) const;
    void overlapMerge(const AlignmentRecord& alignment, std::string& dna, std::string& qualities, std::string& nucigar, int& c_pos1, int& c_pos2, int& q_pos1, int& q_pos2, int& ref_pos, int i, int j) const;
    void getMergedDnaSequence(const BamTools::BamAlignment& alignment);
    /** combines to reads belonging to a paired end to one Alignment Record. They are merged if they overlap. */
    void pairWith(const BamTools::BamAlignment& alignment);

    unsigned int getRecordNr() const;
	int getPhredSum1() const;
	int getPhredSum2() const;
	
	/** Returns probability that alignment pair is correct based on alignment scores alone. */
	double getProbability() const;
    
	/** Returns start position of interval associated with this alignment record.
	 *  In case of a single end read, interval corresponds to the alignment;
	 *  in case of a paired end read, interval corresponds to the whole fragment, i.e. first alignment,
	 *  internal segment, and second alignment.
	 */
	unsigned int getIntervalStart() const;
	/** Returns end position of interval associated with this alignment record, see getIntervalStart(). */
	unsigned int getIntervalEnd() const;
	/** Returns length of intersection between two alignments with respect to the intervals given 
	  * by getIntervalStart() and getIntervalEnd(). */
	size_t intersectionLength(const AlignmentRecord& ap) const;

	/** Returns length of intersection between two alignments with respect to the intervals given 
	  * by getInsertStart() and getInsertEnd(). */
	size_t internalSegmentIntersectionLength(const AlignmentRecord& ap) const;

    /** Returns a map containing the reference positions which are covered by a read.  */
    std::vector<AlignmentRecord::mapValue> coveredPositions() const;

	unsigned int getEnd1() const;
	unsigned int getEnd2() const;
	std::string getName() const;
	unsigned int getStart1() const;
	unsigned int getStart2() const;
	const std::vector<BamTools::CigarOp>& getCigar1() const;
	const std::vector<BamTools::CigarOp>& getCigar2() const;
	const ShortDnaSequence& getSequence1() const;
	const ShortDnaSequence& getSequence2() const;
	int getReadGroup() const;
	unsigned int getInsertStart() const;
	unsigned int getInsertEnd() const;
	unsigned int getInsertLength() const;
	alignment_id_t getID() const;
	void setID(alignment_id_t id);
	bool isSingleEnd() const;
	bool isPairedEnd() const;
    const std::set<int>& getReadNamesSet() const{
        return readNames;
    }
	std::vector<std::string> getReadNames() const;
    const std::vector<char>& getCigar1Unrolled() const;
    const std::vector<char>& getCigar2Unrolled() const;
	int getLengthInclDeletions1() const;
	int getLengthInclDeletions2() const;
	int getLengthInclLongDeletions1() const;
	int getLengthInclLongDeletions2() const;
    const std::vector<mapValue>& getCovmap() const {
        return cov_pos;
    }

    unsigned int getReadCount() const { return readNames.size(); }

    friend double setProbabilities(std::deque<AlignmentRecord*>& reads);
    friend void printReads(std::ostream& output, std::deque<AlignmentRecord*>&, int doc_haplotypes);
};

#endif /* ALIGNMENTRECORD_H_ */
