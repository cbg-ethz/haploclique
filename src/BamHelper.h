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


#ifndef BAMHELPER_H
#define BAMHELPER_H

#include <ostream>
#include <vector>
#include <memory>

#include <boost/unordered_set.hpp>

#include <api/BamAux.h>
#include <api/BamAlignment.h>
#include <api/BamWriter.h>
#include <api/BamReader.h>

#include "IndelLengthDistribution.h"
#include "PositionSet.h"
#include "VariationIndex.h"
#include "ShortDnaSequence.h"

/** Collection of (static) helper functions to handle BAM files, especially
 *  BAM files created by laser-core with YA tags for alternative cigar
 *  strings.
 */
class BamHelper {
public:
	/** A pairing of two alignments to form an alignment pair. */
	typedef struct pairing_t {
		// index of alignment of first/second read end
		short aln_idx1;
		short aln_idx2;
		// index of cigar string (aka "subalignment") of first/second alignment 
		short cigar_idx1;
		short cigar_idx2;
		pairing_t() : aln_idx1(0), aln_idx2(0), cigar_idx1(0), cigar_idx2(0) {}
		pairing_t(short aln_idx1, short aln_idx2) : aln_idx1(aln_idx1), aln_idx2(aln_idx2), cigar_idx1(0), cigar_idx2(0) {}
		pairing_t(short aln_idx1, short aln_idx2, short cigar_idx1, short cigar_idx2) : aln_idx1(aln_idx1), aln_idx2(aln_idx2), cigar_idx1(cigar_idx1), cigar_idx2(cigar_idx2) {}
	} pairing_t;

	/** Details on the alignment such as cigar string and phred scores. */
	typedef struct subalignment_t {
		std::vector<BamTools::CigarOp> cigar;
		// phred score of this alignment, including indel costs
		int phred_score;
		// phred score only counting mismatches, but no indels
		int mismatch_phred_score;
		subalignment_t() : cigar(), phred_score(-1) {}
		subalignment_t(const std::vector<BamTools::CigarOp>& cigar, int phred_score, int mismatch_phred_score) : cigar(cigar.begin(), cigar.end()), phred_score(phred_score), mismatch_phred_score(mismatch_phred_score) {}
		int indel_count() const {
			int result = 0;
			for (size_t i=0; i<cigar.size(); ++i){
				if ((cigar[i].Type == 'I') || (cigar[i].Type == 'D')) result += 1;
			}
			return result;
		}
	} subalignment_t;

	/** Data type representing a set of alignments for a paired end read. */
	typedef struct read_t {
		// rerefences to the original BAM records
		const std::vector<BamTools::BamAlignment*>& alignments1;
		const std::vector<BamTools::BamAlignment*>& alignments2;
		// alternative alignments extracted from the BAM records
		std::vector<std::vector<subalignment_t> > subalignments1;
		std::vector<std::vector<subalignment_t> > subalignments2;
		read_t(const std::vector<BamTools::BamAlignment*>& alignments1, const std::vector<BamTools::BamAlignment*>& alignments2) : alignments1(alignments1), alignments2(alignments2) {}
	} read_t;

	typedef struct alignment_coordinate_t {
		int32_t ref_id;
		int32_t start;
		int32_t end;
		bool operator==(const alignment_coordinate_t& c) const { return (ref_id == c.ref_id) && (start == c.start) && (end == c.end); }
		alignment_coordinate_t() : ref_id(-1), start(0), end(0) {}
		alignment_coordinate_t(const BamTools::BamAlignment& aln) : ref_id(aln.RefID), start(aln.Position), end(aln.GetEndPosition()) {}
	} alignment_coordinates;

	typedef struct aln_pair_t {
		BamTools::BamAlignment* first;
		BamTools::BamAlignment* second;
		aln_pair_t() : first(0), second(0) {}
	} aln_pair_t;

private:
	typedef struct {
		bool operator()(const subalignment_t& s1, const subalignment_t& s2) {
			// 1. Criterion: PHRED score
			if (s1.phred_score != s2.phred_score) return s1.phred_score < s2.phred_score;
			// 2. Criterion number of indels
			int indel_count1 = s1.indel_count();
			int indel_count2 = s1.indel_count();
			if (indel_count1 != indel_count2) return indel_count1 < indel_count2;
			// 3. Prefer leftmost indels
			int i1 = 0;
			int i2 = 0;
			int pos1 = 0;
			int pos2 = 0;
			for (int j=0; j<indel_count1; ++j) {
				while ((s1.cigar[i1].Type != 'D') && (s1.cigar[i1].Type != 'I')) {
					if ((s1.cigar[i1].Type == 'D') || (s1.cigar[i1].Type == 'M') || (s1.cigar[i1].Type == 'X') || (s1.cigar[i1].Type == '=')) {
						pos1 += s1.cigar[i1].Length;
					}
					i1 += 1;
				}
				while ((s2.cigar[i2].Type != 'D') && (s2.cigar[i2].Type != 'I')) {
					if ((s2.cigar[i2].Type == 'D') || (s2.cigar[i2].Type == 'M') || (s2.cigar[i2].Type == 'X') || (s2.cigar[i2].Type == '=')) {
						pos2 += s2.cigar[i2].Length;
					}
					i2 += 1;
				}
				if (pos1 != pos2) return pos1 < pos2;
			}
			// 4. Prefer shorter CIGAR strings
			if (s1.cigar.size() != s2.cigar.size()) return s1.cigar.size() < s2.cigar.size();
			// 5. Last ressort tie breaker: lexicographical order
			for (size_t i=0; i<s1.cigar.size(); ++i) {
				if (s1.cigar[i].Length != s2.cigar[i].Length) return s1.cigar[i].Length < s2.cigar[i].Length;
				if (s1.cigar[i].Type != s2.cigar[i].Type) return s1.cigar[i].Type < s2.cigar[i].Type;
			}
			return false;
		}
	} subalignment_comparator_t;
	
	static void write_alignment_record(BamTools::BamWriter& bam_writer, const std::vector<BamTools::BamAlignment*>& alignments, const std::vector<std::vector<subalignment_t> >& subalignments, int aln_idx, int cigar_idx, const BamTools::BamAlignment* mate_aln, double probability, bool is_primary, bool retain_alternative_cigars, bool reduce_cigar,  bool readgroups_from_names);
	static std::unique_ptr<std::vector<double> > compute_alignment_distribution(const std::vector<std::vector<subalignment_t> >& subalignments, int* best);
public:
	/** Parse a CIGAR string and append result to the given target vector. */
	static void parseCigar(const std::string& cigar, std::vector<BamTools::CigarOp>* target);

	/** Reads CIGAR string, AS, YM, and YA tags and writes all alternative alignments present
	 *  in the given alignment record to the target vector. The "main" alignment given by CIGAR 
	 *  string is written first.
	 */
	static void getSubalignments(const BamTools::BamAlignment& aln, std::vector<subalignment_t>* target);
	
	/** Extracts all subalignments and creates a record for a given read pair. */
	static std::unique_ptr<read_t> createReadRecord(const std::vector<BamTools::BamAlignment*>& alignments1, const std::vector<BamTools::BamAlignment*>& alignments2);

	/** Recompute PHRED score based on insertion/deletion length distributions, CIGAR string, and known PHRED scores of mismatches. */
	static void recalibratePhredScore(subalignment_t* subalignment, const BamTools::BamAlignment& aln, const IndelLengthDistribution& insertion_costs, const IndelLengthDistribution& deletion_costs, PositionSet* snp_set = 0, int phred_offset = 33, VariationIndex* variation_set = 0);

	/** Recalibrates a list of subalignments and sorts them by phred score. */
	static void recalibratePhredScores(std::vector<subalignment_t>* subalignments, const BamTools::BamAlignment& aln, const IndelLengthDistribution& insertion_costs, const IndelLengthDistribution& deletion_costs, PositionSet* snp_set = 0, int phred_offset = 33, VariationIndex* variation_set = 0);

	/** Computes the best pairing of alignments based on a given internal segment size distribution. 
	 *  @param variation_set If given, known indels the lie in the internal segment will be taken into account
	 *                       when computing probability.
	 */
	static pairing_t computeBestPairing(const read_t& read, const HistogramBasedDistribution& internal_segment_size_dist, VariationIndex* variation_set = 0);
	
	/** Writes alignments for a read pair to a BAM file. */
	static void writeAlignments(BamTools::BamWriter& bam_writer, const read_t& read, const pairing_t& best_pairing, bool retain_suboptimal = true, bool retain_alternative_cigars = true, bool reduce_cigar = false, bool readgroups_from_names = false);

	/** Turns all 'X' and '=' events in a CIGAR string to 'M' events. */
	static void cigarReduceMismatches(std::vector<BamTools::CigarOp>* cigar);
	
	/** Expands the XA tag in a given alignment (if present) and appends corresponding records to "target". */
	static void expandXA(const BamTools::BamReader& bam_reader, const BamTools::BamAlignment& aln, std::vector<BamTools::BamAlignment*>* target, boost::unordered_set<BamHelper::alignment_coordinate_t>* coordinates_set = 0, long long* skipped_counter = 0);

	/** Returns variations, i.e. insertions and deletions, present in the given alignments. */
	static std::unique_ptr<std::vector<Variation> > variationsFromAlignment(const BamTools::RefVector& bam_ref_data, const BamTools::BamAlignment& aln);

	/** Reads all alignment pairs that lie in the given region and appends them to
	 *  the given vector "target". Ownership is transferred, i.e. the alignments 
	 *  must be deleted by the caller. Only full pairs, i.e. pair where both read ends are known, are reported.
	 *  Also, only primary alignments are returned. */
	static void readRegion(BamTools::BamReader& bam_reader, int chromosome_id, int start, int end, std::vector<aln_pair_t>* target);

	friend std::size_t hash_value(const BamHelper::alignment_coordinate_t& c);
};

std::ostream& operator<<(std::ostream& os, const BamHelper::alignment_coordinate_t& c);
std::ostream& operator<<(std::ostream& os, const std::vector<BamTools::CigarOp>& cigar);

#endif // BAMHELPER_H
