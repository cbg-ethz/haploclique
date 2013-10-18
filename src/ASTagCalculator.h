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

#ifndef ASTAGCALCULATOR_H_
#define ASTAGCALCULATOR_H_

#include <iostream>
#include <boost/unordered_map.hpp>
#include <bamtools/api/BamReader.h>

#include "NamedDnaSequence.h"

class ASTagCalculator {
private:
	boost::unordered_map<std::string,NamedDnaSequence*>* reference_sequences;
	int bad_score_threshold;
	int phred_offset;
public:
	typedef struct stats_t {
		long long mapped_reads;
		long long reference_missing;
		long long bad_alignments;
		long long overhanging_alignments;
		void operator+=(const stats_t& s) {
			mapped_reads += s.mapped_reads;
			reference_missing += s.reference_missing;
			bad_alignments += s.bad_alignments;
			overhanging_alignments += s.overhanging_alignments;
		}
		stats_t() : mapped_reads(0), reference_missing(0), bad_alignments(0), overhanging_alignments(0) {}
	} stats_t;
	
	ASTagCalculator(const std::string& reference_filename, int phred_offset, int bad_score_threshold);
	virtual ~ASTagCalculator();

	/** Returns false if no tag could be computed because the reference is unknown. */
	bool computeTag(BamTools::BamAlignment& read_aln, const BamTools::RefVector& bam_ref_data, stats_t* stats) const;
	void printWarnings(std::ostream& os, const stats_t& stats) const;
};

#endif /* ASTAGCALCULATOR_H_ */
