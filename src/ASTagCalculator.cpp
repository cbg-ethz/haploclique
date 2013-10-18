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

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "FastaReader.h"

#include "ASTagCalculator.h"

using namespace std;
namespace io = boost::iostreams;

ASTagCalculator::ASTagCalculator(const string& reference_filename, int phred_offset, int bad_score_threshold) {
	this->bad_score_threshold = bad_score_threshold;
	this->phred_offset = phred_offset;

	ifstream reference_istream(reference_filename.c_str());
	bool input_zipped = reference_filename.substr(reference_filename.size()-3,3).compare(".gz") == 0;
	if (reference_istream.fail()) {
		throw std::runtime_error("Error opening file \"" + reference_filename + "\".");
	}
	io::filtering_istream in;
	if (input_zipped) {
		in.push(io::gzip_decompressor());
	}
	in.push(reference_istream);
	reference_sequences = FastaReader::parse(in).release();
	if (reference_sequences->size() == 0) {
		if (input_zipped) {
			throw std::runtime_error("Error: references sequences empty or not properly gzipped.");
		} else {
			throw std::runtime_error("Error: references sequences empty.");
		}
	}
}

ASTagCalculator::~ASTagCalculator() {
	typedef boost::unordered_map<string,NamedDnaSequence*>::const_iterator cit_t;
	for (cit_t it = reference_sequences->begin(); it!=reference_sequences->end(); ++it) {
		delete it->second;
	}
	delete reference_sequences;
}

bool ASTagCalculator::computeTag(BamTools::BamAlignment& read_aln, const BamTools::RefVector& bam_ref_data, stats_t* stats) const {
//	cerr << "Processing read " << read_aln.Name << endl;
	read_aln.RemoveTag("AS");
	if (!read_aln.IsMapped()) return true;
	if (stats != 0) stats->mapped_reads += 1;
	const string& ref_name = bam_ref_data[read_aln.RefID].RefName;
	typedef boost::unordered_map<string,NamedDnaSequence*>::const_iterator cit_t;
	cit_t it = reference_sequences->find(ref_name);
	if (it == reference_sequences->end()) {
		if (stats != 0) stats->reference_missing += 1;
		return false;
	}
	const NamedDnaSequence& ref_sequence = *it->second;
	size_t p_read = 0;
	size_t p_ref = read_aln.Position;
	uint32_t as_tag = 0;
	for (size_t i=0; i<read_aln.CigarData.size(); ++i) {
		const BamTools::CigarOp& cigar_op = read_aln.CigarData[i];
		switch (cigar_op.Type) {
		case 'M':
			for (size_t j=0; j<cigar_op.Length; ++j) {
				if ((p_ref >= ref_sequence.size()) || (read_aln.QueryBases[p_read] != ref_sequence[p_ref])) {
					as_tag += read_aln.Qualities[p_read]-33;
				}
				p_read += 1;
				p_ref += 1;
			}
			break;
		case 'I':
			p_read += cigar_op.Length;
			break;
		case 'D':
			p_ref += cigar_op.Length;
			break;
		case 'S':
			p_read += cigar_op.Length;
			break;
		case 'H':
		case 'P':
		case 'N':
			break;
		default:
			assert(false);
		}
	}
	if (p_ref > ref_sequence.size()) {
		if (stats != 0) stats->overhanging_alignments += 1;
	} else {
		if (as_tag > (uint32_t)bad_score_threshold) {
			if (stats != 0) stats->bad_alignments += 1;
		}
	}
	if (!read_aln.AddTag("AS", "I", as_tag)) {
		assert(false);
	}
	return true;
}

void ASTagCalculator::printWarnings(std::ostream& os, const stats_t& stats) const {
	if (stats.bad_alignments > 0) {
		os << "Warning: there were " << stats.bad_alignments <<" bad alignments (score>" << bad_score_threshold << ")." << endl;
	}
	if (stats.reference_missing > 0) {
		os << "Warning: for " << stats.reference_missing << " out of " << stats.mapped_reads << " alignments, no AS tag could be computed, because the reference sequence is unknown." << endl;
	}
	if (stats.overhanging_alignments > 0) {
		os << "Warning: " << stats.overhanging_alignments << " out of " << stats.mapped_reads << " alignments were skipped, because they stretch beyond the end of the reference." << endl;
	}
}
