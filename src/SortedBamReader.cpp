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

#include <sstream>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include "ShortDnaSequence.h"

#include "SortedBamReader.h"
#include "BamHelper.h"

using namespace std;

SortedBamReader::SortedBamReader(const std::string& filename, bool paired, int max_span, bool expand_xa) {
	if (!paired) throw std::runtime_error("SortedBamReader: not implemented for single reads, yet.");
	this->counter = 0;
	this->progress_messages_os = 0;
	this->progress_message_frequency = -1;
	this->skipped_duplicates = 0;
	this->not_paired_count = 0;
	this->paired = paired;
	this->previous_chromosome = -1;
	this->previous_position = -1;
	this->max_span = max_span;
	this->extra_records = 0;
	this->expand_xa = expand_xa;
	if (!bam_reader.Open(filename)) {
		ostringstream oss;
		oss << "Could not read BAM input from \"" << filename << "\".";
		throw std::runtime_error(oss.str());
	}
	processNextRead();
}

SortedBamReader::~SortedBamReader() {
	aln_by_readnames_t::const_iterator it = aln_by_readnames.begin();
	for (; it != aln_by_readnames.end(); ++it) {
		delete it->second;
	}
}

long long SortedBamReader::hashReadStart(const BamTools::BamAlignment& aln) {
	assert(sizeof(long long) >= 2*sizeof(int32_t));
	long long result = aln.Position;
	result |= ((long long)aln.RefID) << (8*sizeof(int32_t));
	return result;
}

void SortedBamReader::parseCigar(const string& cigar, vector<BamTools::CigarOp>* target) {
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer_t;
	boost::char_separator<char> cigar_separator("", "MIDNSHP=X");
	tokenizer_t tokenizer(cigar, cigar_separator);
	tokenizer_t::const_iterator token_it = tokenizer.begin();
	target->clear();
	while (token_it != tokenizer.end()) {
		int length = boost::lexical_cast<int>(*token_it);
		++token_it;
		assert(token_it != tokenizer.end());
		assert(token_it->size() == 1);
		char type = token_it->at(0);
		target->push_back(BamTools::CigarOp(type,length));
		++token_it;
	}
}

void SortedBamReader::read_bam_record(BamTools::BamAlignment* read_aln, read_record_t* record) {
	assert(paired);
	uint32_t x0 = 0;
	uint32_t x1 = 0;
	if (read_aln->GetTag("X0", x0) && read_aln->GetTag("X1", x1)) {
		string xa = "";
		if (!read_aln->GetTag("XA", xa) && (x0+x1>1)) {
			record->multiplyMappedButNoXA = true;
		}
	}
	vector<BamTools::BamAlignment*>& alignments = read_aln->IsFirstMate()?record->alignments1:record->alignments2;
	boost::unordered_set<BamHelper::alignment_coordinate_t>& coordinates = read_aln->IsFirstMate()?record->coordinates1:record->coordinates2;
	if (alignments.size() > 0) {
		// cerr << "Warning: Skipped extra alignment for \"" << read_aln->Name << "\"." << endl;
		return;
	}
	alignments.push_back(read_aln);
	if (!read_aln->IsMapped()) {
		if (read_aln->IsFirstMate()) record->unmapped1 = true;
		else record->unmapped2 = true;
		return;
	}
	coordinates.insert(BamHelper::alignment_coordinate_t(*read_aln));
	if (expand_xa) {
		BamHelper::expandXA(bam_reader, *read_aln, &alignments, &coordinates, &skipped_duplicates);
	}
}

void SortedBamReader::processNextRead() {
	while (true) {
		auto_ptr<BamTools::BamAlignment> read_aln(new BamTools::BamAlignment());
		if (!bam_reader.GetNextAlignment(*read_aln)){
			next_read = auto_ptr<read_record_t>(0);
			return;
		}
		if (paired) {
			if (read_aln->IsFirstMate() == read_aln->IsSecondMate()) {
				ostringstream oss;
				oss << "Invalid flags for read \"" << read_aln->Name << "\".";
				throw std::runtime_error(oss.str());
			}
		}
		counter += 1;
		if ((progress_messages_os != 0) && (counter % progress_message_frequency == 0)) {
			*progress_messages_os << "Having processed " << counter << " read alignments" << endl;
		}
		if (read_aln->IsFailedQC() || read_aln->IsDuplicate()) {
			// cerr << "Warning: FailedQC or duplicated flag for \"" << read_aln->Name << "\" --> skipped." << endl;
			continue;
		}
		if (previous_chromosome == read_aln->RefID) {
			if (previous_position > read_aln->Position) {
				ostringstream oss;
				oss << "BAM file not sorted. Offending read: \"" << read_aln->Name << "\".";
				throw std::runtime_error(oss.str());
			}
		} else {
			// new chromosome. clear old data.
			assert(aln_by_readnames.size() == aln_by_pos.size());
			not_paired_count += aln_by_readnames.size();
			for (aln_by_pos_t::const_iterator it = aln_by_pos.begin(); it != aln_by_pos.end(); ++it) {
				delete it->second;
			}
			aln_by_pos.clear();
			aln_by_readnames.clear();
		}
		previous_chromosome = read_aln->RefID;
		previous_position = read_aln->Position;
		// discard reads that could not be paired.
		aln_by_pos_t::iterator end_it = aln_by_pos.lower_bound(read_aln->Position - max_span);
		for (aln_by_pos_t::const_iterator it = aln_by_pos.begin(); it != end_it; ++it) {
			assert(it->second != 0);
			aln_by_readnames.erase(it->second->Name);
			delete it->second;
			not_paired_count += 1;
		}
		aln_by_pos.erase(aln_by_pos.begin(), end_it);
		// check whether current read could be found
		aln_by_readnames_t::iterator it = aln_by_readnames.find(read_aln->Name);
		if (it == aln_by_readnames.end()) {
			// if not found, store this alignment for later use
			BamTools::BamAlignment* a = read_aln.release();
			aln_by_readnames[a->Name] = a;
			aln_by_pos.insert(make_pair(a->Position,a));
		} else {
			// check whether a read of this name has been processed before
			if (finished_read_names.find(read_aln->Name) != finished_read_names.end()) {
				if (extra_records <= 1000) {
					cerr << "Warning: extra entry for read: \"" << read_aln->Name << "\" --> skipping" << endl;
				}
				if (extra_records == 1000) {
					cerr << "Limit of 1000 such warnings reached: suppressing further warnings" << endl;
				}
				extra_records += 1;
				next_read = auto_ptr<read_record_t>(0);
				continue;
			}
			// store the first million read names to check whether read names are unique
			if (finished_read_names.size() < 1000000) {
				finished_read_names.insert(read_aln->Name);
			}
			// delete this read from indices
			BamTools::BamAlignment* a = it->second;
			aln_by_readnames.erase(it);
			assert(a != 0);
			pair<aln_by_pos_t::iterator,aln_by_pos_t::iterator> range = aln_by_pos.equal_range(a->Position);
			for (aln_by_pos_t::iterator jt = range.first; jt != range.second; ++jt) {
				if (jt->second == a) {
					aln_by_pos.erase(jt);
					break;
				}
			}
			// create new entry
			next_read = auto_ptr<read_record_t>(new read_record_t());
			read_bam_record(read_aln.release(), next_read.get());
			read_bam_record(a, next_read.get());
			if ((next_read->alignments1.size() == 0) || (next_read->alignments2.size() == 0)) {
				next_read = auto_ptr<read_record_t>(0);
				not_paired_count += 1;
				// cerr << "Skipping strange read" << endl;
				continue;
			}
			return;
		}
	}
}

void SortedBamReader::advance() {
	assert(next_read.get() != 0);
	current_read = next_read;
	assert(next_read.get() == 0);
	processNextRead();
}

bool SortedBamReader::hasNext() const {
	return next_read.get() != 0;
}

const std::string& SortedBamReader::getReadName() const {
	assert(current_read.get() != 0);
	assert(current_read->alignments1.size() > 0);
	return current_read->alignments1[0]->Name;
}

const std::vector<BamTools::BamAlignment*>& SortedBamReader::getAlignments() const {
	assert(!paired);
	assert(current_read.get() != 0);
	return current_read->alignments1;
}

const vector<BamTools::BamAlignment*>& SortedBamReader::getAlignmentsFirst() const {
	assert(paired);
	assert(current_read.get() != 0);
	return current_read->alignments1;
}

const vector<BamTools::BamAlignment*>& SortedBamReader::getAlignmentsSecond() const {
	assert(paired);
	assert(current_read.get() != 0);
	return current_read->alignments2;
}

auto_ptr<vector<BamTools::BamAlignment*> >  SortedBamReader::releaseAlignments() {
	assert(!paired);
	assert(current_read.get() != 0);
	auto_ptr<vector<BamTools::BamAlignment*> > result(new vector<BamTools::BamAlignment*>(current_read->alignments1));
	current_read->alignments1.clear();
	return result;
}

auto_ptr<vector<BamTools::BamAlignment*> >  SortedBamReader::releaseAlignmentsFirst() {
	assert(paired);
	assert(current_read.get() != 0);
	auto_ptr<vector<BamTools::BamAlignment*> > result(new vector<BamTools::BamAlignment*>(current_read->alignments1));
	current_read->alignments1.clear();
	return result;
}

auto_ptr<vector<BamTools::BamAlignment*> > SortedBamReader::releaseAlignmentsSecond() {
	assert(paired);
	assert(current_read.get() != 0);
	auto_ptr<vector<BamTools::BamAlignment*> > result(new vector<BamTools::BamAlignment*>(current_read->alignments2));
	current_read->alignments2.clear();
	return result;
}

bool SortedBamReader::isUnmapped() const {
	assert(!paired);
	assert(current_read.get() != 0);
	return current_read->unmapped1;
}

bool SortedBamReader::isFirstUnmapped() const {
	assert(paired);
	assert(current_read.get() != 0);
	return current_read->unmapped1;
}

bool SortedBamReader::isSecondUnmapped() const {
	assert(paired);
	assert(current_read.get() != 0);
	return current_read->unmapped2;
}

bool SortedBamReader::hasMultipleMappingsButNoXA() const {
	assert(current_read.get() != 0);
	return current_read->multiplyMappedButNoXA;
}

void SortedBamReader::enableProgressMessages(std::ostream& os, int frequency) {
	assert(frequency > 0);
	this->progress_messages_os = &os;
	this->progress_message_frequency = frequency;
}

long long SortedBamReader::getSkippedDuplicates() const {
	return skipped_duplicates;
}

long long SortedBamReader::getNonPairedCount() const {
	return not_paired_count;
}
