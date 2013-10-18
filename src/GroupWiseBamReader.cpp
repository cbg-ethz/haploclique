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

#include "GroupWiseBamReader.h"

using namespace std;

GroupWiseBamReader::GroupWiseBamReader(const std::string& filename, bool paired, bool skip_duplicates, bool expand_xa) : read_name("") {
	this->counter = 0;
	this->progress_messages_os = 0;
	this->progress_message_frequency = -1;
	this->skip_duplicates = skip_duplicates;
	this->skipped_duplicates = 0;
	this->paired = paired;
	this->unmapped1 = true;
	this->unmapped2 = true;
	this->inconsistentXTags = false;
	this->inconsistentXTagsMessage = "";
	this->multiplyMappedButNoXA = false;
	this->multiplyMappedButNoXA1 = false;
	this->multiplyMappedButNoXA2 = false;
	this->expand_xa = expand_xa;
	if (!bam_reader.Open(filename)) {
		ostringstream oss;
		oss << "Could not read BAM input from \"" << filename << "\".";
		throw std::runtime_error(oss.str());
	}
	this->alignments1 = 0;
	this->alignments2 = 0;
	next_read_aln = new BamTools::BamAlignment();
	this->finished = !bam_reader.GetNextAlignment(*next_read_aln);
}

GroupWiseBamReader::~GroupWiseBamReader() {
	if (next_read_aln != 0) delete next_read_aln;
	clearAlignmentVectors();
}

void GroupWiseBamReader::initAlignmentVectors() {
	assert(alignments1 == 0);
	assert(alignments2 == 0);
	alignments1 = new vector<BamTools::BamAlignment*>();
	alignments2 = new vector<BamTools::BamAlignment*>();
}

void GroupWiseBamReader::clearAlignmentVectors() {
	if (alignments1 != 0) {
		for (size_t i = 0; i< alignments1->size(); ++i) {
			assert(alignments1->at(i) != 0);
			delete alignments1->at(i);
		}
		delete alignments1;
		alignments1 = 0;
	}
	if (alignments2 != 0) {
		for (size_t i = 0; i< alignments2->size(); ++i) {
			assert(alignments2->at(i) != 0);
			delete alignments2->at(i);
		}
		delete alignments2;
		alignments2 = 0;
	}
}

void GroupWiseBamReader::advance() {
	assert(!finished);
	read_name = next_read_aln->Name;
	clearAlignmentVectors();
	initAlignmentVectors();
	boost::unordered_set<BamHelper::alignment_coordinate_t> coordinates1;
	boost::unordered_set<BamHelper::alignment_coordinate_t> coordinates2;
	// start positions of the current alignments (stored to check for
	// duplicate alignments)
	boost::unordered_set<long long> alignments1_starts;
	boost::unordered_set<long long> alignments2_starts;
	// number of alignments (including those skipped because they are duplicates)
	int alignments1_count = 0;
	int alignments2_count = 0;
	int n = 0;
	int mapped_count1 = 0;
	int mapped_count2 = 0;
	this->inconsistentXTags = false;
	this->inconsistentXTagsMessage = "";
	this->multiplyMappedButNoXA = false;
	this->multiplyMappedButNoXA1 = false;
	this->multiplyMappedButNoXA2 = false;
	while (true) {
		if (n++ > 0) {
			if (!bam_reader.GetNextAlignment(*next_read_aln)){
				finished = true;
				break;
			}
			if (read_name.compare(next_read_aln->Name) != 0) {
				if (read_names.find(next_read_aln->Name) != read_names.end()) {
					ostringstream oss;
					oss << "Error: Reads not grouped properly. Offending read: \"" << next_read_aln->Name << "\"" << endl;
					throw std::runtime_error(oss.str());
				}
				// store the first million read names to check whether reads are grouped properly.
				if (read_names.size() < 1000000) {
					read_names.insert(next_read_aln->Name);
				}
				break;
			}
		}
		counter += 1;
		if ((progress_messages_os != 0) && (counter % progress_message_frequency == 0)) {
			*progress_messages_os << "Having processed " << counter << " read alignments" << endl;
		}
		if (paired) {
			assert(next_read_aln->IsFirstMate() != next_read_aln->IsSecondMate());
		}
		uint32_t x0 = 0;
		uint32_t x1 = 0;
		if (next_read_aln->GetTag("X0", x0) && next_read_aln->GetTag("X1", x1)) {
			string xa = "";
			if (!next_read_aln->GetTag("XA", xa) && (x0+x1>1)) {
				multiplyMappedButNoXA = true;
				if (!paired || next_read_aln->IsFirstMate()) {
					multiplyMappedButNoXA1 = true;
				} else {
					multiplyMappedButNoXA2 = true;
				}
			}
		}
		if (!paired || next_read_aln->IsFirstMate()) {
			if (next_read_aln->IsMapped()) mapped_count1 += 1;
			addRead(alignments1, coordinates1);
			alignments1_count += 1;
		} else {
			if (next_read_aln->IsMapped()) mapped_count2 += 1;
			addRead(alignments2, coordinates2);
			alignments2_count += 1;
		}
	}
	// Mark reads as unmapped if there are "alignments" that have 
	// the unmapped flag set.
	// TODO: Output some warning if 0 < mapped_count1 < alignments1_count
	this->unmapped1 = mapped_count1 < alignments1_count;
	this->unmapped2 = mapped_count2 < alignments2_count;
	checkX0X1Tags(alignments1,alignments1_count);
	if (paired) {
		checkX0X1Tags(alignments2,alignments2_count);
	}
	if (expand_xa) {
		if (alignments1->size() == 1) {
			BamHelper::expandXA(bam_reader, *(alignments1->at(0)), alignments1, &coordinates1, &skipped_duplicates);
		}
		if (paired && (alignments2->size() == 1)) {
			BamHelper::expandXA(bam_reader, *(alignments2->at(0)), alignments2, &coordinates2, &skipped_duplicates);
		}
	}
}

bool GroupWiseBamReader::hasNext() const {
	return !finished;
}

const std::vector<BamTools::BamAlignment*>& GroupWiseBamReader::getAlignments() const {
	assert(!paired);
	assert(alignments1 != 0);
	return *alignments1;
}

const std::vector<BamTools::BamAlignment*>& GroupWiseBamReader::getAlignmentsFirst() const {
	assert(paired);
	assert(alignments1 != 0);
	return *alignments1;
}

const std::vector<BamTools::BamAlignment*>& GroupWiseBamReader::getAlignmentsSecond() const {
	assert(paired);
	assert(alignments2 != 0);
	return *alignments2;
}

std::auto_ptr<std::vector<BamTools::BamAlignment*> >  GroupWiseBamReader::releaseAlignments() {
	assert(!paired);
	assert(alignments1 != 0);
	auto_ptr<vector<BamTools::BamAlignment*> > result(alignments1);
	alignments1 = 0;
	return result;
}

std::auto_ptr<std::vector<BamTools::BamAlignment*> >  GroupWiseBamReader::releaseAlignmentsFirst() {
	assert(paired);
	assert(alignments1 != 0);
	auto_ptr<vector<BamTools::BamAlignment*> > result(alignments1);
	alignments1 = 0;
	return result;
}

std::auto_ptr<std::vector<BamTools::BamAlignment*> >  GroupWiseBamReader::releaseAlignmentsSecond() {
	assert(paired);
	assert(alignments2 != 0);
	auto_ptr<vector<BamTools::BamAlignment*> > result(alignments2);
	alignments2 = 0;
	return result;
}

bool GroupWiseBamReader::isUnmapped() const {
	assert(!paired);
	return unmapped1;
}

bool GroupWiseBamReader::isFirstUnmapped() const {
	assert(paired);
	return unmapped1;
}

bool GroupWiseBamReader::isSecondUnmapped() const {
	assert(paired);
	return unmapped2;
}

bool GroupWiseBamReader::hasInconsistentXTags() const {
	return inconsistentXTags;
}

bool GroupWiseBamReader::hasMultipleMappingsButNoXA() const {
	return multiplyMappedButNoXA;
}

bool GroupWiseBamReader::hasMultipleMappingsButNoXARead1() const {
	return multiplyMappedButNoXA1;
}

bool GroupWiseBamReader::hasMultipleMappingsButNoXARead2() const {
	return multiplyMappedButNoXA2;
}

const std::string& GroupWiseBamReader::getInconsistentXTagsMessage() const {
	assert(inconsistentXTags);
	return inconsistentXTagsMessage;
}


void GroupWiseBamReader::enableProgressMessages(std::ostream& os, int frequency) {
	assert(frequency > 0);
	this->progress_messages_os = &os;
	this->progress_message_frequency = frequency;
}

long long GroupWiseBamReader::getSkippedDuplicates() const {
	return skipped_duplicates;
}

void GroupWiseBamReader::addRead(vector<BamTools::BamAlignment*>* alignments, boost::unordered_set<BamHelper::alignment_coordinate_t>& coordinates) {
	assert(alignments != 0);
	BamHelper::alignment_coordinate_t c(*next_read_aln);
	// TODO: Do a smarter check for duplicates, not just based on start position!
	if (skip_duplicates && (coordinates.find(c) != coordinates.end())) {
		skipped_duplicates += 1;
		delete next_read_aln;
		// cerr << "Skipping "  << read_aln.Name << endl;
	} else {
		coordinates.insert(c);
		alignments->push_back(next_read_aln);
	}
	next_read_aln = new BamTools::BamAlignment();
}

/** Verifies whether the number of reads complies with the given X0/X1 tags (if given). */
void GroupWiseBamReader::checkX0X1Tags(const vector<BamTools::BamAlignment*>* alignments, int count) {
	vector<BamTools::BamAlignment*>::const_iterator it = alignments->begin();
	for (;it!=alignments->end();++it) {
		uint32_t x0 = 0;
		uint32_t x1 = 0;
		if ((*it)->GetTag("X0", x0) && (*it)->GetTag("X1", x1)) {
			if (count != (int)(x0+x1)) {
				ostringstream oss;
				oss << "Alignment count mismatch for read \"" << (*it)->Name << "\": ";
				oss << "X0=" << x0 << ", X1=" << x1 << ", but found " << count << " alignment(s).";
				inconsistentXTagsMessage = oss.str();
				inconsistentXTags = true;
				return;
			}
		}
	}
}
