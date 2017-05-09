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

#include <vector>
#include <iomanip>

#include "Clique.h"
#include "CliqueFinder.h"

using namespace std;

Clique::Clique(CliqueFinder& parent, size_t index, size_t capacity) : parent(parent) {
	assert(index < capacity);
	alignment_set = new alignment_set_t(capacity);
	alignment_set->set(index, true);
	init();
}

Clique::Clique(CliqueFinder& parent, std::unique_ptr<alignment_set_t>& alignments) : parent(parent) {
	alignment_set = alignments.release();
	init();
}

Clique::~Clique() {
	delete alignment_set;
}

void Clique::init() {
	int n = 0;
	alignment_count = 0;
	for (size_t i=0; i<alignment_set->size(); ++i) {
		if (!alignment_set->test(i)) continue;
		const AlignmentRecord& ap = parent.getAlignmentByIndex(i);
		if (n++==0) {
			leftmost_segment_start = ap.getIntervalStart();
			rightmost_segment_end = ap.getIntervalEnd();
		} else {
			leftmost_segment_start = min(leftmost_segment_start, (size_t)ap.getIntervalStart());
			rightmost_segment_end = max(rightmost_segment_end, (size_t)ap.getIntervalEnd());
		}
        const std::set<int> & readNames = ap.getReadNamesSet();
        cliqueReadNames.insert( readNames.begin(),readNames.end());
		alignment_count += 1;
	}
}

void Clique::computeIntervalIntersection(unsigned int* interval_start, unsigned int* interval_end) {
	int n = 0;
	for (size_t i=0; i<alignment_set->size(); ++i) {
		if (!alignment_set->test(i)) continue;
		const AlignmentRecord& ap = parent.getAlignmentByIndex(i);
		if (n++==0) {
			*interval_start = ap.getIntervalStart();
			*interval_end = ap.getIntervalEnd();
		} else {
			*interval_start = max(*interval_start,ap.getIntervalStart());
			*interval_end = min(*interval_end,ap.getIntervalEnd());
		}
	}
}

void Clique::translate(size_t* translation_table, size_t table_size, size_t new_capacity) {
	alignment_set_t* new_alignment_set = new alignment_set_t(new_capacity);
	for (size_t i=0; i<table_size; ++i) {
		if (alignment_set->test(translation_table[i])) {
			new_alignment_set->set(i,true);
		}
	}
	delete alignment_set;
	alignment_set = new_alignment_set;
}

unique_ptr<alignment_set_t> Clique::intersect(const alignment_set_t& set) const {
	unique_ptr<alignment_set_t> intersection(new alignment_set_t(set));
	(*intersection) &= (*alignment_set);
	return intersection;
}

bool Clique::contains(const Clique& c) {
	return c.alignment_set->is_subset_of(*alignment_set);
}

void Clique::add(size_t index) {
	if (index >= alignment_set->size()) {
		cout << "Clique::add failure from clique: " << *this <<endl;
		cout << "Capacity: " << alignment_set->size() << endl;
	}
	alignment_set->set(index, true);
	const AlignmentRecord& ap = parent.getAlignmentByIndex(index);
	leftmost_segment_start = min(leftmost_segment_start,(size_t)ap.getIntervalStart());
	rightmost_segment_end = max(rightmost_segment_end,(size_t)ap.getIntervalEnd());
    const std::set<int> & readNames = ap.getReadNamesSet();
    cliqueReadNames.insert( readNames.begin(),readNames.end());
    alignment_count += 1;
}

unique_ptr<vector<const AlignmentRecord*> > Clique::getAllAlignments() const {
	unique_ptr<vector<const AlignmentRecord*> > result(new vector<const AlignmentRecord*>());
	for (size_t i=alignment_set->find_first(); i!=alignment_set_t::npos; i=alignment_set->find_next(i)) {
		result->push_back(&parent.getAlignmentByIndex(i));
	}
	return result;
}

void Clique::printSet(std::ostream& os) {
	os << "[";
	int n = 0;
	for (size_t i=alignment_set->find_first(); i!=alignment_set_t::npos; i=alignment_set->find_next(i)) {
		if (n++>0) os << ",";
		os << setw(3) << i;
	}
	os << "]";
}

ostream& operator<<(ostream& ostream, const Clique& clique) {
	vector<alignment_id_t> ids;
	for (size_t i=clique.alignment_set->find_first(); i!=alignment_set_t::npos; i=clique.alignment_set->find_next(i)) {
		ids.push_back(clique.parent.getAlignmentByIndex(i).getID());
	}
	ostream << "[" << clique.leftmost_segment_start << ":" << clique.rightmost_segment_end << "]";
	ostream << "(";
	int n = 0;
	vector<alignment_id_t>::const_iterator it = ids.begin();
	for (; it!=ids.end(); ++it) {
		if (n++>0) ostream << ",";
		ostream << *it;
	}
	ostream << ")";
	return ostream;
}
