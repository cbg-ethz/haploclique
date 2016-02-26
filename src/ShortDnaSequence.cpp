/* Copyright 2012-2014 Tobias Marschall and Armin Töpfer
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

#include <cassert>

#include "ShortDnaSequence.h"
#include <math.h>

using namespace std;

ShortDnaSequence::ShortDnaSequence() : forward((string*)0), forward_qualities((string*)0), backward((string*)0), backward_qualities((string*)0) {
	// cerr << "ShortDnaSequence::ShortDnaSequence()" << endl;
}

ShortDnaSequence::ShortDnaSequence(const std::string& dna, const std::string& qualities) {
	// cerr << "ShortDnaSequence::ShortDnaSequence(" << dna << ", " << qualities << ")" << endl;
	assert(dna.size() == qualities.size());
	this->forward = string_ptr_t(new string(dna.size(),' '));
	this->forward_qualities = string_ptr_t(new string(qualities));
	this->backward = string_ptr_t((string*)0);
	this->backward_qualities = string_ptr_t((string*)0);
	for (size_t i=0; i<dna.size(); ++i) {
		char c = toupper(dna[i]);
		switch (c) {
		case 'A':
		case 'C':
		case 'G':
		case 'T':
			forward->at(i) = c;
			break;
		default:
			forward->at(i) = 'N';
			break;
		}
	}
}

ShortDnaSequence::ShortDnaSequence(string_ptr_t& forward, string_ptr_t& forward_qualities, string_ptr_t& backward, string_ptr_t& backward_qualities) {
	this->forward = forward;
	this->forward_qualities = forward_qualities;
	this->backward = backward;
	this->backward_qualities = backward_qualities;
}

ShortDnaSequence::~ShortDnaSequence() {
	// cerr << "ShortDnaSequence::~ShortDnaSequence()" << endl;
}

ShortDnaSequence ShortDnaSequence::reverseComplement() const {
	// cerr << "ShortDnaSequence::reverseComplement()"  << endl;
	if (backward.get() == 0) {
		backward = string_ptr_t(new string(forward->size(),' '));
		backward_qualities = string_ptr_t(new string(forward->size(),' '));
		int b_pos = forward->size()-1;
		for (size_t i=0; i<forward->size(); ++i, --b_pos) {
			backward_qualities->at(b_pos) = forward_qualities->at(i);
			switch (forward->at(i)) {
			case 'A':
				backward->at(b_pos) = 'T';
				break;
			case 'C':
				backward->at(b_pos) = 'G';
				break;
			case 'G':
				backward->at(b_pos) = 'C';
				break;
			case 'T':
				backward->at(b_pos) = 'A';
				break;
			case 'N':
				backward->at(b_pos) = 'N';
				break;
			default:
				assert(false);
			}
		}
		// cerr << "Computed revcomp: " << *backward << endl;
	}
	return ShortDnaSequence(backward, backward_qualities, forward, forward_qualities);
}

size_t ShortDnaSequence::size() const {
	assert(forward.get() != 0);
	return forward->size();
}

char ShortDnaSequence::operator[](size_t pos) const {
	assert(forward.get() != 0);
    //if(pos >= forward->size()){
    //    int k = 0;
    //}
    assert(pos < forward->size());
	return forward->at(pos);
}

char ShortDnaSequence::qualityChar(size_t pos) const {
	assert(forward_qualities.get() != 0);
    //if(pos>=forward->size()){
    //    int k = 0;
    //}
	assert(pos < forward->size());
	return forward_qualities->at(pos);
}

double ShortDnaSequence::qualityCorrect(size_t pos) const {
    return 1 - pow(10, -(qualityChar(pos)-33)/10.0);
}

double ShortDnaSequence::qualityCorrectLog(size_t pos) const {
    return log10(pow(10, -(qualityChar(pos)-33)/10.0));
}

const std::string& ShortDnaSequence::toString() const {
	assert(forward.get() != 0);
	return *forward;
}

const std::string& ShortDnaSequence::qualityString() const {
	assert(forward_qualities.get() != 0);
	return *forward_qualities;
}

std::ostream& operator<<(std::ostream& os, const ShortDnaSequence& s) {
	assert(s.forward.get() != 0);
	return os << *s.forward;
}
