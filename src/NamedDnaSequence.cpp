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

#include <cassert>
#include <sstream>
#include "NamedDnaSequence.h"

using namespace std;

NamedDnaSequence::NamedDnaSequence(const string& name) : case_sensitive(false), name(name), length(0) {
}

NamedDnaSequence::NamedDnaSequence(const string& name, const string& sequence) : case_sensitive(false), name(name), length(0) {
	append(sequence);
}

NamedDnaSequence::NamedDnaSequence(const string& name, bool case_sensitive) : case_sensitive(case_sensitive), name(name), length(0) {
}

NamedDnaSequence::NamedDnaSequence(const string& name, const string& sequence, bool case_sensitive) : case_sensitive(case_sensitive), name(name), length(0) {
	append(sequence);
}

char NamedDnaSequence::decode(unsigned char c) const {
	switch (c) {
	case 0: return 'A';
	case 1: return 'C';
	case 2: return 'G';
	case 3: return 'T';
	case 4: return 'a';
	case 5: return 'c';
	case 6: return 'g';
	case 7: return 't';
	case 8: return 'n';
	default: return 'N';
	}
}

unsigned char NamedDnaSequence::encode(char c) const {
	if (case_sensitive) {
		switch (c) {
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
		case 'a': return 4;
		case 'c': return 5;
		case 'g': return 6;
		case 't': return 7;
		case 'n': return 8;
		default: return 9;
		}
	} else {
		switch (c) {
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
		case 'a': return 0;
		case 'c': return 1;
		case 'g': return 2;
		case 't': return 3;
		case 'n': return 9;
		default: return 9;
		}
	}
}

void NamedDnaSequence::append(const string& s) {
	string::const_iterator it = s.begin();
	for (;it!=s.end(); ++it) {
		if (length % 2 == 0) {
			sequence.push_back(encode(*it));
		} else {
			sequence[length/2] |= encode(*it)<<4;
		}
		length += 1;
	}
}

char NamedDnaSequence::operator[](size_t pos) const {
	assert(pos < length);
	if (pos % 2 == 0) {
		return decode(sequence[pos>>1] & 0x0f);
	} else {
		return decode(sequence[pos>>1] >>4);
	}
}

std::string NamedDnaSequence::substr(size_t pos, size_t length) const {
	assert(pos+length <= size());
	ostringstream oss;
	for (size_t i=0; i<length; ++i) {
		oss << (*this)[pos+i];
	}
	return oss.str();
}
