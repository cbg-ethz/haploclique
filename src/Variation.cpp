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
#include <stdexcept>
#include <sstream>

#include <boost/functional/hash.hpp>

#include "Variation.h"

using namespace std;

int Variation::getLengthDifference() const {
	switch (type) {
	case Variation::INSERTION:
		return -(int)coord2;
	case Variation::DELETION:
		return ((int)coord2) - ((int)coord1);
	default:
		assert(false);
	}
}

void Variation::canonify(const boost::unordered_map<std::string,NamedDnaSequence*>& references, bool rightmost) {
	typedef boost::unordered_map<string,NamedDnaSequence*>::const_iterator references_it_t;
	references_it_t it = references.find(chromosome);
	if (it == references.end()) {
		ostringstream oss;
		oss << "Cannot find reference " << chromosome;
		throw runtime_error(oss.str());
	}
	canonify(it->second, rightmost);
}

void Variation::canonify(const NamedDnaSequence* reference, bool rightmost) {
	if (type == NONE) return;
	if (type == DELETION) {
		assert(coord1 < coord2);
		assert(coord2 <= reference->size());
		if (rightmost) {
			while (coord2 < reference->size()) {
				if ((*reference)[coord1] != (*reference)[coord2]) break;
				coord1 += 1;
				coord2 += 1;
			}
		} else {
			while (coord1 > 0) {
				if ((*reference)[coord1-1] != (*reference)[coord2-1]) break;
				coord1 -= 1;
				coord2 -= 1;
			}
		}
	}
	if (type == INSERTION) {
		if (sequence.size() == 0) throw runtime_error("Variation::canonify is only applicable to insertions if sequence is known");
		assert(coord1 <= reference->size());
		assert(sequence.size() == coord2);
		if (rightmost) {
			// number of positions by which insertion can be shifted to the right
			int shift = 0;
			while (true) {
				if (coord1 + shift == reference->size()) break;
				if (sequence[shift % coord2] != (*reference)[coord1 + shift]) break;
				shift += 1;
			}
			if (shift == 0) return;
			coord1 += shift;
			// rotate insertion sequence by shift%insertion_length characters
			int j = shift % coord2;
			if (j == 0) return;
			sequence = sequence.substr(j, sequence.size()-j) + sequence.substr(0, j);
		} else {
			// number of positions by which insertion can be shifted to the left
			int shift = 0;
			while (true) {
				if (coord1 - shift == 0) break;
				if (sequence[coord2 - 1 - (shift % coord2)] != (*reference)[coord1 - shift - 1]) break;
				shift += 1;
			}
			if (shift == 0) return;
			coord1 -= shift;
			// rotate insertion sequence by shift%insertion_length characters
			int j = shift % coord2;
			if (j == 0) return;
			sequence = sequence.substr(sequence.size()-j, j) + sequence.substr(0, sequence.size()-j);
		}
	}
}

std::ostream& operator<<(std::ostream& os, const Variation& variation) {
	os << variation.getChromosome() << " " << variation.getCoordinate1() << " ";
	switch (variation.getType()) {
	case Variation::INSERTION:
		os << variation.getCoordinate2() << " INS";
		if (variation.getSequence().size() > 0) {
			os << " " << variation.getSequence();
		}
		break;
	case Variation::DELETION:
		os << (variation.getCoordinate2()-1) << " DEL";
		break;
	default:
		os << variation.getCoordinate2() << " NONE";
		break;
	}
	return os;
}

bool Variation::operator==(const Variation& v) const {
	return (chromosome.compare(v.chromosome) == 0) && (coord1 == v.coord1) && (coord2 == v.coord2) && (type == v.type) && (sequence == v.sequence);
}

std::size_t hash_value(Variation const& v) {
	size_t h = 0;
	boost::hash_combine(h, v.getChromosome());
	boost::hash_combine(h, v.getCoordinate1());
	boost::hash_combine(h, v.getCoordinate2());
	boost::hash_combine(h, v.getType());
	boost::hash_combine(h, v.getSequence());
	return h;
}
