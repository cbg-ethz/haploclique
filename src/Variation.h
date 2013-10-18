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

#ifndef VARIATION_H_
#define VARIATION_H_

#include <iostream>
#include <boost/unordered_map.hpp>

#include "NamedDnaSequence.h"

/** Represents a structural variation. */
class Variation {
public:
	typedef enum { NONE, INSERTION, DELETION } variation_type_t;
private:
	std::string chromosome;
	size_t coord1;
	size_t coord2;
	double p_value;
	std::string sequence;
	variation_type_t type;
public:
	Variation() : chromosome(""), coord1(0), coord2(0), p_value(1.0), sequence(""), type(NONE) {}
	Variation(std::string chromosome, size_t coord1, size_t coord2, double p_value, variation_type_t type) : chromosome(chromosome), coord1(coord1), coord2(coord2), p_value(p_value), sequence(""), type(type) {}
	Variation(std::string chromosome, size_t coord1, size_t coord2, double p_value, const std::string& sequence, variation_type_t type) : chromosome(chromosome), coord1(coord1), coord2(coord2), p_value(p_value), sequence(sequence), type(type) {}

	virtual ~Variation() {}

	std::string getChromosome() const { return chromosome; }

	/** Start position of the event. Refers to a position between two characters, i.e. a value
	 *  of n refers to the position between characters n-1 and n. */
	size_t getCoordinate1() const { return coord1; }

	/** In case of a deletion: end position (see getCoordinate1() for interpretation).
	 *  In case of an insertion: its length. */
	size_t getCoordinate2() const { return coord2; }

	variation_type_t getType() const { return type; }

	/** Returns the sequence of an insertion/deletion. If unknown, returns the empty string. */
	const std::string& getSequence() const { return sequence; }

	double getPValue() const { return p_value; }

	/** Returns the length difference caused by this variation, i.e. returns deviation of
	 *  the expected expected insert size of a read pair that maps over this variation.
	 */
	int getLengthDifference() const;

	/** Removes all information and in particular makes sure that getType() will return NONE. */
	void clear() { *this = Variation(); }

	/** Makes variation canonical, i.e. replaces it by the rightmost equivalent variation. */
	void canonify(const boost::unordered_map<std::string,NamedDnaSequence*>& references, bool rightmost = false);
	void canonify(const NamedDnaSequence* reference, bool rightmost = false);

	bool operator==(const Variation& v) const;
};

std::ostream& operator<<(std::ostream& os, const Variation& variation);

std::size_t hash_value(Variation const& v);

#endif /* VARIATION_H_ */
