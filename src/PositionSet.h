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


#ifndef POSITIONSET_H
#define POSITIONSET_H

#include <boost/unordered_map.hpp>
#include <boost/dynamic_bitset.hpp>

#include <api/BamAux.h>

/** A set of genomic positions represented by bitvectors. */
class PositionSet {
private:
	typedef boost::dynamic_bitset<> chromosome_bitset_t;
	typedef boost::unordered_map<int,chromosome_bitset_t*> position_set_t;
	position_set_t position_set;
public:
	PositionSet();
	virtual ~PositionSet();
	
	/** Reads position set from file and adds these positions to the present set.
	 *  Positions in file are assumed to be 1-based and converted to 0-based positions.
	 *  @param ref_data Data structure providing mapping from reference names to indices.
	 */
	void loadFromFile(std::string filename, const BamTools::RefVector& ref_data);
	
	bool get(int chromosome_id, size_t position);
	void set(int chromosome_id, size_t position, bool value = true);
};

#endif // POSITIONSET_H
