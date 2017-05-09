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

#ifndef EDGECALCULATOR_H_
#define EDGECALCULATOR_H_

#include "AlignmentRecord.h"

class EdgeCalculator {
public:
	virtual ~EdgeCalculator() {}
	
	/** decides whether an edge is to be drawn between two given nodes. */
	virtual bool edgeBetween(const AlignmentRecord& ar1, const AlignmentRecord& ar2) const = 0;

	/** computes a length range. An alignment pair with a length outside this range is
	 *  guaranteed not to have an edge to the given pair ap. */
	virtual void getPartnerLengthRange(const AlignmentRecord& ar, unsigned int* min, unsigned int* max) const = 0;

    virtual void setOverlapCliques(double){}

};

#endif /* EDGECALCULATOR_H_ */
