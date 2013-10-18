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

#ifndef GAUSSIANEDGECALCULATOR_H_
#define GAUSSIANEDGECALCULATOR_H_

#include "AlignmentRecord.h"
#include "EdgeCalculator.h"

class GaussianEdgeCalculator : public EdgeCalculator {
private:
	double significance_level;
	double insert_size_popmean;
	double insert_size_popstddev;
	// if two alignment pairs have a larger insert length differences,
	// then they are incompatible and will not receive an edge.
	int allowable_insert_size_diff;
	double* cached_survival_function;
	static const double sqrt2;
	static const size_t SF_CACHE_SIZE;
	static const double SF_CACHE_FACTOR;
	// survival function (aka complement cumulative distribution function of the standard normal distribution)
	double sf(double x) const;
public:
	GaussianEdgeCalculator(double significance_level, double insert_size_mean, double insert_size_stddev);
	virtual ~GaussianEdgeCalculator();

	/** Decides whether an edge is to be drawn between the two given nodes. */
	virtual bool edgeBetween(const AlignmentRecord& ap1, const AlignmentRecord& ap2) const;

	/** Compute a length range. An alignment pair with a length outside this range is
	 *  guaranteed not to have an edge to the given pair ap. */
	virtual void getPartnerLengthRange(const AlignmentRecord& ap, unsigned int* min, unsigned int* max) const;

};

#endif /* GAUSSIANEDGECALCULATOR_H_ */
