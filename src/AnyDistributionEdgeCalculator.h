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

#ifndef ANYDISTRIBUTIONEDGECALCULATOR_H_
#define ANYDISTRIBUTIONEDGECALCULATOR_H_

#include <vector>

#include "HistogramBasedDistribution.h"
#include "AlignmentRecord.h"
#include "EdgeCalculator.h"

class AnyDistributionEdgeCalculator : public EdgeCalculator {
private:
	double significance_level;
	// if two alignment pairs have a larger insert length differences,
	// then they are incompatible and will not receive an edge.
	int allowable_insert_size_diff;
	std::vector<double>* insert_size_dist;
	// offset of insert_size_dist, i.e. P(insert-size = k) = insert_size_dist[k - offset] 
	int insert_size_dist_offset;
	// right-cumulative distribution of the sum of two independently drawn insert sizes
	std::vector<double>* insert_size_sum_ccdf;
	// offset for insert_size_sum_dist
	int insert_size_sum_ccdf_offset;
	/** Returns the probability that the sum of two (independently drawn) insert sizes is >=k. */
	double insertSizeSumRightTail(int k) const;
	int computeAllowableInsertSizeDiff();
public:
	AnyDistributionEdgeCalculator(double significance_level, const HistogramBasedDistribution& insert_size_dist);
	virtual ~AnyDistributionEdgeCalculator();

	/** Decides whether an edge is to be drawn between the two given nodes. */
	virtual bool edgeBetween(const AlignmentRecord& ap1, const AlignmentRecord& ap2) const;

	/** Compute a length range. An alignment pair with a length outside this range is
	 *  guaranteed not to have an edge to the given pair ap. */
	virtual void getPartnerLengthRange(const AlignmentRecord& ap, unsigned int* min, unsigned int* max) const;

};

#endif /* ANYDISTRIBUTIONEDGECALCULATOR_H_ */
