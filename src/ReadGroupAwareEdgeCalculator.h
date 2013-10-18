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

#ifndef READGROUPAWAREEDGECALCULATOR_H_
#define READGROUPAWAREEDGECALCULATOR_H_

#include <vector>

#include "Types.h"
#include "AlignmentRecord.h"
#include "EdgeCalculator.h"

class ReadGroupAwareEdgeCalculator : public EdgeCalculator {
private:
	double significance_level;
	const std::vector<mean_and_stddev_t>& distribution_params;
	// if two alignment pairs have a larger insert length differences,
	// then they are incompatible and will not receive an edge.
	// when the smaller insert size is from read group rg1 and the larger from read group rg2,
	// then the corresponding allowable_insert_size_diff is stored at 
	// allowable_insert_size_diffs[rg1 * #readgroups + rg2]
	std::vector<double> allowable_insert_size_diffs;
	// for each read group, the minimum/maximum length difference of another alignment pair receiving an edge
	std::vector<int> min_partner_length_diffs;
	std::vector<int> max_partner_length_diffs;
	// stddevs[rg1 * #readgroups + rg2] is the standard deviation of the distribution of the mean insert size
	// of a read pair from rg1 and a read pair from rg2
	std::vector<double> stddevs;
	double* cached_survival_function;
	static const double sqrt2;
	static const size_t SF_CACHE_SIZE;
	static const double SF_CACHE_FACTOR;
	// survival function (aka complement cumulative distribution function of the standard normal distribution)
	double sf(double x) const;
	// decide whether two given alignment pairs are length-compatible
	bool length_compatible(int rg1, int length1, int rg2, int length2) const;
public:
	ReadGroupAwareEdgeCalculator(double significance_level, std::vector<mean_and_stddev_t>& distribution_params);
	virtual ~ReadGroupAwareEdgeCalculator();

	/** Decides whether an edge is to be drawn between the two given nodes. */
	virtual bool edgeBetween(const AlignmentRecord& ap1, const AlignmentRecord& ap2) const;

	/** Compute a length range. An alignment pair with a length outside this range is
	 *  guaranteed not to have an edge to the given pair ap. */
	virtual void getPartnerLengthRange(const AlignmentRecord& ap, unsigned int* min, unsigned int* max) const;

};

#endif /* READGROUPAWAREEDGECALCULATOR_H_ */
