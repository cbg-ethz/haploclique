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

#include <math.h>
#include <limits>
#include <boost/math/distributions/normal.hpp>

#include "ReadGroupAwareEdgeCalculator.h"

using namespace std;

const double ReadGroupAwareEdgeCalculator::sqrt2 = 1.4142135623730951;
const size_t ReadGroupAwareEdgeCalculator::SF_CACHE_SIZE = 38500;
const double ReadGroupAwareEdgeCalculator::SF_CACHE_FACTOR = 0.001;

ReadGroupAwareEdgeCalculator::ReadGroupAwareEdgeCalculator(double significance_level, vector<mean_and_stddev_t>& distribution_params) : distribution_params(distribution_params) {
	this->significance_level = significance_level;
	cached_survival_function = new double[SF_CACHE_SIZE];
	boost::math::normal norm;
	for (size_t x=0; x<SF_CACHE_SIZE; ++x) {
		cached_survival_function[x] = boost::math::cdf(boost::math::complement(norm,x*SF_CACHE_FACTOR));
	}
	int rg_count = distribution_params.size();
	allowable_insert_size_diffs = vector<double>(rg_count*rg_count, 0.0);
	stddevs = vector<double>(rg_count*rg_count, 0.0);
	min_partner_length_diffs = vector<int>(rg_count, numeric_limits<int>::max());
	max_partner_length_diffs = vector<int>(rg_count, numeric_limits<int>::min());
	for (int rg1=0; rg1<rg_count; ++rg1) {
		for (int rg2=0; rg2<rg_count; ++rg2) {
			// standard deviation of the sum of two insert sizes (one from rg1 and one from rg2)
			double stddev = sqrt(distribution_params[rg1].stddev * distribution_params[rg1].stddev + distribution_params[rg2].stddev * distribution_params[rg2].stddev);
			// standard deviation of the mean (instead of the sum), therefore divide by two
			stddevs[rg1*rg_count + rg2] = stddev / 2.0;
			boost::math::normal norm(0,stddev);
			double allowed_diff = -boost::math::quantile(norm, significance_level/2.0);
			allowable_insert_size_diffs[rg1*rg_count + rg2] = allowed_diff;
			double mean_diff = distribution_params[rg2].mean - distribution_params[rg1].mean;
			int min_partner_length_diff = (int)ceil(- allowed_diff + mean_diff);
			int max_partner_length_diff = (int)floor(allowed_diff + mean_diff);
			min_partner_length_diffs[rg1] = min(min_partner_length_diffs[rg1], min_partner_length_diff);
			max_partner_length_diffs[rg1] = max(max_partner_length_diffs[rg1], max_partner_length_diff);
// 			cerr << "rg1: " << rg1 << ", rg2: " << rg2 << endl;
// 			cerr << "   min_partner_length_diff: " << min_partner_length_diff << ", max_partner_length_diff: " << max_partner_length_diff << endl;
// 			cerr << "   range: " << distribution_params[rg1].mean + min_partner_length_diff << " - " << distribution_params[rg1].mean + max_partner_length_diff << endl;
// 			cerr << "   allowed_diff: " << allowed_diff << endl;
		}
// 		cerr << "Extremes for rg " << rg1 << ": min: " << min_partner_length_diffs[rg1] << ", max: " << max_partner_length_diffs[rg1] << endl;
	}
}

ReadGroupAwareEdgeCalculator::~ReadGroupAwareEdgeCalculator() {
	delete [] cached_survival_function;
}

double ReadGroupAwareEdgeCalculator::sf(double x) const {
	if (x>=0.0) {
		size_t i = (size_t)floor(x/SF_CACHE_FACTOR);
		return (i<SF_CACHE_SIZE)?cached_survival_function[i]:0.0;
	} else {
		size_t i = (size_t)ceil(-x/SF_CACHE_FACTOR);
		return (i<SF_CACHE_SIZE)?(1.0-cached_survival_function[i]):1.0;
	}
}

bool ReadGroupAwareEdgeCalculator::length_compatible(int rg1, int length1, int rg2, int length2) const {
	double diff = abs((length1 - distribution_params[rg1].mean) - (length2 - distribution_params[rg2].mean));
	return diff <= allowable_insert_size_diffs[rg1*distribution_params.size() + rg2];
}

void ReadGroupAwareEdgeCalculator::getPartnerLengthRange(const AlignmentRecord& ap, unsigned int *min, unsigned int *max) const {
	assert(ap.getReadGroup() >= 0);
	assert(ap.getReadGroup() < distribution_params.size());
	int min_int = ap.getInsertLength() + min_partner_length_diffs[ap.getReadGroup()];
	if (min_int < 0) *min = 0;
	else *min = (unsigned int)min_int;
	*max = ap.getInsertLength() + max_partner_length_diffs[ap.getReadGroup()];
}

bool ReadGroupAwareEdgeCalculator::edgeBetween(const AlignmentRecord & ap1, const AlignmentRecord & ap2) const {
	int rg1 = ap1.getReadGroup();
	int rg2 = ap2.getReadGroup();
// 	cerr << "    ReadGroupAwareEdgeCalculator::edgeBetween: ";
	if (!length_compatible(rg1, (int)ap1.getInsertLength(), rg2, (int)ap2.getInsertLength())) {
// 		cerr << " NO (length)" << endl;
		return false;
	}
	size_t intersection_length = ap1.intersectionLength(ap2);
	if (intersection_length == 0) {
// 		cerr << " NO (intersection=0)" << endl;
		return false;
	}
	double mean_deviation = ((ap1.getInsertLength() - distribution_params[rg1].mean) + (ap2.getInsertLength() - distribution_params[rg2].mean)) / 2.0;
	int rg_count = distribution_params.size();
	double stddev = stddevs[rg1*rg_count + rg2];
	double interedge = min(1.0, 2.0 * sf((mean_deviation - intersection_length ) / stddev));
	bool result = interedge >= significance_level;
// 	cerr << (result?" YES":"NO (intersection too small)") << endl;
	return result;
}
