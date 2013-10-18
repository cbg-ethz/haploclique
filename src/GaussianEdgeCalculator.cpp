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
#include <boost/math/distributions/normal.hpp>

#include "GaussianEdgeCalculator.h"

using namespace std;

const double GaussianEdgeCalculator::sqrt2 = 1.4142135623730951;
const size_t GaussianEdgeCalculator::SF_CACHE_SIZE = 38500;
const double GaussianEdgeCalculator::SF_CACHE_FACTOR = 0.001;

GaussianEdgeCalculator::GaussianEdgeCalculator(double significance_level, double insert_size_mean, double insert_size_stddev) {
	this->significance_level = significance_level;
	this->insert_size_popmean = insert_size_mean;
	this->insert_size_popstddev = insert_size_stddev;
	cached_survival_function = new double[SF_CACHE_SIZE];
	boost::math::normal norm;
	for (size_t x=0; x<SF_CACHE_SIZE; ++x) {
		cached_survival_function[x] = boost::math::cdf(boost::math::complement(norm,x*SF_CACHE_FACTOR));
	}
	boost::math::normal norm2(0,sqrt2*insert_size_stddev);
	this->allowable_insert_size_diff = (int)floor(-boost::math::quantile(norm2, significance_level/2.0));
	// cout << "allowable insert size difference: " << allowable_insert_size_diff << endl;
}

GaussianEdgeCalculator::~GaussianEdgeCalculator() {
	delete [] cached_survival_function;
}

double GaussianEdgeCalculator::sf(double x) const {
	if (x>=0.0) {
		size_t i = (size_t)floor(x/SF_CACHE_FACTOR);
		return (i<SF_CACHE_SIZE)?cached_survival_function[i]:0.0;
	} else {
		size_t i = (size_t)ceil(-x/SF_CACHE_FACTOR);
		return (i<SF_CACHE_SIZE)?(1.0-cached_survival_function[i]):1.0;
	}
}

void GaussianEdgeCalculator::getPartnerLengthRange(const AlignmentRecord& ap, unsigned int *min, unsigned int *max) const {
	if (allowable_insert_size_diff > ap.getInsertLength()) {
		*min = 0;
	} else {
		*min = ap.getInsertLength() - allowable_insert_size_diff;
	}
	*max = ap.getInsertLength() + allowable_insert_size_diff;
}

bool GaussianEdgeCalculator::edgeBetween(const AlignmentRecord & ap1, const AlignmentRecord & ap2) const {
	if (ap1.isSingleEnd() || ap2.isSingleEnd()) {
		throw runtime_error("Cannot process single-end reads in GaussianEdgeCalculator!");
	}
	double insert_length_diff = abs(((int)ap1.getInsertLength()) - ((int)ap2.getInsertLength()));
	// double meanedge = 2.0 * sf( insert_length_diff / (sqrt2*insert_size_popstddev) );
	if (insert_length_diff > allowable_insert_size_diff) {
		return false;
	}
	size_t intersection_length = ap1.internalSegmentIntersectionLength(ap2);
	if (intersection_length == 0) return false;
	double mean_insert_length = (ap1.getInsertLength() + ap2.getInsertLength()) / 2.0;
	double interedge = min(1.0, 2.0 * sf(sqrt2*(mean_insert_length - intersection_length - insert_size_popmean) / insert_size_popstddev) );
	return interedge >= significance_level;
}


