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

#include <boost/math/distributions/normal.hpp>

#include "ReadSetGroupWiseZTester.h"

using namespace std;

ReadSetGroupWiseZTester::ReadSetGroupWiseZTester(const std::vector<mean_and_stddev_t>& distributions) : distributions(distributions) {
}

double ReadSetGroupWiseZTester::computeSignificance(int insert_size_sum, int count, const std::vector<int>& read_group_counts) const {
// 	double mean_diff = (((double)insert_size_sum) - (count*insert_length_popmean)) / count;
	assert(read_group_counts.size() <= distributions.size());
	// compute mean deviation from mean insert sizes
	double mean_diff = (double)insert_size_sum;
	int k = 0;
	for (size_t rg=0; rg<read_group_counts.size(); ++rg) {
		mean_diff -= read_group_counts[rg]*distributions[rg].mean;
		k += read_group_counts[rg];
	}
	assert(k == count);
	mean_diff /= count;
	// compute standard deviation
	double variance_of_sum = 0.0;
	for (size_t rg=0; rg<read_group_counts.size(); ++rg) {
		variance_of_sum += read_group_counts[rg] * distributions[rg].stddev * distributions[rg].stddev;
	}
	double stddev = sqrt(variance_of_sum) / count;
	double z = mean_diff / stddev;
	boost::math::normal dist;
	return 2 * cdf(complement(dist, fabs(z)));
}
