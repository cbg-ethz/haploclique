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

#include "ReadSetZTester.h"

using namespace std;

ReadSetZTester::ReadSetZTester(double insert_length_popmean, double insert_length_popstddev) : insert_length_popmean(insert_length_popmean), insert_length_popstddev(insert_length_popstddev) {
}

double ReadSetZTester::computeSignificance(int insert_size_sum, int count, const std::vector<int>& read_group_counts) const {
	double mean_diff = (((double)insert_size_sum) - (count*insert_length_popmean)) / count;
	double z = mean_diff * sqrt(count) / insert_length_popstddev;
	boost::math::normal dist;
	return 2 * cdf(complement(dist, fabs(z)));
}
