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


#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <memory>
#include <vector>

class Distributions {
public:
	/** Computes the cumulative distribution, i.e. the left-tail sums. */
	static std::unique_ptr<std::vector<double> > toCDF(const std::vector<double>& distribution);
	/** Computes the complementary cumulative distribution, i.e. the 
	 *  right-tail sums. */
	static std::unique_ptr<std::vector<double> > toCCDF(const std::vector<double>& distribution);
	/** Computes the convolution of two distributions. The size of the
	 *  resulting array is chosen such that all values that are >0.0 in
	 *  double precision are contained. */
	static std::unique_ptr<std::vector<double> > convolve(const std::vector<double>& dist1, const std::vector<double>& dist2, int offset1, int offset2, int* offset_result);
};

#endif // DISTRIBUTIONS_H
