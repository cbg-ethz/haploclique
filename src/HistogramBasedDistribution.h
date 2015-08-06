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

#ifndef HISTOGRAMBASEDDISTRIBUTION_H_
#define HISTOGRAMBASEDDISTRIBUTION_H_

#include <vector>
#include <boost/math/distributions/normal.hpp>

class HistogramBasedDistribution {
public:
	typedef std::pair<int,double> value_t;
private:
	int min;
	int max;
	// pairs of start position and probability
	std::vector<value_t> values;
	double normal_mass(int from, int to, const boost::math::normal& norm) const;
public:
	/** Reads a distribution from a file whose lines have the following format: <start> <end> <probability>
	 * 	where the last column must add up to one.*/
	HistogramBasedDistribution(const std::string& filename);
	HistogramBasedDistribution(int min, int max, const std::vector<value_t>& values);
	virtual ~HistogramBasedDistribution();
	double probability(int value) const;
	/** Minimum value for which a non-zero probability is returned. */
	int minValue() const { return min; }
	/** Maximum value for which a non-zero probability is returned. */
	int maxValue() const { return max; }
	int getMedian();
	
	/** Converts the distribution to a vector of doubles such that result[k] is the
	 *  probability of observing a value of offset+k, where offset is also set
	 *  by this function. A probability mass of (at most) tailCutoff is omitted from 
	 *  each end of the distribution to save space. */
	std::unique_ptr<std::vector<double> > toDistribution(double tailCutoff, int* offset) const;

	/** Prints the distribution along with analogous values obtained from the given normal distribution. */
	void printWithGaussian(std::ostream& os, double mean, double stddev) const;
	/** (Robustly) estimates a Gaussian distribution from the data and writes mean
	 *  and standard deviation to the given variables. */
	void estimateGaussian(double* mean, double* stddev);
	friend std::ostream& operator<<(std::ostream&, const HistogramBasedDistribution& distribution);
};

#endif /* HISTOGRAMBASEDDISTRIBUTION_H_ */
