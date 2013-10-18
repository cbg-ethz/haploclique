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

#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include <memory>
#include <vector>
#include <map>

#include "HistogramBasedDistribution.h"

/** Class to estimate a insert length distribution from observations. */
class Histogram {
public:
	typedef std::map<int,unsigned int> histogram_t;
private:
	long long total_count;
	histogram_t histogram;
public:
	Histogram();
	virtual ~Histogram();

	void add(int insert_length);
	
	size_t size() const;
	
	std::auto_ptr<HistogramBasedDistribution> toDistribution(int min_counts_per_bin) const;
	
	void computeMeanAndStddev(double* mean, double* sd) const;
};

#endif /* HISTOGRAM_H_ */
