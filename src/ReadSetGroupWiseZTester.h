/* Copyright 2013 Tobias Marschall
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

#ifndef READSETGROUPWISEZTESTER_H_
#define READSETGROUPWISEZTESTER_H_

#include <vector>

#include "ReadSetSignificanceTester.h"
#include "Types.h"

/** Performs a significance test on a set of reads from different read groups
  * where each read group's insert size has a Gaussian distribution. */
class ReadSetGroupWiseZTester : public ReadSetSignificanceTester {
private:
	const std::vector<mean_and_stddev_t>& distributions;
public:
	/** Construtor based on mean and standard deviation of Gaussian insert length distribution. */
	ReadSetGroupWiseZTester(const std::vector<mean_and_stddev_t>& distributions);
	/** Test ignores read group information and uses one common Gaussian as given at construction time. */
	virtual double computeSignificance(int insert_size_sum, int count, const std::vector<int>& read_group_counts) const;
	virtual ~ReadSetGroupWiseZTester() {}
};

#endif /* READSETGROUPWISEZTESTER_H_ */
