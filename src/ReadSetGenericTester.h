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

#ifndef READSETGENERICTESTER_H_
#define READSETGENERICTESTER_H_

#include <iostream>
#include <vector>

#include "HistogramBasedDistribution.h"
#include "ReadSetSignificanceTester.h"

/** Performs a significance test on a read group based on an arbitrary distribution. */
class ReadSetGenericTester : public ReadSetSignificanceTester {
private:
	typedef struct sum_distribution_t {
		int offset;
		// probability distribution
		std::vector<double>* distribution;
		// left-cumulative distribution
		std::vector<double>* cdf;
		// right-cumulative distribution
		std::vector<double>* ccdf;
		sum_distribution_t() : offset(0), distribution(0), cdf(0), ccdf(0) {}
		void free() {
			if (distribution != 0) delete distribution;
			if (cdf != 0) delete cdf;
			if (ccdf != 0) delete ccdf;
		}
	} sum_distribution_t;
	// distributions[i] contains the distribution of the sum of insert sizes of i+1 segments
	mutable std::vector<sum_distribution_t> distributions;
	/** Makes sure that cumulative distributions for sums of "count" insert sizes exists. */
	void ensureCumulativeDistribution(int count) const;
public:
	/** Construtor based on an arbitrary insert length distribution. */
	ReadSetGenericTester(const HistogramBasedDistribution& ild);
	/** Test ignores read group information and uses the distribution given at construction time. */
	virtual double computeSignificance(int insert_size_sum, int count, const std::vector<int>& read_group_counts) const;
	virtual ~ReadSetGenericTester();
	/** Caches distributions of insert size sums up to the given count. */
	virtual void cacheDistributions(int maxCount) const;
	virtual void writeCachedDistributions(std::ostream& os) const;
	virtual void readCachedDistributions(std::istream& is) const;
};

#endif /* READSETGENERICTESTER_H_ */
