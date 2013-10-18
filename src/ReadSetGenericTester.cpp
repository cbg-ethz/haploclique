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

#include <cassert>
#include <iostream>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include "Distributions.h"
#include "ReadSetGenericTester.h"

using namespace std;

ReadSetGenericTester::ReadSetGenericTester(const HistogramBasedDistribution& ild) {
	sum_distribution_t entry;
	entry.distribution = ild.toDistribution(1e-4, &entry.offset).release();
	distributions.push_back(entry);
	// ensureDistribution(500);
}

ReadSetGenericTester::~ReadSetGenericTester() {
	vector<sum_distribution_t>::iterator it = distributions.begin();
	for (;it != distributions.end(); ++it) {
		if (it->distribution != 0) delete it->distribution;
		if (it->cdf != 0) delete it->cdf;
		if (it->ccdf != 0) delete it->ccdf;
	}
}

void ReadSetGenericTester::ensureCumulativeDistribution(int count) const {
	cacheDistributions(count);
	if (distributions[count-1].cdf == 0) {
		distributions[count-1].cdf = Distributions::toCDF(*distributions[count-1].distribution).release();
		distributions[count-1].ccdf = Distributions::toCCDF(*distributions[count-1].distribution).release();
		if ((count-1>0) && (count-1<(int)distributions.size()-1)) {
			delete distributions[count-1].distribution;
			distributions[count-1].distribution = 0;
		}
	}
}

void ReadSetGenericTester::cacheDistributions(int maxCount) const {
	assert(maxCount>0);
	while (maxCount-1 >= (int)distributions.size()) {
		int i = distributions.size() - 1;
		assert(distributions[i].distribution != 0);
		sum_distribution_t entry;
		entry.distribution = Distributions::convolve(*distributions[0].distribution, *distributions[i].distribution, distributions[0].offset, distributions[i].offset, &entry.offset).release();
		// cerr << "maxCount: " << (i+1) << ", size: " << entry.distribution->size() << endl;;
		distributions.push_back(entry);
	}
}

void ReadSetGenericTester::writeCachedDistributions(std::ostream& os) const {
	for (size_t i = 0; i<distributions.size(); ++i) {
		os << distributions[i].offset;
		const vector<double>& d = *distributions[i].distribution;
		for (size_t j = 0; j<d.size(); ++j) {
			os << " " << d[j];
		}
		os << endl;
	}
}

void ReadSetGenericTester::readCachedDistributions(std::istream& is) const {
	for (vector<sum_distribution_t>::iterator it = distributions.begin(); it != distributions.end(); ++it) {
		if (it->distribution != 0) delete it->distribution;
		if (it->cdf != 0) delete it->cdf;
		if (it->ccdf != 0) delete it->ccdf;
	}
	distributions.clear();
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer_t;
	boost::char_separator<char> whitespace_separator(" \t");
	assert(!is.fail());
	string line;
	while (getline(is,line)) {
		tokenizer_t tokenizer(line,whitespace_separator);
		sum_distribution_t entry;
		entry.distribution = new vector<double>();
		int n = 0;
		for (tokenizer_t::const_iterator it = tokenizer.begin(); it != tokenizer.end(); ++it) {
			if (n == 0) {
				entry.offset = boost::lexical_cast<int>(*it);
			} else {
				entry.distribution->push_back(boost::lexical_cast<double>(*it));
			}
			n += 1;
		}
		distributions.push_back(entry);
	}
}

double ReadSetGenericTester::computeSignificance(int insert_size_sum, int count, const std::vector<int>& read_group_counts) const {
	ensureCumulativeDistribution(count);
	const sum_distribution_t& entry = distributions[count-1];
	int i = insert_size_sum - entry.offset;
	if ((i<0) || (i>=(int)entry.cdf->size())) {
		return 0.0;
	}
	double one_sided_pvalue = min(entry.cdf->at(i), entry.ccdf->at(i));
	return 2 * one_sided_pvalue;
}
