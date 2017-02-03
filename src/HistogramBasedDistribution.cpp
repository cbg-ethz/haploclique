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

#include <iostream>
#include <fstream>
#include <math.h>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include "HistogramBasedDistribution.h"

using namespace std;

HistogramBasedDistribution::HistogramBasedDistribution(const std::string& filename) {
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer_t;
	boost::char_separator<char> separator(" \t");
	ifstream f(filename.c_str());
	if (f.fail()) {
		ostringstream oss;
		oss << "Could not open file \"" << filename << "\"";
		throw std::runtime_error(oss.str());
	}
	int last = -1;
	double psum = 0.0;
	string line;
	while (getline(f, line)) {
		tokenizer_t tokenizer(line,separator);
		vector<string> tokens(tokenizer.begin(), tokenizer.end());
		if (tokens.size()!=3) {
			throw std::runtime_error("Error parsing insert length distribution.");
		}
		try {
			int start = boost::lexical_cast<int>(tokens[0]);
			int end = boost::lexical_cast<int>(tokens[1]);
			double p = boost::lexical_cast<double>(tokens[2]);
			assert((last==-1) || (start==last+1));
			assert(start <= end);
			if (last==-1) min = start;
			max = end;
			last = end;
			psum += p;
			values.push_back(make_pair(start,p/(end-start+1)));
		} catch(boost::bad_lexical_cast &){
			throw std::runtime_error("Error parsing alignment pair.");
		}
	}
	if (fabs(1.0 - psum) > 1e-6) {
		throw std::runtime_error("Probabilities in insert size distribution do not sum to one");
	}
}

HistogramBasedDistribution::HistogramBasedDistribution(int min, int max, const vector<value_t>& values) : min(min), max(max), values(values) {
}

HistogramBasedDistribution::~HistogramBasedDistribution() {
}

double HistogramBasedDistribution::probability(int value) const {
	if (value < min) return 0.0;
	if (value > max) return 0.0;
	typedef vector<value_t>::const_iterator iter_t;
	iter_t cit = upper_bound(values.begin(), values.end(), make_pair(value,2.0));
	--cit;
	return cit->second;
}

double HistogramBasedDistribution::normal_mass(int from, int to, const boost::math::normal& norm) const {
	return boost::math::cdf(norm,((double)to) + 0.5) - boost::math::cdf(norm,((double)from) - 0.5);
}

void HistogramBasedDistribution::printWithGaussian(ostream& os, double mean, double stddev) const {
	boost::math::normal norm(mean, stddev);
	vector<HistogramBasedDistribution::value_t>::const_iterator cit = values.begin();
	int start = -1;
	double p = -1.0;
	while (true) {
		if (cit != values.begin()) {
			int length = cit->first - start;
			os << start << ' ' << (cit->first-1) << ' ' << (p*length) << ' ' << normal_mass(start,cit->first-1,norm) << endl;
		}
		start = cit->first;
		p = cit->second;
		++cit;
		if (cit == values.end()) {
			int length = max - start + 1;
			os << start << ' ' << max << ' ' << (p*length) << ' ' << normal_mass(start,max,norm) << endl;
			break;
		}
	}
}

std::unique_ptr<std::vector<double> > HistogramBasedDistribution::toDistribution(double tailCutoff, int* offset) const {
	// determine start point
	int min_insert_length = min;
	int min_idx = -1;
	double p_sum = 0.0;
	for (size_t i=0; i<values.size(); ++i) {
		int start = values[i].first;
		int end = (i+1 < values.size())? values[i+1].first-1:max;
		int length = end - start + 1;
		double p = values[i].second * length;
		if (p_sum + p > tailCutoff) {
			min_insert_length = start + int((tailCutoff - p_sum) / values[i].second);
			min_idx = i;
			break;
		}
		p_sum += p;
	}
	// determine end point
	int max_insert_length = max;
	p_sum = 0.0;
	for (size_t i=values.size()-1; i>=0; --i) {
		int start = values[i].first;
		int end = (i+1 < values.size())? values[i+1].first-1:max;
		int length = end - start + 1;
		double p = values[i].second * length;
		if (p_sum + p > tailCutoff) {
			max_insert_length = end - int((tailCutoff - p_sum) / values[i].second);
			break;
		}
		p_sum += p;
	}
	// create result
	unique_ptr<vector<double> > result(new vector<double>());
	int k = min_idx;
	int i = min_insert_length;
	while (i <= max_insert_length) {
		int start = values[k].first;
		int end = (k+1 < values.size())? values[k+1].first-1:max;
		assert(start <= i);
		if (i > end) {
			k += 1;
			if (k >= values.size()) break;
			else continue;
		}
		result->push_back(values[k].second);
		i += 1;
	}
	*offset = min_insert_length;
	return result;
}

int HistogramBasedDistribution::getMedian() {
	double p_sum = 0.0;
	for (size_t i=0; i<values.size(); ++i) {
		int start = values[i].first;
		int end = (i+1 < values.size())? values[i+1].first-1:max;
		int length = end - start + 1;
		double p = values[i].second * length;
		if (p_sum + p >= 0.5) {
			return (int)roundf(start + ((double)length)*(0.5-p_sum)/p);
		}
		p_sum += p;
	}
	assert(false);
}

void HistogramBasedDistribution::estimateGaussian(double* mean, double* stddev) {
	double p_sum = 0.0;
	double q25 = 0;
	bool q25_set = false;
	double q50 = 0;
	bool q50_set = false;
	double q75 = 0;
	for (size_t i=0; i<values.size(); ++i) {
		int start = values[i].first;
		int end = (i+1 < values.size())? values[i+1].first-1:max;
		int length = end - start + 1;
		double p = values[i].second * length;
		if (!q25_set && (p_sum + p >= 0.25)) {
			q25 = start + ((double)length)*(0.25-p_sum)/p;
		}
		if (!q50_set && (p_sum + p >= 0.5)) {
			q50 = start + ((double)length)*(0.5-p_sum)/p;
		}
		if (p_sum + p >= 0.75) {
			q75 = start + ((double)length)*(0.75-p_sum)/p;
			break;
		}
		p_sum += p;
	}
	*mean = q50;
	*stddev = (q75 - q25) / 2 / 0.6744898;
}

std::ostream& operator<<(std::ostream& os, const HistogramBasedDistribution& distribution) {
	vector<HistogramBasedDistribution::value_t>::const_iterator cit = distribution.values.begin();
	int start = -1;
	double p = -1.0;
	while (true) {
		if (cit != distribution.values.begin()) {
			int length = cit->first - start;
			os << start << ' ' << (cit->first-1) << ' ' << (p*length) << endl;
		}
		start = cit->first;
		p = cit->second;
		++cit;
		if (cit == distribution.values.end()) {
			int length = distribution.max - start + 1;
			os << start << ' ' << distribution.max << ' ' << (p*length) << endl;
			break;
		}
	}
	return os;
}
