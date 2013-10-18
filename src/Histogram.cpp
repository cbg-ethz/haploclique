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

#include "Histogram.h"

using namespace std;

Histogram::Histogram() {
	total_count = 0;
}

Histogram::~Histogram() {
}

void Histogram::add(int insert_length) {
	if (histogram.find(insert_length) == histogram.end()) {
		histogram[insert_length] = 1;
	} else {
		histogram[insert_length] += 1;
	}
	total_count += 1;
}

size_t Histogram::size() const {
	return total_count;
}

std::auto_ptr<HistogramBasedDistribution> Histogram::toDistribution(int min_counts_per_bin) const {
	vector<HistogramBasedDistribution::value_t> values;
	int count = 0;
	// sum of counts until now
	int cumulated_count = 0;
	histogram_t::const_iterator it = histogram.begin();
	int interval_start = it->first;
	int min = interval_start;
	int max = min;
	int last = -1;
	for (; it != histogram.end(); ++it) {
		count += it->second;
		cumulated_count += it->second;
		if ((count>min_counts_per_bin) && (total_count-cumulated_count>min_counts_per_bin)) {
			int length = it->first - interval_start + 1;
			values.push_back(HistogramBasedDistribution::value_t(interval_start, double(count)/total_count/length));
			max = it->first;
			count = 0;
			interval_start = it->first + 1;
		}
		last = it->first;
	}
	if (count>0) {
		int length = last - interval_start + 1;
		values.push_back(HistogramBasedDistribution::value_t(interval_start, double(count)/total_count/length));
		max = last;
	}
	return auto_ptr<HistogramBasedDistribution>(new HistogramBasedDistribution(min, max, values));
}

/** Compute mean and stddev of the Gaussian making up the main peak (based on 0.25, 0.5, and 0.75 quantile). */
void Histogram::computeMeanAndStddev(double* mean, double* sd) const {
	// first determine the 0.25, 0.5, and 0.75 quantile
	histogram_t::const_iterator it = histogram.begin();
	long long cumulated_count = 0;
	double q25 = 0;
	bool q25_set = false;
	double q50 = 0;
	bool q50_set = false;
	double q75 = 0;
	for (; it != histogram.end(); ++it) {
		long long new_cumulated_count = cumulated_count + it->second;
		if (!q25_set && (new_cumulated_count >= 0.25*total_count)) {
			q25 = ((double)it->first) - 0.5 + (0.25*total_count - cumulated_count) / it->second;
			q25_set = true;
		}
		if (!q50_set && (new_cumulated_count >= 0.50*total_count)) {
			q50 = ((double)it->first) - 0.5 + (0.50*total_count - cumulated_count) / it->second;
			// q50 = it->first;
			q50_set = true;
		}
		if (new_cumulated_count >= 0.75*total_count) {
			q75 = ((double)it->first) - 0.5 + (0.75*total_count - cumulated_count) / it->second;
			// q75 = it->first;
			break;
		}
		cumulated_count = new_cumulated_count;
	}
	*mean = q50;
	*sd = (q75 - q25) / 2 / 0.6744898;
}
