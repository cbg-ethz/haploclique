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


#include "IndelLengthDistribution.h"

using namespace std;

IndelLengthDistribution::IndelLengthDistribution(const HistogramBasedDistribution& distribution) {
	for (int i=0; i<=distribution.maxValue(); ++i) {
		int phred = (int)round(-log10(distribution.probability(i))*10.0);
		phred_costs.push_back(phred);
	}
}

IndelLengthDistribution::IndelLengthDistribution(const std::vector<int>& phred_costs) : phred_costs(phred_costs) {
}

int IndelLengthDistribution::getPhredCost(int length) const {
	if (length < (int)phred_costs.size()) {
		return phred_costs[length];
	} else {
		// If indel is too long, return cost for longest known indel plus 10
		return phred_costs[phred_costs.size()-1] + 10;
	}
}

std::ostream& operator<<(std::ostream& os, const IndelLengthDistribution& ild) {
	os << '[';
	for (size_t i=0; i<ild.phred_costs.size(); ++i) {
		if (i>0) os << ',';
		os << ild.phred_costs[i];
	}
	os << ']';
}