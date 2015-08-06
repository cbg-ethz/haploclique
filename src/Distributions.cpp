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

#include <algorithm>

#include "Distributions.h"

using namespace std;

unique_ptr<vector<double> > Distributions::toCDF(const vector<double>& distribution) {
	unique_ptr<vector<double> > result(new vector<double>(distribution.size(), 0.0));
	if (distribution.size() == 0) return result;
	result->at(0) = distribution[0];
	for (unsigned int i=1; i<distribution.size(); ++i) {
		result->at(i) = distribution[i] + result->at(i-1);
	}
	return result;
}

unique_ptr<vector<double> > Distributions::toCCDF(const vector<double>& distribution) {
	unique_ptr<vector<double> > result(new vector<double>(distribution.size(), 0.0));
	if (distribution.size() == 0) return result;
	result->at(distribution.size()-1) = distribution[distribution.size()-1];
	for (int i=distribution.size()-2; i>=0; --i) {
		result->at(i) = distribution[i] + result->at(i+1);
	}
	return result;
}

unique_ptr<vector<double> > Distributions::convolve(const vector<double>& dist1, const vector<double>& dist2, int offset1, int offset2, int* offset_result) {
	unique_ptr<vector<double> > result(new vector<double>());
	*offset_result = offset1 + offset2;
	int rightmost_nonzero = -1;
	for (int n = 0; n <= ((int)(dist1.size()+dist2.size()))-2; ++n) {
		double p = 0.0;
		for (int i = max(0, n-((int)dist2.size())+1); i <= min(n, ((int)dist1.size())-1); ++i) {
			p += dist1[i]*dist2[n-i];
		}
		if ((result->size() == 0) && (p == 0.0)) {
			*offset_result += 1;
			continue;
		}
		if (p != 0.0) {
			rightmost_nonzero = result->size();
		}
		result->push_back(p);
	}
	result->resize(rightmost_nonzero+1);
	return result;
}
