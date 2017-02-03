/* 
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

#include <math.h>
#include <boost/math/distributions/normal.hpp>

#include "AnyDistributionEdgeCalculator.h"
#include "Distributions.h"

using namespace std;


AnyDistributionEdgeCalculator::AnyDistributionEdgeCalculator(double significance_level, const HistogramBasedDistribution& ild) : significance_level(significance_level) {
	this->insert_size_dist = ild.toDistribution(1e-4, &insert_size_dist_offset).release();
// 	cerr << "Using internal segment size distribution: " << endl;
// 	for (int i=0; i < (int)insert_size_dist->size(); ++i) {
// 		cerr << (i+insert_size_dist_offset) << " " << insert_size_dist->at(i) << endl;
// 	}
	// since we truncated the tails (a bit), we renormalize back to one
	double sum = 0.0;
	for (size_t i=0; i < insert_size_dist->size(); ++i) sum += insert_size_dist->at(i);
	for (size_t i=0; i < insert_size_dist->size(); ++i) insert_size_dist->at(i) /= sum;
	allowable_insert_size_diff = computeAllowableInsertSizeDiff();
	// cout << "allowable insert size difference: " << allowable_insert_size_diff << endl;
	unique_ptr<vector<double> > insert_size_sum_dist = Distributions::convolve(*insert_size_dist, *insert_size_dist, insert_size_dist_offset, insert_size_dist_offset, &insert_size_sum_ccdf_offset);
	this->insert_size_sum_ccdf = Distributions::toCCDF(*insert_size_sum_dist).release();
}

AnyDistributionEdgeCalculator::~AnyDistributionEdgeCalculator() {
	delete insert_size_dist;
	delete insert_size_sum_ccdf;
}

int AnyDistributionEdgeCalculator::computeAllowableInsertSizeDiff() {
	// first compute distribution of (absolute) differences of two independently sampled insert lengths
	vector<double> diff_dist(2*insert_size_dist->size(), 0.0);
	for (size_t i=0; i<insert_size_dist->size(); ++i) {
		diff_dist[0] += insert_size_dist->at(i) * insert_size_dist->at(i);
		for (size_t j=i+1; j<insert_size_dist->size(); ++j) {
			diff_dist[j-i] += 2 * (insert_size_dist->at(i) * insert_size_dist->at(j));
		}
	}
	// compute right-cumulative distribution
	for (int i=diff_dist.size()-2; i>=0; --i) {
		diff_dist[i] += diff_dist[i+1];
	}
	// determine largest difference k such that P(diff >= k) >= significance_level
	for (int k=1; k < (int)diff_dist.size(); ++k) {
		if (diff_dist[k] < significance_level) return k-1;
	}
	assert(false);
}

double AnyDistributionEdgeCalculator::insertSizeSumRightTail(int k) const {
	if (k<insert_size_sum_ccdf_offset) return 1.0;
	int i = k - insert_size_sum_ccdf_offset;
	if (i < (int)insert_size_sum_ccdf->size()) {
		return insert_size_sum_ccdf->at(i);
	} else {
		return 0.0;
	}
}

void AnyDistributionEdgeCalculator::getPartnerLengthRange(const AlignmentRecord& ap, unsigned int *min, unsigned int *max) const {
	if (allowable_insert_size_diff > (int)ap.getInsertLength()) {
		*min = 0;
	} else {
		*min = ap.getInsertLength() - allowable_insert_size_diff;
	}
	*max = ap.getInsertLength() + allowable_insert_size_diff;
}

bool AnyDistributionEdgeCalculator::edgeBetween(const AlignmentRecord & ap1, const AlignmentRecord & ap2) const {
	double insert_length_diff = abs(((int)ap1.getInsertLength()) - ((int)ap2.getInsertLength()));
// 	// double meanedge = 2.0 * sf( insert_length_diff / (sqrt2*insert_size_popstddev) );
	if (insert_length_diff > allowable_insert_size_diff) {
		return false;
	}
	size_t intersection_length = ap1.intersectionLength(ap2);
	if (intersection_length == 0) return false;
	int insert_length_sum = ap1.getInsertLength() + ap2.getInsertLength() - 2*intersection_length;
	double intersection_pvalue = insertSizeSumRightTail(insert_length_sum);
	return intersection_pvalue >= significance_level;
}
