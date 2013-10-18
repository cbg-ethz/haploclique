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
#include <limits>

#include "VariationCaller.h"

using namespace std;

VariationCaller::VariationCaller(double insert_length_popmean, const ReadSetSignificanceTester& significance_tester) : insert_length_popmean(insert_length_popmean), distributions(0), significance_tester(significance_tester), log_factorial_table(2,0.0) {
}

VariationCaller::VariationCaller(const std::vector<mean_and_stddev_t>* distributions, const ReadSetSignificanceTester& significance_tester) : insert_length_popmean(numeric_limits<double>::quiet_NaN()), distributions(distributions), significance_tester(significance_tester), log_factorial_table(2,0.0) {
}

double VariationCaller::logFactorial(size_t n) const {
	while (n>=log_factorial_table.size()) {
		log_factorial_table.push_back(log_factorial_table[log_factorial_table.size()-1] + log(log_factorial_table.size()));
	}
	return log_factorial_table[n];
}

double VariationCaller::logBinomial(unsigned int n, unsigned int k) const {
	assert(n>=k);
	return logFactorial(n) - logFactorial(k) - logFactorial(n-k);
}

double VariationCaller::logAdd(double log_x, double log_y) const {
	if (boost::math::isnan(log_x)) return log_x;
	if (boost::math::isnan(log_y)) return log_y;
	if (log_y > log_x) {
		if (boost::math::isinf(log_y)) return log_y;
		return log_y + boost::math::log1p(exp(log_x-log_y));
	} else {
		if (boost::math::isinf(log_x)) return log_x;
		return log_x + boost::math::log1p(exp(log_y-log_x));
	}
}

double VariationCaller::log1MinusX(double x) const {
	assert(x>=0.0);
	assert(x<=1.0);
	if (x==1.0) {
		return -numeric_limits<double>::infinity();
	}
	return boost::math::log1p(-x);
}

double VariationCaller::pvalueUpperBound(const vector<const AlignmentRecord*>& alignments) const {
	if (alignments.size() == 0) return 1.0;
	// create a list of alignments along with their individual p-values and sort by p-value
	vector<pair<const AlignmentRecord*,double> > alignments_and_pvalues;
	typedef vector<const AlignmentRecord*>::const_iterator alignment_it_t;
	for (alignment_it_t it = alignments.begin(); it != alignments.end(); ++it) {
		vector<int> read_group_counts;
		if ((*it)->getReadGroup() >= 0) {
			read_group_counts.assign((*it)->getReadGroup(), 0);
			read_group_counts.push_back(1);
		}
		double pvalue = significance_tester.computeSignificance((*it)->getInsertLength(), 1, read_group_counts);
		alignments_and_pvalues.push_back(make_pair(*it, pvalue));
	}
	sort(alignments_and_pvalues.begin(), alignments_and_pvalues.end(), alignment_pvalue_sort_t());
	
	typedef vector<pair<const AlignmentRecord*,double> >::const_iterator alignment_and_pvalue_it_t;
	// create a list of alignments sorted by weights
	vector<const AlignmentRecord*> alignments_by_weight(alignments.begin(), alignments.end());
	sort(alignments_by_weight.begin(), alignments_by_weight.end(), alignment_weight_sort_t());
	// partition the set of alignments into two groups: one with high weights and one with low weights
	size_t high_weight_count = 1;
	// upper bound for the weights in the "high weight" set
	double high_weight_upper = alignments_by_weight[0]->getWeight();
	while ((high_weight_count<alignments_by_weight.size()) && (alignments_by_weight[high_weight_count]->getWeight() >= 0.5 * high_weight_upper)) {
		high_weight_count += 1;
	}
	size_t low_weight_count = alignments.size() - high_weight_count;
	double high_weight_lower = alignments_by_weight[high_weight_count-1]->getWeight();
	double low_weight_upper = 0.0;
	double low_weight_lower = 0.0;
	if (high_weight_count<alignments.size()) {
		low_weight_upper = alignments_by_weight[high_weight_count]->getWeight();
		low_weight_lower = alignments_by_weight[alignments.size()-1]->getWeight();
	}
	alignments_by_weight[high_weight_count-1];
	// cout << "   VC: high-weight: " << high_weight_count << " / " << alignments_by_weight.size() << ", bounds: ";
	// cout << high_weight_upper << ", " << high_weight_lower << ", " << low_weight_upper << ", " << low_weight_lower << endl;
	double log_high_weight_upper = log(high_weight_upper);
	double log_1minus_high_weight_lower = log1MinusX(high_weight_lower);
	double log_low_weight_upper = log(low_weight_upper);
	double log_1minus_low_weight_lower = log1MinusX(low_weight_lower);
	double pvalue = 0.0;
	int insert_size_sum = 0;
	// keep track of the number of read pairs per read group
	vector<int> read_group_counts;
	// the sum of differences from the population mean over the first i alignments
	// double diff_sum = 0.0;
	for (size_t k=0; k<=alignments_and_pvalues.size(); ++k) {
		// compute upper bound for the z-test based p-value: P(H_0 | k most concordant reads are correctly aligned)
		double ztest_pvalue_bound = 1.0;
		if (k>0) {
			pair<const AlignmentRecord*, double>& alignment_and_pvalue = alignments_and_pvalues[k-1];
			insert_size_sum += alignment_and_pvalue.first->getInsertLength();
			int rg = alignment_and_pvalue.first->getReadGroup();
			if (rg >= 0) {
				while (read_group_counts.size() < rg + 1) {
					read_group_counts.push_back(0);
				}
				read_group_counts[rg] += 1;
			}
			ztest_pvalue_bound = significance_tester.computeSignificance(insert_size_sum, k, read_group_counts);
			if (ztest_pvalue_bound < std::numeric_limits<double>::min()) {
				ztest_pvalue_bound = std::numeric_limits<double>::min();
			}
		}
		// compute an upper bound for the (logarithm of the) probability that EXACTLY k alignments are correct
		double log_p_alignments = 0.0;
		if (k==0) {
			// case 1: all alignments are wrong
			for (size_t i=0; i<alignments.size(); ++i) {
				log_p_alignments += log1MinusX(alignments[i]->getWeight());
			}
		} else if (k==alignments.size()) {
			// case 2: all alignments are correct
			for (size_t i=0; i<alignments.size(); ++i) {
				log_p_alignments += log(alignments[i]->getWeight());
			}
		} else {
			// all other cases
			// i gives the number of alignments taken from the high weight set
			log_p_alignments = -numeric_limits<double>::infinity();
			for (size_t i=max(0,((int)k)-(int)low_weight_count); i<=min(k,high_weight_count); ++i) {
				double log_p = 0.0;
				log_p += logBinomial(high_weight_count,i);
				if (i>0) {
					log_p += i*log_high_weight_upper;
				}
				if (high_weight_count-i > 0) {
					log_p += (high_weight_count-i)*log_1minus_high_weight_lower;
				}
				log_p += logBinomial(low_weight_count,k-i);
				if (k-i > 0) {
					log_p += (k-i)*log_low_weight_upper;
				}
				if (low_weight_count-(k-i) > 0) {
					log_p += (low_weight_count-(k-i))*log_1minus_low_weight_lower;
				}
				// cout << "   VC:   i=" << i << ", log_p: " << log_p << endl;
				log_p_alignments = logAdd(log_p_alignments, log_p);
			}
		}
		double p = ztest_pvalue_bound*exp(log_p_alignments);
		// cout << "   VC: k=" << k << ": " << p << endl;
		pvalue += p;
	}
	// cout << "   VC: p-value: " << pvalue << endl;
	if (pvalue < std::numeric_limits<double>::min()) {
		pvalue = std::numeric_limits<double>::min();
	}
	if (pvalue > 1.0) {
		pvalue = 1.0;
	}
	return pvalue;
}
