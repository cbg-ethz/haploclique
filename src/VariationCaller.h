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

#ifndef VARIATIONCALLER_H_
#define VARIATIONCALLER_H_

#include <vector>
#include <algorithm>
#include <limits>

#include <boost/math/distributions/normal.hpp>

#include "Variation.h"
#include "AlignmentRecord.h"
#include "ReadSetSignificanceTester.h"

/** Class providing methods to call a variation given a set of reads supporting
 *  the same phenomenon. */
class VariationCaller {
private:
	double insert_length_popmean;
	const std::vector<mean_and_stddev_t>* distributions;
	const ReadSetSignificanceTester& significance_tester;
	mutable std::vector<double> log_factorial_table;
	// Comparison class to sort alignment pairs by their single p-values
	typedef struct alignment_pvalue_sort_t {
		bool operator()(const std::pair<const AlignmentRecord*, double>& a1, const std::pair<const AlignmentRecord*, double>& a2) {
			return a1.second > a2.second;
		}
	} alignment_pvalue_sort_t;
	// Comparison class to sort alignment pairs by weight (from high to low)
	typedef struct alignment_weight_sort_t {
		bool operator()(const AlignmentRecord* a1, const AlignmentRecord* a2) {
			return a1->getWeight() > a2->getWeight();
		}
		alignment_weight_sort_t() {}
	} alignment_weight_sort_t;
	double logFactorial(size_t n) const;
	/** The logarithm of n choose k. */
	double logBinomial(unsigned int n, unsigned int k) const;
	/** Given log(x) and log(y), return log(x+y). */
	double logAdd(double log_x, double log_y) const;
	double pvalueUpperBound(const std::vector<const AlignmentRecord*>& alignments) const ;
	double log1MinusX(double x) const;
public:
	typedef struct {
		// start, end, and length of the intersection of all insert segments
		unsigned int insert_start;
		unsigned int insert_end;
		int insert_length;
		// difference between weighted mean insert length and population mean
		double diff;
		// sum of weights of alignments
		double total_weight;
	} additional_stats_t;

	VariationCaller(double insert_length_popmean, const ReadSetSignificanceTester& significance_tester);
	VariationCaller(const std::vector<mean_and_stddev_t>* distributions, const ReadSetSignificanceTester& significance_tester);
	virtual ~VariationCaller() {}

	/** Predicts a variation from a set of (compatible) alignment pairs.
	 *  Parameters begin and end must be iterators yielding pointers to
	 *  instances of AlignmentRecord. If "stats" is non-null, additional
	 *  statistics are delivered to that struct. */
	template <typename IteratorType>
	Variation call(IteratorType begin, IteratorType end, additional_stats_t* stats) const {
		additional_stats_t additional_stats;
		if (stats == 0) stats = &additional_stats;
		stats->insert_start = 0;
		stats->insert_end = 0;
		double weighted_mean_diff = 0.0;
		stats->total_weight = 0.0;
		std::string chromosome = "";
		for (IteratorType it(begin); it != end; ++it) {
			const AlignmentRecord& ap = **it;
			assert(!ap.isSingleEnd());
			if (it == begin) {
				stats->insert_start = ap.getInsertStart();
				stats->insert_end = ap.getInsertEnd();
				chromosome = ap.getChromosome();
			} else {
				stats->insert_start = std::max(stats->insert_start, ap.getInsertStart());
				stats->insert_end = std::min(stats->insert_end, ap.getInsertEnd());
				assert(chromosome.compare(ap.getChromosome()) == 0);
			}
			if (distributions == 0) {
				weighted_mean_diff += (ap.getInsertLength() - insert_length_popmean) * ap.getWeight();
			} else {
				int rg = ap.getReadGroup();
				assert(rg < distributions->size());
				weighted_mean_diff += (ap.getInsertLength() - distributions->at(rg).mean) * ap.getWeight();
			}
			stats->total_weight += ap.getWeight();
		}
		weighted_mean_diff /= stats->total_weight;
		stats->diff = weighted_mean_diff;
		stats->insert_length = ((int)stats->insert_end) - stats->insert_start + 1;
		assert(stats->insert_length >= 0);
		std::vector<const AlignmentRecord*> alignments(begin, end);
		double pvalue = pvalueUpperBound(alignments);
		// insertion or deletion?
		int indel_start = 0;
		int indel_end = 0;
		int int_diff = (int)roundf(stats->diff);
		if (int_diff>0) {
			// deletion
			if (stats->insert_length > stats->diff) {
				indel_start = stats->insert_start + (stats->insert_length - int_diff)/2;
				indel_end = indel_start + int_diff;
			} else {
				indel_start = stats->insert_start;
				indel_end = stats->insert_end + 1;
			}
			return Variation(chromosome, indel_start, indel_end, pvalue, Variation::DELETION);
		} else {
			indel_start = ((stats->insert_end+1) + stats->insert_start) / 2;
			return Variation(chromosome, indel_start, -int_diff, pvalue, Variation::INSERTION);
		}
	}
};

#endif /* VARIATIONCALLER_H_ */
