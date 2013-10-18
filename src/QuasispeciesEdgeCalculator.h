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

#ifndef QUASISPECIESEDGECALCULATOR_H_
#define QUASISPECIESEDGECALCULATOR_H_

#include "AlignmentRecord.h"
#include "EdgeCalculator.h"

using namespace std;

class QuasispeciesEdgeCalculator : public EdgeCalculator {
private:
    
    typedef std::map<char, double> inner_map;
    typedef std::map<int, inner_map> outer_map;
    typedef std::map<int, outer_map> outer_outer_map;

    double MIN_OVERLAP;
    static const double FRAME_SHIFT_WEIGHT;
    double Q;
    double EDGE_QUASI_CUTOFF;
    bool FRAMESHIFT_MERGE;
    outer_map ALLEL_FREQUENCIES;
    outer_outer_map ALLEL_FREQUENCIES_INSERTIONS;

    typedef struct overlap_result {
        double probability;
        int hamming;
        int size;
    } overlap_result;
    
    // survival function (aka complement cumulative distribution function of the standard normal distribution)
    double sf(double x) const;
    double computeOverlap(const AlignmentRecord & ap1, const AlignmentRecord & ap2) const;
    QuasispeciesEdgeCalculator::overlap_result singleOverlap(const AlignmentRecord & ap1, const AlignmentRecord & ap2, int strain1, int strain2) const;
    int overlapSize(int e1, int e2, int s1, int s2) const;
public:
    QuasispeciesEdgeCalculator(double Q, double edge_quasi_cutoff, double overlap, bool frameshift_merge, outer_map& allel_frequencies, outer_outer_map& allel_frequencies_insertion_map);
    virtual ~QuasispeciesEdgeCalculator();

    /** Decides whether an edge is to be drawn between the two given nodes. */
    virtual bool edgeBetween(const AlignmentRecord& ap1, const AlignmentRecord& ap2) const;

    /** Compute a length range. An alignment pair with a length outside this range is
     *  guaranteed not to have an edge to the given pair ap. */
    virtual void getPartnerLengthRange(const AlignmentRecord& ap, unsigned int* min, unsigned int* max) const;

};

#endif /* QUASISPECIESEDGECALCULATOR_H_ */
