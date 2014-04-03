/* Copyright 2012-2014 Tobias Marschall and Armin Töpfer
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

    double MIN_OVERLAP_CLIQUES;
    double MIN_OVERLAP_SINGLE;
    static const double FRAME_SHIFT_WEIGHT;
    double Q;
    double EDGE_QUASI_CUTOFF_MIXED;
    double EDGE_QUASI_CUTOFF_SINGLE;
    double EDGE_QUASI_CUTOFF;
    bool FRAMESHIFT_MERGE;
    map<int, double> SIMPSON_MAP;

    typedef struct overlap_result {
        double probability;
        int hamming;
        int size;
    } overlap_result;

    // survival function (aka complement cumulative distribution function of the standard normal distribution)
    double sf(double x) const;
    double computeOverlap(const AlignmentRecord & ap1, const AlignmentRecord & ap2, const double cutoff) const;
    double singleOverlap(const AlignmentRecord & ap1, const AlignmentRecord & ap2, int strain1, int strain2, double MIN_OVERLAP, const double cutoff) const;
    void computeOffsets(int s1, int s2, int e1, int e2, int &offset1, int &offset2, int &overlap) const;
    float overlapSize(int e1, int e2, int s1, int s2) const;
    string tail(std::string const& source, size_t const length) const;
    bool is_disjoint(const std::set<std::string> &set1, const std::set<std::string> &set2) const;
public:
    QuasispeciesEdgeCalculator(double Q, double edge_quasi_cutoff, double overlap, bool frameshift_merge, map<int, double>& simpson_map, double edge_quasi_cutoff_single, double overlap_single, double edge_quasi_cutoff_mixed);
    virtual ~QuasispeciesEdgeCalculator();

    /** Decides whether an edge is to be drawn between the two given nodes. */
    virtual bool edgeBetween(const AlignmentRecord& ap1, const AlignmentRecord& ap2) const;

    /** Compute a length range. An alignment pair with a length outside this range is
     *  guaranteed not to have an edge to the given pair ap. */
    virtual void getPartnerLengthRange(const AlignmentRecord& ap, unsigned int* min, unsigned int* max) const;
};

#endif /* QUASISPECIESEDGECALCULATOR_H_ */
