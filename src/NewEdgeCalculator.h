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

#ifndef NEWEDGECALCULATOR_H_
#define NEWEDGECALCULATOR_H_

#include <set>

#include "AlignmentRecord.h"
#include "EdgeCalculator.h"

using namespace std;

class NewEdgeCalculator : public EdgeCalculator {
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

public:
    NewEdgeCalculator(double Q, double edge_quasi_cutoff, double overlap, bool frameshift_merge, map<int, double>& simpson_map, double edge_quasi_cutoff_single, double overlap_single, double edge_quasi_cutoff_mixed);
    virtual ~NewEdgeCalculator();

    /** Decides whether an edge is to be drawn between the two given nodes. */
    virtual bool edgeBetween(const AlignmentRecord& ap1, const AlignmentRecord& ap2) const;

#endif /* NEWEDGECALCULATOR_H_ */
