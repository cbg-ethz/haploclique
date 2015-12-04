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

#include <math.h>
#include <boost/math/distributions/normal.hpp>
#include <stdlib.h>
#include <map>
#include <set>
#include <algorithm>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

#include "NewEdgeCalculator.h"

 using namespace std;
 using namespace boost;

 const double NewEdgeCalculator::FRAME_SHIFT_WEIGHT = 0.01;
 

 NewEdgeCalculator::NewEdgeCalculator(double Q, double edge_quasi_cutoff, double overlap, bool frameshift_merge, map<int, double>& simpson_map, double edge_quasi_cutoff_single, double overlap_single, double edge_quasi_cutoff_mixed) {
    this->Q = Q;
    this->EDGE_QUASI_CUTOFF = edge_quasi_cutoff;
    this->EDGE_QUASI_CUTOFF_SINGLE = edge_quasi_cutoff_single;
    this->EDGE_QUASI_CUTOFF_MIXED = edge_quasi_cutoff_mixed;
    this->MIN_OVERLAP_CLIQUES = overlap;
    this->MIN_OVERLAP_SINGLE = overlap_single;
    this->FRAMESHIFT_MERGE = frameshift_merge;
    this->SIMPSON_MAP = simpson_map;
}

NewEdgeCalculator::~NewEdgeCalculator() {
}

std::vector<int> commonPositions(const AlignmentRecord::covmap & cov_ap1, const AlignmentRecord::covmap & cov_ap2) const{

}

std::vector<int> tailPositions(const AlignmentRecord::covmap & cov_ap1, const AlignmentRecord::covmap & cov_ap2) const{

}

double calculateProbM(const std::vector<int> & aub, const AlignmentRecord::covmap & cov_ap1, const AlignmentRecord::covmap & cov_ap2) const{

}

double calculateProb0(const std::vector<int> & aub, const AlignmentRecord::covmap & cov_ap1, const AlignmentRecord::covmap & cov_ap2) const{

}


bool NewEdgeCalculator::similarityCriterion(const AlignmentRecord::covmap & cov_ap1, const AlignmentRecord::covmap & cov_ap2) const{
    bool incompatibleGaps = checkGaps(cov_ap1, cov_ap2);
    std::vector<int> aub = commonPositions(cov_ap1, cov_ap2);
    std::vector<int> tail = tailPositions(cov_ap1,cov_ap2);
    double p_m = calculateProbM(aub, cov_ap1, cov_ap2);
    double p_0 = calculateProb0(tail, cov_ap1, cov_ap2);
    double prob = p_m*p_0;
    int test = aub.size()+tail.size();
    double potence = 1.0/test;
    double final_prob = std::pow(prob,potence);
    return false;
}

bool NewEdgeCalculator::insertCriterion(const AlignmentRecord & ap1, const AlignmentRecord & ap1) const{
    return false;
}

bool NewEdgeCalculator::edgeBetween(const AlignmentRecord & ap1, const AlignmentRecord & ap2) const {
    AlignmentRecord::covmap cov_ap1 = ap1.coveredPositions();
    AlignmentRecord::covmap cov_ap2 = ap2.coveredPositions();
    bool insert = true;
    if (ap1.isPairedEnd() && ap2.isPairedEnd()){
        insert = insertCriterion(ap1, ap2);
    }
    bool similarity = similarityCriterion(cov_ap1, cov_ap2);
    return insert && similarity;
}




