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

std::vector<int> commonPositions(const AlignmentRecord::covmap & cov_ap1, const AlignmentRecord::covmap & cov_ap2){
    std::vector<int> res;
    for (auto it = cov_ap1.begin(); it != cov_ap1.end(); ++it){
            auto fit = cov_ap2.find(it->first);
            if (fit!= cov_ap2.end()){
                res.push_back(it->first);
            }
    }
    return res;
}

std::vector<int> tailPositions(const AlignmentRecord::covmap & cov_ap1, const AlignmentRecord::covmap & cov_ap2){
    std::vector<int> res;
    for (auto it = cov_ap1.begin(); it != cov_ap1.end(); ++it){
        auto fit = cov_ap2.find(it->first);
        if (fit == cov_ap2.end()){
            res.push_back(it->first);
        }
    }
    for (auto it = cov_ap2.begin(); it != cov_ap2.end(); ++it){
        auto fit = cov_ap1.find(it->first);
        if (fit == cov_ap1.end()){
            res.push_back(it->first);
        }
    }
    std::sort(res.begin(),res.end());
    return res;
}

double qScore(AlignmentRecord::mapValue& value, char& x){
    if (value.base == x){
        return 1.0 - std::pow(10, (double)(-value.qual-33)/10.0);
    } else {
        return std::pow(10, (double)(-value.qual - 33)/10.0)/3.0;
    }
}

double calculateProbM(const std::vector<int> & aub, const AlignmentRecord::covmap & cov_ap1, const AlignmentRecord::covmap & cov_ap2){
    double res = 0.0;
    double sum = 0.0;
    string bases = "ACTG";
    for (auto i in aub){
        for (char j in bases){
            sum += qScore(cov_ap1[i],j)*qScore(cov_ap2[i],j);
        }
        res *= sum;
    }
    return res;
}

//TO DO: calculate allel frequency distribution beforehand
double calculateProb0(const std::vector<int> & tail, const AlignmentRecord::covmap & cov_ap1, const AlignmentRecord::covmap & cov_ap2){
    double res = 1.0;
    for(i in tail){
        res *= 0.25;
    }
    return res;
}

//TO DO: find out whether gaps are compatible
bool checkGaps(AlignmentRecord::covmap & cov_ap1,AlignmentRecord::covmap & cov_ap2, std::vector<int> & aub){
    bool res = true;
    for (int i = 0; i < aub.size()-1; ++i){
        int ref_diff1 = aub[i+1]-aub[i];
        int ref_diff2 = aub[i+1]-aub[i];
        int pos_diff1 = cov_ap1[aub[i+1]].pir-cov_ap1[aub[i]].pir;
        int pos_diff2 = cov_ap2[aub[i+1]].pir-cov_ap2[aub[i]].pir;
        bool jump1 = cov_ap1[aub[i+1]].read-cov_ap1[aub[i]].read; //=0/1 for no jump/jump
        bool jump2 = cov_ap2[aub[i+1]].read-cov_ap2[aub[i]].read;
    return res;
    }
}

bool similarityCriterion(const AlignmentRecord & a1, const AlignmentRecord::covmap & cov_ap1, const AlignmentRecord & a2, const AlignmentRecord::covmap & cov_ap2, std::vector<int> & aub, std::vector<int> tail){

    //Threshold forprobability that reads were sampled from same haplotype
    double cutoff = 0;
    if (a1.getName().find("Clique") != string::npos && a2.getName().find("Clique") != string::npos) {
        cutoff = this->EDGE_QUASI_CUTOFF;
    } else if (a1.getName().find("Clique") != string::npos || a2.getName().find("Clique") != string::npos) {
        cutoff = this->EDGE_QUASI_CUTOFF_MIXED;
    } else {
        cutoff = this->EDGE_QUASI_CUTOFF_SINGLE;
    }

    //Threshold for Overlap of Read Alignments
    double MIN_OVERLAP = 0;
    if (a1.getName().find("Clique") != string::npos && a2.getName().find("Clique") != string::npos) {
        MIN_OVERLAP = MIN_OVERLAP_CLIQUES;
    } else {
        MIN_OVERLAP = MIN_OVERLAP_SINGLE;
    }

    double p_m = calculateProbM(aub, cov_ap1, cov_ap2);
    double p_0 = calculateProb0(tail, cov_ap1, cov_ap2);
    double prob = p_m*p_0;
    int test = aub.size()+tail.size();
    double potence = 1.0/test;
    double final_prob = std::pow(prob,potence);

    return final_prob >= cutoff;
}

bool NewEdgeCalculator::edgeBetween(const AlignmentRecord & ap1, const AlignmentRecord & ap2) {
    AlignmentRecord::covmap cov_ap1 = ap1.coveredPositions();
    AlignmentRecord::covmap cov_ap2 = ap2.coveredPositions();
    std::vector<int> aub = commonPositions(cov_ap1, cov_ap2);
    if (!checkGaps(cov_ap1, cov_ap2, aub)){
        return false;
    }
    std::vector<int> tail = tailPositions(cov_ap1,cov_ap2);
    return similarityCriterion(ap1, cov_ap1, ap2, cov_ap2, aub, tail);
}




