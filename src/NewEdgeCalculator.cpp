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

std::vector<int> NewEdgeCalculator::commonPositions(const AlignmentRecord::covmap & cov_ap1, const AlignmentRecord::covmap & cov_ap2) const{
    std::vector<int> res;
    for (auto it = cov_ap1.begin(); it != cov_ap1.end(); ++it){
            auto fit = cov_ap2.find(it->first);
            if (fit!= cov_ap2.end()){
                res.push_back(it->first);
            }
    }
    std::sort(res.begin(),res.end());
    return res;
}

double NewEdgeCalculator::qScore(const AlignmentRecord::mapValue& value, char x) const{
    if (value.base == x){
        return 1.0 - value.prob;
    } else {
        return value.prob/3.0;
    }
}

double NewEdgeCalculator::calculateProbM(const std::vector<int> & aub, const AlignmentRecord::covmap & cov_ap1, const AlignmentRecord::covmap & cov_ap2) const{
    double res = 1.0;
    for (auto i : aub){
        auto& v1 = cov_ap1.find(i)->second;
        auto& v2 = cov_ap2.find(i)->second;
        double sum = 0.0;
        sum += qScore(v1,'A')*qScore(v2,'A');
        sum += qScore(v1,'C')*qScore(v2,'C');
        sum += qScore(v1,'T')*qScore(v2,'T');
        sum += qScore(v1,'G')*qScore(v2,'G');
        res *= sum;
    }
    return res;
}

//TO DO: calculate allel frequency distribution beforehand
double NewEdgeCalculator::calculateProb0(const AlignmentRecord::covmap & cov_ap1,const AlignmentRecord::covmap & cov_ap2, int& counter) const{
    double res = 1.0;
    //find positions which are non common positions and look up probability in empirical allele frequency distribution
    for (auto it = cov_ap1.begin(); it != cov_ap1.end(); ++it){
        auto fit = cov_ap2.find(it->first);
        if (fit == cov_ap2.end()){
            counter++;
            auto k = this->SIMPSON_MAP.find(it->first);
            if (k != this->SIMPSON_MAP.end()){
                res *= k->second;
            } else {
                res *= 0.25;
            }
        }
    }
    for (auto it = cov_ap2.begin(); it != cov_ap2.end(); ++it){
        auto fit = cov_ap1.find(it->first);
        if(fit == cov_ap1.end()){
            counter++;
            auto k = this->SIMPSON_MAP.find(it->first);
            if (k != this->SIMPSON_MAP.end()){
                res *= k->second;
            } else {
                res *= 0.25;
            }
        }
    }
    return res;
}

//TO DO: find out whether gaps / insertions are compatible
bool NewEdgeCalculator::checkGaps(const AlignmentRecord::covmap & cov_ap1,const AlignmentRecord::covmap & cov_ap2,const std::vector<int> & aub) const{
    for (int i = 0; i < (signed)aub.size()-1; ++i){
        int ref_diff = aub[i+1]-aub[i];
        auto& cov_ap1_i = cov_ap1.find(aub[i])->second;
        auto& cov_ap1_i1 = cov_ap1.find(aub[i+1])->second;
        auto& cov_ap2_i = cov_ap2.find(aub[i])->second;
        auto& cov_ap2_i1 = cov_ap2.find(aub[i+1])->second;
        int pos_diff1 = cov_ap1_i1.pir-cov_ap1_i.pir; //<0 if jump is given
        int pos_diff2 = cov_ap2_i1.pir-cov_ap2_i.pir; //<0 if jump is given
        bool jump1 = cov_ap1_i1.read-cov_ap1_i.read; //=0/1 for no jump/jump
        bool jump2 = cov_ap2_i1.read-cov_ap2_i.read;
        //insertion
        if(ref_diff == 1 && pos_diff1 != pos_diff2 && (jump1 || jump2) == 0){
            return false;
        }
        //deletion
        else if(ref_diff > 1 && pos_diff1 != pos_diff2 && (jump1 || jump2) == 0){
            return false;
        }
        else if(ref_diff > 1 && pos_diff1 == pos_diff2 && pos_diff1 > 1){
            return false;
        }
    }
    if (aub.empty()) return false;
    return true;
}

bool NewEdgeCalculator::similarityCriterion(const AlignmentRecord & a1, const AlignmentRecord::covmap & cov_ap1, const AlignmentRecord & a2, const AlignmentRecord::covmap & cov_ap2, std::vector<int> & aub) const{

    //Threshold for probability that reads were sampled from same haplotype
    double cutoff = 0;
    int tailLength = 0;
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
    if (aub.size()<=MIN_OVERLAP*std::min(cov_ap1.size(),cov_ap2.size())) return false;
    double p_m = calculateProbM(aub, cov_ap1, cov_ap2);
    double p_0 = calculateProb0(cov_ap1, cov_ap2,tailLength);
    double prob = p_m*p_0;
    int test = aub.size()+tailLength;
    double potence = 1.0/test;
    double final_prob = std::pow(prob,potence);
    //cout << "Final prob: " << final_prob << endl;
    return final_prob >= cutoff;
}

bool NewEdgeCalculator::edgeBetween(const AlignmentRecord & ap1, const AlignmentRecord & ap2) const{
    AlignmentRecord::covmap cov_ap1 = ap1.getCovmap();
    AlignmentRecord::covmap cov_ap2 = ap2.getCovmap();
    std::vector<int> aub = commonPositions(cov_ap1, cov_ap2);
    if (aub.size() == 0 || (!checkGaps(cov_ap1, cov_ap2, aub))){
        return false;
    }
    //cout << "Gaps are compatible" << endl;
    //std::vector<int> tail = tailPositions(cov_ap1,cov_ap2);
    //cout << "Tail positions computed" << endl;
    return similarityCriterion(ap1, cov_ap1, ap2, cov_ap2, aub);
}

void NewEdgeCalculator::getPartnerLengthRange(const AlignmentRecord& ap, unsigned int* min, unsigned int* max) const {
    assert(false);
}




