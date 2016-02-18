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
 

NewEdgeCalculator::NewEdgeCalculator(double Q, double edge_quasi_cutoff, double overlap, bool frameshift_merge, unordered_map<int, double>& simpson_map, double edge_quasi_cutoff_single, double overlap_single, double edge_quasi_cutoff_mixed) {
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

double NewEdgeCalculator::qScore(const AlignmentRecord::mapValue& value, char x) const{
    if (value.base == x){
        return 1.0 - value.prob;
    } else {
        return value.prob/3.0;
    }
}

void NewEdgeCalculator::calculateProbM(const AlignmentRecord::mapValue & val1,const AlignmentRecord::mapValue & val2, double& res) const{
    double sum = 0.0;
    sum += qScore(val1,'A')*qScore(val2,'A');
    sum += qScore(val1,'C')*qScore(val2,'C');
    sum += qScore(val1,'T')*qScore(val2,'T');
    sum += qScore(val1,'G')*qScore(val2,'G');
    res *= sum;
}


void NewEdgeCalculator::calculateProb0(const AlignmentRecord::mapValue & val1, double& res) const{
    //find positions which are non common positions and look up probability in empirical allele frequency distribution
    const auto& k = this->SIMPSON_MAP.find(val1.ref);
    if (k != this->SIMPSON_MAP.end()){
        res *= k->second;
    } else {
        res *= 0.25;
    }
}

//TO DO: find out whether gaps / insertions are compatible
bool NewEdgeCalculator::checkGaps(const std::vector<AlignmentRecord::mapValue> & cov_ap1,const std::vector<AlignmentRecord::mapValue> & cov_ap2,const std::vector<std::pair<int,int>> & aub) const{
    for (int i = 0; i < (signed)aub.size()-1; ++i){
        auto& cov_ap1_i = cov_ap1[aub[i].first];
        auto& cov_ap1_i1 = cov_ap1[aub[i+1].first];
        auto& cov_ap2_i = cov_ap2[aub[i].second];
        auto& cov_ap2_i1 = cov_ap2[aub[i+1].second];

        int ref_diff = cov_ap1_i1.ref-cov_ap1_i.ref;
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
    return true;
}

bool NewEdgeCalculator::similarityCriterion(const AlignmentRecord & a1, const std::vector<AlignmentRecord::mapValue> &cov_ap1, const AlignmentRecord & a2, const std::vector<AlignmentRecord::mapValue> &cov_ap2, double probM, double prob0, int tc, int cc) const{

    //Threshold for probability that reads were sampled from same haplotype
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
    if (cc<=MIN_OVERLAP*std::min(cov_ap1.size(),cov_ap2.size())) return false;
    double prob = probM*prob0;
    int test = cc+tc;
    double potence = 1.0/test;
    double final_prob = std::pow(prob,potence);
    //cout << "Final prob: " << final_prob << endl;
    return final_prob >= cutoff;
}

bool NewEdgeCalculator::edgeBetween(const AlignmentRecord & ap1, const AlignmentRecord & ap2) const{
    const auto& cov_ap1 = ap1.getCovmap();
    const auto& cov_ap2 = ap2.getCovmap();
    //if (ap1.getName().find("Clique_1218") != string::npos && ap2.getName().find("Clique_1219") != string::npos){
    //    int k = 0;
    //}
    /*if(ap1.getName() == "MISEQ-02:83:000000000-A9WYY:1:1101:19916:16016" && ap2.getName()== "MISEQ-02:83:000000000-A9WYY:1:2104:10608:13344"){
            int k = 0;
        }*/
        /*if ((ap1.getName().find("Clique_1174") != string::npos && ap2.getName().find("Clique_1173") != string::npos)){
            cout << ap1.getName() << " " << ap2.getName() << endl;
            cout << ap1.isSingleEnd() << " " << ap2.isSingleEnd() << endl;
            cout << ap1.getSequence1().toString() << endl;
            cout << ap1.getStart1() << " " << ap1.getEnd1() << endl;
            cout << cov_ap1.front().ref << " " << cov_ap1.back().ref << endl;
            cout << ap2.getSequence1().toString() << endl;
            cout << ap2.getStart1() << " " << ap2.getEnd1() << endl;
            cout << cov_ap2.front().ref << " " << cov_ap2.back().ref << endl;
        }*/

    //tail position and common position counter
    int tc = 0;
    int cc = 0;
    double probM = 1.0;
    double prob0 = 1.0;
    int pos1 = 0;
    int pos2 = 0;
    int equalBase = 0;
    std::vector<std::pair<int,int>> aub;

    //iterating the covered positions and computing ProbM and Prob0 simultaneously
    while(pos1<cov_ap1.size() && pos2<cov_ap2.size()){
        auto& v1 = cov_ap1[pos1];
        auto& v2 = cov_ap2[pos2];
        if (v1.ref == v2.ref){
            if(v1.base == v2.base){
                equalBase++;
            }
            calculateProbM(v1,v2,probM);
            aub.push_back(std::make_pair(pos1,pos2));
            pos1++;
            pos2++;
            cc++;
        }
        else if (v1.ref < v2.ref){
            calculateProb0(v1,prob0);
            pos1++;
            tc++;

        }
        else if(v1.ref > v2.ref){
            calculateProb0(v2,prob0);
            pos2++;
            tc++;
        }
    }
    if(equalBase == cov_ap1.size() || equalBase == cov_ap2.size()){
        return true;
    }
    //add remaining entrys but only if they are both on the same read in the case of paired ends
    while(pos1<cov_ap1.size()){
        auto& v1 = cov_ap1[pos1];
        if (v1.read == cov_ap2.back().read){
            calculateProb0(v1,prob0);
            pos1++;
            tc++;
        } else break;
    }
    while(pos2<cov_ap2.size()){
        auto& v2 = cov_ap2[pos2];
        if (v2.read == cov_ap1.back().read){
            calculateProb0(v2,prob0);
            pos2++;
            tc++;
        } else break;
    }

    if (cc == 0 || (!checkGaps(cov_ap1, cov_ap2, aub))){
        return false;
    }
    //cout << "Gaps are compatible" << endl;
    //std::vector<int> tail = tailPositions(cov_ap1,cov_ap2);
    //cout << "Tail positions computed" << endl;
    return similarityCriterion(ap1, cov_ap1, ap2, cov_ap2, probM, prob0, tc, cc);
}

void NewEdgeCalculator::getPartnerLengthRange(const AlignmentRecord& ap, unsigned int* min, unsigned int* max) const {
    assert(false);
}




