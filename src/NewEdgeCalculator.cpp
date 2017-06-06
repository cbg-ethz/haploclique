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
 namespace{
     float phredProb(const char& qual){
         return std::pow(10, (double)(-qual)/10.0);
     }
     std::array<float, 127> compute_error_probs(){
         std::array<float, 127> result;
         for (int i = 33; i < result.size(); i++){
             result[i] = phredProb(i-33);
         }
         return result;
     }
     std::array<std::array<double,127>,127> compute_probM_eqBase_values(){
        std::array<std::array<double, 127>, 127> result;
        std::array<float, 127> error_probs = compute_error_probs();
        for (int i = 33; i < result.size(); i++){
            for(int j = 33; j < result.size(); j++){
                result[i][j] = (1.0-error_probs[i])*(1.0-error_probs[j])+(error_probs[i])*(error_probs[j]/3);
            }
        }
        return result;

     }
     std::array<std::array<double,127>,127> compute_probM_unEqBase_values(){
        std::array<std::array<double, 127>, 127> result;
        std::array<float, 127> error_probs = compute_error_probs();
        for (int i = 33; i < result.size(); i++){
            for(int j = 33; j < result.size(); j++){
                result[i][j] = (1.0-error_probs[i])*(error_probs[j]/3)+(1.0-error_probs[j])*(error_probs[i]/3)+2*(error_probs[i]/3)*(error_probs[j]/3);
            }
        }
        return result;
     }

     std::array<std::array<double,127>,127> probM_eqBase_values = compute_probM_eqBase_values();
     std::array<std::array<double,127>,127> probM_unEqBase_values = compute_probM_unEqBase_values();
 }

NewEdgeCalculator::NewEdgeCalculator(double Q, double edge_quasi_cutoff, double overlap, bool frameshift_merge, unordered_map<int, double>& simpson_map, double edge_quasi_cutoff_single, double overlap_single, double edge_quasi_cutoff_mixed, unsigned int maxPosition, bool noProb0) {
    this->Q = Q;
    this->EDGE_QUASI_CUTOFF = edge_quasi_cutoff;
    this->EDGE_QUASI_CUTOFF_SINGLE = edge_quasi_cutoff_single;
    this->EDGE_QUASI_CUTOFF_MIXED = edge_quasi_cutoff_mixed;
    this->MIN_OVERLAP_CLIQUES = overlap;
    this->MIN_OVERLAP_SINGLE = overlap_single;
    this->FRAMESHIFT_MERGE = frameshift_merge;
    this->SIMPSON_MAP.resize(maxPosition+1);
    for(auto& k_v : simpson_map){
        this->SIMPSON_MAP[k_v.first] = k_v.second;
    }
    this->NOPROB0 = noProb0;
}

NewEdgeCalculator::~NewEdgeCalculator() {
}

void NewEdgeCalculator::calculateProbM(const AlignmentRecord::mapValue & val1,const AlignmentRecord::mapValue & val2, double& res) const{
    double sum = 0.0;
    if(val1.base==val2.base){
        //sum+=(1.0-val1.prob)*(1.0-val2.prob)+(val1.prob)*(val2.prob/3);
        sum+=probM_eqBase_values[val1.qual][val2.qual];
    } else  {
        sum+=probM_unEqBase_values[val1.qual][val2.qual];
    }
    res*=sum;
}



void NewEdgeCalculator::calculateProb0(const AlignmentRecord::mapValue & val1, double& res) const{
    //find positions which are non common positions and look up probability in empirical allele frequency distribution
    //const auto& k = this->SIMPSON_MAP.find(val1.ref);
    auto d = this->SIMPSON_MAP[val1.ref];
    if (d!=0){
        res+=d;
    //if (k != this->SIMPSON_MAP.end()){
    //    res *= k->second;
    } else {
        res +=std::log10(0.25);
    }

}

//TO DO: find out whether gaps / insertions are compatible
bool NewEdgeCalculator::checkGaps(const std::vector<AlignmentRecord::mapValue> & cov_ap1,const std::vector<AlignmentRecord::mapValue> & cov_ap2,const std::vector<std::pair<int,int>> & aub) const{
    //int num_inserts=0;
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
            /*if(num_inserts == 0) num_inserts+= pos_diff1-pos_diff2;
            else*/ return false;
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
        cutoff = std::log10(this->EDGE_QUASI_CUTOFF);
    } else if (a1.getName().find("Clique") != string::npos || a2.getName().find("Clique") != string::npos) {
        cutoff = std::log10(this->EDGE_QUASI_CUTOFF_MIXED);
    } else {
        cutoff = std::log10(this->EDGE_QUASI_CUTOFF_SINGLE);
    }

    //Threshold for Overlap of Read Alignments
    double MIN_OVERLAP = 0;
    if (a1.getName().find("Clique") != string::npos && a2.getName().find("Clique") != string::npos) {
        MIN_OVERLAP = MIN_OVERLAP_CLIQUES;
    } else {
        MIN_OVERLAP = MIN_OVERLAP_SINGLE;
    }
    //double ratio = (float)cc/tc;
    if (cc<=MIN_OVERLAP*std::min(cov_ap1.size(),cov_ap2.size())) return false;
    double prob;
    double potence;
    /*if((double)cc/tc<=0.15 && probM >= 0.8 && MIN_OVERLAP==0.1){
        prob = probM;
        potence = 1.0/(cc);
    } else {*/
    //TEST to see what happens if prob0 is neglected
    if(this->NOPROB0){
        prob = std::log10(probM);
        potence = 1.0/cc;
    } else {
        prob = std::log10(probM) + prob0;
        potence = 1.0/(cc+tc);
    }
    //}
    double final_prob = prob*potence;
    //cout << "Final prob: " << final_prob << endl;
    return final_prob >= cutoff;
}

void NewEdgeCalculator::iterateCovAp(bool pe1, unsigned int& pos2, double& probM, unsigned int& pos1, unsigned int& equalBase, bool pe2, const std::vector<AlignmentRecord::mapValue>& cov_ap2, int& tc, double& prob0, const AlignmentRecord& ap2, const std::vector<AlignmentRecord::mapValue>& cov_ap1, const AlignmentRecord& ap1, int& cc) const
{
    while(pos1<cov_ap1.size() && pos2<cov_ap2.size()){
        auto& v1 = cov_ap1[pos1];
        auto& v2 = cov_ap2[pos2];
        if (v1.ref == v2.ref){
            if(v1.base == v2.base){
                equalBase++;
            }
            calculateProbM(v1,v2,probM);
            aub.emplace_back(pos1,pos2);
            pos1++;
            pos2++;
            cc++;
        }
        else if (v1.ref < v2.ref){
           if (!pe1){
               calculateProb0(v1,prob0);
               tc++;
               pos1++;
           } else {
                if (ap2.getStart1() > ap1.getEnd1() && v1.ref < ap1.getStart2()) pos1++;
                else{
                    calculateProb0(v1,prob0);
                    tc++;
                    pos1++;
                }
           } /*else if (!pe1 && pe2){
                calculateProb0(v1,prob0);
                tc++;
                pos1++;
           }*//*else if (pe1 && pe2){
                if(ap1.getEnd1() < ap2.getStart1() && v1.ref < ap1.getStart2()) pos1++;
                else{
                    calculateProb0(v1,prob0);
                    tc++;
                    pos1++;
                }
           }*/
        }
        else {
               if (!pe2){
                   calculateProb0(v2,prob0);
                   tc++;
                   pos2++;
               } else {
                   if (ap1.getStart1() > ap2.getEnd1() && v2.ref < ap2.getStart2()) pos2++;
                   else{
                       calculateProb0(v2,prob0);
                       tc++;
                       pos2++;
                   }
               } /*else if (pe1 && pe2){
                    if(ap2.getEnd1() < ap1.getStart1() && v1.ref < ap2.getStart2()) pos2++;
                    else{
                        calculateProb0(v2,prob0);
                        tc++;
                        pos2++;
                    }
               }*/
        }
    }
}

void NewEdgeCalculator::iterateRemainingCovAp(int& tc, bool pe1, const std::vector<AlignmentRecord::mapValue>& cov_ap2, bool pe2, double& prob0, const AlignmentRecord& ap1, const AlignmentRecord& ap2, const std::vector<AlignmentRecord::mapValue>& cov_ap1, unsigned int& pos1) const
{
    while(pos1 < cov_ap1.size()){
        auto& v1 = cov_ap1[pos1];
        if(pe1 && pe2){
            if(ap1.getStart2() > ap2.getEnd2() && v1.ref >= ap1.getStart2()) break;
            else {
                calculateProb0(v1,prob0);
                pos1++;
                tc++;
            }
        }
        else if (v1.read == cov_ap2.back().read){
            calculateProb0(v1,prob0);
            pos1++;
            tc++;
        } else break;
    }
}

void NewEdgeCalculator::iterateRemainingCovAp2(const AlignmentRecord& ap2, const std::vector<AlignmentRecord::mapValue>& cov_ap2, bool pe1, const AlignmentRecord& ap1, int& tc, double& prob0, const std::vector<AlignmentRecord::mapValue>& cov_ap1, bool pe2, unsigned int& pos2) const
{
    while(pos2 < cov_ap2.size()){
        auto& v2 = cov_ap2[pos2];
        if(pe1 && pe2){
            if(ap2.getStart2() > ap1.getEnd2() && v2.ref >= ap2.getStart2()) break;
            else {
                calculateProb0(v2,prob0);
                pos2++;
                tc++;
            }
        }
       else if (v2.read == cov_ap1.back().read){
            calculateProb0(v2,prob0);
            pos2++;
            tc++;
        } else break;
    }
}

bool NewEdgeCalculator::edgeBetween(const AlignmentRecord & ap1, const AlignmentRecord & ap2) const{
    const auto& cov_ap1 = ap1.getCovmap();
    const auto& cov_ap2 = ap2.getCovmap();
    // -----    ------
    //                   -----     ----
    if(cov_ap1.back().ref < cov_ap2[0].ref  || cov_ap2.back().ref < cov_ap1[0].ref){
        return false;
    }
    //if (ap1.getName().find("Clique_1218") != string::npos && ap2.getName().find("Clique_1219") != string::npos){
    //    int k = 0;
    //}
    //if(ap1.getName() == "MISEQ-02:63:000000000-A8NP3:1:2102:21435:26782" || ap2.getName() == "MISEQ-02:63:000000000-A8NP3:1:2102:21435:26782"){
    //        int k = 0;
    //}
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
    double prob0 = 0.0;
    unsigned int pos1 = 0;
    unsigned int pos2 = 0;
    unsigned int equalBase = 0;
    bool pe1 = ap1.isPairedEnd();
    bool pe2 = ap2.isPairedEnd();
    aub.clear();

    // special cases of paired end reads for which no edge is allowed
    if(pe1 && pe2){
        unsigned int e11= ap1.getEnd1();
        unsigned int e12= ap1.getEnd2();
        unsigned int s11= ap1.getStart1();
        unsigned int s12= ap1.getStart2();
        unsigned int e21= ap2.getEnd1();
        unsigned int e22= ap2.getEnd2();
        unsigned int s21= ap2.getStart1();
        unsigned int s22= ap2.getStart2();
        //--------   ---------
        //           ---------   -------
        if(e11 < s21 && e12 < s22){
            return false;
        } //         ---------    ---------
        //-------    ---------
        else if(e21 < s11 && e22 < s12){
            return false;
        } // ------------     ----------
        //   -----------                  ---------
        else if(e21 < s12 && e12 < s22){
            return false;
        }//-----------            ----------
        //-----------  ----------
        else if(e11 < s22 && e22 < s12){
            return false;
        }//          --------- ---------
        //----------           ---------
        else if (e21 < s11 && e11 < s22){
            return false;
        }//--------             --------
        //           --------   --------
        else if(e11 < s21 && e21 < s12){
            return false;
        }
    }
    //if(ap1.getName() == "Clique_2268" || ap2.getName() == "Clique_2268"){
    //            int k = 0;
    //}
    //iterating the covered positions and computing ProbM and Prob0 simultaneously
    iterateCovAp(pe1, pos2, probM, pos1, equalBase, pe2, cov_ap2, tc, prob0, ap2, cov_ap1, ap1, cc);

    //given the case that one read is contained in another
    if(equalBase == cov_ap1.size() || equalBase == cov_ap2.size()){
        return true;
    }
    //add remaining entries, but only if they are both in the same read in the case of paired ends
    iterateRemainingCovAp(tc, pe1, cov_ap2, pe2, prob0, ap1, ap2, cov_ap1, pos1);
    iterateRemainingCovAp2(ap2, cov_ap2, pe1, ap1, tc, prob0, cov_ap1, pe2, pos2);
    //no edge if there are no common positions or incompatible gaps
    if (cc == 0 || (!checkGaps(cov_ap1, cov_ap2, aub))){
        return false;
    }
    /*TEST to see how often Prob0 is set to zero
    if(prob0 == 0.0){
        cout << "Prob0 is zero" << endl;
    }*/
    return similarityCriterion(ap1, cov_ap1, ap2, cov_ap2, probM, prob0, tc, cc);
}

void NewEdgeCalculator::getPartnerLengthRange(const AlignmentRecord& ap, unsigned int* min, unsigned int* max) const {
    assert(false);
}

void NewEdgeCalculator::setOverlapCliques(double d){
    this->MIN_OVERLAP_CLIQUES=d;
}




