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

int computeOffset(const std::vector<char>& cigar);
void computeSOffset(const std::vector<char>& cigar,int& c, int& q);
int computeRevOffset(const std::vector<char>& cigar);
int computeRevSOffset(const std::vector<char>& cigar);

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

void NewEdgeCalculator::computeProbM(const char& base1, const char& qual1, const char& base2, const char& qual2, double& res ) const{
    
    double sum = 0.0;
    if(base1 == base2){
        sum+=probM_eqBase_values[qual1][qual2];
    } else {
        sum+=probM_unEqBase_values[qual1][qual2];
    }
    res*=sum;
}

void NewEdgeCalculator::computeProb0(int ref_pos, double& res ) const{
    
    auto d = this->SIMPSON_MAP[ref_pos];
    
    if (d!=0){
        res+=d;
    } else {
        res +=std::log10(0.25);
    }
    
}

bool NewEdgeCalculator::overlapCheckgap(const AlignmentRecord& ap1, const AlignmentRecord& ap2, int& c_pos1, int& c_pos2, int& q_pos1, int& q_pos2, int& ref_pos, double& probM, int& cc, int i, int j) const{
    
    bool flag_gap = false;
    char c1, c2;
    
    const ShortDnaSequence* s1,* s2 = 0;
    if (i == 1){
        c1 = ap1.getCigar1Unrolled()[c_pos1];
        s1 = &ap1.getSequence1();
        if(j == 1){
            c2 = ap2.getCigar1Unrolled()[c_pos2];
            s2 = &ap2.getSequence1();
        } else {
            c2 = ap2.getCigar2Unrolled()[c_pos2];
            s2 = &ap2.getSequence2();
        }
    } else {
        c1 = ap1.getCigar2Unrolled()[c_pos1];
        s1 = &ap1.getSequence2();
        if(j == 1){
            c2 = ap2.getCigar1Unrolled()[c_pos2];
            s2 = &ap2.getSequence1();
        } else {
            c2 = ap2.getCigar2Unrolled()[c_pos2];
            s2 = &ap2.getSequence2();
        }
    }
    
    if((c1 == 'M' && c2 == 'M') || (c1 == 'S' && c2 == 'S') || (c1 == 'I' && c2 == 'I')){
        if (c1 != 'S'){
            computeProbM((*s1)[q_pos1],s1->qualityChar(q_pos1),(*s2)[q_pos2],s2->qualityChar(q_pos2), probM);
            cc++;
        }
        if (c1 != 'I') ref_pos++;
        q_pos1++;
        q_pos2++;
        c_pos1++;
        c_pos2++;
    } else if ((c1 == 'D' && c2 == 'D') || (c1 == 'H' && c2 == 'H') || (c1 == 'D' && c2 == 'H') || (c1 == 'H' && c2 == 'D') || (c1 == 'D' && c2 == 'S') || (c1 == 'S' && c2 == 'D')){
        c_pos1++;
        c_pos2++;
        ref_pos++;
        if (c1 == 'D' || c2 == 'D'){
        }
        if (c1 == 'S'){
            q_pos1++;
        } else if(c2 == 'S'){
            q_pos2++;
        }
    }
    else if ((c1 == 'M' && (c2 == 'H' || c2 == 'S')) || ((c1 == 'H' || c1 == 'S') && c2 == 'M') || (c1 == 'S' && c2 == 'H') || (c1 == 'H' && c2 == 'S')) {
        if (c1 == 'M'){
            ref_pos++;
            c_pos1++;
            c_pos2++;
            q_pos1++;
            if (c2 == 'S') q_pos2++;
        } else if (c2 == 'M'){
            ref_pos++;
            c_pos1++;
            c_pos2++;
            q_pos2++;
            if (c1 == 'S') q_pos1++;
        } else if (c1 == 'S'){
            ref_pos++;
            c_pos1++;
            c_pos2++;
            q_pos1++;
        } else {
            ref_pos++;
            c_pos1++;
            c_pos2++;
            q_pos2++;
        }
    } else if ((c1 == 'I' && (c2 == 'H' || c2 == 'S'))|| ((c1 == 'H' || c1 == 'S') && c2 == 'I')){
        if(c1 == 'I'){
            c_pos1++;
            q_pos1++;
        } else {
            c_pos2++;
            q_pos2++;
        }
    } else{
        flag_gap = true;
        return flag_gap;
    }
    return flag_gap;
}

void NewEdgeCalculator::noOverlapCheckgap(const AlignmentRecord& ap, int& c_pos, int& q_pos, int& ref_pos, double& prob0, int& tc, int i) const{
    char c;
    if (i == 1){
        c = ap.getCigar1Unrolled()[c_pos];
    } else {
        c = ap.getCigar2Unrolled()[c_pos];
    }
    if (c == 'H'){
        ref_pos++;   //reference position
        c_pos++;     //cigar position (index in cigar array)
    } else if (c == 'I') {
        computeProb0(ref_pos, prob0);
        tc++;
        q_pos++;  // index in quality char array
        c_pos++;
    } else if (c == 'D') {
        ref_pos++;
        c_pos++;
    } else if (c == 'S'){
        ref_pos++;
        q_pos++;
        c_pos++;
    } else if (c == 'M'){
        computeProb0(ref_pos, prob0);
        tc++;
        ref_pos++;
        q_pos++;
        c_pos++;
    } else {
        cout << ap.getName() << endl;
        cout << "Cigar string contains inappropriate character: " << c << endl;
        assert(0);
    }
}

bool NewEdgeCalculator::checkGapsSingle(const AlignmentRecord& ap1, const AlignmentRecord& ap2, double& probM, double& prob0, int& cc, int& tc , int i, int j) const{
    
    int offset_f1, offset_f2, offset_b1, offset_b2, ref_s_pos1, ref_e_pos1, ref_s_pos2, ref_e_pos2 = 0;
    //i and j determine which cigar strings / sequences are considered
    if (i == 1){
        //get starting position and ending position according to ref position, paying attention to clipped bases
        //updated ref position including clips
        
        std::vector<char> cigar1 = ap1.getCigar1Unrolled();
        offset_f1 = computeOffset(cigar1);
        offset_b1 = computeRevOffset(cigar1);
        ref_e_pos1 = ap1.getEnd1() + offset_b1;
        ref_s_pos1 = ap1.getStart1() - offset_f1;
        if(j == 1){
            std::vector<char> cigar = ap2.getCigar1Unrolled();
            offset_f2 = computeOffset(cigar);
            offset_b2 = computeRevOffset(cigar);
            ref_s_pos2 = ap2.getStart1() - offset_f2;
            ref_e_pos2 = ap2.getEnd1() + offset_b2;
        } else {
            std::vector<char> cigar = ap2.getCigar2Unrolled();
            offset_f2 = computeOffset(cigar);
            offset_b2 = computeRevOffset(cigar);
            ref_s_pos2 = ap2.getStart2() - offset_f2;
            ref_e_pos2 = ap2.getEnd2() + offset_b2;
        }
    } else {
        //get starting position and ending position according to ref position, paying attention to clipped bases
        //updated ref position including clips
        
        std::vector<char> cigar2 = ap1.getCigar2Unrolled();
        
        offset_f1 = computeOffset(cigar2);
        offset_b1 = computeRevOffset(cigar2);
        ref_s_pos1 = ap1.getStart2()-offset_f1;
        ref_e_pos1 = ap1.getEnd2()+offset_b1;
        
        if(j == 1){
            std::vector<char> cigar = ap2.getCigar1Unrolled();
            offset_f2 = computeOffset(cigar);
            offset_b2 = computeRevOffset(cigar);
            ref_s_pos2 = ap2.getStart1()-offset_f2;
            ref_e_pos2 = ap2.getEnd1()+offset_b2;
        } else {
            std::vector<char> cigar = ap2.getCigar2Unrolled();
            offset_f2 = computeOffset(cigar);
            offset_b2 = computeRevOffset(cigar);
            ref_s_pos2 = ap2.getStart2()-offset_f2;
            ref_e_pos2 = ap2.getEnd2()+offset_b2;
        }
    }
    //position in query!! or quality sequences // phred scores
    int q_pos1 = 0;
    int q_pos2 = 0;
    //position in unrolled cigar vectors
    int c_pos1 = 0;
    int c_pos2 = 0;
    
    bool flag_gap = false;
    //4 cases of different overlaps
    //------------
    //     ------------
    if(ref_s_pos1 <= ref_s_pos2 && ref_e_pos1 <= ref_e_pos2){
        while(ref_s_pos1<ref_s_pos2){
            noOverlapCheckgap(ap1,c_pos1,q_pos1,ref_s_pos1, prob0 , tc, i);
        }
        while(ref_s_pos1<=ref_e_pos1){
            flag_gap = overlapCheckgap(ap1, ap2,c_pos1,c_pos2,q_pos1,q_pos2,ref_s_pos1,  probM , cc, i,j);
            if (flag_gap) {
                return flag_gap;
            }
        }
        while(ref_s_pos1<=ref_e_pos2){ // this should not be needed
            noOverlapCheckgap(ap2, c_pos2, q_pos2, ref_s_pos1, prob0 , tc, j);
        }
    }//------------------------------
    //           ----------
    else if (ref_s_pos1 <= ref_s_pos2 && ref_e_pos1 >= ref_e_pos2){
        while(ref_s_pos1<ref_s_pos2){
            noOverlapCheckgap(ap1,c_pos1,q_pos1,ref_s_pos1, prob0 ,  tc, i);
        }
        while(ref_s_pos1<=ref_e_pos2){
            flag_gap = overlapCheckgap(ap1, ap2,c_pos1,c_pos2,q_pos1,q_pos2,ref_s_pos1, probM , cc, i,j);
            if (flag_gap) {
                return flag_gap;
            }
        }
        while(ref_s_pos1<=ref_e_pos1){
            noOverlapCheckgap(ap1,c_pos1,q_pos1,ref_s_pos1, prob0 ,  tc, i);
        }
        //         ----------
        //--------------------------
    } else if (ref_s_pos1 >= ref_s_pos2 && ref_e_pos1 <= ref_e_pos2){
        while(ref_s_pos2<ref_s_pos1){
            noOverlapCheckgap(ap2, c_pos2, q_pos2, ref_s_pos2, prob0 , tc,  j);
        }
        while(ref_s_pos2<=ref_e_pos1){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_pos1,c_pos2,q_pos1,q_pos2,ref_s_pos2, probM , cc, i,j);
            if (flag_gap) {
                return flag_gap;
            }
        }
        while(ref_s_pos2<=ref_e_pos2){
            noOverlapCheckgap(ap2, c_pos2, q_pos2, ref_s_pos2,prob0 , tc, j);
        }
        //           --------------------
        //---------------------
    } else if (ref_s_pos1 >= ref_s_pos2 && ref_e_pos1 >= ref_e_pos2) {
        while(ref_s_pos2<ref_s_pos1){
            noOverlapCheckgap(ap2, c_pos2, q_pos2, ref_s_pos2, prob0 , tc,  j);
        }
        while(ref_s_pos2<=ref_e_pos2){
            flag_gap = overlapCheckgap(ap1, ap2, c_pos1,c_pos2,q_pos1,q_pos2,ref_s_pos2, probM ,  cc,i,j);
            if (flag_gap) {
                return flag_gap;
            }
        }
        while(ref_s_pos2<=ref_e_pos1){
            noOverlapCheckgap(ap1, c_pos1,q_pos1,ref_s_pos2, prob0 , tc,  i);
        }
    }else{
        cout << " checkGapsSingle " << ap1.getName() << " , " << ap2.getName() << endl << flush;
        assert(0);
    }
    
    return flag_gap;
}

bool NewEdgeCalculator::checkGapsPaired(const AlignmentRecord& ap1, const AlignmentRecord& ap2, double& probM, double& prob0, int& cc, int& tc) const{
    
    bool flag_gap = false;
    //get starting position and ending position according to ref position, paying attention to clipped bases
    
    int offset_f1_c1 = computeOffset(ap1.getCigar1Unrolled());
    int offset_f1_c2 = computeOffset(ap1.getCigar2Unrolled());
    int offset_f2_c1 = computeOffset(ap2.getCigar1Unrolled());
    int offset_f2_c2 = computeOffset(ap2.getCigar2Unrolled());
    
    int offset_b1_c1 = computeRevOffset(ap1.getCigar1Unrolled());
    int offset_b1_c2 = computeRevOffset(ap1.getCigar2Unrolled());
    int offset_b2_c1 = computeRevOffset(ap2.getCigar1Unrolled());
    int offset_b2_c2 = computeRevOffset(ap2.getCigar2Unrolled());
    //updated ref position including clips
    
    int ref_s_pos1_c1 = ap1.getStart1() - offset_f1_c1;
    int ref_e_pos1_c1 = ap1.getEnd1() + offset_b1_c1;
    int ref_s_pos1_c2 = ap1.getStart2() - offset_f1_c2;
    int ref_e_pos1_c2 = ap1.getEnd2() + offset_b1_c2;
    
    int ref_s_pos2_c1 = ap2.getStart1() - offset_f2_c1;
    int ref_e_pos2_c1 = ap2.getEnd1() + offset_b2_c1;
    int ref_s_pos2_c2 = ap2.getStart2() - offset_f2_c2;
    int ref_e_pos2_c2 = ap2.getEnd2() + offset_b2_c2;
    
    //position in query sequences // phred scores
    int q_c1_pos1 = 0;
    int q_c2_pos1 = 0;
    int q_c1_pos2 = 0;
    int q_c2_pos2 = 0;
    //position in unrolled cigar vectors
    int c_c1_pos1 = 0;
    int c_c2_pos1 = 0;
    int c_c1_pos2 = 0;
    int c_c2_pos2 = 0;
    //int i;
    // --------    |  -----------    <-ap1
    //   --------  |     ----------
    if(ap1.getEnd1() < ap2.getStart2() && ap1.getStart2() > ap2.getEnd1()){
        flag_gap = checkGapsSingle(ap1, ap2, probM, prob0,cc, tc, 1,1);
        if (flag_gap) return flag_gap;
        flag_gap = checkGapsSingle(ap1, ap2, probM, prob0, cc, tc, 2,2);
        if (flag_gap) return flag_gap;
        
    }//----------    ------------   <-ap1
    //    ---------------   -----------
    else if(ref_s_pos1_c1 <= ref_s_pos2_c1 && ap1.getEnd1() >= ap2.getStart1() && ap1.getStart2() <= ap2.getEnd1() && ap1.getEnd2() >= ap2.getStart2()){
        while(ref_s_pos1_c1 < ref_s_pos2_c1){
            noOverlapCheckgap(ap1,c_c1_pos1,q_c1_pos1,ref_s_pos1_c1, prob0, tc,  1);
        }
        while(ref_s_pos1_c1 <= ap1.getEnd1()){
            flag_gap = overlapCheckgap(ap1, ap2,c_c1_pos1,c_c1_pos2,q_c1_pos1,q_c1_pos2,ref_s_pos1_c1, probM,  cc,1,1);
            if (flag_gap)
                return flag_gap;
        }
        while(ref_s_pos1_c1< ap1.getStart2()){
            noOverlapCheckgap(ap2,c_c1_pos2,q_c1_pos2,ref_s_pos1_c1, prob0, tc,  1);
        }
        computeSOffset(ap1.getCigar2Unrolled() ,c_c2_pos1,q_c2_pos1);
        while(ref_s_pos1_c1<= ap2.getEnd1()){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c2_pos1,c_c1_pos2,q_c2_pos1,q_c1_pos2,ref_s_pos1_c1, probM, cc, 2,1);
            if (flag_gap)
                return flag_gap;
        }
        while(ref_s_pos1_c1<ap2.getStart2()){
            noOverlapCheckgap(ap1, c_c2_pos1,q_c2_pos1,ref_s_pos1_c1,prob0, tc,  2);
        }
        computeSOffset(ap2.getCigar2Unrolled(),c_c2_pos2,q_c2_pos2);
        while(ref_s_pos1_c1<=ref_e_pos2_c2 && ref_s_pos1_c1 <= ref_e_pos1_c2){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c2_pos1,c_c2_pos2,q_c2_pos1,q_c2_pos2,ref_s_pos1_c1, probM, cc, 2,2);
            if (flag_gap)
                return flag_gap;
        }
        if(ref_s_pos1_c1-1 == ref_e_pos2_c2){
            while(ref_s_pos1_c1<=ref_e_pos1_c2){
                noOverlapCheckgap(ap1,c_c2_pos1,q_c2_pos1,ref_s_pos1_c1, prob0, tc,  2);
            }
        } else if (ref_s_pos1_c1-1 == ref_e_pos1_c2){
            while(ref_s_pos1_c1<=ref_e_pos2_c2){
                noOverlapCheckgap(ap2,c_c2_pos2,q_c2_pos2,ref_s_pos1_c1, prob0,  tc, 2);
            }
        }
    } //-----------       ----------- <-ap1
    //    -----------------               -----------
    else if(ref_s_pos1_c1 <= ref_s_pos2_c1 && ap1.getEnd1() >=ap2.getStart1() && ap1.getStart2() <= ap2.getEnd1() && ap1.getEnd2() < ap2.getStart2()){
        
        while(ref_s_pos1_c1 < ref_s_pos2_c1){
            noOverlapCheckgap(ap1,c_c1_pos1,q_c1_pos1,ref_s_pos1_c1, prob0, tc,  1);
        }
        while(ref_s_pos1_c1 <= ap1.getEnd1() ){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c1_pos1,c_c1_pos2,q_c1_pos1,q_c1_pos2,ref_s_pos1_c1, probM, cc ,1,1);
            if (flag_gap)
                return flag_gap;
        }
        while(ref_s_pos1_c1 < ap1.getStart2()){
            noOverlapCheckgap(ap2 ,c_c1_pos2,q_c1_pos2,ref_s_pos1_c1, prob0, tc,  1);
        }
        computeSOffset(ap1.getCigar2Unrolled() ,c_c2_pos1,q_c2_pos1);
        while(ref_s_pos1_c1<=ref_e_pos2_c1 && ref_s_pos1_c1<=ref_e_pos1_c2){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c2_pos1,c_c1_pos2,q_c2_pos1,q_c1_pos2,ref_s_pos1_c1, probM, cc , 2,1);
            if (flag_gap)
                return flag_gap;
        }
        if(ref_s_pos1_c1-1 == ref_e_pos2_c1 ){
            while(ref_s_pos1_c1<=ref_e_pos1_c2){
                noOverlapCheckgap(ap1 ,c_c2_pos1,q_c2_pos1,ref_s_pos1_c1, prob0 , tc, 2);
            }
        } else if(ref_s_pos1_c1-1 == ref_e_pos1_c2){
            while(ref_s_pos1_c1<=ref_e_pos2_c1){
                noOverlapCheckgap( ap2 ,c_c1_pos2,q_c1_pos2,ref_s_pos1_c1, prob0 ,  tc, 1);
            }
        }
        
    } //----------        ------------  <- ap1
    //                --------    ----------
    else if(ref_s_pos2_c1 <= ref_s_pos1_c2 && ap1.getEnd1() < ap2.getStart1() && ap1.getStart2() <= ap2.getEnd1() && ap1.getEnd2() >= ap2.getStart2()){
        
        while(ref_s_pos2_c1 < ref_s_pos1_c2){
            noOverlapCheckgap(ap2,c_c1_pos2,q_c1_pos2,ref_s_pos2_c1, prob0, tc,  1);
        }
        while(ref_s_pos2_c1 <= ap2.getEnd1()){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c2_pos1,c_c1_pos2,q_c2_pos1,q_c1_pos2,ref_s_pos2_c1, probM ,  cc ,2,1);
            if (flag_gap)
                return flag_gap;
        }
        while(ref_s_pos2_c1 < ap2.getStart2()){
            noOverlapCheckgap( ap1,c_c2_pos1,q_c2_pos1,ref_s_pos2_c1, prob0 , tc,  2);
        }
        computeSOffset(ap2.getCigar2Unrolled(),c_c2_pos2,q_c2_pos2);
        while(ref_s_pos2_c1 <= ref_e_pos2_c2 && ref_s_pos2_c1 <= ref_e_pos1_c2){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c2_pos1,c_c2_pos2,q_c2_pos1,q_c2_pos2,ref_s_pos2_c1, probM,  cc ,2,2);
            if (flag_gap)
                return flag_gap;
        }
        if(ref_s_pos2_c1-1 == ref_e_pos2_c2){
            while(ref_s_pos2_c1<=ref_e_pos1_c2){
                noOverlapCheckgap(ap1,c_c2_pos1,q_c2_pos1,ref_s_pos2_c1, prob0 , tc, 2);
            }
        } else if(ref_s_pos2_c1-1 == ref_e_pos1_c2){
            while(ref_s_pos2_c1<=ref_e_pos2_c2){
                noOverlapCheckgap( ap2 ,c_c2_pos2,q_c2_pos2,ref_s_pos2_c1, prob0, tc,  2);
            }
        }
    } //--------      --------- <-ap1
    //                --  -----------
    else if(ref_s_pos1_c2<=ref_s_pos2_c1 && ap1.getEnd1() < ap2.getStart1() && ap1.getStart2() <= ap2.getEnd1() && ap1.getEnd2() >= ap2.getStart2() ){
        
        while(ref_s_pos1_c2<ref_s_pos2_c1){
            noOverlapCheckgap(ap1,c_c2_pos1,q_c2_pos1,ref_s_pos1_c2, prob0, tc,  2);
        }
        while(ref_s_pos1_c2<=ap2.getEnd1()){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c2_pos1,c_c1_pos2,q_c2_pos1,q_c1_pos2,ref_s_pos1_c2, probM ,  cc ,2,1);
            if (flag_gap)
                return flag_gap;
        }
        while(ref_s_pos1_c2<ap2.getStart2()){
            noOverlapCheckgap( ap1 ,c_c2_pos1,q_c2_pos1,ref_s_pos1_c2, prob0 , tc,  2);
        }
        computeSOffset(ap2.getCigar2Unrolled(),c_c2_pos2,q_c2_pos2);
        while(ref_s_pos1_c2<=ref_e_pos1_c2 && ref_s_pos1_c2<=ref_e_pos2_c2){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c2_pos1,c_c2_pos2,q_c2_pos1,q_c2_pos2,ref_s_pos1_c2, probM, cc , 2,2);
            if (flag_gap)
                return flag_gap;
        }
        if(ref_s_pos1_c2-1==ref_e_pos1_c2){
            while(ref_s_pos1_c2<=ref_e_pos2_c2){
                noOverlapCheckgap( ap2,c_c2_pos2,q_c2_pos2,ref_s_pos1_c2, prob0, tc,  2);
            }
        } else if(ref_s_pos1_c2-1==ref_e_pos2_c2){
            while(ref_s_pos1_c2<=ref_e_pos1_c2){
                noOverlapCheckgap( ap1,c_c2_pos1,q_c2_pos1,ref_s_pos1_c2, prob0, tc,  2);
            }
        }
        
    } // -------------     ----------- <-ap1
    //   ----  ---------------
    else if(ref_s_pos1_c1 <= ref_s_pos2_c1 &&ap1.getStart1() <=ap2.getEnd1() && ap1.getEnd1() >= ap2.getStart2() && ap1.getStart2() <= ap2.getEnd2()){
        
        while(ref_s_pos1_c1<ref_s_pos2_c1){
            noOverlapCheckgap(ap1,c_c1_pos1,q_c1_pos1,ref_s_pos1_c1, prob0, tc,  1);
        }
        while(ref_s_pos1_c1<=ap2.getEnd1()){
            flag_gap = overlapCheckgap(ap1, ap2,c_c1_pos1,c_c1_pos2,q_c1_pos1,q_c1_pos2,ref_s_pos1_c1, probM, cc , 1,1);
            if (flag_gap)
                return flag_gap;
        }
        while(ref_s_pos1_c1<ap2.getStart2()){
            noOverlapCheckgap(ap1,c_c1_pos1,q_c1_pos1,ref_s_pos1_c1, prob0,  tc, 1);
        }
        computeSOffset(ap2.getCigar2Unrolled(),c_c2_pos2,q_c2_pos2);
        while(ref_s_pos1_c1<=ap1.getEnd1()){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c1_pos1,c_c2_pos2,q_c1_pos1,q_c2_pos2,ref_s_pos1_c1, probM , cc , 1,2);
            if (flag_gap)
                return flag_gap;
        }
        while(ref_s_pos1_c1<ap1.getStart2()){
            noOverlapCheckgap(ap2,c_c2_pos2,q_c2_pos2,ref_s_pos1_c1, prob0, tc,  2);
        }
        computeSOffset(ap1.getCigar2Unrolled(),c_c2_pos1,q_c2_pos1);
        while(ref_s_pos1_c1<=ref_e_pos2_c2 && ref_s_pos1_c1<=ref_e_pos1_c2){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c2_pos1,c_c2_pos2,q_c2_pos1,q_c2_pos2,ref_s_pos1_c1, probM ,  cc ,2,2);
            if (flag_gap)
                return flag_gap;
        }
        if(ref_s_pos1_c1-1==ref_e_pos2_c2){
            while(ref_s_pos1_c1<=ref_e_pos1_c2){
                noOverlapCheckgap(ap1,c_c2_pos1,q_c2_pos1,ref_s_pos1_c1, prob0 , tc,  2);
            }
        } else if(ref_s_pos1_c1-1==ref_e_pos1_c2){
            while(ref_s_pos1_c1<=ref_e_pos2_c2){
                noOverlapCheckgap(ap2,c_c2_pos2,q_c2_pos2,ref_s_pos1_c1, prob0,  tc, 2);
            }
        }
    } //----------         -----------  <-ap1
    //  ----   -------
    else if(ref_s_pos1_c1 <= ref_s_pos2_c1 && ap1.getStart2() > ap2.getEnd2() && ap1.getEnd1() >= ap2.getStart2() && ap1.getStart1() <= ap2.getEnd1()){
        
        while(ref_s_pos1_c1 < ref_s_pos2_c1){
            noOverlapCheckgap(ap1,c_c1_pos1,q_c1_pos1,ref_s_pos1_c1, prob0,  tc, 1);
        }
        while(ref_s_pos1_c1 <= ap2.getEnd1()){
            flag_gap = overlapCheckgap(ap1, ap2,c_c1_pos1,c_c1_pos2,q_c1_pos1,q_c1_pos2,ref_s_pos1_c1, probM,  cc ,1,1);
            if (flag_gap)
                return flag_gap;
        }
        while(ref_s_pos1_c1<ap2.getStart2()){
            noOverlapCheckgap(ap1,c_c1_pos1,q_c1_pos1,ref_s_pos1_c1, prob0, tc,  1);
        }
        computeSOffset(ap2.getCigar2Unrolled(),c_c2_pos2,q_c2_pos2);
        
        while(ref_s_pos1_c1<=ref_e_pos1_c1 && ref_s_pos1_c1<=ref_e_pos2_c2){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c1_pos1,c_c2_pos2,q_c2_pos1,q_c2_pos2,ref_s_pos1_c1, probM , cc , 1,2);
            if (flag_gap)
                return flag_gap;
        }
        
        if(ref_s_pos1_c1-1==ref_e_pos1_c1){
            while(ref_s_pos1_c1<=ref_e_pos2_c2){
                noOverlapCheckgap(ap2 ,c_c2_pos2,q_c2_pos2,ref_s_pos1_c1, prob0 , tc,  2);
            }
            
        } else if(ref_s_pos1_c1-1==ref_e_pos2_c2){
            while(ref_s_pos1_c1<=ref_e_pos1_c1){
                noOverlapCheckgap( ap1 ,c_c1_pos1,q_c1_pos1,ref_s_pos1_c1 , prob0 , tc, 1);
            }
        }
    } //   ------   -----------   <-ap1
    //----------------  ----------------
    else if(ref_s_pos2_c1 <= ref_s_pos1_c1 && ap1.getEnd1() >= ap2.getStart1() && ap1.getStart2() <= ap2.getEnd1() && ap1.getEnd2() >= ap2.getStart2()){
        while(ref_s_pos2_c1 < ref_s_pos1_c1){
            noOverlapCheckgap(ap2,c_c1_pos2,q_c1_pos2,ref_s_pos2_c1, prob0, tc,  1);
        }
        while(ref_s_pos2_c1 <=ap1.getEnd1()){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c1_pos1,c_c1_pos2,q_c1_pos1,q_c1_pos2,ref_s_pos2_c1, probM, cc , 1,1);
            if (flag_gap)
                return flag_gap;
        }
        while(ref_s_pos2_c1 < ap1.getStart2()){
            noOverlapCheckgap(ap2 ,c_c1_pos2,q_c1_pos2,ref_s_pos2_c1, prob0 ,  tc, 1 );
        }
        computeSOffset(ap1.getCigar2Unrolled(),c_c2_pos1,q_c2_pos1);
        while(ref_s_pos2_c1<=ap2.getEnd1()){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c2_pos1,c_c1_pos2,q_c2_pos1,q_c1_pos2,ref_s_pos2_c1, probM , cc , 2,1);
            if (flag_gap)
                return flag_gap;
        }
        while(ref_s_pos2_c1<ap2.getStart2()){
            noOverlapCheckgap(ap1 ,c_c2_pos1,q_c2_pos1,ref_s_pos2_c1, prob0 , tc, 2);
        }
        computeSOffset(ap2.getCigar2Unrolled(),c_c2_pos2,q_c2_pos2);
        while(ref_s_pos2_c1<=ref_e_pos1_c2 && ref_s_pos2_c1<=ref_e_pos2_c2){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c2_pos1,c_c2_pos2,q_c2_pos1,q_c2_pos2,ref_s_pos2_c1, probM , cc , 2,2);
            if (flag_gap)
                return flag_gap;
        }
        if(ref_s_pos2_c1-1 == ref_e_pos1_c2){
            while(ref_s_pos2_c1<=ref_e_pos2_c2){
                noOverlapCheckgap( ap2 ,c_c2_pos2,q_c2_pos2,ref_s_pos2_c1, prob0, tc,  2);
            }
        } else if(ref_s_pos2_c1-1 == ref_e_pos2_c2){
            while(ref_s_pos2_c1<=ref_e_pos1_c2){
                noOverlapCheckgap( ap1 ,c_c2_pos1,q_c2_pos1,ref_s_pos2_c1, prob0, tc,  2);
            }
        }
    } //      ----    --------- <-ap1
    //------------------          --------------
    else if(ref_s_pos2_c1 <= ref_s_pos1_c1 && ap1.getEnd1()>=ap2.getStart1() && ap1.getStart2() <=ap2.getEnd1() && ap1.getEnd2() < ap2.getStart2()){
        while(ref_s_pos2_c1<ref_s_pos1_c1){
            noOverlapCheckgap(ap2 ,c_c1_pos2,q_c1_pos2,ref_s_pos2_c1, prob0 , tc,  1);
        }
        while(ref_s_pos2_c1<=ap1.getEnd1()){
            flag_gap = overlapCheckgap( ap1, ap2 ,c_c1_pos1,c_c1_pos2,q_c1_pos1,q_c1_pos2,ref_s_pos2_c1, probM,  cc ,1,1);
            if (flag_gap)
                return flag_gap;
        }
        while(ref_s_pos2_c1<ap1.getStart2()){
            noOverlapCheckgap( ap2,c_c1_pos2,q_c1_pos2,ref_s_pos2_c1, prob0 , tc,  1);
        }
        computeSOffset(ap1.getCigar2Unrolled(),c_c2_pos1,q_c2_pos1);
        while(ref_s_pos2_c1<=ref_e_pos2_c1 && ref_s_pos2_c1<=ref_e_pos1_c2){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c2_pos1,c_c1_pos2,q_c2_pos1,q_c1_pos2,ref_s_pos2_c1, probM ,  cc ,2,1);
            if (flag_gap)
                return flag_gap;
        }
        if(ref_s_pos2_c1-1==ref_e_pos2_c1){
            while(ref_s_pos2_c1<=ref_e_pos1_c2){
                noOverlapCheckgap( ap1,c_c2_pos1,q_c2_pos1,ref_s_pos2_c1, prob0, tc,  2);
            }
        } else if(ref_s_pos2_c1-1 == ref_e_pos1_c2){
            while(ref_s_pos2_c1<=ref_e_pos2_c1){
                noOverlapCheckgap(ap2 ,c_c1_pos2,q_c1_pos2,ref_s_pos2_c1, prob0 , tc,  1);
            }
        }
        
    } //                -------    ------- <-ap1
    //-------------       ------------
    else if(ref_s_pos1_c1 <= ref_s_pos2_c2 && ap1.getStart1() > ap2.getEnd1() && ap1.getEnd1()>=ap2.getStart2() && ap1.getStart2() <= ap2.getEnd2()){
        while(ref_s_pos1_c1<ref_s_pos2_c2){
            noOverlapCheckgap(ap1 ,c_c1_pos1,q_c1_pos1,ref_s_pos1_c1, prob0,  tc, 1);
        }
        while(ref_s_pos1_c1<=ap1.getEnd1()){
            flag_gap = overlapCheckgap( ap1, ap2 ,c_c1_pos1,c_c2_pos2,q_c1_pos1,q_c2_pos2,ref_s_pos1_c1, probM , cc , 1,2);
            if (flag_gap)
                return flag_gap;
        }
        while(ref_s_pos1_c1<ap1.getStart2()){
            noOverlapCheckgap(ap2 ,c_c2_pos2,q_c2_pos2,ref_s_pos1_c1, prob0 , tc,  2);
        }
        computeSOffset(ap1.getCigar2Unrolled(),c_c2_pos1,q_c2_pos1);
        while(ref_s_pos1_c1<= ref_e_pos1_c2 && ref_s_pos1_c1<= ref_e_pos2_c2){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c2_pos1,c_c2_pos2,q_c2_pos1,q_c2_pos2,ref_s_pos1_c1, probM, cc , 2,2);
            if (flag_gap)
                return flag_gap;
        }
        if(ref_s_pos1_c1-1 == ref_e_pos1_c2){
            while(ref_s_pos1_c1<=ref_e_pos2_c2){
                noOverlapCheckgap(ap2 ,c_c2_pos2,q_c2_pos2,ref_s_pos1_c1, prob0, tc, 2);
            }
        } else if(ref_s_pos1_c1-1 == ref_e_pos2_c2){
            while(ref_s_pos1_c1<=ref_e_pos1_c2){
                noOverlapCheckgap(ap1 ,c_c2_pos1,q_c2_pos1,ref_s_pos1_c1, prob0, tc,  2);
            }
        }
        
    } //                   ---   -------- <-ap1
    //-------------     -----------
    else if(ref_s_pos2_c2<=ref_s_pos1_c1 && ap2.getEnd1()<ap1.getStart1() && ap1.getEnd1()>=ap2.getStart2() && ap1.getStart2() <= ap2.getEnd2()){
        while(ref_s_pos2_c2<ref_s_pos1_c1){
            noOverlapCheckgap(ap2 ,c_c2_pos2,q_c2_pos2,ref_s_pos2_c2, prob0 , tc,  2);
        }
        while(ref_s_pos2_c2<=ap1.getEnd1()){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c1_pos1,c_c2_pos2,q_c1_pos1,q_c2_pos2,ref_s_pos2_c2, probM ,  cc ,1,2);
            if (flag_gap)
                return flag_gap;
        }
        while(ref_s_pos2_c2<ap1.getStart2()){
            noOverlapCheckgap(ap2 ,c_c2_pos2,q_c2_pos2,ref_s_pos2_c2, prob0, tc,  2);
        }
        computeSOffset(ap1.getCigar2Unrolled(),c_c2_pos1,q_c2_pos1);
        while(ref_s_pos2_c2<=ref_e_pos1_c2 && ref_s_pos2_c2<=ref_e_pos2_c2){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c2_pos1,c_c2_pos2,q_c2_pos1,q_c2_pos2,ref_s_pos2_c2, probM , cc , 2,2);
            if (flag_gap)
                return flag_gap;
        }
        if(ref_s_pos2_c2-1==ref_e_pos1_c2){
            while(ref_s_pos2_c2<=ref_e_pos2_c2){
                noOverlapCheckgap( ap2 ,c_c2_pos2,q_c2_pos2,ref_s_pos2_c2, prob0, tc,  2);
            }
        } else if(ref_s_pos2_c2-1==ref_e_pos2_c2){
            while(ref_s_pos2_c2<=ref_e_pos1_c2){
                noOverlapCheckgap( ap1 ,c_c2_pos1,q_c2_pos1,ref_s_pos2_c2, prob0, tc,  2);
            }
        }
    }
    //  --------    ------------  <-ap1
    //-----    ----------------
    else if(ref_s_pos2_c1 <= ref_s_pos1_c1 && ap1.getStart1() <= ap2.getEnd1() && ap1.getEnd1() >= ap2.getStart2() && ap1.getStart2() <= ap2.getEnd2()){
        while(ref_s_pos2_c1<ref_s_pos1_c1){
            noOverlapCheckgap(ap2 ,c_c1_pos2,q_c1_pos2,ref_s_pos2_c1, prob0, tc,  1);
        }
        while(ref_s_pos2_c1<=ap2.getEnd1()){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c1_pos1,c_c1_pos2,q_c1_pos1,q_c1_pos2,ref_s_pos2_c1, probM, cc , 1,1);
            if (flag_gap)
                return flag_gap;
        }
        while(ref_s_pos2_c1<ap2.getStart2()){
            noOverlapCheckgap(ap1 ,c_c1_pos1,q_c1_pos1,ref_s_pos2_c1, prob0, tc,  1);
        }
        computeSOffset(ap2.getCigar2Unrolled(),c_c2_pos2,q_c2_pos2);
        while(ref_s_pos2_c1<=ap1.getEnd1()){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c1_pos1,c_c2_pos2,q_c1_pos1,q_c2_pos2,ref_s_pos2_c1, probM, cc , 1,2);
            if (flag_gap)
                return flag_gap;
        }
        while(ref_s_pos2_c1<ap1.getStart2()){
            noOverlapCheckgap(ap2 ,c_c2_pos2,q_c2_pos2,ref_s_pos2_c1, prob0 , tc,  2);
        }
        computeSOffset(ap1.getCigar2Unrolled(),c_c2_pos1,q_c2_pos1);
        while(ref_s_pos2_c1<=ref_e_pos2_c2 && ref_s_pos2_c1<=ref_e_pos1_c2){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c2_pos1,c_c2_pos2,q_c2_pos1,q_c2_pos2,ref_s_pos2_c1, probM, cc , 2,2);
            if (flag_gap)
                return flag_gap;
        }
        if(ref_s_pos2_c1-1==ref_e_pos2_c2){
            while(ref_s_pos2_c1<=ref_e_pos1_c2){
                noOverlapCheckgap(ap1 ,c_c2_pos1,q_c2_pos1,ref_s_pos2_c1, prob0 , tc, 2);
            }
        } else if(ref_s_pos2_c1-1==ref_e_pos1_c2){
            while(ref_s_pos2_c1<=ref_e_pos2_c2){
                noOverlapCheckgap(ap2 ,c_c2_pos2,q_c2_pos2,ref_s_pos2_c1, prob0, tc,  2);
            }
        }
        
    } //   --------------                 -------------- <-ap1
    //----------     -------------
    else if(ref_s_pos2_c1<=ref_s_pos1_c1 && ap1.getStart2() > ap2.getEnd2() && ap1.getStart1() <= ap2.getEnd1() && ap1.getEnd1() >= ap2.getStart2()){
        while(ref_s_pos2_c1<ref_s_pos1_c1){
            noOverlapCheckgap(ap2 ,c_c1_pos2,q_c1_pos2,ref_s_pos2_c1,prob0, tc, 1);
        }
        while(ref_s_pos2_c1<=ap2.getEnd1()){
            flag_gap = overlapCheckgap( ap1, ap2 ,c_c1_pos1,c_c1_pos2,q_c1_pos1,q_c1_pos2,ref_s_pos2_c1, probM ,  cc ,1,1);
            if (flag_gap)
                return flag_gap;
        }
        while(ref_s_pos2_c1<ap2.getStart2()){
            noOverlapCheckgap(ap1 ,c_c1_pos1,q_c1_pos1,ref_s_pos2_c1,prob0,  tc, 1);
        }
        computeSOffset(ap2.getCigar2Unrolled(),c_c2_pos2,q_c2_pos2);
        while(ref_s_pos2_c1<=ref_e_pos1_c1 && ref_s_pos2_c1<=ref_e_pos2_c2){
            flag_gap = overlapCheckgap(ap1, ap2 ,c_c1_pos1,c_c2_pos2,q_c1_pos1,q_c2_pos2,ref_s_pos2_c1, probM , cc , 1,2);
            if (flag_gap)
                return flag_gap;
        }
        if(ref_s_pos2_c1-1==ref_e_pos1_c1){
            while(ref_s_pos2_c1<=ref_e_pos2_c2){
                noOverlapCheckgap(ap2,c_c2_pos2,q_c2_pos2,ref_s_pos2_c1,prob0, tc,  2);
            }
        } else if(ref_s_pos2_c1-1==ref_e_pos2_c2){
            while(ref_s_pos2_c1<=ref_e_pos1_c1){
                noOverlapCheckgap(ap1 ,c_c1_pos1,q_c1_pos1,ref_s_pos2_c1,prob0,  tc, 1);
            }
        }
    }
    else{
        cout << " checkGapsPaired " << ap1.getName() << " , " << ap2.getName() << endl << flush;
        assert(0);
    }
    return flag_gap;
}

bool NewEdgeCalculator::checkGapsMixed(const AlignmentRecord& ap1, const AlignmentRecord& ap2, double& probM, double& prob0, int& cc, int& tc) const{
    bool flag_gap = false;
    if(ap2.isSingleEnd()){
        std::string dna, qualities, nucigar = "";
        int offset_s_f, offset_s_b, offset_p_f1, offset_p_f2, offset_p_b1, offset_p_b2 = 0;
        //get starting position and ending position according to ref position, paying attention to clipped bases
        //updated ref position including clips
        offset_s_f = computeOffset(ap2.getCigar1Unrolled());
        offset_s_b = computeRevOffset(ap2.getCigar1Unrolled());
        
        offset_p_f1 = computeOffset(ap1.getCigar1Unrolled());
        offset_p_f2 = computeOffset(ap1.getCigar2Unrolled());
        offset_p_b1 = computeRevOffset(ap1.getCigar1Unrolled());
        offset_p_b2 = computeRevOffset(ap1.getCigar2Unrolled());
        int ref_s_pos1 = ap2.getStart1()-offset_s_f;
        int ref_e_pos1 = ap2.getEnd1()+offset_s_b;
        int ref_p_s_pos1 = ap1.getStart1() -offset_p_f1;
        int ref_p_e_pos1 = ap1.getEnd1()+offset_p_b1;
        int ref_p_s_pos2 = ap1.getStart2() -offset_p_f2;
        int ref_p_e_pos2 = ap1.getEnd2()+offset_p_b2;
        //position in query sequences // phred scores
        int q_pos1 = 0;
        int q_p_pos1 = 0;
        int q_p_pos2 = 0;
        //position in unrolled cigar vectors
        int c_pos1 = 0;
        int c_p_pos1 = 0;
        int c_p_pos2 = 0;
        // ---------     -------- ->this (second read not changed)
        //----------
        if(ap2.getEnd1() < ap1.getStart2() ){
            flag_gap = checkGapsSingle(ap1, ap2, probM , prob0, cc,tc, 1,1);
            if(flag_gap)
                return flag_gap;
            assert(ap1.getEnd1() < ap1.getStart2() );
        } // -------     -------- ->this (second read not changed)
        //             ----------
        else if (ap2.getStart1() > ap1.getEnd1()){
            
            flag_gap = checkGapsSingle(ap1, ap2 , probM , prob0, cc, tc, 2,1);
            if(flag_gap)
                return flag_gap;
            assert(ap1.getEnd1() < ap1.getStart2() );
        } //----------          -----------   ->this OR  -------       -----------
        //-------------------------------              -----------------------------
        else if(ref_s_pos1 <= ref_p_s_pos1){
            while(ref_s_pos1<ref_p_s_pos1){
                noOverlapCheckgap(ap2 , c_pos1, q_pos1,ref_s_pos1,prob0,  tc, 1);
            }
            while(ref_s_pos1<=ap1.getEnd1()){
                flag_gap = overlapCheckgap(ap1, ap2 ,c_p_pos1,c_pos1,q_p_pos1,q_pos1,ref_s_pos1, probM , cc , 1,1);
                if (flag_gap)
                    return flag_gap;
            }
            while(ref_s_pos1<ap1.getStart2() ){
                noOverlapCheckgap( ap2,c_pos1,q_pos1,ref_s_pos1,prob0,  tc, 1);
            }
            computeSOffset(ap1.getCigar2Unrolled(),c_p_pos2,q_p_pos2);
            while(ref_s_pos1<=ref_p_e_pos2 && ref_s_pos1 <= ref_e_pos1){
                flag_gap = overlapCheckgap(ap1, ap2 ,c_p_pos2,c_pos1,q_p_pos2,q_pos1,ref_s_pos1, probM , cc , 2,1);
                if (flag_gap)
                    return flag_gap;
            }
            if(ref_s_pos1-1 == ref_p_e_pos2){
                while(ref_s_pos1<=ref_e_pos1){
                    noOverlapCheckgap(ap2 ,c_pos1,q_pos1,ref_s_pos1,prob0, tc,  1);
                }
            } else if (ref_s_pos1-1 == ref_e_pos1){
                while(ref_s_pos1<=ref_p_e_pos2){
                    noOverlapCheckgap(ap1 ,c_p_pos2,q_p_pos2,ref_s_pos1,prob0, tc,  2);
                }
            }
        } //----------          ------------ ->this OR ----------       -----------
        //     -------------------------------            ----------------------
        else if(ref_s_pos1 >= ref_p_s_pos1){
            while(ref_p_s_pos1 < ref_s_pos1){
                noOverlapCheckgap(ap1 , c_p_pos1,q_p_pos1,ref_p_s_pos1,prob0,  tc, 1);
            }
            while(ref_p_s_pos1<=ap1.getEnd1()){
                flag_gap = overlapCheckgap(ap1, ap2 ,c_p_pos1,c_pos1,q_p_pos1,q_pos1,ref_p_s_pos1, probM , cc , 1,1);
                if (flag_gap)
                    return flag_gap;
            }
            while(ref_p_s_pos1<ap1.getStart2() ){
                noOverlapCheckgap(ap2 ,c_pos1,q_pos1,ref_p_s_pos1,prob0,  tc, 1);
            }
            computeSOffset(ap1.getCigar2Unrolled(),c_p_pos2,q_p_pos2);
            while(ref_p_s_pos1<=ref_p_e_pos2 && ref_p_s_pos1 <= ref_e_pos1){
                flag_gap = overlapCheckgap(ap1, ap2 ,c_p_pos2,c_pos1,q_p_pos2,q_pos1,ref_p_s_pos1, probM , cc , 2,1);
                if (flag_gap)
                    return flag_gap;
            }
            if(ref_p_s_pos1-1 == ref_p_e_pos2){
                while(ref_p_s_pos1<=ref_e_pos1){
                    noOverlapCheckgap(ap2 ,c_pos1,q_pos1,ref_p_s_pos1,prob0, tc,  1);
                }
            } else if (ref_p_s_pos1-1 == ref_e_pos1){
                while(ref_p_s_pos1<=ref_p_e_pos2){
                    noOverlapCheckgap(ap1 ,c_p_pos2,q_p_pos2,ref_p_s_pos1,prob0, tc,  2);
                }
            }
        }
        else{
            cout << " checkGapsMixed 1 " << ap1.getName() << " , " << ap2.getName() << endl << flush;
            assert(0);
        }
    }
    else if (ap2.isPairedEnd()){
        std::string dna, qualities, nucigar = "";
        int offset_s_f, offset_s_b, offset_p_f1, offset_p_f2, offset_p_b1, offset_p_b2 = 0;
        //get starting position and ending position according to ref position, paying attention to clipped bases
        //updated ref position including clips
        offset_s_f = computeOffset(ap1.getCigar1Unrolled());
        offset_s_b = computeRevOffset(ap1.getCigar1Unrolled());
        offset_p_f1 = computeOffset(ap2.getCigar1Unrolled());
        offset_p_f2 = computeOffset(ap2.getCigar2Unrolled());
        offset_p_b1 = computeRevOffset(ap2.getCigar1Unrolled());
        offset_p_b2 = computeRevOffset(ap2.getCigar2Unrolled());
        int ref_s_pos1 = ap1.getStart1() -offset_s_f;
        int ref_e_pos1 = ap1.getEnd1()+offset_s_b;
        int ref_p_s_pos1 = ap2.getStart1()-offset_p_f1;
        int ref_p_e_pos2 = ap2.getEnd2()+offset_p_b2;
        //position in query sequences // phred scores
        int q_pos1 = 0;
        int q_p_pos1 = 0;
        int q_p_pos2 = 0;
        //position in unrolled cigar vectors
        int c_pos1 = 0;
        int c_p_pos1 = 0;
        int c_p_pos2 = 0;
        
        // ---------  -----------          OR ---------- ------------
        // ---------               ->this                -------------
        if(ap1.getEnd1() < ap2.getStart2()){
            flag_gap = checkGapsSingle(ap1, ap2 , probM , prob0, cc, tc, 1,1);
            if(flag_gap)
                return flag_gap;
            
        } else if (ap1.getStart1()  > ap2.getEnd1()){
            flag_gap = checkGapsSingle(ap1, ap2 , probM , prob0, cc, tc, 1,2);
            if(flag_gap)
                return flag_gap;
            
        }//----------          -----------        OR  -------       -----------
        //----------------------------    <-this    -----------------------------
        else if(ref_s_pos1 <= ref_p_s_pos1){
            
            while(ref_s_pos1 < ref_p_s_pos1){
                noOverlapCheckgap(ap1 , c_pos1, q_pos1,ref_s_pos1,prob0, tc,  1);
            }
            while(ref_s_pos1<=ap2.getEnd1()){
                flag_gap = overlapCheckgap(ap1, ap2 ,c_pos1,c_p_pos1,q_pos1,q_p_pos1,ref_s_pos1, probM , cc , 1,1);
                if (flag_gap){
                    return flag_gap;
                }
            }
            while(ref_s_pos1<ap2.getStart2()){
                noOverlapCheckgap( ap1 ,c_pos1,q_pos1,ref_s_pos1,prob0, tc,  1);
            }
            computeSOffset(ap2.getCigar2Unrolled(),c_p_pos2,q_p_pos2);
            while(ref_s_pos1<=ref_p_e_pos2 && ref_s_pos1 <= ref_e_pos1){
                flag_gap = overlapCheckgap(ap1, ap2 ,c_pos1,c_p_pos2,q_pos1,q_p_pos2,ref_s_pos1, probM , cc , 1,2);
                if (flag_gap)
                    return flag_gap;
            }
            if(ref_s_pos1-1 == ref_p_e_pos2){
                while(ref_s_pos1<=ref_e_pos1){
                    noOverlapCheckgap(ap1 ,c_pos1,q_pos1,ref_s_pos1,prob0,  tc, 1);
                }
            } else if (ref_s_pos1-1 == ref_e_pos1){
                while(ref_s_pos1<=ref_p_e_pos2){
                    noOverlapCheckgap(ap2 ,c_p_pos2,q_p_pos2,ref_s_pos1,prob0,  tc, 2);
                }
            }
            
        }//----------          ------------        OR ----------       -----------
        //     -------------------------------   <-this   ----------------------
        else if(ref_s_pos1 >= ref_p_s_pos1){
            while(ref_p_s_pos1 < ref_s_pos1){
                noOverlapCheckgap(ap2 , c_p_pos1,q_p_pos1,ref_p_s_pos1,prob0,  tc, 1);
            }
            while(ref_p_s_pos1<=ap2.getEnd1()){
                flag_gap = overlapCheckgap(ap1, ap2 ,c_pos1,c_p_pos1,q_pos1,q_p_pos1,ref_p_s_pos1, probM , cc , 1,1);
                if (flag_gap)
                    return flag_gap;
            }
            while(ref_p_s_pos1<ap2.getStart2()){
                noOverlapCheckgap(ap1 ,c_pos1,q_pos1,ref_p_s_pos1,prob0, tc,  1);
            }
            computeSOffset(ap2.getCigar2Unrolled(),c_p_pos2,q_p_pos2);
            while(ref_p_s_pos1<=ref_p_e_pos2 && ref_p_s_pos1 <= ref_e_pos1){
                flag_gap = overlapCheckgap(ap1, ap2 ,c_pos1,c_p_pos2,q_pos1,q_p_pos2,ref_p_s_pos1, probM , cc , 1,2);
                if (flag_gap)
                    return flag_gap;
            }
            if(ref_p_s_pos1-1 == ref_p_e_pos2){
                while(ref_p_s_pos1<=ref_e_pos1){
                    noOverlapCheckgap(ap1 ,c_pos1,q_pos1,ref_p_s_pos1,prob0, tc,  1);
                }
            } else if (ref_p_s_pos1-1 == ref_e_pos1){
                while(ref_p_s_pos1<=ref_p_e_pos2){
                    noOverlapCheckgap(ap2 ,c_p_pos2,q_p_pos2,ref_p_s_pos1,prob0, tc,  2);
                }
            }
        }
        else{
            assert(0);
        }
    }
    else{
        assert(0);
    }
    return flag_gap;
}

bool NewEdgeCalculator::checkGapsCigar(const AlignmentRecord& ap1, const AlignmentRecord& ap2, double& probM, double& prob0, int& cc, int& tc) const{
    
    if (ap1.isSingleEnd() && ap2.isSingleEnd()){
        return (checkGapsSingle(ap1, ap2, probM , prob0 ,cc, tc, 1,1));
    }
    else if (!(ap1.isSingleEnd()) && ap2.isPairedEnd()){
        return (checkGapsPaired(ap1, ap2, probM , prob0, cc, tc ));
    }
    else {
        return(checkGapsMixed(ap1, ap2, probM , prob0 , cc, tc ));
    }
}

//TO DO: find out whether gaps / insertions are compatible
/*bool NewEdgeCalculator::checkGaps(const std::vector<AlignmentRecord::mapValue> & cov_ap1,const std::vector<AlignmentRecord::mapValue> & cov_ap2,const std::vector<std::pair<int,int>> & aub) const{
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
            //if(num_inserts == 0) num_inserts+= pos_diff1-pos_diff2;
            //else
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
}*/

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
    
    if (cc<=MIN_OVERLAP*std::min(cov_ap1.size(),cov_ap2.size())) return false;
    double prob;
    double potence;

    //TEST to see what happens if prob0 is neglected
    if(this->NOPROB0){
        prob = std::log10(probM);
        potence = 1.0/cc;
    } else {
        prob = std::log10(probM) + prob0;
        potence = 1.0/(cc+tc);
    }
    double final_prob = prob*potence;
    return final_prob >= cutoff;
}

/*void NewEdgeCalculator::iterateCovAp(bool pe1, unsigned int& pos2, double& probM, unsigned int& pos1, unsigned int& equalBase, bool pe2, const std::vector<AlignmentRecord::mapValue>& cov_ap2, int& tc, double& prob0, const AlignmentRecord& ap2, const std::vector<AlignmentRecord::mapValue>& cov_ap1, const AlignmentRecord& ap1, int& cc) const
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
           } //else if (!pe1 && pe2){
                //calculateProb0(v1,prob0);
                //tc++;
                //pos1++;
           //}else if (pe1 && pe2){
                //if(ap1.getEnd1() < ap2.getStart1() && v1.ref < ap1.getStart2()) pos1++;
                //else{
                //    calculateProb0(v1,prob0);
                //    tc++;
                //    pos1++;
                //}
           //}
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
               } //else if (pe1 && pe2){
                   // if(ap2.getEnd1() < ap1.getStart1() && v1.ref < ap2.getStart2()) pos2++;
                   // else{
                   //     calculateProb0(v2,prob0);
                   //     tc++;
                   //     pos2++;
                   // }
               //}
        }
    }
}*/

/*void NewEdgeCalculator::iterateRemainingCovAp(int& tc, bool pe1, const std::vector<AlignmentRecord::mapValue>& cov_ap2, bool pe2, double& prob0, const AlignmentRecord& ap1, const AlignmentRecord& ap2, const std::vector<AlignmentRecord::mapValue>& cov_ap1, unsigned int& pos1) const
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
}*/

/*void NewEdgeCalculator::iterateRemainingCovAp2(const AlignmentRecord& ap2, const std::vector<AlignmentRecord::mapValue>& cov_ap2, bool pe1, const AlignmentRecord& ap1, int& tc, double& prob0, const std::vector<AlignmentRecord::mapValue>& cov_ap1, bool pe2, unsigned int& pos2) const
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
}*/

bool NewEdgeCalculator::edgeBetween(const AlignmentRecord & ap1, const AlignmentRecord & ap2) const{
    const auto& cov_ap1 = ap1.getCovmap();
    const auto& cov_ap2 = ap2.getCovmap();
    // -----    ------
    //                   -----     ----
    if(cov_ap1.back().ref < cov_ap2[0].ref  || cov_ap2.back().ref < cov_ap1[0].ref){
        return false;
    }
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
    if (ap1.isSingleEnd() && ap2.isSingleEnd()){
        unsigned int e11= ap1.getEnd1();
        unsigned int s11= ap1.getStart1();
        unsigned int e21= ap2.getEnd1();
        unsigned int s21= ap2.getStart1();
        //-----
        //        -------
        if (e11 < s21){
            return false;
        }
        //        -------
        //-----
        else if (e21 < s11){
            return false;
        }
    }
    
    if(pe1 && pe2){
        unsigned int e11= ap1.getEnd1();
        unsigned int e12= ap1.getEnd2();
        unsigned int s11= ap1.getStart1();
        unsigned int s12= ap1.getStart2();
        unsigned int e21= ap2.getEnd1();
        unsigned int e22= ap2.getEnd2();
        unsigned int s21= ap2.getStart1();
        unsigned int s22= ap2.getStart2();
        //--------   ---------              OR ------   --------
        //           ---------   -------                          -------   --------
        if(e11 < s21 && e12 < s22){
            return false;
        } //         ---------    ---------  OR                    -----     -------
        //-------    ---------                   ----    ------
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
        }//          --------- ---------   OR          ------             ------
        //----------           ---------         ----            ------
        else if (e21 < s11 && e11 < s22){
            return false;
        }//--------             --------          ------          --------
        //           --------   --------   OR             -----              ------
        else if(e11 < s21 && e21 < s12){
            return false;
        }
    }
    
    if (ap1.isSingleEnd() && ap2.isPairedEnd()){
        unsigned int e11= ap1.getEnd1();
        unsigned int s11= ap1.getStart1();
        unsigned int e21= ap2.getEnd1();
        unsigned int s21= ap2.getStart1();
        unsigned int e22= ap2.getEnd2();
        unsigned int s22= ap2.getStart2();
        
        // -------
        //          ------   -------
        if (e11 < s21){
            return false;
        }
        //        -------
        //-----            -------
        else if (e21 < s11 && e11 < s22){
            return false;
        }
        //              -------
        //----   -----
        else if (e22 < s11){
            return false;
        }
    } else if (ap1.isPairedEnd() && ap2.isSingleEnd()){
        unsigned int e11= ap1.getEnd1();
        unsigned int e12= ap1.getEnd2();
        unsigned int s11= ap1.getStart1();
        unsigned int s12= ap1.getStart2();
        unsigned int e21= ap2.getEnd1();
        unsigned int s21= ap2.getStart1();
        // -----   ------
        //                 ------
        if(e12 < s21){
            return false;
        } // ----           ------
        //          ------
        else if (e11 < s21 && e21 < s12){
            return false;
        }
        //           -----    -----
        // ------
        else if (e21 < s11){
            return false;
        }
    }
    
    if (checkGapsCigar(ap1, ap2, probM, prob0, cc, tc)) {
        return false;
    }
    
    if (cc == 0){
        return false;
    }
    
    return similarityCriterion(ap1, cov_ap1, ap2, cov_ap2, probM, prob0, tc, cc);
}

void NewEdgeCalculator::getPartnerLengthRange(const AlignmentRecord& ap, unsigned int* min, unsigned int* max) const {
    assert(false);
}

void NewEdgeCalculator::setOverlapCliques(double d){
    this->MIN_OVERLAP_CLIQUES=d;
}

double NewEdgeCalculator::getOverlapCliques() const {
    return this->MIN_OVERLAP_CLIQUES;
}