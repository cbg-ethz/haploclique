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

#include <vector>
#include <algorithm>
#include <math.h>
#include <ctype.h>
#include <map>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/compare.hpp>

#include "AlignmentRecord.h"
#include "Clique.h"

using namespace std;
using namespace boost;

int phred_sum(const string& phred, char phred_base=33) {
	int result = 0;
	for (size_t i=0; i<phred.size(); ++i) {
		result += phred[i] - phred_base;
	}
	return result;
}

//called by getmergedDnaSequence() to create new CigarData
std::vector<BamTools::CigarOp> createCigar(std::string nucigar){
    std::vector<BamTools::CigarOp> res;
    unsigned int counter = 1;
    unsigned int pos = 0;
    char c;
    while(pos < nucigar.size()){
        c = nucigar[pos];
        pos++;
        while(pos < nucigar.size()){
            if (c == nucigar[pos]){
                counter++;
                pos++;
            } else {
                res.push_back(BamTools::CigarOp(c,counter));
                counter = 1;
                break;
            }
        }
        if(pos == nucigar.size()){
            res.push_back(BamTools::CigarOp(c,counter));
        }
    }
    return res;
}

//helper functions to merge DNA sequences
int agreement(const char& qual1, const char& qual2){
    float prob1 = std::pow(10,(float)-qual1/10);
    float prob2 = std::pow(10,(float)-qual2/10);
    float posterior = (prob1*prob2/3)/(1-prob1-prob2+4*prob1*prob2/3);
    posterior = std::round(-10*log10(posterior));
    return posterior;
}

int disagreement(const char& qual1, const char& qual2){
    float prob1 = std::pow(10,(float)-qual1/10);
    float prob2 = std::pow(10,(float)-qual2/10);
    float posterior = ((prob1*(1-prob2/3))/(prob1+prob2-4*prob1*prob2/3));
    posterior = std::round(-10*log10(posterior));
    return posterior;
}

double phredProb(const char& qual){
    return std::pow(10, (double)(-qual-33)/10.0);
}

std::pair<char,char> computeEntry(const char& base1, const char& qual1, const char& base2, const char& qual2){
    std::pair<char,char> result;

    if (base1==base2){
        result.first = base1;
        result.second = std::min(agreement(qual1-33,qual2-33)+33,126);
    }
    else if (qual1>=qual2) {
        result.first = base1;
        result.second = disagreement(qual1-33,qual2-33)+33;
    } else {
        result.first = base2;
        result.second = disagreement(qual2-33,qual1-33)+33;
    }
    return result;
}


AlignmentRecord::AlignmentRecord(const BamTools::BamAlignment& alignment, int readRef, vector<string>* rnm) : readNameMap(rnm) {
    this->single_end = true;
    this->readNames.insert(readRef);
    this->name = alignment.Name;
    this->start1 = alignment.Position + 1;
    this->end1 = alignment.GetEndPosition();
    this->cigar1 = alignment.CigarData;
    this->sequence1 = ShortDnaSequence(alignment.QueryBases, alignment.Qualities);
    this->phred_sum1 = phred_sum(alignment.Qualities);
    this->length_incl_deletions1 = this->sequence1.size();
    this->length_incl_longdeletions1 = this->sequence1.size();
	for (const auto& it : cigar1) {
        for (unsigned int s = 0; s < it.Length; ++s) {
          	this->cigar1_unrolled.push_back(it.Type);
        }
        if (it.Type == 'D') {
      		this->length_incl_deletions1+=it.Length;
       		if (it.Length > 1) {
       			this->length_incl_longdeletions1+=it.Length;
       		}
        }/*for Debugging
        if (it.Type == 'H'){
            int k = 0;
        }
        if (it.Type == 'N'){
            cout << alignment.Name << endl;
        }
        if (it.Type == 'P'){
           cout << alignment.Name << endl;
        }
        if (it.Type == 'I'){
            int k = 0;
        }*/
    }
    this->cov_pos = this->coveredPositions();
}

AlignmentRecord::AlignmentRecord(unique_ptr<vector<const AlignmentRecord*>>& alignments, unsigned int clique_id) : cigar1_unrolled(), cigar2_unrolled() {
    deque<pair<int, int>> interval;
    vector<ShortDnaSequence> sequences;
    vector<vector<BamTools::CigarOp>> cigars;
    int ct = 0;

    this->readNameMap = (*alignments)[0]->readNameMap;

    for (auto& al : *alignments) {
        if (al->isPairedEnd()) {
            sequences.push_back(al->getSequence1());
            cigars.push_back(al->getCigar1());
            interval.push_back(pair<int, int>(ct, al->getStart1() - 1));
            interval.push_back(pair<int, int>(ct, al->getStart1() - 1 + al->getSequence1().size()));
            ct++;
            sequences.push_back(al->getSequence2());
            cigars.push_back(al->getCigar2());
            interval.push_back(pair<int, int>(ct, al->getStart2() - 1));
            interval.push_back(pair<int, int>(ct, al->getStart2() - 1 + al->getSequence2().size()));
            ct++;
        } else {
            sequences.push_back(al->getSequence1());
            cigars.push_back(al->getCigar1());
            interval.push_back(pair<int, int>(ct, al->getStart1() - 1));
            interval.push_back(pair<int, int>(ct, al->getStart1() - 1 + al->getSequence1().size()));
            ct++;
        }

        this->readNames.insert(al->readNames.begin(), al->readNames.end());
    }

    auto comp = [](pair<int, int> a, pair<int, int> b) { return a.second < b.second; };
    std::sort(interval.begin(), interval.end(), comp);

    //changes this->sequence, this->cigar, this->single_end, this->phredsum, this->cigar, this->sequence
    mergeSequences(interval, sequences, cigars);

    this->name = "Clique_" + to_string(clique_id);

    this->length_incl_deletions1 = this->sequence1.size();
    this->length_incl_longdeletions1 = this->sequence1.size();

    unsigned int length_ct = 0; // DEBUG
	for (const auto& it : cigar1) { //DEBUG
        length_ct += it.Length;        

        for (unsigned int s = 0; s < it.Length; ++s) {
          	this->cigar1_unrolled.push_back(it.Type);
        }
        if (it.Type == 'D') {
      		this->length_incl_deletions1+=it.Length;
       		if (it.Length > 1) {
       			this->length_incl_longdeletions1+=it.Length;
       		}
       	}
    }
    assert(length_ct == (unsigned)this->length_incl_deletions1); //DEBUG


    if(not this->single_end) {
        this->length_incl_deletions2 = this->sequence2.size();
        this->length_incl_longdeletions2 = this->sequence2.size();

        length_ct = 0; //DEBUG
	    for (const auto& it : cigar2) {
            length_ct += it.Length; //DEBUG

            for (unsigned int s = 0; s < it.Length; ++s) {
              	this->cigar2_unrolled.push_back(it.Type);
            }
            if (it.Type == 'D') {
          		this->length_incl_deletions2+=it.Length;
           		if (it.Length > 1) {
           			this->length_incl_longdeletions2+=it.Length;
           		}
           	}
        }
        assert(length_ct == (unsigned)this->length_incl_deletions2); //DEBUG
    }
    this->cov_pos=this->coveredPositions();
}

void AlignmentRecord::mergeSequences(std::deque<std::pair<int, int> > intervals, std::vector<ShortDnaSequence> &to_merge, std::vector<std::vector<BamTools::CigarOp> > &cigars) {

    assert(intervals.size() > 2);

    pair<int, int> p = intervals.front();
    intervals.pop_front();

    map<int, int> overlaps;
    overlaps.insert(p);
    int next = p.second;

    // Find first overlapping interval
    while(overlaps.size() < 2) {
        p = intervals.front();
        intervals.pop_front();
        auto it = overlaps.find(p.first);
        if (it == overlaps.end()) {
            overlaps.insert(p);
        } else {
            overlaps.erase(it);
        }
        next = p.second;
    }

    bool paired = false;

    string seq;
    string qual;
    vector<BamTools::CigarOp> cigar;

    this->start1 = next + 1;
    this->single_end = true;

    // Go over all intervals
    while(not intervals.empty()) {
        if (overlaps.size() < 2) {
            if (paired) {
                this->end2 = next + 1;
                this->phred_sum2 = phred_sum(qual);
                this->sequence2 = ShortDnaSequence(seq, qual);
                this->cigar2 = cigar;
                this->single_end = false;
                break;
            } else {
                this->phred_sum1 = phred_sum(qual);
                this->sequence1 = ShortDnaSequence(seq, qual);
                this->end1 = next + 1;
                this->cigar1 = cigar;
                seq = "";
                qual = "";
                cigar.clear();

                paired = true;

                // Find next overlapping interval
                while(overlaps.size() < 2 and not intervals.empty()) {
                    p = intervals.front();
                    intervals.pop_front();
                    auto it = overlaps.find(p.first);
                    if (it == overlaps.end()) {
                        overlaps.insert(p);
                    } else {
                        overlaps.erase(it);
                    }
                    next = p.second;
                }
                if (intervals.empty()) break;
                this->start2 = next + 1;
            }
        }

        int i = next;
        p = intervals.front();
        intervals.pop_front();
        auto it = overlaps.find(p.first);
        next = p.second;

        getCigarInterval(i, next, cigar, cigars[overlaps.begin()->first], overlaps.begin()->second);

        for (; i < next; i++) {
            map<char, double> basemap;
            
            double max_val = 0.0;
            char max_key = 0;

            for ( auto mapit : overlaps) {
                char c = to_merge[mapit.first][i - mapit.second];
                if (basemap.find(c) == basemap.end()) basemap[c] = 0.0;
                basemap[c] += to_merge[mapit.first].qualityCorrect(i - mapit.second) / overlaps.size();
                if (basemap[c] > max_val){
                    max_key = c;
                    max_val = basemap[c];
                }
            }

            seq.push_back(max_key);
            qual.push_back( (char) round(-10*log10(1-basemap[max_key]) + 33) );           
        }

        if (it == overlaps.end()) {
            overlaps.insert(p);
        } else {
            overlaps.erase(it);
        }
    }
}

//helper functions for merging DNA Sequences to create combined Alignment Record
void AlignmentRecord::noOverlapMerge(std::string& dna, std::string& qualities, std::string& nucigar, int& c_pos, int& q_pos, int& ref_pos){
    char c = this->cigar1_unrolled[c_pos];
    if (c == 'H'){
        nucigar += 'H';
        ref_pos++;
        c_pos++;
    } else if (c == 'I') {
        dna += this->sequence1[q_pos];
        qualities += this->sequence1.qualityChar(q_pos);
        nucigar += 'I';
        q_pos++;
        c_pos++;
    } else if (c == 'D') {
        nucigar += 'D';
        ref_pos++;
        c_pos++;
    } else {
        dna += this->sequence1[q_pos];
        qualities += this->sequence1.qualityChar(q_pos);
        nucigar += c;
        ref_pos++;
        q_pos++;
        c_pos++;
    }
}

void AlignmentRecord::noOverlapMerge(const BamTools::BamAlignment& alignment, std::string& dna, std::string& qualities, std::string& nucigar, std::vector<char>& cigar_temp_unrolled, int& c_pos, int& q_pos, int& ref_pos){
    char c = cigar_temp_unrolled[c_pos];
    if (c == 'H'){
        nucigar += 'H';
        ref_pos++;
        c_pos++;
    } else if (c == 'I') {
        dna += alignment.QueryBases[q_pos];
        qualities += alignment.Qualities[q_pos];
        nucigar += 'I';
        q_pos++;
        c_pos++;
    } else if (c == 'D') {
        nucigar += 'D';
        ref_pos++;
        c_pos++;
    } else {
        dna += alignment.QueryBases[q_pos];
        qualities += alignment.Qualities[q_pos];
        nucigar += c;
        ref_pos++;
        q_pos++;
        c_pos++;
    }
}

void AlignmentRecord::overlapMerge(const BamTools::BamAlignment& alignment, std::string& dna, std::string& qualities, std::string& nucigar, std::vector<char>& cigar_temp_unrolled, int& c_pos1, int& c_pos2, int& q_pos1, int& q_pos2, int& ref_pos){
    char c1 = this->cigar1_unrolled[c_pos1];
    char c2 = cigar_temp_unrolled[c_pos2];
    if((c1 == 'M' && c2 == 'M') || (c1 == 'S' && c2 == 'S') || (c1 == 'I' && c2 == 'I')){
        std::pair<char,char> resPair = computeEntry(this->sequence1[q_pos1],this->sequence1.qualityChar(q_pos1),alignment.QueryBases[q_pos2],alignment.Qualities[q_pos2]);
        dna += resPair.first;
        qualities += resPair.second;
        nucigar += c1;
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
            nucigar += 'D';
        } else {
            nucigar += 'H';
        }
        if (c1 == 'S'){
            q_pos1++;
        } else if(c2 == 'S'){
            q_pos2++;
        }
    } else if ((c1 == 'M' && (c2 == 'D' || c2 == 'H' || c2 == 'S')) || ((c1 == 'D' || c1 == 'H' || c1 == 'S') && c2 == 'M') || (c1 == 'S' && c2 == 'H') || (c1 == 'H' && c2 == 'S')) {
        if (c1 == 'M'){
            nucigar += 'M';
            dna += this->sequence1[q_pos1];
            qualities += this->sequence1.qualityChar(q_pos1);
            ref_pos++;
            c_pos1++;
            c_pos2++;
            q_pos1++;
            if (c2 == 'S') q_pos2++;
        } else if (c2 == 'M'){
            nucigar += 'M';
            dna +=  alignment.QueryBases[q_pos2];
            qualities += alignment.Qualities[q_pos2];
            ref_pos++;
            c_pos1++;
            c_pos2++;
            q_pos2++;
            if (c1 == 'S') q_pos1++;
        } else if (c1 == 'S'){
            nucigar += 'S';
            dna += this->sequence1[q_pos1];
            qualities += this->sequence1.qualityChar(q_pos1);
            ref_pos++;
            c_pos1++;
            c_pos2++;
            q_pos1++;
        } else {
            nucigar += 'S';
            dna +=  alignment.QueryBases[q_pos2];
            qualities += alignment.Qualities[q_pos2];
            ref_pos++;
            c_pos1++;
            c_pos2++;
            q_pos2++;
        }
    } else if (c1 == 'I' || c2 == 'I'){
        if(c1 == 'I'){
            nucigar += 'I';
            dna += this->sequence1[q_pos1];
            qualities += this->sequence1.qualityChar(q_pos1);
            c_pos1++;
            q_pos1++;
        } else {
            nucigar += 'I';
            dna +=  alignment.QueryBases[q_pos2];
            qualities += alignment.Qualities[q_pos2];
            c_pos2++;
            q_pos2++;
        }
    }
}

//creates merged Dna Sequence for overlapping paired end reads, creating new Cigar
void AlignmentRecord::getMergedDnaSequence(const BamTools::BamAlignment& alignment){
        std::string dna = "";
        std::string qualities = "";
        std::string nucigar = "";
        std::vector<char> cigar_temp_unrolled;
        //vector of CigarOp
        for (const auto& it : alignment.CigarData) {
            for (unsigned int s = 0; s < it.Length; ++s) {
                cigar_temp_unrolled.push_back(it.Type);
            }
        }
        //get starting position and ending position according to ref position, paying attention to clipped bases
        int offset_f1 = 0;
        int offset_f2 = 0;
        int offset_b1 = 0;
        int offset_b2 = 0;
        for(auto i : this->cigar1_unrolled){
            if (i == 'S' || i == 'H'){
                offset_f1++;
            } else break;
        }
        for (auto i : cigar_temp_unrolled){
            if (i == 'S' || i == 'H'){
                offset_f2++;
            } else break;
        }
        for (std::vector<char>::reverse_iterator it = this->cigar1_unrolled.rbegin(); it != this->cigar1_unrolled.rend(); ++it){
          if (*it == 'S' || *it == 'H'){
              offset_b1++;
          } else break;
        }
        for (std::vector<char>::reverse_iterator it = cigar_temp_unrolled.rbegin(); it != cigar_temp_unrolled.rend(); ++it){
          if (*it == 'S' || *it == 'H'){
              offset_b2++;
          } else break;
        }
        //updated ref position including clips
        int ref_s_pos1 = this->start1-offset_f1;
        int ref_e_pos1 = this->end1+offset_b1;
        int ref_s_pos2 = alignment.Position+1-offset_f2;
        int ref_e_pos2 = alignment.GetEndPosition()+offset_b2;
        //position in query sequences // phred scores
        int q_pos1 = 0;
        int q_pos2 = 0;
        //position in unrolled cigar vectors
        int c_pos1 = 0;
        int c_pos2 = 0;
        //4 cases of different overlaps
        //------------
        //     ------------
        if(ref_s_pos1 <= ref_s_pos2 && ref_e_pos1 <= ref_e_pos2){
            while(ref_s_pos1<ref_s_pos2){
                noOverlapMerge(dna,qualities,nucigar,c_pos1,q_pos1,ref_s_pos1);
            }
            while(ref_s_pos1<=ref_e_pos1){
                overlapMerge(alignment,dna,qualities,nucigar,cigar_temp_unrolled,c_pos1,c_pos2,q_pos1,q_pos2,ref_s_pos1);
            }
            while(ref_s_pos1<=ref_e_pos2){
                noOverlapMerge(alignment, dna, qualities, nucigar, cigar_temp_unrolled, c_pos2, q_pos2, ref_s_pos1);
            }
        }//------------------------------
            //           ----------
         else if (ref_s_pos1 <= ref_s_pos2 && ref_e_pos1 > ref_e_pos2){
            while(ref_s_pos1<ref_s_pos2){
                noOverlapMerge(dna,qualities,nucigar,c_pos1,q_pos1,ref_s_pos1);
            }
            while(ref_s_pos1<=ref_e_pos2){
                overlapMerge(alignment,dna,qualities,nucigar,cigar_temp_unrolled,c_pos1,c_pos2,q_pos1,q_pos2,ref_s_pos1);
            }
            while(ref_s_pos1<=ref_e_pos1){
                noOverlapMerge(dna,qualities,nucigar,c_pos1,q_pos1,ref_s_pos1);
            }
            //         ----------
            //--------------------------
        } else if (ref_s_pos1 > ref_s_pos2 && ref_e_pos1 <= ref_e_pos2){
            while(ref_s_pos2<ref_s_pos1){
                noOverlapMerge(alignment, dna, qualities, nucigar, cigar_temp_unrolled, c_pos2, q_pos2, ref_s_pos2);
            }
            while(ref_s_pos2<=ref_e_pos1){
                overlapMerge(alignment,dna,qualities,nucigar,cigar_temp_unrolled,c_pos1,c_pos2,q_pos1,q_pos2,ref_s_pos2);
            }
            while(ref_s_pos2<=ref_e_pos2){
                noOverlapMerge(alignment, dna, qualities, nucigar, cigar_temp_unrolled, c_pos2, q_pos2, ref_s_pos2);
            }
            //           --------------------
            //---------------------
        } else {
            assert(ref_s_pos1 > ref_s_pos2 && ref_e_pos1 > ref_e_pos2);
            while(ref_s_pos2<ref_s_pos1){
                noOverlapMerge(alignment, dna, qualities, nucigar, cigar_temp_unrolled, c_pos2, q_pos2, ref_s_pos2);
            }
            while(ref_s_pos2<=ref_e_pos2){
                overlapMerge(alignment,dna,qualities,nucigar,cigar_temp_unrolled,c_pos1,c_pos2,q_pos1,q_pos2,ref_s_pos2);
            }
            while(ref_s_pos2<=ref_e_pos1){
                noOverlapMerge(dna,qualities,nucigar,c_pos1,q_pos1,ref_s_pos2);
            }
        }
        this->start1 = std::min(this->start1,(unsigned int)alignment.Position+1);
        this->end1=std::max((unsigned int)alignment.GetEndPosition(),this->end1);
        this->single_end= true;
        this->cigar1 = createCigar(nucigar);
        this->sequence1=ShortDnaSequence(dna,qualities);
        this->phred_sum1=phred_sum(qualities);
        this->length_incl_deletions1 = this->sequence1.size();
        this->length_incl_longdeletions1 = this->sequence1.size();
        this->cigar1_unrolled.clear();
        for (char i : nucigar){
            this->cigar1_unrolled.push_back(i);
        }
        this->cov_pos = this->coveredPositions();
}

void AlignmentRecord::pairWith(const BamTools::BamAlignment& alignment) {
    if ((unsigned)(alignment.Position+1) > this->end1) {
        this->single_end = false;
        this->start2 = alignment.Position + 1;
        this->end2 = alignment.GetEndPosition();
        if (!(this->end2 > 0)){
            cout << this->name << " end: " << this->end2 << endl;
        }
        this->cigar2 = alignment.CigarData;
        this->sequence2 = ShortDnaSequence(alignment.QueryBases, alignment.Qualities);
        this->phred_sum2 = phred_sum(alignment.Qualities);

        this->length_incl_deletions2 = this->sequence2.size();
        this->length_incl_longdeletions2 = this->sequence2.size();
	    for (const auto& it : cigar2) {
            for (unsigned int s = 0; s < it.Length; ++s) {
              	this->cigar2_unrolled.push_back(it.Type);
            }
            if (it.Type == 'D') {
          		this->length_incl_deletions2+=it.Length;
           		if (it.Length > 1) {
           			this->length_incl_longdeletions2+=it.Length;
           		}
           	}
        }
        this->cov_pos = this->coveredPositions();
    } else if ((unsigned)alignment.GetEndPosition() < this->start1) {
        this->single_end = false;
        this->start2 = this->start1;
        this->end2 = this->end1;
        this->cigar2 = this->cigar1;
        this->sequence2 = this->sequence1;
        this->phred_sum2 = this->phred_sum1;
        this->cigar2_unrolled = this->cigar1_unrolled;
        this->length_incl_deletions2 = this->length_incl_deletions1;
        this->length_incl_longdeletions2 = this->length_incl_longdeletions1;

        this->start1 = alignment.Position + 1;
        this->end1 = alignment.GetEndPosition();
        //if (!(this->end1 > 0)){
        //    cout << this->name << " end: " << this->end1 << endl;
        //}
        this->cigar1 = alignment.CigarData;
        this->sequence1 = ShortDnaSequence(alignment.QueryBases, alignment.Qualities);
        this->phred_sum1 = phred_sum(alignment.Qualities);
        this->length_incl_deletions1 = this->sequence1.size();
        this->length_incl_longdeletions1 = this->sequence1.size();
        this->cigar1_unrolled.clear();
    	for (const auto& it : cigar1) {
            for (unsigned int s = 0; s < it.Length; ++s) {
              	this->cigar1_unrolled.push_back(it.Type);
            }
            if (it.Type == 'D') {
          		this->length_incl_deletions1+=it.Length;
           		if (it.Length > 1) {
           			this->length_incl_longdeletions1+=it.Length;
           		}
           	}
        }
        this->cov_pos = this->coveredPositions();
    }//merging of overlapping paired ends to single end reads
    else {
        this->getMergedDnaSequence(alignment);
    }
}

void AlignmentRecord::getCigarInterval(unsigned int start, unsigned int end, vector<BamTools::CigarOp>& new_cigar, const vector<BamTools::CigarOp>& original_cigar, unsigned int interval_start) {
    assert(start <= end);    
    assert(interval_start <= start);

    if (start == end) return;

    auto it = original_cigar.begin();

    while(interval_start + it->Length < start or it->Type == 'D') {
        if (it->Type != 'D') {   
            interval_start += it->Length;
        }
        it++;
    }

    if (not new_cigar.empty() and it->Type == new_cigar.back().Type) {
        if (interval_start + it->Length >= end) {
            new_cigar.back().Length += end - start;
            return;
        } else {
            new_cigar.back().Length += interval_start + it->Length - start;
        }
        interval_start += it->Length;
        it++;
    }

    while (it != original_cigar.end() and (interval_start + it->Length < end or it->Type == 'D')) {
        
        BamTools::CigarOp op;
        op.Type = it->Type;
        op.Length = min(interval_start + it->Length - start, it->Length);

        new_cigar.push_back(op);

        if (it->Type != 'D') {
            interval_start += it->Length;
        }
        it++;
    }

    if(it != original_cigar.end()) {
        BamTools::CigarOp op;
        op.Type = it->Type;
        op.Length = end - max(interval_start, start);

        new_cigar.push_back(op);
    }
}

size_t AlignmentRecord::intersectionLength(const AlignmentRecord& ap) const {
	assert(single_end == ap.single_end);
	int left = max(getIntervalStart(), ap.getIntervalStart());
	int right = min(getIntervalEnd(), ap.getIntervalEnd()) + 1;
	return max(0, right-left);
}

size_t AlignmentRecord::internalSegmentIntersectionLength(const AlignmentRecord& ap) const {
	int left = max(getInsertStart(), ap.getInsertStart());
	int right = min(getInsertEnd(), ap.getInsertEnd()) + 1;
	return max(0, right-left);
}

//Reads with overlapping paired ends have been merged by pairwith() before this method is called
std::vector<AlignmentRecord::mapValue> AlignmentRecord::coveredPositions() const{
    std::vector<AlignmentRecord::mapValue> cov_positions;
    //position in ref
    int r = this->start1;
    //position in querybases of read
    int q = 0;
    for (unsigned int i = 0; i< this->cigar1_unrolled.size(); ++i){
        char c = this->cigar1_unrolled[i];
        switch(c){
            case 'M': {
                c = this->sequence1[q];
                cov_positions.push_back({r,c,this->sequence1.qualityChar(q),phredProb(c),q,0});
                //char d = this->sequence1[q];
                ++q;
                ++r;
                break;
            }
            case 'D': {
                ++r;
                break;
            }
            case 'S':
            case 'I': {
                ++q;
                break;
            }
            case 'H':
                break;
        }
    }
    //In case of a paired end read
    if (!this->single_end){
        assert(this->start1 <= this->start2);
            //position in ref
            r = this->start2;
            //position in query bases
            q = 0;
            //in case of overlapping paired end it is declared as a single end
            //if (r <= this->end1){
            //   this->single_end = true;
            //}
            for (unsigned int i = 0; i< this->cigar2_unrolled.size(); ++i){
                char c = this->cigar2_unrolled[i];
                switch(c){
                    case 'M': {
                        c = this->sequence2[q];
                        cov_positions.push_back({r,c,this->sequence2.qualityChar(q),phredProb(c),q,1});
                        ++q;
                        ++r;
                        break;
                    }
                    case 'D': {
                        ++r;
                        break;
                    }
                    case 'S':
                    case 'I': {
                        ++q;
                        break;
                    }
                    case 'H':
                        break;
              }
           }
        }
    return cov_positions;
}


int AlignmentRecord::getPhredSum1() const {
	return phred_sum1;
}

int AlignmentRecord::getPhredSum2() const {
	assert(!single_end);
	return phred_sum2;
}

double AlignmentRecord::getProbability() const {
	return probability;
}

unsigned int AlignmentRecord::getEnd1() const {
	return end1;
}

unsigned int AlignmentRecord::getEnd2() const {
	assert(!single_end);
	return end2;
}

std::string AlignmentRecord::getName() const {
	return name;
}

unsigned int AlignmentRecord::getStart1() const {
	return start1;
}

unsigned int AlignmentRecord::getStart2() const {
	assert(!single_end);
	return start2;
}

const std::vector<BamTools::CigarOp>& AlignmentRecord::getCigar1() const {
	return cigar1;
}

const std::vector<BamTools::CigarOp>& AlignmentRecord::getCigar2() const {
	assert(!single_end);
	return cigar2;
}

const ShortDnaSequence& AlignmentRecord::getSequence1() const {
	return sequence1;
}

const ShortDnaSequence& AlignmentRecord::getSequence2() const {
	assert(!single_end);
	return sequence2;
}

unsigned int AlignmentRecord::getIntervalStart() const {
	return start1;
}

unsigned int AlignmentRecord::getIntervalEnd() const {
	if (single_end) {
		return end1;
	} else {
		return end2;
	}
}

unsigned int AlignmentRecord::getInsertStart() const {
	assert(!single_end);
	return end1 + 1;
}

unsigned int AlignmentRecord::getInsertEnd() const {
	assert(!single_end);
	return start2 - 1;
}

unsigned int AlignmentRecord::getInsertLength() const {
	assert(!single_end);
	return start2 - (end1 + 1);
}

alignment_id_t AlignmentRecord::getID() const {
	return id;
}

void AlignmentRecord::setID(alignment_id_t id) {
	this->id = id;
}

bool AlignmentRecord::isSingleEnd() const {
	return single_end;
}

bool AlignmentRecord::isPairedEnd() const {
	return !single_end;
}

std::vector<std::string> AlignmentRecord::getReadNames() const {
    vector<string> rnames;
    for (int i : this->readNames) {
        rnames.push_back((*readNameMap)[i]);
    }
	return rnames;
}

const std::vector<char> AlignmentRecord::getCigar1Unrolled() const {
	return this->cigar1_unrolled;
}
const std::vector<char> AlignmentRecord::getCigar2Unrolled() const {
	return this->cigar2_unrolled;
}
int AlignmentRecord::getLengthInclDeletions1() const {
	return this->length_incl_deletions1;
}
int AlignmentRecord::getLengthInclDeletions2() const {
	return this->length_incl_deletions2;
}
int AlignmentRecord::getLengthInclLongDeletions1() const {
	return this->length_incl_longdeletions1;
}
int AlignmentRecord::getLengthInclLongDeletions2() const {
	return this->length_incl_longdeletions2;
}

double setProbabilities(std::deque<AlignmentRecord*>& reads) {
    double read_usage_ct = 0.0;
    double mean = 1.0 / reads.size();

    for(auto&& r : reads) {
        read_usage_ct += r->getReadCount();
    }

    if (not reads.empty()) {
        read_usage_ct = max(read_usage_ct, (double) reads[0]->readNameMap->size());
    }
    double stdev = 0.0;

    for (auto&& r : reads) {
        r->probability = r->getReadCount() / read_usage_ct;
        
        stdev += (r->probability - mean)*(r->probability - mean);
    }

    return sqrt(1.0 / (reads.size() - 1) * stdev);
}

void printReads(std::ostream& outfile, std::deque<AlignmentRecord*>& reads) {
    auto comp = [](AlignmentRecord* al1, AlignmentRecord* al2) { return al1->probability > al2->probability; };
    std::sort(reads.begin(), reads.end(), comp);

    outfile.precision(5);
    outfile << std::fixed;

    for (auto&& r : reads) {
        outfile << r->name;
        if (not r->single_end) outfile << "|paired";
        outfile << "|ht_freq:" << r->probability;
        outfile << "|start1:" << r->getStart1();
        outfile << "|end1:" << r->getEnd1();
        if (not r->single_end){
            outfile << "|start2:" << r->getStart2();
            outfile << "|end2:" << r->getEnd2();
        }
        outfile << endl;
        outfile << r->sequence1;

        if (not r->single_end) {
            for(unsigned int i = r->end1; i < r->start2; i++) {
                outfile << "N";
            }
            outfile << r->sequence2;
        }
        outfile << endl;
    }
}
