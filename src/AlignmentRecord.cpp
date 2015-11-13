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
       	}
    }
}

AlignmentRecord::AlignmentRecord(unique_ptr<vector<const AlignmentRecord*>>& alignments, unsigned int clique_id) : cigar1_unrolled(), cigar2_unrolled() {
    deque<pair<int, int>> interval;
    vector<ShortDnaSequence> sequences;
    vector<vector<BamTools::CigarOp>> cigars;
    int ct = 0;

    this->readNameMap = (*alignments)[0]->readNameMap;

    for (auto al : *alignments) {
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
    assert(length_ct == this->length_incl_deletions1); //DEBUG


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
        assert(length_ct == this->length_incl_deletions2); //DEBUG
    }
}

void AlignmentRecord::mergeSequences(std::deque<std::pair<int, int>> intervals, vector<ShortDnaSequence> to_merge, vector<vector<BamTools::CigarOp>> cigars) {

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
                if (basemap[c] > max_val) max_key = c;
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

void AlignmentRecord::pairWith(const BamTools::BamAlignment& alignment) {
    this->single_end = false;
    if (alignment.Position > this->start1) {
        this->start2 = alignment.Position + 1;
       this->end2 = alignment.GetEndPosition();
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
    } else {
        this->start2 = this->start1;
        this->end2 = this->end1;
        this->cigar2 = this->cigar1;
        this->sequence2 = this->sequence1;
        this->phred_sum2 = this->phred_sum1;
        this->start1 = alignment.Position + 1;
        this->end1 = alignment.GetEndPosition();
        this->cigar1 = alignment.CigarData;
        this->sequence1 = ShortDnaSequence(alignment.QueryBases, alignment.Qualities);
        this->phred_sum1 = phred_sum(alignment.Qualities);

        this->cigar2_unrolled = this->cigar1_unrolled;
        this->length_incl_deletions2 = this->length_incl_deletions1;
        this->length_incl_longdeletions2 = this->length_incl_longdeletions2;

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
           	}
        }
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

AlignmentRecord::covmap AlignmentRecord::coveredPositions(){
    AlignmentRecord::covmap positions;
    if (single_end){

    }
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
        outfile << "|ht_freq:" << r->probability << endl;

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
