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

#include "QuasispeciesEdgeCalculator.h"

 using namespace std;
 using namespace boost;

 const double QuasispeciesEdgeCalculator::FRAME_SHIFT_WEIGHT = 0.01;
 

 QuasispeciesEdgeCalculator::QuasispeciesEdgeCalculator(double Q, double edge_quasi_cutoff, double overlap, bool frameshift_merge, map<int, double>& simpson_map, double edge_quasi_cutoff_single, double overlap_single, double edge_quasi_cutoff_mixed) {
    this->Q = Q;
    this->EDGE_QUASI_CUTOFF = edge_quasi_cutoff;
    this->EDGE_QUASI_CUTOFF_SINGLE = edge_quasi_cutoff_single;
    this->EDGE_QUASI_CUTOFF_MIXED = edge_quasi_cutoff_mixed;
    this->MIN_OVERLAP_CLIQUES = overlap;
    this->MIN_OVERLAP_SINGLE = overlap_single;
    this->FRAMESHIFT_MERGE = frameshift_merge;
    this->SIMPSON_MAP = simpson_map;
}

QuasispeciesEdgeCalculator::~QuasispeciesEdgeCalculator() {
}

void QuasispeciesEdgeCalculator::getPartnerLengthRange(const AlignmentRecord& ap, unsigned int* min, unsigned int* max) const {
    assert(false);
}

bool QuasispeciesEdgeCalculator::edgeBetween(const AlignmentRecord & ap1, const AlignmentRecord & ap2) const {
    if (ap1.getName().compare(ap2.getName()) == 0) return 1;

    //cerr << ap1.getName() << " " << ap2.getName() ;

    double cutoff = 0;
    if (ap1.getName().find("Clique") != string::npos && ap2.getName().find("Clique") != string::npos) {
        cutoff = EDGE_QUASI_CUTOFF;
    } else if (ap1.getName().find("Clique") != string::npos || ap2.getName().find("Clique") != string::npos) {
        cutoff = this->EDGE_QUASI_CUTOFF_MIXED;
    } else {
        cutoff = EDGE_QUASI_CUTOFF_SINGLE;
    }
    double q = computeOverlap(ap1, ap2, cutoff);
    return q >= cutoff;
}

double QuasispeciesEdgeCalculator::computeOverlap(const AlignmentRecord & ap1, const AlignmentRecord & ap2, const double cutoff) const {

    double MIN_OVERLAP = 0;
    if (ap1.getName().find("Clique") != string::npos && ap2.getName().find("Clique") != string::npos) {
        MIN_OVERLAP = MIN_OVERLAP_CLIQUES;
    } else {
        MIN_OVERLAP = MIN_OVERLAP_SINGLE;
    }

    if (ap1.isSingleEnd() && ap2.isSingleEnd()) {
        int read_size1 = min(ap1.getEnd1() - ap1.getStart1(), ap2.getEnd1() - ap2.getStart1());
        float overlap_size1 = overlapSize(ap1.getEnd1(), ap2.getEnd1(), ap1.getStart1(), ap2.getStart1());
        if (MIN_OVERLAP <= 1) overlap_size1 /= (float) read_size1;
        if (overlap_size1 < MIN_OVERLAP) return 0;
        return singleOverlap(ap1, ap2, 1, 1, MIN_OVERLAP, cutoff);
    } else if (ap1.isSingleEnd() || ap2.isSingleEnd()) {
        if (ap1.isSingleEnd()) {
            int read_size1 = min(ap1.getEnd1() - ap1.getStart1(), ap2.getEnd1() - ap2.getStart1());
            int read_size2 = min(ap1.getEnd1() - ap1.getStart1(), ap2.getEnd2() - ap2.getStart2());
            float overlap_size1 = overlapSize(ap1.getEnd1(), ap2.getEnd1(), ap1.getStart1(), ap2.getStart1());
            float overlap_size2 = overlapSize(ap1.getEnd1(), ap2.getEnd2(), ap1.getStart1(), ap2.getStart2());
            if (MIN_OVERLAP <= 1) {
                overlap_size1 /= (float) read_size1;
                overlap_size2 /= (float) read_size2;
            }
            if ((overlap_size1 >= MIN_OVERLAP && overlap_size2 >= MIN_OVERLAP) || (overlap_size1+overlap_size2 >= MIN_OVERLAP && overlap_size1 > 0 && overlap_size2 > 0)) {
                return singleOverlap(ap1, ap2, 1, 1, MIN_OVERLAP, cutoff)*singleOverlap(ap1, ap2, 1, 2, MIN_OVERLAP, cutoff);
            }
            return 0;
        } else {
            int read_size1 = min(ap1.getEnd1() - ap1.getStart1(), ap2.getEnd1() - ap2.getStart1());
            int read_size2 = min(ap1.getEnd2() - ap1.getStart2(), ap2.getEnd1() - ap2.getStart1());

            float overlap_size1 = overlapSize(ap1.getEnd1(), ap2.getEnd1(), ap1.getStart1(), ap2.getStart1());
            float overlap_size2 = overlapSize(ap1.getEnd2(), ap2.getEnd1(), ap1.getStart2(), ap2.getStart1());
            if (MIN_OVERLAP <= 1) {
                overlap_size1 /= (float) read_size1;
                overlap_size2 /= (float) read_size2;
            }
            if ((overlap_size1 >= MIN_OVERLAP && overlap_size2 >= MIN_OVERLAP) || (overlap_size1+overlap_size2 >= MIN_OVERLAP && overlap_size1 > 0 && overlap_size2 > 0)) {
                return singleOverlap(ap1, ap2, 1, 1, MIN_OVERLAP, cutoff)*singleOverlap(ap1, ap2, 2, 1, MIN_OVERLAP, cutoff);
            }
            return 0;
        }
    } else {
        int read_size1 = min(ap1.getEnd1() - ap1.getStart1(), ap2.getEnd1() - ap2.getStart1());
        float overlap_size1 = overlapSize(ap1.getEnd1(), ap2.getEnd1(), ap1.getStart1(), ap2.getStart1());
        if (MIN_OVERLAP <= 1) overlap_size1 /= (float) read_size1;
        if (overlap_size1 < MIN_OVERLAP) return 0;

        int read_size2 = min(ap1.getEnd2() - ap1.getStart2(), ap2.getEnd2() - ap2.getStart2());
        float overlap_size2 = overlapSize(ap1.getEnd2(), ap2.getEnd2(), ap1.getStart2(), ap2.getStart2());
        if (MIN_OVERLAP <= 1) overlap_size2 /= (float) read_size2;
        if (overlap_size2 < MIN_OVERLAP) return 0;

        return singleOverlap(ap1, ap2, 1, 1, MIN_OVERLAP, cutoff)*singleOverlap(ap1, ap2, 2, 2, MIN_OVERLAP, cutoff);
    }

    return 0;
}

float QuasispeciesEdgeCalculator::overlapSize(int e1, int e2, int s1, int s2) const {
    if (s1 == s2 && e1 == e2) {
        // -----
        // -----
        return e1 - s1;
    } else if (s1 == s2) {
        // ----  AND -----
        // ----- AND ----
        if (e1 < e2) {
            return e1 - s1;
        } else {
            return e2 - s2;
        }
    } else if (e1 == e2) {
        //  ---- AND -----
        // ----- AND  ----
        if (s1 > s2) {
            return e1 - s1;
        } else {
            return e2 - s2;
        }
    } else if (e1 < e2) {
        // ----    AND   --
        //   ----  AND ------
        if (s1 < s2) {
            // ----
            //   ----
            return e1 - s2;
        } else {
            //   --
            // ------
            return e1 - s1;
        }
    } else {
        //   ---- AND ------
        // ----   AND   --
        if (s2 < s1) {
            //   ----
            // ----
            return e2 - s1;
        } else {
            // ------
            //   --
            return e2 - s2;
        }
    }
}

double QuasispeciesEdgeCalculator::singleOverlap(const AlignmentRecord & ap1, const AlignmentRecord & ap2, int strain1, int strain2, double MIN_OVERLAP, const double cutoff) const {
    int e1 = 0;
    int s1 = 0;
    ShortDnaSequence sequence1;
    vector<char> cigar1;

    int e2 = 0;
    int s2 = 0;
    ShortDnaSequence sequence2;
    vector<char> cigar2;

    if (strain1 == 1) {
        s1 = ap1.getStart1();
        e1 = ap1.getEnd1();
        sequence1 = ap1.getSequence1();
        cigar1 = ap1.getCigar1Unrolled();
    } else if (strain1 == 2) {
        s1 = ap1.getStart2();
        e1 = ap1.getEnd2();
        sequence1 = ap1.getSequence2();
        cigar1 = ap1.getCigar2Unrolled();
    }
    if (strain2 == 1) {
        s2 = ap2.getStart1();
        e2 = ap2.getEnd1();
        sequence2 = ap2.getSequence1();
        cigar2 = ap2.getCigar1Unrolled();
    } else if (strain2 == 2) {
        s2 = ap2.getStart2();
        e2 = ap2.getEnd2();
        sequence2 = ap2.getSequence2();
        cigar2 = ap2.getCigar2Unrolled();
    }

    if (s2 < s1) {
        std::swap(s1,s2);
        std::swap(e1,e2);
        std::swap(sequence1,sequence2);
        std::swap(cigar1,cigar2);
    }

    // ====
    // compute overlap
    int offset1 = 0;
    int offset2 = 0;
    int overlap = 0;
    computeOffsets(s1, s2, e1, e2, offset1, offset2, overlap);

    if (offset2 > 0) {
        //cerr << "offset 2: " << offset2 << endl;
    }
    //cerr << "offset 1: " << offset1 << endl;
    //cerr << "offset 2: " << offset2 << endl;
    //cerr << "single overlap: " << overlap << endl;
    if (overlap < MIN_OVERLAP) {
        return 0;
    }

    // ====
    //Offsets for the deletions and insertions that occured
    //upstream of the overlap. They cause index errors if not fixed.
    int offset_deletion1_ = 0;
    int offset_deletion2_ = 0;

    double overlap_probability = 0.0;
    char alphabet[] = {'A', 'C', 'G', 'T'};
    double hamming = 0;
    double total_size = 0;
    if (offset1 >= sequence1.size() || offset2 >= sequence2.size()) {
        return 0;
    }
    ////cerr << offset1 << " offset " << offset2 << endl;
    bool perfect = cutoff == 1.0;

    unsigned long j_seq = 0, j_seq2 = 0;
    unsigned long j_cigar = 0, j_cigar2 = 0;
    int j = 0, j2 = 0;
    for (;;) {
        if (cigar1[j_cigar] != 'S' && cigar2[j_cigar2] != 'S' && cigar1[j_cigar] != 'H' && cigar2[j_cigar2] != 'H') break;
        if (cigar1[j_cigar] == 'S') { j_cigar++; j_seq++; }
        if (cigar2[j_cigar2] == 'S') { j_cigar2++; j_seq2++; }
        if (cigar1[j_cigar] == 'H') { j_cigar++; }
        if (cigar2[j_cigar2] == 'H') { j_cigar2++; }
    }
    for (int j_prefix = 0;;) {
        if (j_prefix == offset1) break;
        switch (cigar1[j_cigar]) {
            case 'M': 
                if (!perfect) {
                    int j_global = s1+j;
                    if (this->SIMPSON_MAP.begin() != this->SIMPSON_MAP.end()
                        && this->SIMPSON_MAP.find(j_global) != this->SIMPSON_MAP.end()) {
                        overlap_probability += log(this->SIMPSON_MAP.at(j_global));
                        //cerr << overlap_probability << endl;
                    }
                }
                j++;
            case 'I':
                j_seq++;
            default: break;
        }
        j_cigar++;
        j_prefix++;
    }
    for (;j_cigar<cigar1.size() && j_cigar2<cigar2.size() && cigar1[j_cigar] != 'S' && cigar2[j_cigar2] != 'S' && j_seq < sequence1.size() && j_seq2 < sequence2.size();) {
        //cerr << cigar1[j_cigar] << cigar2[j_cigar2] << "\t"; 
        if (perfect) {
            if (sequence1[j_seq] != sequence2[j_seq2]) return 0;
            else {
                j_cigar++;
                j_cigar2++;
                j_seq++;
                j_seq2++;
            }
        } else if (cigar1[j_cigar] == cigar2[j_cigar2] && (cigar1[j_cigar] == 'M' || cigar1[j_cigar] == 'I')) {

            //cerr << sequence1[j_seq] << ":" << sequence2[j_seq2];
            double q_x1 = sequence1.qualityCorrect(j_seq);
            // if (q_x1 < 0.001) {
            //     q_x1 = 0.001;
            // }
            double q_x2 = sequence2.qualityCorrect(j_seq2);
            // if (q_x2 < 0.001) {
            //     q_x2 = 0.001;
            // }
            //cerr << "\t" << q_x1 << ":" << q_x2 << "\t";
            double anti_q_x1 = (1.0 - q_x1) / 3.0;
            double anti_q_x2 = (1.0 - q_x2) / 3.0;
            assert(q_x1 <= 1 && q_x2 <= 1);
            double sum = 0.0;
            sum += ((sequence1[j_seq] == alphabet[0] ? q_x1 : anti_q_x1) * (sequence2[j_seq2] == alphabet[0] ? q_x2 : anti_q_x2));
            sum += ((sequence1[j_seq] == alphabet[1] ? q_x1 : anti_q_x1) * (sequence2[j_seq2] == alphabet[1] ? q_x2 : anti_q_x2));
            sum += ((sequence1[j_seq] == alphabet[2] ? q_x1 : anti_q_x1) * (sequence2[j_seq2] == alphabet[2] ? q_x2 : anti_q_x2));
            sum += ((sequence1[j_seq] == alphabet[3] ? q_x1 : anti_q_x1) * (sequence2[j_seq2] == alphabet[3] ? q_x2 : anti_q_x2));
            overlap_probability += log(sum);
            //cerr << "\t" << overlap_probability;
            total_size++;

            //cerr << sum << "\t";
            //cerr << ((sequence1[j_seq] == alphabet[0] ? q_x1 : anti_q_x1)) << ":" << ((sequence2[j_seq2] == alphabet[0] ? q_x2 : anti_q_x2)) << "\t";
            //cerr << ((sequence1[j_seq] == alphabet[1] ? q_x1 : anti_q_x1)) << ":" << ((sequence2[j_seq2] == alphabet[1] ? q_x2 : anti_q_x2)) << "\t";
            //cerr << ((sequence1[j_seq] == alphabet[2] ? q_x1 : anti_q_x1)) << ":" << ((sequence2[j_seq2] == alphabet[2] ? q_x2 : anti_q_x2)) << "\t";
            //cerr << ((sequence1[j_seq] == alphabet[3] ? q_x1 : anti_q_x1)) << ":" << ((sequence2[j_seq2] == alphabet[3] ? q_x2 : anti_q_x2)) << "\t";

            j_cigar++;
            j_cigar2++;
            j_seq++;
            j_seq2++;

            if (j_cigar < cigar1.size() && cigar1[j_cigar] == 'M') {
                j++;
                j2++;
            } 
        } else if (cigar1[j_cigar] == cigar2[j_cigar2] && cigar1[j_cigar] == 'D') {
            j_cigar++;
            j_cigar2++;
        }else if (this->FRAMESHIFT_MERGE) {
            if (cigar2[j_cigar2] == 'M') {
                switch(cigar1[j_cigar]) {
                    case 'D':
                        //cerr << "-:" << sequence2[j_seq2];  
                        j_cigar++;
                        j_cigar2++;
                        j_seq2++;
                        j2++;
                        break;
                    case 'I':
                        //cerr << sequence1[j_seq] << ":-";  
                        j_cigar++;
                        j_seq++;
                        break;
                    default: break;
                }
            } else if (cigar1[j_cigar] == 'M') {
                switch(cigar2[j_cigar2]) {
                    case 'D':
                        //cerr << sequence1[j_seq] << ":-";  
                        j_cigar2++;
                        j_cigar++;
                        j_seq++;
                        j++;
                        break;
                    case 'I':
                        //cerr << "-:" << sequence2[j_seq2];  
                        j_cigar2++;
                        j_seq2++;
                        break;
                }
            } else {
                //cerr << "DOES NOT MATCH" << endl;
                return 0;
            }
        } else {
            return 0;
        }
        //cerr << endl;
    }
    if (!perfect) {
        if (j_cigar < cigar1.size() && cigar1[j_cigar] != 'S') {
            for (;j_cigar<cigar1.size();) {
                if (cigar1[j_cigar] == 'S') break;
                if (cigar1[j_cigar] == 'M') {
                    int j_global = s1+j;
                    if (this->SIMPSON_MAP.begin() != this->SIMPSON_MAP.end()
                        && this->SIMPSON_MAP.find(j_global) != this->SIMPSON_MAP.end()) {
                        overlap_probability += log(this->SIMPSON_MAP.at(j_global));
                        //cerr << overlap_probability << endl;
                    }
                    j++;
                }
                j_cigar++;
            }
        } else if (j_cigar2 < cigar2.size() && cigar2[j_cigar2] != 'S') {
            for (;j_cigar2<cigar2.size();) {
                if (cigar2[j_cigar2] == 'S') break;
                if (cigar2[j_cigar2] == 'M') {
                    int j_global = s2+j2;
                    if (this->SIMPSON_MAP.begin() != this->SIMPSON_MAP.end()
                        && this->SIMPSON_MAP.find(j_global) != this->SIMPSON_MAP.end()) {
                        overlap_probability += log(this->SIMPSON_MAP.at(j_global));
                        //cerr << overlap_probability << endl;
                    }
                    j2++;
                }
                j_cigar2++;
            }
        }
    }
    //cerr << pow(exp(overlap_probability),1.0/total_size) << endl;
    if (perfect) return 1;
    return pow(exp(overlap_probability),1.0/total_size);
}

void QuasispeciesEdgeCalculator::computeOffsets(int s1, int s2, int e1, int e2, int &offset1, int &offset2, int &overlap) const {
    int x_l1 = e1 - s1;
    int x_l2 = e2 - s2;
    if (s1 == s2 && e1 == e2) {
            // -----
            // -----
        overlap = e1 - s1;
        offset1 = 0;
        offset2 = 0;
    } else if (s1 == s2) {
            // ----  AND -----
            // ----- AND ----
        if (e1 < e2) {
            overlap = e1 - s1;
        } else {
            overlap = e2 - s2;
        }
        offset1 = 0;
        offset2 = 0;
    } else if (e1 == e2) {
            //  ---- AND -----
            // ----- AND  ----
        if (s1 > s2) {
            overlap = e1 - s1;
            offset1 = 0;
            offset2 = x_l2 - overlap;
        } else {
            overlap = e2 - s2;
            offset1 = x_l1 - overlap;
            offset2 = 0;
        }
    } else if (e1 < e2) {
            // ----    AND   --
            //   ----  AND ------
        if (s1 < s2) {
                // ----
                //   ----
            overlap = e1 - s2;
            offset1 = x_l1 - overlap;
            offset2 = 0;
        } else {
                //   --
                // ------
            overlap = e1 - s1;
            offset1 = 0;
            offset2 = e2 - (e2 - e1) - overlap - s2;
        }
    } else {
            //   ---- AND ------
            // ----   AND   --
        if (s2 < s1) {
                //   ----
                // ----
                //cout << "#########" << endl;
            overlap = e2 - s1;
            offset1 = 0;
            offset2 = x_l2 - overlap;
        } else {
                // ------
                //   --
            overlap = e2 - s2;
            offset1 = e1 - (e1 - e2) - overlap - s1;
            offset2 = 0;
        }
    }
}
