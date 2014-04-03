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
    if (ap1.getName().compare(ap2.getName()) == 0) {
        return 1;
    }

        // cerr << ap1.getName() << "\t" << ap2.getName() << endl;
    

    double cutoff = 0;
    if (ap1.getName().find("Clique") != string::npos
        && ap2.getName().find("Clique") != string::npos) {
        cutoff = EDGE_QUASI_CUTOFF;
    } else if (ap1.getName().find("Clique") != string::npos
        || ap2.getName().find("Clique") != string::npos) {
        cutoff = this->EDGE_QUASI_CUTOFF_MIXED;
    } else {
        cutoff = EDGE_QUASI_CUTOFF_SINGLE;
    }
    double q = computeOverlap(ap1, ap2, cutoff);
    // if ((ap1.getName().compare(s1) == 0 || ap2.getName().compare(s1) == 0 )
    //     && (ap1.getName().compare(s2) == 0 || ap2.getName().compare(s2) == 0 )) {
        // cerr << endl << "Q: " << q << endl;
    // }
    return q >= cutoff;
}

string QuasispeciesEdgeCalculator::tail(string const& source, size_t const length) const {
  if (length >= source.size()) { return source; }
  return source.substr(source.size() - length);
}

bool QuasispeciesEdgeCalculator::is_disjoint(const std::set<std::string> &set1, const std::set<std::string> &set2) const {

    if(set1.empty() || set2.empty()) return true;

    std::set<std::string>::const_iterator
        it1 = set1.begin(),
        it1End = set1.end();
    std::set<std::string>::const_iterator
        it2 = set2.begin(),
        it2End = set2.end();

    while(it1 != it1End && it2 != it2End)
    {
        int c = (*it1).compare(*it2);
        if(c == 0) { return false; }
        if(c < 0) { it1++; }
        else { it2++; }
    }
    return true;
}

double QuasispeciesEdgeCalculator::computeOverlap(const AlignmentRecord & ap1, const AlignmentRecord & ap2, const double cutoff) const {

    double MIN_OVERLAP = 0;
    if (ap1.getName().find("Clique") != string::npos
        && ap2.getName().find("Clique") != string::npos) {
        MIN_OVERLAP = MIN_OVERLAP_CLIQUES;
    } else {
        MIN_OVERLAP = MIN_OVERLAP_SINGLE;
    }

    if (ap1.isSingleEnd() && ap2.isSingleEnd()) {
        int read_size1 = min(ap1.getEnd1() - ap1.getStart1(), ap2.getEnd1() - ap2.getStart1());
        float overlap_size1 = overlapSize(ap1.getEnd1(), ap2.getEnd1(), ap1.getStart1(), ap2.getStart1());
        if (MIN_OVERLAP <= 1) overlap_size1 /= (float) read_size1;
        if (overlap_size1 < MIN_OVERLAP) return 0;
        // if (!is_disjoint(ap1.getReadNames(),ap2.getReadNames())) return 1;
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
                // if (!is_disjoint(ap1.getReadNames(),ap2.getReadNames())) return 1;
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
                // if (!is_disjoint(ap1.getReadNames(),ap2.getReadNames())) return 1;
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

        // if (!is_disjoint(ap1.getReadNames(),ap2.getReadNames())) return 1;
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

        for (vector<BamTools::CigarOp>::const_iterator it = ap1.getCigar1().begin(); it != ap1.getCigar1().end(); ++it) {
            for (int s = 0; s < it->Length; ++s) cigar1.push_back(it->Type);
        }
    } else if (strain1 == 2) {
        s1 = ap1.getStart2();
        e1 = ap1.getEnd2();
        sequence1 = ap1.getSequence2();

        for (vector<BamTools::CigarOp>::const_iterator it = ap1.getCigar2().begin(); it != ap1.getCigar2().end(); ++it) {
            for (int s = 0; s < it->Length; ++s) cigar1.push_back(it->Type);
        }
    }
    if (strain2 == 1) {
        s2 = ap2.getStart1();
        e2 = ap2.getEnd1();
        sequence2 = ap2.getSequence1();

        for (vector<BamTools::CigarOp>::const_iterator it = ap2.getCigar1().begin(); it != ap2.getCigar1().end(); ++it) {
            for (int s = 0; s < it->Length; ++s) cigar2.push_back(it->Type);
        }
    } else if (strain2 == 2) {
        s2 = ap2.getStart2();
        e2 = ap2.getEnd2();
        sequence2 = ap2.getSequence2();

        for (vector<BamTools::CigarOp>::const_iterator it = ap2.getCigar2().begin(); it != ap2.getCigar2().end(); ++it) {
            for (int s = 0; s < it->Length; ++s) cigar2.push_back(it->Type);
        }
    }

    // ====
    // compute overlap
    int offset1 = 0;
    int offset2 = 0;
    int overlap = 0;
    computeOffsets(s1, s2, e1, e2, offset1, offset2, overlap);

        // cerr << "single overlap: " << overlap << endl;
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
            // cerr << "out of here" << endl;
        return 0;
    }
    /*if (offset1 > 0) {
        offset1 -= 1;
    }
    if (offset2 > 0) {
        offset2 -= 1;
    }*/
    //cerr << offset1 << " offset " << offset2 << endl;
    bool perfect = 0;
    if (ap1.getName().find("Clique") != string::npos
        && ap2.getName().find("Clique") != string::npos
        && cutoff == 1.0) {
        perfect = 1;
    }
    for (int j_compare = 0, j2_compare = 0, prefix = 1, run = 1, run2 = 1, compute_overlap = 0, j_overlap = 0, jm = 0, jm2 = 0,
            shift_ins_prefix = 0, shift_ins_prefix2 = 0,
            shift = 0, shift_del = 0, shift_del_prefix = 0, shift_ins = 0, insertion_index = 0, j = 0, j_cigar = 0,
            shift2 = 0, shift_del2 = 0, shift_del_prefix2 = 0, shift_ins2 = 0, insertion_index2 = 0, j2 = 0, j_cigar2 = 0;;) {
        int j_tmp = j - shift_del_prefix + shift_ins_prefix;
        int j2_tmp = j2 - shift_del_prefix2 + shift_ins_prefix2;
        int jump_single = 0;
        int jump_single2 = 0;
        if ((j - shift == offset1 || offset1 == 0) && prefix) {
            run = 0;
        }
        if ((j2 - shift2 == offset2 || offset2 == 0) && prefix) {
            run2 = 0;
        }
        if (!run && !run2) {
            compute_overlap = 1;
            run = 1;
            run2 = 1;
        }
        if (j_overlap >= overlap) {
            compute_overlap = 0;
        }
        bool skip = 1;
        if (compute_overlap && j < sequence1.size() && j2 < sequence2.size() && j_tmp < sequence1.size() && j2_tmp < sequence2.size()) {
            skip = 0;
            prefix = 0;
            if (cigar1[j_cigar] == cigar2[j_cigar2]) {
                switch (cigar1[j_cigar]) {
                    case 'M':
                    case 'I':
                        if (!perfect) {
                            double q_x1 = sequence1.qualityCorrect(j_tmp);
                            double q_x2 = sequence2.qualityCorrect(j2_tmp);
                            double anti_q_x1 = (1.0 - q_x1) / 3.0;
                            double anti_q_x2 = (1.0 - q_x2) / 3.0;

                            assert(q_x1 <= 1 && q_x2 <= 1);
                            double sum = 0.0;
                            sum += ((sequence1[j_tmp] == alphabet[0] ? q_x1 : anti_q_x1) * (sequence2[j2_tmp] == alphabet[0] ? q_x2 : anti_q_x2));
                            sum += ((sequence1[j_tmp] == alphabet[1] ? q_x1 : anti_q_x1) * (sequence2[j2_tmp] == alphabet[1] ? q_x2 : anti_q_x2));
                            sum += ((sequence1[j_tmp] == alphabet[2] ? q_x1 : anti_q_x1) * (sequence2[j2_tmp] == alphabet[2] ? q_x2 : anti_q_x2));
                            sum += ((sequence1[j_tmp] == alphabet[3] ? q_x1 : anti_q_x1) * (sequence2[j2_tmp] == alphabet[3] ? q_x2 : anti_q_x2));
                            overlap_probability += log(sum);
                        } else {
                            if (sequence1[j_tmp] != sequence2[j2_tmp]) {
                                return 0;
                            }
                        }
                        j_overlap++;
                        total_size++;
                        default: break;
                }
            } else if (this->FRAMESHIFT_MERGE && (cigar1[j_cigar] == 'I' || cigar1[j_cigar] == 'D')  && cigar2[j_cigar2] == 'M' && j_cigar + 1 < cigar1.size() && j_cigar - 1 >= 0 && cigar1[j_cigar - 1] == 'M' && cigar1[j_cigar + 1] == 'M') {
                if (cigar1[j_cigar] == 'I') {
                    jump_single2 = 1;
                }
                //otherwise it's a deletion and just ignore that base
            } else if (this->FRAMESHIFT_MERGE && cigar1[j_cigar] == 'M' && (cigar2[j_cigar2] == 'I' || cigar2[j_cigar2] == 'D') && j_cigar2 + 1 < cigar2.size() && j_cigar2 - 1 >= 0 && cigar2[j_cigar2 - 1] == 'M' && cigar2[j_cigar2 + 1] == 'M') {
                if (cigar2[j_cigar2] == 'I') {
                    jump_single = 1;
                }
                //otherwise it's a deletion and just ignore that base
            } else {
                return 0;
            }
        }
        if (j < sequence1.size() && run && !jump_single) {
            skip = 0;
            int j_global = j - shift - shift_ins + shift_del + s1 - 1;
            //cerr << "jump1 " << j_global << endl;
            if (cigar1[j_cigar] == 'I') {
                if (prefix) shift_ins_prefix++;
                insertion_index++;
                shift_ins++;
                j++;
            } else {
                insertion_index = 0;
                if (cigar1[j_cigar] == 'M') {
                    if (j - shift < offset1 || j - shift > offset1 + overlap) {
                        if (this->SIMPSON_MAP.begin() != this->SIMPSON_MAP.end()
                            && this->SIMPSON_MAP.find(j_global) != this->SIMPSON_MAP.end()) {
                            overlap_probability += log(this->SIMPSON_MAP.at(j_global));
                        }
                        jm++;
                    }
                    j++;
                } else if (cigar1[j_cigar] == 'D') {
                    if (prefix) shift_del_prefix++;
                    shift_del++;
                    if (jm < offset1) offset_deletion1_++;
                } else if (cigar1[j_cigar] == 'S') {
                    shift++;
                    j++;
                }
            }
            j_cigar++;
        }
        if (j2 < sequence2.size() && run2 && !jump_single2) {
            skip = 0;
            int j_global = j2 - shift2 - shift_ins2 + shift_del2 + s2 - 1;
            //cerr << "jump2 " << j_global << endl;
            if (cigar2[j_cigar2] == 'I') {
                if (prefix) shift_ins_prefix2++;
                insertion_index2++;
                shift_ins2++;
                j2++;
            } else {
                insertion_index2 = 0;
                if (cigar2[j_cigar2] == 'M') {
                    if (j2 - shift2 < offset2 || j2 - shift2 > offset2 + overlap) {
                        if (this->SIMPSON_MAP.begin() != this->SIMPSON_MAP.end()
                            && this->SIMPSON_MAP.find(j_global) != this->SIMPSON_MAP.end()) {
                            overlap_probability += log(this->SIMPSON_MAP.at(j_global));
                        }
                        jm2++;
                    }
                    j2++;
                } else if (cigar2[j_cigar2] == 'D') {
                    //cerr << "D2" << endl;
                    if (prefix) shift_del_prefix2++;
                    shift_del2++;
                    if (jm2 < offset2) offset_deletion2_++;
                } else if (cigar2[j_cigar2] == 'S') {
                    shift2++;
                    j2++;
                }
            }
            j_cigar2++;
        }
        if (j == sequence1.size() && j2 == sequence2.size()) {
            break;
        }
        if (skip) {
            break;
        }
    }
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
