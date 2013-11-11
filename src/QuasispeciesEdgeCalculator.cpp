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

QuasispeciesEdgeCalculator::QuasispeciesEdgeCalculator(double Q, double edge_quasi_cutoff, double overlap, bool frameshift_merge, map<int, double>& simpson_map) {
    this->Q = Q;
    this->EDGE_QUASI_CUTOFF = edge_quasi_cutoff;
    this->MIN_OVERLAP = overlap;
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
    //cerr << ap1.getName() << "\t" << ap2.getName() << "\t"; 
    double a = computeOverlap(ap1, ap2);
    //cerr << "done" << endl;
    //    bool x = a > EDGE_QUASI_CUTOFF;
    //    if (a > 0.98)
    //    cerr << a << " " << EDGE_QUASI_CUTOFF << " " << x << endl;

    return a > EDGE_QUASI_CUTOFF;
}

std::string QuasispeciesEdgeCalculator::tail(std::string const& source, size_t const length) const {
  if (length >= source.size()) { return source; }
  return source.substr(source.size() - length);
}

double QuasispeciesEdgeCalculator::computeOverlap(const AlignmentRecord & ap1, const AlignmentRecord & ap2) const {

    double hamming = 0;
    if (ap1.isSingleEnd() && ap2.isSingleEnd()) {
        int read_size1 = min(ap1.getEnd1() - ap1.getStart1(), ap2.getEnd1() - ap2.getStart1());
        float overlap_size1 = -1;
        if (MIN_OVERLAP <= 1) {
            overlap_size1 = (float) overlapSize(ap1.getEnd1(), ap2.getEnd1(), ap1.getStart1(), ap2.getStart1()) / (float) read_size1;
        } else {
            overlap_size1 = (float) overlapSize(ap1.getEnd1(), ap2.getEnd1(), ap1.getStart1(), ap2.getStart1());
        }
        if (overlap_size1 < MIN_OVERLAP) {
            return 0;
        }
        overlap_result r = singleOverlap(ap1, ap2, 1, 1);
        double p = r.probability;
        hamming += r.hamming;
        //        p = pow(p, 1.0 / overlap_size1);
        //if (p > 0.99) cerr << p << "\t" << hamming << endl;
        p = pow(p, 1.0 / (max(ap1.getEnd1(), ap2.getEnd1()) - min(ap1.getStart1(), ap2.getStart1())));
        return p;
    } else if (ap1.isSingleEnd() || ap2.isSingleEnd()) {
        if (ap1.isSingleEnd()) {
            double p = 1.0;
            int read_size1 = min(ap1.getEnd1() - ap1.getStart1(), ap2.getEnd1() - ap2.getStart1());
            int read_size2 = min(ap1.getEnd1() - ap1.getStart1(), ap2.getEnd2() - ap2.getStart2());
            float overlap_size1 = -1;
            float overlap_size2 = -1;
            if (MIN_OVERLAP <= 1) {
                overlap_size1 = (float) overlapSize(ap1.getEnd1(), ap2.getEnd1(), ap1.getStart1(), ap2.getStart1()) / (float) read_size1;
                overlap_size2 = (float) overlapSize(ap1.getEnd1(), ap2.getEnd2(), ap1.getStart1(), ap2.getStart2()) / (float) read_size2;
            } else {
                overlap_size1 = (float) overlapSize(ap1.getEnd1(), ap2.getEnd1(), ap1.getStart1(), ap2.getStart1());
                overlap_size2 = (float) overlapSize(ap1.getEnd1(), ap2.getEnd2(), ap1.getStart1(), ap2.getStart2());
            }
            if (overlap_size1 >= MIN_OVERLAP && overlap_size2 >= MIN_OVERLAP) {
                overlap_result r1 = singleOverlap(ap1, ap2, 1, 1);
                p = r1.probability;
                hamming += r1.hamming;
                overlap_result r2 = singleOverlap(ap1, ap2, 1, 2);
                p *= r2.probability;
                hamming += r2.hamming;
                //                p = pow(p, 1.0 / (overlap_size1 + overlap_size2));
                //if (p > 0.99) cerr << p << "\t" << hamming << endl;
                p = pow(p, 1.0 / (max(ap1.getEnd1(), ap2.getEnd1()) - min(ap1.getStart1(), ap2.getStart1()) + max(ap1.getEnd1(), ap2.getEnd2()) - min(ap1.getStart1(), ap2.getStart2())));
                return p;
            }
            //            else if (overlap_size1 >= MIN_OVERLAP && overlapSize(ap1.getEnd1(), ap2.getEnd1(), ap1.getStart1(), ap2.getStart1()) > 50) {
            //                overlap_result r = singleOverlap(ap1, ap2, 1, 1);
            //                p = r.probability;
            //                p = pow(p, 1.0 / (max(ap1.getEnd1(), ap2.getEnd1()) - min(ap1.getStart1(), ap2.getStart1())));
            //                //                p *= xlr.at(0);
            //                //                hamming += xlr.at(1);
            //                //                p = pow(p, 1.0 / overlap_size1);
            //                return p;
            //            } else if (overlap_size2 >= MIN_OVERLAP && overlapSize(ap1.getEnd1(), ap2.getEnd2(), ap1.getStart1(), ap2.getStart2()) > 50) {
            //                overlap_result r = singleOverlap(ap1, ap2, 1, 2);
            //                p = r.probability;
            //                p = pow(p, 1.0 / (max(ap1.getEnd1(), ap2.getEnd2()) - min(ap1.getStart1(), ap2.getStart2())));
            //                //                p *= xlr.at(0);
            //                //                hamming += xlr.at(1);
            //                //                p = pow(p, 1.0 / overlap_size2);
            //                //if (p > 0.99) cerr << p << "\t" << hamming << endl;
            //                return p;
            //            }
            return 0;
        } else {
            double p = 1.0;
            int read_size1 = min(ap1.getEnd1() - ap1.getStart1(), ap2.getEnd1() - ap2.getStart1());
            int read_size2 = min(ap1.getEnd2() - ap1.getStart2(), ap2.getEnd1() - ap2.getStart1());

            float overlap_size1 = -1;
            float overlap_size2 = -1;
            if (MIN_OVERLAP <= 1) {
                overlap_size1 = (float) overlapSize(ap1.getEnd1(), ap2.getEnd1(), ap1.getStart1(), ap2.getStart1()) / (float) read_size1;
                overlap_size2 = (float) overlapSize(ap1.getEnd2(), ap2.getEnd1(), ap1.getStart2(), ap2.getStart1()) / (float) read_size2;
            } else {
                overlap_size1 = (float) overlapSize(ap1.getEnd1(), ap2.getEnd1(), ap1.getStart1(), ap2.getStart1());
                overlap_size2 = (float) overlapSize(ap1.getEnd2(), ap2.getEnd1(), ap1.getStart2(), ap2.getStart1());
            }

            if (overlap_size1 >= MIN_OVERLAP && overlap_size2 >= MIN_OVERLAP) {
                overlap_result r1 = singleOverlap(ap1, ap2, 1, 1);
                p *= r1.probability;
                hamming += r1.hamming;
                overlap_result r2 = singleOverlap(ap1, ap2, 2, 1);
                p *= r2.probability;
                hamming += r2.hamming;
                //                p = pow(p, 1.0 / (overlap_size1 + overlap_size2));
                //if (p > 0.99) cerr << p << "\t" << hamming << endl;
                p = pow(p, 1.0 / (max(ap1.getEnd1(), ap2.getEnd1()) - min(ap1.getStart1(), ap2.getStart1()) + max(ap1.getEnd2(), ap2.getEnd1()) - min(ap1.getStart2(), ap2.getStart1())));
                return p;
            }
            //            else if (overlap_size1 >= MIN_OVERLAP && overlapSize(ap1.getEnd2(), ap2.getEnd1(), ap1.getStart2(), ap2.getStart1()) > 50 ) {
            //                overlap_result r = singleOverlap(ap1, ap2, 1, 1);
            //                p = r.probability;
            //                p = pow(p, 1.0 / (max(ap1.getEnd1(), ap2.getEnd1()) - min(ap1.getStart1(), ap2.getStart1())));
            //                //                p *= xlr.at(0);
            //                //                hamming += xlr.at(1);
            //                //                p = pow(p, 1.0 / overlap_size1);
            //                //if (p > 0.99) cerr << p << "\t" << hamming << endl;
            //                return p;
            //            } else if (overlapSize(ap1.getEnd1(), ap2.getEnd1(), ap1.getStart1(), ap2.getStart1()) > 50 && overlap_size2 >= MIN_OVERLAP) {
            //                overlap_result r = singleOverlap(ap1, ap2, 2, 1);
            //                p = r.probability;
            //                p = pow(p, 1.0 / (max(ap1.getEnd2(), ap2.getEnd1()) - min(ap1.getStart2(), ap2.getStart1())));
            //                //                xlr = singleOverlap(ap1.getEnd2(), ap2.getEnd1(), ap1.getStart2(), ap2.getStart1(), y_clip_offset1, x_clip_offset2, y_deletions1, x_deletions2, y_insertions1, x_insertions2, ap1.getSequence2(), ap2.getSequence1());
            //                //                p *= xlr.at(0);
            //                //                hamming += xlr.at(1);
            //                //                p = pow(p, 1.0 / overlap_size2);
            //                //if (p > 0.99) cerr << p << "\t" << hamming << endl;
            //                return p;
            //            }
            return 0;
        }
    } else {
        int read_size1 = min(ap1.getEnd1() - ap1.getStart1(), ap2.getEnd1() - ap2.getStart1());
        //        float overlap_size1 = (float) overlapSize(ap1.getEnd1(), ap2.getEnd1(), ap1.getStart1(), ap2.getStart1());float overlap_size1 = -1;
        float overlap_size1 = -1;
        if (MIN_OVERLAP <= 1) {
            overlap_size1 = (float) overlapSize(ap1.getEnd1(), ap2.getEnd1(), ap1.getStart1(), ap2.getStart1()) / (float) read_size1;
        } else {
            overlap_size1 = (float) overlapSize(ap1.getEnd1(), ap2.getEnd1(), ap1.getStart1(), ap2.getStart1());
        }

        if (overlap_size1 < MIN_OVERLAP) {
            return 0;
        }
        int read_size2 = min(ap1.getEnd2() - ap1.getStart2(), ap2.getEnd2() - ap2.getStart2());
        float overlap_size2 = -1;
        if (MIN_OVERLAP <= 1) {
            overlap_size2 = (float) overlapSize(ap1.getEnd2(), ap2.getEnd2(), ap1.getStart2(), ap2.getStart2()) / (float) read_size2;
        } else {
            overlap_size2 = (float) overlapSize(ap1.getEnd2(), ap2.getEnd2(), ap1.getStart2(), ap2.getStart2());
        }
        if (overlap_size2 < MIN_OVERLAP) {
            return 0;
        }


        //        std::vector<std::string> words1;
        //        std::vector<std::string> words2;
        //        string n1 = ap1.getName();
        //        string n2 = ap2.getName();
        //        boost::split(words1, n1, boost::is_any_of(","), boost::token_compress_on);
        //        boost::split(words2, n2, boost::is_any_of(","), boost::token_compress_on);
        //        std::sort(words1.begin(), words1.end());
        //        std::sort(words2.begin(), words2.end());
        //        std::vector<std::string> v(words1.size() + words2.size());
        //
        //        std::vector<std::string>::iterator it;
        //        it = std::set_intersection(words1.begin(), words1.end(), words2.begin(), words2.end(), v.begin());
        //        v.resize(it - v.begin());
        //
        //        if (v.size() >= 5) {
        ////            cerr << "intersection " << v.size() << endl;
        //            return 1;
        //        }

        overlap_result r1 = singleOverlap(ap1, ap2, 1, 1);
        overlap_result r2 = singleOverlap(ap1, ap2, 2, 2);

        double p = r1.probability * r2.probability;

        //        int read_size3 = min(ap1.getEnd2() - ap1.getStart2(), ap2.getEnd1() - ap2.getStart1());
        //        float overlap_size2 = (float) overlapSize(ap1.getEnd2(), ap2.getEnd1(), ap1.getStart2(), ap2.getStart1()) / (float) read_size3;
        //        if (overlap_size3 > 0) {
        //            
        //        }

        p = pow(p, 1.0 / (max(ap1.getEnd1(), ap2.getEnd1()) - min(ap1.getStart1(), ap2.getStart1()) + max(ap1.getEnd2(), ap2.getEnd2()) - min(ap1.getStart2(), ap2.getStart2())));

        //        std::vector<std::string> words1;
        //        std::vector<std::string> words2;
        //        string n1 = ap1.getName();
        //        string n2 = ap2.getName();
        //        boost::split(words1, n1, boost::is_any_of(","), boost::token_compress_on);
        //        boost::split(words2, n2, boost::is_any_of(","), boost::token_compress_on);
        //        cerr << overlap_size1 << " " << words1[1] << " " << r1.probability << " | " << overlap_size2 << " " << words2[1] << " " << r2.probability << " | " << (max(ap1.getEnd1(), ap2.getEnd1()) - min(ap1.getStart1(), ap2.getStart1()) + max(ap1.getEnd2(), ap2.getEnd2()) - min(ap1.getStart2(), ap2.getStart2())) << " : " << r1.size << " : " << r2.size << " = " << p << endl;
        //        cerr << ap1.getEnd1() << " " <<  ap2.getEnd1() << " "  << ap1.getStart1() << " " << ap2.getStart1() << endl;
        hamming = r1.hamming + r2.hamming;
        //        if (p >= EDGE_QUASI_CUTOFF) cerr << (hamming / (overlap_size1 + overlap_size2)) << endl;
        return p;
    }

    return 0;
}

int QuasispeciesEdgeCalculator::overlapSize(int e1, int e2, int s1, int s2) const {
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

QuasispeciesEdgeCalculator::overlap_result QuasispeciesEdgeCalculator::singleOverlap(const AlignmentRecord & ap1, const AlignmentRecord & ap2, int strain1, int strain2) const {

    overlap_result result;
    result.hamming = 0;
    result.probability = 0;

    int e1 = 0;
    int s1 = 0;
    ShortDnaSequence sequence1;
    vector<string> cigar1;

    int e2 = 0;
    int s2 = 0;
    ShortDnaSequence sequence2;
    vector<string> cigar2;

    if (strain1 == 1) {
        s1 = ap1.getStart1();
        e1 = ap1.getEnd1();
        sequence1 = ap1.getSequence1();

        for (vector<BamTools::CigarOp>::const_iterator it = ap1.getCigar1().begin(); it != ap1.getCigar1().end(); ++it) {
            string current;
            switch (it->Type) {
                case 'S': current = "S";
                    break;
                case 'M': current = "M";
                    break;
                case 'I': current = "I";
                    break;
                case 'D': current = "D";
                    break;
                default: break;
            }
            for (int s = 0; s < it->Length; s++) cigar1.push_back(current);
        }
    } else if (strain1 == 2) {
        s1 = ap1.getStart2();
        e1 = ap1.getEnd2();
        sequence1 = ap1.getSequence2();

        for (vector<BamTools::CigarOp>::const_iterator it = ap1.getCigar2().begin(); it != ap1.getCigar2().end(); ++it) {
            string current;
            switch (it->Type) {
                case 'S': current = "S";
                    break;
                case 'M': current = "M";
                    break;
                case 'I': current = "I";
                    break;
                case 'D': current = "D";
                    break;
                default: break;
            }
            for (int s = 0; s < it->Length; s++) cigar1.push_back(current);
        }
    }
    if (strain2 == 1) {
        s2 = ap2.getStart1();
        e2 = ap2.getEnd1();
        sequence2 = ap2.getSequence1();

        for (vector<BamTools::CigarOp>::const_iterator it = ap2.getCigar1().begin(); it != ap2.getCigar1().end(); ++it) {
            string current;
            switch (it->Type) {
                case 'S': current = "S";
                    break;
                case 'M': current = "M";
                    break;
                case 'I': current = "I";
                    break;
                case 'D': current = "D";
                    break;
                default: break;
            }
            for (int s = 0; s < it->Length; s++) cigar2.push_back(current);
        }
    } else if (strain2 == 2) {
        s2 = ap2.getStart2();
        e2 = ap2.getEnd2();
        sequence2 = ap2.getSequence2();

        for (vector<BamTools::CigarOp>::const_iterator it = ap2.getCigar2().begin(); it != ap2.getCigar2().end(); ++it) {
            string current;
            switch (it->Type) {
                case 'S': current = "S";
                    break;
                case 'M': current = "M";
                    break;
                case 'I': current = "I";
                    break;
                case 'D': current = "D";
                    break;
                default: break;
            }
            for (int s = 0; s < it->Length; s++) cigar2.push_back(current);
        }
    }

    // ====
    // compute overlap 
    int x_l1 = e1 - s1;
    int x_l2 = e2 - s2;

    int offset1 = -1;
    int offset2 = -1;

    int overlap = -1;

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
    if (overlap < MIN_OVERLAP) {
        return result;
    }

    // ====
    //Offsets for the deletions and insertions that occured
    //upstream of the overlap. They cause index errors if not fixed.
    int offset_deletion1_ = 0;
    int offset_deletion2_ = 0;

    double overlap_probability = 1.0;
    char alphabet[] = {'A', 'C', 'G', 'T'};
    double hamming = 0;
    int total_size = 0;
    if (offset1 >= sequence1.size() || offset2 >= sequence2.size()) {
        return result;
    }
    for (int prefix = 1, run = 1, run2 = 1, compute_overlap = 0, j_overlap = 0, jm = 0, jm2 = 0, shift = 0, shift_del = 0, shift_ins = 0, insertion_index = 0, j = 0, j_cigar = 0, shift2 = 0, shift_del2 = 0, shift_ins2 = 0, insertion_index2 = 0, j2 = 0, j_cigar2 = 0;;) {
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
        if (compute_overlap && j < sequence1.size() && j2 < sequence2.size()) {
            //            cerr << "overlap" << endl;
            prefix = 0;
            if ((cigar1[j_cigar] == "M" && cigar2[j_cigar2] == "M") || (cigar1[j_cigar] == "I" && cigar2[j_cigar2] == "I")) {
                double q_x1 = sequence1.qualityCorrect(j);
                double q_x2 = sequence2.qualityCorrect(j2);
                assert(q_x1 <= 1 && q_x2 <= 1);
                double sum = 0.0;
                for (int x = 0; x < 4; x++) {
                    //                    cerr << q_x1 << " : " << q_x2 << endl;
                    double p_aj_x = sequence1[j] == alphabet[x] ? q_x1 : (1.0 - q_x1) / 3.0;
                    double p_bj_x = sequence2[j2] == alphabet[x] ? q_x2 : (1.0 - q_x2) / 3.0;
                    sum = sum + (p_aj_x * p_bj_x);
                }
                assert(sum <= 1);
                if (sequence1[j] != sequence2[j2]) {
                    hamming++;
                }
                //                cerr << sum << " ";
                overlap_probability *= sum;
                j_overlap++;
            } else if ((cigar1[j_cigar] != "I" && cigar2[j_cigar2] == "I")
                    || (cigar1[j_cigar] == "I" && cigar2[j_cigar2] != "I")) {
                if (j_cigar + 1 < cigar1.size() && j_cigar2 + 1 < cigar2.size()
                        && j_cigar - 1 >= 0 && j_cigar2 - 1 >= 0) {
                    if (cigar1[j_cigar - 1] == "M" && cigar2[j_cigar2 - 1] == "M"
                            && cigar1[j_cigar + 1] == "M" && cigar2[j_cigar2 + 1] == "M") {
                        //                        cerr << "YO  BITCH";
                        if (cigar1[j_cigar] == "I") {
                            jump_single2 = 1;
                        } else {
                            jump_single = 1;
                        }
                    } else {
                        return result;
                    }
                }
            } else if ((cigar1[j_cigar] != "D" && cigar2[j_cigar2] == "D")
                    || (cigar1[j_cigar] == "D" && cigar2[j_cigar2] != "D")) {
                if (j_cigar + 1 < cigar1.size() && j_cigar2 + 1 < cigar2.size()
                        && j_cigar - 1 >= 0 && j_cigar2 - 1 >= 0) {
                    if (cigar1[j_cigar - 1] == "M" && cigar2[j_cigar2 - 1] == "M"
                            && cigar1[j_cigar + 1] == "M" && cigar2[j_cigar2 + 1] == "M") {
                                               // cerr << "YO  BITCH" << endl;;
                    }
                } else {
                    return result;
                }
            } else if ((cigar1[j_cigar] == "D" && cigar2[j_cigar2] == "D") || cigar1[j_cigar] == "S" || cigar2[j_cigar2] == "S") {
            }
        }
        if (j < sequence1.size() && run && !jump_single) {
            int j_global = j - shift - shift_ins + shift_del + s1 - 1;
            if (cigar1[j_cigar] == "I") {
                insertion_index++;
                shift_ins++;
                j++;
            } else {
                insertion_index = 0;
                if (cigar1[j_cigar] == "M") {
                    if (j - shift < offset1 || j - shift > offset1 + overlap) {
                        if (this->SIMPSON_MAP.begin() != this->SIMPSON_MAP.end()
                                && this->SIMPSON_MAP.find(j_global) != this->SIMPSON_MAP.end()) {
                            overlap_probability *= this->SIMPSON_MAP.at(j_global);
                        }
                        jm++;
                    }
                    j++;
                } else if (cigar1[j_cigar] == "D") {
                    shift_del++;
                    if (jm < offset1) offset_deletion1_++;
                } else if (cigar1[j_cigar] == "S") {
                    shift++;
                    j++;
                }
            }
            j_cigar++;
        }
        if (j2 < sequence2.size() && run2 && !jump_single2) {
            int j_global = j2 - shift2 - shift_ins2 + shift_del2 + s2 - 1;
            if (cigar2[j_cigar2] == "I") {
                insertion_index2++;
                shift_ins2++;
                j2++;
            } else {
                insertion_index2 = 0;
                if (cigar2[j_cigar2] == "M") {
                    if (j2 - shift2 < offset2 || j2 - shift2 > offset2 + overlap) {
                        if (this->SIMPSON_MAP.begin() != this->SIMPSON_MAP.end()
                                && this->SIMPSON_MAP.find(j_global) != this->SIMPSON_MAP.end()) {
                            overlap_probability *= this->SIMPSON_MAP.at(j_global);
                        }
                        jm2++;
                    }
                    j2++;
                } else if (cigar2[j_cigar2] == "D") {
                    shift_del2++;
                    if (jm2 < offset2) offset_deletion2_++;
                } else if (cigar2[j_cigar2] == "S") {
                    shift2++;
                    j2++;
                }
            }
            j_cigar2++;
        }
        if (j == sequence1.size() && j2 == sequence2.size()) {
            break;
        }
        total_size++;
    }

    result.hamming = hamming;
    result.probability = overlap_probability;
    result.size = total_size;
    return result;
}