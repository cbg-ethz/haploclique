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
#include <limits>
#include <map>
#include <boost/unordered_set.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>   
#include "CliqueWriter.h"

using namespace std;
using namespace boost;

CliqueWriter::CliqueWriter(ostream& os, VariationCaller* variation_caller, std::ostream* indel_os, const ReadGroups* read_groups, bool multisample, bool output_all, double fdr_threshold, bool verbose, int min_coverage) : os(os), variation_caller(variation_caller), read_groups(read_groups) {
    this->indel_os = indel_os;
    this->significant_ins_count = -1;
    this->significant_del_count = -1;
    this->total_count = 0;
    this->total_insertion_cliques = 0;
    this->total_deletion_cliques = 0;
    this->output_all = output_all;
    this->fdr_threshold = fdr_threshold;
    this->finished = false;
    this->read_list_os = 0;
    this->verbose = verbose;
    this->multisample = multisample;
    this->clique_count = 0;
    if (multisample) {
        assert(read_groups != 0);
    }
    this->min_coverage = min_coverage;
    this->single_count = 0;
    this->paired_count = 0;
}

CliqueWriter::~CliqueWriter() {
    if (read_list_os != 0) {
        for (size_t i = 0; i < clique_list.size(); ++i) {
            assert(clique_list[i].reads != 0);
            delete clique_list[i].reads;
        }
    }
}

void CliqueWriter::enableReadListOutput(std::ostream& os) {
    assert(!output_all);
    assert(total_count == 0);
    assert(read_list_os == 0);
    read_list_os = &os;
}

void CliqueWriter::callVariation(const vector<const AlignmentRecord*>& pairs, size_t coverage, clique_stats_t* stats) {
    bool DEBUG = 0;
    bool merge = 1;
    bool work = 0;
    bool problem = 0;
    do {
        assert(stats != 0);
        VariationCaller::additional_stats_t vc_stats;
        if (variation_caller != 0) {
            stats->variation = variation_caller->call(pairs.begin(), pairs.end(), &vc_stats);
        } else {
            stats->variation = Variation();
        }
        stats->min_coverage_user = min_coverage;
        stats->clique_size = pairs.size();
        stats->start = vc_stats.insert_start;
        stats->end = vc_stats.insert_end;
        stats->length = vc_stats.insert_length;
        stats->diff = vc_stats.diff;
        stats->total_weight = vc_stats.total_weight;
        stats->coverage = coverage;
        stats->pvalue_corr = min(1.0, stats->variation.getPValue() * pow(2.0, static_cast<int> (stats->coverage)));
        // assume every clique to be not significant until the final FDR control
        stats->is_significant = false;

        stats->readnames = new vector<std::string>();

        vector<std::string> duplicates;
        //Ignore cliques with less than 5 reads.
        //We consider them as trash :)
        //if (pairs.end() - pairs.begin() < min_coverage) return;

        //Compute the borders of both super-reads
        //and save them in stats
        vector<const AlignmentRecord*>::const_iterator it = pairs.begin();
        if (DEBUG) cerr << pairs.end() - pairs.begin() << "\t";

        stats->clique_number = clique_count++;

        int single_end_count = 0;
        int paired_end_count = 0;
        if (pairs.end() - pairs.begin() < min_coverage) {
            return;
        }
        for (; it != pairs.end(); ++it) {
            const AlignmentRecord& ap = **it;
            //        if (std::find(stats->readnames->begin(), stats->readnames->end(), ap.getName())!=stats->readnames->end()) {
            //            duplicates.push_back(ap.getName());
            //            continue;
            //        }
            if (ap.isSingleEnd()) {
                single_end_count += 1;
            } else {
                paired_end_count += 1;
            }
            if (ap.getName().find("-+-") != std::string::npos) {
                string bla = ap.getName();
                vector<string> fields;
                boost::split(fields, bla, is_any_of("-+-"));
                for (int i = 1; i < fields.size(); i++) stats->readnames->push_back(fields[i]);
                //            cerr << fields[0] << "\t" << fields[1] << "\t" << fields[2] << "\t" << fields[3]  << "\t" << fields.size() << "\t" << bla << endl;
            } else {
                //            cerr << "NOOOO" << endl;
                stats->readnames->push_back(ap.getName());
            }
            //cerr << ap.getName() << "\t";
            vector<BamTools::CigarOp>::const_iterator it_cigar = ap.getCigar1().begin();
            int deletions1 = 0;
            int shift1 = 0;
            int shift_end1 = 0;
            int read_begin1 = 1;
            for (; it_cigar != ap.getCigar1().end(); ++it_cigar) {
                if (it_cigar->Type == 'S') {
                    if (read_begin1) {
                        shift1 += it_cigar->Length;
                        read_begin1 = 0;
                    } else {
                        shift_end1 += it_cigar->Length;
                    }
                }
                if (it_cigar->Type == 'D') {
                    deletions1 += it_cigar->Length;
                }
            }
            if (stats->window_start1 == -1 || int(ap.getStart1() + shift1) < stats->window_start1) {
                stats->window_start1 = ap.getStart1() + shift1;
            }
            if (stats->window_end1 == -1 || int(ap.getStart1() + ap.getSequence1().size() + deletions1 - shift_end1) > stats->window_end1) {
                stats->window_end1 = ap.getStart1() + ap.getSequence1().size() + deletions1 - shift_end1;
            }
            if (ap.isPairedEnd()) {
                it_cigar = ap.getCigar2().begin();
                int deletions2 = 0;
                int shift2 = 0;
                int shift_end2 = 0;
                int read_begin2 = 1;
                for (; it_cigar != ap.getCigar2().end(); ++it_cigar) {
                    if (it_cigar->Type == 'S') {
                        if (read_begin1) {
                            shift2 += it_cigar->Length;
                            read_begin2 = 0;
                        } else {
                            shift_end2 += it_cigar->Length;
                        }
                    }
                    if (it_cigar->Type == 'D') {
                        deletions2 += it_cigar->Length;
                    }
                }
                if (stats->window_start2 == -1 || int(ap.getStart2() + shift2) < stats->window_start2) {
                    stats->window_start2 = ap.getStart2() + shift2;
                }
                if (stats->window_end2 == -1 || int(ap.getStart2() + ap.getSequence2().size() + deletions2 - shift_end2) > stats->window_end2) {
                    stats->window_end2 = ap.getStart2() + ap.getSequence2().size() + deletions2 - shift_end2;
                }
            }
        }
        this->single_count += single_end_count;
        this->paired_count += paired_end_count;
        //2D array represents nucleotide distribution for super-read
        int alignment_length1 = stats->window_end1 - stats->window_start1;
        int alignment1[alignment_length1][5];
        double phred1[alignment_length1][5];
        for (int j = 0; j < alignment_length1; j++) {
            for (int k = 0; k < 5; k++) {
                alignment1[j][k] = 0;
                phred1[j][k] = 0;
            }
        }
        if (DEBUG) cerr << "computing alignment, length: " << alignment_length1 << endl;
        //For each read pair, save the first read in alignment
        it = pairs.begin();
        for (; it != pairs.end(); ++it) {
            const AlignmentRecord& ap_alignment = **it;
            //        if (std::find(duplicates.begin(), duplicates.end(), ap_alignment.getName())!=duplicates.end()) {
            //            continue;
            //        }
            vector<BamTools::CigarOp>::const_iterator it_cigar = ap_alignment.getCigar1().begin();

            int alignment_index = ap_alignment.getStart1() - stats->window_start1;
            int sequence_index = 0;

            for (; it_cigar != ap_alignment.getCigar1().end(); ++it_cigar) {
                int cigar_length = it_cigar->Length;
                if (it_cigar->Type == 'S') {
                    sequence_index += cigar_length;
                }
                if (it_cigar->Type == 'D') {
                    for (int k = 0; k < cigar_length; k++) {
                        alignment1[alignment_index++][4]++;
                    }
                }
                if (it_cigar->Type == 'I' || it_cigar->Type == 'M') {
                    if (it_cigar->Type == 'I' && cigar_length == 1 && merge) {
                        //                        cerr << "TRIPPING LAICA BAWS" << endl;
                        alignment_index++;
                    } else {
                        for (int k = 0; k < cigar_length; k++) {
                            int base = shortenBase(ap_alignment.getSequence1()[sequence_index++]);

                            assert(alignment_index < alignment_length1);
                            if (base != -1) {
                                phred1[alignment_index][base] += ap_alignment.getSequence1().qualityCorrectLog(sequence_index - 1);
                                alignment1[alignment_index][base]++;
                            } else {
//                                cerr << "w00000t" << endl;
//                                problem = 1;
                            }
                            alignment_index++;
                        }
                    }
                }
            }
        }
        if (DEBUG) {
            for (int k = 0; k < 5; k++) {
                for (int j = 0; j < alignment_length1; j++) {
                    cerr << alignment1[j][k] << "\t";
                }
                cerr << endl;
            }
            cerr << "\n" << endl;
        }

        //Compute maximum coverage
        stats->maximum_coverage1 = 0;
        int coverage1[alignment_length1];
        for (int j = 0; j < alignment_length1; j++) {
            int sum = 0;
            for (int k = 0; k < 5; k++) {
                sum += alignment1[j][k];
            }
            coverage1[j] = sum;
            if (sum > stats->maximum_coverage1) {
                stats->maximum_coverage1 = sum;
            }
        }

        int min_coverage_local1 = stats->maximum_coverage1 < min_coverage ? 1 : min_coverage;

        //Majority vote for super-read assembly
        bool prefix1 = 1;
        int end1 = 0;
        for (int j = alignment_length1 - 1; j >= 0; j--) {
            if (coverage1[j] >= min_coverage_local1) {
                end1 = j;
                stats->window_end1 -= alignment_length1 - j;
                break;
            }
        }

        for (int j = 0; j < alignment_length1; j++) {
            if (j == end1) break;

            if (prefix1) {
                if (coverage1[j] >= min_coverage_local1) {
                    prefix1 = 0;
                    stats->window_start1 += j;
                } else {
                    continue;
                }
            }

            int local_max = 0;
            int local_base = -1;

            for (int k = 0; k < 5; k++) {
                if (alignment1[j][k] > local_max) {
                    local_max = alignment1[j][k];
                    local_base = k;
                }
            }
            if (local_base < 4) {

                //            if (local_base == -1) {
                //                cerr << "BREAK: " << j << " " << alignment_length1 << endl;
                //
                //                it = pairs.begin();
                //                for (; it != pairs.end(); ++it) {
                //                    const AlignmentRecord& ap_alignment = **it;
                //                    cerr << ap_alignment.getStart1() - stats->window_start1 << " " << ap_alignment.getEnd1() - stats->window_start1 << endl;
                //                }
                //                cerr << "---" << endl;
                //                stats->window_end1 - (alignment_length1 - j - 1);
                //                break;
                //            }
                char base = expandBase(local_base);
                if (base != 'N') {
                    stats->consensus_string1 += base;
                    int p = 33 + round(-10 * log10(pow(10, phred1[j][local_base])));
                    if (p > 126) {
                        p = 126;
                    } else if (p < 33) {
                        p = 126;
                    }
                    stats->phred_string1 += p;

                    stats->coverage_string1 += ('0' + floor(local_max * 10 / stats->maximum_coverage1));
                } else {
                    stats->window_end1--;
                }
            } else if (local_base == 4) {
                stats->window_end1--;
            }
        }

        //same stuff as for first read, this time for second read
        if (paired_end_count > 0) {
            int alignment_length2 = stats->window_end2 - stats->window_start2;
            int alignment2[alignment_length2][5];
            double phred2[alignment_length2][5];
            for (int j = 0; j < alignment_length2; j++) {
                for (int k = 0; k < 5; k++) {
                    alignment2[j][k] = 0;
                    phred2[j][k] = 0;
                }
            }
            if (DEBUG) cerr << "computing alignment, length: " << alignment_length2 << endl;
            it = pairs.begin();
            for (; it != pairs.end(); ++it) {
                const AlignmentRecord& ap_alignment = **it;
                if (ap_alignment.isSingleEnd()) continue;
                vector<BamTools::CigarOp>::const_iterator it_cigar = ap_alignment.getCigar2().begin();

                int alignment_index = ap_alignment.getStart2() - stats->window_start2;
                int sequence_index = 0;

                for (; it_cigar != ap_alignment.getCigar2().end(); ++it_cigar) {
                    if (DEBUG) cerr << "   " << sequence_index << "\t" << alignment_index << "\t" << ap_alignment.getName() << endl;
                    if (DEBUG) cerr << ap_alignment.getSequence2() << endl;
                    if (DEBUG) cerr << it_cigar->Length << it_cigar->Type << endl;
                    int cigar_length = it_cigar->Length;
                    if (it_cigar->Type == 'S') {
                        sequence_index += cigar_length;
                        if (DEBUG) cerr << "S: " << sequence_index << "\t" << alignment_index << endl;
                    }
                    if (it_cigar->Type == 'D') {
                        for (int k = 0; k < cigar_length; k++) {
                            alignment2[alignment_index++][4]++;
                        }
                        if (DEBUG) cerr << "D: " << sequence_index << "\t" << alignment_index << endl;
                    }
                    if (it_cigar->Type == 'I' || it_cigar->Type == 'M') {
                        if (it_cigar->Type == 'I' && cigar_length == 1 && merge) {
                            //                            cerr << "TRIPPING LAICA BAWS" << endl;
                            alignment_index++;
                        } else {
                            for (int k = 0; k < cigar_length; k++) {
                                if (DEBUG) cerr << "MI:" << sequence_index << "\t" << alignment_index << "\t" << alignment_length2 << endl;
                                int base = shortenBase(ap_alignment.getSequence2()[sequence_index++]);
                                if (DEBUG) cerr << base << endl;

                                assert(alignment_index < alignment_length2);
                                if (base != -1) {
                                    phred2[alignment_index][base] += ap_alignment.getSequence2().qualityCorrectLog(sequence_index - 1);
                                    alignment2[alignment_index][base]++;
                                } else {
//                                    cerr << "w000t2" << endl;
//                                    problem = 1;
                                }

                                alignment_index++;
                            }
                            if (DEBUG) cerr << "MI done:" << sequence_index << "\t" << alignment_index << endl;
                        }
                    }
                }
            }
            if (DEBUG) {
                for (int k = 0; k < 5; k++) {
                    for (int j = 0; j < alignment_length2; j++) {
                        cerr << alignment2[j][k] << "\t";
                    }
                    cerr << endl;
                }
                cerr << "\n" << endl;
            }

            //Compute maximum coverage
            stats->maximum_coverage2 = 0;
            int coverage2[alignment_length2];
            for (int j = 0; j < alignment_length2; j++) {
                int sum = 0;
                for (int k = 0; k < 5; k++) {
                    sum += alignment2[j][k];
                }
                coverage2[j] = sum;
                if (sum > stats->maximum_coverage2) {
                    stats->maximum_coverage2 = sum;
                }
            }

            int min_coverage_local2 = stats->maximum_coverage2 < min_coverage ? 1 : min_coverage;

            bool prefix2 = 1;
            int end2 = 0;
            for (int j = alignment_length2 - 1; j >= 0; j--) {
                if (coverage2[j] >= min_coverage_local2) {
                    end2 = j;
                    stats->window_end2 -= alignment_length2 - j;
                    break;
                }
            }

            for (int j = 0; j < alignment_length2; j++) {
                if (j == end2) break;

                if (prefix2) {
                    if (coverage2[j] >= min_coverage_local2) {
                        prefix2 = 0;
                        stats->window_start2 += j;
                    } else {
                        continue;
                    }
                }

                int local_max = 0;
                int local_base = -1;

                for (int k = 0; k < 5; k++) {
                    if (alignment2[j][k] > local_max) {
                        local_max = alignment2[j][k];
                        local_base = k;
                    }
                }
                if (local_base < 4) {
                    //                if (local_base == -1) {
                    //                    stats->window_end2 - (alignment_length2 - j - 1);
                    //                    break;
                    //                }
                    char base = expandBase(local_base);
                    if (base != 'N') {
                        stats->consensus_string2 += base;
                        int p = 33 + round(-10 * log10(pow(10, phred2[j][local_base])));
                        if (p > 126) {
                            p = 126;
                        } else if (p < 33) {
                            p = 126;
                        }
                        stats->phred_string2 += p;
                        stats->coverage_string2 += ('0' + floor(local_max * 10 / stats->maximum_coverage2));
                    } else {
                        stats->window_end2--;
                    }
                } else if (local_base == 4) {
                    stats->window_end2--;
                }
            }
        }
        if (problem || overlapSize(stats)) {
            stats->window_end1 == -1;
            stats->window_start1 == -1;
            stats->window_end2 == -1;
            stats->window_start1 == -1;
            stats->consensus_string1 = "";
            stats->consensus_string2 = "";
            stats->phred_string1 = "";
            stats->phred_string2 = "";
            stats->clique_number = 0;
            this->single_count -= single_end_count;
            this->paired_count -= paired_end_count;
            merge = 0;
            problem = 0;
            cerr << endl << "PROBLEM" << endl;
            continue;
        }
        if (!merge)
            cerr << stats->window_start1 << " " << merge << "|" << endl;
    } while (work);
}

bool CliqueWriter::overlapSize(clique_stats_t* stats) const {
    //    cerr << "Overlap: " << stats->window_start1 << "\t" << stats->window_end1 << "\t" << stats->window_start2 << "\t" << stats->window_end2 << endl;
    if (stats->window_start2 > 0) {
        int s1 = stats->window_start1;
        int s2 = stats->window_start2;
        int e1 = s1 + stats->consensus_string1.size();
        int e2 = s2 + stats->consensus_string2.size();
        if (s2 > e1) {
            //        cerr << "a" << endl;
            // -----
            //       -----
            return 0;
        }
        if (s1 > e2) {
            //        cerr << "b" << endl;
            //       -----
            // -----
            string cons_tmp = stats->consensus_string1;
            string phred_tmp = stats->phred_string1;
            int end_tmp = stats->window_end1;
            int start_tmp = stats->window_start1;

            stats->consensus_string1 = stats->consensus_string2;
            stats->phred_string1 = stats->phred_string2;
            stats->window_end1 = stats->window_end2;
            stats->window_start1 = stats->window_start2;

            stats->consensus_string2 = cons_tmp;
            stats->phred_string2 = phred_tmp;
            stats->window_end2 = end_tmp;
            stats->window_start2 = start_tmp;
            return 0;
        }
        if (s1 == s2 && e1 == e2) {
            // -----
            // -----
        } else if (s1 == s2) {
            // ----  AND -----
            // ----- AND ----
            if (e1 < e2) {
                //            cerr << "c" << endl;
                stats->consensus_string1 = stats->consensus_string2;
                stats->phred_string1 = stats->phred_string2;
                stats->window_end1 = stats->window_end2;
            } else {
            }
        } else if (e1 == e2) {
            //  ---- AND -----
            // ----- AND  ----
            if (s1 > s2) {
                //            cerr << "d" << endl;
                stats->consensus_string1 = stats->consensus_string2;
                stats->phred_string1 = stats->phred_string2;
                stats->window_start1 = stats->window_start2;
            } else {
            }
        } else if (e1 < e2) {
            // ----    AND   --
            //   ----  AND ------
            if (s1 < s2) {
                if (e1 - s2 >= 1) {
                    //                    cerr << "FIXING";
                    string cons = equalStrings(stats->consensus_string1, stats->consensus_string2);
                    if (cons.size() > 0) {
                        for (int i = stats->consensus_string2.size()-(cons.size() - stats->consensus_string1.size()); i < stats->consensus_string2.size(); i++) {
                            stats->phred_string1 += stats->phred_string2.at(i);
                        }
                        stats->consensus_string1 = cons;
                        stats->window_end1 = stats->window_start1 + cons.size();
                        if (stats->phred_string1.size() != cons.size()) {
                            cerr << stats->phred_string1.size() << " " << cons.size() << endl;
                        }
                        assert(stats->phred_string1.size() == cons.size());
                        //                    } else {
                        //                        cerr << "FAILED " << e1 - s2;
                    }
                    //                    cerr << endl;
                }
            } else {
                //            cerr << "f" << endl;
                //   --
                // ------
                stats->consensus_string1 = stats->consensus_string2;
                stats->phred_string1 = stats->phred_string2;
                stats->window_end1 = stats->window_end2;
                stats->window_start1 = stats->window_start2;
            }
        } else {
            //   ---- AND ------
            // ----   AND   --
            if (s2 < s1) {
                //            cerr << "g" << endl;
                //   ----
                // ----

                int offset = 0;
                bool fixing = 1;
                int tries = 0;
                do {
                    if (e1 - s2 >= 1) {
                        if (stats->consensus_string1[e2 - s1 - 1 + offset ] != stats->consensus_string1[stats->consensus_string2.size() - 1]) {
                            //                            cerr << stats->consensus_string2[e1 - s2 - 1 + offset ] << " " << stats->consensus_string2[e1 - s2 - 2 + offset ] << endl;
                            //                            cerr << stats->consensus_string1[stats->consensus_string1.size() - 1] << " " << stats->consensus_string1[stats->consensus_string1.size() - 2] << endl;
                            offset++;
                            fixing = 1;
                        } else {
                            fixing = 0;
                        }
                    } else {
                        fixing = 0;
                    }
                    if (tries > 1) {
                        cerr << "ERROR OVERLAP2" << endl;
                        cerr << stats->consensus_string2 << endl;
                        cerr << stats->consensus_string1 << endl;
                        cerr << e2 - s1 << endl;
                        offset = 0;
                        break;
                    }
                    tries++;
                } while (fixing);

                stats->consensus_string1 = stats->consensus_string2 + stats->consensus_string1.substr(e2 - s1 + offset);
                stats->phred_string1 = stats->phred_string2 + stats->phred_string1.substr(e2 - s1 + offset);
                stats->window_start1 = stats->window_start2;
            } else {
                // ------
                //   --
            }
        }
    }

    stats->window_end2 = 0;
    stats->window_start2 = 0;
    stats->phred_string2 = "";
    stats->consensus_string2 = "";

    int split = 0;
    for (int i = 0; i < stats->consensus_string1.size(); i++) {
        if (stats->consensus_string1.at(i) == 'N') {
            split = 1;
            break;
        }
    }
    if (split) {
        int prefix = 1;
        int suffix = 0;
        string tmp_string1;
        string tmp_phred1;
        int tmp_end = 0;
        for (int i = 0; i < stats->consensus_string1.size(); i++) {
            if (prefix) {
                if (stats->consensus_string1.at(i) == 'N') {
                    split = 1;
                    prefix = 0;
                    tmp_end = stats->window_start1 + i;
                } else {
                    tmp_string1 += stats->consensus_string1.at(i);
                    tmp_phred1 += stats->phred_string1.at(i);
                }
            } else {
                if (prefix == 0 && stats->consensus_string1.at(i) != 'N') {
                    suffix = 1;
                    stats->window_start2 = stats->window_start1 + i;
                    stats->window_end2 = stats->window_end1;
                }
                if (suffix) {
                    stats->consensus_string2 += stats->consensus_string1.at(i);
                    stats->phred_string2 += stats->phred_string1.at(i);

                    if (stats->consensus_string1.at(i) == 'N') {
                        cerr << "FUCK " << i << " " << stats->consensus_string1.size() << endl;
                        return 1;
                    }
                }
            }
        }
        stats->window_end1 = stats->window_start1 + stats->consensus_string1.size();
        stats->consensus_string1 = tmp_string1;

        cerr << "SPLIT" << endl;
    }

    return 0;
}

string CliqueWriter::equalStrings(string s1, string s2) const {
    int e1 = s1.size();
    int e2 = s1.size();
    for (int i = 0; i < e2; i++) {
        bool fit = 1;
        int del = 0;
        for (int j = 0; j <= i; j++) {
            if (s1.at(e1 - 1 - i + j) != s2.at(j)) {
                fit = 0;
                break;
            }
            if (s1.at(e1 - 1 - i + j) == '-') {
                del++;
            }
        }
        if (fit && i - del >= 9) {
            return s1 + s2.substr(i + 1);
        }
    }
    return "";
}

char CliqueWriter::expandBase(int shortBase) {
    switch (shortBase) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        case 4:
            return '-';
        default:
            return 'N';
    }
}

int CliqueWriter::shortenBase(char base) {
    switch (base) {
        case 65:
            return 0;
        case 67:
            return 1;
        case 71:
            return 2;
        case 84:
            return 3;
        case 45:
            return 4;
        case 78:
            return -1;
        default:
            cerr << "Not aware of base: " << base << endl;
            exit(0);
            break;
    }
}

void CliqueWriter::add(std::auto_ptr<Clique> clique) {
    assert(!finished);
    clique_stats_t stats;
    auto_ptr<vector<const AlignmentRecord*> > all_pairs = clique->getAllAlignments();
    assert(all_pairs->size() == clique->size());
    if (!multisample) {
        callVariation(*all_pairs, clique->totalCenterCoverage(), &stats);
    } else {
        // if read groups are available, run tests for all possible combinations of read groups.
        auto_ptr<vector<size_t> > readgroup_wise_coverages = clique->readGroupWiseCoverage();
        assert(readgroup_wise_coverages.get() != 0);
        assert(readgroup_wise_coverages->size() == read_groups->size());
        // create sample wise coverage counts by summing over read groups
        auto_ptr<vector<size_t> > sample_wise_coverages(new vector<size_t>(read_groups->sampleCount(), 0));
        for (size_t i = 0; i < readgroup_wise_coverages->size(); ++i) {
            (*sample_wise_coverages)[read_groups->readGroupIndexToSampleIndex(i)] += readgroup_wise_coverages->at(i);
        }
        // compute number of non-empty subsets of samples
        size_t n = (1 << read_groups->sampleCount()) - 1;
        vector<clique_stats_t> subset_wise_clique_stats(n, clique_stats_t());
        double best_pvalue = numeric_limits<double>::infinity();
        size_t best_pvalue_index = 0;
        std::vector<sample_subset_wise_stats_t> sample_subsset_wise_stats;
        // for every sample, compute probabilty that all alignments are wrong
        vector<double> sample_wrong_probs(read_groups->size(), 1.0);
        for (size_t i = 0; i < all_pairs->size(); ++i) {
            assert(all_pairs->at(i) != 0);
            int rg = all_pairs->at(i)->getReadGroup();
            assert(rg >= 0);
            assert(rg < (int) read_groups->size());
            int sample = read_groups->readGroupIndexToSampleIndex(rg);
            assert(sample >= 0);
            assert(sample < (int) read_groups->sampleCount());
            sample_wrong_probs[sample] *= 1.0 - all_pairs->at(i)->getWeight();
        }
        // i is interpreted as a bitvector telling which read groups are active
        for (size_t i = 1; i <= n; ++i) {
            clique_stats_t* current_stats = &subset_wise_clique_stats[i - 1];
            vector<const AlignmentRecord*> pairs;
            for (size_t j = 0; j < all_pairs->size(); ++j) {
                int rg = all_pairs->at(j)->getReadGroup();
                assert(rg != -1);
                int sample = read_groups->readGroupIndexToSampleIndex(rg);
                assert(sample >= 0);
                assert(sample < (int) read_groups->sampleCount());
                if ((i & (1 << sample)) != 0) {
                    pairs.push_back(all_pairs->at(j));
                }
            }
            if (pairs.size() == 0) {
                sample_subsset_wise_stats.push_back(sample_subset_wise_stats_t(-1.0, 0.0));
                continue;
            }
            double p = 1.0;
            size_t coverage = 0;
            for (size_t j = 0; j < read_groups->sampleCount(); ++j) {
                if ((i & (1 << j)) != 0) {
                    coverage += sample_wise_coverages->at(j);
                    p *= 1.0 - sample_wrong_probs[j];
                } else {
                    p *= sample_wrong_probs[j];
                }
            }
            callVariation(pairs, coverage, current_stats);
            if (current_stats->pvalue_corr < best_pvalue) {
                best_pvalue = current_stats->pvalue_corr;
                best_pvalue_index = i - 1;
            }
            sample_subsset_wise_stats.push_back(sample_subset_wise_stats_t(current_stats->pvalue_corr, p));
        }
        stats = subset_wise_clique_stats[best_pvalue_index];
        stats.sample_subset_wise_stats.assign(sample_subsset_wise_stats.begin(), sample_subsset_wise_stats.end());
        stats.best_sample_combination = best_pvalue_index + 1;
        stats.sample_wise_stats.clear();
        for (size_t i = 0; i < read_groups->sampleCount(); ++i) {
            stats.sample_wise_stats.push_back(sample_wise_stats_t(sample_wise_coverages->at(i)));
        }
        // gather read group information
        for (size_t i = 0; i < all_pairs->size(); ++i) {
            assert(all_pairs->at(i) != 0);
            int rg = all_pairs->at(i)->getReadGroup();
            assert(rg >= 0);
            assert(rg < (int) read_groups->size());
            int sample = read_groups->readGroupIndexToSampleIndex(rg);
            assert(sample >= 0);
            assert(sample < (int) stats.sample_wise_stats.size());
            stats.sample_wise_stats[sample].add(1, all_pairs->at(i)->getWeight());
        }
    }
    switch (stats.variation.getType()) {
        case Variation::INSERTION:
            total_insertion_cliques += 1;
            break;
        case Variation::DELETION:
            total_deletion_cliques += 1;
            break;
        default:
            break;
    }
    total_count += 1;
    bool passed_fdr = stats.pvalue_corr <= fdr_threshold;
    // if read list is to be printed, we need to store the required
    // information
    if ((read_list_os != 0) && passed_fdr) {
        // retrieve all alignments associated the current clique
        auto_ptr<vector<const AlignmentRecord*> > alignments = clique->getAllAlignments();
        stats.reads = new vector<alignment_id_t>();
        for (size_t i = 0; i < alignments->size(); ++i) {
            const AlignmentRecord& ap = *alignments->at(i);
            alignment_id_t as;
            // does readname already exist?
            readname_to_index_bimap_t::left_const_iterator it = readname_to_index.left.find(ap.getName());
            if (it != readname_to_index.left.end()) {
                as.read_name_idx = it->second;
            } else {
                as.read_name_idx = readname_to_index.left.size();
                readname_to_index.insert(readname_to_index_bimap_t::value_type(ap.getName(), as.read_name_idx));
            }
            as.record_nr = ap.getRecordNr();
            stats.reads->push_back(as);
        }
    }
    if (output_all) {
        os << stats << endl;
        if (passed_fdr && (indel_os != 0)) {
            clique_list.push_back(stats);
        }
    } else {
        if (passed_fdr) {
            clique_list.push_back(stats);
        }
    }
}

void CliqueWriter::writeReadlist() {
    assert(read_list_os != 0);
    // determine the set of reads to be written.
    // "first" gives index of clique, "second" gives read index within this clique,
    // i.e. clique_list[first].reads->at(second) gives an alignment
    typedef pair<size_t, size_t> alignment_index_t;
    // Comparator to sort reads according to their name
    typedef map<size_t, vector<alignment_index_t>, readname_comparator_t> read_to_clique_idx_t;
    read_to_clique_idx_t read_to_clique_idx(readname_comparator_t(*this));
    for (size_t i = 0; i < clique_list.size(); ++i) {
        const clique_stats_t& stats = clique_list[i];
        if (!stats.is_significant) continue;
        assert(stats.reads != 0);
        for (size_t j = 0; j < stats.reads->size(); ++j) {
            alignment_id_t& aln_stats = stats.reads->at(j);
            // does read already exist in read_to_clique_idx?
            read_to_clique_idx_t::iterator it = read_to_clique_idx.find(aln_stats.read_name_idx);
            if (it != read_to_clique_idx.end()) {
                // if yes, just add variant index to list
                it->second.push_back(make_pair(i, j));
            } else {
                read_to_clique_idx[aln_stats.read_name_idx] = vector<alignment_index_t>();
                read_to_clique_idx[aln_stats.read_name_idx].push_back(make_pair(i, j));
            }
        }
    }
    read_to_clique_idx_t::const_iterator it = read_to_clique_idx.begin();
    for (; it != read_to_clique_idx.end(); ++it) {
        const string& read_name = readname_to_index.right.at(it->first);
        (*read_list_os) << read_name;
        for (size_t i = 0; i < it->second.size(); ++i) {
            const alignment_index_t& aln_idx = it->second[i];
            const alignment_id_t& as = clique_list[aln_idx.first].reads->at(aln_idx.second);
            (*read_list_os) << "\t" << aln_idx.first << "," << as.record_nr;
        }
        (*read_list_os) << endl;
    }
}

void CliqueWriter::finish() {
    finished = true;
    if (indel_os == 0) return;
    sort(clique_list.begin(), clique_list.end(), clique_stats_comp_t());
    // perform benjamini-hochberg procedure, i.e. determine number significant
    // insertions and deletions
    significant_ins_count = 0;
    significant_del_count = 0;
    size_t insertion_count = 0;
    size_t deletion_count = 0;
    for (size_t i = 0; i < clique_list.size(); ++i) {
        switch (clique_list[i].variation.getType()) {
            case Variation::INSERTION:
                insertion_count += 1;
                clique_list[i].fdr_level = clique_list[i].pvalue_corr * total_insertion_cliques / insertion_count;
                if (clique_list[i].fdr_level <= fdr_threshold) {
                    significant_ins_count = i + 1;
                    clique_list[i].is_significant = true;
                }
                break;
            case Variation::DELETION:
                deletion_count += 1;
                clique_list[i].fdr_level = clique_list[i].pvalue_corr * total_deletion_cliques / deletion_count;
                if (clique_list[i].fdr_level <= fdr_threshold) {
                    significant_del_count = i + 1;
                    clique_list[i].is_significant = true;
                }
                break;
            default:
                assert(false);
        }
    }
    for (size_t i = 0; i < clique_list.size(); ++i) {
        const clique_stats_t& stats = clique_list[i];
        if (stats.is_significant) {
            (*indel_os) << stats.variation << " " << stats.total_weight << " " << stats.clique_size << " " << stats.coverage << " " << stats.start << " " << stats.end << " " << stats.length << " " << stats.diff << " " << stats.variation.getPValue() << " " << stats.pvalue_corr;
            if (stats.fdr_level>=0) {
                (*indel_os) << " " << stats.fdr_level;
            }
            (*indel_os) << endl;
        }
    }
    if (read_list_os != 0) {
        writeReadlist();
    }
}

ostream& operator<<(ostream& os, const CliqueWriter::clique_stats_t& stats) {
    if (stats.readnames->size() >= stats.min_coverage_user) {
        if (stats.consensus_string1.size() > 0) {
            os << "@Clique1_" << stats.clique_number << "-+-" << stats.readnames->at(0);
            for (size_t i = 1; i < stats.readnames->size(); ++i) {
                os << "," << stats.readnames->at(i);
            }
            os << endl;
            os << stats.consensus_string1 << endl;
            //start length size coverage reads
            os << "+" << stats.window_start1 << " " << stats.window_end1 - stats.window_start1 << " " << stats.readnames->size() << " " << stats.coverage_string1 << endl;
            os << stats.phred_string1 << endl;
            if (stats.consensus_string2.size() > 0) {
                os << "@Clique2_" << stats.clique_number << "-+-" << stats.readnames->at(0);
                for (size_t i = 1; i < stats.readnames->size(); ++i) {
                    os << "," << stats.readnames->at(i);
                }
                os << endl;
                os << stats.consensus_string2 << endl;
                os << "+" << stats.window_start2 << " " << stats.window_end2 - stats.window_start2 << " " << stats.readnames->size() << " " << stats.coverage_string2 << endl;
                os << stats.phred_string2 << endl;
            }
        } else if (stats.consensus_string2.size() > 0) {
            os << "@Clique1_" << stats.clique_number << "-+-" << stats.readnames->at(0);
            for (size_t i = 1; i < stats.readnames->size(); ++i) {
                os << "," << stats.readnames->at(i);
            }
            os << endl;
            os << stats.consensus_string2 << endl;
            os << "+" << stats.window_start2 << " " << stats.window_end2 - stats.window_start2 << " " << stats.readnames->size() << " " << stats.coverage_string2 << endl;
            os << stats.phred_string2 << endl;
        }

    }
    return os;
}
