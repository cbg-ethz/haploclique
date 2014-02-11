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
#include <boost/algorithm/string/compare.hpp>
#include "CliqueWriter.h"
#include <sstream>

 using namespace std;
 using namespace boost;

 CliqueWriter::CliqueWriter(ostream& os, VariationCaller* variation_caller, std::ostream* indel_os, const ReadGroups* read_groups, bool multisample, bool output_all, double fdr_threshold, bool verbose, int min_coverage, bool frameshift_merge, std::string suffix) : os(os), variation_caller(variation_caller), read_groups(read_groups) {
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
    this->single_skipped_count = 0;
    this->FRAMESHIFT_MERGE=frameshift_merge;
    this->suffix = suffix;
    if (!this->suffix.empty()) {
        this->suffix += "_";
    }
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
    //cerr << endl << "CALL VARIATION " << pairs.size() << endl;
    //cerr.flush();
    bool DEBUG = 0;
    bool merge = this->FRAMESHIFT_MERGE;
    bool problem = 0;
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
    stats->pvalue_corr = min(1.0, stats->variation.getPValue() * pow(1.2, static_cast<int> (stats->coverage)));
        // assume every clique to be not significant until the final FDR control
    stats->is_significant = false;

    stats->readnames = new vector<std::string>();

    vector<const AlignmentRecord*>::const_iterator it = pairs.begin();

    int single_end_count = 0;
    int paired_end_count = 0;
    ofstream singles;

    for (; it != pairs.end(); ++it) {
        const AlignmentRecord& ap = **it;
        for (int x=0;x<ap.getReadNames().size();++x) {
            stats->readnames->push_back(ap.getReadNames()[x]); 
        }
        stats->clique_size_weighted += ap.getReadCount();
    }
    if (stats->clique_size_weighted < min_coverage) {
        return;
    }
    stats->clique_number = clique_count++;
    bool contains_clique = 0;
    it = pairs.begin();
    for (; it != pairs.end(); ++it) {
        const AlignmentRecord& ap = **it;
        contains_clique |= ap.getName().find("Clique") != std::string::npos;
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
            paired_end_count++;
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
    //2D array represents nucleotide distribution for super-read
    int alignment_length1 = stats->window_end1 - stats->window_start1;
    int alignment1[alignment_length1][5];
    double phred1[alignment_length1][5];

    int alignment_length2 = 0;
    if (paired_end_count > 0) {
        alignment_length2 = stats->window_end2 - stats->window_start2+1;
    }
    int alignment2[alignment_length2][5];
    double phred2[alignment_length2][5];
    int max_length = std::max(alignment_length1,alignment_length2);

    for (int i = 0; i < max_length; ++i) {
        if (i < alignment_length1) {
            std::fill_n(alignment1[i], 5, 0);
            std::fill_n(phred1[i], 5, 0);
        }
        if (i < alignment_length2) {
            std::fill_n(alignment2[i], 5, 0);
            std::fill_n(phred2[i], 5, 0);
        }
    }
    it = pairs.begin();
    for (; it != pairs.end(); ++it) {
        const AlignmentRecord& ap = **it;
        vector<BamTools::CigarOp>::const_iterator it_cigar = ap.getCigar1().begin();

        int alignment_index = ap.getStart1() - stats->window_start1;
        int sequence_index = 0;

        for (; it_cigar != ap.getCigar1().end(); ++it_cigar) {
            int cigar_length = it_cigar->Length;
            switch (it_cigar->Type) {
                case 'S':
                    sequence_index += cigar_length;
                    break;
                case 'D':
                    if (merge && cigar_length == 1) {
                        alignment_index++;
                    } else {
                        for (int k = 0; k < cigar_length; k++) {
                            alignment1[alignment_index++][4]+=ap.getReadCount();
                        }
                    }
                    break;
                case 'I':
                    if (cigar_length == 1 && merge) {
                        sequence_index++;
                        break;
                    }
                case 'M':
                    for (int k = 0; k < cigar_length; k++) {
                        int base = shortenBase(ap.getSequence1()[sequence_index++]);

                        assert(alignment_index < alignment_length1);
                        if (base != -1) {
                            phred1[alignment_index][base] += ap.getSequence1().qualityCorrectLog(sequence_index - 1);
                            alignment1[alignment_index][base]+=ap.getReadCount();
                        }

                        alignment_index++;
                    }
                    break;
            }
        }

        if (ap.isSingleEnd()) continue;
        it_cigar = ap.getCigar2().begin();

        alignment_index = ap.getStart2() - stats->window_start2;
        sequence_index = 0;

        for (; it_cigar != ap.getCigar2().end(); ++it_cigar) {
            int cigar_length = it_cigar->Length;
            switch (it_cigar->Type) {
                case 'S':
                    sequence_index += cigar_length;
                    break;
                case 'D':
                    if (merge && cigar_length == 1) {
                        alignment_index++;
                    } else {
                        for (int k = 0; k < cigar_length; k++) {
                            alignment2[alignment_index++][4]+=ap.getReadCount();
                        }
                    }
                    break;
                case 'I':
                    if (cigar_length == 1 && merge) {
                        sequence_index++;
                        break;
                    }
                case 'M':
                    for (int k = 0; k < cigar_length; k++) {
                        int base = shortenBase(ap.getSequence2()[sequence_index++]);

                        assert(alignment_index < alignment_length2);
                        if (base != -1) {
                            phred2[alignment_index][base] += ap.getSequence2().qualityCorrectLog(sequence_index - 1);
                            alignment2[alignment_index][base]+=ap.getReadCount();
                        }

                        alignment_index++;
                    }
                    break;
            }
        }
    }



    //Compute maximum coverage
    int coverage1[alignment_length1];
    int coverage2[alignment_length2];
    for (int i = 0; i < max_length; ++i) {
        if (i < alignment_length1) {
            int sum = 0;
            for (int k = 0; k < 5; k++) {
                sum += alignment1[i][k];
            }
            coverage1[i] = sum;
        }

        if (i < alignment_length2) {
            int sum = 0;
            for (int k = 0; k < 5; k++) {
                sum += alignment2[i][k];
            }
            coverage2[i] = sum;
        }
    }

    //Majority vote for super-read assembly
    bool prefix1 = 1;
    int end1 = 0;
    for (int j = alignment_length1 - 1; j >= 0; j--) {
        if (coverage1[j] >= min_coverage) {
            end1 = j;
            stats->window_end1 -= alignment_length1 - j;
            break;
        }
    }

    for (int j = 0; j < alignment_length1; j++) {
        if (j > end1) { break; }

        if (prefix1) {
            if (coverage1[j] >= min_coverage) {
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
            } else {
                stats->window_end1--;
            }
        } else if (local_base == 4) {
            stats->window_end1--;
        }
    }

    //same stuff as for first read, this time for second read
    if (paired_end_count > 0) {
        this->paired_count += 1;
        
        bool prefix2 = 1;
        int end2 = 0;
        for (int j = alignment_length2 - 1; j >= 0; j--) {
            if (coverage2[j] >= min_coverage) {
                end2 = j;
                stats->window_end2 -= alignment_length2 - j;
                break;
            }
        }

        for (int j = 0; j < alignment_length2; j++) {
            if (j > end2) break;

            if (prefix2) {
                if (coverage2[j] >= min_coverage) {
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
                } else {
                    stats->window_end2--;
                }
            } else if (local_base == 4) {
                stats->window_end2--;
            }
        }
    } else {
        this->single_count += 1;
    }
    if (problem || overlapSize(stats)) {
        stats->window_end1 = -1;
        stats->window_start1 = -1;
        stats->window_end2 = -1;
        stats->window_start1 = -1;
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
            //continue;
    }
    //if (!merge) cerr << stats->window_start1 << " " << merge << "|" << endl;
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
                    }
                }
            } else {
                stats->consensus_string1 = stats->consensus_string2;
                stats->phred_string1 = stats->phred_string2;
                stats->window_end1 = stats->window_end2;
                stats->window_start1 = stats->window_start2;
            }
        } else {
            //   ---- AND ------
            // ----   AND   --
            if (s2 < s1) {
                string cons = equalStrings(stats->consensus_string2, stats->consensus_string1);
                if (cons.size() > 0) {
                    for (int i = stats->consensus_string1.size()-(cons.size() - stats->consensus_string2.size()); i < stats->consensus_string1.size(); i++) {
                        stats->phred_string2 += stats->phred_string1.at(i);
                    }
                    stats->phred_string1 = stats->phred_string2;
                    stats->consensus_string1 = cons;
                    stats->window_start1 = stats->window_start2;
                    stats->window_end1 = stats->window_start1 + cons.size();
                    if (stats->phred_string1.size() != cons.size()) {
                        cerr << stats->phred_string1.size() << " " << cons.size() << endl;
                    }
                    assert(stats->phred_string1.size() == cons.size());
                }
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

    callVariation(*all_pairs, clique->totalCenterCoverage(), &stats);

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

    if (stats.clique_size_weighted >= stats.min_coverage_user) {
        std::string key;
        key = boost::lexical_cast<std::string>(stats.window_start1)+stats.consensus_string1;
        if (!stats.consensus_string2.empty()) {
            key += boost::lexical_cast<std::string>(stats.window_start2)+stats.consensus_string2;
        }
        if (fastq_map.find(key) == fastq_map.end()) {
            fastq_map[key] = fastq_entry();
            fastq_map[key].name = boost::lexical_cast<std::string>(stats.clique_number);
            fastq_map[key].read_names = new std::set<std::string>();
            for (size_t i = 0; i < stats.readnames->size(); ++i) {
                    fastq_map[key].read_names->insert(stats.readnames->at(i));
            }
            fastq_map[key].pos_1 = stats.window_start1;
            fastq_map[key].seq_1 = stats.consensus_string1;
            fastq_map[key].phreds_1.push_back(stats.phred_string1);
            if (!stats.consensus_string2.empty()) {
                fastq_map[key].pos_2 = stats.window_start2;
                fastq_map[key].seq_2 = stats.consensus_string2;
                fastq_map[key].phreds_2.push_back(stats.phred_string2);
            }
        } else {
            for (size_t i = 0; i < stats.readnames->size(); ++i) {
                fastq_map[key].read_names->insert(stats.readnames->at(i));
            }
            fastq_map[key].phreds_1.push_back(stats.phred_string1);
            if (!stats.consensus_string2.empty()) {
                fastq_map[key].phreds_2.push_back(stats.phred_string2);
            }
        }
    }
    if (stats.window_start1 > 0) {
        printout(stats.window_start1);
    }
    if (output_all) {
        //os << stats << endl;
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

void CliqueWriter::printout(int pos_1) {
    ofstream fr1;
    fr1.open ("data_cliques_paired_R1.fastq", ios::out | ios::app);
    ofstream fr2;
    fr2.open ("data_cliques_paired_R2.fastq", ios::out | ios::app);
    ofstream fs;
    fs.open ("data_cliques_single.fastq", ios::out | ios::app);
    ofstream fc;
    fc.open ("data_clique_to_reads.tsv", ios::out | ios::app);

    std::map<std::string,fastq_entry>::iterator it;
    it = fastq_map.begin();
    for (;it!=fastq_map.end();) {
        fastq_entry f;
        f = it->second;
        if (f.pos_1+10 < pos_1 || pos_1 == -1) {
            if (f.read_names->size() >= this->min_coverage) {
                ofstream* out1;
                if (f.seq_2.empty()) {
                    out1 = &fs;
                } else {
                    out1 = &fr1;
                }
                *out1 << "@Clique_" << this->suffix << f.name << endl;
                std::set<std::string>::iterator set_it;
                set_it = f.read_names->begin();
                fc << "Clique_" << this->suffix << f.name << "\t" << *set_it;
                ++set_it;
                for (;set_it!=f.read_names->end();++set_it) {
                     fc << "," << *set_it;
                }
                fc << endl;
                *out1 << f.seq_1 << endl;
                *out1 << "+" << f.pos_1 << endl;
                *out1 << f.phreds_1[0] << endl;
                if (!f.seq_2.empty()) {
                    fr2 << "@Clique_" << this->suffix << f.name;
                    fr2 << endl;
                    int last = f.seq_2.size()-1;
                    for(int i = last; i >= 0; --i) {
                        switch(f.seq_2[i]) {
                            case 'G':
                                fr2 << 'C';
                                break;
                            case 'A':
                                fr2 << 'T';
                                break;
                            case 'T':
                                fr2 << 'A';
                                break;
                            case 'C':
                                fr2 << 'G';
                                break;
                            case 'N':
                                fr2 << 'N';
                                break;
                            default:
                            cerr << "UNKNOWN BASE " << f.seq_2[i] << endl;
                                break;
                        }
                    }
                    fr2 << endl;
                    fr2 << "+" << f.pos_2 << endl;
                    fr2 << f.phreds_2[0] << endl;
                }
            }
            fastq_map.erase(it++);
        } else {
            ++it;
        }
    }
    fr1.close();
    fr2.close();
    fs.close();
    fc.close();
}

void CliqueWriter::finish() {
    printout(-1);
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
    cerr << "Thou shall not make use of me" << endl;
    return os;
}
