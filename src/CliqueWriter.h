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

#ifndef CLIQUEWRITER_H_
#define CLIQUEWRITER_H_

#include <iostream>
#include <vector>
#include <boost/bimap.hpp>
#include <boost/unordered_map.hpp>
#include <set>

#include "CliqueCollector.h"
#include "Variation.h"
#include "VariationCaller.h"

class CliqueWriter : public CliqueCollector {
private:

    /** Type to store all information needed to uniquely identify an alignment (pair). */
    typedef struct alignment_id_t {
        /** Index of the read name as stored in "readname_to_index". */
        size_t read_name_idx;
        /** Read-specific index of the alignment pair (i.e. if read maps uniquely,
         *  pair_nr is always 0). */
        size_t record_nr;

        alignment_id_t() : read_name_idx(0), record_nr(0) {
        }
    } alignment_id_t;

    /** Information to be stored per sample (if sample information is available). */
    typedef struct sample_wise_stats_t {
        int absolute_support;
        float expected_support;
        int coverage;

        void add(int absolute_support, float expected_support) {
            this->absolute_support += absolute_support;
            this->expected_support += expected_support;
        }

        sample_wise_stats_t() : absolute_support(0), expected_support(0.0), coverage(0) {
        }

        sample_wise_stats_t(int coverage) : absolute_support(0), expected_support(0.0), coverage(coverage) {
        }
    } sample_wise_stats_t;

    /** Information to be stored for all subset of samples (if sample information is available). */
    typedef struct sample_subset_wise_stats_t {
        /** P-value (after local Bonferroni correction) for the present variation when only considering
         *  the alignments for these read groups. */
        double p_value;
        /** Probability for this subset of read groups, that is, the probabilty that at least one alignment
         *  from this group is correct, while all alignments from other read groups are wrong. */
        double probability;

        sample_subset_wise_stats_t() : p_value(-1.0), probability(-1.0) {
        }

        sample_subset_wise_stats_t(double p_value, double probability) : p_value(p_value), probability(probability) {
        }
    } sample_subset_wise_stats_t;

    /** Type to store statistics on a clique. */
    typedef struct clique_stats_t {
        /** The (putative) variation derived from this clique. */
        Variation variation;
        /** Sum of weights of all alignments in the clique. */
        double total_weight;
        /** Number of alignments. */
        size_t clique_size;
        /** Coverage at the center of the intersection of all insert segments.
         *  Refers to the coverage with insert segments (not read alignments)
         *  of all alignments (not only those in the clique). */
        size_t coverage;
        /** Start of intersection of insert segments. */
        unsigned int start;
        /** End of intersection of insert segments. */
        unsigned int end;
        /** Length intersection of insert segments. */
        int length;
        /** Difference of weighted mean of insert segments lengths and
         *  population mean; i.e. negative values indicate insertions
         *  and positive values indicate deletions. */
        double diff;
        /** Locally corrected p-value. */
        double pvalue_corr;
        /** False discovery rate (FDR) achieved when the list of deemed-significant
         *  cliques would end after this one. */
        double fdr_level;
        /** Is this clique significant after the concluding FDR control. */
        bool is_significant;
        /** Bitvector: samples that gave best p-value. */
        int best_sample_combination;
        /** Information on every sample (if sample information is present) */
        std::vector<sample_wise_stats_t> sample_wise_stats;
        /** Information on every subset of samples (if sample information is present) */
        std::vector<sample_subset_wise_stats_t> sample_subset_wise_stats;
        /** If needed, reference to vector of all alignments constitution this clique. */
        std::vector<alignment_id_t>* reads;

        int window_start1;
        int window_end1;
        std::string consensus_string1;
        std::string coverage_string1;
        std::string phred_string1;

        int window_start2;
        int window_end2;
        std::string consensus_string2;
        std::string coverage_string2;
        std::string phred_string2;

        std::vector<std::string>* readnames;
        int clique_number;

        int maximum_coverage1;
        int maximum_coverage2;
        int min_coverage_user;

        int clique_size_weighted;

        clique_stats_t() : variation(), total_weight(0.0), clique_size(0), coverage(0), start(0), end(0), length(0), diff(0), pvalue_corr(0.0), fdr_level(-1.0), is_significant(false), best_sample_combination(-1), reads(0), window_start1(-1), window_end1(-1), window_start2(-1), window_end2(-1), clique_number(0), clique_size_weighted(0) {
        }
    } clique_stats_t;

        //     int window_start[2];
        // int window_end[2];
        // std::string consensus_string[2];
        // std::string coverage_string[2];
        // std::string phred_string[2];

        // std::vector<std::string>* readnames;
        // int clique_number;

        // int maximum_coverage[2];
        // int min_coverage_user;

        // int clique_size_weighted;

        // clique_stats_t() : variation(), total_weight(0.0), clique_size(0), coverage(0), start(0), end(0), length(0), diff(0), pvalue_corr(0.0), fdr_level(-1.0), is_significant(false), best_sample_combination(-1), reads(0), window_start[0](-1), window_end[0](-1), window_start[1](-1), window_end[1](-1), clique_number(0), clique_size_weighted(0) {
        // }


    typedef struct {
        std::vector<int> alignment;
    } insertion;

    /** Comparator type comparing cliques according to their corrected p-value. */
    typedef struct {

        bool operator()(const clique_stats_t& stats0, const clique_stats_t& stats1) {
            return stats0.pvalue_corr < stats1.pvalue_corr;
        }
    } clique_stats_comp_t;

    /** Comparator type comparing read name indices (w.r.t. "readname_to_index")
     *  according to the underlying read names. */
    typedef struct readname_comparator_t {
        CliqueWriter& cw;

        bool operator()(const size_t& index1, const size_t& index2) const {
            const std::string& name1 = cw.readname_to_index.right.at(index1);
            const std::string& name2 = cw.readname_to_index.right.at(index2);
            return name1.compare(name2) < 0;
        }

        readname_comparator_t(CliqueWriter& cw) : cw(cw) {
        }
    } readname_comparator_t;

    typedef struct fastq_entry {
        std::string name;
        std::set<std::string>* read_names;
        std::string seq_1;
        int pos_1;
        std::vector<std::string> phreds_1;
        std::string seq_2;
        int pos_2;
        std::vector<std::string> phreds_2;
    } fastq_entry;

    /** Type for bi-directional map from read names to read name indices. */
    typedef boost::bimap<std::string, size_t> readname_to_index_bimap_t;

    std::ostream& os;
	/** Output stream where indel calls are written to. If null, then no output is made. */
    /** Total number of processed cliques. */
	std::ostream* indel_os;
    long long total_count;
    /** Total number of all insertion-like cliques (i.e. includes non-significant). */
    size_t total_insertion_cliques;
    /** Total number of all deletion-like cliques (i.e. includes non-significant). */
    size_t total_deletion_cliques;
    /** Number of significant insertions (after FDR control). */
    long long significant_ins_count;
    /** Number of significant deletions (after FDR control). */
    long long significant_del_count;
    /** This variation caller is used to perform significance tests on cliques. */
    VariationCaller* variation_caller;
    /** If true, all cliques are output, not only significant ones. */
    bool output_all;
    /** Threshold for the false discovery rate (FDR). */
    double fdr_threshold;
    /** If true, FDR control has already been done and no more cliques 
     *  can be processed. */
    bool finished;
    /** List of cliques that are possibly significant. This list is only populated
     *  if output_all is false. */
    std::vector<clique_stats_t> clique_list;
    /** Stream to write a list of reads associated with significant cliques to.
     *  May be 0. */
    std::ostream* read_list_os;
    /** Bi-directional map that maps read names to indices and back. Only populated
     *  when necessary (e.g. when read list is to be written at the end). Only contains
     *  read that map to at least one (possibly) significant clique. */
    readname_to_index_bimap_t readname_to_index;
    /** If true, additional statistics are output for each written clique. */
    bool verbose;
    /** If true, then dataset is assumed to consist of multiple samples and each clique is evaluated accordingly. */
    bool multisample;
    const ReadGroups* read_groups;
    int min_coverage;
    int clique_count;
    int paired_count;
    int single_count;
    int single_skipped_count;
    bool FRAMESHIFT_MERGE;
    std::map<std::string,fastq_entry> fastq_map;
    std::string suffix;

    void writeReadlist();
    void callVariation(const std::vector<const AlignmentRecord*>& pairs, size_t coverage, clique_stats_t* stats);
    int shortenBase(char base);
    char expandBase(int base);
    bool overlapSize(clique_stats_t* stats) const;
    std::string equalStrings(std::string s1, std::string s2) const;
    void printout(int pos_1);
public:
    CliqueWriter(std::ostream& os, VariationCaller* variation_caller, std::ostream* indel_os, const ReadGroups* read_groups, bool multisample, bool output_all, double fdr_threshold, bool verbose, int min_coverage, bool frameshift_merge, std::string suffix);
    virtual ~CliqueWriter();
    virtual void enableReadListOutput(std::ostream& os);
    virtual void add(std::auto_ptr<Clique> clique);
    virtual void finish();

    virtual long long getTotalCount() {
        return total_count;
    }

    virtual int getPairedCount() {
        return paired_count;
    }

    virtual int getSingleCount() {
        return single_count;
    }

    virtual int getSingleSkippedCount() {
        return single_skipped_count;
    }

    virtual std::map<std::string,fastq_entry> getFastqMap() {
        return fastq_map;
    }

    friend std::ostream& operator<<(std::ostream& os, const clique_stats_t& stats);
};

#endif /* CLIQUEWRITER_H_ */
