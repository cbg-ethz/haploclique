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

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <cassert>
#include <ctime>

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <bamtools/api/BamReader.h>

#include "AlignmentRecord.h"
#include "QuasispeciesEdgeCalculator.h"
#include "CliqueFinder.h"
#include "CliqueWriter.h"
#include "CoverageMonitor.h"
#include "HistogramBasedDistribution.h"
#include "AnyDistributionEdgeCalculator.h"
#include "ReadSetSignificanceTester.h"
#include "ReadSetZTester.h"
#include "ReadSetGenericTester.h"
#include "ReadGroups.h"
#include "ReadGroupAwareEdgeCalculator.h"
#include "ReadSetGroupWiseZTester.h"
#include "GaussianEdgeCalculator.h"

using namespace std;
using namespace boost;
namespace po = boost::program_options;

void usage(const char* name, const po::options_description& options_desc) {
    cerr << "Usage: " << name << " [options]" << endl;
    cerr << endl;
    cerr << "<distribution-file> file with assumed internal segment length distribution." << endl;
    cerr << "                    in default mode, this is a file containing one line with" << endl;
    cerr << "                    mean and standard deviation of the normal distribution" << endl;
    cerr << "                    to be used. Such a file can be generated using option -m" << endl;
    cerr << "                    of insert-length-histogram. Note that the \"internal segment\"" << endl;
    cerr << "                    does NOT include the read (ends), i.e. a fragment disjointly" << endl;
    cerr << "                    consists of two reads (read ends) and an internal segment." << endl;
    cerr << endl;
    cerr << "Reads alignments / alignment pairs from stdin and computes all cliques." << endl;
    cerr << "Format for single-end reads:" << endl;
    cerr << "<read-name> <read-nr> <read-group> <phred-sum> <chromosome> <start> <end> <strand> <cigar> <seq> <qualities> <aln-prob>" << endl;
    cerr << "Format for paired-end reads:" << endl;
    cerr << "<read-name> <pair-nr> <read-group> <phred-sum1> <chromosome1> <start1> <end1> <strand1> <cigar1> <seq1> <qualities1> <phred-sum2> <chromosome2> <start2> <end2> <strand2> <cigar2> <seq2> <qualities2> <aln-pair-prob> <aln-pair-prob-inslength>" << endl;
    cerr << endl;
    cerr << "NOTE: Lines are assumed to be ordered by field 6 (start1)." << endl;
    cerr << endl;
    cerr << options_desc << endl;
    exit(1);
}

bool read_mean_and_sd(const string& filename, double* mean, double* sd) {
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer_t;
    boost::char_separator<char> whitespace_separator(" \t");
    ifstream in(filename.c_str());
    string line;
    if (in.fail() || (!getline(in,line))) {
        return false;
    }
    tokenizer_t tokenizer(line,whitespace_separator);
    vector<string> tokens(tokenizer.begin(), tokenizer.end());
    if (tokens.size() != 2) {
        return false;
    }
    try {
        *mean = boost::lexical_cast<double>(tokens[0]);
        *sd = boost::lexical_cast<double>(tokens[1]);
    } catch(boost::bad_lexical_cast &){
        return false;
    }
    return true;
}

int main(int argc, char* argv[]) {
    // PARAMETERS
    double min_aln_weight;
    double max_insert_length;
    int max_coverage;
    string edge_filename;
    double fdr;
    string reads_output_filename;
    string coverage_output_filename;
    bool verbose = false;
    double edge_quasi_cutoff_cliques;
    double edge_quasi_cutoff_single;
    double Q;
    double overlap_cliques;
    double overlap_single;
    bool frameshift_merge = false;
    int super_read_min_coverage;
    string allel_frequencies_path;
    string mean_and_sd_filename = "";
    string indel_edge_cutoff;
    double indel_edge_sig_level;
    string indel_output_file = "";
    int time_limit;
    bool no_sort;
    string suffix;

    po::options_description options_desc("Allowed options");
    options_desc.add_options()
    ("verbose,v", po::value<bool>(&verbose)->zero_tokens(), "Be verbose: output additional statistics for each variation.")
    ("min_aln_weight,w", po::value<double>(&min_aln_weight)->default_value(0.0016), "Minimum weight of alignment pairs to be considered.")
    ("max_insert_length,l", po::value<double>(&max_insert_length)->default_value(50000), "Maximum insert length of alignments to be considered (0=unlimited).")
    ("max_coverage,c", po::value<int>(&max_coverage)->default_value(200), "Maximum allowed coverage. If exceeded, violating reads are ignored. The number of such ignored reads is printed to stderr (0=unlimited).")
    ("write_edges,e", po::value<string>(&edge_filename)->default_value(""), "Write edges to file of given name.")
    ("fdr,f", po::value<double>(&fdr)->default_value(0.1), "False discovery rate (FDR).")
    ("all,a", "Output all cliques instead of only the significant ones. Cliques are not sorted and last column (FDR) is not computed.")
    ("output_reads,r", po::value<string>(&reads_output_filename)->default_value(""), "Output reads belonging to at least one significant clique to the given filename (along with their assignment to significant cliques.")
    ("output_coverage,C", po::value<string>(&coverage_output_filename)->default_value(""), "Output the coverage with considered insert segments along the chromosome (one line per position) to the given filename.")
    ("edge_quasi_cutoff_cliques,q", po::value<double>(&edge_quasi_cutoff_cliques)->default_value(0.99), "End compatibility probability cutoff between error-corrected reads for quasispecies reconstruction.")
    ("edge_quasi_cutoff_single,g", po::value<double>(&edge_quasi_cutoff_single)->default_value(0.95), "End compatibility probability cutoff between raw<->raw and raw<->error-corrected reads for quasispecies reconstruction.")
    ("random_overlap_probability,Q", po::value<double>(&Q)->default_value(0.9), "Probability that two random reads are equal at the same position.")
    ("frame_shift_merge,m", po::value<bool>(&frameshift_merge)->zero_tokens(), "Reads will be clustered if one has single nucleotide deletions and insertions. Use for PacBio data sets.")
    ("min_overlap_cliques,o", po::value<double>(&overlap_cliques)->default_value(0.9), "Minimum relative overlap between error-corrected reads.")
    ("min_overlap_single,j", po::value<double>(&overlap_single)->default_value(0.6), "Minimum relative overlap between raw<->raw and raw<->error-corrected reads.")
    ("super_read_min_coverage,s", po::value<int>(&super_read_min_coverage)->default_value(2), "Minimum coverage for super-read assembly.")
    ("allel_frequencies,A", po::value<string>(&allel_frequencies_path)->default_value(""), "Minimum coverage for super-read assembly.")
    ("call_indels,I", po::value<string>(&indel_output_file)->default_value(""), "Call indels from cliques. In this mode, the \"classical CLEVER\" edge criterion is used in addition to the new one. Filename to write indels to must be given as parameter.")
    ("mean_and_sd_filename,M", po::value<string>(&mean_and_sd_filename)->default_value(""), "Name of file with mean and standard deviation of insert size distribution (only required if option -I is used).")
    ("indel_edge_sig_level,p", po::value<double>(&indel_edge_sig_level)->default_value(0.2), "Significance level for \"indel\" edges criterion, see option -I (the lower the level, the more edges will be present).")
    ("time_limit,t", po::value<int>(&time_limit)->default_value(10), "Time limit for computation. If exceeded, non processed reads will be written to skipped.")
    ("no_sort,N", po::value<bool>(&no_sort)->zero_tokens(), "Do not sort new clique w.r.t. their bitsets.")
    ("suffix,S", po::value<string>(&suffix)->default_value(""), "Suffix for clique names. Used for parallelization.")
    ;

    if (isatty(fileno(stdin))) {
        usage(argv[0], options_desc);
    }

    po::variables_map options;
    try {
        po::store(po::parse_command_line(argc, argv, options_desc), options);
        po::notify(options);
    } catch (std::exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }

    bool output_all = options.count("all") > 0;
    if (output_all && (reads_output_filename.size() > 0)) {
        cerr << "Error: options -a and -r are mutually exclusive." << endl;
        return 1;
    }
    bool call_indels = indel_output_file.size() > 0;
    if (call_indels && (mean_and_sd_filename.size() == 0)) {
        cerr << "Error: when using option -I, option -M must also be given." << endl;
        return 1;
    }



    std::map<string,string> clique_to_reads_map;
    ifstream tsv_stream("data_clique_to_reads.tsv");
    string tsv_stream_line;
    while (getline(tsv_stream, tsv_stream_line)) {
        std::vector<std::string> words;
        trim_right(tsv_stream_line);
        boost::split(words, tsv_stream_line, boost::is_any_of("\t"), boost::token_compress_on);
        clique_to_reads_map[words[0]] = words[1];
    }
    remove("data_clique_to_reads.tsv");
    
    //read allel frequency distributions
    std::map<int, double> simpson_map;
    //cerr << "PARSE PRIOR";
    cerr.flush();
    if (allel_frequencies_path.size() > 0) {
        ifstream ia(allel_frequencies_path.c_str());
        string ia_line;
        while (getline(ia, ia_line)) {
            std::vector<std::string> words;
            trim_right(ia_line);
            boost::split(words, ia_line, boost::is_any_of("\t"), boost::token_compress_on);

            std::vector<std::string> insertion_words;
            boost::split(insertion_words, words[0], boost::is_any_of("\\."), boost::token_compress_on);
            if (insertion_words.size() > 1) {
            } else {
                simpson_map[atoi(words[0].c_str())] = pow(atof(words[1].c_str()),2)+pow(atof(words[2].c_str()),2)+pow(atof(words[3].c_str()),2)+pow(atof(words[4].c_str()),2)+pow(atof(words[5].c_str()),2);
                //cerr << simpson_map[atoi(words[0].c_str())] << endl;
            }
        }
        ia.close();
    }
    //cerr << "PARSE PRIOR: done" << endl;
    

    clock_t clock_start = clock();
    EdgeCalculator* edge_calculator = 0;
    EdgeCalculator* indel_edge_calculator = 0;
    ReadSetSignificanceTester* significance_tester = 0;
    VariationCaller* variation_caller = 0;
    ReadGroups* read_groups = 0;
    auto_ptr<vector<mean_and_stddev_t> > readgroup_params(0);
    edge_calculator = new QuasispeciesEdgeCalculator(Q, edge_quasi_cutoff_cliques, overlap_cliques, frameshift_merge, simpson_map, edge_quasi_cutoff_single, overlap_single);
    if (call_indels) {
        double insert_mean = -1.0;
        double insert_stddev = -1.0;
        if (!read_mean_and_sd(mean_and_sd_filename, &insert_mean, &insert_stddev)) {
            cerr << "Error reading \"" << mean_and_sd_filename << "\"." << endl;
            return 1;
        }
        cerr << "Null distribution: mean " << insert_mean << ", sd " <<  insert_stddev << endl;
        indel_edge_calculator = new GaussianEdgeCalculator(indel_edge_sig_level,insert_mean,insert_stddev);
        significance_tester = new ReadSetZTester(insert_mean, insert_stddev);
        variation_caller = new VariationCaller(insert_mean, *significance_tester);
    }
    std::ofstream* indel_os = 0;
    if (call_indels) {
        indel_os = new ofstream(indel_output_file.c_str());
    }
    CliqueWriter clique_writer(cout, variation_caller, indel_os, read_groups, false, output_all, fdr, verbose, super_read_min_coverage, frameshift_merge, suffix);
    CliqueFinder clique_finder(*edge_calculator, clique_writer, read_groups, no_sort);
    if (indel_edge_calculator != 0) {
        clique_finder.setSecondEdgeCalculator(indel_edge_calculator);
    }
    EdgeWriter* edge_writer = 0;
    ofstream* edge_ofstream = 0;
    if (edge_filename.size() > 0) {
        edge_ofstream = new ofstream(edge_filename.c_str());
        edge_writer = new EdgeWriter(*edge_ofstream);
        clique_finder.setEdgeWriter(*edge_writer);
    }
    ofstream* reads_ofstream = 0;
    if (reads_output_filename.size() > 0) {
        reads_ofstream = new ofstream(reads_output_filename.c_str());
        clique_writer.enableReadListOutput(*reads_ofstream);
    }

    size_t last_pos = 0;
    int n = 0;
    string line;
    size_t skipped_by_weight = 0;
    size_t skipped_by_length = 0;
    size_t skipped_by_coverage = 0;
    size_t valid_alignments = 0;
    size_t total_alignments = 0;
    cerr << "STATUS";
    cerr.flush();
    while (getline(cin, line)) {
        n += 1;
        total_alignments += 1;
        try {
            AlignmentRecord ap(line, clique_to_reads_map, read_groups);
            if (ap.getIntervalStart() < last_pos) {
                cerr << "Error: Input is not ordered by position (field 6)! Offending line: " << n << endl;
                return 1;
            }
            if (ap.isPairedEnd()) {
                if (ap.getChrom1().compare(ap.getChrom2()) != 0) continue;
                if (ap.getStrand1().compare(ap.getStrand2()) == 0) continue;
            }
            valid_alignments += 1;
            last_pos = ap.getIntervalStart();
            auto_ptr<AlignmentRecord> alignment_autoptr(new AlignmentRecord(ap));
            if (ap.isPairedEnd() && (max_insert_length > 0)) {
                if (alignment_autoptr->getInsertLength() > max_insert_length) {
                //skipped_by_length += 1;
                //continue;
                }
            }
            /*if (alignment_autoptr->getWeight() < min_aln_weight) {
            // cout << "Skipping alignment (weight): "  << alignment_autoptr->getName() << " weight: " << alignment_autoptr->getWeight() << endl;
                skipped << line;
                skipped_by_weight += 1;
                continue;
            }*/
            bool time = ((double) (clock() - clock_start) / CLOCKS_PER_SEC / 60) > time_limit;
            if (max_coverage > 0 || time) {
                if (clique_finder.getCoverageMonitor().probeAlignment(*alignment_autoptr) > (size_t) max_coverage || time) {
                // cout << "Skipping alignment (coverage): "  << alignment_autoptr->getName()  << endl;
                    skipped_by_coverage += 1;
                    continue;
                }
            }
            clique_finder.addAlignment(alignment_autoptr);
        } catch (std::runtime_error&) {
            cerr << "Error parsing input, offending line: " << n << endl;
            return 1;
        }
        cerr << "\rSTATUS: " << total_alignments;
        cerr.flush();
    }
    clique_finder.finish();
    clique_writer.finish();
    cerr << endl;

    if (indel_os != 0) {
        indel_os->close();
        delete indel_os;
    }
    if (edge_calculator != 0) delete edge_calculator;
    if (variation_caller != 0) delete variation_caller;
    if (significance_tester != 0) delete significance_tester;
    if (edge_writer != 0) {
        delete edge_writer;
        delete edge_ofstream;
    }
    if (reads_ofstream != 0) {
        delete reads_ofstream;
    }
    double cpu_time = (double) (clock() - clock_start) / CLOCKS_PER_SEC;
    cerr << "Cliques/Uniques/CPU time:\t" << clique_writer.getPairedCount() << "/" << clique_writer.getSingleCount() << "/" << round(cpu_time) << endl;
    return 0;
}
