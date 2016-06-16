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

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <limits>
#include <cassert>
#include <ctime>
#include <algorithm>
#include <deque>

#include "docopt/docopt.h"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <api/BamReader.h>
#include <api/BamAlignment.h>

#include "AlignmentRecord.h"
#include "QuasispeciesEdgeCalculator.h"
#include "NewEdgeCalculator.h"
#include "CliqueFinder.h"
#include "CLEVER.h"
#include "BronKerbosch.h"
#include "CliqueCollector.h"
#include "AnyDistributionEdgeCalculator.h"
#include "GaussianEdgeCalculator.h"
#include "LogWriter.h"

using namespace std;
using namespace boost;

static const char USAGE[] =
R"(haploclique predicts haplotypes from NGS reads.

Usage:
  haploclique bronkerbosch [options] [--] <bamfile> [<output>]
  haploclique [options] [--] <bamfile> [<output>]

  clever        use the original clever clique finder
  bronkerbosch  use the Bron-Kerbosch based clique finder

Options:
  -q NUM --edge_quasi_cutoff_cliques=NUM  edge calculator option
                                              [default: 0.99]
  -k NUM --edge_quasi_cutoff_mixed=NUM    edge calculator option
                                              [default: 0.97]
  -g NUM --edge_quasi_cutoff_single=NUM   edge calculator option
                                              [default: 0.95]
  -Q NUM --random_overlap_quality=NUM     edge calculator option
                                              [default: 0.9]
  -m --frame_shift_merge                      Reads will be clustered with 
                                              single nucleotide insertions or
                                              deletions. Use for PacBio data.
  -o NUM --min_overlap_cliques=NUM        edge calculator option
                                              [default: 0.9]
  -j NUM --min_overlap_single=NUM         edge calculator option
                                              [default: 0.6]
  -A FILE --allel_frequencies=FILE
  -I FILE --call_indels=FILE              variant calling is not supported
                                              yet.
  -M FILE --mean_and_sd_filename=FILE     Required for option -I
  -p NUM --indel_edge_sig_level=NUM       [default: 0.2]
  -i NUM --iterations=NUM                 Number of iterations.
                                              haploclique will stop if the 
                                              superreads converge.
                                              [default: -1]
  -f NUM --filter=NUM                     Filter out reads with low
                                              frequency at the end.
                                              [default: 0.0]
  -n --no_singletons                          Filter out single read cliques
                                              after first iteration.
  -s NUM --significance=NUM               Filter out reads witch are below
                                              <num> standard deviations.
                                              [default: 3.0]
  -L FILE, --log=FILE                      Write log to <file>.
  -d NUM --doc_haplotypes=NUM              Use for Simulation Study with known
                                           haplotypes to document which reads
                                           contributed to which final cliques (3 or 5).
  -p0 --noProb0                            ignore the tail probabilites during edge
                                           calculation.
  -gff --gff                               Option to create GFF File from output.
  -bam --bam                               Option to create BAM File from output
)";

void usage() {
    cerr << USAGE;
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

deque<AlignmentRecord*>* readBamFile(string filename, vector<string>& readNames, unsigned int& maxPosition) {
    typedef std::unordered_map<std::string, AlignmentRecord*> name_map_t;
    name_map_t names_to_reads;
    deque<AlignmentRecord*>* reads = new deque<AlignmentRecord*>;
    BamTools::BamReader bamreader;
    if (not bamreader.Open(filename)) {
        cerr << bamreader.GetFilename() << endl;
        throw std::runtime_error("Couldn't open Bamfile");
    }
    BamTools::BamAlignment alignment;

    //int readcounter = 0;
    //int singleendcounter = 0;
    //int pairedendcounter = 0;
    //int overlap = 0;

    while (bamreader.GetNextAlignment(alignment)) {
        /*std::size_t found = alignment.QueryBases.find('N');
        if (found!=std::string::npos){
            cout << alignment.Name << endl;
            cout << alignment.QueryBases << endl;
            cout << alignment.Qualities << endl;
            cout << endl;
        }*/
        bool valid = false;
        for (auto& i : alignment.CigarData){
            if (i.Type == 'M' && i.Length > 0) valid = true;
        }
        if(alignment.CigarData.size() > 0 && valid){
            if(names_to_reads.count(alignment.Name) > 0) {
                //cout << alignment.Name << endl;
                names_to_reads[alignment.Name]->pairWith(alignment);
                /*if (names_to_reads[alignment.Name]->isSingleEnd()){
                    ++overlap;
                    ++readcounter;
                } else {
                    pairedendcounter++;
                    ++readcounter;
                }*/
                reads->push_back(names_to_reads[alignment.Name]);
                names_to_reads.erase(alignment.Name);
            } else {
                names_to_reads[alignment.Name] = new AlignmentRecord(alignment, readNames.size(), &readNames);
                readNames.push_back(alignment.Name);
                //++readcounter;
            }
        }
    }

    // Push all single-end reads remaining in names_to_reads into the reads vector. Unmapped reads are filtered out in advance.
    for (const auto& i : names_to_reads) {
        reads->push_back(i.second);
        //readcounter++;
    }

    bamreader.Close();

    auto comp = [](AlignmentRecord* r1, AlignmentRecord* r2) { return r1->getIntervalStart() < r2->getIntervalStart(); };

    sort(reads->begin(), reads->end(), comp);
    cout << "Read BamFile: done" << endl;

    //DEBUGGING
    //check which alignments have no cigar string -> no covered positions -> no edge possible
    /*while(not reads->empty()){
        assert(reads->front() != nullptr);
        unique_ptr<AlignmentRecord> al_ptr(reads->front());
        reads->pop_front();
        if (al_ptr->getCigar1().size() == 0){
            cout << al_ptr->getName() << endl;
        }
    }
    //cout << "Number of Reads: " << readcounter << endl;
    //cout << "Number of Pairs: " << pairedendcounter << endl;
    //cout << "Number of Overlapping Pairs: " << overlap << endl;
    //check number of single ends and paired ends after merging
    int counter = 0;
    while(not reads->empty()){
            assert(reads->front() != nullptr);
            unique_ptr<AlignmentRecord> al_ptr(reads->front());
            reads->pop_front();
            if (al_ptr->isSingleEnd()){
                singleendcounter++;
            } else if (al_ptr->isPairedEnd()){
                pairedendcounter++;
            }
            counter++;
    }
    cout << singleendcounter << " " << pairedendcounter << " " <<  counter << endl;
    */
    maxPosition = (*std::max_element(std::begin(*reads), std::end(*reads),
                              [](const AlignmentRecord* lhs,const AlignmentRecord* rhs){
        return lhs->getIntervalEnd() < rhs->getIntervalEnd();
    }))->getIntervalEnd();


    return reads;
}

int main(int argc, char* argv[]) {
 
    //cout << "Test\n" << endl;

    map<std::string, docopt::value> args
        = docopt::docopt(USAGE,
                         { argv + 1, argv + argc },
                         false);

    for(auto elem : args)
    {
   	cout << elem.first << " " << elem.second << "\n";
    }
    // PARAMETERS
    string bamfile = args["<bamfile>"].asString();
    string outfile = "quasispecies.fasta";
    if (args["<output>"]) outfile = args["<output>"].asString();
    double edge_quasi_cutoff_cliques = stod(args["--edge_quasi_cutoff_cliques"].asString());
    double edge_quasi_cutoff_single = stod(args["--edge_quasi_cutoff_single"].asString());
    double edge_quasi_cutoff_mixed = stod(args["--edge_quasi_cutoff_mixed"].asString());
    double Q = stod(args["--random_overlap_quality"].asString());
    double overlap_cliques = stod(args["--min_overlap_cliques"].asString());
    double overlap_single = stod(args["--min_overlap_single"].asString());
    bool frameshift_merge = args["--frame_shift_merge"].asBool();
    string allel_frequencies_path = "";
    if (args["--allel_frequencies"]) allel_frequencies_path = args["--allel_frequencies"].asString();
    string mean_and_sd_filename = "";
    if (args["--mean_and_sd_filename"]) mean_and_sd_filename = args["--mean_and_sd_filename"].asString();
    string indel_edge_cutoff;
    double indel_edge_sig_level = stod(args["--indel_edge_sig_level"].asString());
    string indel_output_file = "";
    if (args["--call_indels"]) indel_output_file = args["--call_indels"].asString();
    int iterations = stoi(args["--iterations"].asString());
    double filter = stod(args["--filter"].asString());
    double significance = stod(args["--significance"].asString());
    bool filter_singletons = args["--no_singletons"].asBool();
    string logfile = "";
    if (args["--log"]) logfile = args["--log"].asString();
    int doc_haplotypes = 0;
    if (args["--doc_haplotypes"]) doc_haplotypes = stoi(args["--doc_haplotypes"].asString());
    bool noProb0 = args["--noProb0"].asBool();
    //-ot CHAR --output_type=CHAR
    bool bam = args["--bam"].asBool();
    bool gff = args["--gff"].asBool();


    // END PARAMETERS

    bool call_indels = indel_output_file.size() > 0;
    if (call_indels && (mean_and_sd_filename.size() == 0)) {
        cerr << "Error: when using option -I, option -M must also be given." << endl;
        return 1;
    }

    //read allel frequency distributions
    std::unordered_map<int, double> simpson_map;
    unsigned int maxPosition2 = 0;
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
                simpson_map[atoi(words[0].c_str())] = std::log10(pow(atof(words[1].c_str()),2)+pow(atof(words[2].c_str()),2)+pow(atof(words[3].c_str()),2)+pow(atof(words[4].c_str()),2));
                maxPosition2=atoi(words[0].c_str());
                //cerr << simpson_map[atoi(words[0].c_str())] << endl;
            }
        }
        ia.close();
    }
    //cout << "PARSE PRIOR: done" << endl;


    clock_t clock_start = clock();
    vector<string> originalReadNames;
    unsigned int maxPosition1;
    deque<AlignmentRecord*>* reads = readBamFile(bamfile, originalReadNames,maxPosition1);
    EdgeCalculator* edge_calculator = nullptr;
    EdgeCalculator* indel_edge_calculator = nullptr;
    unique_ptr<vector<mean_and_stddev_t> > readgroup_params(nullptr);
    maxPosition1 = (maxPosition1>maxPosition2) ? maxPosition1 : maxPosition2;
    edge_calculator = new NewEdgeCalculator(Q, edge_quasi_cutoff_cliques, overlap_cliques, frameshift_merge, simpson_map, edge_quasi_cutoff_single, overlap_single, edge_quasi_cutoff_mixed, maxPosition1, noProb0);
    if (call_indels) {
        double insert_mean = -1.0;
        double insert_stddev = -1.0;
        if (!read_mean_and_sd(mean_and_sd_filename, &insert_mean, &insert_stddev)) {
            cerr << "Error reading \"" << mean_and_sd_filename << "\"." << endl;
            return 1;
        }
        cerr << "Null distribution: mean " << insert_mean << ", sd " <<  insert_stddev << endl;
        indel_edge_calculator = new GaussianEdgeCalculator(indel_edge_sig_level,insert_mean,insert_stddev);
    }
    std::ofstream* indel_os = nullptr;
    if (call_indels) {
        indel_os = new ofstream(indel_output_file.c_str());
    }
    LogWriter* lw = nullptr;
    if (logfile != "") lw = new LogWriter(logfile);

    CliqueCollector collector(lw);
    CliqueFinder* clique_finder;

    if (args["bronkerbosch"].asBool()) {
        clique_finder = new BronKerbosch(*edge_calculator, collector, lw);
    } else {
        clique_finder = new CLEVER(*edge_calculator, collector, lw);
    }
    if (indel_edge_calculator != 0) {
        clique_finder->setSecondEdgeCalculator(indel_edge_calculator);
    }

    ofstream* reads_ofstream = 0;



    // Main loop
    int ct = 0;
    double stdev = 1.0;
    auto filter_fn = [&](unique_ptr<AlignmentRecord>& read, int size) {
        //if (ct == 8){
        //    cout << read->getName() << " " <<read->getProbability() << " " << 1.0 /(size) << " " << significance*stdev << endl;
        //    bool t =  read->getProbability() < 1.0 /size - significance*stdev;
        //    cout << t << endl;
        //}
        return (ct == 1 and filter_singletons and read->getReadCount() <= 1) or (ct > 1 and significance != 0.0 and read->getProbability() < 1.0 / size - significance*stdev);
    };
    int edgecounter = 0;
    cout << "start: " << originalReadNames.size();
    while (ct != iterations) {
        clique_finder->initialize();
        //cout << "Clique_finder initialized " << ct << endl;
        //if(ct >= 5 && reads->size()<100 ){
        //if(ct== 5){
        //    edge_calculator->setOverlapCliques(0.1);
        //}
        //if(ct==8){
        //    int k = 0;
        //}
        if (lw != nullptr) lw->initialize();
        int size = reads->size();
        while(not reads->empty()) {
            assert(reads->front() != nullptr);

            unique_ptr<AlignmentRecord> al_ptr(reads->front());
            //if(al_ptr->getName() == "Clique_1037"){
            //    int k = 0;
            //}
            reads->pop_front();
            if (filter_fn(al_ptr,size)) continue;

            clique_finder->addAlignment(al_ptr,edgecounter);
            //cout << "addAlignment " << ct << endl;
        }
        cout << "\tedges: " << edgecounter << endl;

        delete reads;
        //cout << "reads deleted " << ct << endl;
        clique_finder->finish();
        //cout << "clique_finder finished " << ct << endl;
        reads = collector.finish();
        //cout << "collector finish " << ct << endl;
        if (lw != nullptr) lw->finish();

        stdev = setProbabilities(*reads);
        //if(ct == 5) break;
        if (clique_finder->hasConverged()) break;
        cout << ct++ << ": " << reads->size();
        edgecounter = 0;
    }

    // Filter superreads according to read probability
    if (filter > 0.0) {
        auto filter_fn = [&](AlignmentRecord* al) { return al->getProbability() < filter;};
        auto new_end_it = std::remove_if(reads->begin(), reads->end(), filter_fn);
        for (auto it = new_end_it; it != reads->end(); it++) {
            delete *it;
            *it = nullptr;
        }
        reads->erase(new_end_it, reads->end());
    }
    ofstream os(outfile + ".fasta", std::ofstream::out);
    setProbabilities(*reads);
    printReads(os, *reads, doc_haplotypes);
    if (gff){
        ofstream os1(outfile + ".gff",std::ofstream::out);
        printGFF(os1,*reads);
    }
    if (bam){
        ofstream os2(outfile + ".bam",std::ofstream::out);
        printBAM(os2,*reads);
    }

    cout << "final: " << reads->size() << endl;

    for (auto&& r : *reads) {
    //    if(r->getName() == "Clique_3370"){
    //        for (auto& i : r->getReadNames()){
    //            cout << i << endl;
    //       }
    //    }
        delete r;
    }

    delete reads;

    if (indel_os != nullptr) {
        indel_os->close();
        delete indel_os;
    }
    if (edge_calculator != nullptr) delete edge_calculator;
    if (lw != nullptr) delete lw;
    if (clique_finder != nullptr) delete clique_finder;
    if (reads_ofstream != nullptr) {
        delete reads_ofstream;
    }
    cout.precision(3);
    cout << std::fixed;
    double cpu_time = (double) (clock() - clock_start) / CLOCKS_PER_SEC;
    cout << "time: " <<  cpu_time << endl;
    return 0;
}
