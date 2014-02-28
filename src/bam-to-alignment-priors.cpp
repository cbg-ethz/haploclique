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
#include <limits>
#include <cassert>
#include <iomanip>
#include <ctime>

#include <boost/program_options.hpp>
#include <boost/unordered_map.hpp>
#include <boost/tokenizer.hpp>
#include <bamtools/api/BamReader.h>
#include <bamtools/api/BamWriter.h>

#include "HistogramBasedDistribution.h"
#include "SortedBamReader.h"
#include "Histogram.h"
#include "ASTagCalculator.h"
#include "ThreadPool.h"
#include "GroupWiseBamReader.h"

using namespace std;
namespace po = boost::program_options;

int phred_base = 33;
int concordant_min = -1;
int concordant_max = -1;
long long skipped_inside_clipping = 0;

typedef boost::unordered_map<std::string, HistogramBasedDistribution*> readgroupwise_distributions_t;

void usage(const char* name, const po::options_description& options_desc) {
	cerr << "Usage: " << name << " [options] <reference.fasta(.gz)> <input.bam>" << endl;
	cerr << endl;
	cerr << "Reads a sorted BAM file (also parsing the BWA-specific XA tags)" << endl;
	cerr << "and outputs a list of plausible alignments (one per line) of read pairs" << endl;
	cerr << "to stdout in the following format:" << endl;
	cerr << endl;
	cerr << "<read-name> <pair-nr> <read-group> <phred-sum1> <chromosome1> <start1> <end1> <strand1> <cigar1> <seq1> <qualities1> <phred-sum2> <chromosome2> <start2> <end2> <strand2> <cigar2> <seq2> <qualities2> <aln-pair-prob> <aln-pair-prob-inslength>" << endl;
	cerr << endl;
	cerr << "Where <aln-pair-prob> is the probability that this alignment pair is " << endl;
	cerr << "correct estimated based on the alignment quality while the probabilities " << endl;
	cerr << "given in <aln-pair-prob-inslength> are estimates based on alignment quality AND" << endl;
	cerr << "insert length." << endl;
	cerr << endl;
	cerr << "<startX> and <endX> coordinates are 1-based and inclusive (i.e. closed intervals)" << endl;
	cerr << "and give the region on the reference the respective read was aligned to." << endl;
	cerr << endl;
	cerr << "In case of single end alignments, the format is changed to:" << endl;
	cerr << "<read-name> <read-nr> <read-group> <phred-sum> <chromosome> <start> <end> <strand> <cigar> <seq> <qualities> <aln-prob>" << endl;
	cerr << endl;
	cerr << options_desc << endl;
	exit(1);
}

// TODO: some duplicate code shared with bam-to-alignment-priors and add-score-tags-to-bam

bool is_clipped(const BamTools::BamAlignment& alignment) {
	vector<BamTools::CigarOp>::const_iterator it = alignment.CigarData.begin();
	for (;it!=alignment.CigarData.end(); ++it) {
		switch (it->Type) {
		case 'S':
		case 'H':
		case 'P':
			return true;
		default:
			continue;
		}
	}
	return false;
}

bool is_unique(const BamTools::BamAlignment& alignment) {
	if (alignment.MapQuality == 0) return false;
	if (!alignment.IsPrimaryAlignment()) return false;
	uint32_t x0_tag = -1;
	uint32_t x1_tag = -1;
	if (alignment.GetTag("X0",x0_tag)) {
		if (x0_tag>1) return false;
	}
	if (alignment.GetTag("X1",x1_tag)) {
		if (x1_tag>0) return false;
	}
	return true;
}

/** Attaches probabilities corresponding with the 'AS' scores of the alignments. */
vector<double>* compute_posteriors(const vector<BamTools::BamAlignment*>& alignments) {
	double probsum = 0.0;
	vector<double>* probs = new vector<double>();
	typedef vector<BamTools::BamAlignment*>::const_iterator aln_iter_t;
	for (aln_iter_t it=alignments.begin(); it!=alignments.end(); ++it) {
		uint32_t as = 0;
		if ((*it)->GetTag("AS", as)) {
			double p = pow(10.0,as/-10.0);
			probs->push_back(p);
			probsum += p;
		} else {
			assert(false);
		}
	}
	if (probsum > 0.0) {
		for (size_t i=0; i<probs->size(); ++i) {
			probs->at(i) /= probsum;
		}
	} else {
		for (size_t i=0; i<probs->size(); ++i) {
			probs->at(i) = 1.0/alignments.size();
		}
	}
	return probs;
}

/** Returns probability "distribution" according to MAPQ fields. Returns 0 if MAPQ is unavailible. */
vector<double>* mapq_probabilities(const vector<BamTools::BamAlignment*>& alignments) {
	typedef vector<BamTools::BamAlignment*>::const_iterator aln_iter_t;
	for (aln_iter_t it=alignments.begin(); it!=alignments.end(); ++it) {
		if ((*it)->MapQuality == 255) return 0;
	}
	double probsum = 0.0;
	vector<double>* probs = new vector<double>();
	for (aln_iter_t it=alignments.begin(); it!=alignments.end(); ++it) {
		if ((*it)->MapQuality == 0) {
			probs->push_back(0.0);
		} else {
			double p = 1.0 - pow(10.0,(*it)->MapQuality/-10.0);
			probs->push_back(p);
			probsum += p;
		}
	}
	if (probsum > 1.0) {
		for (size_t i=0; i<probs->size(); ++i) {
			probs->at(i) /= probsum;
		}
	}
	return probs;
}

int phred_sum(const BamTools::BamAlignment& aln) {
	const string& phred = aln.Qualities;
	int result = 0;
	for (size_t i=0; i<phred.size(); ++i) {
		result += phred[i] - phred_base;
	}
	return result;
}

typedef struct aln_pair_t {
	const BamTools::BamAlignment* aln1;
	const BamTools::BamAlignment* aln2;
	double pair_prob;
	double pair_prob_ins_length;
	aln_pair_t(const BamTools::BamAlignment* aln1, const BamTools::BamAlignment* aln2, double pair_prob, double pair_prob_ins_length) : aln1(aln1), aln2(aln2), pair_prob(pair_prob), pair_prob_ins_length(pair_prob_ins_length) {}
} aln_pair_t;

/** Given an output file (a plain file object) and a list of alignments belonging
to the same read pair, outputs one line per alignment pair. */
void process_read_allpairs(const BamTools::RefVector& ref_vector, const vector<BamTools::BamAlignment*>& alignments1, const vector<BamTools::BamAlignment*>& alignments2, const HistogramBasedDistribution& ild, ostream& out_stream, bool use_mapq) {
	vector<double>* posteriors1;
	vector<double>* posteriors2;
	if (use_mapq) {
		posteriors1 = mapq_probabilities(alignments1);
		if (posteriors1 == 0) return;
		posteriors2 = mapq_probabilities(alignments2);
		if (posteriors2 == 0) {
			delete posteriors1;
			return;
		}
	} else {
		posteriors1 = compute_posteriors(alignments1);
		posteriors2 = compute_posteriors(alignments2);
	}
	assert(posteriors1->size() == alignments1.size());
	assert(posteriors2->size() == alignments2.size());
	// contains an entry for every pair that should be printed
	vector<aln_pair_t> all_pairs;
	double p_sum = 0.0;
	double p_ins_sum = 0.0;
	for (size_t i=0; i<alignments1.size(); ++i) {
		const BamTools::BamAlignment& read_aln1 = *alignments1[i];
		assert(read_aln1.CigarData.size() > 0);
		double posterior1 = posteriors1->at(i);
		assert(read_aln1.IsMapped());
		for (size_t j=0; j<alignments2.size(); ++j) {
			const BamTools::BamAlignment& read_aln2 = *alignments2[j];
			assert(read_aln2.CigarData.size() > 0);
			double posterior2 = posteriors2->at(j);
			assert(read_aln2.IsMapped());
			// skip if strandedness is the same
			if (read_aln1.IsReverseStrand() == read_aln2.IsReverseStrand()) continue;
			// skip if both reads map to different chromosomes
			if (read_aln1.RefID != read_aln2.RefID) continue;
			double p = posterior1*posterior2;
			double p_ins = p;
			const BamTools::BamAlignment& left = (read_aln1.Position<=read_aln2.Position)?read_aln1:read_aln2;
			const BamTools::BamAlignment& right = (read_aln1.Position<=read_aln2.Position)?read_aln2:read_aln1;
//			if (((left.CigarData[left.CigarData.size()-1].Type == 'S') && (left.CigarData[left.CigarData.size()-1].Length > 2)) ||
//				((right.CigarData[0].Type == 'S') && (right.CigarData[0].Length > 2))) {
//				skipped_inside_clipping += 1;
//				delete posteriors1;
//				delete posteriors2;
//				return;
//			}
			int insert_length = right.Position - left.GetEndPosition();
			p_ins *= ild.probability(insert_length);
			p_sum += p;
			p_ins_sum += p_ins;
			// if ((concordant_min!=-1) && (concordant_max!=-1)) {
			// 	if ((concordant_min<=insert_length) && (insert_length<=concordant_max)) continue;
			// }
			all_pairs.push_back(aln_pair_t(&left, &right, p, p_ins));
		}
	}
	string rg = "";
	if (!alignments1[0]->GetTag("RG", rg)) {
		rg = "-1";
	}
	if ((p_sum > 0.0) && (p_ins_sum > 0.0)) {
		int n = 0;
		for (size_t i=0; i<all_pairs.size(); ++i) {
			aln_pair_t& p = all_pairs[i];
			out_stream << p.aln1->Name << " " << (n++) << " " << rg << " "
			     << phred_sum(*(p.aln1)) << " " << ref_vector[p.aln1->RefID].RefName << " " << p.aln1->Position+1 << " " << p.aln1->GetEndPosition() << " " << (p.aln1->IsReverseStrand()?'-':'+') << " " << p.aln1->CigarData << " " << p.aln1->QueryBases << " " << p.aln1->Qualities << " "
			     << phred_sum(*(p.aln2)) << " " << ref_vector[p.aln2->RefID].RefName << " " << p.aln2->Position+1 << " " << p.aln2->GetEndPosition() << " " << (p.aln2->IsReverseStrand()?'-':'+') << " " << p.aln2->CigarData << " " << p.aln2->QueryBases << " " << p.aln2->Qualities << " "
		         << setprecision(16) << p.pair_prob/p_sum << " " << p.pair_prob_ins_length/p_ins_sum << endl;
		}
	}
	delete posteriors1;
	delete posteriors2;
}

void process_read_single_end(const BamTools::RefVector& ref_vector, const vector<BamTools::BamAlignment*>& alignments1, ostream& out_stream, bool use_mapq) {
	vector<double>* posteriors;
	if (use_mapq) {
		posteriors = mapq_probabilities(alignments1);
		if (posteriors == 0) return;
	} else {
		posteriors = compute_posteriors(alignments1);
	}
	assert(posteriors->size() == alignments1.size());
	// contains an entry for every pair that should be printed
	vector<aln_pair_t> all_pairs;
	double p_sum = 0.0;
	for (size_t i=0; i<alignments1.size(); ++i) {
		const BamTools::BamAlignment& aln = *alignments1[i];
		assert(aln.CigarData.size() > 0);
		assert(aln.IsMapped());
		p_sum += posteriors->at(i);
	}
	string rg = "";
	if (!alignments1[0]->GetTag("RG", rg)) {
		rg = "-1";
	}
	if (p_sum > 0.0) {
		for (size_t i=0; i<alignments1.size(); ++i) {
			const BamTools::BamAlignment& aln = *alignments1[i];
			out_stream << aln.Name << " " << i << " " << rg << " "
				<< phred_sum(aln) << " " << ref_vector[aln.RefID].RefName << " " << aln.Position+1 << " " << aln.GetEndPosition() << " " << (aln.IsReverseStrand()?'-':'+') << " " << aln.CigarData << " " << aln.QueryBases << " " << aln.Qualities << " "
				<< setprecision(16) << posteriors->at(i)/p_sum << endl;
		}
	}
	delete posteriors;
}

typedef struct work_package_t {
	vector<BamTools::BamAlignment*>* alignments1;
	vector<BamTools::BamAlignment*>* alignments2;
	const ASTagCalculator& as_tag_calculator;
	const BamTools::RefVector& bam_ref_data;
	const HistogramBasedDistribution& insert_size_dist;
	ASTagCalculator::stats_t as_stats;
	ostringstream output;
	bool use_mapq;

	work_package_t(vector<BamTools::BamAlignment*>* alignments1, vector<BamTools::BamAlignment*>* alignments2, const ASTagCalculator& as_tag_calculator, const BamTools::RefVector& bam_ref_data, const HistogramBasedDistribution& insert_size_dist, bool use_mapq) : alignments1(alignments1), alignments2(alignments2), as_tag_calculator(as_tag_calculator), bam_ref_data(bam_ref_data), insert_size_dist(insert_size_dist), as_stats(), use_mapq(use_mapq) {}

	~work_package_t() {
		free(alignments1);
		if (alignments2 != 0) free(alignments2);
	}

	void free(vector<BamTools::BamAlignment*>* a) {
		assert(a != 0);
		for (size_t i=0; i<a->size(); ++i) {
			assert(a->at(i) != 0);
			delete a->at(i);
		}
		delete a;
	}

	void run() {
		assert(alignments1 != 0);
			// compute alignment scores, i.e. AS tags for all alignments of first read
		for (size_t i=0; i<alignments1->size(); ++i) {
			assert(alignments1->at(i) != 0);
			if (!use_mapq) {
				as_tag_calculator.computeTag(*alignments1->at(i), bam_ref_data, &as_stats);
			}
		}
		if (alignments2 == 0) {
			// single-end case
			process_read_single_end(bam_ref_data, *alignments1, output, use_mapq);
		} else {
			// pairead-end case
			// compute alignment scores, i.e. AS tags for all alignments of second read
			for (size_t i=0; i<alignments2->size(); ++i) {
				assert(alignments2->at(i) != 0);
				if (!use_mapq) {
					as_tag_calculator.computeTag(*alignments2->at(i), bam_ref_data, &as_stats);
				}
			}
			// compute prior probabilities and prepare output
			process_read_allpairs(bam_ref_data, *alignments1, *alignments2, insert_size_dist, output, use_mapq);
		}
	}

} work_package_t;

typedef struct output_writer_t {
	ASTagCalculator::stats_t as_stats;
	BamTools::BamWriter* bam_writer;

	void write(auto_ptr<work_package_t> work) {
		assert(work.get() != 0);
		cout << work->output.str();
		as_stats += work->as_stats;
		if (bam_writer != 0) {
			write_alignments(work->alignments1);
			write_alignments(work->alignments2);
		}
	}

	void write_alignments(vector<BamTools::BamAlignment*>* a) {
		assert(bam_writer != 0);
		assert(a != 0);
		for (size_t i=0; i<a->size(); ++i) {
			assert(a->at(i) != 0);
			bam_writer->SaveAlignment(*(a->at(i)));
		}
	}

	output_writer_t() : bam_writer(0) {}
	output_writer_t(BamTools::BamWriter* bam_writer) : bam_writer(bam_writer) {}
} output_writer_t;

void read_readgroup_list(const string& filename, readgroupwise_distributions_t* result) {
	ifstream is(filename.c_str());
	if (is.fail()) {
		cerr << "Error: could not open \"" << filename << "\"." << endl;
		exit(1);
	}
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer_t;
	boost::char_separator<char> whitespace_separator(" \t");
	string line;
	int linenr = 0;
	while (getline(is,line)) {
		linenr += 1;
		tokenizer_t tokenizer(line, whitespace_separator);
		vector<string> tokens(tokenizer.begin(), tokenizer.end());
		if (tokens.size() != 2) {
			ostringstream oss;
			oss << "Error parsing read group list. Offending line: " << linenr;
			throw std::runtime_error(oss.str());
		}
		const string& readgroup = tokens[0];
		readgroupwise_distributions_t::const_iterator it = result->find(readgroup);
		if (it != result->end()) {
			ostringstream oss;
			oss << "Duplicate readgroup \"" << readgroup << "\" in file \"" << filename << "\"." << endl;
			throw std::runtime_error(oss.str());
		}
		HistogramBasedDistribution* distribution = new HistogramBasedDistribution(tokens[1]);
		(*result)[readgroup] = distribution;
	}
}

int main(int argc, char* argv[]) {
	// PARAMETERS
	bool dont_skip_non_xa = false;
	int max_span;
	int distribution_estimation_count;
	string insert_length_dist_filename = "";
	string rgwise_insert_length_dist_filename = "";
	string mean_sd_filename = "";
	int bad_score_threshold;
	int threads;
	double discard_threshold;
	bool ignore_xa = false;
	bool use_mapq = false;
	bool unsorted = false;
	bool single_end = false;

	po::options_description options_desc("Allowed options");
	options_desc.add_options()
		("dont_skip_non_xa,x", po::value<bool>(&dont_skip_non_xa)->zero_tokens(), "Do not skip reads for which other alignments exist (i.e. X0+X1>1, but no XA tag is present).")
		("ignore_xa", po::value<bool>(&ignore_xa)->zero_tokens(), "Do not use alternative alignments from XA tags.")
		("use_mapq", po::value<bool>(&use_mapq)->zero_tokens(), "Use MAPQ instead of recomputing probability from all alignments.")
		("phred_base,p", po::value<int>(&phred_base)->default_value(33), "Value to substract from ASCII code to get the PHRED quality.")
		("bad_alignment_threshold,b", po::value<int>(&bad_score_threshold)->default_value(1000), "Issue a warning when AS tag is above this value.")
		("max_span,s", po::value<int>(&max_span)->default_value(50000), "Maximal internal segment. Read pairs with larger internal segment will be ignored.")
		("discard_reads,d", po::value<double>(&discard_threshold)->default_value(0.0), "Discard \"concordant\" alignments within the given number of standard deviations (default: disabled).")
		("insert_size_dist,i", po::value<string>(&insert_length_dist_filename)->default_value(""), "Filename of known internal segment size distribution. If not given, this distribution is estimated.")
		("rg_insert_size_dist,r", po::value<string>(&rgwise_insert_length_dist_filename)->default_value(""), "Filename of read-group-wise known internal segment size distributions. Expects two-column text file: <readgroup> <distribution-filename>.")
		("dist_est_count,n", po::value<int>(&distribution_estimation_count)->default_value(5000000), "Number of uniquely mapping reads that are to be used to estimate internal segment size distribution.")
		("mean_and_sd,m", po::value<string>(&mean_sd_filename), "Write (robustly estimated) mean and standard deviation of main peak if internal segment size distribution to given filename.")
		("threads,T", po::value<int>(&threads)->default_value(0), "Number of threads (default: 0 = strictly single-threaded).")
		("unsorted", po::value<bool>(&unsorted)->zero_tokens(), "Don't assume reads to be sorted by position. In this case, alignments must be grouped by read, i.e. alignments of the same read (pair) must be in subsequent lines.")
		("single-end", po::value<bool>(&single_end)->zero_tokens(), "Process single-end reads (instead of pairs).")
	;

	if (argc<3) {
		usage(argv[0], options_desc);
	}
	string reference_filename(argv[argc-2]);
	string bam_input_filename(argv[argc-1]);
	argc -= 2;

	po::variables_map options;
	try {
		po::store(po::parse_command_line(argc, argv, options_desc), options);
		po::notify(options);
	} catch(exception& e) {
		cerr << "error: " << e.what() << "\n";
		return 1;
	}

	if ((discard_threshold != 0.0) && (insert_length_dist_filename.size()>0)) {
		cerr << "Options -i an -d cannot be used together." << endl;
		return 1;
	}

	if (use_mapq && !ignore_xa) {
		cerr << "When using option --use_mapq, then --ignore_xa must also be given" << endl;
		cerr << "Reason: XA tags do not contain MAPQ fields." << endl;
		return 1;
	}

	if ((insert_length_dist_filename.size() > 0) && (rgwise_insert_length_dist_filename.size() > 0)) {
		cerr << "Options -i and -r cannot be used at the same time." << endl;
		return 1;
	}

	if (single_end) {
		if (insert_length_dist_filename.size() > 0) {
			cerr << "Options --single-end and -i cannot be used at the same time." << endl;
			return 1;
		}
		if (rgwise_insert_length_dist_filename.size() > 0) {
			cerr << "Options --single-end and -r cannot be used at the same time." << endl;
			return 1;
		}
		if (mean_sd_filename.size() > 0) {
			cerr << "Options --single-end and -m cannot be used at the same time." << endl;
			return 1;
		}
	}

	clock_t clock_start = clock();

	// Create insert size distribution. Either read it from file (if given) or
	// estimate based on input
	HistogramBasedDistribution* insert_length_distribution = 0;
	readgroupwise_distributions_t* readgroupwise_distributions = 0;
	if (!single_end) {
		if (insert_length_dist_filename.size() > 0) {
			try {
				insert_length_distribution = new HistogramBasedDistribution(insert_length_dist_filename);
			} catch(exception& e) {
				cerr << "Error: " << e.what() << "\n";
				return 1;
			}
		} else if (rgwise_insert_length_dist_filename.size() > 0) {
			try {
				readgroupwise_distributions = new readgroupwise_distributions_t();
				read_readgroup_list(rgwise_insert_length_dist_filename, readgroupwise_distributions);
			} catch(exception& e) {
				cerr << "Error: " << e.what() << "\n";
				return 1;
			}
		} else {
			try {
				cerr << "Estimating internal segment size distribution" << endl;
				auto_ptr<BamReader> bam_reader(0);
				if (unsorted) {
					bam_reader = auto_ptr<BamReader>(new GroupWiseBamReader(bam_input_filename, true, true, !ignore_xa));
				} else {
					bam_reader = auto_ptr<BamReader>(new SortedBamReader(bam_input_filename, true, max_span, !ignore_xa));
				}
				Histogram histogram;
				int i = 0;
				while (bam_reader->hasNext() && (i<distribution_estimation_count)) {
					bam_reader->advance();
					if (bam_reader->isFirstUnmapped() || bam_reader->isSecondUnmapped()) continue;
					if (bam_reader->hasMultipleMappingsButNoXA()) continue;
					const vector<BamTools::BamAlignment*>& alignments1 = bam_reader->getAlignmentsFirst();
					const vector<BamTools::BamAlignment*>& alignments2 = bam_reader->getAlignmentsSecond();
					if (unsorted && ((alignments1.size() == 0) || (alignments2.size() == 0))) {
						cerr << "Error: When using option --unsorted, alignments of both reads in a pair must be" << endl;
						cerr << "       in consecutive lines. This is not the case for read \"" << bam_reader->getReadName() << "\"" << endl;
						return 1;
					}
					assert(alignments1.size() > 0);
					assert(alignments2.size() > 0);
					assert(alignments1[0] != 0);
					assert(alignments2[0] != 0);
					if ((alignments1.size() == 1) && (alignments2.size() == 1))
					if (alignments1[0]->IsReverseStrand() == alignments2[0]->IsReverseStrand()) continue;
					if (alignments1[0]->RefID != alignments2[0]->RefID) continue;
					//if (is_clipped(*alignments1[0]) || is_clipped(*alignments2[0])) continue;
					if (!is_unique(*alignments1[0]) || !is_unique(*alignments2[0])) continue;
					int insert_size = 0;
					if (alignments1[0]->Position <= alignments2[0]->Position) {
						insert_size = alignments2[0]->Position - alignments1[0]->GetEndPosition();
					} else {
						insert_size = alignments1[0]->Position - alignments2[0]->GetEndPosition();
					}
					histogram.add(insert_size);
					i += 1;
				}
				if (i < 1) {
					cerr << "Error: Too few reads estimate internal segment size distribution (" << i << ")." << endl;
					cerr << "Wont proceed with less than 1 uniquely mappable read pairs." << endl;
					return 1;
				}
				if (i < distribution_estimation_count) {
					cerr << "Warning: fewer uniquely mappable read (" << i << ") than asked for (" << distribution_estimation_count << ", see option -n)." << endl;
				}
				double mean;
				double sd;
				histogram.computeMeanAndStddev(&mean, &sd);
				cerr << "Main peak of internal segment length distribution: mean " << mean << ", sd " << sd << endl;
				if (mean_sd_filename.compare("") != 0) {
					ofstream ofs(mean_sd_filename.c_str());
					ofs << mean << " " << sd << endl;
					ofs.close();
				}
				concordant_min = (int)round(mean - discard_threshold * sd);
				concordant_max = (int)round(mean + discard_threshold * sd);
				std::auto_ptr<HistogramBasedDistribution> d = histogram.toDistribution(20);
				assert(d.get() != 0);
				insert_length_distribution = d.release();
			} catch(exception& e) {
				cerr << "Error: " << e.what() << "\n";
				return 1;
			}
		}
	}

	// Read reference genome and initialize objcet that computes alignment scores
	ASTagCalculator tag_calculator(reference_filename,phred_base,bad_score_threshold);
	cerr << "Done reading reference sequences" << endl;

	// count the number of alignments skipped because of missing XA tags
	long long skipped_by_xa = 0;
	BamReader* bam_reader = 0;
	typedef ThreadPool<work_package_t,output_writer_t> thread_pool_t;
	output_writer_t* output_writer = 0;
	try {
		if (unsorted || single_end) {
			bam_reader = new GroupWiseBamReader(bam_input_filename, !single_end, true, !ignore_xa);
		} else {
			bam_reader = new SortedBamReader(bam_input_filename, !single_end, max_span, !ignore_xa);
		}
		bam_reader->enableProgressMessages(cerr, 200000);
		output_writer = new output_writer_t();
		const BamTools::RefVector& bam_ref_data = bam_reader->getReferenceData();
		thread_pool_t thread_pool(threads, 1000, threads, *output_writer);
		while ( bam_reader->hasNext() ) {
			bam_reader->advance();
			if (single_end) {
				if (bam_reader->isUnmapped()) continue;
				if (!dont_skip_non_xa && bam_reader->hasMultipleMappingsButNoXA()) {
					skipped_by_xa += 1;
					continue;
				}
				vector<BamTools::BamAlignment*>* alignments = bam_reader->releaseAlignments().release();
				thread_pool.addTask(auto_ptr<work_package_t>(new work_package_t(alignments, 0, tag_calculator, bam_ref_data, *insert_length_distribution, use_mapq)));
			} else {
				if (bam_reader->isFirstUnmapped() || bam_reader->isSecondUnmapped()) continue;
				if (!dont_skip_non_xa && bam_reader->hasMultipleMappingsButNoXA()) {
					skipped_by_xa += 1;
					continue;
				}
				vector<BamTools::BamAlignment*>* alignments1 = bam_reader->releaseAlignmentsFirst().release();
				vector<BamTools::BamAlignment*>* alignments2 = bam_reader->releaseAlignmentsSecond().release();
				if (readgroupwise_distributions != 0) {
					if (alignments1->size() == 0) continue;
					string rg = "";
					if (!alignments1->at(0)->GetTag("RG", rg)) {
						cerr << "Error: no read group tag present for read " << alignments1->at(0)->Name << endl;
						return 1;
					}
					readgroupwise_distributions_t::const_iterator it = readgroupwise_distributions->find(rg);
					if (it == readgroupwise_distributions->end()) {
						cerr << "Error: read group \"" << rg << "\" of read " << alignments1->at(0)->Name << " is unknown." << endl;
						return 1;
					}
					// cerr << "Found read group specific distribution for read " << alignments1->at(0)->Name << ": p(300) = " << it->second->probability(300) << endl;
					thread_pool.addTask(auto_ptr<work_package_t>(new work_package_t(alignments1, alignments2, tag_calculator, bam_ref_data, *(it->second), use_mapq)));
				} else {
					thread_pool.addTask(auto_ptr<work_package_t>(new work_package_t(alignments1, alignments2, tag_calculator, bam_ref_data, *insert_length_distribution, use_mapq)));
				}
			}
		}
	} catch(exception& e) {
		cerr << "Error: " << e.what() << "\n";
		return 1;
	}
	if (skipped_by_xa > 0) {
		cerr << "Skipped " << skipped_by_xa << " ambiguously mapped reads for which no XA tag was present (to prevent this, use option -x)." << endl;
	}
	if (bam_reader->getSkippedDuplicates() > 0) {
		cerr << "Skipped " << bam_reader->getSkippedDuplicates() << " duplicate alignments." << endl;
	}
	if (bam_reader->getNonPairedCount() > 0) {
		cerr << "Skipped " << bam_reader->getNonPairedCount() << " reads with no matching mate alignment." << endl;
	}
	if (skipped_inside_clipping > 0) {
		cerr << "Skipped " << skipped_inside_clipping << " reads with alignments soft-clipped on the inside of an alignment pair!" << endl;
	}
	assert(output_writer != 0);
	tag_calculator.printWarnings(cerr, output_writer->as_stats);
	delete output_writer;
	if (bam_reader != 0) delete bam_reader;
	if (insert_length_distribution!=0) delete insert_length_distribution;
	if (readgroupwise_distributions != 0) {
		readgroupwise_distributions_t::const_iterator it = readgroupwise_distributions->begin();
		for (; it != readgroupwise_distributions->end(); ++it) delete it->second;
		delete readgroupwise_distributions;
	}
	double cpu_time = (double)(clock() - clock_start) / CLOCKS_PER_SEC;
	cerr << "Total CPU time: " << cpu_time << endl;
	return 0;
}
