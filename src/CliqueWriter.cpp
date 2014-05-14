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

CliqueWriter::CliqueWriter(ostream& os, VariationCaller* variation_caller, std::ostream* indel_os, const ReadGroups* read_groups, bool multisample, bool output_all, double fdr_threshold, bool verbose, int min_coverage, bool frameshift_merge, std::string suffix, int minimal_superread_length) : os(os), variation_caller(variation_caller), read_groups(read_groups) {
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
  this->output_position=0;
  this->minimal_superread_length = minimal_superread_length;
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

void CliqueWriter::callVariation(const vector<const AlignmentRecord*>& pairs, size_t coverage_, clique_stats_t &stats) {
  bool DEBUG = 0;
  bool merge = this->FRAMESHIFT_MERGE;
  bool problem = 0;
  assert(&stats != 0);
  VariationCaller::additional_stats_t vc_stats;
  if (variation_caller != 0) {
    stats.variation = variation_caller->call(pairs.begin(), pairs.end(), &vc_stats);
  } else {
    stats.variation = Variation();
  }
  stats.min_coverage_user = min_coverage;
  stats.clique_size = pairs.size();
  stats.start = vc_stats.insert_start;
  stats.end = vc_stats.insert_end;
  stats.length = vc_stats.insert_length;
  stats.diff = vc_stats.diff;
  stats.total_weight = vc_stats.total_weight;
  stats.coverage = coverage_;
  stats.pvalue_corr = min(1.0, stats.variation.getPValue() * pow(1.2, static_cast<int> (stats.coverage)));
  stats.is_significant = false;
  vector<const AlignmentRecord*>::const_iterator it;
  for (it = pairs.begin(); it != pairs.end(); ++it) {
    const AlignmentRecord& ap = **it;
    // //cerr << ap.getName() << endl;
    std::set<string>::iterator it_s;
    if (!ap.getReadNames().empty()) {
      set<string> from = ap.getReadNames();
      std::copy(from.begin(), from.end(), std::inserter( stats.readnames, stats.readnames.begin()));
    }
    stats.clique_size_weighted += ap.getCount();
    stats.hcount += ap.getHCount();
  }
  //cerr << "WEIGHT: " << stats.clique_size_weighted << endl;

  if (stats.clique_size_weighted < min_coverage) {
    return;
  }

  stats.clique_number = clique_count++;
  bool has_paired_end = 0;
  for (it = pairs.begin(); it != pairs.end(); ++it) {
    const AlignmentRecord& ap = **it;
    vector<BamTools::CigarOp>::const_iterator it_cigar = ap.getCigar1().begin();
    int end = ap.getStart1() + ap.getLengthInclDeletions1();
    if (ap.getStart1() < stats.window_start1 || stats.window_start1 == -1) {
      stats.window_start1 = ap.getStart1();
    }
    if (end > stats.window_end1 || stats.window_end1 == -1) {
      stats.window_end1 = end;
    }
    if (ap.isPairedEnd()) {
      it_cigar = ap.getCigar2().begin();
      int end = ap.getStart2() + ap.getLengthInclDeletions2();
      if (ap.getStart2() < stats.window_start2 || stats.window_start2 == -1) {
        stats.window_start2 = ap.getStart2();
      }
      if (end > stats.window_end2 || stats.window_end2 == -1) {
        stats.window_end2 = end ;
      }
      has_paired_end = 1;
    }
  }


  //2D array represents nucleotide distribution for super-read
  int alignment_length1 = stats.window_end1 - stats.window_start1;
  if (alignment_length1 < this->minimal_superread_length) {
    return;
  }
  int alignment1[alignment_length1][5];
  double phred1[alignment_length1][5];
  char match1[alignment_length1][3];

  int alignment_length2 = 0;
  if (has_paired_end) {
    alignment_length2 = stats.window_end2 - stats.window_start2+1;
    if (alignment_length2 < this->minimal_superread_length) {
      return;
    }
  }
  int alignment2[alignment_length2][5];
  double phred2[alignment_length2][5];
  char match2[alignment_length2][3];
  int max_length = std::max(alignment_length1,alignment_length2);

  for (int i = 0; i < max_length; ++i) {
    if (i < alignment_length1) {
      alignment1[i][0] = 0;
      alignment1[i][1] = 0;
      alignment1[i][2] = 0;
      alignment1[i][3] = 0;
      alignment1[i][4] = 0;

      phred1[i][0] = 0;
      phred1[i][1] = 0;
      phred1[i][2] = 0;
      phred1[i][3] = 0;
      phred1[i][4] = 0;

      match1[i][0] = 0;
      match1[i][1] = 0;
      match1[i][2] = 0;
    }
    if (i < alignment_length2) {
      alignment2[i][0] = 0;
      alignment2[i][1] = 0;
      alignment2[i][2] = 0;
      alignment2[i][3] = 0;
      alignment2[i][4] = 0;

      phred2[i][0] = 0;
      phred2[i][1] = 0;
      phred2[i][2] = 0;
      phred2[i][3] = 0;
      phred2[i][4] = 0;

      match2[i][0] = 0;
      match2[i][1] = 0;
      match2[i][2] = 0;
    }
  }

    // //cerr << endl << "===" << endl;
  for (it = pairs.begin(); it != pairs.end(); ++it) {
    const AlignmentRecord& ap = **it;
    // //cerr << ap.getName() << endl;
    vector<BamTools::CigarOp>::const_iterator it_cigar = ap.getCigar1().begin();
    int alignment_index = ap.getStart1() - stats.window_start1;
    // //cerr << "start: " << ap.getStart1() << "\t" << stats.window_start1 << endl;
    assert(alignment_index >= 0);
    assert(alignment_index < stats.window_end1);
    int sequence_index = 0;

    for (; it_cigar != ap.getCigar1().end(); ++it_cigar) {
      int cigar_length = it_cigar->Length;
      int cigar_type = 0;
      switch (it_cigar->Type) {
        case 'S':
          sequence_index += cigar_length;
          break;
        case 'D':
          cigar_type = 1;
          if (merge && cigar_length == 1) {
            alignment_index++;
          } else {
            for (int k = 0; k < cigar_length; k++) {
              assert(alignment_index < alignment_length1);
              match1[alignment_index][cigar_type]++;
              alignment1[alignment_index++][4]+=ap.getCount();
            }
          }
          break;
        case 'I':
          if (cigar_length == 1 && merge) {
            sequence_index++;
            break;
          }
          cigar_type = 2;
        case 'M':
          for (int k = 0; k < cigar_length; k++) {
            match1[alignment_index][cigar_type]++;
            assert(sequence_index < ap.getSequence1().size());
            assert(alignment_index < alignment_length1);

            int base = shortenBase(ap.getSequence1()[sequence_index++]);
            assert(base <= 4);
            if (base != -1) {
              phred1[alignment_index][base] += ap.getSequence1().qualityCorrectLog(sequence_index - 1);
              if (alignment1[alignment_index][base] + ap.getCount() < 0) {
                //cerr << "GAME OVER" << endl;
                //cerr << "alignment_index: " << alignment_index << "\tbase: " << base << "\tcount: " << ap.getCount() << "\talignment2: " << alignment1[alignment_index][base] << endl;
                exit(0);
              }
              alignment1[alignment_index][base]+=ap.getCount();
            }

            alignment_index++;
          }
          break;
      }
    }

    if (ap.isSingleEnd()) continue;
    it_cigar = ap.getCigar2().begin();

    alignment_index = ap.getStart2() - stats.window_start2;
    sequence_index = 0;

    for (; it_cigar != ap.getCigar2().end(); ++it_cigar) {
      int cigar_length = it_cigar->Length;
      int cigar_type = 0;
      switch (it_cigar->Type) {
        case 'S':
          sequence_index += cigar_length;
          break;
        case 'D':
          cigar_type = 1;
          if (merge && cigar_length == 1) {
            alignment_index++;
          } else {
            for (int k = 0; k < cigar_length; k++) {
              match2[alignment_index][cigar_type]++;
              assert(alignment_index < alignment_length2);
              alignment2[alignment_index++][4]+=ap.getCount();
            }
          }
          break;
        case 'I':
          if (cigar_length == 1 && merge) {
            sequence_index++;
            break;
          }
          cigar_type = 2;
        case 'M':
          for (int k = 0; k < cigar_length; k++) {
            match2[alignment_index][cigar_type]++;
            assert(sequence_index < ap.getSequence2().size());
            assert(alignment_index < alignment_length2);

            int base = shortenBase(ap.getSequence2()[sequence_index++]);
            assert(base <= 4);
            if (base != -1) {
              phred2[alignment_index][base] += ap.getSequence2().qualityCorrectLog(sequence_index - 1);
              if (alignment2[alignment_index][base] + ap.getCount() < 0) {
                //cerr << "GAME OVER" << endl;
                //cerr << "alignment_index: " << alignment_index << "\tbase: " << base << "\tcount: " << ap.getCount() << "\talignment2: " << alignment2[alignment_index][base] << endl;
                exit(0);
              }
              alignment2[alignment_index][base]+=ap.getCount();
            }

            alignment_index++;
          }
          break;
      }
    }
  }

  for (int i = 0; i < alignment_length1; ++i) {
    //cerr << alignment1[i][0] << " " << alignment1[i][1] << " " << alignment1[i][2] << " " << alignment1[i][3] << " " << alignment1[i][4] << endl;
  }

  //Compute maximum coverage
  int coverage1[alignment_length1];
  std::fill_n(coverage1, alignment_length1, 0);
  int max_base1[alignment_length1];
  std::fill_n(max_base1, alignment_length1, 0);
  int max_base_weight1[alignment_length1];
  std::fill_n(max_base_weight1, alignment_length1, 0);
  int max_match1[alignment_length1];
  std::fill_n(max_match1, alignment_length1, 0);
  int max_match_weight1[alignment_length1];
  std::fill_n(max_match_weight1, alignment_length1, 0);
  int coverage2[alignment_length2];
  std::fill_n(coverage2, alignment_length2, 0);
  int max_base2[alignment_length2];
  std::fill_n(max_base2, alignment_length2, 0);
  int max_base_weight2[alignment_length2];
  std::fill_n(max_base_weight2, alignment_length2, 0);
  int max_match2[alignment_length2];
  std::fill_n(max_match2, alignment_length2, 0);
  int max_match_weight2[alignment_length2];
  std::fill_n(max_match_weight2, alignment_length2, 0);
  for (int i = 0; i < max_length; ++i) {
    if (i < alignment_length1) {
      max_base1[i] = -1;
      max_base_weight1[i] = 0;
      int sum = 0;
      for (int k = 0; k < 5; k++) {
        sum += alignment1[i][k];
        if (alignment1[i][k] > max_base_weight1[i] && alignment1[i][k] > 0) {
          max_base_weight1[i] = alignment1[i][k];
          max_base1[i] = k;
        }
      }
      max_match1[i] = -1;
      if (match1[i][0] >= match1[i][1] && match1[i][0] >= match1[i][2] && match1[i][0] > 0) {
        max_match1[i] = 0;
      } else if (match1[i][1] >= match1[i][0] && match1[i][1] >= match1[i][2] && match1[i][1] > 0) {
        max_match1[i] = 1;
      } else if (match1[i][2] >= match1[i][0] && match1[i][2] >= match1[i][1] && match1[i][2] > 0) {
        max_match1[i] = 2;
      }
      max_base_weight1[i] = match1[i][0]+match1[i][1]+match1[i][2];
      coverage1[i] = sum;
    }
    if (i < alignment_length2) {
      max_base2[i] = -1;
      max_base_weight2[i] = 0;
      int sum = 0;
      for (int k = 0; k < 5; k++) {
        sum += alignment2[i][k];
        if (alignment2[i][k] > max_base_weight2[i] && alignment2[i][k] > 0) {
          max_base_weight2[i] = alignment2[i][k];
          max_base2[i] = k;
        }
      }
      max_match2[i] = -1;
      if (match2[i][0] >= match2[i][1] && match2[i][0] >= match2[i][2] && match2[i][0] > 0) {
        max_match2[i] = 0;
      } else if (match2[i][1] >= match2[i][0] && match2[i][1] >= match2[i][2] && match2[i][1] > 0) {
        max_match2[i] = 1;
      } else if (match2[i][2] >= match2[i][0] && match2[i][2] >= match2[i][1] && match2[i][2] > 0) {
        max_match2[i] = 2;
      }
      max_base_weight2[i] = match2[i][0]+match2[i][1]+match2[i][2];
      coverage2[i] = sum;
    }
  }

  for (int i = 0; i < alignment_length1; ++i) {
    //cerr << coverage1[i] << endl;
  }


  bool prefix1 = 1;
  int end1 = 0;
  for (int i = alignment_length1 - 1; i >= 0; i--) {
    if (coverage1[i] >= min_coverage) {
      end1 = i;
      break;
    }
  }

  char current_cigar = 0;
  int current_cigar_count = 0;
  for (int i = 0; i < alignment_length1; i++) {
    if (i > end1) { break; }

    if (prefix1) {
      if (coverage1[i] >= min_coverage) {
        prefix1 = 0;
        stats.window_start1 += i;
      } else {
        continue;
      }
    }

    addCigar(current_cigar, current_cigar_count, max_match1[i], stats, 1, i == end1);

    if (max_match1[i] == 0 || max_match1[i] == 2) {
      if (max_base1[i] < 4) {
        stats.consensus_string1 += expandBase(max_base1[i]);
        int p = 33 + round(-10 * phred1[i][max_base1[i]]);
        if (p > 88) {
          p = 88;
        } else if (p < 33) {
          p = 88;
        }
        stats.phred_string1 += p;
      }
    }
  }

  //same stuff as for first read, this time for second read
  if (has_paired_end) {
    this->paired_count += 1;

    bool prefix2 = 1;
    int end2 = 0;
    current_cigar = 0;
    current_cigar_count = 0;
    for (int i = alignment_length2 - 1; i >= 0; i--) {
      if (coverage2[i] >= min_coverage) {
        end2 = i;
        break;
      }
    }

    for (int i = 0; i < alignment_length2; i++) {
      if (i > end2) break;

      if (prefix2) {
        if (coverage2[i] >= min_coverage) {
          prefix2 = 0;
          stats.window_start2 += i;
        } else {
          continue;
        }
      }
      
      addCigar(current_cigar, current_cigar_count, max_match2[i], stats, 2, i == end2);

      if (max_match2[i] == 0 || max_match2[i] == 2) {
        if (max_base2[i] < 4) {
          stats.consensus_string2 += expandBase(max_base2[i]);
          int p = 33 + round(-10 * phred2[i][max_base2[i]]);
          if (p > 88) {
            p = 88;
          } else if (p < 33) {
            p = 88;
          }
          stats.phred_string2 += p;
        }
      }  
    }
  }

  if (stats.consensus_string1.size() < this->minimal_superread_length) {
    stats.window_end1 = -1;
    stats.window_start1 = -1;
    stats.window_end2 = -1;
    stats.window_start1 = -1;
    stats.consensus_string1 = "";
    stats.consensus_string2 = "";
    stats.phred_string1 = "";
    stats.phred_string2 = "";
    stats.clique_number = 0;
  }
  if (overlapSize(stats)) {
    error("merge");
  }


  //cerr << "XLR" << endl;
  //cerr << stats.window_start1 << " " << stats.cigar_string1 << " " << stats.consensus_string1 << " " << stats.phred_string1 << endl;
  if (has_paired_end) {
    //cerr << stats.window_start2 << " " << stats.cigar_string2 << " " << stats.consensus_string2 << " " << stats.phred_string2 << endl;
  }

  //cerr << endl << endl;
}

bool CliqueWriter::overlapSize(clique_stats_t& stats) const {
  if (stats.window_start2 > 0) {
    int s1 = stats.window_start1;
    int s2 = stats.window_start2;
    int e1 = s1 + stats.consensus_string1.size();
    int e2 = s2 + stats.consensus_string2.size();
    if (s2 > e1) {
      // -----
      //       -----
      return 0;
    }
    if (s1 > e2) {
      //       -----
      // -----
      string cons_tmp = stats.consensus_string1;
      string phred_tmp = stats.phred_string1;
      int end_tmp = stats.window_end1;
      int start_tmp = stats.window_start1;

      stats.consensus_string1 = stats.consensus_string2;
      stats.phred_string1 = stats.phred_string2;
      stats.window_end1 = stats.window_end2;
      stats.window_start1 = stats.window_start2;

      stats.consensus_string2 = cons_tmp;
      stats.phred_string2 = phred_tmp;
      stats.window_end2 = end_tmp;
      stats.window_start2 = start_tmp;
      return 0;
    }
    if (s1 == s2 && e1 == e2) {
      // -----
      // -----
    } else if (s1 == s2) {
      // ----  AND -----
      // ----- AND ----
      if (e1 < e2) {
        stats.consensus_string1 = stats.consensus_string2;
        stats.phred_string1 = stats.phred_string2;
        stats.window_end1 = stats.window_end2;
      } else {
      }
    } else if (e1 == e2) {
      //  ---- AND -----
      // ----- AND  ----
      if (s1 > s2) {
        stats.consensus_string1 = stats.consensus_string2;
        stats.phred_string1 = stats.phred_string2;
        stats.window_start1 = stats.window_start2;
      } else {
      }
    } else if (e1 < e2) {
      // ----    AND   --
      //   ----  AND ------
      if (s1 < s2) {
        if (e1 - s2 >= 1) {
          string cons = equalStrings(stats.consensus_string1, stats.consensus_string2);
          if (cons.size() > 0) {
            for (int i = stats.consensus_string2.size()-(cons.size() - stats.consensus_string1.size()); i < stats.consensus_string2.size(); i++) {
              stats.phred_string1 += stats.phred_string2.at(i);
            }
            stats.consensus_string1 = cons;
            stats.window_end1 = stats.window_start1 + cons.size();
            if (stats.phred_string1.size() != cons.size()) {
              error(boost::lexical_cast<std::string>(stats.phred_string1.size()) + " " + boost::lexical_cast<std::string>(cons.size()));
            }
            assert(stats.phred_string1.size() == cons.size());
          }
        }
      } else {
        stats.consensus_string1 = stats.consensus_string2;
        stats.phred_string1 = stats.phred_string2;
        stats.window_end1 = stats.window_end2;
        stats.window_start1 = stats.window_start2;
      }
    } else {
      //   ---- AND ------
      // ----   AND   --
      if (s2 < s1) {
        string cons = equalStrings(stats.consensus_string2, stats.consensus_string1);
        if (cons.size() > 0) {
          for (int i = stats.consensus_string1.size()-(cons.size() - stats.consensus_string2.size()); i < stats.consensus_string1.size(); i++) {
            stats.phred_string2 += stats.phred_string1.at(i);
          }
          stats.phred_string1 = stats.phred_string2;
          stats.consensus_string1 = cons;
          stats.window_start1 = stats.window_start2;
          stats.window_end1 = stats.window_start1 + cons.size();
          if (stats.phred_string1.size() != cons.size()) {
            error(boost::lexical_cast<std::string>(stats.phred_string1.size()) + " " + boost::lexical_cast<std::string>(cons.size()));
          }
          assert(stats.phred_string1.size() == cons.size());
        }
      } else {
        //     // ------
        //     //   --
      }
    }
  }

  stats.window_end2 = 0;
  stats.window_start2 = 0;
  stats.phred_string2 = "";
  stats.consensus_string2 = "";

  int split = 0;
  for (int i = 0; i < stats.consensus_string1.size(); i++) {
    if (stats.consensus_string1.at(i) == 'N') {
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
    for (int i = 0; i < stats.consensus_string1.size(); i++) {
      if (prefix) {
        if (stats.consensus_string1.at(i) == 'N') {
          split = 1;
          prefix = 0;
          tmp_end = stats.window_start1 + i;
        } else {
          tmp_string1 += stats.consensus_string1.at(i);
          tmp_phred1 += stats.phred_string1.at(i);
        }
      } else {
        if (prefix == 0 && stats.consensus_string1.at(i) != 'N') {
          suffix = 1;
          stats.window_start2 = stats.window_start1 + i;
          stats.window_end2 = stats.window_end1;
        }
        if (suffix) {
          stats.consensus_string2 += stats.consensus_string1.at(i);
          stats.phred_string2 += stats.phred_string1.at(i);

          if (stats.consensus_string1.at(i) == 'N') {
            error(boost::lexical_cast<std::string>(i) + " " + boost::lexical_cast<std::string>(stats.consensus_string1.size()));
            return 1;
          }
        }
      }
    }
    stats.window_end1 = stats.window_start1 + stats.consensus_string1.size();
    stats.consensus_string1 = tmp_string1;

    error("SPLIT");
  }

  return 0;
}

string CliqueWriter::equalStrings(string &s1, string &s2) const {
  int e1 = s1.size();
  int e2 = s2.size();
  int m = min(e2,e1);
  for (int i = 0; i < m; i++) {
    bool fit = 1;
    int del = 0;
    for (int j = 0; j <= i; j++) {
      if (-1-i+j >= 0) { error("NOOO");break; }
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

void CliqueWriter::addCigar(char& current_cigar, int& current_cigar_count, int match, clique_stats_t& stats, int strand, bool last) const {
  char now_cigar_char = 0;
  switch (match) {
    case 0:
     now_cigar_char = 'M';
     break;
    case 1:
     now_cigar_char = 'D';
     break;
    case 2:
     now_cigar_char = 'I';
     break;
  }
  if (current_cigar == 0 && !last) {
    current_cigar = now_cigar_char;
    current_cigar_count = 1;
  } else if (current_cigar == now_cigar_char && !last) {
    current_cigar_count++;
  } else {
    if (last) {
      current_cigar_count++;
    }
    string cigar_string = boost::lexical_cast<std::string>(current_cigar_count) + current_cigar;
    if (strand == 1) {
      stats.cigar_string1 += cigar_string;
    } else{
      stats.cigar_string2 += cigar_string;
    }
    if (!last) {
      current_cigar = now_cigar_char;
      current_cigar_count = 1;
    }
  }
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
    error("Not aware of base: " + boost::lexical_cast<std::string>(base));
    exit(0);
    break;
  }
}

void CliqueWriter::add(std::auto_ptr<Clique> clique) {
  assert(!finished);
  clique_stats_t stats;
  auto_ptr<vector<const AlignmentRecord*> > all_pairs = clique->getAllAlignments();
  assert(all_pairs->size() == clique->size());

  callVariation(*all_pairs, clique->totalCenterCoverage(), stats);

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
    // //cerr << endl << stats.cigar_string1 << endl;
    string key;
    key = stats.consensus_string1+"=";
    if (!stats.consensus_string2.empty()) {
      key += stats.consensus_string2;
    }
    if (fastq_map.find(key) == fastq_map.end()) {
      fastq_map[key] = fastq_entry();
      fastq_map[key].name = boost::lexical_cast<string>(stats.clique_number);
      fastq_map[key].read_names = new set<string>();

      set<string> from = stats.readnames;
      if (!stats.readnames.empty()) copy(from.begin(), from.end(), inserter(*fastq_map[key].read_names, fastq_map[key].read_names->begin()));
      fastq_map[key].pos_1 = stats.window_start1;
      fastq_map[key].seq_1 = stats.consensus_string1;
      fastq_map[key].phreds_1.push_back(stats.phred_string1);
      fastq_map[key].hcount = stats.hcount;
      if (!stats.consensus_string2.empty()) {
        fastq_map[key].pos_2 = stats.window_start2;
        fastq_map[key].seq_2 = stats.consensus_string2;
        fastq_map[key].phreds_2.push_back(stats.phred_string2);
      }
    } else {
      set<string> from = stats.readnames;
      if (!stats.readnames.empty()) copy(from.begin(), from.end(), inserter(*fastq_map[key].read_names, fastq_map[key].read_names->begin()));
      fastq_map[key].phreds_1.push_back(stats.phred_string1);
      if (!stats.consensus_string2.empty()) {
        fastq_map[key].phreds_2.push_back(stats.phred_string2);
      }
      fastq_map[key].hcount += stats.hcount;
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
  // //cerr << "PRINT: " << pos_1 << endl;
  if (fastq_map.empty()) {
    return;
  }
  if (pos_1 < this->output_position+100 && pos_1 != -1) {
    return;
  } else {
    this->output_position=pos_1;
  }
  ofstream fr1;
  fr1.open ("data_cliques_paired_R1.fastq", ios::out | ios::app);
  ofstream fr2;
  fr2.open ("data_cliques_paired_R2.fastq", ios::out | ios::app);
  ofstream fs;
  fs.open ("data_cliques_single.fastq", ios::out | ios::app);
  ofstream fc;
  fc.open ("data_clique_to_reads.tsv", ios::out | ios::app);

  std::map<std::string,fastq_entry>::iterator it;
  for (it = fastq_map.begin();it!=fastq_map.end();) {
    fastq_entry f;
    f = it->second;
    if (f.pos_1+100 < pos_1 || pos_1 == -1) {
      if ((f.read_names->size() + f.hcount) >= this->min_coverage) {
        ofstream* out1;
        if (f.seq_2.empty()) {
          out1 = &fs;
        } else {
          out1 = &fr1;
        }
        *out1 << "@Clique_" << this->suffix << f.name << endl;
        std::set<std::string>::iterator set_it;
        set_it = f.read_names->begin();
        fc << "Clique_" << this->suffix << f.name << "\t";
        if (!f.read_names->empty()) {
          fc << *set_it;
          ++set_it;
          for (;set_it!=f.read_names->end();++set_it) {
            fc << "," << *set_it;
          }
        }
        if (f.hcount != 0) {
          if (!f.read_names->empty()) {
            fc << ",";
          }
          fc << "HCOUNT|" << f.hcount;
        }
        fc << endl;
        *out1 << f.seq_1 << endl;
        //*out1 << "+" << f.pos_1 << endl;
        *out1 << "+" << endl;
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
              error("UNKNOWN BASE " + boost::lexical_cast<std::string>(f.seq_2[i]));
              break;
            }
          }
          fr2 << endl;
          fr2 << "+" << f.pos_2 << endl;
          fr2 << f.phreds_2[0] << endl;
        }
        this->single_count++;
      } else {
        error("Too small " + boost::lexical_cast<std::string>(f.read_names->size()) + " "+ boost::lexical_cast<std::string>(f.hcount));
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

void CliqueWriter::error(std::string s) const {
  std::cerr << "    .o oOOOOOOOo                                            OOOo\n    Ob.OOOOOOOo  OOOo.      oOOo.                      .adOOOOOOO\n    OboO\"\"\"\"\"\"\"\"\"\"\"\".OOo. .oOOOOOo.    OOOo.oOOOOOo..\"\"\"\"\"\"\"\"\"'OO\n    OOP.oOOOOOOOOOOO \"POOOOOOOOOOOo.   `\"OOOOOOOOOP,OOOOOOOOOOOB'\n    `O'OOOO'     `OOOOo\"OOOOOOOOOOO` .adOOOOOOOOO\"oOOO'    `OOOOo\n    .OOOO'            `OOOOOOOOOOOOOOOOOOOOOOOOOO'            `OO\n    OOOOO                 '\"OOOOOOOOOOOOOOOO\"`                oOO\n   oOOOOOba.                .adOOOOOOOOOOba               .adOOOOo.\n  oOOOOOOOOOOOOOba.    .adOOOOOOOOOO@^OOOOOOOba.     .adOOOOOOOOOOOO\n OOOOOOOOOOOOOOOOO.OOOOOOOOOOOOOO\"`  '\"OOOOOOOOOOOOO.OOOOOOOOOOOOOO\n \"OOOO\"       \"YOoOOOOMOIONODOO\"`  .   '\"OOROAOPOEOOOoOY\"     \"OOO\"\n    Y           'OOOOOOOOOOOOOO: .oOOo. :OOOOOOOOOOO?'         :`\n    :            .oO%OOOOOOOOOOo.OOOOOO.oOOOOOOOOOOOO?         .\n    .            oOOP\"%OOOOOOOOoOOOOOOO?oOOOOO?OOOO\"OOo\n                 '%o  OOOO\"%OOOO%\"%OOOOO\"OOOOOO\"OOO':\n                      `$\"  `OOOO' `O\"Y ' `OOOO'  o             .\n    .                  .     OP\"          : o     .\n                              :\n                              .\n\nPlease create an issue at https://github.com/armintoepfer/haploclique/issues\nor write a mail armin.toepfer@gmail.com" << std::endl;
  std::cerr << s<< std::endl;
}
