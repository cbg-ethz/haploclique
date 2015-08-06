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

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include "BamHelper.h"

using namespace std;

void BamHelper::parseCigar(const string& cigar, vector<BamTools::CigarOp>* target) {
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer_t;
	boost::char_separator<char> cigar_separator("", "MIDNSHP=X");
	tokenizer_t tokenizer(cigar, cigar_separator);
	tokenizer_t::const_iterator token_it = tokenizer.begin();
	while (token_it != tokenizer.end()) {
		int length = boost::lexical_cast<int>(*token_it);
		++token_it;
		assert(token_it != tokenizer.end());
		assert(token_it->size() == 1);
		char type = token_it->at(0);
		target->push_back(BamTools::CigarOp(type,length));
		++token_it;
	}
}

void BamHelper::getSubalignments(const BamTools::BamAlignment& aln, vector<subalignment_t>* target) {
	if (!aln.IsMapped()) return;
	uint32_t as_tag = 0;
	if (!aln.GetTag("AS", as_tag)) {
		ostringstream oss;
		oss << "Error: No AS tag present for read \"" << aln.Name << "\"." << endl;
		throw std::runtime_error(oss.str());
	}
	uint32_t ym_tag = 0;
	if (!aln.GetTag("YM", ym_tag)) {
		ostringstream oss;
		oss << "Error: No YM tag present for read \"" << aln.Name << "\"." << endl;
		throw std::runtime_error(oss.str());
	}
	assert(target != 0);
	target->push_back(subalignment_t(aln.CigarData, as_tag, ym_tag));
	string ya_tag = "";
	if (aln.GetTag("YA", ya_tag)) {
		typedef boost::tokenizer<boost::char_separator<char> > tokenizer_t;
		boost::char_separator<char> separator(",;");
		tokenizer_t tokenizer(ya_tag, separator);
		tokenizer_t::const_iterator token_it = tokenizer.begin();
		while (token_it != tokenizer.end()) {
			int phred_score = boost::lexical_cast<int>(*token_it);
			++token_it;
			assert(token_it != tokenizer.end());
			int mismatch_phred_score = boost::lexical_cast<int>(*token_it);
			++token_it;
			assert(token_it != tokenizer.end());
			target->push_back(subalignment_t());
			subalignment_t& subalignment = target->at(target->size()-1);
			subalignment.phred_score = phred_score;
			subalignment.mismatch_phred_score = mismatch_phred_score;
			parseCigar(*token_it, &subalignment.cigar);
			++token_it;
			// if number of subalignments should be excessive, just ignore the rest
			if (target->size() == numeric_limits<unsigned short>::max()) break;
		}
	}
}

unique_ptr<BamHelper::read_t> BamHelper::createReadRecord(const vector<BamTools::BamAlignment*>& alignments1, const vector<BamTools::BamAlignment*>& alignments2) {
	unique_ptr<BamHelper::read_t> result(new read_t(alignments1, alignments2));
	for (size_t i=0; i<alignments1.size(); ++i) {
		result->subalignments1.push_back(vector<subalignment_t>());
		vector<subalignment_t>& v = result->subalignments1[i];
		getSubalignments(*(alignments1[i]), &v);
	}
	for (size_t i=0; i<alignments2.size(); ++i) {
		result->subalignments2.push_back(vector<subalignment_t>());
		vector<subalignment_t>& v = result->subalignments2[i];
		getSubalignments(*(alignments2[i]), &v);
	}
	return result;
}

/** Recompute PHRED score based on insertion/deletion length distributions, CIGAR string, and known PHRED scores of mismatches. */
void BamHelper::recalibratePhredScore(subalignment_t* subalignment, const BamTools::BamAlignment& aln, const IndelLengthDistribution& insertion_costs, const IndelLengthDistribution& deletion_costs, PositionSet* snp_set, int phred_offset, VariationIndex* variation_set) {
	assert(aln.IsMapped());
	assert(subalignment != 0);
	int indel_score = 0;
	int mismatch_score = 0;
	vector<BamTools::CigarOp>::const_iterator it = subalignment->cigar.begin();
	unsigned int ref_pos = aln.Position;
	unsigned int query_pos = 0;
	for (; it != subalignment->cigar.end(); ++it) {
		if (it->Type == 'I') {
			if (variation_set == 0) {
				indel_score += insertion_costs.getPhredCost(it->Length);
			} else {
				unique_ptr<vector<size_t> > variations = variation_set->getInsertions(aln.RefID, ref_pos);
				bool found = false;
				if (variations.get() != 0) {
					string seq = aln.QueryBases.substr(query_pos, it->Length);
					for (size_t j=0; j<variations->size(); ++j) {
						const Variation& v = variation_set->get(variations->at(j));
						assert(v.getType() == Variation::INSERTION);
						assert((unsigned int)v.getCoordinate1() == ref_pos);
						if (seq.compare(v.getSequence()) == 0) {
							found = true;
							break;
						}
					}
				}
				// TODO: account for similar insertions / sequencing errors in insertions
				// only add costs if this insertion has not been found
				// if known, then just add cost of 1 to still give (slight) preference to 
				// alignments without indel
				if (!found) {
					indel_score += insertion_costs.getPhredCost(it->Length);
				} else {
					indel_score += 1;
					// cerr << "FREE insertion" << endl;
				}
			}
			query_pos += it->Length;
		}
		if (it->Type == 'D') {
			if (variation_set == 0) {
				indel_score += deletion_costs.getPhredCost(it->Length);
			} else {
				unique_ptr<vector<size_t> > variations = variation_set->getDeletions(aln.RefID, ref_pos, ref_pos + it->Length);
				// only add costs if no matching variation has been found
				// if known, then just add cost of 1 to still give (slight) preference to 
				// alignments without indel
				if ((variations.get() == 0) || (variations->size() == 0)) {
					indel_score += deletion_costs.getPhredCost(it->Length);
				} else {
					indel_score += 1;
					// cerr << "FREE deletion" << endl;
				}
			}
			ref_pos += it->Length;
		}
		if (it->Type == 'M') {
			if (snp_set != 0) {
				throw std::runtime_error("Encountered 'M' in CIGAR string although 'X'/'=' are necessary to recalibrate against SNP set.");
			}
			ref_pos += it->Length;
			query_pos += it->Length;
		}
		if (it->Type == 'S') {
			query_pos += it->Length;
		}
		if (it->Type == '=') {
			ref_pos += it->Length;
			query_pos += it->Length;
		}
		if (it->Type == 'X') {
			for (size_t i=0; i<it->Length; ++i) {
				if ((snp_set==0) || (!snp_set->get(aln.RefID,ref_pos))) {
					assert(aln.Qualities[query_pos] >= phred_offset); 
					mismatch_score += aln.Qualities[query_pos] - phred_offset;
				}
				ref_pos += 1;
				query_pos += 1;
			}
		}
	}
	// cerr << "Recalibrating " << subalignment->cigar << ": " << subalignment->phred_score << "/" << subalignment->mismatch_phred_score << " --> ";
	// if snp_set is present, then also mismatch scores were reevaluated,
	// if not, we take existing mismatch score and add indel scores.
	if (snp_set == 0) {
		subalignment->phred_score = subalignment->mismatch_phred_score + indel_score;
	} else {
		subalignment->mismatch_phred_score = mismatch_score;
		subalignment->phred_score = mismatch_score + indel_score;
	}
	// cerr << subalignment->phred_score << "/" << subalignment->mismatch_phred_score << endl;
}

void BamHelper::recalibratePhredScores(vector<BamHelper::subalignment_t>* subalignments, const BamTools::BamAlignment& aln, const IndelLengthDistribution& insertion_costs, const IndelLengthDistribution& deletion_costs, PositionSet* snp_set, int phred_offset, VariationIndex* variation_set) {
	assert(subalignments != 0);
	for (size_t i=0; i<subalignments->size(); ++i) {
		recalibratePhredScore(&subalignments->at(i), aln, insertion_costs, deletion_costs, snp_set, phred_offset, variation_set);
	}
	sort(subalignments->begin(), subalignments->end(), subalignment_comparator_t());
}

void BamHelper::write_alignment_record(BamTools::BamWriter& bam_writer, const vector<BamTools::BamAlignment*>& alignments, const vector<vector<BamHelper::subalignment_t> >& subalignments, int aln_idx, int cigar_idx, const BamTools::BamAlignment* mate_aln, double probability, bool is_primary, bool retain_alternative_cigars, bool reduce_cigar, bool readgroups_from_names) {
	assert(aln_idx < alignments.size());
	BamTools::BamAlignment aln(*(alignments[aln_idx]));
	assert(aln_idx < (int)subalignments.size());
	if (aln.IsMapped()) {
		assert(cigar_idx < (int)subalignments[aln_idx].size());
		const BamHelper::subalignment_t& subalignment = subalignments[aln_idx][cigar_idx];
		aln.CigarData.assign(subalignment.cigar.begin(), subalignment.cigar.end());
		uint32_t x = 0;
		float f = 0.0;
		string s = "";
		if (aln.GetTag("AS", x)) aln.RemoveTag("AS");
		if (aln.GetTag("YM", x)) aln.RemoveTag("YM");
		if (aln.GetTag("YA", s)) aln.RemoveTag("YA");
		if (aln.GetTag("YP", f)) aln.RemoveTag("YP");
		if (!aln.AddTag("AS", "I", (uint32_t)subalignment.phred_score)) assert(false);
		if (!aln.AddTag("YM", "I", (uint32_t)subalignment.mismatch_phred_score)) assert(false);
		if (!aln.AddTag("YP", "f", (float)probability)) assert(false);
		if (is_primary) {
			if (!aln.AddTag("X1", "I", (uint32_t)alignments.size())) assert(false);
		}
		if (retain_alternative_cigars && (subalignments[aln_idx].size()>1)) {
			ostringstream oss;
			bool first = true;
			for (size_t i = 0; i < subalignments[aln_idx].size(); ++i) {
				if ((int)i == cigar_idx) continue;
				const BamHelper::subalignment_t& subalignment = subalignments[aln_idx][i];
				if (!first) oss << ';';
				oss << subalignment.phred_score << ',' << subalignment.mismatch_phred_score << ',' << subalignment.cigar;
				first = false;
			}
			if (!aln.AddTag("YA", "Z", oss.str())) assert(false);
		}
		aln.SetIsPrimaryAlignment(is_primary);
		aln.SetIsProperPair(mate_aln != 0);
		if (mate_aln != 0) {
			aln.MateRefID = mate_aln->RefID;
			aln.MatePosition = mate_aln->Position;
			aln.SetIsMateMapped(true);
			aln.SetIsMateReverseStrand(mate_aln->IsReverseStrand());
		} else {
			aln.MateRefID = -1;
			aln.MatePosition = -1;
			aln.SetIsMateMapped(false);
			aln.SetIsMateReverseStrand(false);
		}
		if (probability < 1.0) {
			aln.MapQuality = (uint16_t)round(min(-10.0*log1p(-probability)/log(10.0),100.0));
		} else {
			aln.MapQuality = (uint16_t)100;
		}
		if (reduce_cigar) {
			BamHelper::cigarReduceMismatches(&aln.CigarData);
		}
	} else {
		if (mate_aln != 0) {
			aln.RefID = mate_aln->RefID;
			aln.Position = mate_aln->Position;
			aln.MateRefID = mate_aln->RefID;
			aln.MatePosition = mate_aln->Position;
			aln.SetIsMateMapped(true);
			aln.SetIsMateReverseStrand(mate_aln->IsReverseStrand());
		} else {
			aln.RefID = -1;
			aln.Position = -1;
			aln.MateRefID = -1;
			aln.MatePosition = -1;
			aln.SetIsMateMapped(false);
			aln.SetIsMateReverseStrand(false);
		}
	}
	if (readgroups_from_names) {
		size_t n = aln.Name.find("_");
		if (n == string::npos) {
			ostringstream oss;
			oss << "Error: No underscore present in read name \"" << aln.Name << "\"." << endl;
			throw std::runtime_error(oss.str());
		}
		string readgroup = aln.Name.substr(0, n);
		string s;
		if (aln.GetTag("RG", s)) aln.RemoveTag("RG");
		if (!aln.AddTag("RG", "Z", readgroup)) assert(false);
		aln.Name = aln.Name.substr(n+1, aln.Name.size()-n-1);
	}
	bam_writer.SaveAlignment(aln);
}

BamHelper::pairing_t BamHelper::computeBestPairing(const BamHelper::read_t& read, const HistogramBasedDistribution& internal_segment_size_dist, VariationIndex* variation_set) {
	BamHelper::pairing_t best_pairing(-1,-1,-1,-1);
	double best_prob = 0.0;
	for (size_t i=0; i<read.alignments1.size(); ++i) {
		const BamTools::BamAlignment& aln1 = *(read.alignments1[i]);
		if (!aln1.IsMapped()) continue;
		double p1 = pow(10,-((double)read.subalignments1[i][0].phred_score)/10.0);
		for (size_t j=0; j<read.alignments2.size(); ++j) {
			const BamTools::BamAlignment& aln2 = *(read.alignments2[j]);
			if (!aln2.IsMapped()) continue;
			if (aln1.RefID != aln2.RefID) continue;
			if (aln1.IsReverseStrand() == aln2.IsReverseStrand()) continue;
			if (aln1.Position == aln2.Position) continue;
			double p2 = pow(10,-((double)read.subalignments2[j][0].phred_score)/10.0);
			int left, right;
			if (aln1.Position < aln2.Position) {
				if (aln1.IsReverseStrand()) continue;
				left = aln1.GetEndPosition();
				right = aln2.Position - 1;
			} else {
				if (aln2.IsReverseStrand()) continue;
				left = aln2.GetEndPosition();
				right = aln1.Position - 1;
			}
			int insert_length = right - left + 1;
			double p = p1 * internal_segment_size_dist.probability(insert_length) * p2;
			// if set of known (or putative) variations is given,
			// look for variations that lie inside internal segment and 
			// use best probability
			if (variation_set != 0) {
				unique_ptr<vector<size_t> > variations = variation_set->containedIn(aln1.RefID, left, right+1);
				if (variations.get() != 0) {
					for (size_t k=0; k<variations->size(); ++k) {
						const Variation& v = variation_set->get(variations->at(k));
						double p_new = p1 * internal_segment_size_dist.probability(insert_length - v.getLengthDifference()) * p2;
						if (p_new > p) {
							p = p_new;
						}
					}
				}
			}
			if (p > best_prob) {
				best_prob = p;
				best_pairing = BamHelper::pairing_t(i,j);
			}
		}
	}
	return best_pairing;
}

unique_ptr<vector<double> > BamHelper::compute_alignment_distribution(const vector<vector<BamHelper::subalignment_t> >& subalignments, int* best) {
	unique_ptr<vector<double> > result(new vector<double>());
	if (best != 0) *best = -1;
	// no alignments at all
	if (subalignments.size() == 0) return result;
	// unmapped read
	if ((subalignments.size() == 1) && (subalignments[0].size() == 0)) return result;
	double prob_sum = 0.0;
	double best_prob = 0.0;
	for (size_t i=0; i<subalignments.size(); ++i) {
		double p = pow(10,-((double)subalignments[i][0].phred_score)/10.0);
		result->push_back(p);
		prob_sum += p;
		if ((best != 0) && (p > best_prob)) {
			best_prob = p;
			*best = i;
		}
	}
	assert(prob_sum > 0.0);
	for (size_t i=0; i<subalignments.size(); ++i) {
		result->at(i) /= prob_sum;
	}
	return result;
}

void BamHelper::writeAlignments(BamTools::BamWriter& bam_writer, const BamHelper::read_t& read, const BamHelper::pairing_t& best_pairing, bool retain_suboptimal, bool retain_alternative_cigars, bool reduce_cigar, bool readgroups_from_names) {
	int best1 = -1;
	int best2 = -1;
	unique_ptr<vector<double> > dist1 = compute_alignment_distribution(read.subalignments1, &best1);
	unique_ptr<vector<double> > dist2 = compute_alignment_distribution(read.subalignments2, &best2);
	// Does a (best) pairing exist?
	if (best_pairing.aln_idx1 != -1) {
		write_alignment_record(bam_writer, read.alignments1, read.subalignments1, best_pairing.aln_idx1, best_pairing.cigar_idx1, read.alignments2[best_pairing.aln_idx2], dist1->at(best_pairing.aln_idx1), true, retain_alternative_cigars, reduce_cigar, readgroups_from_names);
		if (retain_suboptimal) {
			for (size_t i = 0; i < read.alignments1.size(); ++i) {
				if (i == best_pairing.aln_idx1) continue;
				write_alignment_record(bam_writer, read.alignments1, read.subalignments1, i, 0, 0, dist1->at(i), false, retain_alternative_cigars, reduce_cigar, readgroups_from_names);
			}
		}
		write_alignment_record(bam_writer, read.alignments2, read.subalignments2, best_pairing.aln_idx2, best_pairing.cigar_idx2, read.alignments1[best_pairing.aln_idx1], dist2->at(best_pairing.aln_idx2), true, retain_alternative_cigars, reduce_cigar, readgroups_from_names);
		if (retain_suboptimal) {
			for (size_t i = 0; i < read.alignments2.size(); ++i) {
				if (i == best_pairing.aln_idx2) continue;
				write_alignment_record(bam_writer, read.alignments2, read.subalignments2, i, 0, 0, dist2->at(i), false, retain_alternative_cigars, reduce_cigar, readgroups_from_names);
			}
		}
	} else {
		for (size_t i = 0; i < read.alignments1.size(); ++i) {
			if ((best1 != -1) && (best1 != i) && (!retain_suboptimal)) continue;
			bool unmapped = !read.alignments1[i]->IsMapped();
			double p = unmapped?0.0:dist1->at(i);
			write_alignment_record(bam_writer, read.alignments1, read.subalignments1, i, 0, (unmapped && (best2 != -1))?read.alignments2[best2]:0, p, i == best1, retain_alternative_cigars, reduce_cigar, readgroups_from_names);
		}
		for (size_t i = 0; i < read.alignments2.size(); ++i) {
			if ((best2 != -1) && (best2 != i) && (!retain_suboptimal)) continue;
			bool unmapped = !read.alignments2[i]->IsMapped();
			double p = unmapped?0.0:dist2->at(i);
			write_alignment_record(bam_writer, read.alignments2, read.subalignments2, i, 0, (unmapped && (best1 != -1))?read.alignments1[best1]:0, p, i == best2, retain_alternative_cigars, reduce_cigar, readgroups_from_names);
		}
	}
}

void BamHelper::cigarReduceMismatches(vector<BamTools::CigarOp>* cigar) {
	assert(cigar != 0);
	vector<BamTools::CigarOp> result;
	size_t j = 0;
	for (size_t i=0; i<cigar->size(); ++i) {
		BamTools::CigarOp& c = cigar->at(i);
		char type = c.Type;
		switch (type) {
			case 'X':
			case '=':
				type = 'M';
				break;
			default:
				break;
		}
		if ((j>0) && (result[j-1].Type == type)) {
			result[j-1].Length += c.Length;
		} else {
			result.push_back(BamTools::CigarOp(type,c.Length));
			j += 1;
		}
	}
	cigar->assign(result.begin(), result.end());
}

void BamHelper::expandXA(const BamTools::BamReader& bam_reader, const BamTools::BamAlignment& aln, std::vector< BamTools::BamAlignment* >* target, boost::unordered_set<BamHelper::alignment_coordinate_t>* coordinates_set, long long* skipped_counter) {
	assert(target != 0);
	string xa = "";
	if (aln.GetTag("XA", xa)) {
		ShortDnaSequence seq(aln.QueryBases, aln.Qualities);
		typedef boost::tokenizer<boost::char_separator<char> > tokenizer_t;
		boost::char_separator<char> separator(",;");
		tokenizer_t tokenizer(xa, separator);
		tokenizer_t::const_iterator token_it = tokenizer.begin();
		while (token_it != tokenizer.end()) {
			unique_ptr<BamTools::BamAlignment> new_aln(new BamTools::BamAlignment());
			new_aln->RefID = bam_reader.GetReferenceID(*token_it);
			++token_it;
			assert(token_it != tokenizer.end());
			int32_t pos = boost::lexical_cast<int>(*token_it);
			new_aln->Position = abs(pos) - 1;
			++token_it;
			assert(token_it != tokenizer.end());
			parseCigar(*token_it, &new_aln->CigarData);
			++token_it;
			assert(token_it != tokenizer.end());
			if (!new_aln->AddTag("NM","i",boost::lexical_cast<int32_t>(*token_it))) assert(false);
			++token_it;
			alignment_coordinate_t coords(*new_aln);
			if ((coordinates_set != 0) && (coordinates_set->find(coords) != coordinates_set->end())) {
				if (skipped_counter != 0) {
					*skipped_counter += 1;
				}
			} else {
				new_aln->Name = aln.Name;
				new_aln->SetIsPaired(aln.IsPaired());
				new_aln->SetIsMapped(true);
				new_aln->SetIsReverseStrand(pos < 0);
				new_aln->SetIsFirstMate(aln.IsFirstMate());
				new_aln->SetIsSecondMate(aln.IsSecondMate());
				new_aln->SetIsPrimaryAlignment(false);
				if (new_aln->IsReverseStrand() == aln.IsReverseStrand()) {
					new_aln->QueryBases = seq.toString();
					new_aln->Qualities = seq.qualityString();
				} else {
					ShortDnaSequence rev = seq.reverseComplement();
					new_aln->QueryBases = rev.toString();
					new_aln->Qualities = rev.qualityString();
				}
				target->push_back(new_aln.release());
				if (coordinates_set != 0) {
					coordinates_set->insert(coords);
				}
			}
		}
	}
}

unique_ptr<vector<Variation> > BamHelper::variationsFromAlignment(const BamTools::RefVector& bam_ref_data, const BamTools::BamAlignment& aln) {
	unique_ptr<vector<Variation> > result(new vector<Variation>());
	if (!aln.IsMapped()) return result;
	int pos = aln.Position;
	for (size_t i = 0; i < aln.CigarData.size(); ++i) {
		const BamTools::CigarOp& op = aln.CigarData[i];
		switch (op.Type) {
			case 'M':
			case 'X':
			case '=':
				pos += op.Length;
				break;
			case 'I':
				result->push_back(Variation(bam_ref_data[aln.RefID].RefName, pos, op.Length, 0.0, Variation::INSERTION));
				break;
			case 'D':
				result->push_back(Variation(bam_ref_data[aln.RefID].RefName, pos, pos+op.Length, 0.0, Variation::DELETION));
				pos += op.Length;
				break;
			default:
				break;
		}
	}
	return result;
}

void BamHelper::readRegion(BamTools::BamReader& bam_reader, int chromosome_id, int start, int end, std::vector<aln_pair_t>* target) {
	bam_reader.SetRegion(chromosome_id, start, chromosome_id, end);
	unique_ptr<BamTools::BamAlignment> aln(new BamTools::BamAlignment());
	typedef boost::unordered_map<string, aln_pair_t> pair_map_t;
	pair_map_t pair_map;
	while (bam_reader.GetNextAlignment(*aln)) {
		if (!aln->IsMapped()) continue;
		if (!aln->IsProperPair()) continue;
		if (!aln->IsPrimaryAlignment()) continue;
		const string& name = aln->Name;
		pair_map_t::iterator it = pair_map.find(name);
		if (it == pair_map.end()) {
			aln_pair_t p;
			if (aln->IsFirstMate()) p.first = aln.release();
			else p.second = aln.release();
			pair_map[name] = p;
		} else {
			if (aln->IsFirstMate()) {
				if (it->second.first == 0) it->second.first = aln.release();
			} else {
				if (it->second.second == 0) it->second.second = aln.release();
			}
		}
		aln = unique_ptr<BamTools::BamAlignment>(new BamTools::BamAlignment());
	}
	assert(target != 0);
	pair_map_t::const_iterator it = pair_map.begin();
	for (; it != pair_map.end(); ++it) {
		if ((it->second.first != 0) && (it->second.second != 0)) {
			target->push_back(it->second);
		} else {
			if (it->second.first != 0) delete it->second.first;
			if (it->second.second != 0) delete it->second.second;
		}
	}
}

std::size_t hash_value(const BamHelper::alignment_coordinate_t& c) {
	size_t h = 0;
	boost::hash_combine(h, c.ref_id);
	boost::hash_combine(h, c.start);
	boost::hash_combine(h, c.end);
	return h;
}

ostream& operator<<(ostream& os, const BamHelper::alignment_coordinate_t& c) {
	os << c.ref_id << ":" << c.start << "-" << c.end;
	return os;
}

ostream& operator<<(ostream& os, const std::vector<BamTools::CigarOp>& cigar) {
	for (size_t i = 0; i<cigar.size(); ++i) {
		os << cigar[i].Length << cigar[i].Type;
	}
	return os;
}
