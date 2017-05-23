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

#include <string.h>
#include <iostream>
#include <vector>
#include <boost/unordered_set.hpp>
#include <boost/dynamic_bitset.hpp>
#include <ctime>
#include <map>
#include "CLEVER.h"

using namespace std;
using namespace boost;

CLEVER::CLEVER(const EdgeCalculator& edge_calculator, CliqueCollector& clique_collector, LogWriter* lw, unsigned int max_cliques, unsigned int limit_clique_size, unsigned int number_of_reads, bool filter_singletons) : CliqueFinder(edge_calculator, clique_collector), lw(lw), max_cliques(max_cliques), limit_clique_size(limit_clique_size), read_in_cliques(number_of_reads), filter_singletons(filter_singletons) {
    capacity = alignment_set_t::bits_per_block;
    alignments = nullptr;
}

CLEVER::~CLEVER() {
	if (cliques!=nullptr) {
		finish();
	}
}

unsigned int CLEVER::getPriorityRead() {
    int min_index = std::min_element(this->read_in_cliques.begin(),this->read_in_cliques.end()) - this->read_in_cliques.begin();
    return min_index;
}

Clique* CLEVER::getLargestClique(unsigned int readref){
    Clique* res_clique = nullptr;
    auto max = std::numeric_limits<std::size_t>::min();
    
    clique_list_t::iterator clique_it = cliques->begin();
    for (;clique_it!=cliques->end(); ++clique_it) {
        auto read_names = (*clique_it)->getCliqueReadNamesSet();
        // first case: clique's set contains read
        if(read_names.find(readref) != read_names.end() && read_names.size() > max){
            max = read_names.size();
            res_clique = *clique_it;
        }
    }
    return res_clique;
}

void CLEVER::initialize() {
    cliques = new clique_list_t();
    alignments = new AlignmentRecord*[capacity];
    capacity = alignment_set_t::bits_per_block;
  	alignment_count = 0;
    next_id = 0;
    alignments_by_length.clear();
    clique_counter = 0;
    converged = true;
    initialized = true;
}

void CLEVER::finish() {
    // entries of read_in_cliques is set to zero at the beginning of each iteration
    std::fill(read_in_cliques.begin(), read_in_cliques.end(), 0);
    if(max_cliques == 0){
        clique_list_t::iterator clique_it = cliques->begin();
        for (;clique_it!=cliques->end(); ++clique_it) {
            Clique* clique = *clique_it;
            clique_collector.add(unique_ptr<Clique>(clique));
        }
    } else {
            while(clique_counter < max_cliques && !cliques->empty()){
                    // get read with largest priority
                    unsigned int readref = getPriorityRead();
                    // get largest Clique to read
                    Clique* clique = getLargestClique(readref);
                    if(clique!=nullptr){
                        if(filter_singletons){
                            if (clique->getCliqueReadCount() > 1){
                                // cout up all entries of reads which are contained in clique
                                for (const auto &r : clique->getCliqueReadNamesSet()){
                                    read_in_cliques[r] += 1;
                                }
                                clique_collector.add(unique_ptr<Clique>(clique));
                                clique_counter++;
                                cliques->remove(clique);
                            // Clique is a singleton and is not added to final clique set
                            } else {
                                read_in_cliques[readref] = std::numeric_limits<int>::max();
                                cliques->remove(clique);
                                delete clique;
                            }
                        // singletons allowed
                        } else {
                            // cout up all entries of reads which are contained in clique
                            for (const auto &r : clique->getCliqueReadNamesSet()){
                                read_in_cliques[r] += 1;
                            }
                            clique_collector.add(unique_ptr<Clique>(clique));
                            clique_counter++;
                            cliques->remove(clique);
                        }
                    } else {
                        read_in_cliques[readref] = std::numeric_limits<int>::max();
                    }
            }
    }

    delete cliques;
    cliques = nullptr;

    for (size_t i=0; i<alignment_count; ++i) {
    	delete alignments[i];
    }
    delete [] alignments;
    alignments = nullptr;

    initialized = false;
}

void CLEVER::reorganize_storage() {
	alignment_set_t set_union(capacity);
    
	clique_list_t::iterator it = cliques->begin();
	for (;it!=cliques->end(); ++it) {
		set_union |= (*it)->getAlignmentSet();
	}
    
    size_t new_alignment_count = set_union.count();
	size_t new_capacity = new_alignment_count * 2;
	// round up to next multiple of block size
	int k = new_capacity % alignment_set_t::bits_per_block;
	if (k > 0) {
		new_capacity += alignment_set_t::bits_per_block - k;
	}
    
	AlignmentRecord** new_alignments = new AlignmentRecord*[new_capacity];
	memset(new_alignments,0,new_capacity*sizeof(AlignmentRecord*));
	// position i in the new tables is equivalent to position "translation_table[i]"
	// in the old table
	size_t* translation_table = new size_t[new_alignment_count];
	size_t j = 0;
	alignments_by_length.clear();
	for (size_t i=0; i<alignment_count; ++i) {
		// test whether alignment i is still in use
		if (set_union[i]) {
            translation_table[j] = i;
			new_alignments[j] = alignments[i];
			// TODO: Once edge criteria are fixed, we can try to (re-)gain some efficiency here...
 			// if (single_end) {
				alignments_by_length.insert(make_pair(0,j));
 			// } else {
 			// 	alignments_by_length.insert(make_pair(new_alignments[j]->getInsertLength(),j));
 			// }
			j += 1;
		} else {
			delete alignments[i];
		}
	}
	// translate bit sets in all active cliques
	it = cliques->begin();
	size_t leftmost_pos = (*it)->leftmostSegmentStart();
	for (;it!=cliques->end(); ++it) {
		(*it)->translate(translation_table, new_alignment_count, new_capacity);
		leftmost_pos = min(leftmost_pos, (*it)->leftmostSegmentStart());
	}

	delete [] translation_table;
	delete [] alignments;
	alignments = new_alignments;
	alignment_count = new_alignment_count;
	capacity = new_capacity;
}

void CLEVER::addAlignment(std::unique_ptr<AlignmentRecord>& alignment_autoptr, int& edgecounter) {
    assert(alignment_autoptr.get() != 0);
	assert(cliques!=0);
    assert(initialized);

	alignment_id_t id = next_id++;
	AlignmentRecord* alignment = alignment_autoptr.release();
	alignment->setID(id);

	// store new alignment
	if (alignment_count==capacity) {
		reorganize_storage();
	}
	size_t index = alignment_count++;
	alignments[index] = alignment;

	// TODO: Once edge criteria are fixed, we can try to (re-)gain some efficiency here...
 	// if (single_end) {
 	// } else {
 	// 	alignments_by_length.insert(make_pair(alignment->getInsertLength(), index));
 	// }
    
	// determine all edges from current alignment pair
	alignment_set_t adjacent(capacity);
	// determine range of insert lengths to be searched
	typedef set<length_and_index_t>::const_iterator iterator_t;
	iterator_t it = alignments_by_length.begin();
	iterator_t end = alignments_by_length.end();
    
	// iterate through all alignments
	for (; it!=end; ++it) {
		const AlignmentRecord* alignment2 = alignments[it->second];
        bool set_edge = edge_calculator.edgeBetween(*alignment, *alignment2);
        if (set_edge && (second_edge_calculator != nullptr)) {
			set_edge = second_edge_calculator->edgeBetween(*alignment, *alignment2);
		}
		if (set_edge) {
            edgecounter++;
			adjacent.set(it->second, true);
			if (lw != nullptr) {
                lw->reportEdge(alignment->getID(), alignment2->getID());
			}
        converged = false;
		}
	}

    alignments_by_length.insert(make_pair(0, index));
	// iterate over all active cliques. output those that lie left of current segment and
	// check intersection with current node for the rest
	clique_list_t::iterator clique_it = cliques->begin();
	// cliques that contain the newly added alignment and
	// therefore need to be checked for subset relations,
	// i.e. if one of these cliques is contained in another,
	// it must be discarded as it is not maximal.
	vector<Clique*> new_cliques;
	while (clique_it!=cliques->end()) {
		Clique* clique = *clique_it;
        if ((clique->rightmostSegmentEnd() < alignment->getIntervalStart() && max_cliques == 0) || ( limit_clique_size!=0 && clique->size() == limit_clique_size)){
            // TO DO !!!! should not be added to clique_collector unless criteria is fulfilled given max_cliques (Does it make things worse if cliques are still contained in cliques of clique_finder?
			clique_it = cliques->erase(clique_it);
			clique_collector.add(unique_ptr<Clique>(clique));
		} else {
			// is there an intersection between nodes adjacent to the new
			// alignment and the currently considered clique?
			unique_ptr<alignment_set_t> intersection = clique->intersect(adjacent);
			if (intersection->any()) {
				// is node adjacent to all nodes in the clique?
				if (intersection->count() == clique->size()) {
					clique->add(index);
					new_cliques.push_back(clique);
					clique_it = cliques->erase(clique_it);
				} else {
					Clique* split_off_clique = new Clique(*this, intersection);
					split_off_clique->add(index);
					new_cliques.push_back(split_off_clique);
					++clique_it;
				}
			} else {
				++clique_it;
			}
		}
	}
	// if current alignment has not been assigned to at least one
	// of the existing cliques, let it form its own singleton clique
	if (new_cliques.size() == 0) {
		new_cliques.push_back(new Clique(*this, index, capacity));
	}
	// check for subset relations and delete cliques that are subsets of others
	for (size_t i=0; i<new_cliques.size(); ++i) {
		if (new_cliques[i]==0) continue;
		Clique* clique_i = new_cliques[i];
		for (size_t j=i+1; j<new_cliques.size(); ++j) {
			if (new_cliques[j]==0) continue;
			Clique* clique_j = new_cliques[j];
			if (clique_i->size()<=clique_j->size()) {
				if (clique_j->contains(*clique_i)) {
					delete clique_i;
					new_cliques[i] = nullptr;
					break;
				}
			} else {
				if (clique_i->contains(*clique_j)) {
					delete clique_j;
					new_cliques[j] = nullptr;
					continue;
				}
			}
		}
	}
	for (size_t i=0; i<new_cliques.size(); ++i) {
		if (new_cliques[i]!=nullptr) {
			cliques->push_back(new_cliques[i]);
		}
	}
}
