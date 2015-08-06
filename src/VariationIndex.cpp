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

#include <boost/unordered_set.hpp>

#include "VariationIndex.h"

using namespace std;

VariationIndex::VariationIndex(const vector<Variation>& variations, bool overlap_queries, const BamTools::RefVector* ref_data) : variations(variations), overlap_queries_allowed(overlap_queries) {
	this->step_size = 5000;
	// PHASE 0: handle chromosome IDs
	query_by_id_allowed = ref_data != 0;
	if (ref_data != 0) {
		for (size_t i=0; i<ref_data->size(); ++i) {
			chromosome_id_map[ref_data->at(i).RefName] = i;
		}
	}
	// PHASE 1: Add add all events to event lists
	for (size_t i=0; i<variations.size(); ++i) {
		const Variation& variation = variations[i];
		int chromosome_id = get_id(variation.getChromosome(), true);
		if (chromosomes.find(chromosome_id) == chromosomes.end()) {
			chromosomes[chromosome_id] = chromosome_record_t();
		}
		chromosome_record_t& chromosome_record = chromosomes[chromosome_id];
		switch (variation.getType()) {
		case Variation::DELETION:
			chromosome_record.events.push_back(event_t(variation.getCoordinate1(), i, DELETION_START));
			chromosome_record.events.push_back(event_t(variation.getCoordinate2(), i, DELETION_END));
			break;
		case Variation::INSERTION:
			chromosome_record.events.push_back(event_t(variation.getCoordinate1(), i, INSERTION));
			break;
		default:
			assert(false);
		}
	}
	// PHASE 2: Sort event lists and create "checkpoints" of active events
	typedef boost::unordered_map<int,chromosome_record_t>::iterator map_it_t;
	for (map_it_t it = chromosomes.begin(); it != chromosomes.end(); ++it) {
		// cout << "================= " << it->first << " =================" << endl;
		vector<event_t>& events = it->second.events;
		vector<checkpoint_t>& checkpoint_list = it->second.checkpoint_list;
		sort(events.begin(), events.end(), event_compare_t());
		typedef vector<event_t>::const_iterator event_it_t;
		size_t i = 0;
		size_t j = 0;
		for (event_it_t event_it=events.begin(); event_it!=events.end(); ++event_it) {
			for (;i*step_size<event_it->position; ++i) {
				checkpoint_list.push_back(checkpoint_t(j));
			}
			j += 1;
		}
		// if overlap queries are requested, populate field "active" in each checkpoint record.
		if (overlap_queries) {
			boost::unordered_set<size_t> active_variations;
			event_it_t event_it=events.begin();
			for (i=0; i<checkpoint_list.size(); ++i) {
				while (event_it!=events.end()) {
					if (event_it->position <= i*step_size) {
						switch (event_it->type) {
						case DELETION_START:
							active_variations.insert(event_it->variation_index);
							break;
						case DELETION_END:
							active_variations.erase(event_it->variation_index);
							break;
						default:
							break;
						}
						++event_it;
					} else {
						break;
					}
				}
				checkpoint_list[i].active.assign(active_variations.begin(), active_variations.end());
			}
		}
	}
}

VariationIndex::~VariationIndex() {
}

const Variation& VariationIndex::get(int variation_id) {
	assert(variation_id >= 0);
	assert(variation_id < variations.size());
	return variations[variation_id];
}

int VariationIndex::get_id(const std::string& chromosome, bool create_new) {
	chromosome_id_map_t::const_iterator it = chromosome_id_map.find(chromosome);
	if (it == chromosome_id_map.end()) {
		if (create_new) {
			int new_id = chromosome_id_map.size();
			chromosome_id_map[chromosome] = new_id;
			return new_id;
		} else {
			return -1;
		}
	} else {
		return it->second;
	}
}

unique_ptr<vector<size_t> > VariationIndex::containedIn(const string& chromosome, size_t start, size_t end) {
	return _containedIn(get_id(chromosome, false), start, end);
}

unique_ptr<vector<size_t> > VariationIndex::containedIn(int chromosome_id, size_t start, size_t end) {
	assert(query_by_id_allowed);
	return _containedIn(chromosome_id, start, end);
}

unique_ptr<vector<size_t> > VariationIndex::_containedIn(int chromosome_id, size_t start, size_t end) {
	unique_ptr<vector<size_t> > result(nullptr);
	if (chromosomes.find(chromosome_id) == chromosomes.end()) {
		return result;
	}
	const chromosome_record_t& chromosome_record = chromosomes[chromosome_id];
	size_t k = (start-1) / step_size;
	if (chromosome_record.checkpoint_list.size() <= k) {
		return result;
	}
	const vector<event_t>& event_list = chromosome_record.events;
	boost::unordered_set<size_t> active_variations;
	for (size_t i=chromosome_record.checkpoint_list[k].index; i<event_list.size(); ++i) {
		const event_t& event = event_list[i];
		if (event.position < start) continue;
		if (event.position > end) break;
		switch (event.type) {
		case INSERTION:
			if (result.get() == 0) {
				result = unique_ptr<vector<size_t> >(new vector<size_t>());
			}
			result->push_back(event.variation_index);
			break;
		case DELETION_START:
			active_variations.insert(event.variation_index);
			break;
		case DELETION_END:
			if (active_variations.find(event.variation_index) != active_variations.end()) {
				if (result.get() == 0) {
					result = unique_ptr<vector<size_t> >(new vector<size_t>());
				}
				result->push_back(event.variation_index);
			}
			break;
		default:
			assert(false);
			break;
		}
	}
	return result;
}

unique_ptr<vector<size_t> > VariationIndex::overlapping(const string& chromosome, size_t start, size_t end) {
	return _overlapping(get_id(chromosome, false), start, end);
}

unique_ptr<vector<size_t> > VariationIndex::overlapping(int chromosome_id, size_t start, size_t end) {
	assert(query_by_id_allowed);
	return _overlapping(chromosome_id, start, end);
}

unique_ptr<vector<size_t> > VariationIndex::_overlapping(int chromosome_id, size_t start, size_t end) {
	if (!overlap_queries_allowed) throw runtime_error("VariationIndex: data structures for overlap queries not initialized.");
	if (start > end) throw runtime_error("VariationIndex::overlapping(): start > end.");;
	if (chromosomes.find(chromosome_id) == chromosomes.end()) {
		return unique_ptr<vector<size_t> >(nullptr);
	}
	const chromosome_record_t& chromosome_record = chromosomes[chromosome_id];
	size_t k = start / step_size;
	if (chromosome_record.checkpoint_list.size() <= k) {
		return unique_ptr<vector<size_t> >(nullptr);
	}
	// step 1: determine set of active events at position "start"
	boost::unordered_set<size_t> result_set(chromosome_record.checkpoint_list[k].active.begin(), chromosome_record.checkpoint_list[k].active.end());
	int i = chromosome_record.checkpoint_list[k].index;
	for (;i<chromosome_record.events.size(); ++i) {
		const event_t& e = chromosome_record.events[i];
		if (e.position > start) break;
		switch (e.type) {
		case DELETION_START:
			result_set.insert(e.variation_index);
			break;
		case DELETION_END:
			result_set.erase(e.variation_index);
			break;
		default:
			break;
		}
	}
	// step 2: determine result by running through all events overlapping the given interval
	for (;i<chromosome_record.events.size(); ++i) {
		const event_t& e = chromosome_record.events[i];
		if (e.position >= end) break;
		switch (e.type) {
		case INSERTION:
		case DELETION_START:
			result_set.insert(e.variation_index);
			break;
		default:
			break;
		}
	}
	if (result_set.size() == 0) return unique_ptr<vector<size_t> >(nullptr);
	else return unique_ptr<vector<size_t> >(new vector<size_t>(result_set.begin(), result_set.end()));
}

std::unique_ptr<std::vector<size_t> > VariationIndex::getDeletions(const std::string& chromosome, size_t start, size_t end) {
	return _getDeletions(get_id(chromosome, false), start, end);
}

std::unique_ptr<std::vector<size_t> > VariationIndex::getDeletions(int chromosome_id, size_t start, size_t end) {
	assert(query_by_id_allowed);
	return _getDeletions(chromosome_id, start, end);
}

std::unique_ptr<std::vector<size_t> > VariationIndex::_getDeletions(int chromosome_id, size_t start, size_t end) {
	if (chromosomes.find(chromosome_id) == chromosomes.end()) {
		return unique_ptr<vector<size_t> >(nullptr);
	}
	const chromosome_record_t& chromosome_record = chromosomes[chromosome_id];
	size_t k = (start-1) / step_size;
	if (chromosome_record.checkpoint_list.size() <= k) {
		return unique_ptr<vector<size_t> >(nullptr);
	}
	const vector<event_t>& event_list = chromosome_record.events;
	boost::unordered_set<size_t> candidates;
	unique_ptr<vector<size_t> > result(nullptr);
	for (size_t i=chromosome_record.checkpoint_list[k].index; i<event_list.size(); ++i) {
		const event_t& event = event_list[i];
		if (event.position > end) break;
		if ((event.position == start) && (event.type == DELETION_START)) {
			candidates.insert(event.variation_index);
		}
		if ((event.position == end) && (event.type == DELETION_END)) {
			if (candidates.find(event.variation_index) != candidates.end()) {
				if (result.get() == 0) {
					result = unique_ptr<vector<size_t> >(new vector<size_t>());
				}
				result->push_back(event.variation_index);
			}
		}
	}
	return result;
}

std::unique_ptr<std::vector<size_t> > VariationIndex::getInsertions(const std::string& chromosome, size_t position) {
	return _getInsertions(get_id(chromosome, false), position);
}

std::unique_ptr<std::vector<size_t> > VariationIndex::getInsertions(int chromosome_id, size_t position) {
	assert(query_by_id_allowed);
	return _getInsertions(chromosome_id, position);
}

std::unique_ptr<std::vector<size_t> > VariationIndex::_getInsertions(int chromosome_id, size_t position) {
	if (chromosomes.find(chromosome_id) == chromosomes.end()) {
		return unique_ptr<vector<size_t> >(nullptr);
	}
	const chromosome_record_t& chromosome_record = chromosomes[chromosome_id];
	size_t k = (position-1) / step_size;
	if (chromosome_record.checkpoint_list.size() <= k) {
		return unique_ptr<vector<size_t> >(nullptr);
	}
	const vector<event_t>& event_list = chromosome_record.events;
	unique_ptr<vector<size_t> > result(nullptr);
	for (size_t i=chromosome_record.checkpoint_list[k].index; i<event_list.size(); ++i) {
		const event_t& event = event_list[i];
		if (event.position > position) break;
		if ((event.position == position) && (event.type == INSERTION)) {
			if (result.get() == 0) {
				result = unique_ptr<vector<size_t> >(new vector<size_t>());
			}
			result->push_back(event.variation_index);
		}
	}
	return result;
}
