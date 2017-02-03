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

#ifndef VARIATIONINDEX_H_
#define VARIATIONINDEX_H_

#include <vector>
#include <boost/unordered_map.hpp>

#include "Variation.h"
#include <api/BamAux.h>

/** Index structure over a set of variations able to quickly answer
 *  "contained-in queries". That is, given an interval, it can (quickly)
 *  return all variations contained in that interval. */
class VariationIndex {
private:
	typedef enum { INSERTION, DELETION_START, DELETION_END } event_type_t;
	typedef struct event_t {
		size_t position;
		size_t variation_index;
		event_type_t type;
		event_t(size_t position, size_t variation_index, event_type_t type) : position(position), variation_index(variation_index), type(type) {}
	} event_t;
	typedef struct {
		bool operator() (const event_t& e1, const event_t& e2) { return (e1.position<e2.position);}
	} event_compare_t;
	typedef struct checkpoint_t {
		// index of the first event after that checkpoint
		size_t index;
		// indices of variations active at the checkpoint position,
		// i.e. all event with positions <= the checkpoint position 
		// have been processed
		std::vector<size_t> active;
		checkpoint_t() : index(0) {}
		checkpoint_t(size_t index) : index(index) {}
	} checkpoint_t;
	// a record of this type is kept for each chromosome
	typedef struct {
		std::vector<event_t> events;
		// in distances of step_size, this vector stores the index of the
		// next event in the event list after the given position, along
		// with all variations active at the checkpoint position;
		// e.g. active_events[n].index gives the index of the first event
		// after position n*step_size
		std::vector<checkpoint_t> checkpoint_list;
	} chromosome_record_t;
	const std::vector<Variation>& variations;
	// Mapping chromosome IDs to sorted lists of events
	boost::unordered_map<int, chromosome_record_t> chromosomes;
	// Mapping chromosome names to IDs
	typedef boost::unordered_map<std::string, int> chromosome_id_map_t;
	chromosome_id_map_t chromosome_id_map;
	bool query_by_id_allowed;
	size_t step_size;
	bool overlap_queries_allowed;
	int get_id(const std::string& chromosome, bool create_new);
	
	std::unique_ptr<std::vector<size_t> > _containedIn(int chromosome_id, size_t start, size_t end);
	std::unique_ptr<std::vector<size_t> > _overlapping(int chromosome_id, size_t start, size_t end);
	std::unique_ptr<std::vector<size_t> > _getDeletions(int chromosome_id, size_t start, size_t end);
	std::unique_ptr<std::vector<size_t> > _getInsertions(int chromosome_id, size_t position);
public:
	/** Indexes the given variations and keeps a reference to them. That is, the given
	 *  reference must be valid throughout the lifetime of this object. 
	 * @param overlap_queries If true, additional space is used to allow for overlap queries
	 *                        as implemented in overlapping().
	 */
	VariationIndex(const std::vector<Variation>& variations, bool overlap_queries = false, const BamTools::RefVector* ref_data = 0);
	virtual ~VariationIndex();

	const Variation& get(int variation_id);

	/** Returns the set of indices (w.r.t. to the vector given at construction time) of
	 *  variations completely contained in the interval from start to end.
	 *  start/end are interpreted as positions between two characters, that is, a value n
	 *  refers to the position between character n and n-1. Returns a null pointer when there
	 *  are no variations in that region. */
	std::unique_ptr<std::vector<size_t> > containedIn(const std::string& chromosome, size_t start, size_t end);
	std::unique_ptr<std::vector<size_t> > containedIn(int chromosome_id, size_t start, size_t end);
	
	/** Returns the set of indices (w.r.t. to the vector given at construction time) of
	 *  variations that overlap with the given interval from start to end.
	 *  start/end are interpreted as positions between two characters, that is, a value n
	 *  refers to the position between character n and n-1. Returns a null pointer when there
	 *  are no variations in that region. */
	std::unique_ptr<std::vector<size_t> > overlapping(const std::string& chromosome, size_t start, size_t end);
	std::unique_ptr<std::vector<size_t> > overlapping(int chromosome_id, size_t start, size_t end);
	
	/** Returns deletions with exactly the given start and end point. */
	std::unique_ptr<std::vector<size_t> > getDeletions(const std::string& chromosome, size_t start, size_t end);
	std::unique_ptr<std::vector<size_t> > getDeletions(int chromosome_id, size_t start, size_t end);
	
	/** Returns deletions with exactly the position. */
	std::unique_ptr<std::vector<size_t> > getInsertions(const std::string& chromosome, size_t position);
	std::unique_ptr<std::vector<size_t> > getInsertions(int chromosome_id, size_t position);
};

#endif /* VARIATIONINDEX_H_ */
