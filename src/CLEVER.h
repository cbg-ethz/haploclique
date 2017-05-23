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

#ifndef CLEVER_H_
#define CLEVER_H_

#include <set>
#include <list>
#include <boost/unordered_map.hpp>
 #include <boost/dynamic_bitset.hpp>

#include "Clique.h"
#include "CliqueFinder.h"
#include "LogWriter.h"

/** Implementation of the Maximal Clique Enumeration algorithm of CLEVER */
class CLEVER : public CliqueFinder {
private:
    size_t capacity;
    AlignmentRecord **alignments;
    LogWriter* lw;
    unsigned int max_cliques;
    unsigned int limit_clique_size;
    unsigned int clique_counter;
    std::vector<unsigned int> read_in_cliques;
    bool filter_singletons;
    typedef std::pair<unsigned int,size_t> length_and_index_t;
    std::set<length_and_index_t> alignments_by_length;
    /** reorganizes the AlignmentRecord data structure to save memory. */
    void reorganize_storage();
public:
    CLEVER(const EdgeCalculator& edge_calculator, CliqueCollector& clique_collector, LogWriter* lw, unsigned int max_cliques, unsigned int limit_clique_size, unsigned int number_of_reads, bool filter_singletons);
    virtual ~CLEVER();
    const AlignmentRecord & getAlignmentByIndex(size_t index) const {
        assert(index<alignment_count);
    	return *(alignments[index]);
    }
    /** returns index of the element with the highest Priority. */
    unsigned int getPriorityRead();
    /** returns the largest clique that contains the read with id equal to "readref". */
    Clique* getLargestClique(unsigned int readref);
    /** adds qualified cliques to clique_collector and filter the rest. */
    void finish();
    void initialize();
    /** constructs the adjacancy bitset for the new alignment and performs clique operations (i.e. create new clique, split a clique) based on the adjacancy bitset of the new Alignment Record. */
    void addAlignment(std::unique_ptr<AlignmentRecord>& ap, int& edgecounter);
};

#endif /* CLEVER_H_ */
