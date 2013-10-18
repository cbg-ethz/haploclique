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

#ifndef CLIQUEFINDER_H_
#define CLIQUEFINDER_H_

#include <set>
#include <list>
#include <boost/unordered_map.hpp>

#include "Clique.h"
#include "CliqueCollector.h"
#include "AlignmentRecord.h"
#include "EdgeCalculator.h"
#include "CoverageMonitor.h"
#include "EdgeWriter.h"
#include "ReadGroups.h"

/** A class to store a set of alignments. */
class CliqueFinder {
private:
    typedef std::pair<unsigned int,size_t> length_and_index_t;
    std::set<length_and_index_t> alignments_by_length;
    AlignmentRecord **alignments;
    size_t alignment_count;
    size_t capacity;
    typedef std::list<Clique*> clique_list_t;
    clique_list_t *cliques;
    CliqueCollector & clique_collector;
    EdgeWriter *edge_writer;
    const EdgeCalculator & edge_calculator;
    const EdgeCalculator *second_edge_calculator;
    CoverageMonitor coverage_monitor;
    alignment_id_t next_id;
    void reorganize_storage();
public:
    CliqueFinder(const EdgeCalculator& edge_calculator, CliqueCollector& clique_collector, const ReadGroups* read_groups);
    virtual ~CliqueFinder();
    void addAlignment(std::auto_ptr<AlignmentRecord> ap);
    void finish();
    CoverageMonitor & getCoverageMonitor() { return coverage_monitor; }

    const AlignmentRecord & getAlignmentByIndex(size_t index) const;

    void setEdgeWriter(EdgeWriter& edge_writer) { this->edge_writer = &edge_writer; }
    
    /** Add a second edge calculator: edges will only be drawn when both edge calculators agree
	 *  that the edge is present.
	 */
    void setSecondEdgeCalculator(const EdgeCalculator* ec) { this->second_edge_calculator = ec; }
};

#endif /* CLIQUEFINDER_H_ */
