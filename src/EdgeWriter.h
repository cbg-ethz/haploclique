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

#ifndef EDGEWRITER_H_
#define EDGEWRITER_H_

#include <iostream>
#include <deque>
#include <boost/unordered_set.hpp>
#include <map>

#include "AlignmentRecord.h"
#include "Types.h"

class EdgeWriter {
	std::ostream& os;
	typedef boost::unordered_set<alignment_id_t> node_set_t;
	std::map<alignment_id_t, node_set_t> nodes;
	bool finished;
	void printNode(alignment_id_t id);
    std::map<alignment_id_t, std::string> vertex_to_read_names;
public:
	EdgeWriter(std::ostream& os);
	virtual ~EdgeWriter();

	/** Sets an edge between the two given nodes. */
	virtual void addEdge(const AlignmentRecord& node1, const AlignmentRecord& node2);

	/** Calling this method indicates that the given node can be written to output
	 *  as no more edges will be added. */
	virtual void setNodeCompleted(const AlignmentRecord& node);

	/** Output all nodes and edges that are still unfinished. */
	virtual void finish();

    virtual void initialize();
};

#endif /* EDGEWRITER_H_ */
