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

#include "EdgeWriter.h"

using namespace std;

EdgeWriter::EdgeWriter(std::ostream& os) : os(os) {
	offset = 0;
	finished = false;
}

EdgeWriter::~EdgeWriter() {
	finish();
}


void EdgeWriter::addEdge(const AlignmentRecord& node1, const AlignmentRecord& node2) {
	assert(!finished);
	assert((node1.getID()>=offset) && (node2.getID()>=offset));
	alignment_id_t max_id = max(node1.getID(), node2.getID());
	while (offset + nodes.size() <= max_id) {
		nodes.push_back(new node_set_t());
	}
	nodes[node1.getID() - offset]->insert(node2.getID());
	nodes[node2.getID() - offset]->insert(node1.getID());
}

void EdgeWriter::printNode(alignment_id_t id) {
	assert(nodes[id-offset]!=0);
	vector<alignment_id_t> adj_nodes(nodes[id-offset]->begin(), nodes[id-offset]->end());
	sort(adj_nodes.begin(), adj_nodes.end());
	os << id << " ->";
	for (size_t i=0; i<adj_nodes.size(); ++i) {
		os << " " << adj_nodes[i];
	}
	os << endl;
}

void EdgeWriter::setNodeCompleted(const AlignmentRecord& node) {
	assert(!finished);
	const alignment_id_t& id = node.getID();
	assert(id>=offset);
	while (offset + nodes.size() <= id) {
		nodes.push_back(new node_set_t());
	}
	printNode(id);
	delete nodes[id-offset];
	nodes[id-offset] = 0;
	while (nodes[0] == 0) {
		nodes.pop_front();
		offset += 1;
	}
}

void EdgeWriter::finish() {
	if (finished) return;
	for (size_t i=0; i<nodes.size(); ++i) {
		if (nodes[i] != 0) {
			printNode(i+offset);
			delete nodes[i];
		}
	}
	nodes.clear();
	finished = true;
}
