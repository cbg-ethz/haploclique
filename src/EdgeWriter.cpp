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

EdgeWriter::EdgeWriter(std::ostream& os) : os(os), vertex_to_read_names() {
	finished = true;
}

EdgeWriter::~EdgeWriter() {
	finish();
}

void EdgeWriter::initialize() {
    nodes.clear();
    vertex_to_read_names.clear();

    finished = false;
}

void EdgeWriter::addEdge(const AlignmentRecord& node1, const AlignmentRecord& node2) {
	assert(!finished);

	nodes[node1.getID()].insert(node2.getID());
	nodes[node2.getID()].insert(node1.getID());

    vertex_to_read_names[node1.getID()] = node1.getName();
    vertex_to_read_names[node2.getID()] = node2.getName();
}

void EdgeWriter::printNode(alignment_id_t id) {
	vector<alignment_id_t> adj_nodes(nodes[id].begin(), nodes[id].end());
	sort(adj_nodes.begin(), adj_nodes.end());
	os << id << ":" << vertex_to_read_names[id] << " ->";
	for (size_t i=0; i<adj_nodes.size(); ++i) {
		os << " " << adj_nodes[i];
	}
	os << endl;
}

void EdgeWriter::setNodeCompleted(const AlignmentRecord& node) {
	assert(!finished);
}

void EdgeWriter::finish() {
	if (finished) return;
	for (auto const &it : nodes) {
		printNode(it.first);
	}
	nodes.clear();
	finished = true;
}
