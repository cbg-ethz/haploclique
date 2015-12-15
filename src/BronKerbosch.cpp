#include "BronKerbosch.h"
#include <boost/unordered_set.hpp>
#include <boost/dynamic_bitset.hpp>
#include <iostream>

// using namespace boost;
using namespace std;

BronKerbosch::BronKerbosch(const EdgeCalculator& edge_calculator, CliqueCollector& clique_collector, LogWriter* lw)
: CliqueFinder(edge_calculator, clique_collector), alignments_(), lw(lw) {
    order_ = nullptr;
    vertices_ = nullptr;
    degree_map_ = nullptr;
    actives_ = nullptr;
    vertices_as_lists_ = nullptr;
}

BronKerbosch::~BronKerbosch() {
    if (cliques != nullptr) {
        finish();
    }
    for (auto&& alignment : alignments_) {
        delete alignment;
    }
    delete order_;
    delete vertices_;
}

void BronKerbosch::initialize() {
    assert(not initialized);

    cliques = new clique_list_t();
    for (auto&& alignment : alignments_) {
        delete alignment;
    }
    alignments_.clear();
    delete order_;
    delete vertices_;
    order_ = nullptr;
    vertices_ = nullptr;
    degree_map_ = new degree_map_t();
    actives_ = new list<adjacency_list_t*>();
    vertices_as_lists_ = new vector<adjacency_list_t*>();

  	alignment_count = 0;
    next_id = 0;

    initialized = true;
    converged = true;
}

void BronKerbosch::finish() {
    assert(initialized);

    degeneracy_order();

    alignment_set_t R(alignment_count);
    alignment_set_t P(alignment_count);
    P.flip();
    alignment_set_t X(alignment_count);

    for (auto&& i : *order_) {

        bronkerbosch(R.set(i), P & (*vertices_)[i], X & (*vertices_)[i]);

        R.reset(i);
        P.reset(i);
        X.set(i);
    }

    for (auto&& clique : *cliques) {
        clique_collector.add(unique_ptr<Clique>(clique));
    }
    delete cliques;
    cliques = nullptr;
    initialized = false;
}

void BronKerbosch::degeneracy_order() {

    // Move all still active vertices into the degree map
    for (auto it = actives_->begin(); it != actives_->end();) {
        list<adjacency_list_t*>::size_type ind = (*it)->second.size();

        if(degree_map_->count(ind) == 0) {           
            degree_map_->emplace(ind, list<adjacency_list_t*>(1, *it) );
	//	degree_map_->insert(std::pair<std::list<adjacency_list_t*>::size_type, std::list<adjacency_list_t*>>(ind,list<adjacency_list_t*>(1, *it)));
	

        } else {
            degree_map_->at(ind).push_back(*it);
        }
        
        it = actives_->erase(it);
    }

    order_ = new list<size_t>();
    vertices_ = new vector<alignment_set_t>(alignment_count, alignment_set_t(alignment_count));

    while (not degree_map_->empty()) {
        auto it = degree_map_->begin();

        if (it->second.empty()) {
            degree_map_->erase(it);
            continue;
        }
        
        adjacency_list_t* list_vertex = it->second.front();
        it->second.pop_front();

        // Update degree map
        for (auto&& index : list_vertex->second) {
            adjacency_list_t* companion = (*vertices_as_lists_)[index];
            list<adjacency_list_t*>::size_type s = companion->second.size();

            (*vertices_)[list_vertex->first].set(companion->first);
            (*vertices_)[companion->first].set(list_vertex->first);

            companion->second.remove(list_vertex->first);

            degree_map_->at(s).remove(companion);

            if (degree_map_->count(s-1) == 0) {
                degree_map_->emplace(s-1, list<adjacency_list_t*>(1, companion) );                
		//degree_map_->insert(std::pair<std::list<adjacency_list_t*>::size_type, std::list<adjacency_list_t*>>(s-1, list<adjacency_list_t*>(1, companion)));
            } else {
                degree_map_->at(s-1).push_front(companion);
            }
        }
        order_->push_back(list_vertex->first);

        delete list_vertex;
    }

    delete degree_map_;
    delete actives_;
    delete vertices_as_lists_;
}

alignment_set_t::size_type BronKerbosch::find_pivot(const alignment_set_t& P, const alignment_set_t& X) {
    alignment_set_t unit = P | X;
    alignment_set_t::size_type max_size = 0;
    alignment_set_t::size_type index = unit.find_first();

    for(auto i = unit.find_next(index); i != alignment_set_t::npos; i = unit.find_next(i)) {
        auto new_size = (P & (*vertices_)[i]).count();        

        if ( new_size > max_size ) {
            max_size = new_size;
            index = i;
        }
    }

    assert(index != alignment_set_t::npos);

    return index;
}

void BronKerbosch::bronkerbosch(alignment_set_t R, alignment_set_t P, alignment_set_t X) {
    assert(cliquesistent(R, P, X));

    if (P.none()) {
        if (X.none()) {
            unique_ptr<alignment_set_t> ptr(new alignment_set_t(R));
            cliques->push_back(new Clique(*this,  ptr));
        }
        return;
    }

    alignment_set_t::size_type pivot = find_pivot(P, X);
    alignment_set_t comp = P - (*vertices_)[pivot];

    for (auto i = comp.find_first(); i != alignment_set_t::npos; i = comp.find_next(i)) {
        bronkerbosch(R.set(i), P & (*vertices_)[i], X & (*vertices_)[i]);

        R.reset(i);
        P.reset(i);
        X.set(i);
    }
}

/*
 * Checks the invariant of the BronKerbosch algorithm
 */
bool BronKerbosch::cliquesistent(const alignment_set_t& R, const alignment_set_t& P, const alignment_set_t& X) {
    auto i = R.find_first();
    alignment_set_t intersect((*vertices_)[i]);
    intersect.set(i);

    for (i=R.find_next(i); i != alignment_set_t::npos; i = R.find_next(i)) {
        alignment_set_t tmp((*vertices_)[i]);
        tmp.set(i);

        intersect &= tmp;
    }
    if (not R.is_subset_of(intersect)) {
        cerr << "R - intersect:";
        printSet(cerr, R - intersect);
        cerr << "vertices in R:" << endl;
        for (auto i=R.find_first(); i != alignment_set_t::npos; i = R.find_next(i)) {
            cerr << i << ":";
            printSet(cerr, (*vertices_)[i]);
        }
        return false;
    }
    if (not (P | X).is_subset_of(intersect)) {
        cerr << "(P | X) - intersect:";
        printSet(cerr, (P | X) - intersect);
        cerr << "vertices in P:" << endl;
        for (auto i=P.find_first(); i != alignment_set_t::npos; i = P.find_next(i)) {
            cerr << i << ":";
            printSet(cerr, (*vertices_)[i]);
        }
        cerr << "vertices in X:" << endl;
        for (auto i=X.find_first(); i != alignment_set_t::npos; i = X.find_next(i)) {
            cerr << i << ":";
            printSet(cerr, (*vertices_)[i]);
        }

        return false;
    }
    return true;
}

void BronKerbosch::printSet(std::ostream& os, alignment_set_t set) {
    for (auto j = set.find_first(); j != alignment_set_t::npos; j = set.find_next(j)) {
        os << " " << j;
    }
    os << endl;
}

void BronKerbosch::printReads(std::ostream& os, alignment_set_t set) {
    auto j = set.find_first();
    os << alignments_[j]->getName();

    for (j = set.find_next(j); j != alignment_set_t::npos; j = set.find_next(j)) {
        os << "," << alignments_[j]->getName();
    }
    os << endl;
}

void BronKerbosch::printEdges(std::string filename) {
    ofstream edgefile(filename, ofstream::out);

    for(unsigned int i = 0; i < alignment_count; i++) {
        alignment_set_t node = (*vertices_)[i];
        edgefile << i << ":" << alignments_[i]->getName() << " ->";
        
        printSet(edgefile, node);
    }
}

void BronKerbosch::addAlignment(std::unique_ptr<AlignmentRecord>& alignment_autoptr) {
	assert(alignment_autoptr.get() != nullptr);
	assert(cliques!=nullptr);
    assert(initialized);

	alignment_id_t id = next_id++;
	AlignmentRecord* alignment = alignment_autoptr.release();
	alignment->setID(id);

	size_t index = alignment_count++;
	alignments_.push_back(alignment);

    adjacency_list_t* vertex = new adjacency_list_t(index, list<size_t>());

    vertices_as_lists_->push_back(vertex);

    for (auto it = actives_->begin(); it != actives_->end();) {
        AlignmentRecord* alignment2 = alignments_[(*it)->first];

        if (alignment->getIntervalStart() > alignment2->getIntervalEnd()) {
            list<adjacency_list_t*>::size_type ind = (*it)->second.size();

            if(degree_map_->count(ind) == 0) {           
                degree_map_->emplace(ind, list<adjacency_list_t*>(1, *it) );
//		degree_map_->insert(std::pair<std::list<adjacency_list_t*>::size_type, std::list<adjacency_list_t*>> (ind, list<adjacency_list_t*>(1, *it)));
            } else {
                degree_map_->at(ind).push_back(*it);
            }
            
            it = actives_->erase(it);
            continue;
        }

        if(edge_calculator.edgeBetween(*alignment, *alignment2) ) {

            if(second_edge_calculator != nullptr) {
                if(not second_edge_calculator->edgeBetween(*alignment, *alignment2) ) {
                    continue;
                }
            }

            // Draw edge between current alignment and old alignment in vert
            get<1>(*vertex).push_back((*it)->first);
            get<1>(**it).push_back(vertex->first);

            if (lw != nullptr) lw->reportEdge(vertex->first, (*it)->first);

            converged = false;
        }
        it++;
    }

    actives_->push_back(vertex);

}
