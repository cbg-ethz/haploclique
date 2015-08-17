#ifndef BRONKERBOSCH_H_
#define BRONKERBOSCH_H_

#include <vector>
#include <map>
#include <list>
#include <utility>

#include "CliqueFinder.h"
#include "LogWriter.h"

class BronKerbosch : public CliqueFinder {
private:
    typedef std::pair<size_t, std::list<size_t>> adjacency_list_t;
    typedef std::map<std::list<adjacency_list_t*>::size_type, std::list<adjacency_list_t*>> degree_map_t;

    std::vector<AlignmentRecord*> alignments_;
    std::list<size_t>* order_;
    std::vector<alignment_set_t>* vertices_;
    degree_map_t* degree_map_;
    std::list<adjacency_list_t*>* actives_;
    std::vector<adjacency_list_t*>* vertices_as_lists_;
    LogWriter* lw;

    void degeneracy_order();
    void bronkerbosch(alignment_set_t R, alignment_set_t P, alignment_set_t X);
    alignment_set_t::size_type find_pivot(const alignment_set_t& P, const alignment_set_t& X);
    void printEdges(std::string filename);
    void printSet(std::ostream& os, alignment_set_t set);
    void printReads(std::ostream& os, alignment_set_t set);
    bool cliquesistent(const alignment_set_t& R, const alignment_set_t& P, const alignment_set_t& X);
public:
    BronKerbosch(const EdgeCalculator& edge_calculator, CliqueCollector& clique_collector, LogWriter* lw);
    virtual ~BronKerbosch();

    virtual const AlignmentRecord & getAlignmentByIndex(size_t index) const {
        assert(index<alignment_count);
    	return *(alignments_[index]);
    }

    void addAlignment(std::unique_ptr<AlignmentRecord>& ap);
    void initialize();
    void finish();
};

#endif // BRONKERBOSCH_H_
