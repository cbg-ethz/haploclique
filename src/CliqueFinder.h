#ifndef CLIQUEFINDER_H_
#define CLIQUEFINDER_H_

#include "CliqueCollector.h"
#include "AlignmentRecord.h"
#include "EdgeCalculator.h"

using namespace std;
using namespace boost;

class CliqueFinder {
protected:
    const EdgeCalculator & edge_calculator;
    CliqueCollector & clique_collector;
    bool initialized;
    bool converged;

    typedef std::list<Clique*> clique_list_t;
    clique_list_t *cliques;    

    size_t alignment_count;
    const EdgeCalculator *second_edge_calculator;
    alignment_id_t next_id;

public:
    CliqueFinder(const EdgeCalculator& edge_calculator, CliqueCollector& clique_collector) : edge_calculator(edge_calculator), clique_collector(clique_collector) {
    	alignment_count = 0;
    	second_edge_calculator = nullptr;
        next_id = 0;
        initialized = false;
        converged = false;        
    }
    
    virtual ~CliqueFinder() {
        
    }

    virtual void addAlignment(std::unique_ptr<AlignmentRecord>& ap, int& edgecounter) = 0;
    virtual void initialize() = 0;
    virtual void finish() = 0;
    virtual const AlignmentRecord & getAlignmentByIndex(size_t index) const = 0;
    /** Add a second edge calculator: edges will only be drawn when both edge calculators agree
	 *  that the edge is present.
	 */
    virtual void setSecondEdgeCalculator(const EdgeCalculator* ec) { this->second_edge_calculator = ec; }
    virtual bool hasConverged() { return converged; };
};

#endif
