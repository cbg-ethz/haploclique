#ifndef CLIQUEFINDER_H_
#define CLIQUEFINDER_H_

#include "CliqueCollector.h"
#include "AlignmentRecord.h"
#include "EdgeCalculator.h"
#include "CoverageMonitor.h"
#include "EdgeWriter.h"
#include "ReadGroups.h"


using namespace std;
using namespace boost;

class CliqueFinder {
protected:
    const EdgeCalculator & edge_calculator;
    CliqueCollector & clique_collector;
    CoverageMonitor coverage_monitor;
    bool initialized;
    bool converged;

    typedef std::list<Clique*> clique_list_t;
    clique_list_t *cliques;    

    size_t alignment_count;
    EdgeWriter *edge_writer;
    const EdgeCalculator *second_edge_calculator;
    alignment_id_t next_id;

public:
    CliqueFinder(const EdgeCalculator& edge_calculator, CliqueCollector& clique_collector, const ReadGroups* read_groups) : edge_calculator(edge_calculator), clique_collector(clique_collector), coverage_monitor(read_groups) {
    	alignment_count = 0;
    	edge_writer = nullptr;
    	second_edge_calculator = nullptr;
        next_id = 0;
        initialized = false;
        converged = false;        
    }
    virtual ~CliqueFinder() {
        
    }

    virtual void addAlignment(std::unique_ptr<AlignmentRecord>& ap) = 0;

    virtual void initialize() = 0;

    virtual void finish() = 0;

    virtual CoverageMonitor & getCoverageMonitor() { return coverage_monitor; }

    virtual const AlignmentRecord & getAlignmentByIndex(size_t index) const = 0;

    virtual void setEdgeWriter(EdgeWriter& edge_writer) { this->edge_writer = &edge_writer; }

    /** Add a second edge calculator: edges will only be drawn when both edge calculators agree
	 *  that the edge is present.
	 */
    virtual void setSecondEdgeCalculator(const EdgeCalculator* ec) { this->second_edge_calculator = ec; }

    virtual bool hasConverged() { return converged; };
};

#endif
