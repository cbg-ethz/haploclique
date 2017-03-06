#include "gtest/gtest.h"

#include <algorithm>
//#include <deque>
//#include <vector>
#include <unordered_map>

#include <boost/algorithm/string.hpp>

#include "AlignmentRecord.h"
#include "NewEdgeCalculator.h"
//#include <unistd.h>

using namespace std;
using namespace boost;

// This test verifies if the algorithm considers an edge between twp identical AlignmentRecords.
TEST(edgeBetweenFunctionTest, edgeBetweenSameAlignmentRecord){
    
    AlignmentRecord alignment;
    EdgeCalculator* edge_calculator = nullptr;
    
    alignment.restoreAlignmentRecord("test/data/simulation/unit_data/alignment_sample.txt");
    
    double Q = 0.9;
    double edge_quasi_cutoff_cliques = 0.99;
    double overlap_cliques = 0.9;
    bool frameshift_merge = false;
    std::unordered_map<int, double> simpson_map;
    double edge_quasi_cutoff_single = 0.95;
    double overlap_single = 0.6;
    double edge_quasi_cutoff_mixed = 0.97;
    unsigned int maxPosition1 = 0;
    bool noProb0 = false;

    edge_calculator = new NewEdgeCalculator(Q, edge_quasi_cutoff_cliques, overlap_cliques, frameshift_merge, simpson_map, edge_quasi_cutoff_single, overlap_single, edge_quasi_cutoff_mixed, maxPosition1, noProb0);
    
    bool set_edge = edge_calculator->edgeBetween(alignment, alignment);
    
    EXPECT_EQ(set_edge, true);
    
    delete edge_calculator;
}


