#include "gtest/gtest.h"

#include <algorithm>
#include <unordered_map>

#include <boost/algorithm/string.hpp>

#include "AlignmentRecord.h"
#include "NewEdgeCalculator.h"


using namespace std;
using namespace boost;

// This test verifies if the algorithm considers an edge between two identical AlignmentRecords.
TEST(edgeBetweenFunctionTest, edgeBetweenSameAlignments){
    
    AlignmentRecord alignment;
    EdgeCalculator* edge_calculator = nullptr;
    
    alignment.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample00.txt");
    
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

// This test verifies if the algorithm considers an edge between two identical AlignmentRecords with soft clip.
TEST(edgeBetweenFunctionTest, edgeBetweenSameAlignmentsSoftClip){
    
    AlignmentRecord alignment;
    EdgeCalculator* edge_calculator = nullptr;
    
    alignment.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample01.txt");
    
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

// This test verifies if the algorithm considers an edge between two similar AlignmentRecords with incompatible gap.
TEST(edgeBetweenFunctionTest, edgeBetweenIncompatibleAlignmentsStartDel){
    
    AlignmentRecord alignment1;
    AlignmentRecord alignment2;
    EdgeCalculator* edge_calculator = nullptr;
    try{
        alignment1.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample02.txt");
    }catch (...){
        cout << "File not found!" << endl;
    }
    
    try{
        alignment2.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample03.txt");
    }catch (...){
        cout << "File 2 not found!" << endl;
    }
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
    
    bool set_edge = edge_calculator->edgeBetween(alignment1, alignment2);
    
    EXPECT_EQ(set_edge, false);
    
    delete edge_calculator;
}

// This test verifies if the algorithm considers an edge between two non_overlapping AlignmentRecords.
TEST(edgeBetweenFunctionTest, edgeBetweenNonOverlappingAlignments){
    
    AlignmentRecord alignment1;
    AlignmentRecord alignment2;
    EdgeCalculator* edge_calculator = nullptr;
    try{
        alignment1.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample04.txt");
    }catch (...){
        cout << "File not found!" << endl;
    }
    try{
        alignment2.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample05.txt");
    }catch (...){
        cout << "File 2 not found!" << endl;
    }
    
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
    
    bool set_edge = edge_calculator->edgeBetween(alignment1, alignment2);
    
    EXPECT_EQ(set_edge, false);
    
    delete edge_calculator;
}

// This test verifies if the algorithm considers an edge between two partially overlapping AlignmentRecords.
// The second pair of the first read is aligned with the first pair of second read.
TEST(edgeBetweenFunctionTest, edgeBetweenOverlappingCrossPairsAlignments){
    
    AlignmentRecord alignment1;
    AlignmentRecord alignment2;
    EdgeCalculator* edge_calculator = nullptr;
    try{
        alignment1.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample06.txt");
    }catch (...){
        cout << "File not found!" << endl;
    }
    
    try{
        alignment2.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample07.txt");
    }catch (...){
        cout << "File 2 not found!" << endl;
    }
    
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
    
    bool set_edge = edge_calculator->edgeBetween(alignment1, alignment2);
    
    EXPECT_EQ(set_edge, false);
    
    delete edge_calculator;
}

// This test verifies if the algorithm considers an edge between two completely different AlignmentRecords.
TEST(edgeBetweenFunctionTest, edgeBetweenDissimilarAlignments){
    
    AlignmentRecord alignment1;
    AlignmentRecord alignment2;
    EdgeCalculator* edge_calculator = nullptr;
    try{
        alignment1.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample00.txt");
    }catch (...){
        cout << "File not found!" << endl;
    }
    
    try{
        alignment2.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample08.txt");
    }catch (...){
        cout << "File 2 not found!" << endl;
    }
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
    
    bool set_edge = edge_calculator->edgeBetween(alignment1, alignment2);
    
    EXPECT_EQ(set_edge, false);
    
    delete edge_calculator;
}

// This test verifies if the algorithm considers an edge between two similar AlignmentRecords that satisfies following condition:
//----  -----
//----          -----
TEST(edgeBetweenFunctionTest, edgeBetweenIncompatibleAlignments1){
    
    AlignmentRecord alignment1;
    AlignmentRecord alignment2;
    EdgeCalculator* edge_calculator = nullptr;
    try{
        alignment1.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample06.txt");
    }catch (...){
        cout << "File not found!" << endl;
    }
    
    try{
        alignment2.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample09.txt");
    }catch (...){
        cout << "File 2 not found!" << endl;
    }
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
    
    bool set_edge = edge_calculator->edgeBetween(alignment1, alignment2);
    
    EXPECT_EQ(set_edge, false);
    
    delete edge_calculator;
}

// This test verifies if the algorithm considers an edge between two similar AlignmentRecords that satisfies following condition:
//----          -----
//----  -----
TEST(edgeBetweenFunctionTest, edgeBetweenIncompatibleAlignments2){
    
    AlignmentRecord alignment1;
    AlignmentRecord alignment2;
    EdgeCalculator* edge_calculator = nullptr;
    try{
        alignment1.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample09.txt");
    }catch (...){
        cout << "File not found!" << endl;
    }
    
    try{
        alignment2.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample06.txt");
    }catch (...){
        cout << "File 2 not found!" << endl;
    }
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
    
    bool set_edge = edge_calculator->edgeBetween(alignment1, alignment2);
    
    EXPECT_EQ(set_edge, false);
    
    delete edge_calculator;
}
