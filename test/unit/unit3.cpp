#include "gtest/gtest.h"

//#include <algorithm>
//#include <unordered_map>

//#include <boost/algorithm/string.hpp>

//#include "AlignmentRecord.h"
//#include "NewEdgeCalculator.h"

#include <vector>
#include <algorithm>
#include <math.h>
#include <ctype.h>
#include <map>
#include <array>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/compare.hpp>

#include "AlignmentRecord.h"
#include "Clique.h"

using namespace std;
using namespace boost;

// This test verifies if the algorithm correctly merge two single AlignmentRecords that satisfies the following condition:
// --------
// --------
TEST(mergeAlignmentRecordsFunctionTest, mergeAlignmentRecordsMatch){
    
    AlignmentRecord ar;
    AlignmentRecord alignment;
    AlignmentRecord trueMergeAlignment;
    int id = 1;
    
    AlignmentRecord* estimatedAlignmentRecord;
    //unique_ptr<vector<const AlignmentRecord*> > result(new vector<const AlignmentRecord*>());
    std::unique_ptr<std::vector<const AlignmentRecord*>> alignments(new vector<const AlignmentRecord*>());
    
    try{
    alignment.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample00.txt");
    ar.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample10.txt");
    }catch (...){
        cout << "File not found!" << endl;
    }
    alignments->push_back(&alignment);
    alignments->push_back(&ar);

    estimatedAlignmentRecord = new AlignmentRecord(alignments, id);
    
    trueMergeAlignment.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample11.txt");
    
    bool mergeRes = ((*estimatedAlignmentRecord)==trueMergeAlignment);
    
    EXPECT_EQ(mergeRes, true);
}

// This test verifies if the algorithm correctly merge two single AlignmentRecords that satisfies the following condition:
// --------
// --------
TEST(mergeAlignmentRecordsFunctionTest, mergeAlignmentRecordsDebug){
    
    AlignmentRecord ar;
    AlignmentRecord alignment;
    AlignmentRecord trueMergeAlignment;
    int id = 1;
    
    AlignmentRecord* estimatedAlignmentRecord;
    std::unique_ptr<std::vector<const AlignmentRecord*>> alignments(new vector<const AlignmentRecord*>());
    
    try{
        alignment.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample12.txt");
        ar.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample13.txt");
    }catch (...){
        cout << "File not found!" << endl;
    }
    
    alignments->push_back(&alignment);
    alignments->push_back(&ar);

    estimatedAlignmentRecord = new AlignmentRecord(alignments, id);
    
    trueMergeAlignment.restoreCompleteAlignmentRecord("test/data/simulation/unit_data/alignment_sample14.txt");
    
    bool mergeRes = ((*estimatedAlignmentRecord)==trueMergeAlignment);
    
    EXPECT_EQ(mergeRes, true);
}
