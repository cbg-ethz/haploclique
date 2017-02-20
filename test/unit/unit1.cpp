#include "gtest/gtest.h"

#include <algorithm>
#include <deque>
#include <vector>
#include <unordered_map>

#include <boost/algorithm/string.hpp>

#include <api/BamReader.h>
#include <api/BamAlignment.h>

#include "AlignmentRecord.h"
#define main HaploMain
#include "haploclique.cpp"
#include <unistd.h> 

using namespace std;
using namespace boost;

// This test verifies if readBamFile function returns correct number of reads for a specific input file.
TEST(readBamFileTest, readBamFileCountSeqs){

    string bamfile = "test/data/simulation/reads_HIV-1_50_01.bam";
    vector<string> originalReadNames;
    unsigned int maxPosition1;
    BamTools::SamHeader header;
    BamTools::RefVector references;
    int readsSize = 1836;
    std::deque<AlignmentRecord*>* reads = readBamFile(bamfile, originalReadNames,maxPosition1,header,references);

    EXPECT_EQ(readsSize, reads->size());

}

// This test verifies if readBamFile function can retrieve any reads from the input file.
TEST(readBamFileExpectTest, readBamFileExpectReadsRetrival){
    
    string bamfile = "test/data/simulation/unit_data/reads_arabis_mosaic_virus_20nc.bam";
    vector<string> originalReadNames;
    unsigned int maxPosition1;
    BamTools::SamHeader header;
    BamTools::RefVector references;
    
    try{
        std::deque<AlignmentRecord*>* reads = readBamFile(bamfile, originalReadNames,maxPosition1,header,references);
        FAIL() << "Expected std::runtime_error" << endl;
    }
    catch(const std::runtime_error& error){
        EXPECT_EQ(error.what(), std::string("No reads could be retrieved from the BamFile. "));
    }
    catch(...){
        FAIL() << "Expected std::runtime_error" << endl;
    }
}

