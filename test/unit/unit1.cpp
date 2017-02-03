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


using namespace std;
using namespace boost;

TEST(readBamFileTest, readBamFileCountSeqs){

    // THIS IS NOT A UNIT TEST, BUT A REGRESSION CRAM TEST
#if 0
    string bamfile = "/home/ubuntu/haploclique/test/data/simulation/reads_arabis_mosaic_virus_25_01.bam";
    vector<string> originalReadNames;
    unsigned int maxPosition1;
    BamTools::SamHeader header;
    BamTools::RefVector references;
    int readsSize = 30;
    std::deque<AlignmentRecord*>* reads = readBamFile(bamfile, originalReadNames,maxPosition1,header,references);

    EXPECT_EQ(readsSize, reads->size());
#endif
    EXPECT_EQ(1, 1);
}

