#include "gtest/gtest.h"

TEST(intAddition, Negative) {
	EXPECT_EQ( -5, (-2 + -3)) << "in case of fail";
}

//double ReadSetZTester::computeSignificance(int insert_size_sum, int count, const std::vector<int>& read_group_counts) const
