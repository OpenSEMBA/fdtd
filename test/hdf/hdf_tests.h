#include <gtest/gtest.h>

extern "C" int test_hdf5_writing_and_reading();

TEST(hdf, writing_and_reading)     {EXPECT_EQ(0, test_hdf5_writing_and_reading()); }
