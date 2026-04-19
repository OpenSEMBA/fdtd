#include <gtest/gtest.h>

extern "C" int test_nfde2sgg_planewave();

TEST(semba_fdtd, nfde2sgg_planewave) { EXPECT_EQ(0, test_nfde2sgg_planewave()); }
