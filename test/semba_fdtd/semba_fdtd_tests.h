#include <gtest/gtest.h>

extern "C" int test_nfde2sgg_planewave();
extern "C" int test_nfde2sgg_planewave_periodic();
extern "C" int test_nfde2sgg_holland();
extern "C" int test_nfde2sgg_sources_voltage();

TEST(semba_fdtd, nfde2sgg_planewave) { EXPECT_EQ(0, test_nfde2sgg_planewave()); }
TEST(semba_fdtd, nfde2sgg_planewave_periodic) { EXPECT_EQ(0, test_nfde2sgg_planewave_periodic()); }
TEST(semba_fdtd, nfde2sgg_holland) { EXPECT_EQ(0, test_nfde2sgg_holland()); }
TEST(semba_fdtd, nfde2sgg_sources_voltage) { EXPECT_EQ(0, test_nfde2sgg_sources_voltage()); }
