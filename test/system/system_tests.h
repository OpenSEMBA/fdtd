#include <gtest/gtest.h>

extern "C" int test_init_solver();
extern "C" int test_rank_remapping();

#ifndef CompileWithMPI
TEST(system, init_solver)     {EXPECT_EQ(0, test_init_solver()); }
TEST(system, rank_remapping)     {EXPECT_EQ(0, test_rank_remapping()); }
#endif
