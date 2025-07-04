#include <gtest/gtest.h>

extern "C" int test_init_solver();

TEST(system, init_solver)     {EXPECT_EQ(0, test_init_solver()); }

