#include <gtest/gtest.h>

extern "C" int test_rotate_spacesteps();
TEST(rotate, rotate_case)    { EXPECT_EQ(0, test_rotate_spacesteps()); }