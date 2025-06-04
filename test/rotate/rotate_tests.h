#include <gtest/gtest.h>

extern "C" int test_simple_rotate();
TEST(rotate, rotate_case)    { EXPECT_EQ(0, test_simple_rotate()); }