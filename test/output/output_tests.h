#include <gtest/gtest.h>

extern "C" int test_initialize();


TEST(output, test_initialize             )    {EXPECT_EQ(0, test_initialize()); }
