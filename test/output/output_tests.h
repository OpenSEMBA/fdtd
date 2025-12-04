#include <gtest/gtest.h>

extern "C" int test_initialize();
extern "C" int test_init_point_probe();

TEST(output, test_initialize)    {EXPECT_EQ(0, test_initialize()); }
TEST(output, test_initialize_point_probe)    {EXPECT_EQ(0, test_init_point_probe()); }
