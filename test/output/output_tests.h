#include <gtest/gtest.h>

extern "C" int test_init_point_probe();
extern "C" int test_update_point_probe();

TEST(output, test_initialize_point_probe)    {EXPECT_EQ(0, test_init_point_probe()); }
TEST(output, test_update_point_probe)    {EXPECT_EQ(0, test_update_point_probe()); }
