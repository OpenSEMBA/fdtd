#include <gtest/gtest.h>

extern "C" int test_allocate_serialize_for_time_domain();
extern "C" int test_allocate_serialize_for_frequency_domain();
extern "C" int test_allocate_current();

TEST(observation, test_allocate_time)    {EXPECT_EQ(0, test_allocate_serialize_for_time_domain()); }
TEST(observation, test_allocate_frequency)    {EXPECT_EQ(0, test_allocate_serialize_for_frequency_domain()); }
TEST(observation, test_allocate_serialize_current)    {EXPECT_EQ(0, test_allocate_current()); }

