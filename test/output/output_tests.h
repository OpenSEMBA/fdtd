#ifdef CompileWithNewOutputModule

#include <gtest/gtest.h>

extern "C" int test_init_point_probe();
extern "C" int test_update_point_probe();
extern "C" int test_flush_point_probe();
extern "C" int test_multiple_flush_point_probe();
extern "C" int test_init_movie_probe();
extern "C" int test_update_movie_probe();
extern "C" int test_flush_movie_probe();
extern "C" int test_init_frequency_slice_probe();
extern "C" int test_update_frequency_slice_probe();
extern "C" int test_flush_frequency_slice_probe();

extern "C" int test_count_required_coords();
extern "C" int test_store_required_coords();
extern "C" int test_is_valid_point_current();
extern "C" int test_is_valid_point_field();


TEST(output, test_initialize_point_probe)    {EXPECT_EQ(0, test_init_point_probe()); }
TEST(output, test_update_point_probe_info)    {EXPECT_EQ(0, test_update_point_probe()); }
TEST(output, test_flush_point_probe_info)    {EXPECT_EQ(0, test_flush_point_probe()); }
TEST(output, test_flush_multiple_point_probe_info)    {EXPECT_EQ(0, test_multiple_flush_point_probe()); }
TEST(output, test_init_movie_probe_for_pec_surface)    {EXPECT_EQ(0, test_init_movie_probe()); }
TEST(output, test_update_movie_probe_for_pec_surface)    {EXPECT_EQ(0, test_update_movie_probe()); }
TEST(output, test_flush_movie_probe_data)    {EXPECT_EQ(0, test_flush_movie_probe()); }
TEST(output, test_init_frequency_slice)    {EXPECT_EQ(0, test_init_frequency_slice_probe()); }
TEST(output, test_update_frequency_slice)    {EXPECT_EQ(0, test_update_frequency_slice_probe()); }
TEST(output, test_flush_frequency_slice)    {EXPECT_EQ(0, test_flush_frequency_slice_probe()); }

TEST(output, test_volumic_utils_count) { EXPECT_EQ(0, test_count_required_coords()); }
TEST(output, test_volumic_utils_store) { EXPECT_EQ(0, test_store_required_coords()); }
TEST(output, test_volumic_utils_valid_current) { EXPECT_EQ(0, test_is_valid_point_current()); }
TEST(output, test_volumic_utils_valid_field) { EXPECT_EQ(0, test_is_valid_point_field()); }

#endif