#include <gtest/gtest.h>

extern "C" int test_allocate_serialize_for_time_domain();
extern "C" int test_allocate_serialize_for_frequency_domain();
extern "C" int test_allocate_current();

extern "C" int test_initial_time_less_than_timestep();
extern "C" int test_timestep_greater_and_mapvtk();
extern "C" int test_timestep_greater_not_mapvtk();
extern "C" int test_freqstep_zero_or_large();
extern "C" int test_volumic_false_true_and_saveall();
extern "C" int test_saveall_branch();
extern "C" int test_final_less_than_initial();
extern "C" int test_huge_cap();

extern "C" int test_init_time_movie_observation();

TEST(observation, test_allocate_time             )    {EXPECT_EQ(0, test_allocate_serialize_for_time_domain()); }
TEST(observation, test_allocate_frequency        )    {EXPECT_EQ(0, test_allocate_serialize_for_frequency_domain()); }
TEST(observation, test_allocate_serialize_current)    {EXPECT_EQ(0, test_allocate_current()); }

TEST(observation, test_preproces_initial_time_less_than_timestep)    {EXPECT_EQ(0, test_initial_time_less_than_timestep()); }
TEST(observation, test_preproces_timestep_greater_and_mapvtk    )    {EXPECT_EQ(0, test_timestep_greater_and_mapvtk()); }
TEST(observation, test_preproces_timestep_greater_not_mapvtk    )    {EXPECT_EQ(0, test_timestep_greater_not_mapvtk()); }
TEST(observation, test_preproces_freqstep_zero_or_large         )    {EXPECT_EQ(0, test_freqstep_zero_or_large()); }
TEST(observation, test_preproces_volumic_false_true_and_saveall )    {EXPECT_EQ(0, test_volumic_false_true_and_saveall()); }
TEST(observation, test_preproces_saveall_branch                 )    {EXPECT_EQ(0, test_saveall_branch()); }
TEST(observation, test_preproces_final_less_than_initial        )    {EXPECT_EQ(0, test_final_less_than_initial()); }
TEST(observation, test_preproces_huge_cap                       )    {EXPECT_EQ(0, test_huge_cap()); }

TEST(observation, test_init_movie_observation                   )    {EXPECT_EQ(0, test_init_time_movie_observation()); }

