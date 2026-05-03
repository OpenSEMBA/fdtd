#include <gtest/gtest.h>

extern "C" int test_searchtag_found();
extern "C" int test_searchtag_notfound();
extern "C" int test_searchtag_empty();

extern "C" int test_check_dielectric_tags_c1p_no_duplicates();
extern "C" int test_check_dielectric_tags_c1p_self_duplicate();
extern "C" int test_check_dielectric_tags_precounting_zero();
extern "C" int test_check_dielectric_tags_c2p();

extern "C" int test_check_animated_tags_c1p_no_duplicates();
extern "C" int test_check_animated_tags_c1p_self_duplicate();
extern "C" int test_check_animated_tags_c2p();

extern "C" int test_check_lossy_tags_no_duplicates();
extern "C" int test_check_lossy_tags_self_duplicate();
extern "C" int test_check_lossy_tags_precounting_zero();
extern "C" int test_check_lossy_tags_with_prev_duplicate();

TEST(preprocess, searchtag_found)    { EXPECT_EQ(0, test_searchtag_found()); }
TEST(preprocess, searchtag_notfound) { EXPECT_EQ(0, test_searchtag_notfound()); }
TEST(preprocess, searchtag_empty)    { EXPECT_EQ(0, test_searchtag_empty()); }

TEST(preprocess, check_dielectric_tags_c1p_no_duplicates) { EXPECT_EQ(0, test_check_dielectric_tags_c1p_no_duplicates()); }
TEST(preprocess, check_dielectric_tags_c1p_self_duplicate){ EXPECT_EQ(0, test_check_dielectric_tags_c1p_self_duplicate()); }
TEST(preprocess, check_dielectric_tags_precounting_zero)  { EXPECT_EQ(0, test_check_dielectric_tags_precounting_zero()); }
TEST(preprocess, check_dielectric_tags_c2p)               { EXPECT_EQ(0, test_check_dielectric_tags_c2p()); }

TEST(preprocess, check_animated_tags_c1p_no_duplicates)   { EXPECT_EQ(0, test_check_animated_tags_c1p_no_duplicates()); }
TEST(preprocess, check_animated_tags_c1p_self_duplicate)  { EXPECT_EQ(0, test_check_animated_tags_c1p_self_duplicate()); }
TEST(preprocess, check_animated_tags_c2p)                 { EXPECT_EQ(0, test_check_animated_tags_c2p()); }

TEST(preprocess, check_lossy_tags_no_duplicates)          { EXPECT_EQ(0, test_check_lossy_tags_no_duplicates()); }
TEST(preprocess, check_lossy_tags_self_duplicate)         { EXPECT_EQ(0, test_check_lossy_tags_self_duplicate()); }
TEST(preprocess, check_lossy_tags_precounting_zero)       { EXPECT_EQ(0, test_check_lossy_tags_precounting_zero()); }
TEST(preprocess, check_lossy_tags_with_prev_duplicate)    { EXPECT_EQ(0, test_check_lossy_tags_with_prev_duplicate()); }
