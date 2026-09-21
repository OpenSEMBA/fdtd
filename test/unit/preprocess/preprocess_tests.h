#ifndef PREPROCESS_TESTS_H
#define PREPROCESS_TESTS_H

#include <gtest/gtest.h>

// External Fortran function declarations (bind(C) from test_preprocess_geom.F90)
extern "C" {
    int test_searchtag();
    int test_searchtag_empty();
    int test_searchtag_single();
    int test_searchtag_special_chars();
    int test_checkDielectricTag_no_dup();
    int test_checkDielectricTag_c2P_duplicate_current();
    int test_checkLossyTag_basic();
    int test_checkLossyTag_duplicate_current();
    int test_checkLossyTag_duplicate_previous();
}

// Test cases following the conformal_tests.h pattern
TEST(preprocess, searchtag) { 
    EXPECT_EQ(0, test_searchtag()); 
}

TEST(preprocess, searchtag_empty) { 
    EXPECT_EQ(0, test_searchtag_empty()); 
}

TEST(preprocess, searchtag_single) { 
    EXPECT_EQ(0, test_searchtag_single()); 
}

TEST(preprocess, searchtag_special_chars) { 
    EXPECT_EQ(0, test_searchtag_special_chars()); 
}

TEST(preprocess, checkDielectricTag_no_dup) { 
    EXPECT_EQ(0, test_checkDielectricTag_no_dup()); 
}

TEST(preprocess, checkDielectricTag_c2P_duplicate_current) {
    EXPECT_EQ(0, test_checkDielectricTag_c2P_duplicate_current());
}

TEST(preprocess, checkLossyTag_basic) { 
    EXPECT_EQ(0, test_checkLossyTag_basic()); 
}

TEST(preprocess, checkLossyTag_duplicate_current) {
    EXPECT_EQ(0, test_checkLossyTag_duplicate_current());
}

TEST(preprocess, checkLossyTag_duplicate_previous) {
    EXPECT_EQ(0, test_checkLossyTag_duplicate_previous());
}

#endif // PREPROCESS_TESTS_H
