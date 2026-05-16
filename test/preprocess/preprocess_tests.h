#ifndef PREPROCESS_TESTS_H
#define PREPROCESS_TESTS_H

#include <gtest/gtest.h>

// External Fortran function declarations (bind(C) from test_preprocess_geom.F90)
extern "C" {
    int test_searchtag();
    int test_searchtag_empty();
    int test_searchtag_single();
    int test_searchtag_special_chars();
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

#endif // PREPROCESS_TESTS_H
