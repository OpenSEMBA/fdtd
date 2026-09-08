#ifndef HEALER_TESTS_H
#define HEALER_TESTS_H

#include <gtest/gtest.h>

extern "C" {
    int test_sort_all_swapped();
    int test_sort_already_ordered();
    int test_sort_x_reversed_yz_ok();
    int test_readjust_grow();
    int test_readjust_shrink();
    int test_readjust_same_size();
}

TEST(healer, sort_all_swapped) {
    EXPECT_EQ(0, test_sort_all_swapped());
}

TEST(healer, sort_already_ordered) {
    EXPECT_EQ(0, test_sort_already_ordered());
}

TEST(healer, sort_x_reversed_yz_ok) {
    EXPECT_EQ(0, test_sort_x_reversed_yz_ok());
}

TEST(healer, readjust_grow) {
    EXPECT_EQ(0, test_readjust_grow());
}

TEST(healer, readjust_shrink) {
    EXPECT_EQ(0, test_readjust_shrink());
}

TEST(healer, readjust_same_size) {
    EXPECT_EQ(0, test_readjust_same_size());
}

#endif // HEALER_TESTS_H
