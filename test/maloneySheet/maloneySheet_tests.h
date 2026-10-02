#ifndef MALONEYSHEET_TESTS_H
#define MALONEYSHEET_TESTS_H

#include <gtest/gtest.h>

extern "C" {
    int test_maloneysheet_e_update();
    int test_maloneysheet_h_correction();
}

TEST(maloneySheet, e_update) {
    EXPECT_EQ(0, test_maloneysheet_e_update());
}

TEST(maloneySheet, h_correction) {
    EXPECT_EQ(0, test_maloneysheet_h_correction());
}

#endif // MALONEYSHEET_TESTS_H
