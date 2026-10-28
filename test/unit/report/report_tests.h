#ifndef REPORT_TESTS_H
#define REPORT_TESTS_H

#include <gtest/gtest.h>

extern "C" {
    int test_log_line_format();
    int test_log_msg_sets_fatal_error();
}

TEST(report, log_line_format) {
    EXPECT_EQ(0, test_log_line_format());
}

TEST(report, log_msg_sets_fatal_error) {
    EXPECT_EQ(0, test_log_msg_sets_fatal_error());
}

#endif // REPORT_TESTS_H
