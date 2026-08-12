#ifndef SGBC_TESTS_H
#define SGBC_TESTS_H

#include <gtest/gtest.h>

extern "C" {
    int test_solve_tridiag_3x3_poisson();
    int test_solve_tridiag_5x5_poisson();
    int test_solve_tridiag_diagonal_system();
}

TEST(sgbc, solve_tridiag_3x3_poisson) {
    EXPECT_EQ(0, test_solve_tridiag_3x3_poisson());
}

TEST(sgbc, solve_tridiag_5x5_poisson) {
    EXPECT_EQ(0, test_solve_tridiag_5x5_poisson());
}

TEST(sgbc, solve_tridiag_diagonal_system) {
    EXPECT_EQ(0, test_solve_tridiag_diagonal_system());
}

#endif // SGBC_TESTS_H
