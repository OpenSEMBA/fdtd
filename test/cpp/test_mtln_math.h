#ifndef TEST_MTLN_MATH_H
#define TEST_MTLN_MATH_H

#include <gtest/gtest.h>
#include <vector>

#include "test_mtln_helpers.h"

namespace utils_m {
std::vector<double> getEigenValues(const std::vector<std::vector<double>>& matrix);
}

TEST(mtln, math_eigvals) {
    using utils_m::getEigenValues;

    const std::vector<std::vector<double>> mat = {
        {0.35, 0.09, -0.44, 0.25},  {0.45, 0.07, -0.33, -0.32},
        {-0.14, -0.54, 0.03, -0.13}, { -0.17, 0.35, 0.17, 0.11}};

    const std::vector<double> ev = getEigenValues(mat);
    int error_cnt = 0;
    if (!mtln_test::checkNear(0.81630361, ev[0], 0.005) || !mtln_test::checkNear(0.0, ev[4], 0.005) ||
        !mtln_test::checkNear(-0.0988341, ev[1], 0.005) ||
        !mtln_test::checkNear(0.41323483, ev[5], 0.005) ||
        !mtln_test::checkNear(-0.0988341, ev[2], 0.005) ||
        !mtln_test::checkNear(-0.41323483, ev[6], 0.005) ||
        !mtln_test::checkNear(-0.05863542, ev[3], 0.005) ||
        !mtln_test::checkNear(0.0, ev[7], 0.005)) {
        error_cnt = 1;
    }
    EXPECT_EQ(error_cnt, 0);
}

TEST(mtln, math_matmul_broadcast) {
    int error_cnt = 0;
    std::vector<std::vector<std::vector<float>>> A(3, std::vector<std::vector<float>>(2, std::vector<float>(2)));
    std::vector<std::vector<std::vector<float>>> B(3, std::vector<std::vector<float>>(2, std::vector<float>(2)));
    std::vector<std::vector<std::vector<float>>> res1(3, std::vector<std::vector<float>>(2, std::vector<float>(2)));
    std::vector<std::vector<std::vector<float>>> res2(3, std::vector<std::vector<float>>(2, std::vector<float>(2)));

    for (int layer = 0; layer < 3; ++layer) {
        A[layer][0][0] = 1.0f;
        A[layer][1][0] = 0.0f;
        A[layer][0][1] = 0.0f;
        A[layer][1][1] = -1.0f;
        B[layer][0][0] = 1.5f;
        B[layer][1][0] = 0.0f;
        B[layer][0][1] = 0.5f;
        B[layer][1][1] = -1.0f;
    }

    for (int i = 0; i < 3; ++i) {
        for (int r = 0; r < 2; ++r) {
            for (int c = 0; c < 2; ++c) {
                float sum = 0.0f;
                for (int k = 0; k < 2; ++k) {
                    sum += A[i][r][k] * B[i][k][c];
                }
                res1[i][r][c] = sum;
            }
        }
        res2[i] = res1[i];
    }

    for (int i = 0; i < 3; ++i) {
        for (int r = 0; r < 2; ++r) {
            for (int c = 0; c < 2; ++c) {
                if (res1[i][r][c] != res2[i][r][c]) {
                    error_cnt += 1;
                }
            }
        }
    }
    EXPECT_EQ(error_cnt, 0);
}

#endif
