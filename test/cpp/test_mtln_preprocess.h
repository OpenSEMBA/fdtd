#ifndef TEST_MTLN_PREPROCESS_H
#define TEST_MTLN_PREPROCESS_H

#include <gtest/gtest.h>
#include <numeric>
#include <vector>

#include "test_mtln_helpers.h"

using mtln_test::buildFourLevelBundle;
using mtln_test::buildSevenCableBundle;

TEST(mtln, preprocess_conductors_before_cable) {
    const auto line_bundle = buildFourLevelBundle();
    int error_cnt = 0;

    if (mtln_preprocess_test::findConductorsBeforeCable("line2", line_bundle.levels[1]) != 0) {
        error_cnt += 1;
    }
    if (mtln_preprocess_test::findConductorsBeforeCable("line3_1", line_bundle.levels[2]) != 0) {
        error_cnt += 1;
    }
    if (mtln_preprocess_test::findConductorsBeforeCable("line3_2", line_bundle.levels[2]) != 2) {
        error_cnt += 1;
    }
    if (mtln_preprocess_test::findConductorsBeforeCable("line4", line_bundle.levels[3]) != 0) {
        error_cnt += 1;
    }
    EXPECT_EQ(error_cnt, 0);
}

TEST(mtln, preprocess_conductors_in_level) {
    const auto line_bundle = buildFourLevelBundle();
    const auto got = mtln_preprocess_test::conductorsInLevel(line_bundle);
    const std::vector<int> expected = {1, 2, 4, 2};
    EXPECT_EQ(got, expected);
}

TEST(mtln, preprocess_zt_conductor_ranges) {
    const auto line_bundle = buildFourLevelBundle();
    const auto conductors_in_level = mtln_preprocess_test::conductorsInLevel(line_bundle);

    const std::vector<int> expected_out = {1, 2, 3, 7};
    const std::vector<std::vector<int>> expected_in = {{2, 3}, {4, 5}, {6, 7}, {8, 9}};

    int error_cnt = 0;
    int cnt = 0;
    for (size_t i = 1; i < line_bundle.levels.size(); ++i) {
        for (const auto& line : line_bundle.levels[i].lines) {
            const int conductor_out = mtln_preprocess_test::findOuterConductorNumber(
                line, line_bundle.levels[i - 1],
                mtln_preprocess_test::sumFirstN(conductors_in_level, static_cast<int>(i) - 1));
            const auto range_in = mtln_preprocess_test::findInnerConductorRange(
                line, line_bundle.levels[i],
                mtln_preprocess_test::sumFirstN(conductors_in_level, static_cast<int>(i)));

            if (expected_out[cnt] != conductor_out) {
                error_cnt += 1;
            }
            if (range_in != expected_in[cnt]) {
                error_cnt += 1;
            }
            cnt += 1;
        }
    }
    EXPECT_EQ(error_cnt, 0);
}

TEST(mtln, preprocess_zt_conductor_ranges_2) {
    const auto line_bundle = buildSevenCableBundle();
    const auto conductors_in_level = mtln_preprocess_test::conductorsInLevel(line_bundle);

    const std::vector<int> expected_out = {1, 2, 3, 4, 7, 8, 9};
    const std::vector<std::vector<int>> expected_in = {
        {2, 3, 4}, {5}, {6, 7}, {8, 9}, {10, 11}, {12}, {13}};

    int error_cnt = 0;
    int cnt = 0;
    for (size_t i = 1; i < line_bundle.levels.size(); ++i) {
        for (const auto& line : line_bundle.levels[i].lines) {
            const int conductor_out = mtln_preprocess_test::findOuterConductorNumber(
                line, line_bundle.levels[i - 1],
                mtln_preprocess_test::sumFirstN(conductors_in_level, static_cast<int>(i) - 1));
            const auto range_in = mtln_preprocess_test::findInnerConductorRange(
                line, line_bundle.levels[i],
                mtln_preprocess_test::sumFirstN(conductors_in_level, static_cast<int>(i)));

            if (expected_out[cnt] != conductor_out) {
                error_cnt += 1;
            }
            if (range_in != expected_in[cnt]) {
                error_cnt += 1;
            }
            cnt += 1;
        }
    }
    EXPECT_EQ(error_cnt, 0);
}

#endif
