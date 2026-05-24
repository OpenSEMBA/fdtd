#ifndef TEST_PLANEWAVE_PW_IN_BOX_H
#define TEST_PLANEWAVE_PW_IN_BOX_H

#include <gtest/gtest.h>

#include "semba_fdtd.h"

#include <filesystem>
#include <string>

namespace pw_in_box_test {

inline std::string casePath(const char* name) {
    return (std::filesystem::path("testData") / "cases" / "planewave" / name).string();
}

} // namespace pw_in_box_test

TEST(PlanewavePwInBox, ShortRunProbeParity_First50Steps) {
    const std::string json = pw_in_box_test::casePath("pw-in-box.fdtd.json");
    ASSERT_TRUE(std::filesystem::exists(json));
    const int err = SEMBA_FDTD_m::SEMBA_FDTD_test::test_run_pw_in_box_probes(
        json,
        pw_in_box_test::casePath("pw-in-box.fdtd_before_Ex_3_3_1.dat"),
        pw_in_box_test::casePath("pw-in-box.fdtd_inbox_Ex_3_3_3.dat"),
        pw_in_box_test::casePath("pw-in-box.fdtd_after_Ex_3_3_5.dat"),
        50);
    EXPECT_EQ(err, 0) << "pw-in-box short probe parity failed with code " << err;
}

TEST(PlanewavePwInBox, MediumRunProbeParity_First90Steps) {
    const std::string json = pw_in_box_test::casePath("pw-in-box.fdtd.json");
    ASSERT_TRUE(std::filesystem::exists(json));
    const int err = SEMBA_FDTD_m::SEMBA_FDTD_test::test_run_pw_in_box_probes(
        json,
        pw_in_box_test::casePath("pw-in-box.fdtd_before_Ex_3_3_1.dat"),
        pw_in_box_test::casePath("pw-in-box.fdtd_inbox_Ex_3_3_3.dat"),
        pw_in_box_test::casePath("pw-in-box.fdtd_after_Ex_3_3_5.dat"),
        90);
    EXPECT_EQ(err, 0) << "pw-in-box 90-step probe parity failed with code " << err;
}

TEST(PlanewavePwInBox, Step100ProbeParity) {
    const std::string json = pw_in_box_test::casePath("pw-in-box.fdtd.json");
    ASSERT_TRUE(std::filesystem::exists(json));
    const int err = SEMBA_FDTD_m::SEMBA_FDTD_test::test_run_pw_in_box_probes(
        json,
        pw_in_box_test::casePath("pw-in-box.fdtd_before_Ex_3_3_1.dat"),
        pw_in_box_test::casePath("pw-in-box.fdtd_inbox_Ex_3_3_3.dat"),
        pw_in_box_test::casePath("pw-in-box.fdtd_after_Ex_3_3_5.dat"),
        100);
    EXPECT_EQ(err, 0) << "pw-in-box 100-step probe parity failed with code " << err;
}

TEST(PlanewavePwInBox, Step120ProbeParity) {
    const std::string json = pw_in_box_test::casePath("pw-in-box.fdtd.json");
    ASSERT_TRUE(std::filesystem::exists(json));
    const int err = SEMBA_FDTD_m::SEMBA_FDTD_test::test_run_pw_in_box_probes(
        json,
        pw_in_box_test::casePath("pw-in-box.fdtd_before_Ex_3_3_1.dat"),
        pw_in_box_test::casePath("pw-in-box.fdtd_inbox_Ex_3_3_3.dat"),
        pw_in_box_test::casePath("pw-in-box.fdtd_after_Ex_3_3_5.dat"),
        120);
    EXPECT_EQ(err, 0) << "pw-in-box 120-step probe parity failed with code " << err;
}

TEST(PlanewavePwInBox, FullRunProbeParity) {
    const std::string json = pw_in_box_test::casePath("pw-in-box.fdtd.json");
    ASSERT_TRUE(std::filesystem::exists(json));
    const int err = SEMBA_FDTD_m::SEMBA_FDTD_test::test_run_pw_in_box_probes(
        json,
        pw_in_box_test::casePath("pw-in-box.fdtd_before_Ex_3_3_1.dat"),
        pw_in_box_test::casePath("pw-in-box.fdtd_inbox_Ex_3_3_3.dat"),
        pw_in_box_test::casePath("pw-in-box.fdtd_after_Ex_3_3_5.dat"),
        -1);
    EXPECT_EQ(err, 0) << "pw-in-box full probe parity failed with code " << err;
}

#endif
