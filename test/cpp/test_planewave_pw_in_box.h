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

inline int runProbeParity(int max_steps) {
    return SEMBA_FDTD_m::SEMBA_FDTD_test::test_run_pw_in_box_probes(
        casePath("pw-in-box.fdtd.json"),
        casePath("pw-in-box.fdtd_before_Ex_3_3_1.dat"),
        casePath("pw-in-box.fdtd_inbox_Ex_3_3_3.dat"),
        casePath("pw-in-box.fdtd_after_Ex_3_3_5.dat"),
        max_steps);
}

inline int runProbeFilesExact(int max_steps) {
    return SEMBA_FDTD_m::SEMBA_FDTD_test::test_run_pw_in_box_probe_files_exact(
        casePath("pw-in-box.fdtd.json"),
        casePath("pw-in-box.fdtd_before_Ex_3_3_1.dat"),
        casePath("pw-in-box.fdtd_inbox_Ex_3_3_3.dat"),
        casePath("pw-in-box.fdtd_after_Ex_3_3_5.dat"),
        max_steps);
}

inline int runPeriodicProbeFilesExact(int max_steps) {
    return SEMBA_FDTD_m::SEMBA_FDTD_test::test_run_pw_in_box_probe_files_exact(
        casePath("pw-with-periodic.fdtd.json"),
        casePath("pw-with-periodic.fdtd_before_Ex_3_3_1.dat"),
        casePath("pw-with-periodic.fdtd_inbox_Ex_3_3_3.dat"),
        casePath("pw-with-periodic.fdtd_after_Ex_3_3_5.dat"),
        max_steps);
}

} // namespace pw_in_box_test

TEST(PlanewavePwInBox, ShortRunProbeParity_First50Steps) {
    const std::string json = pw_in_box_test::casePath("pw-in-box.fdtd.json");
    ASSERT_TRUE(std::filesystem::exists(json));
    const int err = pw_in_box_test::runProbeParity(50);
    EXPECT_EQ(err, 0) << "pw-in-box short probe parity failed with code " << err;
}

TEST(PlanewavePwInBox, MediumRunProbeParity_First90Steps) {
    const std::string json = pw_in_box_test::casePath("pw-in-box.fdtd.json");
    ASSERT_TRUE(std::filesystem::exists(json));
    const int err = pw_in_box_test::runProbeParity(90);
    EXPECT_EQ(err, 0) << "pw-in-box 90-step probe parity failed with code " << err;
}

TEST(PlanewavePwInBox, Step100ProbeParity) {
    const std::string json = pw_in_box_test::casePath("pw-in-box.fdtd.json");
    ASSERT_TRUE(std::filesystem::exists(json));
    const int err = pw_in_box_test::runProbeParity(100);
    EXPECT_EQ(err, 0) << "pw-in-box 100-step probe parity failed with code " << err;
}

TEST(PlanewavePwInBox, Step120ProbeParity) {
    const std::string json = pw_in_box_test::casePath("pw-in-box.fdtd.json");
    ASSERT_TRUE(std::filesystem::exists(json));
    const int err = pw_in_box_test::runProbeParity(120);
    EXPECT_EQ(err, 0) << "pw-in-box 120-step probe parity failed with code " << err;
}

TEST(PlanewavePwInBox, ProbeFilesExact_First120Steps) {
#if !defined(CompileWithRelease)
    GTEST_SKIP() << "Exact probe-file parity is enforced in Release builds.";
#endif
    const std::string json = pw_in_box_test::casePath("pw-in-box.fdtd.json");
    ASSERT_TRUE(std::filesystem::exists(json));
    const int err = pw_in_box_test::runProbeFilesExact(120);
    EXPECT_EQ(err, 0) << "pw-in-box exact probe-file parity through step 120 failed with code " << err;
}

TEST(PlanewavePwInBox, FullRunProbeParity) {
    const std::string json = pw_in_box_test::casePath("pw-in-box.fdtd.json");
    ASSERT_TRUE(std::filesystem::exists(json));
    const int err = pw_in_box_test::runProbeParity(-1);
    EXPECT_EQ(err, 0) << "pw-in-box full probe parity failed with code " << err;
}

TEST(PlanewavePwInBox, ProbeFilesExact_FullRun) {
#if !defined(CompileWithRelease)
    GTEST_SKIP() << "Exact probe-file parity is enforced in Release builds.";
#endif
    const std::string json = pw_in_box_test::casePath("pw-in-box.fdtd.json");
    ASSERT_TRUE(std::filesystem::exists(json));
    const int err = pw_in_box_test::runProbeFilesExact(-1);
    EXPECT_EQ(err, 0) << "pw-in-box exact full probe-file parity failed with code " << err;
}

TEST(PlanewavePwInBox, PeriodicProbeFilesExact_FullRun) {
#if !defined(CompileWithRelease)
    GTEST_SKIP() << "Exact probe-file parity is enforced in Release builds.";
#endif
    const std::string json = pw_in_box_test::casePath("pw-with-periodic.fdtd.json");
    ASSERT_TRUE(std::filesystem::exists(json));
    const int err = pw_in_box_test::runPeriodicProbeFilesExact(-1);
    EXPECT_EQ(err, 0) << "pw-with-periodic exact full probe-file parity failed with code " << err;
}

#endif
