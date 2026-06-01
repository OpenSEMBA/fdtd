#ifndef TEST_PLANEWAVE_TFSF_H
#define TEST_PLANEWAVE_TFSF_H

#include <gtest/gtest.h>

#include "semba_fdtd.h"

#include <cmath>
#include <filesystem>
#include <string>

// +z propagation: the Fortran single-precision golden output is still effectively
// zero at this early time; the stronger parity tests cover the later wave arrival.
TEST(PlanewaveTFSF, InboxExNearZeroBeforeWaveArrival) {
    const std::string json =
        (std::filesystem::path("testData") / "cases" / "planewave" / "pw-in-box.fdtd.json").string();
    ASSERT_TRUE(std::filesystem::exists(json));
    const double ex_inbox = SEMBA_FDTD_m::SEMBA_FDTD_test::test_field_after_tfsf_e_step(json, 0, 3, 3, 3);
    EXPECT_NEAR(ex_inbox, 0.0, 1e-20);
}

TEST(PlanewaveTFSF, BeforeProbeStaysNearZeroAfterFewSteps) {
    const std::string json =
        (std::filesystem::path("testData") / "cases" / "planewave" / "pw-in-box.fdtd.json").string();
    const double ex_before = SEMBA_FDTD_m::SEMBA_FDTD_test::test_field_after_tfsf_e_step(json, 0, 3, 3, 1);
    EXPECT_NEAR(ex_before, 0.0, 1e-3);
}

#endif
