#ifndef TEST_SYSTEM_INIT_SOLVER_H
#define TEST_SYSTEM_INIT_SOLVER_H

#include <gtest/gtest.h>

#include "semba_fdtd.h"

#include <filesystem>
#include <string>

TEST(SystemInitSolver, MatchesFortranHyHzAfterExPulse) {
    const std::string json_path =
        std::filesystem::path("test") / "system" / "init_solver.fdtd.json";
    ASSERT_TRUE(std::filesystem::exists(json_path));
    EXPECT_EQ(SEMBA_FDTD_m::SEMBA_FDTD_test::run_init_solver_test(json_path), 0);
}

#endif
