#ifndef TEST_BORDERSMUR_H
#define TEST_BORDERSMUR_H

#include <gtest/gtest.h>

#include "semba_fdtd.h"

#include <fstream>
#include <filesystem>
#include <limits>
#include <string>

namespace bordersmur_test {

constexpr double C0 = 299792458.0;

inline std::string pulse1dJson() {
    return (std::filesystem::path("testData") / "cases" / "mur" / "pulse-1d-x.fdtd.json").string();
}

inline double expectedMurCx(double dt, double dx) {
#ifdef CompileWithReal8
    const double cnum = dx / (dt * C0);
    return (1.0 - cnum) / (1.0 + cnum);
#else
    const auto one = static_cast<float>(1.0f);
    const auto inv = std::nextafter(one / static_cast<float>(dx),
                                    -std::numeric_limits<float>::infinity());
    const auto cluz = static_cast<float>(C0);
    const auto cnum = static_cast<float>(
        static_cast<double>(one / inv) / (dt * static_cast<double>(cluz)));
    return static_cast<double>((one - cnum) / (one + cnum));
#endif
}

} // namespace bordersmur_test

TEST(BordersMur, MurCxMatchesFortranCalcMurconstants310) {
    const std::string json =
        (std::filesystem::path("testData") / "cases" / "planewave" / "pw-in-box.fdtd.json").string();
    ASSERT_TRUE(std::filesystem::exists(json));
    const auto info = SEMBA_FDTD_m::SEMBA_FDTD_test::test_plane_wave_init(json, 0);
    const double expected = bordersmur_test::expectedMurCx(info.dt, 0.01);
    EXPECT_NEAR(info.murCx, expected, 1e-12);
    EXPECT_NEAR(info.murCy, expected, 1e-12);
    EXPECT_NEAR(info.murCz, expected, 1e-12);
}

TEST(BordersMur, FirstOrderBackHyFace_Fortran1107) {
    const std::string json =
        (std::filesystem::path("testData") / "cases" / "planewave" / "pw-in-box.fdtd.json").string();
    const auto info = SEMBA_FDTD_m::SEMBA_FDTD_test::test_plane_wave_init(json, 0);
    const double expected = 0.4 + info.murCx * (0.5 - 0.2);
    const double got = SEMBA_FDTD_m::SEMBA_FDTD_test::test_mur_apply_back_hy(json);
    EXPECT_NEAR(got, expected, 1e-7);
}

#endif
