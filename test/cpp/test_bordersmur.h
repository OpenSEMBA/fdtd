#ifndef TEST_BORDERSMUR_H
#define TEST_BORDERSMUR_H

#include <gtest/gtest.h>

#include "semba_fdtd.h"

#include <cmath>
#include <fstream>
#include <filesystem>
#include <string>

namespace bordersmur_test {

constexpr double C0 = 299792458.0;

inline std::string pulse1dJson() {
    return (std::filesystem::path("testData") / "cases" / "mur" / "pulse-1d-x.fdtd.json").string();
}

inline double expectedMurCx(double dt, double dx) {
    return (C0 * dt - dx) / (C0 * dt + dx);
}

inline int pulseTransitSteps(int nx, double dx, double dt) {
    return static_cast<int>(std::ceil(2.5 * nx * dx / (C0 * dt)));
}

// Must match timeStep in testData/cases/mur/pulse-1d-x.fdtd.json
constexpr double kPulse1dDt = 1.54066656123717649e-11;

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
    EXPECT_NEAR(got, expected, 1e-12);
}

TEST(BordersMur, PulsePeakDecaysAfterTransit) {
    const std::string json = bordersmur_test::pulse1dJson();
    ASSERT_TRUE(std::filesystem::exists(json));
    const int steps = bordersmur_test::pulseTransitSteps(6, 0.01, bordersmur_test::kPulse1dDt);
    const auto open = SEMBA_FDTD_m::SEMBA_FDTD_test::test_mur_pulse_absorption(
        json, steps, 3, 3, 3, 1.0, false);
    const auto mur = SEMBA_FDTD_m::SEMBA_FDTD_test::test_mur_pulse_absorption(
        json, steps, 3, 3, 3, 1.0, true);
    EXPECT_GT(mur.max_ex_initial, 0.9);
    EXPECT_LT(mur.max_ex_final, 0.25 * open.max_ex_final)
        << "Mur max_ex_final=" << mur.max_ex_final << " open=" << open.max_ex_final
        << " after " << steps << " steps";
}

TEST(BordersMur, CenterProbeNearZeroAfterAbsorption) {
    const std::string json = bordersmur_test::pulse1dJson();
    ASSERT_TRUE(std::filesystem::exists(json));
    const int steps = bordersmur_test::pulseTransitSteps(6, 0.01, bordersmur_test::kPulse1dDt);
    const auto open = SEMBA_FDTD_m::SEMBA_FDTD_test::test_mur_pulse_absorption(
        json, steps, 3, 3, 3, 1.0, false);
    const auto mur = SEMBA_FDTD_m::SEMBA_FDTD_test::test_mur_pulse_absorption(
        json, steps, 3, 3, 3, 1.0, true);
    EXPECT_LE(std::abs(mur.probe_ex_final), std::abs(open.probe_ex_final) + 1e-12)
        << "center Ex mur=" << mur.probe_ex_final << " open=" << open.probe_ex_final;
    EXPECT_LT(mur.max_ex_final, open.max_ex_final);
}

TEST(BordersMur, EnergyDecreasesAfterPeak) {
    const std::string json = bordersmur_test::pulse1dJson();
    ASSERT_TRUE(std::filesystem::exists(json));
    const int steps = bordersmur_test::pulseTransitSteps(6, 0.01, bordersmur_test::kPulse1dDt);
    const auto open = SEMBA_FDTD_m::SEMBA_FDTD_test::test_mur_pulse_absorption(
        json, steps, 3, 3, 3, 1.0, false);
    const auto mur = SEMBA_FDTD_m::SEMBA_FDTD_test::test_mur_pulse_absorption(
        json, steps, 3, 3, 3, 1.0, true);
    EXPECT_LT(mur.energy_final, open.energy_final);
}


TEST(BordersMur, MurReducesPeakVersusOpenBoundary) {
    const std::string json = bordersmur_test::pulse1dJson();
    const int steps = bordersmur_test::pulseTransitSteps(6, 0.01, bordersmur_test::kPulse1dDt);
    const auto open = SEMBA_FDTD_m::SEMBA_FDTD_test::test_mur_pulse_absorption(
        json, steps, 3, 3, 3, 1.0, false);
    const auto mur = SEMBA_FDTD_m::SEMBA_FDTD_test::test_mur_pulse_absorption(
        json, steps, 3, 3, 3, 1.0, true);
    EXPECT_LT(mur.max_ex_final, open.max_ex_final)
        << "Mur maxEx=" << mur.max_ex_final << " open maxEx=" << open.max_ex_final;
    EXPECT_LT(mur.energy_final, open.energy_final);
}

TEST(BordersMur, PulseMatchesFortranProbe) {
    const std::filesystem::path golden_dir =
        std::filesystem::path("testData") / "cases" / "mur";
    const std::filesystem::path golden =
        golden_dir / "pulse-1d-x_probe_Ex_20_1_1.dat";
    if (!std::filesystem::exists(golden)) {
        GTEST_SKIP() << "Fortran golden probe missing: " << golden;
    }
    const std::string json = bordersmur_test::pulse1dJson();
    const int steps = 20;
    const auto r = SEMBA_FDTD_m::SEMBA_FDTD_test::test_mur_pulse_absorption(
        json, steps, 3, 3, 3, 1.0, true);
    std::ifstream in(golden);
    ASSERT_TRUE(in.is_open());
    std::string header;
    std::getline(in, header);
    double t_ref = 0.0, ex_ref = 0.0;
    in >> t_ref >> ex_ref;
    in.close();
    (void)t_ref;
    EXPECT_NEAR(r.probe_ex_final, ex_ref, 5e-4)
        << "C++ probe Ex=" << r.probe_ex_final << " Fortran=" << ex_ref << " at step " << steps;
}

#endif
