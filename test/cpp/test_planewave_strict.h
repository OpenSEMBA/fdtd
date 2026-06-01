#ifndef TEST_PLANEWAVE_STRICT_H
#define TEST_PLANEWAVE_STRICT_H

#include <gtest/gtest.h>

#include "semba_fdtd.h"

#include <cmath>
#include <filesystem>
#include <limits>
#include <string>

namespace planewave_strict_test {

#ifdef CompileWithReal8
using strict_real = double;
#else
using strict_real = float;
#endif

inline std::string pwInBoxJson() {
    return (std::filesystem::path("testData") / "cases" / "planewave" / "pw-in-box.fdtd.json").string();
}

inline std::string pwInBoxPecJson() {
    return (std::filesystem::path("testData") / "cases" / "planewave" / "pw-in-box-pec.fdtd.json").string();
}

inline strict_real r(double value) {
    return static_cast<strict_real>(value);
}

constexpr double kFortranCflDt = 1.5406665526684904e-11;
constexpr double kDeltaEvol = 1.805468626449816074e-13;
constexpr int kNormalSampleNprev = 1350;
constexpr double kSample1350 = 1.425732537018645291e-33;
constexpr double kSample1351 = 1.449875611218499565e-33;

} // namespace planewave_strict_test

TEST(PlanewaveStrict, InitUsesFortranRealKindForCflAndExcitationStep) {
    const std::string json = planewave_strict_test::pwInBoxJson();
    ASSERT_TRUE(std::filesystem::exists(json));
    const auto info = SEMBA_FDTD_m::SEMBA_FDTD_test::test_plane_wave_init(json, 0);

    EXPECT_EQ(info.dt, static_cast<double>(planewave_strict_test::r(planewave_strict_test::kFortranCflDt)));
    EXPECT_EQ(info.deltaevol, static_cast<double>(planewave_strict_test::r(planewave_strict_test::kDeltaEvol)));
    EXPECT_EQ(info.numSamples, 31072);
}

TEST(PlanewaveStrict, SourceBoxUsesFortranSemiOpenIntervals) {
    const std::string json = planewave_strict_test::pwInBoxJson();
    ASSERT_TRUE(std::filesystem::exists(json));
    const auto info = SEMBA_FDTD_m::SEMBA_FDTD_test::test_plane_wave_init(json, 0);

    EXPECT_EQ(info.esqx1, 2);
    EXPECT_EQ(info.esqx2, 4);
    EXPECT_EQ(info.esqy1, 2);
    EXPECT_EQ(info.esqy2, 4);
    EXPECT_EQ(info.esqz1, 2);
    EXPECT_EQ(info.esqz2, 4);
}

TEST(PlanewaveStrict, EvolucionUsesFortranFirstOrderOperationOrder) {
    const std::string json = planewave_strict_test::pwInBoxJson();
    ASSERT_TRUE(std::filesystem::exists(json));

    const auto delta = planewave_strict_test::r(planewave_strict_test::kDeltaEvol);
    const auto nprev = planewave_strict_test::r(planewave_strict_test::kNormalSampleNprev);
    const auto y0 = planewave_strict_test::r(planewave_strict_test::kSample1350);
    const auto y1 = planewave_strict_test::r(planewave_strict_test::kSample1351);
    const auto t_delay = nprev * delta + planewave_strict_test::r(0.5) * delta;
    const auto t_frac = t_delay - nprev * delta;
    const auto expected = ((y1 - y0) / delta) * t_frac + y0;
    const double got = SEMBA_FDTD_m::SEMBA_FDTD_test::test_evolucion(
        json, 0, static_cast<double>(t_delay));

    EXPECT_EQ(got, static_cast<double>(expected));
}

TEST(PlanewaveStrict, IncidentFieldFlushesSinglePrecisionSubnormals) {
    const std::string json = planewave_strict_test::pwInBoxPecJson();
    ASSERT_TRUE(std::filesystem::exists(json));
    const auto info = SEMBA_FDTD_m::SEMBA_FDTD_test::test_plane_wave_init(json, 0);

    const double inc12 = SEMBA_FDTD_m::SEMBA_FDTD_test::test_compute_incid(
        json, 0, 4, 12.0 * info.dt, 3, 3, 1);
    const double inc13 = SEMBA_FDTD_m::SEMBA_FDTD_test::test_compute_incid(
        json, 0, 4, 13.0 * info.dt, 3, 3, 1);

#ifdef CompileWithReal8
    EXPECT_GT(std::abs(inc12), 0.0);
#else
    EXPECT_EQ(inc12, 0.0);
#endif
    EXPECT_GE(std::abs(inc13), static_cast<double>(std::numeric_limits<planewave_strict_test::strict_real>::min()));
}

TEST(PlanewaveStrict, IncidentFieldUsesFortranPlanewaveLightSpeed) {
    const std::string json = planewave_strict_test::pwInBoxPecJson();
    ASSERT_TRUE(std::filesystem::exists(json));
    const auto info = SEMBA_FDTD_m::SEMBA_FDTD_test::test_plane_wave_init(json, 0);

    const double inc = SEMBA_FDTD_m::SEMBA_FDTD_test::test_compute_incid(
        json, 0, 0, 15.5 * info.dt, 3, 3, 2);

#ifdef CompileWithReal8
    EXPECT_GT(std::abs(inc), 0.0);
#else
    EXPECT_EQ(inc, static_cast<double>(planewave_strict_test::r(3.8821768523548882e-35)));
#endif
}

TEST(PlanewaveStrict, GridInverseMatchesFortranSinglePrecisionArrayValue) {
    const std::string json = planewave_strict_test::pwInBoxPecJson();
    ASSERT_TRUE(std::filesystem::exists(json));
    const double idzh = SEMBA_FDTD_m::SEMBA_FDTD_test::test_grid_inverse_z(json, 3);

#ifdef CompileWithReal8
    EXPECT_EQ(idzh, 100.0);
#else
    EXPECT_EQ(idzh, static_cast<double>(planewave_strict_test::r(99.999992370605469)));
#endif
}

#endif
