#ifndef TEST_PLANEWAVE_EVOLUCION_H
#define TEST_PLANEWAVE_EVOLUCION_H

#include <gtest/gtest.h>

#include "semba_fdtd.h"

#include <cmath>
#include <filesystem>
#include <string>

namespace planewave_test {

#ifdef CompileWithReal8
using evol_real = double;
#else
using evol_real = float;
#endif

inline std::string pwInBoxJson() {
    return (std::filesystem::path("testData") / "cases" / "planewave" / "pw-in-box.fdtd.json").string();
}

inline evol_real r(double value) {
    return static_cast<evol_real>(value);
}

// Golden from Fortran planewaves.F90 evolucion L793-798 with gauss_1GHz.exc samples.
constexpr double kDeltaEvol = 1.805468626449816074e-13;
constexpr int kNormalSampleNprev = 1350;
constexpr double kSample1350 = 1.425732537018645291e-33;
constexpr double kSample1351 = 1.449875611218499565e-33;

} // namespace planewave_test

// Fortran planewaves.F90 L793-798: nprev=0 returns 0.
TEST(PlanewaveEvolucion, ZeroWhenNprevBelowOne_Fortran793) {
    const std::string json = planewave_test::pwInBoxJson();
    ASSERT_TRUE(std::filesystem::exists(json));
    const double v = SEMBA_FDTD_m::SEMBA_FDTD_test::test_evolucion(json, 0, 0.5 * planewave_test::kDeltaEvol);
    EXPECT_DOUBLE_EQ(v, 0.0);
}

// Midpoint of a normal-valued interpolation interval (nprev=1350).
TEST(PlanewaveEvolucion, MatchesFortranLines793_798_MidInterval) {
    const std::string json = planewave_test::pwInBoxJson();
    ASSERT_TRUE(std::filesystem::exists(json));
    const auto delta = planewave_test::r(planewave_test::kDeltaEvol);
    const auto nprev = planewave_test::r(planewave_test::kNormalSampleNprev);
    const auto t_delay = nprev * delta + planewave_test::r(0.5) * delta;
    const auto y0 = planewave_test::r(planewave_test::kSample1350);
    const auto y1 = planewave_test::r(planewave_test::kSample1351);
    const auto expected = ((y1 - y0) / delta) * (t_delay - nprev * delta) + y0;
    const double v = SEMBA_FDTD_m::SEMBA_FDTD_test::test_evolucion(json, 0, t_delay);
    EXPECT_EQ(v, static_cast<double>(expected));
}

// Manual Fortran formula using committed sample literals.
TEST(PlanewaveEvolucion, MatchesManualFortranFormula) {
    const auto delta = planewave_test::r(planewave_test::kDeltaEvol);
    const auto nprev = planewave_test::r(planewave_test::kNormalSampleNprev);
    const auto t_delay = nprev * delta + planewave_test::r(0.5) * delta;
    const auto t_frac = t_delay - nprev * delta;
    const auto expected = planewave_test::r(planewave_test::kSample1350) +
        (planewave_test::r(planewave_test::kSample1351) - planewave_test::r(planewave_test::kSample1350)) *
        (t_frac / delta);
    const std::string json = planewave_test::pwInBoxJson();
    const double v = SEMBA_FDTD_m::SEMBA_FDTD_test::test_evolucion(json, 0, t_delay);
    EXPECT_EQ(v, static_cast<double>(expected));
}

// Beyond table: nprev+1 > numSamples -> 0.
TEST(PlanewaveEvolucion, ZeroAfterExcitationTableEnd) {
    const std::string json = planewave_test::pwInBoxJson();
    const auto info = SEMBA_FDTD_m::SEMBA_FDTD_test::test_plane_wave_init(json, 0);
    ASSERT_GT(info.numSamples, 2);
    const double t_delay = static_cast<double>(info.numSamples + 10) * info.deltaevol;
    const double v = SEMBA_FDTD_m::SEMBA_FDTD_test::test_evolucion(json, 0, t_delay);
    EXPECT_DOUBLE_EQ(v, 0.0);
}

#endif
