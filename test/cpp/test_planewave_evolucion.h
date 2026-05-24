#ifndef TEST_PLANEWAVE_EVOLUCION_H
#define TEST_PLANEWAVE_EVOLUCION_H

#include <gtest/gtest.h>

#include "semba_fdtd.h"

#include <cmath>
#include <filesystem>
#include <string>

namespace planewave_test {

inline std::string pwInBoxJson() {
    return (std::filesystem::path("testData") / "cases" / "planewave" / "pw-in-box.fdtd.json").string();
}

// Golden from Fortran planewaves.F90 evolucion L793-798 with gauss_1GHz.exc samples.
constexpr double kDeltaEvol = 1.805468626449816074e-13;
constexpr double kSample0 = 3.720075976020888891e-44;
constexpr double kSample1 = 3.792604493180529565e-44;
constexpr double kEvolMidNprev1 = 3.756340234600708984e-44;

} // namespace planewave_test

// Fortran planewaves.F90 L793-798: nprev=0 returns 0.
TEST(PlanewaveEvolucion, ZeroWhenNprevBelowOne_Fortran793) {
    const std::string json = planewave_test::pwInBoxJson();
    ASSERT_TRUE(std::filesystem::exists(json));
    const double v = SEMBA_FDTD_m::SEMBA_FDTD_test::test_evolucion(json, 0, 0.5 * planewave_test::kDeltaEvol);
    EXPECT_DOUBLE_EQ(v, 0.0);
}

// Midpoint of first interpolation interval (nprev=1).
TEST(PlanewaveEvolucion, MatchesFortranLines793_798_MidInterval) {
    const std::string json = planewave_test::pwInBoxJson();
    ASSERT_TRUE(std::filesystem::exists(json));
    const double t_delay = 1.5 * planewave_test::kDeltaEvol;
    const double v = SEMBA_FDTD_m::SEMBA_FDTD_test::test_evolucion(json, 0, t_delay);
    EXPECT_NEAR(v, planewave_test::kEvolMidNprev1, 1e-55);
}

// Manual Fortran formula using committed sample literals.
TEST(PlanewaveEvolucion, MatchesManualFortranFormula) {
    const int nprev = 1;
    const double t_frac = 0.5 * planewave_test::kDeltaEvol;
    const double expected = planewave_test::kSample0 +
        (planewave_test::kSample1 - planewave_test::kSample0) * (t_frac / planewave_test::kDeltaEvol);
    const std::string json = planewave_test::pwInBoxJson();
    const double v = SEMBA_FDTD_m::SEMBA_FDTD_test::test_evolucion(json, 0, planewave_test::kDeltaEvol + t_frac);
    EXPECT_NEAR(v, expected, 1e-55);
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
