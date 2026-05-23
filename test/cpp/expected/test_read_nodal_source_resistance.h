#ifndef TEST_READ_NODAL_SOURCE_RESISTANCE_H
#define TEST_READ_NODAL_SOURCE_RESISTANCE_H

#include <algorithm>
#include <cmath>
#include "test_smbjson_helpers.h"

// Partial test only — no full Parseador expected builder.

inline constexpr double EXPECTED_CABLE_RESISTANCE_PER_METER = 10000.0;

inline std::string nodalSourceResistanceJsonPath() {
    return TEST_DATA_DIR + "cases/nodalSource/nodalSource.fdtd.json";
}

inline std::string nodalSourceTotalResistanceJsonPath() {
    return TEST_DATA_DIR + "cases/nodalSource/nodalSource_totalResistance.fdtd.json";
}

inline void expectCableResistancePerMeter(const Parseador_t& problem,
                                          double expectedResistance) {
    const double tol = std::max(std::abs(expectedResistance) * 1.0e-6, 1.0e-6);

#ifndef CompileWithMTLN
    EXPECT_EQ(problem.tWires->n_tw, 1);
    EXPECT_EQ(problem.tWires->tw[0].n_twc, 10);
    EXPECT_NEAR(problem.tWires->tw[0].res, expectedResistance, tol);
#else
    ASSERT_TRUE(problem.mtln);
    ASSERT_EQ(problem.mtln->cables.size(), 1u);
    auto* cable = dynamic_cast<mtln_types_m::unshielded_multiwire_t*>(
        problem.mtln->cables[0].ptr.get());
    ASSERT_NE(cable, nullptr) << "Expected an unshielded multiwire cable";
    EXPECT_EQ(cable->step_size.size(), 10u);
    ASSERT_FALSE(cable->resistance_per_meter.empty());
    ASSERT_FALSE(cable->resistance_per_meter[0].empty());
    EXPECT_NEAR(cable->resistance_per_meter[0][0], expectedResistance, tol);
#endif
}

inline void expectTotalResistanceOverride(const Parseador_t& problem,
                                          double expectedResistance) {
    const double tol = std::max(std::abs(expectedResistance) * 1.0e-6, 1.0e-6);

#ifndef CompileWithMTLN
    EXPECT_EQ(problem.tWires->n_tw, 1);
    EXPECT_EQ(problem.tWires->tw[0].n_twc, 10);
    EXPECT_NEAR(problem.tWires->tw[0].res, expectedResistance, tol);
#else
    ASSERT_TRUE(problem.mtln);
    ASSERT_EQ(problem.mtln->cables.size(), 1u);
    auto* cable = dynamic_cast<mtln_types_m::unshielded_multiwire_t*>(
        problem.mtln->cables[0].ptr.get());
    ASSERT_NE(cable, nullptr) << "Expected an unshielded multiwire cable";
    EXPECT_EQ(cable->step_size.size(), 10u);
    ASSERT_FALSE(cable->resistance_per_meter.empty());
    ASSERT_FALSE(cable->resistance_per_meter[0].empty());
    EXPECT_NEAR(cable->resistance_per_meter[0][0], expectedResistance, tol);
#endif
}

#endif // TEST_READ_NODAL_SOURCE_RESISTANCE_H
