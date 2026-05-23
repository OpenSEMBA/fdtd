#ifndef TEST_MTLN_MTL_H
#define TEST_MTLN_MTL_H

#include <gtest/gtest.h>
#include <algorithm>

#include "mtl_m.h"
#include "mtln_types.h"
#include "test_mtln_helpers.h"

using mtln_test::buildLineWithNConductors;
using mtln_test::buildTestBundle;
using mtln_test::checkNear;
using mtln_test::comparePULMatrices;
using mtln_test::emptyTransferImpedance;

TEST(mtln, mtl_wrong_dt) {
    mtl_m::mtl_t line = buildLineWithNConductors(2, "line0", mtln_test::MTL_TYPE_SHIELDED, 1.0);
    EXPECT_NE(line.dt, 1.0);
}

TEST(mtln, mtl_homogeneous) {
    using namespace mtln_types_m;
    using namespace mtl_m;

    const std::vector<std::vector<double>> lpul = {
        {4.4712610E-07, 1.4863653E-07}, {1.4863653E-07, 4.4712610E-07}};
    const std::vector<std::vector<double>> cpul = {
        {2.242e-10, -7.453e-11}, {-7.453e-11, 2.242e-10}};
    const std::vector<std::vector<double>> rpul(2, std::vector<double>(2, 0.0));
    const std::vector<std::vector<double>> gpul(2, std::vector<double>(2, 0.0));
    std::vector<double> step_size(5, 20.0);
    std::vector<segment_t> segments(5);
    for (int i = 0; i < 5; ++i) {
        segments[i].x = i + 1;
        segments[i].y = 1;
        segments[i].z = 1;
        segments[i].orientation = DIRECTION_X_POS;
    }

    int error_cnt = 0;
    const auto zt = emptyTransferImpedance();
    std::vector<multipolar_expansion_t> mE;

    mtl_t line = mtl_shielded(lpul, cpul, rpul, gpul, step_size, "line0", segments, 1e-12, "p", 1,
                              zt);
    comparePULMatrices(error_cnt, line.lpul, lpul);
    comparePULMatrices(error_cnt, line.cpul, cpul);
    comparePULMatrices(error_cnt, line.rpul, rpul);
    comparePULMatrices(error_cnt, line.gpul, gpul);

    line = mtl_unshielded(lpul, cpul, rpul, gpul, step_size, "line0", segments, 1e-12, mE, 0.0);
    comparePULMatrices(error_cnt, line.lpul, lpul);
    comparePULMatrices(error_cnt, line.cpul, cpul);
    comparePULMatrices(error_cnt, line.rpul, rpul);
    comparePULMatrices(error_cnt, line.gpul, gpul);

    EXPECT_EQ(error_cnt, 0);
}

TEST(mtln, mtl_time_step) {
    mtl_m::mtl_t line =
        buildLineWithNConductors(2, "line0", mtln_test::MTL_TYPE_UNSHIELDED, 1e-6);

    const auto phase_velocities = line.getPhaseVelocities();
    const double max_vel = phase_velocities.empty() || phase_velocities[0].empty()
                               ? 0.0
                               : *std::max_element(phase_velocities[0].begin(),
                                                   phase_velocities[0].end());
    const double time_step = line.getMaxTimeStep();

    int error_cnt = 0;
    if (!checkNear(phase_velocities[0][0], 1.05900008e8, 0.01)) {
        error_cnt += 1;
    }
    if (!checkNear(phase_velocities[0][1], 1.05900010e8, 0.01)) {
        error_cnt += 1;
    }
    if (!checkNear(time_step, 1.888573951383424e-07, 0.01)) {
        error_cnt += 1;
    }
    (void)max_vel;
    EXPECT_EQ(error_cnt, 0);
}

TEST(mtln, mtl_bundle_init) {
    using namespace mtln_types_m;
    using namespace mtl_m;

    std::vector<double> step_size(5, 20.0);
    std::vector<segment_t> segments(5);
    for (int i = 0; i < 5; ++i) {
        segments[i].x = i + 1;
        segments[i].y = 1;
        segments[i].z = 1;
        segments[i].orientation = 1;
    }

    const std::vector<std::vector<double>> l1 = {{4.4712610E-07}};
    const std::vector<std::vector<double>> c1 = {{2.242e-10}};
    const std::vector<std::vector<double>> r1 = {{0.0}};
    const std::vector<std::vector<double>> g1 = {{0.0}};
    const auto zt = emptyTransferImpedance();
    std::vector<multipolar_expansion_t> mE;

    mtl_t mtl_in = mtl_shielded(l1, c1, r1, g1, step_size, "line_in", segments, 1e-11,
                                "line_out", 1, zt);
    mtl_t mtl_out = mtl_unshielded(l1, c1, r1, g1, step_size, "line_out", segments, 1e-11, mE, 0.0);

    std::vector<transmission_line_level_t> levels(2);
    levels[0].lines = {mtl_out};
    levels[1].lines = {mtl_in};

    auto bundle = buildTestBundle(levels, "bundle");
    int error_cnt = 0;
    if (bundle.number_of_divisions != 5 || bundle.number_of_conductors != 2) {
        error_cnt += 1;
    }
    const int nc = bundle.number_of_conductors;
    if (mtln_test::test_mtl_bundle_t::at3(bundle.lpul, 0, 0, 0, nc, nc) != mtl_out.lpul[0][0][0]) {
        error_cnt += 1;
    }
    if (mtln_test::test_mtl_bundle_t::at3(bundle.lpul, 0, 1, 1, nc, nc) != mtl_in.lpul[0][0][0]) {
        error_cnt += 1;
    }
    EXPECT_EQ(error_cnt, 0);
}

#endif
