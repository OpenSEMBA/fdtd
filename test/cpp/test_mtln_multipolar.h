#ifndef TEST_MTLN_MULTIPOLAR_H
#define TEST_MTLN_MULTIPOLAR_H

#include <gtest/gtest.h>
#include <cmath>
#include <vector>

#include "mtln_types.h"
#include "multipolar_expansion_m.h"
#include "test_mtln_helpers.h"

using mtln_types_m::RKIND;
using mtln_types_m::box_2d_t;
using mtln_types_m::multipolar_coefficient_t;
using mtln_types_m::multipolar_expansion_t;

TEST(mtln, multipolar_expansion_of_dipole) {
    int error_cnt = 0;

    const std::vector<RKIND> expansionCenter = {0.0, 0.0};
    const RKIND d = 0.1;
    const RKIND r = 1.0;
    std::vector<multipolar_coefficient_t> ab(2);
    ab[0].a = 0.0;
    ab[0].b = 0.0;
    ab[1].a = d;
    ab[1].b = 0.0;

    {
        const std::vector<RKIND> pos = {r, 0.0};
        const RKIND vComputed =
            multipolar_expansion_m::multipolarExpansion2D(pos, ab, expansionCenter);
        const RKIND vExpected = 1.0 / (2.0 * 3.14159265358979323846) * std::log((r + d / 2.0) / (r - d / 2.0));
        if (std::abs(vExpected - vComputed) > 1e-4) {
            error_cnt += 1;
        }
    }
    {
        const std::vector<RKIND> pos = {0.0, r};
        const RKIND vComputed =
            multipolar_expansion_m::multipolarExpansion2D(pos, ab, expansionCenter);
        if (std::abs(vComputed) > 1e-4) {
            error_cnt += 1;
        }
    }
    EXPECT_EQ(error_cnt, 0);
}

TEST(mtln, multipolar_expansion_lansink_two_wires) {
    multipolar_expansion_t mE;
    mE.inner_region.min = {-0.0265, -0.031};
    mE.inner_region.max = {0.0355, 0.031};

    mE.electric.resize(2);
    mE.electric[0].inner_region_average_potential = 0.56086362615993235;
    mE.electric[0].expansion_center = {-0.00497, 0.0};
    mE.electric[0].conductor_potentials = {1.0, 0.59092722039872780};
    mE.electric[0].ab = {{0.94888369862564237, 0.0},
                         {-4.5026111057062017e-07, 0.0},
                         {7.1226466480672610e-05, 7.4826684830590268e-09}};

    mE.electric[1].inner_region_average_potential = 0.80708482435726114;
    mE.electric[1].expansion_center = {0.0099205134400286565, 0.0};
    mE.electric[1].conductor_potentials = {0.84971105674469871, 1.0};
    mE.electric[1].ab = {{1.3644011168458479, 0.0},
                         {1.7915912102364171e-06, 0.0},
                         {1.4620553866347293e-06, -1.4363492844460606e-09}};
    mE.magnetic = mE.electric;

    box_2d_t fdtdCell;
    fdtdCell.min = {-0.100, -0.100};
    fdtdCell.max = {0.100, 0.100};

    int error_cnt = 0;
    const auto computedC = multipolar_expansion_m::getCellCapacitanceOnBox(mE, fdtdCell);
    if (!mtln_test::checkNear(14.08e-12, computedC[0][0], 6e-2)) {
        error_cnt += 1;
    }
    if (!mtln_test::checkNear(43.99e-12, computedC[0][1], 6e-2)) {
        error_cnt += 1;
    }
    if (!mtln_test::checkNear(44.31e-12, computedC[1][0], 6e-2)) {
        error_cnt += 1;
    }
    if (!mtln_test::checkNear(28.79e-12, computedC[1][1], 6e-2)) {
        error_cnt += 1;
    }

    const auto computedL = multipolar_expansion_m::getCellInductanceOnBox(mE, fdtdCell);
    if (!mtln_test::checkNear(791e-9, computedL[0][0], 6e-2)) {
        error_cnt += 1;
    }
    if (!mtln_test::checkNear(253e-9, computedL[0][1], 6e-2)) {
        error_cnt += 1;
    }
    if (!mtln_test::checkNear(251e-9, computedL[1][0], 6e-2)) {
        error_cnt += 1;
    }
    if (!mtln_test::checkNear(387e-9, computedL[1][1], 6e-2)) {
        error_cnt += 1;
    }
    EXPECT_EQ(error_cnt, 0);
}

TEST(mtln, multipolar_expansion_lansink_wire_w_dielectric) {
    multipolar_expansion_t mE;
    mE.inner_region.min = {-0.004, -0.004};
    mE.inner_region.max = {0.004, 0.004};

    mE.electric.resize(1);
    mE.electric[0].inner_region_average_potential = 0.90407844239490087;
    mE.electric[0].expansion_center = {0.0, 0.0};
    mE.electric[0].conductor_potentials = {1.0};
    mE.electric[0].ab = {{0.97667340898489752, 0.0}};

    mE.magnetic.resize(1);
    mE.magnetic[0].inner_region_average_potential = 0.84903792056711014;
    mE.magnetic[0].expansion_center = {0.0, 0.0};
    mE.magnetic[0].conductor_potentials = {1.0};
    mE.magnetic[0].ab = {{0.90929569352397666, 0.0}};

    box_2d_t fdtdCell;
    fdtdCell.min = {-0.0075, -0.0075};
    fdtdCell.max = {0.0075, 0.0075};

    int error_cnt = 0;
    const auto computedC = multipolar_expansion_m::getCellCapacitanceOnBox(mE, fdtdCell);
    if (!mtln_test::checkNear(49e-12, computedC[0][0], 6e-2)) {
        error_cnt += 1;
    }
    const auto computedL = multipolar_expansion_m::getCellInductanceOnBox(mE, fdtdCell);
    if (!mtln_test::checkNear(320e-9, computedL[0][0], 6e-2)) {
        error_cnt += 1;
    }
    EXPECT_EQ(error_cnt, 0);
}

#endif
