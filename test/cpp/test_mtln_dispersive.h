#ifndef TEST_MTLN_DISPERSIVE_H
#define TEST_MTLN_DISPERSIVE_H

#include <gtest/gtest.h>
#include <complex>

#include "dispersive_m.h"
#include "mtl_bundle_m.h"
#include "mtln_types.h"
#include "test_mtln_helpers.h"

namespace {

using mtln_types_m::RKIND;
using mtln_types_m::TRANSFER_IMPEDANCE_DIRECTION_BOTH;
using mtln_types_m::transfer_impedance_per_meter_t;

transfer_impedance_per_meter_t makeTestTransferImpedance() {
    transfer_impedance_per_meter_t zt;
    zt.direction = TRANSFER_IMPEDANCE_DIRECTION_BOTH;
    zt.resistive_term = 1e-2;
    zt.inductive_term = 1e-6;
    zt.poles = {std::complex<RKIND>(-1e6, 1e1)};
    zt.residues = {std::complex<RKIND>(1e4, 1e3)};
    return zt;
}

bool allNonZero(const std::vector<std::vector<std::vector<std::vector<std::complex<RKIND>>>>>& q,
                int div, int c1, int c2) {
    for (const auto& val : q[static_cast<size_t>(div)][static_cast<size_t>(c1)][static_cast<size_t>(c2)]) {
        if (val == std::complex<RKIND>{0.0, 0.0}) {
            return false;
        }
    }
    return true;
}

bool allZero(const std::vector<std::vector<std::vector<std::vector<std::complex<RKIND>>>>>& q,
             int div, int c1, int c2) {
    for (const auto& val : q[static_cast<size_t>(div)][static_cast<size_t>(c1)][static_cast<size_t>(c2)]) {
        if (val != std::complex<RKIND>{0.0, 0.0}) {
            return false;
        }
    }
    return true;
}

void checkDispersiveSizes(const dispersive_m::transfer_impedance_t& ti, int ndiv, int nc, int npoles,
                          int& error_cnt) {
    if (dispersive_m::flatSize4(ti.q1) != static_cast<size_t>(ndiv * nc * nc * npoles)) {
        error_cnt += 1;
    }
    if (dispersive_m::flatSize4(ti.q2) != static_cast<size_t>(ndiv * nc * nc * npoles)) {
        error_cnt += 1;
    }
    if (dispersive_m::flatSize4(ti.q3) != static_cast<size_t>(ndiv * nc * nc * npoles)) {
        error_cnt += 1;
    }
    if (dispersive_m::flatSize3(ti.phi) != static_cast<size_t>(ndiv * nc * npoles)) {
        error_cnt += 1;
    }
    if (dispersive_m::flatSize3Real(ti.d) != static_cast<size_t>(ndiv * nc * nc)) {
        error_cnt += 1;
    }
    if (dispersive_m::flatSize3Real(ti.e) != static_cast<size_t>(ndiv * nc * nc)) {
        error_cnt += 1;
    }
    if (dispersive_m::flatSize3(ti.q1_sum) != static_cast<size_t>(ndiv * nc * nc)) {
        error_cnt += 1;
    }
    if (dispersive_m::flatSize3(ti.q2_sum) != static_cast<size_t>(ndiv * nc * nc)) {
        error_cnt += 1;
    }
    if (dispersive_m::flatSize2(ti.q3_phi) != static_cast<size_t>(ndiv * nc)) {
        error_cnt += 1;
    }
}

void checkCouplingPattern(const dispersive_m::transfer_impedance_t& ti, int& error_cnt) {
    for (int div = 0; div < ti.number_of_divisions; ++div) {
        if (!allNonZero(ti.q1, div, 0, 1)) {
            error_cnt += 1;
        }
        if (!allNonZero(ti.q1, div, 1, 0)) {
            error_cnt += 1;
        }
        if (!allZero(ti.q1, div, 1, 1)) {
            error_cnt += 1;
        }
        if (!allZero(ti.q1, div, 0, 0)) {
            error_cnt += 1;
        }
    }
}

mtl_bundle_m::mtl_bundle_t makeTwoLevelBundle(const transfer_impedance_per_meter_t& zt) {
    mtl_m::transmission_line_level_t level1;
    mtl_m::transmission_line_level_t level2;
    level1.lines.push_back(mtln_test::buildLineWithNConductors(1, "line_out", mtln_test::MTL_TYPE_UNSHIELDED));
    level2.lines.push_back(mtln_test::buildLineWithNConductors(
        1, "line_in", mtln_test::MTL_TYPE_SHIELDED, std::nullopt, "line_out", 1));
    auto bundle = mtl_bundle_m::mtl_bundle_ctor({level1, level2}, "bundle");
    bundle.addTransferImpedance(1, {2}, zt);
    return bundle;
}

} // namespace

TEST(mtln, dispersive_init_1_pole) {
    int error_cnt = 0;
    const auto zt = makeTestTransferImpedance();
    const auto bundle = makeTwoLevelBundle(zt);
    checkDispersiveSizes(bundle.transfer_impedance, 5, 2, 1, error_cnt);
    checkCouplingPattern(bundle.transfer_impedance, error_cnt);
    EXPECT_EQ(error_cnt, 0);
}

TEST(mtln, dispersive_init_2_poles) {
    int error_cnt = 0;
    auto zt = makeTestTransferImpedance();
    zt.poles.push_back(std::complex<RKIND>(-1e6, -1e1));
    zt.residues.push_back(std::complex<RKIND>(1e4, -1e3));
    const auto bundle = makeTwoLevelBundle(zt);
    checkDispersiveSizes(bundle.transfer_impedance, 5, 2, 2, error_cnt);
    EXPECT_EQ(error_cnt, 0);
}

TEST(mtln, dispersive_init_1_pole_3_levels) {
    int error_cnt = 0;
    const auto zt = makeTestTransferImpedance();
    mtl_m::transmission_line_level_t level1;
    mtl_m::transmission_line_level_t level2;
    mtl_m::transmission_line_level_t level3;
    level1.lines.push_back(mtln_test::buildLineWithNConductors(1, "line1", mtln_test::MTL_TYPE_UNSHIELDED));
    level2.lines.push_back(mtln_test::buildLineWithNConductors(
        2, "line2", mtln_test::MTL_TYPE_SHIELDED, std::nullopt, "line1", 1));
    level3.lines.push_back(mtln_test::buildLineWithNConductors(
        1, "line3", mtln_test::MTL_TYPE_SHIELDED, std::nullopt, "line2", 2));
    auto bundle = mtl_bundle_m::mtl_bundle_ctor({level1, level2, level3}, "bundle");
    bundle.addTransferImpedance(2, {3}, zt);
    checkDispersiveSizes(bundle.transfer_impedance, 5, 4, 1, error_cnt);
    EXPECT_EQ(error_cnt, 0);
}

TEST(mtln, dispersive_init_1_pole_lines_with_lumped) {
    int error_cnt = 0;
    const auto zt = makeTestTransferImpedance();
    mtl_m::transmission_line_level_t level1;
    mtl_m::transmission_line_level_t level2;
    auto line_out = mtln_test::buildLineWithNConductors(1, "line_out", mtln_test::MTL_TYPE_UNSHIELDED);
    auto line_in = mtln_test::buildLineWithNConductors(
        1, "line_in", mtln_test::MTL_TYPE_SHIELDED, std::nullopt, "line_out", 1);
    line_out.lumped_elements.addDispersiveLumped(1, 1, zt);
    line_in.lumped_elements.addDispersiveLumped(5, 1, zt);
    level1.lines.push_back(std::move(line_out));
    level2.lines.push_back(std::move(line_in));
    auto bundle = mtl_bundle_m::mtl_bundle_ctor({level1, level2}, "bundle");
    bundle.addTransferImpedance(1, {2}, zt);
    checkDispersiveSizes(bundle.transfer_impedance, 5, 2, 1, error_cnt);
    for (int div = 0; div < bundle.transfer_impedance.number_of_divisions; ++div) {
        if (!allNonZero(bundle.transfer_impedance.q1, div, 0, 1)) {
            error_cnt += 1;
        }
        if (!allNonZero(bundle.transfer_impedance.q1, div, 1, 0)) {
            error_cnt += 1;
        }
    }
    for (int div = 2; div <= 4; ++div) {
        if (!allZero(bundle.transfer_impedance.q1, div - 1, 1, 1)) {
            error_cnt += 1;
        }
        if (!allZero(bundle.transfer_impedance.q1, div - 1, 0, 0)) {
            error_cnt += 1;
        }
    }
    if (!allNonZero(bundle.transfer_impedance.q1, 4, 1, 1)) {
        error_cnt += 1;
    }
    if (!allNonZero(bundle.transfer_impedance.q1, 0, 0, 0)) {
        error_cnt += 1;
    }
    if (!allZero(bundle.transfer_impedance.q1, 4, 0, 0)) {
        error_cnt += 1;
    }
    if (!allZero(bundle.transfer_impedance.q1, 0, 1, 1)) {
        error_cnt += 1;
    }
    EXPECT_EQ(error_cnt, 0);
}

#endif
