#ifndef TEST_MTLN_HELPERS_H
#define TEST_MTLN_HELPERS_H

#include <gtest/gtest.h>
#include <cmath>
#include <optional>
#include <string>
#include <vector>

#include "mtln_types.h"
#include "mtl_m.h"

namespace mtln_test {

static const std::string PATH_TO_TEST_DATA = "testData/";
static const std::string MTL_TYPE_SHIELDED = "shielded";
static const std::string MTL_TYPE_UNSHIELDED = "unshielded";

inline bool checkNear(double target, double number, double rel_tol) {
    const double abs_diff = std::abs(target - number);
    if (abs_diff == 0.0) {
        return true;
    }
    return std::abs(target - number) / target < rel_tol;
}

inline bool checkNearTime(double target, double number, double rel_tol) {
    return checkNear(target, number, rel_tol);
}

inline void comparePULMatrices(int& error_cnt,
                               const std::vector<std::vector<std::vector<double>>>& m_line,
                               const std::vector<std::vector<double>>& m_input) {
    if (m_input.empty() || m_input.size() != m_input[0].size()) {
        error_cnt += 1;
        return;
    }
    for (const auto& slice : m_line) {
        if (slice != m_input) {
            error_cnt += 1;
        }
    }
}

inline mtln_types_m::transfer_impedance_per_meter_t emptyTransferImpedance() {
    mtln_types_m::transfer_impedance_per_meter_t zt;
    zt.inductive_term = 0.0;
    zt.resistive_term = 0.0;
    return zt;
}

inline mtl_m::mtl_t buildLineWithNConductors(
    int n,
    const std::string& name,
    const std::string& type,
    std::optional<double> dt = std::nullopt,
    std::optional<std::string> parent_name = std::nullopt,
    std::optional<int> conductor_in_parent = std::nullopt) {
    using namespace mtln_types_m;
    using namespace mtl_m;

    std::vector<std::vector<double>> lpul(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> cpul(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> rpul(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> gpul(n, std::vector<double>(n, 0.0));
    std::vector<double> step_size(5, 20.0);
    std::vector<segment_t> segments(5);

    for (int i = 0; i < 5; ++i) {
        segments[i].x = 1;
        segments[i].y = i + 1;
        segments[i].z = i + 1;
        segments[i].orientation = DIRECTION_X_POS;
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            rpul[i][j] = 0.0;
            gpul[i][j] = 0.0;
            if (i == j) {
                lpul[i][j] = 4.4712610E-07;
                cpul[i][j] = 2.242e-10;
            } else {
                lpul[i][j] = 1.4863653E-07;
                cpul[i][j] = -7.453e-11;
            }
        }
    }

    const double time_step = dt.value_or(1e-12);
    const transfer_impedance_per_meter_t zt = emptyTransferImpedance();
    std::vector<multipolar_expansion_t> mE;

    if (type == MTL_TYPE_SHIELDED) {
        const std::string parent = parent_name.value_or("p");
        const int conductor = conductor_in_parent.value_or(-1);
        return mtl_shielded(lpul, cpul, rpul, gpul, step_size, name, segments, time_step, parent,
                            conductor, zt);
    }
    if (type == MTL_TYPE_UNSHIELDED) {
        return mtl_unshielded(lpul, cpul, rpul, gpul, step_size, name, segments, time_step, mE, 0.0);
    }
    ADD_FAILURE() << "Unrecognized line type: " << type;
    return mtl_t{};
}

inline mtl_m::transmission_line_bundle_t buildFourLevelBundle() {
    mtl_m::transmission_line_bundle_t bundle;
    bundle.levels.resize(4);
    bundle.levels[0].lines = {buildLineWithNConductors(1, "line1", MTL_TYPE_UNSHIELDED)};
    bundle.levels[1].lines = {buildLineWithNConductors(
        2, "line2", MTL_TYPE_SHIELDED, std::nullopt, "line1", 1)};
    bundle.levels[2].lines = {
        buildLineWithNConductors(2, "line3_1", MTL_TYPE_SHIELDED, std::nullopt, "line2", 1),
        buildLineWithNConductors(2, "line3_2", MTL_TYPE_SHIELDED, std::nullopt, "line2", 2),
    };
    bundle.levels[3].lines = {buildLineWithNConductors(
        2, "line4", MTL_TYPE_SHIELDED, std::nullopt, "line3_2", 2)};
    return bundle;
}

inline mtl_m::transmission_line_bundle_t buildSevenCableBundle() {
    mtl_m::transmission_line_bundle_t bundle;
    bundle.levels.resize(4);
    bundle.levels[0].lines = {buildLineWithNConductors(1, "line1", MTL_TYPE_UNSHIELDED)};
    bundle.levels[1].lines = {buildLineWithNConductors(
        3, "line2", MTL_TYPE_SHIELDED, std::nullopt, "line1", 1)};
    bundle.levels[2].lines = {
        buildLineWithNConductors(1, "line3_1", MTL_TYPE_SHIELDED, std::nullopt, "line2", 1),
        buildLineWithNConductors(2, "line3_2", MTL_TYPE_SHIELDED, std::nullopt, "line2", 2),
        buildLineWithNConductors(2, "line3_3", MTL_TYPE_SHIELDED, std::nullopt, "line2", 3),
    };
    bundle.levels[3].lines = {
        buildLineWithNConductors(2, "line4_1", MTL_TYPE_SHIELDED, std::nullopt, "line3_2", 2),
        buildLineWithNConductors(1, "line4_2", MTL_TYPE_SHIELDED, std::nullopt, "line3_3", 1),
        buildLineWithNConductors(1, "line4_3", MTL_TYPE_SHIELDED, std::nullopt, "line3_3", 2),
    };
    return bundle;
}

struct test_mtl_bundle_t {
    int number_of_conductors = 0;
    int number_of_divisions = 0;
    std::vector<double> lpul;

    static int flat_size(int d1, int d2, int d3) { return d1 * d2 * d3; }

    void initialize_arrays() {
        lpul.assign(flat_size(number_of_divisions, number_of_conductors, number_of_conductors), 0.0);
    }

    static double& at3(std::vector<double>& arr, int i, int j, int k, int d2, int d3) {
        return arr[i * d2 * d3 + j * d3 + k];
    }
};

inline test_mtl_bundle_t buildTestBundle(const std::vector<mtl_m::transmission_line_level_t>& levels,
                                         const std::string& /*name*/) {
    test_mtl_bundle_t bundle;
    int n_cond = 0;
    for (const auto& level : levels) {
        for (const auto& line : level.lines) {
            n_cond += line.number_of_conductors;
        }
    }
    bundle.number_of_conductors = n_cond;
    bundle.number_of_divisions = levels.empty() || levels[0].lines.empty()
                                     ? 0
                                     : static_cast<int>(levels[0].lines[0].step_size.size());
    bundle.initialize_arrays();

    int n_sum = 0;
    const int n_div = bundle.number_of_divisions;
    for (const auto& level : levels) {
        for (const auto& line : level.lines) {
            const int n = line.number_of_conductors;
            for (int seg = 0; seg < n_div; ++seg) {
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        test_mtl_bundle_t::at3(bundle.lpul, seg, n_sum + i, n_sum + j, n_cond, n_cond) =
                            line.lpul[seg][i][j];
                    }
                }
            }
            n_sum += n;
        }
    }
    return bundle;
}

} // namespace mtln_test

namespace mtln_preprocess_test {

inline std::vector<int> conductorsInLevel(const mtl_m::transmission_line_bundle_t& line) {
    std::vector<int> res(line.levels.size(), 0);
    for (size_t i = 0; i < line.levels.size(); ++i) {
        for (const auto& l : line.levels[i].lines) {
            res[i] += l.number_of_conductors;
        }
    }
    return res;
}

inline int findConductorsBeforeCable(const std::string& name,
                                     const mtl_m::transmission_line_level_t& level) {
    int res = 0;
    for (const auto& l : level.lines) {
        if (l.name != name) {
            res += l.number_of_conductors;
        } else {
            return res;
        }
    }
    return res;
}

inline int findOuterConductorNumber(const mtl_m::mtl_t& line,
                                    const mtl_m::transmission_line_level_t& level,
                                    int conductors_in_level) {
    return findConductorsBeforeCable(line.parent_name, level) + conductors_in_level +
           line.conductor_in_parent;
}

inline std::vector<int> findInnerConductorRange(const mtl_m::mtl_t& line,
                                                const mtl_m::transmission_line_level_t& level,
                                                int conductors_in_level) {
    const int start = findConductorsBeforeCable(line.name, level) + conductors_in_level;
    std::vector<int> res(line.number_of_conductors);
    for (int k = 0; k < line.number_of_conductors; ++k) {
        res[k] = start + k + 1;
    }
    return res;
}

inline int sumFirstN(const std::vector<int>& v, int n) {
    int s = 0;
    for (int i = 0; i < n && i < static_cast<int>(v.size()); ++i) {
        s += v[i];
    }
    return s;
}

} // namespace mtln_preprocess_test

#endif
