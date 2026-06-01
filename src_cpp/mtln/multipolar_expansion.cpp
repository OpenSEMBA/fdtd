#include "multipolar_expansion_m.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace multipolar_expansion_m {

namespace {

constexpr RKIND pi = 3.14159265358979323846;
constexpr RKIND EPSILON_VACUUM = 8.854187817e-12;
constexpr RKIND MU_VACUUM = 1.2566370614e-6;

void WarnErrReport(const std::string& msg, bool isFatal) {
    std::cerr << msg << std::endl;
    if (isFatal) {
        throw std::runtime_error(msg);
    }
}

struct integration_grid_t {
    std::vector<RKIND> x;
    std::vector<RKIND> y;
};

integration_grid_t buildIntegrationGridForBox(const box_2d_t& integrationBox,
                                              const box_2d_t& innerRegionBox) {
    const int GRID_INTEGRATION_SAMPLING_POINTS = 100;

    if (integrationBox.min[0] >= innerRegionBox.min[0] ||
        integrationBox.min[1] >= innerRegionBox.min[1] ||
        integrationBox.max[0] <= innerRegionBox.max[0] ||
        integrationBox.max[1] <= innerRegionBox.max[1]) {
        WarnErrReport(
            "Error in mutipolar expansion innerRegion must be fully contained within the integration Box",
            true);
        return {};
    }

    const RKIND innerRegionSize0 = innerRegionBox.max[0] - innerRegionBox.min[0];
    const RKIND innerRegionSize1 = innerRegionBox.max[1] - innerRegionBox.min[1];
    const RKIND integrationBoxSize0 = integrationBox.max[0] - integrationBox.min[0];
    const RKIND integrationBoxSize1 = integrationBox.max[1] - integrationBox.min[1];

    if (integrationBoxSize0 < innerRegionSize0 * 1.25 ||
        integrationBoxSize1 < innerRegionSize1 * 1.25) {
        WarnErrReport(
            "Error in multipolar expansion: integration box is too small for the inner region",
            true);
        return {};
    }

    const int totalPoints = GRID_INTEGRATION_SAMPLING_POINTS * 3 + 1;
    integration_grid_t res;
    res.x.resize(totalPoints);
    res.y.resize(totalPoints);

    for (int x_dim = 0; x_dim < 2; ++x_dim) {
        const std::vector<RKIND> controlPoints = {
            integrationBox.min[x_dim],
            innerRegionBox.min[x_dim],
            innerRegionBox.max[x_dim],
            integrationBox.max[x_dim],
        };

        std::vector<RKIND>& target = (x_dim == 0) ? res.x : res.y;

        for (size_t i = 1; i < controlPoints.size(); ++i) {
            const RKIND minval = controlPoints[i - 1];
            const RKIND maxval = controlPoints[i];
            const RKIND step = (maxval - minval) / GRID_INTEGRATION_SAMPLING_POINTS;
            for (int k = 0; k < GRID_INTEGRATION_SAMPLING_POINTS; ++k) {
                target[(i - 1) * GRID_INTEGRATION_SAMPLING_POINTS + k] = minval + k * step;
            }
        }
        target[totalPoints - 1] = controlPoints.back();
    }

    return res;
}

RKIND boxArea(const box_2d_t& box) {
    return (box.max[0] - box.min[0]) * (box.max[1] - box.min[1]);
}

bool isWithinBox(const box_2d_t& box, const std::vector<RKIND>& point) {
    return point[0] >= box.min[0] && point[1] >= box.min[1] && point[0] <= box.max[0] &&
           point[1] <= box.max[1];
}

} // namespace

RKIND multipolarExpansion2D(const std::vector<RKIND>& position,
                            const std::vector<multipolar_coefficient_t>& ab,
                            const std::vector<RKIND>& expansionCenter) {
    const RKIND rVec0 = position[0] - expansionCenter[0];
    const RKIND rVec1 = position[1] - expansionCenter[1];
    const RKIND r = std::sqrt(rVec0 * rVec0 + rVec1 * rVec1);
    const RKIND phi = std::atan2(rVec1, rVec0);

    RKIND res = 0.0;
    for (size_t n_idx = 0; n_idx < ab.size(); ++n_idx) {
        const int n = static_cast<int>(n_idx) + 1;
        if (n == 1) {
            res -= ab[n_idx].a * std::log(r);
        } else {
            res += (ab[n_idx].a * std::cos((n - 1) * phi) + ab[n_idx].b * std::sin((n - 1) * phi)) /
                   std::pow(r, static_cast<RKIND>(n - 1));
        }
    }

    return res / (2.0 * pi);
}

RKIND getAveragePotential(const field_reconstruction_t& potential,
                          const box_2d_t& innerBox,
                          const box_2d_t& outerBox) {
    const integration_grid_t integrationGrid = buildIntegrationGridForBox(outerBox, innerBox);

    RKIND outerV = 0.0;
    for (size_t m = 1; m < integrationGrid.x.size(); ++m) {
        for (size_t n = 1; n < integrationGrid.y.size(); ++n) {
            const RKIND xMin = integrationGrid.x[m - 1];
            const RKIND xMax = integrationGrid.x[m];
            const RKIND yMin = integrationGrid.y[n - 1];
            const RKIND yMax = integrationGrid.y[n];

            const std::vector<RKIND> midPoint = {0.5 * (xMin + xMax), 0.5 * (yMin + yMax)};
            const RKIND area = (xMax - xMin) * (yMax - yMin);

            if (isWithinBox(innerBox, midPoint)) {
                continue;
            }

            outerV += area * multipolarExpansion2D(midPoint, potential.ab, potential.expansion_center);
        }
    }

    const RKIND innerV = potential.inner_region_average_potential * boxArea(innerBox);
    return (innerV + outerV) / boxArea(outerBox);
}

std::vector<std::vector<RKIND>> getCellCapacitanceOnBox(
    const multipolar_expansion_t& multipolarExpansionParameters,
    const box_2d_t& cellBox) {
    const int N = static_cast<int>(multipolarExpansionParameters.electric.size());
    std::vector<std::vector<RKIND>> res(N, std::vector<RKIND>(N));

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            const RKIND Qj = multipolarExpansionParameters.electric[j].ab[0].a;
            RKIND avVj = getAveragePotential(multipolarExpansionParameters.electric[j],
                                             multipolarExpansionParameters.inner_region,
                                             cellBox);
            const RKIND ViWhenPrescribedVj =
                multipolarExpansionParameters.electric[j].conductor_potentials[i];
            avVj = -avVj + ViWhenPrescribedVj;
            res[i][j] = Qj / avVj * EPSILON_VACUUM;
        }
    }
    return res;
}

std::vector<std::vector<RKIND>> getCellInductanceOnBox(
    const multipolar_expansion_t& multipolarExpansionParameters,
    const box_2d_t& cellBox) {
    const int N = static_cast<int>(multipolarExpansionParameters.magnetic.size());
    std::vector<std::vector<RKIND>> res(N, std::vector<RKIND>(N));

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            const RKIND Ij = multipolarExpansionParameters.magnetic[j].ab[0].a;
            RKIND avAj = getAveragePotential(multipolarExpansionParameters.magnetic[j],
                                             multipolarExpansionParameters.inner_region,
                                             cellBox);
            const RKIND ViWhenPrescribedVj =
                multipolarExpansionParameters.magnetic[j].conductor_potentials[i];
            avAj = -avAj + ViWhenPrescribedVj;
            res[i][j] = avAj / Ij * MU_VACUUM;
        }
    }
    return res;
}

} // namespace multipolar_expansion_m
