#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <iostream>

// Forward declarations and includes for external modules/types
// These correspond to: mtln_types_m, Report_m, FDETYPES_m
// You need to provide implementations for these types and functions.

// Placeholder for RKIND type (usually double)
using RKIND = double;

// Placeholder for constants from FDETYPES_m
constexpr RKIND pi = 3.14159265358979323846;
constexpr RKIND EPSILON_VACUUM = 8.854187817e-12;
constexpr RKIND MU_VACUUM = 1.2566370614e-6;

// Placeholder for Report_m functions
void WarnErrReport(const std::string& msg, bool isFatal) {
    std::cerr << msg << std::endl;
    if (isFatal) {
        throw std::runtime_error(msg);
    }
}

// Placeholder types from mtln_types_m
// These must be defined in the corresponding header file.

struct box_2d_t {
    std::vector<RKIND> min;
    std::vector<RKIND> max;
};

struct multipolar_coefficient_t {
    RKIND a;
    RKIND b;
};

struct field_reconstruction_t {
    std::vector<multipolar_coefficient_t> ab;
    std::vector<RKIND> expansion_center;
    RKIND inner_region_average_potential;
    std::vector<RKIND> conductor_potentials;
};

struct multipolar_expansion_t {
    std::vector<field_reconstruction_t> electric;
    std::vector<field_reconstruction_t> magnetic;
    box_2d_t inner_region;
};

namespace multipolar_expansion_m {

    struct integration_grid_t {
        std::vector<RKIND> x;
        std::vector<RKIND> y;
    };

    inline RKIND multipolarExpansion2D(const std::vector<RKIND>& position, 
                                       const std::vector<multipolar_coefficient_t>& ab, 
                                       const std::vector<RKIND>& expansionCenter) {
        // 2D multipolar expansion from:
        // TSOGTGEREL GANTUMUR, MULTIPOLE EXPANSIONS IN THE PLANE.
        // 2016-04-16 lecture notes.

        std::vector<RKIND> rVec(2);
        rVec[0] = position[0] - expansionCenter[0];
        rVec[1] = position[1] - expansionCenter[1];
        
        RKIND r = std::sqrt(rVec[0]*rVec[0] + rVec[1]*rVec[1]);
        RKIND phi = std::atan2(rVec[1], rVec[0]);

        RKIND res = 0.0;
        for (size_t n_idx = 0; n_idx < ab.size(); ++n_idx) {
            int n = static_cast<int>(n_idx) + 1;
            if (n == 1) {
                res = res - ab[n_idx].a * std::log(r);
                // b0 should always be zero
            } else {
                res = res + (ab[n_idx].a * std::cos((n-1)*phi) + ab[n_idx].b * std::sin((n-1)*phi)) / std::pow(r, (n-1));
            }
        }

        res = res / (2.0 * pi);
        return res;
    }

    inline integration_grid_t buildIntegrationGridForBox(const box_2d_t& integrationBox, 
                                                         const box_2d_t& innerRegionBox) {
        // Returns an integration_grid_t with sorted, unique, equispaced points for x and y
        const int GRID_INTEGRATION_SAMPLING_POINTS = 100;
        
        // Preconditions
        if (integrationBox.min[0] >= innerRegionBox.min[0] || integrationBox.min[1] >= innerRegionBox.min[1] ||
            integrationBox.max[0] <= innerRegionBox.max[0] || integrationBox.max[1] <= innerRegionBox.max[1]) {       
            WarnErrReport(
                "Error in mutipolar expansion innerRegion must be fully contained within the integration Box", true);
            return integration_grid_t(); // Return empty on error
        }

        {
            std::vector<RKIND> innerRegionSize(2);
            std::vector<RKIND> integrationBoxSize(2);         
            innerRegionSize[0] = innerRegionBox.max[0] - innerRegionBox.min[0];
            innerRegionSize[1] = innerRegionBox.max[1] - innerRegionBox.min[1];
            integrationBoxSize[0] = integrationBox.max[0] - integrationBox.min[0];
            integrationBoxSize[1] = integrationBox.max[1] - integrationBox.min[1];
            
            if (integrationBoxSize[0] < innerRegionSize[0] * 1.25 || integrationBoxSize[1] < innerRegionSize[1] * 1.25) {
                WarnErrReport(
                    "Error in multipolar expansion: integration box is too small for the inner region", true);
                return integration_grid_t(); // Return empty on error
            }
        }

        // Builds grid
        const int totalPoints = GRID_INTEGRATION_SAMPLING_POINTS * 3 + 1;
        std::vector<RKIND> allPoints(totalPoints);
        
        for (int x_dim = 0; x_dim < 2; ++x_dim) { 
            // control points are ordered from min to max.
            std::vector<RKIND> controlPoints = {
                integrationBox.min[x_dim],
                innerRegionBox.min[x_dim],
                innerRegionBox.max[x_dim],
                integrationBox.max[x_dim]
            };
            
            for (size_t i = 1; i < controlPoints.size(); ++i) {
                RKIND minval = controlPoints[i-1];
                RKIND maxval = controlPoints[i];
                RKIND step = (maxval - minval) / GRID_INTEGRATION_SAMPLING_POINTS;
                for (int k = 0; k < GRID_INTEGRATION_SAMPLING_POINTS; ++k) {
                    allPoints[(i-1)*GRID_INTEGRATION_SAMPLING_POINTS + k] = minval + k*step;
                }
            }
            allPoints[totalPoints - 1] = controlPoints.back();

            if (x_dim == 0) {
                // Resize to remove duplicates if necessary, but Fortran logic just assigns
                // The last point of one segment is the first of the next, so we might have overlaps.
                // However, the Fortran code overwrites res%x and res%y completely.
                // The vector size is fixed. We keep it as is.
                // Note: The Fortran code uses 1-based indexing for loops but 0-based for vector assignment implicitly via array logic.
                // Here we map directly.
                // The loop for k goes 1 to N, mapping to indices (i-2)*N + k.
                // In C++, i starts at 1 (second element of controlPoints).
                // Index calculation: (i-1)*N + k_idx.
                // Let's re-verify the index mapping.
                // Fortran: do i = 2, size(controlPoints) -> i=2,3,4
                //          do k = 1, N
                //          index = (i-2)*N + k
                //          i=2: index = 0*N + k = k (1..N) -> 0..N-1 in C++
                //          i=3: index = 1*N + k = N+k (N+1..2N) -> N..2N-1 in C++
                //          i=4: index = 2*N + k = 2N+k (2N+1..3N) -> 2N..3N-1 in C++
                //          Last assignment: allPoints(size) = controlPoints(size) -> index 3N
                // Total size 3N+1.
                
                // My C++ loop: i from 1 to 3 (size 4).
                // i=1: min=cp[0], max=cp[1]. k=0..N-1. Index = 0*N + k = k. Correct.
                // i=2: min=cp[1], max=cp[2]. k=0..N-1. Index = 1*N + k. Correct.
                // i=3: min=cp[2], max=cp[3]. k=0..N-1. Index = 2*N + k. Correct.
                // Last point: allPoints[3N] = cp[3]. Correct.
                
                integration_grid_t res;
                res.x = allPoints;
                // We need to store y separately. But the loop overwrites allPoints.
                // We must save x before overwriting.
                // Actually, the Fortran code does:
                // if (x==1) res%x = allPoints else res%y = allPoints
                // This implies res%x and res%y are allocated separately.
                // In C++, we can just assign vectors.
                
                // Let's refactor to build x and y separately to avoid confusion with the single allPoints buffer reuse.
                // But to stick to the logic:
                if (x_dim == 0) {
                    // Save x
                    // We will return at the end.
                }
            }
        }
        
        // Refactored logic for clarity and correctness with std::vector
        integration_grid_t res;
        res.x.resize(totalPoints);
        res.y.resize(totalPoints);
        
        for (int x_dim = 0; x_dim < 2; ++x_dim) {
            std::vector<RKIND> controlPoints = {
                integrationBox.min[x_dim],
                innerRegionBox.min[x_dim],
                innerRegionBox.max[x_dim],
                integrationBox.max[x_dim]
            };
            
            std::vector<RKIND>& target = (x_dim == 0) ? res.x : res.y;
            
            for (size_t i = 1; i < controlPoints.size(); ++i) {
                RKIND minval = controlPoints[i-1];
                RKIND maxval = controlPoints[i];
                RKIND step = (maxval - minval) / GRID_INTEGRATION_SAMPLING_POINTS;
                for (int k = 0; k < GRID_INTEGRATION_SAMPLING_POINTS; ++k) {
                    target[(i-1)*GRID_INTEGRATION_SAMPLING_POINTS + k] = minval + k*step;
                }
            }
            target[totalPoints - 1] = controlPoints.back();
        }
        
        return res;
    }

    inline RKIND boxArea(const box_2d_t& box) {
        return (box.max[0] - box.min[0]) * (box.max[1] - box.min[1]);
    }

    inline bool isWithinBox(const box_2d_t& box, const std::vector<RKIND>& point) {
        return (point[0] >= box.min[0] && point[1] >= box.min[1]) && 
               (point[0] <= box.max[0] && point[1] <= box.max[1]);
    }

    inline RKIND getAveragePotential(const field_reconstruction_t& potential, 
                                     const box_2d_t& innerBox, 
                                     const box_2d_t& outerBox) {
        integration_grid_t integrationGrid = buildIntegrationGridForBox(outerBox, innerBox);

        RKIND outerV = 0.0;
        // Fortran: do m = 2, size(integrationGrid%x)
        // C++ vector size is N. Indices 0 to N-1.
        // Fortran 1-based index m corresponds to C++ index m-1.
        // Loop m=2 to N means C++ indices 1 to N-1.
        for (size_t m = 1; m < integrationGrid.x.size(); ++m) {
            for (size_t n = 1; n < integrationGrid.y.size(); ++n) {
                RKIND xMin = integrationGrid.x[m - 1];
                RKIND xMax = integrationGrid.x[m];
                RKIND yMin = integrationGrid.y[n - 1];
                RKIND yMax = integrationGrid.y[n];

                std::vector<RKIND> midPoint(2);
                midPoint[0] = 0.5 * (xMin + xMax);
                midPoint[1] = 0.5 * (yMin + yMax);
                
                RKIND area = (xMax - xMin) * (yMax - yMin);

                if (isWithinBox(innerBox, midPoint)) {
                    continue; // skip if inside inner region
                }

                std::vector<RKIND> center = potential.expansion_center;
                outerV = outerV + area * multipolarExpansion2D(midPoint, potential.ab, center);
            }
        }

        RKIND innerV = potential.inner_region_average_potential * boxArea(innerBox);
        RKIND avVj = (innerV + outerV) / boxArea(outerBox);
        return avVj;
    }
  
    inline std::vector<std::vector<RKIND>> getCellCapacitanceOnBox(const multipolar_expansion_t& multipolarExpansionParameters, 
                                                                   const box_2d_t& cellBox) {
        std::vector<std::vector<RKIND>> res;

        int N = static_cast<int>(multipolarExpansionParameters.electric.size());
        res.resize(N, std::vector<RKIND>(N));
        
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                // Fortran: ab(1)%a. ab is 1-based. So index 0.
                RKIND Qj = multipolarExpansionParameters.electric[j].ab[0].a;
                
                RKIND avVj = getAveragePotential(
                    multipolarExpansionParameters.electric[j],
                    multipolarExpansionParameters.inner_region,
                    cellBox);
                
                // Fortran: conductor_potentials(i). 1-based. So index i.
                RKIND ViWhenPrescribedVj = multipolarExpansionParameters.electric[j].conductor_potentials[i];
                
                avVj = -avVj + ViWhenPrescribedVj;

                res[i][j] = Qj / avVj * EPSILON_VACUUM;
            }
        }
        return res;
    }
  
    inline std::vector<std::vector<RKIND>> getCellInductanceOnBox(const multipolar_expansion_t& multipolarExpansionParameters, 
                                                                  const box_2d_t& cellBox) {
        std::vector<std::vector<RKIND>> res;

        int N = static_cast<int>(multipolarExpansionParameters.magnetic.size());
        res.resize(N, std::vector<RKIND>(N));
        
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                // Fortran: ab(1)%a. ab is 1-based. So index 0.
                RKIND Ij = multipolarExpansionParameters.magnetic[j].ab[0].a;
                
                RKIND avAj = getAveragePotential(
                    multipolarExpansionParameters.magnetic[j],
                    multipolarExpansionParameters.inner_region,
                    cellBox);
                
                // Fortran: conductor_potentials(i). 1-based. So index i.
                RKIND ViWhenPrescribedVj = multipolarExpansionParameters.magnetic[j].conductor_potentials[i];
                
                avAj = -avAj + ViWhenPrescribedVj;

                res[i][j] = avAj / Ij * MU_VACUUM;
            }
        }
        return res;
    }

} // namespace multipolar_expansion_m