```cpp
#include <vector>
#include <cstdint>

// Forward declarations for types defined in other modules (e.g., FDETYPES_m)
// These structs should be defined in their respective headers.

struct SGGFDTDINFO_t {
    struct {
        bool IsBackPeriodic;
        bool IsFrontPeriodic;
        bool IsLeftPeriodic;
        bool IsRightPeriodic;
        bool IsUpPeriodic;
        bool IsDownPeriodic;
        bool IsBackPMC;
        bool IsFrontPMC;
        bool IsLeftPMC;
        bool IsRightPMC;
        bool IsUpPMC;
        bool IsDownPMC;
        bool IsBackPEC;
        bool IsFrontPEC;
        bool IsLeftPEC;
        bool IsRightPEC;
        bool IsUpPEC;
        bool IsDownPEC;
    } Border;
};

struct logic_control_t {
    bool PeriodicBorders;
    bool PMCBorders;
    bool PECBorders;
};

struct XYZlimit_t {
    int32_t XI;
    int32_t XE;
    int32_t YI;
    int32_t YE;
    int32_t ZI;
    int32_t ZE;
};

struct Border_t {
    bool IsDownPMC;
    bool IsUpPMC;
    bool IsLeftPMC;
    bool IsRightPMC;
    bool IsBackPMC;
    bool IsFrontPMC;
    bool IsDownPeriodic;
    bool IsUpPeriodic;
    bool IsLeftPeriodic;
    bool IsRightPeriodic;
    bool IsBackPeriodic;
    bool IsFrontPeriodic;
};

// Constants typically defined in FDETYPES_m
// Assuming RKIND is double (real(8))
using RKIND = double;

// Indices for cell types, typically defined in FDETYPES_m or similar
// These are usually integer constants.
extern const int32_t iHx;
extern const int32_t iHy;
extern const int32_t iHz;

namespace BORDERS_other_m {

    void InitOtherBorders(const SGGFDTDINFO_t& sgg, logic_control_t& thereAre) {
        thereAre.PeriodicBorders = false;
        if (sgg.Border.IsBackPeriodic || sgg.Border.IsFrontPeriodic || sgg.Border.IsLeftPeriodic || 
            sgg.Border.IsRightPeriodic || sgg.Border.IsUpPeriodic || sgg.Border.IsDownPeriodic) {
            thereAre.PeriodicBorders = true;
        }

        thereAre.PMCBorders = false;
        if (sgg.Border.IsBackPMC || sgg.Border.IsFrontPMC || sgg.Border.IsLeftPMC || 
            sgg.Border.IsRightPMC || sgg.Border.IsUpPMC || sgg.Border.IsDownPMC) {
            thereAre.PMCBorders = true;
        }

        thereAre.PECBorders = false;
        if (sgg.Border.IsBackPEC || sgg.Border.IsFrontPEC || sgg.Border.IsLeftPEC || 
            sgg.Border.IsRightPEC || sgg.Border.IsUpPEC || sgg.Border.IsDownPEC) {
            thereAre.PECBorders = true;
        }
    }

    void MinusCloneMagneticPMC(
        const std::vector<XYZlimit_t>& sggAlloc,
        std::vector<double>& Hx,
        std::vector<double>& Hy,
        std::vector<double>& Hz,
        std::vector<XYZlimit_t>& c,
        int32_t layoutnumber,
        int32_t num_procs,
        const Border_t& sggBorder
    ) {
        // Note: In the original Fortran, 'c' is passed but used as if it were 'sggAlloc' or a local copy.
        // The Fortran code uses 'c(iHx)%ZI' etc. inside the subroutine, but 'c' is declared as a local variable 
        // of type XYZlimit_t(1:6) without initialization in the snippet provided, OR it's a typo for sggAlloc.
        // However, looking closely at the arguments:
        // type(XYZlimit_t), dimension(1:6), intent(in) :: sggAlloc
        // type(XYZlimit_t), dimension(1:6) :: c
        // The code uses 'c(iHx)%ZI'. Since 'c' is not initialized in the snippet, this is likely a bug in the 
        // original Fortran or 'c' is meant to be 'sggAlloc'. 
        // Given the context of FDTD codes, 'c' usually refers to the current cell limits. 
        // Let's assume 'c' is intended to be 'sggAlloc' based on typical usage patterns where limits are passed in.
        // However, to strictly preserve names and logic, we must use 'c'. 
        // But 'c' is uninitialized. 
        // Let's look at the call site if available? No.
        // Let's assume the Fortran code meant to use 'sggAlloc' for the limits, or 'c' is a typo for 'sggAlloc'.
        // Actually, in many FDTD codes, 'c' might be a local variable holding the limits for the current process.
        // Since we cannot initialize 'c' from nothing, and the Fortran code doesn't show it being set, 
        // we will assume 'c' is a typo for 'sggAlloc' in the usage, OR that 'c' is passed by reference and initialized elsewhere?
        // No, 'c' is declared as local 'type(XYZlimit_t), dimension(1:6) :: c'.
        // This is a critical ambiguity. 
        // However, looking at the variable names: sggAlloc is passed in. c is local.
        // It is highly probable that 'c' should be 'sggAlloc'. 
        // Let's check the second subroutine CloneMagneticPeriodic. Same pattern.
        // I will assume 'c' is a typo for 'sggAlloc' for the translation to be functional, 
        // OR I will map 'c' to 'sggAlloc' in the logic.
        // Let's map 'c' to 'sggAlloc' for the purpose of accessing limits.
        
        const auto& limits = sggAlloc; // Using sggAlloc as the source of limits, assuming 'c' was a typo for 'sggAlloc' in the original code's logic or 'c' was meant to be initialized from 'sggAlloc'.

        // Hx Down
        if (sggBorder.IsDownPMC) {
            if (layoutnumber == 0) {
                // Hx(:, :, C(iHx)%ZI-1) = -Hx(:, :, C(iHx)%ZI)
                // Assuming 1D vector storage for 3D array. 
                // We need to know the stride or dimensions. 
                // Since dimensions are not passed, we assume standard row-major or column-major?
                // Fortran is column-major. C++ is row-major.
                // Without explicit dimensions (Nx, Ny, Nz), we cannot correctly index a 1D vector.
                // However, the problem statement says "Convert Fortran arrays: use std::vector".
                // And "Preserve 1-based indexing where Fortran uses it".
                // This implies we might need to keep the 3D structure or handle indexing carefully.
                // Let's assume the vectors Hx, Hy, Hz are flattened 3D arrays.
                // We need the dimensions to map (x,y,z) to index.
                // The limits sggAlloc(iHx) give XI, XE, YI, YE, ZI, ZE.
                // Let's assume the array is stored in Fortran order (Z varies fastest? No, X varies slowest in Fortran? 
                // Fortran: A(1:Nx, 1:Ny, 1:Nz) -> index = i + (j-1)*Nx + (k-1)*Nx*Ny.
                // Let's assume the vector is sized appropriately.
                
                // To make this compilable and correct, we need the dimensions.
                // Since they are not in the signature, we must assume they are implicit or passed elsewhere.
                // For the sake of translation, I will use a helper lambda or assume a standard mapping.
                // But without dimensions, I cannot write the index calculation.
                // I will assume the vectors are 3D arrays represented as std::vector<double> and provide a placeholder 
                // for the indexing logic, or assume the user provides a 3D array class.
                // Given the constraints, I will use a simple 1D index calculation assuming the limits define the size.
                
                int32_t xi = limits[iHx].XI;
                int32_t xe = limits[iHx].XE;
                int32_t yi = limits[iHx].YI;
                int32_t ye = limits[iHx].YE;
                int32_t zi = limits[iHx].ZI;
                int32_t ze = limits[iHx].ZE;
                
                int32_t nx = xe - xi + 1;
                int32_t ny = ye - yi + 1;
                int32_t nz = ze - zi + 1;
                
                // Fortran index: Hx(i, j, k) -> index = (i-1) + (j-1)*nx + (k-1)*nx*ny
                // But our limits are 1-based relative to the global grid? Or local?
                // Assuming the vector Hx covers the range [xi, xe] x [yi, ye] x [zi, ze].
                // And the vector is contiguous.
                
                // Target: Hx(:, :, zi-1) = -Hx(:, :, zi)
                // Source Z index: zi
                // Dest Z index: zi-1
                
                // We need to copy the slice.
                // Since we are using std::vector, we can iterate.
                for (int32_t i = xi; i <= xe; ++i) {
                    for (int32_t j = yi; j <= ye; ++j) {
                        // Calculate indices
                        // Assuming the vector is indexed such that Hx(i,j,k) is at a specific offset.
                        // This is tricky without knowing the exact memory layout.
                        // Let's assume the vector is large enough and we use a helper function.
                        // For this translation, I will assume a helper function `get_index` exists or use a simple formula.
                        // Let's assume the vector is stored in Fortran order (column-major).
                        
                        // Index for (i, j, k) in a 3D array of size Nx x Ny x Nz starting at (xi, yi, zi):
                        // idx = (i - xi) + (j - yi) * nx + (k - zi) * nx * ny
                        // But wait, the vector might be global.
                        // Let's assume the vector Hx is exactly the size of the domain defined by limits.
                        
                        int32_t idx_src = (i - xi) + (j - yi) * nx + (zi - zi) * nx * ny;
                        int32_t idx_dst = (i - xi) + (j - yi) * nx + ((zi - 1) - zi) * nx * ny;
                        
                        // Check bounds for destination
                        if (idx_dst >= 0 && idx_dst < (int32_t)Hx.size()) {
                            Hx[idx_dst] = -Hx[idx_src];
                        }
                    }
                }
            }
        }

        // Hx Up
        if (sggBorder.IsUpPMC) {
            if (layoutnumber == num_procs - 1) {
                int32_t xi = limits[iHx].XI;
                int32_t xe = limits[iHx].XE;
                int32_t yi = limits[iHx].YI;
                int32_t ye = limits[iHx].YE;
                int32_t zi = limits[iHx].ZI;
                int32_t ze = limits[iHx].ZE;
                
                int32_t nx = xe - xi + 1;
                int32_t ny = ye - yi + 1;
                
                // Hx(:, :, ze+1) = -Hx(:, :, ze)
                for (int32_t i = xi; i <= xe; ++i) {
                    for (int32_t j = yi; j <= ye; ++j) {
                        int32_t idx_src = (i - xi) + (j - yi) * nx + (ze - zi) * nx * ny;
                        int32_t idx_dst = (i - xi) + (j - yi) * nx + ((ze + 1) - zi) * nx * ny;
                        
                        if (idx_dst >= 0 && idx_dst < (int32_t)Hx.size()) {
                            Hx[idx_dst] = -Hx[idx_src];
                        }
                    }
                }
            }
        }

        // Hx Left
        if (sggBorder.IsLeftPMC) {
            int32_t xi = limits[iHx].XI;
            int32_t xe = limits[iHx].XE;
            int32_t yi = limits[iHx].YI;
            int32_t ye = limits[iHx].YE;
            int32_t zi = limits[iHx].ZI;
            int32_t ze = limits[iHx].ZE;
            
            int32_t nx = xe - xi + 1;
            int32_t ny = ye - yi + 1;
            int32_t nz = ze - zi + 1;
            
            // Hx(:, yi-1, :) = -Hx(:, yi, :)
            for (int32_t i = xi; i <= xe; ++i) {
                for (int32_t k = zi; k <= ze; ++k) {
                    int32_t idx_src = (i - xi) + (yi - yi) * nx + (k - zi) * nx * ny;
                    int32_t idx_dst = (i - xi) + ((yi - 1) - yi) * nx + (k - zi) * nx * ny;
                    
                    if (idx_dst >= 0 && idx_dst < (int32_t)Hx.size()) {
                        Hx[idx_dst] = -Hx[idx_src];
                    }
                }
            }
        }

        // Hx Right
        if (sggBorder.IsRightPMC) {
            int32_t xi = limits[iHx].XI;
            int32_t xe = limits[iHx].XE;
            int32_t yi = limits[iHx].YI;
            int32_t ye = limits[iHx].YE;
            int32_t zi = limits[iHx].ZI;
            int32_t ze = limits[iHx].ZE;
            
            int32_t nx = xe - xi + 1;
            int32_t ny = ye - yi + 1;
            int32_t nz = ze - zi + 1;
            
            // Hx(:, ye+1, :) = -Hx(:, ye, :)
            for (int32_t i = xi; i <= xe; ++i) {
                for (int32_t k = zi; k <= ze; ++k) {
                    int32_t idx_src = (i - xi) + (ye - yi) * nx + (k - zi) * nx * ny;
                    int32_t idx_dst = (i - xi) + ((ye + 1) - yi) * nx + (k - zi) * nx * ny;
                    
                    if (idx_dst >= 0 && idx_dst < (int32_t)Hx.size()) {
                        Hx[idx_dst] = -Hx[idx_src];
                    }
                }
            }
        }

        // Hy Back
        if (sggBorder.IsBackPMC) {
            int32_t xi = limits[iHy].XI;
            int32_t xe = limits[iHy].XE;
            int32_t yi = limits[iHy].YI;
            int32_t ye = limits[iHy].YE;
            int32_t zi = limits[iHy].ZI;
            int32_t ze = limits[iHy].ZE;
            
            int32_t nx = xe - xi + 1;
            int32_t ny = ye - yi + 1;
            int32_t nz = ze - zi + 1;
            
            // Hy(xi-1, :, :) = -Hy(xi, :, :)
            for (int32_t j = yi; j <= ye; ++j) {
                for (int32_t k = zi; k <= ze; ++k) {
                    int32_t idx_src = (xi - xi) + (j - yi) * nx + (k - zi) * nx * ny;
                    int32_t idx_dst = ((xi - 1) - xi) + (j - yi) * nx + (k - zi) * nx * ny;
                    
                    if (idx_dst >= 0 && idx_dst < (int32_t)Hy.size()) {
                        Hy[idx_dst] = -Hy[idx_src];
                    }
                }
            }
        }

        // Hy Front
        if (sggBorder.IsFrontPMC) {
            int32_t xi = limits[iHy].XI;
            int32_t xe = limits[iHy].XE;
            int32_t yi = limits[iHy].YI;
            int32_t ye = limits[iHy].YE;
            int32_t zi = limits[iHy].ZI;
            int32_t ze = limits[iHy].ZE;
            
            int32_t nx = xe - xi + 1;
            int32_t ny = ye - yi + 1;
            int32_t nz = ze - zi + 1;
            
            // Hy(xe+1, :, :) = -Hy(xe, :, :)
            for (int32_t j = yi; j <= ye; ++j) {
                for (int32_t k = zi; k <= ze; ++k) {
                    int32_t idx_src = (xe - xi) + (j - yi) * nx + (k - zi) * nx * ny;
                    int32_t idx_dst = ((xe + 1) - xi) + (j - yi) * nx + (k - zi) * nx * ny;
                    
                    if (idx_dst >= 0 && idx_dst < (int32_t)Hy.size()) {
                        Hy[idx_dst] = -Hy[idx_src];
                    }
                }
            }
        }

        // Hy Down
        if (sggBorder.IsDownPMC) {
            if (layoutnumber == 0) {
                int32_t xi = limits[iHy].XI;
                int32_t xe = limits[iHy].XE;
                int32_t yi = limits[iHy].YI;
                int32_t ye = limits[iHy].YE;
                int32_t zi = limits[iHy].ZI;
                int32_t ze = limits[iHy].ZE;
                
                int32_t nx = xe - xi + 1;
                int32_t ny = ye - yi + 1;
                
                // Hy(:, :, zi-1) = -Hy(:, :, zi)
                for (int32_t i = xi; i <= xe; ++i) {
                    for (int32_t j = yi; j <= ye; ++j) {
                        int32_t idx_src = (i - xi) + (j - yi) * nx + (zi - zi) * nx * ny;
                        int32_t idx_dst = (i - xi) + (j - yi) * nx + ((zi - 1) - zi) * nx * ny;
                        
                        if (idx_dst >= 0 && idx_dst < (int32_t)Hy.size()) {
                            Hy[idx_dst] = -Hy[idx_src];
                        }
                    }
                }
            }
        }

        // Hy Up
        if (sggBorder.IsUpPMC) {
            if (layoutnumber == num_procs - 1) {
                int32_t xi = limits[iHy].XI;
                int32_t xe = limits[iHy].XE;
                int32_t yi = limits[iHy].YI;
                int32_t ye = limits[iHy].YE;
                int32_t zi = limits[iHy].ZI;
                int32_t ze = limits[iHy].ZE;
                
                int32_t nx = xe - xi + 1;
                int32_t ny = ye - yi + 1;
                
                // Hy(:, :, ze+1) = -Hy(:, :, ze)
                for (int32_t i = xi; i <= xe; ++i) {
                    for (int32_t j = yi; j <= ye; ++j) {
                        int32_t idx_src = (i - xi) + (j - yi) * nx + (ze - zi) * nx * ny;
                        int32_t idx_dst = (i - xi) + (j - yi) * nx + ((ze + 1) - zi) * nx * ny;
                        
                        if (idx_dst >= 0 && idx_dst < (int32_t)Hy.size()) {
                            Hy[idx_dst] = -Hy[idx_src];
                        }
                    }
                }
            }
        }

        // Hz Back
        if (sggBorder.IsBackPMC) {
            int32_t xi = limits[iHz].XI;
            int32_t xe = limits[iHz].XE;
            int32_t yi = limits[iHz].YI;
            int32_t ye = limits[iHz].YE;
            int32_t zi = limits[iHz].ZI;
            int32_t ze = limits[iHz].ZE;
            
            int32_t nx = xe - xi + 1;
            int32_t ny = ye - yi + 1;
            int32_t nz = ze - zi + 1;
            
            // Hz(xi-1, :, :) = -Hz(xi, :, :)
            for (int32_t j = yi; j <= ye; ++j) {
                for (int32_t k = zi; k <= ze; ++k) {
                    int32_t idx_src = (xi - xi) + (j - yi) * nx + (k - zi) * nx * ny;
                    int32_t idx_dst = ((xi - 1) - xi) + (j - yi) * nx + (k - zi) * nx * ny;
                    
                    if (idx_dst >= 0 && idx_dst < (int32_t)Hz.size()) {
                        Hz[idx_dst] = -Hz[idx_src];
                    }
                }
            }
        }

        // Hz Front
        if (sggBorder.IsFrontPMC) {
            int32_t xi = limits[iHz].XI;
            int32_t xe = limits[iHz].XE;
            int32_t yi = limits[iHz].YI;
            int32_t ye = limits[iHz].YE;
            int32_t zi = limits[iHz].ZI;
            int32_t ze = limits[iHz].ZE;
            
            int32_t nx = xe - xi + 1;
            int32_t ny = ye - yi + 1;
            int32_t nz = ze - zi + 1;
            
            // Hz(xe+1, :, :) = -Hz(xe, :, :)
            for (int32_t j = yi; j <= ye; ++j) {
                for (int32_t k = zi; k <= ze; ++k) {
                    int32_t idx_src = (xe - xi) + (j - yi) * nx + (k - zi) * nx * ny;
                    int32_t idx_dst = ((xe + 1) - xi) + (j - yi) * nx + (k - zi) * nx * ny;
                    
                    if (idx_dst >= 0 && idx_dst < (int32_t)Hz.size()) {
                        Hz[idx_dst] = -Hz[idx_src];
                    }
                }
            }
        }

        // Hz Left
        if (sggBorder.IsLeftPMC) {
            int32_t xi = limits[iHz].XI;
            int32_t xe = limits[iHz].XE;
            int32_t yi = limits[iHz].YI;
            int32_t ye = limits[iHz].YE;
            int32_t zi = limits[iHz].ZI;
            int32_t ze = limits[iHz].ZE;
            
            int32_t nx = xe - xi + 1;
            int32_t ny = ye - yi + 1;
            int32_t nz = ze - zi + 1;
            
            // Hz(:, yi-1, :) = -Hz(:, yi, :)
            for (int32_t i = xi; i <= xe; ++i) {
                for (int32_t k = zi; k <= ze; ++k) {
                    int32_t idx_src = (i - xi) + (yi - yi) * nx + (k - zi) * nx * ny;
                    int32_t idx_dst = (i - xi) + ((yi - 1) - yi) * nx + (k - zi) * nx * ny;
                    
                    if (idx_dst >= 0 && idx_dst < (int32_t)Hz.size()) {
                        Hz[idx_dst] = -Hz[idx_src];
                    }
                }
            }
        }

        // Hz Right
        if (sggBorder.IsRightPMC) {
            int32_t xi = limits[iHz].XI;
            int32_t xe = limits[iHz].XE;
            int32_t yi = limits[iHz].YI;
            int32_t ye = limits[iHz].YE;
            int32_t zi = limits[iHz].ZI;
            int32_t ze = limits[iHz].ZE;
            
            int32_t nx = xe - xi + 1;
            int32_t ny = ye - yi + 1;
            int32_t nz = ze - zi + 1;
            
            // Hz(:, ye+1, :) = -Hz(:, ye, :)
            for (int32_t i = xi; i <= xe; ++i) {
                for (int32_t k = zi; k <= ze; ++k) {
                    int32_t idx_src = (i - xi) + (ye - yi) * nx + (k - zi) * nx * ny;
                    int32_t idx_dst = (i - xi) + ((ye + 1) - yi) * nx + (k - zi) * nx * ny;
                    
                    if (idx_dst >= 0 && idx_dst < (int32_t)Hz.size()) {
                        Hz[idx_dst] = -Hz[idx_src];
                    }
                }
            }
        }
    }

    void CloneMagneticPeriodic(
        const std::vector<XYZlimit_t>& sggAlloc,
        std::vector<double>& Hx,
        std::vector<double>& Hy,
        std::vector<double>& Hz,
        std::vector<XYZlimit_t>& c,
        int32_t layoutnumber,
        int32_t num_procs,
        const Border_t& sggBorder
    ) {
        const auto& limits = sggAlloc;

        // Hx Down
        if (sggBorder.IsDownPeriodic) {
            if (layoutnumber == 0) {
                int32_t xi = limits[iHx].XI;
                int32_t xe = limits[iHx].XE;
                int32_t yi = limits[iHx].YI;
                int32_t ye = limits[iHx].YE;
                int32_t zi = limits[iHx].ZI;
                int32_t ze = limits[iHx].ZE;
                
                int32_t nx = xe - xi + 1;
                int32_t ny = ye - yi + 1;
                
                // Hx(:, :, zi-1) = Hx(:, :, ze)
                for (int32_t i = xi; i <= xe; ++i) {
                    for (int32_t j = yi; j <= ye; ++j) {
                        int32_t idx_src = (i - xi) + (j - yi) * nx + (ze - zi) * nx * ny;
                        int32_t idx_dst = (i - xi) + (j - yi) * nx + ((zi - 1) - zi) * nx * ny;
                        
                        if (idx_dst >= 0 && idx_dst < (int32_t)Hx.size()) {
                            Hx[idx_dst] = Hx[idx_src];
                        }
                    }
                }
            }
        }

        // Hx Up
        if (sggBorder.IsUpPeriodic) {
            if (layoutnumber == num_procs - 1) {
                int32_t xi = limits[iHx].XI;
                int32_t xe = limits[iHx].XE;
                int32_t yi = limits[iHx].YI;
                int32_t ye = limits[iHx].YE;
                int32_t zi = limits[iHx].ZI;
                int32_t ze = limits[iHx].ZE;
                
                int32_t nx = xe - xi + 1;
                int32_t ny = ye - yi + 1;
                
                // Hx(:, :, ze+1) = Hx(:, :, zi)
                for (int32_t i = xi; i <= xe; ++i) {
                    for (int32_t j = yi; j <= ye; ++j) {
                        int32_t idx_src = (i - xi) + (j - yi) * nx + (zi - zi) * nx * ny;
                        int32_t idx_dst = (i - xi) + (j - yi) * nx + ((ze + 1) - zi) * nx * ny;
                        
                        if (idx_dst >= 0 && idx_dst < (int32_t)Hx.size()) {
                            Hx[idx_dst] = Hx[idx_src];
                        }
                    }
                }
            }
        }

        // Hx Left
        if (sggBorder.IsLeftPeriodic) {
            int32_t xi = limits[iHx].XI;
            int32_t xe = limits[iHx].XE;
            int32_t yi = limits[iHx].YI;
            int32_t ye = limits[iHx].YE;
            int32_t zi = limits[iHx].ZI;
            int32_t ze = limits[iHx].ZE;
            
            int32_t nx = xe - xi + 1;
            int32_t ny = ye - yi + 1;
            int32_t nz = ze - zi + 1;
            
            // Hx(:, yi-1, :) = Hx(:, ye, :)
            for (int32_t i = xi; i <= xe; ++i) {
                for (int32_t k = zi; k <= ze; ++k) {
                    int32_t idx_src = (i - xi) + (ye - yi) * nx + (k - zi) * nx * ny;
                    int32_t idx_dst = (i - xi) + ((yi - 1) - yi) * nx + (k - zi) * nx * ny;
                    
                    if (idx_dst >= 0 && idx_dst < (int32_t)Hx.size()) {
                        Hx[idx_dst] = Hx[idx_src];
                    }
                }
            }
        }

        // Hx Right
        if (sggBorder.IsRightPeriodic) {
            int32_t xi = limits[iHx].XI;
            int32_t xe = limits[iHx].XE;
            int32_t yi = limits[iHx].YI;
            int32_t ye = limits[iHx].YE;
            int32_t zi = limits[iHx].ZI;
            int32_t ze = limits[iHx].ZE;
            
            int32_t nx = xe - xi + 1;
            int32_t ny = ye - yi + 1;
            int32_t nz = ze - zi + 1;
            
            // Hx(:, ye+1, :) = Hx(:, yi, :)
            for (int32_t i = xi; i <= xe; ++i) {
                for (int32_t k = zi; k <= ze; ++k) {
                    int32_t idx_src = (i - xi) + (yi - yi) * nx + (k - zi) * nx * ny;
                    int32_t idx_dst = (i - xi) + ((ye + 1) - yi) * nx + (k - zi) * nx * ny;
                    
                    if (idx_dst >= 0 && idx_dst < (int32_t)Hx.size()) {
                        Hx[idx_dst] = Hx[idx_src];
                    }
                }
            }
        }

        // Hy Back
        if (sggBorder.IsBackPeriodic) {
            int32_t xi = limits[iHy].XI;
            int32_t xe = limits[iHy].XE;
            int32_t yi = limits[iHy].YI;
            int32_t ye = limits[iHy].YE;
            int32_t zi = limits[iHy].ZI;
            int32_t ze = limits[iHy].ZE;
            
            int32_t nx = xe - xi + 1;
            int32_t ny = ye - yi + 1;
            int32_t nz = ze - zi + 1;
            
            // Hy(xi-1, :, :) = Hy(xe, :, :)
            for (int32_t j = yi; j <= ye; ++j) {
                for (int32_t k = zi; k <= ze; ++k) {
                    int32_t idx_src = (xe - xi) + (j - yi) * nx + (k - zi) * nx * ny;
                    int32_t idx_dst = ((xi - 1) - xi) + (j - yi) * nx + (k - zi) * nx * ny;
                    
                    if (idx_dst >= 0 && idx_dst < (int32_t)Hy.size()) {
                        Hy[idx_dst] = Hy[idx_src];
                    }
                }
            }
        }

        // Hy Front
        if (sggBorder.IsFrontPeriodic) {
            int32_t xi = limits[iHy].XI;
            int32_t xe = limits[iHy].XE;
            int32_t yi = limits[iHy].YI;
            int32_t ye = limits[iHy].YE;
            int32_t zi = limits[iHy].ZI;
            int32_t ze = limits[iHy].ZE;
            
            int32_t nx = xe - xi + 1;
            int32_t ny = ye - yi + 1;
            int32_t nz = ze - zi + 1;
            
            // Hy(xe+1, :, :) = Hy(xi, :, :)
            for (int32_t j = yi; j <= ye; ++j) {
                for (int32_t k = zi; k <= ze; ++