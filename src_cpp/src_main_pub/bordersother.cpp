#include <vector>
#include <string>
#include <iostream>

// Assuming these types are defined in FDETYPES_m
// We need to forward declare or include them. 
// Since I don't have the content of FDETYPES_m, I will define minimal stubs 
// to make the code compile structurally, preserving names.

namespace FDETYPES_m {
    using RKIND = double;
    
    struct logic_control_t {
        bool PeriodicBorders;
        bool PMCBorders;
        bool PECBorders;
    };

    struct XYZlimit_t {
        int XI, XE;
        int YI, YE;
        int ZI, ZE;
    };

    struct Border_t {
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
    };

    struct SGGFDTDINFO_t {
        Border_t Border;
    };
}

// Constants for indexing
// In Fortran, iHx, iHy, iHz are likely integers representing indices into the sggAlloc array.
// Based on the usage sggalloc(iHx), sggalloc(iHy), sggalloc(iHz), they are indices 0, 1, 2 or similar.
// Usually in FDTD codes:
// iHx = 0, iHy = 1, iHz = 2 is a common convention for field components.
// However, looking at the code:
// sggalloc(iHx)%XI ...
// It implies sggAlloc is an array of XYZlimit_t.
// Let's assume standard indices: iHx=0, iHy=1, iHz=2.
// If these are not defined in FDETYPES_m, we define them here as constants.
const int iHx = 0;
const int iHy = 1;
const int iHz = 2;

namespace BORDERS_other_m {

    void InitOtherBorders(const FDETYPES_m::SGGFDTDINFO_t& sgg, FDETYPES_m::logic_control_t& thereAre) {
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
        const std::vector<FDETYPES_m::XYZlimit_t>& sggAlloc,
        const FDETYPES_m::Border_t& sggBorder,
        std::vector<std::vector<std::vector<double>>>& Hx,
        std::vector<std::vector<std::vector<double>>>& Hy,
        std::vector<std::vector<std::vector<double>>>& Hz,
        const std::vector<FDETYPES_m::XYZlimit_t>& c,
        int layoutnumber,
        int num_procs
    ) {
        // Hx Down
        if (sggBorder.IsDownPMC) {
            if (layoutnumber == 0) {
                int z_idx = c[iHx].ZI - 1;
                int z_ref = c[iHx].ZI;
                // Assuming c is passed correctly and has the same dimensions as sggAlloc for limits
                // The loop bounds are implicit in the vector size or passed limits.
                // In Fortran: Hx( : , : ,C(iHx)%ZI-1)=-Hx( : , : ,C(iHx)%ZI)
                // We iterate over X and Y.
                for (int y = 0; y < Hx[0].size(); ++y) {
                    for (int x = 0; x < Hx[0][0].size(); ++x) {
                        Hx[x][y][z_idx] = -Hx[x][y][z_ref];
                    }
                }
            }
        }
        // Hx Up
        if (sggBorder.IsUpPMC) {
            if (layoutnumber == num_procs - 1) {
                int z_idx = c[iHx].ZE + 1;
                int z_ref = c[iHx].ZE;
                for (int y = 0; y < Hx[0].size(); ++y) {
                    for (int x = 0; x < Hx[0][0].size(); ++x) {
                        Hx[x][y][z_idx] = -Hx[x][y][z_ref];
                    }
                }
            }
        }
        // Hx Left
        if (sggBorder.IsLeftPMC) {
            int y_idx = c[iHx].YI - 1;
            int y_ref = c[iHx].YI;
            for (int z = 0; z < Hx[0][0].size(); ++z) {
                for (int x = 0; x < Hx[0][0].size(); ++x) {
                    Hx[x][y_idx][z] = -Hx[x][y_ref][z];
                }
            }
        }
        // Hx Right
        if (sggBorder.IsRightPMC) {
            int y_idx = c[iHx].YE + 1;
            int y_ref = c[iHx].YE;
            for (int z = 0; z < Hx[0][0].size(); ++z) {
                for (int x = 0; x < Hx[0][0].size(); ++x) {
                    Hx[x][y_idx][z] = -Hx[x][y_ref][z];
                }
            }
        }
        
        // Hy Back
        if (sggBorder.IsBackPMC) {
            int x_idx = c[iHy].XI - 1;
            int x_ref = c[iHy].XI;
            for (int z = 0; z < Hy[0][0].size(); ++z) {
                for (int y = 0; y < Hy[0].size(); ++y) {
                    Hy[x_idx][y][z] = -Hy[x_ref][y][z];
                }
            }
        }
        // Hy Front
        if (sggBorder.IsFrontPMC) {
            int x_idx = c[iHy].XE + 1;
            int x_ref = c[iHy].XE;
            for (int z = 0; z < Hy[0][0].size(); ++z) {
                for (int y = 0; y < Hy[0].size(); ++y) {
                    Hy[x_idx][y][z] = -Hy[x_ref][y][z];
                }
            }
        }
        // Hy Down
        if (sggBorder.IsDownPMC) {
            if (layoutnumber == 0) {
                int z_idx = c[iHy].ZI - 1;
                int z_ref = c[iHy].ZI;
                for (int y = 0; y < Hy[0].size(); ++y) {
                    for (int x = 0; x < Hy[0][0].size(); ++x) {
                        Hy[x][y][z_idx] = -Hy[x][y][z_ref];
                    }
                }
            }
        }
        // Hy Up
        if (sggBorder.IsUpPMC) {
            if (layoutnumber == num_procs - 1) {
                int z_idx = c[iHy].ZE + 1;
                int z_ref = c[iHy].ZE;
                for (int y = 0; y < Hy[0].size(); ++y) {
                    for (int x = 0; x < Hy[0][0].size(); ++x) {
                        Hy[x][y][z_idx] = -Hy[x][y][z_ref];
                    }
                }
            }
        }
        
        // Hz Back
        if (sggBorder.IsBackPMC) {
            int x_idx = c[iHz].XI - 1;
            int x_ref = c[iHz].XI;
            for (int z = 0; z < Hz[0][0].size(); ++z) {
                for (int y = 0; y < Hz[0].size(); ++y) {
                    Hz[x_idx][y][z] = -Hz[x_ref][y][z];
                }
            }
        }
        // Hz Front
        if (sggBorder.IsFrontPMC) {
            int x_idx = c[iHz].XE + 1;
            int x_ref = c[iHz].XE;
            for (int z = 0; z < Hz[0][0].size(); ++z) {
                for (int y = 0; y < Hz[0].size(); ++y) {
                    Hz[x_idx][y][z] = -Hz[x_ref][y][z];
                }
            }
        }
        // Hz Left
        if (sggBorder.IsLeftPMC) {
            int y_idx = c[iHz].YI - 1;
            int y_ref = c[iHz].YI;
            for (int z = 0; z < Hz[0][0].size(); ++z) {
                for (int x = 0; x < Hz[0][0].size(); ++x) {
                    Hz[x][y_idx][z] = -Hz[x][y_ref][z];
                }
            }
        }
        // Hz Right
        if (sggBorder.IsRightPMC) {
            int y_idx = c[iHz].YE + 1;
            int y_ref = c[iHz].YE;
            for (int z = 0; z < Hz[0][0].size(); ++z) {
                for (int x = 0; x < Hz[0][0].size(); ++x) {
                    Hz[x][y_idx][z] = -Hz[x][y_ref][z];
                }
            }
        }
    }

    void CloneMagneticPeriodic(
        const std::vector<FDETYPES_m::XYZlimit_t>& sggAlloc,
        const FDETYPES_m::Border_t& sggBorder,
        std::vector<std::vector<std::vector<double>>>& Hx,
        std::vector<std::vector<std::vector<double>>>& Hy,
        std::vector<std::vector<std::vector<double>>>& Hz,
        const std::vector<FDETYPES_m::XYZlimit_t>& c,
        int layoutnumber,
        int num_procs
    ) {
        // Hx Down
        if (sggBorder.IsDownPeriodic) {
            if (layoutnumber == 0) {
                int z_idx = c[iHx].ZI - 1;
                int z_ref = c[iHx].ZE;
                for (int y = 0; y < Hx[0].size(); ++y) {
                    for (int x = 0; x < Hx[0][0].size(); ++x) {
                        Hx[x][y][z_idx] = Hx[x][y][z_ref];
                    }
                }
            }
        }
        // Hx Up
        if (sggBorder.IsUpPeriodic) {
            if (layoutnumber == num_procs - 1) {
                int z_idx = c[iHx].ZE + 1;
                int z_ref = c[iHx].ZI;
                for (int y = 0; y < Hx[0].size(); ++y) {
                    for (int x = 0; x < Hx[0][0].size(); ++x) {
                        Hx[x][y][z_idx] = Hx[x][y][z_ref];
                    }
                }
            }
        }
        // Hx Left
        if (sggBorder.IsLeftPeriodic) {
            int y_idx = c[iHx].YI - 1;
            int y_ref = c[iHx].YE;
            for (int z = 0; z < Hx[0][0].size(); ++z) {
                for (int x = 0; x < Hx[0][0].size(); ++x) {
                    Hx[x][y_idx][z] = Hx[x][y_ref][z];
                }
            }
        }
        // Hx Right
        if (sggBorder.IsRightPeriodic) {
            int y_idx = c[iHx].YE + 1;
            int y_ref = c[iHx].YI;
            for (int z = 0; z < Hx[0][0].size(); ++z) {
                for (int x = 0; x < Hx[0][0].size(); ++x) {
                    Hx[x][y_idx][z] = Hx[x][y_ref][z];
                }
            }
        }
        
        // Hy Back
        if (sggBorder.IsBackPeriodic) {
            int x_idx = c[iHy].XI - 1;
            int x_ref = c[iHy].XE;
            for (int z = 0; z < Hy[0][0].size(); ++z) {
                for (int y = 0; y < Hy[0].size(); ++y) {
                    Hy[x_idx][y][z] = Hy[x_ref][y][z];
                }
            }
        }
        // Hy Front
        if (sggBorder.IsFrontPeriodic) {
            int x_idx = c[iHy].XE + 1;
            int x_ref = c[iHy].XI;
            for (int z = 0; z < Hy[0][0].size(); ++z) {
                for (int y = 0; y < Hy[0].size(); ++y) {
                    Hy[x_idx][y][z] = Hy[x_ref][y][z];
                }
            }
        }
        // Hy Down
        if (sggBorder.IsDownPeriodic) {
            if (layoutnumber == 0) {
                int z_idx = c[iHy].ZI - 1;
                int z_ref = c[iHy].ZE;
                for (int y = 0; y < Hy[0].size(); ++y) {
                    for (int x = 0; x < Hy[0][0].size(); ++x) {
                        Hy[x][y][z_idx] = Hy[x][y][z_ref];
                    }
                }
            }
        }
        // Hy Up
        if (sggBorder.IsUpPeriodic) {
            if (layoutnumber == num_procs - 1) {
                int z_idx = c[iHy].ZE + 1;
                int z_ref = c[iHy].ZI;
                for (int y = 0; y < Hy[0].size(); ++y) {
                    for (int x = 0; x < Hy[0][0].size(); ++x) {
                        Hy[x][y][z_idx] = Hy[x][y][z_ref];
                    }
                }
            }
        }
        
        // Hz Back
        if (sggBorder.IsBackPeriodic) {
            int x_idx = c[iHz].XI - 1;
            int x_ref = c[iHz].XE;
            for (int z = 0; z < Hz[0][0].size(); ++z) {
                for (int y = 0; y < Hz[0].size(); ++y) {
                    Hz[x_idx][y][z] = Hz[x_ref][y][z];
                }
            }
        }
        // Hz Front
        if (sggBorder.IsFrontPeriodic) {
            int x_idx = c[iHz].XE + 1;
            int x_ref = c[iHz].XI;
            for (int z = 0; z < Hz[0][0].size(); ++z) {
                for (int y = 0; y < Hz[0].size(); ++y) {
                    Hz[x_idx][y][z] = Hz[x_ref][y][z];
                }
            }
        }
        // Hz Left
        if (sggBorder.IsLeftPeriodic) {
            int y_idx = c[iHz].YI - 1;
            int y_ref = c[iHz].YE;
            for (int z = 0; z < Hz[0][0].size(); ++z) {
                for (int x = 0; x < Hz[0][0].size(); ++x) {
                    Hz[x][y_idx][z] = Hz[x][y_ref][z];
                }
            }
        }
        // Hz Right
        if (sggBorder.IsRightPeriodic) {
            int y_idx = c[iHz].YE + 1;
            int y_ref = c[iHz].YI;
            for (int z = 0; z < Hz[0][0].size(); ++z) {
                for (int x = 0; x < Hz[0][0].size(); ++x) {
                    Hz[x][y_idx][z] = Hz[x][y_ref][z];
                }
            }
        }
    }

}