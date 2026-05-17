```cpp
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <memory>
#include <stdexcept>

// Forward declarations for external types used in the module
// These would typically come from other headers like FDETYPES_m, Report_m, etc.
struct SGGFDTDINFO_t;
struct bounds_t;

// Assuming these constants/types are defined in FDETYPES_m
#ifndef RKIND
#define RKIND double
#endif

#ifndef INTEGERSIZEOFMEDIAMATRICES
#define INTEGERSIZEOFMEDIAMATRICES int
#endif

// Assuming these enums/constants are defined in FDETYPES_m or similar
#ifndef iHx
#define iHx 1
#endif
#ifndef iHy
#define iHy 2
#endif
#ifndef iHz
#define iHz 3
#endif

#ifndef iEx
#define iEx 1
#endif
#ifndef iEy
#define iEy 2
#endif
#ifndef iEz
#define iEx 3
#endif

#ifndef left
#define left 1
#endif
#ifndef right
#define right 2
#endif
#ifndef down
#define down 3
#endif
#ifndef up
#define up 4
#endif
#ifndef back
#define back 5
#endif
#ifndef front
#define front 6
#endif

// Placeholder for print11 and SEPARADOR from Report_m
void print11(int unit, const std::string& msg) {
    std::cout << msg << std::endl;
}
const std::string SEPARADOR = "========================================";

namespace BORDERS_MUR_m {

    struct xyzlimit_var_t {
        int XI[6];
        int XE[6];
        int YI[6];
        int YE[6];
        int ZI[6];
        int ZE[6];
    };

    struct LR_t {
        std::vector<std::vector<std::vector<RKIND>>> Past_Hx;
        std::vector<std::vector<std::vector<RKIND>>> Past_Hz;
        std::vector<std::vector<std::vector<RKIND>>> PastPast_Hx;
        std::vector<std::vector<std::vector<RKIND>>> PastPast_Hz;
    };

    struct DU_t {
        std::vector<std::vector<std::vector<RKIND>>> Past_Hy;
        std::vector<std::vector<std::vector<RKIND>>> Past_Hx;
        std::vector<std::vector<std::vector<RKIND>>> PastPast_Hy;
        std::vector<std::vector<std::vector<RKIND>>> PastPast_Hx;
    };

    struct BF_t {
        std::vector<std::vector<std::vector<RKIND>>> Past_Hz;
        std::vector<std::vector<std::vector<RKIND>>> Past_Hy;
        std::vector<std::vector<std::vector<RKIND>>> PastPast_Hz;
        std::vector<std::vector<std::vector<RKIND>>> PastPast_Hy;
    };

    // Global variables
    std::vector<std::vector<std::vector<xyzlimit_var_t>>> MURc; // Dimension 4:6, but indexed 1:6 in logic. Using size 7 for 1-based indexing convenience.
    std::vector<LR_t> regLR; // left:right (1:2)
    std::vector<DU_t> regDU; // down:up (3:4)
    std::vector<BF_t> regBF; // back:front (5:6)

    std::vector<RKIND> back_CAB1, back_CAB3, back_cab4;
    std::vector<RKIND> front_CAB1, front_CAB3, front_cab4;
    std::vector<RKIND> left_CAB1, left_CAB3, left_cab4;
    std::vector<RKIND> right_CAB1, right_CAB3, right_cab4;
    std::vector<RKIND> down_CAB1, down_CAB3, down_cab4;
    std::vector<RKIND> up_CAB1, up_CAB3, up_cab4;

    RKIND cluz = 0.0;
    RKIND eps0 = 0.0;
    RKIND mu0 = 0.0;

    // Helper to resize vectors to 1-based indexing if needed, or just use 0-based and adjust access
    // Here we will use 0-based vectors but access with adjusted indices.
    // However, to preserve logic exactly, we might allocate size N+1 and ignore index 0.

    void InitMURBorders(const SGGFDTDINFO_t& sgg, bool& ThereAreMURBorders, bool resume, 
                        const std::vector<RKIND>& Idxh, const std::vector<RKIND>& Idyh, const std::vector<RKIND>& Idzh, 
                        RKIND eps00, RKIND mu00) {
        
        eps0 = eps00; 
        mu0 = mu00;
        cluz = 1.0 / std::sqrt(eps0 * mu0);

        ThereAreMURBorders = false;
        if (sgg.Border.IsBackMUR || sgg.Border.IsFrontMUR || sgg.Border.IsLeftMUR || 
            sgg.Border.IsRightMUR || sgg.Border.IsUpMUR || sgg.Border.IsDownMUR) {
            ThereAreMURBorders = true;
        }
        if (!ThereAreMURBorders) return;

        int num_media = sgg.NumMedia;
        // Allocate CAB arrays. Fortran: 0 : sgg%NumMedia
        back_CAB1.assign(num_media + 1, 0.0);
        back_CAB3.assign(num_media + 1, 0.0);
        back_cab4.assign(num_media + 1, 0.0);
        front_CAB1.assign(num_media + 1, 0.0);
        front_CAB3.assign(num_media + 1, 0.0);
        front_cab4.assign(num_media + 1, 0.0);
        left_CAB1.assign(num_media + 1, 0.0);
        left_CAB3.assign(num_media + 1, 0.0);
        left_cab4.assign(num_media + 1, 0.0);
        right_CAB1.assign(num_media + 1, 0.0);
        right_CAB3.assign(num_media + 1, 0.0);
        right_cab4.assign(num_media + 1, 0.0);
        down_CAB1.assign(num_media + 1, 0.0);
        down_CAB3.assign(num_media + 1, 0.0);
        down_cab4.assign(num_media + 1, 0.0);
        up_CAB1.assign(num_media + 1, 0.0);
        up_CAB3.assign(num_media + 1, 0.0);
        up_cab4.assign(num_media + 1, 0.0);

        // MURc is dimensioned 4:6 in Fortran, but accessed with 1:6 (left,right,down,up,back,front).
        // Let's make it size 7 to use 1-based indexing directly.
        MURc.resize(7);

        for (int field = iHx; field <= iHz; ++field) {
            // Down
            MURc[field].XI[down] = sgg.Sweep[field].XI;
            MURc[field].XE[down] = sgg.Sweep[field].XE;
            MURc[field].YI[down] = sgg.Sweep[field].YI;
            MURc[field].YE[down] = sgg.Sweep[field].YE;
            MURc[field].ZI[down] = sgg.Sweep[field].ZI - 1;
            MURc[field].ZE[down] = MURc[field].ZI[down] + 1;

            // Up
            MURc[field].XI[up] = sgg.Sweep[field].XI;
            MURc[field].XE[up] = sgg.Sweep[field].XE;
            MURc[field].YI[up] = sgg.Sweep[field].YI;
            MURc[field].YE[up] = sgg.Sweep[field].YE;
            MURc[field].ZI[up] = sgg.Sweep[field].ZE;
            MURc[field].ZE[up] = MURc[field].ZI[up] + 1;

            // Left
            MURc[field].XI[left] = sgg.Sweep[field].XI;
            MURc[field].XE[left] = sgg.Sweep[field].XE;
            MURc[field].YI[left] = sgg.Sweep[field].YI - 1;
            MURc[field].YE[left] = MURc[field].YI[left] + 1;
            MURc[field].ZI[left] = sgg.Sweep[field].ZI;
            MURc[field].ZE[left] = sgg.Sweep[field].ZE;

            // Right
            MURc[field].XI[right] = sgg.Sweep[field].XI;
            MURc[field].XE[right] = sgg.Sweep[field].XE;
            MURc[field].YI[right] = sgg.Sweep[field].YE;
            MURc[field].YE[right] = MURc[field].YI[right] + 1;
            MURc[field].ZI[right] = sgg.Sweep[field].ZI;
            MURc[field].ZE[right] = sgg.Sweep[field].ZE;

            // Back
            MURc[field].XI[back] = sgg.Sweep[field].XI - 1;
            MURc[field].XE[back] = MURc[field].XI[back] + 1;
            MURc[field].YI[back] = sgg.Sweep[field].YI;
            MURc[field].YE[back] = sgg.Sweep[field].YE;
            MURc[field].ZI[back] = sgg.Sweep[field].ZI;
            MURc[field].ZE[back] = sgg.Sweep[field].ZE;

            // Front
            MURc[field].XI[front] = sgg.Sweep[field].XE;
            MURc[field].XE[front] = MURc[field].XI[front] + 1;
            MURc[field].YI[front] = sgg.Sweep[field].YI;
            MURc[field].YE[front] = sgg.Sweep[field].YE;
            MURc[field].ZI[front] = sgg.Sweep[field].ZI;
            MURc[field].ZE[front] = sgg.Sweep[field].ZE;
        }

        // Fake coms and ends
        if (!sgg.Border.IsDownMUR) {
            for (int f = 4; f <= 6; ++f) {
                MURc[f].ZI[down] = MURc[f].ZE[down] + 100;
            }
        }
        if (!sgg.Border.IsUpMUR) {
            for (int f = 4; f <= 6; ++f) {
                MURc[f].ZI[up] = MURc[f].ZE[up] + 100;
            }
        }
        if (!sgg.Border.IsLeftMUR) {
            for (int f = 4; f <= 6; ++f) {
                MURc[f].ZI[left] = MURc[f].ZE[left] + 100;
            }
        }
        if (!sgg.Border.IsRightMUR) {
            for (int f = 4; f <= 6; ++f) {
                MURc[f].ZI[right] = MURc[f].ZE[right] + 100;
            }
        }
        if (!sgg.Border.IsFrontMUR) {
            for (int f = 4; f <= 6; ++f) {
                MURc[f].ZI[front] = MURc[f].ZE[front] + 100;
            }
        }
        if (!sgg.Border.IsBackMUR) {
            for (int f = 4; f <= 6; ++f) {
                MURc[f].ZI[back] = MURc[f].ZE[back] + 100;
            }
        }

        // MUR Field component matrix allocation
        // regLR is dimensioned left:right (1:2)
        regLR.resize(3); // Index 1 and 2 used
        for (int region = left; region <= right; ++region) {
            int xi = MURc[iHx].XI[region];
            int xe = MURc[iHx].XE[region];
            int yi = MURc[iHx].YI[region];
            int ye = MURc[iHx].YE[region];
            int zi = MURc[iHx].ZI[region];
            int ze = MURc[iHx].ZE[region];

            // Allocate Past_Hx
            regLR[region].Past_Hx.resize(ze - zi + 1, std::vector<std::vector<RKIND>>(ye - yi + 1, std::vector<RKIND>(xe - xi + 1, 0.0)));
            // Allocate Past_Hz
            xi = MURc[iHz].XI[region];
            xe = MURc[iHz].XE[region];
            yi = MURc[iHz].YI[region];
            ye = MURc[iHz].YE[region];
            zi = MURc[iHz].ZI[region];
            ze = MURc[iHz].ZE[region];
            regLR[region].Past_Hz.resize(ze - zi + 1, std::vector<std::vector<RKIND>>(ye - yi + 1, std::vector<RKIND>(xe - xi + 1, 0.0)));

            if (!resume) {
                // Zeroed by constructor above
            } else {
                // Read from file 14
                // Note: In C++, we'd need a file stream. Assuming global file stream or passing it.
                // For translation, we'll assume a global ifstream or similar mechanism exists.
                // Here we just show the loop structure.
                std::ifstream file14("restart_data.bin", std::ios::binary); // Placeholder
                if (!file14) {
                     // Handle error
                }

                xi = MURc[iHx].XI[region];
                xe = MURc[iHx].XE[region];
                yi = MURc[iHx].YI[region];
                ye = MURc[iHx].YE[region];
                zi = MURc[iHx].ZI[region];
                ze = MURc[iHx].ZE[region];
                
                for (int k = zi; k <= ze; ++k) {
                    for (int j = yi; j <= ye; ++j) {
                        for (int i = xi; i <= xe; ++i) {
                            file14.read(reinterpret_cast<char*>(&regLR[region].Past_Hx[k - zi][j - yi][i - xi]), sizeof(RKIND));
                        }
                    }
                }

                xi = MURc[iHz].XI[region];
                xe = MURc[iHz].XE[region];
                yi = MURc[iHz].YI[region];
                ye = MURc[iHz].YE[region];
                zi = MURc[iHz].ZI[region];
                ze = MURc[iHz].ZE[region];

                for (int k = zi; k <= ze; ++k) {
                    for (int j = yi; j <= ye; ++j) {
                        for (int i = xi; i <= xe; ++i) {
                            file14.read(reinterpret_cast<char*>(&regLR[region].Past_Hz[k - zi][j - yi][i - xi]), sizeof(RKIND));
                        }
                    }
                }
            }
        }

        // regDU is dimensioned down:up (3:4)
        regDU.resize(5); // Index 3 and 4 used
        for (int region = down; region <= up; ++region) {
            int xi = MURc[iHy].XI[region];
            int xe = MURc[iHy].XE[region];
            int yi = MURc[iHy].YI[region];
            int ye = MURc[iHy].YE[region];
            int zi = MURc[iHy].ZI[region];
            int ze = MURc[iHy].ZE[region];

            regDU[region].Past_Hy.resize(ze - zi + 1, std::vector<std::vector<RKIND>>(ye - yi + 1, std::vector<RKIND>(xe - xi + 1, 0.0)));
            
            xi = MURc[iHx].XI[region];
            xe = MURc[iHx].XE[region];
            yi = MURc[iHx].YI[region];
            ye = MURc[iHx].YE[region];
            zi = MURc[iHx].ZI[region];
            ze = MURc[iHx].ZE[region];
            regDU[region].Past_Hx.resize(ze - zi + 1, std::vector<std::vector<RKIND>>(ye - yi + 1, std::vector<RKIND>(xe - xi + 1, 0.0)));

            if (!resume) {
                // Zeroed
            } else {
                std::ifstream file14("restart_data.bin", std::ios::binary);
                if (!file14) {}

                xi = MURc[iHy].XI[region];
                xe = MURc[iHy].XE[region];
                yi = MURc[iHy].YI[region];
                ye = MURc[iHy].YE[region];
                zi = MURc[iHy].ZI[region];
                ze = MURc[iHy].ZE[region];

                for (int k = zi; k <= ze; ++k) {
                    for (int j = yi; j <= ye; ++j) {
                        for (int i = xi; i <= xe; ++i) {
                            file14.read(reinterpret_cast<char*>(&regDU[region].Past_Hy[k - zi][j - yi][i - xi]), sizeof(RKIND));
                        }
                    }
                }

                xi = MURc[iHx].XI[region];
                xe = MURc[iHx].XE[region];
                yi = MURc[iHx].YI[region];
                ye = MURc[iHx].YE[region];
                zi = MURc[iHx].ZI[region];
                ze = MURc[iHx].ZE[region];

                for (int k = zi; k <= ze; ++k) {
                    for (int j = yi; j <= ye; ++j) {
                        for (int i = xi; i <= xe; ++i) {
                            file14.read(reinterpret_cast<char*>(&regDU[region].Past_Hx[k - zi][j - yi][i - xi]), sizeof(RKIND));
                        }
                    }
                }
            }
        }

        // regBF is dimensioned back:front (5:6)
        regBF.resize(7); // Index 5 and 6 used
        for (int region = back; region <= front; ++region) {
            int xi = MURc[iHz].XI[region];
            int xe = MURc[iHz].XE[region];
            int yi = MURc[iHz].YI[region];
            int ye = MURc[iHz].YE[region];
            int zi = MURc[iHz].ZI[region];
            int ze = MURc[iHz].ZE[region];

            regBF[region].Past_Hz.resize(ze - zi + 1, std::vector<std::vector<RKIND>>(ye - yi + 1, std::vector<RKIND>(xe - xi + 1, 0.0)));

            xi = MURc[iHy].XI[region];
            xe = MURc[iHy].XE[region];
            yi = MURc[iHy].YI[region];
            ye = MURc[iHy].YE[region];
            zi = MURc[iHy].ZI[region];
            ze = MURc[iHy].ZE[region];
            regBF[region].Past_Hy.resize(ze - zi + 1, std::vector<std::vector<RKIND>>(ye - yi + 1, std::vector<RKIND>(xe - xi + 1, 0.0)));

            if (!resume) {
                // Zeroed
            } else {
                std::ifstream file14("restart_data.bin", std::ios::binary);
                if (!file14) {}

                xi = MURc[iHz].XI[region];
                xe = MURc[iHz].XE[region];
                yi = MURc[iHz].YI[region];
                ye = MURc[iHz].YE[region];
                zi = MURc[iHz].ZI[region];
                ze = MURc[iHz].ZE[region];

                for (int k = zi; k <= ze; ++k) {
                    for (int j = yi; j <= ye; ++j) {
                        for (int i = xi; i <= xe; ++i) {
                            file14.read(reinterpret_cast<char*>(&regBF[region].Past_Hz[k - zi][j - yi][i - xi]), sizeof(RKIND));
                        }
                    }
                }

                xi = MURc[iHy].XI[region];
                xe = MURc[iHy].XE[region];
                yi = MURc[iHy].YI[region];
                ye = MURc[iHy].YE[region];
                zi = MURc[iHy].ZI[region];
                ze = MURc[iHy].ZE[region];

                for (int k = zi; k <= ze; ++k) {
                    for (int j = yi; j <= ye; ++j) {
                        for (int i = xi; i <= xe; ++i) {
                            file14.read(reinterpret_cast<char*>(&regBF[region].Past_Hy[k - zi][j - yi][i - xi]), sizeof(RKIND));
                        }
                    }
                }
            }
        }

        // Past Past allocation
        for (int region = left; region <= right; ++region) {
            int xi = MURc[iHx].XI[region];
            int xe = MURc[iHx].XE[region];
            int yi = MURc[iHx].YI[region];
            int ye = MURc[iHx].YE[region];
            int zi = MURc[iHx].ZI[region];
            int ze = MURc[iHx].ZE[region];
            regLR[region].PastPast_Hx.resize(ze - zi + 1, std::vector<std::vector<RKIND>>(ye - yi + 1, std::vector<RKIND>(xe - xi + 1, 0.0)));

            xi = MURc[iHz].XI[region];
            xe = MURc[iHz].XE[region];
            yi = MURc[iHz].YI[region];
            ye = MURc[iHz].YE[region];
            zi = MURc[iHz].ZI[region];
            ze = MURc[iHz].ZE[region];
            regLR[region].PastPast_Hz.resize(ze - zi + 1, std::vector<std::vector<RKIND>>(ye - yi + 1, std::vector<RKIND>(xe - xi + 1, 0.0)));

            if (!resume) {
                // Zeroed
            } else {
                std::ifstream file14("restart_data.bin", std::ios::binary);
                if (!file14) {}

                xi = MURc[iHx].XI[region];
                xe = MURc[iHx].XE[region];
                yi = MURc[iHx].YI[region];
                ye = MURc[iHx].YE[region];
                zi = MURc[iHx].ZI[region];
                ze = MURc[iHx].ZE[region];

                for (int k = zi; k <= ze; ++k) {
                    for (int j = yi; j <= ye; ++j) {
                        for (int i = xi; i <= xe; ++i) {
                            file14.read(reinterpret_cast<char*>(&regLR[region].PastPast_Hx[k - zi][j - yi][i - xi]), sizeof(RKIND));
                        }
                    }
                }

                xi = MURc[iHz].XI[region];
                xe = MURc[iHz].XE[region];
                yi = MURc[iHz].YI[region];
                ye = MURc[iHz].YE[region];
                zi = MURc[iHz].ZI[region];
                ze = MURc[iHz].ZE[region];

                for (int k = zi; k <= ze; ++k) {
                    for (int j = yi; j <= ye; ++j) {
                        for (int i = xi; i <= xe; ++i) {
                            file14.read(reinterpret_cast<char*>(&regLR[region].PastPast_Hz[k - zi][j - yi][i - xi]), sizeof(RKIND));
                        }
                    }
                }
            }
        }

        for (int region = down; region <= up; ++region) {
            int xi = MURc[iHy].XI[region];
            int xe = MURc[iHy].XE[region];
            int yi = MURc[iHy].YI[region];
            int ye = MURc[iHy].YE[region];
            int zi = MURc[iHy].ZI[region];
            int ze = MURc[iHy].ZE[region];
            regDU[region].PastPast_Hy.resize(ze - zi + 1, std::vector<std::vector<RKIND>>(ye - yi + 1, std::vector<RKIND>(xe - xi + 1, 0.0)));

            xi = MURc[iHx].XI[region];
            xe = MURc[iHx].XE[region];
            yi = MURc[iHx].YI[region];
            ye = MURc[iHx].YE[region];
            zi = MURc[iHx].ZI[region];
            ze = MURc[iHx].ZE[region];
            regDU[region].PastPast_Hx.resize(ze - zi + 1, std::vector<std::vector<RKIND>>(ye - yi + 1, std::vector<RKIND>(xe - xi + 1, 0.0)));

            if (!resume) {
                // Zeroed
            } else {
                std::ifstream file14("restart_data.bin", std::ios::binary);
                if (!file14) {}

                xi = MURc[iHy].XI[region];
                xe = MURc[iHy].XE[region];
                yi = MURc[iHy].YI[region];
                ye = MURc[iHy].YE[region];
                zi = MURc[iHy].ZI[region];
                ze = MURc[iHy].ZE[region];

                for (int k = zi; k <= ze; ++k) {
                    for (int j = yi; j <= ye; ++j) {
                        for (int i = xi; i <= xe; ++i) {
                            file14.read(reinterpret_cast<char*>(&regDU[region].PastPast_Hy[k - zi][j - yi][i - xi]), sizeof(RKIND));
                        }
                    }
                }

                xi = MURc[iHx].XI[region];
                xe = MURc[iHx].XE[region];
                yi = MURc[iHx].YI[region];
                ye = MURc[iHx].YE[region];
                zi = MURc[iHx].ZI[region];
                ze = MURc[iHx].ZE[region];

                for (int k = zi; k <= ze; ++k) {
                    for (int j = yi; j <= ye; ++j) {
                        for (int i = xi; i <= xe; ++i) {
                            file14.read(reinterpret_cast<char*>(&regDU[region].PastPast_Hx[k - zi][j - yi][i - xi]), sizeof(RKIND));
                        }
                    }
                }
            }
        }

        for (int region = back; region <= front; ++region) {
            int xi = MURc[iHz].XI[region];
            int xe = MURc[iHz].XE[region];
            int yi = MURc[iHz].YI[region];
            int ye = MURc[iHz].YE[region];
            int zi = MURc[iHz].ZI[region];
            int ze = MURc[iHz].ZE[region];
            regBF[region].PastPast_Hz.resize(ze - zi + 1, std::vector<std::vector<RKIND>>(ye - yi + 1, std::vector<RKIND>(xe - xi + 1, 0.0)));

            xi = MURc[iHy].XI[region];
            xe = MURc[iHy].XE[region];
            yi = MURc[iHy].YI[region];
            ye = MURc[iHy].YE[region];
            zi = MURc[iHy].ZI[region];
            ze = MURc[iHy].ZE[region];
            regBF[region].PastPast_Hy.resize(ze - zi + 1, std::vector<std::vector<RKIND>>(ye - yi + 1, std::vector<RKIND>(xe - xi + 1, 0.0)));

            if (!resume) {
                // Zeroed
            } else {
                std::ifstream file14("restart_data.bin", std::ios::binary);
                if (!file14) {}

                xi = MURc[iHz].XI[region];
                xe = MURc[iHz].XE[region];
                yi = MURc[iHz].YI[region];
                ye = MURc[iHz].YE[region];
                zi = MURc[iHz].ZI[region];
                ze = MURc[iHz].ZE[region];

                for (int k = zi; k <= ze; ++k) {
                    for (int j = yi; j <= ye; ++j) {
                        for (int i = xi; i <= xe; ++i) {
                            file14.read(reinterpret_cast<char*>(&regBF[region].PastPast_Hz[k - zi][j - yi][i - xi]), sizeof(RKIND));
                        }
                    }
                }

                xi = MURc[iHy].XI[region];
                xe = MURc[iHy].XE[region];
                yi = MURc[iHy].YI[region];
                ye = MURc[iHy].YE[region];
                zi = MURc[iHy].ZI[region];
                ze = MURc[iHy].ZE[region];

                for (int k = zi; k <= ze; ++k) {
                    for (int j = yi; j <= ye; ++j) {
                        for (int i = xi; i <= xe; ++i) {
                            file14.read(reinterpret_cast<char*>(&regBF[region].PastPast_Hy[k - zi][j - yi][i - xi]), sizeof(RKIND));
                        }
                    }
                }
            }
        }

        calc_murconstants(sgg, Idxh, Idyh, Idzh, eps0, mu0);
    }

    void calc_murconstants(const SGGFDTDINFO_t& sgg, const std::vector<RKIND>& Idxh, const std::vector<RKIND>& Idyh, const std::vector<RKIND>& Idzh, RKIND eps00, RKIND mu00) {
        eps0 = eps00; 
        mu0 = mu00;
        cluz = 1.0 / std::sqrt(eps0 * mu0);

        int num_media = sgg.NumMedia;
        for (int i1 = 0; i1 <= num_media; ++i1) {
            RKIND cnum;
            
            cnum = (1.0 / Idxh[sgg.ALLOC[iEx].XI]) / (sgg.dt * cluz / std::sqrt(sgg.Med[i1].Epr * sgg.Med[i1].Mur));
            back_CAB1[i1] = (1.0 - cnum) / (1.0 + cnum);
            back_CAB3[i1] = 1.0 / (2.0 * cnum * (1.0 + cnum));
            back_cab4[i1] = (2.0 * cnum / (1.0 + cnum) - 4.0 * (1.0 / (2.0 * cnum * (1.0 + cnum))));

            cnum = (1.0 / Idxh[sgg.ALLOC[iEx].XE]) / (sgg.dt * cluz / std::sqrt(sgg.Med[i1].Epr * sgg.Med[i1].Mur));
            front_CAB1[i1] = (1.0 - cnum) / (1.0 + cnum);
            front_CAB3[i1] = 1.0 / (2.0 * cnum * (1.0 + cnum));
            front_cab4[i1] = (2.0 * cnum / (1.0 + cnum) - 4.0 * (1.0 / (2.0 * cnum * (1.0 + cnum))));

            cnum = (1.0 / Idyh[sgg.ALLOC[iEy].YI]) / (sgg.dt * cluz / std::sqrt(sgg.Med[i1].Epr * sgg.Med[i1].Mur));
            left_CAB1[i1] = (1.0 - cnum) / (1.0 + cnum);
            left_CAB3[i1] = 1.0 / (2.0 * cnum * (1.0 + cnum));
            left_cab4[i1] = (2.0 * cnum / (1.0 + cnum) - 4.0 * (1.0 / (2.0 * cnum * (1.0 + cnum))));

            cnum = (1.0 / Idyh[sgg.ALLOC[iEy].YE]) / (sgg.dt * cluz / std::sqrt(sgg.Med[i1].Epr * sgg.Med[i1].Mur));
            right_CAB1[i1] = (1.0 - cnum) / (1.0 + cnum);
            right_CAB3[i1] = 1.0 / (2.0 * cnum * (1.0 + cnum));
            right_cab4[i1] = (2.0 * cnum / (1.0 + cnum)