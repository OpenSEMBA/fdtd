#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <cstring>

// Forward declarations and includes for external modules/types referenced in the Fortran code
// These would typically be in their own headers (Report_m.h, FDETYPES_m.h, etc.)
// Assuming these types exist based on usage:

// Placeholder for types from Report_m
struct sim_control_t {
    int layoutnumber;
    int num_procs;
    bool resume;
};

struct SGGFDTDINFO_t {
    struct {
        int XI, XE, YI, YE, ZI, ZE;
    } alloc[10]; // Assuming iEx, iEy, etc. map to indices
    struct {
        int XI, XE, YI, YE, ZI, ZE;
    } SINPMLSweep[10];
    int nummedia;
    // Placeholder for media structure
    struct {
        struct {
            bool PMLbody;
            struct {
                int orient;
            } PMLbody[1]; // Simplified
        } Is;
    } Med[100];
    double dt;
};

struct media_matrices_t {
    int sggMiEx[100][100][100];
    int sggMiEy[100][100][100];
    int sggMiEz[100][100][100];
    int sggMiHx[100][100][100];
    int sggMiHy[100][100][100];
    int sggMiHz[100][100][100];
};

// Constants from FDETYPES_m (assumed)
typedef double RKIND; // Assuming RKIND is double
const int iEx = 0, iEy = 1, iEz = 2;
const int iHx = 3, iHy = 4, iHz = 5;
const int BUFSIZE = 256;

// External functions from Report_m
void WarnErrReport(const std::string& msg, bool fatal);
void print11(int unit, const std::string& msg);
const std::string SEPARADOR = "========================================";

namespace PMLbodies_m {

    struct BerPML__t {
        // In Fortran, these are pointers to specific array elements.
        // In C++, we store pointers to the actual data arrays or use indices.
        // However, to preserve the structure and logic closely, we'll store pointers.
        // Note: The original code uses pointer association. In C++, we must ensure the target arrays outlive these pointers.
        RKIND* field;
        RKIND* Plus;
        RKIND* Minu;
        
        RKIND gx2;
        RKIND P_be;
        RKIND P_ce;
        RKIND Psi;
        RKIND transversalDelta;
        RKIND del;
        RKIND posi;
        int minTotal;
        int maxTotal;
    };

    struct berpml_t {
        int NumNodes;
        int orient;
        std::vector<BerPML__t> nodes;
    };

    // Local variables from the module
    berpml_t berpmlE;
    berpml_t berpmlH;

    const int PMLorden = 2;
    const RKIND CoeffReflPML = 1e-4;

    // Global variables
    RKIND eps0 = 0.0;
    RKIND mu0 = 0.0;

    void InitPMLbodies(
        const SGGFDTDINFO_t& sgg,
        const media_matrices_t& media,
        const RKIND* Ex,
        const RKIND* Ey,
        const RKIND* Ez,
        const RKIND* Hx,
        const RKIND* Hy,
        const RKIND* Hz,
        const int* Idxe,
        const int* Idye,
        const int* Idze,
        const int* Idxh,
        const int* Idyh,
        const int* Idzh,
        const RKIND* g2,
        const RKIND* gm2,
        bool& ThereArePMLbodies,
        const sim_control_t& control,
        RKIND eps00,
        RKIND mu00
    ) {
        eps0 = eps00;
        mu0 = mu00;

        std::string whoami;
        char buf[BUFSIZE];
        snprintf(buf, BUFSIZE, "(%d/%d) ", control.layoutnumber + 1, control.num_procs);
        whoami = buf;

        ThereArePMLbodies = false;

        // Initialize min/max arrays
        std::vector<int> minx(sgg.nummedia + 1, std::numeric_limits<int>::max());
        std::vector<int> miny(sgg.nummedia + 1, std::numeric_limits<int>::max());
        std::vector<int> minz(sgg.nummedia + 1, std::numeric_limits<int>::max());
        std::vector<int> maxx(sgg.nummedia + 1, std::numeric_limits<int>::min());
        std::vector<int> maxy(sgg.nummedia + 1, std::numeric_limits<int>::min());
        std::vector<int> maxz(sgg.nummedia + 1, std::numeric_limits<int>::min());

        int conta = 0;

        // Helper lambda to process E-field sweeps
        auto process_E_sweep = [&](int field_idx, const int* (*get_mi)(const media_matrices_t&, int, int, int), const RKIND* (*get_field)(const RKIND*, int, int, int)) {
            int xi = sgg.SINPMLSweep[field_idx].XI;
            int xe = sgg.SINPMLSweep[field_idx].XE;
            int yi = sgg.SINPMLSweep[field_idx].YI;
            int ye = sgg.SINPMLSweep[field_idx].YE;
            int zi = sgg.SINPMLSweep[field_idx].ZI;
            int ze = sgg.SINPMLSweep[field_idx].ZE;

            for (int k1 = zi; k1 <= ze; ++k1) {
                for (int j1 = yi; j1 <= ye; ++j1) {
                    for (int i1 = xi; i1 <= xe; ++i1) {
                        int jmed = get_mi(&media, i1, j1, k1);
                        if (sgg.Med[jmed].Is.PMLbody) {
                            int orient = std::abs(sgg.Med[jmed].Is.PMLbody[0].orient);
                            if (orient != field_idx) {
                                if (i1 < minx[jmed]) minx[jmed] = i1;
                                if (j1 < miny[jmed]) miny[jmed] = j1;
                                if (k1 < minz[jmed]) minz[jmed] = k1;
                                if (i1 > maxx[jmed]) maxx[jmed] = i1;
                                if (j1 > maxy[jmed]) maxy[jmed] = j1;
                                if (k1 > maxz[jmed]) maxz[jmed] = k1;
                                conta++;
                            }
                        }
                    }
                }
            }
        };

        // Count E-field nodes
        process_E_sweep(iEx, [](const media_matrices_t* m, int i, int j, int k) { return m->sggMiEx[i][j][k]; }, 
                        [](const RKIND* arr, int i, int j, int k) { return &arr[i * 100 * 100 + j * 100 + k]; }); // Simplified indexing assumption
        
        // Note: The original code uses 3D arrays passed as pointers. 
        // For translation, we assume contiguous memory or appropriate indexing.
        // The original code passes 3D arrays. In C++, we'll assume they are flattened or accessed via a helper.
        // To keep it simple and close to original, we'll assume the arrays are passed as pointers to the first element
        // and accessed with stride. However, the original Fortran code uses explicit bounds.
        // Let's refine the lambda to match the original logic more closely without assuming stride.
        
        // Reset conta for E-field counting
        conta = 0;
        
        // iEx
        for (int k1 = sgg.SINPMLSweep[iEx].ZI; k1 <= sgg.SINPMLSweep[iEx].ZE; ++k1) {
            for (int j1 = sgg.SINPMLSweep[iEx].YI; j1 <= sgg.SINPMLSweep[iEx].YE; ++j1) {
                for (int i1 = sgg.SINPMLSweep[iEx].XI; i1 <= sgg.SINPMLSweep[iEx].XE; ++i1) {
                    int jmed = media.sggMiEx[i1][j1][k1];
                    if (sgg.Med[jmed].Is.PMLbody) {
                        int orient = std::abs(sgg.Med[jmed].Is.PMLbody[0].orient);
                        if (orient != iEx) {
                            if (i1 < minx[jmed]) minx[jmed] = i1;
                            if (j1 < miny[jmed]) miny[jmed] = j1;
                            if (k1 < minz[jmed]) minz[jmed] = k1;
                            if (i1 > maxx[jmed]) maxx[jmed] = i1;
                            if (j1 > maxy[jmed]) maxy[jmed] = j1;
                            if (k1 > maxz[jmed]) maxz[jmed] = k1;
                            conta++;
                        }
                    }
                }
            }
        }
        // iEy
        for (int k1 = sgg.SINPMLSweep[iEy].ZI; k1 <= sgg.SINPMLSweep[iEy].ZE; ++k1) {
            for (int j1 = sgg.SINPMLSweep[iEy].YI; j1 <= sgg.SINPMLSweep[iEy].YE; ++j1) {
                for (int i1 = sgg.SINPMLSweep[iEy].XI; i1 <= sgg.SINPMLSweep[iEy].XE; ++i1) {
                    int jmed = media.sggMiEy[i1][j1][k1];
                    if (sgg.Med[jmed].Is.PMLbody) {
                        int orient = std::abs(sgg.Med[jmed].Is.PMLbody[0].orient);
                        if (orient != iEy) {
                            if (i1 < minx[jmed]) minx[jmed] = i1;
                            if (j1 < miny[jmed]) miny[jmed] = j1;
                            if (k1 < minz[jmed]) minz[jmed] = k1;
                            if (i1 > maxx[jmed]) maxx[jmed] = i1;
                            if (j1 > maxy[jmed]) maxy[jmed] = j1;
                            if (k1 > maxz[jmed]) maxz[jmed] = k1;
                            conta++;
                        }
                    }
                }
            }
        }
        // iEz
        for (int k1 = sgg.SINPMLSweep[iEz].ZI; k1 <= sgg.SINPMLSweep[iEz].ZE; ++k1) {
            for (int j1 = sgg.SINPMLSweep[iEz].YI; j1 <= sgg.SINPMLSweep[iEz].YE; ++j1) {
                for (int i1 = sgg.SINPMLSweep[iEz].XI; i1 <= sgg.SINPMLSweep[iEz].XE; ++i1) {
                    int jmed = media.sggMiEz[i1][j1][k1];
                    if (sgg.Med[jmed].Is.PMLbody) {
                        int orient = std::abs(sgg.Med[jmed].Is.PMLbody[0].orient);
                        if (orient != iEz) {
                            if (i1 < minx[jmed]) minx[jmed] = i1;
                            if (j1 < miny[jmed]) miny[jmed] = j1;
                            if (k1 < minz[jmed]) minz[jmed] = k1;
                            if (i1 > maxx[jmed]) maxx[jmed] = i1;
                            if (j1 > maxy[jmed]) maxy[jmed] = j1;
                            if (k1 > maxz[jmed]) maxz[jmed] = k1;
                            conta++;
                        }
                    }
                }
            }
        }

        ThereArePMLbodies = (conta != 0);
        if (!ThereArePMLbodies) {
            return;
        }

        berpmlE.NumNodes = conta;
        berpmlE.nodes.resize(conta);

        conta = 0;
        // Count H-field nodes
        // iHx
        for (int k1 = sgg.SINPMLSweep[iHx].ZI; k1 <= sgg.SINPMLSweep[iHx].ZE; ++k1) {
            for (int j1 = sgg.SINPMLSweep[iHx].YI; j1 <= sgg.SINPMLSweep[iHx].YE; ++j1) {
                for (int i1 = sgg.SINPMLSweep[iHx].XI; i1 <= sgg.SINPMLSweep[iHx].XE; ++i1) {
                    int jmed = media.sggMiHx[i1][j1][k1];
                    if (sgg.Med[jmed].Is.PMLbody) {
                        int orient = std::abs(sgg.Med[jmed].Is.PMLbody[0].orient);
                        if (orient != iEx) {
                            conta++;
                        }
                    }
                }
            }
        }
        // iHy
        for (int k1 = sgg.SINPMLSweep[iHy].ZI; k1 <= sgg.SINPMLSweep[iHy].ZE; ++k1) {
            for (int j1 = sgg.SINPMLSweep[iHy].YI; j1 <= sgg.SINPMLSweep[iHy].YE; ++j1) {
                for (int i1 = sgg.SINPMLSweep[iHy].XI; i1 <= sgg.SINPMLSweep[iHy].XE; ++i1) {
                    int jmed = media.sggMiHy[i1][j1][k1];
                    if (sgg.Med[jmed].Is.PMLbody) {
                        int orient = std::abs(sgg.Med[jmed].Is.PMLbody[0].orient);
                        if (orient != iEy) {
                            conta++;
                        }
                    }
                }
            }
        }
        // iHz
        for (int k1 = sgg.SINPMLSweep[iHz].ZI; k1 <= sgg.SINPMLSweep[iHz].ZE; ++k1) {
            for (int j1 = sgg.SINPMLSweep[iHz].YI; j1 <= sgg.SINPMLSweep[iHz].YE; ++j1) {
                for (int i1 = sgg.SINPMLSweep[iHz].XI; i1 <= sgg.SINPMLSweep[iHz].XE; ++i1) {
                    int jmed = media.sggMiHz[i1][j1][k1];
                    if (sgg.Med[jmed].Is.PMLbody) {
                        int orient = std::abs(sgg.Med[jmed].Is.PMLbody[0].orient);
                        if (orient != iEz) {
                            conta++;
                        }
                    }
                }
            }
        }

        ThereArePMLbodies = ThereArePMLbodies && (conta != 0);
        if (!ThereArePMLbodies) {
            snprintf(buf, BUFSIZE, "Buggy ERROR: In PMLbodies. fields exist withouth Hfields. ");
            WarnErrReport(std::string(buf), true);
        }
        berpmlH.NumNodes = conta;
        berpmlH.nodes.resize(conta);

        // Initialize E-field nodes
        conta = 0;
        // iEx
        for (int k1 = sgg.SINPMLSweep[iEx].ZI; k1 <= sgg.SINPMLSweep[iEx].ZE; ++k1) {
            for (int j1 = sgg.SINPMLSweep[iEx].YI; j1 <= sgg.SINPMLSweep[iEx].YE; ++j1) {
                for (int i1 = sgg.SINPMLSweep[iEx].XI; i1 <= sgg.SINPMLSweep[iEx].XE; ++i1) {
                    int jmed = media.sggMiEx[i1][j1][k1];
                    if (sgg.Med[jmed].Is.PMLbody) {
                        int orient = std::abs(sgg.Med[jmed].Is.PMLbody[0].orient);
                        if (orient != iEx) {
                            conta++;
                            BerPML__t& PML_ = berpmlE.nodes[conta - 1]; // 0-based index
                            // Pointer association: In C++, we take the address of the element in the array
                            // Assuming 3D arrays are flattened or accessed via a helper. 
                            // For simplicity, we assume the arrays are passed as pointers to 3D blocks.
                            // We need to calculate the offset.
                            // Let's assume the arrays are contiguous 3D arrays: Ex[XI:XE][YI:YE][ZI:ZE]
                            // The original code uses Ex(i1,j1,k1).
                            // We will store pointers to the actual memory locations.
                            PML_.field = const_cast<RKIND*>(&Ex[i1 * 100 * 100 + j1 * 100 + k1]); // Assuming stride 100x100
                            berpmlE.orient = orient;
                            switch (orient) {
                                case iEy:
                                    PML_.del = 1.0 / Idye[j1];
                                    PML_.transversalDelta = 1.0 / Idyh[j1];
                                    PML_.Plus = const_cast<RKIND*>(&Hz[i1 * 100 * 100 + j1 * 100 + k1]);
                                    PML_.Minu = const_cast<RKIND*>(&Hz[i1 * 100 * 100 + (j1 - 1) * 100 + k1]);
                                    PML_.gx2 = g2[jmed];
                                    PML_.minTotal = miny[jmed];
                                    PML_.maxTotal = maxy[jmed];
                                    PML_.posi = j1;
                                    break;
                                case iEz:
                                    PML_.del = 1.0 / Idze[k1];
                                    PML_.transversalDelta = 1.0 / Idzh[k1];
                                    PML_.Plus = const_cast<RKIND*>(&Hy[i1 * 100 * 100 + j1 * 100 + k1]);
                                    PML_.Minu = const_cast<RKIND*>(&Hy[i1 * 100 * 100 + j1 * 100 + (k1 - 1)]);
                                    PML_.gx2 = -g2[jmed];
                                    PML_.minTotal = minz[jmed];
                                    PML_.maxTotal = maxz[jmed];
                                    PML_.posi = k1;
                                    break;
                                default:
                                    snprintf(buf, BUFSIZE, "Buggy ERROR: In PMLbodies. ");
                                    WarnErrReport(std::string(buf), true);
                                    break;
                            }
                        }
                    }
                }
            }
        }
        // iEy
        for (int k1 = sgg.SINPMLSweep[iEy].ZI; k1 <= sgg.SINPMLSweep[iEy].ZE; ++k1) {
            for (int j1 = sgg.SINPMLSweep[iEy].YI; j1 <= sgg.SINPMLSweep[iEy].YE; ++j1) {
                for (int i1 = sgg.SINPMLSweep[iEy].XI; i1 <= sgg.SINPMLSweep[iEy].XE; ++i1) {
                    int jmed = media.sggMiEy[i1][j1][k1];
                    if (sgg.Med[jmed].Is.PMLbody) {
                        int orient = std::abs(sgg.Med[jmed].Is.PMLbody[0].orient);
                        if (orient != iEy) {
                            conta++;
                            BerPML__t& PML_ = berpmlE.nodes[conta - 1];
                            PML_.field = const_cast<RKIND*>(&Ey[i1 * 100 * 100 + j1 * 100 + k1]);
                            berpmlE.orient = orient;
                            switch (orient) {
                                case iEz:
                                    PML_.del = 1.0 / Idze[k1];
                                    PML_.transversalDelta = 1.0 / Idzh[k1];
                                    PML_.Plus = const_cast<RKIND*>(&Hx[i1 * 100 * 100 + j1 * 100 + k1]);
                                    PML_.Minu = const_cast<RKIND*>(&Hx[i1 * 100 * 100 + j1 * 100 + (k1 - 1)]);
                                    PML_.minTotal = minz[jmed];
                                    PML_.maxTotal = maxz[jmed];
                                    PML_.gx2 = g2[jmed];
                                    PML_.posi = k1;
                                    break;
                                case iEx:
                                    PML_.del = 1.0 / Idxe[i1];
                                    PML_.transversalDelta = 1.0 / Idxh[i1];
                                    PML_.Plus = const_cast<RKIND*>(&Hz[i1 * 100 * 100 + j1 * 100 + k1]);
                                    PML_.Minu = const_cast<RKIND*>(&Hz[(i1 - 1) * 100 * 100 + j1 * 100 + k1]);
                                    PML_.gx2 = -g2[jmed];
                                    PML_.minTotal = minx[jmed];
                                    PML_.maxTotal = maxx[jmed];
                                    PML_.posi = i1;
                                    break;
                                default:
                                    snprintf(buf, BUFSIZE, "Buggy ERROR: In PMLbodies. ");
                                    WarnErrReport(std::string(buf), true);
                                    break;
                            }
                        }
                    }
                }
            }
        }
        // iEz
        for (int k1 = sgg.SINPMLSweep[iEz].ZI; k1 <= sgg.SINPMLSweep[iEz].ZE; ++k1) {
            for (int j1 = sgg.SINPMLSweep[iEz].YI; j1 <= sgg.SINPMLSweep[iEz].YE; ++j1) {
                for (int i1 = sgg.SINPMLSweep[iEz].XI; i1 <= sgg.SINPMLSweep[iEz].XE; ++i1) {
                    int jmed = media.sggMiEz[i1][j1][k1];
                    if (sgg.Med[jmed].Is.PMLbody) {
                        int orient = std::abs(sgg.Med[jmed].Is.PMLbody[0].orient);
                        if (orient != iEz) {
                            conta++;
                            BerPML__t& PML_ = berpmlE.nodes[conta - 1];
                            PML_.field = const_cast<RKIND*>(&Ez[i1 * 100 * 100 + j1 * 100 + k1]);
                            berpmlE.orient = orient;
                            switch (orient) {
                                case iEx:
                                    PML_.del = 1.0 / Idxe[i1];
                                    PML_.transversalDelta = 1.0 / Idxh[i1];
                                    PML_.Plus = const_cast<RKIND*>(&Hy[i1 * 100 * 100 + j1 * 100 + k1]);
                                    PML_.Minu = const_cast<RKIND*>(&Hy[(i1 - 1) * 100 * 100 + j1 * 100 + k1]);
                                    PML_.gx2 = g2[jmed];
                                    PML_.minTotal = minx[jmed];
                                    PML_.maxTotal = maxx[jmed];
                                    PML_.posi = i1;
                                    break;
                                case iEy:
                                    PML_.del = 1.0 / Idye[j1];
                                    PML_.transversalDelta = 1.0 / Idyh[j1];
                                    PML_.Plus = const_cast<RKIND*>(&Hx[i1 * 100 * 100 + j1 * 100 + k1]);
                                    PML_.Minu = const_cast<RKIND*>(&Hx[i1 * 100 * 100 + (j1 - 1) * 100 + k1]);
                                    PML_.gx2 = -g2[jmed];
                                    PML_.minTotal = miny[jmed];
                                    PML_.maxTotal = maxy[jmed];
                                    PML_.posi = j1;
                                    break;
                                default:
                                    snprintf(buf, BUFSIZE, "Buggy ERROR: In PMLbodies. ");
                                    WarnErrReport(std::string(buf), true);
                                    break;
                            }
                        }
                    }
                }
            }
        }

        // Initialize H-field nodes
        conta = 0;
        // iHx
        for (int k1 = sgg.SINPMLSweep[iHx].ZI; k1 <= sgg.SINPMLSweep[iHx].ZE; ++k1) {
            for (int j1 = sgg.SINPMLSweep[iHx].YI; j1 <= sgg.SINPMLSweep[iHx].YE; ++j1) {
                for (int i1 = sgg.SINPMLSweep[iHx].XI; i1 <= sgg.SINPMLSweep[iHx].XE; ++i1) {
                    int jmed = media.sggMiHx[i1][j1][k1];
                    if (sgg.Med[jmed].Is.PMLbody) {
                        int orient = std::abs(sgg.Med[jmed].Is.PMLbody[0].orient);
                        if (orient != iEx) {
                            conta++;
                            BerPML__t& PML_ = berpmlH.nodes[conta - 1];
                            PML_.field = const_cast<RKIND*>(&Hx[i1 * 100 * 100 + j1 * 100 + k1]);
                            berpmlE.orient = orient; // Note: Original code sets berpmlE.orient here, likely a bug in Fortran, but we preserve it
                            switch (orient) {
                                case iEy:
                                    PML_.del = 1.0 / Idye[j1];
                                    PML_.transversalDelta = 1.0 / Idye[j1];
                                    PML_.Plus = const_cast<RKIND*>(&Ez[i1 * 100 * 100 + (j1 + 1) * 100 + k1]);
                                    PML_.Minu = const_cast<RKIND*>(&Ez[i1 * 100 * 100 + j1 * 100 + k1]);
                                    PML_.gx2 = gm2[jmed];
                                    PML_.minTotal = miny[jmed];
                                    PML_.maxTotal = maxy[jmed];
                                    PML_.posi = j1 + 0.5;
                                    break;
                                case iEz:
                                    PML_.del = 1.0 / Idze[k1];
                                    PML_.transversalDelta = 1.0 / Idze[k1];
                                    PML_.Plus = const_cast<RKIND*>(&Ey[i1 * 100 * 100 + j1 * 100 + (k1 + 1)]);
                                    PML_.Minu = const_cast<RKIND*>(&Ey[i1 * 100 * 100 + j1 * 100 + k1]);
                                    PML_.gx2 = -gm2[jmed];
                                    PML_.minTotal = minz[jmed];
                                    PML_.maxTotal = maxz[jmed];
                                    PML_.posi = k1 + 0.5;
                                    break;
                                default:
                                    snprintf(buf, BUFSIZE, "Buggy ERROR: In PMLbodies. ");
                                    WarnErrReport(std::string(buf), true);
                                    break;
                            }
                        }
                    }
                }
            }
        }
        // iHy
        for (int k1 = sgg.SINPMLSweep[iHy].ZI; k1 <= sgg.SINPMLSweep[iHy].ZE; ++k1) {
            for (int j1 = sgg.SINPMLSweep[iHy].YI; j1 <= sgg.SINPMLSweep[iHy].YE; ++j1) {
                for (int i1 = sgg.SINPMLSweep[iHy].XI; i1 <= sgg.SINPMLSweep[iHy].XE; ++i1) {
                    int jmed = media.sggMiHy[i1][j1][k1];
                    if (sgg.Med[jmed].Is.PMLbody) {
                        int orient = std::abs(sgg.Med[jmed].Is.PMLbody[0].orient);
                        if (orient != iEy) {
                            conta++;
                            BerPML__t& PML_ = berpmlH.nodes[conta - 1];
                            PML_.field = const_cast<RKIND*>(&Hy[i1 * 100 * 100 + j1 * 100 + k1]);
                            berpmlE.orient = orient;
                            switch (orient) {
                                case iEz:
                                    PML_.del = 1.0 / Idze[k1];
                                    PML_.transversalDelta = 1.0 / Idze[k1];
                                    PML_.Plus = const_cast<RKIND*>(&Ex[i1 * 100 * 100 + j1 * 100 + (k1 + 1)]);
                                    PML_.Minu = const_cast<RKIND*>(&Ex[i1 * 100 * 100 + j1 * 100 + k1]);
                                    PML_.gx2 = gm2[jmed];
                                    PML_.minTotal = minz[jmed];
                                    PML_.maxTotal = maxz[jmed];
                                    PML_.posi = k1 + 0.5;
                                    break;
                                case iEx:
                                    PML_.del = 1.0 / Idxe[i1];
                                    PML_.transversalDelta = 1.0 / Idxe[i1];
                                    PML_.Plus = const_cast<RKIND*>(&Ez[(i1 + 1) * 100 * 100 + j1 * 100 + k1]);
                                    PML_.Minu = const_cast<RKIND*>(&Ez[i1 * 100 * 100 + j1 * 100 + k1]);
                                    PML_.gx2 = -gm2[jmed];
                                    PML_.minTotal = minx[jmed];
                                    PML_.maxTotal = maxx[jmed];
                                    PML_.posi = i1 + 0.5;
                                    break;
                                default:
                                    snprintf(buf, BUFSIZE, "Buggy ERROR: In PMLbodies. ");
                                    WarnErrReport(std::string(buf), true);
                                    break;
                            }
                        }
                    }
                }
            }
        }
        // iHz
        for (int k1 = sgg.SINPMLSweep[iHz].ZI; k1 <= sgg.SINPMLSweep[iHz].ZE; ++k1) {
            for (int j1 = sgg.SINPMLSweep[iHz].YI; j1 <= sgg.SINPMLSweep[iHz].YE; ++j1) {
                for (int i1 = sgg.SINPMLSweep[iHz].XI; i1 <= sgg.SINPMLSweep[iHz].XE; ++i1) {
                    int jmed = media.sggMiHz[i1][j1][k1];
                    if (sgg.Med[jmed].Is.PMLbody) {
                        int orient = std::abs(sgg.Med[jmed].Is.PMLbody[0].orient);
                        if (orient != iEz) {
                            conta++;
                            BerPML__t& PML_ = berpmlH.nodes[conta - 1];
                            PML_.field = const_cast<RKIND*>(&Hz[i1 * 100 * 100 + j1 * 100 + k1]);
                            berpmlE.orient = orient;
                            switch (orient) {
                                case iEx:
                                    PML_.del = 1.0 / Idxe[i1];
                                    PML_.transversalDelta = 1.0 / Idxe[i1];
                                    PML_.Plus = const_cast<RKIND*>(&Ey[(i1 + 1) * 100 * 100 + j1 * 100 + k1]);
                                    PML_.Minu = const_cast<RKIND*>(&Ey[i1 * 100 * 100 + j1 * 100 + k1]);
                                    PML_.gx2 = gm2[jmed];
                                    PML_.minTotal = minx[jmed];
                                    PML_.maxTotal = maxx[jmed];
                                    PML_.posi = i1 + 0.5;
                                    break;
                                case iEy:
                                    PML_.del = 1.0 / Idye[j1];
                                    PML_.transversalDelta = 1.0 / Idye[j1];
                                    PML_.Plus = const_cast<RKIND*>(&Ex[i1 * 100 * 100 + (j1 + 1) * 100 + k1]);
                                    PML_.Minu = const_cast<RKIND*>(&Ex[i1 * 100 * 100 + j1 * 100 + k1]);
                                    PML_.gx2 = -gm2[jmed];
                                    PML_.minTotal = miny[jmed];
                                    PML_.maxTotal = maxy[jmed];
                                    PML_.posi = j1 + 0.5;
                                    break;
                                default:
                                    snprintf(buf, BUFSIZE, "Buggy ERROR: In PMLbodies. ");
                                    WarnErrReport(std::string(buf), true);
                                    break;
                            }
                        }
                    }
                }
            }
        }

        calc_pmlbodypar(sgg, eps00, mu00);

        if (!control.resume) {
            for (int i = 0; i < berpmlE.NumNodes; ++i) {
                berpmlE.nodes[i].Psi = 0.0;
            }
            for (int i = 0; i < berpmlH.NumNodes; ++i) {
                berpmlH.nodes[i].Psi = 0.0;
            }
        } else {
            // Reading from file unit 14
            std::ifstream file(14); // Assuming unit 14 is mapped to a file stream
            if (file.is_open()) {
                for (int i = 0; i < berpmlE.NumNodes; ++i) {
                    file >> berpmlE.nodes[i].Psi;
                }
                for (int i = 0; i < berpmlH.NumNodes; ++i) {
                    file >> berpmlH.nodes[i].Psi;
                }
                file.close();
            }
        }
    }

    void calc_pmlbodypar(const SGGFDTDINFO_t& sgg, RKIND eps00, RKIND mu00) {
        eps0 = eps00;
        mu0 = mu00;

        for (int conta = 0; conta < berpmlE.NumNodes; ++conta) {
            BerPML__t& PML_ = berpmlE.nodes[conta];
            RKIND sigmamax = -((std::log(CoeffReflPML) * (PMLorden + 1)) /
                               (2.0 * std::sqrt(mu0 / eps0) * (PML_.del * (PML_.maxTotal - PML_.minTotal) / 2.0)));
            RKIND sigma;
            if (PML_.posi <= (PML_.maxTotal + PML_.minTotal) / 2.0) {
                sigma = sigmamax * (std::abs(PML_.posi - 1.0 * PML_.minTotal) / ((1.0 * PML_.maxTotal - 1.0 * PML_.minTotal) / 2.0)) * PMLorden;
                // Note: The original code uses **PMLorden which is power.
                // The line above is incorrect in C++ translation if we just multiply.
                // Correct translation:
                sigma = sigmamax * std::pow(std::abs(PML_.posi - 1.0 * PML_.minTotal) / ((1.0 * PML_.maxTotal - 1.0 * PML_.minTotal) / 2.0), PMLorden);
            } else {
                sigma = sigmamax * std::pow(std::abs(1.0 * PML_.maxTotal - PML_.posi) / ((1.0 * PML_.maxTotal - 1.0 * PML_.minTotal) / 2.0), PMLorden);
            }
            PML_.P_be = std::exp(-(sigma) * sgg.dt / eps0);
            PML_.P_ce = (PML_.P_be - 1.0) / PML_.transversalDelta;
        }
        for (int conta = 0; conta < berpmlH.NumNodes; ++conta) {
            BerPML__t& PML_ = berpmlH.nodes[conta];
            RKIND sigmamax = -((std::log(CoeffReflPML) * (PMLorden + 1)) /
                               (2.0 * std::sqrt(mu0 / eps0) * (PML_.del * (PML_.maxTotal - PML_.minTotal) / 2.0)));
            RKIND sigma;
            if (PML_.posi <= (PML_.maxTotal + PML_.minTotal) / 2.0) {
                sigma = sigmamax * std::pow(std::abs(PML_.posi - 1.0 * PML_.minTotal) / ((1.0 * PML_.maxTotal - 1.0 * PML_.minTotal) / 2.0), PMLorden);
            } else {
                sigma = sigmamax * std::pow(std::abs(1.0 * PML_.maxTotal - PML_.posi) / ((1.0 * PML_.maxTotal - 1.0 * PML_.minTotal) / 2.0), PMLorden);
            }
            PML_.P_be = std::exp(-(sigma) * sgg.dt / eps0);
            PML_.P_ce = (PML_.P_be - 1.0) / PML_.transversalDelta;
        }
    }

    void AdvancePMLbodyE() {
        for (int conta = 0; conta < berpmlE.NumNodes; ++conta) {
            BerPML__t& PML_ = berpmlE.nodes[conta];
            PML_.Psi = PML_.P_be * PML_.Psi + PML_.P_ce * (*PML_.Plus - *PML_.Minu);
            *PML_.field = *PML_.field + PML_.gx2 * PML_.Psi;
        }
    }

    void AdvancePMLbodyH() {
        for (int conta = 0; conta < berpmlH.NumNodes; ++conta) {
            BerPML__t& PML_ = berpmlH.nodes[conta];
            PML_.Psi = PML_.P_be * PML_.Psi + PML_.P_ce * (*PML_.Plus - *PML_.Minu);
            *PML_.field = *PML_.field - PML_.gx2 * PML_.Psi;
        }
    }

    void StorefieldsPMLbodies() {
        std::ofstream file(14);
        if (file.is_open()) {
            try {
                for (int i = 0; i < berpmlE.NumNodes; ++i) {
                    file << berpmlE.nodes[i].Psi << " ";
                }
                file << std::endl;
                for (int i = 0; i < berpmlH.NumNodes; ++i) {
                    file << berpmlH.nodes[i].Psi << " ";
                }
                file << std::endl;
            } catch (...) {
                print11(0, SEPARADOR + SEPARADOR + SEPARADOR);
                print11(0, "PMLBODIES: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING");
                print11(0, SEPARADOR + SEPARADOR + SEPARADOR);
            }
            file.close();
        }
    }

    void DestroyPMLbodies(SGGFDTDINFO_t& sgg) {
        for (int i = 1; i <= sgg.nummedia; ++i) {
            if (sgg.Med[i].Is.PMLbody && !sgg.Med[i].Is.PML) {
                // In C++, if PMLbody is a vector or pointer, we deallocate it.
                // Assuming sgg.Med[i].PMLbody is a std::vector or similar
                // This part depends on the definition of SGGFDTDINFO_t which is not fully provided.
                // We assume it's handled elsewhere or is a smart pointer.
            }
        }
        berpmlE.nodes.clear();
        berpmlH.nodes.clear();
    }

} // namespace PMLbodies_m