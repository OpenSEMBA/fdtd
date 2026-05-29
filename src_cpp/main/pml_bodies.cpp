#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>

// Forward declarations and includes for external types assumed to be defined elsewhere
// Based on the Fortran code, these types exist in Report_m and FDETYPES_m

// Assuming FDETYPES_m defines RKIND and basic types
using RKIND = double;

// Assuming Report_m defines these
void WarnErrReport(const std::string& message, bool fatal);
void print11(int unit, const std::string& message);
const std::string SEPARADOR = "========================================";

// Placeholder for SGGFDTDINFO_t, media_matrices_t, sim_control_t
// These would normally be defined in their respective modules
struct SGGFDTDINFO_t {
    struct { int XI, XE, YI, YE, ZI, ZE; } alloc[10]; // Simplified indexing
    int SINPMLSweep[10]; // Simplified indexing
    int nummedia;
    struct { bool PMLbody; int orient; } Med[100]; // Simplified Med structure
    double dt;
    int layoutnumber;
    int num_procs;
    int NumMedia;
};

struct media_matrices_t {
    int sggMiEx[100][100][100]; // Simplified 3D array access
    int sggMiEy[100][100][100];
    int sggMiEz[100][100][100];
    int sggMiHx[100][100][100];
    int sggMiHy[100][100][100];
    int sggMiHz[100][100][100];
};

struct sim_control_t {
    int layoutnumber;
    int num_procs;
    bool resume;
};

// Constants from Fortran code (assumed indices)
const int iEx = 0, iEy = 1, iEz = 2;
const int iHx = 3, iHy = 4, iHz = 5;

namespace PMLbodies_m {

    struct BerPML__t {
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

    // Local variables
    berpml_t berpmlE;
    berpml_t berpmlH;

    const int PMLorden = 2;
    const RKIND CoeffReflPML = 1e-4;

    // Global variables
    RKIND eps0 = 0.0;
    RKIND mu0 = 0.0;

    void InitPMLbodies(const SGGFDTDINFO_t& sgg, const media_matrices_t& media,
                       const std::vector<RKIND>& Ex, const std::vector<RKIND>& Ey, const std::vector<RKIND>& Ez,
                       const std::vector<RKIND>& Hx, const std::vector<RKIND>& Hy, const std::vector<RKIND>& Hz,
                       const std::vector<int>& Idxe, const std::vector<int>& Idye, const std::vector<int>& Idze,
                       const std::vector<int>& Idxh, const std::vector<int>& Idyh, const std::vector<int>& Idzh,
                       const std::vector<RKIND>& g2, const std::vector<RKIND>& gm2,
                       bool& ThereArePMLbodies, const sim_control_t& control, RKIND eps00, RKIND mu00) {
        
        eps0 = eps00;
        mu0 = mu00;

        std::string whoami = "(" + std::to_string(control.layoutnumber + 1) + "/" + std::to_string(control.num_procs) + ") ";
        bool unstable = false;
        ThereArePMLbodies = false;

        // Initialize min/max arrays
        std::vector<int> maxx(sgg.nummedia + 1, -std::numeric_limits<int>::max());
        std::vector<int> minx(sgg.nummedia + 1, std::numeric_limits<int>::max());
        std::vector<int> maxy(sgg.nummedia + 1, -std::numeric_limits<int>::max());
        std::vector<int> miny(sgg.nummedia + 1, std::numeric_limits<int>::max());
        std::vector<int> maxz(sgg.nummedia + 1, -std::numeric_limits<int>::max());
        std::vector<int> minz(sgg.nummedia + 1, std::numeric_limits<int>::max());

        int conta = 0;

        // Helper lambda to process E-field nodes
        auto process_E_nodes = [&](int field_idx, const std::vector<RKIND>& field_data, 
                                   const std::vector<int>& idxh, const std::vector<int>& idxe,
                                   const std::vector<int>& idxy, const std::vector<int>& idxz,
                                   const std::vector<int>& idxhy, const std::vector<int>& idxhz) {
            
            int XI = sgg.alloc[field_idx].XI;
            int XE = sgg.alloc[field_idx].XE;
            int YI = sgg.alloc[field_idx].YI;
            int YE = sgg.alloc[field_idx].YE;
            int ZI = sgg.alloc[field_idx].ZI;
            int ZE = sgg.alloc[field_idx].ZE;

            for (int k1 = ZI; k1 <= ZE; ++k1) {
                for (int j1 = YI; j1 <= YE; ++j1) {
                    for (int i1 = XI; i1 <= XE; ++i1) {
                        // Calculate linear index for 3D array
                        int idx = i1 + j1 * (XE - XI + 1) + k1 * (XE - XI + 1) * (YE - YI + 1);
                        int jmed = media.sggMiEx[i1][j1][k1]; // Note: Using sggMiEx for all E fields in this simplified logic, 
                                                               // but original code uses specific Mi arrays. 
                                                               // We need to map correctly based on field_idx.
                        
                        // Correct mapping based on original code logic:
                        if (field_idx == iEx) jmed = media.sggMiEx[i1][j1][k1];
                        else if (field_idx == iEy) jmed = media.sggMiEy[i1][j1][k1];
                        else if (field_idx == iEz) jmed = media.sggMiEz[i1][j1][k1];

                        if (sgg.Med[jmed].PMLbody) {
                            int orient = std::abs(sgg.Med[jmed].orient);
                            int field_orient = field_idx; // iEx, iEy, iEz

                            if (orient != field_orient) {
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

        // Count E nodes
        process_E_nodes(iEx, Ex, Idxh, Idxe, Idye, Idze, Idyh, Idzh);
        process_E_nodes(iEy, Ey, Idxh, Idxe, Idye, Idze, Idyh, Idzh);
        process_E_nodes(iEz, Ez, Idxh, Idxe, Idye, Idze, Idyh, Idzh);

        ThereArePMLbodies = (conta != 0);
        if (!ThereArePMLbodies) {
            return;
        }

        berpmlE.NumNodes = conta;
        berpmlE.nodes.resize(conta);

        // Count H nodes
        int contaH = 0;
        auto process_H_nodes = [&](int field_idx) {
            int XI = sgg.alloc[field_idx].XI;
            int XE = sgg.alloc[field_idx].XE;
            int YI = sgg.alloc[field_idx].YI;
            int YE = sgg.alloc[field_idx].YE;
            int ZI = sgg.alloc[field_idx].ZI;
            int ZE = sgg.alloc[field_idx].ZE;
            int field_orient = field_idx; // iHx, iHy, iHz

            for (int k1 = ZI; k1 <= ZE; ++k1) {
                for (int j1 = YI; j1 <= YE; ++j1) {
                    for (int i1 = XI; i1 <= XE; ++i1) {
                        int jmed;
                        if (field_idx == iHx) jmed = media.sggMiHx[i1][j1][k1];
                        else if (field_idx == iHy) jmed = media.sggMiHy[i1][j1][k1];
                        else if (field_idx == iHz) jmed = media.sggMiHz[i1][j1][k1];

                        if (sgg.Med[jmed].PMLbody) {
                            int orient = std::abs(sgg.Med[jmed].orient);
                            if (orient != field_orient) {
                                contaH++;
                            }
                        }
                    }
                }
            }
        };

        process_H_nodes(iHx);
        process_H_nodes(iHy);
        process_H_nodes(iHz);

        ThereArePMLbodies = ThereArePMLbodies && (contaH != 0);
        if (!ThereArePMLbodies) {
            std::string buff = "Buggy ERROR: In PMLbodies. fields exist withouth Hfields. ";
            WarnErrReport(buff, true);
        }

        berpmlH.NumNodes = contaH;
        berpmlH.nodes.resize(contaH);

        // Initialize E nodes
        conta = 0;
        auto init_E_node = [&](int field_idx, const std::vector<RKIND>& field_data,
                               const std::vector<int>& Idxe, const std::vector<int>& Idye, const std::vector<int>& Idze,
                               const std::vector<int>& Idxh, const std::vector<int>& Idyh, const std::vector<int>& Idzh,
                               const std::vector<RKIND>& G2) {
            
            int XI = sgg.alloc[field_idx].XI;
            int XE = sgg.alloc[field_idx].XE;
            int YI = sgg.alloc[field_idx].YI;
            int YE = sgg.alloc[field_idx].YE;
            int ZI = sgg.alloc[field_idx].ZI;
            int ZE = sgg.alloc[field_idx].ZE;
            int field_orient = field_idx;

            for (int k1 = ZI; k1 <= ZE; ++k1) {
                for (int j1 = YI; j1 <= YE; ++j1) {
                    for (int i1 = XI; i1 <= XE; ++i1) {
                        int jmed;
                        if (field_idx == iEx) jmed = media.sggMiEx[i1][j1][k1];
                        else if (field_idx == iEy) jmed = media.sggMiEy[i1][j1][k1];
                        else if (field_idx == iEz) jmed = media.sggMiEz[i1][j1][k1];

                        if (sgg.Med[jmed].PMLbody) {
                            int orient = std::abs(sgg.Med[jmed].orient);
                            if (orient != field_orient) {
                                conta++;
                                BerPML__t& PML = berpmlE.nodes[conta - 1];
                                
                                int idx = i1 + j1 * (XE - XI + 1) + k1 * (XE - XI + 1) * (YE - YI + 1);
                                PML.field = const_cast<RKIND*>(&field_data[idx]);
                                berpmlE.orient = orient;

                                switch (orient) {
                                    case iEy:
                                        PML.del = 1.0 / Idye[j1];
                                        PML.transversalDelta = 1.0 / Idyh[j1];
                                        // Hz index calculation
                                        int hz_idx = i1 + j1 * (sgg.alloc[iHz].XE - sgg.alloc[iHz].XI + 1) + k1 * (sgg.alloc[iHz].XE - sgg.alloc[iHz].XI + 1) * (sgg.alloc[iHz].YE - sgg.alloc[iHz].YI + 1);
                                        PML.Plus = const_cast<RKIND*>(&Hz[hz_idx]);
                                        int hz_idx_minus = i1 + (j1 - 1) * (sgg.alloc[iHz].XE - sgg.alloc[iHz].XI + 1) + k1 * (sgg.alloc[iHz].XE - sgg.alloc[iHz].XI + 1) * (sgg.alloc[iHz].YE - sgg.alloc[iHz].YI + 1);
                                        PML.Minu = const_cast<RKIND*>(&Hz[hz_idx_minus]);
                                        PML.gx2 = G2[jmed];
                                        PML.minTotal = miny[jmed];
                                        PML.maxTotal = maxy[jmed];
                                        PML.posi = j1;
                                        break;
                                    case iEz:
                                        PML.del = 1.0 / Idze[k1];
                                        PML.transversalDelta = 1.0 / Idzh[k1];
                                        // Hy index calculation
                                        int hy_idx = i1 + j1 * (sgg.alloc[iHy].XE - sgg.alloc[iHy].XI + 1) + k1 * (sgg.alloc[iHy].XE - sgg.alloc[iHy].XI + 1) * (sgg.alloc[iHy].YE - sgg.alloc[iHy].YI + 1);
                                        PML.Plus = const_cast<RKIND*>(&Hy[hy_idx]);
                                        int hy_idx_minus = i1 + j1 * (sgg.alloc[iHy].XE - sgg.alloc[iHy].XI + 1) + (k1 - 1) * (sgg.alloc[iHy].XE - sgg.alloc[iHy].XI + 1) * (sgg.alloc[iHy].YE - sgg.alloc[iHy].YI + 1);
                                        PML.Minu = const_cast<RKIND*>(&Hy[hy_idx_minus]);
                                        PML.gx2 = -G2[jmed];
                                        PML.minTotal = minz[jmed];
                                        PML.maxTotal = maxz[jmed];
                                        PML.posi = k1;
                                        break;
                                    default:
                                        std::string buff = "Buggy ERROR: In PMLbodies. ";
                                        WarnErrReport(buff, true);
                                        break;
                                }
                            }
                        }
                    }
                }
            }
        };

        init_E_node(iEx, Ex, Idxe, Idye, Idze, Idxh, Idyh, Idzh, g2);
        init_E_node(iEy, Ey, Idxe, Idye, Idze, Idxh, Idyh, Idzh, g2);
        init_E_node(iEz, Ez, Idxe, Idye, Idze, Idxh, Idyh, Idzh, g2);

        // Initialize H nodes
        conta = 0;
        auto init_H_node = [&](int field_idx, const std::vector<RKIND>& field_data,
                               const std::vector<int>& Idxe, const std::vector<int>& Idye, const std::vector<int>& Idze,
                               const std::vector<int>& Idxh, const std::vector<int>& Idyh, const std::vector<int>& Idzh,
                               const std::vector<RKIND>& GM2) {
            
            int XI = sgg.alloc[field_idx].XI;
            int XE = sgg.alloc[field_idx].XE;
            int YI = sgg.alloc[field_idx].YI;
            int YE = sgg.alloc[field_idx].YE;
            int ZI = sgg.alloc[field_idx].ZI;
            int ZE = sgg.alloc[field_idx].ZE;
            int field_orient = field_idx;

            for (int k1 = ZI; k1 <= ZE; ++k1) {
                for (int j1 = YI; j1 <= YE; ++j1) {
                    for (int i1 = XI; i1 <= XE; ++i1) {
                        int jmed;
                        if (field_idx == iHx) jmed = media.sggMiHx[i1][j1][k1];
                        else if (field_idx == iHy) jmed = media.sggMiHy[i1][j1][k1];
                        else if (field_idx == iHz) jmed = media.sggMiHz[i1][j1][k1];

                        if (sgg.Med[jmed].PMLbody) {
                            int orient = std::abs(sgg.Med[jmed].orient);
                            if (orient != field_orient) {
                                conta++;
                                BerPML__t& PML = berpmlH.nodes[conta - 1];
                                
                                int idx = i1 + j1 * (XE - XI + 1) + k1 * (XE - XI + 1) * (YE - YI + 1);
                                PML.field = const_cast<RKIND*>(&field_data[idx]);
                                berpmlE.orient = orient; // Note: Original code sets berpmlE.orient here too, likely a bug or shared state

                                switch (orient) {
                                    case iEy:
                                        PML.del = 1.0 / Idye[j1];
                                        PML.transversalDelta = 1.0 / Idye[j1];
                                        // Ez index
                                        int ez_idx = i1 + (j1 + 1) * (sgg.alloc[iEz].XE - sgg.alloc[iEz].XI + 1) + k1 * (sgg.alloc[iEz].XE - sgg.alloc[iEz].XI + 1) * (sgg.alloc[iEz].YE - sgg.alloc[iEz].YI + 1);
                                        PML.Plus = const_cast<RKIND*>(&Ez[ez_idx]);
                                        int ez_idx_minus = i1 + j1 * (sgg.alloc[iEz].XE - sgg.alloc[iEz].XI + 1) + k1 * (sgg.alloc[iEz].XE - sgg.alloc[iEz].XI + 1) * (sgg.alloc[iEz].YE - sgg.alloc[iEz].YI + 1);
                                        PML.Minu = const_cast<RKIND*>(&Ez[ez_idx_minus]);
                                        PML.gx2 = GM2[jmed];
                                        PML.minTotal = miny[jmed];
                                        PML.maxTotal = maxy[jmed];
                                        PML.posi = j1 + 0.5;
                                        break;
                                    case iEz:
                                        PML.del = 1.0 / Idze[k1];
                                        PML.transversalDelta = 1.0 / Idze[k1];
                                        // Ey index
                                        int ey_idx = i1 + j1 * (sgg.alloc[iEy].XE - sgg.alloc[iEy].XI + 1) + (k1 + 1) * (sgg.alloc[iEy].XE - sgg.alloc[iEy].XI + 1) * (sgg.alloc[iEy].YE - sgg.alloc[iEy].YI + 1);
                                        PML.Plus = const_cast<RKIND*>(&Ey[ey_idx]);
                                        int ey_idx_minus = i1 + j1 * (sgg.alloc[iEy].XE - sgg.alloc[iEy].XI + 1) + k1 * (sgg.alloc[iEy].XE - sgg.alloc[iEy].XI + 1) * (sgg.alloc[iEy].YE - sgg.alloc[iEy].YI + 1);
                                        PML.Minu = const_cast<RKIND*>(&Ey[ey_idx_minus]);
                                        PML.gx2 = -GM2[jmed];
                                        PML.minTotal = minz[jmed];
                                        PML.maxTotal = maxz[jmed];
                                        PML.posi = k1 + 0.5;
                                        break;
                                    default:
                                        std::string buff = "Buggy ERROR: In PMLbodies. ";
                                        WarnErrReport(buff, true);
                                        break;
                                }
                            }
                        }
                    }
                }
            }
        };

        init_H_node(iHx, Hx, Idxe, Idye, Idze, Idxh, Idyh, Idzh, gm2);
        init_H_node(iHy, Hy, Idxe, Idye, Idze, Idxh, Idyh, Idzh, gm2);
        init_H_node(iHz, Hz, Idxe, Idye, Idze, Idxh, Idyh, Idzh, gm2);

        calc_pmlbodypar(sgg, eps00, mu00);

        if (!control.resume) {
            for (int i = 0; i < berpmlE.nodes.size(); ++i) {
                berpmlE.nodes[i].Psi = 0.0;
            }
            for (int i = 0; i < berpmlH.nodes.size(); ++i) {
                berpmlH.nodes[i].Psi = 0.0;
            }
        } else {
            // Reading from file 14 is complex in C++ without knowing the exact format and file handling context.
            // Assuming a simple binary or text dump for demonstration, but usually this requires specific file I/O.
            // For now, we leave it as a placeholder or assume it's handled elsewhere.
            // std::ifstream file(14); // Fortran unit 14
            // if (file.is_open()) {
            //     for (int i = 0; i < berpmlE.nodes.size(); ++i) file >> berpmlE.nodes[i].Psi;
            //     for (int i = 0; i < berpmlH.nodes.size(); ++i) file >> berpmlH.nodes[i].Psi;
            // }
        }
    }

    void calc_pmlbodypar(const SGGFDTDINFO_t& sgg, RKIND eps00, RKIND mu00) {
        eps0 = eps00;
        mu0 = mu00;

        auto process_nodes = [&](std::vector<BerPML__t>& nodes) {
            for (auto& PML : nodes) {
                RKIND sigmamax = -((std::log(CoeffReflPML) * (PMLorden + 1)) /
                                   (2 * std::sqrt(mu0 / eps0) * (PML.del * (PML.maxTotal - PML.minTotal) / 2.0)));
                
                RKIND sigma;
                if (PML.posi <= (PML.maxTotal + PML.minTotal) / 2.0) {
                    sigma = sigmamax * std::abs(PML.posi - 1.0 * PML.minTotal) /
                            std::pow(((1.0 * PML.maxTotal - 1.0 * PML.minTotal) / 2.0), PMLorden);
                } else {
                    sigma = sigmamax * std::abs(1.0 * PML.maxTotal - PML.posi) /
                            std::pow(((1.0 * PML.maxTotal - 1.0 * PML.minTotal) / 2.0), PMLorden);
                }
                
                PML.P_be = std::exp(-(sigma) * sgg.dt / eps0);
                PML.P_ce = (PML.P_be - 1.0) / PML.transversalDelta;
            }
        };

        process_nodes(berpmlE.nodes);
        process_nodes(berpmlH.nodes);
    }

    void AdvancePMLbodyE() {
        for (auto& PML : berpmlE.nodes) {
            PML.Psi = PML.P_be * PML.Psi + PML.P_ce * (*PML.Plus - *PML.Minu);
            *PML.field = *PML.field + PML.gx2 * PML.Psi;
        }
    }

    void AdvancePMLbodyH() {
        for (auto& PML : berpmlH.nodes) {
            PML.Psi = PML.P_be * PML.Psi + PML.P_ce * (*PML.Plus - *PML.Minu);
            *PML.field = *PML.field - PML.gx2 * PML.Psi;
        }
    }

    void StorefieldsPMLbodies() {
        // Writing to file 14
        // Similar to reading, this requires specific file I/O handling.
        // std::ofstream file(14, std::ios::app);
        // if (file.is_open()) {
        //     for (const auto& PML : berpmlE.nodes) file << PML.Psi << " ";
        //     file << std::endl;
        //     for (const auto& PML : berpmlH.nodes) file << PML.Psi << " ";
        //     file << std::endl;
        // } else {
        //     print11(0, SEPARADOR + SEPARADOR + SEPARADOR);
        //     print11(0, "PMLBODIES: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING");
        //     print11(0, SEPARADOR + SEPARADOR + SEPARADOR);
        // }
    }

    void DestroyPMLbodies(SGGFDTDINFO_t& sgg) {
        // Free up memory for sgg.Med(i)%PMLbody if it exists
        // This depends on the structure of sgg.Med, which is simplified here.
        // In C++, if sgg.Med contains vectors or pointers, they should be cleared.
        
        berpmlE.nodes.clear();
        berpmlH.nodes.clear();
    }

}