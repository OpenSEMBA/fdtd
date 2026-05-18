#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdint>

// Forward declarations for external types used in the module
// These would typically come from FDETYPES_m and Report_m
struct SGGFDTDINFO_t;
struct limit_t;
struct sim_control_t;

// Placeholder for RKIND, assuming double precision for real(kind=RKIND)
using RKIND = double;

// Placeholder constants/enums that are likely defined in FDETYPES_m or similar
// These need to be resolved in the actual translation context
enum Direction { Down, Up, Left, Right, Back, Front };
enum FieldIndex { iEx, iEy, iEz, iHx, iHy, iHz };

namespace BORDERS_CPML_m {

    constexpr RKIND StaticFrequency = 1.0e14;

    struct xyzlimit_var_t {
        int32_t XI[6];
        int32_t XE[6];
        int32_t YI[6];
        int32_t YE[6];
        int32_t ZI[6];
        int32_t ZE[6];
    };

    struct LR_t {
        std::vector<std::vector<std::vector<RKIND>>> Psi_Exy;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Ezy;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Hxy;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Hzy;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Exyvac;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Ezyvac;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Hxyvac;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Hzyvac;
    };

    struct DU_t {
        std::vector<std::vector<std::vector<RKIND>>> Psi_Eyz;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Exz;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Hyz;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Hxz;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Eyzvac;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Exzvac;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Hyzvac;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Hxzvac;
    };

    struct BF_t {
        std::vector<std::vector<std::vector<RKIND>>> Psi_Ezx;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Eyx;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Hzx;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Hyx;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Ezxvac;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Eyxvac;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Hzxvac;
        std::vector<std::vector<std::vector<RKIND>>> Psi_Hyxvac;
    };

    // Global variables from the module
    std::vector<xyzlimit_var_t> PMLc(6);

    std::vector<LR_t> regLR(6); // left:right mapped to 0:5
    std::vector<Du_t> regDU(6); // down:up mapped to 0:5
    std::vector<BF_t> regBF(6); // back:front mapped to 0:5

    std::vector<std::vector<RKIND>> sig_max;
    std::vector<std::vector<RKIND>> aPar_max;
    std::vector<std::vector<RKIND>> kPar_max;

    std::vector<RKIND> P_ce_x, P_ce_y, P_ce_z;
    std::vector<RKIND> P_be_x, P_be_y, P_be_z;
    std::vector<RKIND> P_cm_x, P_cm_y, P_cm_z;
    std::vector<RKIND> P_bm_x, P_bm_y, P_bm_z;

    std::vector<RKIND> ce_x, ce_y, ce_z;
    std::vector<RKIND> cm_x, cm_y, cm_z;
    std::vector<RKIND> Ice_x, Ice_y, Ice_z;
    std::vector<RKIND> Icm_x, Icm_y, Icm_z;

    RKIND zvac = 0.0;
    RKIND eps0 = 0.0;
    RKIND mu0 = 0.0;

    RKIND alphamaxpar = 0.0;
    int32_t alphaOrden = 0;
    RKIND kappamaxpar = 0.0;

    std::vector<limit_t> SINPML_fullsize(6);

    std::vector<RKIND> dxe, dye, dze;
    std::vector<RKIND> dxh, dyh, dzh;

    void InitCPMLBorders(const SGGFDTDINFO_t& sgg,
                         const std::vector<RKIND>& temp_dxe,
                         const std::vector<RKIND>& temp_dye,
                         const std::vector<RKIND>& temp_dze,
                         const std::vector<RKIND>& temp_dxh,
                         const std::vector<RKIND>& temp_dyh,
                         const std::vector<RKIND>& temp_dzh,
                         std::vector<int32_t>& Idxe,
                         std::vector<int32_t>& Idye,
                         std::vector<int32_t>& Idze,
                         std::vector<int32_t>& Idxh,
                         std::vector<int32_t>& Idyh,
                         std::vector<int32_t>& Idzh,
                         const std::vector<limit_t>& temp_SINPML_fullsize,
                         const sim_control_t& control,
                         bool& ThereArePMLBorders,
                         RKIND eps00,
                         RKIND mu00) {
        
        eps0 = eps00;
        mu0 = mu00;
        zvac = std::sqrt(mu0 / eps0);
        SINPML_fullsize = temp_SINPML_fullsize;
        alphamaxpar = control.alphamaxpar;
        alphaOrden = control.alphaOrden;
        kappamaxpar = control.kappamaxpar;

        // Resize vectors based on sgg allocations
        // Note: Fortran arrays are 1-based or custom indexed. 
        // We assume std::vector is resized to hold the range [XI, XE].
        // To simplify, we'll resize to XE - XI + 1 and access with offset.
        
        int32_t hxi = sgg.ALLOC(iHx).XI;
        int32_t hxe = sgg.ALLOC(iHx).XE;
        int32_t hyi = sgg.ALLOC(iHy).YI;
        int32_t hye = sgg.ALLOC(iHy).YE;
        int32_t hzi = sgg.ALLOC(iHz).ZI;
        int32_t hze = sgg.ALLOC(iHz).ZE;
        int32_t exi = sgg.ALLOC(iEx).XI;
        int32_t exe = sgg.ALLOC(iEx).XE;
        int32_t eyi = sgg.ALLOC(iEy).YI;
        int32_t eye = sgg.ALLOC(iEy).YE;
        int32_t ezi = sgg.ALLOC(iEz).ZI;
        int32_t eze = sgg.ALLOC(iEz).ZE;

        dxe.resize(hxe - hxi + 1);
        dye.resize(hye - hyi + 1);
        dze.resize(hze - hzi + 1);
        dxh.resize(exe - exi + 1);
        dyh.resize(eye - eyi + 1);
        dzh.resize(eze - ezi + 1);

        // Copy data
        for (int32_t i = 0; i < dxe.size(); ++i) dxe[i] = temp_dxe[hxi + i];
        for (int32_t i = 0; i < dye.size(); ++i) dye[i] = temp_dye[hyi + i];
        for (int32_t i = 0; i < dze.size(); ++i) dze[i] = temp_dze[hzi + i];
        for (int32_t i = 0; i < dxh.size(); ++i) dxh[i] = temp_dxh[exi + i];
        for (int32_t i = 0; i < dyh.size(); ++i) dyh[i] = temp_dyh[eyi + i];
        for (int32_t i = 0; i < dzh.size(); ++i) dzh[i] = temp_dzh[ezi + i];

        ThereArePMLBorders = false;
        if (sgg.Border.IsBackPML || sgg.Border.IsFrontPML || sgg.Border.IsLeftPML || 
            sgg.Border.IsRightPML || sgg.Border.IsUpPML || sgg.Border.IsDownPML) {
            ThereArePMLBorders = true;
        }

        if (!ThereArePMLBorders) return;

        // Find limits
        for (int32_t field = 0; field < 6; ++field) { // iEx to iHz mapped to 0-5
            // Down
            PMLc[field].XI[Down] = sgg.Sweep(field).XI;
            PMLc[field].XE[Down] = sgg.Sweep(field).XE;
            PMLc[field].YI[Down] = sgg.Sweep(field).YI;
            PMLc[field].YE[Down] = sgg.Sweep(field).YE;
            PMLc[field].ZI[Down] = sgg.Sweep(field).ZI;
            PMLc[field].ZE[Down] = std::min(SINPML_fullsize[field].ZI - 1, sgg.Sweep(field).ZE);

            // Up
            PMLc[field].XI[Up] = sgg.Sweep(field).XI;
            PMLc[field].XE[Up] = sgg.Sweep(field).XE;
            PMLc[field].YI[Up] = sgg.Sweep(field).YI;
            PMLc[field].YE[Up] = sgg.Sweep(field).YE;
            PMLc[field].ZI[Up] = std::max(SINPML_fullsize[field].ZE + 1, sgg.Sweep(field).ZI);
            PMLc[field].ZE[Up] = sgg.Sweep(field).ZE;

            // Left
            PMLc[field].XI[Left] = sgg.Sweep(field).XI;
            PMLc[field].XE[Left] = sgg.Sweep(field).XE;
            PMLc[field].YI[Left] = sgg.Sweep(field).YI;
            PMLc[field].YE[Left] = std::min(SINPML_fullsize[field].YI - 1, sgg.Sweep(field).YE);
            PMLc[field].ZI[Left] = sgg.Sweep(field).ZI;
            PMLc[field].ZE[Left] = sgg.Sweep(field).ZE;

            // Right
            PMLc[field].XI[Right] = sgg.Sweep(field).XI;
            PMLc[field].XE[Right] = sgg.Sweep(field).XE;
            PMLc[field].YI[Right] = std::max(SINPML_fullsize[field].YE + 1, sgg.Sweep(field).YI);
            PMLc[field].YE[Right] = sgg.Sweep(field).YE;
            PMLc[field].ZI[Right] = sgg.Sweep(field).ZI;
            PMLc[field].ZE[Right] = sgg.Sweep(field).ZE;

            // Back
            PMLc[field].XI[Back] = sgg.Sweep(field).XI;
            PMLc[field].XE[Back] = std::min(SINPML_fullsize[field].XI - 1, sgg.Sweep(field).XE);
            PMLc[field].YI[Back] = sgg.Sweep(field).YI;
            PMLc[field].YE[Back] = sgg.Sweep(field).YE;
            PMLc[field].ZI[Back] = sgg.Sweep(field).ZI;
            PMLc[field].ZE[Back] = sgg.Sweep(field).ZE;

            // Front
            PMLc[field].XI[Front] = std::max(SINPML_fullsize[field].XE + 1, sgg.Sweep(field).XI);
            PMLc[field].XE[Front] = sgg.Sweep(field).XE;
            PMLc[field].YI[Front] = sgg.Sweep(field).YI;
            PMLc[field].YE[Front] = sgg.Sweep(field).YE;
            PMLc[field].ZI[Front] = sgg.Sweep(field).ZI;
            PMLc[field].ZE[Front] = sgg.Sweep(field).ZE;
        }

        // Allocate PML constants
        sig_max.resize(3, std::vector<RKIND>(2));
        aPar_max.resize(3, std::vector<RKIND>(2));
        kPar_max.resize(3, std::vector<RKIND>(2));

        // Allocate field pointers
        P_ce_x.resize(hxe - hxi + 1);
        P_ce_y.resize(hye - hyi + 1);
        P_ce_z.resize(hze - hzi + 1);
        P_be_x.resize(hxe - hxi + 1);
        P_be_y.resize(hye - hyi + 1);
        P_be_z.resize(hze - hzi + 1);
        P_cm_x.resize(hxe - hxi + 1);
        P_cm_y.resize(hye - hyi + 1);
        P_cm_z.resize(hze - hzi + 1);
        P_bm_x.resize(hxe - hxi + 1);
        P_bm_y.resize(hye - hyi + 1);
        P_bm_z.resize(hze - hzi + 1);

        ce_x.resize(hxe - hxi + 1);
        ce_y.resize(hye - hyi + 1);
        ce_z.resize(hze - hzi + 1);
        cm_x.resize(hxe - hxi + 1);
        cm_y.resize(hye - hyi + 1);
        cm_z.resize(hze - hzi + 1);

        Ice_x.resize(hxe - hxi + 1);
        Ice_y.resize(hye - hyi + 1);
        Ice_z.resize(hze - hzi + 1);
        Icm_x.resize(hxe - hxi + 1);
        Icm_y.resize(hye - hyi + 1);
        Icm_z.resize(hze - hzi + 1);

        // Initialize to 0
        std::fill(ce_x.begin(), ce_x.end(), 0.0);
        std::fill(ce_y.begin(), ce_y.end(), 0.0);
        std::fill(ce_z.begin(), ce_z.end(), 0.0);
        std::fill(cm_x.begin(), cm_x.end(), 0.0);
        std::fill(cm_y.begin(), cm_y.end(), 0.0);
        std::fill(cm_z.begin(), cm_z.end(), 0.0);
        std::fill(Ice_x.begin(), Ice_x.end(), 0.0);
        std::fill(Ice_y.begin(), Ice_y.end(), 0.0);
        std::fill(Ice_z.begin(), Ice_z.end(), 0.0);
        std::fill(Icm_x.begin(), Icm_x.end(), 0.0);
        std::fill(Icm_y.begin(), Icm_y.end(), 0.0);
        std::fill(Icm_z.begin(), Icm_z.end(), 0.0);

        // Depth information matrices
        for (int32_t i = hxi; i <= hxe; ++i) {
            int32_t idx = i - hxi;
            if (i <= SINPML_fullsize[iHx].XI && sgg.PML.NumLayers(1, 1) != 0) {
                ce_x[idx] = 1.0 * (SINPML_fullsize[iHx].XI - i) / sgg.PML.NumLayers(1, 1);
                Ice_x[idx] = 1.0 * (sgg.PML.NumLayers(1, 1) - (SINPML_fullsize[iHx].XI - i)) / sgg.PML.NumLayers(1, 1);
            } else if (i >= SINPML_fullsize[iHx].XE && sgg.PML.NumLayers(1, 2) != 0) {
                ce_x[idx] = 1.0 * (i - SINPML_fullsize[iHx].XE) / sgg.PML.NumLayers(1, 2);
                Ice_x[idx] = 1.0 * (sgg.PML.NumLayers(1, 2) - (i - SINPML_fullsize[iHx].XE)) / sgg.PML.NumLayers(1, 2);
            } else {
                ce_x[idx] = 0.0;
                Ice_x[idx] = 0.0;
            }
        }

        for (int32_t i = hxi; i <= hxe; ++i) {
            int32_t idx = i - hxi;
            if (i <= SINPML_fullsize[iHx].XI - 1 && sgg.PML.NumLayers(1, 1) != 0) {
                cm_x[idx] = 1.0 * (SINPML_fullsize[iHx].XI - (i + 0.5)) / sgg.PML.NumLayers(1, 1);
                Icm_x[idx] = 1.0 * (sgg.PML.NumLayers(1, 1) - (SINPML_fullsize[iHx].XI - (i + 0.5))) / sgg.PML.NumLayers(1, 1);
            } else if (i >= SINPML_fullsize[iHx].XE && sgg.PML.NumLayers(1, 2) != 0) {
                cm_x[idx] = 1.0 * (i - SINPML_fullsize[iHx].XE + 0.5) / sgg.PML.NumLayers(1, 2);
                Icm_x[idx] = 1.0 * (sgg.PML.NumLayers(1, 2) - (i - SINPML_fullsize[iHx].XE + 0.5)) / sgg.PML.NumLayers(1, 2);
            } else {
                cm_x[idx] = 0.0;
                Icm_x[idx] = 0.0;
            }
        }

        for (int32_t j = hyi; j <= hye; ++j) {
            int32_t idx = j - hyi;
            if (j <= SINPML_fullsize[iHy].YI && sgg.PML.NumLayers(2, 1) != 0) {
                ce_y[idx] = 1.0 * (SINPML_fullsize[iHy].YI - j) / sgg.PML.NumLayers(2, 1);
                Ice_y[idx] = 1.0 * (sgg.PML.NumLayers(2, 1) - (SINPML_fullsize[iHy].YI - j)) / sgg.PML.NumLayers(2, 1);
            } else if (j >= SINPML_fullsize[iHy].YE && sgg.PML.NumLayers(2, 2) != 0) {
                ce_y[idx] = 1.0 * (j - SINPML_fullsize[iHy].YE) / sgg.PML.NumLayers(2, 2);
                Ice_y[idx] = 1.0 * (sgg.PML.NumLayers(2, 2) - (j - SINPML_fullsize[iHy].YE)) / sgg.PML.NumLayers(2, 2);
            } else {
                ce_y[idx] = 0.0;
                Ice_y[idx] = 0.0;
            }
        }

        for (int32_t j = hyi; j <= hye; ++j) {
            int32_t idx = j - hyi;
            if (j <= SINPML_fullsize[iHy].YI - 1 && sgg.PML.NumLayers(2, 1) != 0) {
                cm_y[idx] = 1.0 * (SINPML_fullsize[iHy].YI - (j + 0.5)) / sgg.PML.NumLayers(2, 1);
                Icm_y[idx] = 1.0 * (sgg.PML.NumLayers(2, 1) - (SINPML_fullsize[iHy].YI - (j + 0.5))) / sgg.PML.NumLayers(2, 1);
            } else if (j >= SINPML_fullsize[iHy].YE && sgg.PML.NumLayers(2, 2) != 0) {
                cm_y[idx] = 1.0 * (j - SINPML_fullsize[iHy].YE + 0.5) / sgg.PML.NumLayers(2, 2);
                Icm_y[idx] = 1.0 * (sgg.PML.NumLayers(2, 2) - (j - SINPML_fullsize[iHy].YE + 0.5)) / sgg.PML.NumLayers(2, 2);
            } else {
                cm_y[idx] = 0.0;
                Icm_y[idx] = 0.0;
            }
        }

        for (int32_t k = hzi; k <= hze; ++k) {
            int32_t idx = k - hzi;
            if (k <= SINPML_fullsize[iHz].ZI && sgg.PML.NumLayers(3, 1) != 0) {
                ce_z[idx] = 1.0 * (SINPML_fullsize[iHz].ZI - k) / sgg.PML.NumLayers(3, 1);
                Ice_z[idx] = 1.0 * (sgg.PML.NumLayers(3, 1) - (SINPML_fullsize[iHz].ZI - k)) / sgg.PML.NumLayers(3, 1);
            } else if (k >= SINPML_fullsize[iHz].ZE && sgg.PML.NumLayers(3, 2) != 0) {
                ce_z[idx] = 1.0 * (k - SINPML_fullsize[iHz].ZE) / sgg.PML.NumLayers(3, 2);
                Ice_z[idx] = 1.0 * (sgg.PML.NumLayers(3, 2) - (k - SINPML_fullsize[iHz].ZE)) / sgg.PML.NumLayers(3, 2);
            } else {
                ce_z[idx] = 0.0;
                Ice_z[idx] = 0.0;
            }
        }

        for (int32_t k = hzi; k <= hze; ++k) {
            int32_t idx = k - hzi;
            if (k <= SINPML_fullsize[iHz].ZI - 1 && sgg.PML.NumLayers(3, 1) != 0) {
                cm_z[idx] = 1.0 * (SINPML_fullsize[iHz].ZI - (k + 0.5)) / sgg.PML.NumLayers(3, 1);
                Icm_z[idx] = 1.0 * (sgg.PML.NumLayers(3, 1) - (SINPML_fullsize[iHz].ZI - (k + 0.5))) / sgg.PML.NumLayers(3, 1);
            } else if (k >= SINPML_fullsize[iHz].ZE && sgg.PML.NumLayers(3, 2) != 0) {
                cm_z[idx] = 1.0 * (k - SINPML_fullsize[iHz].ZE + 0.5) / sgg.PML.NumLayers(3, 2);
                Icm_z[idx] = 1.0 * (sgg.PML.NumLayers(3, 2) - (k - SINPML_fullsize[iHz].ZE + 0.5)) / sgg.PML.NumLayers(3, 2);
            } else {
                cm_z[idx] = 0.0;
                Icm_z[idx] = 0.0;
            }
        }

        calc_cpmlconstants(sgg, Idxe, Idye, Idze, Idxh, Idyh, Idzh, eps0, mu0);

        // Fake coms and ends
        if (!sgg.Border.IsDownPML) {
            for (int32_t f = 0; f < 6; ++f) PMLc[f].ZI[Down] = PMLc[f].ZE[Down] + 100;
        }
        if (!sgg.Border.IsUpPML) {
            for (int32_t f = 0; f < 6; ++f) PMLc[f].ZI[Up] = PMLc[f].ZE[Up] + 100;
        }
        if (!sgg.Border.IsLeftPML) {
            for (int32_t f = 0; f < 6; ++f) PMLc[f].ZI[Left] = PMLc[f].ZE[Left] + 100;
        }
        if (!sgg.Border.IsRightPML) {
            for (int32_t f = 0; f < 6; ++f) PMLc[f].ZI[Right] = PMLc[f].ZE[Right] + 100;
        }
        if (!sgg.Border.IsFrontPML) {
            for (int32_t f = 0; f < 6; ++f) PMLc[f].ZI[Front] = PMLc[f].ZE[Front] + 100;
        }
        if (!sgg.Border.IsBackPML) {
            for (int32_t f = 0; f < 6; ++f) PMLc[f].ZI[Back] = PMLc[f].ZE[Back] + 100;
        }

        // PML Field component matrix allocation
        // Mapping regions: left=0, right=1, down=2, up=3, back=4, front=5
        // Note: The original code uses enums left:right etc. We assume they map to 0-5 indices for regLR/DU/BF
        
        // Left/Right regions (regLR)
        for (int32_t region = 0; region < 2; ++region) { // 0:Left, 1:Right
            int32_t xi_ex = PMLc[iEx].XI[region];
            int32_t xe_ex = PMLc[iEx].XE[region];
            int32_t yi_ex = PMLc[iEx].YI[region];
            int32_t ye_ex = PMLc[iEx].YE[region];
            int32_t zi_ex = PMLc[iEx].ZI[region];
            int32_t ze_ex = PMLc[iEx].ZE[region];
            
            regLR[region].Psi_Exy.resize(ze_ex - zi_ex + 1, std::vector<std::vector<RKIND>>(ye_ex - yi_ex + 1, std::vector<RKIND>(xe_ex - xi_ex + 1, 0.0)));
            regLR[region].Psi_Ezy.resize(PMLc[iEz].ZE[region] - PMLc[iEz].ZI[region] + 1, std::vector<std::vector<RKIND>>(PMLc[iEz].YE[region] - PMLc[iEz].YI[region] + 1, std::vector<RKIND>(PMLc[iEz].XE[region] - PMLc[iEz].XI[region] + 1, 0.0)));
            regLR[region].Psi_Hxy.resize(PMLc[iHx].ZE[region] - PMLc[iHx].ZI[region] + 1, std::vector<std::vector<RKIND>>(PMLc[iHx].YE[region] - PMLc[iHx].YI[region] + 1, std::vector<RKIND>(PMLc[iHx].XE[region] - PMLc[iHx].XI[region] + 1, 0.0)));
            regLR[region].Psi_Hzy.resize(PMLc[iHz].ZE[region] - PMLc[iHz].ZI[region] + 1, std::vector<std::vector<RKIND>>(PMLc[iHz].YE[region] - PMLc[iHz].YI[region] + 1, std::vector<RKIND>(PMLc[iHz].XE[region] - PMLc[iHz].XI[region] + 1, 0.0)));

            if (!control.resume) {
                // Already initialized to 0 in resize
            } else {
                // Read from file 14
                // Note: File I/O needs to be handled carefully. Assuming file 14 is open.
                std::ifstream file14("14", std::ios::binary); // Placeholder
                if (!file14) {
                    std::cerr << "Error opening file 14 for reading" << std::endl;
                    return;
                }

                for (int32_t k = zi_ex; k <= ze_ex; ++k) {
                    for (int32_t j = yi_ex; j <= ye_ex; ++j) {
                        for (int32_t i = xi_ex; i <= xe_ex; ++i) {
                            file14.read(reinterpret_cast<char*>(&regLR[region].Psi_Exy[k - zi_ex][j - yi_ex][i - xi_ex]), sizeof(RKIND));
                        }
                    }
                }
                // Similar reads for other fields...
                // Due to complexity and lack of full context for file format, this part is simplified.
                // In a real translation, precise binary read logic matching Fortran WRITE/READ is required.
            }
        }

        // Down/Up regions (regDU)
        for (int32_t region = 2; region < 4; ++region) { // 2:Down, 3:Up
            int32_t xi_ey = PMLc[iEy].XI[region];
            int32_t xe_ey = PMLc[iEy].XE[region];
            int32_t yi_ey = PMLc[iEy].YI[region];
            int32_t ye_ey = PMLc[iEy].YE[region];
            int32_t zi_ey = PMLc[iEy].ZI[region];
            int32_t ze_ey = PMLc[iEy].ZE[region];

            regDU[region].Psi_Eyz.resize(ze_ey - zi_ey + 1, std::vector<std::vector<RKIND>>(ye_ey - yi_ey + 1, std::vector<RKIND>(xe_ey - xi_ey + 1, 0.0)));
            regDU[region].Psi_Exz.resize(PMLc[iEx].ZE[region] - PMLc[iEx].ZI[region] + 1, std::vector<std::vector<RKIND>>(PMLc[iEx].YE[region] - PMLc[iEx].YI[region] + 1, std::vector<RKIND>(PMLc[iEx].XE[region] - PMLc[iEx].XI[region] + 1, 0.0)));
            regDU[region].Psi_Hyz.resize(PMLc[iHy].ZE[region] - PMLc[iHy].ZI[region] + 1, std::vector<std::vector<RKIND>>(PMLc[iHy].YE[region] - PMLc[iHy].YI[region] + 1, std::vector<RKIND>(PMLc[iHy].XE[region] - PMLc[iHy].XI[region] + 1, 0.0)));
            regDU[region].Psi_Hxz.resize(PMLc[iHx].ZE[region] - PMLc[iHx].ZI[region] + 1, std::vector<std::vector<RKIND>>(PMLc[iHx].YE[region] - PMLc[iHx].YI[region] + 1, std::vector<RKIND>(PMLc[iHx].XE[region] - PMLc[iHx].XI[region] + 1, 0.0)));

            if (!control.resume) {
                // Already initialized to 0
            } else {
                // Read from file 14
            }
        }

        // Back/Front regions (regBF)
        for (int32_t region = 4; region < 6; ++region) { // 4:Back, 5:Front
            int32_t xi_ez = PMLc[iEz].XI[region];
            int32_t xe_ez = PMLc[iEz].XE[region];
            int32_t yi_ez = PMLc[iEz].YI[region];
            int32_t ye_ez = PMLc[iEz].YE[region];
            int32_t zi_ez = PMLc[iEz].ZI[region];
            int32_t ze_ez = PMLc[iEz].ZE[region];

            regBF[region].Psi_Ezx.resize(ze_ez - zi_ez + 1, std::vector<std::vector<RKIND>>(ye_ez - yi_ez + 1, std::vector<RKIND>(xe_ez - xi_ez + 1, 0.0)));
            regBF[region].Psi_Eyx.resize(PMLc[iEy].ZE[region] - PMLc[iEy].ZI[region] + 1, std::vector<std::vector<RKIND>>(PMLc[iEy].YE[region] - PMLc[iEy].YI[region] + 1, std::vector<RKIND>(PMLc[iEy].XE[region] - PMLc[iEy].XI[region] + 1, 0.0)));
            regBF[region].Psi_Hzx.resize(PMLc[iHz].ZE[region] - PMLc[iHz].ZI[region] + 1, std::vector<std::vector<RKIND>>(PMLc[iHz].YE[region] - PMLc[iHz].YI[region] + 1, std::vector<RKIND>(PMLc[iHz].XE[region] - PMLc[iHz].XI[region] + 1, 0.0)));
            regBF[region].Psi_Hyx.resize(PMLc[iHy].ZE[region] - PMLc[iHy].ZI[region] + 1, std::vector<std::vector<RKIND>>(PMLc[iHy].YE[region] - PMLc[iHy].YI[region] + 1, std::vector<RKIND>(PMLc[iHy].XE[region] - PMLc[iHy].XI[region] + 1, 0.0)));

            if (!control.resume) {
                // Already initialized to 0
            } else {
                // Read from file 14
            }
        }
    }

    void StoreFieldsCPMLBorders() {
        // Placeholder for file writing logic
        // Similar structure to InitCPMLBorders but writing to file 14
        std::ofstream file14("14", std::ios::binary);
        if (!file14) {
            std::cerr << "Error opening file 14 for writing" << std::endl;
            return;
        }

        // Left/Right regions
        for (int32_t region = 0; region < 2; ++region) {
            int32_t xi_ex = PMLc[iEx].XI[region];
            int32_t xe_ex = PMLc[iEx].XE[region];
            int32_t yi_ex = PMLc[iEx].YI[region];
            int32_t ye_ex = PMLc[iEx].YE[region];
            int32_t zi_ex = PMLc[iEx].ZI[region];
            int32_t ze_ex = PMLc[iEx].ZE[region];

            for (int32_t k = zi_ex; k <= ze_ex; ++k) {
                for (int32_t j = yi_ex; j <= ye_ex; ++j) {
                    for (int32_t i = xi_ex; i <= xe_ex; ++i) {
                        file14.write(reinterpret_cast<const char*>(&regLR[region].Psi_Exy[k - zi_ex][j - yi_ex][i - xi_ex]), sizeof(RKIND));
                    }
                }
            }
            // Similar writes for other fields...
        }

        // Down/Up regions
        for (int32_t region = 2; region < 4; ++region) {
            // ...
        }

        // Back/Front regions
        for (int32_t region = 4; region < 6; ++region) {
            // ...
        }
    }

    void calc_cpmlconstants(const SGGFDTDINFO_t& sgg,
                            const std::vector<int32_t>& Idxe,
                            const std::vector<int32_t>& Idye,
                            const std::vector<int32_t>& Idze,
                            const std::vector<int32_t>& Idxh,
                            const std::vector<int32_t>& Idyh,
                            const std::vector<int32_t>& Idzh,
                            RKIND eps0,
                            RKIND mu0) {
        // Implementation placeholder
    }

    void AdvanceelectricCPML() {
        // Implementation placeholder
    }

    void AdvanceMagneticCPML() {
        // Implementation placeholder
    }

    void DestroyCPMLBorders() {
        // Implementation placeholder
    }

    void AdvanceelectricCPML_freespace() {
        // Implementation placeholder
    }

    void AdvanceMagneticCPML_freespace() {
        // Implementation placeholder
    }

} // namespace BORDERS_CPML_m

}
         }
      }

      goto 635;
634:
      print11(0, SEPARADOR + separador + separador);
      print11(0, "BORDERSPML: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING");
      print11(0, SEPARADOR + separador + separador);
635:
      return;
   } // subroutine StoreFieldsCPMLBorders

   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!  Free-up memory
   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   void DestroyCPMLBorders() {
      int region;
      if (sig_max != nullptr) {
         delete[] sig_max;
         delete[] aPar_max;
         delete[] kPar_max;
         sig_max = nullptr;
         aPar_max = nullptr;
         kPar_max = nullptr;
      }
      if (P_ce_x != nullptr) {
         delete[] P_ce_x;
         delete[] P_ce_y;
         delete[] P_ce_z;
         delete[] P_be_x;
         delete[] P_be_y;
         delete[] P_be_z;
         delete[] P_cm_x;
         delete[] P_cm_y;
         delete[] P_cm_z;
         delete[] P_bm_x;
         delete[] P_bm_y;
         delete[] P_bm_z;
         P_ce_x = nullptr;
         P_ce_y = nullptr;
         P_ce_z = nullptr;
         P_be_x = nullptr;
         P_be_y = nullptr;
         P_be_z = nullptr;
         P_cm_x = nullptr;
         P_cm_y = nullptr;
         P_cm_z = nullptr;
         P_bm_x = nullptr;
         P_bm_y = nullptr;
         P_bm_z = nullptr;
      }
      if (ce_x != nullptr) {
         delete[] ce_x;
         delete[] ce_y;
         delete[] ce_z;
         delete[] cm_x;
         delete[] cm_y;
         delete[] cm_z;
         delete[] Ice_x;
         delete[] Ice_y;
         delete[] Ice_z;
         delete[] Icm_x;
         delete[] Icm_y;
         delete[] Icm_z;
         ce_x = nullptr;
         ce_y = nullptr;
         ce_z = nullptr;
         cm_x = nullptr;
         cm_y = nullptr;
         cm_z = nullptr;
         Ice_x = nullptr;
         Ice_y = nullptr;
         Ice_z = nullptr;
         Icm_x = nullptr;
         Icm_y = nullptr;
         Icm_z = nullptr;
      }

      for (REGION = left; REGION <= right; ++REGION) {
         if (regLR[REGION].Psi_Exy != nullptr) {
            delete[] regLR[REGION].Psi_Exy;
            delete[] regLR[REGION].Psi_Ezy;
            delete[] regLR[REGION].Psi_Hxy;
            delete[] regLR[REGION].Psi_Hzy;
            regLR[REGION].Psi_Exy = nullptr;
            regLR[REGION].Psi_Ezy = nullptr;
            regLR[REGION].Psi_Hxy = nullptr;
            regLR[REGION].Psi_Hzy = nullptr;
         }
      }
      for (REGION = down; REGION <= up; ++REGION) {
         if (regDU[REGION].Psi_Eyz != nullptr) {
            delete[] regDU[REGION].Psi_Eyz;
            delete[] regDU[REGION].Psi_Exz;
            delete[] regDU[REGION].Psi_Hyz;
            delete[] regDU[REGION].Psi_Hxz;
            regDU[REGION].Psi_Eyz = nullptr;
            regDU[REGION].Psi_Exz = nullptr;
            regDU[REGION].Psi_Hyz = nullptr;
            regDU[REGION].Psi_Hxz = nullptr;
         }
      }
      for (REGION = back; REGION <= front; ++REGION) {
         if (regBF[REGION].Psi_Ezx != nullptr) {
            delete[] regBF[REGION].Psi_Ezx;
            delete[] regBF[REGION].Psi_Eyx;
            delete[] regBF[REGION].Psi_Hzx;
            delete[] regBF[REGION].Psi_Hyx;
            regBF[REGION].Psi_Ezx = nullptr;
            regBF[REGION].Psi_Eyx = nullptr;
            regBF[REGION].Psi_Hzx = nullptr;
            regBF[REGION].Psi_Hyx = nullptr;
         }
      }
      if (dxe != nullptr) {
         delete[] dxe;
         delete[] dye;
         delete[] dze;
         delete[] dxh;
         delete[] dyh;
         delete[] dzh;
         dxe = nullptr;
         dye = nullptr;
         dze = nullptr;
         dxh = nullptr;
         dyh = nullptr;
         dzh = nullptr;
      }
      return;
   } // subroutine DestroyCPMLBorders

   // **************************************************************************************************
   void AdvanceelectricCPML(int NumMedia, const bounds_t& b,
                            const std::vector<std::vector<std::vector<int>>>& sggMiEx,
                            const std::vector<std::vector<std::vector<int>>>& sggMiEy,
                            const std::vector<std::vector<std::vector<int>>>& sggMiEz,
                            const std::vector<double>& g2,
                            const std::vector<std::vector<std::vector<double>>>& Hx,
                            const std::vector<std::vector<std::vector<double>>>& Hy,
                            const std::vector<std::vector<std::vector<double>>>& Hz,
                            std::vector<std::vector<std::vector<double>>>& Ex,
                            std::vector<std::vector<std::vector<double>>>& Ey,
                            std::vector<std::vector<std::vector<double>>>& Ez) {
      // ---------------------------> inputs <----------------------------------------------------------
      // integer, intent( IN) :: NumMedia
      // type( bounds_t), intent( IN) :: b
      // --->
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEx%NX-1, 0 :  b%sggMiEx%NY-1, 0 :  b%sggMiEx%NZ-1), intent( IN) :: sggMiEx
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEy%NX-1, 0 :  b%sggMiEy%NY-1, 0 :  b%sggMiEy%NZ-1), intent( IN) :: sggMiEy
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEz%NX-1, 0 :  b%sggMiEz%NY-1, 0 :  b%sggMiEz%NZ-1), intent( IN) :: sggMiEz
      // --->
      // real(kind = RKIND), dimension( 0 :  NumMedia), intent( IN) :: g2
      // --->
      // real(kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( IN) :: Hx
      // real(kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( IN) :: Hy
      // real(kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( IN) :: Hz
      // ---------------------------> inputs/outputs <--------------------------------------------------
      // real(kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( INOUT) :: Ex
      // real(kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( INOUT) :: Ey
      // real(kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( INOUT) :: Ez
      // ---------------------------> variables locales <-----------------------------------------------
      int REGION, i, j, k, medio, i_m, j_m, k_m;
      // ---------------------------> empieza AdvanceelectricCPML <-------------------------------------

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = left;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEx].ZI[REGION]; k <= PMLc[iEx].ZE[REGION]; ++k) {
         k_m = k - b.Ex.ZI;
         for (j = PMLc[iEx].YI[REGION]; j <= PMLc[iEx].YE[REGION]; ++j) {
            j_m = j - b.Ex.YI;
            for (i = PMLc[iEx].XI[REGION]; i <= PMLc[iEx].XE[REGION]; ++i) {
               i_m = i - b.Ex.XI;
               // --->
               medio = sggMiEx[i_m][j_m][k_m];
               regLR[REGION].Psi_Exy[i][j][k] = P_be_y[j] * regLR[REGION].Psi_Exy[i][j][k] +
                                                 (Hz[i_m][j_m][k_m] - Hz[i_m][j_m - 1][k_m]) * P_ce_y[j];
               Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m] + G2[medio] * regLR[REGION].Psi_Exy[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEz].ZI[REGION]; k <= PMLc[iEz].ZE[REGION]; ++k) {
         k_m = k - b.Ez.ZI;
         for (j = PMLc[iEz].YI[REGION]; j <= PMLc[iEz].YE[REGION]; ++j) {
            j_m = j - b.Ez.YI;
            for (i = PMLc[iEz].XI[REGION]; i <= PMLc[iEz].XE[REGION]; ++i) {
               i_m = i - b.Ez.XI;
               medio = sggMiEz[i_m][j_m][k_m];
               regLR[REGION].Psi_Ezy[i][j][k] = P_be_y[j] * regLR[REGION].Psi_Ezy[i][j][k] +
                                                (Hx[i_m][j_m][k_m] - Hx[i_m][j_m - 1][k_m]) * P_ce_y[j];
               Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] - G2[medio] * regLR[REGION].Psi_Ezy[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = right;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEx].ZI[REGION]; k <= PMLc[iEx].ZE[REGION]; ++k) {
         k_m = k - b.Ex.ZI;
         for (j = PMLc[iEx].YI[REGION]; j <= PMLc[iEx].YE[REGION]; ++j) {
            j_m = j - b.Ex.YI;
            for (i = PMLc[iEx].XI[REGION]; i <= PMLc[iEx].XE[REGION]; ++i) {
               i_m = i - b.Ex.XI;
               // --->
               medio = sggMiEx[i_m][j_m][k_m];
               regLR[REGION].Psi_Exy[i][j][k] = P_be_y[j] * regLR[REGION].Psi_Exy[i][j][k] +
                                                (Hz[i_m][j_m][k_m] - Hz[i_m][j_m - 1][k_m]) * P_ce_y[j];
               Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m] + G2[medio] * regLR[REGION].Psi_Exy[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEz].ZI[REGION]; k <= PMLc[iEz].ZE[REGION]; ++k) {
         k_m = k - b.Ez.ZI;
         for (j = PMLc[iEz].YI[REGION]; j <= PMLc[iEz].YE[REGION]; ++j) {
            j_m = j - b.Ez.YI;
            for (i = PMLc[iEz].XI[REGION]; i <= PMLc[iEz].XE[REGION]; ++i) {
               i_m = i - b.Ez.XI;
               medio = sggMiEz[i_m][j_m][k_m];
               regLR[REGION].Psi_Ezy[i][j][k] = P_be_y[j] * regLR[REGION].Psi_Ezy[i][j][k] +
                                                (Hx[i_m][j_m][k_m] - Hx[i_m][j_m - 1][k_m]) * P_ce_y[j];
               Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] - G2[medio] * regLR[REGION].Psi_Ezy[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif



      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = down;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEy].ZI[REGION]; k <= PMLc[iEy].ZE[REGION]; ++k) {
         k_m = k - b.Ey.ZI;
         for (j = PMLc[iEy].YI[REGION]; j <= PMLc[iEy].YE[REGION]; ++j) {
            j_m = j - b.Ey.YI;
            for (i = PMLc[iEy].XI[REGION]; i <= PMLc[iEy].XE[REGION]; ++i) {
               i_m = i - b.Ey.XI;
               medio = sggMiEy[i_m][j_m][k_m];
               regDU[REGION].Psi_Eyz[i][j][k] = P_be_z[k] * regDU[REGION].Psi_Eyz[i][j][k] +
                                                (Hx[i_m][j_m][k_m] - Hx[i_m][j_m][k_m - 1]) * P_ce_z[k];
               Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] + G2[medio] * regDU[REGION].Psi_Eyz[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEx].ZI[REGION]; k <= PMLc[iEx].ZE[REGION]; ++k) {
         k_m = k - b.Ex.ZI;
         for (j = PMLc[iEx].YI[REGION]; j <= PMLc[iEx].YE[REGION]; ++j) {
            j_m = j - b.Ex.YI;
            for (i = PMLc[iEx].XI[REGION]; i <= PMLc[iEx].XE[REGION]; ++i) {
               i_m = i - b.Ex.XI;
               medio = sggMiEx[i_m][j_m][k_m];
               regDU[REGION].Psi_Exz[i][j][k] = P_be_z[k] * regDU[REGION].Psi_Exz[i][j][k] +
                                                (Hy[i_m][j_m][k_m] - Hy[i_m][j_m][k_m - 1]) * P_ce_z[k];
               Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m] - G2[medio] * regDU[REGION].Psi_Exz[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif


      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = up;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEy].ZI[REGION]; k <= PMLc[iEy].ZE[REGION]; ++k) {
         k_m = k - b.Ey.ZI;
         for (j = PMLc[iEy].YI[REGION]; j <= PMLc[iEy].YE[REGION]; ++j) {
            j_m = j - b.Ey.YI;
            for (i = PMLc[iEy].XI[REGION]; i <= PMLc[iEy].XE[REGION]; ++i) {
               i_m = i - b.Ey.XI;
               medio = sggMiEy[i_m][j_m][k_m];
               regDU[REGION].Psi_Eyz[i][j][k] = P_be_z[k] * regDU[REGION].Psi_Eyz[i][j][k] +
                                                (Hx[i_m][j_m][k_m] - Hx[i_m][j_m][k_m - 1]) * P_ce_z[k];
               Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] + G2[medio] * regDU[REGION].Psi_Eyz[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEx].ZI[REGION]; k <= PMLc[iEx].ZE[REGION]; ++k) {
         k_m = k - b.Ex.ZI;
         for (j = PMLc[iEx].YI[REGION]; j <= PMLc[iEx].YE[REGION]; ++j) {
            j_m = j - b.Ex.YI;
            for (i = PMLc[iEx].XI[REGION]; i <= PMLc[iEx].XE[REGION]; ++i) {
               i_m = i - b.Ex.XI;
               medio = sggMiEx[i_m][j_m][k_m];
               regDU[REGION].Psi_Exz[i][j][k] = P_be_z[k] * regDU[REGION].Psi_Exz[i][j][k] +
                                                (Hy[i_m][j_m][k_m] - Hy[i_m][j_m][k_m - 1]) * P_ce_z[k];
               Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m] - G2[medio] * regDU[REGION].Psi_Exz[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif




      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = back;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEz].ZI[REGION]; k <= PMLc[iEz].ZE[REGION]; ++k) {
         k_m = k - b.Ez.ZI;
         for (j = PMLc[iEz].YI[REGION]; j <= PMLc[iEz].YE[REGION]; ++j) {
            j_m = j - b.Ez.YI;
            for (i = PMLc[iEz].XI[REGION]; i <= PMLc[iEz].XE[REGION]; ++i) {
               i_m = i - b.Ez.XI;
               medio = sggMiEz[i_m][j_m][k_m];
               regBF[REGION].Psi_Ezx[i][j][k] = P_be_x[i] * regBF[REGION].Psi_Ezx[i][j][k] +
                                                (Hy[i_m][j_m][k_m] - Hy[i_m - 1][j_m][k_m]) * P_ce_x[i];
               Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] + G2[medio] * regBF[REGION].Psi_Ezx[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEy].ZI[REGION]; k <= PMLc[iEy].ZE[REGION]; ++k) {
         k_m = k - b.Ey.ZI;
         for (j = PMLc[iEy].YI[REGION]; j <= PMLc[iEy].YE[REGION]; ++j) {
            j_m = j - b.Ey.YI;
            for (i = PMLc[iEy].XI[REGION]; i <= PMLc[iEy].XE[REGION]; ++i) {
               i_m = i - b.Ey.XI;
               medio = sggMiEy[i_m][j_m][k_m];
               regBF[REGION].Psi_Eyx[i][j][k] = P_be_x[i] * regBF[REGION].Psi_Eyx[i][j][k] +
                                                (Hz[i_m][j_m][k_m] - Hz[i_m - 1][j_m][k_m]) * P_ce_x[i];
               Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] - G2[medio] * regBF[REGION].Psi_Eyx[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif


      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = front;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEz].ZI[REGION]; k <= PMLc[iEz].ZE[REGION]; ++k) {
         k_m = k - b.Ez.ZI;
         for (j = PMLc[iEz].YI[REGION]; j <= PMLc[iEz].YE[REGION]; ++j) {
            j_m = j - b.Ez.YI;
            for (i = PMLc[iEz].XI[REGION]; i <= PMLc[iEz].XE[REGION]; ++i) {
               i_m = i - b.Ez.XI;
               medio = sggMiEz[i_m][j_m][k_m];
               regBF[REGION].Psi_Ezx[i][j][k] = P_be_x[i] * regBF[REGION].Psi_Ezx[i][j][k] +
                                                (Hy[i_m][j_m][k_m] - Hy[i_m - 1][j_m][k_m]) * P_ce_x[i];
               Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] + G2[medio] * regBF[REGION].Psi_Ezx[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEy].ZI[REGION]; k <= PMLc[iEy].ZE[REGION]; ++k) {
         k_m = k - b.Ey.ZI;
         for (j = PMLc[iEy].YI[REGION]; j <= PMLc[iEy].YE[REGION]; ++j) {
            j_m = j - b.Ey.YI;
            for (i = PMLc[iEy].XI[REGION]; i <= PMLc[iEy].XE[REGION]; ++i) {
               i_m = i - b.Ey.XI;
               medio = sggMiEy[i_m][j_m][k_m];
               regBF[REGION].Psi_Eyx[i][j][k] = P_be_x[i] * regBF[REGION].Psi_Eyx[i][j][k] +
                                                (Hz[i_m][j_m][k_m] - Hz[i_m - 1][j_m][k_m]) * P_ce_x[i];
               Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] - G2[medio] * regBF[REGION].Psi_Eyx[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif


      // ---------------------------> acaba AdvanceelectricCPML <---------------------------------------
      return;
   } // subroutine AdvanceelectricCPML
   // **************************************************************************************************
   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!! Advances the magnetic field in the PML
   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   void AdvanceMagneticCPML(int NumMedia, const bounds_t& b,
                            const std::vector<std::vector<std::vector<int>>>& sggMiHx,
                            const std::vector<std::vector<std::vector<int>>>& sggMiHy,
                            const std::vector<std::vector<std::vector<int>>>& sggMiHz,
                            const std::vector<double>& gm2,
                            const std::vector<std::vector<std::vector<double>>>& Ex,
                            const std::vector<std::vector<std::vector<double>>>& Ey,
                            const std::vector<std::vector<std::vector<double>>>& Ez,
                            std::vector<std::vector<std::vector<double>>>& Hx,
                            std::vector<std::vector<std::vector<double>>>& Hy,
                            std::vector<std::vector<std::vector<double>>>& Hz) {
      // ---------------------------> inputs <----------------------------------------------------------
      // integer, intent( IN) :: NumMedia
      // type( bounds_t), intent( IN) :: b
      // --->
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHx%NX-1, 0 :  b%sggMiHx%NY-1, 0 :  b%sggMiHx%NZ-1), intent( IN) :: sggMiHx
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHy%NX-1, 0 :  b%sggMiHy%NY-1, 0 :  b%sggMiHy%NZ-1), intent( IN) :: sggMiHy
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHz%NX-1, 0 :  b%sggMiHz%NY-1, 0 :  b%sggMiHz%NZ-1), intent( IN) :: sggMiHz
      // --->
      // real(kind = RKIND), dimension( 0 :  NumMedia), intent( IN) :: gm2
      // --->
      // real(kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( IN) :: Ex
      // real(kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( IN) :: Ey
      // real(kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( IN) :: Ez
      // ---------------------------> inputs/outputs <--------------------------------------------------
      // real(kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( INOUT) :: Hx
      // real(kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( INOUT) :: Hy
      // real(kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( INOUT) :: Hz
      // ---------------------------> variables locales <-----------------------------------------------
      int REGION, i, j, k, medio, i_m, j_m, k_m;
      // ---------------------------> empieza AdvanceMagneTicCPML <-------------------------------------
      // Hetic Fields PML Zone
      //
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = left;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iHx].ZI[REGION]; k <= PMLc[iHx].ZE[REGION]; ++k) {
         k_m = k - b.Hx.ZI;
         for (j = PMLc[iHx].YI[REGION]; j <= PMLc[iHx].YE[REGION]; ++j) {
            j_m = j - b.Hx.YI;
            for (i = PMLc[iHx].XI[REGION]; i <= PMLc[iHx].XE[REGION]; ++i) {
               i_m = i - b.Hx.XI;
               // --->
               regLR[REGION].Psi_Hxy[i][j][k] = P_bm_y[j] * regLR[REGION].Psi_Hxy[i][j][k] +
                                                (Ez[i_m][j_m + 1][k_m] - Ez[i_m][j_m][k_m]) * P_cm_y[j];
               medio = sggMiHx[i_m][j_m][k_m];
               Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] - GM2[medio] * regLR[REGION].Psi_Hxy[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iHz].ZI[REGION]; k <= PMLc[iHz].ZE[REGION]; ++k) {
         k_m = k - b.Hz.ZI;
         for (j = PMLc[iHz].YI[REGION]; j <= PMLc[iHz].YE[REGION]; ++j) {
            j_m = j - b.Hz.YI;
            for (i = PMLc[iHz].XI[REGION]; i <= PMLc[iHz].XE[REGION]; ++i) {
               i_m = i - b.Hz.XI;
               // --->
               regLR[REGION].Psi_Hzy[i][j][k] = P_bm_y[j] * regLR[REGION].Psi_Hzy[i][j][k] +
                                                (Ex[i_m][j_m + 1][k_m] - Ex[i_m][j_m][k_m]) * P_cm_y[j];
               medio = sggMiHz[i_m][j_m][k_m];
               Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] + GM2[medio] * regLR[REGION].Psi_Hzy[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      REGION = right;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iHx].ZI[REGION]; k <= PMLc[iHx].ZE[REGION]; ++k) {
         k_m = k - b.Hx.ZI;
         for (j = PMLc[iHx].YI[REGION]; j <= PMLc[iHx].YE[REGION]; ++j) {
            j_m = j - b.Hx.YI;
            for (i = PMLc[iHx].XI[REGION]; i <= PMLc[iHx].XE[REGION]; ++i) {
               i_m = i - b.Hx.XI;
               // --->
               regLR[REGION].Psi_Hxy[i][j][k] = P_bm_y[j] * regLR[REGION].Psi_Hxy[i][j][k] +
                                                (Ez[i_m][j_m + 1][k_m] - Ez[i_m][j_m][k_m]) * P_cm_y[j];
               medio = sggMiHx[i_m][j_m][k_m];
               Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] - GM2[medio] * regLR[REGION].Psi_Hxy[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iHz].ZI[REGION]; k <= PMLc[iHz].ZE[REGION]; ++k) {
         k_m = k - b.Hz.ZI;
         for (j = PMLc[iHz].YI[REGION]; j <= PMLc[iHz].YE[REGION]; ++j) {
            j_m = j - b.Hz.YI;
            for (i = PMLc[iHz].XI[REGION]; i <= PMLc[iHz].XE[REGION]; ++i) {
               i_m = i - b.Hz.XI;
               // --->
               regLR[REGION].Psi_Hzy[i][j][k] = P_bm_y[j] * regLR[REGION].Psi_Hzy[i][j][k] +
                                                (Ex[i_m][j_m + 1][k_m] - Ex[i_m][j_m][k_m]) * P_cm_y[j];
               medio = sggMiHz[i_m][j_m][k_m];
               Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] + GM2[medio] * regLR[REGION].Psi_Hzy[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = down;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iHy].ZI[REGION]; k <= PMLc[iHy].ZE[REGION]; ++k) {
         k_m = k - b.Hy.ZI;
         for (j = PMLc[iHy].YI[REGION]; j <= PMLc[iHy].YE[REGION]; ++j) {
            j_m = j - b.Hy.YI;
            for (i = PMLc[iHy].XI[REGION]; i <= PMLc[iHy].XE[REGION]; ++i) {
               i_m = i - b.Hy.XI;
               // --->
               regDU[REGION].Psi_Hyz[i][j][k] = P_bm_z[k] * regDU[REGION].Psi_Hyz[i][j][k] +
                                                (Ex[i_m][j_m][k_m + 1] - Ex[i_m][j_m][k_m]) * P_cm_z[k];
               medio = sggMiHy[i_m][j_m][k_m];
               Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] - GM2[medio] * regDU[REGION].Psi_Hyz[i][j][k];
            } // bucle i
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iHx].ZI[REGION]; k <= PMLc[iHx].ZE[REGION]; ++k) {
         k_m = k - b.Hx.ZI;
         for (j = PMLc[iHx].YI[REGION]; j <= PMLc[iHx].YE[REGION]; ++j) {
            j_m = j - b.Hx.YI;
            for (i = PMLc[iHx].XI[REGION]; i <= PMLc[iHx].XE[REGION]; ++i) {
               i_m = i - b.Hx.XI;
               // --->
               regDU[REGION].Psi_Hxz[i][j][k] = P_bm_z[k] * regDU[REGION].Psi_Hxz[i][j][k] +
                                                (Ey[i_m][j_m][k_m + 1] - Ey[i_m][j_m][k_m]) * P_cm_z[k];
               medio = sggMiHx[i_m][j_m][k_m];
               Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] + GM2[medio] * regDU[REGION].Psi_Hxz[i][j][k];
            } // bucle i
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = up;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iHy].ZI[REGION]; k <= PMLc[iHy].ZE[REGION]; ++k) {
         k_m = k - b.Hy.ZI;
         for (j = PMLc[iHy].YI[REGION]; j <= PMLc[iHy].YE[REGION]; ++j) {
            j_m = j - b.Hy.YI;
            for (i = PMLc[iHy].XI[REGION]; i <= PMLc[iHy].XE[REGION]; ++i) {
               i_m = i - b.Hy.XI;
               // --->
               regDU[REGION].Psi_Hyz[i][j][k] = P_bm_z[k] * regDU[REGION].Psi_Hyz[i][j][k] +
                                                (Ex[i_m][j_m][k_m + 1] - Ex[i_m][j_m][k_m]) * P_cm_z[k];
               medio = sggMiHy[i_m][j_m][k_m];
               Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] - GM2[medio] * regDU[REGION].Psi_Hyz[i][j][k];
            } // bucle i
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iHx].ZI[REGION]; k <= PMLc[iHx].ZE[REGION]; ++k) {
         k_m = k - b.Hx.ZI;
         for (j = PMLc[iHx].YI[REGION]; j <= PMLc[iHx].YE[REGION]; ++j) {
            j_m = j - b.Hx.YI;
            for (i = PMLc[iHx].XI[REGION]; i <= PMLc[iHx].XE[REGION]; ++i) {
               i_m = i - b.Hx.XI;
               // --->
               regDU[REGION].Psi_Hxz[i][j][k] = P_bm_z[k] * regDU[REGION].Psi_Hxz[i][j][k] +
                                                (Ey[i_m][j_m][k_m + 1] - Ey[i_m][j_m][k_m]) * P_cm_z[k];
               medio = sggMiHx[i_m][j_m][k_m];
               Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] + GM2[medio] * regDU[REGION].Psi_Hxz[i][j][k];
            } // bucle i
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = back;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iHz].ZI[REGION]; k <= PMLc[iHz].ZE[REGION]; ++k) {
         k_m = k - b.Hz.ZI;
         for (j = PMLc[iHz].YI[REGION]; j <= PMLc[iHz].YE[REGION]; ++j) {
            j_m = j - b.Hz.YI;
            for (i = PMLc[iHz].XI[REGION]; i <= PMLc[iHz].XE[REGION]; ++i) {
               i_m = i - b.Hz.XI;
               // --->
               regBF[REGION].Psi_Hzx[i][j][k] = P_bm_x[i] * regBF[REGION].Psi_Hzx[i][j][k] +
                                                (Ey[i_m + 1][j_m][k_m] - Ey[i_m][j_m][k_m]) * P_cm_x[i];
               medio = sggMiHz[i_m][j_m][k_m];
               Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] - GM2[medio] * regBF[REGION].Psi_Hzx[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

// This chunk continues the translation of the FDTD simulation code.
// It includes the end of AdvanceMagneTicCPML and the start of calc_cpmlconstants.
// Assumes previous chunks have defined types like SGGFDTDINFO_t, bounds_t, and global arrays.

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
for (k = PMLc[iHy].ZI(REGION); k <= PMLc[iHy].ZE(REGION); ++k) {
    k_m = k - b.Hy.ZI;
    for (j = PMLc[iHy].YI(REGION); j <= PMLc[iHy].YE(REGION); ++j) {
        j_m = j - b.Hy.YI;
        for (i = PMLc[iHy].XI(REGION); i <= PMLc[iHy].XE(REGION); ++i) {
            i_m = i - b.Hy.XI;
            //--->
            regBF[region].Psi_Hyx(i, j, k) = P_bm_x(i) * regBF[REGION].Psi_Hyx(i, j, k) +
                                              (Ez(i_m + 1, j_m, k_m) - Ez(i_m, j_m, k_m)) * P_cm_x(i);
            medio = sggMiHy(i_m, j_m, k_m);
            Hy(i_m, j_m, k_m) = Hy(i_m, j_m, k_m) + GM2(medio) * regBF[REGION].Psi_Hyx(i, j, k);
        }
    }
}
#ifdef CompileWithOpenMP
#pragma omp endparallel
#endif

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REGION = front;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
for (k = PMLc[iHz].ZI(REGION); k <= PMLc[iHz].ZE(REGION); ++k) {
    k_m = k - b.Hz.ZI;
    for (j = PMLc[iHz].YI(REGION); j <= PMLc[iHz].YE(REGION); ++j) {
        j_m = j - b.Hz.YI;
        for (i = PMLc[iHz].XI(REGION); i <= PMLc[iHz].XE(REGION); ++i) {
            i_m = i - b.Hz.XI;
            //--->
            regBF[REGION].Psi_Hzx(i, j, k) = P_bm_x(i) * regBF[REGION].Psi_Hzx(i, j, k) +
                                              (Ey(i_m + 1, j_m, k_m) - Ey(i_m, j_m, k_m)) * P_cm_x(i);
            medio = sggMiHz(i_m, j_m, k_m);
            Hz(i_m, j_m, k_m) = Hz(i_m, j_m, k_m) - GM2(medio) * regBF[REGION].Psi_Hzx(i, j, k);
        }
    }
}
#ifdef CompileWithOpenMP
#pragma omp endparallel
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
for (k = PMLc[iHy].ZI(REGION); k <= PMLc[iHy].ZE(REGION); ++k) {
    k_m = k - b.Hy.ZI;
    for (j = PMLc[iHy].YI(REGION); j <= PMLc[iHy].YE(REGION); ++j) {
        j_m = j - b.Hy.YI;
        for (i = PMLc[iHy].XI(REGION); i <= PMLc[iHy].XE(REGION); ++i) {
            i_m = i - b.Hy.XI;
            //--->
            regBF[region].Psi_Hyx(i, j, k) = P_bm_x(i) * regBF[REGION].Psi_Hyx(i, j, k) +
                                              (Ez(i_m + 1, j_m, k_m) - Ez(i_m, j_m, k_m)) * P_cm_x(i);
            medio = sggMiHy(i_m, j_m, k_m);
            Hy(i_m, j_m, k_m) = Hy(i_m, j_m, k_m) + GM2(medio) * regBF[REGION].Psi_Hyx(i, j, k);
        }
    }
}
#ifdef CompileWithOpenMP
#pragma omp endparallel
#endif

// ---------------------------> acaba AdvanceMagneTicCPML <---------------------------------------
return;
} // end subroutine AdvanceMagneTicCPML

// !!!
// !!!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!   !!! Advances the magnetic field in the PML
// !!!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!   subroutine FreeSpace_AdvanceMagneticCPML( NumMedia, b, gm2, Hx, Hy, Hz, Ex, Ey, Ez)
// !!!      !---------------------------> inputs <----------------------------------------------------------
// !!!      integer, intent( IN) :: NumMedia
// !!!      type( bounds_t), intent( IN) :: b
// !!!      !--->
// !!!      real(kind = RKIND), dimension( 0 :  NumMedia), intent( IN) :: gm2
// !!!      !--->
// !!!      real(kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( IN) :: Ex
// !!!      real(kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( IN) :: Ey
// !!!      real(kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( IN) :: Ez
// !!!      !---------------------------> inputs/outputs <--------------------------------------------------
// !!!      real(kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( INOUT) :: Hx
// !!!      real(kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( INOUT) :: Hy
// !!!      real(kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( INOUT) :: Hz
// !!!      !---------------------------> variables locales <-----------------------------------------------
// !!!      integer(kind=4) :: REGION, i, j, k, i_m, j_m, k_m
// !!!      real(kind = RKIND) :: GM2_1
// !!!      GM2_1=GM2(1)
// !!!      !---------------------------> empieza AdvanceMagneTicCPML <-------------------------------------
// !!!      !Hetic Fields PML Zone
// !!!      !
// !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!      REGION = left
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
// !!!#endif
// !!!      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
// !!!         k_m = k - b%Hx%ZI
// !!!         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
// !!!            j_m = j - b%Hx%YI
// !!!            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
// !!!               i_m = i - b%Hx%XI
// !!!               !--->
// !!!               regLR( REGION)%Psi_Hxy( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hxy( i, j, k) +  &
// !!!               (Ez( i_m, j_m+1, k_m) - Ez( i_m, j_m, k_m)) * P_cm_y( j)
// !!!               Hx( i_m, j_m, k_m)=Hx( i_m, j_m, k_m)-GM2_1*regLR(REGION)%Psi_Hxy(i,j,k)
// !!!            end do
// !!!         end do
// !!!      end do
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP END PARALLEL DO
// !!!#endif
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
// !!!#endif
// !!!      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
// !!!         k_m = k - b%Hz%ZI
// !!!         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
// !!!            j_m = j - b%Hz%YI
// !!!            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
// !!!               i_m = i - b%Hz%XI
// !!!               !--->
// !!!               regLR( REGION)%Psi_Hzy( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hzy( i, j, k) +  &
// !!!               (Ex( i_m, j_m+1, k_m) - Ex( i_m, j_m, k_m)) * P_cm_y( j)
// !!!               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) + GM2_1 * regLR( REGION)%Psi_Hzy( i, j, k)
// !!!            end do
// !!!         end do
// !!!      end do
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP END PARALLEL DO
// !!!#endif
// !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!
// !!!      REGION = right
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
// !!!#endif
// !!!      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
// !!!         k_m = k - b%Hx%ZI
// !!!         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
// !!!            j_m = j - b%Hx%YI
// !!!            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
// !!!               i_m = i - b%Hx%XI
// !!!               !--->
// !!!               regLR( REGION)%Psi_Hxy( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hxy( i, j, k) +  &
// !!!               (Ez( i_m, j_m+1, k_m) - Ez( i_m, j_m, k_m)) * P_cm_y( j)
// !!!               Hx( i_m, j_m, k_m)=Hx( i_m, j_m, k_m)-GM2_1*regLR(REGION)%Psi_Hxy(i,j,k)
// !!!            end do
// !!!         end do
// !!!      end do
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP END PARALLEL DO
// !!!#endif
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
// !!!#endif
// !!!      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
// !!!         k_m = k - b%Hz%ZI
// !!!         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
// !!!            j_m = j - b%Hz%YI
// !!!            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
// !!!               i_m = i - b%Hz%XI
// !!!               !--->
// !!!               regLR( REGION)%Psi_Hzy( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hzy( i, j, k) +  &
// !!!               (Ex( i_m, j_m+1, k_m) - Ex( i_m, j_m, k_m)) * P_cm_y( j)
// !!!               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) + GM2_1 * regLR( REGION)%Psi_Hzy( i, j, k)
// !!!            end do
// !!!         end do
// !!!      end do
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP END PARALLEL DO
// !!!#endif
// !!!
// !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!      REGION = down
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
// !!!#endif
// !!!      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
// !!!         k_m = k - b%Hy%ZI
// !!!         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
// !!!            j_m = j - b%Hy%YI
// !!!            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
// !!!               i_m = i - b%Hy%XI
// !!!               !--->
// !!!               regDU( REGION)%Psi_Hyz( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hyz( i, j, k) +  &
// !!!               (Ex( i_m, j_m, k_m+1) - Ex( i_m, j_m, k_m)) * P_cm_z( k)
// !!!               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) - GM2_1 * regDU( REGION)%Psi_Hyz( i, j, k)
// !!!            end do !bucle i
// !!!         end do
// !!!      end do
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP END PARALLEL DO
// !!!#endif
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
// !!!#endif
// !!!      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
// !!!         k_m = k - b%Hx%ZI
// !!!         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
// !!!            j_m = j - b%Hx%YI
// !!!            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
// !!!               i_m = i - b%Hx%XI
// !!!               !--->
// !!!               regDU( REGION)%Psi_Hxz( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hxz( i, j, k) +  &
// !!!               (Ey( i_m, j_m, k_m+1) - Ey( i_m, j_m, k_m)) * P_cm_z( k)
// !!!               Hx( i_m, j_m, k_m) = Hx( i_m, j_m, k_m) + GM2_1 * regDU( REGION)%Psi_Hxz( i, j, k)
// !!!            end do !bucle i
// !!!         end do
// !!!      end do
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP END PARALLEL DO
// !!!#endif
// !!!
// !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!      REGION = up
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
// !!!#endif
// !!!      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
// !!!         k_m = k - b%Hy%ZI
// !!!         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
// !!!            j_m = j - b%Hy%YI
// !!!            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
// !!!               i_m = i - b%Hy%XI
// !!!               !--->
// !!!               regDU( REGION)%Psi_Hyz( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hyz( i, j, k) +  &
// !!!               (Ex( i_m, j_m, k_m+1) - Ex( i_m, j_m, k_m)) * P_cm_z( k)
// !!!               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) - GM2_1 * regDU( REGION)%Psi_Hyz( i, j, k)
// !!!            end do !bucle i
// !!!         end do
// !!!      end do
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP END PARALLEL DO
// !!!#endif
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
// !!!#endif
// !!!      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
// !!!         k_m = k - b%Hx%ZI
// !!!         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
// !!!            j_m = j - b%Hx%YI
// !!!            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
// !!!               i_m = i - b%Hx%XI
// !!!               !--->
// !!!               regDU( REGION)%Psi_Hxz( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hxz( i, j, k) +  &
// !!!               (Ey( i_m, j_m, k_m+1) - Ey( i_m, j_m, k_m)) * P_cm_z( k)
// !!!               Hx( i_m, j_m, k_m) = Hx( i_m, j_m, k_m) + GM2_1 * regDU( REGION)%Psi_Hxz( i, j, k)
// !!!            end do !bucle i
// !!!         end do
// !!!      end do
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP END PARALLEL DO
// !!!#endif
// !!!
// !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!      REGION=back
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
// !!!#endif
// !!!      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
// !!!         k_m = k - b%Hz%ZI
// !!!         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
// !!!            j_m = j - b%Hz%YI
// !!!            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
// !!!               i_m = i - b%Hz%XI
// !!!               !--->
// !!!               regBF( REGION)%Psi_Hzx( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hzx( i, j, k) +  &
// !!!               (Ey( i_m+1, j_m, k_m) - Ey( i_m, j_m, k_m)) * P_cm_x( i)
// !!!               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) - GM2_1 * regBF( REGION)%Psi_Hzx( i, j, k)
// !!!            end do
// !!!         end do
// !!!      end do
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP END PARALLEL DO
// !!!#endif
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
// !!!#endif
// !!!      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
// !!!         k_m = k - b%Hy%ZI
// !!!         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
// !!!            j_m = j - b%Hy%YI
// !!!            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
// !!!               i_m = i - b%Hy%XI
// !!!               !--->
// !!!               regBF( region)%Psi_Hyx( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hyx( i, j, k) +  &
// !!!               (Ez( i_m+1, j_m, k_m) - Ez( i_m, j_m, k_m)) * P_cm_x( i)
// !!!               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) + GM2_1 * regBF( REGION)%Psi_Hyx( i, j, k)
// !!!            end do
// !!!         end do
// !!!      end do
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP END PARALLEL DO
// !!!#endif
// !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!      REGION=front
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
// !!!#endif
// !!!      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
// !!!         k_m = k - b%Hz%ZI
// !!!         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
// !!!            j_m = j - b%Hz%YI
// !!!            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
// !!!               i_m = i - b%Hz%XI
// !!!               !--->
// !!!               regBF( REGION)%Psi_Hzx( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hzx( i, j, k) +  &
// !!!               (Ey( i_m+1, j_m, k_m) - Ey( i_m, j_m, k_m)) * P_cm_x( i)
// !!!               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) - GM2_1 * regBF( REGION)%Psi_Hzx( i, j, k)
// !!!            end do
// !!!         end do
// !!!      end do
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP END PARALLEL DO
// !!!#endif
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
// !!!#endif
// !!!      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
// !!!         k_m = k - b%Hy%ZI
// !!!         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
// !!!            j_m = j - b%Hy%YI
// !!!            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
// !!!               i_m = i - b%Hy%XI
// !!!               !--->
// !!!               regBF( region)%Psi_Hyx( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hyx( i, j, k) +  &
// !!!               (Ez( i_m+1, j_m, k_m) - Ez( i_m, j_m, k_m)) * P_cm_x( i)
// !!!               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) + GM2_1 * regBF( REGION)%Psi_Hyx( i, j, k)
// !!!            end do
// !!!         end do
// !!!      end do
// !!!#ifdef CompileWithOpenMP
// !!!!$OMP END PARALLEL DO
// !!!#endif
// !!!
// !!!
// !!!      !---------------------------> acaba AdvanceMagneTicCPML <---------------------------------------
// !!!      return
// !!!   endsubroutine FreeSpace_AdvanceMagneTicCPML

void calc_cpmlconstants(const SGGFDTDINFO_t& sgg, int Idxe[], int Idye[], int Idze[], int Idxh[], int Idyh[], int Idzh[], double eps00, double mu00) {
    double eps0, mu0, zvac;
    double del, sigmae, kpare, apare, sigmam, kparm, aparm;
    // Note: BUFSIZE is assumed to be defined elsewhere or handled via std::string if needed, 
    // but here it's just a local variable not strictly used for logic in this chunk.
    // char buff[BUFSIZE]; 

    eps0 = eps00;
    mu0 = mu00;
    // chapuz para convertir la variables de paso en globales
    zvac = sqrt(mu0 / eps0);

    del = 1.0; // una simple inicializacion para que gfortran no se queje
    
    // Find the maximum conductivity for each direcion o=1,2,3 and for the starting and ending layer p=1,2
    Sig_max = 0.0;
    aPar_max = 0.0;
    kPar_max = 0.0;
    
    for (int o = 1; o <= 3; ++o) {
        for (int p = 1; p <= 2; ++p) {
            if ((o == 1) && (p == 1)) del = dxe(sgg.ALLOC[iHx].XI);
            if ((o == 1) && (p == 2)) del = dxe(sgg.ALLOC[iHx].XE);
            if ((o == 2) && (p == 1)) del = dye(sgg.ALLOC[iHy].YI);
            if ((o == 2) && (p == 2)) del = dye(sgg.ALLOC[iHy].YE);
            if ((o == 3) && (p == 1)) del = dze(sgg.ALLOC[iHz].ZI);
            if ((o == 3) && (p == 2)) del = dze(sgg.ALLOC[iHz].ZE);
            
            if (sgg.PML.NumLayers[o][p] != 0) {
                if ((sgg.PML.orden[o][p] == 10) || (sgg.PML.orden[o][p] == 5)) {
                    // gedney sigma optimo no tiene en cuenta el refle (taflove 3 ed, pag 294)
                    // para 5 celdas es exp(-8) y para 16 celdas es exp(-16)
                    Sig_max[o][p] = 0.8 * (sgg.PML.orden[o][p] + 1) / (sqrt(Mu0 / eps0) * del); // cambio tonto no afecta a nada 260919
                } else {
                    Sig_max[o][p] = -((log(sgg.PML.CoeffReflPML[o][p]) * (sgg.PML.orden[o][p] + 1)) /
                                       (2 * sqrt(Mu0 / eps0) * sgg.PML.NumLayers[o][p] * del));
                    // ojo LO SIGUIENTE  estaba maaaaallllllllll porque NO ES MATERIAL INDEPENDENT
                    //                          (2*sqrt(Mu0/eps0)*sqrt(sgg%Med(jmed)%Epr*sgg%Med(jmed)%Mur)*sgg%PML%NumLayers(o,p)*del))
                    // los multilayer petan- !!! !Viene de Gedney, pero esta maaaalllll!! !corregido 20marzo 2011
                }
            } else {
                Sig_max[o][p] = 1.0e29;
            }
        }
    }

    // readjust relatively the alphamaxpar to the maximum conductivity
    for (int o = 1; o <= 3; ++o) {
        for (int p = 1; p <= 2; ++p) {
            if ((o == 1) && (p == 1)) del = dxe(sgg.ALLOC[iEx].XI);
            if ((o == 1) && (p == 2)) del = dxe(sgg.ALLOC[iEx].XE);
            if ((o == 2) && (p == 1)) del = dye(sgg.ALLOC[iEy].YI);
            if ((o == 2) && (p == 2)) del = dye(sgg.ALLOC[iEy].YE);
            if ((o == 3) && (p == 1)) del = dze(sgg.ALLOC[iEz].ZI);
            if ((o == 3) && (p == 2)) del = dze(sgg.ALLOC[iEz].ZE);
            
            if (sgg.PML.NumLayers[o][p] != 0) {
                aPar_max[o][p] = alphamaxpar * Sig_max[o][p];
                kPar_max[o][p] = kappamaxpar;
            } else {
                aPar_max[o][p] = 0.0;
                kPar_max[o][p] = 1.0;
            }
        }
    }

    P_ce_x = 0.0;
    P_ce_y = 0.0;
    P_ce_z = 0.0;
    P_cm_x = 0.0;
    P_cm_y = 0.0;
    P_cm_z = 0.0;
    P_be_x = 1.0;
    P_be_y = 1.0;
    P_be_z = 1.0;
    P_bm_x = 1.0;
    P_bm_y = 1.0;
    P_bm_z = 1.0;
    
    // Calculate the coefficients (in CPML they are the same for every possible medium OJOOOOOOOOOOOOOOOOO)

    // Default !2011 already done in main
    // Idxh        ( : )=    1.0 / dxh( : )
    // Idyh        ( : )=    1.0 / dyh( : )
    // Idzh        ( : )=    1.0 / dzh( : )
    // Idxe        ( : )=    1.0 / dxe( : )
    // Idye        ( : )=    1.0 / dye( : )
    // Idze        ( : )=    1.0 / dze( : )
    
    // Calculate
    for (i = sgg.ALLOC[iEx].XI; i <= sgg.ALLOC[iEx].XE; ++i) {
        if (i <= SINPML_Fullsize[iHx].XI - 1) { // Back
            if ((sgg.PML.orden[1][1] == 0)) {
                Sigmae = Sig_max[1][1];
                kPare = 1.0 + (kPar_max[1][1] - 1);
            } else {
                Sigmae = Sig_max[1][1] * ce_x[i] ** sgg.PML.orden[1][1];
                kPare = 1.0 + (kPar_max[1][1] - 1) * ce_x[i] ** sgg.PML.orden[1][1];
            }
            aPare = aPar_max[1][1] * Ice_x[i] ** alphaOrden; // **sgg%PML%orden(1,1) !!**1.0 !perfil lineal propuesto por Gedney originalmente !gedney lo escala linealmente
            P_be_x[i] = exp(-(sigmae / kPare + aPare) * sgg.dt / Eps0);
            P_ce_x[i] = (sigmae * (P_be_x[i] - 1.0) / (sigmae + kPare * aPare) / kpare) / dxh[i];
            Idxh[i] = 1.0 / (kPare * dxh[i]);
        } else if (i >= SINPML_Fullsize[iHx].XE + 1) { // Front
            if ((sgg.PML.orden[1][2] == 0)) {
                Sigmae = Sig_max[1][2];
                kPare = 1.0 + (kPar_max[1][2] - 1);
            } else {
                Sigmae = Sig_max[1][2] * ce_x[i] ** sgg.PML.orden[1][2];
                kPare = 1.0 + (kPar_max[1][2] - 1) * ce_x[i] ** sgg.PML.orden[1][2];
            }
            aPare = aPar_max[1][2] * Ice_x[i] ** alphaOrden; // **sgg%PML%orden(1,2) !!**1.0 !perfil lineal propuesto por Gedney originalmente
            P_be_x[i] = exp(-(sigmae / kPare + aPare) * sgg.dt / Eps0);
            P_ce_x[i] = (sigmae * (P_be_x[i] - 1.0) / (sigmae + kPare * aPare) / kpare) / dxh[i];
            Idxh[i] = 1.0 / (kPare * dxh[i]);
        }
    }
    
    for (j = sgg.ALLOC[iEy].YI; j <= sgg.ALLOC[iEy].YE; ++j) {
        if (j <= SINPML_Fullsize[iHy].YI - 1) { // Left
            if ((sgg.PML.orden[2][1] == 0)) {
                Sigmae = Sig_max[2][1];
                kPare = 1.0 + (kPar_max[2][1] - 1);
            } else {
                Sigmae = Sig_max[2][1] * ce_y[j] ** sgg.PML.orden[2][1];
                kPare = 1.0 + (kPar_max[2][1] - 1) * ce_y[j] ** sgg.PML.orden[2][1];
            }
            aPare = aPar_max[2][1] * Ice_y[j] ** alphaOrden; // **sgg%PML%orden(2,1) !!**1.0 !perfil lineal propuesto por Gedney originalmente
            P_be_y[j] = exp(-(sigmae / kPare + aPare) * sgg.dt / Eps0);
            P_ce_y[j] = (sigmae * (P_be_y[j] - 1.0) / (sigmae + kPare * aPare) / kpare) / dyh[j];
            IdYh[j] = 1.0 / (kPare * dyh[j]);
        } else if (j >= SINPML_Fullsize[iHy].YE + 1) { // Right
            if ((sgg.PML.orden[2][2] == 0)) {
                Sigmae = Sig_max[2][2];
                kPare = 1.0 + (kPar_max[2][2] - 1);
            } else {
                Sigmae = Sig_max[2][2] * ce_y[j] ** sgg.PML.orden[2][2];
                kPare = 1.0 + (kPar_max[2][2] - 1) * ce_y[j] ** sgg.PML.orden[2][2];
            }
            aPare = aPar_max[2][2] * Ice_y[j] ** alphaOrden; // **sgg%PML%orden(2,2) !!**1.0 !perfil lineal propuesto por Gedney originalmente
            P_be_y[j] = exp(-(sigmae / kPare + aPare) * sgg.dt / Eps0);
            P_ce_y[j] = (sigmae * (P_be_y[j] - 1.0) / (sigmae + kPare * aPare) / kpare) / dyh[j];
            IdYh[j] = 1.0 / (kPare * dyh[j]);
        }
    }
    
    for (k = sgg.ALLOC[iEz].ZI; k <= sgg.ALLOC[iEz].ZE; ++k) {
        if (k <= SINPML_Fullsize[iHz].ZI - 1) { // Down
            if ((sgg.PML.orden[3][1] == 0)) {
                sigmae = Sig_max[3][1];
                kPare = 1.0 + (kPar_max[3][1] - 1);
            } else {
                sigmae = Sig_max[3][1] * ce_z[k] ** sgg.PML.orden[3][1];
                kPare = 1.0 + (kPar_max[3][1] - 1) * ce_z[k] ** sgg.PML.orden[3][1];
            }
            aPare = aPar_max[3][1] * Ice_z[k] ** alphaOrden; // **sgg%PML%orden(3,1) !!**1.0 !perfil lineal propuesto por Gedney originalmente
            P_be_z[k] = exp(-(sigmae / kPare + aPare) * sgg.dt / Eps0);
            P_ce_z[k] = (sigmae * (P_be_z[k] - 1.0) / (sigmae + kPare * aPare) / kpare) / dzh[k];
            Idzh[k] = 1.0 / (kPare * dzh[k]);
        } else if (k >= SINPML_Fullsize[iHz].ZE + 1) { // Up
            if ((sgg.PML.orden[3][2] == 0)) {
                sigmae = Sig_max[3][2];
                kPare = 1.0 + (kPar_max[3][2] - 1);
            } else {
                sigmae = Sig_max[3][2] * ce_z[k] ** sgg.PML.orden[3][2];
                kPare = 1.0 + (kPar_max[3][2] - 1) * ce_z[k] ** sgg.PML.orden[3][2];
            }
            aPare = aPar_max[3][2] * Ice_z[k] ** alphaOrden; // **sgg%PML%orden(3,2) !!**1.0 !perfil lineal propuesto por Gedney originalmente
            P_be_z[k] = exp(-(sigmae / kPare + aPare) * sgg.dt / Eps0);
            P_ce_z[k] = (sigmae * (P_be_z[k] - 1.0) / (sigmae + kPare * aPare) / kpare) / dzh[k];
            Idzh[k] = 1.0 / (kPare * dzh[k]);
        }
    }
    
    // magnetic
    for (i = sgg.ALLOC[iHx].XI; i <= sgg.ALLOC[iHx].XE; ++i) {
        if (i <= SINPML_Fullsize[iHx].XI - 1) { // back
            if ((sgg.PML.orden[1][1] == 0)) {
                Sigmam = Sig_max[1][1];
                kParm = 1.0 + (kPar_max[1][1] - 1);
            } else {
                Sigmam = Sig_max[1][1] * cm_x[i] ** sgg.PML.orden[1][1];
                kParm = 1.0 + (kPar_max[1][1] - 1) * cm_x[i] ** sgg.PML.orden[1][1];
            }
            aParm = aPar_max[1][1] * Icm_x[i] ** alphaOrden; // **sgg%PML%orden(1,1) !!**1.0 !perfil lineal propuesto por Gedney originalmente
        }
    }
}

P_bm_x[i] = std::exp(-(sigmam / kParm + aParm) * sgg.dt / Eps0);
            P_cm_x[i] = (sigmam * (P_bm_x[i] - 1.0_RKIND) / (sigmam + kParm * aParm) / kparm) / dxe[i];
            Idxe[i] = 1.0_RKIND / (kParm * dxe[i]);
            //
            // Note: Fortran write statement converted to commented out C++ equivalent or removed for brevity in logic flow
            // Original: write(buff,'(a,i4,a,5e9.2e2)') 'back(',i,'+d/2), A,S,FcS,FcA,refleLF=',aparm,sigmam,sigmam/(2.0_RKIND * pi*eps0),aparm/(2.0_RKIND * pi*eps0), &
            //          (sqrt(kparm+sigmam/(aParm+1d-15))-1.0_RKIND)/(sqrt(kparm+sigmam/(aParm+1d-15))+1.0_RKIND)
            // if ((sgg%Border%IsBackPML).and.(i>sgg%ALLOC(iHx)%XI)) call print11 (control%layoutnumber, buff)
         } else if (i >= SINPML_Fullsize(iHx).XE) { // front
            if ((sgg.PML.orden(1, 2) == 0)) {
               Sigmam = Sig_max(1, 2);
               kParm = 1.0_RKIND + (kPar_max(1, 2) - 1);
            } else {
               Sigmam = Sig_max(1, 2) * cm_x[i] ** sgg.PML.orden(1, 2);
               kParm = 1.0_RKIND + (kPar_max(1, 2) - 1) * cm_x[i] ** sgg.PML.orden(1, 2);
            }
            aParm = aPar_max(1, 2) * Icm_x[i] ** alphaOrden; // **sgg%PML%orden(1,2) !!**1.0_RKIND !perfil lineal propuesto por Gedney originalmente
            P_bm_x[i] = std::exp(-(sigmam / kParm + aParm) * sgg.dt / Eps0);
            P_cm_x[i] = (sigmam * (P_bm_x[i] - 1.0_RKIND) / (sigmam + kParm * aParm) / kparm) / dxe[i];
            Idxe[i] = 1.0_RKIND / (kParm * dxe[i]);
            //
            // !write(buff,'(a,i4,a,5e9.2e2)') 'front(',i,'+d/2), A,S,FcS,FcA,refleLF=',aparm,sigmam,sigmam/(2.0_RKIND * pi*eps0),aparm/(2.0_RKIND * pi*eps0), &
            // !(sqrt(kparm+sigmam/(aParm+1d-15))-1.0_RKIND)/(sqrt(kparm+sigmam/(aParm+1d-15))+1.0_RKIND)
            // if ((sgg%Border%IsFrontPML).and.(i<sgg%ALLOC(iHx)%XE-1))  call print11 (control%layoutnumber, buff)
         }
      }
      for (j = sgg.ALLOC(iHy).YI; j <= sgg.ALLOC(iHy).YE; ++j) {
         if (j <= SINPML_Fullsize(iHy).YI - 1) { // Left
            if ((sgg.PML.orden(2, 1) == 0)) {
               Sigmam = Sig_max(2, 1);
               kParm = 1.0_RKIND + (kPar_max(2, 1) - 1);
            } else {
               Sigmam = Sig_max(2, 1) * cm_y[j] ** sgg.PML.orden(2, 1);
               kParm = 1.0_RKIND + (kPar_max(2, 1) - 1) * cm_y[j] ** sgg.PML.orden(2, 1);
            }
            aParm = aPar_max(2, 1) * Icm_y[j] ** alphaOrden; // **sgg%PML%orden(2,1) !!**1.0_RKIND !perfil lineal propuesto por Gedney originalmente
            P_bm_y[j] = std::exp(-(sigmam / kParm + aParm) * sgg.dt / Eps0);
            P_cm_y[j] = (sigmam * (P_bm_y[j] - 1.0_RKIND) / (sigmam + kParm * aParm) / kparm) / dye[j];
            Idye[j] = 1.0_RKIND / (kParm * dye[j]);
            //
            // !write(buff,'(a,i4,a,5e9.2e2)') 'left(',j,'+d/2), A,S,FcS,FcA,refleLF=',aparm,sigmam,sigmam/(2.0_RKIND * pi*eps0),aparm/(2.0_RKIND * pi*eps0), &
            // !(sqrt(kparm+sigmam/(aParm+1d-15))-1.0_RKIND)/(sqrt(kparm+sigmam/(aParm+1d-15))+1.0_RKIND)
            // if ((sgg%Border%IsLeftPML).and.(j>sgg%ALLOC(iHy)%YI)) call print11 (control%layoutnumber, buff)
         } else if (j >= SINPML_Fullsize(iHy).YE) { // Right
            if ((sgg.PML.orden(2, 2) == 0)) {
               Sigmam = Sig_max(2, 2);
               kParm = 1.0_RKIND + (kPar_max(2, 2) - 1);
            } else {
               Sigmam = Sig_max(2, 2) * cm_y[j] ** sgg.PML.orden(2, 2);
               kParm = 1.0_RKIND + (kPar_max(2, 2) - 1) * cm_y[j] ** sgg.PML.orden(2, 2);
            }
            aParm = aPar_max(2, 2) * Icm_y[j] ** alphaOrden; // **sgg%PML%orden(2,2) !!**1.0_RKIND !perfil lineal propuesto por Gedney originalmente
            P_bm_y[j] = std::exp(-(sigmam / kParm + aParm) * sgg.dt / Eps0);
            P_cm_y[j] = (sigmam * (P_bm_y[j] - 1.0_RKIND) / (sigmam + kParm * aParm) / kparm) / dye[j];
            Idye[j] = 1.0_RKIND / (kParm * dye[j]);
            //
            // !write(buff,'(a,i4,a,5e9.2e2)') 'right(',j,'+d/2), A,S,FcS,FcA,refleLF=',aparm,sigmam,sigmam/(2.0_RKIND * pi*eps0),aparm/(2.0_RKIND * pi*eps0), &
            // !(sqrt(kparm+sigmam/(aParm+1d-15))-1.0_RKIND)/(sqrt(kparm+sigmam/(aParm+1d-15))+1.0_RKIND)
            // if ((sgg%Border%IsRightPML).and.(j<sgg%ALLOC(iHy)%YE-1)) call print11 (control%layoutnumber, buff)
         }
      }
      for (k = sgg.ALLOC(iHz).ZI; k <= sgg.ALLOC(iHz).ZE; ++k) {
         if (k <= SINPML_Fullsize(iHz).ZI - 1) { // Down
            if ((sgg.PML.orden(3, 1) == 0)) {
               Sigmam = Sig_max(3, 1);
               kParm = 1.0_RKIND + (kPar_max(3, 1) - 1);
            } else {
               Sigmam = Sig_max(3, 1) * cm_z[k] ** sgg.PML.orden(3, 1);
               kParm = 1.0_RKIND + (kPar_max(3, 1) - 1) * cm_z[k] ** sgg.PML.orden(3, 1);
            }
            aParm = aPar_max(3, 1) * Icm_z[k] ** alphaOrden; // **sgg%PML%orden(3,1) !!**1.0_RKIND !perfil lineal propuesto por Gedney originalmente
            P_bm_z[k] = std::exp(-(sigmam / kParm + aParm) * sgg.dt / Eps0);
            P_cm_z[k] = (sigmam * (P_bm_z[k] - 1.0_RKIND) / (sigmam + kParm * aParm) / kparm) / dze[k];
            Idze[k] = 1.0_RKIND / (kParm * dze[k]);
            //
            // !write(buff,'(a,i4,a,5e9.2e2)') 'down(',k,'+d/2), A,S,FcS,FcA,refleLF=',aparm,sigmam,sigmam/(2.0_RKIND * pi*eps0),aparm/(2.0_RKIND * pi*eps0), &
            // !(sqrt(kparm+sigmam/(aParm+1d-15))-1.0_RKIND)/(sqrt(kparm+sigmam/(aParm+1d-15))+1.0_RKIND)
            // if ((sgg%Border%IsDownPML).and.(k>sgg%ALLOC(iHz)%ZI)) call print11 (control%layoutnumber, buff)
         } else if (k >= SINPML_Fullsize(iHz).ZE) { // Up
            if ((sgg.PML.orden(3, 2) == 0)) {
               Sigmam = Sig_max(3, 2);
               kParm = 1.0_RKIND + (kPar_max(3, 2) - 1);
            } else {
               Sigmam = Sig_max(3, 2) * cm_z[k] ** sgg.PML.orden(3, 2);
               kParm = 1.0_RKIND + (kPar_max(3, 2) - 1) * cm_z[k] ** sgg.PML.orden(3, 2);
            }
            aParm = aPar_max(3, 2) * Icm_z[k] ** alphaOrden; // **sgg%PML%orden(3,2) !!**1.0_RKIND !perfil lineal propuesto por Gedney originalmente
            P_bm_z[k] = std::exp(-(sigmam / kParm + aParm) * sgg.dt / Eps0);
            P_cm_z[k] = (sigmam * (P_bm_z[k] - 1.0_RKIND) / (sigmam + kParm * aParm) / kparm) / dze[k];
            Idze[k] = 1.0_RKIND / (kParm * dze[k]);
            //
            // !write(buff,'(a,i4,a,5e9.2e2)') 'up   (',k,'+d/2), A,S,FcS,FcA,refleLF=',aparm,sigmam,sigmam/(2.0_RKIND * pi*eps0),aparm/(2.0_RKIND * pi*eps0), &
            // !(sqrt(kparm+sigmam/(aParm+1d-15))-1.0_RKIND)/(sqrt(kparm+sigmam/(aParm+1d-15))+1.0_RKIND)
            // if ((sgg%Border%IsUpPML).and.(k<sgg%ALLOC(iHz)%ZE-1)) call print11 (control%layoutnumber, buff)
         }
      }

   } // end subroutine calc_cpmlconstants

   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !**************************************************************************************************
   void AdvanceelectricCPML_freespace(int NumMedia, const bounds_t& b, 
                                      const std::vector<std::vector<std::vector<int>>>& sggMiEx, 
                                      const std::vector<std::vector<std::vector<int>>>& sggMiEy, 
                                      const std::vector<std::vector<std::vector<int>>>& sggMiEz, 
                                      const std::vector<double>& g2, 
                                      const std::vector<std::vector<std::vector<double>>>& Hx, 
                                      const std::vector<std::vector<std::vector<double>>>& Hy, 
                                      const std::vector<std::vector<std::vector<double>>>& Hz, 
                                      std::vector<std::vector<std::vector<double>>>& Ex, 
                                      std::vector<std::vector<std::vector<double>>>& Ey, 
                                      std::vector<std::vector<std::vector<double>>>& Ez) {
      // ---------------------------> inputs <----------------------------------------------------------
      // integer, intent( IN) :: NumMedia
      // type( bounds_t), intent( IN) :: b
      // !--->
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEx%NX-1, 0 :  b%sggMiEx%NY-1, 0 :  b%sggMiEx%NZ-1), intent( IN) :: sggMiEx
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEy%NX-1, 0 :  b%sggMiEy%NY-1, 0 :  b%sggMiEy%NZ-1), intent( IN) :: sggMiEy
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEz%NX-1, 0 :  b%sggMiEz%NY-1, 0 :  b%sggMiEz%NZ-1), intent( IN) :: sggMiEz
      // !--->
      // real(kind = RKIND), dimension( 0 :  NumMedia), intent( IN) :: g2
      // !--->
      // real(kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( IN) :: Hx
      // real(kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( IN) :: Hy
      // real(kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( IN) :: Hz
      // ---------------------------> inputs/outputs <--------------------------------------------------
      // real(kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( INOUT) :: Ex
      // real(kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( INOUT) :: Ey
      // real(kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( INOUT) :: Ez
      // ---------------------------> variables locales <-----------------------------------------------
      int REGION, i, j, k, medio, i_m, j_m, k_m;
      // ---------------------------> empieza AdvanceelectricCPML <-------------------------------------

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = left;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc(iEx).ZI(REGION); k <= PMLc(iEx).ZE(REGION); ++k) {
         k_m = k - b.Ex.ZI;
         for (j = PMLc(iEx).YI(REGION); j <= PMLc(iEx).YE(REGION); ++j) {
            j_m = j - b.Ex.YI;
            for (i = PMLc(iEx).XI(REGION); i <= PMLc(iEx).XE(REGION); ++i) {
               i_m = i - b.Ex.XI;
               // --->
               medio = 1;
               regLR(REGION).Psi_Exyvac(i, j, k) = P_be_y[j] * regLR(REGION).Psi_Exyvac(i, j, k) +
                                                    (Hz[i_m][j_m][k_m] - Hz[i_m][j_m - 1][k_m]) * P_ce_y[j];
               Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m] + G2[medio] * regLR(REGION).Psi_Exyvac(i, j, k);
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc(iEz).ZI(REGION); k <= PMLc(iEz).ZE(REGION); ++k) {
         k_m = k - b.Ez.ZI;
         for (j = PMLc(iEz).YI(REGION); j <= PMLc(iEz).YE(REGION); ++j) {
            j_m = j - b.Ez.YI;
            for (i = PMLc(iEz).XI(REGION); i <= PMLc(iEz).XE(REGION); ++i) {
               i_m = i - b.Ez.XI;
               medio = 1;
               regLR(REGION).Psi_Ezyvac(i, j, k) = P_be_y[j] * regLR(REGION).Psi_Ezyvac(i, j, k) +
                                                    (Hx[i_m][j_m][k_m] - Hx[i_m][j_m - 1][k_m]) * P_ce_y[j];
               Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] - G2[medio] * regLR(REGION).Psi_Ezyvac(i, j, k);
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = right;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc(iEx).ZI(REGION); k <= PMLc(iEx).ZE(REGION); ++k) {
         k_m = k - b.Ex.ZI;
         for (j = PMLc(iEx).YI(REGION); j <= PMLc(iEx).YE(REGION); ++j) {
            j_m = j - b.Ex.YI;
            for (i = PMLc(iEx).XI(REGION); i <= PMLc(iEx).XE(REGION); ++i) {
               i_m = i - b.Ex.XI;
               // --->
               medio = 1;
               regLR(REGION).Psi_Exyvac(i, j, k) = P_be_y[j] * regLR(REGION).Psi_Exyvac(i, j, k) +
                                                    (Hz[i_m][j_m][k_m] - Hz[i_m][j_m - 1][k_m]) * P_ce_y[j];
               Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m] + G2[medio] * regLR(REGION).Psi_Exyvac(i, j, k);
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc(iEz).ZI(REGION); k <= PMLc(iEz).ZE(REGION); ++k) {
         k_m = k - b.Ez.ZI;
         for (j = PMLc(iEz).YI(REGION); j <= PMLc(iEz).YE(REGION); ++j) {
            j_m = j - b.Ez.YI;
            for (i = PMLc(iEz).XI(REGION); i <= PMLc(iEz).XE(REGION); ++i) {
               i_m = i - b.Ez.XI;
               medio = 1;
               regLR(REGION).Psi_Ezyvac(i, j, k) = P_be_y[j] * regLR(REGION).Psi_Ezyvac(i, j, k) +
                                                    (Hx[i_m][j_m][k_m] - Hx[i_m][j_m - 1][k_m]) * P_ce_y[j];
               Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] - G2[medio] * regLR(REGION).Psi_Ezyvac(i, j, k);
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif



      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = down;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc(iEy).ZI(REGION); k <= PMLc(iEy).ZE(REGION); ++k) {
         k_m = k - b.Ey.ZI;
         for (j = PMLc(iEy).YI(REGION); j <= PMLc(iEy).YE(REGION); ++j) {
            j_m = j - b.Ey.YI;
            for (i = PMLc(iEy).XI(REGION); i <= PMLc(iEy).XE(REGION); ++i) {
               i_m = i - b.Ey.XI;
               medio = 1;
               regDU(REGION).Psi_Eyzvac(i, j, k) = P_be_z[k] * regDU(REGION).Psi_Eyzvac(i, j, k) +
                                                    (Hx[i_m][j_m][k_m] - Hx[i_m][j_m][k_m - 1]) * P_ce_z[k];
               Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] + G2[medio] * regDU(REGION).Psi_Eyzvac(i, j, k);
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc(iEx).ZI(REGION); k <= PMLc(iEx).ZE(REGION); ++k) {
         k_m = k - b.Ex.ZI;
         for (j = PMLc(iEx).YI(REGION); j <= PMLc(iEx).YE(REGION); ++j) {
            j_m = j - b.Ex.YI;
            for (i = PMLc(iEx).XI(REGION); i <= PMLc(iEx).XE(REGION); ++i) {
               i_m = i - b.Ex.XI;
               medio = 1;
               regDU(REGION).Psi_Exzvac(i, j, k) = P_be_z[k] * regDU(REGION).Psi_Exzvac(i, j, k) +
                                                    (Hy[i_m][j_m][k_m] - Hy[i_m][j_m][k_m - 1]) * P_ce_z[k];
               Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m] - G2[medio] * regDU(REGION).Psi_Exzvac(i, j, k);
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif


      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = up;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc(iEy).ZI(REGION); k <= PMLc(iEy).ZE(REGION); ++k) {
         k_m = k - b.Ey.ZI;
         for (j = PMLc(iEy).YI(REGION); j <= PMLc(iEy).YE(REGION); ++j) {
            j_m = j - b.Ey.YI;
            for (i = PMLc(iEy).XI(REGION); i <= PMLc(iEy).XE(REGION); ++i) {
               i_m = i - b.Ey.XI;
               medio = 1;
               regDU(REGION).Psi_Eyzvac(i, j, k) = P_be_z[k] * regDU(REGION).Psi_Eyzvac(i, j, k) +
                                                    (Hx[i_m][j_m][k_m] - Hx[i_m][j_m][k_m - 1]) * P_ce_z[k];
               Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] + G2[medio] * regDU(REGION).Psi_Eyzvac(i, j, k);
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc(iEx).ZI(REGION); k <= PMLc(iEx).ZE(REGION); ++k) {
         k_m = k - b.Ex.ZI;
         for (j = PMLc(iEx).YI(REGION); j <= PMLc(iEx).YE(REGION); ++j) {
            j_m = j - b.Ex.YI;
            for (i = PMLc(iEx).XI(REGION); i <= PMLc(iEx).XE(REGION); ++i) {
               i_m = i - b.Ex.XI;
               medio = 1;
               regDU(REGION).Psi_Exzvac(i, j, k) = P_be_z[k] * regDU(REGION).Psi_Exzvac(i, j, k) +
                                                    (Hy[i_m][j_m][k_m] - Hy[i_m][j_m][k_m - 1]) * P_ce_z[k];
               Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m] - G2[medio] * regDU(REGION).Psi_Exzvac(i, j, k);
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif




      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = back;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc(iEz).ZI(REGION); k <= PMLc(iEz).ZE(REGION); ++k) {
         k_m = k - b.Ez.ZI;
         for (j = PMLc(iEz).YI(REGION); j <= PMLc(iEz).YE(REGION); ++j) {
            j_m = j - b.Ez.YI;
            for (i = PMLc(iEz).XI(REGION); i <= PMLc(iEz).XE(REGION); ++i) {
               i_m = i - b.Ez.XI;
               medio = 1;
               regBF(REGION).Psi_Ezxvac(i, j, k) = P_be_x[i] * regBF(REGION).Psi_Ezxvac(i, j, k) +
                                                    (Hy[i_m][j_m][k_m] - Hy[i_m - 1][j_m][k_m]) * P_ce_x[i];
               Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] + G2[medio] * regBF(REGION).Psi_Ezxvac(i, j, k);
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc(iEy).ZI(REGION); k <= PMLc(iEy).ZE(REGION); ++k) {
         k_m = k - b.Ey.ZI;
         for (j = PMLc(iEy).YI(REGION); j <= PMLc(iEy).YE(REGION); ++j) {
            j_m = j - b.Ey.YI;
            for (i = PMLc(iEy).XI(REGION); i <= PMLc(iEy).XE(REGION); ++i) {
               i_m = i - b.Ey.XI;
               medio = 1;
               regBF(REGION).Psi_Eyxvac(i, j, k) = P_be_x[i] * regBF(REGION).Psi_Eyxvac(i, j, k) +
                                                    (Hz[i_m][j_m][k_m] - Hz[i_m - 1][j_m][k_m]) * P_ce_x[i];
               Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] - G2[medio] * regBF(REGION).Psi_Eyxvac(i, j, k);
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif


      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = front;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc(iEz).ZI(REGION); k <= PMLc(iEz).ZE(REGION); ++k) {
         k_m = k - b.Ez.ZI;
         for (j = PMLc(iEz).YI(REGION); j <= PMLc(iEz).YE(REGION); ++j) {
            j_m = j - b.Ez.YI;
            for (i = PMLc(iEz).XI(REGION); i <= PMLc(iEz).XE(REGION); ++i) {
               i_m = i - b.Ez.XI;
               medio = 1;
               regBF(REGION).Psi_Ezxvac(i, j, k) = P_be_x[i] * regBF(REGION).Psi_Ezxvac(i, j, k) +
                                                    (Hy[i_m][j_m][k_m] - Hy[i_m - 1][j_m][k_m]) * P_ce_x[i];
               Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] + G2[medio] * regBF(REGION).Psi_Ezxvac(i, j, k);
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc(iEy).ZI(REGION); k <= PMLc(iEy).ZE(REGION); ++k) {
         k_m = k - b.Ey.ZI;
         for (j = PMLc(iEy).YI(REGION); j <= PMLc(iEy).YE(REGION); ++j) {
            j_m = j - b.Ey.YI;
            for (i = PMLc(iEy).XI(REGION); i <= PMLc(iEy).XE(REGION); ++i) {
               i_m = i - b.Ey.XI;
               medio = 1;
               regBF(REGION).Psi_Eyxvac(i, j, k) = P_be_x[i] * regBF(REGION).Psi_Eyxvac(i, j, k) +
                                                    (Hz[i_m][j_m][k_m] - Hz[i_m - 1][j_m][k_m]) * P_ce_x[i];
               Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] - G2[medio] * regBF(REGION).Psi_Eyxvac(i, j, k);
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif


      // ---------------------------> acaba AdvanceelectricCPML <---------------------------------------
      return;
   } // end subroutine AdvanceelectricCPML_freespace
   // !**************************************************************************************************
   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!! Advances the magnetic field in the PML
   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   void AdvanceMagneticCPML_freespace(int NumMedia, const bounds_t& b, 
                                      const std::vector<std::vector<std::vector<int>>>& sggMiHx, 
                                      const std::vector<std::vector<std::vector<int>>>& sggMiHy, 
                                      const std::vector<std::vector<std::vector<int>>>& sggMiHz, 
                                      const std::vector<double>& gm2, 
                                      const std::vector<std::vector<std::vector<double>>>& Ex, 
                                      const std::vector<std::vector<std::vector<double>>>& Ey, 
                                      const std::vector<std::vector<std::vector<double>>>& Ez, 
                                      std::vector<std::vector<std::vector<double>>>& Hx, 
                                      std::vector<std::vector<std::vector<double>>>& Hy, 
                                      std::vector<std::vector<std::vector<double>>>& Hz) {
      // ---------------------------> inputs <----------------------------------------------------------
      // integer, intent( IN) :: NumMedia
      // type( bounds_t), intent( IN) :: b
      // !--->
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHx%NX-1, 0 :  b%sggMiHx%NY-1, 0 :  b%sggMiHx%NZ-1), intent( IN) :: sggMiHx
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHy%NX-1, 0 :  b%sggMiHy%NY-1, 0 :  b%sggMiHy%NZ-1), intent( IN) :: sggMiHy
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHz%NX-1, 0 :  b%sggMiHz%NY-1, 0 :  b%sggMiHz%NZ-1), intent( IN) :: sggMiHz
      // !--->
      // real(kind = RKIND), dimension( 0 :  NumMedia), intent( IN) :: gm2
      // !--->
      // real(kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( IN) :: Ex
      // real(kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( IN) :: Ey
      // real(kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( IN) :: Ez
      // ---------------------------> inputs/outputs <--------------------------------------------------
      // real(kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( INOUT) :: Hx
      // real(kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( INOUT) :: Hy
      // real(kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( INOUT) :: Hz
      // ---------------------------> variables locales <-----------------------------------------------
      int REGION, i, j, k, medio, i_m, j_m, k_m;
      // ---------------------------> empieza AdvanceMagneTicCPML <-------------------------------------
      // Hetic Fields PML Zone
      //
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = left;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc(iHx).ZI(REGION); k <= PMLc(iHx).ZE(REGION); ++k) {
         k_m = k - b.Hx.ZI;
         for (j = PMLc(iHx).YI(REGION); j <= PMLc(iHx).YE(REGION); ++j) {
            j_m = j - b.Hx.YI;
            for (i = PMLc(iHx).XI(REGION); i <= PMLc(iHx).XE(REGION); ++i) {
               i_m = i - b.Hx.XI;
               // --->
               regLR(REGION).Psi_Hxyvac(i, j, k) = P_bm_y[j] * regLR(REGION).Psi_Hxyvac(i, j, k) +
                                                    (Ez[i_m][j_m + 1][k_m] - Ez[i_m][j_m][k_m]) * P_cm_y[j];
               medio = 1;
               Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] - GM2[medio] * regLR(REGION).Psi_Hxyvac(i, j, k);
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc(iHz).ZI(REGION); k <= PMLc(iHz).ZE(REGION); ++k) {
         k_m = k - b.Hz.ZI;
         for (j = PMLc(iHz).YI(REGION); j <= PMLc(iHz).YE(REGION); ++j) {
            j_m = j - b.Hz.YI;
            for (i = PMLc(iHz).XI(REGION); i <= PMLc(iHz).XE(REGION); ++i) {
               i_m = i - b.Hz.XI;
               // --->
               regLR(REGION).Psi_Hzyvac(i, j, k) = P_bm_y[j] * regLR(REGION).Psi_Hzyvac(i, j, k) +
                                                    (Ex[i_m][j_m + 1][k_m] - Ex[i_m][j_m][k_m]) * P_cm_y[j];
               medio = 1;
               Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] + GM2[medio] * regLR(REGION).Psi_Hzyvac(i, j, k);
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      REGION = right;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc(iHx).ZI(REGION); k <= PMLc(iHx).ZE(REGION); ++k) {
         k_m = k - b.Hx.ZI;
         for (j = PMLc(iHx).YI(REGION); j <= PMLc(iHx).YE(REGION); ++j) {
            j_m = j - b.Hx.YI;
            for (i = PMLc(iHx).XI(REGION); i <= PMLc(iHx).XE(REGION); ++i) {
               i_m = i - b.Hx.XI;
               // --->
               regLR(REGION).Psi_Hxyvac(i, j, k) = P_bm_y[j] * regLR(REGION).Psi_Hxyvac(i, j, k) +
                                                    (Ez[i_m][j_m + 1][k_m] - Ez[i_m][j_m][k_m]) * P_cm_y[j];
               medio = 1;
               Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] - GM2[medio] * regLR(REGION).Psi_Hxyvac(i, j, k);
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc(iHz).ZI(REGION); k <= PMLc(iHz).ZE(REGION); ++k) {
         k_m = k - b.Hz.ZI;
         for (j = PMLc(iHz).YI(REGION); j <= PMLc(iHz).YE(REGION); ++j) {
            j_m = j - b.Hz.YI;
            for (i = PMLc(iHz).XI(REGION); i <= PMLc(iHz).XE(REGION); ++i) {
               i_m = i - b.Hz.XI;
               // --->
               regLR(REGION).Psi_Hzyvac(i, j, k) = P_bm_y[j] * regLR(REGION).Psi_Hzyvac(i, j, k) +
                                                    (Ex[i_m][j_m + 1][k_m] - Ex[i_m][j_m][k_m]) * P_cm_y[j];
               medio = 1;
               Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] + GM2[medio] * regLR(REGION).Psi_Hzyvac(i, j, k);
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = down;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc(iHy).ZI(REGION); k <= PMLc(iHy).ZE(REGION); ++k) {
         k_m = k - b.Hy.ZI;
         for (j = PMLc(iHy).YI(REGION); j <= PMLc(iHy).YE(REGION); ++j) {
            j_m = j - b.Hy.YI;
            for (i = PMLc(iHy).XI(REGION); i <= PMLc(iHy).XE(REGION); ++i) {
               i_m = i - b.Hy.XI;
               // --->
               regLR(REGION).Psi_Hyzvac(i, j, k) = P_bm_z[k] * regLR(REGION).Psi_Hyzvac(i, j, k) +
                                                    (Ex[i_m][j_m][k_m] - Ex[i_m][j_m][k_m - 1]) * P_cm_z[k];
               medio = 1;
               Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] + GM2[medio] * regLR(REGION).Psi_Hyzvac(i, j, k);
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc(iHx).ZI(REGION); k <= PMLc(iHx).ZE(REGION); ++k) {
         k_m = k - b.Hx.ZI;
         for (j = PMLc(iHx).YI(REGION); j <= PMLc(iHx).YE(REGION); ++j) {
            j_m = j - b.Hx.YI;
            for (i = PMLc(iHx).XI(REGION); i <= PMLc(iHx).XE(REGION); ++i) {
               i_m = i - b.Hx.XI;
               // --->
               regLR(REGION).Psi_Hxzvac(i, j, k) = P_bm_z[k] * regLR(REGION).Psi_Hxzvac(i, j, k) +
                                                    (Ey[i_m][j_m][k_m] - Ey[i_m][j_m][k_m - 1]) * P_cm_z[k];
               medio = 1;
               Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] - GM2[medio] * regLR(REGION).Psi_Hxzvac(i, j, k);
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

   } // end subroutine AdvanceMagneticCPML_freespace

#ifdef CompileWithOpenMP
#include <omp.h>
#endif

// Assuming necessary headers are included from previous chunks:
// #include "BORDERS_CPML_m.hpp"
// #include <vector>
// #include <string>

namespace BORDERS_CPML_m {

    void AdvanceMagneTicCPML_freespace(
        int REGION,
        const std::vector<int>& P_bm_z,
        const std::vector<int>& P_cm_z,
        const std::vector<int>& P_bm_x,
        const std::vector<int>& P_cm_x,
        const std::vector<double>& GM2,
        // Assuming these are passed by reference or are members of a class context.
        // Based on Fortran structure, they seem to be global or passed in a derived type.
        // We will assume they are accessible via some context or passed as references.
        // For this translation, we assume they are available in the scope or passed as arguments.
        // Since the Fortran code accesses them directly, we assume they are global or part of a class.
        // Let's assume they are passed as references to a context object or are global.
        // To be safe and preserve names, we'll treat them as external or passed.
        // However, looking at the usage, they look like global arrays or members of a large state struct.
        // We will assume they are passed as const references or pointers if needed, 
        // but since the prompt asks to preserve names and this is a chunk, 
        // we assume the surrounding context provides them.
        // Let's assume they are available in the namespace or global scope for this snippet.
        
        // Arrays accessed: Ex, Ey, Ez, Hy, Hx, Hz
        // Derived types: PMLc, b, regDU, regBF
        
        // We need to define the types or assume they are defined elsewhere.
        // Since we can't see the definitions, we assume they are declared in the header.
        
        // Note: The Fortran code uses 1-based indexing for loops but accesses arrays with calculated indices.
        // i_m, j_m, k_m are calculated relative to b%Hy%ZI etc.
        // We assume the arrays Ex, Ey, Ez, Hy, Hx, Hz are 1-based or adjusted accordingly.
        // In C++, we typically use 0-based. If the Fortran code is 1-based, the C++ code should reflect that
        // or the arrays should be sized appropriately.
        // The translation rule says "Preserve 1-based indexing where Fortran uses it".
        // So we will assume the arrays are accessed with 1-based indices.
        
        // However, std::vector is 0-based. We can either:
        // 1. Use raw arrays with 1-based indexing (allocate size N+1 and ignore index 0).
        // 2. Adjust indices by -1.
        // The rule says "Preserve 1-based indexing", which implies the logic should remain the same.
        // If the Fortran code uses 1-based indexing, the C++ code should too.
        // This usually means the arrays are allocated with size (max_index + 1) and index 0 is unused.
        
        // We will assume the arrays are passed as std::vector<double>& and are 1-based (index 0 unused).
        
        // We also need to assume the types for PMLc, b, regDU, regBF are defined.
        // We will assume they are passed as references.
        
        // Since this is a chunk, we assume the types are defined in the header.
        
        // Let's assume the following types are defined in the header:
        // struct PML_type;
        // struct Block_type;
        // struct RegionDU_type;
        // struct RegionBF_type;
        
        // And the following variables are available:
        // PML_type PMLc[];
        // Block_type b;
        // RegionDU_type regDU[];
        // RegionBF_type regBF[];
        // double Ex[], Ey[], Ez[], Hy[], Hx[], Hz[];
        
        // For the purpose of this translation, we will assume these are passed as references or are global.
        // To make the code compile, we need to declare them.
        // Since we don't have the full context, we will assume they are passed as arguments.
        
        // However, the Fortran code does not pass them as arguments. They are likely global or module-level.
        // In C++, we can use extern or pass them.
        // Given the instruction "Preserve all original names", we will assume they are available in the namespace.
        
        // We will write the code assuming these variables are declared in the header or global scope.
        
        // Note: The Fortran code has a typo: regBF( region)%Psi_Hyxvac vs regBF( REGION)%Psi_Hyxvac.
        // We preserve the case as is.
        
        // OpenMP parallel loops
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
        for (int k = PMLc[iHy].ZI[REGION]; k <= PMLc[iHy].ZE[REGION]; ++k) {
            i_m = 0; // Placeholder, will be set in inner loop
            j_m = 0;
            k_m = k - b.Hy.ZI;
            for (int j = PMLc[iHy].YI[REGION]; j <= PMLc[iHy].YE[REGION]; ++j) {
                j_m = j - b.Hy.YI;
                for (int i = PMLc[iHy].XI[REGION]; i <= PMLc[iHy].XE[REGION]; ++i) {
                    i_m = i - b.Hy.XI;
                    // --->
                    regDU[REGION].Psi_Hyzvac[i][j][k] = P_bm_z[k] * regDU[REGION].Psi_Hyzvac[i][j][k] +
                        (Ex[i_m][j_m][k_m + 1] - Ex[i_m][j_m][k_m]) * P_cm_z[k];
                    medio = 1;
                    Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] - GM2[medio] * regDU[REGION].Psi_Hyzvac[i][j][k];
                } // bucle i
            }
        }

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
        for (int k = PMLc[iHx].ZI[REGION]; k <= PMLc[iHx].ZE[REGION]; ++k) {
            i_m = 0;
            j_m = 0;
            k_m = k - b.Hx.ZI;
            for (int j = PMLc[iHx].YI[REGION]; j <= PMLc[iHx].YE[REGION]; ++j) {
                j_m = j - b.Hx.YI;
                for (int i = PMLc[iHx].XI[REGION]; i <= PMLc[iHx].XE[REGION]; ++i) {
                    i_m = i - b.Hx.XI;
                    // --->
                    regDU[REGION].Psi_Hxzvac[i][j][k] = P_bm_z[k] * regDU[REGION].Psi_Hxzvac[i][j][k] +
                        (Ey[i_m][j_m][k_m + 1] - Ey[i_m][j_m][k_m]) * P_cm_z[k];
                    medio = 1;
                    Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] + GM2[medio] * regDU[REGION].Psi_Hxzvac[i][j][k];
                } // bucle i
            }
        }

        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REGION = up;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
        for (int k = PMLc[iHy].ZI[REGION]; k <= PMLc[iHy].ZE[REGION]; ++k) {
            i_m = 0;
            j_m = 0;
            k_m = k - b.Hy.ZI;
            for (int j = PMLc[iHy].YI[REGION]; j <= PMLc[iHy].YE[REGION]; ++j) {
                j_m = j - b.Hy.YI;
                for (int i = PMLc[iHy].XI[REGION]; i <= PMLc[iHy].XE[REGION]; ++i) {
                    i_m = i - b.Hy.XI;
                    // --->
                    regDU[REGION].Psi_Hyzvac[i][j][k] = P_bm_z[k] * regDU[REGION].Psi_Hyzvac[i][j][k] +
                        (Ex[i_m][j_m][k_m + 1] - Ex[i_m][j_m][k_m]) * P_cm_z[k];
                    medio = 1;
                    Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] - GM2[medio] * regDU[REGION].Psi_Hyzvac[i][j][k];
                } // bucle i
            }
        }

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
        for (int k = PMLc[iHx].ZI[REGION]; k <= PMLc[iHx].ZE[REGION]; ++k) {
            i_m = 0;
            j_m = 0;
            k_m = k - b.Hx.ZI;
            for (int j = PMLc[iHx].YI[REGION]; j <= PMLc[iHx].YE[REGION]; ++j) {
                j_m = j - b.Hx.YI;
                for (int i = PMLc[iHx].XI[REGION]; i <= PMLc[iHx].XE[REGION]; ++i) {
                    i_m = i - b.Hx.XI;
                    // --->
                    regDU[REGION].Psi_Hxzvac[i][j][k] = P_bm_z[k] * regDU[REGION].Psi_Hxzvac[i][j][k] +
                        (Ey[i_m][j_m][k_m + 1] - Ey[i_m][j_m][k_m]) * P_cm_z[k];
                    medio = 1;
                    Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] + GM2[medio] * regDU[REGION].Psi_Hxzvac[i][j][k];
                } // bucle i
            }
        }

        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REGION = back;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
        for (int k = PMLc[iHz].ZI[REGION]; k <= PMLc[iHz].ZE[REGION]; ++k) {
            i_m = 0;
            j_m = 0;
            k_m = k - b.Hz.ZI;
            for (int j = PMLc[iHz].YI[REGION]; j <= PMLc[iHz].YE[REGION]; ++j) {
                j_m = j - b.Hz.YI;
                for (int i = PMLc[iHz].XI[REGION]; i <= PMLc[iHz].XE[REGION]; ++i) {
                    i_m = i - b.Hz.XI;
                    // --->
                    regBF[REGION].Psi_Hzxvac[i][j][k] = P_bm_x[i] * regBF[REGION].Psi_Hzxvac[i][j][k] +
                        (Ey[i_m + 1][j_m][k_m] - Ey[i_m][j_m][k_m]) * P_cm_x[i];
                    medio = 1;
                    Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] - GM2[medio] * regBF[REGION].Psi_Hzxvac[i][j][k];
                }
            }
        }

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
        for (int k = PMLc[iHy].ZI[REGION]; k <= PMLc[iHy].ZE[REGION]; ++k) {
            i_m = 0;
            j_m = 0;
            k_m = k - b.Hy.ZI;
            for (int j = PMLc[iHy].YI[REGION]; j <= PMLc[iHy].YE[REGION]; ++j) {
                j_m = j - b.Hy.YI;
                for (int i = PMLc[iHy].XI[REGION]; i <= PMLc[iHy].XE[REGION]; ++i) {
                    i_m = i - b.Hy.XI;
                    // --->
                    regBF[region].Psi_Hyxvac[i][j][k] = P_bm_x[i] * regBF[REGION].Psi_Hyxvac[i][j][k] +
                        (Ez[i_m + 1][j_m][k_m] - Ez[i_m][j_m][k_m]) * P_cm_x[i];
                    medio = 1;
                    Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] + GM2[medio] * regBF[REGION].Psi_Hyxvac[i][j][k];
                }
            }
        }

        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REGION = front;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
        for (int k = PMLc[iHz].ZI[REGION]; k <= PMLc[iHz].ZE[REGION]; ++k) {
            i_m = 0;
            j_m = 0;
            k_m = k - b.Hz.ZI;
            for (int j = PMLc[iHz].YI[REGION]; j <= PMLc[iHz].YE[REGION]; ++j) {
                j_m = j - b.Hz.YI;
                for (int i = PMLc[iHz].XI[REGION]; i <= PMLc[iHz].XE[REGION]; ++i) {
                    i_m = i - b.Hz.XI;
                    // --->
                    regBF[REGION].Psi_Hzxvac[i][j][k] = P_bm_x[i] * regBF[REGION].Psi_Hzxvac[i][j][k] +
                        (Ey[i_m + 1][j_m][k_m] - Ey[i_m][j_m][k_m]) * P_cm_x[i];
                    medio = 1;
                    Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] - GM2[medio] * regBF[REGION].Psi_Hzxvac[i][j][k];
                }
            }
        }

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
        for (int k = PMLc[iHy].ZI[REGION]; k <= PMLc[iHy].ZE[REGION]; ++k) {
            i_m = 0;
            j_m = 0;
            k_m = k - b.Hy.ZI;
            for (int j = PMLc[iHy].YI[REGION]; j <= PMLc[iHy].YE[REGION]; ++j) {
                j_m = j - b.Hy.YI;
                for (int i = PMLc[iHy].XI[REGION]; i <= PMLc[iHy].XE[REGION]; ++i) {
                    i_m = i - b.Hy.XI;
                    // --->
                    regBF[region].Psi_Hyxvac[i][j][k] = P_bm_x[i] * regBF[REGION].Psi_Hyxvac[i][j][k] +
                        (Ez[i_m + 1][j_m][k_m] - Ez[i_m][j_m][k_m]) * P_cm_x[i];
                    medio = 1;
                    Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] + GM2[medio] * regBF[REGION].Psi_Hyxvac[i][j][k];
                }
            }
        }

        // ---------------------------> acaba AdvanceMagneTicCPML <---------------------------------------
    }

} // namespace BORDERS_CPML_m