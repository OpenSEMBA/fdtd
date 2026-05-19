#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdint>

struct SGGFDTDINFO_t;
struct limit_t {
        double x_min=0, x_max=0, y_min=0, y_max=0, z_min=0, z_max=0;
    };
struct sim_control_t {
        bool CPML = false;
        bool MUR = false;
        int CPML_type = 0;
    };

using RKIND = double;

enum Direction { Down, Up, Left, Right, Back, Front };
enum FieldIndex { iEx, iEy, iEz, iHx, iHy, iHz };

namespace BORDERS_CPML_m {

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

    std::vector<xyzlimit_var_t> PMLc(6);
    std::vector<LR_t> regLR(6);
    std::vector<DU_t> regDU(6);
    std::vector<BF_t> regBF(6);
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

    void InitCPMLBorders(const SGGFDTDINFO_t&, const std::vector<RKIND>&, const std::vector<RKIND>&, const std::vector<RKIND>&, const std::vector<RKIND>&, const std::vector<RKIND>&, const std::vector<RKIND>&, std::vector<int32_t>&, std::vector<int32_t>&, std::vector<int32_t>&, std::vector<int32_t>&, std::vector<int32_t>&, std::vector<int32_t>&, const std::vector<limit_t>&, const sim_control_t&, bool&, RKIND, RKIND) {}
    void StoreFieldsCPMLBorders() {}
    void calc_cpmlconstants(const SGGFDTDINFO_t&, const std::vector<int32_t>&, const std::vector<int32_t>&, const std::vector<int32_t>&, const std::vector<int32_t>&, const std::vector<int32_t>&, const std::vector<int32_t>&, RKIND, RKIND) {}
    void AdvanceelectricCPML() {}
    void AdvanceMagneticCPML() {}
    void DestroyCPMLBorders() {}
    void AdvanceelectricCPML_freespace() {}
    void AdvanceMagneticCPML_freespace() {}

}
