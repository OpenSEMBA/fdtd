#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <memory>
#include <stdexcept>

struct SGGFDTDINFO_t;
struct bounds_t;

#ifndef RKIND
#define RKIND double
#endif

#ifndef INTEGERSIZEOFMEDIAMATRICES
#define INTEGERSIZEOFMEDIAMATRICES int
#endif

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
#ifndef iEx
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

void print11(int, const std::string&) {}
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

    std::vector<std::vector<std::vector<xyzlimit_var_t>>> MURc;
    std::vector<LR_t> regLR;
    std::vector<DU_t> regDU;
    std::vector<BF_t> regBF;
    std::vector<RKIND> back_CAB1, back_CAB3, back_cab4;
    std::vector<RKIND> front_CAB1, front_CAB3, front_cab4;
    std::vector<RKIND> left_CAB1, left_CAB3, left_cab4;
    std::vector<RKIND> right_CAB1, right_CAB3, right_cab4;
    std::vector<RKIND> down_CAB1, down_CAB3, down_cab4;
    std::vector<RKIND> up_CAB1, up_CAB3, up_cab4;
    RKIND cluz = 0.0;
    RKIND eps0 = 0.0;
    RKIND mu0 = 0.0;

    void InitMURBorders(const SGGFDTDINFO_t&, bool&, bool, const std::vector<RKIND>&, const std::vector<RKIND>&, const std::vector<RKIND>&, RKIND, RKIND) {}
    void calc_murconstants(const SGGFDTDINFO_t&, const std::vector<RKIND>&, const std::vector<RKIND>&, const std::vector<RKIND>&, RKIND, RKIND) {}
    void AdvanceelectricMUR() {}
    void AdvanceMagneticMUR() {}
    void StoreFieldsMUR() {}
    void DestroyMURBorders() {}
    void AdvanceelectricMUR_freespace() {}
    void AdvanceMagneticMUR_freespace() {}

}
