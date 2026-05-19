#include <vector>
#include <complex>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <memory>

struct media_matrices_t;
struct constants_t;
struct SGGFDTDINFO_t;

using RKIND = double;
using CKIND = std::complex<double>;
using integer4 = int32_t;
using integer8 = int64_t;
using logical = bool;

constexpr int BUFSIZE = 256;
void WarnErrReport(const std::string&, bool = false);
void StopOnError(int, int, const std::string&);
void print11(int, const std::string&);

constexpr int iEx = 0, iEy = 1, iEz = 2;
constexpr int iHx = 3, iHy = 4, iHz = 5;

namespace SGBC_nostoch_m {

    struct val_t {
        std::vector<CKIND> val;
    };

    struct MDfield_t {
        RKIND* FieldPresent;
        RKIND FieldPrevious;
        std::vector<CKIND> Current;
    };

    struct SGBCSurface_t {
        std::vector<RKIND> E;
        std::vector<RKIND> H;
        std::vector<RKIND> E_past;
        RKIND* Efield;
        RKIND* Ha_Plus;
        RKIND* Ha_Minu;
        RKIND* Hb_Plus;
        RKIND* Hb_Minu;
        std::vector<RKIND> delta_entreEinterno;
        std::vector<RKIND> g1;
        std::vector<RKIND> g2a;
        std::vector<RKIND> g2b;
        std::vector<MDfield_t> EDis;
        integer4 numpolres;
        val_t Beta;
        val_t Kappa;
        val_t G3;
        logical correct_ha;
        logical correct_hb;
        logical es_unfilo_placa;
        integer4 depth;
        integer4 jmed;
        std::vector<integer4> capa;
        std::vector<RKIND> G2_interno;
        std::vector<RKIND> GM2_interno;
        std::vector<RKIND> G1_interno;
        std::vector<RKIND> GM1_interno;
        RKIND GM2_externo;
        RKIND Hyee__left;
        RKIND Hyee_right;
        std::vector<RKIND> a;
        std::vector<RKIND> b;
        std::vector<RKIND> c;
        std::vector<RKIND> rb;
        std::vector<RKIND> rh;
        std::vector<RKIND> rhm1;
        RKIND a1;
        RKIND b1;
        RKIND c1;
        RKIND rb1;
        RKIND rh1;
        RKIND an;
        RKIND bn;
        RKIND cn;
        RKIND rbn;
        RKIND rhn;
        std::vector<RKIND> D;
        logical SGBCCrank;
        RKIND transversalDeltaE;
        RKIND transversalDeltaH;
        RKIND alignedlDeltaH;
        std::vector<integer4> med;
        std::vector<CKIND> a11;
        std::vector<CKIND> c11;
        RKIND& g1_0() { return g1[0]; }
        RKIND& g1_1() { return g1[1]; }
        RKIND& g2a_0() { return g2a[0]; }
        RKIND& g2a_1() { return g2a[1]; }
        RKIND& g2b_0() { return g2b[0]; }
        RKIND& g2b_1() { return g2b[1]; }
        integer4& med_0() { return med[0]; }
        integer4& med_1() { return med[1]; }
    };

    struct MalDisp_t {
        integer4 numpolres;
        std::vector<CKIND> a11;
        std::vector<CKIND> c11;
    };

    struct Malon_t {
        logical SGBCDispersive;
        integer4 NumNodes;
        std::vector<SGBCSurface_t> nodes;
        std::vector<MalDisp_t> mediosDis;
    };

    Malon_t malon;
    RKIND eps0 = 8.854187817e-12;
    RKIND mu0 = 1.25663706212e-6;
    RKIND zvac = 376.730313461;
    RKIND cluz = 299792458.0;
    logical SGBCcrank = false;
    logical SGBCDispersive = false;
    RKIND SGBCFreq = 0.0;
    RKIND SGBCresol = 0.0;
    integer4 SGBCdepth = 0;

    void InitSGBCs(
        SGGFDTDINFO_t&, const media_matrices_t&,
        std::vector<std::vector<std::vector<RKIND>>>&,
        std::vector<std::vector<std::vector<RKIND>>>&,
        std::vector<std::vector<std::vector<RKIND>>>&,
        std::vector<std::vector<std::vector<RKIND>>>&,
        std::vector<std::vector<std::vector<RKIND>>>&,
        std::vector<std::vector<std::vector<RKIND>>>&,
        std::vector<RKIND>&, std::vector<RKIND>&, std::vector<RKIND>&,
        std::vector<RKIND>&, std::vector<RKIND>&, std::vector<RKIND>&,
        integer4, integer4, constants_t&, SGGFDTDINFO_t&,
        logical, logical, RKIND, RKIND, logical,
        RKIND, RKIND, RKIND, integer4, logical, logical&) {}
    void calc_SGBCconstants(SGGFDTDINFO_t&, constants_t&, RKIND, RKIND, logical) {}
    void AdvanceSGBCE(RKIND, logical, logical, logical) {}
    void AdvanceSGBCH() {}
    void StoreFieldsSGBCs(logical) {}
    void DestroySGBCs(SGGFDTDINFO_t&) {}
    void calc_g1g2gm1gm2_compo(SGGFDTDINFO_t&, SGBCSurface_t&, RKIND, RKIND, logical) {}
    void g1g2(RKIND, RKIND, RKIND, RKIND&, RKIND&) {}
    void gm1gm2(RKIND, RKIND, RKIND, RKIND&, RKIND&) {}
    void g1g2_Dispersive(RKIND, RKIND, RKIND, RKIND&, RKIND&, std::vector<CKIND>&, std::vector<CKIND>&, std::vector<CKIND>&, integer4, const std::vector<CKIND>&, const std::vector<CKIND>&) {}
    void depth(SGBCSurface_t&, SGGFDTDINFO_t&, integer4, RKIND, RKIND, integer4) {}
    Malon_t* GetSGBCs() { return &malon; }
    void solve_tridiag_distintos(const std::vector<RKIND>&, const std::vector<RKIND>&, const std::vector<RKIND>&, RKIND, RKIND, RKIND, RKIND, RKIND, RKIND, const std::vector<RKIND>&, std::vector<RKIND>&, integer4) {}
    void solve_tridiag_iguales(RKIND, RKIND, RKIND, RKIND, RKIND, RKIND, RKIND, RKIND, RKIND, const std::vector<RKIND>&, std::vector<RKIND>&, integer4) {}
    void test_stab() {}

}
