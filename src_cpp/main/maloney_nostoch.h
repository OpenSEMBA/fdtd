#ifndef SEMBA_FDTD_MALONEY_NOSTOCH_H
#define SEMBA_FDTD_MALONEY_NOSTOCH_H

#include <complex>
#include <cstdint>
#include <vector>

namespace SGBC_nostoch_m {

#ifdef CompileWithReal8
using SGBCReal = double;
#else
using SGBCReal = float;
#endif

struct SGBCLayerDepthResult {
    int32_t depth = 0;
    std::vector<int32_t> capa;
    std::vector<double> delta_entreEinterno;
};

struct SGBCDepthZeroConstants_t {
    double g1 = 1.0;
    double g2a = 0.0;
    double g2b = 0.0;
    double gm2_externo = 0.0;
};

struct SGBCDepthZeroSurface_t {
    bool correct_ha = false;
    bool correct_hb = false;
    int32_t depth = 0;
    double E = 0.0;
    double E_left = 0.0;
    double E_right = 0.0;
    SGBCDepthZeroConstants_t constants;
};

struct SGBCLayer_t {
    double width = 0.0;
    double relativePermittivity = 1.0;
    double relativePermeability = 1.0;
    double electricConductivity = 0.0;
    double magneticConductivity = 0.0;
};

struct SGBCSurface_t {
    bool correct_ha = false;
    bool correct_hb = false;
    bool es_unfilo_placa = false;
    bool SGBCCrank = false;
    int32_t depth = 0;

    SGBCReal Efield = 0.0;
    SGBCReal gm2_externo = 0.0;
    SGBCReal transversalDeltaE = 0.0;
    SGBCReal transversalDeltaH = 0.0;
    SGBCReal alignedDeltaH = 0.0;

    std::vector<int32_t> capa;
    std::vector<SGBCReal> delta_entreEinterno;
    std::vector<SGBCReal> g1;
    std::vector<SGBCReal> g2a;
    std::vector<SGBCReal> g2b;
    std::vector<SGBCReal> E;
    std::vector<SGBCReal> H;
    std::vector<SGBCReal> E_past;
    std::vector<SGBCReal> G2_interno;
    std::vector<SGBCReal> GM2_interno;
    std::vector<SGBCReal> G1_interno;
    std::vector<SGBCReal> GM1_interno;
    SGBCReal Hyee_left = 0.0;
    SGBCReal Hyee_right = 0.0;

    std::vector<SGBCReal> a;
    std::vector<SGBCReal> b;
    std::vector<SGBCReal> c;
    std::vector<SGBCReal> rb;
    std::vector<SGBCReal> rh;
    std::vector<SGBCReal> rhm1;
    SGBCReal a1 = 0.0;
    SGBCReal b1 = 1.0;
    SGBCReal c1 = 0.0;
    SGBCReal rb1 = 0.0;
    SGBCReal rh1 = 0.0;
    SGBCReal an = 0.0;
    SGBCReal bn = 1.0;
    SGBCReal cn = 0.0;
    SGBCReal rbn = 0.0;
    SGBCReal rhn = 0.0;
    std::vector<SGBCReal> D;
    std::vector<SGBCReal> triA;
    std::vector<SGBCReal> triB;
    std::vector<SGBCReal> triC;
    std::vector<SGBCReal> triCp;
    std::vector<SGBCReal> triDp;
    std::vector<SGBCReal> triInvM;
};

struct SGBCHCorrection_t {
    double ha_plus = 0.0;
    double ha_minus = 0.0;
    double hb_plus = 0.0;
    double hb_minus = 0.0;
};

void g1g2(double dt, double epsilon, double sigma, double& G1, double& G2);
void gm1gm2(double dt, double mu, double sigmam, double& Gm1, double& Gm2);
void g1g2_Dispersive(double dt, double epsilon, double sigma,
                     double& G1, double& G2,
                     std::vector<std::complex<double>>& Beta,
                     std::vector<std::complex<double>>& Kappa,
                     std::vector<std::complex<double>>& G3,
                     int32_t numpolres,
                     const std::vector<std::complex<double>>& a11,
                     const std::vector<std::complex<double>>& c11);

SGBCLayerDepthResult calculate_sgbc_layer_depth(
    const std::vector<double>& widths,
    const std::vector<double>& sigmas,
    const std::vector<double>& eprs,
    double SGBCFreq,
    double SGBCresol,
    int32_t SGBCdepth);

SGBCDepthZeroConstants_t calc_depth_zero_sgbc_constants(
    double dt,
    double eps0,
    double mu0,
    double thickness,
    double relativePermittivity,
    double electricConductivity,
    double transversalDeltaE,
    double transversalDeltaH,
    double alignedDeltaH,
    bool correctHa);

SGBCDepthZeroSurface_t make_depth_zero_sgbc_surface(
    const SGBCDepthZeroConstants_t& constants,
    bool correctHa,
    double initialE);

SGBCSurface_t make_sgbc_surface(
    const std::vector<SGBCLayer_t>& layers,
    double dt,
    double eps0,
    double mu0,
    double SGBCFreq,
    double SGBCresol,
    int32_t SGBCdepth,
    bool SGBCCrank,
    bool correctHa,
    bool esUnfiloPlaca,
    double transversalDeltaE,
    double transversalDeltaH,
    double alignedDeltaH,
    double initialE);

void AdvanceSGBCE(SGBCDepthZeroSurface_t& surface,
                  double haPlus,
                  double haMinus,
                  double hbPlus,
                  double hbMinus);

void AdvanceSGBCE(SGBCSurface_t& surface,
                  double haPlus,
                  double haMinus,
                  double hbPlus,
                  double hbMinus);

SGBCHCorrection_t AdvanceSGBCH(const SGBCDepthZeroSurface_t& surface,
                               double eField);

SGBCHCorrection_t AdvanceSGBCH(const SGBCSurface_t& surface,
                               double eField);

void solve_tridiag_distintos(const std::vector<double>& aa,
                             const std::vector<double>& bb,
                             const std::vector<double>& cc,
                             double a1, double b1, double c1,
                             double an, double bn, double cn,
                             const std::vector<double>& d,
                             std::vector<double>& x,
                             int32_t n);
void solve_tridiag_iguales(double aa, double bb, double cc,
                           double a1, double b1, double c1,
                           double an, double bn, double cn,
                           const std::vector<double>& d,
                           std::vector<double>& x,
                           int32_t n);

} // namespace SGBC_nostoch_m

#endif // SEMBA_FDTD_MALONEY_NOSTOCH_H
