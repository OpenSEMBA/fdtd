#include <vector>
#include <complex>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <memory>

#include "maloney_nostoch.h"

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

    struct LegacySGBCSurface_t {
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
        std::vector<LegacySGBCSurface_t> nodes;
        std::vector<MalDisp_t> mediosDis;
    };

    Malon_t malon;
    RKIND eps0 = 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12;
    RKIND mu0 = 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6;
    RKIND zvac = 376.730313461;
    RKIND cluz = 299792458.0;
    constexpr RKIND Pi = 3.141592653589793238462643383279502884;
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
    void AdvanceSGBCE(RKIND, logical, logical, logical) {
        for (auto& compo : malon.nodes) {
            if (compo.depth != 0 || compo.E.empty() || compo.g1.empty() ||
                compo.g2a.empty() || compo.g2b.empty() ||
                compo.Ha_Plus == nullptr || compo.Ha_Minu == nullptr ||
                compo.Hb_Plus == nullptr || compo.Hb_Minu == nullptr) {
                continue;
            }
            compo.E[0] =
                compo.g1[0] * compo.E[0] +
                compo.g2a[0] * (*compo.Ha_Plus - *compo.Ha_Minu) -
                compo.g2b[0] * (*compo.Hb_Plus - *compo.Hb_Minu);
            if (compo.Efield != nullptr) {
                *compo.Efield = compo.E[0];
            }
        }
    }
    void AdvanceSGBCH() {
        for (auto& compo : malon.nodes) {
            if (compo.depth != 0 || compo.E.empty() ||
                compo.Efield == nullptr || compo.Ha_Plus == nullptr ||
                compo.Ha_Minu == nullptr || compo.Hb_Plus == nullptr ||
                compo.Hb_Minu == nullptr) {
                continue;
            }
            if (compo.correct_ha) {
                *compo.Ha_Plus += compo.GM2_externo * (*compo.Efield - compo.E[0]);
                *compo.Ha_Minu -= compo.GM2_externo * (*compo.Efield - compo.E[0]);
            } else if (compo.correct_hb) {
                *compo.Hb_Plus -= compo.GM2_externo * (*compo.Efield - compo.E[0]);
                *compo.Hb_Minu += compo.GM2_externo * (*compo.Efield - compo.E[0]);
            }
        }
    }
    void StoreFieldsSGBCs(logical) {}
    void DestroySGBCs(SGGFDTDINFO_t&) {}
    void calc_g1g2gm1gm2_compo(SGGFDTDINFO_t&, LegacySGBCSurface_t&, RKIND, RKIND, logical) {}
    void g1g2(RKIND dt, RKIND epsilon, RKIND sigma, RKIND& G1, RKIND& G2) {
        G1 = (1.0 - sigma * dt / (2.0 * epsilon)) /
             (1.0 + sigma * dt / (2.0 * epsilon));
        G2 = dt / epsilon /
             (1.0 + sigma * dt / (2.0 * epsilon));
        if (G1 < 0.0) {
            G1 = std::exp(-sigma * dt / epsilon);
            G2 = (1.0 - G1) / sigma;
        }
    }

    void gm1gm2(RKIND dt, RKIND mu, RKIND sigmam, RKIND& Gm1, RKIND& Gm2) {
        Gm1 = (1.0 - sigmam * dt / (2.0 * mu)) /
              (1.0 + sigmam * dt / (2.0 * mu));
        Gm2 = dt / mu /
              (1.0 + sigmam * dt / (2.0 * mu));
        if (Gm1 < 0.0) {
            Gm1 = std::exp(-sigmam * dt / mu);
            Gm2 = (1.0 - Gm1) / sigmam;
        }
    }

    SGBCDepthZeroConstants_t calc_depth_zero_sgbc_constants(
        RKIND dt,
        RKIND eps0Value,
        RKIND mu0Value,
        RKIND thickness,
        RKIND relativePermittivity,
        RKIND electricConductivity,
        RKIND transversalDeltaE,
        RKIND transversalDeltaH,
        RKIND alignedDeltaH,
        bool correctHa) {
        const RKIND sigma =
            electricConductivity * thickness / transversalDeltaH;
        const RKIND epsilon =
            (eps0Value * (transversalDeltaH - thickness) +
             relativePermittivity * eps0Value * thickness) /
            transversalDeltaH;

        RKIND g1 = 1.0;
        RKIND g2 = 0.0;
        g1g2(dt, epsilon, sigma, g1, g2);

        SGBCDepthZeroConstants_t constants;
        constants.g1 = g1;
        if (correctHa) {
            constants.g2a = g2 / transversalDeltaH;
            constants.g2b = g2 / alignedDeltaH;
        } else {
            constants.g2a = g2 / alignedDeltaH;
            constants.g2b = g2 / transversalDeltaH;
        }
        constants.gm2_externo = (dt / mu0Value) / transversalDeltaE;
        return constants;
    }

    SGBCDepthZeroSurface_t make_depth_zero_sgbc_surface(
        const SGBCDepthZeroConstants_t& constants,
        bool correctHa,
        RKIND initialE) {
        SGBCDepthZeroSurface_t surface;
        surface.correct_ha = correctHa;
        surface.correct_hb = !correctHa;
        surface.depth = 0;
        surface.E = initialE;
        surface.E_left = initialE;
        surface.E_right = initialE;
        surface.constants = constants;
        return surface;
    }

    void AdvanceSGBCE(SGBCDepthZeroSurface_t& surface,
                      RKIND haPlus,
                      RKIND haMinus,
                      RKIND hbPlus,
                      RKIND hbMinus) {
        surface.E =
            surface.constants.g1 * surface.E +
            surface.constants.g2a * (haPlus - haMinus) -
            surface.constants.g2b * (hbPlus - hbMinus);
        surface.E_left = surface.E;
        surface.E_right = surface.E;
    }

    SGBCHCorrection_t AdvanceSGBCH(const SGBCDepthZeroSurface_t& surface,
                                   RKIND eField) {
        SGBCHCorrection_t correction;
        if (surface.correct_ha) {
            correction.ha_plus =
                surface.constants.gm2_externo * (eField - surface.E_right);
            correction.ha_minus =
                -surface.constants.gm2_externo * (eField - surface.E_left);
        } else if (surface.correct_hb) {
            correction.hb_plus =
                -surface.constants.gm2_externo * (eField - surface.E_right);
            correction.hb_minus =
                surface.constants.gm2_externo * (eField - surface.E_left);
        }
        return correction;
    }

    static size_t e_index(const SGBCSurface_t& surface, integer4 i) {
        return static_cast<size_t>(i + surface.depth);
    }

    static size_t h_index(const SGBCSurface_t& surface, integer4 i) {
        return static_cast<size_t>(i + surface.depth);
    }

    static SGBCReal sgbc_real(double value) {
        return static_cast<SGBCReal>(value);
    }

    static void g1g2_sgbc(double dt,
                          SGBCReal epsilon,
                          SGBCReal sigma,
                          SGBCReal& G1,
                          SGBCReal& G2) {
        const double epsilonValue = static_cast<double>(epsilon);
        const double sigmaValue = static_cast<double>(sigma);
        const double g1Value =
            (1.0 - sigmaValue * dt / (2.0 * epsilonValue)) /
            (1.0 + sigmaValue * dt / (2.0 * epsilonValue));
        double g2Value =
            dt / epsilonValue /
            (1.0 + sigmaValue * dt / (2.0 * epsilonValue));
        G1 = sgbc_real(g1Value);
        if (G1 < SGBCReal{0.0}) {
            G1 = sgbc_real(std::exp(-sigmaValue * dt / epsilonValue));
            g2Value = static_cast<double>((SGBCReal{1.0} - G1) / sigma);
        }
        G2 = sgbc_real(g2Value);
    }

    static void gm1gm2_sgbc(double dt,
                            SGBCReal mu,
                            SGBCReal sigmam,
                            SGBCReal& GM1,
                            SGBCReal& GM2) {
        const double muValue = static_cast<double>(mu);
        const double sigmamValue = static_cast<double>(sigmam);
        const double gm1Value =
            (1.0 - sigmamValue * dt / (2.0 * muValue)) /
            (1.0 + sigmamValue * dt / (2.0 * muValue));
        double gm2Value =
            dt / muValue /
            (1.0 + sigmamValue * dt / (2.0 * muValue));
        GM1 = sgbc_real(gm1Value);
        if (GM1 < SGBCReal{0.0}) {
            GM1 = sgbc_real(std::exp(-sigmamValue * dt / muValue));
            gm2Value = static_cast<double>((SGBCReal{1.0} - GM1) / sigmam);
        }
        GM2 = sgbc_real(gm2Value);
    }

    static const SGBCLayer_t& layer_at(const std::vector<SGBCLayer_t>& layers,
                                       integer4 oneBasedIndex) {
        if (oneBasedIndex < 1 ||
            static_cast<size_t>(oneBasedIndex) > layers.size()) {
            throw std::out_of_range("SGBC layer index out of range");
        }
        return layers[static_cast<size_t>(oneBasedIndex - 1)];
    }

    static void precompute_tridiag_sgbc(const std::vector<SGBCReal>& a,
                                        const std::vector<SGBCReal>& b,
                                        const std::vector<SGBCReal>& c,
                                        std::vector<SGBCReal>& cp,
                                        std::vector<SGBCReal>& invM,
                                        integer4 n) {
        if (n <= 0) {
            cp.clear();
            invM.clear();
            return;
        }
        const auto count = static_cast<size_t>(n);
        if (cp.size() != count) cp.resize(count);
        if (invM.size() != count) invM.resize(count);
        invM[0] = SGBCReal{1.0} / b[0];
        cp[0] = c[0] * invM[0];
        for (size_t i = 1; i < count; ++i) {
            const SGBCReal m = b[i] - cp[i - 1] * a[i];
            invM[i] = SGBCReal{1.0} / m;
            cp[i] = c[i] * invM[i];
        }
    }

    static void solve_tridiag_prefactored_sgbc(
        const std::vector<SGBCReal>& a,
        const std::vector<SGBCReal>& cp,
        const std::vector<SGBCReal>& invM,
        const std::vector<SGBCReal>& d,
        std::vector<SGBCReal>& x,
        std::vector<SGBCReal>& dp,
        integer4 n) {
        if (n <= 0) {
            x.clear();
            return;
        }
        const auto count = static_cast<size_t>(n);
        if (dp.size() != count) dp.resize(count);
        if (x.size() != count) x.resize(count);
        dp[0] = d[0] * invM[0];
        for (size_t i = 1; i < count; ++i) {
            dp[i] = (d[i] - dp[i - 1] * a[i]) * invM[i];
        }
        x[count - 1] = dp[count - 1];
        for (size_t i = count - 1; i-- > 0;) {
            x[i] = dp[i] - cp[i] * x[i + 1];
        }
    }

    static void calculate_sgbc_constants(SGBCSurface_t& surface,
                                         const std::vector<SGBCLayer_t>& layers,
                                         RKIND dt,
                                         RKIND eps0Value,
                                         RKIND mu0Value) {
        SGBCReal gm1External = 1.0;
        SGBCReal gm2External = 0.0;
        gm1gm2_sgbc(dt, sgbc_real(mu0Value), SGBCReal{0.0},
                    gm1External, gm2External);
        surface.gm2_externo = gm2External / surface.transversalDeltaE;

        if (surface.depth == 0) {
            const auto& layer = layers.front();
            const SGBCReal width = sgbc_real(layer.width);
            const SGBCReal sigmatemp = sgbc_real(layer.electricConductivity);
            const SGBCReal eprtemp = sgbc_real(layer.relativePermittivity);
            SGBCReal epsilon = 0.0;
            SGBCReal sigma = 0.0;
            if (surface.es_unfilo_placa) {
                epsilon =
                    (sgbc_real(eps0Value) * (surface.transversalDeltaH - width / SGBCReal{2.0}) +
                     eprtemp * sgbc_real(eps0Value) * width / SGBCReal{2.0}) /
                    surface.transversalDeltaH;
                sigma =
                    (sigmatemp * width / SGBCReal{2.0}) / surface.transversalDeltaH;
            } else {
                epsilon =
                    (sgbc_real(eps0Value) * (surface.transversalDeltaH - width) +
                     eprtemp * sgbc_real(eps0Value) * width) /
                    surface.transversalDeltaH;
                sigma = sigmatemp * width / surface.transversalDeltaH;
            }

            SGBCReal g1Value = 1.0;
            SGBCReal g2Value = 0.0;
            g1g2_sgbc(dt, epsilon, sigma, g1Value, g2Value);
            surface.g1[0] = g1Value;
            surface.g1[1] = g1Value;
            if (surface.correct_ha) {
                surface.g2a[0] = g2Value / surface.transversalDeltaH;
                surface.g2b[0] = g2Value / surface.alignedDeltaH;
            } else {
                surface.g2a[0] = g2Value / surface.alignedDeltaH;
                surface.g2b[0] = g2Value / surface.transversalDeltaH;
            }
            surface.g2a[1] = surface.g2a[0];
            surface.g2b[1] = surface.g2b[0];

            g1g2_sgbc(dt, eprtemp * sgbc_real(eps0Value), sigmatemp,
                      g1Value, g2Value);
            surface.G1_interno[0] = g1Value;
            surface.G2_interno[0] = g2Value;
            SGBCReal gm1Value = 1.0;
            SGBCReal gm2Value = 0.0;
            gm1gm2_sgbc(dt, sgbc_real(mu0Value), SGBCReal{0.0},
                        gm1Value, gm2Value);
            surface.GM1_interno[0] = gm1Value;
            surface.GM2_interno[0] = gm2Value;
            return;
        }

        for (integer4 side = 0; side <= 1; ++side) {
            const integer4 layerIndex =
                (side == 0) ? 1 : static_cast<integer4>(layers.size());
            const auto& layer = layer_at(layers, layerIndex);
            const SGBCReal delta =
                (side == 0)
                    ? surface.delta_entreEinterno[h_index(surface, -surface.depth)]
                    : surface.delta_entreEinterno[h_index(surface, surface.depth - 1)];
            SGBCReal epsilon = 0.0;
            SGBCReal sigma = 0.0;
            const SGBCReal eprtemp = sgbc_real(layer.relativePermittivity);
            const SGBCReal sigmatemp = sgbc_real(layer.electricConductivity);
            if (surface.es_unfilo_placa) {
                epsilon =
                    (sgbc_real(eps0Value) * (surface.transversalDeltaH + delta / SGBCReal{2.0}) +
                     eprtemp * sgbc_real(eps0Value) * (delta / SGBCReal{2.0})) /
                    (surface.transversalDeltaH + delta);
                sigma =
                    (sigmatemp * (delta / SGBCReal{2.0})) /
                    (surface.transversalDeltaH + delta);
            } else {
                epsilon =
                    (sgbc_real(eps0Value) * surface.transversalDeltaH +
                     eprtemp * sgbc_real(eps0Value) * delta) /
                    (surface.transversalDeltaH + delta);
                sigma =
                    (sigmatemp * delta) /
                    (surface.transversalDeltaH + delta);
            }

            SGBCReal g1Value = 1.0;
            SGBCReal g2Value = 0.0;
            g1g2_sgbc(dt, epsilon, sigma, g1Value, g2Value);
            surface.g1[static_cast<size_t>(side)] = g1Value;
            if (surface.correct_ha) {
                surface.g2a[static_cast<size_t>(side)] =
                    g2Value / (SGBCReal{0.5} * surface.transversalDeltaH + SGBCReal{0.5} * delta);
                surface.g2b[static_cast<size_t>(side)] =
                    g2Value / surface.alignedDeltaH;
            } else {
                surface.g2a[static_cast<size_t>(side)] =
                    g2Value / surface.alignedDeltaH;
                surface.g2b[static_cast<size_t>(side)] =
                    g2Value / (SGBCReal{0.5} * surface.transversalDeltaH + SGBCReal{0.5} * delta);
            }
        }

        for (integer4 i = -surface.depth + 1; i <= surface.depth - 1; ++i) {
            const auto& layer = layer_at(
                layers, surface.capa[h_index(surface, i)]);
            const auto& adjacentLayer = layer_at(
                layers, surface.capa[h_index(surface, i - 1)]);
            const SGBCReal deltaLeft =
                surface.delta_entreEinterno[h_index(surface, i - 1)];
            const SGBCReal deltaRight =
                surface.delta_entreEinterno[h_index(surface, i)];
            const SGBCReal eprtemp =
                (sgbc_real(adjacentLayer.relativePermittivity) * deltaLeft +
                 sgbc_real(layer.relativePermittivity) * deltaRight) /
                (deltaLeft + deltaRight);
            const SGBCReal sigmatemp =
                (sgbc_real(adjacentLayer.electricConductivity) * deltaLeft +
                 sgbc_real(layer.electricConductivity) * deltaRight) /
                (deltaLeft + deltaRight);
            SGBCReal g1Value = 1.0;
            SGBCReal g2Value = 0.0;
            g1g2_sgbc(dt, eprtemp * sgbc_real(eps0Value), sigmatemp,
                      g1Value, g2Value);
            surface.G1_interno[e_index(surface, i)] = g1Value;
            surface.G2_interno[e_index(surface, i)] =
                g2Value / ((deltaRight + deltaLeft) / SGBCReal{2.0});
        }

        for (integer4 i = -surface.depth; i <= surface.depth - 1; ++i) {
            const auto& layer = layer_at(
                layers, surface.capa[h_index(surface, i)]);
            SGBCReal gm1Value = 1.0;
            SGBCReal gm2Value = 0.0;
            gm1gm2_sgbc(dt,
                        sgbc_real(layer.relativePermeability) * sgbc_real(mu0Value),
                        sgbc_real(layer.magneticConductivity),
                        gm1Value,
                        gm2Value);
            surface.GM1_interno[h_index(surface, i)] = gm1Value;
            surface.GM2_interno[h_index(surface, i)] =
                gm2Value / surface.delta_entreEinterno[h_index(surface, i)];
        }

        const SGBCReal signo = surface.correct_ha ? SGBCReal{1.0} : SGBCReal{-1.0};
        const SGBCReal g1eff_0 = surface.g1[0];
        const SGBCReal g1eff_1 = surface.g1[1];
        const SGBCReal g2eff_0 =
            signo * (surface.correct_ha ? surface.g2a[0] : surface.g2b[0]);
        const SGBCReal g2eff_1 =
            signo * (surface.correct_ha ? surface.g2a[1] : surface.g2b[1]);

        for (integer4 i = -surface.depth; i <= surface.depth - 1; ++i) {
            surface.GM2_interno[h_index(surface, i)] =
                signo * surface.GM2_interno[h_index(surface, i)];
        }
        for (integer4 i = -surface.depth + 1; i <= surface.depth - 1; ++i) {
            surface.G2_interno[e_index(surface, i)] =
                signo * surface.G2_interno[e_index(surface, i)];
        }

        for (integer4 i = -surface.depth + 1; i <= surface.depth - 1; ++i) {
            surface.a[e_index(surface, i)] =
                -surface.G2_interno[e_index(surface, i)] *
                surface.GM2_interno[h_index(surface, i - 1)] / SGBCReal{4.0};
            surface.b[e_index(surface, i)] =
                SGBCReal{1.0} +
                surface.G2_interno[e_index(surface, i)] *
                    surface.GM2_interno[h_index(surface, i - 1)] / SGBCReal{4.0} +
                surface.G2_interno[e_index(surface, i)] *
                    surface.GM2_interno[h_index(surface, i)] / SGBCReal{4.0};
            surface.c[e_index(surface, i)] =
                -surface.G2_interno[e_index(surface, i)] *
                surface.GM2_interno[h_index(surface, i)] / SGBCReal{4.0};
            surface.rb[e_index(surface, i)] =
                surface.G1_interno[e_index(surface, i)] -
                surface.G2_interno[e_index(surface, i)] *
                    surface.GM2_interno[h_index(surface, i - 1)] / SGBCReal{4.0} -
                surface.G2_interno[e_index(surface, i)] *
                    surface.GM2_interno[h_index(surface, i)] / SGBCReal{4.0};
            surface.rh[e_index(surface, i)] =
                (surface.G2_interno[e_index(surface, i)] *
                     surface.GM1_interno[h_index(surface, i)] +
                 surface.G2_interno[e_index(surface, i)]) /
                SGBCReal{2.0};
            surface.rhm1[e_index(surface, i)] =
                (surface.G2_interno[e_index(surface, i)] *
                     surface.GM1_interno[h_index(surface, i - 1)] +
                 surface.G2_interno[e_index(surface, i)]) /
                SGBCReal{2.0};
        }

        integer4 i = -surface.depth;
        surface.a1 = 0.0;
        surface.c1 =
            -g2eff_0 * surface.GM2_interno[h_index(surface, i)] / SGBCReal{4.0};
        surface.b1 =
            SGBCReal{1.0} + g2eff_0 * surface.GM2_interno[h_index(surface, i)] / SGBCReal{4.0};
        surface.rb1 =
            g1eff_0 - g2eff_0 * surface.GM2_interno[h_index(surface, i)] / SGBCReal{4.0};
        surface.rh1 =
            (g2eff_0 * surface.GM1_interno[h_index(surface, i)] + g2eff_0) /
            SGBCReal{2.0};

        i = surface.depth;
        surface.cn = 0.0;
        surface.an =
            -g2eff_1 * surface.GM2_interno[h_index(surface, i - 1)] / SGBCReal{4.0};
        surface.bn =
            SGBCReal{1.0} + g2eff_1 * surface.GM2_interno[h_index(surface, i - 1)] / SGBCReal{4.0};
        surface.rbn =
            g1eff_1 - g2eff_1 * surface.GM2_interno[h_index(surface, i - 1)] / SGBCReal{4.0};
        surface.rhn =
            (g2eff_1 * surface.GM1_interno[h_index(surface, i - 1)] +
             g2eff_1) /
            SGBCReal{2.0};
    }

    SGBCSurface_t make_sgbc_surface(
        const std::vector<SGBCLayer_t>& layers,
        RKIND dtValue,
        RKIND eps0Value,
        RKIND mu0Value,
        RKIND SGBCFreqValue,
        RKIND SGBCresolValue,
        integer4 SGBCdepthValue,
        bool SGBCCrankValue,
        bool correctHa,
        bool esUnfiloPlaca,
        RKIND transversalDeltaE,
        RKIND transversalDeltaH,
        RKIND alignedDeltaH,
        RKIND initialE) {
        if (layers.empty()) {
            throw std::invalid_argument("SGBC surface requires at least one layer");
        }
        std::vector<RKIND> widths;
        std::vector<RKIND> sigmas;
        std::vector<RKIND> eprs;
        widths.reserve(layers.size());
        sigmas.reserve(layers.size());
        eprs.reserve(layers.size());
        for (const auto& layer : layers) {
            widths.push_back(layer.width);
            sigmas.push_back(layer.electricConductivity);
            eprs.push_back(layer.relativePermittivity);
        }

        const auto depthResult = calculate_sgbc_layer_depth(
            widths, sigmas, eprs, SGBCFreqValue, SGBCresolValue,
            SGBCdepthValue);

        SGBCSurface_t surface;
        surface.correct_ha = correctHa;
        surface.correct_hb = !correctHa;
        surface.es_unfilo_placa = esUnfiloPlaca;
        surface.depth = depthResult.depth;
        surface.SGBCCrank = SGBCCrankValue && surface.depth >= 2;
        surface.transversalDeltaE = transversalDeltaE;
        surface.transversalDeltaH = transversalDeltaH;
        surface.alignedDeltaH = alignedDeltaH;
        surface.capa = depthResult.capa;
        surface.delta_entreEinterno.assign(
            depthResult.delta_entreEinterno.begin(),
            depthResult.delta_entreEinterno.end());

        const size_t eCount = static_cast<size_t>(2 * surface.depth + 1);
        const size_t hCount =
            (surface.depth > 0) ? static_cast<size_t>(2 * surface.depth) : 1u;
        surface.g1.assign(2, 1.0);
        surface.g2a.assign(2, 0.0);
        surface.g2b.assign(2, 0.0);
        surface.E.assign(eCount, sgbc_real(initialE));
        surface.H.assign(hCount, SGBCReal{0.0});
        surface.E_past.assign(eCount, sgbc_real(initialE));
        surface.G1_interno.assign(eCount, SGBCReal{1.0});
        surface.G2_interno.assign(eCount, SGBCReal{0.0});
        surface.GM1_interno.assign(hCount, SGBCReal{1.0});
        surface.GM2_interno.assign(hCount, SGBCReal{0.0});
        surface.a.assign(eCount, SGBCReal{0.0});
        surface.b.assign(eCount, SGBCReal{1.0});
        surface.c.assign(eCount, SGBCReal{0.0});
        surface.rb.assign(eCount, SGBCReal{0.0});
        surface.rh.assign(eCount, SGBCReal{0.0});
        surface.rhm1.assign(eCount, SGBCReal{0.0});
        surface.D.assign(eCount, SGBCReal{0.0});
        surface.Efield = sgbc_real(initialE);

        calculate_sgbc_constants(surface, layers, dtValue, eps0Value, mu0Value);
        if (surface.SGBCCrank && surface.depth > 0) {
            surface.triA.assign(eCount, SGBCReal{0.0});
            surface.triB.assign(eCount, SGBCReal{1.0});
            surface.triC.assign(eCount, SGBCReal{0.0});
            for (integer4 i = -surface.depth + 1; i <= surface.depth - 1; ++i) {
                const size_t idx = e_index(surface, i);
                surface.triA[idx] = surface.a[idx];
                surface.triB[idx] = surface.b[idx];
                surface.triC[idx] = surface.c[idx];
            }
            surface.triA[e_index(surface, -surface.depth)] = surface.a1;
            surface.triB[e_index(surface, -surface.depth)] = surface.b1;
            surface.triC[e_index(surface, -surface.depth)] = surface.c1;
            surface.triA[e_index(surface, surface.depth)] = surface.an;
            surface.triB[e_index(surface, surface.depth)] = surface.bn;
            surface.triC[e_index(surface, surface.depth)] = surface.cn;
            surface.triCp.assign(eCount, SGBCReal{0.0});
            surface.triDp.assign(eCount, SGBCReal{0.0});
            surface.triInvM.assign(eCount, SGBCReal{0.0});
            precompute_tridiag_sgbc(
                surface.triA, surface.triB, surface.triC,
                surface.triCp, surface.triInvM,
                static_cast<integer4>(eCount));
        }
        return surface;
    }

    void AdvanceSGBCE(SGBCSurface_t& surface,
                      RKIND haPlus,
                      RKIND haMinus,
                      RKIND hbPlus,
                      RKIND hbMinus) {
        const SGBCReal haPlusValue = sgbc_real(haPlus);
        const SGBCReal haMinusValue = sgbc_real(haMinus);
        const SGBCReal hbPlusValue = sgbc_real(hbPlus);
        const SGBCReal hbMinusValue = sgbc_real(hbMinus);
        const integer4 depth = surface.depth;
        if (depth > 0) {
            if (!surface.SGBCCrank) {
                if (surface.correct_ha) {
                    surface.E[e_index(surface, depth)] =
                        surface.g1[1] * surface.E[e_index(surface, depth)] +
                        surface.g2a[1] * (haPlusValue - surface.Hyee_right) -
                        surface.g2b[1] * (hbPlusValue - hbMinusValue);
                    surface.E[e_index(surface, -depth)] =
                        surface.g1[0] * surface.E[e_index(surface, -depth)] +
                        surface.g2a[0] * (surface.Hyee_left - haMinusValue) -
                        surface.g2b[0] * (hbPlusValue - hbMinusValue);
                } else if (surface.correct_hb) {
                    surface.E[e_index(surface, depth)] =
                        surface.g1[1] * surface.E[e_index(surface, depth)] +
                        surface.g2a[1] * (haPlusValue - haMinusValue) -
                        surface.g2b[1] * (hbPlusValue - surface.Hyee_right);
                    surface.E[e_index(surface, -depth)] =
                        surface.g1[0] * surface.E[e_index(surface, -depth)] +
                        surface.g2a[0] * (haPlusValue - haMinusValue) -
                        surface.g2b[0] * (surface.Hyee_left - hbMinusValue);
                }
            }
        } else {
            surface.E[e_index(surface, 0)] =
                surface.g1[0] * surface.E[e_index(surface, 0)] +
                surface.g2a[0] * (haPlusValue - haMinusValue) -
                surface.g2b[0] * (hbPlusValue - hbMinusValue);
        }

        if (surface.SGBCCrank) {
            const auto& oldE = surface.E;
            auto& newE = surface.E_past;
            if (surface.correct_ha) {
                integer4 i = depth;
                surface.D[e_index(surface, i)] =
                    -surface.an * oldE[e_index(surface, i - 1)] +
                    surface.rbn * oldE[e_index(surface, i)] +
                    surface.g2a[1] * (haPlusValue - surface.Hyee_right) -
                    surface.g2b[1] * (hbPlusValue - hbMinusValue);
                i = -depth;
                surface.D[e_index(surface, i)] =
                    -surface.c1 * oldE[e_index(surface, i + 1)] +
                    surface.rb1 * oldE[e_index(surface, i)] +
                    surface.g2a[0] * (surface.Hyee_left - haMinusValue) -
                    surface.g2b[0] * (hbPlusValue - hbMinusValue);
            } else if (surface.correct_hb) {
                integer4 i = depth;
                surface.D[e_index(surface, i)] =
                    -surface.an * oldE[e_index(surface, i - 1)] +
                    surface.rbn * oldE[e_index(surface, i)] +
                    surface.g2a[1] * (haPlusValue - haMinusValue) -
                    surface.g2b[1] * (hbPlusValue - surface.Hyee_right);
                i = -depth;
                surface.D[e_index(surface, i)] =
                    -surface.c1 * oldE[e_index(surface, i + 1)] +
                    surface.rb1 * oldE[e_index(surface, i)] +
                    surface.g2a[0] * (haPlusValue - haMinusValue) -
                    surface.g2b[0] * (surface.Hyee_left - hbMinusValue);
            }

            for (integer4 i = -depth + 1; i <= depth - 1; ++i) {
                surface.D[e_index(surface, i)] =
                    -surface.a[e_index(surface, i)] *
                        oldE[e_index(surface, i - 1)] -
                    surface.c[e_index(surface, i)] *
                        oldE[e_index(surface, i + 1)] +
                    surface.rb[e_index(surface, i)] *
                        oldE[e_index(surface, i)] +
                    surface.rh[e_index(surface, i)] *
                        surface.H[h_index(surface, i)] -
                    surface.rhm1[e_index(surface, i)] *
                        surface.H[h_index(surface, i - 1)];
            }

            const size_t count = surface.E.size();
            solve_tridiag_prefactored_sgbc(
                surface.triA,
                surface.triCp,
                surface.triInvM,
                surface.D,
                newE,
                surface.triDp,
                static_cast<integer4>(count));
        } else {
            for (integer4 i = -depth + 1; i <= depth - 1; ++i) {
                surface.E[e_index(surface, i)] =
                    surface.G1_interno[e_index(surface, i)] *
                        surface.E[e_index(surface, i)] +
                    surface.G2_interno[e_index(surface, i)] *
                        (surface.H[h_index(surface, i)] -
                         surface.H[h_index(surface, i - 1)]);
            }
        }

        if (surface.SGBCCrank) {
            const auto& oldE = surface.E;
            const auto& newE = surface.E_past;
            for (integer4 i = -depth; i <= depth - 1; ++i) {
                surface.H[h_index(surface, i)] =
                    surface.GM1_interno[h_index(surface, i)] *
                        surface.H[h_index(surface, i)] +
                    surface.GM2_interno[h_index(surface, i)] / SGBCReal{2.0} *
                        (newE[e_index(surface, i + 1)] -
                         newE[e_index(surface, i)] +
                         oldE[e_index(surface, i + 1)] -
                         oldE[e_index(surface, i)]);
            }
            surface.Hyee_left = surface.H[h_index(surface, -depth)];
            surface.Hyee_right = surface.H[h_index(surface, depth - 1)];
            surface.Efield =
                (newE[e_index(surface, -depth)] +
                 newE[e_index(surface, depth)]) /
                SGBCReal{2.0};
            std::swap(surface.E, surface.E_past);
        } else if (depth != 0) {
            for (integer4 i = -depth; i <= depth - 1; ++i) {
                surface.H[h_index(surface, i)] =
                    surface.GM1_interno[h_index(surface, i)] *
                        surface.H[h_index(surface, i)] +
                    surface.GM2_interno[h_index(surface, i)] *
                        (surface.E[e_index(surface, i + 1)] -
                         surface.E[e_index(surface, i)]);
            }
            surface.Hyee_left = surface.H[h_index(surface, -depth)];
            surface.Hyee_right = surface.H[h_index(surface, depth - 1)];
            surface.Efield =
                (surface.E[e_index(surface, -depth)] +
                 surface.E[e_index(surface, depth)]) /
                SGBCReal{2.0};
        } else {
            surface.Efield = surface.E[e_index(surface, 0)];
        }
    }

    SGBCHCorrection_t AdvanceSGBCH(const SGBCSurface_t& surface,
                                   RKIND eField) {
        SGBCHCorrection_t correction;
        const integer4 depth = surface.depth;
        const SGBCReal eFieldValue = sgbc_real(eField);
        if (surface.correct_ha) {
            correction.ha_plus =
                surface.gm2_externo *
                (eFieldValue - surface.E[e_index(surface, depth)]);
            correction.ha_minus =
                -surface.gm2_externo *
                (eFieldValue - surface.E[e_index(surface, -depth)]);
        } else if (surface.correct_hb) {
            correction.hb_plus =
                -surface.gm2_externo *
                (eFieldValue - surface.E[e_index(surface, depth)]);
            correction.hb_minus =
                surface.gm2_externo *
                (eFieldValue - surface.E[e_index(surface, -depth)]);
        }
        return correction;
    }

    void g1g2_Dispersive(RKIND dt, RKIND epsilon, RKIND sigma,
                         RKIND& G1, RKIND& G2,
                         std::vector<CKIND>& Beta,
                         std::vector<CKIND>& Kappa,
                         std::vector<CKIND>& G3,
                         integer4 numpolres,
                         const std::vector<CKIND>& a11,
                         const std::vector<CKIND>& c11) {
        if (numpolres < 0) {
            throw std::invalid_argument("numpolres must be non-negative");
        }
        const auto count = static_cast<size_t>(numpolres);
        if (a11.size() < count || c11.size() < count) {
            throw std::invalid_argument("a11/c11 shorter than numpolres");
        }
        Beta.resize(count);
        Kappa.resize(count);
        G3.resize(count);
        for (size_t i = 0; i < count; ++i) {
            Kappa[i] = (1.0 + a11[i] * dt / 2.0) /
                       (1.0 - a11[i] * dt / 2.0);
            Beta[i] = (c11[i] * dt) /
                      (1.0 - a11[i] * dt / 2.0);
        }
        RKIND tempo = 0.0;
        for (size_t i = 0; i < count; ++i) {
            tempo += std::real(Beta[i]);
        }
        G1 = (2.0 * epsilon + tempo - sigma * dt) /
             (2.0 * epsilon + tempo + sigma * dt);
        G2 = 2.0 * dt / (2.0 * epsilon + tempo + sigma * dt);
        for (size_t i = 0; i < count; ++i) {
            G3[i] = G2 / 2.0 * (1.0 + Kappa[i]);
        }
    }

    SGBCLayerDepthResult calculate_sgbc_layer_depth(
        const std::vector<RKIND>& widths,
        const std::vector<RKIND>& sigmas,
        const std::vector<RKIND>& eprs,
        RKIND SGBCFreqValue,
        RKIND SGBCresolValue,
        integer4 SGBCdepthValue) {
        const size_t numcapas = widths.size();
        if (numcapas == 0 || sigmas.size() != numcapas || eprs.size() != numcapas) {
            throw std::invalid_argument("SGBC layer vectors must be non-empty and same-sized");
        }

        SGBCLayerDepthResult result;
        bool ultimacapamas1 = false;
        for (int precuenta = 0; precuenta <= 1; ++precuenta) {
            int celdafinal = 0;
            if (precuenta == 1) {
                if (result.depth % 2 != 0) {
                    result.depth += 1;
                    ultimacapamas1 = true;
                } else {
                    ultimacapamas1 = false;
                }
                result.depth = static_cast<integer4>(result.depth / 2.0);
                const size_t storageSize =
                    (result.depth > 0) ? static_cast<size_t>(2 * result.depth) : 1u;
                result.capa.assign(storageSize, 0);
                result.delta_entreEinterno.assign(storageSize, 0.0);
                celdafinal = -result.depth - 1;
            }

            for (size_t layer = 0; layer < numcapas; ++layer) {
                const SGBCReal width = sgbc_real(widths[layer]);
                const SGBCReal sigma = sgbc_real(sigmas[layer]);
                const SGBCReal sgbcFreq = sgbc_real(SGBCFreqValue);
                const SGBCReal sgbcResol = sgbc_real(SGBCresolValue);
                const SGBCReal piValue = sgbc_real(Pi);
                const SGBCReal mu0Value = sgbc_real(mu0);
                const SGBCReal epsilon =
                    sgbc_real(eprs[layer]) * sgbc_real(eps0);
                const SGBCReal skin_depth =
                    SGBCReal{1.0} /
                    (std::sqrt(SGBCReal{2.0}) * sgbcFreq * piValue *
                     std::pow(mu0Value * mu0Value *
                                  (SGBCReal{4.0} * epsilon * epsilon +
                                   sigma * sigma /
                                       (sgbcFreq * sgbcFreq * piValue * piValue)),
                              SGBCReal{0.25}) *
                     std::sin(std::atan2(SGBCReal{2.0} * piValue * epsilon * mu0Value,
                                         -(mu0Value * sigma) / sgbcFreq) /
                              SGBCReal{2.0}));

                integer4 anchocapa = 0;
                if (SGBCdepthValue == 0) {
                    if (numcapas > 1) {
                        throw std::invalid_argument(
                            "SGBCDepth=0 and numcapas>1 not compatible");
                    }
                    anchocapa = 1;
                } else if (SGBCdepthValue > 0) {
                    anchocapa = SGBCdepthValue;
                } else {
                    anchocapa = 1 + static_cast<integer4>(
                                        sgbcResol * width / skin_depth);
                }
                if (anchocapa < 2) anchocapa = 2;

                if (precuenta == 0) {
                    if (SGBCdepthValue == 0) {
                        result.depth = 0;
                    } else {
                        result.depth += anchocapa;
                    }
                } else {
                    if (SGBCdepthValue == 0) {
                        result.capa[0] = static_cast<integer4>(layer + 1);
                        result.delta_entreEinterno[0] = width;
                    } else {
                        const int celdainicial = celdafinal + 1;
                        celdafinal = celdainicial + anchocapa - 1;
                        if (layer + 1 == numcapas && ultimacapamas1) {
                            anchocapa += 1;
                            celdafinal += 1;
                        }
                        const SGBCReal delta =
                            width / static_cast<SGBCReal>(anchocapa);
                        for (int cell = celdainicial; cell <= celdafinal; ++cell) {
                            const size_t idx =
                                static_cast<size_t>(cell + result.depth);
                            result.capa[idx] = static_cast<integer4>(layer + 1);
                            result.delta_entreEinterno[idx] =
                                static_cast<RKIND>(delta);
                        }
                    }
                }
            }

            if (precuenta == 1 && result.depth != 0 &&
                celdafinal != result.depth - 1) {
                throw std::runtime_error("SGBC layer rounding mismatch");
            }
        }
        return result;
    }

    void depth(LegacySGBCSurface_t&, SGGFDTDINFO_t&, integer4, RKIND, RKIND, integer4) {}
    Malon_t* GetSGBCs() { return &malon; }
    void solve_tridiag_distintos(const std::vector<RKIND>& aa,
                                 const std::vector<RKIND>& bb,
                                 const std::vector<RKIND>& cc,
                                 RKIND a1, RKIND b1, RKIND c1,
                                 RKIND an, RKIND bn, RKIND cn,
                                 const std::vector<RKIND>& d,
                                 std::vector<RKIND>& x,
                                 integer4 n) {
        if (n <= 0) {
            x.clear();
            return;
        }
        const auto count = static_cast<size_t>(n);
        if (aa.size() < count || bb.size() < count || cc.size() < count ||
            d.size() < count) {
            throw std::invalid_argument("tridiagonal inputs shorter than n");
        }
        std::vector<RKIND> a(count), b(count), c(count), cp(count), dp(count);
        if (count == 1) {
            a[0] = an;
            b[0] = bn;
            c[0] = cn;
        } else {
            a[0] = a1;
            b[0] = b1;
            c[0] = c1;
            a[count - 1] = an;
            b[count - 1] = bn;
            c[count - 1] = cn;
        }
        for (size_t i = 1; i + 1 < count; ++i) {
            a[i] = aa[i];
            b[i] = bb[i];
            c[i] = cc[i];
        }
        cp[0] = c[0] / b[0];
        dp[0] = d[0] / b[0];
        for (size_t i = 1; i < count; ++i) {
            const RKIND m = b[i] - cp[i - 1] * a[i];
            cp[i] = c[i] / m;
            dp[i] = (d[i] - dp[i - 1] * a[i]) / m;
        }
        x.assign(count, 0.0);
        x[count - 1] = dp[count - 1];
        for (size_t i = count - 1; i-- > 0;) {
            x[i] = dp[i] - cp[i] * x[i + 1];
        }
    }

    void solve_tridiag_iguales(RKIND aa, RKIND bb, RKIND cc,
                               RKIND a1, RKIND b1, RKIND c1,
                               RKIND an, RKIND bn, RKIND cn,
                               const std::vector<RKIND>& d,
                               std::vector<RKIND>& x,
                               integer4 n) {
        if (n <= 0) {
            x.clear();
            return;
        }
        const auto count = static_cast<size_t>(n);
        std::vector<RKIND> a(count, aa), b(count, bb), c(count, cc);
        if (count == 1) {
            a[0] = an;
            b[0] = bn;
            c[0] = cn;
        } else {
            a[0] = a1;
            b[0] = b1;
            c[0] = c1;
            a[count - 1] = an;
            b[count - 1] = bn;
            c[count - 1] = cn;
        }
        solve_tridiag_distintos(a, b, c, a1, b1, c1, an, bn, cn, d, x, n);
    }
    void test_stab() {}

}
