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
    void AdvanceSGBCE(RKIND, logical, logical, logical) {}
    void AdvanceSGBCH() {}
    void StoreFieldsSGBCs(logical) {}
    void DestroySGBCs(SGGFDTDINFO_t&) {}
    void calc_g1g2gm1gm2_compo(SGGFDTDINFO_t&, SGBCSurface_t&, RKIND, RKIND, logical) {}
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
                const RKIND width = widths[layer];
                const RKIND sigma = sigmas[layer];
                const RKIND epsilon = eprs[layer] * eps0;
                const RKIND skin_depth =
                    1.0 /
                    (std::sqrt(2.0) * SGBCFreqValue * Pi *
                     std::pow(mu0 * mu0 *
                                  (4.0 * epsilon * epsilon +
                                   sigma * sigma /
                                       (SGBCFreqValue * SGBCFreqValue * Pi * Pi)),
                              0.25) *
                     std::sin(std::atan2(2.0 * Pi * epsilon * mu0,
                                         -(mu0 * sigma) / SGBCFreqValue) /
                              2.0));

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
                                        SGBCresolValue * width / skin_depth);
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
                        const RKIND delta = width / anchocapa;
                        for (int cell = celdainicial; cell <= celdafinal; ++cell) {
                            const size_t idx =
                                static_cast<size_t>(cell + result.depth);
                            result.capa[idx] = static_cast<integer4>(layer + 1);
                            result.delta_entreEinterno[idx] = delta;
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

    void depth(SGBCSurface_t&, SGGFDTDINFO_t&, integer4, RKIND, RKIND, integer4) {}
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
