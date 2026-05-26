#ifndef SEMBA_FDTD_MALONEY_NOSTOCH_H
#define SEMBA_FDTD_MALONEY_NOSTOCH_H

#include <complex>
#include <cstdint>
#include <vector>

namespace SGBC_nostoch_m {

struct SGBCLayerDepthResult {
    int32_t depth = 0;
    std::vector<int32_t> capa;
    std::vector<double> delta_entreEinterno;
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
