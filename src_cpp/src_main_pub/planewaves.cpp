#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdexcept>

namespace FDETYPES_m {
    constexpr int RKIND = 8;
    constexpr int BUFSIZE = 256;
}

namespace Report_m {
    void stoponerror(int, int, const std::string&) {}
    void WarnErrReport(const std::string&) {}
    struct coorsxyzP_t {
        struct PhysCoor_t {
            std::vector<double> x;
            std::vector<double> y;
            std::vector<double> z;
        };
        std::vector<PhysCoor_t> PhysCoor;
    };
}

struct SGGFDTDINFO_t;
struct media_matrices_t;
struct limit_t;

enum FieldIndex { iEx = 0, iEy = 1, iEz = 2, iHx = 3, iHy = 4, iHz = 5 };

namespace ilumina_m {

    double cluz = 0.0;
    double zvac = 0.0;
    double eps0 = 0.0;
    double mu0 = 0.0;
    Report_m::coorsxyzP_t Punto;

    struct ehxyz_t {
        int Ex = -15, Ey = -15, Ez = -15, Hx = -15, Hy = -15, Hz = -15;
    };

    struct tfidaa_t {
        ehxyz_t com, fin, tra, fro, izq, der, aba, arr;
    };

    struct ijk_t {
        tfidaa_t i, j, k;
    };

    std::vector<ijk_t> TrFr;
    std::vector<ijk_t> IzDe;
    std::vector<ijk_t> AbAr;
    std::vector<bool> IluminaTr;
    std::vector<bool> IluminaFr;
    std::vector<bool> IluminaIz;
    std::vector<bool> IluminaDe;
    std::vector<bool> IluminaAr;
    std::vector<bool> IluminaAb;
    std::vector<int> numus;
    std::vector<double> deltaevol;
    std::vector<std::vector<std::vector<double>>> fpw;
    std::vector<std::vector<double>> distanciaInicial;
    std::vector<std::vector<double>> pxpw;
    std::vector<std::vector<double>> pypw;
    std::vector<std::vector<double>> pzpw;
    std::vector<std::vector<double>> INCERT;
    std::vector<std::vector<double>> evol;

    void InitPlaneWave(const SGGFDTDINFO_t&, const media_matrices_t&, int, int, const std::vector<limit_t>&, bool&, bool, double, double) {}
    void AdvancePlaneWave(const SGGFDTDINFO_t&, int, int, int, std::vector<std::vector<std::vector<double>>>&, std::vector<std::vector<std::vector<double>>>&, std::vector<std::vector<std::vector<double>>>&, std::vector<std::vector<std::vector<double>>>&, std::vector<std::vector<std::vector<double>>>&, std::vector<std::vector<std::vector<double>>>&, std::vector<std::vector<std::vector<double>>>&, std::vector<std::vector<std::vector<double>>>&, std::vector<std::vector<std::vector<double>>>&, std::vector<std::vector<std::vector<double>>>&, std::vector<std::vector<std::vector<double>>>&, std::vector<std::vector<std::vector<double>>>&, bool) {}
    void StoreFieldsPlaneWave() {}
    void DestroyPlaneWave(SGGFDTDINFO_t&) {}
    void calc_planewaveconstants(const SGGFDTDINFO_t&, double, double) {}
    void DestroyIlumina(SGGFDTDINFO_t&) {}

}
