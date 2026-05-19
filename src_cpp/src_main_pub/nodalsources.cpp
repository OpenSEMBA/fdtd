#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdint>

struct SGGFDTDINFO_t;
struct NodalSource_t;
struct XYZlimit_t;

void WarnErrReport(const std::string&) {}

constexpr int32_t RKIND = 8;
constexpr int32_t INTEGERSIZEOFMEDIAMATRICES = 4;
constexpr int BUFSIZE = 256;

enum Direction { iEx, iEy, iEz, iHx, iHy, iHz };

namespace nodalsources_m {

    struct xyzlimit_singlescaled_t {
        int32_t XI, XE, YI, YE, ZI, ZE;
        double amplitude;
    };

    struct NodalLocal_t {
        std::vector<double> evol;
        double deltaevol;
        int32_t numus;
        xyzlimit_singlescaled_t punto;
        bool IsInitialValue;
    };

    struct nodsou_t {
        int32_t NumHard = 0;
        int32_t NumSoft = 0;
        std::vector<NodalLocal_t> nodHard;
        std::vector<NodalLocal_t> nodSoft;
    };

    nodsou_t Nodal_Ex;
    nodsou_t Nodal_Ey;
    nodsou_t Nodal_Ez;
    nodsou_t Nodal_Hx;
    nodsou_t Nodal_Hy;
    nodsou_t Nodal_Hz;

    void InitnodalSources(const SGGFDTDINFO_t&, int32_t, int32_t, const std::vector<NodalSource_t>&, const std::vector<XYZlimit_t>&, bool&, bool&) {}
    double evolucion(double, const NodalLocal_t&) { return 0.0; }
    void DestroyNodal(SGGFDTDINFO_t&) {}
    void AdvancenodalE(const SGGFDTDINFO_t&, const std::vector<std::vector<std::vector<int32_t>>>&, const std::vector<std::vector<std::vector<int32_t>>>&, const std::vector<std::vector<std::vector<int32_t>>>&, int32_t, int32_t, const struct bounds_t&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, std::vector<std::vector<std::vector<double>>>&, std::vector<std::vector<std::vector<double>>>&, std::vector<std::vector<std::vector<double>>>&, bool) {}
    void AdvancenodalH(const SGGFDTDINFO_t&, const std::vector<std::vector<std::vector<int32_t>>>&, const std::vector<std::vector<std::vector<int32_t>>>&, const std::vector<std::vector<std::vector<int32_t>>>&, int32_t, int32_t, const struct bounds_t&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, std::vector<std::vector<std::vector<double>>>&, std::vector<std::vector<std::vector<double>>>&, std::vector<std::vector<std::vector<double>>>&, bool) {}

}
