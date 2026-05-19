#include <vector>
#include <cmath>
#include <string>

namespace DMMA_m {

    using RKIND = double;
    extern const int iEz;
    extern const int iEx;
    extern const int iEy;
    extern const double pi;
    extern const double EPS0;
    extern const double MU0;

    void dmma_thin_Slot(
        RKIND, RKIND, RKIND,
        const std::vector<RKIND>&, int, int,
        RKIND, RKIND, RKIND,
        std::vector<std::vector<RKIND>>&,
        std::vector<std::vector<RKIND>>&,
        RKIND, RKIND) {}

}
