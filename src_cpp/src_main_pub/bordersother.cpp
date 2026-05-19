#include <vector>
#include <cstdint>

struct SGGFDTDINFO_t;
struct logic_control_t;
struct XYZlimit_t;
struct Border_t;

#ifndef RKIND
#define RKIND double
#endif

namespace BORDERS_other_m {

    void InitOtherBorders(const SGGFDTDINFO_t&, logic_control_t&) {}
    void MinusCloneMagneticPMC(
        const std::vector<XYZlimit_t>&,
        std::vector<RKIND>&,
        std::vector<RKIND>&,
        std::vector<RKIND>&,
        std::vector<XYZlimit_t>&,
        int, int, const Border_t&) {}
    void CloneMagneticPeriodic(
        const std::vector<XYZlimit_t>&,
        std::vector<std::vector<std::vector<RKIND>>>&,
        std::vector<std::vector<std::vector<RKIND>>>&,
        std::vector<std::vector<std::vector<RKIND>>>&,
        std::vector<XYZlimit_t>&,
        int32_t, int32_t, const Border_t&) {}

}
