#include <vector>
#include <memory>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <array>

struct Thinwires_t {
    std::vector<int> WireTipoMedio;
    int NumDifferentWires = 0;
    struct ChargeNodes_t {
        int indexnode = -1;
        double ChargePresent = 0.0;
        double ChargePast = 0.0;
        bool IsAttachedtoVoltage = false;
        bool IsMur = false;
        bool IsBackDownLeftMur = false;
        bool IsFrontUpRightMur = false;
        bool IsPeriodic = false;
        bool IsPEC = false;
        bool IsLossy = false;
        bool HasIsource = false;
        bool IsHeterogeneousJunction = false;
        bool Exists = false;
        bool proc = false;
        bool Is_LeftEnd = false;
        bool Is_RightEnd = false;
        bool IsInSingleRLCsegment = false;
        int NumCurrentMinus = 0;
        int NumCurrentPlus = 0;
        int i = -1, j = -1, k = -1;
        double cteMur = 0.0, cteProp = 0.0, oRIGctePlain = 0.0, ctePlain = 0.0;
    } NullNode;
    struct CurrentSegments_t {
        double R = 0.0, Resist = 0.0, Resist_devia = 0.0;
        double C = 0.0, L = 0.0, Lintrinsic = 0.0;
        int NumParallel = 1, origindex = 0, indexsegment = 0;
        double currentpast = 0.0, current = 0.0, qplus_qminus = 0.0;
        double current_for_devia = 0.0, qplus_qminus_for_devia = 0.0;
        double Efield_main2wire_for_devia = 0.0;
        double inv_Lind_acum = 0.0, Lind_acum = 0.0, Lind = 0.0, Lind_devia = 0.0;
        double HEUR_safety = 0.0, delta = 0.0;
        double deltaTransv1 = 0.0, deltaTransv2 = 0.0;
        double cte1 = 0.0, cte2 = 0.0, cte3 = 0.0;
        double cte1_for_devia = 0.0, cte2_for_devia = 0.0, cte3_for_devia = 0.0;
        double cte5 = 0.0;
        int ilibre = -1, jlibre = -1, klibre = -1;
        int i = -1, j = -1, k = -1, tipofield = -1;
        bool IsPMC = false, orientadoalreves = false;
        bool HasVsource = false, IsShielded = false;
        bool HasAbsorbing_RightEnd = false, HasAbsorbing_LeftEnd = false;
        bool HasParallel_RightEnd = false, HasParallel_LeftEnd = false;
        bool HasSeries_RightEnd = false, HasSeries_LeftEnd = false;
        bool IsEnd_norLeft_norRight = false;
        bool Is_LeftEnd = false, Is_RightEnd = false;
        ChargeNodes_t* chargePlus = nullptr;
        ChargeNodes_t* chargeMinus = nullptr;
    } NullSegment;
    double olddt = 0.0;
    int NumNeededCurrentUpMPI = 0;
    int NumNeededCurrentDownMPI = 0;
};

struct limit_t { int XI, XE, YI, YE, ZI, ZE; };

struct SGGFDTDINFO_t {
    double dt = 0.0;
    int NumMedia = 0;
    struct { int XI, XE, YI, YE, ZI, ZE; } Alloc[10];
    struct { int ZI, ZE; } Sweep[10];
    struct {
        double Epr = 1.0, Mur = 1.0;
        struct { bool Is = {false}; struct { bool ThinWire = false; int numsegmentos = 0; struct { bool multirabo = false; int i, j, k, ori, origIndex; bool Is_LeftEnd = false, Is_RightEnd = false, IsEnd_norLeft_norRight = false; int ilibre = -1, jlibre = -1, klibre = -1; } segm[100]; struct { bool disp = false, disp_RightEnd = false, disp_LeftEnd = false; } wire[1]; bool HasAbsorbing_RightEnd = false, HasAbsorbing_LeftEnd = false, HasParallel_RightEnd = false, HasParallel_LeftEnd = false, HasSeries_RightEnd = false, HasSeries_LeftEnd = false; } wire[1]; } Med[100];
    } sgg;
};

struct sim_control_t {
    int layoutnumber = 0;
    int num_procs = 1;
    bool strictOLD = false;
    bool connectendings = false;
};

namespace HollandWires_m {

    const double HEUR_RADIUSOVERDELTA = 10.0;
    bool thereAreVsources = false;
    bool thereAreIsources = false;
    bool thereAreMurConditions = false;
    Thinwires_t HWires;
    std::vector<double> InvEps;
    std::vector<double> InvMu;
    std::vector<double> OldInvEps;
    std::vector<double> OldInvMu;
    double eps0 = 0.0;
    double mu0 = 0.0;

    void WarnErrReport(const std::string&, bool = false) {}

    void InitWires(
        SGGFDTDINFO_t&,
        const std::vector<int>&, std::vector<int>&, std::vector<int>&, std::vector<int>&,
        std::vector<int>&, std::vector<int>&, std::vector<int>&,
        std::vector<int>&, std::vector<int>&, std::vector<int>&,
        std::vector<int>&, std::vector<int>&, std::vector<int>&,
        const std::vector<double>&,
        const std::vector<limit_t>&, const std::vector<limit_t>&,
        double&, double, double, const sim_control_t&) {}
    void AdvanceWiresE() {}
    void AdvanceWiresH() {}
    void AdvanceWiresEcrank() {}
    void StoreFieldsWires() {}
    void DestroyWires() {}
    void GetHwires() {}
    void ReportWireJunctions() {}
    void calc_wirehollandconstants() {}

}
