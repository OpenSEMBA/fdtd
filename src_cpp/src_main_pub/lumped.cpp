#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <limits>
#include <memory>

namespace Report_m {
    void WarnErrReport(const std::string&, bool);
    void print11(int, const std::string&);
    const std::string SEPARADOR = "========================================";
}

namespace FDETYPES_m {
    using RKIND = double;
    const int BUFSIZE = 256;
}

struct AllocInfo { int XI, XE, YI, YE, ZI, ZE; };
struct MediaIs { bool Lumped = false; bool PML = false; };
struct LumpedParam { int Orient = 0; bool inductor = false; bool diodo = false; bool resistor = false; bool capacitor = false; double R = 0.0; double L = 0.0; double C = 0.0; double DiodB = 0.0; double DiodIsat = 0.0; double Rtime_on = 0.0; double Rtime_off = 0.0; };
struct MediaStruct { MediaIs Is; std::vector<LumpedParam> Lumped; double epr = 1.0; double sigma = 0.0; bool sigmareasignado = false; };
struct SGGFDTDINFO_t { std::vector<AllocInfo> alloc; std::vector<AllocInfo> SINPMLSweep; std::vector<int> sggMiEx; std::vector<int> sggMiEy; std::vector<int> sggMiEz; int NumMedia = 0; std::vector<MediaStruct> Med; double dt = 0.0; };
struct media_matrices_t {};
struct sim_control_t { int layoutnumber = 0; int num_procs = 1; bool resume = false; bool stochastic = false; };

enum IndexType { iEx, iEy, iEz, iHx, iHy, iHz };

struct Nodes_t {
    double alignedDeltaE = 0.0;
    double transversalDeltaHa = 0.0;
    double transversalDeltaHb = 0.0;
    int Orient = 0;
    int jmed = 0;
    double* Efield = nullptr;
    double* Ha_Plus = nullptr;
    double* Ha_Minu = nullptr;
    double* Hb_Plus = nullptr;
    double* Hb_Minu = nullptr;
    double EfieldPrevPrev = 0.0;
    double EfieldPrev = 0.0;
    double Jcur = 0.0;
    double EfieldPrevPrev_for_devia = 0.0;
    double EfieldPrev_for_devia = 0.0;
    double Jcur_for_devia = 0.0;
    double G1 = 0.0;
    double G2a = 0.0;
    double G2b = 0.0;
    double GJ = 0.0;
    double sigmaEffResistInduct = 0.0;
    double currentCoeff = 0.0;
    double G1_usual = 0.0;
    double G2a_usual = 0.0;
    double G2b_usual = 0.0;
    double diodeB = 0.0;
    double diodepreA = 0.0;
};

struct LumpedElem_t {
    int NumNodes = 0;
    std::vector<Nodes_t> Nodes;
};

LumpedElem_t LumpElem;
double eps0 = 0.0;
double mu0 = 0.0;
double zvac = 0.0;
double cluz = 0.0;

void inject_devialumped(const SGGFDTDINFO_t&, int, bool, bool, Nodes_t&);
void calc_lumped_deviaconsts(const SGGFDTDINFO_t&, Nodes_t&, int, double, double, double, double, double, double, double);

namespace Lumped_m {

    void InitLumped(const SGGFDTDINFO_t&, const media_matrices_t&,
                    const std::vector<std::vector<std::vector<double>>>&,
                    const std::vector<std::vector<std::vector<double>>>&,
                    const std::vector<std::vector<std::vector<double>>>&,
                    const std::vector<std::vector<std::vector<double>>>&,
                    const std::vector<std::vector<std::vector<double>>>&,
                    const std::vector<std::vector<std::vector<double>>>&,
                    const std::vector<int>&, const std::vector<int>&, const std::vector<int>&,
                    const std::vector<int>&, const std::vector<int>&, const std::vector<int>&,
                    bool&, double, double, const sim_control_t&) {}
    void AdvanceLumpedE(const SGGFDTDINFO_t&, int, bool, bool) {}
    void calc_lumpedconstants(const SGGFDTDINFO_t&, double, double) {}
    void StoreFieldsLumpeds(bool) {}
    void DestroyLumped(SGGFDTDINFO_t&) {}
    double newton_raphson(double, double, double) { return 0.0; }
    LumpedElem_t* Getlumped() { return &LumpElem; }

}
