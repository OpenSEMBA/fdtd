#ifndef LUMPED_SLIM_H
#define LUMPED_SLIM_H

#include <vector>

namespace lumped_slim {

struct LumpedMaterial_t {
    bool resistor = false;
    bool inductor = false;
    bool capacitor = false;
    bool diodo = false;
    double R = 0.0;
    double L = 0.0;
    double C = 0.0;
    double Rtime_on = 0.0;
    double Rtime_off = 1.0e30;
    int orient = 1;
    double epr = 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12;
    double sigma = 0.0;
};

struct LumpedNode_t {
    float* Efield = nullptr;
    float* Ha_Plus = nullptr;
    float* Ha_Minu = nullptr;
    float* Hb_Plus = nullptr;
    float* Hb_Minu = nullptr;

    double alignedDeltaE = 0.0;
    double transversalDeltaHa = 0.0;
    double transversalDeltaHb = 0.0;
    int orient = 1;

    LumpedMaterial_t mat;

    double EfieldPrevPrev = 0.0;
    double EfieldPrev = 0.0;
    double Jcur = 0.0;

    double G1 = 0.0;
    double G2a = 0.0;
    double G2b = 0.0;
    double GJ = 0.0;
    double sigmaEffResistInduct = 0.0;
    double currentCoeff = 0.0;
    double G1_usual = 0.0;
    double G2a_usual = 0.0;
    double G2b_usual = 0.0;
};

class LumpedSolver_t {
public:
    std::vector<LumpedNode_t> nodes;

    void clear() { nodes.clear(); }

    void calcConstants(double dt, double eps0, double mu0);

    void advance(int timestep, double dt);
};

} // namespace lumped_slim

#endif
