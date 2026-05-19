#include <string>
#include <vector>
#include <iostream>
struct SGGFDTDINFO_t { int NumberRequest=0; double dt=0.0; struct { bool Volumic=false; int nP=0; int What[10]={}; int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; bool done=false,flushed=false,Begun=false,TimeDomain=false,FreqDomain=false; } observation[100]; struct { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; } Alloc[10]; };
struct sim_control_t { std::string wiresflavor=""; int layoutnumber=0; int num_procs=1; bool resume=false; bool stochastic=false; };
struct media_matrices_t { int NumMed=0; std::vector<int> sggMiNo, sggMiEx, sggMiEy, sggMiEz, sggMiHx, sggMiHy, sggMiHz; };
struct limit_t { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; };
struct thereAre_t { bool Wires=false; };
struct g_t { int g2=0; };
struct sinPMLSweep_t { std::vector<int> sggMiNo; };
class timestepping_t {
public:
    SGGFDTDINFO_t sgg;
    sim_control_t control;
    media_matrices_t media;
    limit_t lim;
    thereAre_t thereAre;
    g_t g;
    sinPMLSweep_t sinPML_fullsize;
    double dtcritico=0.0, eps0=0.0, mu0=0.0;
    double *Ex=nullptr, *Ey=nullptr, *Ez=nullptr;
    double *Hx=nullptr, *Hy=nullptr, *Hz=nullptr;
    int *Idxe=nullptr, *Idye=nullptr, *Idze=nullptr;
    int *Idxh=nullptr, *Idyh=nullptr, *Idzh=nullptr;
    int fullsize=0;
    timestepping_t() {}
    ~timestepping_t() {
        delete[] Ex; delete[] Ey; delete[] Ez;
        delete[] Hx; delete[] Hy; delete[] Hz;
        delete[] Idxe; delete[] Idye; delete[] Idze;
        delete[] Idxh; delete[] Idyh; delete[] Idzh;
    }
    void advanceE() {}
    void advanceH() {}
    void applyBCs() {}
    void sampleProbes() {}
    void writeOutput() {}
    void run(int maxSteps) {
        for(int i=0; i<maxSteps; i++) {
            advanceE();
            advanceH();
            applyBCs();
            sampleProbes();
        }
    }
};
