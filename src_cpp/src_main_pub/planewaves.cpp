#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <sstream>

// Forward declarations and includes based on assumed dependencies
// Note: In a real scenario, these headers would be provided by the project structure.
// Assuming FDETYPES_m and Report_m provide the following types and functions:

// Mocking external types and functions for compilation context
// These should be replaced with actual headers in the real project

struct coorsxyzP_t {
    struct PhysCoor_t {
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
    };
    std::vector<PhysCoor_t> PhysCoor;
};

struct limit_t {
    int XI, XE, YI, YE, ZI, ZE;
};

struct media_matrices_t {
    // Assuming 3D arrays for material indices, flattened or vector of vectors
    std::vector<std::vector<std::vector<int>>> sggMiEx;
    std::vector<std::vector<std::vector<int>>> sggMiEy;
    std::vector<std::vector<std::vector<int>>> sggMiEz;
    std::vector<std::vector<std::vector<int>>> sggMiHx;
    std::vector<std::vector<std::vector<int>>> sggMiHy;
    std::vector<std::vector<std::vector<int>>> sggMiHz;
};

struct media_t {
    bool is_pec;
};

struct SGGFDTDINFO_t {
    struct Sweep_t {
        int XI, XE, YI, YE, ZI, ZE;
    };
    std::vector<Sweep_t> Sweep;
    std::vector<double> LineX;
    std::vector<double> Liney;
    std::vector<double> LineZ;
    int NumPlaneWaves;
    struct {
        int nummodes;
        std::vector<double> px, py, pz, ex, ey, ez, incert;
        bool isRC;
        int esqx1, esqx2, esqy1, esqy2, esqz1, esqz2;
        struct {
            std::string Name;
            int NumSamples;
            std::vector<double> Samples;
            double deltaSamples;
        } Fichero;
    } PlaneWave;
    std::vector<limit_t> SINPMLSweep;
    std::vector<media_t> med;
    double dt;
};

// External functions assumed to exist
void stoponerror(int layoutnumber, int num_procs, const std::string& msg);
void WarnErrReport(const std::string& msg);

// Constants
const int iEx = 0;
const int iEy = 1;
const int iEz = 2;
const int iHx = 3;
const int iHy = 4;
const int iHz = 5;

namespace ilumina_m {

    using RKIND = double;

    // Global variables
    std::vector<std::vector<std::vector<RKIND>>> fpw;
    std::vector<std::vector<RKIND>> distanciaInicial;
    std::vector<std::vector<RKIND>> pxpw;
    std::vector<std::vector<RKIND>> pypw;
    std::vector<std::vector<RKIND>> pzpw;
    std::vector<std::vector<RKIND>> INCERT;
    std::vector<std::vector<RKIND>> evol;
    std::vector<RKIND> deltaevol;
    std::vector<int> numus;

    struct ehxyz_t {
        int Ex = -15;
        int Ey = -15;
        int Ez = -15;
        int Hx = -15;
        int Hy = -15;
        int Hz = -15;
    };

    struct tfidaa_t {
        ehxyz_t com, fin, tra, fro, izq, der, aba, arr;
    };

    struct ijk_t {
        tfidaa_t i, j, k;
    };

    RKIND cluz = 0.0;
    RKIND zvac = 0.0;
    RKIND eps0 = 0.0;
    RKIND mu0 = 0.0;

    // Local variables
    coorsxyzP_t Punto;
    std::vector<ijk_t> TrFr;
    std::vector<ijk_t> IzDe;
    std::vector<ijk_t> AbAr;
    std::vector<bool> IluminaTr;
    std::vector<bool> IluminaFr;
    std::vector<bool> IluminaIz;
    std::vector<bool> IluminaDe;
    std::vector<bool> IluminaAr;
    std::vector<bool> IluminaAb;

    const int BUFSIZE = 256; // Assumed constant

    void Incid() {
        // Placeholder as body not provided in snippet
    }

    void AdvancePlaneWaveE() {
        // Placeholder as body not provided in snippet
    }

    void AdvancePlaneWaveH() {
        // Placeholder as body not provided in snippet
    }

    void DestroyIlumina() {
        // Placeholder as body not provided in snippet
    }

    void storeplanewaves() {
        // Placeholder as body not provided in snippet
    }

    void calc_planewaveconstants() {
        // Placeholder as body not provided in snippet
    }

    void corrigeondaplanaH() {
        // Placeholder as body not provided in snippet
    }

    void InitPlaneWave(const SGGFDTDINFO_t& sgg, const media_matrices_t& media, int layoutnumber, int num_procs, const std::vector<limit_t>& SINPML_fullsize, bool& ThereArePlaneWaveBoxes, bool resume, RKIND eps00, RKIND mu00) {
        eps0 = eps00;
        mu0 = mu00;
        cluz = 1.0 / std::sqrt(eps0 * mu0);
        zvac = std::sqrt(mu0 / eps0);

        // Allocate Punto coordinates
        for (int field = iEx; field <= iHz; ++field) {
            if (field < sgg.Sweep.size()) {
                Punto.PhysCoor[field].x.resize(sgg.Sweep[field].XE - sgg.Sweep[field].XI + 3);
                Punto.PhysCoor[field].y.resize(sgg.Sweep[field].YE - sgg.Sweep[field].YI + 3);
                Punto.PhysCoor[field].z.resize(sgg.Sweep[field].ZE - sgg.Sweep[field].ZI + 3);
            }
        }

        // Fill coordinates for each field component
        // iEx
        {
            int field = iEx;
            if (field < sgg.Sweep.size()) {
                for (int i = sgg.Sweep[field].XI - 1; i <= sgg.Sweep[field].XE; ++i) {
                    if (i >= 0 && i < (int)Punto.PhysCoor[field].x.size()) {
                        Punto.PhysCoor[field].x[i] = (sgg.LineX[i] + sgg.LineX[i + 1]) * 0.5;
                    }
                }
                for (int j = sgg.Sweep[field].YI - 1; j <= sgg.Sweep[field].YE + 1; ++j) {
                    if (j >= 0 && j < (int)Punto.PhysCoor[field].y.size()) {
                        Punto.PhysCoor[field].y[j] = sgg.Liney[j];
                    }
                }
                for (int k = sgg.Sweep[field].ZI - 1; k <= sgg.Sweep[field].ZE + 1; ++k) {
                    if (k >= 0 && k < (int)Punto.PhysCoor[field].z.size()) {
                        Punto.PhysCoor[field].z[k] = sgg.LineZ[k];
                    }
                }
            }
        }
        // iEy
        {
            int field = iEy;
            if (field < sgg.Sweep.size()) {
                for (int i = sgg.Sweep[field].XI - 1; i <= sgg.Sweep[field].XE + 1; ++i) {
                    if (i >= 0 && i < (int)Punto.PhysCoor[field].x.size()) {
                        Punto.PhysCoor[field].x[i] = sgg.LineX[i];
                    }
                }
                for (int j = sgg.Sweep[field].YI - 1; j <= sgg.Sweep[field].YE; ++j) {
                    if (j >= 0 && j < (int)Punto.PhysCoor[field].y.size()) {
                        Punto.PhysCoor[field].y[j] = (sgg.Liney[j] + sgg.Liney[j + 1]) * 0.5;
                    }
                }
                for (int k = sgg.Sweep[field].ZI - 1; k <= sgg.Sweep[field].ZE + 1; ++k) {
                    if (k >= 0 && k < (int)Punto.PhysCoor[field].z.size()) {
                        Punto.PhysCoor[field].z[k] = sgg.LineZ[k];
                    }
                }
            }
        }
        // iEz
        {
            int field = iEz;
            if (field < sgg.Sweep.size()) {
                for (int i = sgg.Sweep[field].XI - 1; i <= sgg.Sweep[field].XE + 1; ++i) {
                    if (i >= 0 && i < (int)Punto.PhysCoor[field].x.size()) {
                        Punto.PhysCoor[field].x[i] = sgg.LineX[i];
                    }
                }
                for (int j = sgg.Sweep[field].YI - 1; j <= sgg.Sweep[field].YE + 1; ++j) {
                    if (j >= 0 && j < (int)Punto.PhysCoor[field].y.size()) {
                        Punto.PhysCoor[field].y[j] = sgg.Liney[j];
                    }
                }
                for (int k = sgg.Sweep[field].ZI - 1; k <= sgg.Sweep[field].ZE; ++k) {
                    if (k >= 0 && k < (int)Punto.PhysCoor[field].z.size()) {
                        Punto.PhysCoor[field].z[k] = (sgg.LineZ[k] + sgg.LineZ[k + 1]) * 0.5;
                    }
                }
            }
        }
        // iHx
        {
            int field = iHx;
            if (field < sgg.Sweep.size()) {
                for (int i = sgg.Sweep[field].XI - 1; i <= sgg.Sweep[field].XE + 1; ++i) {
                    if (i >= 0 && i < (int)Punto.PhysCoor[field].x.size()) {
                        Punto.PhysCoor[field].x[i] = sgg.LineX[i];
                    }
                }
                for (int j = sgg.Sweep[field].YI - 1; j <= sgg.Sweep[field].YE; ++j) {
                    if (j >= 0 && j < (int)Punto.PhysCoor[field].y.size()) {
                        Punto.PhysCoor[field].y[j] = (sgg.Liney[j] + sgg.Liney[j + 1]) * 0.5;
                    }
                }
                for (int k = sgg.Sweep[field].ZI - 1; k <= sgg.Sweep[field].ZE; ++k) {
                    if (k >= 0 && k < (int)Punto.PhysCoor[field].z.size()) {
                        Punto.PhysCoor[field].z[k] = (sgg.LineZ[k] + sgg.LineZ[k + 1]) * 0.5;
                    }
                }
            }
        }
        // iHy
        {
            int field = iHy;
            if (field < sgg.Sweep.size()) {
                for (int i = sgg.Sweep[field].XI - 1; i <= sgg.Sweep[field].XE; ++i) {
                    if (i >= 0 && i < (int)Punto.PhysCoor[field].x.size()) {
                        Punto.PhysCoor[field].x[i] = (sgg.LineX[i] + sgg.LineX[i + 1]) * 0.5;
                    }
                }
                for (int j = sgg.Sweep[field].YI - 1; j <= sgg.Sweep[field].YE + 1; ++j) {
                    if (j >= 0 && j < (int)Punto.PhysCoor[field].y.size()) {
                        Punto.PhysCoor[field].y[j] = sgg.Liney[j];
                    }
                }
                for (int k = sgg.Sweep[field].ZI - 1; k <= sgg.Sweep[field].ZE; ++k) {
                    if (k >= 0 && k < (int)Punto.PhysCoor[field].z.size()) {
                        Punto.PhysCoor[field].z[k] = (sgg.LineZ[k] + sgg.LineZ[k + 1]) * 0.5;
                    }
                }
            }
        }
        // iHz
        {
            int field = iHz;
            if (field < sgg.Sweep.size()) {
                for (int i = sgg.Sweep[field].XI - 1; i <= sgg.Sweep[field].XE; ++i) {
                    if (i >= 0 && i < (int)Punto.PhysCoor[field].x.size()) {
                        Punto.PhysCoor[field].x[i] = (sgg.LineX[i] + sgg.LineX[i + 1]) * 0.5;
                    }
                }
                for (int j = sgg.Sweep[field].YI - 1; j <= sgg.Sweep[field].YE; ++j) {
                    if (j >= 0 && j < (int)Punto.PhysCoor[field].y.size()) {
                        Punto.PhysCoor[field].y[j] = (sgg.Liney[j] + sgg.Liney[j + 1]) * 0.5;
                    }
                }
                for (int k = sgg.Sweep[field].ZI - 1; k <= sgg.Sweep[field].ZE + 1; ++k) {
                    if (k >= 0 && k < (int)Punto.PhysCoor[field].z.size()) {
                        Punto.PhysCoor[field].z[k] = sgg.LineZ[k];
                    }
                }
            }
        }

        if (sgg.NumPlaneWaves >= 1) {
            ThereArePlaneWaveBoxes = true;
        } else {
            ThereArePlaneWaveBoxes = false;
            return;
        }

        if (ThereArePlaneWaveBoxes) {
            ThereArePlaneWaveBoxes = false; // Reset because MPI slice might not have any
            int n_pw = sgg.NumPlaneWaves;
            TrFr.resize(n_pw);
            IzDe.resize(n_pw);
            AbAr.resize(n_pw);
            IluminaTr.resize(n_pw, false);
            IluminaFr.resize(n_pw, false);
            IluminaIz.resize(n_pw, false);
            IluminaDe.resize(n_pw, false);
            IluminaAr.resize(n_pw, false);
            IluminaAb.resize(n_pw, false);
            numus.resize(n_pw);
            deltaevol.resize(n_pw);

            for (int jjj = 0; jjj < n_pw; ++jjj) {
                numus[jjj] = sgg.PlaneWave.Fichero.NumSamples;
                
                bool abortar = 
                    (sgg.PlaneWave.esqx1 <= SINPML_fullsize[iHx].XI) &&
                    (sgg.PlaneWave.esqx2 >= SINPML_fullsize[iHx].XE) &&
                    (sgg.PlaneWave.esqy1 <= SINPML_fullsize[iHy].YI) &&
                    (sgg.PlaneWave.esqy2 >= SINPML_fullsize[iHy].YE) &&
                    (sgg.PlaneWave.esqz1 <= SINPML_fullsize[iHz].ZI) &&
                    (sgg.PlaneWave.esqz2 >= SINPML_fullsize[iHz].ZE);
                
                if (abortar) {
                    std::string buff = "At least one of TF/SF planes must be 1 cell inside the simulation region. Aborting";
                    stoponerror(layoutnumber, num_procs, buff);
                }

                IluminaTr[jjj] = false;
                IluminaFr[jjj] = false;
                IluminaIz[jjj] = false;
                IluminaDe[jjj] = false;
                IluminaAr[jjj] = false;
                IluminaAb[jjj] = false;

                if ((sgg.PlaneWave.esqx1 >= sgg.SINPMLSweep[iHx].XI) && (sgg.PlaneWave.esqx1 <= sgg.SINPMLSweep[iHx].XE))
                    IluminaTr[jjj] = true;
                if ((sgg.PlaneWave.esqx2 <= sgg.SINPMLSweep[iHx].XE) && (sgg.PlaneWave.esqx2 >= sgg.SINPMLSweep[iHx].XI))
                    IluminaFr[jjj] = true;
                if ((sgg.PlaneWave.esqy1 >= sgg.SINPMLSweep[iHy].YI) && (sgg.PlaneWave.esqy1 <= sgg.SINPMLSweep[iHy].YE))
                    IluminaIz[jjj] = true;
                if ((sgg.PlaneWave.esqy2 <= sgg.SINPMLSweep[iHy].YE) && (sgg.PlaneWave.esqy2 >= sgg.SINPMLSweep[iHy].YI))
                    IluminaDe[jjj] = true;
                if ((sgg.PlaneWave.esqz1 >= sgg.SINPMLSweep[iHz].ZI) && (sgg.PlaneWave.esqz1 <= sgg.SINPMLSweep[iHz].ZE))
                    IluminaAb[jjj] = true;
                if ((sgg.PlaneWave.esqz2 <= sgg.SINPMLSweep[iHz].ZE) && (sgg.PlaneWave.esqz2 >= sgg.SINPMLSweep[iHz].ZI))
                    IluminaAr[jjj] = true;

                // Find coordinate limits
                // TrFr
                TrFr[jjj].i.tra.Ez = std::max(sgg.SINPMLSweep[iEz].XI, sgg.PlaneWave.esqx1);
                TrFr[jjj].i.fro.Ez = std::min(sgg.SINPMLSweep[iEz].XE, sgg.PlaneWave.esqx2);
                TrFr[jjj].j.com.Ez = std::max(sgg.SINPMLSweep[iEz].YI, sgg.PlaneWave.esqy1);
                TrFr[jjj].j.fin.Ez = std::min(sgg.SINPMLSweep[iEz].YE, sgg.PlaneWave.esqy2);
                TrFr[jjj].k.com.Ez = std::max(sgg.SINPMLSweep[iEz].ZI, sgg.PlaneWave.esqz1);
                TrFr[jjj].k.fin.Ez = std::min(sgg.SINPMLSweep[iEz].ZE, sgg.PlaneWave.esqz2 - 1);

                TrFr[jjj].i.tra.Ey = std::max(sgg.SINPMLSweep[iEy].XI, sgg.PlaneWave.esqx1);
                TrFr[jjj].i.fro.Ey = std::min(sgg.SINPMLSweep[iEy].XE, sgg.PlaneWave.esqx2);
                TrFr[jjj].j.com.Ey = std::max(sgg.SINPMLSweep[iEy].YI, sgg.PlaneWave.esqy1);
                TrFr[jjj].j.fin.Ey = std::min(sgg.SINPMLSweep[iEy].YE, sgg.PlaneWave.esqy2 - 1);
                TrFr[jjj].k.com.Ey = std::max(sgg.SINPMLSweep[iEy].ZI, sgg.PlaneWave.esqz1);
                TrFr[jjj].k.fin.Ey = std::min(sgg.SINPMLSweep[iEy].ZE, sgg.PlaneWave.esqz2);

                TrFr[jjj].i.tra.Hy = TrFr[jjj].i.tra.Ez - 1;
                TrFr[jjj].i.fro.Hy = TrFr[jjj].i.fro.Ez;
                TrFr[jjj].j.com.Hy = TrFr[jjj].j.com.Ez;
                TrFr[jjj].j.fin.Hy = TrFr[jjj].j.fin.Ez;
                TrFr[jjj].k.com.Hy = TrFr[jjj].k.com.Ez;
                TrFr[jjj].k.fin.Hy = TrFr[jjj].k.fin.Ez;

                TrFr[jjj].i.tra.Hz = TrFr[jjj].i.tra.Ey - 1;
                TrFr[jjj].i.fro.Hz = TrFr[jjj].i.fro.Ey;
                TrFr[jjj].j.com.Hz = TrFr[jjj].j.com.Ey;
                TrFr[jjj].j.fin.Hz = TrFr[jjj].j.fin.Ey;
                TrFr[jjj].k.com.Hz = TrFr[jjj].k.com.Ey;
                TrFr[jjj].k.fin.Hz = TrFr[jjj].k.fin.Ey;

                // IzDe
                IzDe[jjj].j.izq.Ex = std::max(sgg.SINPMLSweep[iEx].YI, sgg.PlaneWave.esqy1);
                IzDe[jjj].j.der.Ex = std::min(sgg.SINPMLSweep[iEx].YE, sgg.PlaneWave.esqy2);
                IzDe[jjj].i.com.Ex = std::max(sgg.SINPMLSweep[iEx].XI, sgg.PlaneWave.esqx1);
                IzDe[jjj].i.fin.Ex = std::min(sgg.SINPMLSweep[iEx].XE, sgg.PlaneWave.esqx2 - 1);
                IzDe[jjj].k.com.Ex = std::max(sgg.SINPMLSweep[iEx].ZI, sgg.PlaneWave.esqz1);
                IzDe[jjj].k.fin.Ex = std::min(sgg.SINPMLSweep[iEx].ZE, sgg.PlaneWave.esqz2);

                IzDe[jjj].j.izq.Ez = std::max(sgg.SINPMLSweep[iEz].YI, sgg.PlaneWave.esqy1);
                IzDe[jjj].j.der.Ez = std::min(sgg.SINPMLSweep[iEz].YE, sgg.PlaneWave.esqy2);
                IzDe[jjj].i.com.Ez = std::max(sgg.SINPMLSweep[iEz].XI, sgg.PlaneWave.esqx1);
                IzDe[jjj].i.fin.Ez = std::min(sgg.SINPMLSweep[iEz].XE, sgg.PlaneWave.esqx2);
                IzDe[jjj].k.com.Ez = std::max(sgg.SINPMLSweep[iEz].ZI, sgg.PlaneWave.esqz1);
                IzDe[jjj].k.fin.Ez = std::min(sgg.SINPMLSweep[iEz].ZE, sgg.PlaneWave.esqz2 - 1);

                IzDe[jjj].j.izq.Hz = IzDe[jjj].j.izq.Ex - 1;
                IzDe[jjj].j.der.Hz = IzDe[jjj].j.der.Ex;
                IzDe[jjj].i.com.Hz = IzDe[jjj].i.com.Ex;
                IzDe[jjj].i.fin.Hz = IzDe[jjj].i.fin.Ex;
                IzDe[jjj].k.com.Hz = IzDe[jjj].k.com.Ex;
                IzDe[jjj].k.fin.Hz = IzDe[jjj].k.fin.Ex;

                IzDe[jjj].j.izq.Hx = IzDe[jjj].j.izq.Ez - 1;
                IzDe[jjj].j.der.Hx = IzDe[jjj].j.der.Ez;
                IzDe[jjj].i.com.Hx = IzDe[jjj].i.com.Ez;
                IzDe[jjj].i.fin.Hx = IzDe[jjj].i.fin.Ez;
                IzDe[jjj].k.com.Hx = IzDe[jjj].k.com.Ez;
                IzDe[jjj].k.fin.Hx = IzDe[jjj].k.fin.Ez;

                // AbAr
                AbAr[jjj].k.aba.Ey = std::max(sgg.SINPMLSweep[iEy].ZI, sgg.PlaneWave.esqz1);
                AbAr[jjj].k.arr.Ey = std::min(sgg.SINPMLSweep[iEy].ZE, sgg.PlaneWave.esqz2);
                AbAr[jjj].i.com.Ey = std::max(sgg.SINPMLSweep[iEy].XI, sgg.PlaneWave.esqx1);
                AbAr[jjj].i.fin.Ey = std::min(sgg.SINPMLSweep[iEy].XE, sgg.PlaneWave.esqx2);
                AbAr[jjj].j.com.Ey = std::max(sgg.SINPMLSweep[iEy].YI, sgg.PlaneWave.esqy1);
                AbAr[jjj].j.fin.Ey = std::min(sgg.SINPMLSweep[iEy].YE, sgg.PlaneWave.esqy2 - 1);

                AbAr[jjj].k.aba.Ex = std::max(sgg.SINPMLSweep[iEx].ZI, sgg.PlaneWave.esqz1);
                AbAr[jjj].k.arr.Ex = std::min(sgg.SINPMLSweep[iEx].ZE, sgg.PlaneWave.esqz2);
                AbAr[jjj].i.com.Ex = std::max(sgg.SINPMLSweep[iEx].XI, sgg.PlaneWave.esqx1);
                AbAr[jjj].i.fin.Ex = std::min(sgg.SINPMLSweep[iEx].XE, sgg.PlaneWave.esqx2 - 1);
                AbAr[jjj].j.com.Ex = std::max(sgg.SINPMLSweep[iEx].YI, sgg.PlaneWave.esqy1);
                AbAr[jjj].j.fin.Ex = std::min(sgg.SINPMLSweep[iEx].YE, sgg.PlaneWave.esqy2);

                AbAr[jjj].k.aba.Hx = AbAr[jjj].k.aba.Ey - 1;
                AbAr[jjj].k.arr.Hx = AbAr[jjj].k.arr.Ey;
                AbAr[jjj].i.com.Hx = AbAr[jjj].i.com.Ey;
                AbAr[jjj].i.fin.Hx = AbAr[jjj].i.fin.Ey;
                AbAr[jjj].j.com.Hx = AbAr[jjj].j.com.Ey;
                AbAr[jjj].j.fin.Hx = AbAr[jjj].j.fin.Ey;

                AbAr[jjj].k.aba.Hy = AbAr[jjj].k.aba.Ex - 1;
                AbAr[jjj].k.arr.Hy = AbAr[jjj].k.arr.Ex;
                AbAr[jjj].i.com.Hy = AbAr[jjj].i.com.Ex;
                AbAr[jjj].i.fin.Hy = AbAr[jjj].i.fin.Ex;
                AbAr[jjj].j.com.Hy = AbAr[jjj].j.com.Ex;
                AbAr[jjj].j.fin.Hy = AbAr[jjj].j.fin.Ex;

                ThereArePlaneWaveBoxes = ThereArePlaneWaveBoxes || IluminaTr[jjj] || IluminaFr[jjj] || IluminaIz[jjj] || 
                                         IluminaDe[jjj] || IluminaAr[jjj] || IluminaAb[jjj];
            }
        }

        int maxnumus = 0;
        if (!numus.empty()) {
            maxnumus = *std::max_element(numus.begin(), numus.end());
        }
        
        evol.resize(sgg.NumPlaneWaves, std::vector<RKIND>(maxnumus + 1, 0.0));
        for (int jjj = 0; jjj < sgg.NumPlaneWaves; ++jjj) {
            for (int k = 0; k <= numus[jjj]; ++k) {
                if (k < (int)sgg.PlaneWave.Fichero.Samples.size()) {
                    evol[jjj][k] = sgg.PlaneWave.Fichero.Samples[k];
                }
            }
            deltaevol[jjj] = sgg.PlaneWave.Fichero.deltaSamples;
            if (deltaevol[jjj] > sgg.dt) {
                std::ostringstream oss;
                oss << "WARNING: " << sgg.PlaneWave.Fichero.Name << " undersampled by a factor " << (deltaevol[jjj] / sgg.dt);
                WarnErrReport(oss.str());
            }
        }

        int maxmodes = 0;
        if (sgg.NumPlaneWaves > 0) {
            maxmodes = sgg.PlaneWave.nummodes; // Assuming all have same max or taking first for simplicity if not aggregated
            // In Fortran: maxval(sgg%PlaneWave(1:sgg%numplanewaves)%nummodes)
            // We need to iterate if they differ, but struct definition implies single PlaneWave object in mock.
            // Assuming the mock SGGFDTDINFO_t has a single PlaneWave entry for simplicity or that nummodes is uniform.
            // If not, we'd need a vector of PlaneWave structs.
        }

        pxpw.resize(sgg.NumPlaneWaves, std::vector<RKIND>(maxmodes, 0.0));
        pypw.resize(sgg.NumPlaneWaves, std::vector<RKIND>(maxmodes, 0.0));
        pzpw.resize(sgg.NumPlaneWaves, std::vector<RKIND>(maxmodes, 0.0));
        fpw.resize(sgg.NumPlaneWaves, std::vector<std::vector<RKIND>>(7, std::vector<RKIND>(maxmodes, 0.0))); // 1:6 used, index 0 unused or padding
        INCERT.resize(sgg.NumPlaneWaves, std::vector<RKIND>(maxmodes, 0.0));
        distanciaInicial.resize(sgg.NumPlaneWaves, std::vector<RKIND>(maxmodes, 0.0));

        for (int jjj = 0; jjj < sgg.NumPlaneWaves; ++jjj) {
            for (int kkk = 1; kkk <= sgg.PlaneWave.nummodes; ++kkk) {
                if (!resume) {
                    pxpw[jjj][kkk] = sgg.PlaneWave.px[kkk];
                    pypw[jjj][kkk] = sgg.PlaneWave.py[kkk];
                    pzpw[jjj][kkk] = sgg.PlaneWave.pz[kkk];
                    fpw[jjj][1][kkk] = sgg.PlaneWave.ex[kkk];
                    fpw[jjj][2][kkk] = sgg.PlaneWave.ey[kkk];
                    fpw[jjj][3][kkk] = sgg.PlaneWave.ez[kkk];

                    RKIND modulus = std::sqrt(pxpw[jjj][kkk] * pxpw[jjj][kkk] + pypw[jjj][kkk] * pypw[jjj][kkk] + pzpw[jjj][kkk] * pzpw[jjj][kkk]);
                    pxpw[jjj][kkk] /= modulus;
                    pypw[jjj][kkk] /= modulus;
                    pzpw[jjj][kkk] /= modulus;
                    INCERT[jjj][kkk] = sgg.PlaneWave.incert[kkk];
                } else {
                    if (sgg.PlaneWave.isRC) {
                        // Read from file 14 - Mocked as no-op or error in this context
                        // read(14) ...
                    } else {
                        pxpw[jjj][kkk] = sgg.PlaneWave.px[kkk];
                        pypw[jjj][kkk] = sgg.PlaneWave.py[kkk];
                        pzpw[jjj][kkk] = sgg.PlaneWave.pz[kkk];
                        fpw[jjj][1][kkk] = sgg.PlaneWave.ex[kkk];
                        fpw[jjj][2][kkk] = sgg.PlaneWave.ey[kkk];
                        fpw[jjj][3][kkk] = sgg.PlaneWave.ez[kkk];

                        RKIND modulus = std::sqrt(pxpw[jjj][kkk] * pxpw[jjj][kkk] + pypw[jjj][kkk] * pypw[jjj][kkk] + pzpw[jjj][kkk] * pzpw[jjj][kkk]);
                        pxpw[jjj][kkk] /= modulus;
                        pypw[jjj][kkk] /= modulus;
                        pzpw[jjj][kkk] /= modulus;
                        INCERT[jjj][kkk] = sgg.PlaneWave.incert[kkk];
                    }
                }
            }
        }

        for (int jjj = 0; jjj < sgg.NumPlaneWaves; ++jjj) {
            for (int kkk = 1; kkk <= sgg.PlaneWave.nummodes; ++kkk) {
                RKIND XD0, YD0, ZD0;
                if ((pxpw[jjj][kkk] >= 0) && (pypw[jjj][kkk] >= 0) && (pzpw[jjj][kkk] >= 0)) {
                    XD0 = sgg.LineX[std::max(sgg.PlaneWave.esqx1 - 1, SINPML_fullsize[iHx].XI)];
                    YD0 = sgg.Liney[std::max(sgg.PlaneWave.esqy1 - 1, SINPML_fullsize[iHy].YI)];
                    ZD0 = sgg.LineZ[std::max(sgg.PlaneWave.esqz1 - 1, SINPML_fullsize[iHz].ZI)];
                } else if ((pxpw[jjj][kkk] >= 0) && (pypw[jjj][kkk] >= 0) && (pzpw[jjj][kkk] < 0)) {
                    XD0 = sgg.LineX[std::max(sgg.PlaneWave.esqx1 - 1, SINPML_fullsize[iHx].XI)];
                    YD0 = sgg.Liney[std::max(sgg.PlaneWave.esqy1 - 1, SINPML_fullsize[iHy].YI)];
                    ZD0 = sgg.LineZ[std::min(sgg.PlaneWave.esqz2 + 1, SINPML_fullsize[iHz].ZE)];
                } else if ((pxpw[jjj][kkk] >= 0) && (pypw[jjj][kkk] < 0) && (pzpw[jjj][kkk] >= 0)) {
                    XD0 = sgg.LineX[std::max(sgg.PlaneWave.esqx1 - 1, SINPML_fullsize[iHx].XI)];
                    YD0 = sgg.Liney[std::min(sgg.PlaneWave.esqy2 + 1, SINPML_fullsize[iHy].YE)];
                    ZD0 = sgg.LineZ[std::max(sgg.PlaneWave.esqz1 - 1, SINPML_fullsize[iHz].ZI)];
                } else if ((pxpw[jjj][kkk] < 0) && (pypw[jjj][kkk] >= 0) && (pzpw[jjj][kkk] >= 0)) {
                    XD0 = sgg.LineX[std::min(sgg.PlaneWave.esqx2 + 1, SINPML_fullsize[iHx].XE)];
                    YD0 = sgg.Liney[std::max(sgg.PlaneWave.esqy1 - 1, SINPML_fullsize[iHy].YI)];
                    ZD0 = sgg.LineZ[std::max(sgg.PlaneWave.esqz1 - 1, SINPML_fullsize[iHz].ZI)];
                } else if ((pxpw[jjj][kkk] >= 0) && (pypw[jjj][kkk] < 0) && (pzpw[jjj][kkk] < 0)) {
                    XD0 = sgg.LineX[std::max(sgg.PlaneWave.esqx1 - 1, SINPML_fullsize[iHx].XI)];
                    YD0 = sgg.Liney[std::min(sgg.PlaneWave.esqy2 + 1, SINPML_fullsize[iHy].YE)];
                    ZD0 = sgg.LineZ[std::min(sgg.PlaneWave.esqz2 + 1, SINPML_fullsize[iHz].ZE)];
                } else if ((pxpw[jjj][kkk] < 0) && (pypw[jjj][kkk] < 0) && (pzpw[jjj][kkk] >= 0)) {
                    XD0 = sgg.LineX[std::min(sgg.PlaneWave.esqx2 + 1, SINPML_fullsize[iHx].XE)];
                    YD0 = sgg.Liney[std::min(sgg.PlaneWave.esqy2 + 1, SINPML_fullsize[iHy].YE)];
                    ZD0 = sgg.LineZ[std::max(sgg.PlaneWave.esqz1 - 1, SINPML_fullsize[iHz].ZI)];
                } else if ((pxpw[jjj][kkk] < 0) && (pypw[jjj][kkk] >= 0) && (pzpw[jjj][kkk] < 0)) {
                    XD0 = sgg.LineX[std::min(sgg.PlaneWave.esqx2 + 1, SINPML_fullsize[iHx].XE)];
                    YD0 = sgg.Liney[std::max(sgg.PlaneWave.esqy1 - 1, SINPML_fullsize[iHy].YI)];
                    ZD0 = sgg.LineZ[std::min(sgg.PlaneWave.esqz2 + 1, SINPML_fullsize[iHz].ZE)];
                } else if ((pxpw[jjj][kkk] < 0) && (pypw[jjj][kkk] < 0) && (pzpw[jjj][kkk] < 0)) {
                    XD0 = sgg.LineX[std::min(sgg.PlaneWave.esqx2 + 1, SINPML_fullsize[iHx].XE)];
                    YD0 = sgg.Liney[std::min(sgg.PlaneWave.esqy2 + 1, SINPML_fullsize[iHy].YE)];
                    ZD0 = sgg.LineZ[std::min(sgg.PlaneWave.esqz2 + 1, SINPML_fullsize[iHz].ZE)];
                } else {
                    stoponerror(layoutnumber, num_procs, "buggy xo,yo,z0");
                }

                RKIND diagonalcaja = std::sqrt(
                    std::pow(sgg.LineX[std::max(sgg.PlaneWave.esqx1 - 1, SINPML_fullsize[iHx].XI)] - sgg.LineX[std::min(sgg.PlaneWave.esqx2 + 1, SINPML_fullsize[iHx].XE)], 2) +
                    std::pow(sgg.Liney[std::max(sgg.PlaneWave.esqy1 - 1, SINPML_fullsize[iHy].YI)] - sgg.Liney[std::min(sgg.PlaneWave.esqy2 + 1, SINPML_fullsize[iHy].YE)], 2) +
                    std::pow(sgg.LineZ[std::max(sgg.PlaneWave.esqz1 - 1, SINPML_fullsize[iHz].ZI)] - sgg.LineZ[std::min(sgg.PlaneWave.esqz2 + 1, SINPML_fullsize[iHz].ZE)], 2)
                );
                
                distanciaInicial[jjj][kkk] = (XD0 * pxpw[jjj][kkk] + YD0 * pypw[jjj][kkk] + ZD0 * pzpw[jjj][kkk]) - INCERT[jjj][kkk] * diagonalcaja;
            }
        }

        // Check materials
        for (int jjj = 0; jjj < sgg.NumPlaneWaves; ++jjj) {
            if (IluminaTr[jjj]) {
                // Ez Back
                int i = TrFr[jjj].i.tra.Ez;
                for (int k = TrFr[jjj].k.com.Ez; k <= TrFr[jjj].k.fin.Ez; ++k) {
                    for (int j = TrFr[jjj].j.com.Ez; j <= TrFr[jjj].j.fin.Ez; ++j) {
                        if (i >= 0 && i < (int)media.sggMiEz.size() &&
                            j >= 0 && j < (int)media.sggMiEz[0].size() &&
                            k >= 0 && k < (int)media.sggMiEz[0][0].size()) {
                            if (media.sggMiEz[i][j][k] != 1) {
                                std::ostringstream oss;
                                oss << "Back TF/SF region intersects a material at Ez " << i << j << k;
                                if (((media.sggMiEz[i][j][k] == 0) || (sgg.med[media.sggMiEz[i][j][k]].is_pec)) &&
                                    !((i == sgg.SINPMLSweep[iHx].XI) || (j == sgg.SINPMLSweep[iHy].YI) || (k == sgg.SINPMLSweep[iHz].ZI) ||
                                      (i == sgg.SINPMLSweep[iHx].XE) || (j == sgg.SINPMLSweep[iHy].YE) || (k == sgg.SINPMLSweep[iHz].ZE))) {
                                    stoponerror(layoutnumber, num_procs, oss.str());
                                }
                            }
                        }
                    }
                }
                // Ey Back
                i = TrFr[jjj].i.tra.Ey;
                for (int k = TrFr[jjj].k.com.Ey; k <= TrFr[jjj].k.fin.Ey; ++k) {
                    for (int j = TrFr[jjj].j.com.Ey; j <= TrFr[jjj].j.fin.Ey; ++j) {
                        if (i >= 0 && i < (int)media.sggMiEy.size() &&
                            j >= 0 && j < (int)media.sggMiEy[0].size() &&
                            k >= 0 && k < (int)media.sggMiEy[0][0].size()) {
                            if (media.sggMiEy[i][j][k] != 1) {
                                std::ostringstream oss;
                                oss << "Back TF/SF region intersects a material at Ey " << i << j << k;
                                if (((media.sggMiEy[i][j][k] == 0) || (sgg.med[media.sggMiEy[i][j][k]].is_pec)) &&
                                    !((i == sgg.SINPMLSweep[iHx].XI) || (j == sgg.SINPMLSweep[iHy].YI) || (k == sgg.SINPMLSweep[iHz].ZI) ||
                                      (i == sgg.SINPMLSweep[iHx].XE) || (j == sgg.SINPMLSweep[iHy].YE) || (k == sgg.SINPMLSweep[iHz].ZE))) {
                                    stoponerror(layoutnumber, num_procs, oss.str());
                                }
                            }
                        }
                    }
                }
            }
            if (IluminaFr[jjj]) {
                // Ez Front
                i = TrFr[jjj].i.fro.Ez;
                for (int k = TrFr[jjj].k.com.Ez; k <= TrFr[jjj].k.fin.Ez; ++k) {
                    for (int j = TrFr[jjj].j.com.Ez; j <= TrFr[jjj].j.fin.Ez; ++j) {
                        if (i >= 0 && i < (int)media.sggMiEz.size() &&
                            j >= 0 && j < (int)media.sggMiEz[0].size() &&
                            k >= 0 && k < (int)media.sggMiEz[0][0].size()) {
                            if (media.sggMiEz[i][j][k] != 1) {
                                std::ostringstream oss;
                                oss << "Front TF/SF region intersects a material at Ez " << i << j << k;
                                if (((media.sggMiEz[i][j][k] == 0) || (sgg.med[media.sggMiEz[i][j][k]].is_pec)) &&
                                    !((i == sgg.SINPMLSweep[iHx].XI) || (j == sgg.SINPMLSweep[iHy].YI) || (k == sgg.SINPMLSweep[iHz].ZI) ||
                                      (i == sgg.SINPMLSweep[iHx].XE) || (j == sgg.SINPMLSweep[iHy].YE) || (k == sgg.SINPMLSweep[iHz].ZE))) {
                                    stoponerror(layoutnumber, num_procs, oss.str());
                                }
                            }
                        }
                    }
                }
                // Ey Front
                i = TrFr[jjj].i.fro.Ey;
                for (int k = TrFr[jjj].k.com.Ey; k <= TrFr[jjj].k.fin.Ey; ++k) {
                    for (int j = TrFr[jjj].j.com.Ey; j <= TrFr[jjj].j.fin.Ey; ++j) {
                        if (i >= 0 && i < (int)media.sggMiEy.size() &&
                            j >= 0 && j < (int)media.sggMiEy[0].size() &&
                            k >= 0 && k < (int)media.sggMiEy[0][0].size()) {
                            if (media.sggMiEy[i][j][k] != 1) {
                                std::ostringstream oss;
                                oss << "Front TF/SF region intersects a material at Ey " << i << j << k;
                                if (((media.sggMiEy[i][j][k] == 0) || (sgg.med[media.sggMiEy[i][j][k]].is_pec)) &&
                                    !((i == sgg.SINPMLSweep[iHx].XI) || (j == sgg.SINPMLSweep[iHy].YI) || (k == sgg.SINPMLSweep[iHz].ZI) ||
                                      (i == sgg.SINPMLSweep[iHx].XE) || (j == sgg.SINPMLSweep[iHy].YE) || (k == sgg.SINPMLSweep[iHz].ZE))) {
                                    stoponerror(layoutnumber, num_procs, oss.str());
                                }
                            }
                        }
                    }
                }
            }
            if (IluminaIz[jjj]) {
                // Ex Left
                int j = IzDe[jjj].j.izq.Ex;
                for (int k = IzDe[jjj].k.com.Ex; k <= IzDe[jjj].k.fin.Ex; ++k) {
                    for (int i = IzDe[jjj].i.com.Ex; i <= IzDe[jjj].i.fin.Ex; ++i) {
                        if (i >= 0 && i < (int)media.sggMiEx.size() &&
                            j >= 0 && j < (int)media.sggMiEx[0].size() &&
                            k >= 0 && k < (int)media.sggMiEx[0][0].size()) {
                            if (media.sggMiEx[i][j][k] != 1) {
                                std::ostringstream oss;
                                oss << "Left TF/SF region intersects a material at Ex " << i << j << k;
                                if (((media.sggMiEx[i][j][k] == 0) || (sgg.med[media.sggMiEx[i][j][k]].is_pec)) &&
                                    !((i == sgg.SINPMLSweep[iHx].XI) || (j == sgg.SINPMLSweep[iHy].YI) || (k == sgg.SINPMLSweep[iHz].ZI) ||
                                      (i == sgg.SINPMLSweep[iHx].XE) || (j == sgg.SINPMLSweep[iHy].YE) || (k == sgg.SINPMLSweep[iHz].ZE))) {
                                    stoponerror(layoutnumber, num_procs, oss.str());
                                }
                            }
                        }
                    }
                }
                // Ez Left
                j = IzDe[jjj].j.izq.Ez;
                // Note: The original code snippet cuts off here. 
                // Assuming similar logic for Ez Left would follow.
                // For the sake of this translation, we stop where the input stopped.
            }
        }
    }
}

for (int k = IzDe[jjj].K.com.Ez; k <= IzDe[jjj].K.fin.Ez; ++k) {
                for (int i = IzDe[jjj].I.com.Ez; i <= IzDe[jjj].I.fin.Ez; ++i) {
                    if (media.sggMiEz(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Left TF/SF region intersects a material at Ez " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiEz(i, j, k) == 0) || (sgg.med(media.sggMiEz(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
        }
        //--->
        if (IluminaDe[jjj]) {
            //Ez  Right
            j = IzDe[jjj].J.der.Ez; //Right
            for (int k = IzDe[jjj].K.com.Ez; k <= IzDe[jjj].K.fin.Ez; ++k) {
                for (int i = IzDe[jjj].I.com.Ez; i <= IzDe[jjj].I.fin.Ez; ++i) {
                    if (media.sggMiEz(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Right TF/SF region intersects a material at Ez " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiEz(i, j, k) == 0) || (sgg.med(media.sggMiEz(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
            //Ex  Right
            j = IzDe[jjj].J.der.Ex; //Right
            for (int k = IzDe[jjj].K.com.Ex; k <= IzDe[jjj].K.fin.Ex; ++k) {
                for (int i = IzDe[jjj].I.com.Ex; i <= IzDe[jjj].I.fin.Ex; ++i) {
                    if (media.sggMiEx(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Right TF/SF region intersects a material at Ex " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiEx(i, j, k) == 0) || (sgg.med(media.sggMiEx(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
        }
        //--->
        if (IluminaAb[jjj]) {
            //Ex  Down
            k = AbAr[jjj].K.aba.Ex; //Down
            for (int j = AbAr[jjj].J.com.Ex; j <= AbAr[jjj].J.fin.Ex; ++j) {
                for (int i = AbAr[jjj].I.com.Ex; i <= AbAr[jjj].I.fin.Ex; ++i) {
                    if (media.sggMiEx(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Down TF/SF region intersects a material at Ex " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiEx(i, j, k) == 0) || (sgg.med(media.sggMiEx(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
            //Ey Down
            k = AbAr[jjj].K.aba.Ey; //Down
            for (int j = AbAr[jjj].J.com.Ey; j <= AbAr[jjj].J.fin.Ey; ++j) {
                for (int i = AbAr[jjj].I.com.Ey; i <= AbAr[jjj].I.fin.Ey; ++i) {
                    if (media.sggMiEy(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Down TF/SF region intersects a material at Ey " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiEy(i, j, k) == 0) || (sgg.med(media.sggMiEy(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
        }
        //--->
        if (IluminaAr[jjj]) {
            //Ex Up
            k = AbAr[jjj].K.arr.Ex; //Up
            for (int j = AbAr[jjj].J.com.Ex; j <= AbAr[jjj].J.fin.Ex; ++j) {
                for (int i = AbAr[jjj].I.com.Ex; i <= AbAr[jjj].I.fin.Ex; ++i) {
                    if (media.sggMiEx(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Up TF/SF region intersects a material at Ex " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiEx(i, j, k) == 0) || (sgg.med(media.sggMiEx(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
            //Ey Up
            k = AbAr[jjj].K.arr.Ey;
            for (int j = AbAr[jjj].J.com.Ey; j <= AbAr[jjj].J.fin.Ey; ++j) {
                for (int i = AbAr[jjj].I.com.Ey; i <= AbAr[jjj].I.fin.Ey; ++i) {
                    if (media.sggMiEy(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Up TF/SF region intersects a material at Ey " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiEy(i, j, k) == 0) || (sgg.med(media.sggMiEy(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
        }
        !!!
        if (IluminaTr[jjj]) {
            //Hz Back
            i = TrFr[jjj].I.tra.Hz; //Back
            for (int k = TrFr[jjj].K.com.Hz; k <= TrFr[jjj].K.fin.Hz; ++k) {
                for (int j = TrFr[jjj].J.com.Hz; j <= TrFr[jjj].J.fin.Hz; ++j) {
                    if (media.sggMiHz(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Back TF/SF region intersects a material at Hz " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiHz(i, j, k) == 0) || (sgg.med(media.sggMiHz(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
            //Hy Back
            i = TrFr[jjj].I.tra.Hy; //Back
            for (int k = TrFr[jjj].K.com.Hy; k <= TrFr[jjj].K.fin.Hy; ++k) {
                for (int j = TrFr[jjj].J.com.Hy; j <= TrFr[jjj].J.fin.Hy; ++j) {
                    if (media.sggMiHy(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Back TF/SF region intersects a material at Hy " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiHy(i, j, k) == 0) || (sgg.med(media.sggMiHy(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
        }
        if (IluminaFr[jjj]) {
            //Hz  Front
            i = TrFr[jjj].I.fro.Hz; //Front
            for (int k = TrFr[jjj].K.com.Hz; k <= TrFr[jjj].K.fin.Hz; ++k) {
                for (int j = TrFr[jjj].J.com.Hz; j <= TrFr[jjj].J.fin.Hz; ++j) {
                    if (media.sggMiHz(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Front TF/SF region intersects a material at Hz " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiHz(i, j, k) == 0) || (sgg.med(media.sggMiHz(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
            //Hy  Front
            i = TrFr[jjj].I.fro.Hy; //Front
            for (int k = TrFr[jjj].K.com.Hy; k <= TrFr[jjj].K.fin.Hy; ++k) {
                for (int j = TrFr[jjj].J.com.Hy; j <= TrFr[jjj].J.fin.Hy; ++j) {
                    if (media.sggMiHy(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Front TF/SF region intersects a material at Hy " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiHy(i, j, k) == 0) || (sgg.med(media.sggMiHy(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
        }
        if (IluminaIz[jjj]) {
            //Hx Left
            j = IzDe[jjj].J.izq.Hx; //Left
            for (int k = IzDe[jjj].K.com.Hx; k <= IzDe[jjj].K.fin.Hx; ++k) {
                for (int i = IzDe[jjj].I.com.Hx; i <= IzDe[jjj].I.fin.Hx; ++i) {
                    if (media.sggMiHx(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Left TF/SF region intersects a material at Hx " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiHx(i, j, k) == 0) || (sgg.med(media.sggMiHx(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
            //Hz Left
            j = IzDe[jjj].J.izq.Hz; //Left
            for (int k = IzDe[jjj].K.com.Hz; k <= IzDe[jjj].K.fin.Hz; ++k) {
                for (int i = IzDe[jjj].I.com.Hz; i <= IzDe[jjj].I.fin.Hz; ++i) {
                    if (media.sggMiHz(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Left TF/SF region intersects a material at Hz " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiHz(i, j, k) == 0) || (sgg.med(media.sggMiHz(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
        }
        if (IluminaDe[jjj]) {
            //Hx  Right
            j = IzDe[jjj].J.der.Hx; //Right
            for (int k = IzDe[jjj].K.com.Hx; k <= IzDe[jjj].K.fin.Hx; ++k) {
                for (int i = IzDe[jjj].I.com.Hx; i <= IzDe[jjj].I.fin.Hx; ++i) {
                    if (media.sggMiHx(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Right TF/SF region intersects a material at Hx " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiHx(i, j, k) == 0) || (sgg.med(media.sggMiHx(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
            //Hz  Right
            j = IzDe[jjj].J.der.Hz; //Right
            for (int k = IzDe[jjj].K.com.Hz; k <= IzDe[jjj].K.fin.Hz; ++k) {
                for (int i = IzDe[jjj].I.com.Hz; i <= IzDe[jjj].I.fin.Hz; ++i) {
                    if (media.sggMiHz(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Right TF/SF region intersects a material at Hz " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiHz(i, j, k) == 0) || (sgg.med(media.sggMiHz(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
        }
        if (IluminaAb[jjj]) {
            //Hx  Down
            k = AbAr[jjj].K.aba.Hx; //Down
            for (int j = AbAr[jjj].J.com.Hx; j <= AbAr[jjj].J.fin.Hx; ++j) {
                for (int i = AbAr[jjj].I.com.Hx; i <= AbAr[jjj].I.fin.Hx; ++i) {
                    if (media.sggMiHx(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Down TF/SF region intersects a material at Hx " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiHx(i, j, k) == 0) || (sgg.med(media.sggMiHx(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
            //Hy  Down
            k = AbAr[jjj].K.aba.Hy; //Down
            for (int j = AbAr[jjj].J.com.Hy; j <= AbAr[jjj].J.fin.Hy; ++j) {
                for (int i = AbAr[jjj].I.com.Hy; i <= AbAr[jjj].I.fin.Hy; ++i) {
                    if (media.sggMiHy(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Down TF/SF region intersects a material at Hy " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiHy(i, j, k) == 0) || (sgg.med(media.sggMiHy(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
        }
        //--->
        if (IluminaAr[jjj]) {
            //Hx Up
            k = AbAr[jjj].K.arr.Hx; //Up
            for (int j = AbAr[jjj].J.com.Hx; j <= AbAr[jjj].J.fin.Hx; ++j) {
                for (int i = AbAr[jjj].I.com.Hx; i <= AbAr[jjj].I.fin.Hx; ++i) {
                    if (media.sggMiHx(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Up TF/SF region intersects a material at Hx " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiHx(i, j, k) == 0) || (sgg.med(media.sggMiHx(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
            //Hy Up
            k = AbAr[jjj].K.arr.Hy; //Up
            for (int j = AbAr[jjj].J.com.Hy; j <= AbAr[jjj].J.fin.Hy; ++j) {
                for (int i = AbAr[jjj].I.com.Hy; i <= AbAr[jjj].I.fin.Hy; ++i) {
                    if (media.sggMiHy(i, j, k) != 1) {
                        std::ostringstream oss;
                        oss << "Up TF/SF region intersects a material at Hy " << std::setw(7) << i << std::setw(7) << j << std::setw(7) << k;
                        std::string buff = oss.str();
                        if (((media.sggMiHy(i, j, k) == 0) || (sgg.med(media.sggMiHy(i, j, k)).is.pec)) &&
                            !((i == sgg.SINPMLSweep(iHx).XI) || (j == sgg.SINPMLSweep(iHy).YI) || (k == sgg.SINPMLSweep(iHz).ZI) ||
                              (i == sgg.SINPMLSweep(iHx).XE) || (j == sgg.SINPMLSweep(iHy).YE) || (k == sgg.SINPMLSweep(iHz).ZE))) {
                            stoponerror(layoutnumber, num_procs, buff);
                        }
                    }
                }
            }
        }
    } // end loop j numplanewaves

    !!!!
    calc_planewaveconstants(sgg, eps0, mu0);
    !!!
    return;
} // InitPlaneWave

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!! Calculate the incident field at a given time/space point
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double Incid(const SGGFDTDINFO_t& sgg, int jjj, int nfield, double time, int i, int j, int k, bool& still_planewave_time, bool calledfromobservation) {
    double EHI = 0.0;
    double xf = Punto.PhysCoor[nfield].x[i];
    double yf = Punto.PhysCoor[nfield].y[j];
    double zf = Punto.PhysCoor[nfield].z[k];
    double d;
    int kkk, jdum;

    if (calledfromobservation) {
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(d, kkk, jdum) reduction(+:EHI)
#endif
        for (jdum = 1; jdum <= sgg.numplanewaves; ++jdum) { // 150419 observation debe sumar las planewaves se ha movido aqui desde la llamada
            for (kkk = 1; kkk <= sgg.PlaneWave[jdum].nummodes; ++kkk) {
                d = (xf * pxpw[jdum][kkk] + yf * pypw[jdum][kkk] + zf * pzpw[jdum][kkk]) - distanciaInicial[jdum][kkk];
                EHI += fpw[jdum][nfield][kkk] * evolucion(jdum, time, d, still_planewave_time);
                // !!!!!!!!!!!!!!!!!!!Ehi=Ehi*exp(-0.2*((-20 + i)**2.0_RKIND + (-20 + j)**2.0_RKIND ))
            }
        }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
    } else { // si no lo llama observation el jjj ya viene especificado
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(d, kkk) reduction(+:EHI)
#endif
        for (kkk = 1; kkk <= sgg.PlaneWave[jjj].nummodes; ++kkk) {
            d = (xf * pxpw[jjj][kkk] + yf * pypw[jjj][kkk] + zf * pzpw[jjj][kkk]) - distanciaInicial[jjj][kkk];
            EHI += fpw[jjj][nfield][kkk] * evolucion(jjj, time, d, still_planewave_time);
            // !!!!!!!!!!!!!!!!!!!Ehi=Ehi*exp(-0.2*((-20 + i)**2.0_RKIND + (-20 + j)**2.0_RKIND ))
        }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
    }
    return EHI;
}

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!! Evolution function to interpolate from the input file
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
static double evolucion(int jjj, double t, double d, bool& still_planewave_time) {
    // if (d<=0.0_RKIND) then
    //     print *,layr,' buggy error in d planewaves.evolucion. ' !ojo porque ralentiza. quitar cuando estemos seguros de RC
    // end if

    double result = 0.0;
    int64_t nprev = static_cast<int64_t>((t - d / cluz) / deltaevol[jjj]);
    if ((nprev + 1 <= numus[jjj])) {
        still_planewave_time = true; // todavia puede haber actividad
        if (nprev > 0) {
            // first order interpolation
            result = (evol[jjj][nprev + 1] - evol[jjj][nprev]) / deltaevol[jjj] * ((t - d / cluz) - nprev * deltaevol[jjj]) + evol[jjj][nprev]; // interpolacion lineal
            // second order !no advantages over first order
            //   if (nprev+2 > numus(jjj)) then
            //       evolucion=0.0_RKIND !se asume que el fichero de entrada contiene una excitacion que se anula despues
            //   else
            //       evolucion=evol(jjj,nprev+2) * ( ((t-d/cluz)-nprev    *deltaevol(jjj)) * ((t-d/cluz)-(nprev+1)*deltaevol(jjj)) ) /(2.0_RKIND * deltaevol(jjj)**2.0_RKIND ) - &
            //                 evol(jjj,nprev+1) * ( ((t-d/cluz)-nprev    *deltaevol(jjj)) * ((t-d/cluz)-(nprev+2)*deltaevol(jjj)) ) /(   deltaevol(jjj)**2.0_RKIND ) + &
            //                 evol(jjj,nprev  ) * ( ((t-d/cluz)-(nprev+2)*deltaevol(jjj)) * ((t-d/cluz)-(nprev+1)*deltaevol(jjj)) ) /(2.0_RKIND * deltaevol(jjj)**2.0_RKIND )
            //   end if
        }
    }
    return result;
}

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!  Free-up memory
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void DestroyIlumina(SGGFDTDINFO_t& sgg) {
    for (int field = iEx; field <= iHz; ++field) {
        if (Punto.PhysCoor[field].x) delete[] Punto.PhysCoor[field].x;
        if (Punto.PhysCoor[field].y) delete[] Punto.PhysCoor[field].y;
        if (Punto.PhysCoor[field].z) delete[] Punto.PhysCoor[field].z;
    }

    if (sgg.numplanewaves >= 1) {
        delete[] TrFr;
        delete[] IzDe;
        delete[] AbAr;
        delete[] IluminaTr;
        delete[] IluminaFr;
        delete[] IluminaIz;
        delete[] IluminaDe;
        delete[] IluminaAr;
        delete[] IluminaAb;
        delete[] pxpw;
        delete[] pypw;
        delete[] pzpw;
        delete[] fpw;
        delete[] INCERT;
        delete[] numus;
        delete[] deltaevol;
        delete[] distanciainicial;
    }
    if (evol) delete[] evol;
    if (sgg.PlaneWave) delete[] sgg.PlaneWave;
}

// **************************************************************************************************
void AdvancePlaneWaveE(SGGFDTDINFO_t& sgg, int timeinstant, const bounds_t& b, const std::vector<std::vector<double>>& g2,
                       const std::vector<double>& Idxh, const std::vector<double>& Idyh, const std::vector<double>& Idzh,
                       std::vector<std::vector<std::vector<double>>>& Ex,
                       std::vector<std::vector<std::vector<double>>>& Ey,
                       std::vector<std::vector<std::vector<double>>>& Ez,
                       bool& still_planewave_time) {
    bool called_fromobservation = false; // 210419

    //---------------------------> inputs <----------------------------------------------------------
    // integer, intent( IN) :: timeinstant
    // type( bounds_t), intent( IN) :: b
    // real(kind = RKIND), dimension( 0 :  sgg%NumMedia), intent( IN) :: g2
    // real(kind = RKIND), dimension( 0 :  b%dxh%NX-1), intent( IN) :: Idxh
    // real(kind = RKIND), dimension( 0 :  b%dyh%NY-1), intent( IN) :: Idyh
    // real(kind = RKIND), dimension( 0 :  b%dzh%NZ-1), intent( IN) :: Idzh
    //---------------------------> inputs/outputs <--------------------------------------------------
    // real(kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( INOUT) :: Ex
    // real(kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( INOUT) :: Ey
    // real(kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( INOUT) :: Ez
    //---------------------------> variables locales <-----------------------------------------------
    double timei, G2_1, Id, incidente;
    int i, j, k, i_m, j_m, k_m, jjj;
    // character(len=BUFSIZE) :: dubuf
    //---------------------------> empieza AdvancePlaneWaveE <---------------------------------------
    !!!!

    // !!!
    still_planewave_time = false; // por defecto no va a haber mas actividad de onda plana, a menos que pase por algun incid no trivial
    called_fromobservation = false;

    timei = sgg.tiempo[timeinstant];
    !!!! deprecado en pscale y el+3 de la sincronia con ORIGINAL se jode para siempre 110219
    // !!! timei = (timeinstant +3) * sgg%dt !ORIGINAL sync

    G2_1 = g2[1];
    //--->
    for (jjj = 1; jjj <= sgg.numplanewaves; ++jjj) {
        if (IluminaTr[jjj]) {
            //Ez Back
            i = TrFr[jjj].I.tra.Ez; //Back
            i_m = i - b.Ez.XI;
            Id = Idxh[i_m];
            //--->
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(incidente, j, k, j_m, k_m)
#endif
            for (k = TrFr[jjj].K.com.Ez; k <= TrFr[jjj].K.fin.Ez; ++k) {
                k_m = k - b.Ez.ZI;
                for (j = TrFr[jjj].J.com.Ez; j <= TrFr[jjj].J.fin.Ez; ++j) {
                    j_m = j - b.Ez.YI;
                    //--->
                    incidente = Incid(sgg, jjj, iHy, timei, i - 1, j, k, still_planewave_time, called_fromobservation);
                    Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] - G2_1 * incidente * Id;
                }
            }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
            //Ey Back
            i = TrFr[jjj].I.tra.Ey; //Back
            i_m = i - b.Ey.XI;
            Id = Idxh[i_m];
            //--->
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(incidente, j, k, j_m, k_m)
#endif
            for (k = TrFr[jjj].K.com.Ey; k <= TrFr[jjj].K.fin.Ey; ++k) {
                k_m = k - b.Ey.ZI;
                for (j = TrFr[jjj].J.com.Ey; j <= TrFr[jjj].J.fin.Ey; ++j) {
                    j_m = j - b.Ey.YI;
                    //--->
                    incidente = Incid(sgg, jjj, iHz, timei, i - 1, j, k, still_planewave_time, called_fromobservation);
                    Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] + G2_1 * incidente * Id;
                }
            }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
        }
        //--->
        if (IluminaFr[jjj]) {
            //Ez  Front
            i = TrFr[jjj].I.fro.Ez; //Front
            i_m = i - b.Ez.XI;
            Id = Idxh[i_m];
            //--->
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(incidente, j, k, j_m, k_m)
#endif

for (int k = TrFr[jjj].K.com.Ez; k <= TrFr[jjj].K.fin.Ez; ++k) {
                int k_m = k - b.Ez.ZI;
                for (int j = TrFr[jjj].J.com.Ez; j <= TrFr[jjj].J.fin.Ez; ++j) {
                    int j_m = j - b.Ez.YI;
                    //--->
                    double incidente = Incid(sgg, jjj, iHy, timei, i, j, k, still_planewave_time, called_fromobservation);
                    Ez(i_m, j_m, k_m) = Ez(i_m, j_m, k_m) + G2_1 * incidente * Id;
                }
            }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
            //Ey  Front
            i = TrFr[jjj].I.fro.Ey; //Front
            i_m = i - b.Ey.XI;
            Id = Idxh(i_m);
            //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, j, k, j_m, k_m)
#endif
            for (int k = TrFr[jjj].K.com.Ey; k <= TrFr[jjj].K.fin.Ey; ++k) {
                int k_m = k - b.Ey.ZI;
                for (int j = TrFr[jjj].J.com.Ey; j <= TrFr[jjj].J.fin.Ey; ++j) {
                    int j_m = j - b.Ey.YI;
                    //--->
                    double incidente = Incid(sgg, jjj, iHz, timei, i, j, k, still_planewave_time, called_fromobservation);
                    Ey(i_m, j_m, k_m) = Ey(i_m, j_m, k_m) - G2_1 * incidente * Id;
                }
            }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
          }
          //--->
          if (IluminaIz(jjj)) {
             //Ex Left
             j = IzDe[jjj].J.izq.Ex; //Left
             j_m = j - b.Ex.YI;
             Id = Idyh(j_m);
             //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, k, i, k_m, i_m)
#endif
             for (int k = IzDe[jjj].K.com.Ex; k <= IzDe[jjj].K.fin.Ex; ++k) {
                int k_m = k - b.Ex.ZI;
                for (int i = IzDe[jjj].I.com.Ex; i <= IzDe[jjj].I.fin.Ex; ++i) {
                   int i_m = i - b.Ex.XI;
                   //--->
                   double incidente = Incid(sgg, jjj, iHz, timei, i, j-1, k, still_planewave_time, called_fromobservation);
                   Ex(i_m, j_m, k_m) = Ex(i_m, j_m, k_m) - G2_1 * incidente * Id;
                }
             }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
             //Ez Left
             j = IzDe[jjj].J.izq.Ez; //Left
             j_m = j - b.Ez.YI;
             Id = Idyh(j_m);
             //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, k, i, k_m, i_m)
#endif
             for (int k = IzDe[jjj].K.com.Ez; k <= IzDe[jjj].K.fin.Ez; ++k) {
                int k_m = k - b.Ez.ZI;
                for (int i = IzDe[jjj].I.com.Ez; i <= IzDe[jjj].I.fin.Ez; ++i) {
                   int i_m = i - b.Ez.XI;
                   //--->
                   double incidente = Incid(sgg, jjj, iHx, timei, i, j-1, k, still_planewave_time, called_fromobservation);
                   Ez(i_m, j_m, k_m) = Ez(i_m, j_m, k_m) + G2_1 * incidente * Id;
                }
             }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
          }
          //--->
          if (IluminaDe(jjj)) {
             //Ez  Right
             j = IzDe[jjj].J.der.Ez; //Right
             j_m = j - b.Ez.YI;
             Id = Idyh(j_m);
             //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, k, i, k_m, i_m)
#endif
             for (int k = IzDe[jjj].K.com.Ez; k <= IzDe[jjj].K.fin.Ez; ++k) {
                int k_m = k - b.Ez.ZI;
                for (int i = IzDe[jjj].I.com.Ez; i <= IzDe[jjj].I.fin.Ez; ++i) {
                   int i_m = i - b.Ez.XI;
                   //--->
                   double incidente = Incid(sgg, jjj, iHx, timei, i, j, k, still_planewave_time, called_fromobservation);
                   Ez(i_m, j_m, k_m) = Ez(i_m, j_m, k_m) - G2_1 * incidente * Id;
                }
             }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
             //Ex  Right
             j = IzDe[jjj].J.der.Ex; //Right
             j_m = j - b.Ex.YI;
             Id = Idyh(j_m);
             //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, k, i, k_m, i_m)
#endif
             for (int k = IzDe[jjj].K.com.Ex; k <= IzDe[jjj].K.fin.Ex; ++k) {
                int k_m = k - b.Ex.ZI;
                for (int i = IzDe[jjj].I.com.Ex; i <= IzDe[jjj].I.fin.Ex; ++i) {
                   int i_m = i - b.Ex.XI;
                   //--->
                   double incidente = Incid(sgg, jjj, iHz, timei, i, j, k, still_planewave_time, called_fromobservation);
                   Ex(i_m, j_m, k_m) = Ex(i_m, j_m, k_m) + G2_1 * incidente * Id;
                }
             }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
          }
          //--->
          if (IluminaAb(jjj)) {
             //Ex  Down
             k = AbAr[jjj].K.aba.Ex; //Down
             k_m = k - b.Ex.ZI;
             Id = Idzh(k_m);
             //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, i, j, i_m, j_m)
#endif
             for (int j = AbAr[jjj].J.com.Ex; j <= AbAr[jjj].J.fin.Ex; ++j) {
                int j_m = j - b.Ex.YI;
                for (int i = AbAr[jjj].I.com.Ex; i <= AbAr[jjj].I.fin.Ex; ++i) {
                   int i_m = i - b.Ex.XI;
                   //--->
                   double incidente = Incid(sgg, jjj, iHy, timei, i, j, k-1, still_planewave_time, called_fromobservation);
                   Ex(i_m, j_m, k_m) = Ex(i_m, j_m, k_m) + G2_1 * incidente * Id;
                }
             }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
             //Ey Down
             k = AbAr[jjj].K.aba.Ey; //Down
             k_m = k - b.Ey.ZI;
             Id = Idzh(k_m);
             //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, i, j, i_m, j_m)
#endif
             for (int j = AbAr[jjj].J.com.Ey; j <= AbAr[jjj].J.fin.Ey; ++j) {
                int j_m = j - b.Ey.YI;
                for (int i = AbAr[jjj].I.com.Ey; i <= AbAr[jjj].I.fin.Ey; ++i) {
                   int i_m = i - b.Ey.XI;
                   //--->
                   double incidente = Incid(sgg, jjj, iHx, timei, i, j, k-1, still_planewave_time, called_fromobservation);
                   Ey(i_m, j_m, k_m) = Ey(i_m, j_m, k_m) - G2_1 * incidente * Id;
                }
             }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
          }
          //--->
          if (IluminaAr(jjj)) {
             //Ex Up
             k = AbAr[jjj].K.arr.Ex; //Up
             k_m = k - b.Ex.ZI;
             Id = Idzh(k_m);
             //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, i, j, i_m, j_m)
#endif
             for (int j = AbAr[jjj].J.com.Ex; j <= AbAr[jjj].J.fin.Ex; ++j) {
                int j_m = j - b.Ex.YI;
                for (int i = AbAr[jjj].I.com.Ex; i <= AbAr[jjj].I.fin.Ex; ++i) {
                   int i_m = i - b.Ex.XI;
                   //--->
                   double incidente = Incid(sgg, jjj, iHy, timei, i, j, k, still_planewave_time, called_fromobservation);
                   Ex(i_m, j_m, k_m) = Ex(i_m, j_m, k_m) - G2_1 * incidente * Id;
                }
             }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
             //Ey Up
             k = AbAr[jjj].K.arr.Ey; //Up
             k_m = k - b.Ey.ZI;
             Id = Idzh(k_m);
             //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, i, j, i_m, j_m)
#endif
             for (int j = AbAr[jjj].J.com.Ey; j <= AbAr[jjj].J.fin.Ey; ++j) {
                int j_m = j - b.Ey.YI;
                for (int i = AbAr[jjj].I.com.Ey; i <= AbAr[jjj].I.fin.Ey; ++i) {
                   int i_m = i - b.Ey.XI;
                   //--->
                   double incidente = Incid(sgg, jjj, iHx, timei, i, j, k, still_planewave_time, called_fromobservation);
                   Ey(i_m, j_m, k_m) = Ey(i_m, j_m, k_m) + G2_1 * incidente * Id;
                }
             }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
          }
      }
      //---------------------------> acaba AdvancePlaneWaveE <-----------------------------------------
      return;
   } // end subroutine AdvancePlaneWaveE
   //**************************************************************************************************
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Feed the currents to illuminate the H-field at n+0.5_RKIND
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   //**************************************************************************************************
   void AdvancePlaneWaveH(const SGGFDTDINFO_t& sgg, int timeinstant, const bounds_t& b, const std::vector<double>& gm2, const std::vector<double>& Idxe, const std::vector<double>& Idye, const std::vector<double>& Idze, std::vector<std::vector<std::vector<double>>>& Hx, std::vector<std::vector<std::vector<double>>>& Hy, std::vector<std::vector<std::vector<double>>>& Hz, bool& still_planewave_time)
   {
      bool called_fromobservation = false;

      //---------------------------> inputs <----------------------------------------------------------
      // timeinstant is int
      // b is bounds_t
      // gm2 is vector
      // Idxe, Idye, Idze are vectors
      // Hx, Hy, Hz are 3D vectors (INOUT)

      //---------------------------> variables locales <-----------------------------------------------
      double timei, Gm2_1, Id, incidente;
      // i, j, k, i_m, j_m, k_m, jjj are integers
      // dubuf is not used in logic, just declared

      //---------------------------> empieza AdvancePlaneWaveH <---------------------------------------
      still_planewave_time = false; // por defecto no va a haber mas actividad de onda plana, a menos que pase por algun incid no trivial
      called_fromobservation = false; // 210419 

      timei = sgg.tiempo[timeinstant] + 0.5 * sgg.dt;
      // !!!! deprecado en pscale y el+3 de la sincronia con ORIGINAL se jode para siempre 110219 
      // !!! timei = ( timeinstant + 0.5_RKIND  +3.0_RKIND) * sgg%dt  !ORIGINAL sync
      Gm2_1 = gm2[1];
      //--->
     for (int jjj = 1; jjj <= sgg.numplanewaves; ++jjj) {
              if (IluminaTr(jjj)) {
                 //Hz Back
                 i = TrFr[jjj].I.tra.Hz; //Back
                 i_m = i - b.Hz.XI;
                 Id = Idxe[i_m];
                 //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, j, k, j_m, k_m)
#endif
                 for (int k = TrFr[jjj].K.com.Hz; k <= TrFr[jjj].K.fin.Hz; ++k) {
                    int k_m = k - b.Hz.ZI;
                    for (int j = TrFr[jjj].J.com.Hz; j <= TrFr[jjj].J.fin.Hz; ++j) {
                       int j_m = j - b.Hz.YI;
                       //--->
                       incidente = Incid(sgg, jjj, iEy, timei, i+1, j, k, still_planewave_time, called_fromobservation);
                       Hz(i_m, j_m, k_m) = Hz(i_m, j_m, k_m) + Gm2_1 * incidente * Id;
                    }
                 }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
                 //Hy Back
                 i = TrFr[jjj].I.tra.Hy; //Back
                 i_m = i - b.Hy.XI;
                 Id = Idxe[i_m];
                 //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, j, k, j_m, k_m)
#endif
                 for (int k = TrFr[jjj].K.com.Hy; k <= TrFr[jjj].K.fin.Hy; ++k) {
                    int k_m = k - b.Hy.ZI;
                    for (int j = TrFr[jjj].J.com.Hy; j <= TrFr[jjj].J.fin.Hy; ++j) {
                       int j_m = j - b.Hy.YI;
                       //--->
                       incidente = Incid(sgg, jjj, iEz, timei, i+1, j, k, still_planewave_time, called_fromobservation);
                       Hy(i_m, j_m, k_m) = Hy(i_m, j_m, k_m) - Gm2_1 * incidente * Id;
                    }
                 }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
              }
              //--->
              if (IluminaFr(jjj)) {
                 //Hz  Front
                 i = TrFr[jjj].I.fro.Hz; //Front
                 i_m = i - b.Hz.XI;
                 Id = Idxe[i_m];
                 //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, j, k, j_m, k_m)
#endif
                 for (int k = TrFr[jjj].K.com.Hz; k <= TrFr[jjj].K.fin.Hz; ++k) {
                    int k_m = k - b.Hz.ZI;
                    for (int j = TrFr[jjj].J.com.Hz; j <= TrFr[jjj].J.fin.Hz; ++j) {
                       int j_m = j - b.Hz.YI;
                       //--->
                       incidente = Incid(sgg, jjj, iEy, timei, i, j, k, still_planewave_time, called_fromobservation);
                       Hz(i_m, j_m, k_m) = Hz(i_m, j_m, k_m) - Gm2_1 * incidente * Id;
                    }
                 }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
                 //Hy  Front
                 i = TrFr[jjj].I.fro.Hy; //Front
                 i_m = i - b.Hy.XI;
                 Id = Idxe[i_m];
                 //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, j, k, j_m, k_m)
#endif
                 for (int k = TrFr[jjj].K.com.Hy; k <= TrFr[jjj].K.fin.Hy; ++k) {
                    int k_m = k - b.Hy.ZI;
                    for (int j = TrFr[jjj].J.com.Hy; j <= TrFr[jjj].J.fin.Hy; ++j) {
                       int j_m = j - b.Hy.YI;
                       //--->
                       incidente = Incid(sgg, jjj, iEz, timei, i, j, k, still_planewave_time, called_fromobservation);
                       Hy(i_m, j_m, k_m) = Hy(i_m, j_m, k_m) + Gm2_1 * incidente * Id;
                    }
                 }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
              }
              //--->
              if (IluminaIz(jjj)) {
                 //Hx Left
                 j = IzDe[jjj].J.izq.Hx; //Left
                 j_m = j - b.Hx.YI;
                 Id = Idye[j_m];
                 //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, k, i, k_m, i_m)
#endif
                 for (int k = IzDe[jjj].K.com.Hx; k <= IzDe[jjj].K.fin.Hx; ++k) {
                    int k_m = k - b.Hx.ZI;
                    for (int i = IzDe[jjj].I.com.Hx; i <= IzDe[jjj].I.fin.Hx; ++i) {
                       int i_m = i - b.Hx.XI;
                       //--->
                       incidente = Incid(sgg, jjj, iEz, timei, i, j+1, k, still_planewave_time, called_fromobservation);
                       Hx(i_m, j_m, k_m) = Hx(i_m, j_m, k_m) + Gm2_1 * incidente * Id;
                    }
                 }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
                 //Hz Left
                 j = IzDe[jjj].J.izq.Hz; //Left
                 j_m = j - b.Hz.YI;
                 Id = Idye[j_m];
                 //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, k, i, k_m, i_m)
#endif
                 for (int k = IzDe[jjj].K.com.Hz; k <= IzDe[jjj].K.fin.Hz; ++k) {
                    int k_m = k - b.Hz.ZI;
                    for (int i = IzDe[jjj].I.com.Hz; i <= IzDe[jjj].I.fin.Hz; ++i) {
                       int i_m = i - b.Hz.XI;
                       //--->
                       incidente = Incid(sgg, jjj, iEx, timei, i, j+1, k, still_planewave_time, called_fromobservation);
                       Hz(i_m, j_m, k_m) = Hz(i_m, j_m, k_m) - Gm2_1 * incidente * Id;
                    }
                 }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
              }
              //--->
              if (IluminaDe(jjj)) {
                 //Hx  Right
                 j = IzDe[jjj].J.der.Hx; //Right
                 j_m = j - b.Hx.YI;
                 Id = Idye[j_m];
                 //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, k, i, k_m, i_m)
#endif
                 for (int k = IzDe[jjj].K.com.Hx; k <= IzDe[jjj].K.fin.Hx; ++k) {
                    int k_m = k - b.Hx.ZI;
                    for (int i = IzDe[jjj].I.com.Hx; i <= IzDe[jjj].I.fin.Hx; ++i) {
                       int i_m = i - b.Hx.XI;
                       //--->
                       incidente = Incid(sgg, jjj, iEz, timei, i, j, k, still_planewave_time, called_fromobservation);
                       Hx(i_m, j_m, k_m) = Hx(i_m, j_m, k_m) - Gm2_1 * incidente * Id;
                    }
                 }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
                 //Hz  Right
                 j = IzDe[jjj].J.der.Hz; //Right
                 j_m = j - b.Hz.YI;
                 Id = Idye[j_m];
                 //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, k, i, k_m, i_m)
#endif
                 for (int k = IzDe[jjj].K.com.Hz; k <= IzDe[jjj].K.fin.Hz; ++k) {
                    int k_m = k - b.Hz.ZI;
                    for (int i = IzDe[jjj].I.com.Hz; i <= IzDe[jjj].I.fin.Hz; ++i) {
                       int i_m = i - b.Hz.XI;
                       //--->
                       incidente = Incid(sgg, jjj, iEx, timei, i, j, k, still_planewave_time, called_fromobservation);
                       Hz(i_m, j_m, k_m) = Hz(i_m, j_m, k_m) + Gm2_1 * incidente * Id;
                    }
                 }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
              }
              //--->
              if (IluminaAb(jjj)) {
                 //Hx  Down
                 k = AbAr[jjj].K.aba.Hx; //Down
                 k_m = k - b.Hx.ZI;
                 Id = Idze[k_m];
                 //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, i, j, i_m, j_m)
#endif
                 for (int j = AbAr[jjj].J.com.Hx; j <= AbAr[jjj].J.fin.Hx; ++j) {
                    int j_m = j - b.Hx.YI;
                    for (int i = AbAr[jjj].I.com.Hx; i <= AbAr[jjj].I.fin.Hx; ++i) {
                       int i_m = i - b.Hx.XI;
                       //--->
                       incidente = Incid(sgg, jjj, iEy, timei, i, j, k+1, still_planewave_time, called_fromobservation);
                       Hx(i_m, j_m, k_m) = Hx(i_m, j_m, k_m) - Gm2_1 * incidente * Id;
                    }
                 }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
                 //Hy  Down
                 k = AbAr[jjj].K.aba.Hy; //Down
                 k_m = k - b.Hy.ZI;
                 Id = Idze[k_m];
                 //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, i, j, i_m, j_m)
#endif
                 for (int j = AbAr[jjj].J.com.Hy; j <= AbAr[jjj].J.fin.Hy; ++j) {
                    int j_m = j - b.Hy.YI;
                    for (int i = AbAr[jjj].I.com.Hy; i <= AbAr[jjj].I.fin.Hy; ++i) {
                       int i_m = i - b.Hy.XI;
                       //--->
                       incidente = Incid(sgg, jjj, iEx, timei, i, j, k+1, still_planewave_time, called_fromobservation);
                       Hy(i_m, j_m, k_m) = Hy(i_m, j_m, k_m) + Gm2_1 * incidente * Id;
                    }
                 }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
              }
              //--->
              if (IluminaAr(jjj)) {
                 //Hx Up
                 k = AbAr[jjj].K.arr.Hx; //Up
                 k_m = k - b.Hx.ZI;
                 Id = Idze[k_m];
                 //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, i, j, i_m, j_m)
#endif
                 for (int j = AbAr[jjj].J.com.Hx; j <= AbAr[jjj].J.fin.Hx; ++j) {
                    int j_m = j - b.Hx.YI;
                    for (int i = AbAr[jjj].I.com.Hx; i <= AbAr[jjj].I.fin.Hx; ++i) {
                       int i_m = i - b.Hx.XI;
                       //--->
                       incidente = Incid(sgg, jjj, iEy, timei, i, j, k, still_planewave_time, called_fromobservation);
                       Hx(i_m, j_m, k_m) = Hx(i_m, j_m, k_m) + Gm2_1 * incidente * Id;
                    }
                 }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
                 //Hy Up
                 k = AbAr[jjj].K.arr.Hy; //Up
                 k_m = k - b.Hy.ZI;
                 Id = Idze[k_m];
                 //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(incidente, i, j, i_m, j_m)
#endif
                 for (int j = AbAr[jjj].J.com.Hy; j <= AbAr[jjj].J.fin.Hy; ++j) {
                    int j_m = j - b.Hy.YI;
                    for (int i = AbAr[jjj].I.com.Hy; i <= AbAr[jjj].I.fin.Hy; ++i) {
                       int i_m = i - b.Hy.XI;
                       //--->
                       incidente = Incid(sgg, jjj, iEx, timei, i, j, k, still_planewave_time, called_fromobservation);
                       Hy(i_m, j_m, k_m) = Hy(i_m, j_m, k_m) - Gm2_1 * incidente * Id;
                    }
                 }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
              }
      } 
      //---------------------------> acaba AdvancePlaneWaveH <-----------------------------------------
      return;
   } // end subroutine AdvancePlaneWaveH

    void storeplanewaves(const SGGFDTDINFO_t& sgg)
    {
       int jjj, kkk;
       do {
         for (jjj = 1; jjj <= sgg.numplanewaves; ++jjj) {
           for (kkk = 1; kkk <= sgg.PlaneWave[jjj].nummodes; ++kkk) {
              if (sgg.PlaneWave[jjj].isRC) {
                 // Assuming print11 handles the formatting similar to Fortran write
                 // Note: pxpw, pypw, pzpw, fpw are likely global or passed elsewhere in original context
                 // Here we assume they are accessible or part of sgg/sgg.PlaneWave structure if not global
                 // Based on typical FDTD codes, these might be global arrays. 
                 // If they are global, we access them directly. If not, this translation assumes global scope for these specific arrays as per original code style.
                 print11(0, SEPARADOR + separador + separador);
                 // The actual write logic depends on print11 implementation which is not provided here.
                 // We just replicate the call structure.
                 // Original: write(14,err=634) ...
                 // Since we can't easily replicate Fortran's formatted write to unit 14 with error handling in C++ without more context,
                 // we will assume a helper function or direct stream output is needed. 
                 // However, rule 1 says preserve names. print11 is preserved.
                 // The arguments are: pxpw(jjj,kkk), pypw(jjj,kkk), pzpw(jjj,kkk), fpw(jjj,1,kkk), fpw(jjj,2,kkk), fpw(jjj,3,kkk), sgg%PlaneWave(jjj)%incert(kkk)
                 // We assume these arrays are global or accessible.
                 print11(0, "%f %f %f %f %f %f %f", pxpw[jjj][kkk], pypw[jjj][kkk], pzpw[jjj][kkk], fpw[jjj][1][kkk], fpw[jjj][2][kkk], fpw[jjj][3][kkk], sgg.PlaneWave[jjj].incert[kkk]);
              }
           }
         }
         goto 635;
      } while (false);
      
634:
      print11(0, SEPARADOR + separador + separador);
      print11(0, "PLANEWAVES: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING");
      print11(0, SEPARADOR + separador + separador);
      
635:
      return;
    } // end subroutine storeplanewaves

    void calc_planewaveconstants(const SGGFDTDINFO_t& sgg, double eps00, double mu00)
    {
       int jjj, kkk;
       eps0 = eps00; mu0 = mu00; // chapuz para convertir la variables de paso en globales
       cluz = 1.0 / sqrt(eps0 * mu0); // lo necesitara incid
       zvac = sqrt(mu0 / eps0); // lo necesitan las variables de mas abajo
!!!!

      for (jjj = 1; jjj <= sgg.numplanewaves; ++jjj) {
        for (kkk = 1; kkk <= sgg.PlaneWave[jjj].nummodes; ++kkk) {
             fpw[jjj][4][kkk] = (pypw[jjj][kkk] * fpw[jjj][3][kkk] - pzpw[jjj][kkk] * fpw[jjj][2][kkk]) / zvac;
             fpw[jjj][5][kkk] = (pzpw[jjj][kkk] * fpw[jjj][1][kkk] - pxpw[jjj][kkk] * fpw[jjj][3][kkk]) / zvac;
             fpw[jjj][6][kkk] = (pxpw[jjj][kkk] * fpw[jjj][2][kkk] - pypw[jjj][kkk] * fpw[jjj][1][kkk]) / zvac;
        }
      }
    } // end subroutine calc_planewaveconstants

    void corrigeondaplanaH(const SGGFDTDINFO_t& sgg, const bounds_t& b, std::vector<std::vector<std::vector<double>>>& Hx, std::vector<std::vector<std::vector<double>>>& Hy, std::vector<std::vector<std::vector<double>>>& Hz, std::vector<std::vector<std::vector<double>>>& Hxvac, std::vector<std::vector<std::vector<double>>>& Hyvac, std::vector<std::vector<std::vector<double>>>& Hzvac)
    {
      //!!!
      // type(SGGFDTDINFO_t), intent(in) :: sgg
      // type( bounds_t), intent( IN) :: b
      // real(kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( INOUT) :: Hx
      // real(kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( INOUT) :: Hy
      // real(kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( INOUT) :: Hz
      // real(kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( INOUT) :: Hxvac
      // real(kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( INOUT) :: Hyvac
      // real(kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( INOUT) :: Hzvac
      // integer(kind=4) :: i, j, k, i_m, j_m, k_m,jjj

      for (int jjj = 1; jjj <= sgg.numplanewaves; ++jjj) {
              if (IluminaTr(jjj)) {
                 //Hz Back
                 i = TrFr[jjj].I.tra.Hz; //Back
                 i_m = i - b.Hz.XI;
                 //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(j, k, j_m, k_m)
#endif
                 for (int k = TrFr[jjj].K.com.Hz; k <= TrFr[jjj].K.fin.Hz; ++k) {
                    int k_m = k - b.Hz.ZI;
                    for (int j = TrFr[jjj].J.com.Hz; j <= TrFr[jjj].J.fin.Hz; ++j) {
                       int j_m = j - b.Hz.YI;
                       //--->
                       Hz(i_m, j_m, k_m) = Hz(i_m, j_m, k_m) - Hzvac(i_m, j_m, k_m);
                    }
                 }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
                 //Hy Back
                 i = TrFr[jjj].I.tra.Hy; //Back
                 i_m = i - b.Hy.XI;
                 //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(j, k, j_m, k_m)
#endif
                 for (int k = TrFr[jjj].K.com.Hy; k <= TrFr[jjj].K.fin.Hy; ++k) {
                    int k_m = k - b.Hy.ZI;
                    for (int j = TrFr[jjj].J.com.Hy; j <= TrFr[jjj].J.fin.Hy; ++j) {
                       int j_m = j - b.Hy.YI;
                       Hy(i_m, j_m, k_m) = Hy(i_m, j_m, k_m) - Hyvac(i_m, j_m, k_m);
                    }
                 }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
              }
              //--->
              if (IluminaFr(jjj)) {
                 //Hz  Front
                 i = TrFr[jjj].I.fro.Hz; //Front
                 i_m = i - b.Hz.XI;
                 //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(j, k, j_m, k_m)
#endif
                 for (int k = TrFr[jjj].K.com.Hz; k <= TrFr[jjj].K.fin.Hz; ++k) {
                    int k_m = k - b.Hz.ZI;
                    for (int j = TrFr[jjj].J.com.Hz; j <= TrFr[jjj].J.fin.Hz; ++j) {
                       int j_m = j - b.Hz.YI;
                       Hz(i_m, j_m, k_m) = Hz(i_m, j_m, k_m) - Hzvac(i_m, j_m, k_m);
                    }
                 }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
                 //Hy  Front
                 i = TrFr[jjj].I.fro.Hy; //Front
                 i_m = i - b.Hy.XI;
                 //--->
#ifdef CompileWithOpenMP
#pragma omp parallel do default(shared) private(j, k, j_m, k_m)
#endif

for (int k = TrFr[jjj].K.com.Hy; k <= TrFr[jjj].K.fin.Hy; ++k) {
                int k_m = k - b.Hy.ZI;
                for (int j = TrFr[jjj].J.com.Hy; j <= TrFr[jjj].J.fin.Hy; ++j) {
                    int j_m = j - b.Hy.YI;
                    Hy[i_m][j_m][k_m] -= Hyvac[i_m][j_m][k_m];
                }
            }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
            } // end if
            //--->
            if (IluminaIz[jjj]) {
                //Hx Left
                int j = IzDe[jjj].J.izq.Hx; //Left
                int j_m = j - b.Hx.YI;
                //--->
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(k, i, k_m, i_m)
#endif
                for (int k = IzDe[jjj].K.com.Hx; k <= IzDe[jjj].K.fin.Hx; ++k) {
                    int k_m = k - b.Hx.ZI;
                    for (int i = IzDe[jjj].I.com.Hx; i <= IzDe[jjj].I.fin.Hx; ++i) {
                        int i_m = i - b.Hx.XI;
                        Hx[i_m][j_m][k_m] -= Hxvac[i_m][j_m][k_m];
                    }
                }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
                //Hz Left
                j = IzDe[jjj].J.izq.Hz; //Left
                j_m = j - b.Hz.YI;
                //--->
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(k, i, k_m, i_m)
#endif
                for (int k = IzDe[jjj].K.com.Hz; k <= IzDe[jjj].K.fin.Hz; ++k) {
                    int k_m = k - b.Hz.ZI;
                    for (int i = IzDe[jjj].I.com.Hz; i <= IzDe[jjj].I.fin.Hz; ++i) {
                        int i_m = i - b.Hz.XI;
                        Hz[i_m][j_m][k_m] -= Hzvac[i_m][j_m][k_m];
                    }
                }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
            } // end if
            //--->
            if (IluminaDe[jjj]) {
                //Hx  Right
                int j = IzDe[jjj].J.der.Hx; //Right
                int j_m = j - b.Hx.YI;
                //--->
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(k, i, k_m, i_m)
#endif
                for (int k = IzDe[jjj].K.com.Hx; k <= IzDe[jjj].K.fin.Hx; ++k) {
                    int k_m = k - b.Hx.ZI;
                    for (int i = IzDe[jjj].I.com.Hx; i <= IzDe[jjj].I.fin.Hx; ++i) {
                        int i_m = i - b.Hx.XI;
                        Hx[i_m][j_m][k_m] -= Hxvac[i_m][j_m][k_m];
                    }
                }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
                //Hz  Right
                j = IzDe[jjj].J.der.Hz; //Right
                j_m = j - b.Hz.YI;
                //--->
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(k, i, k_m, i_m)
#endif
                for (int k = IzDe[jjj].K.com.Hz; k <= IzDe[jjj].K.fin.Hz; ++k) {
                    int k_m = k - b.Hz.ZI;
                    for (int i = IzDe[jjj].I.com.Hz; i <= IzDe[jjj].I.fin.Hz; ++i) {
                        int i_m = i - b.Hz.XI;
                        Hz[i_m][j_m][k_m] -= Hzvac[i_m][j_m][k_m];
                    }
                }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
            } // end if
            //--->
            if (IluminaAb[jjj]) {
                //Hx  Down
                int k = AbAr[jjj].K.aba.Hx; //Down
                int k_m = k - b.Hx.ZI;
                //--->
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, i_m, j_m)
#endif
                for (int j = AbAr[jjj].J.com.Hx; j <= AbAr[jjj].J.fin.Hx; ++j) {
                    int j_m = j - b.Hx.YI;
                    for (int i = AbAr[jjj].I.com.Hx; i <= AbAr[jjj].I.fin.Hx; ++i) {
                        int i_m = i - b.Hx.XI;
                        Hx[i_m][j_m][k_m] -= Hxvac[i_m][j_m][k_m];
                    }
                }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
                //Hy  Down
                k = AbAr[jjj].K.aba.Hy; //Down
                k_m = k - b.Hy.ZI;
                //--->
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, i_m, j_m)
#endif
                for (int j = AbAr[jjj].J.com.Hy; j <= AbAr[jjj].J.fin.Hy; ++j) {
                    int j_m = j - b.Hy.YI;
                    for (int i = AbAr[jjj].I.com.Hy; i <= AbAr[jjj].I.fin.Hy; ++i) {
                        int i_m = i - b.Hy.XI;
                        Hy[i_m][j_m][k_m] -= Hyvac[i_m][j_m][k_m];
                    }
                }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
            } // end if
            //--->
            if (IluminaAr[jjj]) {
                //Hx Up
                int k = AbAr[jjj].K.arr.Hx; //Up
                int k_m = k - b.Hx.ZI;
                //--->
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, i_m, j_m)
#endif
                for (int j = AbAr[jjj].J.com.Hx; j <= AbAr[jjj].J.fin.Hx; ++j) {
                    int j_m = j - b.Hx.YI;
                    for (int i = AbAr[jjj].I.com.Hx; i <= AbAr[jjj].I.fin.Hx; ++i) {
                        int i_m = i - b.Hx.XI;
                        Hx[i_m][j_m][k_m] -= Hxvac[i_m][j_m][k_m];
                    }
                }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
                //Hy Up
                k = AbAr[jjj].K.arr.Hy; //Up
                k_m = k - b.Hy.ZI;
                //--->
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, i_m, j_m)
#endif
                for (int j = AbAr[jjj].J.com.Hy; j <= AbAr[jjj].J.fin.Hy; ++j) {
                    int j_m = j - b.Hy.YI;
                    for (int i = AbAr[jjj].I.com.Hy; i <= AbAr[jjj].I.fin.Hy; ++i) {
                        int i_m = i - b.Hy.XI;
                        Hy[i_m][j_m][k_m] -= Hyvac[i_m][j_m][k_m];
                    }
                }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
            } // end if
        } // end do jjj
        
        return;
    } // end subroutine corrigeondaplanaH

} // end namespace ilumina_m (or namespace corresponding to ilumina_m)