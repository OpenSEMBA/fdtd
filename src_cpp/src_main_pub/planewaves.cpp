#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdexcept>

// Forward declarations and includes for external types/modules
// Assuming FDETYPES_m defines RKIND, BUFSIZE, and basic types
// Assuming Report_m defines stoponerror, WarnErrReport, and coorsxyzP_t
// Assuming these are available in the global namespace or a specific namespace

// Placeholder for types defined in FDETYPES_m and Report_m
// In a real translation, these would be included from headers
namespace FDETYPES_m {
    constexpr int RKIND = 8; // Assuming double precision
    constexpr int BUFSIZE = 256;
}

namespace Report_m {
    // Placeholder for external functions
    void stoponerror(int layoutnumber, int num_procs, const std::string& msg) {
        std::cerr << "ERROR: " << msg << std::endl;
        throw std::runtime_error(msg);
    }

    void WarnErrReport(const std::string& msg) {
        std::cerr << "WARNING: " << msg << std::endl;
    }

    // Placeholder for coorsxyzP_t
    struct coorsxyzP_t {
        // Structure content depends on FDETYPES_m/Report_m
        // Assuming it contains PhysCoor which is an array of coordinate structures
        struct PhysCoor_t {
            std::vector<double> x;
            std::vector<double> y;
            std::vector<double> z;
        };
        std::vector<PhysCoor_t> PhysCoor;
    };
}

// Placeholder for SGGFDTDINFO_t, media_matrices_t, limit_t
// These are complex types likely defined in other modules
struct SGGFDTDINFO_t {
    int NumPlaneWaves;
    int numplanewaves; // Alias or duplicate
    struct {
        int XI, XE, YI, YE, ZI, ZE;
    } Sweep[6]; // Assuming iEx..iHz map to indices 0..5 or similar
    std::vector<double> LineX;
    std::vector<double> Liney;
    std::vector<double> LineZ;
    std::vector<double> Linex; // Alias or duplicate
    std::vector<double> Liney; // Alias or duplicate
    std::vector<double> Linez; // Alias or duplicate
    struct {
        int XI, XE, YI, YE, ZI, ZE;
    } SINPMLSweep[6];
    double dt;
    struct {
        int NumSamples;
        std::vector<double> Samples;
        double deltaSamples;
        std::string Name;
    } Fichero;
    struct {
        int nummodes;
        std::vector<double> px;
        std::vector<double> py;
        std::vector<double> pz;
        std::vector<double> ex;
        std::vector<double> ey;
        std::vector<double> ez;
        std::vector<double> incert;
        bool isRC;
        int esqx1, esqx2, esqy1, esqy2, esqz1, esqz2;
    } PlaneWave;
};

struct media_matrices_t {
    // Placeholder for media data
    int sggMiEz(int i, int j, int k) const { return 1; }
    int sggMiEy(int i, int j, int k) const { return 1; }
    int sggMiEx(int i, int j, int k) const { return 1; }
    struct {
        bool is_pec() const { return false; }
    } med(int idx) const { return {}; }
};

struct limit_t {
    int XI, XE, YI, YE, ZI, ZE;
};

// Enumerations for field indices
enum FieldIndex {
    iEx = 0,
    iEy = 1,
    iEz = 2,
    iHx = 3,
    iHy = 4,
    iHz = 5
};

namespace ilumina_m {

    // Global variables
    double cluz = 0.0;
    double zvac = 0.0;
    double eps0 = 0.0;
    double mu0 = 0.0;

    // Local variables (static-like persistence)
    Report_m::coorsxyzP_t Punto;

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

    void InitPlaneWave(const SGGFDTDINFO_t& sgg, const media_matrices_t& media, int layoutnumber, int num_procs, const std::vector<limit_t>& SINPML_fullsize, bool& ThereArePlaneWaveBoxes, bool resume, double eps00, double mu00) {
        eps0 = eps00;
        mu0 = mu00;
        cluz = 1.0 / std::sqrt(eps0 * mu0);
        zvac = std::sqrt(mu0 / eps0);

        // Allocate Punto coordinates
        // Note: Fortran arrays are 1-based, C++ vectors are 0-based.
        // We will resize vectors to accommodate 1-based indexing if needed, or adjust indices.
        // The code uses explicit indices like sgg%Sweep(field)%XI-1.
        // We assume sgg indices are compatible with vector indexing or we adjust.
        // For simplicity, we'll resize to max index + 1 and ignore index 0 if 1-based.
        
        for (int field = iEx; field <= iHz; ++field) {
            int xi = sgg.Sweep[field].XI - 1;
            int xe = sgg.Sweep[field].XE + 1;
            int yi = sgg.Sweep[field].YI - 1;
            int ye = sgg.Sweep[field].YE + 1;
            int zi = sgg.Sweep[field].ZI - 1;
            int ze = sgg.Sweep[field].ZE + 1;

            // Ensure vectors are large enough
            if (Punto.PhysCoor.size() <= field) {
                Punto.PhysCoor.resize(field + 1);
            }
            
            Punto.PhysCoor[field].x.resize(xe + 1);
            Punto.PhysCoor[field].y.resize(ye + 1);
            Punto.PhysCoor[field].z.resize(ze + 1);
        }

        // Fill coordinates for each field component
        // iEx
        {
            int field = iEx;
            int xi = sgg.Sweep[field].XI - 1;
            int xe = sgg.Sweep[field].XE;
            for (int i = xi; i <= xe; ++i) {
                Punto.PhysCoor[field].x[i] = 0.5 * (sgg.LineX[i] + sgg.LineX[i + 1]);
            }
            int yi = sgg.Sweep[field].YI - 1;
            int ye = sgg.Sweep[field].YE + 1;
            for (int j = yi; j <= ye; ++j) {
                Punto.PhysCoor[field].y[j] = sgg.Liney[j];
            }
            int zi = sgg.Sweep[field].ZI - 1;
            int ze = sgg.Sweep[field].ZE + 1;
            for (int k = zi; k <= ze; ++k) {
                Punto.PhysCoor[field].z[k] = sgg.LineZ[k];
            }
        }

        // iEy
        {
            int field = iEy;
            int xi = sgg.Sweep[field].XI - 1;
            int xe = sgg.Sweep[field].XE + 1;
            for (int i = xi; i <= xe; ++i) {
                Punto.PhysCoor[field].x[i] = sgg.LineX[i];
            }
            int yi = sgg.Sweep[field].YI - 1;
            int ye = sgg.Sweep[field].YE;
            for (int j = yi; j <= ye; ++j) {
                Punto.PhysCoor[field].y[j] = 0.5 * (sgg.Liney[j] + sgg.Liney[j + 1]);
            }
            int zi = sgg.Sweep[field].ZI - 1;
            int ze = sgg.Sweep[field].ZE + 1;
            for (int k = zi; k <= ze; ++k) {
                Punto.PhysCoor[field].z[k] = sgg.LineZ[k];
            }
        }

        // iEz
        {
            int field = iEz;
            int xi = sgg.Sweep[field].XI - 1;
            int xe = sgg.Sweep[field].XE + 1;
            for (int i = xi; i <= xe; ++i) {
                Punto.PhysCoor[field].x[i] = sgg.LineX[i];
            }
            int yi = sgg.Sweep[field].YI - 1;
            int ye = sgg.Sweep[field].YE + 1;
            for (int j = yi; j <= ye; ++j) {
                Punto.PhysCoor[field].y[j] = sgg.Liney[j];
            }
            int zi = sgg.Sweep[field].ZI - 1;
            int ze = sgg.Sweep[field].ZE;
            for (int k = zi; k <= ze; ++k) {
                Punto.PhysCoor[field].z[k] = 0.5 * (sgg.LineZ[k] + sgg.LineZ[k + 1]);
            }
        }

        // iHx
        {
            int field = iHx;
            int xi = sgg.Sweep[field].XI - 1;
            int xe = sgg.Sweep[field].XE + 1;
            for (int i = xi; i <= xe; ++i) {
                Punto.PhysCoor[field].x[i] = sgg.LineX[i];
            }
            int yi = sgg.Sweep[field].YI - 1;
            int ye = sgg.Sweep[field].YE;
            for (int j = yi; j <= ye; ++j) {
                Punto.PhysCoor[field].y[j] = 0.5 * (sgg.Liney[j] + sgg.Liney[j + 1]);
            }
            int zi = sgg.Sweep[field].ZI - 1;
            int ze = sgg.Sweep[field].ZE;
            for (int k = zi; k <= ze; ++k) {
                Punto.PhysCoor[field].z[k] = 0.5 * (sgg.LineZ[k] + sgg.LineZ[k + 1]);
            }
        }

        // iHy
        {
            int field = iHy;
            int xi = sgg.Sweep[field].XI - 1;
            int xe = sgg.Sweep[field].XE;
            for (int i = xi; i <= xe; ++i) {
                Punto.PhysCoor[field].x[i] = 0.5 * (sgg.LineX[i] + sgg.LineX[i + 1]);
            }
            int yi = sgg.Sweep[field].YI - 1;
            int ye = sgg.Sweep[field].YE + 1;
            for (int j = yi; j <= ye; ++j) {
                Punto.PhysCoor[field].y[j] = sgg.Liney[j];
            }
            int zi = sgg.Sweep[field].ZI - 1;
            int ze = sgg.Sweep[field].ZE;
            for (int k = zi; k <= ze; ++k) {
                Punto.PhysCoor[field].z[k] = 0.5 * (sgg.LineZ[k] + sgg.LineZ[k + 1]);
            }
        }

        // iHz
        {
            int field = iHz;
            int xi = sgg.Sweep[field].XI - 1;
            int xe = sgg.Sweep[field].XE;
            for (int i = xi; i <= xe; ++i) {
                Punto.PhysCoor[field].x[i] = 0.5 * (sgg.LineX[i] + sgg.LineX[i + 1]);
            }
            int yi = sgg.Sweep[field].YI - 1;
            int ye = sgg.Sweep[field].YE;
            for (int j = yi; j <= ye; ++j) {
                Punto.PhysCoor[field].y[j] = 0.5 * (sgg.Liney[j] + sgg.Liney[j + 1]);
            }
            int zi = sgg.Sweep[field].ZI - 1;
            int ze = sgg.Sweep[field].ZE + 1;
            for (int k = zi; k <= ze; ++k) {
                Punto.PhysCoor[field].z[k] = sgg.LineZ[k];
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
            int n_pw = sgg.numplanewaves;
            TrFr.resize(n_pw + 1);
            IzDe.resize(n_pw + 1);
            AbAr.resize(n_pw + 1);
            IluminaTr.resize(n_pw + 1, false);
            IluminaFr.resize(n_pw + 1, false);
            IluminaIz.resize(n_pw + 1, false);
            IluminaDe.resize(n_pw + 1, false);
            IluminaAr.resize(n_pw + 1, false);
            IluminaAb.resize(n_pw + 1, false);
            numus.resize(n_pw + 1);
            deltaevol.resize(n_pw + 1);

            for (int jjj = 1; jjj <= sgg.NumPlaneWaves; ++jjj) {
                // Note: Accessing sgg.PlaneWave(jjj) implies sgg.PlaneWave is an array/vector
                // Assuming sgg.PlaneWave is 1-based or we adjust index
                // For this translation, assuming sgg.PlaneWave is a vector and we access jjj-1 if 0-based
                // However, Fortran code uses jjj directly. We will assume sgg.PlaneWave is resized to n_pw+1
                // or we adjust access. Let's assume sgg.PlaneWave is a vector and we use jjj-1 for 0-based access
                // BUT the code uses sgg.PlaneWave(jjj). If it's a vector, it should be sgg.PlaneWave[jjj-1].
                // To preserve logic, let's assume the structure supports 1-based indexing or we map it.
                // Given the complexity, I will assume sgg.PlaneWave is a vector and use jjj-1.
                
                // Correction: The code uses sgg.PlaneWave(jjj). If this is a derived type array, 
                // in C++ it's likely a std::vector. Fortran 1-based indexing means jjj=1 is index 0 in C++.
                // I will use jjj-1 for access to sgg.PlaneWave.
                
                int idx = jjj - 1;
                numus[jjj] = sgg.PlaneWave[idx].NumSamples;

                bool abortar = 
                    (sgg.PlaneWave[idx].esqx1 <= SINPML_fullsize[iHx].XI) &&
                    (sgg.PlaneWave[idx].esqx2 >= SINPML_fullsize[iHx].XE) &&
                    (sgg.PlaneWave[idx].esqy1 <= SINPML_fullsize[iHy].YI) &&
                    (sgg.PlaneWave[idx].esqy2 >= SINPML_fullsize[iHy].YE) &&
                    (sgg.PlaneWave[idx].esqz1 <= SINPML_fullsize[iHz].ZI) &&
                    (sgg.PlaneWave[idx].esqz2 >= SINPML_fullsize[iHz].ZE);

                if (abortar) {
                    std::string buff = "At least one of TF/SF planes must be 1 cell inside the simulation region. Aborting";
                    Report_m::stoponerror(layoutnumber, num_procs, buff);
                }

                IluminaTr[jjj] = false;
                IluminaFr[jjj] = false;
                IluminaIz[jjj] = false;
                IluminaDe[jjj] = false;
                IluminaAr[jjj] = false;
                IluminaAb[jjj] = false;

                if ((sgg.PlaneWave[idx].esqx1 >= sgg.SINPMLSweep[iHx].XI) && (sgg.PlaneWave[idx].esqx1 <= sgg.SINPMLSweep[iHx].XE)) {
                    IluminaTr[jjj] = true;
                }
                if ((sgg.PlaneWave[idx].esqx2 <= sgg.SINPMLSweep[iHx].XE) && (sgg.PlaneWave[idx].esqx2 >= sgg.SINPMLSweep[iHx].XI)) {
                    IluminaFr[jjj] = true;
                }
                if ((sgg.PlaneWave[idx].esqy1 >= sgg.SINPMLSweep[iHy].YI) && (sgg.PlaneWave[idx].esqy1 <= sgg.SINPMLSweep[iHy].YE)) {
                    IluminaIz[jjj] = true;
                }
                if ((sgg.PlaneWave[idx].esqy2 <= sgg.SINPMLSweep[iHy].YE) && (sgg.PlaneWave[idx].esqy2 >= sgg.SINPMLSweep[iHy].YI)) {
                    IluminaDe[jjj] = true;
                }
                if ((sgg.PlaneWave[idx].esqz1 >= sgg.SINPMLSweep[iHz].ZI) && (sgg.PlaneWave[idx].esqz1 <= sgg.SINPMLSweep[iHz].ZE)) {
                    IluminaAb[jjj] = true;
                }
                if ((sgg.PlaneWave[idx].esqz2 <= sgg.SINPMLSweep[iHz].ZE) && (sgg.PlaneWave[idx].esqz2 >= sgg.SINPMLSweep[iHz].ZI)) {
                    IluminaAr[jjj] = true;
                }

                // Find coordinate limits
                TrFr[jjj].i.tra.Ez = std::max(sgg.SINPMLSweep[iEz].XI, sgg.PlaneWave[idx].esqx1);
                TrFr[jjj].i.fro.Ez = std::min(sgg.SINPMLSweep[iEz].XE, sgg.PlaneWave[idx].esqx2);
                TrFr[jjj].j.com.Ez = std::max(sgg.SINPMLSweep[iEz].YI, sgg.PlaneWave[idx].esqy1);
                TrFr[jjj].j.fin.Ez = std::min(sgg.SINPMLSweep[iEz].YE, sgg.PlaneWave[idx].esqy2);
                TrFr[jjj].k.com.Ez = std::max(sgg.SINPMLSweep[iEz].ZI, sgg.PlaneWave[idx].esqz1);
                TrFr[jjj].k.fin.Ez = std::min(sgg.SINPMLSweep[iEz].ZE, sgg.PlaneWave[idx].esqz2 - 1);

                TrFr[jjj].i.tra.Ey = std::max(sgg.SINPMLSweep[iEy].XI, sgg.PlaneWave[idx].esqx1);
                TrFr[jjj].i.fro.Ey = std::min(sgg.SINPMLSweep[iEy].XE, sgg.PlaneWave[idx].esqx2);
                TrFr[jjj].j.com.Ey = std::max(sgg.SINPMLSweep[iEy].YI, sgg.PlaneWave[idx].esqy1);
                TrFr[jjj].j.fin.Ey = std::min(sgg.SINPMLSweep[iEy].YE, sgg.PlaneWave[idx].esqy2 - 1);
                TrFr[jjj].k.com.Ey = std::max(sgg.SINPMLSweep[iEy].ZI, sgg.PlaneWave[idx].esqz1);
                TrFr[jjj].k.fin.Ey = std::min(sgg.SINPMLSweep[iEy].ZE, sgg.PlaneWave[idx].esqz2);

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

                IzDe[jjj].j.izq.Ex = std::max(sgg.SINPMLSweep[iEx].YI, sgg.PlaneWave[idx].esqy1);
                IzDe[jjj].j.der.Ex = std::min(sgg.SINPMLSweep[iEx].YE, sgg.PlaneWave[idx].esqy2);
                IzDe[jjj].i.com.Ex = std::max(sgg.SINPMLSweep[iEx].XI, sgg.PlaneWave[idx].esqx1);
                IzDe[jjj].i.fin.Ex = std::min(sgg.SINPMLSweep[iEx].XE, sgg.PlaneWave[idx].esqx2 - 1);
                IzDe[jjj].k.com.Ex = std::max(sgg.SINPMLSweep[iEx].ZI, sgg.PlaneWave[idx].esqz1);
                IzDe[jjj].k.fin.Ex = std::min(sgg.SINPMLSweep[iEx].ZE, sgg.PlaneWave[idx].esqz2);

                IzDe[jjj].j.izq.Ez = std::max(sgg.SINPMLSweep[iEz].YI, sgg.PlaneWave[idx].esqy1);
                IzDe[jjj].j.der.Ez = std::min(sgg.SINPMLSweep[iEz].YE, sgg.PlaneWave[idx].esqy2);
                IzDe[jjj].i.com.Ez = std::max(sgg.SINPMLSweep[iEz].XI, sgg.PlaneWave[idx].esqx1);
                IzDe[jjj].i.fin.Ez = std::min(sgg.SINPMLSweep[iEz].XE, sgg.PlaneWave[idx].esqx2);
                IzDe[jjj].k.com.Ez = std::max(sgg.SINPMLSweep[iEz].ZI, sgg.PlaneWave[idx].esqz1);
                IzDe[jjj].k.fin.Ez = std::min(sgg.SINPMLSweep[iEz].ZE, sgg.PlaneWave[idx].esqz2 - 1);

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

                AbAr[jjj].k.aba.Ey = std::max(sgg.SINPMLSweep[iEy].ZI, sgg.PlaneWave[idx].esqz1);
                AbAr[jjj].k.arr.Ey = std::min(sgg.SINPMLSweep[iEy].ZE, sgg.PlaneWave[idx].esqz2);
                AbAr[jjj].i.com.Ey = std::max(sgg.SINPMLSweep[iEy].XI, sgg.PlaneWave[idx].esqx1);
                AbAr[jjj].i.fin.Ey = std::min(sgg.SINPMLSweep[iEy].XE, sgg.PlaneWave[idx].esqx2);
                AbAr[jjj].j.com.Ey = std::max(sgg.SINPMLSweep[iEy].YI, sgg.PlaneWave[idx].esqy1);
                AbAr[jjj].j.fin.Ey = std::min(sgg.SINPMLSweep[iEy].YE, sgg.PlaneWave[idx].esqy2 - 1);

                AbAr[jjj].k.aba.Ex = std::max(sgg.SINPMLSweep[iEx].ZI, sgg.PlaneWave[idx].esqz1);
                AbAr[jjj].k.arr.Ex = std::min(sgg.SINPMLSweep[iEx].ZE, sgg.PlaneWave[idx].esqz2);
                AbAr[jjj].i.com.Ex = std::max(sgg.SINPMLSweep[iEx].XI, sgg.PlaneWave[idx].esqx1);
                AbAr[jjj].i.fin.Ex = std::min(sgg.SINPMLSweep[iEx].XE, sgg.PlaneWave[idx].esqx2 - 1);
                AbAr[jjj].j.com.Ex = std::max(sgg.SINPMLSweep[iEx].YI, sgg.PlaneWave[idx].esqy1);
                AbAr[jjj].j.fin.Ex = std::min(sgg.SINPMLSweep[iEx].YE, sgg.PlaneWave[idx].esqy2);

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
        for (int jjj = 1; jjj <= sgg.numplanewaves; ++jjj) {
            if (numus[jjj] > maxnumus) maxnumus = numus[jjj];
        }

        evol.resize(sgg.numplanewaves + 1, std::vector<double>(maxnumus + 1, 0.0));
        for (int jjj = 1; jjj <= sgg.numplanewaves; ++jjj) {
            int idx = jjj - 1;
            for (int k = 0; k <= numus[jjj]; ++k) {
                evol[jjj][k] = sgg.PlaneWave[idx].Samples[k];
            }
            deltaevol[jjj] = sgg.PlaneWave[idx].deltaSamples;
            if (deltaevol[jjj] > sgg.dt) {
                std::string buff = "WARNING: " + sgg.PlaneWave[idx].Name + " undersampled by a factor " + std::to_string(deltaevol[jjj] / sgg.dt);
                Report_m::WarnErrReport(buff);
            }
        }

        int maxmodes = 0;
        for (int jjj = 1; jjj <= sgg.numplanewaves; ++jjj) {
            int idx = jjj - 1;
            if (sgg.PlaneWave[idx].nummodes > maxmodes) maxmodes = sgg.PlaneWave[idx].nummodes;
        }

        pxpw.resize(sgg.numplanewaves + 1, std::vector<double>(maxmodes + 1, 0.0));
        pypw.resize(sgg.numplanewaves + 1, std::vector<double>(maxmodes + 1, 0.0));
        pzpw.resize(sgg.numplanewaves + 1, std::vector<double>(maxmodes + 1, 0.0));
        fpw.resize(sgg.numplanewaves + 1, std::vector<std::vector<double>>(7, std::vector<double>(maxmodes + 1, 0.0))); // 1:6 indices
        INCERT.resize(sgg.numplanewaves + 1, std::vector<double>(maxmodes + 1, 0.0));
        distanciaInicial.resize(sgg.numplanewaves + 1, std::vector<double>(maxmodes + 1, 0.0));

        for (int jjj = 1; jjj <= sgg.numplanewaves; ++jjj) {
            int idx = jjj - 1;
            for (int kkk = 1; kkk <= sgg.PlaneWave[idx].nummodes; ++kkk) {
                if (!resume) {
                    pxpw[jjj][kkk] = sgg.PlaneWave[idx].px[kkk-1]; // Assuming 0-based vector
                    pypw[jjj][kkk] = sgg.PlaneWave[idx].py[kkk-1];
                    pzpw[jjj][kkk] = sgg.PlaneWave[idx].pz[kkk-1];
                    fpw[jjj][1][kkk] = sgg.PlaneWave[idx].ex[kkk-1];
                    fpw[jjj][2][kkk] = sgg.PlaneWave[idx].ey[kkk-1];
                    fpw[jjj][3][kkk] = sgg.PlaneWave[idx].ez[kkk-1];

                    double modulus = std::sqrt(pxpw[jjj][kkk] * pxpw[jjj][kkk] + pypw[jjj][kkk] * pypw[jjj][kkk] + pzpw[jjj][kkk] * pzpw[jjj][kkk]);
                    pxpw[jjj][kkk] /= modulus;
                    pypw[jjj][kkk] /= modulus;
                    pzpw[jjj][kkk] /= modulus;
                    INCERT[jjj][kkk] = sgg.PlaneWave[idx].incert[kkk-1];
                } else {
                    if (sgg.PlaneWave[idx].isRC) {
                        // Read from file 14
                        // In C++, this would require file I/O. Placeholder.
                        // double px, py, pz, ex, ey, ez, inc;
                        // std::cin >> px >> py >> pz >> ex >> ey >> ez >> inc; // Placeholder
                        // pxpw[jjj][kkk] = px; ...
                    } else {
                        pxpw[jjj][kkk] = sgg.PlaneWave[idx].px[kkk-1];
                        pypw[jjj][kkk] = sgg.PlaneWave[idx].py[kkk-1];
                        pzpw[jjj][kkk] = sgg.PlaneWave[idx].pz[kkk-1];
                        fpw[jjj][1][kkk] = sgg.PlaneWave[idx].ex[kkk-1];
                        fpw[jjj][2][kkk] = sgg.PlaneWave[idx].ey[kkk-1];
                        fpw[jjj][3][kkk] = sgg.PlaneWave[idx].ez[kkk-1];

                        double modulus = std::sqrt(pxpw[jjj][kkk] * pxpw[jjj][kkk] + pypw[jjj][kkk] * pypw[jjj][kkk] + pzpw[jjj][kkk] * pzpw[jjj][kkk]);
                        pxpw[jjj][kkk] /= modulus;
                        pypw[jjj][kkk] /= modulus;
                        pzpw[jjj][kkk] /= modulus;
                        INCERT[jjj][kkk] = sgg.PlaneWave[idx].incert[kkk-1];
                    }
                }
            }
        }

        for (int jjj = 1; jjj <= sgg.numplanewaves; ++jjj) {
            int idx = jjj - 1;
            for (int kkk = 1; kkk <= sgg.PlaneWave[idx].nummodes; ++kkk) {
                double XD0, YD0, ZD0;
                if ((pxpw[jjj][kkk] >= 0) && (pypw[jjj][kkk] >= 0) && (pzpw[jjj][kkk] >= 0)) {
                    XD0 = sgg.Linex[std::max(sgg.PlaneWave[idx].esqx1 - 1, SINPML_fullsize[iHx].XI)];
                    YD0 = sgg.Liney[std::max(sgg.PlaneWave[idx].esqy1 - 1, SINPML_fullsize[iHy].YI)];
                    ZD0 = sgg.Linez[std::max(sgg.PlaneWave[idx].esqz1 - 1, SINPML_fullsize[iHz].ZI)];
                } else if ((pxpw[jjj][kkk] >= 0) && (pypw[jjj][kkk] >= 0) && (pzpw[jjj][kkk] < 0)) {
                    XD0 = sgg.Linex[std::max(sgg.PlaneWave[idx].esqx1 - 1, SINPML_fullsize[iHx].XI)];
                    YD0 = sgg.Liney[std::max(sgg.PlaneWave[idx].esqy1 - 1, SINPML_fullsize[iHy].YI)];
                    ZD0 = sgg.Linez[std::min(sgg.PlaneWave[idx].esqz2 + 1, SINPML_fullsize[iHz].ZE)];
                } else if ((pxpw[jjj][kkk] >= 0) && (pypw[jjj][kkk] < 0) && (pzpw[jjj][kkk] >= 0)) {
                    XD0 = sgg.Linex[std::max(sgg.PlaneWave[idx].esqx1 - 1, SINPML_fullsize[iHx].XI)];
                    YD0 = sgg.Liney[std::min(sgg.PlaneWave[idx].esqy2 + 1, SINPML_fullsize[iHy].YE)];
                    ZD0 = sgg.Linez[std::max(sgg.PlaneWave[idx].esqz1 - 1, SINPML_fullsize[iHz].ZI)];
                } else if ((pxpw[jjj][kkk] < 0) && (pypw[jjj][kkk] >= 0) && (pzpw[jjj][kkk] >= 0)) {
                    XD0 = sgg.Linex[std::min(sgg.PlaneWave[idx].esqx2 + 1, SINPML_fullsize[iHx].XE)];
                    YD0 = sgg.Liney[std::max(sgg.PlaneWave[idx].esqy1 - 1, SINPML_fullsize[iHy].YI)];
                    ZD0 = sgg.Linez[std::max(sgg.PlaneWave[idx].esqz1 - 1, SINPML_fullsize[iHz].ZI)];
                } else if ((pxpw[jjj][kkk] >= 0) && (pypw[jjj][kkk] < 0) && (pzpw[jjj][kkk] < 0)) {
                    XD0 = sgg.Linex[std::max(sgg.PlaneWave[idx].esqx1 - 1, SINPML_fullsize[iHx].XI)];
                    YD0 = sgg.Liney[std::min(sgg.PlaneWave[idx].esqy2 + 1, SINPML_fullsize[iHy].YE)];
                    ZD0 = sgg.Linez[std::min(sgg.PlaneWave[idx].esqz2 + 1, SINPML_fullsize[iHz].ZE)];
                } else if ((pxpw[jjj][kkk] < 0) && (pypw[jjj][kkk] < 0) && (pzpw[jjj][kkk] >= 0)) {
                    XD0 = sgg.Linex[std::min(sgg.PlaneWave[idx].esqx2 + 1, SINPML_fullsize[iHx].XE)];
                    YD0 = sgg.Liney[std::min(sgg.PlaneWave[idx].esqy2 + 1, SINPML_fullsize[iHy].YE)];
                    ZD0 = sgg.Linez[std::max(sgg.PlaneWave[idx].esqz1 - 1, SINPML_fullsize[iHz].ZI)];
                } else if ((pxpw[jjj][kkk] < 0) && (pypw[jjj][kkk] >= 0) && (pzpw[jjj][kkk] < 0)) {
                    XD0 = sgg.Linex[std::min(sgg.PlaneWave[idx].esqx2 + 1, SINPML_fullsize[iHx].XE)];
                    YD0 = sgg.Liney[std::max(sgg.PlaneWave[idx].esqy1 - 1, SINPML_fullsize[iHy].YI)];
                    ZD0 = sgg.Linez[std::min(sgg.PlaneWave[idx].esqz2 + 1, SINPML_fullsize[iHz].ZE)];
                } else if ((pxpw[jjj][kkk] < 0) && (pypw[jjj][kkk] < 0) && (pzpw[jjj][kkk] < 0)) {
                    XD0 = sgg.Linex[std::min(sgg.PlaneWave[idx].esqx2 + 1, SINPML_fullsize[iHx].XE)];
                    YD0 = sgg.Liney[std::min(sgg.PlaneWave[idx].esqy2 + 1, SINPML_fullsize[iHy].YE)];
                    ZD0 = sgg.Linez[std::min(sgg.PlaneWave[idx].esqz2 + 1, SINPML_fullsize[iHz].ZE)];
                } else {
                    Report_m::stoponerror(layoutnumber, num_procs, "buggy xo,yo,z0");
                }

                double diagonalcaja = std::sqrt(
                    std::pow(sgg.Linex[std::max(sgg.PlaneWave[idx].esqx1 - 1, SINPML_fullsize[iHx].XI)] - sgg.Linex[std::min(sgg.PlaneWave[idx].esqx2 + 1, SINPML_fullsize[iHx].XE)], 2) +
                    std::pow(sgg.Liney[std::max(sgg.PlaneWave[idx].esqy1 - 1, SINPML_fullsize[iHy].YI)] - sgg.Liney[std::min(sgg.PlaneWave[idx].esqy2 + 1, SINPML_fullsize[iHy].YE)], 2) +
                    std::pow(sgg.Linez[std::max(sgg.PlaneWave[idx].esqz1 - 1, SINPML_fullsize[iHz].ZI)] - sgg.Linez[std::min(sgg.PlaneWave[idx].esqz2 + 1, SINPML_fullsize[iHz].ZE)], 2)
                );

                distanciaInicial[jjj][kkk] = (XD0 * pxpw[jjj][kkk] + YD0 * pypw[jjj][kkk] + ZD0 * pzpw[jjj][kkk]) - INCERT[jjj][kkk] * diagonalcaja;
            }
        }

        // Check materials
        for (int jjj = 1; jjj <= sgg.numplanewaves; ++jjj) {
            if (IluminaTr[jjj]) {
                // Ez Back
                int i = TrFr[jjj].i.tra.Ez;
                for (int k = TrFr[jjj].k.com.Ez; k <= TrFr[jjj].k.fin.Ez; ++k) {
                    for (int j = TrFr[jjj].j.com.Ez; j <= TrFr[jjj].j.fin.Ez; ++j) {
                        if (media.sggMiEz(i, j, k) != 1) {
                            std::string buff = "Back TF/SF region intersects a material at Ez " + std::to_string(i) + " " + std::to_string(j) + " " + std::to_string(k);
                            if (((media.sggMiEz(i, j, k) == 0) || (media.med(media.sggMiEz(i, j, k)).is_pec())) && 
                                !((i == sgg.SINPMLSweep[iHx].XI) || (j == sgg.SINPMLSweep[iHy].YI) || (k == sgg.SINPMLSweep[iHz].ZI) || 
                                  (i == sgg.SINPMLSweep[iHx].XE) || (j == sgg.SINPMLSweep[iHy].YE) || (k == sgg.SINPMLSweep[iHz].ZE))) {
                                Report_m::stoponerror(layoutnumber, num_procs, buff);
                            }
                        }
                    }
                }
                // Ey Back
                i = TrFr[jjj].i.tra.Ey;
                for (int k = TrFr[jjj].k.com.Ey; k <= TrFr[jjj].k.fin.Ey; ++k) {
                    for (int j = TrFr[jjj].j.com.Ey; j <= TrFr[jjj].j.fin.Ey; ++j) {
                        if (media.sggMiEy(i, j, k) != 1) {
                            std::string buff = "Back TF/SF region intersects a material at Ey " + std::to_string(i) + " " + std::to_string(j) + " " + std::to_string(k);
                            if (((media.sggMiEy(i, j, k) == 0) || (media.med(media.sggMiEy(i, j, k)).is_pec())) && 
                                !((i == sgg.SINPMLSweep[iHx].XI) || (j == sgg.SINPMLSweep[iHy].YI) || (k == sgg.SINPMLSweep[iHz].ZI) || 
                                  (i == sgg.SINPMLSweep[iHx].XE) || (j == sgg.SINPMLSweep[iHy].YE) || (k == sgg.SINPMLSweep[iHz].ZE))) {
                                Report_m::stoponerror(layoutnumber, num_procs, buff);
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
                        if (media.sggMiEz(i, j, k) != 1) {
                            std::string buff = "Front TF/SF region intersects a material at Ez " + std::to_string(i) + " " + std::to_string(j) + " " + std::to_string(k);
                            if (((media.sggMiEz(i, j, k) == 0) || (media.med(media.sggMiEz(i, j, k)).is_pec())) && 
                                !((i == sgg.SINPMLSweep[iHx].XI) || (j == sgg.SINPMLSweep[iHy].YI) || (k == sgg.SINPMLSweep[iHz].ZI) || 
                                  (i == sgg.SINPMLSweep[iHx].XE) || (j == sgg.SINPMLSweep[iHy].YE) || (k == sgg.SINPMLSweep[iHz].ZE))) {
                                Report_m::stoponerror(layoutnumber, num_procs, buff);
                            }
                        }
                    }
                }
                // Ey Front
                i = TrFr[jjj].i.fro.Ey;
                for (int k = TrFr[jjj].k.com.Ey; k <= TrFr[jjj].k.fin.Ey; ++k) {
                    for (int j = TrFr[jjj].j.com.Ey; j <= TrFr[jjj].j.fin.Ey; ++j) {
                        if (media.sggMiEy(i, j, k) != 1) {
                            std::string buff = "Front TF/SF region intersects a material at Ey " + std::to_string(i) + " " + std::to_string(j) + " " + std::to_string(k);
                            if (((media.sggMiEy(i, j, k) == 0) || (media.med(media.sggMiEy(i, j, k)).is_pec())) && 
                                !((i == sgg.SINPMLSweep[iHx].XI) || (j == sgg.SINPMLSweep[iHy].YI) || (k == sgg.SINPMLSweep[iHz].ZI) || 
                                  (i == sgg.SINPMLSweep[iHx].XE) || (j == sgg.SINPMLSweep[iHy].YE) || (k == sgg.SINPMLSweep[iHz].ZE))) {
                                Report_m::stoponerror(layoutnumber, num_procs, buff);
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
                        if (media.sggMiEx(i, j, k) != 1) {
                            std::string buff = "Left TF/SF region intersects a material at Ex " + std::to_string(i) + " " + std::to_string(j) + " " + std::to_string(k);
                            if (((media.sggMiEx(i, j, k) == 0) || (media.med(media.sggMiEx(i, j, k)).is_pec())) && 
                                !((i == sgg.SINPMLSweep[iHx].XI) || (j == sgg.SINPMLSweep[iHy].YI) || (k == sgg.SINPMLSweep[iHz].ZI) || 
                                  (i == sgg.SINPMLSweep[iHx].XE) || (j == sgg.SINPMLSweep[iHy].YE) || (k == sgg.SINPMLSweep[iHz].ZE))) {
                                Report_m::stoponerror(layoutnumber, num_procs, buff);
                            }
                        }
                    }
                }
                // Ez Left
                j = IzDe[jjj].j.izq.Ez;
                // ... (Code cuts off here in the prompt, so we stop here)
            }
        }
    }
}

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdexcept>

// Forward declarations and includes for external types/functions used in this chunk
// These would typically be in a header file generated from previous chunks
struct SGGFDTDINFO_t;
struct bounds_t;
struct PhysCoor_t;

extern PhysCoor_t Punto;
extern std::vector<std::vector<double>> pxpw;
extern std::vector<std::vector<double>> pypw;
extern std::vector<std::vector<double>> pzpw;
extern std::vector<std::vector<std::vector<double>>> fpw;
extern std::vector<double> distanciaInicial;
extern std::vector<double> INCERT;
extern std::vector<int> numus;
extern std::vector<double> deltaevol;
extern std::vector<std::vector<double>> evol;

// Assuming these are defined in previous chunks or headers
// TrFr, IzDe, AbAr, IluminaTr, etc. are likely global or member arrays/vectors
extern std::vector<PlaneWaveInfo_t> TrFr; 
extern std::vector<PlaneWaveInfo_t> IzDe;
extern std::vector<PlaneWaveInfo_t> AbAr;
extern std::vector<bool> IluminaTr;
extern std::vector<bool> IluminaFr;
extern std::vector<bool> IluminaIz;
extern std::vector<bool> IluminaDe;
extern std::vector<bool> IluminaAr;
extern std::vector<bool> IluminaAb;

// Helper function declarations
void stoponerror(int layoutnumber, int num_procs, const std::string& buff);
double Incid(const SGGFDTDINFO_t& sgg, int jjj, int nfield, double time, int i, int j, int k, bool& still_planewave_time, bool calledfromobservation);
void calc_planewaveconstants(const SGGFDTDINFO_t& sgg, double eps0, double mu0);
void DestroyIlumina(SGGFDTDINFO_t& sgg);

// Constants
const int iEx = 1; // Assuming indices based on context, usually 1,2,3 for Ex,Ey,Ez
const int iEy = 2;
const int iHz = 3;
const int iHx = 1;
const int iHy = 2;
const int iHz_H = 3; // Avoid conflict if iHz is used for index and field
// Note: The Fortran code uses iEx, iEy, iHz, iHx, iHy, iHz as integer constants for field indices.
// We assume they are defined elsewhere or define them here if not present.
// Based on typical FDTD: Ex=1, Ey=2, Ez=3, Hx=1, Hy=2, Hz=3? Or distinct?
// The code uses iHy in Incid call for Ez update, iHz for Ey update.
// Let's assume standard mapping or that these are global constants.
// For translation safety, we assume they are available in scope or defined as constexpr.
constexpr int iEx_const = 1;
constexpr int iEy_const = 2;
constexpr int iEz_const = 3;
constexpr int iHx_const = 1;
constexpr int iHy_const = 2;
constexpr int iHz_const = 3;

// RKIND is typically double
using RKIND = double;
const RKIND RKIND_VAL = 1.0; // Marker for kind

// cluz is likely speed of light or similar constant
extern double cluz;

// BUFSIZE for character strings
constexpr int BUFSIZE = 1024;

// Placeholder for PlaneWaveInfo_t structure members accessed in the code
// This structure would have been defined in previous chunks
struct PlaneWaveInfo_t {
    struct {
        struct {
            int com;
            int fin;
            int der;
            int izq;
            int tra;
            int fro;
            int aba;
            int arr;
        } I, J, K;
    } %; // Simulating Fortran % access
    int nummodes;
};

// Note: The Fortran code uses % for structure access. In C++, we use .
// The struct PlaneWaveInfo_t above is a placeholder. The actual struct should match the Fortran derived type.

void InitPlaneWave_impl(int numplanewaves, int layoutnumber, int num_procs, const SGGFDTDINFO_t& sgg, double eps0, double mu0) {
    // This function encapsulates the logic from the first part of the chunk
    // Note: The Fortran code is inside a subroutine InitPlaneWave. 
    // We are translating the body.
    
    // Variables from context
    // media, sgg, IzDe, IluminaDe, etc. are assumed to be accessible
    
    // Loop over planewaves
    for (int jjj = 1; jjj <= numplanewaves; ++jjj) {
        // Left TF/SF region intersects a material at Ez
        if (IluminaIz[jjj]) { // Assuming IluminaIz is 1-based or adjusted
            int j = IzDe[jjj].J.izq.Ez; // Left
            for (int k = IzDe[jjj].K.com.Ez; k <= IzDe[jjj].K.fin.Ez; ++k) {
                for (int i = IzDe[jjj].I.com.Ez; i <= IzDe[jjj].I.fin.Ez; ++i) {
                    if (media.sggMiEz(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Left TF/SF region intersects a material at Ez %7d%7d%7d", i, j, k);
                        if (((media.sggMiEz(i, j, k) == 0) || (sgg.med(media.sggMiEz(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
        }
        
        // Right TF/SF region intersects a material at Ez
        if (IluminaDe[jjj]) {
            int j = IzDe[jjj].J.der.Ez; // Right
            for (int k = IzDe[jjj].K.com.Ez; k <= IzDe[jjj].K.fin.Ez; ++k) {
                for (int i = IzDe[jjj].I.com.Ez; i <= IzDe[jjj].I.fin.Ez; ++i) {
                    if (media.sggMiEz(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Right TF/SF region intersects a material at Ez %7d%7d%7d", i, j, k);
                        if (((media.sggMiEz(i, j, k) == 0) || (sgg.med(media.sggMiEz(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
            // Ex Right
            j = IzDe[jjj].J.der.Ex; // Right
            for (int k = IzDe[jjj].K.com.Ex; k <= IzDe[jjj].K.fin.Ex; ++k) {
                for (int i = IzDe[jjj].I.com.Ex; i <= IzDe[jjj].I.fin.Ex; ++i) {
                    if (media.sggMiEx(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Right TF/SF region intersects a material at Ex %7d%7d%7d", i, j, k);
                        if (((media.sggMiEx(i, j, k) == 0) || (sgg.med(media.sggMiEx(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
        }
        
        // Down TF/SF region intersects a material at Ex
        if (IluminaAb[jjj]) {
            int k = AbAr[jjj].K.aba.Ex; // Down
            for (int j = AbAr[jjj].J.com.Ex; j <= AbAr[jjj].J.fin.Ex; ++j) {
                for (int i = AbAr[jjj].I.com.Ex; i <= AbAr[jjj].I.fin.Ex; ++i) {
                    if (media.sggMiEx(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Down TF/SF region intersects a material at Ex %7d%7d%7d", i, j, k);
                        if (((media.sggMiEx(i, j, k) == 0) || (sgg.med(media.sggMiEx(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
            // Ey Down
            k = AbAr[jjj].K.aba.Ey; // Down
            for (int j = AbAr[jjj].J.com.Ey; j <= AbAr[jjj].J.fin.Ey; ++j) {
                for (int i = AbAr[jjj].I.com.Ey; i <= AbAr[jjj].I.fin.Ey; ++i) {
                    if (media.sggMiEy(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Down TF/SF region intersects a material at Ey %7d%7d%7d", i, j, k);
                        if (((media.sggMiEy(i, j, k) == 0) || (sgg.med(media.sggMiEy(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
        }
        
        // Up TF/SF region intersects a material at Ex
        if (IluminaAr[jjj]) {
            int k = AbAr[jjj].K.arr.Ex; // Up
            for (int j = AbAr[jjj].J.com.Ex; j <= AbAr[jjj].J.fin.Ex; ++j) {
                for (int i = AbAr[jjj].I.com.Ex; i <= AbAr[jjj].I.fin.Ex; ++i) {
                    if (media.sggMiEx(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Up TF/SF region intersects a material at Ex %7d%7d%7d", i, j, k);
                        if (((media.sggMiEx(i, j, k) == 0) || (sgg.med(media.sggMiEx(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
            // Ey Up
            k = AbAr[jjj].K.arr.Ey;
            for (int j = AbAr[jjj].J.com.Ey; j <= AbAr[jjj].J.fin.Ey; ++j) {
                for (int i = AbAr[jjj].I.com.Ey; i <= AbAr[jjj].I.fin.Ey; ++i) {
                    if (media.sggMiEy(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Up TF/SF region intersects a material at Ey %7d%7d%7d", i, j, k);
                        if (((media.sggMiEy(i, j, k) == 0) || (sgg.med(media.sggMiEy(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
        }
        
        // Back TF/SF region intersects a material at Hz
        if (IluminaTr[jjj]) {
            int i = TrFr[jjj].I.tra.Hz; // Back
            for (int k = TrFr[jjj].K.com.Hz; k <= TrFr[jjj].K.fin.Hz; ++k) {
                for (int j = TrFr[jjj].J.com.Hz; j <= TrFr[jjj].J.fin.Hz; ++j) {
                    if (media.sggMiHz(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Back TF/SF region intersects a material at Hz %7d%7d%7d", i, j, k);
                        if (((media.sggMiHz(i, j, k) == 0) || (sgg.med(media.sggMiHz(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
            // Hy Back
            i = TrFr[jjj].I.tra.Hy; // Back
            for (int k = TrFr[jjj].K.com.Hy; k <= TrFr[jjj].K.fin.Hy; ++k) {
                for (int j = TrFr[jjj].J.com.Hy; j <= TrFr[jjj].J.fin.Hy; ++j) {
                    if (media.sggMiHy(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Back TF/SF region intersects a material at Hy %7d%7d%7d", i, j, k);
                        if (((media.sggMiHy(i, j, k) == 0) || (sgg.med(media.sggMiHy(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
        }
        
        // Front TF/SF region intersects a material at Hz
        if (IluminaFr[jjj]) {
            int i = TrFr[jjj].I.fro.Hz; // Front
            for (int k = TrFr[jjj].K.com.Hz; k <= TrFr[jjj].K.fin.Hz; ++k) {
                for (int j = TrFr[jjj].J.com.Hz; j <= TrFr[jjj].J.fin.Hz; ++j) {
                    if (media.sggMiHz(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Front TF/SF region intersects a material at Hz %7d%7d%7d", i, j, k);
                        if (((media.sggMiHz(i, j, k) == 0) || (sgg.med(media.sggMiHz(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
            // Hy Front
            i = TrFr[jjj].I.fro.Hy; // Front
            for (int k = TrFr[jjj].K.com.Hy; k <= TrFr[jjj].K.fin.Hy; ++k) {
                for (int j = TrFr[jjj].J.com.Hy; j <= TrFr[jjj].J.fin.Hy; ++j) {
                    if (media.sggMiHy(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Front TF/SF region intersects a material at Hy %7d%7d%7d", i, j, k);
                        if (((media.sggMiHy(i, j, k) == 0) || (sgg.med(media.sggMiHy(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
        }
        
        // Left TF/SF region intersects a material at Hx
        if (IluminaIz[jjj]) {
            int j = IzDe[jjj].J.izq.Hx; // Left
            for (int k = IzDe[jjj].K.com.Hx; k <= IzDe[jjj].K.fin.Hx; ++k) {
                for (int i = IzDe[jjj].I.com.Hx; i <= IzDe[jjj].I.fin.Hx; ++i) {
                    if (media.sggMiHx(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Left TF/SF region intersects a material at Hx %7d%7d%7d", i, j, k);
                        if (((media.sggMiHx(i, j, k) == 0) || (sgg.med(media.sggMiHx(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
            // Hz Left
            j = IzDe[jjj].J.izq.Hz; // Left
            for (int k = IzDe[jjj].K.com.Hz; k <= IzDe[jjj].K.fin.Hz; ++k) {
                for (int i = IzDe[jjj].I.com.Hz; i <= IzDe[jjj].I.fin.Hz; ++i) {
                    if (media.sggMiHz(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Left TF/SF region intersects a material at Hz %7d%7d%7d", i, j, k);
                        if (((media.sggMiHz(i, j, k) == 0) || (sgg.med(media.sggMiHz(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
        }
        
        // Right TF/SF region intersects a material at Hx
        if (IluminaDe[jjj]) {
            int j = IzDe[jjj].J.der.Hx; // Right
            for (int k = IzDe[jjj].K.com.Hx; k <= IzDe[jjj].K.fin.Hx; ++k) {
                for (int i = IzDe[jjj].I.com.Hx; i <= IzDe[jjj].I.fin.Hx; ++i) {
                    if (media.sggMiHx(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Right TF/SF region intersects a material at Hx %7d%7d%7d", i, j, k);
                        if (((media.sggMiHx(i, j, k) == 0) || (sgg.med(media.sggMiHx(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
            // Hz Right
            j = IzDe[jjj].J.der.Hz; // Right
            for (int k = IzDe[jjj].K.com.Hz; k <= IzDe[jjj].K.fin.Hz; ++k) {
                for (int i = IzDe[jjj].I.com.Hz; i <= IzDe[jjj].I.fin.Hz; ++i) {
                    if (media.sggMiHz(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Right TF/SF region intersects a material at Hz %7d%7d%7d", i, j, k);
                        if (((media.sggMiHz(i, j, k) == 0) || (sgg.med(media.sggMiHz(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
        }
        
        // Down TF/SF region intersects a material at Hx
        if (IluminaAb[jjj]) {
            int k = AbAr[jjj].K.aba.Hx; // Down
            for (int j = AbAr[jjj].J.com.Hx; j <= AbAr[jjj].J.fin.Hx; ++j) {
                for (int i = AbAr[jjj].I.com.Hx; i <= AbAr[jjj].I.fin.Hx; ++i) {
                    if (media.sggMiHx(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Down TF/SF region intersects a material at Hx %7d%7d%7d", i, j, k);
                        if (((media.sggMiHx(i, j, k) == 0) || (sgg.med(media.sggMiHx(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
            // Hy Down
            k = AbAr[jjj].K.aba.Hy; // Down
            for (int j = AbAr[jjj].J.com.Hy; j <= AbAr[jjj].J.fin.Hy; ++j) {
                for (int i = AbAr[jjj].I.com.Hy; i <= AbAr[jjj].I.fin.Hy; ++i) {
                    if (media.sggMiHy(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Down TF/SF region intersects a material at Hy %7d%7d%7d", i, j, k);
                        if (((media.sggMiHy(i, j, k) == 0) || (sgg.med(media.sggMiHy(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
        }
        
        // Up TF/SF region intersects a material at Hx
        if (IluminaAr[jjj]) {
            int k = AbAr[jjj].K.arr.Hx; // Up
            for (int j = AbAr[jjj].J.com.Hx; j <= AbAr[jjj].J.fin.Hx; ++j) {
                for (int i = AbAr[jjj].I.com.Hx; i <= AbAr[jjj].I.fin.Hx; ++i) {
                    if (media.sggMiHx(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Up TF/SF region intersects a material at Hx %7d%7d%7d", i, j, k);
                        if (((media.sggMiHx(i, j, k) == 0) || (sgg.med(media.sggMiHx(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
            // Hy Up
            k = AbAr[jjj].K.arr.Hy; // Up
            for (int j = AbAr[jjj].J.com.Hy; j <= AbAr[jjj].J.fin.Hy; ++j) {
                for (int i = AbAr[jjj].I.com.Hy; i <= AbAr[jjj].I.fin.Hy; ++i) {
                    if (media.sggMiHy(i, j, k) != 1) {
                        char buff[BUFSIZE];
                        snprintf(buff, BUFSIZE, "Up TF/SF region intersects a material at Hy %7d%7d%7d", i, j, k);
                        if (((media.sggMiHy(i, j, k) == 0) || (sgg.med(media.sggMiHy(i, j, k)).is.pec)) && 
                            !((i == sgg.SINPMLSweep(iHx_const).XI) || (j == sgg.SINPMLSweep(iHy_const).YI) || (k == sgg.SINPMLSweep(iHz_const).ZI) ||
                              (i == sgg.SINPMLSweep(iHx_const).XE) || (j == sgg.SINPMLSweep(iHy_const).YE) || (k == sgg.SINPMLSweep(iHz_const).ZE))) {
                            stoponerror(layoutnumber, num_procs, std::string(buff));
                        }
                    }
                }
            }
        }
    } // end do jjj

    calc_planewaveconstants(sgg, eps0, mu0);
}

// Function Incid
RKIND Incid(const SGGFDTDINFO_t& sgg, int jjj, int nfield, double time, int i, int j, int k, bool& still_planewave_time, bool calledfromobservation) {
    RKIND EhI = 0.0;
    double xf = Punto.PhysCoor[nfield].x[i];
    double yf = Punto.PhysCoor[nfield].y[j];
    double zf = Punto.PhysCoor[nfield].z[k];

    if (calledfromobservation) {
        // Parallel region would be handled by OpenMP pragmas in original
        // For sequential translation, we just loop
        for (int jdum = 1; jdum <= sgg.numplanewaves; ++jdum) {
            for (int kkk = 1; kkk <= sgg.PlaneWave[jdum].nummodes; ++kkk) {
                double d = (xf * pxpw[jdum][kkk] + yf * pypw[jdum][kkk] + zf * pzpw[jdum][kkk]) - distanciaInicial[jdum][kkk];
                EhI += fpw[jdum][nfield][kkk] * evolucion(jdum, time, d, still_planewave_time);
            }
        }
    } else {
        for (int kkk = 1; kkk <= sgg.PlaneWave[jjj].nummodes; ++kkk) {
            double d = (xf * pxpw[jjj][kkk] + yf * pypw[jjj][kkk] + zf * pzpw[jjj][kkk]) - distanciaInicial[jjj][kkk];
            EhI += fpw[jjj][nfield][kkk] * evolucion(jjj, time, d, still_planewave_time);
        }
    }
    return EhI;
}

// Helper function for Incid: evolucion
RKIND evolucion(int jjj, double t, double d, bool& still_planewave_time) {
    RKIND result = 0.0;
    long long nprev = static_cast<long long>((t - d / cluz) / deltaevol[jjj]);
    
    if ((nprev + 1 <= numus[jjj])) {
        still_planewave_time = true;
        if (nprev > 0) {
            // First order interpolation
            result = (evol[jjj][nprev + 1] - evol[jjj][nprev]) / deltaevol[jjj] * ((t - d / cluz) - nprev * deltaevol[jjj]) + evol[jjj][nprev];
        }
    }
    return result;
}

// DestroyIlumina
void DestroyIlumina(SGGFDTDINFO_t& sgg) {
    for (int field = iEx_const; field <= iHz_const; ++field) {
        if (Punto.PhysCoor[field].x) delete[] Punto.PhysCoor[field].x;
        if (Punto.PhysCoor[field].y) delete[] Punto.PhysCoor[field].y;
        if (Punto.PhysCoor[field].z) delete[] Punto.PhysCoor[field].z;
    }

    if (sgg.numplanewaves >= 1) {
        // Deallocate global arrays
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
        delete[] distanciaInicial;
    }
    if (evol) delete[] evol;
    if (sgg.PlaneWave) delete[] sgg.PlaneWave;
}

// AdvancePlaneWaveE
void AdvancePlaneWaveE(const SGGFDTDINFO_t& sgg, int timeinstant, const std::vector<double>& b_dxh, const std::vector<double>& b_dyh, const std::vector<double>& b_dzh, 
                       std::vector<std::vector<std::vector<RKIND>>>& Ex, std::vector<std::vector<std::vector<RKIND>>>& Ey, std::vector<std::vector<std::vector<RKIND>>>& Ez,
                       bool& still_planewave_time) {
    bool called_fromobservation = false;
    still_planewave_time = false;
    
    double timei = sgg.tiempo[timeinstant];
    RKIND G2_1 = sgg.G2[1]; // Assuming G2 is 0-based or 1-based, Fortran is 1-based. Code uses G2(1).
    
    // Note: The Fortran code passes b as bounds_t. We assume b has members like dxh.NX, Ez.XI, etc.
    // For this translation, we assume the structure of b is accessible.
    
    for (int jjj = 1; jjj <= sgg.numplanewaves; ++jjj) {
        if (IluminaTr[jjj]) {
            // Ez Back
            int i = TrFr[jjj].I.tra.Ez;
            int i_m = i - b.Ez.XI;
            RKIND Id = b_dxh[i_m];
            
            // Parallel region
            for (int k = TrFr[jjj].K.com.Ez; k <= TrFr[jjj].K.fin.Ez; ++k) {
                int k_m = k - b.Ez.ZI;
                for (int j = TrFr[jjj].J.com.Ez; j <= TrFr[jjj].J.fin.Ez; ++j) {
                    int j_m = j - b.Ez.YI;
                    RKIND incidente = Incid(sgg, jjj, iHy_const, timei, i - 1, j, k, still_planewave_time, called_fromobservation);
                    Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] - G2_1 * incidente * Id;
                }
            }
            
            // Ey Back
            i = TrFr[jjj].I.tra.Ey;
            i_m = i - b.Ey.XI;
            Id = b_dxh[i_m];
            
            for (int k = TrFr[jjj].K.com.Ey; k <= TrFr[jjj].K.fin.Ey; ++k) {
                int k_m = k - b.Ey.ZI;
                for (int j = TrFr[jjj].J.com.Ey; j <= TrFr[jjj].J.fin.Ey; ++j) {
                    int j_m = j - b.Ey.YI;
                    RKIND incidente = Incid(sgg, jjj, iHz_const, timei, i - 1, j, k, still_planewave_time, called_fromobservation);
                    Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] + G2_1 * incidente * Id;
                }
            }
        }
        
        if (IluminaFr[jjj]) {
            // Ez Front
            int i = TrFr[jjj].I.fro.Ez;
            int i_m = i - b.Ez.XI;
            RKIND Id = b_dxh[i_m];
            
            // Parallel region
            for (int k = TrFr[jjj].K.com.Ez; k <= TrFr[jjj].K.fin.Ez; ++k) {
                int k_m = k - b.Ez.ZI;
                for (int j = TrFr[jjj].J.com.Ez; j <= TrFr[jjj].J.fin.Ez; ++j) {
                    int j_m = j - b.Ez.YI;
                    RKIND incidente = Incid(sgg, jjj, iHy_const, timei, i + 1, j, k, still_planewave_time, called_fromobservation);
                    Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] + G2_1 * incidente * Id;
                }
            }
            
            // Ey Front
            i = TrFr[jjj].I.fro.Ey;
            i_m = i - b.Ey.XI;
            Id = b_dxh[i_m];
            
            for (int k = TrFr[jjj].K.com.Ey; k <= TrFr[jjj].K.fin.Ey; ++k) {
                int k_m = k - b.Ey.ZI;
                for (int j = TrFr[jjj].J.com.Ey; j <= TrFr[jjj].J.fin.Ey; ++j) {
                    int j_m = j - b.Ey.YI;
                    RKIND incidente = Incid(sgg, jjj, iHz_const, timei, i + 1, j, k, still_planewave_time, called_fromobservation);
                    Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] - G2_1 * incidente * Id;
                }
            }
        }
    }
}

// This is a continuation chunk. Previous chunks defined types like SGGFDTDINFO_t, bounds_t, TrFr_t, IzDe_t, AbAr_t, and functions like Incid, print11, etc.
// Global variables like pxpw, pypw, pzpw, fpw, eps0, mu0, cluz, zvac, SEPARADOR, separador are assumed to be defined elsewhere or in a global namespace.
// RKIND is assumed to be defined (e.g., using double).

void AdvancePlaneWaveH(const SGGFDTDINFO_t& sgg, int timeinstant, const bounds_t& b, double gm2, const std::vector<double>& Idxe, const std::vector<double>& Idye, const std::vector<double>& Idze, std::vector<std::vector<std::vector<double>>>& Hx, std::vector<std::vector<std::vector<double>>>& Hy, std::vector<std::vector<std::vector<double>>>& Hz, bool& still_planewave_time) {
    bool called_fromobservation = false;

    double timei = sgg.tiempo[timeinstant] + 0.5 * sgg.dt;
    double Gm2_1 = gm2[1];

    for (int jjj = 1; jjj <= sgg.numplanewaves; ++jjj) {
        if (IluminaTr(jjj)) {
            // Hz Back
            int i = TrFr[jjj].I.tra.Hz;
            int i_m = i - b.Hx.XI;
            double Id = Idxe[i_m];

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(j, k, j_m, k_m)
#endif
            for (int k = TrFr[jjj].K.com.Hz; k <= TrFr[jjj].K.fin.Hz; ++k) {
                int k_m = k - b.Hx.ZI;
                for (int j = TrFr[jjj].J.com.Hz; j <= TrFr[jjj].J.fin.Hz; ++j) {
                    int j_m = j - b.Hx.YI;
                    double incidente = Incid(sgg, jjj, iEy, timei, i + 1, j, k, still_planewave_time, called_fromobservation);
                    Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] + Gm2_1 * incidente * Id;
                }
            }

            // Hy Back
            i = TrFr[jjj].I.tra.Hy;
            i_m = i - b.Hy.XI;
            Id = Idxe[i_m];

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(j, k, j_m, k_m)
#endif
            for (int k = TrFr[jjj].K.com.Hy; k <= TrFr[jjj].K.fin.Hy; ++k) {
                int k_m = k - b.Hy.ZI;
                for (int j = TrFr[jjj].J.com.Hy; j <= TrFr[jjj].J.fin.Hy; ++j) {
                    int j_m = j - b.Hy.YI;
                    double incidente = Incid(sgg, jjj, iEz, timei, i + 1, j, k, still_planewave_time, called_fromobservation);
                    Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] - Gm2_1 * incidente * Id;
                }
            }
        }

        if (IluminaFr(jjj)) {
            // Hz Front
            int i = TrFr[jjj].I.fro.Hz;
            int i_m = i - b.Hx.XI;
            double Id = Idxe[i_m];

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(j, k, j_m, k_m)
#endif
            for (int k = TrFr[jjj].K.com.Hz; k <= TrFr[jjj].K.fin.Hz; ++k) {
                int k_m = k - b.Hx.ZI;
                for (int j = TrFr[jjj].J.com.Hz; j <= TrFr[jjj].J.fin.Hz; ++j) {
                    int j_m = j - b.Hx.YI;
                    double incidente = Incid(sgg, jjj, iEy, timei, i, j, k, still_planewave_time, called_fromobservation);
                    Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] - Gm2_1 * incidente * Id;
                }
            }

            // Hy Front
            i = TrFr[jjj].I.fro.Hy;
            i_m = i - b.Hy.XI;
            Id = Idxe[i_m];

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(j, k, j_m, k_m)
#endif
            for (int k = TrFr[jjj].K.com.Hy; k <= TrFr[jjj].K.fin.Hy; ++k) {
                int k_m = k - b.Hy.ZI;
                for (int j = TrFr[jjj].J.com.Hy; j <= TrFr[jjj].J.fin.Hy; ++j) {
                    int j_m = j - b.Hy.YI;
                    double incidente = Incid(sgg, jjj, iEz, timei, i, j, k, still_planewave_time, called_fromobservation);
                    Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] + Gm2_1 * incidente * Id;
                }
            }
        }

        if (IluminaIz(jjj)) {
            // Hx Left
            int j = IzDe[jjj].J.izq.Hx;
            int j_m = j - b.Hx.YI;
            double Id = Idye[j_m];

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(k, i, k_m, i_m)
#endif
            for (int k = IzDe[jjj].K.com.Hx; k <= IzDe[jjj].K.fin.Hx; ++k) {
                int k_m = k - b.Hx.ZI;
                for (int i = IzDe[jjj].I.com.Hx; i <= IzDe[jjj].I.fin.Hx; ++i) {
                    int i_m = i - b.Hx.XI;
                    double incidente = Incid(sgg, jjj, iEz, timei, i, j + 1, k, still_planewave_time, called_fromobservation);
                    Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] + Gm2_1 * incidente * Id;
                }
            }

            // Hz Left
            j = IzDe[jjj].J.izq.Hz;
            j_m = j - b.Hx.YI;
            Id = Idye[j_m];

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(k, i, k_m, i_m)
#endif
            for (int k = IzDe[jjj].K.com.Hz; k <= IzDe[jjj].K.fin.Hz; ++k) {
                int k_m = k - b.Hx.ZI;
                for (int i = IzDe[jjj].I.com.Hz; i <= IzDe[jjj].I.fin.Hz; ++i) {
                    int i_m = i - b.Hx.XI;
                    double incidente = Incid(sgg, jjj, iEx, timei, i, j + 1, k, still_planewave_time, called_fromobservation);
                    Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] - Gm2_1 * incidente * Id;
                }
            }
        }

        if (IluminaDe(jjj)) {
            // Hx Right
            int j = IzDe[jjj].J.der.Hx;
            int j_m = j - b.Hx.YI;
            double Id = Idye[j_m];

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(k, i, k_m, i_m)
#endif
            for (int k = IzDe[jjj].K.com.Hx; k <= IzDe[jjj].K.fin.Hx; ++k) {
                int k_m = k - b.Hx.ZI;
                for (int i = IzDe[jjj].I.com.Hx; i <= IzDe[jjj].I.fin.Hx; ++i) {
                    int i_m = i - b.Hx.XI;
                    double incidente = Incid(sgg, jjj, iEz, timei, i, j, k, still_planewave_time, called_fromobservation);
                    Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] - Gm2_1 * incidente * Id;
                }
            }

            // Hz Right
            j = IzDe[jjj].J.der.Hz;
            j_m = j - b.Hx.YI;
            Id = Idye[j_m];

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(k, i, k_m, i_m)
#endif
            for (int k = IzDe[jjj].K.com.Hz; k <= IzDe[jjj].K.fin.Hz; ++k) {
                int k_m = k - b.Hx.ZI;
                for (int i = IzDe[jjj].I.com.Hz; i <= IzDe[jjj].I.fin.Hz; ++i) {
                    int i_m = i - b.Hx.XI;
                    double incidente = Incid(sgg, jjj, iEx, timei, i, j, k, still_planewave_time, called_fromobservation);
                    Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] + Gm2_1 * incidente * Id;
                }
            }
        }

        if (IluminaAb(jjj)) {
            // Hx Down
            int k = AbAr[jjj].K.aba.Hx;
            int k_m = k - b.Hx.ZI;
            double Id = Idze[k_m];

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, i_m, j_m)
#endif
            for (int j = AbAr[jjj].J.com.Hx; j <= AbAr[jjj].J.fin.Hx; ++j) {
                int j_m = j - b.Hx.YI;
                for (int i = AbAr[jjj].I.com.Hx; i <= AbAr[jjj].I.fin.Hx; ++i) {
                    int i_m = i - b.Hx.XI;
                    double incidente = Incid(sgg, jjj, iEy, timei, i, j, k + 1, still_planewave_time, called_fromobservation);
                    Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] - Gm2_1 * incidente * Id;
                }
            }

            // Hy Down
            k = AbAr[jjj].K.aba.Hy;
            k_m = k - b.Hy.ZI;
            Id = Idze[k_m];

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, i_m, j_m)
#endif
            for (int j = AbAr[jjj].J.com.Hy; j <= AbAr[jjj].J.fin.Hy; ++j) {
                int j_m = j - b.Hy.YI;
                for (int i = AbAr[jjj].I.com.Hy; i <= AbAr[jjj].I.fin.Hy; ++i) {
                    int i_m = i - b.Hy.XI;
                    double incidente = Incid(sgg, jjj, iEx, timei, i, j, k + 1, still_planewave_time, called_fromobservation);
                    Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] + Gm2_1 * incidente * Id;
                }
            }
        }

        if (IluminaAr(jjj)) {
            // Hx Up
            int k = AbAr[jjj].K.arr.Hx;
            int k_m = k - b.Hx.ZI;
            double Id = Idze[k_m];

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, i_m, j_m)
#endif
            for (int j = AbAr[jjj].J.com.Hx; j <= AbAr[jjj].J.fin.Hx; ++j) {
                int j_m = j - b.Hx.YI;
                for (int i = AbAr[jjj].I.com.Hx; i <= AbAr[jjj].I.fin.Hx; ++i) {
                    int i_m = i - b.Hx.XI;
                    double incidente = Incid(sgg, jjj, iEy, timei, i, j, k, still_planewave_time, called_fromobservation);
                    Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] + Gm2_1 * incidente * Id;
                }
            }

            // Hy Up
            k = AbAr[jjj].K.arr.Hy;
            k_m = k - b.Hy.ZI;
            Id = Idze[k_m];

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, i_m, j_m)
#endif
            for (int j = AbAr[jjj].J.com.Hy; j <= AbAr[jjj].J.fin.Hy; ++j) {
                int j_m = j - b.Hy.YI;
                for (int i = AbAr[jjj].I.com.Hy; i <= AbAr[jjj].I.fin.Hy; ++i) {
                    int i_m = i - b.Hy.XI;
                    double incidente = Incid(sgg, jjj, iEx, timei, i, j, k, still_planewave_time, called_fromobservation);
                    Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] - Gm2_1 * incidente * Id;
                }
            }
        }
    }
}

void storeplanewaves(const SGGFDTDINFO_t& sgg) {
    for (int jjj = 1; jjj <= sgg.numplanewaves; ++jjj) {
        for (int kkk = 1; kkk <= sgg.PlaneWave[jjj].nummodes; ++kkk) {
            if (sgg.PlaneWave[jjj].isRC) {
                // Assuming print11 handles the formatting similar to Fortran write
                // Note: In C++, we might need to construct the string or use fprintf directly if print11 is a wrapper
                // For now, assuming print11 signature matches the intent
                print11(0, SEPARADOR + separador + separador);
                // The actual write content: pxpw(jjj,kkk), pypw(jjj,kkk), pzpw(jjj,kkk), fpw(jjj,1,kkk), fpw(jjj,2,kkk), fpw(jjj,3,kkk), sgg%PlaneWave(jjj)%incert(kkk)
                // Since print11 is a custom function, we pass the arguments. 
                // The Fortran write(14,err=634) implies writing to unit 14. 
                // We assume print11(0, ...) is the error handler or a specific print function. 
                // The actual data writing logic needs to be mapped. 
                // Given the error handling structure, let's assume a helper or direct file write.
                // However, the prompt says "Translate ONLY this chunk". 
                // The Fortran code: write(14,err=634) ...
                // C++ equivalent would likely be fprintf or a custom write function.
                // Since print11 is called in the error block, it's likely a logging function.
                // Let's assume there's a global file pointer or a write function.
                // For simplicity and preserving names, we'll call a hypothetical write function or use fprintf if 14 is a file descriptor.
                // But to be safe with "preserve names", we just translate the logic.
                
                // Re-reading: write(14,err=634) ... goto 635; 634 call print11...
                // This is a standard Fortran error handling pattern.
                // In C++, we might use a try-catch or check a return value.
                // Let's assume a function `write_planewave_data` or similar exists, or we use fprintf.
                // Given the constraints, I will use fprintf to unit 14 (assuming it's a FILE*).
                // If 14 is not a FILE*, this translation assumes a global `FILE* unit14`.
                
                // To strictly follow "preserve names" and not invent globals, I will assume the existence of a mechanism to write to unit 14.
                // However, without context, I'll use a placeholder or standard C++ file I/O if 14 is known.
                // Let's assume `unit14` is a global FILE*.
                fprintf(unit14, "%g %g %g %g %g %g %g\n", 
                       pxpw[jjj][kkk], pypw[jjj][kkk], pzpw[jjj][kkk], 
                       fpw[jjj][1][kkk], fpw[jjj][2][kkk], fpw[jjj][3][kkk], 
                       sgg.PlaneWave[jjj].incert[kkk]);
            }
        }
    }
    return;
    
634:
    print11(0, SEPARADOR + separador + separador);
    print11(0, "PLANEWAVES: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING");
    print11(0, SEPARADOR + separador + separador);
635:
    return;
}

void calc_planewaveconstants(const SGGFDTDINFO_t& sgg, double eps00, double mu00) {
    eps0 = eps00;
    mu0 = mu00;
    cluz = 1.0 / sqrt(eps0 * mu0);
    zvac = sqrt(mu0 / eps0);

    for (int jjj = 1; jjj <= sgg.numplanewaves; ++jjj) {
        for (int kkk = 1; kkk <= sgg.PlaneWave[jjj].nummodes; ++kkk) {
            fpw[jjj][4][kkk] = (pypw[jjj][kkk] * fpw[jjj][3][kkk] - pzpw[jjj][kkk] * fpw[jjj][2][kkk]) / zvac;
            fpw[jjj][5][kkk] = (pzpw[jjj][kkk] * fpw[jjj][1][kkk] - pxpw[jjj][kkk] * fpw[jjj][3][kkk]) / zvac;
            fpw[jjj][6][kkk] = (pxpw[jjj][kkk] * fpw[jjj][2][kkk] - pypw[jjj][kkk] * fpw[jjj][1][kkk]) / zvac;
        }
    }
}

void corrigeondaplanaH(const SGGFDTDINFO_t& sgg, const bounds_t& b, std::vector<std::vector<std::vector<double>>>& Hx, std::vector<std::vector<std::vector<double>>>& Hy, std::vector<std::vector<std::vector<double>>>& Hz, std::vector<std::vector<std::vector<double>>>& Hxvac, std::vector<std::vector<std::vector<double>>>& Hyvac, std::vector<std::vector<std::vector<double>>>& Hzvac) {
    for (int jjj = 1; jjj <= sgg.numplanewaves; ++jjj) {
        if (IluminaTr(jjj)) {
            // Hz Back
            int i = TrFr[jjj].I.tra.Hz;
            int i_m = i - b.Hx.XI;

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(j, k, j_m, k_m)
#endif
            for (int k = TrFr[jjj].K.com.Hz; k <= TrFr[jjj].K.fin.Hz; ++k) {
                int k_m = k - b.Hx.ZI;
                for (int j = TrFr[jjj].J.com.Hz; j <= TrFr[jjj].J.fin.Hz; ++j) {
                    int j_m = j - b.Hx.YI;
                    Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] - Hzvac[i_m][j_m][k_m];
                }
            }

            // Hy Back
            i = TrFr[jjj].I.tra.Hy;
            i_m = i - b.Hy.XI;

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(j, k, j_m, k_m)
#endif
            for (int k = TrFr[jjj].K.com.Hy; k <= TrFr[jjj].K.fin.Hy; ++k) {
                int k_m = k - b.Hy.ZI;
                for (int j = TrFr[jjj].J.com.Hy; j <= TrFr[jjj].J.fin.Hy; ++j) {
                    int j_m = j - b.Hy.YI;
                    Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] - Hyvac[i_m][j_m][k_m];
                }
            }
        }

        if (IluminaFr(jjj)) {
            // Hz Front
            int i = TrFr[jjj].I.fro.Hz;
            int i_m = i - b.Hx.XI;

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(j, k, j_m, k_m)
#endif
            for (int k = TrFr[jjj].K.com.Hz; k <= TrFr[jjj].K.fin.Hz; ++k) {
                int k_m = k - b.Hx.ZI;
                for (int j = TrFr[jjj].J.com.Hz; j <= TrFr[jjj].J.fin.Hz; ++j) {
                    int j_m = j - b.Hx.YI;
                    Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] - Hzvac[i_m][j_m][k_m];
                }
            }

            // Hy Front
            i = TrFr[jjj].I.fro.Hy;
            i_m = i - b.Hy.XI;

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(j, k, j_m, k_m)
#endif
            for (int k = TrFr[jjj].K.com.Hy; k <= TrFr[jjj].K.fin.Hy; ++k) {
                int k_m = k - b.Hy.ZI;
                for (int j = TrFr[jjj].J.com.Hy; j <= TrFr[jjj].J.fin.Hy; ++j) {
                    int j_m = j - b.Hy.YI;
                    Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] - Hyvac[i_m][j_m][k_m];
                }
            }
        }
    }
}

#ifdef CompileWithOpenMP
#include <omp.h>
#endif

// Assuming previous chunks have defined the necessary types and structures:
// - TrFr, IzDe, AbAr structures with nested com/fin/izq/der/aba/arr members for I, J, K
// - b structure with XI, YI, ZI members for Hx, Hy, Hz
// - Hx, Hy, Hz, Hxvac, Hyvac arrays (likely 3D arrays or flattened)
// - IluminaIz, IluminaDe, IluminaAb, IluminaAr functions or arrays
// - corrigeondaplanaH subroutine context

// Note: The original code uses 1-based indexing for loops but accesses arrays with adjusted indices.
// We preserve the logic exactly. Array access assumes 0-based or 1-based depending on previous definitions.
// Since Fortran arrays are often 1-based, we assume the arrays Hx, Hy, Hz, etc., are accessed with
// indices that might need adjustment if C++ uses 0-based. However, the code calculates i_m, j_m, k_m
// which are likely the actual indices into the C++ arrays. We will assume the arrays are sized appropriately
// and i_m, j_m, k_m are valid indices.

void ilumina_m::corrigeondaplanaH() {
    // Assuming 'jjj' is a loop variable from an outer scope not shown in this chunk.
    // Based on the structure, this is likely inside a loop over 'jjj'.
    // Since this is a continuation, we assume 'jjj' is defined in the calling context or this function
    // takes 'jjj' as an argument. However, the end of the chunk shows 'end do' and 'return',
    // suggesting this is the end of a subroutine. The 'do k = ...' loops are inside an 'if' block
    // which is inside a loop over 'jjj'.
    // Let's assume this code is part of a larger function where 'jjj' is iterated.
    // But the prompt says "Translate ONLY this chunk". The chunk starts with a do loop for k.
    // This implies the outer loop over 'jjj' was in a previous chunk.
    // We will translate the code as is, assuming 'jjj' is available in the scope.

    // However, looking at the end: "end do" matches the start "do k = ..."? No.
    // The structure is:
    // do k = ...
    //    ...
    // end do
    // #ifdef ...
    // end if
    // if (IluminaIz(jjj)) { ... }
    // ...
    // end do  <-- This likely closes a loop over 'jjj' that started in a previous chunk.
    // return
    // end subroutine

    // Since I cannot see the start of the 'jjj' loop, I will translate the visible code.
    // I will assume 'jjj' is a local variable or parameter.

    // To make this compile, we need to assume the types.
    // Let's assume the following types are defined in the header or previous chunks:
    // struct IndexRange { int com, fin; };
    // struct TrFrType { struct { IndexRange I, J, K; } com, fin; };
    // Similar for IzDe and AbAr.
    // struct Bounds { int XI, YI, ZI; };
    // struct bType { Bounds Hx, Hy, Hz; };
    // Arrays Hx, Hy, Hz, Hxvac, Hyvac are likely std::vector or raw arrays.

    // We will write the code assuming these types exist.

    // Note: The original code has OpenMP directives. We preserve them.

    // The 'do k' loop at the beginning:
    for (int k = TrFr[jjj].K.com.Hy; k <= TrFr[jjj].K.fin.Hy; ++k) {
        int k_m = k - b.Hy.ZI;
        for (int j = TrFr[jjj].J.com.Hy; j <= TrFr[jjj].J.fin.Hy; ++j) {
            int j_m = j - b.Hy.YI;
            Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] - Hyvac[i_m][j_m][k_m];
        }
    }

#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif

    if (IluminaIz(jjj)) {
        // Hx Left
        int j = IzDe[jjj].J.izq.Hx; // Left
        int j_m = j - b.Hx.YI;

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(k, i, k_m, i_m)
#endif
        for (int k = IzDe[jjj].K.com.Hx; k <= IzDe[jjj].K.fin.Hx; ++k) {
            int k_m = k - b.Hx.ZI;
            for (int i = IzDe[jjj].I.com.Hx; i <= IzDe[jjj].I.fin.Hx; ++i) {
                int i_m = i - b.Hx.XI;
                Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] - Hxvac[i_m][j_m][k_m];
            }
        }

#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif

        // Hz Left
        j = IzDe[jjj].J.izq.Hz; // Left
        j_m = j - b.Hz.YI;

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(k, i, k_m, i_m)
#endif
        for (int k = IzDe[jjj].K.com.Hz; k <= IzDe[jjj].K.fin.Hz; ++k) {
            int k_m = k - b.Hz.ZI;
            for (int i = IzDe[jjj].I.com.Hz; i <= IzDe[jjj].I.fin.Hz; ++i) {
                int i_m = i - b.Hx.XI; // Note: Original code uses b%Hx%XI here, likely a typo in Fortran or intentional
                Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] - Hzvac[i_m][j_m][k_m];
            }
        }

#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
    }

    if (IluminaDe(jjj)) {
        // Hx Right
        int j = IzDe[jjj].J.der.Hx; // Right
        j_m = j - b.Hx.YI;

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(k, i, k_m, i_m)
#endif
        for (int k = IzDe[jjj].K.com.Hx; k <= IzDe[jjj].K.fin.Hx; ++k) {
            int k_m = k - b.Hx.ZI;
            for (int i = IzDe[jjj].I.com.Hx; i <= IzDe[jjj].I.fin.Hx; ++i) {
                int i_m = i - b.Hx.XI;
                Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] - Hxvac[i_m][j_m][k_m];
            }
        }

#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif

        // Hz Right
        j = IzDe[jjj].J.der.Hz; // Right
        j_m = j - b.Hz.YI;

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(k, i, k_m, i_m)
#endif
        for (int k = IzDe[jjj].K.com.Hz; k <= IzDe[jjj].K.fin.Hz; ++k) {
            int k_m = k - b.Hz.ZI;
            for (int i = IzDe[jjj].I.com.Hz; i <= IzDe[jjj].I.fin.Hz; ++i) {
                int i_m = i - b.Hx.XI; // Note: Original code uses b%Hx%XI here
                Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] - Hzvac[i_m][j_m][k_m];
            }
        }

#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
    }

    if (IluminaAb(jjj)) {
        // Hx Down
        int k = AbAr[jjj].K.aba.Hx; // Down
        int k_m = k - b.Hx.ZI;

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, i_m, j_m)
#endif
        for (int j = AbAr[jjj].J.com.Hx; j <= AbAr[jjj].J.fin.Hx; ++j) {
            int j_m = j - b.Hx.YI;
            for (int i = AbAr[jjj].I.com.Hx; i <= AbAr[jjj].I.fin.Hx; ++i) {
                int i_m = i - b.Hx.XI;
                Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] - Hxvac[i_m][j_m][k_m];
            }
        }

#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif

        // Hy Down
        k = AbAr[jjj].K.aba.Hy; // Down
        k_m = k - b.Hy.ZI;

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, i_m, j_m)
#endif
        for (int j = AbAr[jjj].J.com.Hy; j <= AbAr[jjj].J.fin.Hy; ++j) {
            int j_m = j - b.Hy.YI;
            for (int i = AbAr[jjj].I.com.Hy; i <= AbAr[jjj].I.fin.Hy; ++i) {
                int i_m = i - b.Hy.XI;
                Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] - Hyvac[i_m][j_m][k_m];
            }
        }

#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
    }

    if (IluminaAr(jjj)) {
        // Hx Up
        int k = AbAr[jjj].K.arr.Hx; // Up
        int k_m = k - b.Hx.ZI;

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, i_m, j_m)
#endif
        for (int j = AbAr[jjj].J.com.Hx; j <= AbAr[jjj].J.fin.Hx; ++j) {
            int j_m = j - b.Hx.YI;
            for (int i = AbAr[jjj].I.com.Hx; i <= AbAr[jjj].I.fin.Hx; ++i) {
                int i_m = i - b.Hx.XI;
                Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] - Hxvac[i_m][j_m][k_m];
            }
        }

#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif

        // Hy Up
        int k = AbAr[jjj].K.arr.Hy; // Up
        k_m = k - b.Hy.ZI;

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, i_m, j_m)
#endif
        for (int j = AbAr[jjj].J.com.Hy; j <= AbAr[jjj].J.fin.Hy; ++j) {
            int j_m = j - b.Hy.YI;
            for (int i = AbAr[jjj].I.com.Hy; i <= AbAr[jjj].I.fin.Hy; ++i) {
                int i_m = i - b.Hy.XI;
                Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] - Hyvac[i_m][j_m][k_m];
            }
        }

#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
    }
}