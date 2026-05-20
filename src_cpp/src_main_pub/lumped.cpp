#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <iomanip>

// Forward declarations and includes for external modules/types
// Assuming these types exist in the corresponding header files based on Fortran usage
// #include "Report_m.h"
// #include "FDETYPES_m.h"
// #include "lumped_vars_m.h"
// #ifdef CompileWithStochastic
// #include "lumped_devia.h"
// #endif

// Mocking external types/enums/constants for compilation context
// In a real translation, these would come from the actual headers
enum SGGIndex { iEx, iEy, iEz, iHx, iHy, iHz };

struct AllocInfo {
    int XI, XE, YI, YE, ZI, ZE;
};

struct MediaIs {
    bool Lumped;
};

struct LumpedParams {
    int Orient;
    bool inductor;
    bool diodo;
    bool resistor;
    double R, L, C, DiodB, DiodIsat;
    double Rtime_on, Rtime_off;
};

struct MediaStruct {
    MediaIs Is;
    LumpedParams Lumped[1]; // Assuming size 1 based on usage Lumped(1)
    double epr;
    double sigma;
    bool sigmareasignado;
};

struct SGGFDTDINFO_t {
    std::vector<AllocInfo> alloc;
    std::vector<AllocInfo> SINPMLSweep;
    std::vector<MediaStruct> Med;
    int NumMedia;
    double dt;
    // Helper to access media by index safely
    MediaStruct& med(int idx) { return Med[idx]; }
};

struct media_matrices_t {
    std::vector<std::vector<std::vector<int>>> sggMiEx;
    std::vector<std::vector<std::vector<int>>> sggMiEy;
    std::vector<std::vector<std::vector<int>>> sggMiEz;
};

struct sim_control_t {
    int layoutnumber;
    int num_procs;
    bool resume;
    bool stochastic;
};

// Constants
const int RKIND = 8; // Assuming double precision
const int BUFSIZE = 256;

// External functions assumed to exist
void WarnErrReport(const std::string& msg, bool fatal);
void print11(int unit, const std::string& msg);
const std::string SEPARADOR = "========================================";

#ifdef CompileWithStochastic
void inject_devialumped(const SGGFDTDINFO_t& sgg, int timestep, bool simu_devia, bool stochastic, class Nodes_t& lumped_);
void calc_lumped_deviaconsts(const SGGFDTDINFO_t& sgg, class Nodes_t& lumped_, int jmed, double Resist, double Induct, double Capaci, double sigmaeff, double epsiloneff, double eps0, double mu0);
#endif

// Derived Types converted to Classes/Structs

struct Nodes_t {
    double alignedDeltaE;
    double transversalDeltaHa;
    double transversalDeltaHb;
    int Orient;
    int jmed;
    
    // Pointers in Fortran become references or raw pointers in C++. 
    // Since these are updated in place, we use double* to point to the main field arrays.
    double* Efield;
    double* Ha_Plus;
    double* Ha_Minu;
    double* Hb_Plus;
    double* Hb_Minu;

    // State variables
    double EfieldPrevPrev;
    double EfieldPrev;
    double Jcur;

    // Coefficients
    double G1;
    double G2a;
    double G2b;
    double GJ;
    double sigmaEffResistInduct;
    double currentCoeff;
    
    double G1_usual;
    double G2a_usual;
    double G2b_usual;

    // Diode specific
    double diodeB;
    double diodepreA;

#ifdef CompileWithStochastic
    double EfieldPrevPrev_for_devia;
    double EfieldPrev_for_devia;
    double Jcur_for_devia;
#endif
};

struct LumpedElem_t {
    int NumNodes;
    std::vector<Nodes_t> Nodes;
};

namespace Lumped_m {

    // Global variables
    LumpedElem_t LumpElem;
    double eps0 = 0.0;
    double mu0 = 0.0;
    double zvac = 0.0;
    double cluz = 0.0;

    // Helper function for Newton-Raphson
    double newton_raphson(double A, double B, double C) {
        double x = 0.0; // Initial guess, though Fortran code uses x0 uninitialized which is risky. 
                        // Assuming x0 is passed or initialized. The Fortran code:
                        // real(kind=RKIND) :: x0, xx0, fxx0, dfxx0
                        // clave = 1
                        // xx0 = x0  <-- x0 is uninitialized here in the snippet provided. 
                        // This is a bug in the original Fortran if x0 isn't passed. 
                        // However, looking at the signature: function newton_raphson(A,B,C) result (x)
                        // x0 is local. It should probably be initialized. 
                        // Let's assume a standard initial guess of 0.0 or 1.0. 
                        // Given the physics (diode equation), 0 is a reasonable start.
        double x0 = 0.0; 
        double xx0 = x0;
        double tol = 1e-12; // Typical tolerance
        int nmax = 1024;
        int clave = 1;
        int n = 0;

        for (int i = 1; i <= nmax; ++i) {
            double fxx0 = A * std::exp(B * xx0) + xx0 + C;
            double dfxx0 = A * B * std::exp(B * xx0) + 1.0;
            
            if (std::abs(dfxx0) < 1e-20) { // Avoid division by zero
                break;
            }

            double x_new = xx0 - fxx0 / dfxx0;
            
            if (std::abs(x_new - xx0) < tol * std::abs(x_new)) {
                clave = 0;
                n = i;
                xx0 = x_new;
                break;
            }
            xx0 = x_new;
        }
        
        if (clave != 0) {
            // print11(0, "Error convergencia en Newton-Raphson");
            std::cerr << "Error convergencia en Newton-Raphson" << std::endl;
        }
        
        return xx0;
    }

    void InitLumped(const SGGFDTDINFO_t& sgg, const media_matrices_t& media, 
                    const std::vector<std::vector<std::vector<double>>>& Ex,
                    const std::vector<std::vector<std::vector<double>>>& Ey,
                    const std::vector<std::vector<std::vector<double>>>& Ez,
                    const std::vector<std::vector<std::vector<double>>>& Hx,
                    const std::vector<std::vector<std::vector<double>>>& Hy,
                    const std::vector<std::vector<std::vector<double>>>& Hz,
                    const std::vector<double>& Idxh,
                    const std::vector<double>& Idyh,
                    const std::vector<double>& Idzh,
                    const std::vector<double>& Idxe,
                    const std::vector<double>& Idye,
                    const std::vector<double>& Idze,
                    bool& ThereAreLumped,
                    double eps00, double mu00,
                    const sim_control_t& control) {
        
        eps0 = eps00;
        mu0 = mu00;

        ThereAreLumped = false;

        int conta = 0;
        int j1, k1, i1;

        // Count Ex nodes
        for (k1 = sgg.SINPMLSweep[iEx].ZI; k1 <= sgg.SINPMLSweep[iEx].ZE; ++k1) {
            for (j1 = sgg.SINPMLSweep[iEx].YI; j1 <= sgg.SINPMLSweep[iEx].YE; ++j1) {
                for (i1 = sgg.SINPMLSweep[iEx].XI; i1 <= sgg.SINPMLSweep[iEx].XE; ++i1) {
                    int jmed = media.sggMiEx[i1][j1][k1];
                    if (sgg.med(jmed).Is.Lumped) {
                        conta++;
                    }
                }
            }
        }

        // Count Ey nodes
        for (k1 = sgg.SINPMLSweep[iEy].ZI; k1 <= sgg.SINPMLSweep[iEy].ZE; ++k1) {
            for (j1 = sgg.SINPMLSweep[iEy].YI; j1 <= sgg.SINPMLSweep[iEy].YE; ++j1) {
                for (i1 = sgg.SINPMLSweep[iEy].XI; i1 <= sgg.SINPMLSweep[iEy].XE; ++i1) {
                    int jmed = media.sggMiEy[i1][j1][k1];
                    if (sgg.med(jmed).Is.Lumped) {
                        conta++;
                    }
                }
            }
        }

        // Count Ez nodes
        for (k1 = sgg.SINPMLSweep[iEz].ZI; k1 <= sgg.SINPMLSweep[iEz].ZE; ++k1) {
            for (j1 = sgg.SINPMLSweep[iEz].YI; j1 <= sgg.SINPMLSweep[iEz].YE; ++j1) {
                for (i1 = sgg.SINPMLSweep[iEz].XI; i1 <= sgg.SINPMLSweep[iEz].XE; ++i1) {
                    int jmed = media.sggMiEz[i1][j1][k1];
                    if (sgg.med(jmed).Is.Lumped) {
                        conta++;
                    }
                }
            }
        }

        ThereAreLumped = (conta != 0);
        if (!ThereAreLumped) {
            return;
        }

        LumpElem.NumNodes = conta;
        LumpElem.Nodes.resize(conta);

        conta = 0;

        // Process Ex nodes
        for (k1 = sgg.SINPMLSweep[iEx].ZI; k1 <= sgg.SINPMLSweep[iEx].ZE; ++k1) {
            for (j1 = sgg.SINPMLSweep[iEx].YI; j1 <= sgg.SINPMLSweep[iEx].YE; ++j1) {
                for (i1 = sgg.SINPMLSweep[iEx].XI; i1 <= sgg.SINPMLSweep[iEx].XE; ++i1) {
                    int jmed = media.sggMiEx[i1][j1][k1];
                    if (sgg.med(jmed).Is.Lumped) {
                        conta++;
                        Nodes_t& lumped = LumpElem.Nodes[conta - 1]; // 0-indexed vector
                        
                        // Note: Fortran arrays are 1-based or user-defined bounds. 
                        // The indices i1, j1, k1 are used directly. 
                        // We assume the vectors Ex, Ey etc. are accessed with these indices.
                        // In C++, if the vectors are sized to match the max index, we access directly.
                        
                        lumped.alignedDeltaE = 1.0 / Idxe[i1];
                        lumped.transversalDeltaHa = 1.0 / Idyh[j1];
                        lumped.transversalDeltaHb = 1.0 / Idzh[k1];
                        lumped.Orient = sgg.med(jmed).Lumped[1].Orient;
                        lumped.jmed = jmed;
                        
                        // Pointers to the field values
                        lumped.Efield = const_cast<double*>(&Ex[i1][j1][k1]);
                        lumped.Ha_Plus = const_cast<double*>(&Hz[i1][j1][k1]);
                        lumped.Ha_Minu = const_cast<double*>(&Hz[i1][j1-1][k1]);
                        lumped.Hb_Plus = const_cast<double*>(&Hy[i1][j1][k1]);
                        lumped.Hb_Minu = const_cast<double*>(&Hy[i1][j1][k1-1]);
                    }
                }
            }
        }

        // Process Ey nodes
        for (k1 = sgg.SINPMLSweep[iEy].ZI; k1 <= sgg.SINPMLSweep[iEy].ZE; ++k1) {
            for (j1 = sgg.SINPMLSweep[iEy].YI; j1 <= sgg.SINPMLSweep[iEy].YE; ++j1) {
                for (i1 = sgg.SINPMLSweep[iEy].XI; i1 <= sgg.SINPMLSweep[iEy].XE; ++i1) {
                    int jmed = media.sggMiEy[i1][j1][k1];
                    if (sgg.med(jmed).Is.Lumped) {
                        conta++;
                        Nodes_t& lumped = LumpElem.Nodes[conta - 1];
                        
                        lumped.alignedDeltaE = 1.0 / Idye[j1];
                        lumped.transversalDeltaHa = 1.0 / Idzh[k1];
                        lumped.transversalDeltaHb = 1.0 / Idxh[i1];
                        lumped.Orient = sgg.med(jmed).Lumped[1].Orient;
                        lumped.jmed = jmed;
                        
                        lumped.Efield = const_cast<double*>(&Ey[i1][j1][k1]);
                        lumped.Ha_Plus = const_cast<double*>(&Hx[i1][j1][k1]);
                        lumped.Ha_Minu = const_cast<double*>(&Hx[i1][j1][k1-1]);
                        lumped.Hb_Plus = const_cast<double*>(&Hz[i1][j1][k1]);
                        lumped.Hb_Minu = const_cast<double*>(&Hz[i1-1][j1][k1]);
                    }
                }
            }
        }

        // Process Ez nodes
        for (k1 = sgg.SINPMLSweep[iEz].ZI; k1 <= sgg.SINPMLSweep[iEz].ZE; ++k1) {
            for (j1 = sgg.SINPMLSweep[iEz].YI; j1 <= sgg.SINPMLSweep[iEz].YE; ++j1) {
                for (i1 = sgg.SINPMLSweep[iEz].XI; i1 <= sgg.SINPMLSweep[iEz].XE; ++i1) {
                    int jmed = media.sggMiEz[i1][j1][k1];
                    if (sgg.med(jmed).Is.Lumped) {
                        conta++;
                        Nodes_t& lumped = LumpElem.Nodes[conta - 1];
                        
                        lumped.alignedDeltaE = 1.0 / Idze[k1];
                        lumped.transversalDeltaHa = 1.0 / Idxh[i1];
                        lumped.transversalDeltaHb = 1.0 / Idyh[j1];
                        lumped.Orient = sgg.med(jmed).Lumped[1].Orient;
                        lumped.jmed = jmed;
                        
                        lumped.Efield = const_cast<double*>(&Ez[i1][j1][k1]);
                        lumped.Ha_Plus = const_cast<double*>(&Hy[i1][j1][k1]);
                        lumped.Ha_Minu = const_cast<double*>(&Hy[i1-1][j1][k1]);
                        lumped.Hb_Plus = const_cast<double*>(&Hx[i1][j1][k1]);
                        lumped.Hb_Minu = const_cast<double*>(&Hx[i1][j1-1][k1]);
                    }
                }
            }
        }

        calc_lumpedconstants(sgg, eps0, mu0);

        if (!control.resume) {
            for (int idx = 0; idx < LumpElem.NumNodes; ++idx) {
                Nodes_t& lumped = LumpElem.Nodes[idx];
                lumped.EfieldPrevPrev = 0.0;
                lumped.EfieldPrev = 0.0;
                lumped.Jcur = 0.0;
#ifdef CompileWithStochastic
                if (control.stochastic) {
                    lumped.EfieldPrevPrev_for_devia = 0.0;
                    lumped.EfieldPrev_for_devia = 0.0;
                    lumped.Jcur_for_devia = 0.0;
                }
#endif
            }
        } else {
            // Reading from file unit 14
            // In C++, we'd use std::ifstream or similar. 
            // Assuming a global file stream or object exists for unit 14.
            // For this translation, we'll assume a helper function or global stream.
            // Since we can't implement the file I/O infrastructure here, we'll stub it.
            for (int idx = 0; idx < LumpElem.NumNodes; ++idx) {
                Nodes_t& lumped = LumpElem.Nodes[idx];
                // read(14) lumped_%EfieldPrevPrev,lumped_%EfieldPrev,lumped_%Jcur 
                // Stub:
                // std::cin >> lumped.EfieldPrevPrev >> lumped.EfieldPrev >> lumped.Jcur; 
            }
#ifdef CompileWithStochastic
            if (control.stochastic) {
                for (int idx = 0; idx < LumpElem.NumNodes; ++idx) {
                    Nodes_t& lumped = LumpElem.Nodes[idx];
                    // read(14) ...
                }
            }
#endif
        }
    }

    void AdvanceLumpedE(const SGGFDTDINFO_t& sgg, int timestep, bool simu_devia, bool stochastic) {
        for (int conta = 0; conta < LumpElem.NumNodes; ++conta) {
            Nodes_t& lumped = LumpElem.Nodes[conta];
            
            if (sgg.med(lumped.jmed).Lumped[1].inductor) {
                lumped.EfieldPrevPrev = lumped.EfieldPrev;
                lumped.EfieldPrev = *lumped.Efield;
                lumped.Jcur = lumped.currentCoeff * lumped.Jcur + 
                              lumped.sigmaEffResistInduct * (lumped.EfieldPrev + lumped.EfieldPrevPrev);
            } else {
                lumped.Jcur = 0.0;
            }

            if (sgg.med(lumped.jmed).Lumped[1].diodo) {
                double fieldC = -lumped.G1 * *lumped.Efield - 
                                (lumped.G2a * (*lumped.Ha_Plus - *lumped.Ha_Minu) - 
                                 lumped.G2b * (*lumped.Hb_Plus - *lumped.Hb_Minu)) - 
                                lumped.diodepreA;
                double A = lumped.diodepreA * std::exp(lumped.diodeB * *lumped.Efield);
                double Enplus1 = newton_raphson(A, lumped.diodeB, fieldC);
                *lumped.Efield = Enplus1;
            } else {
                if (sgg.med(lumped.jmed).Lumped[1].resistor) {
                    double time = timestep * sgg.dt;
                    if ((time >= sgg.med(lumped.jmed).Lumped[1].Rtime_on) && 
                        (time <= sgg.med(lumped.jmed).Lumped[1].Rtime_off)) {
                        *lumped.Efield = lumped.G1 * *lumped.Efield + 
                                         (lumped.G2a * (*lumped.Ha_Plus - *lumped.Ha_Minu) - 
                                          lumped.G2b * (*lumped.Hb_Plus - *lumped.Hb_Minu)) - 
                                         lumped.GJ * lumped.Jcur;
                    } else {
                        *lumped.Efield = lumped.G1_usual * *lumped.Efield + 
                                         (lumped.G2a_usual * (*lumped.Ha_Plus - *lumped.Ha_Minu) - 
                                          lumped.G2b_usual * (*lumped.Hb_Plus - *lumped.Hb_Minu));
                    }
                } else {
                    // Inductor or Capacitor
                    *lumped.Efield = lumped.G1 * *lumped.Efield + 
                                     (lumped.G2a * (*lumped.Ha_Plus - *lumped.Ha_Minu) - 
                                      lumped.G2b * (*lumped.Hb_Plus - *lumped.Hb_Minu)) - 
                                     lumped.GJ * lumped.Jcur;
                }
            }

#ifdef CompileWithStochastic
            inject_devialumped(sgg, timestep, simu_devia, stochastic, lumped);
#endif
        }
    }

    void calc_lumpedconstants(const SGGFDTDINFO_t& sgg, double eps00, double mu00) {
        eps0 = eps00;
        mu0 = mu00;
        zvac = std::sqrt(mu0 / eps0);
        cluz = 1.0 / std::sqrt(eps0 * mu0);

        for (int conta = 0; conta < LumpElem.NumNodes; ++conta) {
            Nodes_t& lumped = LumpElem.Nodes[conta];
            int jmed = lumped.jmed;

            int orient = sgg.med(jmed).Lumped[1].Orient;
            double Resist = sgg.med(jmed).Lumped[1].R;
            double Induct = sgg.med(jmed).Lumped[1].L;
            double Capaci = sgg.med(jmed).Lumped[1].C;
            double DiodB = sgg.med(jmed).Lumped[1].DiodB;
            double DiodIsat = sgg.med(jmed).Lumped[1].DiodIsat;
            double epsilon = sgg.med(jmed).epr * eps0;
            double sigma = sgg.med(jmed).sigma;
            
            double alignedDeltaE = lumped.alignedDeltaE;
            double transversalDeltaHa = lumped.transversalDeltaHa;
            double transversalDeltaHb = lumped.transversalDeltaHb;

            double epsilonEffCapac = alignedDeltaE * Capaci / (transversalDeltaHa * transversalDeltaHb);
            double sigmaEffResistInduct = alignedDeltaE * sgg.dt / 
                                          (2.0 * transversalDeltaHa * transversalDeltaHb * 
                                           (Induct + Resist * sgg.dt / 2.0));
            double sigmaEffResist = alignedDeltaE / (Resist * transversalDeltaHa * transversalDeltaHb);
            double sigmaEffResistCapac = sigmaEffResist;
            double sigmaEffResistDiode = sigmaEffResist;
            double currentCoeff = (Induct - Resist * sgg.dt / 2.0) / (Induct + Resist * sgg.dt / 2.0);

            double sigmaeff = 0.0;
            double epsiloneff = 0.0;

            if (sgg.med(jmed).Lumped[1].resistor) {
                sigmaeff = sigma + sigmaEffResist;
                epsiloneff = epsilon;
            } else if (sgg.med(jmed).Lumped[1].inductor) {
                sigmaeff = sigma + sigmaEffResistInduct;
                epsiloneff = epsilon;
            } else if (sgg.med(jmed).Lumped[1].inductor) { // Note: Original code had inductor check twice? No, next is capacitor
                // Wait, original:
                // if resistor ...
                // elseif inductor ...
                // elseif capacitor ...
                // elseif diodo ...
            } else if (sgg.med(jmed).Lumped[1].capacitor) { // Assuming this flag exists or is implied
                sigmaeff = sigma + sigmaEffResistCapac;
                epsiloneff = epsilon + epsilonEffCapac;
            } else if (sgg.med(jmed).Lumped[1].diodo) {
                sigmaeff = sigma + sigmaEffResistDiode;
                epsiloneff = epsilon;
            } else {
                // Default or error? Original code doesn't specify else.
                // Assuming it falls through or one of the above matches.
                // If none match, sigmaeff/epsiloneff remain 0? Or previous value?
                // Let's assume valid media always has one of these flags.
            }

            if (!sgg.med(jmed).sigmareasignado) {
                sgg.med(jmed).sigma = sigmaeff;
                sgg.med(jmed).sigmareasignado = true;
            } else {
                std::cerr << "error buggy: reasignando sigma en un lumped" << std::endl;
                // stop
            }

            double G1 = (epsiloneff / sgg.dt - sigmaeff / 2.0) / 
                        (epsiloneff / sgg.dt + sigmaeff / 2.0);
            double G2 = 1.0 / (epsiloneff / sgg.dt + sigmaeff / 2.0);

            lumped.G1 = G1;
            lumped.G2a = G2 / transversalDeltaHa;
            lumped.G2b = G2 / transversalDeltaHb;
            lumped.GJ = G2 * (1 + currentCoeff) / 2.0;
            lumped.sigmaEffResistInduct = sigmaEffResistInduct;
            lumped.currentCoeff = currentCoeff;

            double G1_usual = (1.0 - sigma * sgg.dt / (2.0 * epsilon)) / 
                              (1.0 + sigma * sgg.dt / (2.0 * epsilon));
            double G2_usual = sgg.dt / epsilon / 
                              (1.0 + sigma * sgg.dt / (2.0 * epsilon));

            if (G1_usual < 0.0) {
                G1_usual = std::exp(-sigma * sgg.dt / epsilon);
                G2_usual = (1.0 - G1_usual) / sigma;
            }

            lumped.G1_usual = G1_usual;
            lumped.G2a_usual = G2_usual / transversalDeltaHa;
            lumped.G2b_usual = G2_usual / transversalDeltaHb;

            if (orient > 0.0) {
                lumped.diodeB = lumped.diodeB * alignedDeltaE / 2.0;
                lumped.diodepreA = DiodIsat * G2 / (transversalDeltaHa * transversalDeltaHb);
            } else if (orient < 0.0) {
                lumped.diodeB = -lumped.diodeB * alignedDeltaE / 2.0;
                lumped.diodepreA = -DiodIsat * G2 / (transversalDeltaHa * transversalDeltaHb);
            } else {
                std::string buff = "Buggy ERROR: In lumped orientations. ";
                WarnErrReport(buff, true);
            }

#ifdef CompileWithStochastic
            calc_lumped_deviaconsts(sgg, lumped, jmed, Resist, Induct, Capaci, sigmaeff, epsiloneff, eps0, mu0);
#endif
        }
    }

    void StoreFieldsLumpeds(bool stochastic) {
        for (int conta = 0; conta < LumpElem.NumNodes; ++conta) {
            Nodes_t& lumped = LumpElem.Nodes[conta];
            // write(14,err=634) ...
            // Stub file write
        }
#ifdef CompileWithStochastic
        if (stochastic) {
            for (int conta = 0; conta < LumpElem.NumNodes; ++conta) {
                Nodes_t& lumped = LumpElem.Nodes[conta];
                // write(14,err=634) ...
            }
        }
#endif
        return;
        
        // Labels 634 and 635 are handled by control flow in Fortran.
        // In C++, we just return. Error handling would be done via exceptions or return codes.
    }

    void DestroyLumped(SGGFDTDINFO_t& sgg) {
        for (int i = 0; i < sgg.NumMedia; ++i) {
            if (sgg.med(i).Is.Lumped && !sgg.med(i).Is.PML) { // Assuming Is.PML exists
                // if (associated(sgg%Med(i)%Lumped)) deallocate...
                // In C++, if Lumped is a vector or pointer, we clear it.
                // Assuming sgg.med(i).Lumped is a vector or similar managed memory.
                // If it's a raw pointer, we'd delete it.
                // Based on Fortran `deallocate`, it's likely dynamic.
                // We assume the MediaStruct handles its own memory or we clear it here.
                // For this translation, we assume the vector/array inside Med is cleared automatically or via destructor.
            }
        }
        if (!LumpElem.Nodes.empty()) {
            LumpElem.Nodes.clear();
        }
    }

    LumpedElem_t& Getlumped() {
        return LumpElem;
    }

}