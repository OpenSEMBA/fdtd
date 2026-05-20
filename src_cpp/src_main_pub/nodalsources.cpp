#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <memory>

// Assuming these types are defined in other headers based on the Fortran 'use' statements
// We need forward declarations or includes for:
// FDETYPES_m, Report_m, SGGFDTDINFO_t, NodalSource_t, XYZlimit_t, bounds_t

// Placeholder definitions for external types to make the code compile conceptually
// In a real scenario, these would come from "FDETYPES_m.hpp" and "Report_m.hpp"

using RKIND = double;
using INTEGERSIZEOFMEDIAMATRICES = int;
enum { BUFSIZE = 256 };

struct XYZlimit_t {
    int XI, XE, YI, YE, ZI, ZE;
};

struct NodalSourcePoint_t {
    double xc, yc, zc;
    int XI, XE, YI, YE, ZI, ZE;
};

struct FileData_t {
    std::string Name;
    int NumSamples;
    double deltaSamples;
    std::vector<double> Samples;
};

struct NodalSource_t {
    bool IsElec;
    bool IsHard;
    bool IsInitialValue;
    int numpuntos;
    std::vector<NodalSourcePoint_t> punto;
    FileData_t fichero;
};

struct SGGFDTDINFO_t {
    double dt;
    std::vector<double> tiempo;
    // Assuming Med is a vector of media objects
    struct Media_t {
        struct Is_t {
            bool PEC;
        } Is;
    };
    std::vector<Media_t> Med;
    std::vector<NodalSource_t>* NodalSource;
};

struct bounds_t {
    struct {
        int NX, NY, NZ, XI, XE, YI, YE, ZI, ZE;
    } sggMiEx, sggMiEy, sggMiEz, sggMiHx, sggMiHy, sggMiHz, Ex, Ey, Ez, Hx, Hy, Hz, dxh, dyh, dzh;
};

// Forward declaration for warning function
void WarnErrReport(const std::string& msg);

// Constants for sweep indices (assumed to be defined in FDETYPES_m or similar)
enum { iEx = 0, iEy = 1, iEz = 2, iHx = 3, iHy = 4, iHz = 5 };

namespace nodalsources_m {

    struct xyzlimit_singlescaled_t {
        int XI, XE, YI, YE, ZI, ZE;
        double amplitude;
    };

    struct NodalLocal_t {
        std::vector<double> evol;
        double deltaevol;
        int numus;
        xyzlimit_singlescaled_t punto;
        bool IsInitialValue;
    };

    struct nodsou_t {
        int NumHard = 0;
        int NumSoft = 0;
        std::vector<NodalLocal_t> nodHard;
        std::vector<NodalLocal_t> nodSoft;
    };

    // Global variables, previously save, target
    nodsou_t Nodal_Ex;
    nodsou_t Nodal_Ey;
    nodsou_t Nodal_Ez;
    nodsou_t Nodal_Hx;
    nodsou_t Nodal_Hy;
    nodsou_t Nodal_Hz;

    void InitnodalSources(const SGGFDTDINFO_t& sgg, int layoutnumber, int NumNodalSources, 
                          const std::vector<NodalSource_t>& sggNodalSource, 
                          const std::vector<XYZlimit_t>& sggSweep, 
                          bool& ThereAreNodalE, bool& ThereAreNodalH) {
        
        int numNodalSoft_Ex = 0;
        int numNodalSoft_Ey = 0;
        int numNodalSoft_Ez = 0;
        int numNodalSoft_Hx = 0;
        int numNodalSoft_Hy = 0;
        int numNodalSoft_Hz = 0;
        int numNodalHard_Ex = 0;
        int numNodalHard_Ey = 0;
        int numNodalHard_Ez = 0;
        int numNodalHard_Hx = 0;
        int numNodalHard_Hy = 0;
        int numNodalHard_Hz = 0;

        ThereAreNodalE = false;
        ThereAreNodalH = false;

        for (int j = 0; j < NumNodalSources; ++j) {
            if (sggNodalSource[j].IsElec) {
                for (int i = 0; i < sggNodalSource[j].numpuntos; ++i) {
                    if (sggNodalSource[j].punto[i].xc != 0.0) {
                        if (sggNodalSource[j].IsHard) {
                            numNodalHard_Ex++;
                        } else {
                            numNodalSoft_Ex++;
                        }
                    }
                    if (sggNodalSource[j].punto[i].yc != 0.0) {
                        if (sggNodalSource[j].IsHard) {
                            numNodalHard_Ey++;
                        } else {
                            numNodalSoft_Ey++;
                        }
                    }
                    if (sggNodalSource[j].punto[i].zc != 0.0) {
                        if (sggNodalSource[j].IsHard) {
                            numNodalHard_Ez++;
                        } else {
                            numNodalSoft_Ez++;
                        }
                    }
                }
            } else {
                for (int i = 0; i < sggNodalSource[j].numpuntos; ++i) {
                    if (sggNodalSource[j].punto[i].xc != 0.0) {
                        if (sggNodalSource[j].IsHard) {
                            numNodalHard_Hx++;
                        } else {
                            numNodalSoft_Hx++;
                        }
                    }
                    if (sggNodalSource[j].punto[i].yc != 0.0) {
                        if (sggNodalSource[j].IsHard) {
                            numNodalHard_Hy++;
                        } else {
                            numNodalSoft_Hy++;
                        }
                    }
                    if (sggNodalSource[j].punto[i].zc != 0.0) {
                        if (sggNodalSource[j].IsHard) {
                            numNodalHard_Hz++;
                        } else {
                            numNodalSoft_Hz++;
                        }
                    }
                }
            }
        }

        if (numNodalSoft_Ex + numNodalSoft_Ey + numNodalSoft_Ez != 0) {
            ThereAreNodalE = true;
            Nodal_Ex.nodSoft.resize(numNodalSoft_Ex);
            Nodal_Ey.nodSoft.resize(numNodalSoft_Ey);
            Nodal_Ez.nodSoft.resize(numNodalSoft_Ez);
        }
        if (numNodalHard_Ex + numNodalHard_Ey + numNodalHard_Ez != 0) {
            ThereAreNodalE = true;
            Nodal_Ex.nodHard.resize(numNodalHard_Ex);
            Nodal_Ey.nodHard.resize(numNodalHard_Ey);
            Nodal_Ez.nodHard.resize(numNodalHard_Ez);
        }
        if (numNodalSoft_Hx + numNodalSoft_Hy + numNodalSoft_Hz != 0) {
            ThereAreNodalH = true;
            Nodal_Hx.nodSoft.resize(numNodalSoft_Hx);
            Nodal_Hy.nodSoft.resize(numNodalSoft_Hy);
            Nodal_Hz.nodSoft.resize(numNodalSoft_Hz);
        }
        if (numNodalHard_Hx + numNodalHard_Hy + numNodalHard_Hz != 0) {
            ThereAreNodalH = true;
            Nodal_Hx.nodHard.resize(numNodalHard_Hx);
            Nodal_Hy.nodHard.resize(numNodalHard_Hy);
            Nodal_Hz.nodHard.resize(numNodalHard_Hz);
        }

        Nodal_Ex.NumHard = 0;
        Nodal_Ey.NumHard = 0;
        Nodal_Ez.NumHard = 0;
        Nodal_Hx.NumHard = 0;
        Nodal_Hy.NumHard = 0;
        Nodal_Hz.NumHard = 0;

        Nodal_Ex.NumSoft = 0;
        Nodal_Ey.NumSoft = 0;
        Nodal_Ez.NumSoft = 0;
        Nodal_Hx.NumSoft = 0;
        Nodal_Hy.NumSoft = 0;
        Nodal_Hz.NumSoft = 0;

        // Helper lambda to create nodal sources
        auto CreateNodal = [&](int layout, nodsou_t& dummy, const NodalSource_t& sggdummy, const XYZlimit_t& sweep, int index, double amplit) {
            if (sggdummy.IsHard) {
                dummy.NumHard++;
                dummy.nodHard[dummy.NumHard - 1].IsInitialValue = sggdummy.IsInitialValue;
                dummy.nodHard[dummy.NumHard - 1].punto.XI = std::max(sggdummy.punto[index].XI, sweep.XI);
                dummy.nodHard[dummy.NumHard - 1].punto.XE = std::min(sggdummy.punto[index].XE, sweep.XE);
                dummy.nodHard[dummy.NumHard - 1].punto.YI = std::max(sggdummy.punto[index].YI, sweep.YI);
                dummy.nodHard[dummy.NumHard - 1].punto.YE = std::min(sggdummy.punto[index].YE, sweep.YE);
                dummy.nodHard[dummy.NumHard - 1].punto.ZI = std::max(sggdummy.punto[index].ZI, sweep.ZI);
                dummy.nodHard[dummy.NumHard - 1].punto.ZE = std::min(sggdummy.punto[index].ZE, sweep.ZE);
                dummy.nodHard[dummy.NumHard - 1].punto.amplitude = amplit;
                dummy.nodHard[dummy.NumHard - 1].deltaevol = sggdummy.fichero.deltaSamples;
                
                if (dummy.nodHard[dummy.NumHard - 1].deltaevol > sgg.dt) {
                    std::string buff = "WARNING: " + sggdummy.fichero.Name + " undersampled by a factor " + std::to_string(dummy.nodHard[dummy.NumHard - 1].deltaevol / sgg.dt);
                    WarnErrReport(buff);
                }
                dummy.nodHard[dummy.NumHard - 1].numus = sggdummy.fichero.NumSamples;
                dummy.nodHard[dummy.NumHard - 1].evol = sggdummy.fichero.Samples;
            } else {
                dummy.NumSoft++;
                dummy.nodSoft[dummy.NumSoft - 1].IsInitialValue = sggdummy.IsInitialValue;
                dummy.nodSoft[dummy.NumSoft - 1].punto.XI = std::max(sggdummy.punto[index].XI, sweep.XI);
                dummy.nodSoft[dummy.NumSoft - 1].punto.XE = std::min(sggdummy.punto[index].XE, sweep.XE);
                dummy.nodSoft[dummy.NumSoft - 1].punto.YI = std::max(sggdummy.punto[index].YI, sweep.YI);
                dummy.nodSoft[dummy.NumSoft - 1].punto.YE = std::min(sggdummy.punto[index].YE, sweep.YE);
                dummy.nodSoft[dummy.NumSoft - 1].punto.ZI = std::max(sggdummy.punto[index].ZI, sweep.ZI);
                dummy.nodSoft[dummy.NumSoft - 1].punto.ZE = std::min(sggdummy.punto[index].ZE, sweep.ZE);
                dummy.nodSoft[dummy.NumSoft - 1].punto.amplitude = amplit;
                dummy.nodSoft[dummy.NumSoft - 1].deltaevol = sggdummy.fichero.deltaSamples;

                if (dummy.nodSoft[dummy.NumSoft - 1].deltaevol > sgg.dt) {
                    std::string buff = "WARNING: " + sggdummy.fichero.Name + " undersampled by a factor " + std::to_string(dummy.nodSoft[dummy.NumSoft - 1].deltaevol / sgg.dt);
                    WarnErrReport(buff);
                }
                dummy.nodSoft[dummy.NumSoft - 1].numus = sggdummy.fichero.NumSamples;
                dummy.nodSoft[dummy.NumSoft - 1].evol = sggdummy.fichero.Samples;
            }
        };

        for (int j = 0; j < NumNodalSources; ++j) {
            if (sggNodalSource[j].IsElec) {
                for (int i = 0; i < sggNodalSource[j].numpuntos; ++i) {
                    double amplit = sggNodalSource[j].punto[i].xc;
                    if (amplit != 0.0) {
                        CreateNodal(layoutnumber, Nodal_Ex, sggNodalSource[j], sggSweep[iEx], i, amplit);
                    }
                    amplit = sggNodalSource[j].punto[i].yc;
                    if (amplit != 0.0) {
                        CreateNodal(layoutnumber, Nodal_Ey, sggNodalSource[j], sggSweep[iEy], i, amplit);
                    }
                    amplit = sggNodalSource[j].punto[i].zc;
                    if (amplit != 0.0) {
                        CreateNodal(layoutnumber, Nodal_Ez, sggNodalSource[j], sggSweep[iEz], i, amplit);
                    }
                }
            } else {
                for (int i = 0; i < sggNodalSource[j].numpuntos; ++i) {
                    double amplit = sggNodalSource[j].punto[i].xc;
                    if (amplit != 0.0) {
                        CreateNodal(layoutnumber, Nodal_Hx, sggNodalSource[j], sggSweep[iHx], i, amplit);
                    }
                    amplit = sggNodalSource[j].punto[i].yc;
                    if (amplit != 0.0) {
                        CreateNodal(layoutnumber, Nodal_Hy, sggNodalSource[j], sggSweep[iHy], i, amplit);
                    }
                    amplit = sggNodalSource[j].punto[i].zc;
                    if (amplit != 0.0) {
                        CreateNodal(layoutnumber, Nodal_Hz, sggNodalSource[j], sggSweep[iHz], i, amplit);
                    }
                }
            }
        }
    }

    double evolucion(double t, const NodalLocal_t& dummy) {
        if (dummy.IsInitialValue) {
            if (static_cast<int>(t / dummy.deltaevol) != 0) {
                std::cout << "error en initial values. " << std::endl;
                // In C++, we might throw an exception or return a specific error code
                // For now, we'll just return 0 or handle it as per original intent
                return 0.0; 
            }
            return dummy.evol[0];
        }

        double deltaevol = dummy.deltaevol;
        const std::vector<double>& evol = dummy.evol;
        int numus = dummy.numus;

        double result = 0.0;
        int nprev = static_cast<int>(t / deltaevol);

        if ((nprev + 1 > numus) || (nprev + 1 <= 0)) {
            result = 0.0;
        } else {
            result = (evol[nprev + 1] - evol[nprev]) / deltaevol * (t - nprev * deltaevol) + evol[nprev];
        }
        return result;
    }

    void DestroyNodal(SGGFDTDINFO_t& sgg) {
        // Vectors are automatically cleared when resized to 0 or destroyed
        // But we need to clear the global variables
        Nodal_Ex.nodHard.clear();
        Nodal_Ex.nodSoft.clear();
        Nodal_Ey.nodHard.clear();
        Nodal_Ey.nodSoft.clear();
        Nodal_Ez.nodHard.clear();
        Nodal_Ez.nodSoft.clear();
        Nodal_Hx.nodHard.clear();
        Nodal_Hx.nodSoft.clear();
        Nodal_Hy.nodHard.clear();
        Nodal_Hy.nodSoft.clear();
        Nodal_Hz.nodHard.clear();
        Nodal_Hz.nodSoft.clear();

        Nodal_Ex.NumHard = 0;
        Nodal_Ex.NumSoft = 0;
        Nodal_Ey.NumHard = 0;
        Nodal_Ey.NumSoft = 0;
        Nodal_Ez.NumHard = 0;
        Nodal_Ez.NumSoft = 0;
        Nodal_Hx.NumHard = 0;
        Nodal_Hx.NumSoft = 0;
        Nodal_Hy.NumHard = 0;
        Nodal_Hy.NumSoft = 0;
        Nodal_Hz.NumHard = 0;
        Nodal_Hz.NumSoft = 0;

        if (sgg.NodalSource) {
            delete sgg.NodalSource;
            sgg.NodalSource = nullptr;
        }
    }

    void AdvancenodalE(const SGGFDTDINFO_t& sgg, 
                       const std::vector<std::vector<std::vector<IntegERSIZEOFMEDIAMATRICES>>>& sggMiEx,
                       const std::vector<std::vector<std::vector<IntegERSIZEOFMEDIAMATRICES>>>& sggMiEy,
                       const std::vector<std::vector<std::vector<IntegERSIZEOFMEDIAMATRICES>>>& sggMiEz,
                       int NumMedia, int timeinstant, 
                       const bounds_t& b,
                       const std::vector<double>& g2,
                       const std::vector<double>& Idxh,
                       const std::vector<double>& Idyh,
                       const std::vector<double>& Idzh,
                       std::vector<std::vector<std::vector<double>>>& Ex,
                       std::vector<std::vector<std::vector<double>>>& Ey,
                       std::vector<std::vector<std::vector<double>>>& Ez,
                       bool simu_devia) {
        
        double timei = sgg.tiempo[timeinstant];

        // Process Ex
        for (int ii = 0; ii < Nodal_Ex.NumHard; ++ii) {
            if (Nodal_Ex.nodHard[ii].IsInitialValue && (timeinstant != 0)) {
                continue;
            }
            double amp = Nodal_Ex.nodHard[ii].punto.amplitude;
            for (int k = Nodal_Ex.nodHard[ii].punto.ZI; k <= Nodal_Ex.nodHard[ii].punto.ZE; ++k) {
                int k_m = k - b.Ex.ZI;
                for (int j = Nodal_Ex.nodHard[ii].punto.YI; j <= Nodal_Ex.nodHard[ii].punto.YE; ++j) {
                    int j_m = j - b.Ex.YI;
                    for (int i = Nodal_Ex.nodHard[ii].punto.XI; i <= Nodal_Ex.nodHard[ii].punto.XE; ++i) {
                        int i_m = i - b.Ex.XI;
                        int medio = sggMiEx[i_m][j_m][k_m];
                        if (!simu_devia) {
                            if (!sgg.Med[medio].Is.PEC) {
                                Ex[i_m][j_m][k_m] = amp * evolucion(timei, Nodal_Ex.nodHard[ii]);
                            }
                        } else {
                            if (!sgg.Med[medio].Is.PEC) {
                                Ex[i_m][j_m][k_m] = 0.0;
                            }
                        }
                    }
                }
            }
        }

        for (int ii = 0; ii < Nodal_Ex.NumSoft; ++ii) {
            if (Nodal_Ex.nodSoft[ii].IsInitialValue && (timeinstant != 0)) {
                continue;
            }
            double amp = Nodal_Ex.nodSoft[ii].punto.amplitude;
            for (int k = Nodal_Ex.nodSoft[ii].punto.ZI; k <= Nodal_Ex.nodSoft[ii].punto.ZE; ++k) {
                int k_m = k - b.Ex.ZI;
                for (int j = Nodal_Ex.nodSoft[ii].punto.YI; j <= Nodal_Ex.nodSoft[ii].punto.YE; ++j) {
                    int j_m = j - b.Ex.YI;
                    for (int i = Nodal_Ex.nodSoft[ii].punto.XI; i <= Nodal_Ex.nodSoft[ii].punto.XE; ++i) {
                        int i_m = i - b.Ex.XI;
                        int medio = sggMiEx[i_m][j_m][k_m];
                        if (!simu_devia) {
                            if (!sgg.Med[medio].Is.PEC) {
                                Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m] - g2[medio] * Idyh[j_m] * Idzh[k_m] * amp * evolucion(timei, Nodal_Ex.nodSoft[ii]);
                            }
                        } else {
                            if (!sgg.Med[medio].Is.PEC) {
                                Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m];
                            }
                        }
                    }
                }
            }
        }

        // Process Ey
        for (int ii = 0; ii < Nodal_Ey.NumHard; ++ii) {
            if (Nodal_Ey.nodHard[ii].IsInitialValue && (timeinstant != 0)) {
                continue;
            }
            double amp = Nodal_Ey.nodHard[ii].punto.amplitude;
            for (int k = Nodal_Ey.nodHard[ii].punto.ZI; k <= Nodal_Ey.nodHard[ii].punto.ZE; ++k) {
                int k_m = k - b.Ey.ZI;
                for (int j = Nodal_Ey.nodHard[ii].punto.YI; j <= Nodal_Ey.nodHard[ii].punto.YE; ++j) {
                    int j_m = j - b.Ey.YI;
                    for (int i = Nodal_Ey.nodHard[ii].punto.XI; i <= Nodal_Ey.nodHard[ii].punto.XE; ++i) {
                        int i_m = i - b.Ey.XI;
                        int medio = sggMiEy[i_m][j_m][k_m];
                        if (!simu_devia) {
                            if (!sgg.Med[medio].Is.PEC) {
                                Ey[i_m][j_m][k_m] = amp * evolucion(timei, Nodal_Ey.nodHard[ii]);
                            }
                        } else {
                            if (!sgg.Med[medio].Is.PEC) {
                                Ey[i_m][j_m][k_m] = 0.0;
                            }
                        }
                    }
                }
            }
        }

        for (int ii = 0; ii < Nodal_Ey.NumSoft; ++ii) {
            if (Nodal_Ey.nodSoft[ii].IsInitialValue && (timeinstant != 0)) {
                continue;
            }
            double amp = Nodal_Ey.nodSoft[ii].punto.amplitude;
            for (int k = Nodal_Ey.nodSoft[ii].punto.ZI; k <= Nodal_Ey.nodSoft[ii].punto.ZE; ++k) {
                int k_m = k - b.Ey.ZI;
                for (int j = Nodal_Ey.nodSoft[ii].punto.YI; j <= Nodal_Ey.nodSoft[ii].punto.YE; ++j) {
                    int j_m = j - b.Ey.YI;
                    for (int i = Nodal_Ey.nodSoft[ii].punto.XI; i <= Nodal_Ey.nodSoft[ii].punto.XE; ++i) {
                        int i_m = i - b.Ey.XI;
                        int medio = sggMiEy[i_m][j_m][k_m];
                        if (!simu_devia) {
                            if (!sgg.Med[medio].Is.PEC) {
                                Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] - g2[medio] * Idxh[i_m] * Idzh[k_m] * amp * evolucion(timei, Nodal_Ey.nodSoft[ii]);
                            }
                        } else {
                            if (!sgg.Med[medio].Is.PEC) {
                                Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m];
                            }
                        }
                    }
                }
            }
        }

        // Process Ez
        for (int ii = 0; ii < Nodal_Ez.NumHard; ++ii) {
            if (Nodal_Ez.nodHard[ii].IsInitialValue && (timeinstant != 0)) {
                continue;
            }
            double amp = Nodal_Ez.nodHard[ii].punto.amplitude;
            for (int k = Nodal_Ez.nodHard[ii].punto.ZI; k <= Nodal_Ez.nodHard[ii].punto.ZE; ++k) {
                int k_m = k - b.Ez.ZI;
                for (int j = Nodal_Ez.nodHard[ii].punto.YI; j <= Nodal_Ez.nodHard[ii].punto.YE; ++j) {
                    int j_m = j - b.Ez.YI;
                    for (int i = Nodal_Ez.nodHard[ii].punto.XI; i <= Nodal_Ez.nodHard[ii].punto.XE; ++i) {
                        int i_m = i - b.Ez.XI;
                        int medio = sggMiEz[i_m][j_m][k_m];
                        if (!simu_devia) {
                            if (!sgg.Med[medio].Is.PEC) {
                                Ez[i_m][j_m][k_m] = amp * evolucion(timei, Nodal_Ez.nodHard[ii]);
                            }
                        } else {
                            if (!sgg.Med[medio].Is.PEC) {
                                Ez[i_m][j_m][k_m] = 0.0;
                            }
                        }
                    }
                }
            }
        }

        for (int ii = 0; ii < Nodal_Ez.NumSoft; ++ii) {
            if (Nodal_Ez.nodSoft[ii].IsInitialValue && (timeinstant != 0)) {
                continue;
            }
            double amp = Nodal_Ez.nodSoft[ii].punto.amplitude;
            for (int k = Nodal_Ez.nodSoft[ii].punto.ZI; k <= Nodal_Ez.nodSoft[ii].punto.ZE; ++k) {
                int k_m = k - b.Ez.ZI;
                for (int j = Nodal_Ez.nodSoft[ii].punto.YI; j <= Nodal_Ez.nodSoft[ii].punto.YE; ++j) {
                    int j_m = j - b.Ez.YI;
                    for (int i = Nodal_Ez.nodSoft[ii].punto.XI; i <= Nodal_Ez.nodSoft[ii].punto.XE; ++i) {
                        int i_m = i - b.Ez.XI;
                        int medio = sggMiEz[i_m][j_m][k_m];
                        if (!simu_devia) {
                            if (!sgg.Med[medio].Is.PEC) {
                                Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] - g2[medio] * Idyh[j_m] * Idxh[i_m] * amp * evolucion(timei, Nodal_Ez.nodSoft[ii]);
                            }
                        } else {
                            if (!sgg.Med[medio].Is.PEC) {
                                Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m];
                            }
                        }
                    }
                }
            }
        }
    }

    void AdvancenodalH(const SGGFDTDINFO_t& sgg,
                       const std::vector<std::vector<std::vector<IntegERSIZEOFMEDIAMATRICES>>>& sggMiHx,
                       const std::vector<std::vector<std::vector<IntegERSIZEOFMEDIAMATRICES>>>& sggMiHy,
                       const std::vector<std::vector<std::vector<IntegERSIZEOFMEDIAMATRICES>>>& sggMiHz,
                       int NumMedia, int timeinstant,
                       const bounds_t& b,
                       const std::vector<double>& gm2,
                       const std::vector<double>& Idxe,
                       const std::vector<double>& Idye,
                       const std::vector<double>& Idze,
                       std::vector<std::vector<std::vector<double>>>& Hx,
                       std::vector<std::vector<std::vector<double>>>& Hy,
                       std::vector<std::vector<std::vector<double>>>& Hz,
                       bool simu_devia) {
        
        double timei = sgg.tiempo[timeinstant];
        double GM2_1;

        if (simu_devia) {
            // The original code cuts off here, so we leave it empty or add a placeholder
            // Based on the pattern, it likely does similar loops for H fields
        }
    }

} // namespace nodalsources_m

std::cout << "Devia H nodal/field sources untested. Aborting" << std::endl;
            std::exit(1);
        }
        GM2_1 = GM2[0];
        
        // ---------------------------> empieza AdvancenodalH <---------------------------------------
        
        timei = sgg->tiempo(timeinstant) + 0.5_RKIND * sgg->dt;
        // !!!! deprecado en pscale y el+3 de la sincronia con ORIGINAL se jode para siempre 110219 
        // !!! timei = ( timeinstant + 0.5_RKIND  +3.0_RKIND) * sgg->dt  !ORIGINAL sync

        for (int ii = 0; ii < Nodal_Hx.numHard; ++ii) {
            if (Nodal_Hx.nodHard[ii].IsInitialValue && (timeinstant != 0)) {
                continue;
            }
            //     
            amp = Nodal_Hx.nodHard[ii].punto.amplitude;
            for (int k = Nodal_Hx.nodHard[ii].punto.zi; k <= Nodal_Hx.nodHard[ii].punto.ze; ++k) {
                k_m = k - b.Hx.ZI;
                for (int j = Nodal_Hx.nodHard[ii].punto.yi; j <= Nodal_Hx.nodHard[ii].punto.ye; ++j) {
                    j_m = j - b.Hx.YI;
                    for (int i = Nodal_Hx.nodHard[ii].punto.xi; i <= Nodal_Hx.nodHard[ii].punto.xe; ++i) {
                        i_m = i - b.Hx.XI;
                        medio = sggMiHx(i_m, j_m, k_m);
                        if (!sgg->Med(medio).Is.PMC) {
                            Hx(i_m, j_m, k_m) = amp * evolucion(timei, Nodal_Hx.nodHard[ii]);
                        }
                    }
                }
            }
        }
        //
        for (int ii = 0; ii < Nodal_Hx.numSoft; ++ii) {
            if (Nodal_Hx.nodSoft[ii].IsInitialValue && (timeinstant != 0)) {
                continue;
            }
            //     
            amp = Nodal_Hx.nodSoft[ii].punto.amplitude;
            for (int k = Nodal_Hx.nodSoft[ii].punto.zi; k <= Nodal_Hx.nodSoft[ii].punto.ze; ++k) {
                k_m = k - b.Hx.ZI;
                for (int j = Nodal_Hx.nodSoft[ii].punto.yi; j <= Nodal_Hx.nodSoft[ii].punto.ye; ++j) {
                    j_m = j - b.Hx.YI;
                    for (int i = Nodal_Hx.nodSoft[ii].punto.xi; i <= Nodal_Hx.nodSoft[ii].punto.xe; ++i) {
                        i_m = i - b.Hx.XI;
                        medio = sggMiHx(i_m, j_m, k_m);
                        if (!sgg->Med(medio).Is.PMC) {
                            Hx(i_m, j_m, k_m) = Hx(i_m, j_m, k_m) - Gm2(medio) * Idye(j_m) * Idze(k_m) * amp * evolucion(timei, Nodal_Hx.nodSoft[ii]);
                        }
                    }
                }
            }
        }
        //
        //
        for (int ii = 0; ii < Nodal_Hy.numHard; ++ii) {
            if (Nodal_Hy.nodHard[ii].IsInitialValue && (timeinstant != 0)) {
                continue;
            }
            //
            amp = Nodal_Hy.nodHard[ii].punto.amplitude;
            for (int k = Nodal_Hy.nodHard[ii].punto.zi; k <= Nodal_Hy.nodHard[ii].punto.ze; ++k) {
                k_m = k - b.Hy.ZI;
                for (int j = Nodal_Hy.nodHard[ii].punto.yi; j <= Nodal_Hy.nodHard[ii].punto.ye; ++j) {
                    j_m = j - b.Hy.YI;
                    for (int i = Nodal_Hy.nodHard[ii].punto.xi; i <= Nodal_Hy.nodHard[ii].punto.xe; ++i) {
                        i_m = i - b.Hy.XI;
                        medio = sggMiHx(i_m, j_m, k_m);
                        if (!sgg->Med(medio).Is.PMC) {
                            Hy(i_m, j_m, k_m) = amp * evolucion(timei, Nodal_Hy.nodHard[ii]);
                        }
                    }
                }
            }
        }
        //
        for (int ii = 0; ii < Nodal_Hy.numSoft; ++ii) {
            if (Nodal_Hy.nodSoft[ii].IsInitialValue && (timeinstant != 0)) {
                continue;
            }
            //     
            amp = Nodal_Hy.nodSoft[ii].punto.amplitude;
            for (int k = Nodal_Hy.nodSoft[ii].punto.zi; k <= Nodal_Hy.nodSoft[ii].punto.ze; ++k) {
                k_m = k - b.Hy.ZI;
                for (int j = Nodal_Hy.nodSoft[ii].punto.yi; j <= Nodal_Hy.nodSoft[ii].punto.ye; ++j) {
                    j_m = j - b.Hy.YI;
                    for (int i = Nodal_Hy.nodSoft[ii].punto.xi; i <= Nodal_Hy.nodSoft[ii].punto.xe; ++i) {
                        i_m = i - b.Hy.XI;
                        medio = sggMiHy(i_m, j_m, k_m);
                        if (!sgg->Med(medio).Is.PMC) {
                            Hy(i_m, j_m, k_m) = Hy(i_m, j_m, k_m) - Gm2(medio) * Idxe(i_m) * Idze(k_m) * amp * evolucion(timei, Nodal_Hy.nodSoft[ii]);
                        }
                    }
                }
            }
        }

        for (int ii = 0; ii < Nodal_Hz.numHard; ++ii) {
            if (Nodal_Hz.nodHard[ii].IsInitialValue && (timeinstant != 0)) {
                continue;
            }
            //
            amp = Nodal_Hz.nodHard[ii].punto.amplitude;
            for (int k = Nodal_Hz.nodHard[ii].punto.zi; k <= Nodal_Hz.nodHard[ii].punto.ze; ++k) {
                k_m = k - b.Hz.ZI;
                for (int j = Nodal_Hz.nodHard[ii].punto.yi; j <= Nodal_Hz.nodHard[ii].punto.ye; ++j) {
                    j_m = j - b.Hz.YI;
                    for (int i = Nodal_Hz.nodHard[ii].punto.xi; i <= Nodal_Hz.nodHard[ii].punto.xe; ++i) {
                        i_m = i - b.Hz.XI;
                        medio = sggMiHx(i_m, j_m, k_m);
                        if (!sgg->Med(medio).Is.PMC) {
                            Hz(i_m, j_m, k_m) = amp * evolucion(timei, Nodal_Hz.nodHard[ii]);
                        }
                    }
                }
            }
        }
        //
        for (int ii = 0; ii < Nodal_Hz.numSoft; ++ii) {
            if (Nodal_Hz.nodSoft[ii].IsInitialValue && (timeinstant != 0)) {
                continue;
            }
            //     
            amp = Nodal_Hz.nodSoft[ii].punto.amplitude;
            for (int k = Nodal_Hz.nodSoft[ii].punto.zi; k <= Nodal_Hz.nodSoft[ii].punto.ze; ++k) {
                k_m = k - b.Hz.ZI;
                for (int j = Nodal_Hz.nodSoft[ii].punto.yi; j <= Nodal_Hz.nodSoft[ii].punto.ye; ++j) {
                    j_m = j - b.Hz.YI;
                    for (int i = Nodal_Hz.nodSoft[ii].punto.xi; i <= Nodal_Hz.nodSoft[ii].punto.xe; ++i) {
                        i_m = i - b.Hz.XI;
                        medio = sggMiHz(i_m, j_m, k_m);
                        if (!sgg->Med(medio).Is.PMC) {
                            Hz(i_m, j_m, k_m) = Hz(i_m, j_m, k_m) - Gm2(medio) * Idye(j_m) * Idxe(i_m) * amp * evolucion(timei, Nodal_Hz.nodSoft[ii]);
                        }
                    }
                }
            }
        }

        return;
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Function to publish the private output data (used in postprocess)
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!

    void getnodal(nodsou_t*& rNodal_Ex, nodsou_t*& rNodal_Ey, nodsou_t*& rNodal_Ez, nodsou_t*& rNodal_Hx, nodsou_t*& rNodal_Hy, nodsou_t*& rNodal_Hz) {
        rNodal_Ex = Nodal_Ex;
        rNodal_Ey = Nodal_Ey;
        rNodal_Ez = Nodal_Ez;
        rNodal_Hx = Nodal_Hx;
        rNodal_Hy = Nodal_Hy;
        rNodal_Hz = Nodal_Hz;
    }

} // namespace nodalsources_m