#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdint>

// Forward declarations and includes for external types used in the Fortran code
// These would typically come from FDETYPES_m and Report_m
// Assuming these types exist in the global namespace or a specific namespace

// Placeholder for types from FDETYPES_m
struct SGGFDTDINFO_t;
struct NodalSource_t;
struct XYZlimit_t;

// Placeholder for types from Report_m
void WarnErrReport(const std::string& msg);

// Constants assumed from FDETYPES_m
constexpr int32_t RKIND = 8; // Assuming double precision
constexpr int32_t INTEGERSIZEOFMEDIAMATRICES = 4; // Assuming 32-bit integers
constexpr int BUFSIZE = 256; // Assumed buffer size

// Enumerations assumed from context
enum Direction { iEx, iEy, iEz, iHx, iHy, iHz };

namespace nodalsources_m {

    struct xyzlimit_singlescaled_t {
        int32_t XI, XE, YI, YE, ZI, ZE;
        double amplitude;
    };

    struct NodalLocal_t {
        std::vector<double> evol;
        double deltaevol;
        int32_t numus;
        xyzlimit_singlescaled_t punto;
        bool IsInitialValue;
    };

    struct nodsou_t {
        int32_t NumHard = 0;
        int32_t NumSoft = 0;
        std::vector<NodalLocal_t> nodHard;
        std::vector<NodalLocal_t> nodSoft;
    };

    // Global variables
    nodsou_t Nodal_Ex;
    nodsou_t Nodal_Ey;
    nodsou_t Nodal_Ez;
    nodsou_t Nodal_Hx;
    nodsou_t Nodal_Hy;
    nodsou_t Nodal_Hz;

    // Helper function for CreateNodal logic (extracted to avoid deep nesting)
    void CreateNodalHelper(nodsou_t& dummy, const NodalSource_t& sggdummy, const XYZlimit_t& sggSweep, int32_t index, double amplit, bool isHard, SGGFDTDINFO_t& sgg) {
        if (isHard) {
            dummy.nodHard.push_back(NodalLocal_t());
            int32_t idx = dummy.nodHard.size() - 1;
            dummy.NumHard++;
            
            dummy.nodHard[idx].IsInitialValue = sggdummy.IsInitialValue;
            dummy.nodHard[idx].punto.XI = std::max(sggdummy.punto[index].XI, sggSweep.XI);
            dummy.nodHard[idx].punto.XE = std::min(sggdummy.punto[index].XE, sggSweep.XE);
            dummy.nodHard[idx].punto.YI = std::max(sggdummy.punto[index].YI, sggSweep.YI);
            dummy.nodHard[idx].punto.YE = std::min(sggdummy.punto[index].YE, sggSweep.YE);
            dummy.nodHard[idx].punto.ZI = std::max(sggdummy.punto[index].ZI, sggSweep.ZI);
            dummy.nodHard[idx].punto.ZE = std::min(sggdummy.punto[index].ZE, sggSweep.ZE);
            
            dummy.nodHard[idx].punto.amplitude = amplit;
            dummy.nodHard[idx].deltaevol = sggdummy.fichero.deltaSamples;
            
            if (dummy.nodHard[idx].deltaevol > sgg.dt) {
                std::string buff = "WARNING: " + sggdummy.Fichero.Name + " undersampled by a factor " + std::to_string(dummy.nodHard[idx].deltaevol / sgg.dt);
                WarnErrReport(buff);
            }
            
            dummy.nodHard[idx].numus = sggdummy.Fichero.NumSamples;
            dummy.nodHard[idx].evol = sggdummy.fichero.Samples;
        } else {
            dummy.nodSoft.push_back(NodalLocal_t());
            int32_t idx = dummy.nodSoft.size() - 1;
            dummy.NumSoft++;
            
            dummy.nodSoft[idx].IsInitialValue = sggdummy.IsInitialValue;
            dummy.nodSoft[idx].punto.XI = std::max(sggdummy.punto[index].XI, sggSweep.XI);
            dummy.nodSoft[idx].punto.XE = std::min(sggdummy.punto[index].XE, sggSweep.XE);
            dummy.nodSoft[idx].punto.YI = std::max(sggdummy.punto[index].YI, sggSweep.YI);
            dummy.nodSoft[idx].punto.YE = std::min(sggdummy.punto[index].YE, sggSweep.YE);
            dummy.nodSoft[idx].punto.ZI = std::max(sggdummy.punto[index].ZI, sggSweep.ZI);
            dummy.nodSoft[idx].punto.ZE = std::min(sggdummy.punto[index].ZE, sggSweep.ZE);
            
            dummy.nodSoft[idx].punto.amplitude = amplit;
            dummy.nodSoft[idx].deltaevol = sggdummy.fichero.deltaSamples;
            
            if (dummy.nodSoft[idx].deltaevol > sgg.dt) {
                std::string buff = "WARNING: " + sggdummy.Fichero.Name + " undersampled by a factor " + std::to_string(dummy.nodSoft[idx].deltaevol / sgg.dt);
                WarnErrReport(buff);
            }
            
            dummy.nodSoft[idx].numus = sggdummy.Fichero.NumSamples;
            dummy.nodSoft[idx].evol = sggdummy.fichero.Samples;
        }
    }

    void InitnodalSources(const SGGFDTDINFO_t& sgg, int32_t layoutnumber, int32_t NumNodalSources, const std::vector<NodalSource_t>& sggNodalSource, const std::vector<XYZlimit_t>& sggSweep, bool& ThereArenodalE, bool& ThereArenodalH) {
        int32_t numNodalSoft_Ex = 0, numNodalSoft_Ey = 0, numNodalSoft_Ez = 0;
        int32_t numNodalSoft_Hx = 0, numNodalSoft_Hy = 0, numNodalSoft_Hz = 0;
        int32_t numNodalHard_Ex = 0, numNodalHard_Ey = 0, numNodalHard_Ez = 0;
        int32_t numNodalHard_Hx = 0, numNodalHard_Hy = 0, numNodalHard_Hz = 0;

        ThereArenodalE = false;
        ThereArenodalH = false;

        for (int32_t j = 0; j < NumNodalSources; ++j) {
            if (sggNodalSource[j].IsElec) {
                for (int32_t i = 0; i < sggNodalSource[j].numpuntos; ++i) {
                    if (sggNodalSource[j].punto[i].xc != 0.0) {
                        if (sggNodalSource[j].IsHard) numNodalHard_Ex++;
                        else numNodalSoft_Ex++;
                    }
                    if (sggNodalSource[j].punto[i].yc != 0.0) {
                        if (sggNodalSource[j].IsHard) numNodalHard_Ey++;
                        else numNodalSoft_Ey++;
                    }
                    if (sggNodalSource[j].punto[i].zc != 0.0) {
                        if (sggNodalSource[j].IsHard) numNodalHard_Ez++;
                        else numNodalSoft_Ez++;
                    }
                }
            } else {
                for (int32_t i = 0; i < sggNodalSource[j].numpuntos; ++i) {
                    if (sggNodalSource[j].punto[i].xc != 0.0) {
                        if (sggNodalSource[j].IsHard) numNodalHard_Hx++;
                        else numNodalSoft_Hx++;
                    }
                    if (sggNodalSource[j].punto[i].yc != 0.0) {
                        if (sggNodalSource[j].IsHard) numNodalHard_Hy++;
                        else numNodalSoft_Hy++;
                    }
                    if (sggNodalSource[j].punto[i].zc != 0.0) {
                        if (sggNodalSource[j].IsHard) numNodalHard_Hz++;
                        else numNodalSoft_Hz++;
                    }
                }
            }
        }

        if (numNodalSoft_Ex + numNodalSoft_Ey + numNodalSoft_Ez != 0) {
            ThereArenodalE = true;
            Nodal_Ex.nodSoft.resize(numNodalSoft_Ex);
            Nodal_Ey.nodSoft.resize(numNodalSoft_Ey);
            Nodal_Ez.nodSoft.resize(numNodalSoft_Ez);
        }
        if (numNodalHard_Ex + numNodalHard_Ey + numNodalHard_Ez != 0) {
            ThereArenodalE = true;
            Nodal_Ex.nodHard.resize(numNodalHard_Ex);
            Nodal_Ey.nodHard.resize(numNodalHard_Ey);
            Nodal_Ez.nodHard.resize(numNodalHard_Ez);
        }
        if (numNodalSoft_Hx + numNodalSoft_Hy + numNodalSoft_Hz != 0) {
            ThereArenodalH = true;
            Nodal_Hx.nodSoft.resize(numNodalSoft_Hx);
            Nodal_Hy.nodSoft.resize(numNodalSoft_Hy);
            Nodal_Hz.nodSoft.resize(numNodalSoft_Hz);
        }
        if (numNodalHard_Hx + numNodalHard_Hy + numNodalHard_Hz != 0) {
            ThereArenodalH = true;
            Nodal_Hx.nodHard.resize(numNodalHard_Hx);
            Nodal_Hy.nodHard.resize(numNodalHard_Hy);
            Nodal_Hz.nodHard.resize(numNodalHard_Hz);
        }

        Nodal_Ex.NumHard = 0; Nodal_Ey.NumHard = 0; Nodal_Ez.NumHard = 0;
        Nodal_Hx.NumHard = 0; Nodal_Hy.NumHard = 0; Nodal_Hz.NumHard = 0;
        Nodal_Ex.NumSoft = 0; Nodal_Ey.NumSoft = 0; Nodal_Ez.NumSoft = 0;
        Nodal_Hx.NumSoft = 0; Nodal_Hy.NumSoft = 0; Nodal_Hz.NumSoft = 0;

        for (int32_t j = 0; j < NumNodalSources; ++j) {
            if (sggNodalSource[j].IsElec) {
                for (int32_t i = 0; i < sggNodalSource[j].numpuntos; ++i) {
                    double amplit = sggNodalSource[j].punto[i].xc;
                    if (amplit != 0.0) {
                        CreateNodalHelper(Nodal_Ex, sggNodalSource[j], sggSweep[iEx], i, amplit, sggNodalSource[j].IsHard, const_cast<SGGFDTDINFO_t&>(sgg));
                    }
                    amplit = sggNodalSource[j].punto[i].yc;
                    if (amplit != 0.0) {
                        CreateNodalHelper(Nodal_Ey, sggNodalSource[j], sggSweep[iEy], i, amplit, sggNodalSource[j].IsHard, const_cast<SGGFDTDINFO_t&>(sgg));
                    }
                    amplit = sggNodalSource[j].punto[i].zc;
                    if (amplit != 0.0) {
                        CreateNodalHelper(Nodal_Ez, sggNodalSource[j], sggSweep[iEz], i, amplit, sggNodalSource[j].IsHard, const_cast<SGGFDTDINFO_t&>(sgg));
                    }
                }
            } else {
                for (int32_t i = 0; i < sggNodalSource[j].numpuntos; ++i) {
                    double amplit = sggNodalSource[j].punto[i].xc;
                    if (amplit != 0.0) {
                        CreateNodalHelper(Nodal_Hx, sggNodalSource[j], sggSweep[iHx], i, amplit, sggNodalSource[j].IsHard, const_cast<SGGFDTDINFO_t&>(sgg));
                    }
                    amplit = sggNodalSource[j].punto[i].yc;
                    if (amplit != 0.0) {
                        CreateNodalHelper(Nodal_Hy, sggNodalSource[j], sggSweep[iHy], i, amplit, sggNodalSource[j].IsHard, const_cast<SGGFDTDINFO_t&>(sgg));
                    }
                    amplit = sggNodalSource[j].punto[i].zc;
                    if (amplit != 0.0) {
                        CreateNodalHelper(Nodal_Hz, sggNodalSource[j], sggSweep[iHz], i, amplit, sggNodalSource[j].IsHard, const_cast<SGGFDTDINFO_t&>(sgg));
                    }
                }
            }
        }
    }

    double evolucion(double t, const NodalLocal_t& dummy) {
        if (dummy.IsInitialValue) {
            if (static_cast<int32_t>(t / dummy.deltaevol) != 0) {
                std::cerr << "error en initial values." << std::endl;
                // In C++, we might throw an exception or handle differently, but preserving 'stop' behavior via exit or assert
                // For translation purposes, we'll just return 0 or handle as error
                return 0.0; 
            }
            return dummy.evol[0];
        }

        double deltaevol = dummy.deltaevol;
        const std::vector<double>& evol = dummy.evol;
        int32_t numus = dummy.numus;

        double evolucion_val = 0.0;
        int32_t nprev = static_cast<int32_t>(t / deltaevol);

        if ((nprev + 1 > numus) || (nprev + 1 <= 0)) {
            evolucion_val = 0.0;
        } else {
            evolucion_val = (evol[nprev + 1] - evol[nprev]) / deltaevol * (t - nprev * deltaevol) + evol[nprev];
        }

        return evolucion_val;
    }

    void DestroyNodal(SGGFDTDINFO_t& sgg) {
        if (Nodal_Ex.NumSoft + Nodal_Ey.NumSoft + Nodal_Ez.NumSoft != 0) {
            Nodal_Ex.nodSoft.clear();
            Nodal_Ey.nodSoft.clear();
            Nodal_Ez.nodSoft.clear();
        }
        if (Nodal_Ex.NumHard + Nodal_Ey.NumHard + Nodal_Ez.NumHard != 0) {
            Nodal_Ex.nodHard.clear();
            Nodal_Ey.nodHard.clear();
            Nodal_Ez.nodHard.clear();
        }
        if (Nodal_Hx.NumSoft + Nodal_Hy.NumSoft + Nodal_Hz.NumSoft != 0) {
            Nodal_Hx.nodSoft.clear();
            Nodal_Hy.nodSoft.clear();
            Nodal_Hz.nodSoft.clear();
        }
        if (Nodal_Hx.NumHard + Nodal_Hy.NumHard + Nodal_Hz.NumHard != 0) {
            Nodal_Hx.nodHard.clear();
            Nodal_Hy.nodHard.clear();
            Nodal_Hz.nodHard.clear();
        }

        if (!sgg.NodalSource.empty()) {
            sgg.NodalSource.clear();
        }
    }

    template <typename T>
    void AdvancenodalE_impl(const SGGFDTDINFO_t& sgg, 
                            const std::vector<std::vector<std::vector<T>>>& sggMiEx,
                            const std::vector<std::vector<std::vector<T>>>& sggMiEy,
                            const std::vector<std::vector<std::vector<T>>>& sggMiEz,
                            int32_t NumMedia, int32_t timeinstant, 
                            const struct bounds_t& b,
                            const std::vector<double>& g2,
                            const std::vector<double>& Idyh,
                            const std::vector<double>& Idzh,
                            std::vector<std::vector<std::vector<double>>>& Ex,
                            std::vector<std::vector<std::vector<double>>>& Ey,
                            std::vector<std::vector<std::vector<double>>>& Ez,
                            bool simu_devia) {
        
        double timei = sgg.tiempo[timeinstant];
        double amp;
        int32_t i, j, k, i_m, j_m, k_m, ii, medio;

        // Ex Hard
        for (ii = 0; ii < Nodal_Ex.NumHard; ++ii) {
            if (Nodal_Ex.nodHard[ii].IsInitialValue && timeinstant != 0) continue;
            amp = Nodal_Ex.nodHard[ii].punto.amplitude;
            for (k = Nodal_Ex.nodHard[ii].punto.ZI; k <= Nodal_Ex.nodHard[ii].punto.ZE; ++k) {
                k_m = k - b.Ex.ZI;
                for (j = Nodal_Ex.nodHard[ii].punto.YI; j <= Nodal_Ex.nodHard[ii].punto.YE; ++j) {
                    j_m = j - b.Ex.YI;
                    for (i = Nodal_Ex.nodHard[ii].punto.XI; i <= Nodal_Ex.nodHard[ii].punto.XE; ++i) {
                        i_m = i - b.Ex.XI;
                        medio = sggMiEx[i_m][j_m][k_m];
                        if (!simu_devia) {
                            if (!sgg.Med(medio).Is.PEC) {
                                Ex[i_m][j_m][k_m] = amp * evolucion(timei, Nodal_Ex.nodHard[ii]);
                            }
                        } else {
                            if (!sgg.Med(medio).Is.PEC) {
                                Ex[i_m][j_m][k_m] = 0.0;
                            }
                        }
                    }
                }
            }
        }

        // Ex Soft
        for (ii = 0; ii < Nodal_Ex.NumSoft; ++ii) {
            if (Nodal_Ex.nodSoft[ii].IsInitialValue && timeinstant != 0) continue;
            amp = Nodal_Ex.nodSoft[ii].punto.amplitude;
            for (k = Nodal_Ex.nodSoft[ii].punto.ZI; k <= Nodal_Ex.nodSoft[ii].punto.ZE; ++k) {
                k_m = k - b.Ex.ZI;
                for (j = Nodal_Ex.nodSoft[ii].punto.YI; j <= Nodal_Ex.nodSoft[ii].punto.YE; ++j) {
                    j_m = j - b.Ex.YI;
                    for (i = Nodal_Ex.nodSoft[ii].punto.XI; i <= Nodal_Ex.nodSoft[ii].punto.XE; ++i) {
                        i_m = i - b.Ex.XI;
                        medio = sggMiEx[i_m][j_m][k_m];
                        if (!simu_devia) {
                            if (!sgg.Med(medio).Is.PEC) {
                                Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m] - g2[medio] * Idyh[j_m] * Idzh[k_m] * amp * evolucion(timei, Nodal_Ex.nodSoft[ii]);
                            }
                        } else {
                            if (!sgg.Med(medio).Is.PEC) {
                                Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m];
                            }
                        }
                    }
                }
            }
        }

        // Ey Hard
        for (ii = 0; ii < Nodal_Ey.NumHard; ++ii) {
            if (Nodal_Ey.nodHard[ii].IsInitialValue && timeinstant != 0) continue;
            amp = Nodal_Ey.nodHard[ii].punto.amplitude;
            for (k = Nodal_Ey.nodHard[ii].punto.ZI; k <= Nodal_Ey.nodHard[ii].punto.ZE; ++k) {
                k_m = k - b.Ey.ZI;
                for (j = Nodal_Ey.nodHard[ii].punto.YI; j <= Nodal_Ey.nodHard[ii].punto.YE; ++j) {
                    j_m = j - b.Ey.YI;
                    for (i = Nodal_Ey.nodHard[ii].punto.XI; i <= Nodal_Ey.nodHard[ii].punto.XE; ++i) {
                        i_m = i - b.Ey.XI;
                        medio = sggMiEy[i_m][j_m][k_m];
                        if (!simu_devia) {
                            if (!sgg.Med(medio).Is.PEC) {
                                Ey[i_m][j_m][k_m] = amp * evolucion(timei, Nodal_Ey.nodHard[ii]);
                            }
                        } else {
                            if (!sgg.Med(medio).Is.PEC) {
                                Ey[i_m][j_m][k_m] = 0.0;
                            }
                        }
                    }
                }
            }
        }

        // Ey Soft
        for (ii = 0; ii < Nodal_Ey.NumSoft; ++ii) {
            if (Nodal_Ey.nodSoft[ii].IsInitialValue && timeinstant != 0) continue;
            amp = Nodal_Ey.nodSoft[ii].punto.amplitude;
            for (k = Nodal_Ey.nodSoft[ii].punto.ZI; k <= Nodal_Ey.nodSoft[ii].punto.ZE; ++k) {
                k_m = k - b.Ey.ZI;
                for (j = Nodal_Ey.nodSoft[ii].punto.YI; j <= Nodal_Ey.nodSoft[ii].punto.YE; ++j) {
                    j_m = j - b.Ey.YI;
                    for (i = Nodal_Ey.nodSoft[ii].punto.XI; i <= Nodal_Ey.nodSoft[ii].punto.XE; ++i) {
                        i_m = i - b.Ey.XI;
                        medio = sggMiEy[i_m][j_m][k_m];
                        if (!simu_devia) {
                            if (!sgg.Med(medio).Is.PEC) {
                                Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] - g2[medio] * Idyh[j_m] * Idzh[k_m] * amp * evolucion(timei, Nodal_Ey.nodSoft[ii]);
                            }
                        } else {
                            if (!sgg.Med(medio).Is.PEC) {
                                Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m];
                            }
                        }
                    }
                }
            }
        }

        // Ez Hard
        for (ii = 0; ii < Nodal_Ez.NumHard; ++ii) {
            if (Nodal_Ez.nodHard[ii].IsInitialValue && timeinstant != 0) continue;
            amp = Nodal_Ez.nodHard[ii].punto.amplitude;
            for (k = Nodal_Ez.nodHard[ii].punto.ZI; k <= Nodal_Ez.nodHard[ii].punto.ZE; ++k) {
                k_m = k - b.Ez.ZI;
                for (j = Nodal_Ez.nodHard[ii].punto.YI; j <= Nodal_Ez.nodHard[ii].punto.YE; ++j) {
                    j_m = j - b.Ez.YI;
                    for (i = Nodal_Ez.nodHard[ii].punto.XI; i <= Nodal_Ez.nodHard[ii].punto.XE; ++i) {
                        i_m = i - b.Ez.XI;
                        medio = sggMiEz[i_m][j_m][k_m];
                        if (!simu_devia) {
                            if (!sgg.Med(medio).Is.PEC) {
                                Ez[i_m][j_m][k_m] = amp * evolucion(timei, Nodal_Ez.nodHard[ii]);
                            }
                        } else {
                            if (!sgg.Med(medio).Is.PEC) {
                                Ez[i_m][j_m][k_m] = 0.0;
                            }
                        }
                    }
                }
            }
        }

        // Ez Soft
        for (ii = 0; ii < Nodal_Ez.NumSoft; ++ii) {
            if (Nodal_Ez.nodSoft[ii].IsInitialValue && timeinstant != 0) continue;
            amp = Nodal_Ez.nodSoft[ii].punto.amplitude;
            for (k = Nodal_Ez.nodSoft[ii].punto.ZI; k <= Nodal_Ez.nodSoft[ii].punto.ZE; ++k) {
                k_m = k - b.Ez.ZI;
                for (j = Nodal_Ez.nodSoft[ii].punto.YI; j <= Nodal_Ez.nodSoft[ii].punto.YE; ++j) {
                    j_m = j - b.Ez.YI;
                    for (i = Nodal_Ez.nodSoft[ii].punto.XI; i <= Nodal_Ez.nodSoft[ii].punto.XE; ++i) {
                        i_m = i - b.Ez.XI;
                        medio = sggMiEz[i_m][j_m][k_m];
                        if (!simu_devia) {
                            if (!sgg.Med(medio).Is.PEC) {
                                Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] - g2[medio] * Idyh[j_m] * Idzh[k_m] * amp * evolucion(timei, Nodal_Ez.nodSoft[ii]);
                            }
                        } else {
                            if (!sgg.Med(medio).Is.PEC) {
                                Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m];
                            }
                        }
                    }
                }
            }
        }
    }

    void AdvancenodalE(const SGGFDTDINFO_t& sgg, 
                       const std::vector<std::vector<std::vector<int32_t>>>& sggMiEx,
                       const std::vector<std::vector<std::vector<int32_t>>>& sggMiEy,
                       const std::vector<std::vector<std::vector<int32_t>>>& sggMiEz,
                       int32_t NumMedia, int32_t timeinstant, 
                       const struct bounds_t& b,
                       const std::vector<double>& g2,
                       const std::vector<double>& Idyh,
                       const std::vector<double>& Idzh,
                       std::vector<std::vector<std::vector<double>>>& Ex,
                       std::vector<std::vector<std::vector<double>>>& Ey,
                       std::vector<std::vector<std::vector<double>>>& Ez,
                       bool simu_devia) {
        AdvancenodalE_impl(sgg, sggMiEx, sggMiEy, sggMiEz, NumMedia, timeinstant, b, g2, Idyh, Idzh, Ex, Ey, Ez, simu_devia);
    }

    void AdvancenodalH(const SGGFDTDINFO_t& sgg,
                       const std::vector<std::vector<std::vector<int32_t>>>& sggMiHx,
                       const std::vector<std::vector<std::vector<int32_t>>>& sggMiHy,
                       const std::vector<std::vector<std::vector<int32_t>>>& sggMiHz,
                       int32_t NumMedia, int32_t timeinstant,
                       const struct bounds_t& b,
                       const std::vector<double>& gm2,
                       const std::vector<double>& Idxe,
                       const std::vector<double>& Idye,
                       const std::vector<double>& Idze,
                       std::vector<std::vector<std::vector<double>>>& Hx,
                       std::vector<std::vector<std::vector<double>>>& Hy,
                       std::vector<std::vector<std::vector<double>>>& Hz,
                       bool simu_devia) {
        // Implementation truncated in source, providing empty stub or similar structure
        // Note: The source code cuts off inside this function.
        // Assuming similar logic to AdvancenodalE but for H fields.
        // Since the source is incomplete, we cannot provide a full translation.
        // This is a placeholder to satisfy compilation requirements based on the provided snippet.
    }

} // namespace nodalsources_m

#include <iostream>
#include <string>
#include <vector>
#include <memory>

// Forward declarations and includes for types used in this chunk
// Assuming these are defined in previous chunks or headers
// struct sgg_t;
// struct b_t;
// struct nodsou_t;
// std::vector<sgg_t> sgg;
// b_t b;
// nodsou_t Nodal_Ex, Nodal_Ey, Nodal_Ez, Nodal_Hx, Nodal_Hy, Nodal_Hz;
// std::vector<double> Gm2;
// std::vector<double> Idye;
// std::vector<double> Idze;
// std::vector<double> Idxe;
// std::vector<int> sggMiHx;
// std::vector<int> sggMiHy;
// std::vector<int> sggMiHz;
// double Hx[], Hy[], Hz[]; // Assuming these are global or passed contextually
// double evolucion(double time, const nodsou_t& source);
// int RKIND; // Assuming this is a constant or type indicator

namespace nodalsources_m {

    void AdvancenodalH(int timeinstant, int N, const sgg_t& sgg, const b_t& b, 
                       double Hx[], double Hy[], double Hz[], 
                       const std::vector<double>& Gm2, 
                       const std::vector<double>& Idye, 
                       const std::vector<double>& Idze, 
                       const std::vector<double>& Idxe,
                       const std::vector<int>& sggMiHx,
                       const std::vector<int>& sggMiHy,
                       const std::vector<int>& sggMiHz) {
        
        // Note: In the original Fortran, Hx, Hy, Hz, Gm2, Idye, etc., seem to be 
        // accessed via global variables or a common block. Here we pass them or assume 
        // they are accessible. The translation assumes a context where these are available.
        // For strict translation, we might need to pass them or access a global state object.
        // Given the "preserve names" rule, we assume these variables exist in the scope 
        // or are members of a class/module. Let's assume they are passed or global for this snippet.
        
        // Re-defining variables based on Fortran context if not passed:
        // Since this is a chunk, we assume the necessary context (sgg, b, arrays) is available.
        // We will use references to the global/module state if this were a class method, 
        // or pass them. Let's assume a class context or global state for simplicity in translation 
        // of this specific block, but strictly speaking, we should translate the signature.
        // However, the prompt asks to translate the chunk. The chunk calls a subroutine.
        // The subroutine signature in Fortran is implicit or derived from context.
        // Let's assume the variables Hx, Hy, Hz, Gm2, Idye, Idze, Idxe, sggMiHx, sggMiHy, sggMiHz
        // are accessible in the current scope (e.g., global or class members).

        double GM2_1;
        if (Nodal_Hx.numHard == 0) {
            std::cout << "Devia H nodal/field sources untested. Aborting" << std::endl;
            // In C++, we might throw an exception or exit. 
            // Fortran 'stop' terminates the program.
            std::exit(1); 
        }
        GM2_1 = Gm2[1]; // Assuming 1-based indexing in Fortran, so index 1 is the second element? 
                         // Or is it 1-based array? Fortran arrays are often 1-based.
                         // If Gm2 is std::vector, and Fortran uses Gm2(1), it usually means index 1.
                         // If the vector is 0-based, we might need Gm2[0]. 
                         // However, rule 4 says "Preserve 1-based indexing where Fortran uses it".
                         // This usually implies the logic uses 1-based indices. 
                         // If the underlying C++ container is 0-based, we adjust access.
                         // But if we are preserving names, we might have a wrapper or just use 0-based 
                         // if the Fortran code was actually 0-based in implementation despite syntax.
                         // Let's assume standard Fortran 1-based indexing for Gm2(1) -> Gm2[0] if vector is 0-based.
                         // BUT, if the vector is sized to N+1 and we use index 1, we keep index 1.
                         // Let's assume the C++ vector is 0-based and Fortran index 1 maps to C++ index 0.
                         // Wait, "Preserve 1-based indexing" usually means the *logic* remains 1-based, 
                         // so if Fortran does Gm2(1), and we use std::vector, we do Gm2[0].
                         // However, if the array is allocated with size N+1, we might keep index 1.
                         // Let's stick to the most common translation: Fortran 1-based -> C++ 0-based.
                         // So Gm2(1) -> Gm2[0].
                         // Let's re-read rule 4: "Preserve 1-based indexing where Fortran uses it."
                         // This is ambiguous. It could mean "keep the indices 1, 2, 3..." in the code.
                         // If I change Gm2(1) to Gm2[0], I am changing the index.
                         // If I change Gm2(1) to Gm2[1], I am preserving the index, but the vector must be 1-based.
                         // To preserve names and logic, I will assume the C++ arrays/vectors are sized to accommodate 1-based indexing 
                         // (i.e., size N+1) and accessed with 1-based indices.
        
        // timei calculation
        double timei = sgg.tiempo[timeinstant] + 0.5 * sgg.dt; 
        // Note: _RKIND is a kind parameter. In C++, we use double.

        // barridonodalhardHx
        for (int ii = 1; ii <= Nodal_Hx.numHard; ++ii) {
            if (Nodal_Hx.nodHard[ii].IsInitialValue && (timeinstant != 0)) {
                continue;
            }
            double amp = Nodal_Hx.nodHard[ii].punto.amplitude;
            for (int k = Nodal_Hx.nodHard[ii].punto.zi; k <= Nodal_Hx.nodHard[ii].punto.ze; ++k) {
                int k_m = k - b.Hx.ZI;
                for (int j = Nodal_Hx.nodHard[ii].punto.yi; j <= Nodal_Hx.nodHard[ii].punto.ye; ++j) {
                    int j_m = j - b.Hx.YI;
                    for (int i = Nodal_Hx.nodHard[ii].punto.xi; i <= Nodal_Hx.nodHard[ii].punto.xe; ++i) {
                        int i_m = i - b.Hx.XI;
                        int medio = sggMiHx[i_m]; // Assuming sggMiHx is 3D flattened or accessed differently?
                                                  // Fortran: sggMiHx(i_m,j_m,k_m). 
                                                  // If it's a 3D array in Fortran, it's likely contiguous.
                                                  // In C++, if it's std::vector, we need to map 3D to 1D.
                                                  // Or if it's a 3D vector<vector<vector<int>>>.
                                                  // Let's assume a helper or direct access if structure is known.
                                                  // For translation, we'll keep the syntax similar if possible, 
                                                  // but C++ doesn't have 3D arrays with 3 indices natively in std::vector.
                                                  // We'll assume sggMiHx is accessed via a function or flattened index.
                                                  // Let's assume a flattened access for now: sggMiHx(i_m, j_m, k_m) -> sggMiHx[i_m * dimY * dimZ + j_m * dimZ + k_m]
                                                  // But we don't have dims. Let's assume it's a 3D vector or we have a macro/function.
                                                  // To preserve names, let's assume sggMiHx is a 3D array-like object.
                        medio = sggMiHx(i_m, j_m, k_m); // Assuming operator() or similar
                        if (!sgg.Med(medio).Is.PMC) {
                            Hx[i_m][j_m][k_m] = amp * evolucion(timei, Nodal_Hx.nodHard[ii]);
                        }
                    }
                }
            }
        }

        // barridonodalsoftHx
        for (int ii = 1; ii <= Nodal_Hx.numSoft; ++ii) {
            if (Nodal_Hx.nodSoft[ii].IsInitialValue && (timeinstant != 0)) {
                continue;
            }
            double amp = Nodal_Hx.nodSoft[ii].punto.amplitude;
            for (int k = Nodal_Hx.nodSoft[ii].punto.zi; k <= Nodal_Hx.nodSoft[ii].punto.ze; ++k) {
                int k_m = k - b.Hx.ZI;
                for (int j = Nodal_Hx.nodSoft[ii].punto.yi; j <= Nodal_Hx.nodSoft[ii].punto.ye; ++j) {
                    int j_m = j - b.Hx.YI;
                    for (int i = Nodal_Hx.nodSoft[ii].punto.xi; i <= Nodal_Hx.nodSoft[ii].punto.xe; ++i) {
                        int i_m = i - b.Hx.XI;
                        int medio = sggMiHx(i_m, j_m, k_m);
                        if (!sgg.Med(medio).Is.PMC) {
                            Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] - Gm2(medio) * Idye[j_m] * Idze[k_m] * amp * evolucion(timei, Nodal_Hx.nodSoft[ii]);
                        }
                    }
                }
            }
        }

        // barridonodalhardHy
        for (int ii = 1; ii <= Nodal_Hy.numHard; ++ii) {
            if (Nodal_Hy.nodHard[ii].IsInitialValue && (timeinstant != 0)) {
                continue;
            }
            double amp = Nodal_Hy.nodHard[ii].punto.amplitude;
            for (int k = Nodal_Hy.nodHard[ii].punto.zi; k <= Nodal_Hy.nodHard[ii].punto.ze; ++k) {
                int k_m = k - b.Hy.ZI;
                for (int j = Nodal_Hy.nodHard[ii].punto.yi; j <= Nodal_Hy.nodHard[ii].punto.ye; ++j) {
                    int j_m = j - b.Hy.YI;
                    for (int i = Nodal_Hy.nodHard[ii].punto.xi; i <= Nodal_Hy.nodHard[ii].punto.xe; ++i) {
                        int i_m = i - b.Hy.XI;
                        int medio = sggMiHx(i_m, j_m, k_m);
                        if (!sgg.Med(medio).Is.PMC) {
                            Hy[i_m][j_m][k_m] = amp * evolucion(timei, Nodal_Hy.nodHard[ii]);
                        }
                    }
                }
            }
        }

        // barridonodalsoftHy
        for (int ii = 1; ii <= Nodal_Hy.numSoft; ++ii) {
            if (Nodal_Hy.nodSoft[ii].IsInitialValue && (timeinstant != 0)) {
                continue;
            }
            double amp = Nodal_Hy.nodSoft[ii].punto.amplitude;
            for (int k = Nodal_Hy.nodSoft[ii].punto.zi; k <= Nodal_Hy.nodSoft[ii].punto.ze; ++k) {
                int k_m = k - b.Hy.ZI;
                for (int j = Nodal_Hy.nodSoft[ii].punto.yi; j <= Nodal_Hy.nodSoft[ii].punto.ye; ++j) {
                    int j_m = j - b.Hy.YI;
                    for (int i = Nodal_Hy.nodSoft[ii].punto.xi; i <= Nodal_Hy.nodSoft[ii].punto.xe; ++i) {
                        int i_m = i - b.Hy.XI;
                        int medio = sggMiHy(i_m, j_m, k_m);
                        if (!sgg.Med(medio).Is.PMC) {
                            Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] - Gm2(medio) * Idxe[i_m] * Idze[k_m] * amp * evolucion(timei, Nodal_Hy.nodSoft[ii]);
                        }
                    }
                }
            }
        }

        // barridonodalhardHz
        for (int ii = 1; ii <= Nodal_Hz.numHard; ++ii) {
            if (Nodal_Hz.nodHard[ii].IsInitialValue && (timeinstant != 0)) {
                continue;
            }
            double amp = Nodal_Hz.nodHard[ii].punto.amplitude;
            for (int k = Nodal_Hz.nodHard[ii].punto.zi; k <= Nodal_Hz.nodHard[ii].punto.ze; ++k) {
                int k_m = k - b.Hz.ZI;
                for (int j = Nodal_Hz.nodHard[ii].punto.yi; j <= Nodal_Hz.nodHard[ii].punto.ye; ++j) {
                    int j_m = j - b.Hz.YI;
                    for (int i = Nodal_Hz.nodHard[ii].punto.xi; i <= Nodal_Hz.nodHard[ii].punto.xe; ++i) {
                        int i_m = i - b.Hz.XI;
                        int medio = sggMiHx(i_m, j_m, k_m);
                        if (!sgg.Med(medio).Is.PMC) {
                            Hz[i_m][j_m][k_m] = amp * evolucion(timei, Nodal_Hz.nodHard[ii]);
                        }
                    }
                }
            }
        }

        // barridonodalsoftHz
        for (int ii = 1; ii <= Nodal_Hz.numSoft; ++ii) {
            if (Nodal_Hz.nodSoft[ii].IsInitialValue && (timeinstant != 0)) {
                continue;
            }
            double amp = Nodal_Hz.nodSoft[ii].punto.amplitude;
            for (int k = Nodal_Hz.nodSoft[ii].punto.zi; k <= Nodal_Hz.nodSoft[ii].punto.ze; ++k) {
                int k_m = k - b.Hz.ZI;
                for (int j = Nodal_Hz.nodSoft[ii].punto.yi; j <= Nodal_Hz.nodSoft[ii].punto.ye; ++j) {
                    int j_m = j - b.Hz.YI;
                    for (int i = Nodal_Hz.nodSoft[ii].punto.xi; i <= Nodal_Hz.nodSoft[ii].punto.xe; ++i) {
                        int i_m = i - b.Hz.XI;
                        int medio = sggMiHz(i_m, j_m, k_m);
                        if (!sgg.Med(medio).Is.PMC) {
                            Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] - Gm2(medio) * Idye[j_m] * Idxe[i_m] * amp * evolucion(timei, Nodal_Hz.nodSoft[ii]);
                        }
                    }
                }
            }
        }

        return;
    }

    void getnodal(nodsou_t*& rNodal_Ex, nodsou_t*& rNodal_Ey, nodsou_t*& rNodal_Ez, 
                  nodsou_t*& rNodal_Hx, nodsou_t*& rNodal_Hy, nodsou_t*& rNodal_Hz) {
        
        rNodal_Ex = &Nodal_Ex;
        rNodal_Ey = &Nodal_Ey;
        rNodal_Ez = &Nodal_Ez;
        rNodal_Hx = &Nodal_Hx;
        rNodal_Hy = &Nodal_Hy;
        rNodal_Hz = &Nodal_Hz;

        return;
    }

} // namespace nodalsources_m