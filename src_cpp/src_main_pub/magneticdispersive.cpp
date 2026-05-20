#include <vector>
#include <complex>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

// Assuming these types are defined in other headers based on the Fortran 'use' statements
// We need to forward declare or include them. Since I don't have the source for FDETYPES_m and Report_m,
// I will assume standard definitions or provide stubs if necessary. 
// However, to strictly follow "Preserve ALL names", I will assume the existence of these types.

// Stub definitions for external types to make the code compile conceptually.
// In a real scenario, these would be included from their respective headers.

struct SGGFDTDINFO_t {
    int NumMedia;
    // Assuming these structures exist based on usage
    struct {
        int XI, XE, YI, YE, ZI, ZE;
    } Alloc[10]; // Placeholder size, assuming indices like iHx are constants
    struct {
        int XI, XE, YI, YE, ZI, ZE;
    } Sweep[10];
    
    struct {
        bool IsMdispersive;
        bool IsMdispersiveANIS;
        bool IsPML;
        // Placeholder for Mdispersive data
        struct {
            int numpolres11;
            std::vector<double> c11;
            std::vector<double> a11;
            double Sigmam11;
            double mu11;
        } Mdispersive[10]; // Placeholder size
    } Med[10]; // Placeholder size

    double dt;
};

// Constants for indices
const int iHx = 0;
const int iHy = 1;
const int iHz = 2;

// Types from FDETYPES_m
typedef double RKIND;
typedef std::complex<double> CKIND;

struct media_matrices_t {
    std::vector<std::vector<std::vector<int>>> sggMiHx;
    std::vector<std::vector<std::vector<int>>> sggMiHy;
    std::vector<std::vector<std::vector<int>>> sggMiHz;
};

// Types from Report_m
void print11(int unit, const std::string& msg) {
    std::cout << msg << std::endl;
}

const std::string SEPARADOR = "========================================";

namespace Mdispersives_m {

    struct field_t {
        int i, j, k;
        int WhatField;
        double* FieldPresent; // Pointer to background field
        double FieldPrevious;
        std::vector<std::complex<double>> Current;
    };

    struct Mdispersive_t {
        int indexmed;
        int numnodesHx, numnodesHy, numnodesHz;
        std::vector<std::complex<double>> Beta;
        std::vector<std::complex<double>> Kappa;
        std::vector<std::complex<double>> GM3;
        std::vector<field_t> NodesHx;
        std::vector<field_t> NodesHy;
        std::vector<field_t> NodesHz;
    };

    struct Mdispersive2_t {
        int NumMdispersives;
        std::vector<Mdispersive_t> Medium;
    };

    // LOCAL VARIABLES
    Mdispersive2_t MDutton;

    void InitMdispersives(const SGGFDTDINFO_t& sgg, const media_matrices_t& media, 
                          std::vector<double>& GM1, std::vector<double>& GM2,
                          std::vector<std::vector<std::vector<double>>>& Hx,
                          std::vector<std::vector<std::vector<double>>>& Hy,
                          std::vector<std::vector<std::vector<double>>>& Hz,
                          bool& ThereAreMdispersives, bool resume) {
        
        MDutton.Medium.clear();
        MDutton.NumMdispersives = 0;

        ThereAreMdispersives = false;
        int conta = 0;
        for (int jmed = 1; jmed <= sgg.NumMedia; ++jmed) {
            if (sgg.Med[jmed].IsMdispersive && !sgg.Med[jmed].IsMdispersiveANIS) {
                conta++;
            }
        }

        MDutton.NumMdispersives = conta;
        MDutton.Medium.resize(MDutton.NumMdispersives);

        conta = 0;
        for (int jmed = 1; jmed <= sgg.NumMedia; ++jmed) {
            if (sgg.Med[jmed].IsMdispersive && !sgg.Med[jmed].IsMdispersiveANIS) {
                conta++;
                int idx = conta - 1; // 0-based index for vector
                MDutton.Medium[idx].indexmed = jmed;
                int numpolres = sgg.Med[jmed].Mdispersive[1].numpolres11;
                MDutton.Medium[idx].numpolres11 = numpolres;

                MDutton.Medium[idx].Beta.resize(numpolres, 0.0);
                MDutton.Medium[idx].Kappa.resize(numpolres, 0.0);
                MDutton.Medium[idx].GM3.resize(numpolres, 0.0);

                for (int i1 = 1; i1 <= numpolres; ++i1) {
                    int i1_idx = i1 - 1;
                    double a11_val = sgg.Med[jmed].Mdispersive[1].a11[i1_idx];
                    double c11_val = sgg.Med[jmed].Mdispersive[1].c11[i1_idx];
                    
                    MDutton.Medium[idx].Kappa[i1_idx] = (1.0 + a11_val * sgg.dt / 2.0) / 
                                                        (1.0 - a11_val * sgg.dt / 2.0);
                    MDutton.Medium[idx].Beta[i1_idx] = (c11_val * sgg.dt) / 
                                                       (1.0 - a11_val * sgg.dt / 2.0);
                }
            }
        }

        // Calculate coefficients
        for (int jmed = 1; jmed <= MDutton.NumMdispersives; ++jmed) {
            int tempindex = MDutton.Medium[jmed-1].indexmed;
            int numpolres = sgg.Med[tempindex].Mdispersive[1].numpolres11;
            double tempo = 0.0;
            for (int i1 = 1; i1 <= numpolres; ++i1) {
                tempo += std::real(MDutton.Medium[jmed-1].Beta[i1-1]);
            }
            
            double mu11 = sgg.Med[tempindex].Mdispersive[1].mu11;
            double Sigmam11 = sgg.Med[tempindex].Mdispersive[1].Sigmam11;
            
            GM1[tempindex] = (2.0 * mu11 + tempo - Sigmam11 * sgg.dt) / 
                             (2.0 * mu11 + tempo + Sigmam11 * sgg.dt);
            GM2[tempindex] = (2.0 * sgg.dt) / 
                             (2.0 * mu11 + tempo + Sigmam11 * sgg.dt);
            
            for (int i1 = 1; i1 <= numpolres; ++i1) {
                MDutton.Medium[jmed-1].GM3[i1-1] = GM2[tempindex] / 2.0 * (1.0 + MDutton.Medium[jmed-1].Kappa[i1-1]);
            }
        }

        for (int jmed = 1; jmed <= MDutton.NumMdispersives; ++jmed) {
            int tempindex = MDutton.Medium[jmed-1].indexmed;
            
            // Hx
            int conta_hx = 0;
            for (int k1 = sgg.Sweep[iHx].ZI; k1 <= sgg.Sweep[iHx].ZE; ++k1) {
                for (int j1 = sgg.Sweep[iHx].YI; j1 <= sgg.Sweep[iHx].YE; ++j1) {
                    for (int i1 = sgg.Sweep[iHx].XI; i1 <= sgg.Sweep[iHx].XE; ++i1) {
                        if (media.sggMiHx[i1][j1][k1] == tempindex) {
                            conta_hx++;
                        }
                    }
                }
            }
            ThereAreMdispersives = ThereAreMdispersives || (conta_hx != 0);
            MDutton.Medium[jmed-1].numnodesHx = conta_hx;
            MDutton.Medium[jmed-1].NodesHx.resize(conta_hx);
            
            for (int i1 = 0; i1 < conta_hx; ++i1) {
                MDutton.Medium[jmed-1].NodesHx[i1].Current.resize(sgg.Med[tempindex].Mdispersive[1].numpolres11);
            }
            
            conta_hx = 0;
            for (int k1 = sgg.Sweep[iHx].ZI; k1 <= sgg.Sweep[iHx].ZE; ++k1) {
                for (int j1 = sgg.Sweep[iHx].YI; j1 <= sgg.Sweep[iHx].YE; ++j1) {
                    for (int i1 = sgg.Sweep[iHx].XI; i1 <= sgg.Sweep[iHx].XE; ++i1) {
                        if (media.sggMiHx[i1][j1][k1] == tempindex) {
                            conta_hx++;
                            int idx = conta_hx - 1;
                            MDutton.Medium[jmed-1].NodesHx[idx].i = i1;
                            MDutton.Medium[jmed-1].NodesHx[idx].j = j1;
                            MDutton.Medium[jmed-1].NodesHx[idx].k = k1;
                            MDutton.Medium[jmed-1].NodesHx[idx].WhatField = iHx;
                            MDutton.Medium[jmed-1].NodesHx[idx].FieldPresent = &Hx[i1][j1][k1];
                        }
                    }
                }
            }

            // Hy
            int conta_hy = 0;
            for (int k1 = sgg.Sweep[iHy].ZI; k1 <= sgg.Sweep[iHy].ZE; ++k1) {
                for (int j1 = sgg.Sweep[iHy].YI; j1 <= sgg.Sweep[iHy].YE; ++j1) {
                    for (int i1 = sgg.Sweep[iHy].XI; i1 <= sgg.Sweep[iHy].XE; ++i1) {
                        if (media.sggMiHy[i1][j1][k1] == tempindex) {
                            conta_hy++;
                        }
                    }
                }
            }
            ThereAreMdispersives = ThereAreMdispersives || (conta_hy != 0);
            MDutton.Medium[jmed-1].numnodesHy = conta_hy;
            MDutton.Medium[jmed-1].NodesHy.resize(conta_hy);
            
            for (int i1 = 0; i1 < conta_hy; ++i1) {
                MDutton.Medium[jmed-1].NodesHy[i1].Current.resize(sgg.Med[tempindex].Mdispersive[1].numpolres11);
            }
            
            conta_hy = 0;
            for (int k1 = sgg.Sweep[iHy].ZI; k1 <= sgg.Sweep[iHy].ZE; ++k1) {
                for (int j1 = sgg.Sweep[iHy].YI; j1 <= sgg.Sweep[iHy].YE; ++j1) {
                    for (int i1 = sgg.Sweep[iHy].XI; i1 <= sgg.Sweep[iHy].XE; ++i1) {
                        if (media.sggMiHy[i1][j1][k1] == tempindex) {
                            conta_hy++;
                            int idx = conta_hy - 1;
                            MDutton.Medium[jmed-1].NodesHy[idx].i = i1;
                            MDutton.Medium[jmed-1].NodesHy[idx].j = j1;
                            MDutton.Medium[jmed-1].NodesHy[idx].k = k1;
                            MDutton.Medium[jmed-1].NodesHy[idx].WhatField = iHy;
                            MDutton.Medium[jmed-1].NodesHy[idx].FieldPresent = &Hy[i1][j1][k1];
                        }
                    }
                }
            }

            // Hz
            int conta_hz = 0;
            for (int k1 = sgg.Sweep[iHz].ZI; k1 <= sgg.Sweep[iHz].ZE; ++k1) {
                for (int j1 = sgg.Sweep[iHz].YI; j1 <= sgg.Sweep[iHz].YE; ++j1) {
                    for (int i1 = sgg.Sweep[iHz].XI; i1 <= sgg.Sweep[iHz].XE; ++i1) {
                        if (media.sggMiHz[i1][j1][k1] == tempindex) {
                            conta_hz++;
                        }
                    }
                }
            }
            ThereAreMdispersives = ThereAreMdispersives || (conta_hz != 0);
            MDutton.Medium[jmed-1].numnodesHz = conta_hz;
            MDutton.Medium[jmed-1].NodesHz.resize(conta_hz);
            
            for (int i1 = 0; i1 < conta_hz; ++i1) {
                MDutton.Medium[jmed-1].NodesHz[i1].Current.resize(sgg.Med[tempindex].Mdispersive[1].numpolres11);
            }
            
            conta_hz = 0;
            for (int k1 = sgg.Sweep[iHz].ZI; k1 <= sgg.Sweep[iHz].ZE; ++k1) {
                for (int j1 = sgg.Sweep[iHz].YI; j1 <= sgg.Sweep[iHz].YE; ++j1) {
                    for (int i1 = sgg.Sweep[iHz].XI; i1 <= sgg.Sweep[iHz].XE; ++i1) {
                        if (media.sggMiHz[i1][j1][k1] == tempindex) {
                            conta_hz++;
                            int idx = conta_hz - 1;
                            MDutton.Medium[jmed-1].NodesHz[idx].i = i1;
                            MDutton.Medium[jmed-1].NodesHz[idx].j = j1;
                            MDutton.Medium[jmed-1].NodesHz[idx].k = k1;
                            MDutton.Medium[jmed-1].NodesHz[idx].WhatField = iHz;
                            MDutton.Medium[jmed-1].NodesHz[idx].FieldPresent = &Hz[i1][j1][k1];
                        }
                    }
                }
            }
        }

        // Resume or start
        for (int jmed = 1; jmed <= MDutton.NumMdispersives; ++jmed) {
            int numpolres = sgg.Med[MDutton.Medium[jmed-1].indexmed].Mdispersive[1].numpolres11;
            if (!resume) {
                // Hx, Jx
                for (int i1 = 0; i1 < MDutton.Medium[jmed-1].numnodesHx; ++i1) {
                    MDutton.Medium[jmed-1].NodesHx[i1].FieldPrevious = 0.0;
                    for (int k1 = 0; k1 < numpolres; ++k1) {
                        MDutton.Medium[jmed-1].NodesHx[i1].Current[k1] = 0.0;
                    }
                }
                // Hy, Jy
                for (int i1 = 0; i1 < MDutton.Medium[jmed-1].numnodesHy; ++i1) {
                    MDutton.Medium[jmed-1].NodesHy[i1].FieldPrevious = 0.0;
                    for (int k1 = 0; k1 < numpolres; ++k1) {
                        MDutton.Medium[jmed-1].NodesHy[i1].Current[k1] = 0.0;
                    }
                }
                // Hz, Jz
                for (int i1 = 0; i1 < MDutton.Medium[jmed-1].numnodesHz; ++i1) {
                    MDutton.Medium[jmed-1].NodesHz[i1].FieldPrevious = 0.0;
                    for (int k1 = 0; k1 < numpolres; ++k1) {
                        MDutton.Medium[jmed-1].NodesHz[i1].Current[k1] = 0.0;
                    }
                }
            } else {
                // Hx, Jx
                for (int i1 = 0; i1 < MDutton.Medium[jmed-1].numnodesHx; ++i1) {
                    std::cin >> MDutton.Medium[jmed-1].NodesHx[i1].FieldPrevious; // Assuming unit 14 is stdin or a specific file stream
                    for (int k1 = 0; k1 < numpolres; ++k1) {
                        std::cin >> MDutton.Medium[jmed-1].NodesHx[i1].Current[k1];
                    }
                }
                // Hy, Jy
                for (int i1 = 0; i1 < MDutton.Medium[jmed-1].numnodesHy; ++i1) {
                    std::cin >> MDutton.Medium[jmed-1].NodesHy[i1].FieldPrevious;
                    for (int k1 = 0; k1 < numpolres; ++k1) {
                        std::cin >> MDutton.Medium[jmed-1].NodesHy[i1].Current[k1];
                    }
                }
                // Hz, Jz
                for (int i1 = 0; i1 < MDutton.Medium[jmed-1].numnodesHz; ++i1) {
                    std::cin >> MDutton.Medium[jmed-1].NodesHz[i1].FieldPrevious;
                    for (int k1 = 0; k1 < numpolres; ++k1) {
                        std::cin >> MDutton.Medium[jmed-1].NodesHz[i1].Current[k1];
                    }
                }
            }
        }
    }

    void AdvanceMdispersiveH(const SGGFDTDINFO_t& sgg) {
        for (int jmed = 1; jmed <= MDutton.NumMdispersives; ++jmed) {
            int numpolres = MDutton.Medium[jmed-1].numpolres11;
            
            // Hx, Jx
            for (int i1 = 0; i1 < MDutton.Medium[jmed-1].numnodesHx; ++i1) {
                field_t& tempnode = MDutton.Medium[jmed-1].NodesHx[i1];
                for (int k1 = 0; k1 < numpolres; ++k1) {
                    tempnode.FieldPresent = tempnode.FieldPresent - std::real(MDutton.Medium[jmed-1].GM3[k1] * tempnode.Current[k1]);
                }
                for (int k1 = 0; k1 < numpolres; ++k1) {
                    tempnode.Current[k1] = MDutton.Medium[jmed-1].Kappa[k1] * tempnode.Current[k1] + 
                                           MDutton.Medium[jmed-1].Beta[k1] / sgg.dt * (tempnode.FieldPresent - tempnode.FieldPrevious);
                }
                tempnode.FieldPrevious = tempnode.FieldPresent;
            }
            
            // Hy, Jy
            for (int i1 = 0; i1 < MDutton.Medium[jmed-1].numnodesHy; ++i1) {
                field_t& tempnode = MDutton.Medium[jmed-1].NodesHy[i1];
                for (int k1 = 0; k1 < numpolres; ++k1) {
                    tempnode.FieldPresent = tempnode.FieldPresent - std::real(MDutton.Medium[jmed-1].GM3[k1] * tempnode.Current[k1]);
                }
                for (int k1 = 0; k1 < numpolres; ++k1) {
                    tempnode.Current[k1] = MDutton.Medium[jmed-1].Kappa[k1] * tempnode.Current[k1] + 
                                           MDutton.Medium[jmed-1].Beta[k1] / sgg.dt * (tempnode.FieldPresent - tempnode.FieldPrevious);
                }
                tempnode.FieldPrevious = tempnode.FieldPresent;
            }

            // Hz, Jz
            for (int i1 = 0; i1 < MDutton.Medium[jmed-1].numnodesHz; ++i1) {
                field_t& tempnode = MDutton.Medium[jmed-1].NodesHz[i1];
                for (int k1 = 0; k1 < numpolres; ++k1) {
                    tempnode.FieldPresent = tempnode.FieldPresent - std::real(MDutton.Medium[jmed-1].GM3[k1] * tempnode.Current[k1]);
                }
                for (int k1 = 0; k1 < numpolres; ++k1) {
                    tempnode.Current[k1] = MDutton.Medium[jmed-1].Kappa[k1] * tempnode.Current[k1] + 
                                           MDutton.Medium[jmed-1].Beta[k1] / sgg.dt * (tempnode.FieldPresent - tempnode.FieldPrevious);
                }
                tempnode.FieldPrevious = tempnode.FieldPresent;
            }
        }
    }

    void StoreFieldsMdispersives() {
        for (int jmed = 1; jmed <= MDutton.NumMdispersives; ++jmed) {
            int numpolres = MDutton.Medium[jmed-1].numpolres11;
            
            // Hx, Jx
            for (int i1 = 0; i1 < MDutton.Medium[jmed-1].numnodesHx; ++i1) {
                std::cout << MDutton.Medium[jmed-1].NodesHx[i1].FieldPrevious << std::endl;
                for (int k1 = 0; k1 < numpolres; ++k1) {
                    std::cout << MDutton.Medium[jmed-1].NodesHx[i1].Current[k1] << std::endl;
                }
            }
            // Hy, Jy
            for (int i1 = 0; i1 < MDutton.Medium[jmed-1].numnodesHy; ++i1) {
                std::cout << MDutton.Medium[jmed-1].NodesHy[i1].FieldPrevious << std::endl;
                for (int k1 = 0; k1 < numpolres; ++k1) {
                    std::cout << MDutton.Medium[jmed-1].NodesHy[i1].Current[k1] << std::endl;
                }
            }
            // Hz, Jz
            for (int i1 = 0; i1 < MDutton.Medium[jmed-1].numnodesHz; ++i1) {
                std::cout << MDutton.Medium[jmed-1].NodesHz[i1].FieldPrevious << std::endl;
                for (int k1 = 0; k1 < numpolres; ++k1) {
                    std::cout << MDutton.Medium[jmed-1].NodesHz[i1].Current[k1] << std::endl;
                }
            }
        }
    }

    void DestroyMdispersives(SGGFDTDINFO_t& sgg) {
        // Free up memory for sgg.Med
        for (int i = 1; i <= sgg.NumMedia; ++i) {
            if (sgg.Med[i].IsMdispersive && !sgg.Med[i].IsPML && !sgg.Med[i].IsMdispersiveANIS) {
                sgg.Med[i].Mdispersive[1].c11.clear();
                sgg.Med[i].Mdispersive[1].a11.clear();
            }
        }
        for (int i = 1; i <= sgg.NumMedia; ++i) {
            if (sgg.Med[i].IsMdispersive && !sgg.Med[i].IsPML && !sgg.Med[i].IsMdispersiveANIS) {
                // In C++, if Mdispersive is a vector or array of structs, we just clear it.
                // Assuming it was dynamically allocated in Fortran, we simulate deallocation by clearing.
                sgg.Med[i].Mdispersive[1].c11.clear();
                sgg.Med[i].Mdispersive[1].a11.clear();
            }
        }

        for (int jmed = 1; jmed <= MDutton.NumMdispersives; ++jmed) {
            MDutton.Medium[jmed-1].Beta.clear();
            MDutton.Medium[jmed-1].Kappa.clear();
            MDutton.Medium[jmed-1].GM3.clear();

            for (int i1 = 0; i1 < MDutton.Medium[jmed-1].numnodesHx; ++i1) {
                MDutton.Medium[jmed-1].NodesHx[i1].Current.clear();
            }
            MDutton.Medium[jmed-1].NodesHx.clear();
            
            for (int i1 = 0; i1 < MDutton.Medium[jmed-1].numnodesHy; ++i1) {
                MDutton.Medium[jmed-1].NodesHy[i1].Current.clear();
            }
            MDutton.Medium[jmed-1].NodesHy.clear();
            
            for (int i1 = 0; i1 < MDutton.Medium[jmed-1].numnodesHz; ++i1) {
                MDutton.Medium[jmed-1].NodesHz[i1].Current.clear();
            }
            MDutton.Medium[jmed-1].NodesHz.clear();
        }
        MDutton.Medium.clear();
    }

}