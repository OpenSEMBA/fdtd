#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <string>

// Forward declarations for types defined in other modules
// These would typically be in separate headers like FDETYPES_m.h and Report_m.h

// Assuming FDETYPES_m provides:
// RKIND, CKIND, SGGFDTDINFO_t, media_matrices_t, iHx, iHy, iHz, SEPARADOR, print11
// We need to mock these or assume they are available.
// For the purpose of this translation, we assume the existence of these types/functions
// as they are used but not defined in this snippet.

// Mocking necessary external types/functions to make the code compile conceptually
// In a real scenario, these would be included from their respective headers.

// extern double RKIND; // Placeholder for kind parameter logic, usually handled by typedefs
// extern std::complex<double> CKIND; // Placeholder

// Assuming these types are defined elsewhere:
// struct SGGFDTDINFO_t;
// struct media_matrices_t;
// enum { iHx, iHy, iHz };
// std::string SEPARADOR;
// void print11(int, const std::string&);

// To make this a standalone compilable unit for demonstration, we define minimal stubs
// if they are not provided. However, per instructions, we preserve names.
// We will assume the following includes/types exist in the environment:
// #include "FDETYPES_m.h"
// #include "Report_m.h"

// If those headers are not available, the following stubs are necessary for compilation:
/*
namespace FDETYPES_m {
    using RKIND = double;
    using CKIND = std::complex<double>;
    
    struct Alloc_t {
        int XI, XE, YI, YE, ZI, ZE;
    };
    
    struct Sweep_t {
        int XI, XE, YI, YE, ZI, ZE;
    };
    
    struct Is_t {
        bool Mdispersive;
        bool MdispersiveANIS;
        bool PML;
    };
    
    struct Mdispersive_t {
        int numpolres11;
        std::vector<double> a11;
        std::vector<double> c11;
        double Sigmam11;
        double mu11;
    };
    
    struct Media_t {
        Is_t Is;
        std::vector<Mdispersive_t> Mdispersive;
    };
    
    struct SGGFDTDINFO_t {
        int NumMedia;
        std::vector<Media_t> Med;
        std::vector<Alloc_t> Alloc;
        std::vector<Sweep_t> Sweep;
        double dt;
    };
}

namespace Report_m {
    extern std::string SEPARADOR;
    void print11(int, const std::string&);
}

// Mock implementation for print11
void Report_m::print11(int, const std::string& s) {
    std::cout << s << std::endl;
}
*/

// However, to strictly follow "Translate this Fortran code", we assume the headers are present.
// We will write the C++ code assuming the necessary headers are included.

// Note: The Fortran code uses 1-based indexing for arrays in many places.
// C++ uses 0-based indexing. We will adjust vector sizes and loops accordingly,
// or use 1-based indexing by allocating size N+1 and ignoring index 0.
// Given the complexity and the instruction to preserve names and logic, 
// we will use std::vector and adjust indices to 0-based where appropriate for C++ idioms,
// BUT since the Fortran code explicitly uses 1-based loops (1:N), we will keep the logic
// consistent. If the Fortran code accesses index 1 to N, the C++ vector should be sized N+1
// and we use indices 1 to N, or we resize to N and use 0 to N-1.
// Looking at `allocate (MDutton%Medium(1 : MDutton%NumMdispersives))`, it uses 1-based.
// We will use 1-based indexing for vectors by allocating size N+1 and using indices 1..N.

namespace Mdispersives_m {

    // Assuming these are defined in FDETYPES_m
    // using RKIND = double;
    // using CKIND = std::complex<double>;

    struct field_t {
        int32_t i, j, k;
        int32_t WhatField;
        double* FieldPresent; // Pointer to the background field
        double FieldPrevious;
        std::vector<std::complex<double>> Current;
    };

    struct Mdispersive_t {
        int32_t indexmed;
        int32_t numnodesHx, numnodesHy, numnodesHz;
        int32_t numpolres11;
        std::vector<std::complex<double>> Beta;
        std::vector<std::complex<double>> Kappa;
        std::vector<std::complex<double>> GM3;
        std::vector<field_t> NodesHx;
        std::vector<field_t> NodesHy;
        std::vector<field_t> NodesHz;
    };

    struct Mdispersive2_t {
        int32_t NumMdispersives;
        std::vector<Mdispersive_t> Medium;
    };

    // LOCAL VARIABLES
    Mdispersive2_t MDutton;

    void InitMdispersives(FDETYPES_m::SGGFDTDINFO_t& sgg, 
                          media_matrices_t& media, 
                          std::vector<double>& GM1, 
                          std::vector<double>& GM2, 
                          std::vector<std::vector<std::vector<double>>>& Hx, 
                          std::vector<std::vector<std::vector<double>>>& Hy, 
                          std::vector<std::vector<std::vector<double>>>& Hz,
                          bool& ThereAreMdispersives, 
                          bool resume) {
        
        MDutton.Medium.clear();
        MDutton.Medium.shrink_to_fit(); // Ensure empty state

        ThereAreMdispersives = false;
        int conta = 0;
        for (int jmed = 1; jmed <= sgg.NumMedia; ++jmed) {
            if (sgg.Med[jmed].Is.Mdispersive && !sgg.Med[jmed].Is.MdispersiveANIS) {
                conta++;
            }
        }

        MDutton.NumMdispersives = conta;
        // Resize to NumMdispersives + 1 to accommodate 1-based indexing
        MDutton.Medium.resize(MDutton.NumMdispersives + 1);
        
        conta = 0;
        for (int jmed = 1; jmed <= sgg.NumMedia; ++jmed) {
            if (sgg.Med[jmed].Is.Mdispersive && !sgg.Med[jmed].Is.MdispersiveANIS) {
                conta++;
                MDutton.Medium[conta].indexmed = jmed;
                MDutton.Medium[conta].numpolres11 = sgg.Med[jmed].Mdispersive[1].numpolres11;
                
                int n = MDutton.Medium[conta].numpolres11;
                MDutton.Medium[conta].Beta.resize(n + 1, 0.0);
                MDutton.Medium[conta].Kappa.resize(n + 1, 0.0);
                MDutton.Medium[conta].GM3.resize(n + 1, 0.0);
                
                for (int i1 = 1; i1 <= n; ++i1) {
                    MDutton.Medium[conta].Kappa[i1] = (1.0 + sgg.Med[jmed].Mdispersive[1].a11[i1] * sgg.dt / 2.0) /
                                                      (1.0 - sgg.Med[jmed].Mdispersive[1].a11[i1] * sgg.dt / 2.0);
                    MDutton.Medium[conta].Beta[i1] = (sgg.Med[jmed].Mdispersive[1].c11[i1] * sgg.dt) /
                                                     (1.0 - sgg.Med[jmed].Mdispersive[1].a11[i1] * sgg.dt / 2.0);
                }
            }
        }

        // Calculate coefficients
        for (int jmed = 1; jmed <= MDutton.NumMdispersives; ++jmed) {
            int tempindex = MDutton.Medium[jmed].indexmed;
            int numpolres = sgg.Med[tempindex].Mdispersive[1].numpolres11;
            double tempo = 0.0;
            for (int i1 = 1; i1 <= numpolres; ++i1) { // Note: Fortran used NumPolRes, likely a typo for numpolres or global. Assuming numpolres based on context.
                 // If NumPolRes is a global constant, it should be defined. Here we assume it refers to the number of poles for this medium.
                 // The Fortran code says `Do i1=1,NumPolRes`. In the previous block, `numpolres` was defined.
                 // It is highly likely `NumPolRes` is a typo for `numpolres` or a global. 
                 // Given `numpolres` is local, we use `numpolres`.
                 tempo += std::real(MDutton.Medium[jmed].Beta[i1]);
            }
            
            GM1[tempindex] = (2.0 * sgg.Med[tempindex].Mdispersive[1].mu11 + tempo - sgg.Med[tempindex].Mdispersive[1].Sigmam11 * sgg.dt) /
                             (2.0 * sgg.Med[tempindex].Mdispersive[1].mu11 + tempo + sgg.Med[tempindex].Mdispersive[1].Sigmam11 * sgg.dt);
            GM2[tempindex] = (2.0 * sgg.dt) /
                             (2.0 * sgg.Med[tempindex].Mdispersive[1].mu11 + tempo + sgg.Med[tempindex].Mdispersive[1].Sigmam11 * sgg.dt);
            
            for (int i1 = 1; i1 <= numpolres; ++i1) {
                MDutton.Medium[jmed].GM3[i1] = GM2[tempindex] / 2.0 * (1.0 + MDutton.Medium[jmed].Kappa[i1]);
            }
        }

        for (int jmed = 1; jmed <= MDutton.NumMdispersives; ++jmed) {
            int tempindex = MDutton.Medium[jmed].indexmed;
            
            // !!!Hx
            conta = 0;
            for (int k1 = sgg.Sweep[iHx].ZI; k1 <= sgg.Sweep[iHx].ZE; ++k1) {
                for (int j1 = sgg.Sweep[iHx].YI; j1 <= sgg.Sweep[iHx].YE; ++j1) {
                    for (int i1 = sgg.Sweep[iHx].XI; i1 <= sgg.Sweep[iHx].XE; ++i1) {
                        if (media.sggMiHx[i1][j1][k1] == tempindex) {
                            conta++;
                        }
                    }
                }
            }

            ThereAreMdispersives = ThereAreMdispersives || (conta != 0);
            MDutton.Medium[jmed].NumNodesHx = conta;
            MDutton.Medium[jmed].NodesHx.resize(conta + 1);
            for (int i1 = 1; i1 <= conta; ++i1) {
                int n = sgg.Med[tempindex].Mdispersive[1].numpolres11;
                MDutton.Medium[jmed].NodesHx[i1].Current.resize(n + 1);
            }
            
            conta = 0;
            for (int k1 = sgg.Sweep[iHx].ZI; k1 <= sgg.Sweep[iHx].ZE; ++k1) {
                for (int j1 = sgg.Sweep[iHx].YI; j1 <= sgg.Sweep[iHx].YE; ++j1) {
                    for (int i1 = sgg.Sweep[iHx].XI; i1 <= sgg.Sweep[iHx].XE; ++i1) {
                        if (media.sggMiHx[i1][j1][k1] == tempindex) {
                            conta++;
                            MDutton.Medium[jmed].NodesHx[conta].i = i1;
                            MDutton.Medium[jmed].NodesHx[conta].j = j1;
                            MDutton.Medium[jmed].NodesHx[conta].k = k1;
                            MDutton.Medium[jmed].NodesHx[conta].WhatField = iHx;
                            MDutton.Medium[jmed].NodesHx[conta].FieldPresent = &Hx[i1][j1][k1];
                        }
                    }
                }
            }

            // !!!Hy
            conta = 0;
            for (int k1 = sgg.Sweep[iHy].ZI; k1 <= sgg.Sweep[iHy].ZE; ++k1) {
                for (int j1 = sgg.Sweep[iHy].YI; j1 <= sgg.Sweep[iHy].YE; ++j1) {
                    for (int i1 = sgg.Sweep[iHy].XI; i1 <= sgg.Sweep[iHy].XE; ++i1) {
                        if (media.sggMiHy[i1][j1][k1] == tempindex) {
                            conta++;
                        }
                    }
                }
            }

            ThereAreMdispersives = ThereAreMdispersives || (conta != 0);
            MDutton.Medium[jmed].NumNodesHy = conta;
            MDutton.Medium[jmed].NodesHy.resize(conta + 1);
            for (int i1 = 1; i1 <= conta; ++i1) {
                int n = sgg.Med[tempindex].Mdispersive[1].numpolres11;
                MDutton.Medium[jmed].NodesHy[i1].Current.resize(n + 1);
            }
            
            conta = 0;
            for (int k1 = sgg.Sweep[iHy].ZI; k1 <= sgg.Sweep[iHy].ZE; ++k1) {
                for (int j1 = sgg.Sweep[iHy].YI; j1 <= sgg.Sweep[iHy].YE; ++j1) {
                    for (int i1 = sgg.Sweep[iHy].XI; i1 <= sgg.Sweep[iHy].XE; ++i1) {
                        if (media.sggMiHy[i1][j1][k1] == tempindex) {
                            conta++;
                            MDutton.Medium[jmed].NodesHy[conta].i = i1;
                            MDutton.Medium[jmed].NodesHy[conta].j = j1;
                            MDutton.Medium[jmed].NodesHy[conta].k = k1;
                            MDutton.Medium[jmed].NodesHy[conta].WhatField = iHy;
                            MDutton.Medium[jmed].NodesHy[conta].FieldPresent = &Hy[i1][j1][k1];
                        }
                    }
                }
            }

            // !!!Hz
            conta = 0;
            for (int k1 = sgg.Sweep[iHz].ZI; k1 <= sgg.Sweep[iHz].ZE; ++k1) {
                for (int j1 = sgg.Sweep[iHz].YI; j1 <= sgg.Sweep[iHz].YE; ++j1) {
                    for (int i1 = sgg.Sweep[iHz].XI; i1 <= sgg.Sweep[iHz].XE; ++i1) {
                        if (media.sggMiHz[i1][j1][k1] == tempindex) {
                            conta++;
                        }
                    }
                }
            }

            ThereAreMdispersives = ThereAreMdispersives || (conta != 0);
            MDutton.Medium[jmed].NumNodesHz = conta;
            MDutton.Medium[jmed].NodesHz.resize(conta + 1);
            for (int i1 = 1; i1 <= conta; ++i1) {
                int n = sgg.Med[tempindex].Mdispersive[1].numpolres11;
                MDutton.Medium[jmed].NodesHz[i1].Current.resize(n + 1);
            }
            
            conta = 0;
            for (int k1 = sgg.Sweep[iHz].ZI; k1 <= sgg.Sweep[iHz].ZE; ++k1) {
                for (int j1 = sgg.Sweep[iHz].YI; j1 <= sgg.Sweep[iHz].YE; ++j1) {
                    for (int i1 = sgg.Sweep[iHz].XI; i1 <= sgg.Sweep[iHz].XE; ++i1) {
                        if (media.sggMiHz[i1][j1][k1] == tempindex) {
                            conta++;
                            MDutton.Medium[jmed].NodesHz[conta].i = i1;
                            MDutton.Medium[jmed].NodesHz[conta].j = j1;
                            MDutton.Medium[jmed].NodesHz[conta].k = k1;
                            MDutton.Medium[jmed].NodesHz[conta].WhatField = iHz;
                            MDutton.Medium[jmed].NodesHz[conta].FieldPresent = &Hz[i1][j1][k1];
                        }
                    }
                }
            }
        }

        // Resume or start
        for (int jmed = 1; jmed <= MDutton.NumMdispersives; ++jmed) {
            int numpolres = sgg.Med[MDutton.Medium[jmed].indexmed].Mdispersive[1].numpolres11;
            if (!resume) {
                // Hx, Jx
                for (int i1 = 1; i1 <= MDutton.Medium[jmed].NumNodesHx; ++i1) {
                    MDutton.Medium[jmed].NodesHx[i1].FieldPrevious = 0.0;
                    for (int k1 = 1; k1 <= numpolres; ++k1) {
                        MDutton.Medium[jmed].NodesHx[i1].Current[k1] = 0.0;
                    }
                }
                // Hy, Jy
                for (int i1 = 1; i1 <= MDutton.Medium[jmed].NumNodesHy; ++i1) {
                    MDutton.Medium[jmed].NodesHy[i1].FieldPrevious = 0.0;
                    for (int k1 = 1; k1 <= numpolres; ++k1) {
                        MDutton.Medium[jmed].NodesHy[i1].Current[k1] = 0.0;
                    }
                }
                // Hz, Jz
                for (int i1 = 1; i1 <= MDutton.Medium[jmed].NumNodesHz; ++i1) {
                    MDutton.Medium[jmed].NodesHz[i1].FieldPrevious = 0.0;
                    for (int k1 = 1; k1 <= numpolres; ++k1) {
                        MDutton.Medium[jmed].NodesHz[i1].Current[k1] = 0.0;
                    }
                }
            } else {
                // Hx, Jx
                for (int i1 = 1; i1 <= MDutton.Medium[jmed].NumNodesHx; ++i1) {
                    std::cin >> MDutton.Medium[jmed].NodesHx[i1].FieldPrevious; // Assuming unit 14 is stdin or a file stream
                    for (int k1 = 1; k1 <= numpolres; ++k1) {
                        std::cin >> MDutton.Medium[jmed].NodesHx[i1].Current[k1];
                    }
                }
                // Hy, Jy
                for (int i1 = 1; i1 <= MDutton.Medium[jmed].NumNodesHy; ++i1) {
                    std::cin >> MDutton.Medium[jmed].NodesHy[i1].FieldPrevious;
                    for (int k1 = 1; k1 <= numpolres; ++k1) {
                        std::cin >> MDutton.Medium[jmed].NodesHy[i1].Current[k1];
                    }
                }
                // Hz, Jz
                for (int i1 = 1; i1 <= MDutton.Medium[jmed].NumNodesHz; ++i1) {
                    std::cin >> MDutton.Medium[jmed].NodesHz[i1].FieldPrevious;
                    for (int k1 = 1; k1 <= numpolres; ++k1) {
                        std::cin >> MDutton.Medium[jmed].NodesHz[i1].Current[k1];
                    }
                }
            }
        }
    }

    void AdvanceMdispersiveH(FDETYPES_m::SGGFDTDINFO_t& sgg) {
        for (int jmed = 1; jmed <= MDutton.NumMdispersives; ++jmed) {
            int numpolres = MDutton.Medium[jmed].numpolres11;
            
            // Hx, Jx
            for (int i1 = 1; i1 <= MDutton.Medium[jmed].NumNodesHx; ++i1) {
                field_t& tempnode = MDutton.Medium[jmed].NodesHx[i1];
                for (int k1 = 1; k1 <= numpolres; ++k1) {
                    tempnode.FieldPresent = tempnode.FieldPresent - std::real(MDutton.Medium[jmed].GM3[k1] * tempnode.Current[k1]);
                }
                for (int k1 = 1; k1 <= numpolres; ++k1) {
                    tempnode.Current[k1] = MDutton.Medium[jmed].Kappa[k1] * tempnode.Current[k1] + 
                                           MDutton.Medium[jmed].Beta[k1] / sgg.dt * (tempnode.FieldPresent - tempnode.FieldPrevious);
                }
                tempnode.FieldPrevious = tempnode.FieldPresent;
            }
            
            // Hy, Jy
            for (int i1 = 1; i1 <= MDutton.Medium[jmed].NumNodesHy; ++i1) {
                field_t& tempnode = MDutton.Medium[jmed].NodesHy[i1];
                for (int k1 = 1; k1 <= numpolres; ++k1) {
                    tempnode.FieldPresent = tempnode.FieldPresent - std::real(MDutton.Medium[jmed].GM3[k1] * tempnode.Current[k1]);
                }
                for (int k1 = 1; k1 <= numpolres; ++k1) {
                    tempnode.Current[k1] = MDutton.Medium[jmed].Kappa[k1] * tempnode.Current[k1] + 
                                           MDutton.Medium[jmed].Beta[k1] / sgg.dt * (tempnode.FieldPresent - tempnode.FieldPrevious);
                }
                tempnode.FieldPrevious = tempnode.FieldPresent;
            }

            // Hz, Jz
            for (int i1 = 1; i1 <= MDutton.Medium[jmed].NumNodesHz; ++i1) {
                field_t& tempnode = MDutton.Medium[jmed].NodesHz[i1];
                for (int k1 = 1; k1 <= numpolres; ++k1) {
                    tempnode.FieldPresent = tempnode.FieldPresent - std::real(MDutton.Medium[jmed].GM3[k1] * tempnode.Current[k1]);
                }
                for (int k1 = 1; k1 <= numpolres; ++k1) {
                    tempnode.Current[k1] = MDutton.Medium[jmed].Kappa[k1] * tempnode.Current[k1] + 
                                           MDutton.Medium[jmed].Beta[k1] / sgg.dt * (tempnode.FieldPresent - tempnode.FieldPrevious);
                }
                tempnode.FieldPrevious = tempnode.FieldPresent;
            }
        }
    }

    void StoreFieldsMdispersives() {
        for (int jmed = 1; jmed <= MDutton.NumMdispersives; ++jmed) {
            int numpolres = MDutton.Medium[jmed].numpolres11;
            
            // Hx, Jx
            for (int i1 = 1; i1 <= MDutton.Medium[jmed].NumNodesHx; ++i1) {
                std::cout << MDutton.Medium[jmed].NodesHx[i1].FieldPrevious << std::endl;
                for (int k1 = 1; k1 <= numpolres; ++k1) {
                    std::cout << MDutton.Medium[jmed].NodesHx[i1].Current[k1] << std::endl;
                }
            }
            // Hy, Jy
            for (int i1 = 1; i1 <= MDutton.Medium[jmed].NumNodesHy; ++i1) {
                std::cout << MDutton.Medium[jmed].NodesHy[i1].FieldPrevious << std::endl;
                for (int k1 = 1; k1 <= numpolres; ++k1) {
                    std::cout << MDutton.Medium[jmed].NodesHy[i1].Current[k1] << std::endl;
                }
            }
            // Hz, Jz
            for (int i1 = 1; i1 <= MDutton.Medium[jmed].NumNodesHz; ++i1) {
                std::cout << MDutton.Medium[jmed].NodesHz[i1].FieldPrevious << std::endl;
                for (int k1 = 1; k1 <= numpolres; ++k1) {
                    std::cout << MDutton.Medium[jmed].NodesHz[i1].Current[k1] << std::endl;
                }
            }
        }
    }

    void DestroyMdispersives(FDETYPES_m::SGGFDTDINFO_t& sgg) {
        // Free up memory
        for (int i = 1; i <= sgg.NumMedia; ++i) {
            if (sgg.Med[i].Is.Mdispersive && !sgg.Med[i].Is.PML && !sgg.Med[i].Is.MdispersiveANIS) {
                // Assuming Mdispersive is a vector, index 1 is the first element
                if (sgg.Med[i].Mdispersive.size() > 1) {
                    sgg.Med[i].Mdispersive[1].c11.clear();
                    sgg.Med[i].Mdispersive[1].a11.clear();
                }
            }
        }
        for (int i = 1; i <= sgg.NumMedia; ++i) {
            if (sgg.Med[i].Is.Mdispersive && !sgg.Med[i].Is.PML && !sgg.Med[i].Is.MdispersiveANIS) {
                sgg.Med[i].Mdispersive.clear();
            }
        }

        for (int jmed = 1; jmed <= MDutton.NumMdispersives; ++jmed) {
            MDutton.Medium[jmed].Beta.clear();
            MDutton.Medium[jmed].Kappa.clear();
            MDutton.Medium[jmed].GM3.clear();

            for (int i1 = 1; i1 <= MDutton.Medium[jmed].NumNodesHx; ++i1) {
                MDutton.Medium[jmed].NodesHx[i1].Current.clear();
            }
            MDutton.Medium[jmed].NodesHx.clear();
            
            for (int i1 = 1; i1 <= MDutton.Medium[jmed].NumNodesHy; ++i1) {
                MDutton.Medium[jmed].NodesHy[i1].Current.clear();
            }
            MDutton.Medium[jmed].NodesHy.clear();
            
            for (int i1 = 1; i1 <= MDutton.Medium[jmed].NumNodesHz; ++i1) {
                MDutton.Medium[jmed].NodesHz[i1].Current.clear();
            }
            MDutton.Medium[jmed].NodesHz.clear();
        }

        MDutton.Medium.clear();
    }

}