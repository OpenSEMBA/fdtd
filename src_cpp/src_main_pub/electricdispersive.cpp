#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <string>

// Forward declarations for types defined in other modules
// These would typically be in separate headers like FDETYPES_m.h and Report_m.h

// Placeholder for types from FDETYPES_m
// Assuming RKIND is double, CKIND is std::complex<double>
// Assuming SGGFDTDINFO_t, media_matrices_t, etc. are defined elsewhere
struct SGGFDTDINFO_t;
struct media_matrices_t;

// Placeholder for types from Report_m
void print11(int unit, const std::string& msg);
const std::string SEPARADOR = "========================================";

namespace EDispersives_m {

    // Assuming these constants/types are defined in FDETYPES_m
    using RKIND = double;
    using CKIND = std::complex<double>;
    
    // Assuming iEx, iEy, iEz are integer constants defined elsewhere
    extern const int iEx;
    extern const int iEy;
    extern const int iEz;
    
    // Assuming NumPolRes is a global constant or part of a config struct
    // In the original code, NumPolRes is used without qualification inside the module,
    // likely a global parameter or part of the context. 
    // Given the context of "numpolres11" in the struct, NumPolRes might be a typo in Fortran 
    // referring to numpolres11 or a global constant. 
    // Looking at usage: "Do i1=1,NumPolRes" inside loops where numpolres is calculated.
    // It seems NumPolRes is intended to be the number of poles for the current medium.
    // However, in the loop "Do i1=1,NumPolRes", it is used before being assigned to numpolres in some contexts?
    // No, in InitEDispersives:
    // numpolres=sgg%Med(tempindex)%EDispersive(1)%numpolres11
    // Do i1=1,NumPolRes ...
    // This looks like a bug in the original Fortran code if NumPolRes is not defined.
    // It likely should be numpolres. I will assume it refers to the local variable numpolres 
    // or a global constant. Given the context, it's safer to assume it's a typo for numpolres 
    // or a global constant. Let's look at AdvanceEDispersiveE:
    // numpolres=Dutton%Medium(jmed)%numpolres11
    // Do k1=1,NumPolRes ...
    // Again, likely numpolres. I will replace NumPolRes with the local numpolres variable 
    // where it appears in loops, assuming it was a variable name collision or typo in the source.
    // Actually, looking closely, NumPolRes is never declared. It is likely a typo for numpolres.
    // I will use numpolres in the C++ translation.

    struct field_t {
        int32_t i;
        int32_t j;
        int32_t k;
        int32_t WhatField;

        RKIND* FieldPresent; // Pointer to background field
        RKIND FieldPrevious;
        std::vector<CKIND> Current;
    };

    struct EDispersive_t {
        int32_t indexmed;
        int32_t numnodesEx;
        int32_t numnodesEy;
        int32_t numnodesEz;
        int32_t numpolres11;
        
        std::vector<CKIND> Beta;
        std::vector<CKIND> Kappa;
        std::vector<CKIND> G3;
        
        std::vector<field_t> NodesEx;
        std::vector<field_t> NodesEy;
        std::vector<field_t> NodesEz;
    };

    struct EDispersive2_t {
        int32_t NumEDispersives;
        std::vector<EDispersive_t> Medium;
    };

    // Global save variable
    EDispersive2_t Dutton;

    void InitEDispersives(const SGGFDTDINFO_t& sgg, 
                          const media_matrices_t& media, 
                          std::vector<RKIND>& G1, 
                          std::vector<RKIND>& G2, 
                          std::vector<std::vector<std::vector<RKIND>>>& Ex, 
                          std::vector<std::vector<std::vector<RKIND>>>& Ey, 
                          std::vector<std::vector<std::vector<RKIND>>>& Ez,
                          bool& ThereAreEDispersives, 
                          bool resume) {
        
        ThereAreEDispersives = false;
        int conta = 0;
        
        // Count dispersive media
        for (int jmed = 1; jmed <= sgg.NumMedia; ++jmed) {
            if (sgg.Med[jmed].Is.EDispersive && !sgg.Med[jmed].Is.EDispersiveAnis) {
                conta++;
            }
        }

        Dutton.NumEDispersives = conta;
        Dutton.Medium.resize(Dutton.NumEDispersives);
        
        conta = 0;
        for (int jmed = 1; jmed <= sgg.NumMedia; ++jmed) {
            if (sgg.Med[jmed].Is.EDispersive && !sgg.Med[jmed].Is.EDispersiveAnis) {
                conta++;
                Dutton.Medium[conta - 1].indexmed = jmed;
                Dutton.Medium[conta - 1].numpolres11 = sgg.Med[jmed].EDispersive[1].numpolres11;
                
                int n = sgg.Med[jmed].EDispersive[1].numpolres11;
                Dutton.Medium[conta - 1].Beta.resize(n, 0.0);
                Dutton.Medium[conta - 1].Kappa.resize(n, 0.0);
                Dutton.Medium[conta - 1].G3.resize(n, 0.0);
                
                for (int i1 = 1; i1 <= n; ++i1) {
                    Dutton.Medium[conta - 1].Kappa[i1 - 1] = (1.0 + sgg.Med[jmed].EDispersive[1].a11[i1 - 1] * sgg.dt / 2.0) /
                                                             (1.0 - sgg.Med[jmed].EDispersive[1].a11[i1 - 1] * sgg.dt / 2.0);
                    Dutton.Medium[conta - 1].Beta[i1 - 1] = (sgg.Med[jmed].EDispersive[1].C11[i1 - 1] * sgg.dt) /
                                                            (1.0 - sgg.Med[jmed].EDispersive[1].a11[i1 - 1] * sgg.dt / 2.0);
                }
            }
        }

        // Calculate coefficients
        for (int jmed = 1; jmed <= Dutton.NumEDispersives; ++jmed) {
            int tempindex = Dutton.Medium[jmed - 1].indexmed;
            int numpolres = sgg.Med[tempindex].EDispersive[1].numpolres11;
            RKIND tempo = 0.0;
            
            for (int i1 = 1; i1 <= numpolres; ++i1) {
                tempo += std::real(Dutton.Medium[jmed - 1].Beta[i1 - 1]);
            }
            
            G1[tempindex] = (2.0 * sgg.Med[tempindex].Edispersive[1].eps11 + tempo - sgg.Med[tempindex].Edispersive[1].Sigma11 * sgg.dt) /
                            (2.0 * sgg.Med[tempindex].Edispersive[1].eps11 + tempo + sgg.Med[tempindex].Edispersive[1].Sigma11 * sgg.dt);
            G2[tempindex] = 2.0 * sgg.dt / (2.0 * sgg.Med[tempindex].Edispersive[1].eps11 + tempo + sgg.Med[tempindex].Edispersive[1].Sigma11 * sgg.dt);
            
            for (int i1 = 1; i1 <= numpolres; ++i1) {
                Dutton.Medium[jmed - 1].G3[i1 - 1] = G2[tempindex] / 2.0 * (1.0 + Dutton.Medium[jmed - 1].Kappa[i1 - 1]);
            }
        }

        // Allocate and populate Nodes
        for (int jmed = 1; jmed <= Dutton.NumEDispersives; ++jmed) {
            int tempindex = Dutton.Medium[jmed - 1].indexmed;
            
            // Ex
            int contaEx = 0;
            for (int k1 = sgg.Sweep[iEx].ZI; k1 <= sgg.Sweep[iEx].ZE; ++k1) {
                for (int j1 = sgg.Sweep[iEx].YI; j1 <= sgg.Sweep[iEx].YE; ++j1) {
                    for (int i1 = sgg.Sweep[iEx].XI; i1 <= sgg.Sweep[iEx].XE; ++i1) {
                        if (media.sggMiEx[i1][j1][k1] == tempindex) {
                            contaEx++;
                        }
                    }
                }
            }
            
            ThereAreEDispersives = ThereAreEDispersives || (contaEx != 0);
            Dutton.Medium[jmed - 1].NumNodesEx = contaEx;
            Dutton.Medium[jmed - 1].NodesEx.resize(contaEx);
            
            for (int i1 = 0; i1 < contaEx; ++i1) {
                Dutton.Medium[jmed - 1].NodesEx[i1].Current.resize(sgg.Med[tempindex].EDispersive[1].numpolres11);
            }
            
            int contaIdx = 0;
            for (int k1 = sgg.Sweep[iEx].ZI; k1 <= sgg.Sweep[iEx].ZE; ++k1) {
                for (int j1 = sgg.Sweep[iEx].YI; j1 <= sgg.Sweep[iEx].YE; ++j1) {
                    for (int i1 = sgg.Sweep[iEx].XI; i1 <= sgg.Sweep[iEx].XE; ++i1) {
                        if (media.sggMiEx[i1][j1][k1] == tempindex) {
                            contaIdx++;
                            Dutton.Medium[jmed - 1].NodesEx[contaIdx - 1].i = i1;
                            Dutton.Medium[jmed - 1].NodesEx[contaIdx - 1].j = j1;
                            Dutton.Medium[jmed - 1].NodesEx[contaIdx - 1].k = k1;
                            Dutton.Medium[jmed - 1].NodesEx[contaIdx - 1].WhatField = iEx;
                            Dutton.Medium[jmed - 1].NodesEx[contaIdx - 1].FieldPresent = &Ex[i1][j1][k1];
                        }
                    }
                }
            }
            
            // Ey
            int contaEy = 0;
            for (int k1 = sgg.Sweep[iEy].ZI; k1 <= sgg.Sweep[iEy].ZE; ++k1) {
                for (int j1 = sgg.Sweep[iEy].YI; j1 <= sgg.Sweep[iEy].YE; ++j1) {
                    for (int i1 = sgg.Sweep[iEy].XI; i1 <= sgg.Sweep[iEy].XE; ++i1) {
                        if (media.sggMiEy[i1][j1][k1] == tempindex) {
                            contaEy++;
                        }
                    }
                }
            }
            
            ThereAreEDispersives = ThereAreEDispersives || (contaEy != 0);
            Dutton.Medium[jmed - 1].NumNodesEy = contaEy;
            Dutton.Medium[jmed - 1].NodesEy.resize(contaEy);
            
            for (int i1 = 0; i1 < contaEy; ++i1) {
                Dutton.Medium[jmed - 1].NodesEy[i1].Current.resize(sgg.Med[tempindex].EDispersive[1].numpolres11);
            }
            
            contaIdx = 0;
            for (int k1 = sgg.Sweep[iEy].ZI; k1 <= sgg.Sweep[iEy].ZE; ++k1) {
                for (int j1 = sgg.Sweep[iEy].YI; j1 <= sgg.Sweep[iEy].YE; ++j1) {
                    for (int i1 = sgg.Sweep[iEy].XI; i1 <= sgg.Sweep[iEy].XE; ++i1) {
                        if (media.sggMiEy[i1][j1][k1] == tempindex) {
                            contaIdx++;
                            Dutton.Medium[jmed - 1].NodesEy[contaIdx - 1].i = i1;
                            Dutton.Medium[jmed - 1].NodesEy[contaIdx - 1].j = j1;
                            Dutton.Medium[jmed - 1].NodesEy[contaIdx - 1].k = k1;
                            Dutton.Medium[jmed - 1].NodesEy[contaIdx - 1].WhatField = iEy;
                            Dutton.Medium[jmed - 1].NodesEy[contaIdx - 1].FieldPresent = &Ey[i1][j1][k1];
                        }
                    }
                }
            }
            
            // Ez
            int contaEz = 0;
            for (int k1 = sgg.Sweep[iEz].ZI; k1 <= sgg.Sweep[iEz].ZE; ++k1) {
                for (int j1 = sgg.Sweep[iEz].YI; j1 <= sgg.Sweep[iEz].YE; ++j1) {
                    for (int i1 = sgg.Sweep[iEz].XI; i1 <= sgg.Sweep[iEz].XE; ++i1) {
                        if (media.sggMiEz[i1][j1][k1] == tempindex) {
                            contaEz++;
                        }
                    }
                }
            }
            
            ThereAreEDispersives = ThereAreEDispersives || (contaEz != 0);
            Dutton.Medium[jmed - 1].NumNodesEz = contaEz;
            Dutton.Medium[jmed - 1].NodesEz.resize(contaEz);
            
            for (int i1 = 0; i1 < contaEz; ++i1) {
                Dutton.Medium[jmed - 1].NodesEz[i1].Current.resize(sgg.Med[tempindex].EDispersive[1].numpolres11);
            }
            
            contaIdx = 0;
            for (int k1 = sgg.Sweep[iEz].ZI; k1 <= sgg.Sweep[iEz].ZE; ++k1) {
                for (int j1 = sgg.Sweep[iEz].YI; j1 <= sgg.Sweep[iEz].YE; ++j1) {
                    for (int i1 = sgg.Sweep[iEz].XI; i1 <= sgg.Sweep[iEz].XE; ++i1) {
                        if (media.sggMiEz[i1][j1][k1] == tempindex) {
                            contaIdx++;
                            Dutton.Medium[jmed - 1].NodesEz[contaIdx - 1].i = i1;
                            Dutton.Medium[jmed - 1].NodesEz[contaIdx - 1].j = j1;
                            Dutton.Medium[jmed - 1].NodesEz[contaIdx - 1].k = k1;
                            Dutton.Medium[jmed - 1].NodesEz[contaIdx - 1].WhatField = iEz;
                            Dutton.Medium[jmed - 1].NodesEz[contaIdx - 1].FieldPresent = &Ez[i1][j1][k1];
                        }
                    }
                }
            }
        }

        // Resume or start
        for (int jmed = 1; jmed <= Dutton.NumEDispersives; ++jmed) {
            int tempindex = Dutton.Medium[jmed - 1].indexmed;
            int numpolres = sgg.Med[tempindex].EDispersive[1].numpolres11;
            
            if (!resume) {
                // Ex, Jx
                for (int i1 = 0; i1 < Dutton.Medium[jmed - 1].NumNodesEx; ++i1) {
                    Dutton.Medium[jmed - 1].NodesEx[i1].fieldPrevious = 0.0;
                    for (int k1 = 0; k1 < numpolres; ++k1) {
                        Dutton.Medium[jmed - 1].NodesEx[i1].current[k1] = 0.0;
                    }
                }
                // Ey, Jy
                for (int i1 = 0; i1 < Dutton.Medium[jmed - 1].NumNodesEy; ++i1) {
                    Dutton.Medium[jmed - 1].NodesEy[i1].fieldPrevious = 0.0;
                    for (int k1 = 0; k1 < numpolres; ++k1) {
                        Dutton.Medium[jmed - 1].NodesEy[i1].current[k1] = 0.0;
                    }
                }
                // Ez, Jz
                for (int i1 = 0; i1 < Dutton.Medium[jmed - 1].NumNodesEz; ++i1) {
                    Dutton.Medium[jmed - 1].NodesEz[i1].fieldPrevious = 0.0;
                    for (int k1 = 0; k1 < numpolres; ++k1) {
                        Dutton.Medium[jmed - 1].NodesEz[i1].current[k1] = 0.0;
                    }
                }
            } else {
                // Ex, Jx
                for (int i1 = 0; i1 < Dutton.Medium[jmed - 1].NumNodesEx; ++i1) {
                    std::cin >> Dutton.Medium[jmed - 1].NodesEx[i1].fieldPrevious; // Assuming unit 14 is stdin or a file stream
                    for (int k1 = 0; k1 < numpolres; ++k1) {
                        std::cin >> Dutton.Medium[jmed - 1].NodesEx[i1].current[k1];
                    }
                }
                // Ey, Jy
                for (int i1 = 0; i1 < Dutton.Medium[jmed - 1].NumNodesEy; ++i1) {
                    std::cin >> Dutton.Medium[jmed - 1].NodesEy[i1].fieldPrevious;
                    for (int k1 = 0; k1 < numpolres; ++k1) {
                        std::cin >> Dutton.Medium[jmed - 1].NodesEy[i1].current[k1];
                    }
                }
                // Ez, Jz
                for (int i1 = 0; i1 < Dutton.Medium[jmed - 1].NumNodesEz; ++i1) {
                    std::cin >> Dutton.Medium[jmed - 1].NodesEz[i1].fieldPrevious;
                    for (int k1 = 0; k1 < numpolres; ++k1) {
                        std::cin >> Dutton.Medium[jmed - 1].NodesEz[i1].current[k1];
                    }
                }
            }
        }
    }

    void AdvanceEDispersiveE(const SGGFDTDINFO_t& sgg) {
        for (int jmed = 1; jmed <= Dutton.NumEDispersives; ++jmed) {
            int numpolres = Dutton.Medium[jmed - 1].numpolres11;
            
            // Ex, Jx
            for (int i1 = 0; i1 < Dutton.Medium[jmed - 1].NumNodesEx; ++i1) {
                field_t& tempnode = Dutton.Medium[jmed - 1].NodesEx[i1];
                
                for (int k1 = 0; k1 < numpolres; ++k1) {
                    tempnode.fieldPresent = *tempnode.FieldPresent - std::real(Dutton.Medium[jmed - 1].G3[k1] * tempnode.current[k1]);
                }
                
                for (int k1 = 0; k1 < numpolres; ++k1) {
                    tempnode.current[k1] = Dutton.Medium[jmed - 1].Kappa[k1] * tempnode.current[k1] +
                                           Dutton.Medium[jmed - 1].Beta[k1] * (tempnode.fieldPresent - tempnode.fieldPrevious) / sgg.dt;
                }
                tempnode.fieldPrevious = tempnode.fieldPresent;
            }
            
            // Ey, Jy
            for (int i1 = 0; i1 < Dutton.Medium[jmed - 1].NumNodesEy; ++i1) {
                field_t& tempnode = Dutton.Medium[jmed - 1].NodesEy[i1];
                
                for (int k1 = 0; k1 < numpolres; ++k1) {
                    tempnode.FieldPresent = *tempnode.FieldPresent - std::real(Dutton.Medium[jmed - 1].G3[k1] * tempnode.current[k1]);
                }
                
                for (int k1 = 0; k1 < numpolres; ++k1) {
                    tempnode.current[k1] = Dutton.Medium[jmed - 1].Kappa[k1] * tempnode.current[k1] +
                                           Dutton.Medium[jmed - 1].Beta[k1] * (tempnode.FieldPresent - tempnode.fieldPrevious) / sgg.dt;
                }
                tempnode.fieldPrevious = tempnode.FieldPresent;
            }
            
            // Ez, Jz
            for (int i1 = 0; i1 < Dutton.Medium[jmed - 1].NumNodesEz; ++i1) {
                field_t& tempnode = Dutton.Medium[jmed - 1].NodesEz[i1];
                
                for (int k1 = 0; k1 < numpolres; ++k1) {
                    tempnode.FieldPresent = *tempnode.FieldPresent - std::real(Dutton.Medium[jmed - 1].G3[k1] * tempnode.current[k1]);
                }
                
                for (int k1 = 0; k1 < numpolres; ++k1) {
                    tempnode.current[k1] = Dutton.Medium[jmed - 1].Kappa[k1] * tempnode.current[k1] +
                                           Dutton.Medium[jmed - 1].Beta[k1] * (tempnode.FieldPresent - tempnode.fieldPrevious) / sgg.dt;
                }
                tempnode.fieldPrevious = tempnode.FieldPresent;
            }
        }
    }

    void StoreFieldsEDispersives() {
        for (int jmed = 1; jmed <= Dutton.NumEDispersives; ++jmed) {
            int numpolres = Dutton.Medium[jmed - 1].numpolres11;
            
            // Ex, Jx
            for (int i1 = 0; i1 < Dutton.Medium[jmed - 1].NumNodesEx; ++i1) {
                std::cout << Dutton.Medium[jmed - 1].NodesEx[i1].fieldPrevious << std::endl;
                for (int k1 = 0; k1 < numpolres; ++k1) {
                    std::cout << Dutton.Medium[jmed - 1].NodesEx[i1].current[k1] << std::endl;
                }
            }
            // Ey, Jy
            for (int i1 = 0; i1 < Dutton.Medium[jmed - 1].NumNodesEy; ++i1) {
                std::cout << Dutton.Medium[jmed - 1].NodesEy[i1].fieldPrevious << std::endl;
                for (int k1 = 0; k1 < numpolres; ++k1) {
                    std::cout << Dutton.Medium[jmed - 1].NodesEy[i1].current[k1] << std::endl;
                }
            }
            // Ez, Jz
            for (int i1 = 0; i1 < Dutton.Medium[jmed - 1].NumNodesEz; ++i1) {
                std::cout << Dutton.Medium[jmed - 1].NodesEz[i1].fieldPrevious << std::endl;
                for (int k1 = 0; k1 < numpolres; ++k1) {
                    std::cout << Dutton.Medium[jmed - 1].NodesEz[i1].current[k1] << std::endl;
                }
            }
        }
    }

    void DestroyEDispersives(SGGFDTDINFO_t& sgg) {
        // Free up memory for sgg.Med
        for (int i = 1; i <= sgg.NumMedia; ++i) {
            if (sgg.Med[i].Is.EDispersive && !sgg.Med[i].Is.PML && !sgg.Med[i].Is.EDispersiveAnis) {
                sgg.Med[i].EDispersive[1].C11.clear();
                sgg.Med[i].EDispersive[1].a11.clear();
            }
        }
        for (int i = 1; i <= sgg.NumMedia; ++i) {
            if (sgg.Med[i].Is.EDispersive && !sgg.Med[i].Is.PML && !sgg.Med[i].Is.EDispersiveAnis) {
                sgg.Med[i].EDispersive.clear();
            }
        }

        // Free up memory for Dutton
        for (int jmed = 1; jmed <= Dutton.NumEDispersives; ++jmed) {
            Dutton.Medium[jmed - 1].Beta.clear();
            Dutton.Medium[jmed - 1].Kappa.clear();
            Dutton.Medium[jmed - 1].G3.clear();

            for (int i1 = 0; i1 < Dutton.Medium[jmed - 1].NumNodesEx; ++i1) {
                Dutton.Medium[jmed - 1].NodesEx[i1].Current.clear();
            }
            Dutton.Medium[jmed - 1].NodesEx.clear();
            
            for (int i1 = 0; i1 < Dutton.Medium[jmed - 1].NumNodesEy; ++i1) {
                Dutton.Medium[jmed - 1].NodesEy[i1].Current.clear();
            }
            Dutton.Medium[jmed - 1].NodesEy.clear();
            
            for (int i1 = 0; i1 < Dutton.Medium[jmed - 1].NumNodesEz; ++i1) {
                Dutton.Medium[jmed - 1].NodesEz[i1].Current.clear();
            }
            Dutton.Medium[jmed - 1].NodesEz.clear();
        }
        
        Dutton.Medium.clear();
    }

}