#include <vector>
#include <memory>
#include <string>
#include <iostream>

// Assuming FDETYPES_m provides these types/constants. 
// Since the full definition isn't provided, we define placeholders or assume standard equivalents.
// In a real scenario, these would come from the translated FDETYPES_m.

using rkind = double;
using RKIND = double;

enum FieldIndex { iEx, iEy, iEz, iHx, iHy, iHz };

struct SGGFDTDINFO_t {
    int NumMedia;
    // Placeholder for Media structure. 
    // The original code accesses sgg%Med(jmed)%Is%Anisotropic, etc.
    // We assume a structure that mimics this hierarchy.
    struct MediaInfo {
        struct IsFlags {
            bool Anisotropic;
            bool ThinSlot;
        } Is;
        struct AnisotropicData {
            double sigma[3][3];
            double sigmam[3][3];
            double mur[3][3];
            double epr[3][3];
        } Anisotropic[1]; // Assuming size 1 based on usage Anisotropic(1)
        
        // Placeholder for other members needed by SGG
        struct SINPMLSweep_t {
            int XI, XE, YI, YE, ZI, ZE;
        } SINPMLSweep[6]; // Assuming 6 fields: Ex, Ey, Ez, Hx, Hy, Hz
    };
    std::vector<MediaInfo> Med;
    
    struct SharedElement {
        int i, j, k;
        int Field; // FieldIndex
        int Times;
        int SharedMed; // Assuming integer index or ID
    };
    
    struct SharedList {
        int conta;
        std::vector<SharedElement> elem;
    } Eshared, Hshared;
};

struct media_matrices_t {
    // Placeholder for media matrices. 
    // Original code accesses media%sggMiEx(i1,j1,k1) etc.
    // Assuming 3D arrays or flattened vectors.
    std::vector<std::vector<std::vector<int>>> sggMiEx;
    std::vector<std::vector<std::vector<int>>> sggMiEy;
    std::vector<std::vector<std::vector<int>>> sggMiEz;
    std::vector<std::vector<std::vector<int>>> sggMiHx;
    std::vector<std::vector<std::vector<int>>> sggMiHy;
    std::vector<std::vector<std::vector<int>>> sggMiHz;
};

// Forward declarations for types used in Anisotropic_m
struct Coeff_t {
    double eexx, eexy, eexz, eeyx, eeyy, eeyz, eezx, eezy, eezz;
    double ehxx, ehxy, ehxz, ehyx, ehyy, ehyz, ehzx, ehzy, ehzz;
    double hexx, hexy, hexz, heyx, heyy, heyz, hezx, hezy, hezz;
    double hhxx, hhxy, hhxz, hhyx, hhyy, hhyz, hhzx, hhzy, hhzz;
};

struct LocalSharedElement_t {
    int times;
    std::vector<int> SharedMed;
    Coeff_t coeff;
};

struct Anisotropicinfo_t {
    int indexmed;
    int numnodesEx, numnodesEy, numnodesEz;
    int numnodesHx, numnodesHy, numnodesHz;
    
    std::vector<int> Ex_i, Ey_i, Ez_i;
    std::vector<int> Hx_i, Hy_i, Hz_i;
    
    std::vector<int> Ex_j, Ey_j, Ez_j;
    std::vector<int> Hx_j, Hy_j, Hz_j;
    
    std::vector<int> Ex_k, Ey_k, Ez_k;
    std::vector<int> Hx_k, Hy_k, Hz_k;
    
    std::vector<double> Ex_value, Ey_value, Ez_value;
    std::vector<double> Hx_value, Hy_value, Hz_value;
    
    std::vector<LocalSharedElement_t> Ex_Shared, Ey_Shared, Ez_Shared;
    std::vector<LocalSharedElement_t> Hx_Shared, Hy_Shared, Hz_Shared;
    
    Coeff_t coeff;
    bool IsOnlyThinSlot;
    
    double sigma[3][3];
    double epr[3][3];
    double mur[3][3];
    double sigmaM[3][3]; // sigmam in Fortran
};

struct AnisotropicMed_t {
    int NumMed;
    std::vector<Anisotropicinfo_t> info;
};

namespace Anisotropic_m {

    // Global variables
    AnisotropicMed_t AniMed;
    double eps0 = 0.0;
    double mu0 = 0.0;
    double cluz = 0.0;
    double zvac = 0.0;

    // Helper to access 3D array element (assuming row-major for C++)
    // Note: The original Fortran code uses 1-based indexing for loops but array access might be 0 or 1 based depending on declaration.
    // Assuming standard C++ 0-based indexing for the vectors defined below.
    // If the input data is 1-based, adjustments might be needed in the caller.
    // Here we assume the input vectors sggMi... are sized appropriately and accessed with 0-based indices 
    // corresponding to i1-1, j1-1, k1-1 if the Fortran code was 1-based.
    // However, to preserve logic exactly, we will assume the input arrays passed are already aligned or 
    // we access them directly. Given the loop `Do i1=...`, if Fortran arrays are 1-based, 
    // we must subtract 1. Let's assume standard FDTD grid where indices map directly to vector indices if 0-based.
    // If the Fortran code uses 1-based indexing for `sggMiEx(i1,j1,k1)`, then in C++ we use `sggMiEx[i1-1][j1-1][k1-1]`.
    // Since I don't have the definition of `media_matrices_t`'s internal storage, I will assume it's a 3D vector 
    // and the indices `i1, j1, k1` from the Fortran loop (which likely start at 1) need to be adjusted to 0-based.
    
    inline int get_media_index(const std::vector<std::vector<std::vector<int>>>& arr, int i, int j, int k) {
        // Assuming 0-based indexing in C++ vector. 
        // If Fortran was 1-based, i1, j1, k1 are 1-based.
        return arr[i-1][j-1][k-1];
    }

    void InitAnisotropic(const SGGFDTDINFO_t& sgg, const media_matrices_t& media, bool& ThereAreAnisotropic, bool& ThereAreThinSlot, double eps00, double mu00) {
        eps0 = eps00;
        mu0 = mu00;

        ThereAreAnisotropic = false;
        ThereAreThinSlot = false;
        int conta = 0;
        
        // Count anisotropic media
        for (int jmed = 1; jmed <= sgg.NumMedia; ++jmed) {
            if (sgg.Med[jmed-1].Is.Anisotropic) {
                conta++;
            }
        }

        AniMed.NumMed = conta;
        AniMed.info.resize(conta);
        
        conta = 0;
        for (int jmed = 1; jmed <= sgg.NumMedia; ++jmed) {
            if (sgg.Med[jmed-1].Is.Anisotropic) {
                conta++;
                AniMed.info[conta-1].indexmed = jmed; // 1-based index stored
                if (sgg.Med[jmed-1].Is.ThinSlot) {
                    AniMed.info[conta-1].IsOnlyThinSlot = true;
                } else {
                    AniMed.info[conta-1].IsOnlyThinSlot = false;
                }
            }
        }

        // Copy coefficients
        for (int jmed = 1; jmed <= AniMed.NumMed; ++jmed) {
            int idx = AniMed.info[jmed-1].indexmed;
            AniMed.info[jmed-1].sigma[0][0] = sgg.Med[idx-1].Anisotropic[0].sigma[0][0];
            // ... copying all sigma elements would be tedious. 
            // Assuming a copy function or loop. For brevity, I'll use a loop or memcpy if types match.
            // Since sigma is double[3][3], we can copy manually.
            for(int r=0; r<3; ++r)
                for(int c=0; c<3; ++c) {
                    AniMed.info[jmed-1].sigma[r][c] = sgg.Med[idx-1].Anisotropic[0].sigma[r][c];
                    AniMed.info[jmed-1].sigmam[r][c] = sgg.Med[idx-1].Anisotropic[0].sigmam[r][c];
                    AniMed.info[jmed-1].mur[r][c] = sgg.Med[idx-1].Anisotropic[0].mur[r][c];
                    AniMed.info[jmed-1].epr[r][c] = sgg.Med[idx-1].Anisotropic[0].epr[r][c];
                }
        }

        // Process each anisotropic medium
        for (int jmed = 1; jmed <= AniMed.NumMed; ++jmed) {
            int tempindex = AniMed.info[jmed-1].indexmed;
            
            // --- Ex ---
            conta = 0;
            // Assuming SINPMLSweep is indexed 0..5 corresponding to iEx..iHz
            // iEx is likely 0, iEy 1, etc.
            const auto& sweepEx = sgg.SINPMLSweep[0];
            for (int k1 = sweepEx.ZI; k1 <= sweepEx.ZE; ++k1) {
                for (int j1 = sweepEx.YI; j1 <= sweepEx.YE; ++j1) {
                    for (int i1 = sweepEx.XI; i1 <= sweepEx.XE; ++i1) {
                        if (get_media_index(media.sggMiEx, i1, j1, k1) == tempindex) {
                            conta++;
                        }
                    }
                }
            }
            ThereAreAnisotropic = ThereAreAnisotropic || (conta != 0);
            ThereAreThinSlot = ThereAreAnisotropic && AniMed.info[jmed-1].IsOnlyThinSlot;
            AniMed.info[jmed-1].numnodesEx = conta;
            
            AniMed.info[jmed-1].Ex_i.resize(conta);
            AniMed.info[jmed-1].Ex_j.resize(conta);
            AniMed.info[jmed-1].Ex_k.resize(conta);
            AniMed.info[jmed-1].Ex_value.resize(conta, 0.0);
            AniMed.info[jmed-1].Ex_Shared.resize(conta);
            
            conta = 0;
            for (int k1 = sweepEx.ZI; k1 <= sweepEx.ZE; ++k1) {
                for (int j1 = sweepEx.YI; j1 <= sweepEx.YE; ++j1) {
                    for (int i1 = sweepEx.XI; i1 <= sweepEx.XE; ++i1) {
                        if (get_media_index(media.sggMiEx, i1, j1, k1) == tempindex) {
                            conta++;
                            AniMed.info[jmed-1].Ex_Shared[conta-1].times = 1;
                            AniMed.info[jmed-1].Ex_i[conta-1] = i1;
                            AniMed.info[jmed-1].Ex_j[conta-1] = j1;
                            AniMed.info[jmed-1].Ex_k[conta-1] = k1;
                        }
                    }
                }
            }

            // --- Ey ---
            conta = 0;
            const auto& sweepEy = sgg.SINPMLSweep[1];
            for (int k1 = sweepEy.ZI; k1 <= sweepEy.ZE; ++k1) {
                for (int j1 = sweepEy.YI; j1 <= sweepEy.YE; ++j1) {
                    for (int i1 = sweepEy.XI; i1 <= sweepEy.XE; ++i1) {
                        if (get_media_index(media.sggMiEy, i1, j1, k1) == tempindex) {
                            conta++;
                        }
                    }
                }
            }
            ThereAreAnisotropic = ThereAreAnisotropic || (conta != 0);
            ThereAreThinSlot = ThereAreAnisotropic && AniMed.info[jmed-1].IsOnlyThinSlot;
            AniMed.info[jmed-1].numnodesEy = conta;
            
            AniMed.info[jmed-1].Ey_i.resize(conta);
            AniMed.info[jmed-1].Ey_j.resize(conta);
            AniMed.info[jmed-1].Ey_k.resize(conta);
            AniMed.info[jmed-1].Ey_value.resize(conta, 0.0);
            AniMed.info[jmed-1].Ey_Shared.resize(conta);
            
            conta = 0;
            for (int k1 = sweepEy.ZI; k1 <= sweepEy.ZE; ++k1) {
                for (int j1 = sweepEy.YI; j1 <= sweepEy.YE; ++j1) {
                    for (int i1 = sweepEy.XI; i1 <= sweepEy.XE; ++i1) {
                        if (get_media_index(media.sggMiEy, i1, j1, k1) == tempindex) {
                            conta++;
                            AniMed.info[jmed-1].Ey_Shared[conta-1].times = 1;
                            AniMed.info[jmed-1].Ey_i[conta-1] = i1;
                            AniMed.info[jmed-1].Ey_j[conta-1] = j1;
                            AniMed.info[jmed-1].Ey_k[conta-1] = k1;
                        }
                    }
                }
            }

            // --- Ez ---
            conta = 0;
            const auto& sweepEz = sgg.SINPMLSweep[2];
            for (int k1 = sweepEz.ZI; k1 <= sweepEz.ZE; ++k1) {
                for (int j1 = sweepEz.YI; j1 <= sweepEz.YE; ++j1) {
                    for (int i1 = sweepEz.XI; i1 <= sweepEz.XE; ++i1) {
                        if (get_media_index(media.sggMiEz, i1, j1, k1) == tempindex) {
                            conta++;
                        }
                    }
                }
            }
            ThereAreAnisotropic = ThereAreAnisotropic || (conta != 0);
            ThereAreThinSlot = ThereAreAnisotropic && AniMed.info[jmed-1].IsOnlyThinSlot;
            AniMed.info[jmed-1].numnodesEz = conta;
            
            AniMed.info[jmed-1].Ez_i.resize(conta);
            AniMed.info[jmed-1].Ez_j.resize(conta);
            AniMed.info[jmed-1].Ez_k.resize(conta);
            AniMed.info[jmed-1].Ez_value.resize(conta, 0.0);
            AniMed.info[jmed-1].Ez_Shared.resize(conta);
            
            conta = 0;
            for (int k1 = sweepEz.ZI; k1 <= sweepEz.ZE; ++k1) {
                for (int j1 = sweepEz.YI; j1 <= sweepEz.YE; ++j1) {
                    for (int i1 = sweepEz.XI; i1 <= sweepEz.XE; ++i1) {
                        if (get_media_index(media.sggMiEz, i1, j1, k1) == tempindex) {
                            conta++;
                            AniMed.info[jmed-1].Ez_Shared[conta-1].times = 1;
                            AniMed.info[jmed-1].Ez_i[conta-1] = i1;
                            AniMed.info[jmed-1].Ez_j[conta-1] = j1;
                            AniMed.info[jmed-1].Ez_k[conta-1] = k1;
                        }
                    }
                }
            }

            // --- Hx ---
            conta = 0;
            const auto& sweepHx = sgg.SINPMLSweep[3];
            for (int k1 = sweepHx.ZI; k1 <= sweepHx.ZE; ++k1) {
                for (int j1 = sweepHx.YI; j1 <= sweepHx.YE; ++j1) {
                    for (int i1 = sweepHx.XI; i1 <= sweepHx.XE; ++i1) {
                        if (get_media_index(media.sggMiHx, i1, j1, k1) == tempindex) {
                            conta++;
                        }
                    }
                }
            }
            ThereAreAnisotropic = ThereAreAnisotropic || (conta != 0);
            ThereAreThinSlot = ThereAreAnisotropic && AniMed.info[jmed-1].IsOnlyThinSlot;
            AniMed.info[jmed-1].numnodesHx = conta;
            
            AniMed.info[jmed-1].Hx_i.resize(conta);
            AniMed.info[jmed-1].Hx_j.resize(conta);
            AniMed.info[jmed-1].Hx_k.resize(conta);
            AniMed.info[jmed-1].Hx_value.resize(conta, 0.0);
            AniMed.info[jmed-1].Hx_Shared.resize(conta);
            
            conta = 0;
            for (int k1 = sweepHx.ZI; k1 <= sweepHx.ZE; ++k1) {
                for (int j1 = sweepHx.YI; j1 <= sweepHx.YE; ++j1) {
                    for (int i1 = sweepHx.XI; i1 <= sweepHx.XE; ++i1) {
                        if (get_media_index(media.sggMiHx, i1, j1, k1) == tempindex) {
                            conta++;
                            AniMed.info[jmed-1].Hx_Shared[conta-1].times = 1;
                            AniMed.info[jmed-1].Hx_i[conta-1] = i1;
                            AniMed.info[jmed-1].Hx_j[conta-1] = j1;
                            AniMed.info[jmed-1].Hx_k[conta-1] = k1;
                        }
                    }
                }
            }

            // --- Hy ---
            conta = 0;
            const auto& sweepHy = sgg.SINPMLSweep[4];
            for (int k1 = sweepHy.ZI; k1 <= sweepHy.ZE; ++k1) {
                for (int j1 = sweepHy.YI; j1 <= sweepHy.YE; ++j1) {
                    for (int i1 = sweepHy.XI; i1 <= sweepHy.XE; ++i1) {
                        if (get_media_index(media.sggMiHy, i1, j1, k1) == tempindex) {
                            conta++;
                        }
                    }
                }
            }
            ThereAreAnisotropic = ThereAreAnisotropic || (conta != 0);
            ThereAreThinSlot = ThereAreAnisotropic && AniMed.info[jmed-1].IsOnlyThinSlot;
            AniMed.info[jmed-1].numnodesHy = conta;
            
            AniMed.info[jmed-1].Hy_i.resize(conta);
            AniMed.info[jmed-1].Hy_j.resize(conta);
            AniMed.info[jmed-1].Hy_k.resize(conta);
            AniMed.info[jmed-1].Hy_value.resize(conta, 0.0);
            AniMed.info[jmed-1].Hy_Shared.resize(conta);
            
            conta = 0;
            for (int k1 = sweepHy.ZI; k1 <= sweepHy.ZE; ++k1) {
                for (int j1 = sweepHy.YI; j1 <= sweepHy.YE; ++j1) {
                    for (int i1 = sweepHy.XI; i1 <= sweepHy.XE; ++i1) {
                        if (get_media_index(media.sggMiHy, i1, j1, k1) == tempindex) {
                            conta++;
                            AniMed.info[jmed-1].Hy_Shared[conta-1].times = 1;
                            AniMed.info[jmed-1].Hy_i[conta-1] = i1;
                            AniMed.info[jmed-1].Hy_j[conta-1] = j1;
                            AniMed.info[jmed-1].Hy_k[conta-1] = k1;
                        }
                    }
                }
            }

            // --- Hz ---
            conta = 0;
            const auto& sweepHz = sgg.SINPMLSweep[5];
            for (int k1 = sweepHz.ZI; k1 <= sweepHz.ZE; ++k1) {
                for (int j1 = sweepHz.YI; j1 <= sweepHz.YE; ++j1) {
                    for (int i1 = sweepHz.XI; i1 <= sweepHz.XE; ++i1) {
                        if (get_media_index(media.sggMiHz, i1, j1, k1) == tempindex) {
                            conta++;
                        }
                    }
                }
            }
            ThereAreAnisotropic = ThereAreAnisotropic || (conta != 0);
            ThereAreThinSlot = ThereAreAnisotropic && AniMed.info[jmed-1].IsOnlyThinSlot;
            AniMed.info[jmed-1].numnodesHz = conta;
            
            AniMed.info[jmed-1].Hz_i.resize(conta);
            AniMed.info[jmed-1].Hz_j.resize(conta);
            AniMed.info[jmed-1].Hz_k.resize(conta);
            AniMed.info[jmed-1].Hz_value.resize(conta, 0.0);
            AniMed.info[jmed-1].Hz_Shared.resize(conta);
            
            conta = 0;
            for (int k1 = sweepHz.ZI; k1 <= sweepHz.ZE; ++k1) {
                for (int j1 = sweepHz.YI; j1 <= sweepHz.YE; ++j1) {
                    for (int i1 = sweepHz.XI; i1 <= sweepHz.XE; ++i1) {
                        if (get_media_index(media.sggMiHz, i1, j1, k1) == tempindex) {
                            conta++;
                            AniMed.info[jmed-1].Hz_Shared[conta-1].times = 1;
                            AniMed.info[jmed-1].Hz_i[conta-1] = i1;
                            AniMed.info[jmed-1].Hz_j[conta-1] = j1;
                            AniMed.info[jmed-1].Hz_k[conta-1] = k1;
                        }
                    }
                }
            }
        }

        // Update shared times and allocate SharedMed
        for (int j1 = 1; j1 <= sgg.Eshared.conta; ++j1) {
            bool found = false;
            for (int jmed = 1; jmed <= AniMed.NumMed && !found; ++jmed) {
                for (int i1 = 1; i1 <= AniMed.info[jmed-1].numnodesEx; ++i1) {
                    if ((sgg.Eshared.elem[j1-1].i == AniMed.info[jmed-1].Ex_i[i1-1]) &&
                        (sgg.Eshared.elem[j1-1].j == AniMed.info[jmed-1].Ex_j[i1-1]) &&
                        (sgg.Eshared.elem[j1-1].k == AniMed.info[jmed-1].Ex_k[i1-1]) &&
                        (sgg.Eshared.elem[j1-1].Field == iEx)) {
                        
                        AniMed.info[jmed-1].Ex_Shared[i1-1].times = sgg.Eshared.elem[j1-1].Times;
                        if (sgg.Eshared.elem[j1-1].Times > 1) {
                            AniMed.info[jmed-1].Ex_Shared[i1-1].SharedMed.resize(sgg.Eshared.elem[j1-1].Times);
                        }
                        found = true;
                        break;
                    }
                }
            }
            
            if (!found) {
                for (int jmed = 1; jmed <= AniMed.NumMed && !found; ++jmed) {
                    for (int i1 = 1; i1 <= AniMed.info[jmed-1].numnodesEy; ++i1) {
                        if ((sgg.Eshared.elem[j1-1].i == AniMed.info[jmed-1].Ey_i[i1-1]) &&
                            (sgg.Eshared.elem[j1-1].j == AniMed.info[jmed-1].Ey_j[i1-1]) &&
                            (sgg.Eshared.elem[j1-1].k == AniMed.info[jmed-1].Ey_k[i1-1]) &&
                            (sgg.Eshared.elem[j1-1].Field == iEy)) {
                            
                            AniMed.info[jmed-1].Ey_Shared[i1-1].times = sgg.Eshared.elem[j1-1].Times;
                            if (sgg.Eshared.elem[j1-1].Times > 1) {
                                AniMed.info[jmed-1].Ey_Shared[i1-1].SharedMed.resize(sgg.Eshared.elem[j1-1].Times);
                            }
                            found = true;
                            break;
                        }
                    }
                }
            }
            
            if (!found) {
                for (int jmed = 1; jmed <= AniMed.NumMed && !found; ++jmed) {
                    for (int i1 = 1; i1 <= AniMed.info[jmed-1].numnodesEz; ++i1) {
                        if ((sgg.Eshared.elem[j1-1].i == AniMed.info[jmed-1].Ez_i[i1-1]) &&
                            (sgg.Eshared.elem[j1-1].j == AniMed.info[jmed-1].Ez_j[i1-1]) &&
                            (sgg.Eshared.elem[j1-1].k == AniMed.info[jmed-1].Ez_k[i1-1]) &&
                            (sgg.Eshared.elem[j1-1].Field == iEz)) {
                            
                            AniMed.info[jmed-1].Ez_Shared[i1-1].times = sgg.Eshared.elem[j1-1].Times;
                            if (sgg.Eshared.elem[j1-1].Times > 1) {
                                AniMed.info[jmed-1].Ez_Shared[i1-1].SharedMed.resize(sgg.Eshared.elem[j1-1].Times);
                            }
                            found = true;
                            break;
                        }
                    }
                }
            }
        }

        for (int j1 = 1; j1 <= sgg.Hshared.conta; ++j1) {
            bool found = false;
            for (int jmed = 1; jmed <= AniMed.NumMed && !found; ++jmed) {
                for (int i1 = 1; i1 <= AniMed.info[jmed-1].numnodesHx; ++i1) {
                    if ((sgg.Hshared.elem[j1-1].i == AniMed.info[jmed-1].Hx_i[i1-1]) &&
                        (sgg.Hshared.elem[j1-1].j == AniMed.info[jmed-1].Hx_j[i1-1]) &&
                        (sgg.Hshared.elem[j1-1].k == AniMed.info[jmed-1].Hx_k[i1-1]) &&
                        (sgg.Hshared.elem[j1-1].Field == iHx)) {
                        
                        AniMed.info[jmed-1].Hx_Shared[i1-1].times = sgg.Hshared.elem[j1-1].Times;
                        if (sgg.Hshared.elem[j1-1].Times > 1) {
                            AniMed.info[jmed-1].Hx_Shared[i1-1].SharedMed.resize(sgg.Hshared.elem[j1-1].Times);
                        }
                        found = true;
                        break;
                    }
                }
            }
            
            if (!found) {
                for (int jmed = 1; jmed <= AniMed.NumMed && !found; ++jmed) {
                    for (int i1 = 1; i1 <= AniMed.info[jmed-1].numnodesHy; ++i1) {
                        if ((sgg.Hshared.elem[j1-1].i == AniMed.info[jmed-1].Hy_i[i1-1]) &&
                            (sgg.Hshared.elem[j1-1].j == AniMed.info[jmed-1].Hy_j[i1-1]) &&
                            (sgg.Hshared.elem[j1-1].k == AniMed.info[jmed-1].Hy_k[i1-1]) &&
                            (sgg.Hshared.elem[j1-1].Field == iHy)) {
                            
                            AniMed.info[jmed-1].Hy_Shared[i1-1].times = sgg.Hshared.elem[j1-1].Times;
                            if (sgg.Hshared.elem[j1-1].Times > 1) {
                                AniMed.info[jmed-1].Hy_Shared[i1-1].SharedMed.resize(sgg.Hshared.elem[j1-1].Times);
                            }
                            found = true;
                            break;
                        }
                    }
                }
            }
            
            if (!found) {
                for (int jmed = 1; jmed <= AniMed.NumMed && !found; ++jmed) {
                    for (int i1 = 1; i1 <= AniMed.info[jmed-1].numnodesHz; ++i1) {
                        if ((sgg.Hshared.elem[j1-1].i == AniMed.info[jmed-1].Hz_i[i1-1]) &&
                            (sgg.Hshared.elem[j1-1].j == AniMed.info[jmed-1].Hz_j[i1-1]) &&
                            (sgg.Hshared.elem[j1-1].k == AniMed.info[jmed-1].Hz_k[i1-1]) &&
                            (sgg.Hshared.elem[j1-1].Field == iHz)) {
                            
                            AniMed.info[jmed-1].Hz_Shared[i1-1].times = sgg.Hshared.elem[j1-1].Times;
                            if (sgg.Hshared.elem[j1-1].Times > 1) {
                                AniMed.info[jmed-1].Hz_Shared[i1-1].SharedMed.resize(sgg.Hshared.elem[j1-1].Times);
                            }
                            found = true;
                            break;
                        }
                    }
                }
            }
        }

        // Store indexes of shared media
        for (int jmed = 1; jmed <= AniMed.NumMed; ++jmed) {
            for (int i1 = 1; i1 <= AniMed.info[jmed-1].numnodesEx; ++i1) {
                int conta = 0;
                for (int j1 = 1; j1 <= sgg.Eshared.conta; ++j1) {
                    if ((sgg.Eshared.elem[j1-1].i == AniMed.info[jmed-1].Ex_i[i1-1]) &&
                        (sgg.Eshared.elem[j1-1].j == AniMed.info[jmed-1].Ex_j[i1-1]) &&
                        (sgg.Eshared.elem[j1-1].k == AniMed.info[jmed-1].Ex_k[i1-1]) &&
                        (sgg.Eshared.elem[j1-1].Field == iEx)) {
                        conta++;
                        AniMed.info[jmed-1].Ex_Shared[i1-1].SharedMed[conta-1] = sgg.Eshared.elem[j1-1].SharedMed;
                    }
                }
            }
        }
        
        for (int jmed = 1; jmed <= AniMed.NumMed; ++jmed) {
            for (int i1 = 1; i1 <= AniMed.info[jmed-1].numnodesEy; ++i1) {
                int conta = 0;
                for (int j1 = 1; j1 <= sgg.Eshared.conta; ++j1) {
                    if ((sgg.Eshared.elem[j1-1].i == AniMed.info[jmed-1].Ey_i[i1-1]) &&
                        (sgg.Eshared.elem[j1-1].j == AniMed.info[jmed-1].Ey_j[i1-1]) &&
                        (sgg.Eshared.elem[j1-1].k == AniMed.info[jmed-1].Ey_k[i1-1]) &&
                        (sgg.Eshared.elem[j1-1].Field == iEy)) {
                        conta++;
                        AniMed.info[jmed-1].Ey_Shared[i1-1].SharedMed[conta-1] = sgg.Eshared.elem[j1-1].SharedMed;
                    }
                }
            }
        }
        
        for (int jmed = 1; jmed <= AniMed.NumMed; ++jmed) {
            for (int i1 = 1; i1 <= AniMed.info[jmed-1].numnodesEz; ++i1) {
                int conta = 0;
                for (int j1 = 1; j1 <= sgg.Eshared.conta; ++j1) {
                    if ((sgg.Eshared.elem[j1-1].i == AniMed.info[jmed-1].Ez_i[i1-1]) &&
                        (sgg.Eshared.elem[j1-1].j == AniMed.info[jmed-1].Ez_j[i1-1]) &&
                        (sgg.Eshared.elem[j1-1].k == AniMed.info[jmed-1].Ez_k[i1-1]) &&
                        (sgg.Eshared.elem[j1-1].Field == iEz)) {
                        conta++;
                        AniMed.info[jmed-1].Ez_Shared[i1-1].SharedMed[conta-1] = sgg.Eshared.elem[j1-1].SharedMed;
                    }
                }
            }
        }
        
        for (int jmed = 1; jmed <= AniMed.NumMed; ++jmed) {
            for (int i1 = 1; i1 <= AniMed.info[jmed-1].numnodesHx; ++i1) {
                int conta = 0;
                for (int j1 = 1; j1 <= sgg.Hshared.conta; ++j1) {
                    if ((sgg.Hshared.elem[j1-1].i == AniMed.info[jmed-1].Hx_i[i1-1]) &&
                        (sgg.Hshared.elem[j1-1].j == AniMed.info[jmed-1].Hx_j[i1-1]) &&
                        (sgg.Hshared.elem[j1-1].k == AniMed.info[jmed-1].Hx_k[i1-1]) &&
                        (sgg.Hshared.elem[j1-1].Field == iHx)) {
                        conta++;
                        AniMed.info[jmed-1].Hx_Shared[i1-1].SharedMed[conta-1] = sgg.Hshared.elem[j1-1].SharedMed;
                    }
                }
            }
        }
        
        for (int jmed = 1; jmed <= AniMed.NumMed; ++jmed) {
            for (int i1 = 1; i1 <= AniMed.info[jmed-1].numnodesHy; ++i1) {
                int conta = 0;
                for (int j1 = 1; j1 <= sgg.Hshared.conta; ++j1) {
                    if ((sgg.Hshared.elem[j1-1].i == AniMed.info[jmed-1].Hy_i[i1-1]) &&
                        (sgg.Hshared.elem[j1-1].j == AniMed.info[jmed-1].Hy_j[i1-1]) &&
                        (sgg.Hshared.elem[j1-1].k == AniMed.info[jmed-1].Hy_k[i1-1]) &&
                        (sgg.Hshared.elem[j1-1].Field == iHy)) {
                        conta++;
                        AniMed.info[jmed-1].Hy_Shared[i1-1].SharedMed[conta-1] = sgg.Hshared.elem[j1-1].SharedMed;
                    }
                }
            }
        }
        
        for (int jmed = 1; jmed <= AniMed.NumMed; ++jmed) {
            for (int i1 = 1; i1 <= AniMed.info[jmed-1].numnodesHz; ++i1) {
                int conta = 0;
                for (int j1 = 1; j1 <= sgg.Hshared.conta; ++j1) {
                    if ((sgg.Hshared.elem[j1-1].i == AniMed.info[jmed-1].Hz_i[i1-1]) &&
                        (sgg.Hshared.elem[j1-1].j == AniMed.info[jmed-1].Hz_j[i1-1]) &&
                        (sgg.Hshared.elem[j1-1].k == AniMed.info[jmed-1].Hz_k[i1-1]) &&
                        (sgg.Hshared.elem[j1-1].Field == iHz)) {
                        conta++;
                        AniMed.info[jmed-1].Hz_Shared[i1-1].SharedMed[conta-1] = sgg.Hshared.elem[j1-1].SharedMed;
                    }
                }
            }
        }

        // Create matrices (calculate averaged coefficients)
        for (int jmed = 1; jmed <= AniMed.NumMed; ++jmed) {
            int idx = AniMed.info[jmed-1].indexmed;
            // dummyAnisProp points to sgg%med(indexmed)%Anisotropic(1)
            // We use a reference to the local copy in sgg
            const auto& dummyAnisProp = sgg.Med[idx-1].Anisotropic[0];
            
            for (int i1 = 1; i1 <= AniMed.info[jmed-1].numnodesEx; ++i1) {
                int conta = AniMed.info[jmed-1].Ex_Shared[i1-1].times;
                if (conta > 1) {
                    // Reset to base / conta
                    for(int r=0; r<3; ++r)
                        for(int c=0; c<3; ++c) {
                            AniMed.info[jmed-1].sigma[r][c] = dummyAnisProp.sigma[r][c] / conta;
                            AniMed.info[jmed-1].sigmam[r][c] = dummyAnisProp.sigmam[r][c] / conta;
                            AniMed.info[jmed-1].mur[r][c] = dummyAnisProp.mur[r][c] / conta;
                            AniMed.info[jmed-1].epr[r][c] = dummyAnisProp.epr[r][c] / conta;
                        }
                    
                    for (int k1 = 1; k1 < conta; ++k1) {
                        int sharedIdx = AniMed.info[jmed-1].Ex_Shared[i1-1].SharedMed[k1-1];
                        // dummyAnisShared points to sgg%med(sharedIdx)%Anisotropic(1)
                        // Note: sharedIdx is an index into sgg.Med. 
                        // We assume SharedMed stores the media index directly.
                        const auto& dummyAnisShared = sgg.Med[sharedIdx-1].Anisotropic[0];
                        
                        for(int r=0; r<3; ++r)
                            for(int c=0; c<3; ++c) {
                                AniMed.info[jmed-1].sigma[r][c] += dummyAnisShared.sigma[r][c] / conta;
                                AniMed.info[jmed-1].sigmam[r][c] += dummyAnisShared.sigmam[r][c] / conta;
                                AniMed.info[jmed-1].mur[r][c] += dummyAnisShared.mur[r][c] / conta;
                                AniMed.info[jmed-1].epr[r][c] += dummyAnisShared.epr[r][c] / conta;
                            }
                    }
                }
            }

            for (int i1 = 1; i1 <= AniMed.info[jmed-1].numnodesEy; ++i1) {
                int conta = AniMed.info[jmed-1].Ey_Shared[i1-1].times;
                if (conta > 1) {
                    for(int r=0; r<3; ++r)
                        for(int c=0; c<3; ++c) {
                            AniMed.info[jmed-1].sigma[r][c] = dummyAnisProp.sigma[r][c] / conta;
                            AniMed.info[jmed-1].sigmam[r][c] = dummyAnisProp.sigmam[r][c] / conta;
                            AniMed.info[jmed-1].mur[r][c] = dummyAnisProp.mur[r][c] / conta;
                            AniMed.info[jmed-1].epr[r][c] = dummyAnisProp.epr[r][c] / conta;
                        }
                    
                    for (int k1 = 1; k1 < conta; ++k1) {
                        int sharedIdx = AniMed.info[jmed-1].Ey_Shared[i1-1].SharedMed[k1-1];
                        const auto& dummyAnisShared = sgg.Med[sharedIdx-1].Anisotropic[0];
                        
                        for(int r=0; r<3; ++r)
                            for(int c=0; c<3; ++c) {
                                AniMed.info[jmed-1].sigma[r][c] += dummyAnisShared.sigma[r][c] / conta;
                                AniMed.info[jmed-1].sigmam[r][c] += dummyAnisShared.sigmam[r][c] / conta;
                                AniMed.info[jmed-1].mur[r][c] += dummyAnisShared.mur[r][c] / conta;
                                AniMed.info[jmed-1].epr[r][c] += dummyAnisShared.epr[r][c] / conta;
                            }
                    }
                }
            }

            for (int i1 = 1; i1 <= AniMed.info[jmed-1].numnodesEz; ++i1) {
                int conta = AniMed.info[jmed-1].Ez_Shared[i1-1].times;
                if (conta > 1) {
                    for(int r=0; r<3; ++r)
                        for(int c=0; c<3; ++c) {
                            AniMed.info[jmed-1].sigma[r][c] = dummyAnisProp.sigma[r][c] / conta;
                            AniMed.info[jmed-1].sigmam[r][c] = dummyAnisProp.sigmam[r][c] / conta;
                            AniMed.info[jmed-1].mur[r][c] = dummyAnisProp.mur[r][c] / conta;
                            AniMed.info[jmed-1].epr[r][c] = dummyAnisProp.epr[r][c] / conta;
                        }
                    
                    for (int k1 = 1; k1 < conta; ++k1) {
                        int sharedIdx = AniMed.info[jmed-1].Ez_Shared[i1-1].SharedMed[k1-1];
                        const auto& dummyAnisShared = sgg.Med[sharedIdx-1].Anisotropic[0];
                        
                        for(int r=0; r<3; ++r)
                            for(int c=0; c<3; ++c) {
                                AniMed.info[jmed-1].sigma[r][c] += dummyAnisShared.sigma[r][c] / conta;
                                AniMed.info[jmed-1].sigmam[r][c] += dummyAnisShared.sigmam[r][c] / conta;
                                AniMed.info[jmed-1].mur[r][c] += dummyAnisShared.mur[r][c] / conta;
                                AniMed.info[jmed-1].epr[r][c] += dummyAnisShared.epr[r][c] / conta;
                            }
                    }
                }
            }
            
            // Note: The original code snippet cuts off here. 
            // Assuming similar logic for Hx, Hy, Hz follows.
            // Since the input code was truncated, I will stop here to match the input.
        }
    }

    // Stubs for other required functions to satisfy "Convert ALL subroutines"
    // The input code did not provide implementations for these, but they were declared public.
    void AdvanceAnisotropicE() {
        // Implementation not provided in source
    }

    void AdvanceAnisotropich() {
        // Implementation not provided in source
    }

    void DestroyAnisotropic() {
        AniMed.info.clear();
        AniMed.NumMed = 0;
    }

    void calc_anisotropicconstants() {
        // Implementation not provided in source
    }

} // namespace Anisotropic_m

AniMed.info[jmed].mur += dummyAnisShared.mur / conta;
            AniMed.info[jmed].epr += dummyAnisShared.epr / conta;
         }
      }
   }

   for (int i1 = 1; i1 <= AniMed.info[jmed].NumNodesHx; ++i1) {
      conta = AniMed.info[jmed].Hx_Shared[i1].times;
      if (conta > 1) {
         AniMed.info[jmed].sigma = dummyAnisProp.sigma / conta;
         AniMed.info[jmed].sigmam = dummyAnisProp.sigmam / conta;
         AniMed.info[jmed].mur = dummyAnisProp.mur / conta;
         AniMed.info[jmed].epr = dummyAnisProp.epr / conta;
         for (int k1 = 1; k1 <= conta - 1; ++k1) {
            dummyAnisShared = sgg.med[AniMed.info[jmed].Hx_Shared[i1].SharedMed[k1]].Anisotropic[1];
            AniMed.info[jmed].sigma += dummyAnisShared.sigma / conta;
            AniMed.info[jmed].sigmam += dummyAnisShared.sigmam / conta;
            AniMed.info[jmed].mur += dummyAnisShared.mur / conta;
            AniMed.info[jmed].epr += dummyAnisShared.epr / conta;
         }
      }
   }

   for (int i1 = 1; i1 <= AniMed.info[jmed].NumNodesHy; ++i1) {
      conta = AniMed.info[jmed].Hy_Shared[i1].times;
      if (conta > 1) {
         AniMed.info[jmed].sigma = dummyAnisProp.sigma / conta;
         AniMed.info[jmed].sigmam = dummyAnisProp.sigmam / conta;
         AniMed.info[jmed].mur = dummyAnisProp.mur / conta;
         AniMed.info[jmed].epr = dummyAnisProp.epr / conta;
         for (int k1 = 1; k1 <= conta - 1; ++k1) {
            dummyAnisShared = sgg.med[AniMed.info[jmed].Hy_Shared[i1].SharedMed[k1]].Anisotropic[1];
            AniMed.info[jmed].sigma += dummyAnisShared.sigma / conta;
            AniMed.info[jmed].sigmam += dummyAnisShared.sigmam / conta;
            AniMed.info[jmed].mur += dummyAnisShared.mur / conta;
            AniMed.info[jmed].epr += dummyAnisShared.epr / conta;
         }
      }
   }

   for (int i1 = 1; i1 <= AniMed.info[jmed].NumNodesHz; ++i1) {
      conta = AniMed.info[jmed].Hz_Shared[i1].times;
      if (conta > 1) {
         AniMed.info[jmed].sigma = dummyAnisProp.sigma / conta;
         AniMed.info[jmed].sigmam = dummyAnisProp.sigmam / conta;
         AniMed.info[jmed].mur = dummyAnisProp.mur / conta;
         AniMed.info[jmed].epr = dummyAnisProp.epr / conta;
         for (int k1 = 1; k1 <= conta - 1; ++k1) {
            dummyAnisShared = sgg.med[AniMed.info[jmed].Hz_Shared[i1].SharedMed[k1]].Anisotropic[1];
            AniMed.info[jmed].sigma += dummyAnisShared.sigma / conta;
            AniMed.info[jmed].sigmam += dummyAnisShared.sigmam / conta;
            AniMed.info[jmed].mur += dummyAnisShared.mur / conta;
            AniMed.info[jmed].epr += dummyAnisShared.epr / conta;
         }
      }
   }
   }
   }

   calc_anisotropicconstants(sgg, eps0, mu0);

   return;

} // End of InitAnisotropic

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// subroutine to advance the E field in the Anisotropic (no need to advance the magnetic field)
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

void AdvanceAnisotropicE(SGGFDTDINFO& sggAlloc,
                         std::vector<std::vector<std::vector<double>>>& Ex,
                         std::vector<std::vector<std::vector<double>>>& Ey,
                         std::vector<std::vector<std::vector<double>>>& Ez,
                         std::vector<std::vector<std::vector<double>>>& Hx,
                         std::vector<std::vector<std::vector<double>>>& Hy,
                         std::vector<std::vector<std::vector<double>>>& Hz,
                         const std::vector<double>& Idxe,
                         const std::vector<double>& Idye,
                         const std::vector<double>& Idze,
                         const std::vector<double>& Idxh,
                         const std::vector<double>& Idyh,
                         const std::vector<double>& Idzh) {

   Coeff* coeff = nullptr;
   Anisotropicinfo* dummy = nullptr;
   int jmed, i1, i, j, k;

   for (jmed = 1; jmed <= AniMed.NumMed; ++jmed) {
      dummy = &AniMed.info[jmed];
      
      for (i1 = 1; i1 <= dummy->NumNodesEx; ++i1) {
         if (dummy->Ex_Shared[i1].times == 1) {
            coeff = &AniMed.info[jmed].coeff;
         } else {
            coeff = &AniMed.info[jmed].Ex_Shared[i1].Coeff;
         }

         i = dummy->Ex_i[i1];
         j = dummy->Ex_j[i1];
         k = dummy->Ex_k[i1];
         
         if (dummy->IsOnlyThinSlot) {
            dummy->Ex_value[i1] = ex(i, j, k) - coeff->ehxx * (hz(i, j - 1, k) - hz(i, j, k)) * Idyh[j] +
                                  coeff->ehxx * (hy(i, j, k - 1) - hy(i, j, k)) * Idzh[k];
         } else {
            dummy->Ex_value[i1] = 
               (4.0 * coeff->eexx * ex(i, j, k) + coeff->eexy * ey(i, j - 1, k) + coeff->eexy * ey(i, j, k) +
                coeff->eexy * ey(i + 1, j - 1, k) + coeff->eexy * ey(i + 1, j, k) + coeff->eexz * ez(i, j, k - 1) +
                coeff->eexz * ez(i, j, k) + coeff->eexz * ez(i + 1, j, k - 1) + coeff->eexz * ez(i + 1, j, k)) / 4.0 -
               ((coeff->ehxz * hy(i - 1, j, k - 1) + coeff->ehxz * hy(i - 1, j, k) - coeff->ehxz * hy(i + 1, j, k - 1) -
                 coeff->ehxz * hy(i + 1, j, k) - coeff->ehxy * hz(i - 1, j - 1, k) - coeff->ehxy * hz(i - 1, j, k) +
                 coeff->ehxy * hz(i + 1, j - 1, k) + coeff->ehxy * hz(i + 1, j, k)) * Idxe[i]) / 4.0 +
               ((coeff->ehxz * hx(i, j - 1, k - 1) + coeff->ehxz * hx(i, j - 1, k) - coeff->ehxz * hx(i, j, k - 1) -
                 coeff->ehxz * hx(i, j, k) + coeff->ehxz * hx(i + 1, j - 1, k - 1) +
                 coeff->ehxz * hx(i + 1, j - 1, k) - coeff->ehxz * hx(i + 1, j, k - 1) - coeff->ehxz * hx(i + 1, j, k) -
                 4.0 * coeff->ehxx * hz(i, j - 1, k) + 4.0 * coeff->ehxx * hz(i, j, k)) * Idyh[j]) / 4.0 -
               ((coeff->ehxy * hx(i, j - 1, k - 1) - coeff->ehxy * hx(i, j - 1, k) + coeff->ehxy * hx(i, j, k - 1) -
                 coeff->ehxy * hx(i, j, k) + coeff->ehxy * hx(i + 1, j - 1, k - 1) -
                 coeff->ehxy * hx(i + 1, j - 1, k) + coeff->ehxy * hx(i + 1, j, k - 1) - coeff->ehxy * hx(i + 1, j, k) -
                 4.0 * coeff->ehxx * hy(i, j, k - 1) + 4.0 * coeff->ehxx * hy(i, j, k)) * Idzh[k]) / 4.0;
         }
      }

      for (i1 = 1; i1 <= dummy->NumNodesEy; ++i1) {
         if (dummy->Ey_Shared[i1].times == 1) {
            coeff = &AniMed.info[jmed].coeff;
         } else {
            coeff = &AniMed.info[jmed].Ey_Shared[i1].Coeff;
         }
         
         i = dummy->Ey_i[i1];
         j = dummy->Ey_j[i1];
         k = dummy->Ey_k[i1];
         
         if (dummy->IsOnlyThinSlot) {
            dummy->Ey_value[i1] = ey(i, j, k) + coeff->ehyy * (hz(i - 1, j, k) - hz(i, j, k)) * Idxh[i] -
                                  coeff->ehyy * (hx(i, j, k - 1) - hx(i, j, k)) * Idzh[k];
         } else {
            dummy->Ey_value[i1] = 
               (coeff->eeyx * ex(i - 1, j, k) + coeff->eeyx * ex(i - 1, j + 1, k) + coeff->eeyx * ex(i, j, k) +
                coeff->eeyx * ex(i, j + 1, k) + 4.0 * coeff->eeyy * ey(i, j, k) + coeff->eeyz * ez(i, j, k - 1) +
                coeff->eeyz * ez(i, j, k) + coeff->eeyz * ez(i, j + 1, k - 1) + coeff->eeyz * ez(i, j + 1, k)) / 4.0 -
               ((coeff->ehyz * hy(i - 1, j, k - 1) + coeff->ehyz * hy(i - 1, j, k) +
                 coeff->ehyz * hy(i - 1, j + 1, k - 1) + coeff->ehyz * hy(i - 1, j + 1, k) -
                 coeff->ehyz * hy(i, j, k - 1) - coeff->ehyz * hy(i, j, k) - coeff->ehyz * hy(i, j + 1, k - 1) -
                 coeff->ehyz * hy(i, j + 1, k) - 4.0 * coeff->ehyy * hz(i - 1, j, k) + 4.0 * coeff->ehyy * hz(i, j, k)) * Idxh[i]) / 4.0 +
               ((coeff->ehyz * hx(i, j - 1, k - 1) + coeff->ehyz * hx(i, j - 1, k) -
                 coeff->ehyz * hx(i, j + 1, k - 1) - coeff->ehyz * hx(i, j + 1, k) -
                 coeff->ehyx * hz(i - 1, j - 1, k) + coeff->ehyx * hz(i - 1, j + 1, k) -
                 coeff->ehyx * hz(i, j - 1, k) + coeff->ehyx * hz(i, j + 1, k)) * Idye[j]) / 4.0 -
               ((4.0 * coeff->ehyy * hx(i, j, k - 1) - 4.0 * coeff->ehyy * hx(i, j, k) - coeff->ehyx * hy(i - 1, j, k - 1) +
                 coeff->ehyx * hy(i - 1, j, k) - coeff->ehyx * hy(i - 1, j + 1, k - 1) +
                 coeff->ehyx * hy(i - 1, j + 1, k) - coeff->ehyx * hy(i, j, k - 1) + coeff->ehyx * hy(i, j, k) -
                 coeff->ehyx * hy(i, j + 1, k - 1) + coeff->ehyx * hy(i, j + 1, k)) * Idzh[k]) / 4.0;
         }
      }

      for (i1 = 1; i1 <= dummy->NumNodesEz; ++i1) {
         if (dummy->Ez_Shared[i1].times == 1) {
            coeff = &AniMed.info[jmed].coeff;
         } else {
            coeff = &AniMed.info[jmed].Ez_Shared[i1].Coeff;
         }
         
         i = dummy->Ez_i[i1];
         j = dummy->Ez_j[i1];
         k = dummy->Ez_k[i1];
         
         if (dummy->IsOnlyThinSlot) {
            dummy->Ez_value[i1] = ez(i, j, k) - coeff->ehzz * (hy(i - 1, j, k) - hy(i, j, k)) * Idxh[i] +
                                  coeff->ehzz * (hx(i, j - 1, k) - hx(i, j, k)) * Idyh[j];
         } else {
            dummy->Ez_value[i1] = 
               (coeff->eezx * ex(i - 1, j, k) + coeff->eezx * ex(i - 1, j, k + 1) + coeff->eezx * ex(i, j, k) +
                coeff->eezx * ex(i, j, k + 1) + coeff->eezy * ey(i, j - 1, k) + coeff->eezy * ey(i, j - 1, k + 1) +
                coeff->eezy * ey(i, j, k) + coeff->eezy * ey(i, j, k + 1) + 4.0 * coeff->eezz * ez(i, j, k)) / 4.0 -
               ((4.0 * coeff->ehzz * hy(i - 1, j, k) - 4.0 * coeff->ehzz * hy(i, j, k) - coeff->ehzy * hz(i - 1, j - 1, k) -
                 coeff->ehzy * hz(i - 1, j - 1, k + 1) - coeff->ehzy * hz(i - 1, j, k) -
                 coeff->ehzy * hz(i - 1, j, k + 1) + coeff->ehzy * hz(i, j - 1, k) +
                 coeff->ehzy * hz(i, j - 1, k + 1) + coeff->ehzy * hz(i, j, k) + coeff->ehzy * hz(i, j, k + 1)) * Idxh[i]) / 4.0 +
               ((4.0 * coeff->ehzz * hx(i, j - 1, k) - 4.0 * coeff->ehzz * hx(i, j, k) - coeff->ehzx * hz(i - 1, j - 1, k) -
                 coeff->ehzx * hz(i - 1, j - 1, k + 1) + coeff->ehzx * hz(i - 1, j, k) +
                 coeff->ehzx * hz(i - 1, j, k + 1) - coeff->ehzx * hz(i, j - 1, k) -
                 coeff->ehzx * hz(i, j - 1, k + 1) + coeff->ehzx * hz(i, j, k) + coeff->ehzx * hz(i, j, k + 1)) * Idyh[j]) / 4.0 -
               ((coeff->ehzy * hx(i, j - 1, k - 1) - coeff->ehzy * hx(i, j - 1, k + 1) +
                 coeff->ehzy * hx(i, j, k - 1) - coeff->ehzy * hx(i, j, k + 1) - coeff->ehzx * hy(i - 1, j, k - 1) +
                 coeff->ehzx * hy(i - 1, j, k + 1) - coeff->ehzx * hy(i, j, k - 1) + coeff->ehzx * hy(i, j, k + 1)) * Idze[k]) / 4.0;
         }
      }
   }

   // store it
   for (jmed = 1; jmed <= AniMed.NumMed; ++jmed) {
      dummy = &AniMed.info[jmed];
      for (i1 = 1; i1 <= dummy->NumNodesEx; ++i1) {
         i = dummy->Ex_i[i1];
         j = dummy->Ex_j[i1];
         k = dummy->Ex_k[i1];
         ex(i, j, k) = dummy->Ex_value[i1];
      }
      for (i1 = 1; i1 <= dummy->NumNodesEy; ++i1) {
         i = dummy->Ey_i[i1];
         j = dummy->Ey_j[i1];
         k = dummy->Ey_k[i1];
         ey(i, j, k) = dummy->Ey_value[i1];
      }
      for (i1 = 1; i1 <= dummy->NumNodesEz; ++i1) {
         i = dummy->Ez_i[i1];
         j = dummy->Ez_j[i1];
         k = dummy->Ez_k[i1];
         ez(i, j, k) = dummy->Ez_value[i1];
      }
   }

   return;
} // End of AdvanceAnisotropicE

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// subroutine to advance the H field in the Anisotropic (no need to advance the electric field)
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

void AdvanceAnisotropicH(SGGFDTDINFO& sggAlloc,
                         std::vector<std::vector<std::vector<double>>>& Ex,
                         std::vector<std::vector<std::vector<double>>>& Ey,
                         std::vector<std::vector<std::vector<double>>>& Ez,
                         std::vector<std::vector<std::vector<double>>>& Hx,
                         std::vector<std::vector<std::vector<double>>>& Hy,
                         std::vector<std::vector<std::vector<double>>>& Hz,
                         const std::vector<double>& Idxe,
                         const std::vector<double>& Idye,
                         const std::vector<double>& Idze,
                         const std::vector<double>& Idxh,
                         const std::vector<double>& Idyh,
                         const std::vector<double>& Idzh) {

   Anisotropicinfo* dummy = nullptr;
   Coeff* coeff = nullptr;
   int jmed, i1, i, j, k;

   for (jmed = 1; jmed <= AniMed.NumMed; ++jmed) {
      dummy = &AniMed.info[jmed];
      
      for (i1 = 1; i1 <= dummy->NumNodesHx; ++i1) {
         if (dummy->Hx_Shared[i1].times == 1) {
            coeff = &AniMed.info[jmed].coeff;
         } else {
            coeff = &AniMed.info[jmed].Hx_Shared[i1].Coeff;
         }
         
         i = dummy->Hx_i[i1];
         j = dummy->Hx_j[i1];
         k = dummy->Hx_k[i1];
         
         if (dummy->IsOnlyThinSlot) {
            dummy->Hx_value[i1] = hx(i, j, k) + coeff->hexx * (ez(i, j, k) - ez(i, j + 1, k)) * Idye[j] -
                                  coeff->hexx * (ey(i, j, k) - ey(i, j, k + 1)) * Idze[k];
         } else {
            dummy->Hx_value[i1] = 
               (4.0 * coeff->hhxx * hx(i, j, k) + coeff->hhxy * hy(i - 1, j, k) + coeff->hhxy * hy(i - 1, j + 1, k) +
                coeff->hhxy * hy(i, j, k) + coeff->hhxy * hy(i, j + 1, k) + coeff->hhxz * hz(i - 1, j, k) +
                coeff->hhxz * hz(i - 1, j, k + 1) + coeff->hhxz * hz(i, j, k) + coeff->hhxz * hz(i, j, k + 1)) / 4.0 +
               ((coeff->hexz * ey(i - 1, j, k) + coeff->hexz * ey(i - 1, j, k + 1) - coeff->hexz * ey(i + 1, j, k) -
                 coeff->hexz * ey(i + 1, j, k + 1) - coeff->hexy * ez(i - 1, j, k) - coeff->hexy * ez(i - 1, j + 1, k) +
                 coeff->hexy * ez(i + 1, j, k) + coeff->hexy * ez(i + 1, j + 1, k)) * Idxh[i]) / 4.0 -
               ((coeff->hexz * ex(i - 1, j, k) + coeff->hexz * ex(i - 1, j, k + 1) - coeff->hexz * ex(i - 1, j + 1, k) -
                 coeff->hexz * ex(i - 1, j + 1, k + 1) + coeff->hexz * ex(i, j, k) + coeff->hexz * ex(i, j, k + 1) -
                 coeff->hexz * ex(i, j + 1, k) - coeff->hexz * ex(i, j + 1, k + 1) - 4.0 * coeff->hexx * ez(i, j, k) +
                 4.0 * coeff->hexx * ez(i, j + 1, k)) * Idye[j]) / 4.0 +
               ((coeff->hexy * ex(i - 1, j, k) - coeff->hexy * ex(i - 1, j, k + 1) + coeff->hexy * ex(i - 1, j + 1, k) -
                 coeff->hexy * ex(i - 1, j + 1, k + 1) + coeff->hexy * ex(i, j, k) - coeff->hexy * ex(i, j, k + 1) +
                 coeff->hexy * ex(i, j + 1, k) - coeff->hexy * ex(i, j + 1, k + 1) - 4.0 * coeff->hexx * ey(i, j, k) +
                 4.0 * coeff->hexx * ey(i, j, k + 1)) * Idze[k]) / 4.0;
         }
      }

      for (i1 = 1; i1 <= dummy->NumNodesHy; ++i1) {
         if (dummy->Hy_Shared[i1].times == 1) {
            coeff = &AniMed.info[jmed].coeff;
         } else {
            coeff = &AniMed.info[jmed].Hy_Shared[i1].Coeff;
         }
         
         i = dummy->Hy_i[i1];
         j = dummy->Hy_j[i1];
         k = dummy->Hy_k[i1];
         
         if (dummy->IsOnlyThinSlot) {
            dummy->Hy_value[i1] = hy(i, j, k) - coeff->heyy * (ez(i, j, k) - ez(i + 1, j, k)) * Idxe[i] +
                                  coeff->heyy * (ex(i, j, k) - ex(i, j, k + 1)) * Idze[k];
         } else {
            dummy->Hy_value[i1] = 
               (coeff->hhyx * hx(i, j - 1, k) + coeff->hhyx * hx(i, j, k) + coeff->hhyx * hx(i + 1, j - 1, k) +
                coeff->hhyx * hx(i + 1, j, k) + 4.0 * coeff->hhyy * hy(i, j, k) + coeff->hhyz * hz(i, j - 1, k) +
                coeff->hhyz * hz(i, j - 1, k + 1) + coeff->hhyz * hz(i, j, k) + coeff->hhyz * hz(i, j, k + 1)) / 4.0 +
               ((coeff->heyz * ey(i, j - 1, k) + coeff->heyz * ey(i, j - 1, k + 1) + coeff->heyz * ey(i, j, k) +
                 coeff->heyz * ey(i, j, k + 1) - coeff->heyz * ey(i + 1, j - 1, k) -
                 coeff->heyz * ey(i + 1, j - 1, k + 1) - coeff->heyz * ey(i + 1, j, k) -
                 coeff->heyz * ey(i + 1, j, k + 1) - 4.0 * coeff->heyy * ez(i, j, k) + 4.0 * coeff->heyy * ez(i + 1, j, k)) * Idxe[i]) / 4.0 -
               ((coeff->heyz * ex(i, j - 1, k) + coeff->heyz * ex(i, j - 1, k + 1) - coeff->heyz * ex(i, j + 1, k) -
                 coeff->heyz * ex(i, j + 1, k + 1) - coeff->heyx * ez(i, j - 1, k) +
                 coeff->heyx * ez(i, j + 1, k) - coeff->heyx * ez(i + 1, j - 1, k) + coeff->heyx * ez(i + 1, j + 1, k)) * Idyh[j]) / 4.0 +
               ((4.0 * coeff->heyy * ex(i, j, k) - 4.0 * coeff->heyy * ex(i, j, k + 1) - coeff->heyx * ey(i, j - 1, k) +
                 coeff->heyx * ey(i, j - 1, k + 1) - coeff->heyx * ey(i, j, k) + coeff->heyx * ey(i, j, k + 1) -
                 coeff->heyx * ey(i + 1, j - 1, k) + coeff->heyx * ey(i + 1, j - 1, k + 1) -
                 coeff->heyx * ey(i + 1, j, k) + coeff->heyx * ey(i + 1, j, k + 1)) * Idze[k]) / 4.0;
         }
      }

      for (i1 = 1; i1 <= dummy->NumNodesHz; ++i1) {
         if (dummy->Hz_Shared[i1].times == 1) {
            coeff = &AniMed.info[jmed].coeff;
         } else {
            coeff = &AniMed.info[jmed].Hz_Shared[i1].Coeff;
         }
         
         i = dummy->Hz_i[i1];
         j = dummy->Hz_j[i1];
         k = dummy->Hz_k[i1];
         
         if (dummy->IsOnlyThinSlot) {
            dummy->Hz_value[i1] = hz(i, j, k) + coeff->hezz * (ey(i, j, k) - ey(i + 1, j, k)) * Idxe[i] -
                                  coeff->hezz * (ex(i, j, k) - ex(i, j + 1, k)) * Idye[j];
         } else {
            dummy->Hz_value[i1] = 
               (coeff->hhzx * hx(i, j, k - 1) + coeff->hhzx * hx(i, j, k) + coeff->hhzx * hx(i + 1, j, k - 1) +
                coeff->hhzx * hx(i + 1, j, k) + coeff->hhzy * hy(i, j, k - 1) + coeff->hhzy * hy(i, j, k) +
                coeff->hhzy * hy(i, j + 1, k - 1) + coeff->hhzy * hy(i, j + 1, k) + 4.0 * coeff->hhzz * hz(i, j, k)) / 4.0 +
               ((4.0 * coeff->hezz * ey(i, j, k) - 4.0 * coeff->hezz * ey(i + 1, j, k) - coeff->hezy * ez(i, j, k - 1) -
                 coeff->hezy * ez(i, j, k) - coeff->hezy * ez(i, j + 1, k - 1) - coeff->hezy * ez(i, j + 1, k) +
                 coeff->hezy * ez(i + 1, j, k - 1) + coeff->hezy * ez(i + 1, j, k) +
                 coeff->hezy * ez(i + 1, j + 1, k - 1) + coeff->hezy * ez(i + 1, j + 1, k)) * Idxe[i]) / 4.0 -
               ((4.0 * coeff->hezz * ex(i, j, k) - 4.0 * coeff->hezz * ex(i, j + 1, k) - coeff->hezx * ez(i, j, k - 1) -
                 coeff->hezx * ez(i, j, k) + coeff->hezx * ez(i, j + 1, k - 1) + coeff->hezx * ez(i, j + 1, k) -
                 coeff->hezx * ez(i + 1, j, k - 1) - coeff->hezx * ez(i + 1, j, k) +
                 coeff->hezx * ez(i + 1, j + 1, k - 1) + coeff->hezx * ez(i + 1, j + 1, k)) * Idye[j]) / 4.0 +
               ((coeff->hezy * ex(i, j, k - 1) - coeff->hezy * ex(i, j, k + 1) + coeff->hezy * ex(i, j + 1, k - 1) -
                 coeff->hezy * ex(i, j + 1, k + 1) - coeff->hezx * ey(i, j, k - 1) + coeff->hezx * ey(i, j, k + 1) -
                 coeff->hezx * ey(i + 1, j, k - 1) + coeff->hezx * ey(i + 1, j, k + 1)) * Idzh[k]) / 4.0;
         }
      }
   }

   // store it
   for (jmed = 1; jmed <= AniMed.NumMed; ++jmed) {
      dummy = &AniMed.info[jmed];
      for (i1 = 1; i1 <= dummy->NumNodesHx; ++i1) {
         i = dummy->Hx_i[i1];
         j = dummy->Hx_j[i1];
         k = dummy->Hx_k[i1];
         Hx(i, j, k) = dummy->Hx_value[i1];
      }
      for (i1 = 1; i1 <= dummy->NumNodesHy; ++i1) {
         i = dummy->Hy_i[i1];
         j = dummy->Hy_j[i1];
         k = dummy->Hy_k[i1];
         Hy(i, j, k) = dummy->Hy_value[i1];
      }
      for (i1 = 1; i1 <= dummy->NumNodesHz; ++i1) {
         i = dummy->Hz_i[i1];
         j = dummy->Hz_j[i1];
         k = dummy->Hz_k[i1];
         Hz(i, j, k) = dummy->Hz_value[i1];
      }
   }

   return;
} // End of AdvanceAnisotropicH

void DestroyAnisotropic(SGGFDTDINFO& sgg) {
   int jmed, i;
   
   // free up memory
   for (jmed = 1; jmed <= AniMed.NumMed; ++jmed) {
      AniMed.info[jmed].Ex_i.clear();
      AniMed.info[jmed].Ex_j.clear();
      AniMed.info[jmed].Ex_k.clear();
      AniMed.info[jmed].Ey_i.clear();
      AniMed.info[jmed].Ey_j.clear();
      AniMed.info[jmed].Ey_k.clear();
      AniMed.info[jmed].Ez_i.clear();
      AniMed.info[jmed].Ez_j.clear();
      AniMed.info[jmed].Ez_k.clear();
      AniMed.info[jmed].Hx_i.clear();
      AniMed.info[jmed].Hx_j.clear();
      AniMed.info[jmed].Hx_k.clear();
      AniMed.info[jmed].Hy_i.clear();
      AniMed.info[jmed].Hy_j.clear();
      AniMed.info[jmed].Hy_k.clear();
      AniMed.info[jmed].Hz_i.clear();
      AniMed.info[jmed].Hz_j.clear();
      AniMed.info[jmed].Hz_k.clear();
      AniMed.info[jmed].Ex_value.clear();
      AniMed.info[jmed].Ex_Shared.clear();
      AniMed.info[jmed].Ey_value.clear();
      AniMed.info[jmed].Ey_Shared.clear();
      AniMed.info[jmed].Ez_value.clear();
      AniMed.info[jmed].Ez_Shared.clear();
      AniMed.info[jmed].Hx_value.clear();
      AniMed.info[jmed].Hx_Shared.clear();
      AniMed.info[jmed].Hy_value.clear();
      AniMed.info[jmed].Hy_Shared.clear();
      AniMed.info[jmed].Hz_value.clear();
      AniMed.info[jmed].Hz_Shared.clear();
   }
   
   for (i = 1; i <= sgg.NumMedia; ++i) {
      if (sgg.Med[i].Is.Anisotropic) {
         sgg.Med[i].Anisotropic.clear();
      }
   }
   AniMed.info.clear();
}

// function para publicar el valor de Med
AnisotropicMed* GetMed() {
   return &AniMed;
}

// found by the mathematica notebook
void calc_anisotropicconstants(SGGFDTDINFO& sgg, double& eps00, double& mu00) {
   double sigma[3][3], epr[3][3], mur[3][3], sigmaM[3][3];
   int jmed;
   
   eps0 = eps00;
   mu0 = mu00; // chapuz para convertir la variables de paso en globales
   zvac = sqrt(mu0 / eps0);
   cluz = 1.0 / sqrt(eps0 * mu0);

   // calculate the coefficients
   for (jmed = 1; jmed <= AniMed.NumMed; ++jmed) {
      // Note: In C++, we need to copy the 3x3 arrays or pass references if the struct holds them.
      // Assuming AniMed.info[jmed].sigma is a 3x3 array or similar structure accessible.
      // For translation purposes, assuming direct access or copy is needed.
      // Since Fortran arrays are contiguous, we might need to handle this carefully.
      // Here we assume simple assignment/copy for the 3x3 blocks.
      
      // Copy sigma
      for(int r=0; r<3; ++r)
         for(int c=0; c<3; ++c)
            sigma[r][c] = AniMed.info[jmed].sigma[r][c]; // Assuming 0-indexed internal storage or mapped access
         
      // Copy sigmam
      for(int r=0; r<3; ++r)
         for(int c=0; c<3; ++c)
            sigmaM[r][c] = AniMed.info[jmed].sigmam[r][c];
         
      // Copy mur
      for(int r=0; r<3; ++r)
         for(int c=0; c<3; ++c)
            mur[r][c] = AniMed.info[jmed].mur[r][c];
         
      // Copy epr
      for(int r=0; r<3; ++r)
         for(int c=0; c<3; ++c)
            epr[r][c] = AniMed.info[jmed].epr[r][c];

      // copied from hirf_anis1.nb
      CalculateCoeff(epr, mur, sigma, sigmaM, sgg.dt, AniMed.info[jmed].coeff);
   }
}

return;
}

void CalculateCoeff(coeff_t& coeff, const std::array<std::array<double, 3>, 3>& epr, const std::array<std::array<double, 3>, 3>& mur, const std::array<std::array<double, 3>, 3>& sigma, const std::array<std::array<double, 3>, 3>& sigmaM, double dt) {
    // Note: mur is not used in the calculation logic provided in the Fortran snippet, 
    // but is passed as an argument. We keep it for signature compatibility.
    // sigmaM is also passed but not used in the visible assignments.

    // Helper lambda for common denominator terms to reduce repetition if needed, 
    // but direct translation is safer for exact numerical equivalence.
    
    // Common terms appearing in denominators and numerators
    // Let's define some intermediate variables to make the C++ code readable and closer to the Fortran structure
    
    // Denominator base term (appears in many places)
    // D_base = -((2 * eps0 * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) + (2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(2,3) + dt * sigma(2,3))
    // However, the full denominator is a sum of products. Let's compute the full denominator first.
    
    // Denominator components
    // Term 1: -((2 * eps0 * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) + (2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(2,3) + dt * sigma(2,3))
    double term1_denom = -((2.0 * eps0 * epr[0][2] + dt * sigma[0][2]) * (2.0 * eps0 * epr[1][1] + dt * sigma[1][1])) + (2.0 * eps0 * epr[0][1] + dt * sigma[0][1]) * (2.0 * eps0 * epr[1][2] + dt * sigma[1][2]);
    
    // Term 2: -((2 * eps0 * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) + (2 * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr(2,3) + dt * sigma(2,3))
    double term2_denom = -((2.0 * eps0 * epr[0][2] + dt * sigma[0][2]) * (2.0 * eps0 * epr[1][0] + dt * sigma[1][0])) + (2.0 * eps0 * epr[0][0] + dt * sigma[0][0]) * (2.0 * eps0 * epr[1][2] + dt * sigma[1][2]);
    
    // Term 3: -((2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) + (2 * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))
    double term3_denom = -((2.0 * eps0 * epr[0][1] + dt * sigma[0][1]) * (2.0 * eps0 * epr[1][0] + dt * sigma[1][0])) + (2.0 * eps0 * epr[0][0] + dt * sigma[0][0]) * (2.0 * eps0 * epr[1][1] + dt * sigma[1][1]);
    
    // Denominator parts involving sigma/dt
    double denom_part1 = (eps0 * epr[2][0]) / dt + sigma[2][0] / 2.0;
    double denom_part2 = (eps0 * epr[2][1]) / dt + sigma[2][1] / 2.0;
    double denom_part3 = (eps0 * epr[2][2]) / dt + sigma[2][2] / 2.0;
    
    double denominator = term1_denom * denom_part1 - term2_denom * denom_part2 + term3_denom * denom_part3;

    // --- coeff%eexx ---
    // Numerator parts
    double num_eexx_part1 = -((2.0 * eps0 * epr[0][2] + dt * sigma[0][2]) * (2.0 * eps0 * epr[1][1] + dt * sigma[1][1])) + (2.0 * eps0 * epr[0][1] + dt * sigma[0][1]) * (2.0 * eps0 * epr[1][2] + dt * sigma[1][2]);
    double num_eexx_part2 = (eps0 * epr[2][0]) / dt - sigma[2][0] / 2.0;
    double num_eexx_part3 = (eps0 * epr[1][0]) / dt - sigma[1][0] / 2.0;
    double num_eexx_part4 = (2.0 * eps0 * epr[0][2] + dt * sigma[0][2]) * (2.0 * eps0 * epr[2][1] + dt * sigma[2][1]) - (2.0 * eps0 * epr[0][1] + dt * sigma[0][1]) * (2.0 * eps0 * epr[2][2] + dt * sigma[2][2]);
    double num_eexx_part5 = (eps0 * epr[0][0]) / dt - sigma[0][0] / 2.0;
    double num_eexx_part6 = -((2.0 * eps0 * epr[1][2] + dt * sigma[1][2]) * (2.0 * eps0 * epr[2][1] + dt * sigma[2][1])) + (2.0 * eps0 * epr[1][1] + dt * sigma[1][1]) * (2.0 * eps0 * epr[2][2] + dt * sigma[2][2]);
    
    coeff.eexx = (num_eexx_part1 * num_eexx_part2 + num_eexx_part3 * num_eexx_part4 + num_eexx_part5 * num_eexx_part6) / denominator;

    // --- coeff%eexy ---
    // This is a complex expression. We will translate it directly.
    // Numerator
    double num_eexy_term1 = (epr[0][2] * epr[2][1] - epr[0][1] * epr[2][2]) * sigma[1][1] + epr[1][2] * (-(epr[2][1] * sigma[0][1]) + epr[0][1] * sigma[2][1]) + epr[1][1] * (epr[2][2] * sigma[0][1] - epr[0][2] * sigma[2][1]);
    double num_eexy_term2 = epr[2][1] * (sigma[0][2] * sigma[1][1] - sigma[0][1] * sigma[1][2]) + epr[1][1] * (-(sigma[0][2] * sigma[2][1]) + sigma[0][1] * sigma[2][2]) + epr[0][1] * (sigma[1][2] * sigma[2][1] - sigma[1][1] * sigma[2][2]);
    double num_eexy_num = 4.0 * dt * eps0 * (2.0 * eps0 * num_eexy_term1 + dt * num_eexy_term2);
    
    // Denominator part 1 (same as base denominator structure but with different sigma/dt parts? No, looking at Fortran:
    // Denom is: 8 * eps0**3 * (det_epr) + ...
    // Let's calculate the determinant of epr part:
    double det_epr_part = epr[0][2] * (epr[1][1] * epr[2][0] - epr[1][0] * epr[2][1]) + epr[0][1] * (-(epr[1][2] * epr[2][0]) + epr[1][0] * epr[2][2]) + epr[0][0] * (epr[1][2] * epr[2][1] - epr[1][1] * epr[2][2]);
    
    double num_eexy_term3 = epr[2][0] * sigma[0][2] * sigma[1][1] + epr[2][2] * (sigma[0][1] * sigma[1][0] - sigma[0][0] * sigma[1][1]) - epr[2][0] * sigma[0][1] * sigma[1][2] + epr[2][1] * (-(sigma[0][2] * sigma[1][0]) + sigma[0][0] * sigma[1][2]) - epr[1][2] * sigma[0][1] * sigma[2][0] + epr[1][1] * sigma[0][2] * sigma[2][0] + epr[0][2] * sigma[1][1] * sigma[2][0] - epr[0][1] * sigma[1][2] * sigma[2][0] + epr[1][2] * sigma[0][0] * sigma[2][1] - epr[1][0] * sigma[0][2] * sigma[2][1] - epr[0][2] * sigma[1][0] * sigma[2][1] + epr[0][0] * sigma[1][2] * sigma[2][1] + (-(epr[1][1] * sigma[0][0]) + epr[1][0] * sigma[0][1] + epr[0][1] * sigma[1][0] - epr[0][0] * sigma[1][1]) * sigma[2][2];
    
    double num_eexy_term4 = epr[1][0] * epr[2][2] * sigma[0][1] - epr[1][0] * epr[2][1] * sigma[0][2] - epr[0][2] * epr[2][1] * sigma[1][0] + epr[0][1] * epr[2][2] * sigma[1][0] + epr[0][2] * epr[2][0] * sigma[1][1] - epr[0][0] * epr[2][2] * sigma[1][1] - epr[0][1] * epr[2][0] * sigma[1][2] + epr[0][0] * epr[2][1] * sigma[1][2] - epr[0][2] * epr[1][0] * sigma[2][1] + epr[1][2] * (epr[2][1] * sigma[0][0] - epr[2][0] * sigma[0][1] - epr[0][1] * sigma[2][0] + epr[0][0] * sigma[2][1]) + epr[0][1] * epr[1][0] * sigma[2][2] + epr[1][1] * (-(epr[2][2] * sigma[0][0]) + epr[2][0] * sigma[0][2] + epr[0][2] * sigma[2][0] - epr[0][0] * sigma[2][2]);
    
    double num_eexy_term5 = sigma[0][2] * (sigma[1][1] * sigma[2][0] - sigma[1][0] * sigma[2][1]) + sigma[0][1] * (-(sigma[1][2] * sigma[2][0]) + sigma[1][0] * sigma[2][2]) + sigma[0][0] * (sigma[1][2] * sigma[2][1] - sigma[1][1] * sigma[2][2]);
    
    double num_eexy = num_eexy_num / (8.0 * std::pow(eps0, 3) * det_epr) + 2.0 * dt * dt * eps0 * num_eexy_term3 + 4.0 * dt * eps0 * eps0 * num_eexy_term4 + dt * dt * dt * num_eexy_term5;
    
    coeff.eexy = num_eexy;

    // --- coeff%eexz ---
    // Very similar to eexy but with index swaps for z (index 2)
    double num_eexz_term1 = (epr[0][2] * epr[2][1] - epr[0][1] * epr[2][2]) * sigma[1][2] + epr[1][2] * (-(epr[2][1] * sigma[0][2]) + epr[0][1] * sigma[2][2]) + epr[1][1] * (epr[2][2] * sigma[0][2] - epr[0][2] * sigma[2][2]);
    double num_eexz_term2 = epr[2][2] * (sigma[0][2] * sigma[1][1] - sigma[0][1] * sigma[1][2]) + epr[1][2] * (-(sigma[0][2] * sigma[2][1]) + sigma[0][1] * sigma[2][2]) + epr[0][2] * (sigma[1][2] * sigma[2][1] - sigma[1][1] * sigma[2][2]);
    double num_eexz_num = 4.0 * dt * eps0 * (2.0 * eps0 * num_eexz_term1 + dt * num_eexz_term2);
    
    double num_eexz_term3 = epr[2][0] * sigma[0][2] * sigma[1][1] + epr[2][2] * (sigma[0][1] * sigma[1][0] - sigma[0][0] * sigma[1][1]) - epr[2][0] * sigma[0][1] * sigma[1][2] + epr[2][1] * (-(sigma[0][2] * sigma[1][0]) + sigma[0][0] * sigma[1][2]) - epr[1][2] * sigma[0][1] * sigma[2][0] + epr[1][1] * sigma[0][2] * sigma[2][0] + epr[0][2] * sigma[1][1] * sigma[2][0] - epr[0][1] * sigma[1][2] * sigma[2][0] + epr[1][2] * sigma[0][0] * sigma[2][1] - epr[1][0] * sigma[0][2] * sigma[2][1] - epr[0][2] * sigma[1][0] * sigma[2][1] + epr[0][0] * sigma[1][2] * sigma[2][1] + (-(epr[1][1] * sigma[0][0]) + epr[1][0] * sigma[0][1] + epr[0][1] * sigma[1][0] - epr[0][0] * sigma[1][1]) * sigma[2][2];
    
    double num_eexz_term4 = epr[1][0] * epr[2][2] * sigma[0][1] - epr[1][0] * epr[2][1] * sigma[0][2] - epr[0][2] * epr[2][1] * sigma[1][0] + epr[0][1] * epr[2][2] * sigma[1][0] + epr[0][2] * epr[2][0] * sigma[1][1] - epr[0][0] * epr[2][2] * sigma[1][1] - epr[0][1] * epr[2][0] * sigma[1][2] + epr[0][0] * epr[2][1] * sigma[1][2] - epr[0][2] * epr[1][0] * sigma[2][1] + epr[1][2] * (epr[2][1] * sigma[0][0] - epr[2][0] * sigma[0][1] - epr[0][1] * sigma[2][0] + epr[0][0] * sigma[2][1]) + epr[0][1] * epr[1][0] * sigma[2][2] + epr[1][1] * (-(epr[2][2] * sigma[0][0]) + epr[2][0] * sigma[0][2] + epr[0][2] * sigma[2][0] - epr[0][0] * sigma[2][2]);
    
    double num_eexz_term5 = sigma[0][2] * (sigma[1][1] * sigma[2][0] - sigma[1][0] * sigma[2][1]) + sigma[0][1] * (-(sigma[1][2] * sigma[2][0]) + sigma[1][0] * sigma[2][2]) + sigma[0][0] * (sigma[1][2] * sigma[2][1] - sigma[1][1] * sigma[2][2]);
    
    coeff.eexz = num_eexz_num / (8.0 * std::pow(eps0, 3) * det_epr) + 2.0 * dt * dt * eps0 * num_eexz_term3 + 4.0 * dt * eps0 * eps0 * num_eexz_term4 + dt * dt * dt * num_eexz_term5;

    // --- coeff%eeyx ---
    double num_eeyx_term1 = (epr[0][2] * epr[2][0] - epr[0][0] * epr[2][2]) * sigma[1][0] + epr[1][2] * (-(epr[2][0] * sigma[0][0]) + epr[0][0] * sigma[2][0]) + epr[1][0] * (epr[2][2] * sigma[0][0] - epr[0][2] * sigma[2][0]);
    double num_eeyx_term2 = epr[2][0] * (sigma[0][2] * sigma[1][0] - sigma[0][0] * sigma[1][2]) + epr[1][0] * (-(sigma[0][2] * sigma[2][0]) + sigma[0][0] * sigma[2][2]) + epr[0][0] * (sigma[1][2] * sigma[2][0] - sigma[1][0] * sigma[2][2]);
    double num_eeyx_num = 4.0 * dt * eps0 * (2.0 * eps0 * num_eeyx_term1 + dt * num_eeyx_term2);
    
    double det_epr_part2 = epr[0][2] * (-(epr[1][1] * epr[2][0]) + epr[1][0] * epr[2][1]) + epr[0][1] * (epr[1][2] * epr[2][0] - epr[1][0] * epr[2][2]) + epr[0][0] * (-(epr[1][2] * epr[2][1]) + epr[1][1] * epr[2][2]);
    
    double num_eeyx_term3 = -(epr[2][0] * sigma[0][2] * sigma[1][1]) + epr[2][2] * (-(sigma[0][1] * sigma[1][0]) + sigma[0][0] * sigma[1][1]) + epr[2][0] * sigma[0][1] * sigma[1][2] + epr[2][1] * (sigma[0][2] * sigma[1][0] - sigma[0][0] * sigma[1][2]) + epr[1][2] * sigma[0][1] * sigma[2][0] - epr[1][1] * sigma[0][2] * sigma[2][0] - epr[0][2] * sigma[1][1] * sigma[2][0] + epr[0][1] * sigma[1][2] * sigma[2][0] - epr[1][2] * sigma[0][0] * sigma[2][1] + epr[1][0] * sigma[0][2] * sigma[2][1] + epr[0][2] * sigma[1][0] * sigma[2][1] - epr[0][0] * sigma[1][2] * sigma[2][1] + (epr[1][1] * sigma[0][0] - epr[1][0] * sigma[0][1] - epr[0][1] * sigma[1][0] + epr[0][0] * sigma[1][1]) * sigma[2][2];
    
    double num_eeyx_term4 = -(epr[1][0] * epr[2][2] * sigma[0][1]) + epr[1][0] * epr[2][1] * sigma[0][2] + epr[0][2] * epr[2][1] * sigma[1][0] - epr[0][1] * epr[2][2] * sigma[1][0] - epr[0][2] * epr[2][0] * sigma[1][1] + epr[0][0] * epr[2][2] * sigma[1][1] + epr[0][1] * epr[2][0] * sigma[1][2] - epr[0][0] * epr[2][1] * sigma[1][2] + epr[0][2] * epr[1][0] * sigma[2][1] + epr[1][2] * (-(epr[2][1] * sigma[0][0]) + epr[2][0] * sigma[0][1] + epr[0][1] * sigma[2][0] - epr[0][0] * sigma[2][1]) - epr[0][1] * epr[1][0] * sigma[2][2] + epr[1][1] * (epr[2][2] * sigma[0][0] - epr[2][0] * sigma[0][2] - epr[0][2] * sigma[2][0] + epr[0][0] * sigma[2][2]);
    
    double num_eeyx_term5 = sigma[0][2] * (-(sigma[1][1] * sigma[2][0]) + sigma[1][0] * sigma[2][1]) + sigma[0][1] * (sigma[1][2] * sigma[2][0] - sigma[1][0] * sigma[2][2]) + sigma[0][0] * (-(sigma[1][2] * sigma[2][1]) + sigma[1][1] * sigma[2][2]);
    
    coeff.eeyx = num_eeyx_num / (8.0 * std::pow(eps0, 3) * det_epr_part2) + 2.0 * dt * dt * eps0 * num_eeyx_term3 + 4.0 * dt * eps0 * eps0 * num_eeyx_term4 + dt * dt * dt * num_eeyx_term5;

    // --- coeff%eeyy ---
    double num_eeyy_part1 = (2.0 * eps0 * epr[0][2] + dt * sigma[0][2]) * (2.0 * eps0 * epr[1][0] + dt * sigma[1][0]) - (2.0 * eps0 * epr[0][0] + dt * sigma[0][0]) * (2.0 * eps0 * epr[1][2] + dt * sigma[1][2]);
    double num_eeyy_part2 = (eps0 * epr[2][1]) / dt - sigma[2][1] / 2.0;
    double num_eeyy_part3 = (eps0 * epr[1][1]) / dt - sigma[1][1] / 2.0;
    double num_eeyy_part4 = -((2.0 * eps0 * epr[0][2] + dt * sigma[0][2]) * (2.0 * eps0 * epr[2][0] + dt * sigma[2][0])) + (2.0 * eps0 * epr[0][0] + dt * sigma[0][0]) * (2.0 * eps0 * epr[2][2] + dt * sigma[2][2]);
    double num_eeyy_part5 = (eps0 * epr[0][1]) / dt - sigma[0][1] / 2.0;
    double num_eeyy_part6 = (2.0 * eps0 * epr[1][2] + dt * sigma[1][2]) * (2.0 * eps0 * epr[2][0] + dt * sigma[2][0]) - (2.0 * eps0 * epr[1][0] + dt * sigma[1][0]) * (2.0 * eps0 * epr[2][2] + dt * sigma[2][2]);
    
    coeff.eeyy = (num_eeyy_part1 * num_eeyy_part2 + num_eeyy_part3 * num_eeyy_part4 + num_eeyy_part5 * num_eeyy_part6) / denominator;

    // --- coeff%eeyz ---
    double num_eeyz_term1 = (-(epr[0][2] * epr[2][0]) + epr[0][0] * epr[2][2]) * sigma[1][2] + epr[1][2] * (epr[2][0] * sigma[0][2] - epr[0][0] * sigma[2][2]) + epr[1][0] * (-(epr[2][2] * sigma[0][2]) + epr[0][2] * sigma[2][2]);
    double num_eeyz_term2 = epr[2][2] * (-(sigma[0][2] * sigma[1][0]) + sigma[0][0] * sigma[1][2]) + epr[1][2] * (sigma[0][2] * sigma[2][0] - sigma[0][0] * sigma[2][2]) + epr[0][2] * (-(sigma[1][2] * sigma[2][0]) + sigma[1][0] * sigma[2][2]);
    double num_eeyz_num = 4.0 * dt * eps0 * (2.0 * eps0 * num_eeyz_term1 + dt * num_eeyz_term2);
    
    double num_eeyz_term3 = epr[2][0] * sigma[0][2] * sigma[1][1] + epr[2][2] * (sigma[0][1] * sigma[1][0] - sigma[0][0] * sigma[1][1]) - epr[2][0] * sigma[0][1] * sigma[1][2] + epr[2][1] * (-(sigma[0][2] * sigma[1][0]) + sigma[0][0] * sigma[1][2]) - epr[1][2] * sigma[0][1] * sigma[2][0] + epr[1][1] * sigma[0][2] * sigma[2][0] + epr[0][2] * sigma[1][1] * sigma[2][0] - epr[0][1] * sigma[1][2] * sigma[2][0] + epr[1][2] * sigma[0][0] * sigma[2][1] - epr[1][0] * sigma[0][2] * sigma[2][1] - epr[0][2] * sigma[1][0] * sigma[2][1] + epr[0][0] * sigma[1][2] * sigma[2][1] + (-(epr[1][1] * sigma[0][0]) + epr[1][0] * sigma[0][1] + epr[0][1] * sigma[1][0] - epr[0][0] * sigma[1][1]) * sigma[2][2];
    
    double num_eeyz_term4 = epr[1][0] * epr[2][2] * sigma[0][1] - epr[1][0] * epr[2][1] * sigma[0][2] - epr[0][2] * epr[2][1] * sigma[1][0] + epr[0][1] * epr[2][2] * sigma[1][0] + epr[0][2] * epr[2][0] * sigma[1][1] - epr[0][0] * epr[2][2] * sigma[1][1] - epr[0][1] * epr[2][0] * sigma[1][2] + epr[0][0] * epr[2][1] * sigma[1][2] - epr[0][2] * epr[1][0] * sigma[2][1] + epr[1][2] * (epr[2][1] * sigma[0][0] - epr[2][0] * sigma[0][1] - epr[0][1] * sigma[2][0] + epr[0][0] * sigma[2][1]) + epr[0][1] * epr[1][0] * sigma[2][2] + epr[1][1] * (-(epr[2][2] * sigma[0][0]) + epr[2][0] * sigma[0][2] + epr[0][2] * sigma[2][0] - epr[0][0] * sigma[2][2]);
    
    double num_eeyz_term5 = sigma[0][2] * (sigma[1][1] * sigma[2][0] - sigma[1][0] * sigma[2][1]) + sigma[0][1] * (-(sigma[1][2] * sigma[2][0]) + sigma[1][0] * sigma[2][2]) + sigma[0][0] * (sigma[1][2] * sigma[2][1] - sigma[1][1] * sigma[2][2]);
    
    coeff.eeyz = num_eeyz_num / (8.0 * std::pow(eps0, 3) * det_epr) + 2.0 * dt * dt * eps0 * num_eeyz_term3 + 4.0 * dt * eps0 * eps0 * num_eeyz_term4 + dt * dt * dt * num_eeyz_term5;

    // --- coeff%eezx ---
    double num_eezx_term1 = (epr[0][1] * epr[2][0] - epr[0][0] * epr[2][1]) * sigma[1][0] + epr[1][1] * (-(epr[2][0] * sigma[0][0]) + epr[0][0] * sigma[2][0]) + epr[1][0] * (epr[2][1] * sigma[0][0] - epr[0][1] * sigma[2][0]);
    double num_eezx_term2 = epr[2][0] * (sigma[0][1] * sigma[1][0] - sigma[0][0] * sigma[1][1]) + epr[1][0] * (-(sigma[0][1] * sigma[2][0]) + sigma[0][0] * sigma[2][1]) + epr[0][0] * (sigma[1][1] * sigma[2][0] - sigma[1][0] * sigma[2][1]);
    double num_eezx_num = 4.0 * dt * eps0 * (2.0 * eps0 * num_eezx_term1 + dt * num_eezx_term2);
    
    double num_eezx_term3 = epr[2][0] * sigma[0][2] * sigma[1][1] + epr[2][2] * (sigma[0][1] * sigma[1][0] - sigma[0][0] * sigma[1][1]) - epr[2][0] * sigma[0][1] * sigma[1][2] + epr[2][1] * (-(sigma[0][2] * sigma[1][0]) + sigma[0][0] * sigma[1][2]) - epr[1][2] * sigma[0][1] * sigma[2][0] + epr[1][1] * sigma[0][2] * sigma[2][0] + epr[0][2] * sigma[1][1] * sigma[2][0] - epr[0][1] * sigma[1][2] * sigma[2][0] + epr[1][2] * sigma[0][0] * sigma[2][1] - epr[1][0] * sigma[0][2] * sigma[2][1] - epr[0][2] * sigma[1][0] * sigma[2][1] + epr[0][0] * sigma[1][2] * sigma[2][1] + (-(epr[1][1] * sigma[0][0]) + epr[1][0] * sigma[0][1] + epr[0][1] * sigma[1][0] - epr[0][0] * sigma[1][1]) * sigma[2][2];
    
    double num_eezx_term4 = epr[1][0] * epr[2][2] * sigma[0][1] - epr[1][0] * epr[2][1] * sigma[0][2] - epr[0][2] * epr[2][1] * sigma[1][0] + epr[0][1] * epr[2][2] * sigma[1][0] + epr[0][2] * epr[2][0] * sigma[1][1] - epr[0][0] * epr[2][2] * sigma[1][1] - epr[0][1] * epr[2][0] * sigma[1][2] + epr[0][0] * epr[2][1] * sigma[1][2] - epr[0][2] * epr[1][0] * sigma[2][1] + epr[1][2] * (epr[2][1] * sigma[0][0] - epr[2][0] * sigma[0][1] - epr[0][1] * sigma[2][0] + epr[0][0] * sigma[2][1]) + epr[0][1] * epr[1][0] * sigma[2][2] + epr[1][1] * (-(epr[2][2] * sigma[0][0]) + epr[2][0] * sigma[0][2] + epr[0][2] * sigma[2][0] - epr[0][0] * sigma[2][2]);
    
    double num_eezx_term5 = sigma[0][2] * (sigma[1][1] * sigma[2][0] - sigma[1][0] * sigma[2][1]) + sigma[0][1] * (-(sigma[1][2] * sigma[2][0]) + sigma[1][0] * sigma[2][2]) + sigma[0][0] * (sigma[1][2] * sigma[2][1] - sigma[1][1] * sigma[2][2]);
    
    coeff.eezx = num_eezx_num / (8.0 * std::pow(eps0, 3) * det_epr) + 2.0 * dt * dt * eps0 * num_eezx_term3 + 4.0 * dt * eps0 * eps0 * num_eezx_term4 + dt * dt * dt * num_eezx_term5;

    // --- coeff%eezy ---
    double num_eezy_term1 = (-(epr[0][1] * epr[2][0]) + epr[0][0] * epr[2][1]) * sigma[1][1] + epr[1][1] * (epr[2][0] * sigma[0][1] - epr[0][0] * sigma[2][1]) + epr[1][0] * (-(epr[2][1] * sigma[0][1]) + epr[0][1] * sigma[2][1]);
    double num_eezy_term2 = epr[2][1] * (-(sigma[0][1] * sigma[1][0]) + sigma[0][0] * sigma[1][1]) + epr[1][1] * (sigma[0][1] * sigma[2][0] - sigma[0][0] * sigma[2][1]) + epr[0][1] * (-(sigma[1][1] * sigma[2][0]) + sigma[1][0] * sigma[2][1]);
    double num_eezy_num = 4.0 * dt * eps0 * (2.0 * eps0 * num_eezy_term1 + dt * num_eezy_term2);
    
    double det_epr_part3 = epr[0][2] * (-(epr[1][1] * epr[2][0]) + epr[1][0] * epr[2][1]) + epr[0][1] * (epr[1][2] * epr[2][0] - epr[1][0] * epr[2][2]) + epr[0][0] * (-(epr[1][2] * epr[2][1]) + epr[1][1] * epr[2][2]);
    
    double num_eezy_term3 = -(epr[2][0] * sigma[0][2] * sigma[1][1]) + epr[2][2] * (-(sigma[0][1] * sigma[1][0]) + sigma[0][0] * sigma[1][1]) + epr[2][0] * sigma[0][1] * sigma[1][2] + epr[2][1] * (sigma[0][2] * sigma[1][0] - sigma[0][0] * sigma[1][2]) + epr[1][2] * sigma[0][1] * sigma[2][0] - epr[1][1] * sigma[0][2] * sigma[2][0] - epr[0][2] * sigma[1][1] * sigma[2][0] + epr[0][1] * sigma[1][2] * sigma[2][0] - epr[1][2] * sigma[0][0] * sigma[2][1] + epr[1][0] * sigma[0][2] * sigma[2][1] + epr[0][2] * sigma[1][0] * sigma[2][1] - epr[0][0] * sigma[1][2] * sigma[2][1] + (epr[1][1] * sigma[0][0] - epr[1][0] * sigma[0][1] - epr[0][1] * sigma[1][0] + epr[0][0] * sigma[1][1]) * sigma[2][2];
    
    double num_eezy_term4 = -(epr[1][0] * epr[2][2] * sigma[0][1]) + epr[1][0] * epr[2][1] * sigma[0][2] + epr[0][2] * epr[2][1] * sigma[1][0] - epr[0][1] * epr[2][2] * sigma[1][0] - epr[0][2] * epr[2][0] * sigma[1][1] + epr[0][0] * epr[2][2] * sigma[1][1] + epr[0][1] * epr[2][0] * sigma[1][2] - epr[0][0] * epr[2][1] * sigma[1][2] + epr[0][2] * epr[1][0] * sigma[2][1] + epr[1][2] * (-(epr[2][1] * sigma[0][0]) + epr[2][0] * sigma[0][1] + epr[0][1] * sigma[2][0] - epr[0][0] * sigma[2][1]) - epr[0][1] * epr[1][0] * sigma[2][2] + epr[1][1] * (epr[2][2] * sigma[0][0] - epr[2][0] * sigma[0][2] - epr[0][2] * sigma[2][0] + epr[0][0] * sigma[2][2]);
    
    double num_eezy_term5 = sigma[0][2] * (-(sigma[1][1] * sigma[2][0]) + sigma[1][0] * sigma[2][1]) + sigma[0][1] * (sigma[1][2] * sigma[2][0] - sigma[1][0] * sigma[2][2]) + sigma[0][0] * (-(sigma[1][2] * sigma[2][1]) + sigma[1][1] * sigma[2][2]);
    
    coeff.eezy = num_eezy_num / (8.0 * std::pow(eps0, 3) * det_epr_part3) + 2.0 * dt * dt * eps0 * num_eezy_term3 + 4.0 * dt * eps0 * eps0 * num_eezy_term4 + dt * dt * dt * num_eezy_term5;

    // --- coeff%eezz ---
    double num_eezz_part1 = (eps0 * epr[1][2]) / dt - sigma[1][2] / 2.0;
    double num_eezz_part2 = (2.0 * eps0 * epr[0][1] + dt * sigma[0][1]) * (2.0 * eps0 * epr[2][0] + dt * sigma[2][0]) - (2.0 * eps0 * epr[0][0] + dt * sigma[0][0]) * (2.0 * eps0 * epr[2][1] + dt * sigma[2][1]);
    double num_eezz_part3 = (eps0 * epr[0][2]) / dt - sigma[0][2] / 2.0;
    double num_eezz_part4 = -((2.0 * eps0 * epr[1][1] + dt * sigma[1][1]) * (2.0 * eps0 * epr[2][0] + dt * sigma[2][0])) + (2.0 * eps0 * epr[1][0] + dt * sigma[1][0]) * (2.0 * eps0 * epr[2][1] + dt * sigma[2][1]);
    double num_eezz_part5 = -((2.0 * eps0 * epr[0][1] + dt * sigma[0][1]) * (2.0 * eps0 * epr[1][0] + dt * sigma[1][0])) + (2.0 * eps0 * epr[0][0] + dt * sigma[0][0]) * (2.0 * eps0 * epr[1][1] + dt * sigma[1][1]);
    double num_eezz_part6 = (eps0 * epr[2][2]) / dt - sigma[2][2] / 2.0;
    
    coeff.eezz = (num_eezz_part1 * num_eezz_part2 + num_eezz_part3 * num_eezz_part4 + num_eezz_part5 * num_eezz_part6) / denominator;

    // --- coeff%ehxx ---
    double num_ehxx = -((2.0 * eps0 * epr[1][2] + dt * sigma[1][2]) * (2.0 * eps0 * epr[2][1] + dt * sigma[2][1])) + (2.0 * eps0 * epr[1][1] + dt * sigma[1][1]) * (2.0 * eps0 * epr[2][2] + dt * sigma[2][2]);
    coeff.ehxx = num_ehxx / denominator;

    // --- coeff%ehxy ---
    double num_ehxy = (2.0 * eps0 * epr[0][2] + dt * sigma[0][2]) * (2.0 * eps0 * epr[2][1] + dt * sigma[2][1]) - (2.0 * eps0 * epr[0][1] + dt * sigma[0][1]) * (2.0 * eps0 * epr[2][2] + dt * sigma[2][2]);
    coeff.ehxy = num_ehxy / denominator;

    // --- coeff%ehxz ---
    double num_ehzx = -((2.0 * eps0 * epr[0][2] + dt * sigma[0][2]) * (2.0 * eps0 * epr[1][1] + dt * sigma[1][1])) + (2.0 * eps0 * epr[0][1] + dt * sigma[0][1]) * (2.0 * eps0 * epr[1][2] + dt * sigma[1][2]);
    coeff.ehxz = num_ehzx / denominator;

    // --- coeff%ehyx ---
    double num_ehyx = (2.0 * eps0 * epr[1][2] + dt * sigma[1][2]) * (2.0 * eps0 * epr[2][0] + dt * sigma[2][0]) - (2.0 * eps0 * epr[1][0] + dt * sigma[1][0]) * (2.0 * eps0 * epr[2][2] + dt * sigma[2][2]);
    coeff.ehyx = num_ehyx / denominator;

    // --- coeff%ehyy ---
    double num_ehyy = -((2.0 * eps0 * epr[0][2] + dt * sigma[0][2]) * (2.0 * eps0 * epr[2][0] + dt * sigma[2][0])) + (2.0 * eps0 * epr[0][0] + dt * sigma[0][0]) * (2.0 * eps0 * epr[2][2] + dt * sigma[2][2]);
    coeff.ehyy = num_ehyy / denominator;

    // --- coeff%ehyz ---
    double num_ehyz = (2.0 * eps0 * epr[0][2] + dt * sigma[0][2]) * (2.0 * eps0 * epr[1][0] + dt * sigma[1][0]) - (2.0 * eps0 * epr[0][0] + dt * sigma[0][0]) * (2.0 * eps0 * epr[1][2] + dt * sigma[1][2]);
    coeff.ehyz = num_ehyz / denominator;

    // --- coeff%ehzx ---
    double num_ehzx2 = -((2.0 * eps0 * epr[1][1] + dt * sigma[1][1]) * (2.0 * eps0 * epr[2][0] + dt * sigma[2][0])) + (2.0 * eps0 * epr[1][0] + dt * sigma[1][0]) * (2.0 * eps0 * epr[2][1] + dt * sigma[2][1]);
    coeff.ehzx = num_ehzx2 / denominator;

    // --- coeff%ehzy ---
    // The Fortran code cuts off here, but based on the pattern, it would be:
    // coeff%ehzy = ...
    // Since the input is truncated, we stop here.
}

coeff.ehzz = ((-(2 * eps0 * epr(1, 2) + dt * sigma(1, 2)) * (2 * eps0 * epr(2, 1) + dt * sigma(2, 1))) + (2 * eps0 * epr(1, 1) + dt * sigma(1, 1)) * (2 * eps0 * epr(2, 2) + dt * sigma(2, 2))) / ((-(2 * eps0 * epr(1, 3) + dt * sigma(1, 3)) * (2 * eps0 * epr(2, 2) + dt * sigma(2, 2))) + (2 * eps0 * epr(1, 2) + dt * sigma(1, 2)) * (2 * eps0 * epr(2, 3) + dt * sigma(2, 3))) * ((eps0 * epr(3, 1)) / dt + sigma(3, 1) / 2.0) - ((-(2 * eps0 * epr(1, 3) + dt * sigma(1, 3)) * (2 * eps0 * epr(2, 1) + dt * sigma(2, 1))) + (2 * eps0 * epr(1, 1) + dt * sigma(1, 1)) * (2 * eps0 * epr(2, 3) + dt * sigma(2, 3))) * ((eps0 * epr(3, 2)) / dt + sigma(3, 2) / 2.0) + ((-(2 * eps0 * epr(1, 2) + dt * sigma(1, 2)) * (2 * eps0 * epr(2, 1) + dt * sigma(2, 1))) + (2 * eps0 * epr(1, 1) + dt * sigma(1, 1)) * (2 * eps0 * epr(2, 2) + dt * sigma(2, 2))) * ((eps0 * epr(3, 3)) / dt + sigma(3, 3) / 2.0);

        coeff.hhxx = ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) + (2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 1)) / dt - sigmam(3, 1) / 2.0) + ((mu0 * mur(2, 1)) / dt - sigmam(2, 1) / 2.0) * ((2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(3, 2) + dt * sigmam(3, 2)) - (2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(3, 3) + dt * sigmam(3, 3))) + ((mu0 * mur(1, 1)) / dt - sigmam(1, 1) / 2.0) * ((-(2 * mu0 * mur(2, 3) + dt * sigmam(2, 3)) * (2 * mu0 * mur(3, 2) + dt * sigmam(3, 2))) + (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2)) * (2 * mu0 * mur(3, 3) + dt * sigmam(3, 3))) / ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) + (2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 1)) / dt + sigmam(3, 1) / 2.0) - ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 2)) / dt + sigmam(3, 2) / 2.0) + ((-(2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) * ((mu0 * mur(3, 3)) / dt + sigmam(3, 3) / 2.0);

        coeff.hhxy = (4 * dt * mu0 * (2 * mu0 * ((mur(1, 3) * mur(3, 2) - mur(1, 2) * mur(3, 3)) * sigmam(2, 2) + mur(2, 3) * (-(mur(3, 2) * sigmam(1, 2)) + mur(1, 2) * sigmam(3, 2)) + mur(2, 2) * (mur(3, 3) * sigmam(1, 2) - mur(1, 3) * sigmam(3, 2))) + dt * (mur(3, 2) * (sigmam(1, 3) * sigmam(2, 2) - sigmam(1, 2) * sigmam(2, 3)) + mur(2, 2) * (-(sigmam(1, 3) * sigmam(3, 2)) + sigmam(1, 2) * sigmam(3, 3)) + mur(1, 2) * (sigmam(2, 3) * sigmam(3, 2) - sigmam(2, 2) * sigmam(3, 3))))) / (8 * mu0 * mu0 * mu0 * (mur(1, 3) * (mur(2, 2) * mur(3, 1) - mur(2, 1) * mur(3, 2)) + mur(1, 2) * (-(mur(2, 3) * mur(3, 1)) + mur(2, 1) * mur(3, 3)) + mur(1, 1) * (mur(2, 3) * mur(3, 2) - mur(2, 2) * mur(3, 3))) + 2 * dt * dt * mu0 * (mur(3, 1) * sigmam(1, 3) * sigmam(2, 2) + mur(3, 3) * (sigmam(1, 2) * sigmam(2, 1) - sigmam(1, 1) * sigmam(2, 2)) - mur(3, 1) * sigmam(1, 2) * sigmam(2, 3) + mur(3, 2) * (-(sigmam(1, 3) * sigmam(2, 1)) + sigmam(1, 1) * sigmam(2, 3)) - mur(2, 3) * sigmam(1, 2) * sigmam(3, 1) + mur(2, 2) * sigmam(1, 3) * sigmam(3, 1) + mur(1, 3) * sigmam(2, 2) * sigmam(3, 1) - mur(1, 2) * sigmam(2, 3) * sigmam(3, 1) + mur(2, 3) * sigmam(1, 1) * sigmam(3, 2) - mur(2, 1) * sigmam(1, 3) * sigmam(3, 2) - mur(1, 3) * sigmam(2, 1) * sigmam(3, 2) + mur(1, 1) * sigmam(2, 3) * sigmam(3, 2) + (-(mur(2, 2) * sigmam(1, 1)) + mur(2, 1) * sigmam(1, 2) + mur(1, 2) * sigmam(2, 1) - mur(1, 1) * sigmam(2, 2)) * sigmam(3, 3)) + 4 * dt * mu0 * mu0 * (mur(2, 1) * mur(3, 3) * sigmam(1, 2) - mur(2, 1) * mur(3, 2) * sigmam(1, 3) - mur(1, 3) * mur(3, 2) * sigmam(2, 1) + mur(1, 2) * mur(3, 3) * sigmam(2, 1) + mur(1, 3) * mur(3, 1) * sigmam(2, 2) - mur(1, 1) * mur(3, 3) * sigmam(2, 2) - mur(1, 2) * mur(3, 1) * sigmam(2, 3) + mur(1, 1) * mur(3, 2) * sigmam(2, 3) - mur(1, 3) * mur(2, 1) * sigmam(3, 2) + mur(2, 3) * (mur(3, 2) * sigmam(1, 1) - mur(3, 1) * sigmam(1, 2) - mur(1, 2) * sigmam(3, 1) + mur(1, 1) * sigmam(3, 2)) + mur(1, 2) * mur(2, 1) * sigmam(3, 3) + mur(2, 2) * (-(mur(3, 3) * sigmam(1, 1)) + mur(3, 1) * sigmam(1, 3) + mur(1, 3) * sigmam(3, 1) - mur(1, 1) * sigmam(3, 3))) + dt * dt * dt * (sigmam(1, 3) * (sigmam(2, 2) * sigmam(3, 1) - sigmam(2, 1) * sigmam(3, 2)) + sigmam(1, 2) * (-(sigmam(2, 3) * sigmam(3, 1)) + sigmam(2, 1) * sigmam(3, 3)) + sigmam(1, 1) * (sigmam(2, 3) * sigmam(3, 2) - sigmam(2, 2) * sigmam(3, 3))));

        coeff.hhxz = (4 * dt * mu0 * (2 * mu0 * ((mur(1, 3) * mur(3, 2) - mur(1, 2) * mur(3, 3)) * sigmam(2, 3) + mur(2, 3) * (-(mur(3, 2) * sigmam(1, 3)) + mur(1, 2) * sigmam(3, 3)) + mur(2, 2) * (mur(3, 3) * sigmam(1, 3) - mur(1, 3) * sigmam(3, 3))) + dt * (mur(3, 3) * (sigmam(1, 3) * sigmam(2, 2) - sigmam(1, 2) * sigmam(2, 3)) + mur(2, 3) * (-(sigmam(1, 3) * sigmam(3, 2)) + sigmam(1, 2) * sigmam(3, 3)) + mur(1, 3) * (sigmam(2, 3) * sigmam(3, 2) - sigmam(2, 2) * sigmam(3, 3))))) / (8 * mu0 * mu0 * mu0 * (mur(1, 3) * (mur(2, 2) * mur(3, 1) - mur(2, 1) * mur(3, 2)) + mur(1, 2) * (-(mur(2, 3) * mur(3, 1)) + mur(2, 1) * mur(3, 3)) + mur(1, 1) * (mur(2, 3) * mur(3, 2) - mur(2, 2) * mur(3, 3))) + 2 * dt * dt * mu0 * (mur(3, 1) * sigmam(1, 3) * sigmam(2, 2) + mur(3, 3) * (sigmam(1, 2) * sigmam(2, 1) - sigmam(1, 1) * sigmam(2, 2)) - mur(3, 1) * sigmam(1, 2) * sigmam(2, 3) + mur(3, 2) * (-(sigmam(1, 3) * sigmam(2, 1)) + sigmam(1, 1) * sigmam(2, 3)) - mur(2, 3) * sigmam(1, 2) * sigmam(3, 1) + mur(2, 2) * sigmam(1, 3) * sigmam(3, 1) + mur(1, 3) * sigmam(2, 2) * sigmam(3, 1) - mur(1, 2) * sigmam(2, 3) * sigmam(3, 1) + mur(2, 3) * sigmam(1, 1) * sigmam(3, 2) - mur(2, 1) * sigmam(1, 3) * sigmam(3, 2) - mur(1, 3) * sigmam(2, 1) * sigmam(3, 2) + mur(1, 1) * sigmam(2, 3) * sigmam(3, 2) + (-(mur(2, 2) * sigmam(1, 1)) + mur(2, 1) * sigmam(1, 2) + mur(1, 2) * sigmam(2, 1) - mur(1, 1) * sigmam(2, 2)) * sigmam(3, 3)) + 4 * dt * mu0 * mu0 * (mur(2, 1) * mur(3, 3) * sigmam(1, 2) - mur(2, 1) * mur(3, 2) * sigmam(1, 3) - mur(1, 3) * mur(3, 2) * sigmam(2, 1) + mur(1, 2) * mur(3, 3) * sigmam(2, 1) + mur(1, 3) * mur(3, 1) * sigmam(2, 2) - mur(1, 1) * mur(3, 3) * sigmam(2, 2) - mur(1, 2) * mur(3, 1) * sigmam(2, 3) + mur(1, 1) * mur(3, 2) * sigmam(2, 3) - mur(1, 3) * mur(2, 1) * sigmam(3, 2) + mur(2, 3) * (mur(3, 2) * sigmam(1, 1) - mur(3, 1) * sigmam(1, 2) - mur(1, 2) * sigmam(3, 1) + mur(1, 1) * sigmam(3, 2)) + mur(1, 2) * mur(2, 1) * sigmam(3, 3) + mur(2, 2) * (-(mur(3, 3) * sigmam(1, 1)) + mur(3, 1) * sigmam(1, 3) + mur(1, 3) * sigmam(3, 1) - mur(1, 1) * sigmam(3, 3))) + dt * dt * dt * (sigmam(1, 3) * (sigmam(2, 2) * sigmam(3, 1) - sigmam(2, 1) * sigmam(3, 2)) + sigmam(1, 2) * (-(sigmam(2, 3) * sigmam(3, 1)) + sigmam(2, 1) * sigmam(3, 3)) + sigmam(1, 1) * (sigmam(2, 3) * sigmam(3, 2) - sigmam(2, 2) * sigmam(3, 3))));

        coeff.hhyx = (4 * dt * mu0 * (2 * mu0 * ((mur(1, 3) * mur(3, 1) - mur(1, 1) * mur(3, 3)) * sigmam(2, 1) + mur(2, 3) * (-(mur(3, 1) * sigmam(1, 1)) + mur(1, 1) * sigmam(3, 1)) + mur(2, 1) * (mur(3, 3) * sigmam(1, 1) - mur(1, 3) * sigmam(3, 1))) + dt * (mur(3, 1) * (sigmam(1, 3) * sigmam(2, 1) - sigmam(1, 1) * sigmam(2, 3)) + mur(2, 1) * (-(sigmam(1, 3) * sigmam(3, 1)) + sigmam(1, 1) * sigmam(3, 3)) + mur(1, 1) * (sigmam(2, 3) * sigmam(3, 1) - sigmam(2, 1) * sigmam(3, 3))))) / (8 * mu0 * mu0 * mu0 * (mur(1, 3) * (-(mur(2, 2) * mur(3, 1)) + mur(2, 1) * mur(3, 2)) + mur(1, 2) * (mur(2, 3) * mur(3, 1) - mur(2, 1) * mur(3, 3)) + mur(1, 1) * (-(mur(2, 3) * mur(3, 2)) + mur(2, 2) * mur(3, 3))) + 2 * dt * dt * mu0 * (-(mur(3, 1) * sigmam(1, 3) * sigmam(2, 2)) + mur(3, 3) * (-(sigmam(1, 2) * sigmam(2, 1)) + sigmam(1, 1) * sigmam(2, 2)) + mur(3, 1) * sigmam(1, 2) * sigmam(2, 3) + mur(3, 2) * (sigmam(1, 3) * sigmam(2, 1) - sigmam(1, 1) * sigmam(2, 3)) + mur(2, 3) * sigmam(1, 2) * sigmam(3, 1) - mur(2, 2) * sigmam(1, 3) * sigmam(3, 1) - mur(1, 3) * sigmam(2, 2) * sigmam(3, 1) + mur(1, 2) * sigmam(2, 3) * sigmam(3, 1) - mur(2, 3) * sigmam(1, 1) * sigmam(3, 2) + mur(2, 1) * sigmam(1, 3) * sigmam(3, 2) + mur(1, 3) * sigmam(2, 1) * sigmam(3, 2) - mur(1, 1) * sigmam(2, 3) * sigmam(3, 2) + (mur(2, 2) * sigmam(1, 1) - mur(2, 1) * sigmam(1, 2) - mur(1, 2) * sigmam(2, 1) + mur(1, 1) * sigmam(2, 2)) * sigmam(3, 3)) + 4 * dt * mu0 * mu0 * (-(mur(2, 1) * mur(3, 3) * sigmam(1, 2)) + mur(2, 1) * mur(3, 2) * sigmam(1, 3) + mur(1, 3) * mur(3, 2) * sigmam(2, 1) - mur(1, 2) * mur(3, 3) * sigmam(2, 1) - mur(1, 3) * mur(3, 1) * sigmam(2, 2) + mur(1, 1) * mur(3, 3) * sigmam(2, 2) + mur(1, 2) * mur(3, 1) * sigmam(2, 3) - mur(1, 1) * mur(3, 2) * sigmam(2, 3) + mur(1, 3) * mur(2, 1) * sigmam(3, 2) + mur(2, 3) * (-(mur(3, 2) * sigmam(1, 1)) + mur(3, 1) * sigmam(1, 2) + mur(1, 2) * sigmam(3, 1) - mur(1, 1) * sigmam(3, 2)) - mur(1, 2) * mur(2, 1) * sigmam(3, 3) + mur(2, 2) * (mur(3, 3) * sigmam(1, 1) - mur(3, 1) * sigmam(1, 3) - mur(1, 3) * sigmam(3, 1) + mur(1, 1) * sigmam(3, 3))) + dt * dt * dt * (sigmam(1, 3) * (-(sigmam(2, 2) * sigmam(3, 1)) + sigmam(2, 1) * sigmam(3, 2)) + sigmam(1, 2) * (sigmam(2, 3) * sigmam(3, 1) - sigmam(2, 1) * sigmam(3, 3)) + sigmam(1, 1) * (-(sigmam(2, 3) * sigmam(3, 2)) + sigmam(2, 2) * sigmam(3, 3))));

        coeff.hhyy = (((2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1)) - (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 2)) / dt - sigmam(3, 2) / 2.0) + ((mu0 * mur(2, 2)) / dt - sigmam(2, 2) / 2.0) * ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(3, 1) + dt * sigmam(3, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(3, 3) + dt * sigmam(3, 3))) + ((mu0 * mur(1, 2)) / dt - sigmam(1, 2) / 2.0) * ((2 * mu0 * mur(2, 3) + dt * sigmam(2, 3)) * (2 * mu0 * mur(3, 1) + dt * sigmam(3, 1)) - (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1)) * (2 * mu0 * mur(3, 3) + dt * sigmam(3, 3)))) / ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) + (2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 1)) / dt + sigmam(3, 1) / 2.0) - ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 2)) / dt + sigmam(3, 2) / 2.0) + ((-(2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) * ((mu0 * mur(3, 3)) / dt + sigmam(3, 3) / 2.0);

        coeff.hhyz = (4 * dt * mu0 * (2 * mu0 * ((-(mur(1, 3) * mur(3, 1)) + mur(1, 1) * mur(3, 3)) * sigmam(2, 3) + mur(2, 3) * (mur(3, 1) * sigmam(1, 3) - mur(1, 1) * sigmam(3, 3)) + mur(2, 1) * (-(mur(3, 3) * sigmam(1, 3)) + mur(1, 3) * sigmam(3, 3))) + dt * (mur(3, 3) * (-(sigmam(1, 3) * sigmam(2, 1)) + sigmam(1, 1) * sigmam(2, 3)) + mur(2, 3) * (sigmam(1, 3) * sigmam(3, 1) - sigmam(1, 1) * sigmam(3, 3)) + mur(1, 3) * (-(sigmam(2, 3) * sigmam(3, 1)) + sigmam(2, 1) * sigmam(3, 3))))) / (8 * mu0 * mu0 * mu0 * (mur(1, 3) * (mur(2, 2) * mur(3, 1) - mur(2, 1) * mur(3, 2)) + mur(1, 2) * (-(mur(2, 3) * mur(3, 1)) + mur(2, 1) * mur(3, 3)) + mur(1, 1) * (mur(2, 3) * mur(3, 2) - mur(2, 2) * mur(3, 3))) + 2 * dt * dt * mu0 * (mur(3, 1) * sigmam(1, 3) * sigmam(2, 2) + mur(3, 3) * (sigmam(1, 2) * sigmam(2, 1) - sigmam(1, 1) * sigmam(2, 2)) - mur(3, 1) * sigmam(1, 2) * sigmam(2, 3) + mur(3, 2) * (-(sigmam(1, 3) * sigmam(2, 1)) + sigmam(1, 1) * sigmam(2, 3)) - mur(2, 3) * sigmam(1, 2) * sigmam(3, 1) + mur(2, 2) * sigmam(1, 3) * sigmam(3, 1) + mur(1, 3) * sigmam(2, 2) * sigmam(3, 1) - mur(1, 2) * sigmam(2, 3) * sigmam(3, 1) + mur(2, 3) * sigmam(1, 1) * sigmam(3, 2) - mur(2, 1) * sigmam(1, 3) * sigmam(3, 2) - mur(1, 3) * sigmam(2, 1) * sigmam(3, 2) + mur(1, 1) * sigmam(2, 3) * sigmam(3, 2) + (-(mur(2, 2) * sigmam(1, 1)) + mur(2, 1) * sigmam(1, 2) + mur(1, 2) * sigmam(2, 1) - mur(1, 1) * sigmam(2, 2)) * sigmam(3, 3)) + 4 * dt * mu0 * mu0 * (mur(2, 1) * mur(3, 3) * sigmam(1, 2) - mur(2, 1) * mur(3, 2) * sigmam(1, 3) - mur(1, 3) * mur(3, 2) * sigmam(2, 1) + mur(1, 2) * mur(3, 3) * sigmam(2, 1) + mur(1, 3) * mur(3, 1) * sigmam(2, 2) - mur(1, 1) * mur(3, 3) * sigmam(2, 2) - mur(1, 2) * mur(3, 1) * sigmam(2, 3) + mur(1, 1) * mur(3, 2) * sigmam(2, 3) - mur(1, 3) * mur(2, 1) * sigmam(3, 2) + mur(2, 3) * (mur(3, 2) * sigmam(1, 1) - mur(3, 1) * sigmam(1, 2) - mur(1, 2) * sigmam(3, 1) + mur(1, 1) * sigmam(3, 2)) + mur(1, 2) * mur(2, 1) * sigmam(3, 3) + mur(2, 2) * (-(mur(3, 3) * sigmam(1, 1)) + mur(3, 1) * sigmam(1, 3) + mur(1, 3) * sigmam(3, 1) - mur(1, 1) * sigmam(3, 3))) + dt * dt * dt * (sigmam(1, 3) * (sigmam(2, 2) * sigmam(3, 1) - sigmam(2, 1) * sigmam(3, 2)) + sigmam(1, 2) * (-(sigmam(2, 3) * sigmam(3, 1)) + sigmam(2, 1) * sigmam(3, 3)) + sigmam(1, 1) * (sigmam(2, 3) * sigmam(3, 2) - sigmam(2, 2) * sigmam(3, 3))));

        coeff.hhzx = (4 * dt * mu0 * (2 * mu0 * ((mur(1, 2) * mur(3, 1) - mur(1, 1) * mur(3, 2)) * sigmam(2, 1) + mur(2, 2) * (-(mur(3, 1) * sigmam(1, 1)) + mur(1, 1) * sigmam(3, 1)) + mur(2, 1) * (mur(3, 2) * sigmam(1, 1) - mur(1, 2) * sigmam(3, 1))) + dt * (mur(3, 1) * (sigmam(1, 2) * sigmam(2, 1) - sigmam(1, 1) * sigmam(2, 2)) + mur(2, 1) * (-(sigmam(1, 2) * sigmam(3, 1)) + sigmam(1, 1) * sigmam(3, 2)) + mur(1, 1) * (sigmam(2, 2) * sigmam(3, 1) - sigmam(2, 1) * sigmam(3, 2))))) / (8 * mu0 * mu0 * mu0 * (mur(1, 3) * (mur(2, 2) * mur(3, 1) - mur(2, 1) * mur(3, 2)) + mur(1, 2) * (-(mur(2, 3) * mur(3, 1)) + mur(2, 1) * mur(3, 3)) + mur(1, 1) * (mur(2, 3) * mur(3, 2) - mur(2, 2) * mur(3, 3))) + 2 * dt * dt * mu0 * (mur(3, 1) * sigmam(1, 3) * sigmam(2, 2) + mur(3, 3) * (sigmam(1, 2) * sigmam(2, 1) - sigmam(1, 1) * sigmam(2, 2)) - mur(3, 1) * sigmam(1, 2) * sigmam(2, 3) + mur(3, 2) * (-(sigmam(1, 3) * sigmam(2, 1)) + sigmam(1, 1) * sigmam(2, 3)) - mur(2, 3) * sigmam(1, 2) * sigmam(3, 1) + mur(2, 2) * sigmam(1, 3) * sigmam(3, 1) + mur(1, 3) * sigmam(2, 2) * sigmam(3, 1) - mur(1, 2) * sigmam(2, 3) * sigmam(3, 1) + mur(2, 3) * sigmam(1, 1) * sigmam(3, 2) - mur(2, 1) * sigmam(1, 3) * sigmam(3, 2) - mur(1, 3) * sigmam(2, 1) * sigmam(3, 2) + mur(1, 1) * sigmam(2, 3) * sigmam(3, 2) + (-(mur(2, 2) * sigmam(1, 1)) + mur(2, 1) * sigmam(1, 2) + mur(1, 2) * sigmam(2, 1) - mur(1, 1) * sigmam(2, 2)) * sigmam(3, 3)) + 4 * dt * mu0 * mu0 * (mur(2, 1) * mur(3, 3) * sigmam(1, 2) - mur(2, 1) * mur(3, 2) * sigmam(1, 3) - mur(1, 3) * mur(3, 2) * sigmam(2, 1) + mur(1, 2) * mur(3, 3) * sigmam(2, 1) + mur(1, 3) * mur(3, 1) * sigmam(2, 2) - mur(1, 1) * mur(3, 3) * sigmam(2, 2) - mur(1, 2) * mur(3, 1) * sigmam(2, 3) + mur(1, 1) * mur(3, 2) * sigmam(2, 3) - mur(1, 3) * mur(2, 1) * sigmam(3, 2) + mur(2, 3) * (mur(3, 2) * sigmam(1, 1) - mur(3, 1) * sigmam(1, 2) - mur(1, 2) * sigmam(3, 1) + mur(1, 1) * sigmam(3, 2)) + mur(1, 2) * mur(2, 1) * sigmam(3, 3) + mur(2, 2) * (-(mur(3, 3) * sigmam(1, 1)) + mur(3, 1) * sigmam(1, 3) + mur(1, 3) * sigmam(3, 1) - mur(1, 1) * sigmam(3, 3))) + dt * dt * dt * (sigmam(1, 3) * (sigmam(2, 2) * sigmam(3, 1) - sigmam(2, 1) * sigmam(3, 2)) + sigmam(1, 2) * (-(sigmam(2, 3) * sigmam(3, 1)) + sigmam(2, 1) * sigmam(3, 3)) + sigmam(1, 1) * (sigmam(2, 3) * sigmam(3, 2) - sigmam(2, 2) * sigmam(3, 3))));

        coeff.hhzy = (4 * dt * mu0 * (2 * mu0 * ((-(mur(1, 2) * mur(3, 1)) + mur(1, 1) * mur(3, 2)) * sigmam(2, 2) + mur(2, 2) * (mur(3, 1) * sigmam(1, 2) - mur(1, 1) * sigmam(3, 2)) + mur(2, 1) * (-(mur(3, 2) * sigmam(1, 2)) + mur(1, 2) * sigmam(3, 2))) + dt * (mur(3, 2) * (-(sigmam(1, 2) * sigmam(2, 1)) + sigmam(1, 1) * sigmam(2, 2)) + mur(2, 2) * (sigmam(1, 2) * sigmam(3, 1) - sigmam(1, 1) * sigmam(3, 2)) + mur(1, 2) * (-(sigmam(2, 2) * sigmam(3, 1)) + sigmam(2, 1) * sigmam(3, 2))))) / (8 * mu0 * mu0 * mu0 * (mur(1, 3) * (-(mur(2, 2) * mur(3, 1)) + mur(2, 1) * mur(3, 2)) + mur(1, 2) * (mur(2, 3) * mur(3, 1) - mur(2, 1) * mur(3, 3)) + mur(1, 1) * (-(mur(2, 3) * mur(3, 2)) + mur(2, 2) * mur(3, 3))) + 2 * dt * dt * mu0 * (-(mur(3, 1) * sigmam(1, 3) * sigmam(2, 2)) + mur(3, 3) * (-(sigmam(1, 2) * sigmam(2, 1)) + sigmam(1, 1) * sigmam(2, 2)) + mur(3, 1) * sigmam(1, 2) * sigmam(2, 3) + mur(3, 2) * (sigmam(1, 3) * sigmam(2, 1) - sigmam(1, 1) * sigmam(2, 3)) + mur(2, 3) * sigmam(1, 2) * sigmam(3, 1) - mur(2, 2) * sigmam(1, 3) * sigmam(3, 1) - mur(1, 3) * sigmam(2, 2) * sigmam(3, 1) + mur(1, 2) * sigmam(2, 3) * sigmam(3, 1) - mur(2, 3) * sigmam(1, 1) * sigmam(3, 2) + mur(2, 1) * sigmam(1, 3) * sigmam(3, 2) + mur(1, 3) * sigmam(2, 1) * sigmam(3, 2) - mur(1, 1) * sigmam(2, 3) * sigmam(3, 2) + (mur(2, 2) * sigmam(1, 1) - mur(2, 1) * sigmam(1, 2) - mur(1, 2) * sigmam(2, 1) + mur(1, 1) * sigmam(2, 2)) * sigmam(3, 3)) + 4 * dt * mu0 * mu0 * (-(mur(2, 1) * mur(3, 3) * sigmam(1, 2)) + mur(2, 1) * mur(3, 2) * sigmam(1, 3) + mur(1, 3) * mur(3, 2) * sigmam(2, 1) - mur(1, 2) * mur(3, 3) * sigmam(2, 1) - mur(1, 3) * mur(3, 1) * sigmam(2, 2) + mur(1, 1) * mur(3, 3) * sigmam(2, 2) + mur(1, 2) * mur(3, 1) * sigmam(2, 3) - mur(1, 1) * mur(3, 2) * sigmam(2, 3) + mur(1, 3) * mur(2, 1) * sigmam(3, 2) + mur(2, 3) * (-(mur(3, 2) * sigmam(1, 1)) + mur(3, 1) * sigmam(1, 2) + mur(1, 2) * sigmam(3, 1) - mur(1, 1) * sigmam(3, 2)) - mur(1, 2) * mur(2, 1) * sigmam(3, 3) + mur(2, 2) * (mur(3, 3) * sigmam(1, 1) - mur(3, 1) * sigmam(1, 3) - mur(1, 3) * sigmam(3, 1) + mur(1, 1) * sigmam(3, 3))) + dt * dt * dt * (sigmam(1, 3) * (-(sigmam(2, 2) * sigmam(3, 1)) + sigmam(2, 1) * sigmam(3, 2)) + sigmam(1, 2) * (sigmam(2, 3) * sigmam(3, 1) - sigmam(2, 1) * sigmam(3, 3)) + sigmam(1, 1) * (-(sigmam(2, 3) * sigmam(3, 2)) + sigmam(2, 2) * sigmam(3, 3))));

        coeff.hhzz = (((mu0 * mur(2, 3)) / dt - sigmam(2, 3) / 2.0) * ((2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(3, 1) + dt * sigmam(3, 1)) - (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(3, 2) + dt * sigmam(3, 2))) + ((mu0 * mur(1, 3)) / dt - sigmam(1, 3) / 2.0) * ((-(2 * mu0 * mur(2, 2) + dt * sigmam(2, 2)) * (2 * mu0 * mur(3, 1) + dt * sigmam(3, 1))) + (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1)) * (2 * mu0 * mur(3, 2) + dt * sigmam(3, 2))) + ((-(2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) * ((mu0 * mur(3, 3)) / dt - sigmam(3, 3) / 2.0)) / ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) + (2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 1)) / dt + sigmam(3, 1) / 2.0) - ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 2)) / dt + sigmam(3, 2) / 2.0) + ((-(2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) * ((mu0 * mur(3, 3)) / dt + sigmam(3, 3) / 2.0);

        coeff.hexx = ((-(2 * mu0 * mur(2, 3) + dt * sigmam(2, 3)) * (2 * mu0 * mur(3, 2) + dt * sigmam(3, 2))) + (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2)) * (2 * mu0 * mur(3, 3) + dt * sigmam(3, 3))) / ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) + (2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 1)) / dt + sigmam(3, 1) / 2.0) - ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 2)) / dt + sigmam(3, 2) / 2.0) + ((-(2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) * ((mu0 * mur(3, 3)) / dt + sigmam(3, 3) / 2.0);

        coeff.hexy = ((2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(3, 2) + dt * sigmam(3, 2)) - (2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(3, 3) + dt * sigmam(3, 3))) / ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) + (2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 1)) / dt + sigmam(3, 1) / 2.0) - ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 2)) / dt + sigmam(3, 2) / 2.0) + ((-(2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) * ((mu0 * mur(3, 3)) / dt + sigmam(3, 3) / 2.0);

        coeff.hexz = ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) + (2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) / ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) + (2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 1)) / dt + sigmam(3, 1) / 2.0) - ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 2)) / dt + sigmam(3, 2) / 2.0) + ((-(2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) * ((mu0 * mur(3, 3)) / dt + sigmam(3, 3) / 2.0);

        coeff.heyx = ((2 * mu0 * mur(2, 3) + dt * sigmam(2, 3)) * (2 * mu0 * mur(3, 1) + dt * sigmam(3, 1)) - (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1)) * (2 * mu0 * mur(3, 3) + dt * sigmam(3, 3))) / ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) + (2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 1)) / dt + sigmam(3, 1) / 2.0) - ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 2)) / dt + sigmam(3, 2) / 2.0) + ((-(2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) * ((mu0 * mur(3, 3)) / dt + sigmam(3, 3) / 2.0);

        coeff.heyy = ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(3, 1) + dt * sigmam(3, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(3, 3) + dt * sigmam(3, 3))) / ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) + (2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 1)) / dt + sigmam(3, 1) / 2.0) - ((-(2 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 2)) / dt + sigmam(3, 2) / 2.0) + ((-(2 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) * ((mu0 * mur(3, 3)) / dt + sigmam(3, 3) / 2.0);

coeff.heyz = ((2.0 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2.0 * mu0 * mur(2, 1) + dt * sigmam(2, 1)) - (2.0 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2.0 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) / ((-((2.0 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2.0 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) + (2.0 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2.0 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 1)) / dt + sigmam(3, 1) / 2.0) - (-((2.0 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2.0 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2.0 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2.0 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 2)) / dt + sigmam(3, 2) / 2.0) + (-((2.0 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2.0 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2.0 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2.0 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) * ((mu0 * mur(3, 3)) / dt + sigmam(3, 3) / 2.0));

        coeff.hezx = (-((2.0 * mu0 * mur(2, 2) + dt * sigmam(2, 2)) * (2.0 * mu0 * mur(3, 1) + dt * sigmam(3, 1))) + (2.0 * mu0 * mur(2, 1) + dt * sigmam(2, 1)) * (2.0 * mu0 * mur(3, 2) + dt * sigmam(3, 2))) / ((-((2.0 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2.0 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) + (2.0 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2.0 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 1)) / dt + sigmam(3, 1) / 2.0) - (-((2.0 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2.0 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2.0 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2.0 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 2)) / dt + sigmam(3, 2) / 2.0) + (-((2.0 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2.0 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2.0 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2.0 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) * ((mu0 * mur(3, 3)) / dt + sigmam(3, 3) / 2.0));

        coeff.hezy = ((2.0 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2.0 * mu0 * mur(3, 1) + dt * sigmam(3, 1)) - (2.0 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2.0 * mu0 * mur(3, 2) + dt * sigmam(3, 2))) / ((-((2.0 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2.0 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) + (2.0 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2.0 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 1)) / dt + sigmam(3, 1) / 2.0) - (-((2.0 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2.0 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2.0 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2.0 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 2)) / dt + sigmam(3, 2) / 2.0) + (-((2.0 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2.0 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2.0 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2.0 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) * ((mu0 * mur(3, 3)) / dt + sigmam(3, 3) / 2.0));

        coeff.hezz = (-((2.0 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2.0 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2.0 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2.0 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) / ((-((2.0 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2.0 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) + (2.0 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2.0 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 1)) / dt + sigmam(3, 1) / 2.0) - (-((2.0 * mu0 * mur(1, 3) + dt * sigmam(1, 3)) * (2.0 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2.0 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2.0 * mu0 * mur(2, 3) + dt * sigmam(2, 3))) * ((mu0 * mur(3, 2)) / dt + sigmam(3, 2) / 2.0) + (-((2.0 * mu0 * mur(1, 2) + dt * sigmam(1, 2)) * (2.0 * mu0 * mur(2, 1) + dt * sigmam(2, 1))) + (2.0 * mu0 * mur(1, 1) + dt * sigmam(1, 1)) * (2.0 * mu0 * mur(2, 2) + dt * sigmam(2, 2))) * ((mu0 * mur(3, 3)) / dt + sigmam(3, 3) / 2.0));

        return;
    }

} // namespace Anisotropic_m