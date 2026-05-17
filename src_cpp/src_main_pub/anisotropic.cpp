```cpp
#include <vector>
#include <cmath>
#include <memory>

// Forward declarations for types defined in other modules
struct SGGFDTDINFO_t;
struct media_matrices_t;
struct Anisotropic_t;
struct XYZlimit_t;

// Assuming FDETYPES_m defines these
using rkind = double;
using RKIND = double;
using integer_4 = int32_t;
using integer_tiempo = int32_t;
using real_kind_tiempo = double;

// Placeholder for iEx, iEy, etc. usually enums or constants from FDETYPES_m
enum FieldIndex { iEx, iEy, iEz, iHx, iHy, iHz };

namespace Anisotropic_m {

    struct Coeff_t {
        double eexx, eexy, eexz, eeyx, eeyy, eeyz, eezx, eezy, eezz;
        double ehxx, ehxy, ehxz, ehyx, ehyy, ehyz, ehzx, ehzy, ehzz;
        double hexx, hexy, hexz, heyx, heyy, heyz, hezx, hezy, hezz;
        double hhxx, hhxy, hhxz, hhyx, hhyy, hhyz, hhzx, hhzy, hhzz;
    };

    struct LocalSharedElement_t {
        integer_4 times;
        std::vector<integer_4> SharedMed;
        Coeff_t coeff;
    };

    struct Anisotropicinfo_t {
        integer_4 indexmed;
        integer_4 numnodesEx, numnodesEy, numnodesEz;
        integer_4 numnodesHx, numnodesHy, numnodesHz;
        std::vector<integer_4> Ex_i, Ey_i, Ez_i, Hx_i, Hy_i, Hz_i;
        std::vector<integer_4> Ex_j, Ey_j, Ez_j, Hx_j, Hy_j, Hz_j;
        std::vector<integer_4> Ex_k, Ey_k, Ez_k, Hx_k, Hy_k, Hz_k;
        std::vector<double> Ex_value, Ey_value, Ez_value, Hx_value, Hy_value, Hz_value;
        std::vector<LocalSharedElement_t> Ex_Shared, Ey_Shared, Ez_Shared;
        std::vector<LocalSharedElement_t> Hx_Shared, Hy_Shared, Hz_Shared;
        
        Coeff_t coeff;
        bool IsOnlyThinSlot;

        double sigma[3][3];
        double epr[3][3];
        double mur[3][3];
        double sigmaM[3][3];
    };

    struct AnisotropicMed_t {
        integer_4 NumMed;
        std::vector<Anisotropicinfo_t> info;
    };

    // Global variables
    AnisotropicMed_t AniMed;
    double eps0 = 0.0;
    double mu0 = 0.0;
    double cluz = 0.0;
    double zvac = 0.0;

    void AdvanceAnisotropicE(const std::vector<XYZlimit_t>& sggAlloc,
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
                             const std::vector<double>& Idzh);

    void AdvanceAnisotropicH(const std::vector<XYZlimit_t>& sggAlloc,
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
                             const std::vector<double>& Idzh);

    void InitAnisotropic(const SGGFDTDINFO_t& sgg, const media_matrices_t& media, bool& ThereAreAnisotropic, bool& ThereAreThinSlot, double eps00, double mu00);
    void DestroyAnisotropic(SGGFDTDINFO_t& sgg);
    AnisotropicMed_t* GetMed();
    void calc_anisotropicconstants(const SGGFDTDINFO_t& sgg, double& eps00, double& mu00);
    void CalculateCoeff(Coeff_t& coeff, const double sigma[3][3], const double epr[3][3], const double mur[3][3], const double sigmam[3][3], double dt);

    // Implementation details

    void InitAnisotropic(const SGGFDTDINFO_t& sgg, const media_matrices_t& media, bool& ThereAreAnisotropic, bool& ThereAreThinSlot, double eps00, double mu00) {
        eps0 = eps00;
        mu0 = mu00;

        ThereAreAnisotropic = false;
        ThereAreThinSlot = false;
        int conta = 0;
        for (int jmed = 1; jmed <= sgg.NumMedia; ++jmed) {
            if (sgg.Med[jmed].Is.Anisotropic) {
                conta++;
            }
        }

        AniMed.NumMed = conta;
        AniMed.info.resize(AniMed.NumMed + 1); // 1-based indexing simulation
        conta = 0;
        for (int jmed = 1; jmed <= sgg.NumMedia; ++jmed) {
            if (sgg.Med[jmed].Is.Anisotropic) {
                conta++;
                AniMed.info[conta].indexmed = jmed;
                if (sgg.Med(jmed).Is.ThinSlot) {
                    AniMed.info[conta].IsOnlyThinSlot = true;
                } else {
                    AniMed.info[conta].IsOnlyThinSlot = false;
                }
            }
        }

        for (int jmed = 1; jmed <= AniMed.NumMed; ++jmed) {
            AniMed.info[jmed].sigma[0][0] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).sigma[0][0];
            AniMed.info[jmed].sigma[0][1] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).sigma[0][1];
            AniMed.info[jmed].sigma[0][2] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).sigma[0][2];
            AniMed.info[jmed].sigma[1][0] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).sigma[1][0];
            AniMed.info[jmed].sigma[1][1] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).sigma[1][1];
            AniMed.info[jmed].sigma[1][2] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).sigma[1][2];
            AniMed.info[jmed].sigma[2][0] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).sigma[2][0];
            AniMed.info[jmed].sigma[2][1] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).sigma[2][1];
            AniMed.info[jmed].sigma[2][2] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).sigma[2][2];

            AniMed.info[jmed].sigmam[0][0] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).sigmam[0][0];
            AniMed.info[jmed].sigmam[0][1] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).sigmam[0][1];
            AniMed.info[jmed].sigmam[0][2] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).sigmam[0][2];
            AniMed.info[jmed].sigmam[1][0] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).sigmam[1][0];
            AniMed.info[jmed].sigmam[1][1] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).sigmam[1][1];
            AniMed.info[jmed].sigmam[1][2] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).sigmam[1][2];
            AniMed.info[jmed].sigmam[2][0] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).sigmam[2][0];
            AniMed.info[jmed].sigmam[2][1] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).sigmam[2][1];
            AniMed.info[jmed].sigmam[2][2] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).sigmam[2][2];

            AniMed.info[jmed].mur[0][0] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).mur[0][0];
            AniMed.info[jmed].mur[0][1] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).mur[0][1];
            AniMed.info[jmed].mur[0][2] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).mur[0][2];
            AniMed.info[jmed].mur[1][0] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).mur[1][0];
            AniMed.info[jmed].mur[1][1] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).mur[1][1];
            AniMed.info[jmed].mur[1][2] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).mur[1][2];
            AniMed.info[jmed].mur[2][0] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).mur[2][0];
            AniMed.info[jmed].mur[2][1] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).mur[2][1];
            AniMed.info[jmed].mur[2][2] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).mur[2][2];

            AniMed.info[jmed].epr[0][0] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).epr[0][0];
            AniMed.info[jmed].epr[0][1] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).epr[0][1];
            AniMed.info[jmed].epr[0][2] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).epr[0][2];
            AniMed.info[jmed].epr[1][0] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).epr[1][0];
            AniMed.info[jmed].epr[1][1] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).epr[1][1];
            AniMed.info[jmed].epr[1][2] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).epr[1][2];
            AniMed.info[jmed].epr[2][0] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).epr[2][0];
            AniMed.info[jmed].epr[2][1] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).epr[2][1];
            AniMed.info[jmed].epr[2][2] = sgg.med(AniMed.info[jmed].indexmed).Anisotropic(1).epr[2][2];
        }

        for (int jmed = 1; jmed <= AniMed.NumMed; ++jmed) {
            int tempindex = AniMed.info[jmed].indexmed;
            
            // Ex
            conta = 0;
            for (int k1 = sgg.SINPMLSweep(iEx).ZI; k1 <= sgg.SINPMLSweep(iEx).ZE; ++k1) {
                for (int j1 = sgg.SINPMLSweep(iEx).YI; j1 <= sgg.SINPMLSweep(iEx).YE; ++j1) {
                    for (int i1 = sgg.SINPMLSweep(iEx).XI; i1 <= sgg.SINPMLSweep(iEx).XE; ++i1) {
                        if (media.sggMiEx(i1, j1, k1) == tempindex) conta++;
                    }
                }
            }
            ThereAreAnisotropic = ThereAreAnisotropic || (conta != 0);
            ThereAreThinSlot = ThereAreAnisotropic && AniMed.info[jmed].IsOnlyThinSlot;
            AniMed.info[jmed].numnodesEx = conta;
            AniMed.info[jmed].Ex_i.resize(conta + 1);
            AniMed.info[jmed].Ex_j.resize(conta + 1);
            AniMed.info[jmed].Ex_k.resize(conta + 1);
            AniMed.info[jmed].Ex_value.resize(conta + 1, 0.0);
            AniMed.info[jmed].Ex_Shared.resize(conta + 1);
            
            conta = 0;
            for (int k1 = sgg.SINPMLSweep(iEx).ZI; k1 <= sgg.SINPMLSweep(iEx).ZE; ++k1) {
                for (int j1 = sgg.SINPMLSweep(iEx).YI; j1 <= sgg.SINPMLSweep(iEx).YE; ++j1) {
                    for (int i1 = sgg.SINPMLSweep(iEx).XI; i1 <= sgg.SINPMLSweep(iEx).XE; ++i1) {
                        if (media.sggMiEx(i1, j1, k1) == tempindex) {
                            conta++;
                            AniMed.info[jmed].Ex_Shared[conta].times = 1;
                            AniMed.info[jmed].Ex_i[conta] = i1;
                            AniMed.info[jmed].Ex_j[conta] = j1;
                            AniMed.info[jmed].Ex_k[conta] = k1;
                        }
                    }
                }
            }

            // Ey
            conta = 0;
            for (int k1 = sgg.SINPMLSweep(iEy).ZI; k1 <= sgg.SINPMLSweep(iEy).ZE; ++k1) {
                for (int j1 = sgg.SINPMLSweep(iEy).YI; j1 <= sgg.SINPMLSweep(iEy).YE; ++j1) {
                    for (int i1 = sgg.SINPMLSweep(iEy).XI; i1 <= sgg.SINPMLSweep(iEy).XE; ++i1) {
                        if (media.sggMiEy(i1, j1, k1) == tempindex) conta++;
                    }
                }
            }
            ThereAreAnisotropic = ThereAreAnisotropic || (conta != 0);
            ThereAreThinSlot = ThereAreAnisotropic && AniMed.info[jmed].IsOnlyThinSlot;
            AniMed.info[jmed].numnodesEy = conta;
            AniMed.info[jmed].Ey_i.resize(conta + 1);
            AniMed.info[jmed].Ey_j.resize(conta + 1);
            AniMed.info[jmed].Ey_k.resize(conta + 1);
            AniMed.info[jmed].Ey_value.resize(conta + 1, 0.0);
            AniMed.info[jmed].Ey_Shared.resize(conta + 1);
            
            conta = 0;
            for (int k1 = sgg.SINPMLSweep(iEy).ZI; k1 <= sgg.SINPMLSweep(iEy).ZE; ++k1) {
                for (int j1 = sgg.SINPMLSweep(iEy).YI; j1 <= sgg.SINPMLSweep(iEy).YE; ++j1) {
                    for (int i1 = sgg.SINPMLSweep(iEy).XI; i1 <= sgg.SINPMLSweep(iEy).XE; ++i1) {
                        if (media.sggMiEy(i1, j1, k1) == tempindex) {
                            conta++;
                            AniMed.info[jmed].Ey_Shared[conta].times = 1;
                            AniMed.info[jmed].Ey_i[conta] = i1;
                            AniMed.info[jmed].Ey_j[conta] = j1;
                            AniMed.info[jmed].Ey_k[conta] = k1;
                        }
                    }
                }
            }

            // Ez
            conta = 0;
            for (int k1 = sgg.SINPMLSweep(iEz).ZI; k1 <= sgg.SINPMLSweep(iEz).ZE; ++k1) {
                for (int j1 = sgg.SINPMLSweep(iEz).YI; j1 <= sgg.SINPMLSweep(iEz).YE; ++j1) {
                    for (int i1 = sgg.SINPMLSweep(iEz).XI; i1 <= sgg.SINPMLSweep(iEz).XE; ++i1) {
                        if (media.sggMiEz(i1, j1, k1) == tempindex) conta++;
                    }
                }
            }
            ThereAreAnisotropic = ThereAreAnisotropic || (conta != 0);
            ThereAreThinSlot = ThereAreAnisotropic && AniMed.info[jmed].IsOnlyThinSlot;
            AniMed.info[jmed].numnodesEz = conta;
            AniMed.info[jmed].Ez_i.resize(conta + 1);
            AniMed.info[jmed].Ez_j.resize(conta + 1);
            AniMed.info[jmed].Ez_k.resize(conta + 1);
            AniMed.info[jmed].Ez_value.resize(conta + 1, 0.0);
            AniMed.info[jmed].Ez_Shared.resize(conta + 1);
            
            conta = 0;
            for (int k1 = sgg.SINPMLSweep(iEz).ZI; k1 <= sgg.SINPMLSweep(iEz).ZE; ++k1) {
                for (int j1 = sgg.SINPMLSweep(iEz).YI; j1 <= sgg.SINPMLSweep(iEz).YE; ++j1) {
                    for (int i1 = sgg.SINPMLSweep(iEz).XI; i1 <= sgg.SINPMLSweep(iEz).XE; ++i1) {
                        if (media.sggMiEz(i1, j1, k1) == tempindex) {
                            conta++;
                            AniMed.info[jmed].Ez_Shared[conta].times = 1;
                            AniMed.info[jmed].Ez_i[conta] = i1;
                            AniMed.info[jmed].Ez_j[conta] = j1;
                            AniMed.info[jmed].Ez_k[conta] = k1;
                        }
                    }
                }
            }

            // Hx
            conta = 0;
            for (int k1 = sgg.SINPMLSweep(iHx).ZI; k1 <= sgg.SINPMLSweep(iHx).ZE; ++k1) {
                for (int j1 = sgg.SINPMLSweep(iHx).YI; j1 <= sgg.SINPMLSweep(iHx).YE; ++j1) {
                    for (int i1 = sgg.SINPMLSweep(iHx).XI; i1 <= sgg.SINPMLSweep(iHx).XE; ++i1) {
                        if (media.sggMiHx(i1, j1, k1) == tempindex) conta++;
                    }
                }
            }
            ThereAreAnisotropic = ThereAreAnisotropic || (conta != 0);
            ThereAreThinSlot = ThereAreAnisotropic && AniMed.info[jmed].IsOnlyThinSlot;
            AniMed.info[jmed].numnodesHx = conta;
            AniMed.info[jmed].Hx_i.resize(conta + 1);
            AniMed.info[jmed].Hx_j.resize(conta + 1);
            AniMed.info[jmed].Hx_k.resize(conta + 1);
            AniMed.info[jmed].Hx_value.resize(conta + 1, 0.0);
            AniMed.info[jmed].Hx_Shared.resize(conta + 1);
            
            conta = 0;
            for (int k1 = sgg.SINPMLSweep(iHx).ZI; k1 <= sgg.SINPMLSweep(iHx).ZE; ++k1) {
                for (int j1 = sgg.SINPMLSweep(iHx).YI; j1 <= sgg.SINPMLSweep(iHx).YE; ++j1) {
                    for (int i1 = sgg.SINPMLSweep(iHx).XI; i1 <= sgg.SINPMLSweep(iHx).XE; ++i1) {
                        if (media.sggMiHx(i1, j1, k1) == tempindex) {
                            conta++;
                            AniMed.info[jmed].Hx_Shared[conta].times = 1;
                            AniMed.info[jmed].Hx_i[conta] = i1;
                            AniMed.info[jmed].Hx_j[conta] = j1;
                            AniMed.info[jmed].Hx_k[conta] = k1;
                        }
                    }
                }
            }

            // Hy
            conta = 0;
            for (int k1 = sgg.SINPMLSweep(iHy).ZI; k1 <= sgg.SINPMLSweep(iHy).ZE; ++k1) {
                for (int j1 = sgg.SINPMLSweep(iHy).YI; j1 <= sgg.SINPMLSweep(iHy).YE; ++j1) {
                    for (int i1 = sgg.SINPMLSweep(iHy).XI; i1 <= sgg.SINPMLSweep(iHy).XE; ++i1) {
                        if (media.sggMiHy(i1, j1, k1) == tempindex) conta++;
                    }
                }
            }
            ThereAreAnisotropic = ThereAreAnisotropic || (conta != 0);
            ThereAreThinSlot = ThereAreAnisotropic && AniMed.info[jmed].IsOnlyThinSlot;
            AniMed.info[jmed].numnodesHy = conta;
            AniMed.info[jmed].Hy_i.resize(conta + 1);
            AniMed.info[jmed].Hy_j.resize(conta + 1);
            AniMed.info[jmed].Hy_k.resize(conta + 1);
            AniMed.info[jmed].Hy_value.resize(conta + 1, 0.0);
            AniMed.info[jmed].Hy_Shared.resize(conta + 1);
            
            conta = 0;
            for (int k1 = sgg.SINPMLSweep(iHy).ZI; k1 <= sgg.SINPMLSweep(iHy).ZE; ++k1) {
                for (int j1 = sgg.SINPMLSweep(iHy).YI; j1 <= sgg.SINPMLSweep(iHy).YE; ++j1) {
                    for (int i1 = sgg.SINPMLSweep(iHy).XI; i1 <= sgg.SINPMLSweep(iHy).XE; ++i1) {
                        if (media.sggMiHy(i1, j1, k1) == tempindex) {
                            conta++;
                            AniMed.info[jmed].Hy_Shared[conta].times = 1;
                            AniMed.info[jmed].Hy_i[conta] = i1;
                            AniMed.info[jmed].Hy_j[conta] = j1;
                            AniMed.info[jmed].Hy_k[conta] = k1;
                        }
                    }
                }
            }

            // Hz
            conta = 0;
            for (int k1 = sgg.SINPMLSweep(iHz).ZI; k1 <= sgg.SINPMLSweep(iHz).ZE; ++k1) {
                for (int j1 = sgg.SINPMLSweep(iHz).YI; j1 <= sgg.SINPMLSweep(iHz).YE; ++j1) {
                    for (int i1 = sgg.SINPMLSweep(iHz).XI; i1 <= sgg.SINPMLSweep(iHz).XE; ++i1) {
                        if (media.sggMiHz(i1, j1, k1) == tempindex) conta++;
                    }
                }
            }
            ThereAreAnisotropic = ThereAreAnisotropic || (conta != 0);
            ThereAreThinSlot = ThereAreAnisotropic && AniMed.info[jmed].IsOnlyThinSlot;
            AniMed.info[jmed].numnodesHz = conta;
            AniMed.info[jmed].Hz_i.resize(conta + 1);
            AniMed.info[jmed].Hz_j.resize(conta + 1);
            AniMed.info[jmed].Hz_k.resize(conta + 1);
            AniMed.info[jmed].Hz_value.resize(conta + 1, 0.0);
            AniMed.info[jmed].Hz_Shared.resize(conta + 1);
            
            conta = 0;
            for (int k1 = sgg.SINPMLSweep(iHz).ZI; k1 <= sgg.SINPMLSweep(iHz).ZE; ++k1) {
                for (int j1 = sgg.SINPMLSweep(iHz).YI; j1 <= sgg.SINPMLSweep(iHz).YE; ++j1) {
                    for (int i1 = sgg.SINPMLSweep(iHz).XI; i1 <= sgg.SINPMLSweep(iHz).XE; ++i1) {
                        if (media.sggMiHz(i1, j1, k1) == tempindex) {
                            conta++;
                            AniMed.info[jmed].Hz_Shared[conta].times = 1;
                            AniMed.info[jmed].Hz_i[conta] = i1;
                            AniMed.info[jmed].Hz_j[conta] = j1;
                            AniMed.info[jmed].Hz_k[conta] = k1;
                        }
                    }
                }
            }
        }

        // Update shared times
        for (int j1 = 1; j1 <= sgg.Eshared.conta; ++j1) {
            bool found = false;
            for (int jmed = 1; jmed <= AniMed.NumMed && !found; ++jmed) {
                for (int i1 = 1; i1 <= AniMed.info[jmed].numnodesEx; ++i1) {
                    if (sgg.Eshared.elem[j1].i == AniMed.info[jmed].Ex_i[i1] &&
                        sgg.Eshared.elem[j1].j == AniMed.info[jmed].Ex_j[i1] &&
                        sgg.Eshared.elem[j1].k == AniMed.info[jmed].Ex_k[i1] &&
                        sgg.Eshared.elem[j1].Field == iEx) {
                        AniMed.info[jmed].Ex_Shared[i1].times = sgg.Eshared.elem[j1].Times;
                        if (sgg.Eshared.elem[j1].Times > 1) {
                            AniMed.info[jmed].Ex_Shared[i1].SharedMed.resize(sgg.Eshared.elem[j1].Times + 1);
                        }
                        found = true;
                    }
                }
            }
            if (!found) {
                for (int jmed = 1; jmed <= AniMed.NumMed && !found; ++jmed) {
                    for (int i1 = 1; i1 <= AniMed.info[jmed].numnodesEy; ++i1) {
                        if (sgg.Eshared.elem[j1].i == AniMed.info[jmed].Ey_i[i1] &&
                            sgg.Eshared.elem[j1].j == AniMed.info[jmed].Ey_j[i1] &&
                            sgg.Eshared.elem[j1].k == AniMed.info[jmed].Ey_k[i1] &&
                            sgg.Eshared.elem[j1].Field == iEy) {
                            AniMed.info[jmed].Ey_Shared[i1].times = sgg.Eshared.elem[j1].Times;
                            if (sgg.Eshared.elem[j1].Times > 1) {
                                AniMed.info[jmed].Ey_Shared[i1].SharedMed.resize(sgg.Eshared.elem[j1].Times + 1);
                            }
                            found = true;
                        }
                    }
                }
            }
            if (!found) {
                for (int jmed = 1; jmed <= AniMed.NumMed && !found; ++jmed) {
                    for (int i1 = 1; i1 <= AniMed.info[jmed].numnodesEz; ++i1) {
                        if (sgg.Eshared.elem[j1].i == AniMed.info[jmed].Ez_i[i1] &&
                            sgg.Eshared.elem[j1].j == AniMed.info[jmed].Ez_j[i1] &&
                            sgg.Eshared.elem[j1].k == AniMed.info[jmed].Ez_k[i1] &&
                            sgg.Eshared.elem[j1].Field == iEz) {
                            AniMed.info[jmed].Ez_Shared[i1].times = sgg.Eshared.elem[j1].Times;
                            if (sgg.Eshared.elem[j1].Times > 1) {
                                AniMed.info[jmed].Ez_Shared[i1].SharedMed.resize(sgg.Eshared.elem[j1].Times + 1);
                            }
                            found = true;
                        }
                    }
                }
            }
        }

        for (int j1 = 1; j1 <= sgg.Hshared.conta; ++j1) {
            bool found = false;
            for (int jmed = 1; jmed <= AniMed.NumMed && !found; ++jmed) {
                for (int i1 = 1; i1 <= AniMed.info[jmed].numnodesHx; ++i1) {
                    if (sgg.Hshared.elem[j1].i == AniMed.info[jmed].Hx_i[i1] &&
                        sgg.Hshared.elem[j1].j == AniMed.info[jmed].Hx_j[i1] &&
                        sgg.Hshared.elem[j1].k == AniMed.info[jmed].Hx_k[i1] &&
                        sgg.Hshared.elem[j1].Field == iHx) {
                        AniMed.info[jmed].Hx_Shared[i1].times = sgg.Hshared.elem[j1].Times;
                        if (sgg.Hshared.elem[j1].Times > 1) {
                            AniMed.info[jmed].Hx_Shared[i1].SharedMed.resize(sgg.Hshared.elem[j1].Times + 1);
                        }
                        found = true;
                    }
                }
            }
            if (!found) {
                for (int jmed = 1; jmed <= AniMed.NumMed && !found; ++jmed) {
                    for (int i1 = 1; i1 <= AniMed.info[jmed].numnodesHy; ++i1) {
                        if (sgg.Hshared.elem[j1].i == AniMed.info[jmed].Hy_i[i1] &&
                            sgg.Hshared.elem[j1].j == AniMed.info[jmed].Hy_j[i1] &&
                            sgg.Hshared.elem[j1].k == AniMed.info[jmed].Hy_k[i1] &&
                            sgg.Hshared.elem[j1].Field == iHy) {
                            AniMed.info[jmed].Hy_Shared[i1].times = sgg.Hshared.elem[j1].Times;
                            if (sgg.Hshared.elem[j1].Times > 1) {
                                AniMed.info[jmed].Hy_Shared[i1].SharedMed.resize(sgg.Hshared.elem[j1].Times + 1);
                            }
                            found = true;
                        }
                    }
                }
            }
            if (!found) {
                for (int jmed = 1; jmed <= AniMed.NumMed && !found; ++jmed) {
                    for (int i1 = 1; i1 <= AniMed.info[jmed].numnodesHz; ++i1) {
                        if (sgg.Hshared.elem[j1].i == AniMed.info[jmed].Hz_i[i1] &&
                            sgg.Hshared.elem[j1].j == AniMed.info[jmed].Hz_j[i1] &&
                            sgg.Hshared.elem[j1].k == AniMed.info[jmed].Hz_k[i1] &&
                            sgg.Hshared.elem[j1].Field == iHz) {
                            AniMed.info[jmed].Hz_Shared[i1].times = sgg.Hshared.elem[j1].Times;
                            if (sgg.Hshared.elem[j1].Times > 1) {
                                AniMed.info[jmed].Hz_Shared[i1].SharedMed.resize(sgg.Hshared.elem[j1].Times + 1);
                            }
                            found = true;
                        }
                    }
                }
            }
        }

        // Store shared media indices
        for (int jmed = 1; jmed <= AniMed.NumMed; ++jmed) {
            for (int i1 = 1; i1 <= AniMed.info[jmed].numnodesEx; ++i1) {
                int conta = 0;
                for (int j1 = 1; j1 <= sgg.Eshared.conta; ++j1) {
                    if (sgg.Eshared.elem[j1].i == AniMed.info[jmed].Ex_i[i1] &&
                        sgg.Eshared.elem[j1].j == AniMed.info[jmed].Ex_j[i1] &&