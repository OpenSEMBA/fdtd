#include <vector>
#include <memory>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>

// Assuming these headers exist based on the Fortran 'use' statements
// #include "Report_m.h"
// #include "FDETYPES_m.h"
// #include "wiresHolland_constants_m.h"
// #ifdef CompileWithStochastic
// #include "wiresHolland_devia.h"
// #endif
// #ifdef CompileWithThickWires
// #include "Thick_m.h"
// #endif

// Placeholder types to match Fortran derived types
// In a real translation, these would be defined in their respective header files

struct Thinwires_t {
    std::vector<int> WireTipoMedio;
    int NumDifferentWires = 0;
    
    // Null Node and Segment placeholders
    struct ChargeNodes_t {
        int indexnode = -1;
        double ChargePresent = 0.0;
        double ChargePast = 0.0;
        bool IsAttachedtoVoltage = false;
        bool IsMur = false;
        bool IsBackDownLeftMur = false;
        bool IsFrontUpRightMur = false;
        bool IsPeriodic = false;
        bool IsPEC = false;
        bool IsLossy = false;
        bool HasIsource = false;
        bool IsHeterogeneousJunction = false;
        bool Exists = false;
        bool proc = false;
        bool Is_LeftEnd = false;
        bool Is_RightEnd = false;
        bool IsInSingleRLCsegment = false;
        int NumCurrentMinus = 0;
        int NumCurrentPlus = 0;
        int i = -1;
        int j = -1;
        int k = -1;
        double cteMur = 0.0;
        double cteProp = 0.0;
        double oRIGctePlain = 0.0;
        double ctePlain = 0.0;
        // Pointers to fields would be handled differently in C++, simplified here
    } NullNode;

    struct CurrentSegments_t {
        double R = 0.0;
        double Resist = 0.0;
        double Resist_devia = 0.0;
        double C = 0.0;
        double L = 0.0;
        double Lintrinsic = 0.0;
        int NumParallel = 1;
        int origindex = 0;
        int indexsegment = 0;
        double currentpast = 0.0;
        double current = 0.0;
        double qplus_qminus = 0.0;
        double current_for_devia = 0.0;
        double qplus_qminus_for_devia = 0.0;
        double Efield_main2wire_for_devia = 0.0;
        double inv_Lind_acum = 0.0;
        double Lind_acum = 0.0;
        double Lind = 0.0;
        double Lind_devia = 0.0;
        double HEUR_safety = 0.0;
        double delta = 0.0;
        double deltaTransv1 = 0.0;
        double deltaTransv2 = 0.0;
        double cte1 = 0.0;
        double cte2 = 0.0;
        double cte3 = 0.0;
        double cte1_for_devia = 0.0;
        double cte2_for_devia = 0.0;
        double cte3_for_devia = 0.0;
        double cte5 = 0.0;
        int ilibre = -1;
        int jlibre = -1;
        int klibre = -1;
        int i = -1;
        int j = -1;
        int k = -1;
        int tipofield = -1;
        bool IsPMC = false;
        bool orientadoalreves = false;
        bool HasVsource = false;
        bool IsShielded = false;
        bool HasAbsorbing_RightEnd = false;
        bool HasAbsorbing_LeftEnd = false;
        bool HasParallel_RightEnd = false;
        bool HasParallel_LeftEnd = false;
        bool HasSeries_RightEnd = false;
        bool HasSeries_LeftEnd = false;
        bool IsEnd_norLeft_norRight = false;
        bool Is_LeftEnd = false;
        bool Is_RightEnd = false;
        // Pointers to ChargeNodes_t
        ChargeNodes_t* chargePlus = nullptr;
        ChargeNodes_t* chargeMinus = nullptr;
    } NullSegment;

    double olddt = 0.0;
    int NumNeededCurrentUpMPI = 0;
    int NumNeededCurrentDownMPI = 0;
};

struct limit_t {
    int XI, XE, YI, YE, ZI, ZE;
};

struct SGGFDTDINFO_t {
    double dt = 0.0;
    int NumMedia = 0;
    struct {
        int XI, XE, YI, YE, ZI, ZE;
    } Alloc[10]; // Assuming indices like iEx, iEy etc map to integers
    
    struct {
        int ZI, ZE;
    } Sweep[10];

    struct {
        double Epr = 1.0;
        double Mur = 1.0;
        struct {
            bool Is = {false};
            struct {
                bool ThinWire = false;
                int numsegmentos = 0;
                struct {
                    bool multirabo = false;
                    int i, j, k, ori, origIndex;
                    bool Is_LeftEnd = false;
                    bool Is_RightEnd = false;
                    bool IsEnd_norLeft_norRight = false;
                    int ilibre = -1;
                    int jlibre = -1;
                    int klibre = -1;
                    // Dispersive elements would be here
                } segm[100]; // Placeholder size
                struct {
                    bool disp = false;
                    bool disp_RightEnd = false;
                    bool disp_LeftEnd = false;
                } wire[1];
                bool HasAbsorbing_RightEnd = false;
                bool HasAbsorbing_LeftEnd = false;
                bool HasParallel_RightEnd = false;
                bool HasParallel_LeftEnd = false;
                bool HasSeries_RightEnd = false;
                bool HasSeries_LeftEnd = false;
            } wire[1];
        } Med[100]; // Placeholder size
    } sgg;
};

struct sim_control_t {
    int layoutnumber = 0;
    int num_procs = 1;
    bool strictOLD = false;
    bool connectendings = false;
};

// Global variables from the module
namespace HollandWires_m {

    const double HEUR_RADIUSOVERDELTA = 10.0;

    bool thereAreVsources = false;
    bool thereAreIsources = false;
    bool thereAreMurConditions = false;

    Thinwires_t HWires;
    
    // Pointers converted to vectors for memory management, 
    // though in C++ we might use raw pointers or unique_ptr depending on ownership.
    // Here we use std::vector to simulate the dynamic allocation from Fortran.
    std::vector<double> InvEps;
    std::vector<double> InvMu;
    std::vector<double> OldInvEps;
    std::vector<double> OldInvMu;

    double eps0 = 0.0;
    double mu0 = 0.0;

    // Helper function to simulate WarnErrReport
    void WarnErrReport(const std::string& msg, bool is_error = false) {
        std::cout << "[HollandWires] " << msg << std::endl;
        if (is_error) {
            std::cerr << "Error: " << msg << std::endl;
            // In Fortran this might stop execution, here we just print
        }
    }

    void InitWires(
        SGGFDTDINFO_t& sgg,
        const std::vector<int>& sggMiNo,
        std::vector<int>& sggMiEx,
        std::vector<int>& sggMiEy,
        std::vector<int>& sggMiEz,
        std::vector<int>& sggMiHx,
        std::vector<int>& sggMiHy,
        std::vector<int>& sggMiHz,
        std::vector<int>& Idxe,
        std::vector<int>& Idye,
        std::vector<int>& Idze,
        std::vector<int>& Idxh,
        std::vector<int>& Idyh,
        std::vector<int>& Idzh,
        const std::vector<double>& G2,
        const std::vector<limit_t>& SINPML_fullsize,
        const std::vector<limit_t>& fullsize,
        double& dtcritico,
        double eps00,
        double mu00,
        const sim_control_t& control
    ) {
        // Dummy variables
        double eps000 = eps00;
        double mu000 = mu00;

        // Initialize globals
        eps0 = eps00;
        mu0 = mu00;
        HWires.olddt = sgg.dt;

        // Reset pointers (vectors)
        HWires.WireTipoMedio.clear();
        // CurrentSegment, ChargeNode are internal structures, no need to clear vectors unless they hold dynamic data
        
        dtcritico = sgg.dt;

        // Direction strings (not strictly needed for logic but part of original)
        // std::string dir[10]; 
        // dir[iEx] = " X "; ...

        thereAreVsources = false;
        thereAreIsources = false;
        thereAreMurConditions = false;

        // MPI Barrier placeholder
        // #ifdef CompileWithMPI
        // MPI_Barrier(...);
        // #endif

        // Initialize Adjacency placeholder
        // adj.YESsegment(1:2) = -1;

        // Allocate InvEps and InvMu
        int num_media = sgg.NumMedia;
        InvEps.resize(num_media + 1);
        InvMu.resize(num_media + 1);
        
        for (int m = 0; m <= num_media; ++m) {
            InvEps[m] = 1.0 / (eps0 * sgg.sgg.Med[m].Epr);
            InvMu[m] = 1.0 / (mu0 * sgg.sgg.Med[m].Mur);
        }

        OldInvEps = InvEps;
        OldInvMu = InvMu;

        // MPI counters
        HWires.NumNeededCurrentUpMPI = 0;
        HWires.NumNeededCurrentDownMPI = 0;

        // Initialize NullNode
        HWires.NullNode.indexnode = -1;
        HWires.NullNode.ChargePresent = 0.0;
        HWires.NullNode.ChargePast = 0.0;
        HWires.NullNode.IsAttachedtoVoltage = false;
        HWires.NullNode.IsMur = false;
        HWires.NullNode.IsBackDownLeftMur = false;
        HWires.NullNode.IsFrontUpRightMur = false;
        HWires.NullNode.IsPeriodic = false;
        HWires.NullNode.IsPEC = false;
        HWires.NullNode.IsLossy = false;
        HWires.NullNode.HasIsource = false;
        HWires.NullNode.IsHeterogeneousJunction = false;
        HWires.NullNode.Exists = false;
        HWires.NullNode.proc = false;
        HWires.NullNode.Is_LeftEnd = false;
        HWires.NullNode.Is_RightEnd = false;
        HWires.NullNode.IsInSingleRLCsegment = false;
        HWires.NullNode.NumCurrentMinus = 0;
        HWires.NullNode.NumCurrentPlus = 0;
        HWires.NullNode.i = -1;
        HWires.NullNode.j = -1;
        HWires.NullNode.k = -1;
        HWires.NullNode.cteMur = 0.0;
        HWires.NullNode.cteProp = 0.0;
        HWires.NullNode.oRIGctePlain = 0.0;
        HWires.NullNode.ctePlain = 0.0;

        // Initialize NullSegment
        HWires.NullSegment.R = 0.0;
        HWires.NullSegment.Resist = 0.0;
        HWires.NullSegment.Resist_devia = 0.0;
        HWires.NullSegment.C = 0.0;
        HWires.NullSegment.L = 0.0;
        HWires.NullSegment.Lintrinsic = 0.0;
        HWires.NullSegment.NumParallel = 1;
        HWires.NullSegment.origindex = 0;
        HWires.NullSegment.indexsegment = 0;
        HWires.NullSegment.currentpast = 0.0;
        HWires.NullSegment.current = 0.0;
        HWires.NullSegment.qplus_qminus = 0.0;
        HWires.NullSegment.current_for_devia = 0.0;
        HWires.NullSegment.qplus_qminus_for_devia = 0.0;
        HWires.NullSegment.Efield_main2wire_for_devia = 0.0;
        HWires.NullSegment.inv_Lind_acum = 0.0;
        HWires.NullSegment.Lind_acum = 0.0;
        HWires.NullSegment.Lind = 0.0;
        HWires.NullSegment.Lind_devia = 0.0;
        HWires.NullSegment.HEUR_safety = 0.0;
        HWires.NullSegment.delta = 0.0;
        HWires.NullSegment.deltaTransv1 = 0.0;
        HWires.NullSegment.deltaTransv2 = 0.0;
        HWires.NullSegment.cte1 = 0.0;
        HWires.NullSegment.cte2 = 0.0;
        HWires.NullSegment.cte3 = 0.0;
        HWires.NullSegment.cte1_for_devia = 0.0;
        HWires.NullSegment.cte2_for_devia = 0.0;
        HWires.NullSegment.cte3_for_devia = 0.0;
        HWires.NullSegment.cte5 = 0.0;
        HWires.NullSegment.ilibre = -1;
        HWires.NullSegment.jlibre = -1;
        HWires.NullSegment.klibre = -1;
        HWires.NullSegment.i = -1;
        HWires.NullSegment.j = -1;
        HWires.NullSegment.k = -1;
        HWires.NullSegment.tipofield = -1;
        HWires.NullSegment.IsPMC = false;
        HWires.NullSegment.orientadoalreves = false;
        HWires.NullSegment.HasVsource = false;
        HWires.NullSegment.IsShielded = false;
        HWires.NullSegment.HasAbsorbing_RightEnd = false;
        HWires.NullSegment.HasAbsorbing_LeftEnd = false;
        HWires.NullSegment.HasParallel_RightEnd = false;
        HWires.NullSegment.HasParallel_LeftEnd = false;
        HWires.NullSegment.HasSeries_RightEnd = false;
        HWires.NullSegment.HasSeries_LeftEnd = false;
        HWires.NullSegment.IsEnd_norLeft_norRight = false;
        HWires.NullSegment.Is_LeftEnd = false;
        HWires.NullSegment.Is_RightEnd = false;
        HWires.NullSegment.chargePlus = &HWires.NullNode;
        HWires.NullSegment.chargeMinus = &HWires.NullNode;

        bool ThereAreWires = false;

        // Detect thin wires
        int conta = 0;
        for (int jmed = 1; jmed <= sgg.NumMedia; ++jmed) {
            if (sgg.sgg.Med[jmed].wire[1].Is.ThinWire) {
                conta++;
            }
        }

        HWires.NumDifferentWires = conta;
        HWires.WireTipoMedio.resize(conta);
        
        conta = 0;
        for (int jmed = 1; jmed <= sgg.NumMedia; ++jmed) {
            if (sgg.sgg.Med[jmed].wire[1].Is.ThinWire) {
                ThereAreWires = true;
                conta++;
                HWires.WireTipoMedio[conta - 1] = jmed;
            }
        }

        // Check for dispersive wires
        for (int iwi = 0; iwi < HWires.NumDifferentWires; ++iwi) {
            int med_idx = HWires.WireTipoMedio[iwi];
            if (sgg.sgg.Med[med_idx].wire[1].disp || 
                sgg.sgg.Med[med_idx].wire[1].disp_RightEnd || 
                sgg.sgg.Med[med_idx].wire[1].disp_LeftEnd) {
                WarnErrReport("Dispersive wire or connectors unsupported in Holland wires", true);
            }
        }

        // Reset multirabo
        for (int iwi = 0; iwi < HWires.NumDifferentWires; ++iwi) {
            int med_idx = HWires.WireTipoMedio[iwi];
            for (int iwj = 0; iwj < sgg.sgg.Med[med_idx].wire[1].numsegmentos; ++iwj) {
                sgg.sgg.Med[med_idx].wire[1].segm[iwj].multirabo = false;
            }
        }

        if (!ThereAreWires) return;

        WarnErrReport("----------------------------------------------------------------");

        // Process segments
        for (int iwi = 0; iwi < HWires.NumDifferentWires; ++iwi) {
            int med_idx = HWires.WireTipoMedio[iwi];
            for (int iwj = 0; iwj < sgg.sgg.Med[med_idx].wire[1].numsegmentos; ++iwj) {
                if (!control.strictOLD) {
                    sgg.sgg.Med[med_idx].wire[1].segm[iwj].ilibre = -1;
                    sgg.sgg.Med[med_idx].wire[1].segm[iwj].jlibre = -1;
                    sgg.sgg.Med[med_idx].wire[1].segm[iwj].klibre = -1;

                    if (control.connectendings) {
                        if (sgg.sgg.Med[med_idx].wire[1].numsegmentos == 1) {
                            sgg.sgg.Med[med_idx].wire[1].segm[iwj].Is_LeftEnd = true;
                            sgg.sgg.Med[med_idx].wire[1].segm[iwj].Is_RightEnd = true;
                            
                            if (sgg.sgg.Med[med_idx].wire[1].HasAbsorbing_RightEnd) 
                                sgg.sgg.Med[med_idx].wire[1].HasAbsorbing_LeftEnd = true;
                            if (sgg.sgg.Med[med_idx].wire[1].HasAbsorbing_LeftEnd) 
                                sgg.sgg.Med[med_idx].wire[1].HasAbsorbing_RightEnd = true;
                            if (sgg.sgg.Med[med_idx].wire[1].HasParallel_RightEnd) 
                                sgg.sgg.Med[med_idx].wire[1].HasParallel_LeftEnd = true;
                            if (sgg.sgg.Med[med_idx].wire[1].HasParallel_LeftEnd) 
                                sgg.sgg.Med[med_idx].wire[1].HasParallel_RightEnd = true;
                            if (sgg.sgg.Med[med_idx].wire[1].HasSeries_RightEnd) 
                                sgg.sgg.Med[med_idx].wire[1].HasSeries_LeftEnd = true;
                            if (sgg.sgg.Med[med_idx].wire[1].HasSeries_LeftEnd) 
                                sgg.sgg.Med[med_idx].wire[1].HasSeries_RightEnd = true;
                        }

                        sgg.sgg.Med[med_idx].wire[1].segm[iwj].Is_LeftEnd = 
                            sgg.sgg.Med[med_idx].wire[1].segm[iwj].Is_LeftEnd && 
                            (sgg.sgg.Med[med_idx].wire[1].HasParallel_LeftEnd || 
                             sgg.sgg.Med[med_idx].wire[1].HasSeries_LeftEnd || 
                             sgg.sgg.Med[med_idx].wire[1].HasAbsorbing_LeftEnd);
                             
                        sgg.sgg.Med[med_idx].wire[1].segm[iwj].Is_RightEnd = 
                            sgg.sgg.Med[med_idx].wire[1].segm[iwj].Is_RightEnd && 
                            (sgg.sgg.Med[med_idx].wire[1].HasParallel_RightEnd || 
                             sgg.sgg.Med[med_idx].wire[1].HasSeries_RightEnd || 
                             sgg.sgg.Med[med_idx].wire[1].HasAbsorbing_RightEnd);
                    }
                    
                    sgg.sgg.Med[med_idx].wire[1].segm[iwj].IsEnd_norLeft_norRight = 
                        !(sgg.sgg.Med[med_idx].wire[1].segm[iwj].Is_LeftEnd || 
                          sgg.sgg.Med[med_idx].wire[1].segm[iwj].Is_RightEnd);
                }
            }
        }
    }

    // Other functions would be implemented similarly
    void AdvanceWiresE(...) {}
    void AdvanceWiresH(...) {}
    void AdvanceWiresEcrank(...) {}
    void StoreFieldsWires(...) {}
    void DestroyWires(...) {}
    void GetHwires(...) {}
    void ReportWireJunctions(...) {}
    void calc_wirehollandconstants(...) {}

}

conectado1 = conectado1 || ((i1 == i2) && (j1 == j2 + 1) && (k1 == k2));
                    //
                    conectado2 = conectado2 || ((i1 + 1 == i2) && (j1 == j2) && (k1 == k2));
                    conectado2 = conectado2 || ((i1 + 1 == i2) && (j1 == j2 + 1) && (k1 == k2));
                    case (iEz):
                    conectado1 = conectado1 || ((i1 == i2) && (j1 == j2) && (k1 == k2));
                    conectado1 = conectado1 || ((i1 == i2) && (j1 == j2) && (k1 == k2 + 1));
                    //
                    conectado2 = conectado2 || ((i1 + 1 == i2) && (j1 == j2) && (k1 == k2));
                    conectado2 = conectado2 || ((i1 + 1 == i2) && (j1 == j2) && (k1 == k2 + 1));
                    break;
                    }
                    case (iEy):
                    switch (whatfield2) {
                    case (iEy):
                    conectado1 = conectado1 || ((j1 == j2 + 1) && (k1 == k2) && (i1 == i2));
                    conectado2 = conectado2 || ((j1 + 1 == j2) && (k1 == k2) && (i1 == i2));
                    break;
                    case (iEz):
                    conectado1 = conectado1 || ((j1 == j2) && (k1 == k2) && (i1 == i2));
                    conectado1 = conectado1 || ((j1 == j2) && (k1 == k2 + 1) && (i1 == i2));
                    //
                    conectado2 = conectado2 || ((j1 + 1 == j2) && (k1 == k2) && (i1 == i2));
                    conectado2 = conectado2 || ((j1 + 1 == j2) && (k1 == k2 + 1) && (i1 == i2));
                    break;
                    case (iEx):
                    conectado1 = conectado1 || ((j1 == j2) && (k1 == k2) && (i1 == i2));
                    conectado1 = conectado1 || ((j1 == j2) && (k1 == k2) && (i1 == i2 + 1));
                    //
                    conectado2 = conectado2 || ((j1 + 1 == j2) && (k1 == k2) && (i1 == i2));
                    conectado2 = conectado2 || ((j1 + 1 == j2) && (k1 == k2) && (i1 == i2 + 1));
                    break;
                    }
                    break;
                    case (iEz):
                    switch (whatfield2) {
                    case (iEz):
                    conectado1 = conectado1 || ((k1 == k2 + 1) && (i1 == i2) && (j1 == j2));
                    conectado2 = conectado2 || ((k1 + 1 == k2) && (i1 == i2) && (j1 == j2));
                    break;
                    case (iEx):
                    conectado1 = conectado1 || ((k1 == k2) && (i1 == i2) && (j1 == j2));
                    conectado1 = conectado1 || ((k1 == k2) && (i1 == i2 + 1) && (j1 == j2));
                    //
                    conectado2 = conectado2 || ((k1 + 1 == k2) && (i1 == i2) && (j1 == j2));
                    conectado2 = conectado2 || ((k1 + 1 == k2) && (i1 == i2 + 1) && (j1 == j2));
                    break;
                    case (iEy):
                    conectado1 = conectado1 || ((k1 == k2) && (i1 == i2) && (j1 == j2));
                    conectado1 = conectado1 || ((k1 == k2) && (i1 == i2) && (j1 == j2 + 1));
                    //
                    conectado2 = conectado2 || ((k1 + 1 == k2) && (i1 == i2) && (j1 == j2));
                    conectado2 = conectado2 || ((k1 + 1 == k2) && (i1 == i2) && (j1 == j2 + 1));
                    break;
                    }
                    break;
                    }

                    conectado = conectado1 && conectado2;
                    if (conectado) {
                        exit buskk;
                    }
                    }
                    } while (buskk);

                    if (control.connectendings) {
                        if ((!sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_LeftEnd) &&
                            (!sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_RightEnd)) {
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].IsEnd_norLeft_norRight = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].IsEnd_norLeft_norRight ||
                                ((!conectado) && conectado2) &&
                                (!(sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].HasParallel_LeftEnd || sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].HasSeries_LeftEnd ||
                                    sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].HasAbsorbing_LeftEnd));
                        }
                        if ((!sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_LeftEnd) &&
                            (!sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_RightEnd)) {
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].IsEnd_norLeft_norRight = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].IsEnd_norLeft_norRight ||
                                ((!conectado) && conectado1) &&
                                (!(sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].HasParallel_RightEnd || sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].HasSeries_RightEnd ||
                                    sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].HasAbsorbing_RightEnd));
                        }
                    }
                    sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].IsEnd_norLeft_norRight =
                        (sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].IsEnd_norLeft_norRight) &&
                        (!conectado) && (conectado1 || conectado2) &&
                        (!(sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_LeftEnd ||
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_RightEnd)); // si hay mas de uno este se pone a .true.
                    // detecta cual es el extremo libre
                    if ((!conectado) && conectado1) {
                        switch (whatfield) {
                        case (iEx):
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ilibre = i1 + 1;
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].jlibre = j1;
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].klibre = k1;
                            break;
                        case (iEy):
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ilibre = i1;
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].jlibre = j1 + 1;
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].klibre = k1;
                            break;
                        case (iEz):
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ilibre = i1;
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].jlibre = j1;
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].klibre = k1 + 1;
                            break;
                        }
                    }
                    else if ((!conectado) && conectado2) {
                        sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ilibre = i1;
                        sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].jlibre = j1;
                        sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].klibre = k1;
                    }

                    // caso especial
                    if (sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].numsegmentos == 1) {
                        sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].IsEnd_norLeft_norRight = false;
                        sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_LeftEnd = true;
                        sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_RightEnd = true;
                        switch (whatfield) {
                        case (iEx):
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ilibre = i1 + 1;
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].jlibre = j1;
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].klibre = k1;
                            break;
                        case (iEy):
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ilibre = i1;
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].jlibre = j1 + 1;
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].klibre = k1;
                            break;
                        case (iEz):
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ilibre = i1;
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].jlibre = j1;
                            sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].klibre = k1 + 1;
                            break;
                        }
                    }
                    // check for intermediate RLC error

                    if ((sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ilibre == -1) ||
                        (sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].jlibre == -1) ||
                        (sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].klibre == -1)) {
                        if (!conectado) {
                            sprintf(buff, "wir0_BUGGYERROR: Non-Intermediate multi-segment WIRE. %7d%7d%7d%7d %s", origIndex, i1, j1, k1, dir(whatfield).c_str());
                            if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff, true);
                        }
                        if (((sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_RightEnd) &&
                            (sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].HasParallel_RightEnd || sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].HasSeries_RightEnd)) ||
                            ((sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_LeftEnd) &&
                            (sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].HasParallel_LeftEnd || sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].HasSeries_LeftEnd ||
                                sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].HasAbsorbing_LeftEnd))) {
                            if (conectado) {
                                sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_RightEnd = false;
                                sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_LeftEnd = false;
                                sprintf(buff, "wir0_WARNING: Intermediate segment with RLC. Neglecting RLC %7d%7d%7d%7d %s", origIndex, i1, j1, k1, dir(whatfield).c_str());
                                if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff);
                            }
                            else {
                                sprintf(buff, "wir0_BUGGYERROR: Non-Intermediate multi-segment WIRE with RLC. %7d%7d%7d%7d %s", origIndex, i1, j1, k1, dir(whatfield).c_str());
                                if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff, true);
                            }
                        }
                    }
                    //
                    }
                    else { // del strictOLD

                        sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ilibre = -1;
                        sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].jlibre = -1;
                        sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].klibre = -1;
                        sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].IsEnd_norLeft_norRight = false; // irrelevante en strictOLD
                        i1 = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].i;
                        j1 = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].j;
                        k1 = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].k;
                        whatfield = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ori;
                        origindex = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].origIndex;
                        //
                        if (sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_LeftEnd) {
                            dummy1 = 1;
                            dummyfin = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].numsegmentos;
                        }
                        else if (sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_RightEnd) {
                            dummy1 = -1;
                            dummyfin = 1;
                        }
                        if ((sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_LeftEnd) || (sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_RightEnd)) {
                            if (sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].numsegmentos != 1) {
                                dummy2 = -1;
                                buscakk2: for (iwjjj = iwj + dummy1; dummy1 > 0 ? iwjjj <= dummyfin : iwjjj >= dummyfin; iwjjj += dummy1) { // atras o adelante
                                    i2 = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].i;
                                    j2 = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].j;
                                    k2 = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].k;
                                    whatfield2 = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].ori;
                                    if ((i1 == i2) && (j1 == j2) && (k1 == k2) && (whatfield == whatfield2)) {
                                        dummy2 = -dummy2; // detecta numero de rabitos para o impar aunque yo luego en las uniones solo trato 2 rabitos como mucho
                                    }
                                    else {
                                        continue;
                                        break buscakk2;
                                    }
                                }
                                if (dummy2 == -1) {
                                    dummy3 = 0;
                                }
                                else {
                                    dummy3 = 1;
                                }

                                //
                                sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ilibre = i1;
                                sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].jlibre = j1;
                                sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].klibre = k1;
                                if (whatfield2 == whatfield) {
                                    if (abs(i1 - i2) + abs(j1 - j2) + abs(k1 - k2) > 1) {
                                        sprintf(buff, "wir0_ERROR: strictOLD LeftEnd/RightEnd segment disconnected. %7d%7d%7d%7d %s", origIndex, i1, j1, k1, dir(whatfield).c_str());
                                        if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff, true);
                                    }
                                    if (i1 > i2) {
                                        sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ilibre = i1 + 1 - dummy3;
                                    }
                                    else if (i1 < i2) {
                                        sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ilibre = i1 + dummy3;
                                    }
                                    if (j1 > j2) {
                                        sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].jlibre = j1 + 1 - dummy3;
                                    }
                                    else if (j1 < j2) {
                                        sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].jlibre = j1 + dummy3;
                                    }
                                    if (k1 > k2) {
                                        sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].klibre = k1 + 1 - dummy3;
                                    }
                                    else if (k1 < k2) {
                                        sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].klibre = k1 + dummy3;
                                    }
                                }
                                else {
                                    if (abs(i1 - i2) + abs(j1 - j2) + abs(k1 - k2) > 2) {
                                        sprintf(buff, "wir0_ERROR: strictOLD LeftEnd/RightEnd segment disconnected. %7d%7d%7d%7d %s", origIndex, i1, j1, k1, dir(whatfield).c_str());
                                        if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff, true);
                                    }
                                    switch (whatfield) {
                                    case (iEx):
                                        if (i2 == i1) sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ilibre = i1 + 1;
                                        break;
                                    case (iEy):
                                        if (j2 == j1) sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].jlibre = j1 + 1;
                                        break;
                                    case (iEz):
                                        if (k2 == k1) sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].klibre = k1 + 1;
                                        break;
                                    }
                                }
                            }
                            else { // DEL NUMERO SEGMENTOS '2014 NO PORTADO A !CHECK
                                switch (whatfield) {
                                case (iEx):
                                    sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ilibre = i1 + 1;
                                    sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].jlibre = j1;
                                    sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].klibre = k1;
                                    break;
                                case (iEy):
                                    sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ilibre = i1;
                                    sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].jlibre = j1 + 1;
                                    sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].klibre = k1;
                                    break;
                                case (iEz):
                                    sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ilibre = i1;
                                    sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].jlibre = j1;
                                    sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].klibre = k1 + 1;
                                    break;
                                }
                            } // DEL NUMERO SEGMENTOS '2014 NO PORTADO A !CHECK
                        } // Del Left_End Right_End
                        // !!!!!!!!!!!!!!!!!!!
                    } // del strictOLD

                    }
                    } while (1);
                    }
                    }

                    // preprocesa para eliminar multirabos luego se utiliza en repetido

                    if (control.strictOLD) {
                        for (iwi = 1; iwi <= HWires.NumDifferentWires; iwi++) {
                            for (iwj = 1; iwj <= sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].numsegmentos; iwj++) {
                                if (!sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].multirabo) {
                                    multirabos = 1;
                                    i1 = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].i;
                                    j1 = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].j;
                                    k1 = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].k;
                                    whatfield = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ori;
                                    ORIGINDEX = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].origindex;
                                    // precontaje
                                    buscarabos: for (iwjjj = iwj + 1; iwjjj <= sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].numsegmentos; iwjjj++) { // el Right_End aunque no se tape si debe detectarse a efectos par/impar
                                        i2 = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].i;
                                        j2 = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].j;
                                        k2 = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].k;
                                        whatfield2 = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].ori;
                                        ORIGINDEX2 = sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].origindex;
                                        if ((i1 == i2) && (j1 == j2) && (k1 == k2) && (whatfield == whatfield2)) {
                                            multirabos = multirabos + 1;
                                        }
                                        else {
                                            primernorabo = origindex2;
                                            Jprimernorabo = iwjjj;
                                            break buscarabos;
                                        }
                                    }
                                    // machaca rabos
                                    if (multirabos != 1) {
                                        //
                                        taparabos: for (iwjjj = iwj + (2 - multirabos % 2); iwjjj <= sgg->Med(HWires.WireTipoMedio(iwi)).wire[1].numsegmentos - 1; iwjjj++) { // el Right_End no debe taparse ya se tapa el de dentro en el otro buc

}
            le
                           i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%i
                           j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%j
                           k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%k
                           whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%ori
                           ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex
                           if ((i1==i2).and.(j1==j2).and.(k1==k2).and.(whatfield==whatfield2)) then
                              write (buff,'(a,i7,3I7,a,a,i7)')  'wir0_WARNING: strictOLD Redundannt zig-zag rabito, will be eliminated to Mod(2)', origIndex2,i2,j2,k2,dir(whatfield2), &
                              'by segment ',primernorabo
                              if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff)
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%multirabo=.true.
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%multiraboDE=primernorabo
                           else
                              exit taparabos
                           end if
                        end do taparabos
                     end if !del multirabos no nulo
                  end if
               end do
            end do
            !ahora la vuelta
            do iwi=1,HWires%NumDifferentWires
               do iwj=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos,1,-1
                  if (.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo) then
                     multirabos=1
                     i1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
                     j1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
                     k1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
                     whatfield=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
                     ORIGINDEX=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origindex
                     buscarabos2: do iwjjj=iwj-1,1,-1 !el Left_End aunque no se tape si debe detectarse a efectos par/impar
                        i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%i
                        j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%j
                        k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%k
                        whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%ori
                        ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex
                        if ((i1==i2).and.(j1==j2).and.(k1==k2).and.(whatfield==whatfield2)) then
                           multirabos=multirabos+1
                        else
                           primernorabo=origindex2
                           Jprimernorabo=iwjjj
                           exit buscarabos2
                        end if
                     end do buscarabos2

                     !machaca rabos
                     if (multirabos/=1) then
                        !
                        taparabos2: do iwjjj=iwj-(2-mod(multirabos,2)),2,-1 !el Left_End no debe taparse
                           i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%i
                           j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%j
                           k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%k
                           whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%ori
                           ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex
                           if ((i1==i2).and.(j1==j2).and.(k1==k2).and.(whatfield==whatfield2)) then
                              !faltaria eliminar sondas en multirabos 7/2/14
                              write (buff,'(a,i7,3I7,a,a,i7)')  'wir0_WARNING: strictOLD Redundannt zig-zag rabito, will be eliminated to Mod(2)', origIndex2,i2,j2,k2,dir(whatfield2), &
                              'by segment ',primernorabo
                              if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff)
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%multirabo=.true.
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%multiraboDE=primernorabo
                           else
                              exit taparabos2
                           end if
                        end do taparabos2
                     end if !del multirabos no nulo
                  end if
               end do
            end do


            !!!!!!!!!!!!!!!!
            !segunda pasada para procesar TAPARRABOS
            !!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!

            if (control%TAPARRABOS) then
               do iwi=1,HWires%NumDifferentWires
                  do iwj=1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
                     if (.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo) then
                        multirabos=1
                        i1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
                        j1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
                        k1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
                        whatfield=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
                        ORIGINDEX=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origindex
                        !precontaje
                        buscarabos6: do iwjjj=iwj+1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos !el Right_End aunque no se tape si debe detectarse a efectos par/impar
                           i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%i
                           j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%j
                           k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%k
                           whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%ori
                           ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex
                           if ((i1==i2).and.(j1==j2).and.(k1==k2).and.(whatfield==whatfield2)) then
                              multirabos=multirabos+1
                           else
                              primernorabo=origindex2
                              Jprimernorabo=IWjjj
                              exit buscarabos6
                           end if
                        end do buscarabos6
                        !machaca rabos
                        if ((mod(multirabos,2)/=1).and.(Jprimernorabo /= sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos)) then !no al RightEnd
                           if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj          )%Is_LeftEnd) then
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj          )%Is_LeftEnd=.false.
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%Is_LeftEnd=.true. !pasa caracter Left_End al primernorabo
                              !!!tocado esto por el problema de gra_powerline_simple.nfde 190916 el rabito se quedaba con el libre mal computado. Los ilibre,jlbre,klibre para Left_End son el primer punto directamente
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%ilibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%i  !!!!sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ilibre
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%jlibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%j  !!!!sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%jlibre
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%klibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%k  !!!!sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%klibre
!!!y retocado 200117 por el problema del ntc1 de la demo de getafe que tambien calculaba mal el libre (no tiene que ser por guevos el primer punto, habra que ver como se conecta con el siguiente!!!)
                                 if ((sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%i)==(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo+1)%i).and. &
                                     (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%j)==(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo+1)%j).and. &
                                     (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%k)==(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo+1)%k)) then
                                     select case (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%ori)
                                     case (1)
                                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%ilibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%ilibre+1
                                     case (2)
                                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%jlibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%jlibre+1
                                     case (3)
                                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%klibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%klibre+1
                                     end select
                                end if
                           end if
                           !
                           if ((.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_LeftEnd).AND.(.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_RightEnd)) then
                              i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
                              j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
                              k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
                              whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
                              ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origindex
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo=.true.
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multiraboDE=primernorabo !el ultimo cazado es el primer NOrabo
                              write (buff,'(a,i7,3I7,a,a,i7)')  'wir0_WARNING: strictOLD and taparrabos redundannt zig-zag rabito, will be ALSO eliminated ', origIndex2,i2,j2,k2,dir(whatfield2), &
                              'by segment ',primernorabo
                              if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff)
                           end if
                           !
                           if  ((.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%Is_LeftEnd).AND.(.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%Is_RightEnd)) then
                              i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%i
                              j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%j
                              k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%k
                              whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%ori
                              ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%origindex
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%multirabo=.true.
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%multiraboDE=primernorabo !el ultimo cazado es el primer NOrabo
                              write (buff,'(a,i7,3I7,a,a,i7)')  'wir0_WARNING: strictOLD and taparrabos redundannt zig-zag rabito, will be ALSO eliminated ', origIndex2,i2,j2,k2,dir(whatfield2), &
                              'by segment ',primernorabo
                              if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff)
                           end if
                        end if !del multirabos no nulo
                     end if
                  end do
               end do
               !ahora la vuelta
               do iwi=1,HWires%NumDifferentWires
                  do iwj=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos,1,-1
                     if (.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo) then
                        multirabos=1
                        i1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
                        j1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
                        k1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
                        whatfield=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
                        ORIGINDEX=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origindex
                        buscarabos7: do iwjjj=iwj-1,1,-1 !el Left_End aunque no se tape si debe detectarse a efectos par/impar
                           i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%i
                           j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%j
                           k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%k
                           whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%ori
                           ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex
                           if ((i1==i2).and.(j1==j2).and.(k1==k2).and.(whatfield==whatfield2)) then
                              multirabos=multirabos+1
                           else
                              primernorabo=origindex2
                              Jprimernorabo=iwjjj
                              exit buscarabos7
                           end if
                        end do buscarabos7

                        !machaca rabos
                        if ((mod(multirabos,2)/=1).and.(Jprimernorabo /= 1)) then   !no el Left_End
                           if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj          )%Is_RightEnd) then
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj          )%Is_RightEnd=.false.
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%Is_RightEnd=.true. !pasa caracter RightEnd al primernorabo
                              !!!tocado esto por el problema de gra_powerline_simple.nfde 190916 el rabito se quedaba con el libre mal computado. Los ilibre,jlbre,klibre para Left_End son el primer punto directamente
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%ilibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%i  !!!!sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ilibre
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%jlibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%j  !!!!sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%jlibre
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%klibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%k  !!!!sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%klibre
!!!y retocado 200117 por el problema del ntc1 de la demo de getafe que tambien calculaba mal el libre (no tiene que ser por guevos el primer punto, habra que ver como se conecta con el siguiente!!!)
                                 if ((sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%i)==(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo-1)%i).and. &
                                     (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%j)==(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo-1)%j).and. &
                                     (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%k)==(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo-1)%k)) then
                                     select case (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%ori)
                                     case (1)
                                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%ilibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%ilibre+1
                                     case (2)
                                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%jlibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%jlibre+1
                                     case (3)
                                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%klibre=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(Jprimernorabo)%klibre+1
                                     end select
                                end if
                           
                           end if
                           !
                           if ((.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_LeftEnd).AND.(.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_RightEnd)) then
                              i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
                              j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
                              k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
                              whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
                              ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origindex
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo=.true.
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multiraboDE=primernorabo !el ultimo cazado es el primer NOrabo
                              write (buff,'(a,i7,3I7,a,a,i7)')  'wir0_WARNING: strictOLD and taparrabos redundannt zig-zag rabito, will be ALSO eliminated ', origIndex2,i2,j2,k2,dir(whatfield2), &
                              'by segment ',primernorabo
                              if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff)
                           end if
                           !
                           if  ((.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%Is_LeftEnd).AND.(.NOT.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%Is_RightEnd)) then
                              i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%i
                              j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%j
                              k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%k
                              whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%ori
                              ORIGINDEX2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%origindex
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%multirabo=.true.
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%multiraboDE=primernorabo !el ultimo cazado es el primer NOrabo
                              write (buff,'(a,i7,3I7,a,a,i7)')  'wir0_WARNING: strictOLD and taparrabos redundannt zig-zag rabito, will be ALSO eliminated ', origIndex2,i2,j2,k2,dir(whatfield2), &
                              'by segment ',primernorabo
                              if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff)
                           end if
                        end if !del multirabos no nulo
                     end if
                  end do
               end do
            end if !control%TAPARRABOS

}  'wir0_WARNING: strictOLD and taparrabos redundannt zig-zag rabito, will be ALSO eliminated ', origIndex2,i2,j2,k2,dir(whatfield2), &
                              'by segment ',primernorabo
                              if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff)
                           end if
                        end if !del multirabos no nulo
                     end if
                  end do
               end do
            end if !del TAPARRABOS
            !
         end if !del if strictOLD

         !chequeo de errores buggy multirabo
         do iwi=1,HWires%NumDifferentWires
            do iwj=1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
               if ((sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%multirabo).and. &
               ((sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_LeftEnd).or.(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_RightEnd))) then
                  i1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
                  j1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
                  k1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
                  whatfield=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
                  ORIGINDEX=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origindex
                  write (buff,'(a,i7,3I7,a)')  'wir0_BuggyERROR: strictOLD LeftEnd/RightEnd cannot be multirabo ', origIndex,i1,j1,k1,' '//dir(whatfield)
                  if ((k1 >= ZI).and.(k1 <= ZE)) call WarnErrReport(buff,.true.)
               end if
            end do
         end do





         !!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !report them
         do iwi=1,HWires%NumDifferentWires
            do iwj=1,      sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
               i1=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
               j1=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
               k1=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
               i1libre=   sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ilibre
               j1libre=   sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%jlibre
               k1libre=   sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%klibre
               whatfield= sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
               origindex= sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origIndex
               LeftEnd_index=  sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%LeftEnd
               RightEnd_index=  sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%RightEnd
               !
               if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_LeftEnd.and.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_RightEnd) then
                  write (buff,'(a,4I7,a,3I7,a,2i7)')  'wir0_INFO: Ending segment (LeftEnd_and_Right)',origIndex,i1,j1,k1,'-', &
                  i1libre,j1libre,k1libre,' '//dir(whatfield),LeftEnd_index,RightEnd_index
                  if ((k1 >= ZI).and.(k1 <= ZE).and.control%verbose) call WarnErrReport(buff)
               elseif (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_LeftEnd) then
                  write (buff,'(a,4I7,a,3I7,a,i7)')  'wir0_INFO: Ending segment (LeftEnd   )',origIndex,i1,j1,k1,'-', &
                  i1libre,j1libre,k1libre,' '//dir(whatfield),LeftEnd_index
                  if ((k1 >= ZI).and.(k1 <= ZE).and.control%verbose) call WarnErrReport(buff)
               elseif (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_RightEnd) then
                  write (buff,'(a,4I7,a,3I7,a,i7)')  'wir0_INFO: Ending segment (RightEnd   )',origIndex,i1,j1,k1,'-', &
                  i1libre,j1libre,k1libre,' '//dir(whatfield),RightEnd_index
                  if ((k1 >= ZI).and.(k1 <= ZE).and.control%verbose) call WarnErrReport(buff)
               end if
               if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%IsEnd_norLeft_norRight) then
                  if (control%connectendings) then
                     write (buff,'(a,4I7,a,3I7,a)')  'wir0_INFO: Ending segment (other )',origIndex,i1,j1,k1,'-', &
                     i1libre,j1libre,k1libre,' '//dir(whatfield)
                     if ((k1 >= ZI).and.(k1 <= ZE).and.control%verbose) call WarnErrReport(buff)
                  else
                     sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%IsEnd_norLeft_norRight =.false.
                     write (buff,'(a,4I7,a,3I7,a)')  'wir0_WARNING: Resetting Ending segment (other ) to NON-ENDING',origIndex,i1,j1,k1,'-', &
                     i1libre,j1libre,k1libre,' '//dir(whatfield)
                     if ((k1 >  ZI).and.(k1 <= ZE).and.(whatfield /= iEz)) call WarnErrReport(buff)
                     if ((k1 >= ZI).and.(k1 <= ZE).and.(whatfield == iEz)) call WarnErrReport(buff)
                  end if
               end if
            end do
         end do

         !Segment pre-counting incluyendo deteccion de repetidos

         do iwi=1,HWires%NumDifferentWires
            do iwj=1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
               sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%repetido=.false.
            end do
         end do
         !
         conta=0
         do iwi=1,HWires%NumDifferentWires
            do iwj=1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
               if (((.not.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%repetido).or.control%strictOLD).and.(.not.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo)) then
                  i1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
                  j1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
                  k1=       sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
                  whatfield=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
                  !
                  do iwjjj=iwj+1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
                     i2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%i
                     j2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%j
                     k2=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%k
                     whatfield2=sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%ori
                     repetido = (i1 == i2).and.(j1 == j2).and.(k1 == k2).and.(whatfield == whatfield2)
                     if (repetido) then
                        if     (.not.(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%Is_LeftEnd .or. &
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%Is_RightEnd )) then

                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%repetido=repetido.or. &
                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%repetido
                        elseif (.not.(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_LeftEnd .or. &
                        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_RightEnd )) then

                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%repetido=repetido.or. &
                           sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%repetido
                        else
                           !aviso pero tomo una decision. md 260213 a veces lo duplica en principio y final!!!!!!

                           if ( ((abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_R_RightEnd) < 1.0e-12_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_R_RightEnd  ) < 1.0e-12_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_C_RightEnd) < 1.0e-12_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_C_RightEnd  ) > 1.0e7_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_L_RightEnd) < 1.0e-12_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_L_RightEnd  ) < 1.0e-12_RKIND_wires)).and. &
                           !
                           ((abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_R_LeftEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_R_LeftEnd  ) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_C_LeftEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_C_LeftEnd  ) <= 1.0e7_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_L_LeftEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_L_LeftEnd  ) >= 1.0e-12_RKIND_wires)) ) then
                              if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%Is_RightEnd) then
                                 sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%repetido=.true.
                                 if (control%strictOLD) then
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Keeping both', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 else
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Removing the second one (no RightEnd RLC)', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 end if
                                 if ((k1 >  ZI).and.(k1 <= ZE).and.(whatfield /= iEz).and.control%verbose) call WarnErrReport(buff)
                                 if ((k1 >= ZI).and.(k1 <= ZE).and.(whatfield == iEz).and.control%verbose) call WarnErrReport(buff)
                              elseif (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_RightEnd) then
                                 sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%repetido=.true.
                                 if (control%strictOLD) then
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Keeping both', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 else
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Removing the first one (no RightEnd RLC)', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 end if
                                 if ((k1 >  ZI).and.(k1 <= ZE).and.(whatfield /= iEz).and.control%verbose) call WarnErrReport(buff)
                                 if ((k1 >= ZI).and.(k1 <= ZE).and.(whatfield == iEz).and.control%verbose) call WarnErrReport(buff)
                              end if
                           end if

                           if ( ((abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_R_RightEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_R_RightEnd  ) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_C_RightEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_C_RightEnd  ) <= 1.0e7_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_L_RightEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_L_RightEnd  ) >= 1.0e-12_RKIND_wires)).and. &
                           !
                           ((abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_R_LeftEnd) < 1.0e-12_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_R_LeftEnd  ) < 1.0e-12_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_C_LeftEnd) < 1.0e-12_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_C_LeftEnd  ) > 1.0e7_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_L_LeftEnd) < 1.0e-12_RKIND_wires).and. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_L_LeftEnd  ) < 1.0e-12_RKIND_wires)) ) then
                              if (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%Is_LeftEnd) then
                                 sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%repetido=.true.
                                 if (control%strictOLD) then
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Keeping both', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 else
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Removing the second one (NO LeftEnd RLC, Only RightEnd RLC)', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 end if
                                 if ((k1 >  ZI).and.(k1 <= ZE).and.(whatfield /= iEz).and.control%verbose) call WarnErrReport(buff)
                                 if ((k1 >= ZI).and.(k1 <= ZE).and.(whatfield == iEz).and.control%verbose) call WarnErrReport(buff)
                              elseif (sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%Is_LeftEnd) then
                                 sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%repetido=.true.
                                 if (control%strictOLD) then
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Keeping both', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 else
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Removing the first one (NO LeftEnd RLC, Only RightEnd RLC)', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 end if
                                 if ((k1 >  ZI).and.(k1 <= ZE).and.(whatfield /= iEz).and.control%verbose) call WarnErrReport(buff)
                                 if ((k1 >= ZI).and.(k1 <= ZE).and.(whatfield == iEz).and.control%verbose) call WarnErrReport(buff)
                              end if
                           end if
                           !
                           if ( ((abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_R_LeftEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_R_LeftEnd  ) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_C_LeftEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_C_LeftEnd  ) <= 1.0e7_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_L_LeftEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_L_LeftEnd  ) >= 1.0e-12_RKIND_wires)).AND. &
                           !
                           ((abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_R_RightEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_R_RightEnd  ) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_C_RightEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_C_RightEnd  ) <= 1.0e7_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_L_RightEnd) >= 1.0e-12_RKIND_wires).or. &
                           (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_L_RightEnd  ) >= 1.0e-12_RKIND_wires)) ) then

                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%repetido=repetido.or. &
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%repetido

                              if (  (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_R_RightEnd -    &
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_R_LeftEnd) < 1.0e-12_RKIND_wires).and. &
                              (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_R_RightEnd -    &
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_R_LeftEnd  ) < 1.0e-12_RKIND_wires).and. &
                              (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_C_RightEnd -    &
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_C_LeftEnd) < 1.0e-12_RKIND_wires).and. &
                              (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_C_RightEnd -    &
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_C_LeftEnd  ) < 1.0e-12_RKIND_wires).and. &
                              (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_L_RightEnd -    &
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Parallel_L_LeftEnd) < 1.0e-12_RKIND_wires).and. &
                              (abs(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_L_RightEnd -    &
                              sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%Series_L_LeftEnd  ) < 1.0e-12_RKIND_wires)) then
                                 if (control%strictOLD) then
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment with the same RLC. Keeping both', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 else
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment with the same RLC. Removing the second one', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 end if
                                 if ((k1 >  ZI).and.(k1 <= ZE).and.(whatfield /= iEz).and.control%verbose) call WarnErrReport(buff)
                                 if ((k1 >= ZI).and.(k1 <= ZE).and.(whatfield == iEz).and.control%verbose) call WarnErrReport(buff)
                              end if
                           end if
                        end if
                     end if
                  end do
               end if
            end do
         end do

ff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment with the same RLC. Will remove the second one', &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                    sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 end if
                                 if ((k1 >  ZI).and.(k1 <= ZE).and.(whatfield /= iEz).and.control%verbose) call WarnErrReport(buff)
                                 if ((k1 >= ZI).and.(k1 <= ZE).and.(whatfield == iEz).and.control%verbose) call WarnErrReport(buff)
                              else
                                 write (buff,'(a,2i7,3i7)')  'wir0_ERROR: Duplicate terminal LeftEnd and RightEnd Parallel segment with non-null different RLC', &
                                 sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj  )%origindex, &
                                 sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwjjj)%origindex,i1,j1,k1
                                 if ((k1 >  ZI).and.(k1 <= ZE).and.(whatfield /= iEz)) call WarnErrReport(buff,.true.)
                                 if ((k1 >= ZI).and.(k1 <= ZE).and.(whatfield == iEz)) call WarnErrReport(buff,.true.)
                              end if
                           end if
                        end if
                     end if
                  end do
                  !
                  !clipping: in case of direct .nfde reading the PREPROCESSor has not clipped this data
                  if ((i1 >= sgg%Alloc(whatfield)%XI).and. &
                  (i1 <= sgg%Alloc(whatfield)%XE).and. &
                  (j1 >= sgg%Alloc(whatfield)%YI).and. &
                  (j1 <= sgg%Alloc(whatfield)%YE).and. &
                  (k1 >= sgg%Alloc(whatfield)%ZI).and. &
                  (k1 <= sgg%Alloc(whatfield)%ZE)) then
                     conta=conta+1
                  end if
               end if
            end do
         end do
         if (conta==0) therearewires=.false.
         if (.not.therearewires) return
         !Report duplicated segments
         do iwi=1,HWires%NumDifferentWires
            do iwj=1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
               i1=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
               j1=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
               k1=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
               whatfield =sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
               origindex= sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origIndex
               if ((sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%repetido).and.(.not.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo)) then
                  if (control%strictOLD) then
                     write (buff,'(a,4i7,a)')  'wir0_WARNING: Keeping duplicate (Parallel) intra-WIRE segment', &
                     origindex,i1,j1,k1,' '//dir(whatfield)
                  else
                     write (buff,'(a,4i7,a)')  'wir0_WARNING: Removing duplicate (Parallel) intra-WIRE segment and voiding ASSOCIATED probes ', &
                     origindex,i1,j1,k1,' '//dir(whatfield)
                  end if
                  if ((k1 >= ZI).and.(k1 <= ZE).and.control%verbose) call WarnErrReport(buff)
               end if
            end do
         end do

         HWires%NumCurrentSegments=conta
         !inicializa ctes segmentos
         allocate (HWires%CurrentSegment(1 : HWires%NumCurrentSegments))
         do i1=1,HWires%NumCurrentSegments
            nullify (HWires%CurrentSegment(i1)%ChargePlus,HWires%CurrentSegment(i1)%ChargeMinus,  &
            HWires%CurrentSegment(i1)%TipoWire, &
            HWires%CurrentSegment(i1)%Efield_main2wire,HWires%CurrentSegment(i1)%Efield_wire2main)
            !
            HWires%CurrentSegment(i1)%R      = 0.0_RKIND_wires
            HWires%CurrentSegment(i1)%Resist      = 0.0_RKIND_wires
            HWires%CurrentSegment(i1)%Resist_devia      = 0.0_RKIND_wires
            HWires%CurrentSegment(i1)%C      = 0.0_RKIND_wires
            HWires%CurrentSegment(i1)%L      = 0.0_RKIND_wires
            HWires%CurrentSegment(i1)%Lintrinsic    = 0.0_RKIND_wires
            HWires%CurrentSegment(i1)%proc   =.false.
            HWires%CurrentSegment(i1)%NumParallel    =1
            HWires%CurrentSegment(i1)%origindex       =i1
            HWires%CurrentSegment(i1)%indexsegment    =i1
            HWires%CurrentSegment(i1)%currentpast         =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%current         =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%inv_Lind_acum            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%Lind_acum            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%HEUR_safety =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%Lind            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%Lind_devia            =0.0_RKIND_wires
            !!! HWires%CurrentSegment(i1)%logRoverR0      =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%delta           =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%deltaTransv1    =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%deltaTransv2    =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%cte1            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%cte2            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%cte3            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%cte1_for_devia            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%cte2_for_devia            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%cte3_for_devia            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%cte5            =0.0_RKIND_wires
            HWires%CurrentSegment(i1)%ilibre               =-1
            HWires%CurrentSegment(i1)%jlibre               =-1
            HWires%CurrentSegment(i1)%klibre               =-1
            HWires%CurrentSegment(i1)%i               =-1
            HWires%CurrentSegment(i1)%j               =-1
            HWires%CurrentSegment(i1)%k               =-1
            HWires%CurrentSegment(i1)%tipofield       =-1
            HWires%CurrentSegment(i1)%IsPMC           =.false.
            HWires%CurrentSegment(i1)%orientadoalreves =.false.
            HWires%CurrentSegment(i1)%HasVsource      =.false.
            HWires%CurrentSegment(i1)%IsShielded      =.false.
            HWires%CurrentSegment(i1)%HasAbsorbing_RightEnd =.false.
            HWires%CurrentSegment(i1)%HasAbsorbing_LeftEnd =.false.
            HWires%CurrentSegment(i1)%HasParallel_RightEnd =.false.
            HWires%CurrentSegment(i1)%HasParallel_LeftEnd =.false.
            HWires%CurrentSegment(i1)%HasSeries_RightEnd   =.false.
            HWires%CurrentSegment(i1)%HasSeries_LeftEnd   =.false.
            HWires%CurrentSegment(i1)%IsEnd_norLeft_norRight   =.false.
            HWires%CurrentSegment(i1)%Is_LeftEnd   =.false.
            HWires%CurrentSegment(i1)%Is_RightEnd   =.false.
         end do

         !assign segment info
         conta=0
         do iwi=1,HWires%NumDifferentWires
            do iwj=1,sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos
               i1libre=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ilibre
               j1libre=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%jlibre
               k1libre=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%klibre
               i1=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%i
               j1=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%j
               k1=        sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%k
               whatfield= sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%ori
               origindex= sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origIndex
               IsEnd_norLeft_norRight=     sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%IsEnd_norLeft_norRight
               Is_LeftEnd=     sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_LeftEnd
               Is_RightEnd=     sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%Is_RightEnd
               if (((.not.(sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%repetido)).or.control%strictOLD).and.(.not.sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo)) then
                  !clipping: in case of direct .nfde reading the PREPROCESSor has not clipped this data
                  if ((i1 >= sgg%Alloc(whatfield)%XI).and. &
                  (i1 <= sgg%Alloc(whatfield)%XE).and. &
                  (j1 >= sgg%Alloc(whatfield)%YI).and. &
                  (j1 <= sgg%Alloc(whatfield)%YE).and. &
                  (k1 >= sgg%Alloc(whatfield)%ZI).and. &
                  (k1 <= sgg%Alloc(whatfield)%ZE)) then
                     conta=conta+1
                     HWires%CurrentSegment(conta)%IsEnd_norLeft_norRight=IsEnd_norLeft_norRight
                     HWires%CurrentSegment(conta)%Is_LeftEnd =Is_LeftEnd
                     HWires%CurrentSegment(conta)%Is_RightEnd =Is_RightEnd
                     HWires%CurrentSegment(conta)%origindex=origindex
                     HWires%CurrentSegment(conta)%tipofield=whatfield
                     HWires%CurrentSegment(conta)%ilibre=i1libre
                     HWires%CurrentSegment(conta)%jlibre=j1libre
                     HWires%CurrentSegment(conta)%klibre=k1libre
                     HWires%CurrentSegment(conta)%i=i1
                     HWires%CurrentSegment(conta)%j=j1
                     HWires%CurrentSegment(conta)%k=k1
                     HWires%CurrentSegment(conta)%ie=i1
                     HWires%CurrentSegment(conta)%je=j1
                     HWires%CurrentSegment(conta)%ke=k1
                     HWires%CurrentSegment(conta)%indexmed=HWires%WireTipoMedio(iwi)
                     HWires%CurrentSegment(conta)%TipoWire=>sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)
                     !!only for the observation sign to match (not used in this routine)
                     if (iwj < sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%numsegmentos) then
                        if (.not.control%strictOLD) then
                           HWires%CurrentSegment(conta)%orientadoalreves = &
                           (i1 > sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%i).or. &
                           (j1 > sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%j).or. &
                           (k1 > sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj+1)%k)
                        else
                           HWires%CurrentSegment(conta)%orientadoalreves =.false. !later corrected
                        end if
                     elseif     (iwj > 1 ) then
                        !only for the observation sign to match (not used in this routine)
                        if (.not.control%strictOLD) then
                           HWires%CurrentSegment(conta)%orientadoalreves = &
                           (i1 < sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%i).or. &
                           (j1 < sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%j).or. &
                           (k1 < sgg%Med(HWires%WireTipoMedio(iwi))%wire(1)%segm(iwj-1)%k)
                        else
                           HWires%CurrentSegment(conta)%orientadoalreves =.false. !later corrected
                        end if
                     end if
                     !
                     select case (HWires%CurrentSegment(conta)%tipofield)
                     case (iEx)
                        ! default
                        HWires%CurrentSegment(conta)%Efield_wire2main => Ex(i1,j1,k1) 
                        HWires%CurrentSegment(conta)%Efield_main2wire => Ex(i1,j1,k1) 
!!al final se reapuntan los del thickness
                        HWires%CurrentSegment(conta)%delta=1.0_RKIND_wires / Idxe(i1)    !ojo esto de los delta habra que corregirlo  para uniones
                        HWires%CurrentSegment(conta)%deltaTransv1=1.0_RKIND_wires / Idyh(j1)
                        if (k1 <= sgg%ALLOC(iEz)%ZE) then !esta corriente en el limite de los alloc nunca se precisa
                           !updateo por exceso para que haya una celda de solapamiento y la union topologica de
                           !hilos no se pierda por culpa de la puta particion MPI
                           HWires%CurrentSegment(conta)%deltaTransv2=1.0_RKIND_wires / Idzh(k1)
                        else !no se comete error alguno
                           !es solo para que el indice no reviente 2012
                           HWires%CurrentSegment(conta)%deltaTransv2=1.0_RKIND_wires / Idzh(k1-1)
                        end if
                        !dama
                        HWires%CurrentSegment(conta)%ie     = i1+1
                        HWires%CurrentSegment(conta)%x      = j1+0.25_RKIND_wires
                        HWires%CurrentSegment(conta)%y      = k1+0.25_RKIND_wires 
                        sggmiE = sggmiEx(i1,j1,k1); call deembed_peclossyconformal_segments(sggmiE); sggmiEx(i1,j1,k1)= sggmiE !por si se ha modificado !ojo agresivo 180220 !ojo cambiado aqui 170323 de sitio pq no se habian puesto los deltatrans y salia division por cero
                        !fin dama
                      case (iEy)    
                        ! default
                        HWires%CurrentSegment(conta)%Efield_wire2main => Ey(i1,j1,k1) 
                        HWires%CurrentSegment(conta)%Efield_main2wire => Ey(i1,j1,k1)      
                        !
                        HWires%CurrentSegment(conta)%delta=1.0_RKIND_wires / Idye(j1)
                        if (k1 <= sgg%ALLOC(iEz)%ZE) then !esta corriente en el limite de los alloc nunca se precisa
                           !updateo por exceso para que haya una celda de solapamiento y la union topologica de
                           !hilos no se pierda por culpa de la puta particion MPI
                           HWires%CurrentSegment(conta)%deltaTransv1=1.0_RKIND_wires / Idzh(k1)
                        else
                           HWires%CurrentSegment(conta)%deltaTransv1=1.0_RKIND_wires / Idzh(k1-1) !no se comete error alguno
                           !es solo para que el indice no reviente 2012
                        end if
                        HWires%CurrentSegment(conta)%deltaTransv2=1.0_RKIND_wires / Idxh(i1)
                        !dama
                        HWires%CurrentSegment(conta)%je     = j1+1
                        HWires%CurrentSegment(conta)%x      = k1+0.25_RKIND_wires
                        HWires%CurrentSegment(conta)%y      = i1+0.25_RKIND_wires
                        sggmiE = sggmiEy(i1,j1,k1); call deembed_peclossyconformal_segments(sggmiE); sggmiEy(i1,j1,k1)= sggmiE !por si se ha modificado !ojo agresivo 180220 !ojo cambiado aqui 170323 de sitio pq no se habian puesto los deltatrans y salia division por cero
                        !fin dama
                      case (iEz)  
                        ! default
                        HWires%CurrentSegment(conta)%Efield_wire2main => Ez(i1,j1,k1) 
                        HWires%CurrentSegment(conta)%Efield_main2wire => Ez(i1,j1,k1)  
                        !
                        HWires%CurrentSegment(conta)%delta=1.0_RKIND_wires / Idze(k1)
                        HWires%CurrentSegment(conta)%deltaTransv1=1.0_RKIND_wires / Idxh(i1)
                        HWires%CurrentSegment(conta)%deltaTransv2=1.0_RKIND_wires / Idyh(j1)
                        !dama
                        HWires%CurrentSegment(conta)%ke     = k1+1
                        HWires%CurrentSegment(conta)%x      = i1+0.25_RKIND_wires
                        HWires%CurrentSegment(conta)%y      = j1+0.25_RKIND_wires
                        sggmiE = sggmiEz(i1,j1,k1); call deembed_peclossyconformal_segments(sggmiE); sggmiEz(i1,j1,k1)= sggmiE !por si se ha modificado !ojo agresivo 180220 !ojo cambiado aqui 170323 de sitio pq no se habian puesto los deltatrans y salia division por cero
                        !fin dama
                     end select
                  end if
               end if !del repetido
            end do
         end do
 
!!fin niapa 171216
      !
      HWires%NumCurrentSegments=conta
      !

      !hacer agujeros

      do i1=1,HWires%NumCurrentSegments
         segmento=>HWires%CurrentSegment(i1)
         i=segmento%i
         j=segmento%j
         k=segmento%k
         whatfield= segmento%tipofield
         IsEnd_norLeft_norRight=segmento%IsEnd_norLeft_norRight
         Is_LeftEnd=segmento%Is_LeftEnd
         Is_RightEnd=segmento%Is_RightEnd
         if ((i > SINPML_fullsize(whatfield)%XI).and. &
         (i < SINPML_fullsize(whatfield)%XE).and. &
         (j > SINPML_fullsize(whatfield)%YI).and. &
         (j < SINPML_fullsize(whatfield)%YE).and. &
         (k > SINPML_fullsize(whatfield)%ZI).and. &
         (k < SINPML_fullsize(whatfield)%ZE)) then
            if (control%makeholes.and.(.not.IsEnd_norLeft_norRight).and.(.not.Is_LeftEnd).and.(.not.Is_RightEnd)) then
                if (control%num_procs==0) then 
                    call stoponerror(control%layoutnumber,control%num_procs,'Makeholes not available for MPI. Stoppping. ')
                end if
               select case (whatfield)
                case (iEx)
                  sggmiHx(i  ,j  ,k      ) = 1
                  sggmiHx(i  ,j-1,k      ) = 1
                  sggmiHx(i  ,j  ,k -1   ) = 1
                  sggmiHx(i  ,j-1,k -1   ) = 1
                  sggmiHx(i+1,j  ,k      ) = 1
                  sggmiHx(i+1,j-1,k      ) = 1
                  sggmiHx(i+1,j  ,k -1   ) = 1
                  sggmiHx(i+1,j-1,k -1   ) = 1
                  if (.not.sgg%med(sggmiEx(i  ,j  ,k      ))%Is%ThinWire) then
                     sggmiEx(i  ,j  ,k      ) = 1
                     write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                     i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                     call WarnErrReport(buff)
                  end if
                  if (.not.sgg%med(sggmiEy(i  ,j  ,k      ))%Is%ThinWire)  then
                     sggmiEy(i  ,j  ,k      ) = 1
                     write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                     i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                     call WarnErrReport(buff)
                  end if
                  if (.not.sgg%med(sggmiEy(i  ,j-1,k      ))%Is%ThinWire)  then
                     sggmiEy(i  ,j-1,k      ) = 1
                     write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                     i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                     call WarnErrReport(buff)
                  end if
                  if       ((k <=  sgg%alloc(iEz)%ZE).and.(k >= sgg%alloc(iEz)%ZI)) then
                     if (.not.sgg%med(sggmiEz(i  ,j  ,k      ))%Is%ThinWire)  then
                        sggmiEz(i  ,j  ,k      ) = 1
                        write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', &
                        i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
                        call WarnErrReport(buff)
                     end if
                  end if
                  if   ((k-1 <=  sgg%alloc(iEz)%ZE).and.(k -1 >=sgg%alloc(iEz)%ZI)) then
                     if (.not.sgg%med(sggmiEz(i  ,j  ,k-1    ))%Is%ThinWire

} else {
                        sggmiEz(i, j, k - 1) = 1;
                        sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                                i, j, k, segmento.origIndex, dir(whatfield).c_str());
                        WarnErrReport(buff);
                    }
                }
            }
            if (!sgg.med(sggmiEy(i + 1, j, k)).Is.ThinWire) {
                sggmiEy(i + 1, j, k) = 1;
                sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                        i, j, k, segmento.origIndex, dir(whatfield).c_str());
                WarnErrReport(buff);
            }
            if (!sgg.med(sggmiEy(i + 1, j - 1, k)).Is.ThinWire) {
                sggmiEy(i + 1, j - 1, k) = 1;
                sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                        i, j, k, segmento.origIndex, dir(whatfield).c_str());
                WarnErrReport(buff);
            }
            if ((k <= sgg.alloc(iEz).ZE) && (k >= sgg.alloc(iEz).ZI)) {
                if (!sgg.med(sggmiEz(i + 1, j, k)).Is.ThinWire) {
                    sggmiEz(i + 1, j, k) = 1;
                    sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                            i, j, k, segmento.origIndex, dir(whatfield).c_str());
                    WarnErrReport(buff);
                }
            }
            if ((k - 1 <= sgg.alloc(iEz).ZE) && (k - 1 >= sgg.alloc(iEz).ZI)) {
                if (!sgg.med(sggmiEz(i + 1, j, k - 1)).Is.ThinWire) {
                    sggmiEz(i + 1, j, k - 1) = 1;
                    sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                            i, j, k, segmento.origIndex, dir(whatfield).c_str());
                    WarnErrReport(buff);
                }
            }
            break;
        case iEy:
            sggmiHy(i, j, k) = 1;
            sggmiHy(i, j, k - 1) = 1;
            sggmiHy(i - 1, j, k) = 1;
            sggmiHy(i - 1, j, k - 1) = 1;
            sggmiHy(i, j + 1, k) = 1;
            sggmiHy(i, j + 1, k - 1) = 1;
            sggmiHy(i - 1, j + 1, k) = 1;
            sggmiHy(i - 1, j + 1, k - 1) = 1;
            if (!sgg.med(sggmiEy(i, j, k)).Is.ThinWire) {
                sggmiEy(i, j, k) = 1;
                sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                        i, j, k, segmento.origIndex, dir(whatfield).c_str());
                WarnErrReport(buff);
            }
            if ((k <= sgg.alloc(iEz).ZE) && (k >= sgg.alloc(iEz).ZI)) {
                if (!sgg.med(sggmiEz(i, j, k)).Is.ThinWire) {
                    sggmiEz(i, j, k) = 1;
                    sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                            i, j, k, segmento.origIndex, dir(whatfield).c_str());
                    WarnErrReport(buff);
                }
            }
            if ((k - 1 <= sgg.alloc(iEz).ZE) && (k - 1 >= sgg.alloc(iEz).ZI)) {
                if (!sgg.med(sggmiEz(i, j, k - 1)).Is.ThinWire) {
                    sggmiEz(i, j, k - 1) = 1;
                    sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                            i, j, k, segmento.origIndex, dir(whatfield).c_str());
                    WarnErrReport(buff);
                }
            }
            if (!sgg.med(sggmiEx(i, j, k)).Is.ThinWire) {
                sggmiEx(i, j, k) = 1;
                sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                        i, j, k, segmento.origIndex, dir(whatfield).c_str());
                WarnErrReport(buff);
            }
            if (!sgg.med(sggmiEx(i - 1, j, k)).Is.ThinWire) {
                sggmiEx(i - 1, j, k) = 1;
                sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                        i, j, k, segmento.origIndex, dir(whatfield).c_str());
                WarnErrReport(buff);
            }
            if ((k <= sgg.alloc(iEz).ZE) && (k >= sgg.alloc(iEz).ZI)) {
                if (!sgg.med(sggmiEz(i, j + 1, k)).Is.ThinWire) {
                    sggmiEz(i, j + 1, k) = 1;
                    sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                            i, j, k, segmento.origIndex, dir(whatfield).c_str());
                    WarnErrReport(buff);
                }
            }
            if ((k - 1 <= sgg.alloc(iEz).ZE) && (k - 1 >= sgg.alloc(iEz).ZI)) {
                if (!sgg.med(sggmiEz(i, j + 1, k - 1)).Is.ThinWire) {
                    sggmiEz(i, j + 1, k - 1) = 1;
                    sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                            i, j, k, segmento.origIndex, dir(whatfield).c_str());
                    WarnErrReport(buff);
                }
            }
            if (!sgg.med(sggmiEx(i, j + 1, k)).Is.ThinWire) {
                sggmiEx(i, j + 1, k) = 1;
                sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                        i, j, k, segmento.origIndex, dir(whatfield).c_str());
                WarnErrReport(buff);
            }
            if (!sgg.med(sggmiEx(i - 1, j + 1, k)).Is.ThinWire) {
                sggmiEx(i - 1, j + 1, k) = 1;
                sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                        i, j, k, segmento.origIndex, dir(whatfield).c_str());
                WarnErrReport(buff);
            }
            break;
        case iEz:
            sggmiHz(i, j, k) = 1;
            sggmiHz(i - 1, j, k) = 1;
            sggmiHz(i, j - 1, k) = 1;
            sggmiHz(i - 1, j - 1, k) = 1;
            sggmiHz(i, j, k + 1) = 1;
            sggmiHz(i - 1, j, k + 1) = 1;
            sggmiHz(i, j - 1, k + 1) = 1;
            sggmiHz(i - 1, j - 1, k + 1) = 1;
            if (!sgg.med(sggmiEz(i, j, k)).Is.ThinWire) {
                sggmiEz(i, j, k) = 1;
                sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                        i, j, k, segmento.origIndex, dir(whatfield).c_str());
                WarnErrReport(buff);
            }
            if ((k <= sgg.alloc(iEx).ZE) && (k >= sgg.alloc(iEx).ZI)) {
                if (!sgg.med(sggmiEx(i, j, k)).Is.ThinWire) {
                    sggmiEx(i, j, k) = 1;
                    sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                            i, j, k, segmento.origIndex, dir(whatfield).c_str());
                    WarnErrReport(buff);
                }
                if (!sgg.med(sggmiEx(i - 1, j, k)).Is.ThinWire) {
                    sggmiEx(i - 1, j, k) = 1;
                    sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                            i, j, k, segmento.origIndex, dir(whatfield).c_str());
                    WarnErrReport(buff);
                }
            }
            if ((k <= sgg.alloc(iEy).ZE) && (k >= sgg.alloc(iEy).ZI)) {
                if (!sgg.med(sggmiEy(i, j, k)).Is.ThinWire) {
                    sggmiEy(i, j, k) = 1;
                    sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                            i, j, k, segmento.origIndex, dir(whatfield).c_str());
                    WarnErrReport(buff);
                }
                if (!sgg.med(sggmiEy(i, j - 1, k)).Is.ThinWire) {
                    sggmiEy(i, j - 1, k) = 1;
                    sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                            i, j, k, segmento.origIndex, dir(whatfield).c_str());
                    WarnErrReport(buff);
                }
            }
            if ((k + 1 <= sgg.alloc(iEx).ZE) && (k + 1 >= sgg.alloc(iEx).ZI)) {
                if (!sgg.med(sggmiEx(i, j, k + 1)).Is.ThinWire) {
                    sggmiEx(i, j, k + 1) = 1;
                    sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                            i, j, k, segmento.origIndex, dir(whatfield).c_str());
                    WarnErrReport(buff);
                }
                if (!sgg.med(sggmiEx(i - 1, j, k + 1)).Is.ThinWire) {
                    sggmiEx(i - 1, j, k + 1) = 1;
                    sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                            i, j, k, segmento.origIndex, dir(whatfield).c_str());
                    WarnErrReport(buff);
                }
            }
            if ((k + 1 <= sgg.alloc(iEy).ZE) && (k + 1 >= sgg.alloc(iEy).ZI)) {
                if (!sgg.med(sggmiEy(i, j, k + 1)).Is.ThinWire) {
                    sggmiEy(i, j, k + 1) = 1;
                    sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                            i, j, k, segmento.origIndex, dir(whatfield).c_str());
                    WarnErrReport(buff);
                }
                if (!sgg.med(sggmiEy(i, j - 1, k + 1)).Is.ThinWire) {
                    sggmiEy(i, j - 1, k + 1) = 1;
                    sprintf(buff, "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at %7d%7d%7d for WIRE-segment %7d %s",
                            i, j, k, segmento.origIndex, dir(whatfield).c_str());
                    WarnErrReport(buff);
                }
            }
            break;
        }
    }
}

// !!!calculo y gestion autoinduccion
for (int i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
    int jmed = HWires.CurrentSegment(i1).indexmed;
    double desp = HWires.CurrentSegment(i1).delta;
    double despT1 = HWires.CurrentSegment(i1).deltaTransv1;
    double despT2 = HWires.CurrentSegment(i1).deltaTransv2;
    
    // Esto no debe ser preciso. la hipotesis es que vuelve promediando a una celda
    if (control.wirethickness > 1) {        // solo tiene impacto en la %cte2 de acople y en la fuente de carga para que sea de voltaje
        despT1 = control.wirethickness * despT1;
        despT2 = control.wirethickness * despT2;
    }
    double r0 = HWires.CurrentSegment(i1).TipoWire.Radius;   // CON SU RADIO

    if (r0 < 1e-30) {
        sprintf(buff, "wir0_ERROR: ire radius cannot be null");
        if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff, true);
    }
    if (r0 < 1e-9 * desp) {
        sprintf(buff, "wir0_WARNING: WIRE radius too small %.2e", r0);
        WarnErrReport(buff);
    }
    if ((r0 > 0.5_RKIND_wires * despT1) || (r0 > 0.5_RKIND_wires * despT2)) {
        sprintf(buff, "wir0_WARNING: WIRE radius greater than half a space-step. Reduced to this limit %.2e",
                std::min(0.5_RKIND_wires * despT1, 0.5_RKIND_wires * despT2));
        if ((HWires.CurrentSegment(i1).k >= ZI) && (HWires.CurrentSegment(i1).k <= ZE)) WarnErrReport(buff);
    }

    std::string ind_model = control.inductance_model;
    // trim and adjustl equivalent
    while (!ind_model.empty() && ind_model.front() == ' ') ind_model.erase(ind_model.begin());
    while (!ind_model.empty() && ind_model.back() == ' ') ind_model.pop_back();

    if (ind_model == "berenger") {
        //---------------------------------------------------------------------------------
        // Slanted=Berenger
        // The second best one for the edelvik PEC box (Boutayeb's PPT is the best one for this case)
        // menos ruidoso que el de Boutayeb en uniones de hilos paralelos
        HWires.CurrentSegment(i1).Lind =
            (1.0_RKIND_wires / (4.0_RKIND_wires * pi * InvMu(jmed))) *
            (log((despT1 * despT1 + despT2 * despT2) / (4.0_RKIND_wires * r0 * r0)) +
             despT1 / despT2 * atan(despT2 / despT1) +
             despT2 / despT1 * atan(despT1 / despT2) +
             pi * r0 * r0 / (despT2 * despT1) - 3.0_RKIND_wires);
        //---------------------------------------------------------------------------------
        // Slanted corrected for wires of radius > 0.3_RKIND_wires Delta
        // Untested
        // just divides by a correction factor equal to that used by Boutayeb in his correction
        // proposed by Grando
        if ((r0 > 0.3_RKIND_wires * despT1) || (r0 > 0.3_RKIND_wires * despT2)) {
            HWires.CurrentSegment(i1).Lind = HWires.CurrentSegment(i1).Lind /
                                             (1.0_RKIND_wires - pi * r0 * r0 / (despT1 * despT2));
        }
    } else if (ind_model == "ledfelt") {
        //---------------------------------------------------------------------------------
        // Ledfelt thesis square symbols
        // The best for the antenna dipole, but it is not the best for edelvik's cavity
        // Dec'11 me ha dado muy bien al revisar el TC2_OF-SWAB-A2-B
        // warning: there is a mistake in formula 8.15 parenthesis
        HWires.CurrentSegment(i1).Lind =
            (1.0_RKIND_wires / (4.0_RKIND_wires * pi * InvMu(jmed))) *
            (log((despT1 * despT1 + despT2 * despT2) / (r0 * r0)) +
             despT1 / despT2 * atan(despT2 / despT1) +
             despT2 / despT1 * atan(despT1 / despT2) +
             pi * r0 * r0 / (16.0_RKIND_wires * despT2 * despT1) - 3.0_RKIND_wires);
        // Slanted makes this correction for radius>0.3_RKIND_wires Delta.  a>0.3_RKIND_wires Delta
        // I use it also in Ledleft !2012
        // never tested
        if ((r0 > 0.3_RKIND_wires * despT1) || (r0 > 0.3_RKIND_wires * despT2)) {
            HWires.CurrentSegment(i1).Lind = HWires.CurrentSegment(i1).Lind /
                                             (1.0_RKIND_wires - pi * r0 * r0 / (despT1 * despT2));
        }
    } else if (ind_model == "boutayeb") {
        //!!---------------------------------------------------------------------------------
        // Boutayeb's PPT. For radius>0.3_RKIND_wires Delta coincides with that of Berenger
        // adds a correction for radiuos <0.3_RKIND_wires Delta
        // and is equal to Guiffault for radius > 0.5_RKIND_wires Delta
        // The BEST ONE for edelvik_box and for HIRF NTC2.0_RKIND_wires
        HWires.CurrentSegment(i1).Lind =
            (1.0_RKIND_wires / (4.0_RKIND_wires * pi * InvMu(jmed))) *
            (log((despT1 * despT1 + despT2 * despT2) / (4.0_RKIND_wires * r0 * r0)) +
             despT1 / despT2 * atan(despT2 / despT1) +
             despT2 / despT1 * atan(despT1 / despT2) +
             pi * r0 * r0 / (despT2 * despT1) - 3.0_RKIND_wires);

        if ((r0 < 0.3_RKIND_wires * despT1) || (r0 < 0.3_RKIND_wires * despT2)) {
            HWires.CurrentSegment(i1).Lind = HWires.CurrentSegment(i1).Lind - 0.57_RKIND_wires / (4.0_RKIND_wires * pi * InvMu(jmed));
        }
        //---------------------------------------------------------------------------------
        // Boutayeb's PPT for radius>0.3_RKIND_wires Delta.  a>0.3_RKIND_wires Delta.
        // (Slanted corrects for 0.3_RKIND_wires Delta while boutayeb does it for 0.5_RKIND_wires Delta, I take Slanted's)
        // Untested
        // just divides by a correction factor (warning becomes negative for r0/delta >0.56)
        if ((r0 > 0.3_RKIND_wires * despT1) || (r0 > 0.3_RKIND_wires * despT2)) {
            HWires.CurrentSegment(i1).Lind = HWires.CurrentSegment(i1).Lind /
                                             (1.0_RKIND_wires - pi * r0 * r0 / (despT1 * despT2));
        }
        //!!---------------------------------------------------------------------------------
    } else {
        buff = "wir0_ERROR: Incorrect inductance model";
        WarnErrReport(buff, true);
    }
    if (HWires.CurrentSegment(i1).Lind < 0.0_RKIND_wires) {
        buff = "wir0_ERROR: Wrong self-inductance. ";
        WarnErrReport(buff, true);
    }
    // !!!!for junctions after holland
    // !!!    HWires%CurrentSegment(i1)%logRoverR0=log(((desp*despT1*despT2)**(1.0_RKIND_wires/3.0_RKIND_wires))/r0)  !HWires%CurrentSegment(i1)%Lind ! !holland pp90 (17)(18)
    // !!!!22/07/15 el cambio de radio es equivalente a un cambio de delta.... uso algo mejor
    // !!!  HWires%CurrentSegment(i1)%logRoverR0=HWires%CurrentSegment(i1)%Lind * InvMu(jmed) * InvEPS(jmed)
    // !!!!
}

//
LindProb.resize(HWires.NumCurrentSegments + 1);
for (int i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
    LindProb(i1) = true;
    HWires.CurrentSegment(i1).inv_Lind_acum = 1.0_RKIND_wires / HWires.CurrentSegment(i1).Lind;
    HWires.CurrentSegment(i1).HEUR_safety = (sgg.dt * sgg.dt) / (eps0 * HWires.CurrentSegment(i1).deltaTransv1 * HWires.CurrentSegment(i1).deltaTransv2);
}

for (int i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
    if (LindProb(i1)) {
        WireSegment* org = &HWires.CurrentSegment(i1);
        org->NumParallel = 1;
        org->Lind_acum = org->Lind;
        for (int j1 = i1 + 1; j1 <= HWires.NumCurrentSegments; ++j1) {
            WireSegment* fin = &HWires.CurrentSegment(j1);
            if ((org->i == fin->i) && (org->j == fin->j) && (org->k == fin->k) && (org->tipofield == fin->tipofield)) {

org->NumParallel = org->NumParallel + 1;
            if (control->stableradholland) org->Lind_acum = org->Lind_acum + fin->Lind;
         }
      }
      for (int j1 = i1 + 1; j1 < HWires->NumCurrentSegments; ++j1) {
         fin = HWires->CurrentSegment[j1];
         if ((org->i == fin->i) && (org->j == fin->j) && (org->k == fin->k) && (org->tipofield == fin->tipofield)) {
            fin->NumParallel = org->NumParallel;
            fin->Lind_acum = org->Lind_acum;
            LindProb[j1] = false;
         }
      }
   }
}
delete[] LindProb;

// !!!!solo para pruebas
// !!!open (629,file='param.txt',form='formatted')
// !!!    read (629,*) deltadummy
// !!!close (629)
// !!!    print *,'---------deltadummy=',deltadummy
// !!!do i1=1,HWires%NumCurrentSegments
// !!!    HWires%CurrentSegment(i1)%inv_Lind_acum=deltadummy*mu0
// !!!end do
// !!!!!!!!
// !!!!!!!!!!!!!!!!!\E7\E7\E7\E7\E7\E7\E7\E7\E7\E7
// !!!!!!!! hago 120715 mi criterior
dtcritico = sgg->dt;
for (int i1 = 1; i1 <= HWires->NumCurrentSegments; ++i1) {
   dummy = HWires->CurrentSegment[i1];
   jmed = dummy->indexmed;
   desp = dummy->delta;
   despT1 = dummy->deltaTransv1;
   despT2 = dummy->deltaTransv2;
   r0 = dummy->TipoWire->Radius;
   if (r0 < 1e-30) {
      sprintf(buff, "wir0_ERROR: wire radius cannot be null");
      if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff, true);
   }
   // !!!!!!!!!!!!correccion bruta 13/07/15
   deltadummy = dummy->inv_Lind_acum * dummy->HEUR_safety / 0.9_RKIND_wires; // EL 0.9 ES POR MI TRANQUILIDAD
   if (deltadummy > 1.0_RKIND_wires) {
      if (control->stableradholland) {
         if (strcmp(adjustl(control->inductance_model), "boutayeb") == 0) {
            b = -4.0 * pi * InvMu[jmed] * (dummy->Lind * deltadummy) + log(despT1 * despT1 + despT2 * despT2) + despT1 / despT2 * atan(despT2 / despT1) + despT2 / despT1 * atan(despT1 / despT2) - 3.0_RKIND_wires;
            a = pi / (despT2 * despT1);
            if ((r0 < 0.3_RKIND_wires * despT1) || (r0 < 0.3_RKIND_wires * despT2)) {
               B = B - 0.57;
            }
            newr0 = sqrt(-Lambert(-A * exp(b) / 4.0_RKIND_wires) / A);
            // !!!doublechecking
            b = -4.0 * pi * InvMu[jmed] * (dummy->Lind) + log(despT1 * despT1 + despT2 * despT2) + despT1 / despT2 * atan(despT2 / despT1) + despT2 / despT1 * atan(despT1 / despT2) - 3.0_RKIND_wires;
            a = pi / (despT2 * despT1);
            if ((r0 < 0.3_RKIND_wires * despT1) || (r0 < 0.3_RKIND_wires * despT2)) {
               B = B - 0.57;
            }
            OLDR0 = sqrt(-Lambert(-A * exp(b) / 4.0_RKIND_wires) / A);
            // !!!!!!!!!!
            sprintf(buff, "wir0_WARNING: AUTOMATIC CORRECTION OF L/mu0=%e for r0=%e%e %i wires at %i%i%i to L/mu0=%e for newr0=%e",
                    dummy->Lind / mu0, r0, oldr0, dummy->NumParallel, dummy->i, dummy->j, dummy->k,
                    dummy->Lind * deltadummy / mu0, newr0);
         } else {
            sprintf(buff, "wir0_WARNING: AUTOMATIC CORRECTION OF L/mu0=%e %i wires at %i%i%i to L/mu0=%e",
                    dummy->Lind / mu0, dummy->NumParallel, dummy->i, dummy->j, dummy->k,
                    dummy->Lind * deltadummy / mu0);
         }
         if ((dummy->k > ZI) && (dummy->k <= ZE)) WarnErrReport(buff);
         dummy->Lind = dummy->Lind * deltadummy; // bajo repartiendo proporcialmente
      } else {
         sprintf(buff, "wir0_SEVEREWARNING: L/mu0=%e in %i wires at %i%i%i smaller (posibly unstable) than L/mu0=%e",
                 dummy->Lind / mu0, dummy->NumParallel, dummy->i, dummy->j, dummy->k,
                 dummy->Lind * deltadummy / mu0);
         if ((dummy->k > ZI) && (dummy->k <= ZE)) WarnErrReport(buff);
         dtcritico = min(sgg->dt / sqrt(deltadummy), dtcritico);
      }
   }
}
// !!if (dtcritico<sgg%dt) then
// !!            write(buff,'(a,e9.2e2,a,e9.2e2)') &
// !!            &    'wir0_ERROR: UNSTABLE sgg%dt, decrease wire radius, number of Parallel WIREs, or make sgg%dt < ',dtcritico
// !!            call WarnErrReport(buff,.true.)
// !!end if

// !!!!!!!!!!!fin !mi criterio 13/07/15

// Grounding R_RightEnd and R_LeftEnd resistances info
for (int i1 = 1; i1 <= HWires->NumCurrentSegments; ++i1) {
   segmento = HWires->CurrentSegment[i1];
   //
   if (segmento->TipoWire->HasAbsorbing_LeftEnd) {
      if (segmento->Is_LeftEnd) {
         segmento->HasAbsorbing_LeftEnd = true;
         sprintf(buff, "wir1_INFO: Absorbing conditions in terminal LeftEnd segment %7i%7i%7i%7i%7i",
                 segmento->origIndex, segmento->i, segmento->j, segmento->k, segmento->tipofield);
         if ((segmento->k >= ZI) && (segmento->k <= ZE) && control->verbose) WarnErrReport(buff);
      }
   }
   if (segmento->TipoWire->HasAbsorbing_RightEnd) {
      if (segmento->Is_RightEnd) {
         segmento->HasAbsorbing_RightEnd = true;
         sprintf(buff, "wir1_WARNING: Absorbing conditions  in terminal RightEnd segment %7i%7i%7i%7i%7i",
                 segmento->origIndex, segmento->i, segmento->j, segmento->k, segmento->tipofield);
         if ((segmento->k >= ZI) && (segmento->k <= ZE)) WarnErrReport(buff);
      }
   }
   //
   if (segmento->TipoWire->HasParallel_LeftEnd) {
      if (segmento->Is_LeftEnd) {
         segmento->HasParallel_LeftEnd = true;
         sprintf(buff, "wir1_WARNING: Parallel RLC in terminal LeftEnd segment %7i%7i%7i%7i%7i",
                 segmento->origIndex, segmento->i, segmento->j, segmento->k, segmento->tipofield);
         if ((segmento->k >= ZI) && (segmento->k <= ZE)) WarnErrReport(buff);
      }
   }
   if (segmento->TipoWire->HasParallel_RightEnd) {
      if (segmento->Is_RightEnd) {
         segmento->HasParallel_RightEnd = true;
         sprintf(buff, "wir1_WARNING: Parallel RLC in terminal RightEnd segment %7i%7i%7i%7i%7i",
                 segmento->origIndex, segmento->i, segmento->j, segmento->k, segmento->tipofield);
         if ((segmento->k >= ZI) && (segmento->k <= ZE)) WarnErrReport(buff);
      }
   }

   if (segmento->TipoWire->HasSeries_LeftEnd) {
      if (segmento->Is_LeftEnd) {
         segmento->HasSeries_LeftEnd = true;
         sprintf(buff, "wir1_WARNING: Series RLC in terminal LeftEnd segment %7i%7i%7i%7i%7i",
                 segmento->origIndex, segmento->i, segmento->j, segmento->k, segmento->tipofield);
         if ((segmento->k >= ZI) && (segmento->k <= ZE)) WarnErrReport(buff);
      }
   }
   if (segmento->TipoWire->HasSeries_RightEnd) {
      if (segmento->Is_RightEnd) {
         segmento->HasSeries_RightEnd = true;
         sprintf(buff, "wir1_WARNING: Series RLC in terminal RightEnd segment %7i%7i%7i%7i%7i",
                 segmento->origIndex, segmento->i, segmento->j, segmento->k, segmento->tipofield);
         if ((segmento->k >= ZI) && (segmento->k <= ZE)) WarnErrReport(buff);
      }
   }
}
//
// Create the final update constants for the advance of the currents
// It takes into account the extra inductance and resistance per unit length specified in ORIGINAL
// It also takes into account the Series/Parallel Grounding Inductance at and the end segments TR and TL !untested
// Junctions do no affect to these constants (later taken into account by means of the fractionplus and
// fractionminus constants)
for (int i1 = 1; i1 <= HWires->NumCurrentSegments; ++i1) {
   // constantes de actualizacion
   dummy = HWires->CurrentSegment[i1];
   resist = 0.0_RKIND_wires;

   // !!!for lossy groundings
   i = dummy->i;
   j = dummy->j;
   k = dummy->k;
   whatfield = dummy->tipofield;
   //
   rlossy = 0.0_RKIND_wires;
   sigt = 0.0_RKIND_wires;
   sigtPlus = 0.0_RKIND_wires;
   sigtMinu = 0.0_RKIND_wires;
   IsLossy = false; IsPEC = false;
   IsLossyPlus = false; IsPECPlus = false;
   IsLossyMinu = false; IsPECMinu = false;
   //
   switch (whatfield) {
       case iEx:
          esPML = sgg->med[sggmiEx(i, j, k)]->is->PML;
          break;
       case iEy:
          esPML = sgg->med[sggmiEy(i, j, k)]->is->PML;
          break;
       case iEz:
          esPML = sgg->med[sggmiEz(i, j, k)]->is->PML;
          break;
   }
   if (esPML) {
      continue;
      // !!\C7  dummy%IsShielded=.true.
   } else {
       if ((k <= sgg->alloc[iEZ]->ZE) && (k >= sgg->alloc[iEZ]->ZI)) {
          kmenos1 = k - 1;
          kmas1 = k + 1;
          //
          // esta informacion solo se utiliza si realmente luego hay un nodo terminal y se suma la resistencia. En cualquier otro caso no sirver para nada
          // de todos modos hay un bug en la deteccion. sgg 110815
          // se usal la informacion nodal (lo que sigue algun d\EDa se mover\E1 a la rutina de generacion nodal que se creo en preprocess y se podra dejar solo lo que sigue 110815
          // !!!bug 270224 gg Cuando se sobrepasa el MPI  kmenos1=k o kmas1=k da un buggy error. Pero no se toma decision alguna. se puede ignorar con -ignoreerrors
          switch (whatfield) {
           case iEx:
             med[0] = sggMiEx[i + 1][j][k];
             med[1] = sggMiEy[i + 1][j][k];
             med[2] = sggMiEy[i + 1][j - 1][k];
             med[3] = sggMiEz[i + 1][j][k];
             if (kmenos1 < sgg->alloc[iEz]->ZI) { kmenos1 = k; }
             med[4] = sggMiEz[i + 1][j][kmenos1];
             med[5] = sggMiNo[i + 1][j][k];
             //
             med[6] = sggMiEx[i - 1][j][k];
             med[7] = sggMiEy[i][j][k];
             med[8] = sggMiEy[i][j - 1][k];
             med[9] = sggMiEz[i][j][k];
             if (kmenos1 < sgg->alloc[iEz]->ZI) { kmenos1 = k; }
             med[10] = sggMiEz[i][j][kmenos1];
             med[11] = sggMiNo[i][j][k];
             break;
           case iEy:
             med[0] = sggMiEy[i][j + 1][k];
             med[1] = sggMiEz[i][j + 1][k];
             if (kmenos1 < sgg->alloc[iEz]->ZI) { kmenos1 = k; }
             med[2] = sggMiEz[i][j + 1][kmenos1];
             med[3] = sggMiEx[i][j + 1][k];
             med[4] = sggMiEx[i - 1][j + 1][k];
             med[5] = sggMiNo[i][j + 1][k];
             //
             med[6] = sggMiEy[i][j - 1][k];
             med[7] = sggMiEz[i][j][k];
             if (kmenos1 < sgg->alloc[iEz]->ZI) { kmenos1 = k; }
             med[8] = sggMiEz[i][j][kmenos1];
             med[9] = sggMiEx[i][j][k];
             med[10] = sggMiEx[i - 1][j][k];
             med[11] = sggMiNo[i][j][k];
             break;
           case iEz:
             // !!!ojooooo 270224 esta logica esta mal porque machaco las variables kmas1 y kmenos1.... corregir.... no tiene impacto pq el bucle es informativo y no se toman decisiones. en todo caso solo afectaria a MPI!
             // !!!pero esta no es la razon 27024 por la que sucede el error gg   (    8/   40) wir1_BUGGYERROR:  Lossy, pec,  2.           329         180         160   15000.0000000000      F F
             if (kmas1 > sgg->alloc[iEz]->ZE) { kmas1 = k; }
             med[0] = sggMiEz[i][j][kmas1];
             if (kmas1 > sgg->alloc[iEx]->ZE) { kmas1 = k; }
             med[1] = sggMiEx[i][j][kmas1];
             if (kmas1 > sgg->alloc[iEx]->ZE) { kmas1 = k; }
             med[2] = sggMiEx[i - 1][j][kmas1];
             if (kmas1 > sgg->alloc[iEy]->ZE) { kmas1 = k; }
             med[3] = sggMiEy[i][j][kmas1];
             if (kmas1 > sgg->alloc[iEy]->ZE) { kmas1 = k; }
             med[4] = sggMiEy[i][j - 1][kmas1];
             if (kmas1 > sgg->alloc[iHz]->ZE) { kmas1 = k; }
             med[5] = sggMiNo[i][j][kmas1];
             //
             if (kmenos1 < sgg->alloc[iEz]->ZI) { kmenos1 = k; }
             med[6] = sggMiEz[i][j][kmenos1];
             med[7] = sggMiEx[i][j][k];
             med[8] = sggMiEx[i - 1][j][k];
             med[9] = sggMiEy[i][j][k];
             med[10] = sggMiEy[i][j - 1][k];
             med[11] = sggMiNo[i][j][k];
             break;
          }
          // hay que tratar cada extremo por separado
          for (int nm = 0; nm <= 5; ++nm) {
             IsLossyPlus = IsLossyPlus || sgg->Med[med[nm]]->Is->Lossy;
             IsPECPlus = IsPECPlus || sgg->med[med[nm]]->is->PEC || (med[nm] == 0);
             if (!IsPECPlus) {
                sigtPlus = max(sigtPlus, sgg->Med[med[nm]]->sigma);
             }
          }
          for (int nm = 6; nm <= 11; ++nm) {
             IsLossyMinu = IsLossyMinu || sgg->Med[med[nm]]->Is->Lossy;
             IsPECMinu = IsPECMinu || sgg->med[med[nm]]->is->PEC || (med[nm] == 0);
             if (!IsPECMInu) {
                sigtMinu = max(sigtMinu, sgg->Med[med[nm]]->sigma);
             }
          }
          // !!!telaranias quitadas 060215
          if (IsPECPlus) { // domina el pec en un nodo con multiples medios
             IsLossyPlus = false;
          }
          if (IsPECMinu) { // domina el pec en un nodo con multiples medios
             IsLossyMinu = false;
          }
          sigt = max(sigtplus, sigtMinu); // con que uno de los dos extremos sea Lossy hay que corregir la resistencia de contacto
          IsLossy = IsLossyPlus || IsLossyMinu; // con que uno de los dos extremos sea Lossy hay que corregir la resistencia de contacto
          IsPEC = IsPECPlus && IsPECMinu; // !!!importante para ser PEC tienen que serlo los dos (es un caso trivial de un segemento unido a pec por los dos sitios. pero se da en el siva)
          // !!!checking de coherencia
          if (isPEC && IsLossy) { // es decir: los dos extremos PEC y alguno lossy es que hay un error
             sprintf(buff, "wir1_BUGGYERROR:  Lossy, pec, 1.  %i %i %i %e %i %i", i, j, k, sigt, ispec, isLossy);
             WarnErrReport(buff, true);
          }
          // !!!
          if ((!ispec) && (!isLossy)) { // algun extremo no pec y ambos extremos no lossy debe haber error si la conductividad no es nula
             if (abs(sigt) > 1.0e-19_RKIND_wires) {
                // 270224 creo que el bug gg 270224 es simplemente que conviven pec y lossy en plus o minus siendo el otro minus/plus vacio. Se voidea
                // 270224 el islosyplus/minus a false por ser pec. el ispec es false porque uno es vacio y el islossy tambien es vacio
                // 270224 pero se ha almacenado alguna sigt no nula. El arreglo es simplemente convertir esto en un warning de casuistica de conexionado en vez de un buggyerror
                // !!comentado a 010324 por los motivos anteriores    write (buff,*)  'wir1_BUGGYERROR:  Lossy, pec,  2.  ', i,j,k,sigt,ispec,isLossy
                // !!comentado a 010324 por los motivos anteriores    call WarnErrReport(buff,.true.)
             }
             sprintf(buff, "wir1_WARNING:  Lossy, pec,  2.: A wire is connected at some ending both to a lossy and to a PEC edge. Assuming it PEC  %i %i %i %e %i %i %i %i",
                     i, j, k, sigt, ispec, isLossy, IsPECPlus, IsPECMinu);
             WarnErrReport(buff);
          }
          if (isLossy) { // alguno extremo lossy con conductividad desconocida (multiports de ss)
             if (abs(sigt) < 1.0e-19_RKIND_wires) {
                sigt = 1e4; // asignale una resistencia de contacto por defecto si es nula !tipico de los composites de ss !habra algun dia que afinar esto
                sprintf(buff, "wir1_WARNING:  A Lossy segment with unknown conductivity. Assuming a STANDARD value of 1e4 S/m %i %i %i %e %i %i %i %i",
                        i, j, k, sigt, ispec, isLossy, IsLossyPlus, IsLossyMinu);
                WarnErrReport(buff);
             }
          }

          // !!!!!!!!!!!!!!!!!!!!hasta aqui casuistica. Aniade ahora la resistencia si procede con la formula de tercero de fisicas
          if (isLossy) {
             // rlossy=rlossy + 1.0_RKIND_wires/(2.0_RKIND_wires * pi*dummy%TipoWire%radius*sigt)/dummy%delta   !p.u.l.     !pruebas distintas 0223
             //  rlossy=rlossy + 1.0_RKIND_wires/(2.0_RKIND_wires * pi*dummy%DELTA          *sigt)/dummy%delta   !p.u.l.     !pruebas distintas 0223
             rlossy = rlossy; // + 1.0_RKIND_wires/(2.0_RKIND_wires * pi * (control%factordelta*dummy%DELTA +control%factorradius*dummy%TipoWire%radius) *sigt)/dummy%delta   !p.u.l. !comentado a 290323 porque no me fio de esto hasta que no se valide bien agb
          }
       }
   }
   //
   resist = dummy->TipoWire->R;
   givenautoin = HWires->CurrentSegment[i1]->TipoWire->L;
   dummy->givenautoin = givenautoin;
   dummy->resist = resist;
//
   resist_devia = dummy->TipoWire->R_devia;
   givenautoin_devia = HWires->CurrentSegment[i1]->TipoWire->L_devia;
   dummy->givenautoin_devia = givenautoin_devia;
   dummy->resist_devia = resist_devia;

   // !!bug'inest OLD 020413 !\E7 ahora solo se tratan resistencias y se aniaden a los segementos finales

   if ((dummy->Is_LeftEnd) && (isLossy || isLossy)) {
      // no tengo en cuenta el caso particularisimo de un solo segmento conectado a lossy por los dos extremos !habria que sumarle la resistencia dos veces pero la casuistica se enfollona !\E7
      if ((!dummy->HasParallel_LeftEnd) && (!dummy->HasSeries_LeftEnd) && (!dummy->HasAbsorbing_LeftEnd)) {
         resist = resist + rlossy;
         sprintf(buff, "wir1_INFO: Adding Lossy material resistence to LeftEnd segment in contact with lossy without a terminal RLC %7i%7i%7i%7i %s",
                 dummy->origIndex, dummy->i, dummy->j, dummy->k, dir(dummy->tipofield));
      } else {
         resist = resist + rlossy;
         sprintf(buff, "wir1_INFO: Adding Lossy material resistence to LeftEnd segment grounded through RLC %7i%7i%7i%7i %s",
                 dummy->origIndex, dummy->i, dummy->j, dummy->k, dir(dummy->tipofield));
      }
      if ((dummy->k > ZI) && (dummy->k <= ZE) && (dummy->tipofield != iEz) && control->verbose) WarnErrReport(buff);
      if ((dummy->k >= ZI) && (dummy->k <= ZE) && (dummy->tipofield == iEz) && control->verbose) WarnErrReport(buff);
   }
   if ((dummy->Is_RightEnd) && (isLossy || isLossy)) {
      // no tengo en cuenta el caso particularisimo de un solo segmento conectado a lossy por los dos extremos !habria que sumarle la resistencia dos veces pero la casuistica se enfollona !\E7
      if ((!dummy->HasParallel_RightEnd) && (!dummy->HasSeries_RightEnd) && (!dummy->HasAbsorbing_Ri

} else {
                resist = resist + rlossy;
                write(buff, '(a,4i7,a)') = 'wir1_INFO: Adding Lossy material resistence to RightEnd segment grounded through RLC ' + 
                    dummy->origIndex + dummy->i + dummy->j + dummy->k + ' ' + dir(dummy->tipofield);
            }
            if ((dummy->k > ZI) && (dummy->k <= ZE) && (dummy->tipofield != iEz) && control->verbose) {
                WarnErrReport(buff);
            }
            if ((dummy->k >= ZI) && (dummy->k <= ZE) && (dummy->tipofield == iEz) && control->verbose) {
                WarnErrReport(buff);
            }
        }
        // lOS QUE NO TENGAN RESITENCIAS Y ESTEN EN ABIERTO NO LOS CONECTO A LOSSY. sI SE QUIEREN HACER CONEXIONES A LOSSY
        // HAY QUE ESPECIFICAR UNA RESISTENCIA LUMPED
        // lUEGO SI SE CONECTARAN A pec DIRECTAMENTE SI LA TOPOLOGIA LO MANDA
        //
        if ((dummy->IsEnd_norLeft_norRight) && (isLossy || isLossy)) {
            //no tengo en cuenta el caso particularisimo de un solo segmento conectado a lossy por los dos extremos !habria que sumarle la resistencia dos veces pero la casuistica se enfollona !\E7
            if ((!dummy->HasParallel_LeftEnd) && (!dummy->HasSeries_LeftEnd) && (!dummy->HasParallel_RightEnd) && (!dummy->HasSeries_RightEnd) && 
                (!dummy->HasAbsorbing_RightEnd)) {
                resist = resist + rlossy;
                write(buff, '(a,4i7,a)') = 'wir1_INFO: Adding Lossy material resistence to Ending segment (other) segment in contact with lossy without a terminal RLC ' + 
                    dummy->origIndex + dummy->i + dummy->j + dummy->k + ' ' + dir(dummy->tipofield);
                if ((dummy->k > ZI) && (dummy->k <= ZE) && (dummy->tipofield != iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
                if ((dummy->k >= ZI) && (dummy->k <= ZE) && (dummy->tipofield == iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
            } else {
                resist = resist + rlossy;
                write(buff, '(a,4i7,a)') = 'wir1_BUGGYERROR:  Lossy material resistence to Ending (other) segment grounded through RLC () ' + 
                    dummy->origIndex + dummy->i + dummy->j + dummy->k + ' ' + dir(dummy->tipofield);
                if ((dummy->k > ZI) && (dummy->k <= ZE) && (dummy->tipofield != iEz)) {
                    WarnErrReport(buff, true);
                }
                if ((dummy->k >= ZI) && (dummy->k <= ZE) && (dummy->tipofield == iEz)) {
                    WarnErrReport(buff, true);
                }
            }
        }

        if (dummy->HasParallel_RightEnd) {
            givenautoin = givenautoin + dummy->TipoWire->Parallel_L_RightEnd / dummy->delta; //se le suma la autoinduccion !2011 \E7 untested
            dummy->givenautoin = givenautoin;
            //stoch
            givenautoin_devia = givenautoin_devia + dummy->TipoWire->Parallel_L_RightEnd_devia / dummy->delta; //se le suma la autoinduccion !2011 \E7 untested
            dummy->givenautoin_devia = givenautoin_devia;

            //!!bug'inest OLD 020413 !\E7 ahora solo se tratan resistencias y se aniaden a los segmentos finales

            resist = resist + dummy->TipoWire->Parallel_R_RightEnd / dummy->delta; //se le suma la autoinduccion !2011 \E7 untested
            resist_devia = resist_devia + dummy->TipoWire->Parallel_R_RightEnd_devia / dummy->delta; //se le suma la autoinduccion !2011 \E7 untested
            if (dummy->TipoWire->Parallel_R_RightEnd != 0.0_RKIND_wires) {
                write(buff, '(a,4i7,a)') = 'wir1_INFO: Adding Parallel RightEnd Resistance in segment ' + 
                    dummy->origIndex + dummy->i + dummy->j + dummy->k + ' ' + dir(dummy->tipofield);
                if ((dummy->k > ZI) && (dummy->k <= ZE) && (dummy->tipofield != iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
                if ((dummy->k >= ZI) && (dummy->k <= ZE) && (dummy->tipofield == iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
            } else {
                write(buff, '(a,4i7,a)') = 'wir1_INFO: Adding Parallel RightEnd null-Resistance in segment ' + 
                    dummy->origIndex + dummy->i + dummy->j + dummy->k + ' ' + dir(dummy->tipofield);
                if ((dummy->k > ZI) && (dummy->k <= ZE) && (dummy->tipofield != iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
                if ((dummy->k >= ZI) && (dummy->k <= ZE) && (dummy->tipofield == iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
            }

            //(ojo que es per unit length la intrinsea)
            if (dummy->TipoWire->Parallel_L_RightEnd != 0.0_RKIND_wires) {
                write(buff, '(a,4i7,a)') = 'wir1_INFO: Adding Parallel RightEnd Inductance in segment ' + 
                    dummy->origIndex + dummy->i + dummy->j + dummy->k + ' ' + dir(dummy->tipofield);
                if ((dummy->k > ZI) && (dummy->k <= ZE) && (dummy->tipofield != iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
                if ((dummy->k >= ZI) && (dummy->k <= ZE) && (dummy->tipofield == iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
            }
            //aniado tambien al ultimo segmento la resistencia y peto si hay capacitancias
            if (dummy->TipoWire->Parallel_C_RightEnd >= 1.0e-12_RKIND_wires) {
                write(buff, '(a,4i7,a)') = 'wir1_ERROR: (Currently unsupported)  Capacitances in Parallel RightEnd at segment ' + 
                    dummy->origIndex + dummy->i + dummy->j + dummy->k + ' ' + dir(dummy->tipofield);
                if ((dummy->k > ZI) && (dummy->k <= ZE) && (dummy->tipofield != iEz)) {
                    WarnErrReport(buff, true);
                }
                if ((dummy->k >= ZI) && (dummy->k <= ZE) && (dummy->tipofield == iEz)) {
                    WarnErrReport(buff, true);
                }
            } else {
                dummy->TipoWire->Parallel_C_RightEnd = 0.0_RKIND_wires;
            }
        }
        if (dummy->HasParallel_LeftEnd) {
            givenautoin = givenautoin + dummy->TipoWire->Parallel_L_LeftEnd / dummy->delta; //se le suma la autoinduccion !2011 \E7 untested
            dummy->givenautoin = givenautoin;
            //
            givenautoin_devia = givenautoin_devia + dummy->TipoWire->Parallel_L_LeftEnd_devia / dummy->delta; //se le suma la autoinduccion !2011 \E7 untested
            dummy->givenautoin_devia = givenautoin_devia;
            //!!bug'inest OLD 020413 !\E7 ahora solo se tratan resistencias y se aniaden a los segementos finales

            resist = resist + dummy->TipoWire->Parallel_R_LeftEnd / dummy->delta; //se le suma la autoinduccion !2011 \E7 untested
            resist_devia = resist_devia + dummy->TipoWire->Parallel_R_LeftEnd_devia / dummy->delta; //se le suma la autoinduccion !2011 \E7 untested
            if (dummy->TipoWire->Parallel_R_LeftEnd != 0.0_RKIND_wires) {
                write(buff, '(a,4i7,a)') = 'wir1_INFO: Adding Parallel LeftEnd Resistance in segment ' + 
                    dummy->origIndex + dummy->i + dummy->j + dummy->k + ' ' + dir(dummy->tipofield);
                if ((dummy->k > ZI) && (dummy->k <= ZE) && (dummy->tipofield != iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
                if ((dummy->k >= ZI) && (dummy->k <= ZE) && (dummy->tipofield == iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
            } else {
                write(buff, '(a,4i7,a)') = 'wir1_INFO: Adding Parallel LeftEnd null-Resistance in segment ' + 
                    dummy->origIndex + dummy->i + dummy->j + dummy->k + ' ' + dir(dummy->tipofield);
                if ((dummy->k > ZI) && (dummy->k <= ZE) && (dummy->tipofield != iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
                if ((dummy->k >= ZI) && (dummy->k <= ZE) && (dummy->tipofield == iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
            }

            if (dummy->TipoWire->Parallel_L_LeftEnd != 0.0_RKIND_wires) {
                write(buff, '(a,4i7,a)') = 'wir1_INFO: Adding Parallel LeftEnd Inductance in segment ' + 
                    dummy->origIndex + dummy->i + dummy->j + dummy->k + ' ' + dir(dummy->tipofield);
                if ((dummy->k > ZI) && (dummy->k <= ZE) && (dummy->tipofield != iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
                if ((dummy->k >= ZI) && (dummy->k <= ZE) && (dummy->tipofield == iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
            }
            if (dummy->TipoWire->Parallel_C_LeftEnd >= 1.0e-12_RKIND_wires) {
                write(buff, '(a,4i7,a)') = 'wir1_ERROR: (Currently unsupported)  Capacitances in Parallel LeftEnd at segment ' + 
                    dummy->origIndex + dummy->i + dummy->j + dummy->k + ' ' + dir(dummy->tipofield);
                if ((dummy->k > ZI) && (dummy->k <= ZE) && (dummy->tipofield != iEz)) {
                    WarnErrReport(buff, true);
                }
                if ((dummy->k >= ZI) && (dummy->k <= ZE) && (dummy->tipofield == iEz)) {
                    WarnErrReport(buff, true);
                }
            } else {
                dummy->TipoWire->Parallel_C_LeftEnd = 0.0_RKIND_wires;
            }
        }
        //
        if (dummy->HasSeries_RightEnd) {
            givenautoin = givenautoin + dummy->TipoWire->Series_L_RightEnd / dummy->delta; //se le suma la autoinduccion !2011 \E7 untested
            resist = resist + dummy->TipoWire->Series_R_RightEnd / dummy->delta;
            dummy->givenautoin = givenautoin;
            dummy->resist = resist;
            //
            givenautoin_devia = givenautoin_devia + dummy->TipoWire->Series_L_RightEnd_devia / dummy->delta; //se le suma la autoinduccion !2011 \E7 untested
            resist_devia = resist_devia + dummy->TipoWire->Series_R_RightEnd_devia / dummy->delta;
            dummy->givenautoin_devia = givenautoin_devia;
            dummy->resist_devia = resist_devia;
            //(ojo que es per unit length la intrinsea)
            if (dummy->TipoWire->Series_L_RightEnd != 0.0_RKIND_wires) {
                write(buff, '(a,4i7,a)') = 'wir1_INFO: Adding Series RightEnd Inductance in segment ' + 
                    dummy->origIndex + dummy->i + dummy->j + dummy->k + ' ' + dir(dummy->tipofield);
                if ((dummy->k > ZI) && (dummy->k <= ZE) && (dummy->tipofield != iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
                if ((dummy->k >= ZI) && (dummy->k <= ZE) && (dummy->tipofield == iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
            }
            if (dummy->TipoWire->Series_R_RightEnd != 0.0_RKIND_wires) {
                write(buff, '(a,4i7,a)') = 'wir1_INFO: Adding Series RightEnd Resistance in segment ' + 
                    dummy->origIndex + dummy->i + dummy->j + dummy->k + ' ' + dir(dummy->tipofield);
                if ((dummy->k > ZI) && (dummy->k <= ZE) && (dummy->tipofield != iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
                if ((dummy->k >= ZI) && (dummy->k <= ZE) && (dummy->tipofield == iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
            }
            if (dummy->TipoWire->Series_C_RightEnd <= 1.0e7_RKIND_wires) {
                write(buff, '(a,4i7,a)') = 'wir1_ERROR: (Currently unsupported)  Capacitances smaller than 1.0e7_RKIND_wires (inf) in Series RightEnd at segment ' + 
                    dummy->origIndex + dummy->i + dummy->j + dummy->k + ' ' + dir(dummy->tipofield);
                if ((dummy->k > ZI) && (dummy->k <= ZE) && (dummy->tipofield != iEz)) {
                    WarnErrReport(buff, true);
                }
                if ((dummy->k >= ZI) && (dummy->k <= ZE) && (dummy->tipofield == iEz)) {
                    WarnErrReport(buff, true);
                }
            } else {
                dummy->TipoWire->Series_C_RightEnd = 2.0e7_RKIND_wires;
            }

        }
        if (dummy->HasSeries_LeftEnd) {
            givenautoin = givenautoin + dummy->TipoWire->Series_L_LeftEnd / dummy->delta; //se le suma la autoinduccion !2011 \E7 untested
            resist = resist + dummy->TipoWire->Series_R_LeftEnd / dummy->delta; //se le suma la autoinduccion !2011 \E7 untested
            dummy->givenautoin = givenautoin;
            dummy->resist = resist;
            //
            givenautoin_devia = givenautoin_devia + dummy->TipoWire->Series_L_LeftEnd_devia / dummy->delta; //se le suma la autoinduccion !2011 \E7 untested
            resist_devia = resist_devia + dummy->TipoWire->Series_R_LeftEnd_devia / dummy->delta; //se le suma la autoinduccion !2011 \E7 untested
            dummy->givenautoin_devia = givenautoin_devia;
            dummy->resist_devia = resist_devia;
            if (dummy->TipoWire->Series_L_LeftEnd != 0.0_RKIND_wires) {
                write(buff, '(a,4i7,a)') = 'wir1_INFO: Adding Series LeftEnd Inductance in segment ' + 
                    dummy->origIndex + dummy->i + dummy->j + dummy->k + ' ' + dir(dummy->tipofield);
                if ((dummy->k > ZI) && (dummy->k <= ZE) && (dummy->tipofield != iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
                if ((dummy->k >= ZI) && (dummy->k <= ZE) && (dummy->tipofield == iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
            }
            if (dummy->TipoWire->Series_R_LeftEnd != 0.0_RKIND_wires) {
                write(buff, '(a,4i7,a)') = 'wir1_INFO: Adding Series LeftEnd Resistance in segment ' + 
                    dummy->origIndex + dummy->i + dummy->j + dummy->k + ' ' + dir(dummy->tipofield);
                if ((dummy->k > ZI) && (dummy->k <= ZE) && (dummy->tipofield != iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
                if ((dummy->k >= ZI) && (dummy->k <= ZE) && (dummy->tipofield == iEz) && control->verbose) {
                    WarnErrReport(buff);
                }
            }
            if (dummy->TipoWire->Series_C_LeftEnd <= 1.0e7_RKIND_wires) {
                write(buff, '(a,4i7,a)') = 'wir1_ERROR: (Currently unsupported)  Capacitances  smaller than 1.0e7_RKIND_wires (inf) in Series LeftEnd atn segment ' + 
                    dummy->origIndex + dummy->i + dummy->j + dummy->k + ' ' + dir(dummy->tipofield);
                if ((dummy->k > ZI) && (dummy->k <= ZE) && (dummy->tipofield != iEz)) {
                    WarnErrReport(buff, true);
                }
                if ((dummy->k >= ZI) && (dummy->k <= ZE) && (dummy->tipofield == iEz)) {
                    WarnErrReport(buff, true);
                }
            } else {
                dummy->TipoWire->Series_C_LeftEnd = 2.0e7_RKIND_wires;
            }

        }
    
        dummy->R = resist;
        dummy->L = givenautoin;
        dummy->givenautoin = givenautoin;
        dummy->resist = resist;
        //
        dummy->R_devia = resist_devia;
        dummy->L_devia = givenautoin_devia;
        dummy->givenautoin_devia = givenautoin_devia;
        dummy->resist_devia = resist_devia;
        //!!ojooooo  acumulo en %lind toda la autoinduccion para que los calculos de capacidad la tengan en cuenta completa    
        if (!control->fieldtotl) {
            dummy->Lind = dummy->Lind + givenautoin;
            //
            dummy->Lind_devia = dummy->Lind_devia + givenautoin_devia;
        } else {
            if (givenautoin <= tiny(1.0e-1)) {
                StopOnError(0, 0, 'Fieldtotl not compatible with null given self inductance. ');
            }
            //
            dummy->Lind = givenautoin;
            //
            dummy->Lind_devia = givenautoin_devia;
        } //para el fieldtotl no tento en cuenta mas que la autoin y no devuelvo nada !100517

        wiresconstantes(control->fieldtotl, dummy, G2, sgg);

    }

    //!!!
    //!!copiado 23/04/2014
    //
    if (trim(adjustl(control->wiresflavor)) == 'transition') {
        //contar multilines
        NumMultilines = 0;
        for (is1 = 1; is1 <= HWires->NumCurrentSegments; is1++) {
            if (!HWires->CurrentSegment[is1]->proc) {
                NumMultilines++;
                for (is2 = is1 + 1; is2 <= HWires->NumCurrentSegments; is2++) {
                    if ((HWires->CurrentSegment[is1]->i == HWires->CurrentSegment[is2]->i) &&
                        (HWires->CurrentSegment[is1]->j == HWires->CurrentSegment[is2]->j) &&
                        (HWires->CurrentSegment[is1]->k == HWires->CurrentSegment[is2]->k) &&
                        (HWires->CurrentSegment[is1]->tipofield == HWires->CurrentSegment[is2]->tipofield) &&
                        !HWires->CurrentSegment[is1]->proc) {

                        HWires->CurrentSegment[is2]->proc = true;
                    }
                }
            }
        }

        for (is1 = 1; is1 <= HWires->NumCurrentSegments; is1++) {
            HWires->CurrentSegment[is1]->proc = false;
        }

        HWires->NumMultilines = NumMultilines;
        write(buff, *) = 'wir1_INFO: Numero de multilineas ' + NumMultilines;
        if (control->verbose) {
            WarnErrReport(buff);
        }
        allocate(HWires->Multilines, 1, NumMultilines);
        //contar paralelos
        contmtln = 0;
        for (is1 = 1; is1 <= HWires->NumCurrentSegments; is1++) {
            if (!HWires->CurrentSegment[is1]->proc) {
                contmtln++;
                if (contmtln > NumMultilines) {
                    write(buff, *) = 'wir0_BUGGYERROR: Demasiados multihilos';
                    if ((k1 >= ZI) && (k1 <= ZE)) {
                        WarnErrReport(buff, true);
                    }
                    ThereAreWires = false;
                }
                NumParallel = 1;
                for (is2 = is1 + 1; is2 <= HWires->NumCurrentSegments; is2++) {
                    if ((HWires->CurrentSegment[is1]->i == HWires->CurrentSegment[is2]->i) &&
                        (HWires->CurrentSegment[is1]->j == HWires->CurrentSegment[is2]->j) &&
                        (HWires->CurrentSegment[is1]->k == HWires->CurrentSegment[is2]->k) &&
                        (HWires->CurrentSegment[is1]->tipofield == HWires->CurrentSegment[is2]->tipofield)) {

                        NumParallel++;
                        HWires->CurrentSegment[is2]->proc = true;
                    }
                }
                HWires->Multilines[contmtln]->NumParallel = NumParallel;
                allocate(HWires->Multilines[contmtln]->Segments, 1, NumParallel);
            }
        }

        for (is1 = 1; is1 <= HWires->NumCurrentSegments; is1++) {
            HWires->CurrentSegment[is1]->proc = false;
        }

        //asignar multilines
        contmtln = 0;
        for (is1 = 1; is1 <= HWires->NumCurrentSegments; is1++) {
            if (!HWires->CurrentSegment[is1]->proc) {
                contmtln++;
                if (contmtln > NumMultilines) {
                    write(buff, *) = 'wir0_BUGGYERROR: Demasiados multihilos';
                    if ((k1 >= ZI) && (k1 <= ZE)) {
                        WarnErrReport(buff, true);
                    }
                    ThereAreWires = false;
                }
                contprll = 0;
                for (is2 = is1; is2 <= HWires->NumCurrentSegments; is2++) {
                    if ((HWires->CurrentSegment[is1]->i == HWires->CurrentSegment[is2]->i) &&
                        (HWires->CurrentSegment[is1]->j == HWires->CurrentSegment[is2]->j) &&
                        (HWires->CurrentSegment[is1]->k == HWires->CurrentSegment[is2]->k) &&
                        (HWires->CurrentSegment[is1]->tipofield == HWires->CurrentSegment[is2]->tipofield)) {

                        contprll++;
                        if (contprll > HWires->Multilines[contmtln]->NumParallel) {
                            write(buff, *) = 'wir0_BUGGYERROR: Demasiados hilos paralelos';
                            if ((k1 >= ZI) && (k1 <= ZE)) {
                                WarnErrReport(buff, true);
                            }
                            ThereAreWires = false;
                        }
                        HWires->CurrentSegment[is2]->proc = true;
                        HWires->Multilines[contmtln]->Segments[contprll]->ptr = HWires->CurrentSegment[is2];
                    }
                }
            }
        }

        //asignamos posiciones en celda
        for (iw1 = 1; iw1 <= NumMultilines; iw1++) {
            N = HWires->Multilines[iw1]->NumParallel;
            for (i1 = 2; i1 <= N; i1++) {
                HWires->Multilines[iw1]->Segments[i1]->ptr->x = HWires->Multilines[iw1]->Segments[i1]->ptr->x + cos(static_cast<RKIND_wires>(i1 - 2) / static_cast<RKIND_wires>(N - 1) * 2.0_RKIND_wires * pi) * 0.25_RKIND_wires;
                HWires->Multilines[iw1]->Segments[i1]->ptr->y = HWires->Multilines[iw1]->Segments[i1]->ptr->y + sin(static_cast<RKIND_wires>(i1 - 2) / static_cast<RKIND_wires>(N - 1) * 2.0_RKIND_wires * pi) * 0.25_RKIND_wires;
            }
        }

        //asignar autoinducciones
        //asignar constantes de evolucion

        for (iw1 = 1; iw1 <= NumMultilines; iw1++) {

NumParallel = HWires.Multilines[iw1].NumParallel;
            HWires.Multilines[iw1].C.resize(NumParallel, std::vector<double>(NumParallel));
            HWires.Multilines[iw1].R.resize(NumParallel, std::vector<double>(NumParallel));
            HWires.Multilines[iw1].L.resize(NumParallel, std::vector<double>(NumParallel));
            
            for (int is1 = 1; is1 <= NumParallel; ++is1) {
                dl = HWires.Multilines[iw1].Segments[is1].ptr->delta;
                dx1 = HWires.Multilines[iw1].Segments[is1].ptr->deltaTransv1;
                dx2 = HWires.Multilines[iw1].Segments[is1].ptr->deltaTransv2;
                r0 = HWires.Multilines[iw1].Segments[is1].ptr->TipoWire->radius;

                if (r0 < 1e-30) {
                    std::ostringstream buff_stream;
                    buff_stream << "wir0_ERROR: Wire radius cannot be null";
                    if ((k1 >= ZI) && (k1 <= ZE)) {
                        WarnErrReport(buff_stream.str(), true);
                    }
                }

                imed = HWires.Multilines[iw1].Segments[is1].ptr->indexmed;
                
                // Calculate Lintrinsic
                double log_term = std::log((std::pow(dx1, 2.0) + std::pow(dx2, 2.0)) / (4.0 * r0 * r0));
                double atan_term1 = (dx1 / dx2) * std::atan(dx2 / dx1);
                double atan_term2 = (dx2 / dx1) * std::atan(dx1 / dx2);
                double area_term = (M_PI * r0 * r0) / (dx2 * dx1);
                
                HWires.Multilines[iw1].Segments[is1].ptr->Lintrinsic = 
                    (1.0 / (4.0 * M_PI * InvMu(imed))) * (log_term + atan_term1 + atan_term2 + area_term - 3.0);

                if ((r0 < 0.3 * dx1) || (r0 < 0.3 * dx2)) {
                    HWires.Multilines[iw1].Segments[is1].ptr->Lintrinsic -= 0.57 / (4.0 * M_PI * InvMu(imed));
                }

                if ((r0 > 0.3 * dx1) || (r0 > 0.3 * dx2)) {
                    HWires.Multilines[iw1].Segments[is1].ptr->Lintrinsic /= (1.0 - (M_PI * r0 * r0) / (dx2 * dx1));
                }

                HWires.Multilines[iw1].Segments[is1].ptr->L += HWires.Multilines[iw1].Segments[is1].ptr->Lintrinsic;
                HWires.Multilines[iw1].Segments[is1].ptr->C += 1.0 / (InvMu(imed) * InvEps(imed) * HWires.Multilines[iw1].Segments[is1].ptr->Lintrinsic);

                for (int is2 = 1; is2 <= NumParallel; ++is2) {
                    if (is1 == is2) {
                        HWires.Multilines[iw1].C[is1-1][is2-1] = HWires.Multilines[iw1].Segments[is1].ptr->C;
                        HWires.Multilines[iw1].R[is1-1][is2-1] = HWires.Multilines[iw1].Segments[is1].ptr->R;
                        HWires.Multilines[iw1].L[is1-1][is2-1] = HWires.Multilines[iw1].Segments[is1].ptr->L;
                    } else {
                        HWires.Multilines[iw1].C[is1-1][is2-1] = 0.0;
                        HWires.Multilines[iw1].R[is1-1][is2-1] = 0.0;

                        double dx = (HWires.Multilines[iw1].Segments[is1].ptr->x - HWires.Multilines[iw1].Segments[is2].ptr->x) * HWires.Multilines[iw1].Segments[is2].ptr->deltaTransv1;
                        double dy = (HWires.Multilines[iw1].Segments[is1].ptr->y - HWires.Multilines[iw1].Segments[is2].ptr->y) * HWires.Multilines[iw1].Segments[is2].ptr->deltaTransv2;
                        dist = std::sqrt(std::pow(dx, 2.0) + std::pow(dy, 2.0));
                        
                        phi = std::atan2((HWires.Multilines[iw1].Segments[is1].ptr->y - HWires.Multilines[iw1].Segments[is2].ptr->y), 
                                         (HWires.Multilines[iw1].Segments[is1].ptr->x - HWires.Multilines[iw1].Segments[is2].ptr->x));
                        
                        HWires.Multilines[iw1].L[is1-1][is2-1] = (1.0 / (2.0 * M_PI * InvMu(imed))) * 
                            F(HWires.Multilines[iw1].Segments[is1].ptr->deltaTransv1 / 2.0, 
                              HWires.Multilines[iw1].Segments[is1].ptr->deltaTransv2 / 2.0,
                              HWires.Multilines[iw1].Segments[is1].ptr->TipoWire->radius, 
                              HWires.Multilines[iw1].Segments[is2].ptr->TipoWire->radius,
                              dist, phi);
                    }
                }
            }

            std::ostringstream buff_stream2;
            buff_stream2 << "wir1_INFO: Multihilo " << iw1 << " con numero de hilos paralelos " << NumParallel;
            if (control.verbose) {
                WarnErrReport(buff_stream2.str());
            }

            HWires.Multilines[iw1].b1I.resize(NumParallel, std::vector<double>(NumParallel));
            HWires.Multilines[iw1].b2I.resize(NumParallel, std::vector<double>(NumParallel));
            HWires.Multilines[iw1].b3I.resize(NumParallel, std::vector<double>(NumParallel));
            
            std::vector<std::vector<double>> Den(NumParallel, std::vector<double>(NumParallel));
            for(int r=0; r<NumParallel; ++r) {
                for(int c=0; c<NumParallel; ++c) {
                    Den[r][c] = HWires.Multilines[iw1].L[r][c] + HWires.Multilines[iw1].R[r][c] * sgg.dt / 2.0;
                }
            }
            
            MatInv(NumParallel, Den);

            // MatMul(Den, HWires.Multilines[iw1].L - HWires.Multilines[iw1].R * sgg.dt / 2.0)
            std::vector<std::vector<double>> RHS(NumParallel, std::vector<double>(NumParallel));
            for(int r=0; r<NumParallel; ++r) {
                for(int c=0; c<NumParallel; ++c) {
                    RHS[r][c] = HWires.Multilines[iw1].L[r][c] - HWires.Multilines[iw1].R[r][c] * sgg.dt / 2.0;
                }
            }
            HWires.Multilines[iw1].b1I = MatMul(Den, RHS);

            // yo lo hago con cargas y no hay que multiplicar por C
            dl = HWires.Multilines[iw1].Segments[1].ptr->delta; // tomo el de 1 !dama no lo tenia
            imed = HWires.Multilines[iw1].Segments[1].ptr->indexmed; // lo aniado yo pq no estaba definido sgg 141118 !tomo el de 1 !dama no lo tenia
            
            // HWires.Multilines[iw1].b2I = -sgg.dt/dl * InvMu(imed) * InvEps(imed) * MatMul(Den, HWires.Multilines[iw1].L)
            std::vector<std::vector<double>> MulDenL = MatMul(Den, HWires.Multilines[iw1].L);
            for(int r=0; r<NumParallel; ++r) {
                for(int c=0; c<NumParallel; ++c) {
                    HWires.Multilines[iw1].b2I[r][c] = -sgg.dt / dl * InvMu(imed) * InvEps(imed) * MulDenL[r][c];
                }
            }
            
            // HWires.Multilines[iw1].b3I = sgg.dt * Den
            for(int r=0; r<NumParallel; ++r) {
                for(int c=0; c<NumParallel; ++c) {
                    HWires.Multilines[iw1].b3I[r][c] = sgg.dt * Den[r][c];
                }
            }

            // Den is local, no need to explicitly deallocate in C++ if using vector, but logically it goes out of scope.

            for (int is1 = 1; is1 <= NumParallel; ++is1) {
                HWires.Multilines[iw1].Segments[is1].ptr->bI = HWires.Multilines[iw1].b2I[is1-1][is1-1];
            }
        }
    } // end if !DEL FLAVOR transition
    // !!!!!!!!!!fin dama
    // !!!

    // voltage sources
    for (int i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
        Segmento* segmento = &HWires.CurrentSegment[i1];
        if (segmento->TipoWire->VsourceExists) {
            for (int k1 = 1; k1 <= segmento->Tipowire->numvoltagesources; ++k1) {
                if ((segmento->i == segmento->Tipowire->Vsource[k1].I) &&
                    (segmento->j == segmento->Tipowire->Vsource[k1].J) &&
                    (segmento->k == segmento->Tipowire->Vsource[k1].K)) {
                    
                    segmento->Vsource = &segmento->TipoWire->Vsource[k1];
                    segmento->HasVsource = true;
                    thereareVsources = true;

                    segmento->resist += segmento->TipoWire->vsource[k1].resistance / segmento->delta;

                    if ((segmento->HasSeries_RightEnd) || (segmento->HasSeries_LeftEnd)) {
                        std::ostringstream buff_stream3;
                        buff_stream3 << "wir1_INFO: Voltage source with Series RL (C neglected if present) impedance in segment "
                                     << segmento->origIndex << " " << segmento->i << " " << segmento->j << " " << segmento->k << " " << segmento->tipofield;
                        if ((segmento->k >= ZI) && (segmento->k <= ZE) && control.verbose) {
                            WarnErrReport(buff_stream3.str());
                        }
                    } else if ((segmento->HasParallel_RightEnd) || (segmento->HasParallel_LeftEnd)) {
                        std::ostringstream buff_stream4;
                        buff_stream4 << "wir1_WARNING: Voltage source with Parallel RLC  in segment "
                                     << segmento->origIndex << " " << segmento->i << " " << segmento->j << " " << segmento->k << " " << segmento->tipofield;
                        if ((segmento->k >= ZI) && (segmento->k <= ZE)) {
                            WarnErrReport(buff_stream4.str());
                        }
                    } else if ((segmento->HasAbsorbing_RightEnd) || (segmento->HasAbsorbing_LeftEnd)) {
                        std::ostringstream buff_stream5;
                        buff_stream5 << "wir1_WARNING: Voltage source with Absorbing  in segment "
                                     << segmento->origIndex << " " << segmento->i << " " << segmento->j << " " << segmento->k << " " << segmento->tipofield;
                        if ((segmento->k >= ZI) && (segmento->k <= ZE)) {
                            WarnErrReport(buff_stream5.str());
                        }
                    } else {
                        std::ostringstream buff_stream6;
                        buff_stream6 << "wir1_INFO: Voltage source with null internal resistence in segment "
                                     << segmento->origIndex << " " << segmento->i << " " << segmento->j << " " << segmento->k << " " << segmento->tipofield;
                        if ((segmento->k >= ZI) && (segmento->k <= ZE) && control.verbose) {
                            WarnErrReport(buff_stream6.str());
                        }
                    }

                    if (segmento->Vsource->Fichero.DeltaSamples > sgg.dt) {
                        std::ostringstream buff_stream7;
                        buff_stream7 << "wir1_WARNING: " << segmento->Vsource->Fichero.Name 
                                     << " undersampled by a factor " 
                                     << segmento->Vsource->Fichero.DeltaSamples / sgg.dt;
                        if ((segmento->k >= ZI) && (segmento->k <= ZE)) {
                            WarnErrReport(buff_stream7.str());
                        }
                    }
                }
            }
        }
    }

    // END OF THE PROCESSING OF CURRENT SEGMENT INFORMATION

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    // START OF CHARGE NODAL PROCESSING
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    
    HWires.NumChargeNodes = 2 * HWires.NumCurrentSegments;
    // There can be up to twice the ammount of segments (this allocation could be reduced later -I don't do it since memory is enough-)
    HWires.ChargeNode.resize(HWires.NumChargeNodes);
    
    for (int i1 = 1; i1 <= HWires.NumChargeNodes; ++i1) {
        // Nullify pointers
        HWires.ChargeNode[i1-1].CurrentPlus_1 = nullptr;
        HWires.ChargeNode[i1-1].CurrentPlus_2 = nullptr;
        HWires.ChargeNode[i1-1].CurrentPlus_3 = nullptr;
        HWires.ChargeNode[i1-1].CurrentPlus_4 = nullptr;
        HWires.ChargeNode[i1-1].CurrentPlus_5 = nullptr;
        HWires.ChargeNode[i1-1].CurrentPlus_6 = nullptr;
        HWires.ChargeNode[i1-1].CurrentPlus_7 = nullptr;
        HWires.ChargeNode[i1-1].CurrentPlus_8 = nullptr;
        HWires.ChargeNode[i1-1].CurrentPlus_9 = nullptr;
        
        HWires.ChargeNode[i1-1].CurrentMinus_1 = nullptr;
        HWires.ChargeNode[i1-1].CurrentMinus_2 = nullptr;
        HWires.ChargeNode[i1-1].CurrentMinus_3 = nullptr;
        HWires.ChargeNode[i1-1].CurrentMinus_4 = nullptr;
        HWires.ChargeNode[i1-1].CurrentMinus_5 = nullptr;
        HWires.ChargeNode[i1-1].CurrentMinus_6 = nullptr;
        HWires.ChargeNode[i1-1].CurrentMinus_7 = nullptr;
        HWires.ChargeNode[i1-1].CurrentMinus_8 = nullptr;
        HWires.ChargeNode[i1-1].CurrentMinus_9 = nullptr;
        
        HWires.ChargeNode[i1-1].NodeInside = nullptr;

        HWires.ChargeNode[i1-1].indexnode = -1;
        HWires.ChargeNode[i1-1].ChargePresent = 0.0;
        HWires.ChargeNode[i1-1].ChargePast = 0.0;
        HWires.ChargeNode[i1-1].IsAttachedtoVoltage = false;
        HWires.ChargeNode[i1-1].IsMur = false;
        HWires.ChargeNode[i1-1].IsBackDownLeftMur = false;
        HWires.ChargeNode[i1-1].IsFrontUpRightMur = false;
        HWires.ChargeNode[i1-1].IsPeriodic = false;
        HWires.ChargeNode[i1-1].IsPEC = false;
        
        HWires.ChargeNode[i1-1].already_YEEadvanced_byconformal_changedtoPECfield1 = nullptr;
        HWires.ChargeNode[i1-1].already_YEEadvanced_byconformal_changedtoPECfield2 = nullptr;
        HWires.ChargeNode[i1-1].already_YEEadvanced_byconformal_changedtoPECfield3 = nullptr;
        HWires.ChargeNode[i1-1].already_YEEadvanced_byconformal_changedtoPECfield4 = nullptr;
        HWires.ChargeNode[i1-1].already_YEEadvanced_byconformal_changedtoPECfield5 = nullptr;
        HWires.ChargeNode[i1-1].already_YEEadvanced_byconformal_changedtoPECfield6 = nullptr;
        
        HWires.ChargeNode[i1-1].IsLossy = false;
        HWires.ChargeNode[i1-1].HasIsource = false;
        HWires.ChargeNode[i1-1].IsHeterogeneousJunction = false;
        
        // just for informative !not implemented in MPI unsure behaviour under mpi 2011 \E7
        HWires.ChargeNode[i1-1].Exists = false;
        HWires.ChargeNode[i1-1].Is_LeftEnd = false;
        HWires.ChargeNode[i1-1].Is_RightEnd = false;
        HWires.ChargeNode[i1-1].IsInSingleRLCsegment = false;
        HWires.ChargeNode[i1-1].NumCurrentMinus = 0;
        HWires.ChargeNode[i1-1].NumCurrentPlus = 0;
        HWires.ChargeNode[i1-1].i = -1;
        HWires.ChargeNode[i1-1].j = -1;
        HWires.ChargeNode[i1-1].k = -1;
        HWires.ChargeNode[i1-1].cteMur = 0.0;
        HWires.ChargeNode[i1-1].cteProp = 0.0;
        HWires.ChargeNode[i1-1].oRIGctePlain = 0.0;
        HWires.ChargeNode[i1-1].ctePlain = 0.0;
        
        for (int j1 = 1; j1 <= 2 * MaxNumCurrentMinusPlus; ++j1) {
            HWires.ChargeNode[i1-1].YESsegment[j1-1] = -j1; // default !servira luego para contar nodos repetidos
        }
    }

    std::ostringstream buff_stream8;
    buff_stream8 << "----------------------------------------------------------------";
    WarnErrReport(buff_stream8.str());
    
    // Detect nodes which are Between two segments (posible duplicities are removed)
    conta = 0;
    NUMESEG = HWires.NumCurrentSegments;
    // runtina de adjacencias
    for (int i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
        Segmento* org = &HWires.CurrentSegment[i1];
        Segmento* orgmenos1 = nullptr;
        Segmento* orgmas1 = nullptr;
        
        if (i1 > 1) {
            if (HWires.CurrentSegment[i1-1].indexmed == HWires.CurrentSegment[i1].indexmed) {
                orgmenos1 = &HWires.CurrentSegment[i1-1];
            }
        }
        if (i1 < NUMESEG) {
            if (HWires.CurrentSegment[i1+1].indexmed == HWires.CurrentSegment[i1].indexmed) {
                orgmas1 = &HWires.CurrentSegment[i1+1];
            }
        }

        for (int j1 = 1; j1 <= HWires.NumCurrentSegments; ++j1) {
            Segmento* fin = &HWires.CurrentSegment[j1];
            Segmento* FINmenos1 = nullptr;
            Segmento* FINmas1 = nullptr;
            
            if (j1 > 1) {
                FINmenos1 = &HWires.CurrentSegment[j1-1];
            }
            if (j1 < NUMESEG) {
                FINmas1 = &HWires.CurrentSegment[j1+1];
            }

            if (i1 != j1) {
                // calls the adjacent routine to see wether the org and fin segments are connected
                AdjacencyResult Adj = TestAdjacency(org, i1, fin, j1, control.connectendings, control.isolategroupgroups, control.strictOLD, ZI, ZE, NUMESEG, orgmenos1, orgmas1, FINmenos1, FINmas1, control.verbose);

                if (Adj.Is) {
                    bool NodeExists = false;
                    for (int nn = 1; nn <= conta; ++nn) {
                        if ((HWires.ChargeNode[nn-1].i == Adj.i) &&
                            (HWires.ChargeNode[nn-1].j == Adj.j) &&
                            (HWires.ChargeNode[nn-1].k == Adj.k)) {
                            
                            if ((Adj.IsHeterogeneousJunction && HWires.ChargeNode[nn-1].IsHeterogeneousJunction) ||
                                ((!Adj.IsHeterogeneousJunction) && (!HWires.ChargeNode[nn-1].IsHeterogeneousJunction))) {
                                
                                if (Adj.YESsegment[0] != -1) {
                                    bool proceed1 = true;
                                    bool proceed2 = true;
                                    
                                    for (int j2 = 1; j2 <= 2 * MaxNumCurrentMinusPlus; ++j2) {
                                        if (HWires.ChargeNode[nn-1].YESsegment[j2-1] == Adj.YESsegment[0]) {
                                            proceed1 = false;
                                            break;
                                        }
                                    }
                                    
                                    for (int j2 = 1; j2 <= 2 * MaxNumCurrentMinusPlus; ++j2) {
                                        if (HWires.ChargeNode[nn-1].YESsegment[j2-1] == Adj.YESsegment[1]) {
                                            proceed2 = false;
                                            break;
                                        }
                                    }
                                    
                                    if ((!proceed1) || (!proceed2) || (!(proceed1 && proceed2))) {
                                        // This node was already present. There is connectivity and the node actually exists
                                        NodeExists = true;
                                        if (proceed1) {
                                            for (int j2 = 1; j2 <= 2 * MaxNumCurrentMinusPlus; ++j2) {
                                                if (HWires.ChargeNode[nn-1].YESsegment[j2-1] < 0) {
                                                    HWires.ChargeNode[nn-1].YESsegment[j2-1] = Adj.YESsegment[0];
                                                    break;
                                                }
                                            }
                                        }
                                        if (proceed2) {
                                            for (int j2 = 1; j2 <= 2 * MaxNumCurrentMinusPlus; ++j2) {
                                                if (HWires.ChargeNode[nn-1].YESsegment[j2-1] < 0) {
                                                    HWires.ChargeNode[nn-1].YESsegment[j2-1] = Adj.YESsegment[1];
                                                    break;
                                                }
                                            }
                                        }
                                        goto bprin_end;
                                    } else {
                                        NodeExists = false;
                                    }
                                }
                            }
                        }
                    }
                    bprin_end:;
                    
                    if (!NodeExists) {
                        conta++;
                        if (Adj.YESsegment[0] != -1) {
                            for (int j2 = 1; j2 <= 2 * MaxNumCurrentMinusPlus; ++j2) {
                                if (HWires.ChargeNode[conta-1].YESsegment[j2-1] < 0) {
                                    HWires.ChargeNode[conta-1].YESsegment[j2-1] = Adj.YESsegment[0];
                                    break;
                                }
                            }
                            for (int j2 = 1; j2 <= 2 * MaxNumCurrentMinusPlus; ++j2) {
                                if (HWires.ChargeNode[conta-1].YESsegment[j2-1] < 0) {
                                    HWires.ChargeNode[conta-1].YESsegment[j2-1] = Adj.YESsegment[1];
                                    break;
                                }
                            }
                        }
                        
                        HWires.ChargeNode[conta-1].IsHeterogeneousJunction = Adj.IsHeterogeneousJunction;
                        HWires.ChargeNode[conta-1].Exists = true;
                        HWires.ChargeNode[conta-1].i = Adj.i;
                        HWires.ChargeNode[conta-1].j = Adj.j;
                        HWires.ChargeNode[conta-1].k = Adj.k;
                    }
                }

                if (Adj.Is && Adj.BothEndingsConnected) { // busca el otro extremo y crealo tambien !SOLO NECESARIO EN WIRES DE UN SOLO SEGMENTO
                    NodeExists = false;
                    
                    if (Adj.i != org->i) {
                        Adj.i = org->i;
                    } else if (Adj.i != org->ilibre) {
                        Adj.i = org->ilibre;
                    } else {
                        Adj.i = org->i;
                    }
                    
                    if (Adj.J != org->J) {
                        Adj.J = org->J;
                    } else if (Adj.J != org->Jlibre) {
                        Adj.J = org->Jlibre;
                    } else {
                        Adj.J = org->J;
                    }
                    
                    if (Adj.K != org->K) {
                        Adj.K = org->K;
                    } else if (Adj.K != org->Klibre) {
                        Adj.K = org->Klibre;
                    } else {
                        Adj.K = org->K;
                    }

                    for (int nn = 1; nn <= conta; ++nn) {
                        if ((HWires.ChargeNode[nn-1].i == Adj.i) &&
                            (HWires.ChargeNode[nn-1].j == Adj.j) &&
                            (HWires.ChargeNode[nn-1].k == Adj.k)) {
                            
                            if ((Adj.IsHeterogeneousJunction && HWires.ChargeNode[nn-1].IsHeterogeneousJunction) ||
                                ((!Adj.IsHeterogeneousJunction) && (!HWires.ChargeNode[nn-1].IsHeterogeneousJunction))) {
                                
                                if (Adj.YESsegment[0] != -1) {
                                    bool proceed1 = true;
                                    bool proceed2 = true;
                                    
                                    for (int j2 = 1; j2 <= 2 * MaxNumCurrentMinusPlus; ++j2) {
                                        if (HWires.ChargeNode[nn-1].YESsegment[j2-1] == Adj.YESsegment[0]) {
                                            proceed1 = false;
                                            break;
                                        }
                                    }

}
                        }
                        for (int j2 = 1; j2 <= 2 * MaxNumCurrentMinusPlus; ++j2) {
                            if (HWires.ChargeNode(nn).YESsegment(j2) == adj.YESsegment(2)) {
                                proceed2 = false;
                                goto b112_exit;
                            }
                        }
                        b112_exit:;
                        if ((!proceed1) || (!proceed2) || (!proceed1 || !proceed2)) {
                            // This node was already present. There is connectivity and the node actually exists
                            NodeExists = true;
                            if (proceed1) {
                                for (int j2 = 1; j2 <= 2 * MaxNumCurrentMinusPlus; ++j2) {
                                    if (HWires.ChargeNode(nn).YESsegment(j2) < 0) {
                                        HWires.ChargeNode(nn).YESsegment(j2) = adj.YESsegment(1);
                                        goto b22_exit;
                                    }
                                }
                                b22_exit:;
                            }
                            if (PROCEED2) {
                                for (int j2 = 1; j2 <= 2 * MaxNumCurrentMinusPlus; ++j2) {
                                    if (HWires.ChargeNode(nn).YESsegment(j2) < 0) {
                                        HWires.ChargeNode(nn).YESsegment(j2) = adj.YESsegment(2);
                                        goto b42_exit;
                                    }
                                }
                                b42_exit:;
                            }
                            goto bprin2_exit;
                        } else {
                            NodeExists = false;
                        }
                    }
                }
            }
        }
        bprin2_exit:;
        if (!NodeExists) {
            conta = conta + 1;
            if (adj.YESsegment(1) != -1) {
                for (int j2 = 1; j2 <= 2 * MaxNumCurrentMinusPlus; ++j2) {
                    if (HWires.ChargeNode(conta).YESsegment(j2) < 0) {
                        HWires.ChargeNode(conta).YESsegment(j2) = adj.YESsegment(1);
                        goto b62_exit;
                    }
                }
                b62_exit:;
                for (int j2 = 1; j2 <= 2 * MaxNumCurrentMinusPlus; ++j2) {
                    if (HWires.ChargeNode(conta).YESsegment(j2) < 0) {
                        HWires.ChargeNode(conta).YESsegment(j2) = adj.YESsegment(2);
                        goto b82_exit;
                    }
                }
                b82_exit:;
            }
            //
            HWires.ChargeNode(conta).IsHeterogeneousJunction = adj.IsHeterogeneousJunction;
            HWires.ChargeNode(conta).Exists = true;
            HWires.ChargeNode(conta).i = adj.i;
            HWires.ChargeNode(conta).j = adj.j;
            HWires.ChargeNode(conta).k = adj.k;
        }
    }
} // del i1/=ji1
//
    }
}
HWires.NumChargeNodes = conta;
// reactualize the number of nodes


std::string buff = "----------------------------------------------------------------";
WarnErrReport(buff);

// Now I remove duplicate nodes, using info of the shared segments
// More than one of such node may exist if there are commond adjacencies
// They are not sequentially detected
for (int nn = 1; nn <= HWires.NumChargeNodes; ++nn) {
    for (int nnn = nn + 1; nnn <= HWires.NumChargeNodes; ++nnn) {
        if ((HWires.ChargeNode(nn).i == HWires.ChargeNode(nnn).i) &&
            (HWires.ChargeNode(nn).j == HWires.ChargeNode(nnn).j) &&
            (HWires.ChargeNode(nn).k == HWires.ChargeNode(nnn).k)) {
            conta = 0;
            for (int j1 = 1; j1 <= 2 * MaxNumCurrentMinusPlus; ++j1) {
                for (int k1 = 1; k1 <= 2 * MaxNumCurrentMinusPlus; ++k1) {
                    if (HWires.ChargeNode(nn).YESsegment(k1) == HWires.ChargeNode(nnn).YESsegment(j1)) {
                        conta = conta + 1;
                    }
                }
            }
            if (conta == 2 * MaxNumCurrentMinusPlus) { // They coincide
                HWires.ChargeNode(nn).Exists = false; // voids it
            }
        }
    }
}

// compress the authentic nodes
conta = HWires.NumChargeNodes;
for (int nn = 1; nn <= HWires.NumChargeNodes; ++nn) {
    if (!HWires.ChargeNode(nn).Exists) {
        conta = conta - 1;
        for (int nnn = nn; nnn <= HWires.NumChargeNodes - 1; ++nnn) {
            HWires.ChargeNode(nnn) = HWires.ChargeNode(nnn + 1);
        }
        continue;
    }
}
// voids the rest 
for (int i1 = conta + 1; i1 <= HWires.NumChargeNodes; ++i1) {

    HWires.ChargeNode(i1).CurrentPlus_1 = nullptr; HWires.ChargeNode(i1).CurrentPlus_2 = nullptr; HWires.ChargeNode(i1).CurrentPlus_3 = nullptr;
    HWires.ChargeNode(i1).CurrentPlus_4 = nullptr; HWires.ChargeNode(i1).CurrentPlus_5 = nullptr; HWires.ChargeNode(i1).CurrentPlus_6 = nullptr;
    HWires.ChargeNode(i1).CurrentPlus_7 = nullptr; HWires.ChargeNode(i1).CurrentPlus_8 = nullptr; HWires.ChargeNode(i1).CurrentPlus_9 = nullptr;
    HWires.ChargeNode(i1).CurrentMinus_1 = nullptr; HWires.ChargeNode(i1).CurrentMinus_2 = nullptr; HWires.ChargeNode(i1).CurrentMinus_3 = nullptr;
    HWires.ChargeNode(i1).CurrentMinus_4 = nullptr; HWires.ChargeNode(i1).CurrentMinus_5 = nullptr; HWires.ChargeNode(i1).CurrentMinus_6 = nullptr;
    HWires.ChargeNode(i1).CurrentMinus_7 = nullptr; HWires.ChargeNode(i1).CurrentMinus_8 = nullptr; HWires.ChargeNode(i1).CurrentMinus_9 = nullptr;
    HWires.ChargeNode(i1).NodeInside = nullptr;
    HWires.ChargeNode(i1).ChargePresent = 0.0_RKIND_wires;
    HWires.ChargeNode(i1).ChargePast = 0.0_RKIND_wires;
    HWires.ChargeNode(i1).IsAttachedtoVoltage = false;
    HWires.ChargeNode(i1).IsMur = false;
    HWires.ChargeNode(i1).IsBackDownLeftMur = false;
    HWires.ChargeNode(i1).IsFrontUpRightMur = false;
    HWires.ChargeNode(i1).Isperiodic = false;
    HWires.ChargeNode(i1).IsPEC = false;
    HWires.ChargeNode(i1).already_YEEadvanced_byconformal_changedtoPECfield1 = nullptr;
    HWires.ChargeNode(i1).already_YEEadvanced_byconformal_changedtoPECfield2 = nullptr;
    HWires.ChargeNode(i1).already_YEEadvanced_byconformal_changedtoPECfield3 = nullptr;
    HWires.ChargeNode(i1).already_YEEadvanced_byconformal_changedtoPECfield4 = nullptr;
    HWires.ChargeNode(i1).already_YEEadvanced_byconformal_changedtoPECfield5 = nullptr;
    HWires.ChargeNode(i1).already_YEEadvanced_byconformal_changedtoPECfield6 = nullptr;
    HWires.ChargeNode(i1).IsLossy = false;
    HWires.ChargeNode(i1).HasIsource = false;
    HWires.ChargeNode(i1).IsHeterogeneousJunction = false;
    HWires.ChargeNode(i1).Exists = false;
    HWires.ChargeNode(i1).Is_LeftEnd = false;
    HWires.ChargeNode(i1).Is_RightEnd = false;
    HWires.ChargeNode(i1).NumCurrentMinus = 0;
    HWires.ChargeNode(i1).NumCurrentPlus = 0;
    HWires.ChargeNode(i1).i = -1;
    HWires.ChargeNode(i1).j = -1;
    HWires.ChargeNode(i1).k = -1;
    HWires.ChargeNode(i1).cteMur = 0.0_RKIND_wires;
    HWires.ChargeNode(i1).cteProp = 0.0_RKIND_wires;
    HWires.ChargeNode(i1).oRIGctePlain = 0.0_RKIND_wires;
    HWires.ChargeNode(i1).ctePlain = 0.0_RKIND_wires;
    for (int j1 = 1; j1 <= 2 * MaxNumCurrentMinusPlus; ++j1) {
        HWires.ChargeNode(i1).YESsegment(j1) = -j1; // default !used for duplicate nodes
    }
}
//
HWires.NumChargeNodes = conta; // reactualize the number of authentic nodes


// I point segments nodal info to the nodes connected to each segment
for (int i1 = 1; i1 <= HWires.NumChargeNodes; ++i1) {
    ChargeNode& nodo = HWires.ChargeNode(i1);
    for (int j1 = 1; j1 <= HWires.NumCurrentSegments; ++j1) {
        CurrentSegment& segmento = HWires.CurrentSegment(j1);
        bool proceed = false;
        for (int k1 = 1; k1 <= 2 * MaxNumCurrentMinusPlus; ++k1) {
            if (j1 == nodo.YESsegment(k1)) {
                proceed = true;
                goto total_exit;
            }
        }
        total_exit:;
        if (proceed) {
            if ((segmento.i == nodo.i) && (segmento.j == nodo.j) && (segmento.k == nodo.k)) {
                segmento.ChargeMinus = &HWires.ChargeNode(i1);
            }
            if ((segmento.i + 1 == nodo.i) && (segmento.j == nodo.j) && (segmento.k == nodo.k) && (segmento.tipofield == iEx)) {
                segmento.ChargePlus = &HWires.ChargeNode(i1);
            }
            if ((segmento.i == nodo.i) && (segmento.j + 1 == nodo.j) && (segmento.k == nodo.k) && (segmento.tipofield == iEy)) {
                segmento.ChargePlus = &HWires.ChargeNode(i1);
            }
            if ((segmento.i == nodo.i) && (segmento.j == nodo.j) && (segmento.k + 1 == nodo.k) && (segmento.tipofield == iEz)) {
                segmento.ChargePlus = &HWires.ChargeNode(i1);
            }
        } else {
            continue;
        }
    }
}

// I point nodes to their segments, detecting possible junctions (I set the IsJunction Flag)
// Beware that I set the IsJunction flag even in the case of TWO segments which are both in the plus
// or in the minus direction (later I set another flag for the authentic junctions (up to 2*MaxNumCurrentMinusPlus) segments
for (int j1 = 1; j1 <= HWires.NumCurrentSegments; ++j1) {
    CurrentSegment& segmento = HWires.CurrentSegment(j1);
    if (segmento.ChargePlus != nullptr) {
        segmento.ChargePlus->NumCurrentMinus = segmento.ChargePlus->NumCurrentMinus + 1;
        if (segmento.ChargePlus->NumCurrentMinus > 9) {
            std::string BUFF = "wir1_ERROR: More than 9 Minus WIREs joined";
            WarnErrReport(BUFF, true);
        }
        switch (segmento.ChargePlus->NumCurrentMinus) {
            case 1: segmento.ChargePlus->CurrentMinus_1 = &HWires.CurrentSegment(j1); break;
            case 2: segmento.ChargePlus->CurrentMinus_2 = &HWires.CurrentSegment(j1); break;
            case 3: segmento.ChargePlus->CurrentMinus_3 = &HWires.CurrentSegment(j1); break;
            case 4: segmento.ChargePlus->CurrentMinus_4 = &HWires.CurrentSegment(j1); break;
            case 5: segmento.ChargePlus->CurrentMinus_5 = &HWires.CurrentSegment(j1); break;
            case 6: segmento.ChargePlus->CurrentMinus_6 = &HWires.CurrentSegment(j1); break;
            case 7: segmento.ChargePlus->CurrentMinus_7 = &HWires.CurrentSegment(j1); break;
            case 8: segmento.ChargePlus->CurrentMinus_8 = &HWires.CurrentSegment(j1); break;
            case 9: segmento.ChargePlus->CurrentMinus_9 = &HWires.CurrentSegment(j1); break;
        }
    }
    if (segmento.ChargeMinus != nullptr) {
        segmento.ChargeMinus->NumCurrentPlus = segmento.ChargeMinus->NumCurrentPlus + 1;
        if (segmento.ChargeMinus->NumCurrentPlus > 9) {
            std::string BUFF = "wir1_ERROR: More than 9 Plus WIREs joined";
            WarnErrReport(BUFF, true);
        }
        switch (segmento.ChargeMinus->NumCurrentPlus) {
            case 1: segmento.ChargeMinus->CurrentPlus_1 = &HWires.CurrentSegment(j1); break;
            case 2: segmento.ChargeMinus->CurrentPlus_2 = &HWires.CurrentSegment(j1); break;
            case 3: segmento.ChargeMinus->CurrentPlus_3 = &HWires.CurrentSegment(j1); break;
            case 4: segmento.ChargeMinus->CurrentPlus_4 = &HWires.CurrentSegment(j1); break;
            case 5: segmento.ChargeMinus->CurrentPlus_5 = &HWires.CurrentSegment(j1); break;
            case 6: segmento.ChargeMinus->CurrentPlus_6 = &HWires.CurrentSegment(j1); break;
            case 7: segmento.ChargeMinus->CurrentPlus_7 = &HWires.CurrentSegment(j1); break;
            case 8: segmento.ChargeMinus->CurrentPlus_8 = &HWires.CurrentSegment(j1); break;
            case 9: segmento.ChargeMinus->CurrentPlus_9 = &HWires.CurrentSegment(j1); break;
        }
    }
}

// Now I detect hanging nodes (at wire terminations) and increment the node info
for (int i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
    if (HWires.CurrentSegment(i1).ChargePlus == nullptr) {
        // segments must end in a charge node
        conta = conta + 1;
        HWires.CurrentSegment(i1).ChargePlus = &HWires.ChargeNode(conta);
        HWires.ChargeNode(conta).CurrentMinus_1 = &HWires.CurrentSegment(i1);
        HWires.ChargeNode(conta).NumCurrentMinus = 1;
        HWires.ChargeNode(conta).i = HWires.ChargeNode(conta).CurrentMinus_1->i;
        HWires.ChargeNode(conta).j = HWires.ChargeNode(conta).CurrentMinus_1->j;
        HWires.ChargeNode(conta).k = HWires.ChargeNode(conta).CurrentMinus_1->k;
        HWires.ChargeNode(conta).Exists = true;
        if (HWires.ChargeNode(conta).CurrentMinus_1->tipofield == iEx) {
            HWires.ChargeNode(conta).i = HWires.ChargeNode(conta).CurrentMinus_1->i + 1;
        }
        if (HWires.ChargeNode(conta).CurrentMinus_1->tipofield == iEy) {
            HWires.ChargeNode(conta).j = HWires.ChargeNode(conta).CurrentMinus_1->j + 1;
        }
        if (HWires.ChargeNode(conta).CurrentMinus_1->tipofield == iEz) {
            HWires.ChargeNode(conta).k = HWires.ChargeNode(conta).CurrentMinus_1->k + 1;
        }
    }
    if (HWires.CurrentSegment(i1).ChargeMinus == nullptr) {
        // segments must end in a charge node
        conta = conta + 1;
        HWires.CurrentSegment(i1).ChargeMinus = &HWires.ChargeNode(conta);
        HWires.ChargeNode(conta).CurrentPlus_1 = &HWires.CurrentSegment(i1);
        HWires.ChargeNode(conta).NumCurrentPlus = 1;
        HWires.ChargeNode(conta).i = HWires.ChargeNode(conta).CurrentPlus_1->i;
        HWires.ChargeNode(conta).j = HWires.ChargeNode(conta).CurrentPlus_1->j;
        HWires.ChargeNode(conta).k = HWires.ChargeNode(conta).CurrentPlus_1->k;
        HWires.ChargeNode(conta).Exists = true;
        if (HWires.ChargeNode(conta).CurrentPlus_1->tipofield == iEx) {
            HWires.ChargeNode(conta).i = HWires.ChargeNode(conta).CurrentPlus_1->i;
        }
        if (HWires.ChargeNode(conta).CurrentPlus_1->tipofield == iEy) {
            HWires.ChargeNode(conta).j = HWires.ChargeNode(conta).CurrentPlus_1->j;
        }
        if (HWires.ChargeNode(conta).CurrentPlus_1->tipofield == iEz) {
            HWires.ChargeNode(conta).k = HWires.ChargeNode(conta).CurrentPlus_1->k;
        }
    }
}
HWires.NumChargeNodes = conta; // reactualize POR ULTIMA VEZ the number of charge nodes
//!!!!!!!!!!


for (int i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
    CurrentSegment& segmento = HWires.CurrentSegment(i1);
    bool asignado = false;
    if (segmento.Is_LeftEnd || segmento.IsEnd_norLeft_norRight) {
        if ((segmento.chargeplus->i == segmento.ilibre) &&
            (segmento.chargeplus->j == segmento.jlibre) &&
            (segmento.chargeplus->k == segmento.klibre)) {
            segmento.chargeplus->Is_LeftEnd = true;
            asignado = true;
        }
        if ((segmento.chargeMinus->i == segmento.ilibre) &&
            (segmento.chargeMinus->j == segmento.jlibre) &&
            (segmento.chargeMinus->k == segmento.klibre)) {
            segmento.chargeMinus->Is_LeftEnd = true;
            asignado = true;
        }
    }
    if (segmento.Is_RightEnd) {
        if ((segmento.chargeplus->i == segmento.ilibre) &&
            (segmento.chargeplus->j == segmento.jlibre) &&
            (segmento.chargeplus->k == segmento.klibre)) {
            segmento.chargeplus->Is_RightEnd = true;
            asignado = true;
        }
        if ((segmento.chargeMinus->i == segmento.ilibre) &&
            (segmento.chargeMinus->j == segmento.jlibre) &&
            (segmento.chargeMinus->k == segmento.klibre)) {
            segmento.chargeMinus->Is_RightEnd = true;
            asignado = true;
        }
    }
    // caso particular hilos de 1  segmento
    if (segmento.Is_LeftEnd && segmento.Is_RightEnd) {
        segmento.chargeplus->Is_LeftEnd = true;
        segmento.chargeMinus->Is_RightEnd = true;
        asignado = true;
        // bug OLD 060313 hilos de un solo segmento con RLC. Correcion en caso de formatos antiguos solo
        if (control.connectendings) {
            segmento.chargeplus->Is_LeftEnd = true;
            segmento.chargeplus->Is_RightEnd = true;
            segmento.chargeMinus->Is_RightEnd = true;
            segmento.chargeMinus->Is_LeftEnd = true;
        }
        segmento.chargeplus->IsInSingleRLCsegment = (segmento.HasParallel_LeftEnd || segmento.HasSeries_LeftEnd) ||
            (segmento.HasParallel_RightEnd || segmento.HasSeries_RightEnd);
        segmento.chargeminus->IsInSingleRLCsegment = (segmento.HasParallel_LeftEnd || segmento.HasSeries_LeftEnd) ||
            (segmento.HasParallel_RightEnd || segmento.HasSeries_RightEnd);
    }
    if (!asignado) {
        if (segmento.HasParallel_LeftEnd || segmento.HasSeries_LeftEnd) {
            if (segmento.chargeplus->isPEC || segmento.chargeplus->isLossy || segmento.chargeplus->isLossy || segmento.chargeplus->isLossy) {
                if ((!segmento.chargeplus->Is_LeftEnd) && (!segmento.chargeplus->Is_RightEnd)) {
                    segmento.chargeplus->Is_LeftEnd = true;
                    std::ostringstream buff_stream;
                    buff_stream << "wir1_INFO: Forcing non-terminal node to LeftEnd to attach RLC "
                                << segmento.chargeplus->i << segmento.chargeplus->j << segmento.chargeplus->k;
                    std::string buff = buff_stream.str();
                    if ((segmento.chargeplus->k > ZI) && (segmento.chargeplus->k <= ZE) && control.verbose) WarnErrReport(buff);
                }
            } else if (segmento.chargeminus->isPEC || segmento.chargeminus->isLossy || segmento.chargeminus->isLossy || segmento.chargeminus->isLossy) {
                if ((!segmento.chargeminus->Is_LeftEnd) && (!segmento.chargeminus->Is_RightEnd)) {
                    segmento.chargeminus->Is_LeftEnd = true;
                    std::ostringstream buff_stream;
                    buff_stream << "wir1_INFO: Forcing non-terminal node to LeftEnd to attach RLC "
                                << segmento.chargeminus->i << segmento.chargeminus->j << segmento.chargeminus->k;
                    std::string buff = buff_stream.str();
                    if ((segmento.chargeminus->k > ZI) && (segmento.chargeminus->k <= ZE) && control.verbose) WarnErrReport(buff);
                }
            } else {
                if (segmento.TipoWire->Series_C_LeftEnd < 1.0e7_RKIND_wires) {
                    std::ostringstream buff_stream;
                    buff_stream << "wir1_ERROR: Series LeftEnd Capacitance in INTERMEDIATE segment smaller than 1e7 (inf) "
                                << segmento.origIndex << segmento.i << segmento.j << segmento.k << dir(segmento.tipofield);
                    std::string buff = buff_stream.str();
                    if ((segmento.k > ZI) && (segmento.k <= ZE) && (segmento.tipofield != iEz)) WarnErrReport(buff, true);
                    if ((segmento.k >= ZI) && (segmento.k <= ZE) && (segmento.tipofield == iEz)) WarnErrReport(buff, true);
                }
                if (segmento.TipoWire->Parallel_C_LeftEnd != 0.0_RKIND_wires) {
                    std::ostringstream buff_stream;
                    buff_stream << "wir1_ERROR: Parallel LeftEnd Capacitance in INTERMEDIATE segment "
                                << segmento.origIndex << segmento.i << segme

if ((segmento->k > ZI) && (segmento->k <= ZE) && (segmento->tipofield != iEz)) WarnErrReport(buff, true);
                    if ((segmento->k >= ZI) && (segmento->k <= ZE) && (segmento->tipofield == iEz)) WarnErrReport(buff, true);
                }
                if (segmento->TipoWire->Parallel_R_LeftEnd != 0.0_RKIND_wires) {
                    sprintf(buff, "wir1_ERROR:  Parallel LeftEnd Resistance in INTERMEDIATE segment %7i%7i%7i%7i",
                            segmento->origIndex, segmento->i, segmento->j, segmento->k);
                    if ((segmento->k > ZI) && (segmento->k <= ZE) && (segmento->tipofield != iEz)) WarnErrReport(buff, true);
                    if ((segmento->k >= ZI) && (segmento->k <= ZE) && (segmento->tipofield == iEz)) WarnErrReport(buff, true);
                }
            }
        }
        if (segmento->HasParallel_RightEnd || segmento->HasSeries_RightEnd) {
            if (segmento->chargeplus->isPEC || segmento->chargeplus->isLossy || segmento->chargeplus->isLossy || segmento->chargeplus->isLossy) {
                if ((!segmento->chargeplus->Is_LeftEnd) && (!segmento->chargeplus->Is_RightEnd)) {
                    segmento->chargeplus->Is_RightEnd = true;
                    sprintf(buff, "wir1_INFO: Forcing non-terminal node to RightEnd to attach RLC %7i%7i%7i",
                            segmento->chargeplus->i, segmento->chargeplus->j, segmento->chargeplus->k);
                    if ((segmento->chargeplus->k > ZI) && (segmento->chargeplus->k <= ZE) && control->verbose) WarnErrReport(buff);
                }
            } else if (segmento->chargeminus->isPEC || segmento->chargeminus->isLossy || segmento->chargeminus->isLossy || segmento->chargeminus->isLossy) {
                if ((!segmento->chargeminus->Is_LeftEnd) && (!segmento->chargeminus->Is_RightEnd)) {
                    segmento->chargeminus->Is_RightEnd = true;
                    sprintf(buff, "wir1_INFO: Forcing non-terminal node to RightEnd to attach RLC %7i%7i%7i",
                            segmento->chargeminus->i, segmento->chargeminus->j, segmento->chargeminus->k);
                    if ((segmento->chargeminus->k > ZI) && (segmento->chargeminus->k <= ZE) && control->verbose) WarnErrReport(buff);
                }
            } else {
                if (segmento->TipoWire->Series_C_RightEnd < 1.0e7_RKIND_wires) {
                    sprintf(buff, "wir1_ERROR: Series RightEnd Capacitance in INTERMEDIATE segment  smaller than 1e7 (inf) %7i%7i%7i%7i",
                            segmento->origIndex, segmento->i, segmento->j, segmento->k);
                    if ((segmento->k > ZI) && (segmento->k <= ZE) && (segmento->tipofield != iEz)) WarnErrReport(buff, true);
                    if ((segmento->k >= ZI) && (segmento->k <= ZE) && (segmento->tipofield == iEz)) WarnErrReport(buff, true);
                }
                if (segmento->TipoWire->Parallel_C_RightEnd != 0.0_RKIND_wires) {
                    sprintf(buff, "wir1_ERROR: Parallel RightEnd Capacitance in INTERMEDIATE segment %7i%7i%7i%7i",
                            segmento->origIndex, segmento->i, segmento->j, segmento->k);
                    if ((segmento->k > ZI) && (segmento->k <= ZE) && (segmento->tipofield != iEz)) WarnErrReport(buff, true);
                    if ((segmento->k >= ZI) && (segmento->k <= ZE) && (segmento->tipofield == iEz)) WarnErrReport(buff, true);
                }
                if (segmento->TipoWire->Parallel_R_RightEnd != 0.0_RKIND_wires) {
                    sprintf(buff, "wir1_ERROR: Series RightEnd Resistance in INTERMEDIATE segment %7i%7i%7i%7i",
                            segmento->origIndex, segmento->i, segmento->j, segmento->k);
                    if ((segmento->k > ZI) && (segmento->k <= ZE) && (segmento->tipofield != iEz)) WarnErrReport(buff, true);
                    if ((segmento->k >= ZI) && (segmento->k <= ZE) && (segmento->tipofield == iEz)) WarnErrReport(buff, true);
                }
            }
        }
    }
}

detect_peclossyconformal_nodes();

// END of geometrical PREPROCESSING of junctions and adjacencies

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Adjust constants and sort of nodes

// find PMC segments
for (i1 = 1; i1 <= HWires->NumCurrentSegments; i1++) {
    segmento = HWires->CurrentSegment(i1);
    switch (segmento->tipofield) {
        case iEx:
            if ((segmento->i + 1 == sgg->Alloc(segmento->tipofield)->XE) && (sgg->Border->IsFrontPMC)) {
                segmento->IsPMC = true;
            }
            if ((segmento->i == sgg->Alloc(segmento->tipofield)->XI) && (sgg->Border->IsBackPMC)) {
                segmento->IsPMC = true;
            }
            break;
        case iEy:
            if ((segmento->j + 1 == sgg->Alloc(segmento->tipofield)->YE) && (sgg->Border->IsRightPMC)) {
                segmento->IsPMC = true;
            }
            if ((segmento->j == sgg->Alloc(segmento->tipofield)->YI) && (sgg->Border->IsLeftPMC)) {
                segmento->IsPMC = true;
            }
            break;
        case iEz:
            if ((segmento->k + 1 == sgg->Alloc(segmento->tipofield)->ZE) && (sgg->Border->IsUpPMC)) {
                segmento->IsPMC = true;
            }
            if ((segmento->k == sgg->Alloc(segmento->tipofield)->ZI) && (sgg->Border->IsDownPMC)) {
                segmento->IsPMC = true;
            }
            break;
    }
    if (segmento->isPMC) {
        sprintf(buff, "wir1_WARNING: PMC endings in WIREs are UNTESTED");
        if ((segmento->k >= ZI) && (segmento->k <= ZE)) WarnErrReport(buff);
    }
}

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// find nodes at the boundaries
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// ojo que lo del shielding siguiente no es valido si hay branching. habria que shieldear canda segmento del junction !\E7\E7
for (i1 = 1; i1 <= HWires->NumChargeNodes; i1++) {
    nodo = HWires->Chargenode(i1);
    if (nodo->i == SINPML_FULLSIZE(iHx)->XI) {
        if (sgg->Border->IsBackPML || sgg->Border->IsBackMur) {
            thereAreMurConditions = true;
            nodo->IsMur = true;
            nodo->IsBackDownLeftMur = true;
            nodo->IsFrontUpRightMur = false;
            dummy = nodo->CurrentPlus_1;
            // adjust constants
            HWires->Chargenode(i1)->NodeInside = HWires->Chargenode(i1)->CurrentPlus_1->ChargePlus;
            // !!dummy->IsShielded                                                  = true;
            // !!dummy->ChargePlus->CurrentPlus_1->IsShielded                         = true;
            // !!dummy->ChargePlus->CurrentPlus_1->ChargePlus->CurrentPlus_1->IsShielded = true;
            nodo->cteMur = (sqrt(InvMu(dummy->indexmed) * InvEps(dummy->indexmed)) * sgg->dt - dummy->delta) /
                           (sqrt(InvMu(dummy->indexmed) * InvEps(dummy->indexmed)) * sgg->dt + dummy->delta);

            //
            sprintf(buff, "wir1_WARNING: Node %7i%7i%7i with BACK Mur conditions", nodo->I, nodo->J, nodo->k);
            if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
        } else if (sgg->Border->IsBackPeriodic) {
            nodo->IsPeriodic = true;
            sprintf(buff, "wir1_WARNING: Node %7i%7i%7i with BACK Periodic conditions", nodo->I, nodo->J, nodo->k);
            if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
        } else if (sgg->Border->IsBackPEC) {
            nodo->IsPec = true;
            sprintf(buff, "wir1_WARNING: Node %7i%7i%7i with BACK PEC boundary conditions", nodo->I, nodo->J, nodo->k);
            if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
        }
    }
    if (nodo->i == SINPML_FULLSIZe(iHx)->XE) {
        if (sgg->Border->IsFrontPML || sgg->Border->IsFrontMur) {
            thereAreMurConditions = true;
            nodo->IsMur = true;
            nodo->IsBackDownLeftMur = false;
            nodo->IsFrontUpRightMur = true;
            dummy = nodo->CurrentMInus_1;
            // adjust constantas shielding 3 levels of segments (no junctions permitted) !must be tested!!!!
            HWires->Chargenode(i1)->NodeInside = HWires->Chargenode(i1)->CurrentMinus_1->ChargeMinus;
            // !!dummy->IsShielded                                                      = true;  !One current back
            // !!dummy->ChargeMinus->CurrentMinus_1->IsShielded                           = true;  !Two currents back
            // !!dummy->ChargeMinus->CurrentMinus_1->ChargeMinus->CurrentMinus_1->IsShielded = true;  !Three current back
            nodo->cteMur = (sqrt(InvMu(dummy->indexmed) * InvEps(dummy->indexmed)) * sgg->dt - dummy->delta) /
                           (sqrt(InvMu(dummy->indexmed) * InvEps(dummy->indexmed)) * sgg->dt + dummy->delta);
            //
            sprintf(buff, "wir1_WARNING: Node %7i%7i%7i with FRONT Mur conditions", nodo->I, nodo->J, nodo->k);
            if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
        } else if (sgg->Border->IsFrontPeriodic) {
            nodo->IsPeriodic = true;
            sprintf(buff, "wir1_WARNING: Node %7i%7i%7i with Front Periodic conditions", nodo->I, nodo->J, nodo->k);
            if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
        } else if (sgg->Border->IsFrontPEC) {
            nodo->IsPec = true;
            sprintf(buff, "wir1_WARNING: Node %7i%7i%7i with FRONT PEC boundary conditions", nodo->I, nodo->J, nodo->k);
            if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
        }
    }
    //
    if (nodo->j == SINPML_FULLSIZe(iHy)->YI) {
        if (sgg->Border->IsLeftPML || sgg->Border->IsLeftMur) {
            thereAreMurConditions = true;
            nodo->IsMur = true;
            nodo->IsBackDownLeftMur = true;
            nodo->IsFrontUpRightMur = false;
            dummy = nodo->CurrentPlus_1;
            // adjust constants
            HWires->Chargenode(i1)->NodeInside = HWires->Chargenode(i1)->CurrentPlus_1->ChargePlus;
            // !!dummy->IsShielded                                                  = true;
            // !!dummy->ChargePlus->CurrentPlus_1->IsShielded                         = true;
            // !!dummy->ChargePlus->CurrentPlus_1->ChargePlus->CurrentPlus_1->IsShielded = true;
            nodo->cteMur = (sqrt(InvMu(dummy->indexmed) * InvEps(dummy->indexmed)) * sgg->dt - dummy->delta) /
                           (sqrt(InvMu(dummy->indexmed) * InvEps(dummy->indexmed)) * sgg->dt + dummy->delta);

            //
            sprintf(buff, "wir1_WARNING: Node %7i%7i%7i with LEFT Mur conditions", nodo->I, nodo->J, nodo->k);
            if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
        } else if (sgg->Border->IsLeftPeriodic) {
            nodo->IsPeriodic = true;
            sprintf(buff, "wir1_WARNING: Node %7i%7i%7i with Left Periodic conditions", nodo->I, nodo->J, nodo->k);
            if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
        } else if (sgg->Border->IsLeftPEC) {
            nodo->IsPec = true;
            sprintf(buff, "wir1_WARNING: Node %7i%7i%7i with LEFT PEC boundary conditions", nodo->I, nodo->J, nodo->k);
            if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
        }
    }
    if (nodo->j == SINPML_FULLSIZe(iHy)->YE) {
        if (sgg->Border->IsRightPML || sgg->Border->IsRightMur) {
            thereAreMurConditions = true;
            nodo->IsMur = true;
            nodo->IsBackDownLeftMur = false;
            nodo->IsFrontUpRightMur = true;
            dummy = nodo->CurrentMInus_1;
            // adjust constantas shielding 3 levels of segments (no junctions permitted) !must be tested!!!!
            HWires->Chargenode(i1)->NodeInside = HWires->Chargenode(i1)->CurrentMinus_1->ChargeMinus;
            // !!dummy->IsShielded                                                      = true;  !One current back
            // !!dummy->ChargeMinus->CurrentMinus_1->IsShielded                           = true;  !Two currents back
            // !!dummy->ChargeMinus->CurrentMinus_1->ChargeMinus->CurrentMinus_1->IsShielded = true;  !Three current back
            nodo->cteMur = (sqrt(InvMu(dummy->indexmed) * InvEps(dummy->indexmed)) * sgg->dt - dummy->delta) /
                           (sqrt(InvMu(dummy->indexmed) * InvEps(dummy->indexmed)) * sgg->dt + dummy->delta);
            //
            sprintf(buff, "wir1_WARNING: Node %7i%7i%7i with RIGHT Mur conditions", nodo->I, nodo->J, nodo->k);
            if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
        } else if (sgg->Border->IsRightPeriodic) {
            nodo->IsPeriodic = true;
            sprintf(buff, "wir1_WARNING: Node %7i%7i%7i with Right Periodic conditions", nodo->I, nodo->J, nodo->k);
            if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
        } else if (sgg->Border->IsRightPEC) {
            nodo->IsPec = true;
            sprintf(buff, "wir1_WARNING: Node %7i%7i%7i with RIGHT PEC boundary conditions", nodo->I, nodo->J, nodo->k);
            if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
        }
    }
    //
    if (nodo->k == SINPML_FULLSIZe(iHz)->ZI) {
        if (sgg->Border->IsDownPML || sgg->Border->IsDownMur) {
            thereAreMurConditions = true;
            nodo->IsMur = true;
            nodo->IsBackDownLeftMur = true;
            nodo->IsFrontUpRightMur = false;
            dummy = nodo->CurrentPlus_1;
            // adjust constants
            HWires->Chargenode(i1)->NodeInside = HWires->Chargenode(i1)->CurrentPlus_1->ChargePlus;
            // !!dummy->IsShielded                                                  = true;
            // !!dummy->ChargePlus->CurrentPlus_1->IsShielded                         = true;
            // !!dummy->ChargePlus->CurrentPlus_1->ChargePlus->CurrentPlus_1->IsShielded = true;
            nodo->cteMur = (sqrt(InvMu(dummy->indexmed) * InvEps(dummy->indexmed)) * sgg->dt - dummy->delta) /
                           (sqrt(InvMu(dummy->indexmed) * InvEps(dummy->indexmed)) * sgg->dt + dummy->delta);

            //
            sprintf(buff, "wir1_WARNING: Node %7i%7i%7i with DOWN Mur conditions", nodo->I, nodo->J, nodo->k);
            if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
        } else if (sgg->Border->IsDownPeriodic) {
            nodo->IsPeriodic = true;
            sprintf(buff, "wir1_WARNING: Node %7i%7i%7i with Down Periodic conditions", nodo->I, nodo->J, nodo->k);
            if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
        } else if (sgg->Border->IsDownPEC) {
            nodo->IsPec = true;
            sprintf(buff, "wir1_WARNING: Node %7i%7i%7i with DOWN PEC boundary conditions", nodo->I, nodo->J, nodo->k);
            if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
        }
    }
    if (nodo->k == SINPML_FULLSIZe(iHz)->ZE) {
        if (sgg->Border->IsUpPML || sgg->Border->IsUpMur) {
            thereAreMurConditions = true;
            nodo->IsMur = true;
            nodo->IsBackDownLeftMur = false;
            nodo->IsFrontUpRightMur = true;
            dummy = nodo->CurrentMInus_1;
            // adjust constantas shielding 3 levels of segments (no junctions permitted) !must be tested!!!!
            HWires->Chargenode(i1)->NodeInside = HWires->Chargenode(i1)->CurrentMinus_1->ChargeMinus;
            // !!dummy->IsShielded                                                      = true;  !One current back
            // !!dummy->ChargeMinus->CurrentMinus_1->IsShielded                           = true;  !Two currents back
            // !!dummy->ChargeMinus->CurrentMinus_1->ChargeMinus->CurrentMinus_1->IsShielded = true;  !Three current back
            nodo->cteMur = (sqrt(InvMu(dummy->indexmed) * InvEps(dummy->indexmed)) * sgg->dt - dummy->delta) /
                           (sqrt(InvMu(dummy->indexmed) * InvEps(dummy->indexmed)) * sgg->dt + dummy->delta);
            //
            sprintf(buff, "wir1_WARNING: Node %7i%7i%7i with UP Mur conditions", nodo->I, nodo->J, nodo->k);
            if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
        } else if (sgg->Border->IsUpPeriodic) {
            nodo->IsPeriodic = true;
            sprintf(buff, "wir1_WARNING: Node %7i%7i%7i with Up Periodic conditions", nodo->I, nodo->J, nodo->k);
            if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
        } else if (sgg->Border->IsUpPEC) {
            nodo->IsPec = true;
            sprintf(buff, "wir1_WARNING: Node %7i%7i%7i with UP PEC boundary conditions", nodo->I, nodo->J, nodo->k);
            if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
        }
    }
    if ((nodo->isMur) && (nodo->NumCurrentPlus + nodo->NumCurrentMinus > 1)) {
        sprintf(buff, "wir1_WARNING: Node %7i%7i%7i is both a non-open and a Mur", nodo->I, nodo->J, nodo->k);
        if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
    }
}

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Constantes de actualizacion de nodos
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

for (i1 = 1; i1 <= HWires->NumChargeNodes; i1++) {
    nodo = HWires->ChargeNode(i1);
    // asign constants taking junctions into account
    deltadummy = 0.0_RKIND_wires;
    indexnode = 0;
    for (j1 = 1; j1 <= nodo->NumCurrentPlus + nodo->NumCurrentMinus; j1++) {
        if (j1 <= nodo->NumCurrentPlus) {
            switch (j1) {
                case 1: segmento = nodo->CurrentPlus_1; break;
                case 2: segmento = nodo->CurrentPlus_2; break;
                case 3: segmento = nodo->CurrentPlus_3; break;
                case 4: segmento = nodo->CurrentPlus_4; break;
                case 5: segmento = nodo->CurrentPlus_5; break;
                case 6: segmento = nodo->CurrentPlus_6; break;
                case 7: segmento = nodo->CurrentPlus_7; break;
                case 8: segmento = nodo->CurrentPlus_8; break;
                case 9: segmento = nodo->CurrentPlus_9; break;
            }
        } else {
            switch (j1 - nodo->NumCurrentPlus) {
                case 1: segmento = nodo->CurrentMinus_1; break;
                case 2: segmento = nodo->CurrentMinus_2; break;
                case 3: segmento = nodo->CurrentMinus_3; break;
                case 4: segmento = nodo->CurrentMinus_4; break;
                case 5: segmento = nodo->CurrentMinus_5; break;
                case 6: segmento = nodo->CurrentMinus_6; break;
                case 7: segmento = nodo->CurrentMinus_7; break;
                case 8: segmento = nodo->CurrentMinus_8; break;
                case 9: segmento = nodo->CurrentMinus_9; break;
            }
        }
        deltadummy = deltadummy + segmento->delta / 2.0_RKIND_wires; // holland pag. 92
        indexnode = indexnode + segmento->OrigIndex;
        //
    }
    //
    nodo->indexnode = indexnode; // just for normalization of numbers
}

// and debugging mpi/nonmpi
        nodo.CteProp = 1.0_RKIND_wires;
        nodo.CtePlain = sgg.dt / deltadummy;
        if (nodo.NumCurrentPlus + nodo.NumCurrentMinus == 1) {
            nodo.CtePlain = sgg.dt / (2.0_RKIND_wires * deltadummy); // correct the ending in case of open terminal
        }
        NODO.oRIGctePlain = nodo.CtePlain; // lo salvo para que al ungroundear retome esta constante
    }

    // !!!!!!!!!
    // corrige cteplain para PEC lossy nodes
    for (i1 = 1; i1 <= HWires.NumChargeNodes; ++i1) {
        nodo = HWires.ChargeNode(i1);
        if (nodo.IsLossy || nodo.isPEC) {
            nodo.CtePlain = 0.0_RKIND_wires;
            nodo.cteprop = 0.0_RKIND_wires;
        }
    }

    // !!!!FINAL reporting of unGROUNDING OF PEC/LOSSY NODES

    for (i1 = 1; i1 <= HWires.NumChargeNodes; ++i1) {
        nodo = Hwires.ChargeNode(i1);
        if ((nodo.k >= sgg.Sweep(iHz).ZI) && (nodo.k <= sgg.Sweep(iHz).ZE)) { // check that the node is inside the layout (just for MPI)
            if (nodo.IsPec || nodo.IsLossy) {
                if ((nodo.Is_LeftEnd) || (nodo.Is_RightEnd) || (nodo.NumCurrentPlus + nodo.NumCurrentMinus < 2)) { // es un nodo terminal
                    continue;
                } else { // no es un nodo terminal LeftEnd o RightEnd o abierto pero Si es un nodo material
                    if (!control.groundwires) {
                        if (nodo.IsPec) {
                            sprintf(buff, 'wir1_BUGGYERROR: NON-terminal node lying on PEC ()%i7%3i7%2i3%a%2i3%a', 
                                'wir1_BUGGYERROR: NON-terminal node lying on PEC ()', 
                                nodo.indexnode, nodo.i, nodo.j, nodo.k, 
                                ' (', nodo.numcurrentminus, nodo.numcurrentplus, ')');
                        }
                        if (nodo.IsLossy) {
                            sprintf(buff, 'wir1_BUGGYERROR: NON-terminal node lying on Lossy ()%i7%3i7%2i3%a%2i3%a', 
                                'wir1_BUGGYERROR: NON-terminal node lying on Lossy ()', 
                                nodo.indexnode, nodo.i, nodo.j, nodo.k, 
                                ' (', nodo.numcurrentminus, nodo.numcurrentplus, ')');
                        }
                        if ((nodo.k >= ZI) && (nodo.k <= ZE)) {
                            WarnErrReport(buff, true);
                        }
                        // if (nodo.isPEC)   nodo.isPEC = false;
                        // if (nodo.isLossy) nodo.isLossy = false;
                        // !!!restaura las constantes
                        // nodo.CteProp = 1.0_RKIND_wires;
                        // nodo.Cteplain = nodo.OrigCteplain;
                    } else {
                        // lo deja como esta
                        if (nodo.IsPec) {
                            sprintf(buff, 'wir1_INFO: Leaving grounded a  NON-Terminal node lying on PEC (LAUNCHED WITH -groundwires) %i7%3i7%2i3%a%2i3%a', 
                                'wir1_INFO: Leaving grounded a  NON-Terminal node lying on PEC (LAUNCHED WITH -groundwires) ', 
                                nodo.indexnode, nodo.i, nodo.j, nodo.k, 
                                ' (', nodo.numcurrentminus, nodo.numcurrentplus, ')');
                        }
                        if (nodo.IsLossy) {
                            sprintf(buff, 'wir1_INFO: Leaving grounded a  NON-Terminal node lying on Lossy (LAUNCHED WITH -groundwires) %i7%3i7%2i3%a%2i3%a', 
                                'wir1_INFO: Leaving grounded a  NON-Terminal node lying on Lossy (LAUNCHED WITH -groundwires) ', 
                                nodo.indexnode, nodo.i, nodo.j, nodo.k, 
                                ' (', nodo.numcurrentminus, nodo.numcurrentplus, ')');
                        }
                        if ((nodo.k >= ZI) && (nodo.k <= ZE) && control.verbose) {
                            WarnErrReport(buff);
                        }
                    }
                }
            }
        }
    }

    // gestion condiciones mur en nodos frontera

    for (i1 = 1; i1 <= HWires.NumChargeNodes; ++i1) {
        nodo = HWires.ChargeNode(i1);

        for (j1 = 1; j1 <= nodo.NumCurrentPlus + nodo.NumCurrentMinus; ++j1) {
            if (j1 <= nodo.NumCurrentPlus) {
                switch (j1) {
                    case 1: segmento = nodo.CurrentPlus_1; break;
                    case 2: segmento = nodo.CurrentPlus_2; break;
                    case 3: segmento = nodo.CurrentPlus_3; break;
                    case 4: segmento = nodo.CurrentPlus_4; break;
                    case 5: segmento = nodo.CurrentPlus_5; break;
                    case 6: segmento = nodo.CurrentPlus_6; break;
                    case 7: segmento = nodo.CurrentPlus_7; break;
                    case 8: segmento = nodo.CurrentPlus_8; break;
                    case 9: segmento = nodo.CurrentPlus_9; break;
                }
            } else {
                switch (j1 - nodo.NumCurrentPlus) {
                    case 1: segmento = nodo.CurrentMinus_1; break;
                    case 2: segmento = nodo.CurrentMinus_2; break;
                    case 3: segmento = nodo.CurrentMinus_3; break;
                    case 4: segmento = nodo.CurrentMinus_4; break;
                    case 5: segmento = nodo.CurrentMinus_5; break;
                    case 6: segmento = nodo.CurrentMinus_6; break;
                    case 7: segmento = nodo.CurrentMinus_7; break;
                    case 8: segmento = nodo.CurrentMinus_8; break;
                    case 9: segmento = nodo.CurrentMinus_9; break;
                }
            }

            // HOLD 251019
            if ((segmento.HasAbsorbing_RightEnd) && (nodo.Is_RightEnd)) {
                thereAreMurConditions = true;
                nodo.IsMur = true;
                if (nodo.CurrentPlus_1 is associated) { 
                    nodo.NodeInside = nodo.CurrentPlus_1.ChargePlus;
                    dummy = nodo.CurrentPlus_1;
                } else if (nodo.CurrentMinus_1 is associated) {  
                    dummy = nodo.CurrentMinus_1;
                    nodo.NodeInside = nodo.CurrentMinus_1.ChargeMinus;
                } else {
                    sprintf(buff, 'wir0_BUGGYERROR: Mur1 on wires wrong. ');
                    WarnErrReport(buff, true);
                }
                nodo.cteMur = (sqrt(InvMu(dummy.indexmed) * InvEps(dummy.indexmed)) * sgg.dt - dummy.delta) / 
                              (sqrt(InvMu(dummy.indexmed) * InvEps(dummy.indexmed)) * sgg.dt + dummy.delta);

            } else if ((segmento.HasAbsorbing_LeftEnd) && (nodo.Is_LeftEnd)) {
                thereAreMurConditions = true;
                nodo.IsMur = true;
                if (nodo.CurrentPlus_1 is associated) { 
                    nodo.NodeInside = nodo.CurrentPlus_1.ChargePlus;
                    dummy = nodo.CurrentPlus_1;
                } else if (nodo.CurrentMinus_1 is associated) {  
                    dummy = nodo.CurrentMinus_1;
                    nodo.NodeInside = nodo.CurrentMinus_1.ChargeMinus;
                } else {
                    sprintf(buff, 'wir0_BUGGYERROR: Mur1 on wires wrong. ');
                    WarnErrReport(buff, true);
                }
                nodo.cteMur = (sqrt(InvMu(dummy.indexmed) * InvEps(dummy.indexmed)) * sgg.dt - dummy.delta) / 
                              (sqrt(InvMu(dummy.indexmed) * InvEps(dummy.indexmed)) * sgg.dt + dummy.delta);
            }
            // !!!
        } // del barrido de nodos !movido aqui 07/02/14 bug esn.nfde solo en reporte
    }

    // I create the FranctionPlus and FractionMinus actualization constants needed to update the currents in junctions
    for (i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
        segmento = HWires.CurrentSegment(i1);
        nodo = segmento.ChargeMinus;
        segmento.FractionMinus = -1.0e30_RKIND_wires; // valor absurdo tiene que entrar siempre en algun case
        // Junctions and plain nodes require this correction in case of two distinct radius segments meet
        DenominatorFractionMinusDummy = 0.0_RKIND_wires;
        deltadummy1 = 0.0_RKIND_wires;
        if (nodo.NumCurrentPlus <= 0) {
            StopOnError(0, 0, 'Bug fractionminuss. ');
        }
        for (j1 = 1; j1 <= nodo.NumCurrentPlus; ++j1) {
            switch (j1) {
                case 1: dummy = nodo.CurrentPlus_1; break;
                case 2: dummy = nodo.CurrentPlus_2; break;
                case 3: dummy = nodo.CurrentPlus_3; break;
                case 4: dummy = nodo.CurrentPlus_4; break;
                case 5: dummy = nodo.CurrentPlus_5; break;
                case 6: dummy = nodo.CurrentPlus_6; break;
                case 7: dummy = nodo.CurrentPlus_7; break;
                case 8: dummy = nodo.CurrentPlus_8; break;
                case 9: dummy = nodo.CurrentPlus_9; break;
            }
            deltadummy1 = deltadummy1 + dummy.delta;
            DenominatorFractionMinusDummy = DenominatorFractionMinusDummy + dummy.delta / (dummy.Lind * InvMu(dummy.indexmed) * InvEPS(dummy.indexmed));
        }
        for (j1 = 1; j1 <= nodo.NumCurrentMinus; ++j1) {
            switch (j1) {
                case 1: dummy = nodo.CurrentMinus_1; break;
                case 2: dummy = nodo.CurrentMinus_2; break;
                case 3: dummy = nodo.CurrentMinus_3; break;
                case 4: dummy = nodo.CurrentMinus_4; break;
                case 5: dummy = nodo.CurrentMinus_5; break;
                case 6: dummy = nodo.CurrentMinus_6; break;
                case 7: dummy = nodo.CurrentMinus_7; break;
                case 8: dummy = nodo.CurrentMinus_8; break;
                case 9: dummy = nodo.CurrentMinus_9; break;
            }
            deltadummy1 = deltadummy1 + dummy.delta;
            DenominatorFractionMinusDummy = DenominatorFractionMinusDummy + dummy.delta / (dummy.Lind * InvMu(dummy.indexmed) * InvEPS(dummy.indexmed));
        }
        segmento.FractionMinus = (deltadummy1 / (segmento.Lind * InvMu(segmento.indexmed) * InvEPS(segmento.indexmed))) / DenominatorFractionMinusDummy; // Hollond paper on wires'81 (page 91) <--- esta mal. lo de berenger en -wiresflavor new es lo correcto.. yo lo hago ahora 22/07/15
        //
        nodo = segmento.ChargePlus;
        segmento.FractionPlus = -1.0e30_RKIND_wires; // valor absurdo tiene que entrar siempre en algun case
        DenominatorFractionPlusDummy = 0.0_RKIND_wires;
        deltadummy2 = 0.0_RKIND_wires;
        if (nodo.NumCurrentMinus <= 0) {
            StopOnError(0, 0, 'Bug fractionplus. ');
        }
        for (j1 = 1; j1 <= nodo.NumCurrentMinus; ++j1) {
            switch (j1) {
                case 1: dummy = nodo.CurrentMinus_1; break;
                case 2: dummy = nodo.CurrentMinus_2; break;
                case 3: dummy = nodo.CurrentMinus_3; break;
                case 4: dummy = nodo.CurrentMinus_4; break;
                case 5: dummy = nodo.CurrentMinus_5; break;
                case 6: dummy = nodo.CurrentMinus_6; break;
                case 7: dummy = nodo.CurrentMinus_7; break;
                case 8: dummy = nodo.CurrentMinus_8; break;
                case 9: dummy = nodo.CurrentMinus_9; break;
            }
            deltadummy2 = deltadummy2 + dummy.delta;
            DenominatorFractionPlusDummy = DenominatorFractionPlusDummy + dummy.delta / (dummy.Lind * InvMu(dummy.indexmed) * InvEPS(dummy.indexmed));
        }
        for (j1 = 1; j1 <= nodo.NumCurrentPlus; ++j1) {
            switch (j1) {
                case 1: dummy = nodo.CurrentPlus_1; break;
                case 2: dummy = nodo.CurrentPlus_2; break;
                case 3: dummy = nodo.CurrentPlus_3; break;
                case 4: dummy = nodo.CurrentPlus_4; break;
                case 5: dummy = nodo.CurrentPlus_5; break;
                case 6: dummy = nodo.CurrentPlus_6; break;
                case 7: dummy = nodo.CurrentPlus_7; break;
                case 8: dummy = nodo.CurrentPlus_8; break;
                case 9: dummy = nodo.CurrentPlus_9; break;
            }
            deltadummy2 = deltadummy2 + dummy.delta;
            DenominatorFractionPlusDummy = DenominatorFractionPlusDummy + dummy.delta / (dummy.Lind * InvMu(dummy.indexmed) * InvEPS(dummy.indexmed));
        }
        segmento.FractionPlus = (deltadummy2 / (segmento.Lind * InvMu(segmento.indexmed) * InvEPS(segmento.indexmed))) / DenominatorFractionPlusDummy;
        continue;
    }

    // End of the adjusting of constants and sort of nodes

    // detect inverse oriented segments bug OLD segmentos al reves 11/03/15

    for (conta = 1; conta <= HWires.NumCurrentSegments; ++conta) {
        // only for the observation sign to match (not used in this routine)
        if (HWires.CurrentSegment(conta).chargeplus.currentplus_1 is associated) {
            if ((HWires.CurrentSegment(conta).chargeplus.currentplus_1.origindex < HWires.CurrentSegment(conta).origindex) && 
                (HWires.CurrentSegment(conta).chargeplus.currentplus_1.indexmed == HWires.CurrentSegment(conta).indexmed)) { // tienen que estar en el mismo hilo
                HWires.CurrentSegment(conta).orientadoalreves = true; // relies on ORIGINAL orientation
            }
        }
        if (HWires.CurrentSegment(conta).chargeplus.currentplus_2 is associated) {
            if ((HWires.CurrentSegment(conta).chargeplus.currentplus_2.origindex < HWires.CurrentSegment(conta).origindex) && 
                (HWires.CurrentSegment(conta).chargeplus.currentplus_2.indexmed == HWires.CurrentSegment(conta).indexmed)) { // tienen que estar en el mismo hilo
                HWires.CurrentSegment(conta).orientadoalreves = true; // relies on ORIGINAL orientation
            }
        }
        if (HWires.CurrentSegment(conta).chargeplus.currentminus_1 is associated) { // puede que sea yo mismo, por eso miro tambien el currentminus2
            if ((HWires.CurrentSegment(conta).chargeplus.currentminus_1.origindex < HWires.CurrentSegment(conta).origindex) && 
                (HWires.CurrentSegment(conta).chargeplus.currentminus_1.indexmed == HWires.CurrentSegment(conta).indexmed)) { // tienen que estar en el mismo hilo
                HWires.CurrentSegment(conta).orientadoalreves = true; // relies on ORIGINAL orientation
            }
        }
        if (HWires.CurrentSegment(conta).chargeplus.currentminus_2 is associated) {
            if ((HWires.CurrentSegment(conta).chargeplus.currentminus_2.origindex < HWires.CurrentSegment(conta).origindex) && 
                (HWires.CurrentSegment(conta).chargeplus.currentminus_2.indexmed == HWires.CurrentSegment(conta).indexmed)) { // tienen que estar en el mismo hilo
                HWires.CurrentSegment(conta).orientadoalreves = true; // relies on ORIGINAL orientation
            }
        }
        // !!!!!!!!!
        if (HWires.CurrentSegment(conta).chargeminus.currentminus_1 is associated) {
            if ((HWires.CurrentSegment(conta).chargeminus.currentminus_1.origindex > HWires.CurrentSegment(conta).origindex) && 
                (HWires.CurrentSegment(conta).chargeminus.currentminus_1.indexmed == HWires.CurrentSegment(conta).indexmed)) {
                HWires.CurrentSegment(conta).orientadoalreves = true; // relies on ORIGINAL orientation
            }
        }
        if (HWires.CurrentSegment(conta).chargeminus.currentminus_2 is associated) {
            if ((HWires.CurrentSegment(conta).chargeminus.currentminus_2.origindex > HWires.CurrentSegment(conta).origindex) && 
                (HWires.CurrentSegment(conta).chargeminus.currentminus_2.indexmed == HWires.CurrentSegment(conta).indexmed)) {
                HWires.CurrentSegment(conta).orientadoalreves = true; // relies on ORIGINAL orientation
            }
        }
        if (HWires.CurrentSegment(conta).chargeminus.currentplus_1 is associated) { // puede que sea yo mismo, por eso miro tambien el currentplus2
            if ((HWires.CurrentSegment(conta).chargeminus.currentplus_1.origindex > HWires.CurrentSegment(conta).origindex) && 
                (HWires.CurrentSegment(conta).chargeminus.currentplus_1.indexmed == HWires.CurrentSegment(conta).indexmed)) {
                HWires.CurrentSegment(conta).orientadoalreves = true; // relies on ORIGINAL orientation
            }
        }
        if (HWires.CurrentSegment(conta).chargeminus.currentplus_2 is associated) {
            if ((HWires.CurrentSegment(conta).chargeminus.currentplus_2.origindex > HWires.CurrentSegment(conta).origindex) && 
                (HWires.CurrentSegment(conta).chargeminus.currentplus_2.indexmed == HWires.CurrentSegment(conta).indexmed)) {
                HWires.CurrentSegment(conta).orientadoalreves = true; // relies on ORIGINAL orientation
            }
        }
    }

    // Current sources (only permitted in Plain nodes -no hanging-)
    // Read the file with the time evolution of the source
    // I find ambiguity in the geometrical position (-untested with OLD-) which is given in .nfde
    // at at segment index (i,j,k). So  I assume that the current source
    // is located at the node wich is the minus direction of the segment. That is I detect if the the CurrentPlus_1
    // is there some current source

    for (i1 = 1; i1 <= HWires.NumChargeNodes; ++i1) {
        if (HWires.ChargeNode(i1).CurrentPlus_1 is associated) {
            dummy = HWires.ChargeNode(i1).CurrentPlus_1;
            for (k1 = 1; k1 <= dummy.TipoWire.numcurrentsources; ++k1) {
                if (dummy.TipoWire.IsourceExists) {
                    if ((dummy.i == dummy.TipoWire.Isource(k1).I) && 
                        (dummy.j == dummy.TipoWire.Isource(k1).J) && 
                        (dummy.k == dummy.TipoWire.Isource(k1).K)) {
                        if (HWires.ChargeNode(i1).NumCurrentPlus + HWires.ChargeNode(i1).NumCurrentMinus > 2) {
                            sprintf(buff, 'wir1_ERROR: Current source at %i7%a', 
                                'wir1_ERROR: Current source at ', HWIREs.ChargeNode(i1).indexnode, ' in junctions forbidden');
                            WarnErrReport(buff, true);
                        }
                        HWires.ChargeNode(i1).Isource = dummy.TipoWire.Isource(k1);
                        HWires.ChargeNode(i1).HasIsource = true;

                        thereareIsources = true;
                        sprintf(buff, 'wir1_INFO: Current source at node %i7%3i7',  
                            'wir1_INFO: Current source at node ', 
                            nodo.indexnode, nodo.i, nodo.j, nodo.k);
                        if ((HWires.ChargeNode(i1).k > ZI) && (HWires.ChargeNode(i1).k <= ZE) && control.verbose) {
                            WarnErrReport(buff);
                        }
                        //
                        if (HWires.ChargeNode(i1).Isource.Fichero.DeltaSamples > sgg.dt) {
                            sprintf(buff, 'wir1_WARNING: %s%e15.4e3', 
                                'wir1_WARNING: ', trim(adjustl(HWires.ChargeNode(i1).Isource.Fichero.Name)), 
                                ' undersampled by a factor ', HWires.ChargeNode(i1).Isource.Fichero.DeltaSamples / sgg.dt);
                            if ((HWires.ChargeNode(i1).k > ZI) && (HWires.ChargeNode(i1).k <= ZE)) {
                                WarnErrReport(buff);
                            }
                        }
                    }
                }
            }
        }
    }

    // End of geometrical PREPROCESSING
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    resume_casuistics();
    
    // !!!!!!!!!find crank-nicolson coefficients
    
    if (control.wirecrank) {
        init_wirecrank();
    }

    // Se arregla lo del resuming permit scaling 221118 pero realmente no se por que. deberia ser lo mismo. Es redundante, pero....
    eps000 = eps0; mu000 = mu0; // niapa para evitar lo del rkind 020719
    calc_wirehollandconstants(sgg, G2, control.fieldtotl, control.wiresflavor, mu000, eps000, control.simu_devia);
    
#ifdef CompileWithThickWires
    // init extra thick stuff !hay que hacerlo aqui al final con toda la info
    // ojo si permittivity scaling y si mpi.... habra que rehacer
    for (conta = 1; conta <= HWires.NumCurrentSegments; ++conta) {
        if (control.wirethickness > 1) {
            init_thick(sgg, eps000, mu000, Ex, Ey, Ez, Hx, Hy, Hz, Idxe, Idye, Idze, Idxh, Idyh, Idzh, HWires.CurrentSegment(conta), control.wirethickness);
        }
    }
#endif      
    return;

    contains

    // !!!!!!!!!
    
    void deembed_peclossyconformal_segments(int sggMiE) {
        
        // primero los conformal 130220 %Is%split_and_useless
        if ((sgg.Med(sggmiE).Is.split_and_useless) && 
             !(IsEnd_norLeft_norRight || Is_LeftEnd || Is_RightEnd)) {   // NO NO NO ES UN TERMINAL
            d();
        }
    }

eembed_segment;
        std::ostringstream buff_stream;
        buff_stream << "wir0_WARNING: YES de-embedding a NON-TERMINAL conformal split_and_useless WIRE segment: " 
                    << sggmiE << " " 
                    << HWires.CurrentSegment(conta).origIndex << " " 
                    << HWires.CurrentSegment(conta).i << " " 
                    << HWires.CurrentSegment(conta).j << " " 
                    << HWires.CurrentSegment(conta).k << " " 
                    << HWires.CurrentSegment(conta).tipofield;
        std::string buff = buff_stream.str();
        if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff);

    } else if ((sgg.Med(sggmiE).Is.split_and_useless) && 
               (IsEnd_norLeft_norRight || Is_LeftEnd || Is_RightEnd)) {   // SI SI SI ES UN TERMINAL
        deembed_segment();
        std::ostringstream buff_stream2;
        buff_stream2 << "wir0_ERROR: YES-TERMINAL WIRE SEGMENT IN A CONFORMAL split_and_useless SURFACE (): " 
                     << sggmiE << " " 
                     << HWires.CurrentSegment(conta).origIndex << " " 
                     << HWires.CurrentSegment(conta).i << " " 
                     << HWires.CurrentSegment(conta).j << " " 
                     << HWires.CurrentSegment(conta).k << " " 
                     << HWires.CurrentSegment(conta).tipofield;
        std::string buff2 = buff_stream2.str();
        if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff2, true);
    } else if ((sgg.Med(sggmiE).Is.already_YEEadvanced_byconformal) &&  // !!!!!!!!!!!!already_YEEadvanced_byconformal
               !(IsEnd_norLeft_norRight || Is_LeftEnd || Is_RightEnd)) {   // NO NO NO ES UN TERMINAL                    
        if (!control.fieldtotl) {
            HWires.CurrentSegment(conta).cte5 = sgg.dt / eps0 / (HWires.CurrentSegment(conta).deltaTransv1 * HWires.CurrentSegment(conta).deltaTransv2);
            // HWires.CurrentSegment(conta).cte5 = HWires.CurrentSegment(conta).cte5 / 3.5;
        } else {
            HWires.CurrentSegment(conta).cte5 = 0.0_RKIND_wires;
        }
        std::ostringstream buff_stream3;
        buff_stream3 << "wir0_WARNING: NO de-embedding a NON-TERMINAL conformal already_YEEadvanced_byconformal  WIRE segment: " 
                     << sggmiE << " " 
                     << HWires.CurrentSegment(conta).origIndex << " " 
                     << HWires.CurrentSegment(conta).i << " " 
                     << HWires.CurrentSegment(conta).j << " " 
                     << HWires.CurrentSegment(conta).k << " " 
                     << HWires.CurrentSegment(conta).tipofield;
        std::string buff3 = buff_stream3.str();
        if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff3);
    } else if ((sgg.Med(sggmiE).Is.already_YEEadvanced_byconformal) &&
               (IsEnd_norLeft_norRight || Is_LeftEnd || Is_RightEnd)) {  // SI SI SI ES UN TERMINAL                    
        if (!control.fieldtotl) {
            HWires.CurrentSegment(conta).cte5 = sgg.dt / eps0 / (HWires.CurrentSegment(conta).deltaTransv1 * HWires.CurrentSegment(conta).deltaTransv2);
            // HWires.CurrentSegment(conta).cte5 = HWires.CurrentSegment(conta).cte5 / 3.5;
        } else {
            HWires.CurrentSegment(conta).cte5 = 0.0_RKIND_wires;
        }
        std::ostringstream buff_stream4;
        buff_stream4 << "wir0_WARNING: NO de-embedding a YES-TERMINAL conformal already_YEEadvanced_byconformal  WIRE segment: " 
                     << sggmiE << " " 
                     << HWires.CurrentSegment(conta).origIndex << " " 
                     << HWires.CurrentSegment(conta).i << " " 
                     << HWires.CurrentSegment(conta).j << " " 
                     << HWires.CurrentSegment(conta).k << " " 
                     << HWires.CurrentSegment(conta).tipofield;
        std::string buff4 = buff_stream4.str();
        if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff4);
    } else if ((sggmiE == 0) || (sgg.med(sggmiE).is.pec) || 
               (std::abs(sgg.Med(sggmiE).sigma) >= 1.0e-15_RKIND_wires) || 
               (std::abs(sgg.Med(sggmiE).sigmam) >= 1.0e-15_RKIND_wires) || 
               sgg.Med(sggmiE).Is.Lossy) {
        deembed_segment();
        //
        if (!(IsEnd_norLeft_norRight || Is_LeftEnd || Is_RightEnd)) { // NO NO NO ES UN TERMINAL
            if ((sggmiE == 0) || (sgg.med(sggmiE).is.pec)) {
                std::ostringstream buff_stream5;
                buff_stream5 << "wir0_WARNING: YES De-embedding a NON-TERMINAL struct segment from PEC " 
                             << sggmiE << " " 
                             << HWires.CurrentSegment(conta).origIndex << " " 
                             << HWires.CurrentSegment(conta).i << " " 
                             << HWires.CurrentSegment(conta).j << " " 
                             << HWires.CurrentSegment(conta).k << " " 
                             << HWires.CurrentSegment(conta).tipofield;
                std::string buff5 = buff_stream5.str();
                if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff5);
            } else {
                std::ostringstream buff_stream6;
                buff_stream6 << "wir0_WARNING: YES De-embedding a NON-TERMINAL struct segment from Lossy " 
                             << sggmiE << " " 
                             << HWires.CurrentSegment(conta).origIndex << " " 
                             << HWires.CurrentSegment(conta).i << " " 
                             << HWires.CurrentSegment(conta).j << " " 
                             << HWires.CurrentSegment(conta).k << " " 
                             << HWires.CurrentSegment(conta).tipofield;
                std::string buff6 = buff_stream6.str();
                if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff6);
            }
        } else { // SI SI SI ES UN TERMINAL
            if ((sggmiE == 0) || (sgg.med(sggmiE).is.pec)) {
                deembed_segment(); // ojoo esto no estaba a 290323 y lo he aniadido porque parece que tiene sentido deembed si lo dice
                std::ostringstream buff_stream7;
                buff_stream7 << "wir0_SEVEREWARNING: YES de-embedding a YES-TERMINAL struct SEGMENT from PEC  " 
                             << sggmiE << " " 
                             << HWires.CurrentSegment(conta).origIndex << " " 
                             << HWires.CurrentSegment(conta).i << " " 
                             << HWires.CurrentSegment(conta).j << " " 
                             << HWires.CurrentSegment(conta).k << " " 
                             << HWires.CurrentSegment(conta).tipofield;
                std::string buff7 = buff_stream7.str();
                if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff7);
            } else {
                std::ostringstream buff_stream8;
                buff_stream8 << "wir0_SEVEREWARNING: YES de-embedding a YES-TERMINAL struct SEGMENT from Lossy  " 
                             << sggmiE << " " 
                             << HWires.CurrentSegment(conta).origIndex << " " 
                             << HWires.CurrentSegment(conta).i << " " 
                             << HWires.CurrentSegment(conta).j << " " 
                             << HWires.CurrentSegment(conta).k << " " 
                             << HWires.CurrentSegment(conta).tipofield;
                std::string buff8 = buff_stream8.str();
                if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff8);
            }
        }
        if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff); // Note: 'buff' here refers to the last assigned buffer in the previous block, or potentially undefined if logic flow differs, but preserving original variable usage.
    }
    // luego normales
    else {
        if ((sggmiE != 1) && (!sgg.Med(sggmiE).Is.ThinWire)) {
            std::ostringstream buff_stream9;
            buff_stream9 << "wir0_WARNING: NO de-embedding a terminal/non-terminal segment in Lossless medium " 
                         << sggmiE << " " 
                         << HWires.CurrentSegment(conta).origIndex << " " 
                         << HWires.CurrentSegment(conta).i << " " 
                         << HWires.CurrentSegment(conta).j << " " 
                         << HWires.CurrentSegment(conta).k << " " 
                         << HWires.CurrentSegment(conta).tipofield;
            std::string buff9 = buff_stream9.str();
            WarnErrReport(buff9);
        }
    }

    } // end of deembed_peclossyconformal_segments

    // !!!

    void deembed_segment() {
        // AJUSTAR LA CONSTANTE !no se precisa realmente porque luego se iran su campos a cero (ver nota antes) 180220
                
        if (!control.fieldtotl) {
            HWires.CurrentSegment(conta).cte5 = sgg.dt / eps0 / (HWires.CurrentSegment(conta).deltaTransv1 * HWires.CurrentSegment(conta).deltaTransv2);
            // HWires.CurrentSegment(conta).cte5 = HWires.CurrentSegment(conta).cte5 / 3.5;
        } else {
            HWires.CurrentSegment(conta).cte5 = 0.0_RKIND_wires;
        }
                    
        HWires.CurrentSegment(conta).Efield_wire2main = HWires.null_field; // YES de-embedding
        HWires.CurrentSegment(conta).Efield_main2wire = HWires.null_field; // YES de-embedding
    }

    void detect_peclossyconformal_nodes() {
        bool embed;

        // detect PEC and lossy nodes
        for (int i1 = 1; i1 <= HWires.NumChargeNodes; ++i1) {
            nodo = HWires.ChargeNode(i1);
            i = nodo->i;
            j = nodo->j;
            k = nodo->k;

            if ((i > SINPML_fullsize(iEx).XI) && 
                (i <= SINPML_fullsize(iEx).XE) && 
                (j > SINPML_fullsize(iEy).YI) && 
                (j <= SINPML_fullsize(iEy).YE) && 
                (k > SINPML_fullsize(iEz).ZI) && 
                (k <= SINPML_fullsize(iEz).ZE)) {

                if ((k <= sgg.alloc(iEZ).ZE) && (k >= sgg.alloc(iEZ).ZI)) {
                    kmenos1 = k - 1;
                    if (k - 1 < sgg.alloc(iEz).ZI) kmenos1 = k;
                    
                    nodo->IsLossy = 
                        sgg.Med(sggMino(i, j, k)).Is.Lossy || 
                        sgg.Med(sggMiEx(i, j, k)).Is.Lossy || 
                        sgg.Med(sggMiEy(i, j, k)).Is.Lossy || 
                        sgg.Med(sggMiEz(i, j, k)).Is.Lossy || 
                        sgg.Med(sggMiEx(i - 1, j, k)).Is.Lossy || 
                        sgg.Med(sggMiEy(i, j - 1, k)).Is.Lossy || 
                        sgg.Med(sggMiEz(i, j, kmenos1)).Is.Lossy;

                    sigt = MAX(MAX(MAX(sgg.Med(sggMiEx(i, j, k)).sigma, 
                                       sgg.Med(sggMiEy(i, j, k)).sigma),
                                       MAX(sgg.Med(sggMiEz(i, j, k)).sigma,  
                                       sgg.Med(sggMiEx(i - 1, j, k)).sigma)),
                                       MAX(sgg.Med(sggMiEy(i, j - 1, k)).sigma,  
                                       sgg.Med(sggMiEz(i, j, kmenos1)).sigma));
                    nodo->islossy = nodo->islossy || (std::abs(sigt) > 1.0e-19_RKIND_wires);

                    // 091024 para poner bien uniones hilo-lossy sgbc conformal !ojooo es bruto !bug FU50_50mm_conf    
                    // fallaria en un hipotetico conformal dielectrico con hilos unidos... habria que distinguir bien ojooooo y propagar este cambio al resto de sabores de hilos
                    // !!!!! if( sgg.Med(sggMino(i, j, k)).Is.already_YEEadvanced_byconformal || ... ) ...

                    nodo->ispec = 
                        (sggMiNo(i, j, k) == 0) || (sgg.med(sggMiNo(i, j, k)).is.pec) || 
                        (sggMiEx(i, j, k) == 0) || (sgg.med(sggMiEx(i, j, k)).is.pec) || 
                        (sggMiEy(i, j, k) == 0) || (sgg.med(sggMiEy(i, j, k)).is.pec) || 
                        (sggMiEz(i, j, k) == 0) || (sgg.med(sggMiEz(i, j, k)).is.pec) || 
                        (sggMiEx(i - 1, j, k) == 0) || (sgg.med(sggMiEx(i - 1, j, k)).is.pec) || 
                        (sggMiEy(i, j - 1, k) == 0) || (sgg.med(sggMiEy(i, j - 1, k)).is.pec) || 
                        (sggMiEz(i, j, kmenos1) == 0) || (sgg.med(sggMiEz(i, j, kmenos1)).is.pec);

                    if ((!((nodo->Is_RightEnd || nodo->Is_LeftEnd))) && (nodo->ispec || nodo->islossy)) {
                        if (nodo->ispec) {
                            std::ostringstream buff_stream10;
                            buff_stream10 << "wir1_WARNING: Un-grounding NON-TERMINAL node from PEC " << i << " " << j << " " << k;
                            std::string buff10 = buff_stream10.str();
                            if ((k > ZI) && (k <= ZE)) WarnErrReport(buff10);
                        }
                        if (nodo->islossy) {
                            std::ostringstream buff_stream11;
                            buff_stream11 << "wir1_WARNING: Un-grounding NON-TERMINAL node from Lossy " << i << " " << j << " " << k;
                            std::string buff11 = buff_stream11.str();
                            if ((k > ZI) && (k <= ZE)) WarnErrReport(buff11);
                        }
                        nodo->ispec = false;
                        nodo->islossy = false;
                    }

                    // pedazo de niapa para poner los nodos conformal a voltage nulo y sus segmentos conformal tambien si son already_YEEadvanced_byconformal (old notouch o no_touch) y que Dios reparta suerte 140220
                    if ((nodo->Is_RightEnd || nodo->Is_LeftEnd)) { // !!!los busy_nodes se han puesto a pec cuando en realidad estan unidos a already_YEEadvanced_byconformal .and.(.not.(nodo%ispec.or.nodo%islossy))) then !solo si es un extremo y no estaba ya puesto a pec            
                        medio1 = sggMiEx(i, j, k);
                        medio1m = sggMiEx(i - 1, j, k);
                        medio2 = sggMiEy(i, j, k);
                        medio2m = sggMiEy(i, j - 1, k);
                        medio3 = sggMiEz(i, j, k);
                        medio3m = sggMiEz(i, j, kmenos1);

                        if (sgg.med(medio1).is.split_and_useless || sgg.med(medio2).is.split_and_useless || sgg.med(medio3).is.split_and_useless || 
                            sgg.med(medio1m).is.split_and_useless || sgg.med(medio2m).is.split_and_useless || sgg.med(medio3m).is.split_and_useless) {
                            std::ostringstream buff_stream12;
                            buff_stream12 << "wir1_BUGGYERROR: Conformal node CONNECTED TO at least one split_and_useless conformal edge that cannot be safely set to PEC () " << i << " " << j << " " << k;
                            std::string buff12 = buff_stream12.str();
                            if ((k > ZI) && (k <= ZE)) WarnErrReport(buff12, true);
                        }

                        if (sgg.med(medio1).is.already_YEEadvanced_byconformal || sgg.med(medio2).is.already_YEEadvanced_byconformal || sgg.med(medio3).is.already_YEEadvanced_byconformal || 
                            sgg.med(medio1m).is.already_YEEadvanced_byconformal || sgg.med(medio2m).is.already_YEEadvanced_byconformal || sgg.med(medio3m).is.already_YEEadvanced_byconformal) {
                            nodo->ispec = true; // luego se pondra nodo%cteplain = 0 para todos los pec
                            std::ostringstream buff_stream13;
                            buff_stream13 << "wir1_WARNING: Conformal already_YEEadvanced_byconformal node changed to PEC grounded node at " << i << " " << j << " " << k;
                            std::string buff13 = buff_stream13.str();
                            if ((k > ZI) && (k <= ZE)) WarnErrReport(buff13);
                            
                            // ademas anotar y poner a cero el efield correspondiente already_YEEadvanced_byconformal overrideando el conformal_advance_E() que ya se ha hecho antes a partir de esta version   
                            if (sgg.med(medio1).is.already_YEEadvanced_byconformal) {
                                check_embed(embed, iEx, i, j, k);
                                if (!embed) {
                                    // sggmiEx(i,j,k)=0;!ojoo quitar luego solo para visualiz
                                    nodo->already_YEEadvanced_byconformal_changedtoPECfield1 = Ex(i, j, k);
                                    std::ostringstream buff_stream14;
                                    buff_stream14 << "wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, also changed to PEC line at Ex " << i << " " << j << " " << k;
                                    std::string buff14 = buff_stream14.str();
                                    if ((k > ZI) && (k <= ZE)) WarnErrReport(buff14);
                                } else {
                                    std::ostringstream buff_stream15;
                                    buff_stream15 << "wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, cannot be also changed to PEC line for embedding a wire at Ex " << i << " " << j << " " << k;
                                    std::string buff15 = buff_stream15.str();
                                    if ((k > ZI) && (k <= ZE)) WarnErrReport(buff15);
                                }
                            }
                            if (sgg.med(medio1m).is.already_YEEadvanced_byconformal) {
                                check_embed(embed, iEx, i - 1, j, k);
                                if (!embed) {
                                    // sggmiEx(i-1,j,k)=0;!ojoo quitar luego solo para visualiz
                                    nodo->already_YEEadvanced_byconformal_changedtoPECfield2 = Ex(i - 1, j, k);
                                    std::ostringstream buff_stream16;
                                    buff_stream16 << "wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, also changed to PEC line at mEx " << i - 1 << " " << j << " " << k;
                                    std::string buff16 = buff_stream16.str();
                                    if ((k > ZI) && (k <= ZE)) WarnErrReport(buff16);
                                } else {
                                    std::ostringstream buff_stream17;
                                    buff_stream17 << "wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, cannot be also changed to PEC line for embedding a wire at mEx " << i - 1 << " " << j << " " << k;
                                    std::string buff17 = buff_stream17.str();
                                    if ((k > ZI) && (k <= ZE)) WarnErrReport(buff17);
                                }
                            }
                            if (sgg.med(medio2).is.already_YEEadvanced_byconformal) {
                                check_embed(embed, iEy, i, j, k);
                                if (!embed) {
                                    // sggmiEy(i,j,k)=0;!ojoo quitar luego solo para visualiz
                                    nodo->already_YEEadvanced_byconformal_changedtoPECfield3 = Ey(i, j, k);
                                    std::ostringstream buff_stream18;
                                    buff_stream18 << "wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, also changed to PEC line at Ey " << i << " " << j << " " << k;
                                    std::string buff18 = buff_stream18.str();
                                    if ((k > ZI) && (k <= ZE)) WarnErrReport(buff18);
                                } else {
                                    std::ostringstream buff_stream19;
                                    buff_stream19 << "wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, cannot be also changed to PEC line for embedding a wire at Ey " << i << " " << j << " " << k;
                                    std::string buff19 = buff_stream19.str();
                                    if ((k > ZI) && (k <= ZE)) WarnErrReport(buff19);
                                }
                            }
                            if (sgg.med(medio2m).is.already_YEEadvanced_byconformal) {
                                check_embed(embed, iEy, i, j - 1, k);
                                if (!embed) {
                                    // sggmiEy(i,j-1,k)=0;!ojoo quitar luego solo para visualiz
                                    nodo->already_YEEadvanced_byconformal_changedtoPECfield4 = Ey(i, j - 1, k);
                                    std::ostringstream buff_stream20;
                                    buff_stream20 << "wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, also changed to PEC line at mEy " << i << " " << j - 1 << " " << k;
                                    std::string buff20 = buff_stream20.str();
                                    if ((k > ZI) && (k <= ZE)) WarnErrReport(buff20);
                                } else {
                                    std::ostringstream buff_stream21;
                                    buff_stream21 << "wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, cannot be also changed to PEC line for embedding a wire at mEy " << i << " " << j - 1 << " " << k;
                                    std::string buff21 = buff_stream21.str();
                                    if ((k > ZI) && (k <= ZE)) WarnErrReport(buff21);
                                }
                            }
                            if (sgg.med(medio3).is.already_YEEadvanced_byconformal) {
                                check_embed(embed, iEz, i, j, k);
                                if (!embed) {
                                    // sggmiEz(i,j,k)=0;!ojoo quitar luego solo para visualiz
                                    nodo->already_YEEadvanced_byconformal_changedtoPECfield5 = Ez(i, j, k);
                                    std::ostringstream buff_stream22;
                                    buff_stream22 << "wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, also changed to PEC line at Ez " << i << " " << j << " " << k;
                                    std::string buff22 = buff_stream22.str();
                                    if ((k > ZI) && (k <= ZE)) WarnErrReport(buff22);
                                } else {
                                    std::ostringstream buff_stream23;
                                    buff_stream23 << "wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, cannot be also changed to PEC line for embedding a wire at Ez " << i << " " << j << " " << k;
                                    std::string buff23 = buff_stream23.str();
                                    if ((k > ZI) && (k <= ZE)) WarnErrReport(buff23);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

}
        }
        if (sgg.med(medio3m).is.already_YEEadvanced_byconformal) {
            check_embed(embed, iEz, i, j, kmenos1);
            if (!embed) {
                //         sggmiEz(i,j,kmenos1)=0;!ojoo quitar luego solo para visualiz
                nodo.already_YEEadvanced_byconformal_changedtoPECfield6 = &Ez(i, j, kmenos1);
                std::cout << "wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, also changed to PEC line at mEz " << i << " " << j << " " << kmenos1 << std::endl;
                if ((k > ZI) && (k <= ZE)) WarnErrReport(buff);
            } else {
                std::cout << "wir1_WARNING: Conformal already_YEEadvanced_byconformal segment atached to PEC grounded node, cannot be also changed to PEC line for embedding a wire at mEz " << i << " " << j << " " << kmenos1 << std::endl;
                if ((k > ZI) && (k <= ZE)) WarnErrReport(buff);
            }
        }

        else if (sgg.med(medio1).is.split_and_useless && sgg.med(medio2).is.split_and_useless && sgg.med(medio3).is.split_and_useless &&
                 sgg.med(medio1m).is.split_and_useless && sgg.med(medio2m).is.split_and_useless && sgg.med(medio3m).is.split_and_useless) {
            std::cout << "wir1_ERROR: Conformal split_and_useless node NOT changed (IMPOSIBLE) to PEC grounded node at " << i << " " << j << " " << k << std::endl;
            if ((k > ZI) && (k <= ZE)) WarnErrReport(buff, true);
            //ojo no tendria arreglo porque aunque se updatee a 0, este campo luego esta partido en dos en conformal y es split_and_useless
        }
        //esta a PEC/lossy o al punietero aire
        //!!!lo he sacado del if. solo para reporte y que coincida con lo que da estructurado  280323
        if (nodo.ispec) {
            std::cout << "wir1_INFO: (SHOULD BE REDUNDANT) Terminal Node grounded to PEC " << i << " " << j << " " << k << std::endl;
            if ((k > ZI) && (k <= ZE)) WarnErrReport(buff);
        } else if (nodo.islossy) {
            std::cout << "wir1_INFO: (SHOULD BE REDUNDANT) Terminal Node grounded to Lossy " << i << " " << j << " " << k << std::endl;
            if ((k > ZI) && (k <= ZE)) WarnErrReport(buff);
        } else if (nodo.IsHeterogeneousJunction) {
            std::cout << "wir1_INFO: (SHOULD BE REDUNDANT) Terminal Node is heteoreneous junction " << i << " " << j << " " << k << std::endl;
            if ((k > ZI) && (k <= ZE)) WarnErrReport(buff);
        } else {
            std::cout << "wir1_INFO: (SHOULD BE REDUNDANT) Terminal Node embedded in air " << i << " " << j << " " << k << std::endl;
            if ((k > ZI) && (k <= ZE)) WarnErrReport(buff);
        }
        if ((k > ZI) && (k <= ZE)) WarnErrReport(buff);
    }
    //!!!!fin pedazo de niapa

    //!!!!! ojo cambio AGRESIVO: sgg 08092016 a peticion OLD. si un nodo es heterogeneousjunction automaticamente no es ni pec ni lossy !esto es delicado !validado con gra_simple_conectado179_162_173.nfde
    if (nodo.isPEC && nodo.IsHeterogeneousJunction) {
        nodo.isPEC = false;
        std::cout << "wir1_ERROR: (Deprecated 170220-010324: no longer an error, though I stop. If sure relaunch with -ignoreerrors) ENL/ENR PEC grounded node detached from ground-material for being a wire junction" << i << " " << j << " " << k << " " << nodo.ispec << " " << nodo.islossy << std::endl;
        if ((k > ZI) && (k <= ZE)) WarnErrReport(buff, true);
    }
    if (nodo.isLossy && nodo.IsHeterogeneousJunction) {
        nodo.isLossy = false;
        std::cout << "wir1_ERROR: (Deprecated 170220-010324: no longer an error, though I stop. If sure relaunch with -ignoreerrors) ENL/ENR Lossy grounded node detached from ground-material for being a wire junction" << i << " " << j << " " << k << " " << nodo.ispec << " " << nodo.islossy << std::endl;
        if ((k > ZI) && (k <= ZE)) WarnErrReport(buff, true);
    }
    //!!!!!!!!! fin cambio sgg 080912016

    if (nodo.isPec) {
        nodo.isLossy = false; //tocado 030615 para evitar conflictos en nodos frontera entre Lossy y pec
    }
    //
}
} //del fullsize
}

} // detect_peclossyconformal_nodes

void check_embed(bool& embed, int tipofieldo, int io, int jo, int ko) {
    bool embed;
    int ib, jb, kb, tipofieldb;
    int io, jo, ko, tipofieldo;
    CurrentSegments_t* dummy = nullptr;
    //!!!!!!!!!!!!!!!!!!!!!!!! embed=.true.; return !!!!ojoooo sgg tocado a mano para ver bug conformal 220323

    embed = false;
    if (nodo.CurrentMinus_1 != nullptr) {
        dummy = nodo.CurrentMinus_1;
        auxem(embed, tipofieldo, io, jo, ko, dummy);
    }
    if (nodo.CurrentMinus_2 != nullptr) {
        dummy = nodo.CurrentMinus_2;
        auxem(embed, tipofieldo, io, jo, ko, dummy);
    }
    if (nodo.CurrentMinus_3 != nullptr) {
        dummy = nodo.CurrentMinus_3;
        auxem(embed, tipofieldo, io, jo, ko, dummy);
    }
    if (nodo.CurrentMinus_4 != nullptr) {
        dummy = nodo.CurrentMinus_4;
        auxem(embed, tipofieldo, io, jo, ko, dummy);
    }
    if (nodo.CurrentMinus_5 != nullptr) {
        dummy = nodo.CurrentMinus_5;
        auxem(embed, tipofieldo, io, jo, ko, dummy);
    }
    if (nodo.CurrentMinus_6 != nullptr) {
        dummy = nodo.CurrentMinus_6;
        auxem(embed, tipofieldo, io, jo, ko, dummy);
    }
    if (nodo.CurrentMinus_7 != nullptr) {
        dummy = nodo.CurrentMinus_7;
        auxem(embed, tipofieldo, io, jo, ko, dummy);
    }
    if (nodo.CurrentMinus_8 != nullptr) {
        dummy = nodo.CurrentMinus_8;
        auxem(embed, tipofieldo, io, jo, ko, dummy);
    }
    if (nodo.CurrentMinus_9 != nullptr) {
        dummy = nodo.CurrentMinus_9;
        auxem(embed, tipofieldo, io, jo, ko, dummy);
    }
    //
    if (nodo.CurrentPlus_1 != nullptr) {
        dummy = nodo.CurrentPlus_1;
        auxem(embed, tipofieldo, io, jo, ko, dummy);
    }
    if (nodo.CurrentPlus_2 != nullptr) {
        dummy = nodo.CurrentPlus_2;
        auxem(embed, tipofieldo, io, jo, ko, dummy);
    }
    if (nodo.CurrentPlus_3 != nullptr) {
        dummy = nodo.CurrentPlus_3;
        auxem(embed, tipofieldo, io, jo, ko, dummy);
    }
    if (nodo.CurrentPlus_4 != nullptr) {
        dummy = nodo.CurrentPlus_4;
        auxem(embed, tipofieldo, io, jo, ko, dummy);
    }
    if (nodo.CurrentPlus_5 != nullptr) {
        dummy = nodo.CurrentPlus_5;
        auxem(embed, tipofieldo, io, jo, ko, dummy);
    }
    if (nodo.CurrentPlus_6 != nullptr) {
        dummy = nodo.CurrentPlus_6;
        auxem(embed, tipofieldo, io, jo, ko, dummy);
    }
    if (nodo.CurrentPlus_7 != nullptr) {
        dummy = nodo.CurrentPlus_7;
        auxem(embed, tipofieldo, io, jo, ko, dummy);
    }
    if (nodo.CurrentPlus_8 != nullptr) {
        dummy = nodo.CurrentPlus_8;
        auxem(embed, tipofieldo, io, jo, ko, dummy);
    }
    if (nodo.CurrentPlus_9 != nullptr) {
        dummy = nodo.CurrentPlus_9;
        auxem(embed, tipofieldo, io, jo, ko, dummy);
    }

    return;
}

void auxem(bool& embed, int tipofieldo, int io, int jo, int ko, CurrentSegments_t* dummy) {
    bool embed;
    int ib, jb, kb, tipofieldb;
    int io, jo, ko, tipofieldo;
    CurrentSegments_t* dummy;
    ib = dummy->i;
    jb = dummy->j;
    kb = dummy->k;
    tipofieldb = dummy->tipofield;
    embed = embed || ((ib == io) && (jb == jo) && (kb == ko) && (tipofieldo == tipofieldb));
    return;
}

void resume_casuistics() {
    //In case of resuming a problem, re-starting currents and charges must be read instead of initialized
    //resuming
    if (!control.resume) {
        for (int i1 = 1; i1 <= HWires.NumChargeNodes; ++i1) {
            HWires.ChargeNode[i1].ChargePresent = 0.0_RKIND_wires;
            HWires.ChargeNode[i1].ChargePast = 0.0_RKIND_wires;
        }
        for (int i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
            HWires.CurrentSegment[i1].Current = 0.0_RKIND_wires;
            HWires.CurrentSegment[i1].qplus_qminus = 0.0_RKIND_wires;
            HWires.CurrentSegment[i1].current_for_devia = 0.0_RKIND_wires;
            HWires.CurrentSegment[i1].qplus_qminus_for_devia = 0.0_RKIND_wires;
            HWires.CurrentSegment[i1].Efield_main2wire_for_devia = 0.0_RKIND_wires;
        }
#ifdef CompileWithMPI
        //Initialize MPI counters (later written to disk)
        Hwires.NumNeededCurrentUpMPI = 0;
        Hwires.NumNeededCurrentDownMPI = 0;
#endif
    } else {
        for (int i1 = 1; i1 <= HWires.NumChargeNodes; ++i1) {
            READ(14, HWires.ChargeNode[i1].ChargePresent);
            READ(14, HWires.ChargeNode[i1].ChargePast);
        }
        for (int i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
            READ(14, HWires.CurrentSegment[i1].Current);
            READ(14, HWires.CurrentSegment[i1].qplus_qminus);
            READ(14, HWires.CurrentSegment[i1].current_for_devia);
            READ(14, HWires.CurrentSegment[i1].qplus_qminus_for_devia);
            READ(14, HWires.CurrentSegment[i1].Efield_main2wire_for_devia);
        }
#ifdef CompileWithMPI
        READ(14, Hwires.NumNeededCurrentUpMPI, Hwires.NumNeededCurrentDownMPI);
        Hwires.MPIUpNeededCurrentSegment = new CurrentSegments_t[Hwires.NumNeededCurrentUpMPI + 1];
        Hwires.MPIDownNeededCurrentSegment = new CurrentSegments_t[Hwires.NumNeededCurrentDownMPI + 1];
        for (int i1 = 1; i1 <= Hwires.NumNeededCurrentUpMPI; ++i1) {
            READ(14, HWires.MPIUpNeededCurrentSegment[i1].Current);
        }
        for (int i1 = 1; i1 <= Hwires.NumNeededCurrentDownMPI; ++i1) {
            READ(14, HWires.MPIDownNeededCurrentSegment[i1].Current);
        }
        //for new wires MPI march'12 sgg bug en multiwires mpi

#endif
    }
}

void init_wirecrank() {
    if (control.layoutnumber != 0) {
        buff = "Unsupported crank in MPI";
        WarnErrReport(buff, true);
    }
    for (int n = 1; n <= HWires.NumCurrentSegments; ++n) {
        Segmento = &HWires.CurrentSegment[n];
        jmed = Segmento->indexmed;
        Resist = Segmento->TipoWire.R;
        Autoin = Segmento->Lind; //!!!+ Segmento->TipoWire.L !he comentado esto a 110517 porque con lo del fieldtotl ya he metido en %lind tambien la tipowire%L
        Capaci = 1.0_RKIND_wires / (InvMu(jmed) * InvEps(jmed) * Autoin);
        deltax = Segmento->delta;
        !!!!
        segmento->rightCHminus = 1.0_RKIND_wires / (deltax * Capaci);
        segmento->rightCHplus = -segmento->rightCHminus;
        segmento->rightCUminus = (sgg.dt) / (4.0_RKIND_wires * deltax * deltax * Capaci);
        segmento->rightCUplus = segmento->rightCUminus;
        !!!!!!
        segmento->upperdiag = -sgg.dt / (4.0_RKIND_wires * deltax * deltax * Capaci); //comun a todos
        segmento->lowerdiag = segmento->upperdiag; //comun a todos
        !!!
        if ((Segmento->ChargeMinus.NumCurrentMinus == 1) &&
            (Segmento->ChargeMinus.NumCurrentPlus == 1) &&
            (Segmento->ChargePlus.NumCurrentMinus == 1) &&
            (Segmento->ChargePlus.NumCurrentPlus == 1)) { //es un intermedio
            segmento->diag = (sgg.dt / (deltax * deltax * RKIND_wires * Capaci) + (2_RKIND_wires * Autoin) / sgg.dt + resist) / 2.0_RKIND_wires;
            segmento->rightCU = (-sgg.dt / (2.0_RKIND_wires * deltax * deltax * RKIND_wires * Capaci) + Autoin / sgg.dt - resist / 2.0_RKIND_wires);
        } else {
            segmento->diag = ((3_RKIND_wires * sgg.dt) / (2.0_RKIND_wires * deltax * deltax * Capaci) + (2_RKIND_wires * Autoin) / sgg.dt + resist) / 2.0_RKIND_wires; //comun a todos
            segmento->rightCU = ((-3_RKIND_wires * sgg.dt) / (4.0_RKIND_wires * deltax * deltax * Capaci) + Autoin / sgg.dt - resist / 2.0_RKIND_wires);
            if ((Segmento->ChargeMinus.NumCurrentMinus + Segmento->ChargeMinus.NumCurrentPlus == 1) &&
                (Segmento->ChargePlus.NumCurrentMinus == 1) &&
                (Segmento->ChargePlus.NumCurrentPlus == 1)) {
                segmento->lowerdiag = -1e30; //lo voideo pq no debe tener efecto
            } else if ((Segmento->ChargeMinus.NumCurrentMinus == 1) &&
                       (Segmento->ChargeMinus.NumCurrentPlus == 1) &&
                       (Segmento->ChargePlus.NumCurrentMinus + Segmento->ChargePlus.NumCurrentPlus == 1)) {
                segmento->upperdiag = -1e30_RKIND_wires; //lo voideo pq no debe tener efecto
            } else {
                buff = "Unsupported crank";
                WarnErrReport(buff, true);
            }
        }
    }
}

} // InitWires
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//end of initialization subroutine
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//!!!!141118 Permittivity scaling

void calc_wirehollandconstants(SGGFDTDINFO_t& sgg, std::vector<double>& G2, bool fieldtotl, const std::string& wiresflavor, double mu00, double eps00, bool simu_devia) {
    bool simu_devia;
    SGGFDTDINFO_t& sgg;
    std::vector<double>& G2;
    std::string wiresflavor;
    bool fieldtotl;
    double eps00, mu00;
    double dl;
    std::vector<std::vector<double>>* Den = nullptr;
    int n, jmed, layoutnumber, iw1, is1, is2, i1, NumParallel, imed;

    double df1, df3, df2, Ddf1, Ddf3, Ddf2, vf1, vf3, vf2, runit;

    CurrentSegments_t* dummy = nullptr;
    ChargeNodes_t* Nodo = nullptr;
    TMultiline_t* Multiline = nullptr;

    eps0 = eps00;
    mu0 = mu00; //chapuz para convertir la variables de paso en globales

    OldInvEps.assign(OldInvEps.begin(), OldInvEps.begin() + sgg.NumMedia + 1);
    for (int i = 0; i <= sgg.NumMedia; ++i) {
        OldInvEps[i] = InvEps[i];
    }
    OldInvMu.assign(OldInvMu.begin(), OldInvMu.begin() + sgg.NumMedia + 1);
    for (int i = 0; i <= sgg.NumMedia; ++i) {
        OldInvMu[i] = InvMu[i];
    }
    InvEps.assign(sgg.NumMedia + 1, 0.0);
    for (int i = 0; i <= sgg.NumMedia; ++i) {
        InvEps[i] = 1.0_RKIND_wires / (Eps0 * sgg.Med[i].Epr);
    }
    InvMu.assign(sgg.NumMedia + 1, 0.0);
    for (int i = 0; i <= sgg.NumMedia; ++i) {
        InvMu[i] = 1.0_RKIND_wires / (Mu0 * sgg.Med[i].Mur);
    }

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Constantes nodales
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for (int i1 = 1; i1 <= HWires.NumChargeNodes; ++i1) {
        nodo = &HWires.ChargeNode[i1];
        nodo->CtePlain = nodo->CtePlain / Hwires.olddt * sgg.dt;
    }

    //ctes Mur
    for (int i1 = 1; i1 <= HWires.NumChargeNodes; ++i1) {
        nodo = &HWires.Chargenode[i1];
        if (nodo->IsBackDownLeftMur) {
            dummy = nodo->CurrentPlus_1;
            nodo->cteMur = (sqrt(InvMu(dummy->indexmed) * InvEps(dummy->indexmed)) * sgg.dt - dummy->delta) /
                           (sqrt(InvMu(dummy->indexmed) * InvEps(dummy->indexmed)) * sgg.dt + dummy->delta);
        } else if (nodo->IsFrontUpRightMur) {
            dummy = nodo->CurrentMInus_1;
            nodo->cteMur = (sqrt(InvMu(dummy->indexmed) * InvEps(dummy->indexmed)) * sgg.dt - dummy->delta) /
                           (sqrt(InvMu(dummy->indexmed) * InvEps(dummy->indexmed)) * sgg.dt + dummy->delta);
        }
    }

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Constantes de segmentos
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for (int i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
        //constantes de actualizacion
        dummy = &HWires.CurrentSegment[i1];
        jmed = dummy->indexmed;
        dummy->Lind = dummy->Lind * OldInvMu(jmed) / InvMu(jmed); //escalo la autoinduccion por invmu (aunque el mu no se escala nunca)
        dummy->Lind_devia = dummy->Lind_devia * OldInvMu(jmed) / InvMu(jmed); //escalo la autoinduccion por invmu (aunque el mu no se escala nunca)
        wiresconstantes(fieldtotl, dummy, G2, sgg);
    }
    //!!!!!!!!!!!!!!!!!!!!
    //Junctions segmento%FractionMinus & segmento%FractionPlus not affected by permittivity scaling

    //Solapes con conformal (niapa 171216) no modificados en %cte5

    //!!!!!!!!!!!!!!! Berenger transition !untested
    //!!!copiado 23/04/2014
    if (wiresflavor == "transition") {
        for (iw1 = 1; iw1 <= HWires.NumMultilines; ++iw1) {
            NumParallel = HWires.Multilines[iw1].NumParallel;
            for (is1 = 1; is1 <= NumParallel; ++is1) {
                imed = HWires.Multilines[iw1].Segments[is1].ptr->indexmed;
                HWires.Multilines[iw1].Segments[is1].ptr->Lintrinsic = HWires.Multilines[iw1].Segments[is1].ptr->Lintrinsic * OldInvMu(imed) / InvMu(imed);
                HWires.Multilines[iw1].Segments[is1].ptr->L = HWires.Multilines[iw1].Segments[is1].ptr->L * OldInvMu(imed) / InvMu(imed);
                HWires.Multilines[iw1].Segments[is1].ptr->C = HWires.Multilines[iw1].Segments[is1].ptr->C * OldInvEps(imed) / InvEps(imed);
                for (is2 = 1; is2 <= NumParallel; ++is2) {
                    HWires.Multilines[iw1].C[is1][is2] = HWires.Multilines[iw1].Segments[is1].ptr->C;
                    HWires.Multilines[iw1].L[is1][is2] = HWires.Multilines[iw1].Segments[is1].ptr->L;
                }
            }
            Den = new std::vector<std::vector<double>>(NumParallel + 1, std::vector<double>(NumParallel + 1));
            for (int r = 1; r <= NumParallel; ++r) {
                for (int c = 1; c <= NumParallel; ++c) {
                    (*Den)[r][c] = HWires.Multilines[iw1].L[r][c] + HWires.Multilines[iw1].R[r][c] * sgg.dt / 2.0_RKIND_wires;
                }
            }
            MatInv(NumParallel, *Den);

            HWires.Multilines[iw1].b1I = MatMul(*Den, HWires.Multilines[iw1].L - HWires.Multilines[iw1].R * sgg.dt / 2.0_RKIND_wires);
            dl = HWires.Multilines[iw1].Segments[1].ptr->delta; //tomo el de 1 !dama no lo tenia
            imed = HWires.Multilines[iw1].Segments[1].ptr->indexmed; //lo aniado yo pq no estaba definido sgg 141118 !tomo el de 1 !dama no lo tenia
            HWires.Multilines[iw1].b2I = -sgg.dt / dl * InvMu(imed) * InvEps(imed) * MatMul(*Den, HWires.Multilines[iw1].L);
            HWires.Multilines[iw1].b3I = sgg.dt * *Den;

            delete Den;
            Den = nullptr;

            for (is1 = 1; is1 <= NumParallel; ++is1) {
                HWires.Multilines[iw1].Segments[is1].ptr->bI = HWires.Multilines[iw1].b2I[is1][is1];
            }
        }
    }

    //!!!guarda el dt anterior

    Hwires.olddt = sgg.dt;

    return;
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//Returns info on the adjacency of two segments (first and second)
//all the wires (Parallel) are collapsed and everything is ohmicly connected if adjacent
//If -connectendings wires are not collapsed, but the endings are connected even if not specified by the .nfde
//Parallel coincident wires with coincident nodes are
//different wires with different nodes (later will be handled with the Berenger BerengerMTLN routine).
//segments are collapsed into ONLY 1 -any- and reported accordingly in the warnings file. Furthermore coincident
//nodes are connected ohmicly -are collapsed also- in this latter case.
//

bool TestAdjacency(CurrentSegments_t& first, int numfirst, CurrentSegments_t& second, int numsecond, bool connectendings, bool isolategroupgroups, bool strictOLD, int ZI, int ZE, int NUMESEG, CurrentSegments_t& firstmenos1, CurrentSegments_t& FIRSTMAS1, CurrentSegments_t& secondmenos1, CurrentSegments_t& secondmas1, bool verbose) {
    CurrentSegments_t& first, second, firstmen;
    // Function body not provided in the chunk, returning placeholder
    return false;
}

os1, firstmas1, secondmenos1, secondmas1;
    CurrentSegments_t* firstprevio = nullptr;
    CurrentSegments_t* secondprevio = nullptr;
    adyacc_t adj;
    bool verbose;
    bool conexionados, RRConnected, RLConnected, LLConnected, LRConnected, EndingsConnected, connectendings, isolategroupgroups, RequestedConnection, strictOLD, success, entro1, entro2, entro3;
    int numfirst, numsecond, ZI, ZE, void_var, offx, offy, offz, NUMESEG;
    std::array<std::string, 3> DIR;
    std::string buff;
    RequestedConnection = false;

    DIR[iEx] = " X ";
    DIR[iEy] = " Y ";
    DIR[iEz] = " Z ";
    adj.YESsegment[0] = -1;
    adj.YESsegment[1] = -1;
    adj.Is = false;
    adj.Parallel = false;
    adj.IsHeterogeneousJunction = false;
    adj.BothEndingsConnected = false;

    if (numfirst == numsecond) return; // trivial case: equal segment

    if (first->tipofield == second->tipofield) {
        if ((first->j == second->j) && (first->k == second->k) && (first->tipofield == iEx)) {
            adj.j = first->j;
            adj.k = first->k;
            if (first->i == second->i - 1) {
                adj.Is = true;
                adj.i = second->i;
            }
            if (first->i == second->i + 1) {
                adj.Is = true;
                adj.i = first->i;
            }
            if (first->i == second->i) {
                adj.Parallel = true;
                adj.i = first->i;
            }
        } else if ((first->i == second->i) && (first->k == second->k) && (first->tipofield == iEy)) {
            adj.i = first->i;
            adj.k = first->k;
            if (first->j == second->j - 1) {
                adj.Is = true;
                adj.j = second->j;
            }
            if (first->j == second->j + 1) {
                adj.Is = true;
                adj.j = first->j;
            }
            if (first->j == second->j) {
                adj.Parallel = true;
                adj.j = first->j;
            }
        } else if ((first->i == second->i) && (first->j == second->j) && (first->tipofield == iEz)) {
            adj.i = first->i;
            adj.j = first->j;
            if (first->k == second->k - 1) {
                adj.Is = true;
                adj.k = second->k;
            }
            if (first->k == second->k + 1) {
                adj.Is = true;
                adj.k = first->k;
            }
            if (first->k == second->k) {
                adj.Parallel = true;
                adj.k = first->k;
            }
        }
    } else if (((first->tipofield == iEx) && (second->tipofield == iEy)) &&
               (((first->i == second->i) && (first->j == second->j) && (first->k == second->k)) ||
                ((first->i == second->i) && (first->j == second->j + 1) && (first->k == second->k)))) {
        adj.Is = true;
        adj.i = first->i;
        adj.j = first->j;
        adj.k = first->k;
    } else if (((first->tipofield == iEx) && (second->tipofield == iEy)) &&
               (((first->i == second->i - 1) && (first->j == second->j) && (first->k == second->k)) ||
                ((first->i == second->i - 1) && (first->j == second->j + 1) && (first->k == second->k)))) {
        adj.Is = true;
        adj.i = first->i + 1;
        adj.j = first->j;
        adj.k = first->k;
    } else if (((first->tipofield == iEy) && (second->tipofield == iEz)) &&
               (((first->j == second->j) && (first->k == second->k) && (first->i == second->i)) ||
                ((first->j == second->j) && (first->k == second->k + 1) && (first->i == second->i)))) {
        adj.Is = true;
        adj.j = first->j;
        adj.k = first->k;
        adj.i = first->i;
    } else if (((first->tipofield == iEy) && (second->tipofield == iEz)) &&
               (((first->j == second->j - 1) && (first->k == second->k) && (first->i == second->i)) ||
                ((first->j == second->j - 1) && (first->k == second->k + 1) && (first->i == second->i)))) {
        adj.Is = true;
        adj.j = first->j + 1;
        adj.k = first->k;
        adj.i = first->i;
    } else if (((first->tipofield == iEz) && (second->tipofield == iEx)) &&
               (((first->k == second->k) && (first->i == second->i) && (first->j == second->j)) ||
                ((first->k == second->k) && (first->i == second->i + 1) && (first->j == second->j)))) {
        adj.Is = true;
        adj.k = first->k;
        adj.i = first->i;
        adj.j = first->j;
    } else if (((first->tipofield == iEz) && (second->tipofield == iEx)) &&
               (((first->k == second->k - 1) && (first->i == second->i) && (first->j == second->j)) ||
                ((first->k == second->k - 1) && (first->i == second->i + 1) && (first->j == second->j)))) {
        adj.Is = true;
        adj.k = first->k + 1;
        adj.i = first->i;
        adj.j = first->j;
    } else if (((first->tipofield == iEy) && (second->tipofield == iEx)) &&
               (((first->i == second->i) && (first->j == second->j) && (first->k == second->k)) ||
                ((first->i == second->i + 1) && (first->j == second->j) && (first->k == second->k)))) {
        adj.Is = true;
        adj.i = first->i;
        adj.j = first->j;
        adj.k = first->k;
    } else if (((first->tipofield == iEy) && (second->tipofield == iEx)) &&
               (((first->i == second->i) && (first->j == second->j - 1) && (first->k == second->k)) ||
                ((first->i == second->i + 1) && (first->j == second->j - 1) && (first->k == second->k)))) {
        adj.Is = true;
        adj.i = first->i;
        adj.j = first->j + 1;
        adj.k = first->k;
    } else if (((first->tipofield == iEx) && (second->tipofield == iEz)) &&
               (((first->k == second->k) && (first->i == second->i) && (first->j == second->j)) ||
                ((first->k == second->k + 1) && (first->i == second->i) && (first->j == second->j)))) {
        adj.Is = true;
        adj.k = first->k;
        adj.i = first->i;
        adj.j = first->j;
    } else if (((first->tipofield == iEx) && (second->tipofield == iEz)) &&
               (((first->k == second->k) && (first->i == second->i - 1) && (first->j == second->j)) ||
                ((first->k == second->k + 1) && (first->i == second->i - 1) && (first->j == second->j)))) {
        adj.Is = true;
        adj.k = first->k;
        adj.i = first->i + 1;
        adj.j = first->j;
    } else if (((first->tipofield == iEz) && (second->tipofield == iEy)) &&
               (((first->j == second->j) && (first->k == second->k) && (first->i == second->i)) ||
                ((first->j == second->j + 1) && (first->k == second->k) && (first->i == second->i)))) {
        adj.Is = true;
        adj.j = first->j;
        adj.k = first->k;
        adj.i = first->i;
    } else if (((first->tipofield == iEz) && (second->tipofield == iEy)) &&
               (((first->j == second->j) && (first->k == second->k - 1) && (first->i == second->i)) ||
                ((first->j == second->j + 1) && (first->k == second->k - 1) && (first->i == second->i)))) {
        adj.Is = true;
        adj.j = first->j;
        adj.k = first->k + 1;
        adj.i = first->i;
    }

    // added July'12 to connect requested endings
    if (adj.Parallel) adj.Is = true;
    // RightConnected adyacency foreseen in the .nfde file      // added to solve aNTh-THW_PW bug

    // aniadido 6/2/14 para tratar bien ESN.nfde (solo sirve para dar warnings a OLD con sus rabillos)
    if ((first->Is_LeftEnd) && (second->Is_LeftEnd)) {
        RequestedConnection = (!connectendings) && (first->indexmed != second->indexmed) &&
                              ((first->tipowire->LeftEnd == second->tipowire->LeftEnd));
    }
    if ((first->Is_LeftEnd) && (second->Is_RightEnd)) {
        RequestedConnection = (!connectendings) && (first->indexmed != second->indexmed) &&
                              ((first->tipowire->LeftEnd == second->tipowire->RightEnd));
    }
    if ((first->Is_RightEnd) && (second->Is_RightEnd)) {
        RequestedConnection = (!connectendings) && (first->indexmed != second->indexmed) &&
                              ((first->tipowire->RightEnd == second->tipowire->RightEnd));
    }
    if ((first->Is_RightEnd) && (second->Is_LeftEnd)) {
        RequestedConnection = (!connectendings) && (first->indexmed != second->indexmed) &&
                              ((first->tipowire->RightEnd == second->tipowire->LeftEnd));
    }

    // esto lo he quitado el 6/2/14 pq en ESN.nfde ambos extremos estaban adjancentes pero solo se pedia conexion en uno y por lo tanto hay que discernir con la condicional
    // de mas arriba
    // !!!!     RequestedConnection = (.not.connectendings).and.(first%indexmed /= second%indexmed).and. &
    // !!!!                                    ((first%tipowire%LeftEnd == second%tipowire%LeftEnd).or. &
    // !!!!                                     (first%tipowire%LeftEnd == second%tipowire%RightEnd).or. &
    // !!!!                                     (first%tipowire%RightEnd == second%tipowire%LeftEnd).or. &
    // !!!!                                     (first%tipowire%RightEnd == second%tipowire%RightEnd))  !solo para reporte bug OLD 12/09/13 Model_unidos en .nfde

    // quito cualquier referencia a LextremoI,j,k porque yo calculo bien los extremos libres tanto si es strictOLD como si no

    LLConnected = (first->indexmed != second->indexmed) &&
                  ((second->Is_LeftEnd && first->Is_LeftEnd) ||
                   (first->IsEnd_norLeft_norRight && second->IsEnd_norLeft_norRight) ||
                   (first->IsEnd_norLeft_norRight && second->Is_LeftEnd) ||
                   (first->Is_LeftEnd && second->IsEnd_norLeft_norRight));

    LRConnected = (first->indexmed != second->indexmed) &&
                  ((second->Is_RightEnd && first->Is_LeftEnd) ||
                   (first->IsEnd_norLeft_norRight && second->IsEnd_norLeft_norRight));

    RLConnected = (first->indexmed != second->indexmed) &&
                  ((second->Is_LeftEnd && first->Is_RightEnd) ||
                   (first->IsEnd_norLeft_norRight && second->IsEnd_norLeft_norRight));

    RRConnected = (first->indexmed != second->indexmed) &&
                  ((second->Is_RightEnd && first->Is_RightEnd) ||
                   (first->IsEnd_norLeft_norRight && second->IsEnd_norLeft_norRight) ||
                   (first->IsEnd_norLeft_norRight && second->Is_RightEnd) ||
                   (first->Is_RightEnd && second->IsEnd_norLeft_norRight));

    if ((first->Is_LeftEnd && first->Is_RightEnd) && (!(second->Is_LeftEnd && second->Is_RightEnd))) {
        conexionados =
            (((second->ilibre == first->ilibre) &&
              (second->jlibre == first->jlibre) &&
              (second->klibre == first->klibre)) ||
             ((second->ilibre == first->i) &&
              (second->jlibre == first->j) &&
              (second->klibre == first->k)));

        LLConnected = LLConnected && conexionados;
        LRConnected = LRConnected && conexionados;
        RLConnected = RLConnected && conexionados;
        RRConnected = RRConnected && conexionados;
    } else if ((second->Is_LeftEnd && second->Is_RightEnd) && (!(first->Is_LeftEnd && first->Is_RightEnd))) {
        conexionados =
            (((second->ilibre == first->ilibre) &&
              (second->jlibre == first->jlibre) &&
              (second->klibre == first->klibre)) ||
             ((second->i == first->ilibre) &&
              (second->j == first->jlibre) &&
              (second->k == first->klibre)));

        LLConnected = LLConnected && conexionados;
        LRConnected = LRConnected && conexionados;
        RLConnected = RLConnected && conexionados;
        RRConnected = RRConnected && conexionados;
    } else if ((second->Is_LeftEnd && second->Is_RightEnd) && (first->Is_LeftEnd && first->Is_RightEnd)) {
        conexionados =
            (((second->ilibre == first->ilibre) &&
              (second->jlibre == first->jlibre) &&
              (second->klibre == first->klibre)) ||
             ((second->i == first->ilibre) &&
              (second->j == first->jlibre) &&
              (second->k == first->klibre)) ||
             ((second->ilibre == first->i) &&
              (second->jlibre == first->j) &&
              (second->klibre == first->k)));

        LLConnected = LLConnected && conexionados;
        LRConnected = LRConnected && conexionados;
        RLConnected = RLConnected && conexionados;
        RRConnected = RRConnected && conexionados;
    } else {
        conexionados =
            (second->ilibre == first->ilibre) &&
            (second->jlibre == first->jlibre) &&
            (second->klibre == first->klibre);

        LLConnected = LLConnected && conexionados;
        LRConnected = LRConnected && conexionados;
        RLConnected = RLConnected && conexionados;
        RRConnected = RRConnected && conexionados;
    }

    if (!connectendings) LLConnected = LLConnected && (first->TIPOWIRE->LeftEnd == second->TIPOWIRE->LeftEnd);
    if (!connectendings) LRConnected = LRConnected && (first->TIPOWIRE->LeftEnd == second->TIPOWIRE->RightEnd);
    if (!connectendings) RLConnected = RLConnected && (first->TIPOWIRE->RightEnd == second->TIPOWIRE->LeftEnd);
    if (!connectendings) RRConnected = RRConnected && (first->TIPOWIRE->RightEnd == second->TIPOWIRE->RightEnd);

    EndingsConnected = (LLConnected || LRConnected || RLConnected || RRConnected);
    // casos especiales con un hilo de un segmento

    if (adj.Is && connectendings) {
        void_var = 0;
        if (second->Is_LeftEnd) void_var++;
        if (second->Is_RightEnd) void_var++;
        if (second->IsEnd_norLeft_norRight) void_var++;
        if ((first->Is_LeftEnd && first->Is_RightEnd) && (void_var == 1)) {
            if (
                ((first->i == second->ilibre) && (first->j == second->jlibre) && (first->k == second->klibre)) ||
                ((first->ilibre == second->ilibre) && (first->jlibre == second->jlibre) && (first->klibre == second->klibre))) {
                EndingsConnected = true;
                adj.i = second->ilibre;
                adj.j = second->jlibre;
                adj.k = second->klibre;
            }
        }
        void_var = 0;
        if (first->Is_LeftEnd) void_var++;
        if (first->Is_RightEnd) void_var++;
        if (first->IsEnd_norLeft_norRight) void_var++;
        if ((void_var == 1) && (second->Is_LeftEnd && second->Is_RightEnd)) {
            if (
                ((first->ilibre == second->i) && (first->jlibre == second->j) && (first->klibre == second->k)) ||
                ((first->ilibre == second->ilibre) && (first->jlibre == second->jlibre) && (first->klibre == second->klibre))) {
                EndingsConnected = true;
                adj.i = first->ilibre;
                adj.j = first->jlibre;
                adj.k = first->klibre;
            }
        }
    }

    // los cuatro conectados (puede pasar en amelet o en nfde)
    if ((first->Is_LeftEnd && first->Is_RightEnd) && (second->Is_LeftEnd && second->Is_RightEnd)) {
        if (connectendings) {
            if (
                ((first->i == second->ilibre) && (first->j == second->jlibre) && (first->k == second->klibre)) ||
                ((first->ilibre == second->ilibre) && (first->jlibre == second->jlibre) && (first->klibre == second->klibre))) {
                EndingsConnected = true;
                adj.i = second->ilibre;
                adj.j = second->jlibre;
                adj.k = second->klibre;
            } else if (
                ((first->ilibre == second->i) && (first->jlibre == second->j) && (first->klibre == second->k)) ||
                ((first->ilibre == second->ilibre) && (first->jlibre == second->jlibre) && (first->klibre == second->klibre))) {
                EndingsConnected = true;
                adj.i = first->ilibre;
                adj.j = first->jlibre;
                adj.k = first->klibre;
            } else if (
                ((first->i == second->i) && (first->j == second->j) && (first->k == second->k))) {
                EndingsConnected = true;
                adj.i = first->i;
                adj.j = first->j;
                adj.k = first->k;
            } else if (
                ((first->ilibre == second->ilibre) && (first->jlibre == second->jlibre) && (first->klibre == second->klibre))) {
                EndingsConnected = true;
                adj.i = first->ilibre;
                adj.j = first->jlibre;
                adj.k = first->klibre;
            }
            //
            adj.BothEndingsConnected = adj.Parallel;
        } else {
            EndingsConnected = (LLConnected || LRConnected || RLConnected || RRConnected); // redundante ya hecho antes
            adj.BothEndingsConnected = ((RRConnected && LLConnected) || (RLConnected && LRConnected));
        }
    }

    // fin casos especiales
    if (adj.Is) {
        if (first->indexmed != second->indexmed) {
            if (!EndingsConnected) {
                ADJ->IS = false;
                adj.IsHeterogeneousJunction = false;
                if (adj.Parallel) {
                    if ((second->Is_RightEnd || second->Is_LeftEnd || second->IsEnd_norLeft_norRight) &&
                        (first->Is_RightEnd || first->Is_LeftEnd || first->IsEnd_norLeft_norRight)) {
                        if (RequestedConnection) {
                            buff = "wir2_ERROR: Requested connection on non-connected Parallel Adjacent ENDING segments from multiWIREs:  ";
                            if ((first->k >= ZI) && (first->k <= ZE)) WarnErrReport(buff, true);
                            std::stringstream ss;
                            ss << first->origindex << " " << first->i << " " << first->j << " " << first->k << " " << DIR[first->tipofield] << " "
                               << second->origindex << " " << second->i << " " << second->j << " " << second->k << " " << DIR[second->tipofield];
                            buff = ss.str();
                            if ((first->k >= ZI) && (first->k <= ZE)) WarnErrReport(buff, true);
                        } else {
                            buff = "wir2_WARNING: DISCONNECTING Parallel Adjacent ENDING segments from multiWIREs:  ";
                            if ((first->k >= ZI) && (first->k <= ZE)) WarnErrReport(buff);
                            std::stringstream ss;
                            ss << first->origindex << " " << first->i << " " << first->j << " " << first->k << " " << DIR[first->tipofield] << " "
                               << second->origindex << " " << second->i << " " << second->j << " " << second->k << " " << DIR[second->tipofield];
                            buff = ss.str();
                            if ((first->k >= ZI) && (first->k <= ZE)) WarnErrReport(buff);
                        }
                    } else {
                        buff = "wir2_INFO: DISCONNECTING Parallel Adjacent intermediate segments from multiWIREs:  ";
                        if ((first->k >= ZI) && (first->k <= ZE) && verbose) WarnErrReport(buff);
                        std::stringstream ss;
                        ss << first->origindex << " " << first->i << " " << first->j << " " << first->k << " " << DIR[first->tipofield] << " "
                           << second->origindex << " " << second->i << " " << second->j << " " << second->k << " " << DIR[second->tipofield];
                        buff = ss.str();
                        if ((first->k >= ZI) && (first->k <= ZE) && verbose) WarnErrReport(buff);
                    }
                } else {
                    if ((second->Is_RightEnd || second->Is_LeftEnd || second->IsEnd_norLeft_norRight) &&
                        (first->Is_RightEnd || first->Is_LeftEnd || first->IsEnd_norLeft_norRight)) {
                        if (RequestedConnection) {
                            buff = "wir2_ERROR: Requested connection on non-connected Non-Parallel Adjacent ENDING segments from multiWIREs:  ";
                            if ((first->k >= ZI) && (first->k <= ZE)) WarnErrReport(buff);
                            std::stringstream ss;
                            ss << first->origindex << " " << first->i << " " << first->j << " " << first->k << " " << DIR[first->tipofield] << " "
                               << second->origindex << " " << second->i << " " << second->j << " " << second->k << " " << DIR[second->tipofield];
                            buff = ss.str();
                            if ((first->k >= ZI) && (first->k <= ZE)) WarnErrReport(buff, true);
                        } else {
                            buff = "wir2_WARNING: DISCONNECTING NON-Parallel Adjacent ENDING segments from multiWIREs:  ";
                            if ((first->k >= ZI) && (first->k <= ZE)) WarnErrReport(buff);
                            std::stringstream ss;
                            ss << first->origindex << " " << first->i << " " << first->j << " " << first->k << " " << DIR[first->tipofield] << " "
                               << second->origindex << " " << second->i << " " << second->j << " " << second->k << " " << DIR[second->tipofield];
                            buff = ss.str();
                            if ((first->k >= ZI) && (first->k <= ZE)) WarnErrReport(buff);
                        }
                    }
                }
            }
        }
    }

pofield),
                           second->origindex, second->i, second->j, second->k, " " + dir(second->tipofield));
                           if ((first->k >= ZI) && (first->k <= ZE)) WarnErrReport(buff);
                        }
                     }
                     else {
                        sprintf(buff, "wir2_INFO: DISCONNECTING Non-Parallel Adjacent intermediate segments from multiWIREs:");
                        if ((first->k >= ZI) && (first->k <= ZE) && verbose) WarnErrReport(buff);
                        sprintf(buff, "%7d%7d%7d%7d %s%7d%7d%7d%7d %s", first->origindex, first->i, first->j, first->k, dir(first->tipofield).c_str(),
                           second->origindex, second->i, second->j, second->k, dir(second->tipofield).c_str());
                        if ((first->k >= ZI) && (first->k <= ZE) && verbose) WarnErrReport(buff);
                     }
                  }
               }
               else if (EndingsConnected) {
                  if (adj->Parallel) {
                     sprintf(buff, "wir2_INFO: CONNECTING Parallel Adjacent ENDING segments from multiWIREs:  ");
                  }
                  else {
                     sprintf(buff, "wir2_INFO: CONNECTING Non-Parallel Adjacent ENDING segments from multiWIREs:");
                  }
                  if ((first->k >= ZI) && (first->k <= ZE) && verbose) WarnErrReport(buff);
                  sprintf(buff, "%7d%7d%7d%7d %s%7d%7d%7d%7d %s", first->origindex, first->i, first->j, first->k, dir(first->tipofield).c_str(),
                     second->origindex, second->i, second->j, second->k, dir(second->tipofield).c_str());
                  if ((first->k >= ZI) && (first->k <= ZE) && verbose) WarnErrReport(buff);
                  sprintf(buff, "           AT :%7d%7d%7d", adj->i, adj->j, adj->k);
                  if ((first->k >= ZI) && (first->k <= ZE) && verbose) WarnErrReport(buff);

                  adj->IsHeterogeneousJunction = true;
                  ADJ->IS = true;
                  adj->YESsegment(1) = numfirst;
                  adj->YESsegment(2) = numsecond;
                  if (LLConnected || LRConnected || RLConnected || RRConnected) {
                     // esto daba problemas con edel_wuking2.nfde y en general con el ultimo cable si solo tiene un solo segmento no unido, poniendo nodo%heterogeneousjunction=.true. al nodo libre en la vuelta (numfirst>numsecond)
                     // comentado, pues, el 141218 y aniadido codigo 
                     // adj->i = first->ilibre;
                     // adj->j = first->jlibre;
                     // adj->k = first->klibre;
                     entro1 = false; entro2 = false; entro3 = false;
                     // lo de antes descomentado no daba problemas OLD. En correo 150219 entro el primer buggy de abajo. Por tanto aniado mas
                     // casuistica. Ahora puede que falle el edel_wuking2.nfde por este cambio, pero prefiero no descomentar lo anterior. Arreglar cuando pase !!
                     if (
                        ((first->i == second->ilibre) && (first->j == second->jlibre) && (first->k == second->klibre))) {
                        adj->i = second->ilibre;
                        adj->j = second->jlibre;
                        adj->k = second->klibre;
                        entro1 = true;
                     }
                     if (
                        ((first->ilibre == second->i) && (first->jlibre == second->j) && (first->klibre == second->k))) {
                        adj->i = first->ilibre;
                        adj->j = first->jlibre;
                        adj->k = first->klibre;
                        entro2 = true;
                     }
                     // esta es la casuistica que aniado 150219                   
                     if (
                        ((first->ilibre == second->ilibre) && (first->jlibre == second->jlibre) && (first->klibre == second->klibre))) {
                        adj->i = second->ilibre;
                        adj->j = second->jlibre;
                        adj->k = second->klibre;
                        entro3 = true;
                     }
                     // fin casuistica que anidado 150219
                     if ((!entro1) && (!entro2) && (!entro3)) {
                        sprintf(buff, "wir2_BUGGYERROR: bug in adjacencies entro. %d%d%d%d%d%d", first->ilibre, first->jlibre, first->klibre, second->ilibre, second->jlibre, second->klibre);
                        WarnErrReport(buff, true);
                     }
                     if (entro1 && entro2) {
                        if (!entro3) {
                           sprintf(buff, "wir2_BUGGYERROR: bug in adjacencies entro. %d%d%d%d%d%d", first->ilibre, first->jlibre, first->klibre, second->ilibre, second->jlibre, second->klibre);
                           WarnErrReport(buff, true);
                        }
                     }
                     // fina aniadido 141218
                  }
               }
            }
            ELSEIF (First->indexmed == second->indexmed) {
               if (adj->Parallel) { // paralelos del mismo hilo
                  if (!strictOLD) {
                     adj->is = false;
                     // NUNCA DEBERIA ENTRAR AQUI porque los paralelos se han colapsado en la version no estricta
                     sprintf(buff, "wir2_BUGGYERROR: DISCONNECTING Parallel segments from the same WIRE:");
                     WarnErrReport(buff, true);
                     sprintf(buff, "%7d%7d%7d%7d %s%7d%7d%7d%7d %s", first->origindex, first->i, first->j, first->k, dir(first->tipofield).c_str(),
                        second->origindex, second->i, second->j, second->k, dir(second->tipofield).c_str());
                     WarnErrReport(buff, true);
                  }
                  else {
                     // !!!!!!!!!!aqui es donde viene el meollo porque NO he quitado segmentos repetidos
                     // voy a suponer que forman una cadena ordenada

                     if (abs(numfirst - numsecond) > 1) {
                        adj->is = false;
                        adj->IsHeterogeneousJunction = false;
                        sprintf(buff, "wir2_INFO: DISCONNECTING NON-CORRELATIVE Parallel segments from the same WIRE:");
                        if (verbose) WarnErrReport(buff);
                        sprintf(buff, "%7d%7d%7d%7d %s%7d%7d%7d%7d %s", first->origindex, first->i, first->j, first->k, dir(first->tipofield).c_str(),
                           second->origindex, second->i, second->j, second->k, dir(second->tipofield).c_str());
                        if (verbose) WarnErrReport(buff);
                     }
                     else {
                        sucCess = false;
                        firstprevio = nullptr;
                        secondprevio = nullptr;
                        if (firstmenos1 != nullptr) {
                           if (!((firstmenos1->i == first->i) && (firstmenos1->j == first->j) && (firstmenos1->k == first->k) && (firstmenos1->tipofield == first->tipofield))) {
                              firstprevio = firstmenos1;
                           }
                           else {
                              if (firstmas1 != nullptr) {
                                 if (!((firstmas1->i == first->i) && (firstmas1->j == first->j) && (firstmas1->k == first->k) && (firstmas1->tipofield == first->tipofield))) {
                                    firstprevio = firstmas1;
                                 }
                              }
                           }
                        }
                        else {
                           if (firstmas1 != nullptr) {
                              if (!((firstmas1->i == first->i) && (firstmas1->j == first->j) && (firstmas1->k == first->k) && (firstmas1->tipofield == first->tipofield))) {
                                 firstprevio = firstmas1;
                              }
                           }
                        }
                        if (secondmenos1 != nullptr) {
                           if (!((secondmenos1->i == second->i) && (secondmenos1->j == second->j) && (secondmenos1->k == second->k) && (secondmenos1->tipofield == second->tipofield))) {
                              secondprevio = secondmenos1;
                           }
                           else {
                              if (secondmas1 != nullptr) {
                                 if (!((secondmas1->i == second->i) && (secondmas1->j == second->j) && (secondmas1->k == second->k) && (secondmas1->tipofield == second->tipofield))) {
                                    secondprevio = secondmas1;
                                 }
                              }
                           }
                        }
                        else {
                           if (secondmas1 != nullptr) {
                              if (!((secondmas1->i == second->i) && (secondmas1->j == second->j) && (secondmas1->k == second->k) && (secondmas1->tipofield == second->tipofield))) {
                                 secondprevio = secondmas1;
                              }
                           }
                        }
                        //
                        offx = 0; offy = 0; offz = 0;
                        if ((first->tipofield == second->tipofield) && (first->tipofield == iEx)) offx = 1;
                        if ((first->tipofield == second->tipofield) && (first->tipofield == iEy)) offy = 1;
                        if ((first->tipofield == second->tipofield) && (first->tipofield == iEz)) offz = 1;
                        if (firstprevio != nullptr) {
                           switch (first->tipofield) {
                           case iEx:
                              if (firstprevio->i == first->i) {
                                 success = true;
                                 adj->i = first->i + 1;
                                 adj->j = first->j;
                                 adj->k = first->k;
                              }
                              else if ((firstprevio->i - first->i) == 1) {
                                 success = true;
                                 adj->i = first->i;
                                 adj->j = first->j;
                                 adj->k = first->k;
                              }
                              else if ((firstprevio->i - first->i) == -1) {
                                 success = true;
                                 adj->i = first->i + offx;
                                 adj->j = first->j;
                                 adj->k = first->k;
                              }
                              break;
                           case iEy:
                              if (firstprevio->j == first->j) {
                                 success = true;
                                 adj->i = first->i;
                                 adj->j = first->j + 1;
                                 adj->k = first->k;
                              }
                              else if ((firstprevio->j - first->j) == 1) {
                                 success = true;
                                 adj->i = first->i;
                                 adj->j = first->j;
                                 adj->k = first->k;
                              }
                              else if ((firstprevio->j - first->j) == -1) {
                                 success = true;
                                 adj->i = first->i;
                                 adj->j = first->j + offy;
                                 adj->k = first->k;
                              }
                              break;
                           case iEz:
                              if (firstprevio->k == first->k) {
                                 success = true;
                                 adj->i = first->i;
                                 adj->j = first->j;
                                 adj->k = first->k + 1;
                              }
                              else if ((firstprevio->k - first->k) == 1) {
                                 success = true;
                                 adj->i = first->i;
                                 adj->j = first->j;
                                 adj->k = first->k;
                              }
                              else if ((firstprevio->k - first->k) == -1) {
                                 success = true;
                                 adj->i = first->i;
                                 adj->j = first->j;
                                 adj->k = first->k + offz;
                              }
                              break;
                           }
                        }
                        else if (secondprevio != nullptr) {
                           switch (second->tipofield) {
                           case iEx:
                              if (secondprevio->i == second->i) {
                                 success = true;
                                 adj->i = second->i + 1;
                                 adj->j = second->j;
                                 adj->k = second->k;
                              }
                              else if ((secondprevio->i - second->i) == 1) {
                                 success = true;
                                 adj->i = second->i;
                                 adj->j = second->j;
                                 adj->k = second->k;
                              }
                              else if ((secondprevio->i - second->i) == -1) {
                                 success = true;
                                 adj->i = second->i + offx;
                                 adj->j = second->j;
                                 adj->k = second->k;
                              }
                              break;
                           case iEy:
                              if (secondprevio->j == second->j) {
                                 success = true;
                                 adj->i = second->i;
                                 adj->j = second->j + 1;
                                 adj->k = second->k;
                              }
                              else if ((secondprevio->j - second->j) == 1) {
                                 success = true;
                                 adj->i = second->i;
                                 adj->j = second->j;
                                 adj->k = second->k;
                              }
                              else if ((secondprevio->j - second->j) == -1) {
                                 success = true;
                                 adj->i = second->i;
                                 adj->j = second->j + offy;
                                 adj->k = second->k;
                              }
                              break;
                           case iEz:
                              if (secondprevio->k == second->k) {
                                 success = true;
                                 adj->i = second->i;
                                 adj->j = second->j;
                                 adj->k = second->k + 1;
                              }
                              else if ((secondprevio->k - second->k) == 1) {
                                 success = true;
                                 adj->i = second->i;
                                 adj->j = second->j;
                                 adj->k = second->k;
                              }
                              else if ((secondprevio->k - second->k) == -1) {
                                 success = true;
                                 adj->i = second->i;
                                 adj->j = second->j;
                                 adj->k = second->k + offz;
                              }
                              break;
                           }
                        }
                        else {
                           if ((first->Is_LeftEnd || first->Is_RightEnd) && (second->Is_LeftEnd || second->Is_RightEnd)) {
                              success = true;
                              adj->i = second->i + offx;
                              adj->j = second->j + offy;
                              adj->k = second->k + offz;
                           }
                        }
                        if (success) {
                           adj->IsHeterogeneousJunction = false;
                           ADJ->IS = true;
                           adj->YESsegment(1) = numfirst;
                           adj->YESsegment(2) = numsecond;
                           sprintf(buff, "wir2_INFO: CONNECTING CORRELATIVE Parallel segments from the same WIRE (rabitos):");
                           if (verbose) WarnErrReport(buff);
                           sprintf(buff, "%7d%7d%7d%7d %s%7d%7d%7d%7d %s", first->origindex, first->i, first->j, first->k, dir(first->tipofield).c_str(),
                              second->origindex, second->i, second->j, second->k, dir(second->tipofield).c_str());
                           if (verbose) WarnErrReport(buff);
                           sprintf(buff, "           AT :%7d%7d%7d", adj->i, adj->j, adj->k);
                           if (verbose) WarnErrReport(buff);
                        }
                        else {
                           adj->IsHeterogeneousJunction = false;
                           ADJ->IS = false;
                           sprintf(buff, "wir2_BUGGYERROR:  Cannot determine point of contact of Parallel intra-WIRE segment connection (mas de dos rabitos doblados?). ");
                           WarnErrReport(buff, true);
                           sprintf(buff, "%7d%7d%7d%7d %s%7d%7d%7d%7d %s", first->origindex, first->i, first->j, first->k, dir(first->tipofield).c_str(),
                              second->origindex, second->i, second->j, second->k, dir(second->tipofield).c_str());
                           WarnErrReport(buff, true);
                        }
                        //
                     }
                  } // DEL NO STRICTO
               }
               else { // no son paralelos pero si estan en el mismo hilo
                  if (!strictOLD) {
                     adj->IsHeterogeneousJunction = false;
                     ADJ->IS = true;
                     adj->YESsegment(1) = numfirst;
                     adj->YESsegment(2) = numsecond;
                     if (abs(numfirst - numsecond) > 1) {
                        sprintf(buff, "wir2_INFO: CONNECTING NON-CORRELATIVE Non-Parallel segments from the same WIRE:");
                        if (verbose) WarnErrReport(buff);
                     }
                     else {
                        sprintf(buff, "wir2_INFO: CONNECTING CORRELATIVE Non-Parallel segments from the same WIRE:");
                        if (verbose) WarnErrReport(buff);
                     }
                     sprintf(buff, "%7d%7d%7d%7d %s%7d%7d%7d%7d %s", first->origindex, first->i, first->j, first->k, dir(first->tipofield).c_str(),
                        second->origindex, second->i, second->j, second->k, dir(second->tipofield).c_str());
                     WarnErrReport(buff);
                     sprintf(buff, "           AT :%7d%7d%7d", adj->i, adj->j, adj->k);
                     WarnErrReport(buff);
                  }
                  else {
                     if (abs(numfirst - numsecond) > 1) { // solo si estan leidos contiguamente se toman como adyacentes
                        adj->IsHeterogeneousJunction = false;
                        ADJ->IS = false;
                        sprintf(buff, "wir2_INFO: DISCONNECTING NON-CORRELATIVE Non-Parallel segments from the same WIRE:");
                        if (verbose) WarnErrReport(buff); // demasiado verbose
                        sprintf(buff, "%7d%7d%7d%7d %s%7d%7d%7d%7d %s", first->origindex, first->i, first->j, first->k, dir(first->tipofield).c_str(),
                           second->origindex, second->i, second->j, second->k, dir(second->tipofield).c_str());
                        if (verbose) WarnErrReport(buff);
                     }
                     else {
                        adj->IsHeterogeneousJunction = false;
                        ADJ->IS = true;
                        adj->YESsegment(1) = numfirst;
                        adj->YESsegment(2) = numsecond;
                        sprintf(buff, "wir2_INFO: CONNECTING CORRELATIVE Non-Parallel segments from the same WIRE:");
                        if (verbose) WarnErrReport(buff); // demasiado verbose
                        sprintf(buff, "%7d%7d%7d%7d %s%7d%7d%7d%7d %s", first->origindex, first->i, first->j, first->k, dir(first->tipofield).c_str(),
                           second->origindex, second->i, second->j, second->k, dir(second->tipofield).c_str());
                        if (verbose) WarnErrReport(buff);
                        sprintf(buff, "           AT :%7d%7d%7d", adj->i, adj->j, adj->k);
                        if (verbose) WarnErrReport(buff);
                     }
                  }
               }
            }
         }
      }
      //
      // !!!!!!!!!isolategroupsgroups
      if (isolategroupgroups) {
         if ((adj->is) && (first->tipowire->LeftEnd != second->tipowire->LeftEnd)) {
            sprintf(buff, "wir2_WARNING: DISCONNECTING PREV

IOUSLY Connected segments from multiWIREs being in DIFFERENT GROUPGROUPS:  '
            if ((first%k >= ZI).and.(first%k <= ZE)) call WarnErrReport(buff)
            write (buff,'(i7,3i7,a,i7,3i7,a)') first%tipowire%LeftEnd ,first%origindex,first%i,first%j,first%k,' '//dir(first%tipofield),&
            second%tipowire%LeftEnd,second%origindex,second%i,second%j,second%k,' '//dir(second%tipofield)
            if ((first%k >= ZI).and.(first%k <= ZE)) call WarnErrReport(buff)
            adj%is=.false.
            adj%BothEndingsConnected=.false.
         end if
      end if
      !!!!!!!!!!
      return
   end function TestAdjacency



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Advancing charge and current routine
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   subroutine AdvanceWiresE(sgg,timeinstant, layoutnumber,wiresflavor,simu_devia,stochastic,experimentalVideal,wirethickness,eps0,mu0)
                    
      real(kind=RKIND), intent(in) :: eps0,mu0
      integer, intent(in) :: wirethickness
      logical :: simu_devia,stochastic,experimentalVideal
      type(SGGFDTDINFO_t), intent(in) :: sgg

      integer(kind=4) :: n,jmed,layoutnumber,iw1,is1,is2

      integer(kind=4), intent(in) :: timeinstant
      real(kind=RKIND_wires) :: Iplus,IMinus,Qplus,QMinus,timei
      real(kind=RKIND_wires) :: Vincid,Iincid
      type(CurrentSegments_t), pointer  :: Segmento, Segmento2
      type(ChargeNodes_t), pointer  :: Nodo
      type(TMultiline_t), pointer                      :: Multiline
      character(len=*), INTENT(in) :: wiresflavor
      timei = sgg%tiempo(timeinstant) 
      !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !FIRST ADVANCE THE CHARGE from n+1.0_RKIND_wires / 2 to n+3/2 using the current known at n+1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !it is important that charge is first advanced for time-stepping coherency
      !Crucial also for MPI correctness!!!,
      !Beware than for MPI some currents outside the layout are also calculated with wrong information from charges, but
      !the MPI exchange call will overwrite them with the correct ones coming from the adjacent layers
      !WARNING: MPI does not handle correctly PERIODIC MIRRORING !\E7 20sept11 in currents !to dooooooooo

!!!!!!!!!!!!!!!051223 no pongo a cero los already_YEEadvanced_byconformal_changedtoPECfield !No son cero, sino lo que se updatee mas lo que se inyecte
!!!#ifdef CompileWithOpenMP
!!!!$OMP PARALLEL do DEFAULT(SHARED)  private(Nodo)
!!!#endif
!!!      do n=1,HWires%NumChargeNodes
!!!         Nodo => HWires%ChargeNode(n)
!!!!!!!140220 pon a PEC los viejos notouch=already_YEEadvanced_byconformal_changedtoPECfield que tenga conectado un nodo
!!!         !debe ser lo primero que se hace para overridear el call conformal_advance_E
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield1)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield1=0.0_RKIND
!!!         end if
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield2)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield2=0.0_RKIND
!!!         end if
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield3)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield3=0.0_RKIND
!!!         end if
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield4)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield4=0.0_RKIND
!!!         end if
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield5)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield5=0.0_RKIND
!!!         end if
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield6)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield6=0.0_RKIND
!!!         end if
!!!      end do
!!!#ifdef CompileWithOpenMP
!!!!$OMP END PARALLEL DO
!!!#endif
!
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED)  private(Iplus,IMinus,Nodo)
#endif
      do n=1,HWires%NumChargeNodes
         Nodo => HWires%ChargeNode(n)
!
         if (nodo%exists) then
            Nodo%ChargePast=Nodo%ChargePresent
            !I sum up the plus and the minus currents accordingly in advance
            Iplus = 0.0_RKIND_wires
            IMinus = 0.0_RKIND_wires
            If (Nodo%NumCurrentPlus>=1) Iplus = Nodo%CurrentPlus_1%Current
            If (Nodo%NumCurrentPlus>=2) Iplus = Iplus+Nodo%CurrentPlus_2%Current
            If (Nodo%NumCurrentPlus>=3) Iplus = Iplus+Nodo%CurrentPlus_3%Current
            If (Nodo%NumCurrentPlus>=4) Iplus = Iplus+Nodo%CurrentPlus_4%Current
            If (Nodo%NumCurrentPlus>=5) Iplus = Iplus+Nodo%CurrentPlus_5%Current
            If (Nodo%NumCurrentPlus>=6) Iplus = Iplus+Nodo%CurrentPlus_6%Current
            If (Nodo%NumCurrentPlus>=7) Iplus = Iplus+Nodo%CurrentPlus_7%Current
            If (Nodo%NumCurrentPlus>=8) Iplus = Iplus+Nodo%CurrentPlus_8%Current
            If (Nodo%NumCurrentPlus>=9) Iplus = Iplus+Nodo%CurrentPlus_9%Current
            !
            If (Nodo%NumCurrentMinus>=1) IMinus = Nodo%CurrentMinus_1%Current
            If (Nodo%NumCurrentMinus>=2) IMinus = IMinus+Nodo%CurrentMinus_2%Current
            If (Nodo%NumCurrentMinus>=3) IMinus = IMinus+Nodo%CurrentMinus_3%Current
            If (Nodo%NumCurrentMinus>=4) IMinus = IMinus+Nodo%CurrentMinus_4%Current
            If (Nodo%NumCurrentMinus>=5) IMinus = IMinus+Nodo%CurrentMinus_5%Current
            If (Nodo%NumCurrentMinus>=6) IMinus = IMinus+Nodo%CurrentMinus_6%Current
            If (Nodo%NumCurrentMinus>=7) IMinus = IMinus+Nodo%CurrentMinus_7%Current
            If (Nodo%NumCurrentMinus>=8) IMinus = IMinus+Nodo%CurrentMinus_8%Current
            If (Nodo%NumCurrentMinus>=9) IMinus = IMinus+Nodo%CurrentMinus_9%Current

            if ((Nodo%NumCurrentMinus == 1).and.(Nodo%NumCurrentPlus == 0)) then
                if (Nodo%IsPeriodic) then
                    Iplus = + Iminus
                else
            !The node is a true terminal one and the mirror of the current is employed for updating (Edelvik's treatment similar to PMC)
                    Iplus = -Iminus
                end if
            end if
            if ((Nodo%NumCurrentMinus == 0).and.(Nodo%NumCurrentPlus == 1)) then
                if (Nodo%IsPeriodic) then
                    IMinus = +Iplus
                else
                    IMinus = -Iplus
                end if
            end if
            !Algoritmo comun a toda la casuistica !feb 13
            if (.not.nodo%IsMur)   then
               Nodo%ChargePresent =  Nodo%CteProp*Nodo%ChargePast   - Nodo%CtePlain*(Iplus-IMinus)
            else
                continue
            end if
            !        if (nodo%IsAttachedtoVoltage) Nodo%ChargePresent =  0.0_RKIND_wires

         end if !del nodo%exists !voideo algunos cuando spliteo
      end do
!!!no toco lo que esta comentado a 230323 para fuentes duras/blandas
      !!!las convierto en duras mas abajo. correos jag rayos junio'15
      !!!Transparent current source feeding in the charge nodes if necessary
      !!!if (.not.simu_devia) then              
      !!!if (thereAreIsources) then
      !!    do n=1,HWires%NumChargeNodes
      !!        if (HWires%ChargeNode(n)%exists) then
      !!        If  (HWires%ChargeNode(n)%HasIsource) then
      !!            Nodo => HWires%ChargeNode(n)
      !!            Iincid=evolucion(timei-unmedio*sgg%dt,Nodo%Isource%Fichero%Samples, &
      !!                              Nodo%Isource%Fichero%DeltaSamples,Nodo%Isource%Fichero%NumSamples)
      !!            Nodo%ChargePresent                         =Nodo%ChargePresent                          +  Nodo%CtePlain                          * Iincid
      !!!!!!!!no funciona esta idea previa de fuente de corriente  quizas porque violo kirchhoff
      !!            Nodo%CurrentPlus_1%Chargeplus%ChargePresent=Nodo%CurrentPlus_1%Chargeplus%ChargePresent -  Nodo%CurrentPlus_1%Chargeplus%CtePlain * Iincid
      !!        end if
      !!        end if
      !!    end do
      !!!end if
      !!!end if
      
!aniado eso a 230323 para fuentes de corrientes duras blandas. Me inspiro en lo que esta comentado antes   
      if (.not.simu_devia) then              
          if (thereAreIsources) then
              do n=1,HWires%NumChargeNodes
                  if (HWires%ChargeNode(n)%exists) then
                     If  (HWires%ChargeNode(n)%HasIsource) then
                        Nodo => HWires%ChargeNode(n)
                        Iincid=evolucion(timei-unmedio*sgg%dt,Nodo%Isource%Fichero%Samples, &
                                        Nodo%Isource%Fichero%DeltaSamples,Nodo%Isource%Fichero%NumSamples)
                        Nodo%ChargePresent = Nodo%ChargePresent    +  Nodo%CtePlain     * Iincid
                     end if
                  end if
              end do
          end if
      end if

      !Absorbing Mur boundary conditions if necessary in the charges
      if (thereAreMurConditions) then
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED)  private(Nodo)
#endif
         do n=1,HWires%NumChargeNodes
            Nodo => HWires%ChargeNode(n)
            if (nodo%exists.and.nodo%IsMur) then
               Nodo%ChargePresent=Nodo%NodeInside%ChargePast + Nodo%cteMur*(Nodo%NodeInside%ChargePresent - Nodo%ChargePast)
            end if
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !END OF CHARGE ADVANCING from n+1.0_RKIND_wires / 2 to n+3/2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Inject the current at n+1 into the electring FDTD field previously calculated at n+1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !actualizo bien el campo con los globales (hilos paralelos)
      !I sum up all the currents (in case of unshielded segment)
      !this is the correct approach for Parallel wires (later to be corrected for BerengerMTLN)
      !

      !voids the null_field vale for embedded and paralel segments (if requested at run time with -groundwires and without the -stableradholland)
      HWires%null_field = 0.0_RKIND_wires
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do n=1,HWires%NumCurrentSegments
         Segmento => HWires%CurrentSegment(n)
         if (.not.Segmento%IsShielded) then
!171216quitado            Segmento%Efield_wire2main_past = real(Segmento%Efield_wire2main,kind=RKIND_wires)
#ifdef CompileWithThickWires
             if (wirethickness>1) then
                call Advance_Thick_Efield_wire2main(sgg,segmento,eps0,mu0,wirethickness)
             end if 
#endif                      
             if (wirethickness==1) then
                segmento%efield_wire2main=real(segmento%efield_wire2main,kind=rkind_wires) &
                    - segmento%cte5 * segmento%current
             end if
         end if
      end do

      !revoids the null_field value (unvoided above) for embedded and paralel segments (if requested at run time with -groundwires)
      HWires%null_field = 0.0_RKIND_wires
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !CURRENT ADVANCING from n+1 to n+2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !

      if (trim(adjustl(wiresflavor))=='transition') then
         !recuerdo intensidades
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED)   private(Segmento)
#endif
         do is1 = 1,HWires%NumCurrentSegments
            Segmento => HWires%CurrentSegment(is1)
            Segmento%CurrentPast = Segmento%Current
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif

#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED)   private(Multiline,Segmento,Segmento2,Qplus,QMinus)
#endif
         do iw1 = 1,HWires%NumMultilines
            Multiline => HWires%Multilines(iw1)
            do is1 = 1,Multiline%NumParallel
               Segmento => Multiline%Segments(is1)%ptr
               Segmento%Current = 0.0_RKIND_wires
               do is2 = 1,Multiline%NumParallel
                  Segmento2 => Multiline%Segments(is2)%ptr
                  Qplus  = Segmento2%ChargePlus%ChargePresent
                  QMinus = Segmento2%ChargeMinus%ChargePresent
                  Segmento2 => Multiline%Segments(is2)%ptr
                  Segmento%Current = Segmento%Current                                         + &
                  Multiline%b1I(is1,is2)* Segmento2%CurrentPast            + &
                  Multiline%b2I(is1,is2)*(Segmento2%fractionPlus*Qplus-Segmento2%fractionMinus*QMinus)
                  if(.not.(Segmento%IsShielded.and.Segmento2%IsShielded)) then
                      !!!lo he descomentado a 300523 porque creo que estaba mal. ojo si algun dia se usa el flavor transition. pongo un stop para avisar
                     Segmento%Current = Segmento%Current + Multiline%b3I(is1,is2)*real(Segmento2%Efield_main2wire,kind=RKIND_wires)
                     stop
                  end if
               end do
            end do
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         !!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!
      elseif (trim(adjustl(wiresflavor))=='holland') then
         !!!!!!!!!!!!!!!!!!!
         !!!!!!MY FLAVOR!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED)   private(Segmento,Qplus,QMinus,jmed)
#endif
         do n=1,HWires%NumCurrentSegments
            Segmento => HWires%CurrentSegment(n)
            Segmento%CurrentPast = Segmento%Current !save this data only for the probes observation right time positioning
            !--->
            if (Segmento%IsPMC) then
               Segmento%Current=0.0_RKIND_wires
            else
               Qplus  = Segmento%ChargePlus%ChargePresent
               QMinus = Segmento%ChargeMinus%ChargePresent
               Segmento%qplus_qminus=Segmento%fractionPlus*Qplus-Segmento%fractionMinus*QMinus
               Segmento%Current=Segmento%cte1*Segmento%Current - Segmento%cte3*(Segmento%qplus_qminus)
               if (.not.Segmento%IsShielded) then
#ifdef CompileWithThickWires
                   if (wirethickness>1) then
                       call Advance_Thick_Efield_main2wire(sgg,segmento,eps0,mu0,wirethickness)
                   end if 
#endif                            
                   if (wirethickness==1) then
                     Segmento%Current = Segmento%Current + &
                         Segmento%cte2*real(Segmento%Efield_main2wire,kind=RKIND_wires)
                   end if           
               end if
            end if
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         !!!!!!!!!!!!
         !!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!
      else
         call WarnErrReport('wir0_BUGGYERROR: this flavor not permitted in sgg routine',.true.)
      end if !wiresflavor
      !igual en ambos
      if (.not.simu_devia) then             
          if (thereAreVsources) then
             !Feed the transparent voltage sources as appropriate charges in the current segments (Like in Edelvik's model)
             do n=1,HWires%NumCurrentSegments
                If  (HWires%CurrentSegment(n)%HasVsource) then
                   Segmento => HWires%CurrentSegment(n)
                   Vincid=evolucion(timei,Segmento%Vsource%Fichero%Samples, &
                   Segmento%Vsource%Fichero%DeltaSamples,Segmento%Vsource%Fichero%NumSamples)
                   !!!if (experimentalVideal) then
                   !!!    if ((.not.Segmento%ChargePlus%isPEC).and.Segmento%ChargeMinus%isPEC) then
                   !!!         Segmento%ChargePlus%ChargePresent = Vincid  / (Segmento%Lind * InvMu(Segmento%indexmed)*InvEps(Segmento%indexmed))
                   !!!    elseif ((.not.Segmento%ChargeMinus%isPEC).and.Segmento%ChargePlus%isPEC) then
                   !!!         Segmento%ChargeMinus%ChargePresent =  Vincid  / (Segmento%Lind * InvMu(Segmento%indexmed)*InvEps(Segmento%indexmed))
                   !!!    else
                   !!!        print *,'error en experimentalVideal 200621'
                   !!!    end if
                   !!!else
                  !lo de siempre. aniado lo anterior para ver lo de las fuentes duras !C=Q/V-> C=c^-2/L-> Q=V*C=V/(L c^2)
                  !fuentes de voltaje. para que sean de carga hay que comentar el Lind
                  Segmento%Current = Segmento%Current + &
                     Segmento%cte3 * Vincid  / (Segmento%Lind * InvMu(Segmento%indexmed)*InvEps(Segmento%indexmed))
                  !I use the capacitance to find the incident charge
                  !assuming that the evolution file contains a voltage, not a charge
                   !!!end if
                end if
             end do
          end if
      end if
      
!!!230323 comento lo que sigue porque ya manejo fuentes duras y blandas antes, pero esto funciono con jag rayos junio'15
      !!!!!!080615  uso una fuente dura de corriente correos jag simulacion rayos Junio'15
      !!!if (.not.simu_devia) then             
      !!!    if (thereAreIsources) then
      !!!       do n=1,HWires%NumChargeNodes
      !!!          if (HWires%ChargeNode(n)%exists) then
      !!!             If  (HWires%ChargeNode(n)%HasIsource) then
      !!!                Nodo => HWires%ChargeNode(n)
      !!!                Iincid=evolucion(timei,Nodo%Isource%Fichero%Samples, &
      !!!                Nodo%Isource%Fichero%DeltaSamples,Nodo%Isource%Fichero%NumSamples)
      !!!                Nodo%CurrentPlus_1%Current = Iincid
      !!!             end if
      !!!          end if
      !!!       end do
      !!!    end if
      !!!end if
      !!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !END OF CURRENT ADVANCING from n+1 to n+2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!
!stochastic
#ifdef CompileWithStochastic
      if (stochastic.and.simu_devia) call inject_deviasources(layoutnumber,Hwires) !solo son los segmentos los afectados
#endif
 
!!!!!!!!!!!!!!!051223 no pongo a cero los already_YEEadvanced_byconformal_changedtoPECfield !No son cero, sino lo que se updatee mas lo que se inyecte
!!!machaca el campo que haya metido el wires

!!!#ifdef CompileWithOpenMP
!!!!$OMP PARALLEL do DEFAULT(SHARED)  private(Nodo)
!!!#endif
!!!      do n=1,HWires%NumChargeNodes
!!!         Nodo => HWires%ChargeNode(n)
!!!!!!!140220 pon a PEC los viejos notouch=already_YEEadvanced_byconformal_changedtoPECfield que tenga conectado un nodo
!!!         !debe ser lo primero que se hace para overridear el call conformal_advance_E 
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield1)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield1=0.0_RKIND
!!!         end if
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield2)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield2=0.0_RKIND
!!!         end if
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield3)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield3=0.0_RKIND
!!!         end if
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield4)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield4=0.0_RKIND
!!!         end if
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield5)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield5=0.0_RKIND
!!!         end if
!!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield6)) then
!!!             nodo%already_YEEadvanced_byconformal_changedtoP

ECfield6 = 0.0_RKIND;
    //         end if
    //      end do
    // #ifdef CompileWithOpenMP
    // !$OMP END PARALLEL DO
    // #endif
    //
      
    return;

    } // end subroutine AdvanceWiresE

    // !!!
    // !!!!!!!!!!!!!!!!
    // !!!!!!


    void AdvanceWiresH(SGGFDTDINFO_t& sgg, int timeinstant, int layoutnumber, const std::string& wiresflavor, bool simu_devia, bool stochastic, bool experimentalVideal, int wirethickness, double eps0, double mu0)
    {
                    
        double eps0_in = eps0;
        double mu0_in = mu0;
        int wirethickness_in = wirethickness;
        bool simu_devia_in = simu_devia;
        bool stochastic_in = stochastic;
        bool experimentalVideal_in = experimentalVideal;
        SGGFDTDINFO_t& sgg_in = sgg;

        int n, jmed, layoutnumber_in, iw1, is1, is2;

        int timeinstant_in = timeinstant;
        double Iplus, IMinus, Qplus, QMinus, timei;
        double Vincid, Iincid;
        CurrentSegments_t* Segmento, * Segmento2;
        ChargeNodes_t* Nodo;
        TMultiline_t* Multiline;
        const std::string& wiresflavor_in = wiresflavor;
        timei = sgg.tiempo(timeinstant); 
   
        //

        for (n = 1; n <= HWires.NumCurrentSegments; ++n) {
            Segmento = &HWires.CurrentSegment(n);
#ifdef CompileWithThickWires
            if (wirethickness > 1) {
                Advance_Thick_Hfield_wire2main(sgg, *Segmento, eps0, mu0);
            }
#endif                      
            if (wirethickness == 1) {
                continue;
            }
        }

      
        return;

    } // end subroutine AdvanceWiresH

    // !!!


    // !!!!!!!
    // !!!!!!!

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !Advancing charge and current routine
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    void AdvanceWiresEcrank(SGGFDTDINFO_t& sgg, int timeinstant, int layoutnumber, const std::string& wiresflavor, bool simu_devia, bool stochastic)
    {
        SGGFDTDINFO_t& sgg_in = sgg;
        bool simu_devia_in = simu_devia;
        bool stochastic_in = stochastic;

        // !!!

        int n, jmed, layoutnumber_in, iw1, is1, is2;

        int timeinstant_in = timeinstant;
        double Iplus, IMinus, IplusPast, IMinusPast, source, timei;
        double Vincid, Iincid;
        CurrentSegments_t* Segmento, * Segmento2;
        ChargeNodes_t* Nodo;
        TMultiline_t* Multiline;
        const std::string& wiresflavor_in = wiresflavor;
        std::vector<double> a(HWires.NumCurrentSegments), b(HWires.NumCurrentSegments), c(HWires.NumCurrentSegments), d(HWires.NumCurrentSegments), x(HWires.NumCurrentSegments);
      
        timei = sgg.tiempo(timeinstant); 
        // !!!
        Iplus = -1.0; IMinus = -1.0;
      
        // !!!! deprecado en pscale  110219 

      
        // !!!mal a 0624 !comentado en todos los sabores de wires
        // !!!#ifdef CompileWithOpenMP
        // !!!!$OMP PARALLEL do DEFAULT(SHARED)  private(Nodo)
        // !!!#endif
        // !!!      do n=1,HWires%NumChargeNodes
        // !!!         Nodo => HWires%ChargeNode(n)
        // !!!!!!!140220 pon a PEC los viejos notouch=already_YEEadvanced_byconformal_changedtoPECfield que tenga conectado un nodo
        // !!!         !debe ser lo primero que se hace para overridear el call conformal_advance_E 
        // !!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield1)) then
        // !!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield1=0.0_RKIND
        // !!!         end if
        // !!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield2)) then
        // !!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield2=0.0_RKIND
        // !!!         end if
        // !!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield3)) then
        // !!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield3=0.0_RKIND
        // !!!         end if
        // !!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield4)) then
        // !!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield4=0.0_RKIND
        // !!!         end if
        // !!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield5)) then
        // !!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield5=0.0_RKIND
        // !!!         end if
        // !!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield6)) then
        // !!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield6=0.0_RKIND
        // !!!         end if
        // !!!      end do
        // !!!#ifdef CompileWithOpenMP
        // !!!!$OMP END PARALLEL DO
        // !!!#endif
        // !

        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // !Inject the current at n+1 into the electring FDTD field previously calculated at n+1
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        for (n = 1; n <= HWires.NumCurrentSegments; ++n) {
            Segmento = &HWires.CurrentSegment(n);
            //171216quitado         Segmento%Efield_wire2main_past = Segmento%Efield_wire2main
            Segmento->Efield_wire2main = static_cast<double>(Segmento->Efield_wire2main) - Segmento->cte5 * Segmento->Current;
        }


        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // !CURRENT ADVANCING from n+1 to n+2
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // !

        for (n = 1; n <= HWires.NumCurrentSegments; ++n) {
            Segmento = &HWires.CurrentSegment(n);
            b[n - 1] = segmento->diag;
            c[n - 1] = segmento->upperdiag;
            a[n - 1] = segmento->lowerdiag;
            d[n - 1] =        segmento->rightCU       * Segmento->Current;
            d[n - 1] = d[n - 1] + segmento->rightCHminus  * Segmento->ChargeMinus->ChargePresent;
            d[n - 1] = d[n - 1] + segmento->rightCHplus   * Segmento->ChargePlus->ChargePresent;
            d[n - 1] = d[n - 1] + static_cast<double>(Segmento->Efield_wire2main);
            if (Segmento->ChargeMinus->NumCurrentMinus == 1) d[n - 1] = d[n - 1] + segmento->rightCUminus  * Segmento->ChargeMinus->CurrentMinus_1->Current;
            if (Segmento->ChargePlus->NumCurrentPlus == 1) d[n - 1] = d[n - 1] + segmento->rightCUplus   * Segmento->ChargePlus->CurrentPlus_1->Current;
            if (!simu_devia) {             
                If (Segmento->HasVsource) {
                    source = evolucion(timei, Segmento->Vsource->Fichero->Samples, &
                    Segmento->Vsource->Fichero->DeltaSamples, Segmento->Vsource->Fichero->NumSamples);
                    d[n - 1] = d[n - 1] + source;
                }
            }
            
        }
        n = HWires.NumCurrentSegments;
        solve_tridiag_wires(a, b, c, d, x, n);
        //  a - sub-diagonal (means it is the diagonal below the main diagonal)
        //  b - the main diagonal
        //  c - sup-diagonal (means it is the diagonal above the main diagonal)
        //  d - right part
        //  x - the answer
        //  n - number of equations

#ifdef CompileWithOpenMP
        #pragma omp parallel for default(shared) private(Segmento, jmed)
#endif
        for (n = 1; n <= HWires.NumCurrentSegments; ++n) {
            Segmento = &HWires.CurrentSegment(n);
            Segmento->CurrentPast = Segmento->Current; //save this data only for the probes observation right time positioning
            Segmento->Current = x[n - 1];
        }
#ifdef CompileWithOpenMP
        #pragma omp parallel for end
#endif



        //
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // !Next ADVANCE THE CHARGE from n-1.0_RKIND_wires / 2 to n+1/2 using the current known at n+1/2
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CompileWithOpenMP
        #pragma omp parallel for default(shared) private(Iplus, IMinus, Nodo)
#endif
        for (n = 1; n <= HWires.NumChargeNodes; ++n) {
            Nodo = &HWires.ChargeNode(n);
            if (nodo->exists) {
                Nodo->ChargePast = Nodo->ChargePresent;
                If (Nodo->NumCurrentPlus == 1) Iplus = Nodo->CurrentPlus_1->Current;
                If (Nodo->NumCurrentPlus == 1) IplusPast = Nodo->CurrentPlus_1->CurrentPast;
                //
                If (Nodo->NumCurrentMinus == 1) IMinus = Nodo->CurrentMinus_1->Current;
                If (Nodo->NumCurrentMinus == 1) IMinusPast = Nodo->CurrentMinus_1->CurrentPast;
                //
                if ((Nodo->NumCurrentMinus == 1) && (Nodo->NumCurrentPlus == 0)) Iplus = -IMinus;
                if ((Nodo->NumCurrentMinus == 0) && (Nodo->NumCurrentPlus == 1)) IMinus = -Iplus;
                if ((Nodo->NumCurrentMinus == 1) && (Nodo->NumCurrentPlus == 0)) IplusPast = -IMinusPast;
                if ((Nodo->NumCurrentMinus == 0) && (Nodo->NumCurrentPlus == 1)) IMinusPast = -IplusPast;

                Nodo->ChargePresent = Nodo->CteProp * Nodo->ChargePast - Nodo->CtePlain * ((Iplus + IplusPast) / 2.0_RKIND_wires - (IMinus + IminusPast) / 2.0_RKIND_wires);

            }
        }
#ifdef CompileWithOpenMP
        #pragma omp parallel for end
#endif

        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // !END OF CHARGE ADVANCING from n+1.0_RKIND_wires / 2 to n+3/2
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
        // !!!mal a 1024 !comentado en todos los sabores de wires
        // !!!#ifdef CompileWithOpenMP
        // !!!!$OMP PARALLEL do DEFAULT(SHARED)  private(Nodo)
        // !!!#endif
        // !!!      do n=1,HWires%NumChargeNodes
        // !!!         Nodo => HWires%ChargeNode(n)
        // !!!!!!!140220 pon a PEC los viejos notouch=already_YEEadvanced_byconformal_changedtoPECfield que tenga conectado un nodo
        // !!!         !debe ser lo primero que se hace para overridear el call conformal_advance_E 
        // !!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield1)) then
        // !!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield1=0.0_RKIND
        // !!!         end if
        // !!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield2)) then
        // !!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield2=0.0_RKIND
        // !!!         end if
        // !!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield3)) then
        // !!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield3=0.0_RKIND
        // !!!         end if
        // !!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield4)) then
        // !!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield4=0.0_RKIND
        // !!!         end if
        // !!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield5)) then
        // !!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield5=0.0_RKIND
        // !!!         end if
        // !!!         if (associated(nodo%already_YEEadvanced_byconformal_changedtoPECfield6)) then
        // !!!             nodo%already_YEEadvanced_byconformal_changedtoPECfield6=0.0_RKIND
        // !!!         end if
        // !!!      end do
        // !!!#ifdef CompileWithOpenMP
        // !!!!$OMP END PARALLEL DO
        // !!!#endif
        // !

        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // !END OF CURRENT ADVANCING from n+1 to n+2
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // !!!
        return;
    } // end subroutine AdvanceWiresEcrank


    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!! Function to interpolate the evolution files at the desired time
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    double evolucion(double t, const std::vector<double>& evol, double deltaevol, int numus)
    {
        int nprev;
        // integer(kind=8) :: nprev
        // real(kind=RKIND_wires) , dimension(0 : numus) :: evol
        // real(kind=RKIND_wires) :: deltaevol, t
        // !
        nprev = static_cast<int>((t) / deltaevol);
        if ((nprev + 1 > numus) || (NPREV + 1 <= 0)) { //SI NPREV<0 ES PORQUE SE HA DESBORADO EL ENTERO !BUG MIGEL 130614
            return 0.0_RKIND_wires; //if running out of samples, assume they are null
        }
        else {
            return (evol[nprev + 1] - evol[nprev]) / deltaevol * ((t) - nprev * deltaevol) + evol[nprev]; //linear interpolation
        }

        return 0.0; // Dummy return to satisfy compiler
    } // end function evolucion


    // !!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !Routine to store the fields for resumed problems
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!

    void StoreFieldsWires()
    {
        int i1;

        //store data for resuming
        for (i1 = 1; i1 <= HWires.NumChargeNodes; ++i1) {
            write(14, err = 634) HWires.ChargeNode(i1).ChargePresent;
            write(14, err = 634) HWires.ChargeNode(i1).ChargePast;
        }
        for (i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
            write(14, err = 634) HWires.CurrentSegment(i1).Current;
            write(14, err = 634) HWires.CurrentSegment(i1).qplus_qminus;
            write(14, err = 634) HWires.CurrentSegment(i1).current_for_devia;
            write(14, err = 634) HWires.CurrentSegment(i1).qplus_qminus_for_devia;
            write(14, err = 634) HWires.CurrentSegment(i1).Efield_main2wire_for_devia;
        }
        //
#ifdef CompileWithMPI
        write(14, err = 634) Hwires.NumNeededCurrentUpMPI, Hwires.NumNeededCurrentDownMPI;
        for (i1 = 1; i1 <= Hwires.NumNeededCurrentUpMPI; ++i1) {
            write(14, err = 634) HWires.MPIUpNeededCurrentSegment(i1).Current;
        }
        for (i1 = 1; i1 <= Hwires.NumNeededCurrentDownMPI; ++i1) {
            write(14, err = 634) HWires.MPIDownNeededCurrentSegment(i1).Current;
        }
#endif

        goto 635;
    634:   print11(0, SEPARADOR + separador + separador);
        print11(0, "WIRES: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING");
        print11(0, SEPARADOR + separador + separador);          
    635:   return;
    } // end subroutine StoreFieldsWires

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !Routine to free-up memory upon termination
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!

    void DestroyWires(SGGFDTDINFO_t& sgg)
    {
        SGGFDTDINFO_t& sgg_in = sgg;
        int i;

        //free up memory !ojo no se como hacerlo
        for (i = 1; i <= sgg.NumMedia; ++i) {
            if (sgg.Med(i).Is->ThinWire) {
                if (sgg.Med(i).wire(1).Vsource) {
                    delete sgg.Med(i).wire(1).Vsource;
                    sgg.Med(i).wire(1).Vsource = nullptr;
                }
                if (sgg.Med(i).wire(1).Isource) {
                    delete sgg.Med(i).wire(1).Isource;
                    sgg.Med(i).wire(1).Isource = nullptr;
                }
                if (sgg.Med(i).wire) {
                    delete[] sgg.Med(i).wire;
                    sgg.Med(i).wire = nullptr;
                }
            }
        }

        if (HWires.WireTipoMedio) {
            delete HWires.WireTipoMedio;
            HWires.WireTipoMedio = nullptr;
        }
        if (HWires.CurrentSegment) {
            delete[] HWires.CurrentSegment;
            HWires.CurrentSegment = nullptr;
        }
        if (HWires.ChargeNode) {
            delete[] HWires.ChargeNode;
            HWires.ChargeNode = nullptr;
        }
#ifdef CompileWithMPI
        if (Hwires.NumNeededCurrentUpMPI > 0) {
            delete[] HWires.MPIUpNeededCurrentSegment;
            HWires.MPIUpNeededCurrentSegment = nullptr;
        }
        if (Hwires.NumNeededCurrentDownMPI > 0) {
            delete[] HWires.MPIDownNeededCurrentSegment;
            HWires.MPIDownNeededCurrentSegment = nullptr;
        }
#endif
    } // end subroutine

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !Function to publish the private wire data (used in observation and in MPI)
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!

    Thinwires_t* GetHwires()
    {
        Thinwires_t* r;
        r = &Hwires;
        return r;
    } // end function GetHwires




    // !!!!!!!!!!!!!!!!!!!Just for reporting
    // !!! It also set the IsHeterogeneousJunction flag correctly (not needed for anything else than reporting)
    void ReportWireJunctions(int layoutnumber, int num_procs, bool therearewires, int ZI, int ZE, bool groundwires, bool strictOLD, bool verbose)
    {

        bool paralelos, groundwires_in, therearewires_in, Terminal, IsHeterogeneousJunction, paraErr, strictOLD_in, verbose_in;
        std::string buff;
        int i1, j1, layoutnumber_in, zi, ze, ierr, num_procs_in, indio;
        int mini = 1000000000, minj = 1000000000, mink = 1000000000, maxi = -1000000000, maxj = -1000000000, maxk = -1000000000;
        CurrentSegments_t* org, * fin;
        std::string DIR[3];
        std::string ig;
        ChargeNodes_t* nodo;
        struct nodosopentoair_t {
            int i, j, k, indexnode;
        };
      
        std::vector<nodosopentoair_t> nodosopentoair;


        DIR[iEx] = " X ";
        DIR[iEy] = " Y ";
        DIR[iEz] = " Z ";




#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, ierr);
#endif
        if (therearewires) {
            buff = "----------------------------------------------------------------";
            WarnErrReport(buff);
        }

        for (i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
            if (HWires.CurrentSegment(i1).i < mini) mini = HWires.CurrentSegment(i1).i;
            if (HWires.CurrentSegment(i1).j < minj) minj = HWires.CurrentSegment(i1).j;
            if (HWires.CurrentSegment(i1).k < mink) mink = HWires.CurrentSegment(i1).k;
            if (HWires.CurrentSegment(i1).i > maxi) maxi = HWires.CurrentSegment(i1).i;
            if (HWires.CurrentSegment(i1).j > maxj) maxj = HWires.CurrentSegment(i1).j;
            if (HWires.CurrentSegment(i1).k > maxk) maxk = HWires.CurrentSegment(i1).k;
            //
            if ((HWires.CurrentSegment(i1).ChargeMinus->NumCurrentPlus > MaxNumCurrentMinusPlus) ||
                (HWires.CurrentSegment(i1).ChargeMinus->NumCurrentMinus > MaxNumCurrentMinusPlus)) {
                buff = "wir3_BUGGYERROR: More than " + std::to_string(MaxNumCurrentMinusPlus) + " plus/minus junctions at " + std::to_string(HWires.CurrentSegment(i1).origindex) + "," + std::to_string(HWires.CurrentSegment(i1).i) + "," + std::to_string(HWires.CurrentSegment(i1).j) + "," + std::to_string(HWires.CurrentSegment(i1).k);
                WarnErrReport(buff, true);
                buff = "Contact with  to increase this limit.";
                WarnErrReport(buff, true);
            }
        }
        if ((maxi >= mini) && (maxj >= minj) && (maxk >= mink)) buff = "wir3_INFO: BBOX for Holland WIREs " + std::to_string(mini) + " " + std::to_string(minj) + " " + std::to_string(mink) + " " + std::to_string(maxi) + " " + std::to_string(maxj) + " " + std::to_string(maxk);
        if (HWires.NumCurrentSegments != 0) WarnErrReport(buff);


        for (i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
            org = &HWires.CurrentSegment(i1);
            for (j1 = i1 + 1; j1 <= HWires.NumCurrentSegments; ++j1) {
                fin = &HWires.CurrentSegment(j1);
                paralelos = (org->i == fin->i) && (org->j == fin->j) && (org->k == fin->k) && (org->tipofield == fin->tipofield);
                if (paralelos) {
                    if (org->indexmed != fin->indexmed) {
                        buff = "wir3_INFO: Parallel segments from different wires (multiWIRE) " + std::to_string(org->origindex) + " " + std::to_string(fin->origindex) + " at " + std::to_string(org->i) + " " + std::to_string(org->j) + " " + std::to_string(org->k);

                        if ((org->k >= ZI) && (org->k <= ZE) && verbose) WarnErrReport(buff);
                    }
                    else {
                        if (!strictOLD) {
                            buff = "wir3_BUGGYERROR: Parallel segments from the same WIRE. " + std::to_string(org->origindex) + " " + std::to_string(fin->origindex) + " at " + std::to_string(org->i) + " " + std::to_string(org->j) + " " + std::to_string(org->k);
                            WarnErrReport(buff, true);
                        }
                        else {
                            buff = "wir3_INFO: Parallel segments from the same wire (multiWIRE) " + std::to_string(org->origindex) + " " + std::to_string(fin->origindex) + " at " + std::to_string(org->i) + " " + std::to_string(org->j) + " " + std::to_string(org->k);

                            if ((org->k >= ZI) && (org->k <= ZE) && verbose) WarnErrReport(buff);
                        }
                    }
                }
            }
        }
        //

#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, ierr);
#endif

        for (i1 = 1; i1 <= HWires.NumChargeNodes; ++i1) {
            nodo = &HWires.Chargenode(i1);
            if (nodo->exists) {
                if ((nodo->NumCurrentMinus == 1) && (nodo->NumCurrentPlus == 1)) {
                    //special cases of 1plus and 1minus
                    IsHeterogeneousJunction = (nodo->CurrentMinus_1->indexmed != nodo->CurrentPlus_1->indexmed);
                }
                else if (nodo->NumCurrentMinus + nodo->NumCurrentPlus == 1) {
                    IsHeterogeneousJunction = false;
                }
                else {
                    IsHeterogeneousJunction = true;
                    if (nodo->NumCurrentMInus >= 2)
                        IsHeterogeneousJunction = IsHeterogeneousJunction &&
                        (nodo->CurrentMinus_1->indexmed != nodo->CurrentMinus_2->indexmed);
                    if (nodo->NumCurrentMInus >= 3)
                        IsHeterogeneousJunction = IsHeterogeneousJunction &&
                        (nodo->CurrentMinus_1->indexmed != nodo->CurrentMinus_3->indexmed) &&
                        (nodo->CurrentMinus_2->indexmed != nodo->CurrentMinus_3->indexmed);
                    if (nodo->NumCurrentMInus >= 4)
                        IsHeterogeneousJunction = IsHeterogeneousJunction &&
                        (nodo->CurrentMinus_1->indexmed != nodo->CurrentMinus_4->indexmed) &&
                        (nodo->CurrentMinus_2->indexmed != nodo->CurrentMinus_4->indexmed) &&
                        (nodo->CurrentMinus_3->indexmed != nodo->CurrentMinus_4->indexmed);
                    if (nodo->NumCurrentMInus >= 5)
                        IsHeterogeneousJunction = IsHeterogeneousJunction &&
                        (nodo->CurrentMinus_1->indexmed != nodo->CurrentMinus_5->indexmed) &&
                        (nodo->CurrentMinus_2->indexmed != nodo->CurrentMinus_5->indexmed) &&
                        (nodo->CurrentMinus_3->indexmed != nodo->CurrentMinus_5->indexmed) &&
                        (nodo->CurrentMinus_4->indexmed != nodo->CurrentMinus_5->indexmed);
                    if (nodo->NumCurrentMInus >= 6)
                        IsHeterogeneousJunction = IsHeterogeneousJunction &&
                        (nodo->CurrentMinus_1->indexmed != nodo->CurrentMinus_6->indexmed) &&
                        (nodo->CurrentMinus

} else if (
                (nodo.CurrentMinus_6.indexmed != nodo.CurrentMinus_1.indexmed) &&
                (nodo.CurrentMinus_6.indexmed != nodo.CurrentMinus_2.indexmed) &&
                (nodo.CurrentMinus_6.indexmed != nodo.CurrentMinus_3.indexmed) &&
                (nodo.CurrentMinus_6.indexmed != nodo.CurrentMinus_4.indexmed) &&
                (nodo.CurrentMinus_6.indexmed != nodo.CurrentMinus_5.indexmed)
            ) {
                if (nodo.NumCurrentMInus >= 7) {
                    IsHeterogeneousJunction = IsHeterogeneousJunction &&
                        (nodo.CurrentMinus_1.indexmed != nodo.CurrentMinus_7.indexmed) &&
                        (nodo.CurrentMinus_2.indexmed != nodo.CurrentMinus_7.indexmed) &&
                        (nodo.CurrentMinus_3.indexmed != nodo.CurrentMinus_7.indexmed) &&
                        (nodo.CurrentMinus_4.indexmed != nodo.CurrentMinus_7.indexmed) &&
                        (nodo.CurrentMinus_5.indexmed != nodo.CurrentMinus_7.indexmed) &&
                        (nodo.CurrentMinus_6.indexmed != nodo.CurrentMinus_7.indexmed);
                }
                if (nodo.NumCurrentMInus >= 8) {
                    IsHeterogeneousJunction = IsHeterogeneousJunction &&
                        (nodo.CurrentMinus_1.indexmed != nodo.CurrentMinus_8.indexmed) &&
                        (nodo.CurrentMinus_2.indexmed != nodo.CurrentMinus_8.indexmed) &&
                        (nodo.CurrentMinus_3.indexmed != nodo.CurrentMinus_8.indexmed) &&
                        (nodo.CurrentMinus_4.indexmed != nodo.CurrentMinus_8.indexmed) &&
                        (nodo.CurrentMinus_5.indexmed != nodo.CurrentMinus_8.indexmed) &&
                        (nodo.CurrentMinus_6.indexmed != nodo.CurrentMinus_8.indexmed) &&
                        (nodo.CurrentMinus_7.indexmed != nodo.CurrentMinus_8.indexmed);
                }
                if (nodo.NumCurrentMInus >= 9) {
                    IsHeterogeneousJunction = IsHeterogeneousJunction &&
                        (nodo.CurrentMinus_1.indexmed != nodo.CurrentMinus_9.indexmed) &&
                        (nodo.CurrentMinus_2.indexmed != nodo.CurrentMinus_9.indexmed) &&
                        (nodo.CurrentMinus_3.indexmed != nodo.CurrentMinus_9.indexmed) &&
                        (nodo.CurrentMinus_4.indexmed != nodo.CurrentMinus_9.indexmed) &&
                        (nodo.CurrentMinus_5.indexmed != nodo.CurrentMinus_9.indexmed) &&
                        (nodo.CurrentMinus_6.indexmed != nodo.CurrentMinus_9.indexmed) &&
                        (nodo.CurrentMinus_7.indexmed != nodo.CurrentMinus_9.indexmed) &&
                        (nodo.CurrentMinus_8.indexmed != nodo.CurrentMinus_9.indexmed);
                }
            }

            if (nodo.NumCurrentPlus >= 2) {
                IsHeterogeneousJunction = IsHeterogeneousJunction &&
                    (nodo.CurrentPlus_1.indexmed != nodo.CurrentPlus_2.indexmed);
            }
            if (nodo.NumCurrentPlus >= 3) {
                IsHeterogeneousJunction = IsHeterogeneousJunction &&
                    (nodo.CurrentPlus_1.indexmed != nodo.CurrentPlus_3.indexmed) &&
                    (nodo.CurrentPlus_2.indexmed != nodo.CurrentPlus_3.indexmed);
            }
            if (nodo.NumCurrentPlus >= 4) {
                IsHeterogeneousJunction = IsHeterogeneousJunction &&
                    (nodo.CurrentPlus_1.indexmed != nodo.CurrentPlus_4.indexmed) &&
                    (nodo.CurrentPlus_2.indexmed != nodo.CurrentPlus_4.indexmed) &&
                    (nodo.CurrentPlus_3.indexmed != nodo.CurrentPlus_4.indexmed);
            }
            if (nodo.NumCurrentPlus >= 5) {
                IsHeterogeneousJunction = IsHeterogeneousJunction &&
                    (nodo.CurrentPlus_1.indexmed != nodo.CurrentPlus_5.indexmed) &&
                    (nodo.CurrentPlus_2.indexmed != nodo.CurrentPlus_5.indexmed) &&
                    (nodo.CurrentPlus_3.indexmed != nodo.CurrentPlus_5.indexmed) &&
                    (nodo.CurrentPlus_4.indexmed != nodo.CurrentPlus_5.indexmed);
            }
            if (nodo.NumCurrentPlus >= 6) {
                IsHeterogeneousJunction = IsHeterogeneousJunction &&
                    (nodo.CurrentPlus_1.indexmed != nodo.CurrentPlus_6.indexmed) &&
                    (nodo.CurrentPlus_2.indexmed != nodo.CurrentPlus_6.indexmed) &&
                    (nodo.CurrentPlus_3.indexmed != nodo.CurrentPlus_6.indexmed) &&
                    (nodo.CurrentPlus_4.indexmed != nodo.CurrentPlus_6.indexmed) &&
                    (nodo.CurrentPlus_5.indexmed != nodo.CurrentPlus_6.indexmed);
            }
            if (nodo.NumCurrentPlus >= 7) {
                IsHeterogeneousJunction = IsHeterogeneousJunction &&
                    (nodo.CurrentPlus_1.indexmed != nodo.CurrentPlus_7.indexmed) &&
                    (nodo.CurrentPlus_2.indexmed != nodo.CurrentPlus_7.indexmed) &&
                    (nodo.CurrentPlus_3.indexmed != nodo.CurrentPlus_7.indexmed) &&
                    (nodo.CurrentPlus_4.indexmed != nodo.CurrentPlus_7.indexmed) &&
                    (nodo.CurrentPlus_5.indexmed != nodo.CurrentPlus_7.indexmed) &&
                    (nodo.CurrentPlus_6.indexmed != nodo.CurrentPlus_7.indexmed);
            }
            if (nodo.NumCurrentPlus >= 8) {
                IsHeterogeneousJunction = IsHeterogeneousJunction &&
                    (nodo.CurrentPlus_1.indexmed != nodo.CurrentPlus_8.indexmed) &&
                    (nodo.CurrentPlus_2.indexmed != nodo.CurrentPlus_8.indexmed) &&
                    (nodo.CurrentPlus_3.indexmed != nodo.CurrentPlus_8.indexmed) &&
                    (nodo.CurrentPlus_4.indexmed != nodo.CurrentPlus_8.indexmed) &&
                    (nodo.CurrentPlus_5.indexmed != nodo.CurrentPlus_8.indexmed) &&
                    (nodo.CurrentPlus_6.indexmed != nodo.CurrentPlus_8.indexmed) &&
                    (nodo.CurrentPlus_7.indexmed != nodo.CurrentPlus_8.indexmed);
            }
            if (nodo.NumCurrentPlus >= 9) {
                IsHeterogeneousJunction = IsHeterogeneousJunction &&
                    (nodo.CurrentPlus_1.indexmed != nodo.CurrentPlus_9.indexmed) &&
                    (nodo.CurrentPlus_2.indexmed != nodo.CurrentPlus_9.indexmed) &&
                    (nodo.CurrentPlus_3.indexmed != nodo.CurrentPlus_9.indexmed) &&
                    (nodo.CurrentPlus_4.indexmed != nodo.CurrentPlus_9.indexmed) &&
                    (nodo.CurrentPlus_5.indexmed != nodo.CurrentPlus_9.indexmed) &&
                    (nodo.CurrentPlus_6.indexmed != nodo.CurrentPlus_9.indexmed) &&
                    (nodo.CurrentPlus_7.indexmed != nodo.CurrentPlus_9.indexmed) &&
                    (nodo.CurrentPlus_8.indexmed != nodo.CurrentPlus_9.indexmed);
            }
        }

        if ((IsHeterogeneousJunction && nodo.IsHeterogeneousJunction) ||
            (!IsHeterogeneousJunction && !nodo.IsHeterogeneousJunction)) {
            continue;
        } else {
            std::ostringstream buff;
            buff << "wir3_BUGGYERROR: Heterogeneous Junctions mismatch "
                 << std::setw(7) << std::setfill(' ') << nodo.indexnode
                 << std::setw(7) << std::setfill(' ') << nodo.i
                 << std::setw(7) << std::setfill(' ') << nodo.j
                 << std::setw(7) << std::setfill(' ') << nodo.k
                 << " ("
                 << std::setw(2) << std::setfill(' ') << nodo.numcurrentminus
                 << std::setw(2) << std::setfill(' ') << nodo.numcurrentplus
                 << "). ";
            if ((nodo.k >= ZI) && (nodo.k <= ZE)) {
                WarnErrReport(buff.str(), true);
            }
        }
    } // end if exist

} // do i1

#ifdef CompileWithMPI
    MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif

// !!!!!!!!!!!!!!!!
for (int i1 = 1; i1 <= HWires.NumChargeNodes; ++i1) {
    nodo = Hwires.ChargeNode(i1);
    if (nodo.exists) {
        if (nodo.isLossy || nodo.isPEC) {
            std::string ig = "";
            if (nodo.isLossy) ig = " to Lossy";
            if (nodo.isPEC) ig = " to PEC";

            if ((nodo.Is_LeftEnd) || (nodo.Is_RightEnd)) {
                std::ostringstream buff;
                buff << "wir3_INFO: Terminal (LeftEnd/RightEnd) node directly GROUNDED  "
                     << std::setw(7) << std::setfill(' ') << nodo.indexnode
                     << std::setw(7) << std::setfill(' ') << nodo.i
                     << std::setw(7) << std::setfill(' ') << nodo.j
                     << std::setw(7) << std::setfill(' ') << nodo.k
                     << " ("
                     << std::setw(2) << std::setfill(' ') << nodo.numcurrentminus
                     << std::setw(2) << std::setfill(' ') << nodo.numcurrentplus
                     << ")" << ig
                     << std::setw(12) << std::setprecision(2) << std::scientific << Nodo.CteProp
                     << std::setw(12) << std::setprecision(2) << std::scientific << Nodo.CtePlain;
                if ((nodo.k >= ZI) && (nodo.k <= ZE) && verbose) {
                    WarnErrReport(buff.str());
                }
            } else if ((nodo.NumCurrentPlus + nodo.NumCurrentMinus < 2)) {
                std::ostringstream buff;
                buff << "wir3_WARNING: Terminal (other) node  direcly GROUNDED  "
                     << std::setw(7) << std::setfill(' ') << nodo.indexnode
                     << std::setw(7) << std::setfill(' ') << nodo.i
                     << std::setw(7) << std::setfill(' ') << nodo.j
                     << std::setw(7) << std::setfill(' ') << nodo.k
                     << " ("
                     << std::setw(2) << std::setfill(' ') << nodo.numcurrentminus
                     << std::setw(2) << std::setfill(' ') << nodo.numcurrentplus
                     << ")" << ig
                     << std::setw(12) << std::setprecision(2) << std::scientific << Nodo.CteProp
                     << std::setw(12) << std::setprecision(2) << std::scientific << Nodo.CtePlain;
                if ((nodo.k >= ZI) && (nodo.k <= ZE)) {
                    WarnErrReport(buff.str());
                }
            } else {
                if (groundwires) {
                    std::ostringstream buff;
                    buff << "wir3_INFO: NON-terminal node directly GROUNDED (-groundWIREs)"
                         << std::setw(7) << std::setfill(' ') << nodo.indexnode
                         << std::setw(7) << std::setfill(' ') << nodo.i
                         << std::setw(7) << std::setfill(' ') << nodo.j
                         << std::setw(7) << std::setfill(' ') << nodo.k
                         << " ("
                         << std::setw(2) << std::setfill(' ') << nodo.numcurrentminus
                         << std::setw(2) << std::setfill(' ') << nodo.numcurrentplus
                         << ")" << ig
                         << std::setw(12) << std::setprecision(2) << std::scientific << Nodo.CteProp
                         << std::setw(12) << std::setprecision(2) << std::scientific << Nodo.CtePlain;
                    if ((nodo.k >= ZI) && (nodo.k <= ZE) && verbose) {
                        WarnErrReport(buff.str());
                    }
                } else {
                    std::ostringstream buff;
                    buff << "wir3_WARNING: NON-terminal node directly GROUNDED "
                         << std::setw(7) << std::setfill(' ') << nodo.indexnode
                         << std::setw(7) << std::setfill(' ') << nodo.i
                         << std::setw(7) << std::setfill(' ') << nodo.j
                         << std::setw(7) << std::setfill(' ') << nodo.k
                         << " ("
                         << std::setw(2) << std::setfill(' ') << nodo.numcurrentminus
                         << std::setw(2) << std::setfill(' ') << nodo.numcurrentplus
                         << ")" << ig
                         << std::setw(12) << std::setprecision(2) << std::scientific << Nodo.CteProp
                         << std::setw(12) << std::setprecision(2) << std::scientific << Nodo.CtePlain;
                    if ((nodo.k >= ZI) && (nodo.k <= ZE)) {
                        WarnErrReport(buff.str());
                    }
                }
            }
        }
    }
}

for (int i1 = 1; i1 <= HWires.NumChargeNodes; ++i1) {
    nodo = Hwires.ChargeNode(i1);
    if (nodo.exists) {
        if (nodo.IsHeterogeneousJunction) {
            std::ostringstream buff;
            buff << "wir3_INFO: MultiWIRE Junction made at node "
                 << std::setw(7) << std::setfill(' ') << nodo.indexnode
                 << std::setw(7) << std::setfill(' ') << nodo.i
                 << std::setw(7) << std::setfill(' ') << nodo.j
                 << std::setw(7) << std::setfill(' ') << nodo.k
                 << " ("
                 << std::setw(2) << std::setfill(' ') << nodo.numcurrentminus
                 << std::setw(2) << std::setfill(' ') << nodo.numcurrentplus
                 << ")";
            if ((nodo.k >= ZI) && (nodo.k <= ZE) && verbose) {
                WarnErrReport(buff.str());
            }
            bool Terminal = false; // 7/2/14 esti estaba mal

            if (nodo.NumCurrentMinus >= 1) {
                Terminal = Terminal ||
                    (nodo.CurrentMinus_1.Is_LeftEnd || nodo.CurrentMinus_1.Is_RightEnd || nodo.CurrentMinus_1.IsEnd_norLeft_norRight);
                std::ostringstream buff_seg;
                buff_seg << "         Segment -1: "
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_1.origindex
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_1.i
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_1.j
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_1.k
                         << std::setw(7) << std::setfill(' ') << dir(nodo.CurrentMinus_1.tipofield)
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_1.indexmed;
                if ((nodo.k >= ZI) && (nodo.k <= ZE) && verbose) {
                    WarnErrReport(buff_seg.str());
                }
            }
            if (nodo.NumCurrentMinus >= 2) {
                Terminal = Terminal ||
                    (nodo.CurrentMinus_2.Is_LeftEnd || nodo.CurrentMinus_2.Is_RightEnd || nodo.CurrentMinus_2.IsEnd_norLeft_norRight);
                std::ostringstream buff_seg;
                buff_seg << "         Segment -2: "
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_2.origindex
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_2.i
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_2.j
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_2.k
                         << std::setw(7) << std::setfill(' ') << dir(nodo.CurrentMinus_2.tipofield)
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_2.indexmed;
                if ((nodo.k >= ZI) && (nodo.k <= ZE) && verbose) {
                    WarnErrReport(buff_seg.str());
                }
            }
            if (nodo.NumCurrentMinus >= 3) {
                Terminal = Terminal ||
                    (nodo.CurrentMinus_3.Is_LeftEnd || nodo.CurrentMinus_3.Is_RightEnd || nodo.CurrentMinus_3.IsEnd_norLeft_norRight);
                std::ostringstream buff_seg;
                buff_seg << "         Segment -3: "
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_3.origindex
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_3.i
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_3.j
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_3.k
                         << std::setw(7) << std::setfill(' ') << dir(nodo.CurrentMinus_3.tipofield)
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_3.indexmed;
                if ((nodo.k >= ZI) && (nodo.k <= ZE) && verbose) {
                    WarnErrReport(buff_seg.str());
                }
            }
            if (nodo.NumCurrentMinus >= 4) {
                Terminal = Terminal ||
                    (nodo.CurrentMinus_4.Is_LeftEnd || nodo.CurrentMinus_4.Is_RightEnd || nodo.CurrentMinus_4.IsEnd_norLeft_norRight);
                std::ostringstream buff_seg;
                buff_seg << "         Segment -4: "
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_4.origindex
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_4.i
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_4.j
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_4.k
                         << std::setw(7) << std::setfill(' ') << dir(nodo.CurrentMinus_4.tipofield)
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_4.indexmed;
                if ((nodo.k >= ZI) && (nodo.k <= ZE) && verbose) {
                    WarnErrReport(buff_seg.str());
                }
            }
            if (nodo.NumCurrentMinus >= 5) {
                Terminal = Terminal ||
                    (nodo.CurrentMinus_5.Is_LeftEnd || nodo.CurrentMinus_5.Is_RightEnd || nodo.CurrentMinus_5.IsEnd_norLeft_norRight);
                std::ostringstream buff_seg;
                buff_seg << "         Segment -5: "
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_5.origindex
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_5.i
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_5.j
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_5.k
                         << std::setw(7) << std::setfill(' ') << dir(nodo.CurrentMinus_5.tipofield)
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_5.indexmed;
                if ((nodo.k >= ZI) && (nodo.k <= ZE) && verbose) {
                    WarnErrReport(buff_seg.str());
                }
            }
            if (nodo.NumCurrentMinus >= 6) {
                Terminal = Terminal ||
                    (nodo.CurrentMinus_6.Is_LeftEnd || nodo.CurrentMinus_6.Is_RightEnd || nodo.CurrentMinus_6.IsEnd_norLeft_norRight);
                std::ostringstream buff_seg;
                buff_seg << "         Segment -6: "
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_6.origindex
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_6.i
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_6.j
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_6.k
                         << std::setw(7) << std::setfill(' ') << dir(nodo.CurrentMinus_6.tipofield)
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_6.indexmed;
                if ((nodo.k >= ZI) && (nodo.k <= ZE) && verbose) {
                    WarnErrReport(buff_seg.str());
                }
            }
            if (nodo.NumCurrentMinus >= 7) {
                Terminal = Terminal ||
                    (nodo.CurrentMinus_7.Is_LeftEnd || nodo.CurrentMinus_7.Is_RightEnd || nodo.CurrentMinus_7.IsEnd_norLeft_norRight);
                std::ostringstream buff_seg;
                buff_seg << "         Segment -7: "
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_7.origindex
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_7.i
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_7.j
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_7.k
                         << std::setw(7) << std::setfill(' ') << dir(nodo.CurrentMinus_7.tipofield)
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_7.indexmed;
                if ((nodo.k >= ZI) && (nodo.k <= ZE) && verbose) {
                    WarnErrReport(buff_seg.str());
                }
            }
            if (nodo.NumCurrentMinus >= 8) {
                Terminal = Terminal ||
                    (nodo.CurrentMinus_8.Is_LeftEnd || nodo.CurrentMinus_8.Is_RightEnd || nodo.CurrentMinus_8.IsEnd_norLeft_norRight);
                std::ostringstream buff_seg;
                buff_seg << "         Segment -8: "
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_8.origindex
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_8.i
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_8.j
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_8.k
                         << std::setw(7) << std::setfill(' ') << dir(nodo.CurrentMinus_8.tipofield)
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_8.indexmed;
                if ((nodo.k >= ZI) && (nodo.k <= ZE) && verbose) {
                    WarnErrReport(buff_seg.str());
                }
            }
            if (nodo.NumCurrentMinus >= 9) {
                Terminal = Terminal ||
                    (nodo.CurrentMinus_9.Is_LeftEnd || nodo.CurrentMinus_9.Is_RightEnd || nodo.CurrentMinus_9.IsEnd_norLeft_norRight);
                std::ostringstream buff_seg;
                buff_seg << "         Segment -9: "
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_9.origindex
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_9.i
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_9.j
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_9.k
                         << std::setw(7) << std::setfill(' ') << dir(nodo.CurrentMinus_9.tipofield)
                         << std::setw(7) << std::setfill(' ') << nodo.CurrentMinus_9.indexmed;
                if ((nodo.k >= ZI) && (nodo.k <= ZE)) {
                    WarnErrReport(buff_seg.str());
                }
            }

            if (nodo.NumCurrentplus >= 1) {
                Terminal = Terminal ||
                    (nodo.CurrentPlus_1.Is_LeftEnd || nodo.CurrentPlus_1.Is_RightEnd || nodo.CurrentPlus_1.IsEnd_norLeft_norRight);
                std::ostringstream buff_seg;
                buff_seg << "         Segment +1: "
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_1.origindex
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_1.i
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_1.j
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_1.k
                         << std::setw(7) << std::setfill(' ') << dir(nodo.Currentplus_1.tipofield)
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_1.indexmed;
                if ((nodo.k >= ZI) && (nodo.k <= ZE) && verbose) {
                    WarnErrReport(buff_seg.str());
                }
            }
            if (nodo.NumCurrentplus >= 2) {
                Terminal = Terminal ||
                    (nodo.CurrentPlus_2.Is_LeftEnd || nodo.CurrentPlus_2.Is_RightEnd || nodo.CurrentPlus_2.IsEnd_norLeft_norRight);
                std::ostringstream buff_seg;
                buff_seg << "         Segment +2: "
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_2.origindex
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_2.i
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_2.j
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_2.k
                         << std::setw(7) << std::setfill(' ') << dir(nodo.Currentplus_2.tipofield)
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_2.indexmed;
                if ((nodo.k >= ZI) && (nodo.k <= ZE) && verbose) {
                    WarnErrReport(buff_seg.str());
                }
            }
            if (nodo.NumCurrentplus >= 3) {
                Terminal = Terminal ||
                    (nodo.CurrentPlus_3.Is_LeftEnd || nodo.CurrentPlus_3.Is_RightEnd || nodo.CurrentPlus_3.IsEnd_norLeft_norRight);
                std::ostringstream buff_seg;
                buff_seg << "         Segment +3: "
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_3.origindex
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_3.i
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_3.j
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_3.k
                         << std::setw(7) << std::setfill(' ') << dir(nodo.Currentplus_3.tipofield)
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_3.indexmed;
                if ((nodo.k >= ZI) && (nodo.k <= ZE)) {
                    WarnErrReport(buff_seg.str());
                }
            }
            if (nodo.NumCurrentplus >= 4) {
                Terminal = Terminal ||
                    (nodo.CurrentPlus_4.Is_LeftEnd || nodo.CurrentPlus_4.Is_RightEnd || nodo.CurrentPlus_4.IsEnd_norLeft_norRight);
                std::ostringstream buff_seg;
                buff_seg << "         Segment +4: "
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_4.origindex
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_4.i
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_4.j
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_4.k
                         << std::setw(7) << std::setfill(' ') << dir(nodo.Currentplus_4.tipofield)
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_4.indexmed;
                if ((nodo.k >= ZI) && (nodo.k <= ZE) && verbose) {
                    WarnErrReport(buff_seg.str());
                }
            }
            if (nodo.NumCurrentplus >= 5) {
                Terminal = Terminal ||
                    (nodo.CurrentPlus_5.Is_LeftEnd || nodo.CurrentPlus_5.Is_RightEnd || nodo.CurrentPlus_5.IsEnd_norLeft_norRight);
                std::ostringstream buff_seg;
                buff_seg << "         Segment +5: "
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_5.origindex
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_5.i
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_5.j
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_5.k
                         << std::setw(7) << std::setfill(' ') << dir(nodo.Currentplus_5.tipofield)
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_5.indexmed;
                if ((nodo.k >= ZI) && (nodo.k <= ZE) && verbose) {
                    WarnErrReport(buff_seg.str());
                }
            }
            if (nodo.NumCurrentplus >= 6) {
                Terminal = Terminal ||
                    (nodo.CurrentPlus_6.Is_LeftEnd || nodo.CurrentPlus_6.Is_RightEnd || nodo.CurrentPlus_6.IsEnd_norLeft_norRight);
                std::ostringstream buff_seg;
                buff_seg << "         Segment +6: "
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_6.origindex
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_6.i
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_6.j
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_6.k
                         << std::setw(7) << std::setfill(' ') << dir(nodo.Currentplus_6.tipofield)
                         << std::setw(7) << std::setfill(' ') << nodo.Currentplus_6.indexmed;
                if ((nodo.k >= ZI) && (nodo.k <= ZE) && verbose) {
                    WarnErrReport(buff_seg.str());
                }
            }
        }
    }
}

if (nodo->NumCurrentplus >= 6) {
                    Terminal = Terminal ||
                               (nodo->CurrentPlus_6->Is_LeftEnd || nodo->CurrentPlus_6->Is_RightEnd || nodo->CurrentPlus_6->IsEnd_norLeft_norRight);
                    sprintf(buff, "         Segment +6: %7i%7i%7i%7i%s%7i",
                            "         Segment +6: ",
                            nodo->Currentplus_6->origindex,
                            nodo->Currentplus_6->i,
                            nodo->Currentplus_6->j,
                            nodo->Currentplus_6->k,
                            dir(nodo->Currentplus_6->tipofield),
                            nodo->Currentplus_6->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE) && verbose) WarnErrReport(buff);
                }
                if (nodo->NumCurrentplus >= 7) {
                    Terminal = Terminal ||
                               (nodo->CurrentPlus_7->Is_LeftEnd || nodo->CurrentPlus_7->Is_RightEnd || nodo->CurrentPlus_7->IsEnd_norLeft_norRight);
                    sprintf(buff, "         Segment +7: %7i%7i%7i%7i%s%7i",
                            "         Segment +7: ",
                            nodo->Currentplus_7->origindex,
                            nodo->Currentplus_7->i,
                            nodo->Currentplus_7->j,
                            nodo->Currentplus_7->k,
                            dir(nodo->Currentplus_7->tipofield),
                            nodo->Currentplus_7->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE) && verbose) WarnErrReport(buff);
                }
                if (nodo->NumCurrentplus >= 8) {
                    Terminal = Terminal ||
                               (nodo->CurrentPlus_8->Is_LeftEnd || nodo->CurrentPlus_8->Is_RightEnd || nodo->CurrentPlus_8->IsEnd_norLeft_norRight);
                    sprintf(buff, "         Segment +8: %7i%7i%7i%7i%s%7i",
                            "         Segment +8: ",
                            nodo->Currentplus_8->origindex,
                            nodo->Currentplus_8->i,
                            nodo->Currentplus_8->j,
                            nodo->Currentplus_8->k,
                            dir(nodo->Currentplus_8->tipofield),
                            nodo->Currentplus_8->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE) && verbose) WarnErrReport(buff);
                }
                if (nodo->NumCurrentplus >= 9) {
                    Terminal = Terminal ||
                               (nodo->CurrentPlus_9->Is_LeftEnd || nodo->CurrentPlus_9->Is_RightEnd || nodo->CurrentPlus_9->IsEnd_norLeft_norRight);
                    sprintf(buff, "         Segment +9: %7i%7i%7i%7i%s%7i",
                            "         Segment +9: ",
                            nodo->Currentplus_9->origindex,
                            nodo->Currentplus_9->i,
                            nodo->Currentplus_9->j,
                            nodo->Currentplus_9->k,
                            dir(nodo->Currentplus_9->tipofield),
                            nodo->Currentplus_9->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE) && verbose) WarnErrReport(buff);
                }
                if (!Terminal) {
                    sprintf(buff, "wir3_BUGGYERROR: Some of the segments are not terminal");
                    WarnErrReport(buff, true);
                }
            }
        }
    }

    // !!!!!!!!!!!!!!!!

    for (i1 = 1; i1 <= HWires->NumChargeNodes; i1++) {
        nodo = Hwires->ChargeNode(i1);
        if (nodo->exists) {
            if ((!nodo->IsHeterogeneousJunction) &&
                ((nodo->NumCurrentMinus + nodo->NumCurrentPlus) >= 3)) { // homogeneo y verdadera union
                if (!strictOLD) {
                    sprintf(buff, "wir3_INFO: Intra-WIRE Junction made between: ");
                    paraerr = false;
                    if ((nodo->k >= ZI) && (nodo->k <= ZE) && verbose) WarnErrReport(buff, paraerr);
                } else {
                    sprintf(buff, "wir3_ERROR: Intra-WIRE Junction of more than 2 segments is forbidden: ");
                    paraerr = true;
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
                sprintf(buff, "         At node: %7i%7i%7i",
                        "         At node: ",
                        nodo->i,
                        nodo->j,
                        nodo->k);
                if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                if (nodo->NumCurrentMinus >= 1) {
                    sprintf(buff, "         Segment -1: %7i%7i%7i%7i%s%7i",
                            "         Segment -1: ",
                            nodo->CurrentMinus_1->origindex,
                            nodo->CurrentMinus_1->i,
                            nodo->CurrentMinus_1->j,
                            nodo->CurrentMinus_1->k,
                            dir(nodo->CurrentMinus_1->tipofield),
                            nodo->CurrentMinus_1->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
                if (nodo->NumCurrentMinus >= 2) {
                    sprintf(buff, "         Segment -2: %7i%7i%7i%7i%s%7i",
                            "         Segment -2: ",
                            nodo->CurrentMinus_2->origindex,
                            nodo->CurrentMinus_2->i,
                            nodo->CurrentMinus_2->j,
                            nodo->CurrentMinus_2->k,
                            dir(nodo->CurrentMinus_2->tipofield),
                            nodo->CurrentMinus_2->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
                if (nodo->NumCurrentMinus >= 3) {
                    sprintf(buff, "         Segment -3: %7i%7i%7i%7i%s%7i",
                            "         Segment -3: ",
                            nodo->CurrentMinus_3->origindex,
                            nodo->CurrentMinus_3->i,
                            nodo->CurrentMinus_3->j,
                            nodo->CurrentMinus_3->k,
                            dir(nodo->CurrentMinus_3->tipofield),
                            nodo->CurrentMinus_3->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
                if (nodo->NumCurrentMinus >= 4) {
                    sprintf(buff, "         Segment -4: %7i%7i%7i%7i%s%7i",
                            "         Segment -4: ",
                            nodo->CurrentMinus_4->origindex,
                            nodo->CurrentMinus_4->i,
                            nodo->CurrentMinus_4->j,
                            nodo->CurrentMinus_4->k,
                            dir(nodo->CurrentMinus_4->tipofield),
                            nodo->CurrentMinus_4->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
                if (nodo->NumCurrentMinus >= 5) {
                    sprintf(buff, "         Segment -5: %7i%7i%7i%7i%s%7i",
                            "         Segment -5: ",
                            nodo->CurrentMinus_5->origindex,
                            nodo->CurrentMinus_5->i,
                            nodo->CurrentMinus_5->j,
                            nodo->CurrentMinus_5->k,
                            dir(nodo->CurrentMinus_5->tipofield),
                            nodo->CurrentMinus_5->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
                if (nodo->NumCurrentMinus >= 6) {
                    sprintf(buff, "         Segment -6: %7i%7i%7i%7i%s%7i",
                            "         Segment -6: ",
                            nodo->CurrentMinus_6->origindex,
                            nodo->CurrentMinus_6->i,
                            nodo->CurrentMinus_6->j,
                            nodo->CurrentMinus_6->k,
                            dir(nodo->CurrentMinus_6->tipofield),
                            nodo->CurrentMinus_6->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
                if (nodo->NumCurrentMinus >= 7) {
                    sprintf(buff, "         Segment -7: %7i%7i%7i%7i%s%7i",
                            "         Segment -7: ",
                            nodo->CurrentMinus_7->origindex,
                            nodo->CurrentMinus_7->i,
                            nodo->CurrentMinus_7->j,
                            nodo->CurrentMinus_7->k,
                            dir(nodo->CurrentMinus_7->tipofield),
                            nodo->CurrentMinus_7->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
                if (nodo->NumCurrentMinus >= 8) {
                    sprintf(buff, "         Segment -8: %7i%7i%7i%7i%s%7i",
                            "         Segment -8: ",
                            nodo->CurrentMinus_8->origindex,
                            nodo->CurrentMinus_8->i,
                            nodo->CurrentMinus_8->j,
                            nodo->CurrentMinus_8->k,
                            dir(nodo->CurrentMinus_8->tipofield),
                            nodo->CurrentMinus_8->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
                if (nodo->NumCurrentMinus >= 9) {
                    sprintf(buff, "         Segment -9: %7i%7i%7i%7i%s%7i",
                            "         Segment -9: ",
                            nodo->CurrentMinus_9->origindex,
                            nodo->CurrentMinus_9->i,
                            nodo->CurrentMinus_9->j,
                            nodo->CurrentMinus_9->k,
                            dir(nodo->CurrentMinus_9->tipofield),
                            nodo->CurrentMinus_9->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
                //
                if (nodo->NumCurrentplus >= 1) {
                    sprintf(buff, "         Segment +1: %7i%7i%7i%7i%s%7i",
                            "         Segment +1: ",
                            nodo->Currentplus_1->origindex,
                            nodo->Currentplus_1->i,
                            nodo->Currentplus_1->j,
                            nodo->Currentplus_1->k,
                            dir(nodo->Currentplus_1->tipofield),
                            nodo->Currentplus_1->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
                if (nodo->NumCurrentplus >= 2) {
                    sprintf(buff, "         Segment +2: %7i%7i%7i%7i%s%7i",
                            "         Segment +2: ",
                            nodo->Currentplus_2->origindex,
                            nodo->Currentplus_2->i,
                            nodo->Currentplus_2->j,
                            nodo->Currentplus_2->k,
                            dir(nodo->Currentplus_2->tipofield),
                            nodo->Currentplus_2->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
                if (nodo->NumCurrentplus >= 3) {
                    sprintf(buff, "         Segment +3: %7i%7i%7i%7i%s%7i",
                            "         Segment +3: ",
                            nodo->Currentplus_3->origindex,
                            nodo->Currentplus_3->i,
                            nodo->Currentplus_3->j,
                            nodo->Currentplus_3->k,
                            dir(nodo->Currentplus_3->tipofield),
                            nodo->Currentplus_3->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
                if (nodo->NumCurrentplus >= 4) {
                    sprintf(buff, "         Segment +4: %7i%7i%7i%7i%s%7i",
                            "         Segment +4: ",
                            nodo->Currentplus_4->origindex,
                            nodo->Currentplus_4->i,
                            nodo->Currentplus_4->j,
                            nodo->Currentplus_4->k,
                            dir(nodo->Currentplus_4->tipofield),
                            nodo->Currentplus_4->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
                if (nodo->NumCurrentplus >= 5) {
                    sprintf(buff, "         Segment +5: %7i%7i%7i%7i%s%7i",
                            "         Segment +5: ",
                            nodo->Currentplus_5->origindex,
                            nodo->Currentplus_5->i,
                            nodo->Currentplus_5->j,
                            nodo->Currentplus_5->k,
                            dir(nodo->Currentplus_5->tipofield),
                            nodo->Currentplus_5->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
                if (nodo->NumCurrentplus >= 6) {
                    sprintf(buff, "         Segment +6: %7i%7i%7i%7i%s%7i",
                            "         Segment +6: ",
                            nodo->Currentplus_6->origindex,
                            nodo->Currentplus_6->i,
                            nodo->Currentplus_6->j,
                            nodo->Currentplus_6->k,
                            dir(nodo->Currentplus_6->tipofield),
                            nodo->Currentplus_6->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
                if (nodo->NumCurrentplus >= 7) {
                    sprintf(buff, "         Segment +7: %7i%7i%7i%7i%s%7i",
                            "         Segment +7: ",
                            nodo->Currentplus_7->origindex,
                            nodo->Currentplus_7->i,
                            nodo->Currentplus_7->j,
                            nodo->Currentplus_7->k,
                            dir(nodo->Currentplus_7->tipofield),
                            nodo->Currentplus_7->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
                if (nodo->NumCurrentplus >= 8) {
                    sprintf(buff, "         Segment +8: %7i%7i%7i%7i%s%7i",
                            "         Segment +8: ",
                            nodo->Currentplus_8->origindex,
                            nodo->Currentplus_8->i,
                            nodo->Currentplus_8->j,
                            nodo->Currentplus_8->k,
                            dir(nodo->Currentplus_8->tipofield),
                            nodo->Currentplus_8->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
                if (nodo->NumCurrentplus >= 9) {
                    sprintf(buff, "         Segment +9: %7i%7i%7i%7i%s%7i",
                            "         Segment +9: ",
                            nodo->Currentplus_9->origindex,
                            nodo->Currentplus_9->i,
                            nodo->Currentplus_9->j,
                            nodo->Currentplus_9->k,
                            dir(nodo->Currentplus_9->tipofield),
                            nodo->Currentplus_9->indexmed);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff, paraerr);
                }
            }
        }
    }

#ifdef CompileWithMPI
    MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif

    // !!!!!!!!!!!!!!!!

    nodosopentoair = new NodeData[HWires->NumChargeNodes + 1];
    indio = 0;
    for (i1 = 1; i1 <= HWires->NumChargeNodes; i1++) {
        nodo = Hwires->ChargeNode(i1);
        if (nodo->exists) {
            if ((nodo->numcurrentPlus + nodo->numcurrentMinus) < 2) {
                if (!(nodo->IsPec || nodo->IsLossy)) {
                    indio = indio + 1;
                    nodosopentoair[indio].indexnode = nodo->indexnode;
                    nodosopentoair[indio].i = nodo->i;
                    nodosopentoair[indio].j = nodo->j;
                    nodosopentoair[indio].k = nodo->k;
                    sprintf(buff, "wir3_WARNING: Node open to air %7i%7i%7i%7i",
                            "wir3_WARNING: Node open to air ",
                            nodo->indexnode,
                            nodo->i,
                            nodo->j,
                            nodo->k);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE)) WarnErrReport(buff);
                } else {
                    sprintf(buff, "wir3_INFO: NON-JUNCTION Node  GROUNDED %7i%7i%7i%7i",
                            "wir3_INFO: NON-JUNCTION Node  GROUNDED ",
                            nodo->indexnode,
                            nodo->i,
                            nodo->j,
                            nodo->k);
                    if ((nodo->k >= ZI) && (nodo->k <= ZE) && verbose) WarnErrReport(buff);
                }
            }
        }
    }

    for (i1 = 1; i1 <= indio; i1++) {
        for (j1 = i1 + 1; j1 <= indio; j1++) {
            if ((nodosopentoair[i1].i == nodosopentoair[j1].i) &&
                (nodosopentoair[i1].j == nodosopentoair[j1].j) &&
                (nodosopentoair[i1].k == nodosopentoair[j1].k)) {
                sprintf(buff, "wir3_ERROR: TWO nodes at the same location open to air. Should be connected? %7i%7i%7i%7i%7i%7i%7i%7i%7i",
                        "wir3_ERROR: TWO nodes at the same location open to air. Should be connected? ",
                        indio, i1, j1,
                        nodosopentoair[i1].indexnode,
                        nodosopentoair[j1].indexnode,
                        nodosopentoair[i1].i,
                        nodosopentoair[i1].j,
                        nodosopentoair[i1].k);
                if ((nodosopentoair[i1].k >= ZI) && (nodosopentoair[i1].k <= ZE)) WarnErrReport(buff, true);
            }
        }
    }
    delete[] nodosopentoair;

#ifdef CompileWithMPI
    MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
    // more testing

    for (i1 = 1; i1 <= HWires->NumCurrentSegments; i1++) {
        if (!(HWires->CurrentSegment(i1)->chargeplus->exists && HWires->CurrentSegment(i1)->chargeminus->exists)) {
            sprintf(buff, "wir3_BUGGYERROR: Bug in WIRE node assignment. %7i%7i%7i%7i%s%3i%3i%s",
                    "wir3_BUGGYERROR: Bug in WIRE node assignment. ",
                    nodo->indexnode,
                    nodo->i,
                    nodo->j,
                    nodo->k,
                    " (",
                    nodo->numcurrentminus,
                    nodo->numcurrentplus,
                    ")");
            WarnErrReport(buff, true);
        }
    }

    // !!!!writes the lines in a DXF file if requested with -map
    // !!!    do i1=1,HWires%NumCurrentSegments
    // !!!        write(dxfbuff,'(a)') 'LINE'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '8'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%indexmed+20
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '62'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%indexmed+20
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!    select case(HWires%CurrentSegment(i1)%tipofield)
    // !!!    case(iEx)
    // !!!        write(dxfbuff,'(a)') '10'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%i
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '20'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%J
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '30'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%K
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '11'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%i+1
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '21'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%J
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '31'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%K
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '0'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!    case(iEy)
    // !!!        write(dxfbuff,'(a)') '10'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%i
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '20'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%J
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '30'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%K
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '11'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%i
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '21'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%J+1
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '31'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%K
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '0'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!    case(iEz)
    // !!!        write(dxfbuff,'(a)') '10'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%i
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '20'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%J
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '30'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%K
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '11'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%i
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '21'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%J
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '31'
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(i7)') HWires%CurrentSegment(i1)%K+1
    // !!!        call DXFWRITE(DXFBUFF)
    // !!!        write(dxfbuff,'(a)') '0'

#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, ierr);
#endif

      return;
   } // end subroutine ReportWireJunctions

   // Functions to get the values of the non diagonal elements of the autoinduction
   double F(double A, double B, double a1, double a2, double d, double phi)
   {
      int i, j;
      double frac;

      frac = 1.0 / (8.0 * A * B);
      double F_val = pi * a2 * a2 * (1.0 - 2.0 * log(a2 / d)) - 3.0 / (2.0 * frac);
      for (i = 1; i <= 2; i++) {
         for (j = 1; j <= 2; j++) {
            F_val += Gkl(i, j, A, B, d, phi) + Hkl(i, j, A, B, d, phi);
         }
      }

      F_val *= frac;
      return F_val;
   }

   double Gkl(int k, int l, double A, double B, double d, double phi)
   {
      double Ak, Bl;
      Ak = Ai(k, A, d, phi);
      Bl = Bi(l, B, d, phi);
      return Ak * Bl * log((Ak * Ak + Bl * Bl) / (d * d));
   }

   double Hkl(int k, int l, double A, double B, double d, double phi)
   {
      double Ak, Bl;
      Ak = Ai(k, A, d, phi);
      Bl = Bi(l, B, d, phi);
      return Ak * Ak * atan(Bl / Ak) + Bl * Bl * atan(Ak / Bl);
   }

   double Ai(int i, double A, double d, double phi)
   {
      double ai = -1.0;
      switch (i) {
         case 1:
            ai = A - d * cos(phi);
            break;
         case 2:
            ai = A + d * cos(phi);
            break;
      }
      return ai;
   }

   double Bi(int i, double B, double d, double phi)
   {
      double bi = -1.0;
      switch (i) {
         case 1:
            bi = B - d * sin(phi);
            break;
         case 2:
            bi = B + d * sin(phi);
            break;
      }
      return bi;
   }

   // Matrix inversion routine
   void MatInv(int N, std::vector<std::vector<double>>& M)
   {
      std::vector<double> P(N);
      std::vector<double> y(N);
      std::vector<std::vector<double>> B(N, std::vector<double>(N));
      std::vector<std::vector<double>> eye(N, std::vector<double>(N));

      for (int i = 0; i < N; i++) {
         for (int j = 0; j < N; j++) {
            eye[i][j] = 0.0;
         }
         eye[i][i] = 1.0;
      }

      for (int i = 0; i < N; i++) {
         P[i] = i;
      }

      for (int k = 0; k < N; k++) {
         double val = 0.0;
         int pivot = -1;
         for (int i = k; i < N; i++) {
            if (abs(B[i][k]) > val) {
               val = abs(B[i][k]);
               pivot = i;
            }
         }

         if (val == 0.0) {
            std::string buff = "WIR1_ERROR: Inversion de matriz fallida";
            WarnErrReport(buff, true);
            return;
         }

         int tmpi = P[k];
         P[k] = P[pivot];
         P[pivot] = tmpi;

         for (int i = 0; i < N; i++) {
            double tmpr = B[k][i];
            B[k][i] = B[pivot][i];
            B[pivot][i] = tmpr;
         }

         for (int i = k + 1; i < N; i++) {
            B[i][k] = B[i][k] / B[k][k];
            for (int j = k + 1; j < N; j++) {
               B[i][j] = B[i][j] - B[i][k] * B[k][j];
            }
         }
      }

      for (int k = 0; k < N; k++) {
         for (int i = 0; i < N; i++) {
            y[i] = eye[P[i]][k];
            for (int j = 0; j < i; j++) {
               y[i] = y[i] - B[i][j] * y[j];
            }
         }
         for (int i = N - 1; i >= 0; i--) {
            M[i][k] = y[i];
            for (int j = i + 1; j < N; j++) {
               M[i][k] = M[i][k] - B[i][j] * M[j][k];
            }
            M[i][k] = M[i][k] / B[i][i];
         }
      }
   }

   double Lambert(double z)
   {
      double x = z;
      double expx, newx;
      double eps = abs(x / 100000.0);
      int maxn = 1 << 12;
      int n = 1;

      if (z < -1.0 / 3.0) {
         return -1.0; // no va a converger nunca
      }

      while ((abs(x * exp(x) - z) > eps) && (n < maxn)) {
         expx = exp(x);
         // Newton method alternative
         newx = x - (x * expx - z) / (expx * (x + 1.0) - (x + 2.0) * (x * expx - z) / (2.0 * x + 2.0));
         x = newx;
         n++;
      }

      expx = exp(x);
      if (n >= maxn) {
         return -1.0;
      }
      return x;
   }

   // Tridiagonal solver
   void solve_tridiag_wires(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d, std::vector<double>& x, int n)
   {
      std::vector<double> cp(n);
      std::vector<double> dp(n);
      double m;

      // initialize c-prime and d-prime
      cp[0] = c[0] / b[0];
      dp[0] = d[0] / b[0];

      // solve for vectors c-prime and d-prime
      for (int i = 1; i < n; i++) {
         m = b[i] - cp[i - 1] * a[i];
         cp[i] = c[i] / m;
         dp[i] = (d[i] - dp[i - 1] * a[i]) / m;
      }

      // initialize x
      x[n - 1] = dp[n - 1];

      // solve for x from the vectors c-prime and d-prime
      for (int i = n - 2; i >= 0; i--) {
         x[i] = dp[i] - cp[i] * x[i + 1];
      }
   }

   void wiresconstantes(bool fieldtotl, CurrentSegments_t* dummy, const std::vector<double>& G2, const SGGFDTDINFO_t& sgg)
   {
      if (!fieldtotl) {
         dummy->cte5 = G2(dummy->indexmed) / (dummy->deltaTransv1 * dummy->deltaTransv2);
         // dummy->cte5 = dummy->cte5 / 3.5;
      } else {
         dummy->cte5 = 0.0; // esta es la cte de acoplo wire2main
      }

      // Lind ya contiene el givenautoin y a los autoin LeftEnd y RightEnd, que entiendo que de escalarse el mu tambien afectaria a todos ellos !ojo
      dummy->cte1 = ((dummy->Lind) / sgg.dt - dummy->Resist / 2.0) / ((dummy->Lind) / sgg.dt + dummy->Resist / 2.0);
      dummy->cte3 = InvMu(dummy->indexmed) * InvEps(dummy->indexmed) / (dummy->delta) * (dummy->Lind) / ((dummy->Lind) / sgg.dt + dummy->Resist / 2.0);
      // dummy->cte3 = dummy->cte3 / 3.5;
      dummy->cte2 = 1.0 / ((dummy->Lind) / sgg.dt + dummy->Resist / 2.0);

#ifdef CompileWithStochastic
      calc_wirehollandconstants_for_devia(sgg, dummy, InvEps, InvMu);
#endif
   }