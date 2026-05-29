#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <algorithm>

// Assuming necessary headers and types are defined in the environment
// This is a translation of the Fortran module HollandWires_m

namespace HollandWires_m {

    // Constants and types assumed from other modules
    // Using double for real(kind=RKIND_wires) and int for integer
    // Using std::vector for arrays
    
    // Placeholder types to match Fortran structures
    struct limit_t {
        int XI, XE, YI, YE, ZI, ZE;
    };

    struct SGGFDTDINFO_t {
        int Alloc[10]; // Simplified, assuming indices like iHx, iHy etc map to specific slots
        double dt;
        struct {
            int ZI, ZE;
        } Sweep[10];
        struct {
            double Epr;
            double Mur;
        } Med[100]; // Assuming max NumMedia
        int NumMedia;
    };

    struct sim_control_t {
        int layoutnumber;
        int num_procs;
    };

    struct Thinwires_t {
        int olddt;
        int NumNeededCurrentUpMPI;
        int NumNeededCurrentDownMPI;
        // Other members would be defined here based on usage
        struct {
            int indexnode;
        } NullNode;
    };

    struct CurrentSegments_t {
        // Members not fully specified in snippet
    };

    struct ChargeNodes_t {
        // Members not fully specified in snippet
    };

    struct adyacc_t {
        int YESsegment[2];
    };

    // Global variables from the module
    bool thereAreVsources = false;
    bool thereAreIsources = false;
    bool thereAreMurConditions = false;

    Thinwires_t HWires;
    
    std::vector<double> InvEps;
    std::vector<double> InvMu;
    std::vector<double> OldInvEps;
    std::vector<double> OldInvMu;

    double eps0 = 0.0;
    double mu0 = 0.0;

    // Function declarations
    void InitWires(
        SGGFDTDINFO_t& sgg,
        const std::vector<int>& sggMiNo,
        const std::vector<int>& sggMiEx,
        const std::vector<int>& sggMiEy,
        const std::vector<int>& sggMiEz,
        const std::vector<int>& sggMiHx,
        const std::vector<int>& sggMiHy,
        const std::vector<int>& sggMiHz,
        bool& ThereAreWires,
        const std::vector<double>& Ex,
        const std::vector<double>& Ey,
        const std::vector<double>& Ez,
        const std::vector<double>& Hx,
        const std::vector<double>& Hy,
        const std::vector<double>& Hz,
        const std::vector<int>& Idxe,
        const std::vector<int>& Idye,
        const std::vector<int>& Idze,
        const std::vector<int>& Idxh,
        const std::vector<int>& Idyh,
        const std::vector<int>& Idzh,
        const std::vector<double>& G2,
        const std::vector<limit_t>& SINPML_fullsize,
        const std::vector<limit_t>& fullsize,
        double& dtcritico,
        double eps00,
        double mu00,
        const sim_control_t& control
    );

    void AdvanceWiresE(...); // Placeholder
    void AdvanceWiresH(...); // Placeholder
    void AdvanceWiresEcrank(...); // Placeholder
    void StoreFieldsWires(...); // Placeholder
    void DestroyWires(...); // Placeholder
    void GetHwires(...); // Placeholder
    void ReportWireJunctions(...); // Placeholder
    void calc_wirehollandconstants(...); // Placeholder

    // Implementation of InitWires
    void InitWires(
        SGGFDTDINFO_t& sgg,
        const std::vector<int>& sggMiNo,
        const std::vector<int>& sggMiEx,
        const std::vector<int>& sggMiEy,
        const std::vector<int>& sggMiEz,
        const std::vector<int>& sggMiHx,
        const std::vector<int>& sggMiHy,
        const std::vector<int>& sggMiHz,
        bool& ThereAreWires,
        const std::vector<double>& Ex,
        const std::vector<double>& Ey,
        const std::vector<double>& Ez,
        const std::vector<double>& Hx,
        const std::vector<double>& Hy,
        const std::vector<double>& Hz,
        const std::vector<int>& Idxe,
        const std::vector<int>& Idye,
        const std::vector<int>& Idze,
        const std::vector<int>& Idxh,
        const std::vector<int>& Idyh,
        const std::vector<int>& Idzh,
        const std::vector<double>& G2,
        const std::vector<limit_t>& SINPML_fullsize,
        const std::vector<limit_t>& fullsize,
        double& dtcritico,
        double eps00,
        double mu00,
        const sim_control_t& control
    ) {
        double eps000 = eps00; // Dummy variables
        double mu000 = mu00;   // Dummy variables

        // Note: In Fortran, these are passed by reference and modified/used.
        // In C++, we assume the vectors are sized appropriately or we access them with bounds checking.
        // The Fortran code uses pointer association for some arrays, which is complex to map directly.
        // For this translation, we assume the vectors passed contain the data in the correct order.

        ThereAreWires = false;

        // chapuz para convertir la variables de paso en globales
        eps0 = eps00;
        mu0 = mu00;

        // guarda el dt original (para permit scaling)
        HWires.olddt = sgg.dt;

        int ZI = sgg.Sweep[2].ZI; // Assuming iHz maps to index 2 or similar, need exact mapping
        int ZE = sgg.Sweep[2].ZE;

        // dir array initialization
        std::string dir[10];
        dir[0] = " X "; // iEx
        dir[1] = " Y "; // iEy
        dir[2] = " Z "; // iEz

        thereAreVsources = false;
        thereAreIsources = false;
        thereAreMurConditions = false;

        // whoami string construction
        std::string whoami = "(" + std::to_string(control.layoutnumber + 1) + "/" + std::to_string(control.num_procs) + ") ";

#ifdef CompileWithMPI
        int ierr = 0;
        // MPI_Barrier call would go here
        // call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

        HWires.NullNode.indexnode = -1;

        // Initialize adjacency
        HWires.NullNode.indexnode = -1; // Already set above
        
        // Reset pointers (simulated by nulling out or clearing vectors if they were pointers)
        // HWires.WireTipoMedio => null()
        // HWires.CurrentSegment => null()
        // HWires.ChargeNode => null()

        dtcritico = sgg.dt;

        // Initialize InvEps and InvMu
        // allocate (InvEps(0 : sgg%NumMedia),InvMu(0 : sgg%NumMedia))
        InvEps.assign(sgg.NumMedia + 1, 0.0);
        InvMu.assign(sgg.NumMedia + 1, 0.0);

        for (int i = 0; i <= sgg.NumMedia; ++i) {
            InvEps[i] = 1.0 / (eps0 * sgg.Med[i].Epr);
            InvMu[i] = 1.0 / (mu0 * sgg.Med[i].Mur);
        }

        // Initialize OldInvEps and OldInvMu
        OldInvEps = InvEps;
        OldInvMu = InvMu;

#ifdef CompileWithMPI
        Hwires.NumNeededCurrentUpMPI = 0;
        Hwires.NumNeededCurrentDownMPI = 0;
#endif

        // Trivial initialization
        HWires.NullNode.indexnode = -1;
    }

} // namespace HollandWires_m

HWires.NullNode.ChargePresent = 0.0;
        HWires.NullNode.ChargePast = 0.0;
        HWires.NullNode.IsAttachedtoVoltage = false;
        HWires.NullNode.IsMur = false;
        HWires.NullNode.IsBackDownLeftMur = false;
        HWires.NullNode.IsFrontUpRightMur = false;
        HWires.NullNode.IsPeriodic = false;
        HWires.NullNode.IsPEC = false;
        HWires.NullNode.already_YEEadvanced_byconformal_changedtoPECfield1 = nullptr;
        HWires.NullNode.already_YEEadvanced_byconformal_changedtoPECfield2 = nullptr;
        HWires.NullNode.already_YEEadvanced_byconformal_changedtoPECfield3 = nullptr;
        HWires.NullNode.already_YEEadvanced_byconformal_changedtoPECfield4 = nullptr;
        HWires.NullNode.already_YEEadvanced_byconformal_changedtoPECfield5 = nullptr;
        HWires.NullNode.already_YEEadvanced_byconformal_changedtoPECfield6 = nullptr;
        HWires.NullNode.IsLossy = false;
        HWires.NullNode.HasIsource = false;
        HWires.NullNode.IsHeterogeneousJunction = false;
        //just for informative !not implemented in MPI unsure behaviour under mpi 2011 \E7
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

        HWires.NullSegment.R = 0.0;
        HWires.NullSegment.Resist = 0.0;
        HWires.NullSegment.Resist_devia = 0.0;
        HWires.NullSegment.C = 0.0;
        HWires.NullSegment.L = 0.0;
        HWires.NullSegment.Lintrinsic = 0.0;
        HWires.NullSegment.NumParallel = 1;
        HWires.NullSegment.origindex = i1;
        HWires.NullSegment.indexsegment = i1;
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
        //!!! HWires.NullSegment.logRoverR0 = 0.0;
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

        ThereAreWires = false;

        //detect thin wires : same radius implies same medium independently of its orientation
        conta = 0;
        for (jmed = 1; jmed <= sgg.NumMedia; ++jmed) {
            if (sgg.Med[jmed].Is.ThinWire) conta++;
        }

        HWires.NumDifferentWires = conta;
        HWires.WireTipoMedio.resize(HWires.NumDifferentWires);
        conta = 0;
        for (jmed = 1; jmed <= sgg.NumMedia; ++jmed) {
            if (sgg.Med[jmed].Is.ThinWire) {
                ThereAreWires = true;
                conta++;
                HWires.WireTipoMedio[conta - 1] = jmed;
            }
        }

        //!!!chequeo previo 210323 para que no haya conectores dispersivos
        for (iwi = 1; iwi <= HWires.NumDifferentWires; ++iwi) {
            if ((sgg.Med[HWires.WireTipoMedio[iwi - 1]].wire[1].disp.size() > 0) ||
                (sgg.Med[HWires.WireTipoMedio[iwi - 1]].wire[1].disp_RightEnd.size() > 0) ||
                (sgg.Med[HWires.WireTipoMedio[iwi - 1]].wire[1].disp_LeftEnd.size() > 0)) {
                std::string buff = "Dispersive wire or connectors unsupported in Holland wires";
                WarnErrReport(buff, true);
            }
        }

        for (iwi = 1; iwi <= HWires.NumDifferentWires; ++iwi) {
            for (iwj = 1; iwj <= sgg.Med[HWires.WireTipoMedio[iwi - 1]].wire[1].numsegmentos; ++iwj) {
                sgg.Med[HWires.WireTipoMedio[iwi - 1]].wire[1].segm[iwj - 1].multirabo = false;
            }
        }

        if (!therearewires) return;

        std::string buff = "----------------------------------------------------------------";
        WarnErrReport(buff);

        // it directly reads the segments specified in the .nfde file

        // detects endings and set ending=.true.
        // esto implica que podra haber LeftEnd y RightEnd declarados y ademas ending

        for (iwi = 1; iwi <= HWires.NumDifferentWires; ++iwi) {
            for (iwj = 1; iwj <= sgg.Med[HWires.WireTipoMedio[iwi - 1]].wire[1].numsegmentos; ++iwj) {
                if (!control.strictOLD) {
                    sgg.Med[HWires.WireTipoMedio[iwi - 1]].wire[1].segm[iwj - 1].ilibre = -1;
                    sgg.Med[HWires.WireTipoMedio[iwi - 1]].wire[1].segm[iwj - 1].jlibre = -1;
                    sgg.Med[HWires.WireTipoMedio[iwi - 1]].wire[1].segm[iwj - 1].klibre = -1;
                    // (en formatos antiguas esta info solo viene bien si hay cargas (md me las pone al final y al principio) pero en general puede estar mal y hay que resetearla
                    if (control.connectendings) {
                        if (sgg.Med[HWires.WireTipoMedio[iwi - 1]].wire[1].numsegmentos == 1) {
                            sgg.Med[HWires.WireTipoMedio[iwi - 1]].wire[1].segm[iwj - 1].Is_LeftEnd = true;
                            sgg.Med[HWires.WireTipoMedio[iwi - 1]].wire[1].segm[iwj - 1].Is_RightEnd = true;
                            // ojooo esto esta bien? parecen los TL y TR intercambiados sgg 251019 pero no lo toco
                            if (sgg.Med[HWires.WireTipoMedio[iwi - 1]].wire[1].HasAbsorbing_RightEnd) {
                                sgg.Med[HWires.WireTipoMedio[iwi - 1]].wire[1].HasAbsorbing_LeftEnd = true;
                            }
                        }
                    }
                }
            }
        }

if (sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].HasAbsorbing_LeftEnd) 
                            sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].HasAbsorbing_RightEnd = true;
                        if (sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].HasParallel_RightEnd) 
                            sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].HasParallel_LeftEnd = true;
                        if (sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].HasParallel_LeftEnd) 
                            sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].HasParallel_RightEnd = true;
                        if (sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].HasSeries_RightEnd) 
                            sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].HasSeries_LeftEnd = true;
                        if (sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].HasSeries_LeftEnd) 
                            sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].HasSeries_RightEnd = true;
                    }
                    
                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_LeftEnd = 
                        sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_LeftEnd &&
                        (sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].HasParallel_LeftEnd || 
                         sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].HasSeries_LeftEnd ||
                         sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].HasAbsorbing_LeftEnd);
                    
                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_RightEnd = 
                        sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_RightEnd &&
                        (sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].HasParallel_RightEnd || 
                         sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].HasSeries_RightEnd ||
                         sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].HasAbsorbing_RightEnd);
                }
                
                sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].IsEnd_norLeft_norRight = 
                    !(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_LeftEnd || 
                      sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_RightEnd);
                
                i1 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].i;
                j1 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].j;
                k1 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].k;
                whatfield = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ori;
                origindex = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].origIndex;
                
                conectado1 = false;
                conectado2 = false;
                
                for (int iwjjj = 1; iwjjj <= sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].numsegmentos; ++iwjjj) {
                    i2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].i;
                    j2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].j;
                    k2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].k;
                    whatfield2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].ori;
                    origindex2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].origIndex;
                    
                    switch (whatfield) {
                        case iEx:
                            switch (whatfield2) {
                                case iEx:
                                    conectado1 = conectado1 || ((i1 == i2 + 1) && (j1 == j2) && (k1 == k2));
                                    conectado2 = conectado2 || ((i1 + 1 == i2) && (j1 == j2) && (k1 == k2));
                                    break;
                                case iEy:
                                    conectado1 = conectado1 || ((i1 == i2) && (j1 == j2) && (k1 == k2));
                                    conectado1 = conectado1 || ((i1 == i2) && (j1 == j2 + 1) && (k1 == k2));
                                    conectado2 = conectado2 || ((i1 + 1 == i2) && (j1 == j2) && (k1 == k2));
                                    conectado2 = conectado2 || ((i1 + 1 == i2) && (j1 == j2 + 1) && (k1 == k2));
                                    break;
                                case iEz:
                                    conectado1 = conectado1 || ((i1 == i2) && (j1 == j2) && (k1 == k2));
                                    conectado1 = conectado1 || ((i1 == i2) && (j1 == j2) && (k1 == k2 + 1));
                                    conectado2 = conectado2 || ((i1 + 1 == i2) && (j1 == j2) && (k1 == k2));
                                    conectado2 = conectado2 || ((i1 + 1 == i2) && (j1 == j2) && (k1 == k2 + 1));
                                    break;
                            }
                            break;
                        case iEy:
                            switch (whatfield2) {
                                case iEy:
                                    conectado1 = conectado1 || ((j1 == j2 + 1) && (k1 == k2) && (i1 == i2));
                                    conectado2 = conectado2 || ((j1 + 1 == j2) && (k1 == k2) && (i1 == i2));
                                    break;
                                case iEz:
                                    conectado1 = conectado1 || ((j1 == j2) && (k1 == k2) && (i1 == i2));
                                    conectado1 = conectado1 || ((j1 == j2) && (k1 == k2 + 1) && (i1 == i2));
                                    conectado2 = conectado2 || ((j1 + 1 == j2) && (k1 == k2) && (i1 == i2));
                                    conectado2 = conectado2 || ((j1 + 1 == j2) && (k1 == k2 + 1) && (i1 == i2));
                                    break;
                                case iEx:
                                    conectado1 = conectado1 || ((j1 == j2) && (k1 == k2) && (i1 == i2));
                                    conectado1 = conectado1 || ((j1 == j2) && (k1 == k2) && (i1 == i2 + 1));
                                    conectado2 = conectado2 || ((j1 + 1 == j2) && (k1 == k2) && (i1 == i2));
                                    conectado2 = conectado2 || ((j1 + 1 == j2) && (k1 == k2) && (i1 == i2 + 1));
                                    break;
                            }
                            break;
                        case iEz:
                            switch (whatfield2) {
                                case iEz:
                                    conectado1 = conectado1 || ((k1 == k2 + 1) && (i1 == i2) && (j1 == j2));
                                    conectado2 = conectado2 || ((k1 + 1 == k2) && (i1 == i2) && (j1 == j2));
                                    break;
                                case iEx:
                                    conectado1 = conectado1 || ((k1 == k2) && (i1 == i2) && (j1 == j2));
                                    conectado1 = conectado1 || ((k1 == k2) && (i1 == i2 + 1) && (j1 == j2));
                                    conectado2 = conectado2 || ((k1 + 1 == k2) && (i1 == i2) && (j1 == j2));
                                    conectado2 = conectado2 || ((k1 + 1 == k2) && (i1 == i2 + 1) && (j1 == j2));
                                    break;
                                case iEy:
                                    conectado1 = conectado1 || ((k1 == k2) && (i1 == i2) && (j1 == j2));
                                    conectado1 = conectado1 || ((k1 == k2) && (i1 == i2) && (j1 == j2 + 1));
                                    conectado2 = conectado2 || ((k1 + 1 == k2) && (i1 == i2) && (j1 == j2));
                                    conectado2 = conectado2 || ((k1 + 1 == k2) && (i1 == i2) && (j1 == j2 + 1));
                                    break;
                            }
                            break;
                    }
                    
                    conectado = conectado1 && conectado2;
                    if (conectado) {
                        break;
                    }
                }
                
                if (control.connectendings) {
                    if ((!sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_LeftEnd) &&
                        (!sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_RightEnd)) {
                        sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].IsEnd_norLeft_norRight = 
                            sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].IsEnd_norLeft_norRight ||
                            ((!conectado) && conectado2) &&

} else if (!(sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].HasParallel_LeftEnd || sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].HasSeries_LeftEnd || sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].HasAbsorbing_LeftEnd)) {
                        }
                        if ((!sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_LeftEnd) && (!sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_RightEnd)) {
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].IsEnd_norLeft_norRight = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].IsEnd_norLeft_norRight || (((!conectado) && conectado1) && (!sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].HasParallel_RightEnd || sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].HasSeries_RightEnd || sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].HasAbsorbing_RightEnd));
                        }
                    }
                    sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].IsEnd_norLeft_norRight = (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].IsEnd_norLeft_norRight) && (!conectado) && (conectado1 || conectado2) && (!sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_LeftEnd || sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_RightEnd); //si hay mas de uno este se pone a .true.
                    //detecta cual es el extremo libre
                    if ((!conectado) && conectado1) {
                        switch (whatfield) {
                        case iEx:
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ilibre = i1 + 1;
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].jlibre = j1;
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].klibre = k1;
                            break;
                        case iEy:
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ilibre = i1;
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].jlibre = j1 + 1;
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].klibre = k1;
                            break;
                        case iEz:
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ilibre = i1;
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].jlibre = j1;
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].klibre = k1 + 1;
                            break;
                        }
                    } else if ((!conectado) && conectado2) {
                        sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ilibre = i1;
                        sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].jlibre = j1;
                        sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].klibre = k1;
                    }

                    //caso especial
                    if (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].numsegmentos == 1) {
                        sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].IsEnd_norLeft_norRight = false;
                        sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_LeftEnd = true;
                        sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_RightEnd = true;
                        switch (whatfield) {
                        case iEx:
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ilibre = i1 + 1;
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].jlibre = j1;
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].klibre = k1;
                            break;
                        case iEy:
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ilibre = i1;
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].jlibre = j1 + 1;
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].klibre = k1;
                            break;
                        case iEz:
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ilibre = i1;
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].jlibre = j1;
                            sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].klibre = k1 + 1;
                            break;
                        }
                    }
                    //check for intermediate RLC error

                    if ((sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ilibre == -1) || (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].jlibre == -1) || (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].klibre == -1)) {
                        if (!conectado) {
                            sprintf(buff, "wir0_BUGGYERROR: Non-Intermediate multi-segment WIRE. %7d%7d%7d%7d%s", origIndex, i1, j1, k1, dir(whatfield).c_str());
                            if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff, true);
                        }
                        if (((sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_RightEnd) && (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].HasParallel_RightEnd || sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].HasSeries_RightEnd)) || ((sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_LeftEnd) && (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].HasParallel_LeftEnd || sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].HasSeries_LeftEnd || sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].HasAbsorbing_LeftEnd))) {
                            if (conectado) {
                                sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_RightEnd = false;
                                sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_LeftEnd = false;
                                sprintf(buff, "wir0_WARNING: Intermediate segment with RLC. Neglecting RLC %7d%7d%7d%7d%s", origIndex, i1, j1, k1, dir(whatfield).c_str());
                                if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff);
                            } else {
                                sprintf(buff, "wir0_BUGGYERROR: Non-Intermediate multi-segment WIRE with RLC. %7d%7d%7d%7d%s", origIndex, i1, j1, k1, dir(whatfield).c_str());
                                if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff, true);
                            }
                        }
                    }
                    //
                } else { //del strictOLD

                    sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ilibre = -1;
                    sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].jlibre = -1;
                    sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].klibre = -1;
                    sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].IsEnd_norLeft_norRight = false; //irrelevante en strictOLD
                    i1 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].i;
                    j1 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].j;
                    k1 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].k;
                    whatfield = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ori;
                    origindex = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].origIndex;
                    //
                    if (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_LeftEnd) {
                        dummy1 = 1;
                        dummyfin = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].numsegmentos;
                    } else if (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_RightEnd) {
                        dummy1 = -1;
                        dummyfin = 1;
                    }
                }

if ( ((sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_LeftEnd) || (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_RightEnd)) ) {
                     if (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].numsegmentos != 1) {
                        dummy2 = -1;
                        for (int iwjjj = iwj + dummy1; iwjjj <= dummyfin; iwjjj += dummy1) { //atras o adelante
                           i2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].i;
                           j2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].j;
                           k2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].k;
                           whatfield2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].ori;
                           if ((i1 == i2) && (j1 == j2) && (k1 == k2) && (whatfield == whatfield2)) {
                              dummy2 = -dummy2; //detecta numero de rabitos para o impar aunque yo luego en las uniones solo trato 2 rabitos como mucho
                           } else {
                              continue;
                              break;
                           }
                        }
                        if (dummy2 == -1) {
                           dummy3 = 0;
                        } else {
                           dummy3 = 1;
                        }

                        //
                        sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ilibre = i1;
                        sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].jlibre = j1;
                        sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].klibre = k1;
                        if (whatfield2 == whatfield) {
                           if (abs(i1 - i2) + abs(j1 - j2) + abs(k1 - k2) > 1) {
                              sprintf(buff, "wir0_ERROR: strictOLD LeftEnd/RightEnd segment disconnected. %7d%7d%7d%7d %s", origIndex, i1, j1, k1, dir(whatfield).c_str());
                              if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff, true);
                           }
                           if (i1 > i2) {
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ilibre = i1 + 1 - dummy3;
                           } else if (i1 < i2) {
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ilibre = i1 + dummy3;
                           }
                           if (j1 > j2) {
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].jlibre = j1 + 1 - dummy3;
                           } else if (j1 < j2) {
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].jlibre = j1 + dummy3;
                           }
                           if (k1 > k2) {
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].klibre = k1 + 1 - dummy3;
                           } else if (k1 < k2) {
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].klibre = k1 + dummy3;
                           }
                        } else {
                           if (abs(i1 - i2) + abs(j1 - j2) + abs(k1 - k2) > 2) {
                              sprintf(buff, "wir0_ERROR: strictOLD LeftEnd/RightEnd segment disconnected. %7d%7d%7d%7d %s", origIndex, i1, j1, k1, dir(whatfield).c_str());
                              if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff, true);
                           }
                           switch (whatfield) {
                            case iEx:
                              if (i2 == i1) sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ilibre = i1 + 1;
                              break;
                            case iEy:
                              if (j2 == j1) sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].jlibre = j1 + 1;
                              break;
                            case iEz:
                              if (k2 == k1) sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].klibre = k1 + 1;
                              break;
                           }
                        }
                     } else { //DEL NUMERO SEGMENTOS '2014 NO PORTADO A !CHECK
                        switch (whatfield) {
                         case iEx:
                           sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ilibre = i1 + 1;
                           sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].jlibre = j1;
                           sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].klibre = k1;
                           break;
                         case iEy:
                           sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ilibre = i1;
                           sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].jlibre = j1 + 1;
                           sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].klibre = k1;
                           break;
                         case iEz:
                           sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ilibre = i1;
                           sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].jlibre = j1;
                           sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].klibre = k1 + 1;
                           break;
                        }
                     } //DEL NUMERO SEGMENTOS '2014 NO PORTADO A !CHECK
                  } //Del Left_End Right_End
                  //!!!!!!!!!!!!!!!!!!
               } // del strictOLD

            } // end do
         } // end do


         //preprocesa para eliminar multirabos luego se utiliza en repetido

         if (control.strictOLD) {
            for (int iwi = 1; iwi <= HWires.NumDifferentWires; iwi++) {
               for (int iwj = 1; iwj <= sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].numsegmentos; iwj++) {
                  if (!sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].multirabo) {
                     multirabos = 1;
                     i1 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].i;
                     j1 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].j;
                     k1 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].k;
                     whatfield = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ori;
                     ORIGINDEX = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].origindex;
                     //precontaje
                     for (int iwjjj = iwj + 1; iwjjj <= sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].numsegmentos; iwjjj++) { //el Right_End aunque no se tape si debe detectarse a efectos par/impar
                        i2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].i;
                        j2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].j;
                        k2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].k;
                        whatfield2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].ori;
                        ORIGINDEX2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].origindex;
                        if ((i1 == i2) && (j1 == j2) && (k1 == k2) && (whatfield == whatfield2)) {
                           multirabos = multirabos + 1;
                        } else {
                           primernorabo = origindex2;
                           Jprimernorabo = iwjjj;
                           break;
                        }
                     }
                     //machaca rabos
                     if (multirabos != 1) {
                        //
                        for (int iwjjj = iwj + (2 - mod(multirabos, 2)); iwjjj <= sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].numsegmentos - 1; iwjjj++) { //el Right_End no debe taparse ya se tapa el de dentro en el otro bucle
                           i2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].i;
                           j2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].j;
                           k2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].k;
                           whatfield2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].ori;
                           ORIGINDEX2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].origindex;

if ((i1 == i2) && (j1 == j2) && (k1 == k2) && (whatfield == whatfield2)) {
                                std::ostringstream buff_stream;
                                buff_stream << "wir0_WARNING: strictOLD Redundannt zig-zag rabito, will be eliminated to Mod(2)" << std::setw(7) << origIndex2 
                                            << std::setw(7) << i2 << std::setw(7) << j2 << std::setw(7) << k2 
                                            << dir(whatfield2) << "by segment " << primernorabo;
                                std::string buff = buff_stream.str();
                                if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff);
                                sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].multirabo = true;
                                sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].multiraboDE = primernorabo;
                            } else {
                                goto end_taparabos;
                            }
                        }
                        end_taparabos:;
                    }
                }
            }
            //ahora la vuelta
            for (iwi = 1; iwi <= HWires.NumDifferentWires; ++iwi) {
                for (iwj = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].numsegmentos; iwj >= 1; --iwj) {
                    if (!sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].multirabo) {
                        multirabos = 1;
                        i1 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].i;
                        j1 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].j;
                        k1 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].k;
                        whatfield = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ori;
                        ORIGINDEX = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].origindex;
                        for (iwjjj = iwj - 1; iwjjj >= 1; --iwjjj) { //el Left_End aunque no se tape si debe detectarse a efectos par/impar
                            i2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].i;
                            j2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].j;
                            k2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].k;
                            whatfield2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].ori;
                            ORIGINDEX2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].origindex;
                            if ((i1 == i2) && (j1 == j2) && (k1 == k2) && (whatfield == whatfield2)) {
                                multirabos = multirabos + 1;
                            } else {
                                primernorabo = origindex2;
                                Jprimernorabo = iwjjj;
                                goto end_buscarabos2;
                            }
                        }
                        end_buscarabos2:;

                        //machaca rabos
                        if (multirabos != 1) {
                            //
                            for (iwjjj = iwj - (2 - mod(multirabos, 2)); iwjjj >= 2; --iwjjj) { //el Left_End no debe taparse
                                i2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].i;
                                j2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].j;
                                k2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].k;
                                whatfield2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].ori;
                                ORIGINDEX2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].origindex;
                                if ((i1 == i2) && (j1 == j2) && (k1 == k2) && (whatfield == whatfield2)) {
                                    //faltaria eliminar sondas en multirabos 7/2/14
                                    std::ostringstream buff_stream2;
                                    buff_stream2 << "wir0_WARNING: strictOLD Redundannt zig-zag rabito, will be eliminated to Mod(2)" << std::setw(7) << origIndex2 
                                                 << std::setw(7) << i2 << std::setw(7) << j2 << std::setw(7) << k2 
                                                 << dir(whatfield2) << "by segment " << primernorabo;
                                    std::string buff2 = buff_stream2.str();
                                    if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff2);
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].multirabo = true;
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].multiraboDE = primernorabo;
                                } else {
                                    goto end_taparabos2;
                                }
                            }
                            end_taparabos2:;
                        }
                    }
                }
            }


            //!!!!!!!!!!!!!!!!
            //segunda pasada para procesar TAPARRABOS
            //!!!!!!!!!!!!!!!!!!
            //!!!!!!!!!!!

            if (control.TAPARRABOS) {
                for (iwi = 1; iwi <= HWires.NumDifferentWires; ++iwi) {
                    for (iwj = 1; iwj <= sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].numsegmentos; ++iwj) {
                        if (!sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].multirabo) {
                            multirabos = 1;
                            i1 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].i;
                            j1 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].j;
                            k1 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].k;
                            whatfield = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ori;
                            ORIGINDEX = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].origindex;
                            //precontaje
                            for (iwjjj = iwj + 1; iwjjj <= sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].numsegmentos; ++iwjjj) { //el Right_End aunque no se tape si debe detectarse a efectos par/impar
                                i2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].i;
                                j2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].j;
                                k2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].k;
                                whatfield2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].ori;
                                ORIGINDEX2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].origindex;
                                if ((i1 == i2) && (j1 == j2) && (k1 == k2) && (whatfield == whatfield2)) {
                                    multirabos = multirabos + 1;
                                } else {
                                    primernorabo = origindex2;
                                    Jprimernorabo = IWjjj;
                                    goto end_buscarabos6;
                                }
                            }
                            end_buscarabos6:;
                            //machaca rabos
                            if ((mod(multirabos, 2) != 1) && (Jprimernorabo != sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].numsegmentos)) { //no al RightEnd
                                if (sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_LeftEnd) {
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_LeftEnd = false;
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[Jprimernorabo].Is_LeftEnd = true; //pasa caracter Left_End al primernorabo
                                    //!!!tocado esto por el problema de gra_powerline_simple.nfde 190916 el rabito se quedaba con el libre mal computado. Los ilibre,jlbre,klibre para Left_End son el primer punto directamente
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[Jprimernorabo].ilibre = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[Jprimernorabo].i; //!!!!sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].ilibre
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[Jprimernorabo].jlibre = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[Jprimernorabo].j; //!!!!sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].jlibre
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[Jprimernorabo].klibre = sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[Jprimernorabo].k; //!!!!sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].klibre
                                    //y retocado 200117 por el problema del ntc1 de la demo de getafe que tambien calculaba mal el libre (no tiene que ser por guevos el primer punto, habra que ver como se conecta con el siguiente!!!)
                                }
                            }
                        }
                    }
                }
            }

if ((sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].i) == (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo + 1].i) &&
                                     (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].j) == (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo + 1].j) &&
                                     (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].k) == (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo + 1].k)) {
                                     switch (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].ori) {
                                     case 1:
                                        sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].ilibre = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].ilibre + 1;
                                        break;
                                     case 2:
                                        sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].jlibre = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].jlibre + 1;
                                        break;
                                     case 3:
                                        sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].klibre = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].klibre + 1;
                                        break;
                                     }
                                }
                           }
                           //
                           if ((!sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_LeftEnd) && (!sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_RightEnd)) {
                              i2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].i;
                              j2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].j;
                              k2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].k;
                              whatfield2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ori;
                              ORIGINDEX2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].origindex;
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].multirabo = true;
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].multiraboDE = primernorabo; //el ultimo cazado es el primer NOrabo
                              sprintf(buff, "wir0_WARNING: strictOLD and taparrabos redundannt zig-zag rabito, will be ALSO eliminated %7d%7d%7d%7d%s%s%7d", origIndex2, i2, j2, k2, dir(whatfield2).c_str(), "by segment ", primernorabo);
                              if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff);
                           }
                           //
                           if ((!sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj + 1].Is_LeftEnd) && (!sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj + 1].Is_RightEnd)) {
                              i2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj + 1].i;
                              j2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj + 1].j;
                              k2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj + 1].k;
                              whatfield2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj + 1].ori;
                              ORIGINDEX2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj + 1].origindex;
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj + 1].multirabo = true;
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj + 1].multiraboDE = primernorabo; //el ultimo cazado es el primer NOrabo
                              sprintf(buff, "wir0_WARNING: strictOLD and taparrabos redundannt zig-zag rabito, will be ALSO eliminated %7d%7d%7d%7d%s%s%7d", origIndex2, i2, j2, k2, dir(whatfield2).c_str(), "by segment ", primernorabo);
                              if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff);
                           }
                        } //del multirabos no nulo
                     }
                  }
               }
               //ahora la vuelta
               for (iwi = 0; iwi < HWires.NumDifferentWires; ++iwi) {
                  for (int iwj = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].numsegmentos - 1; iwj >= 0; --iwj) {
                     if (!sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].multirabo) {
                        multirabos = 1;
                        i1 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].i;
                        j1 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].j;
                        k1 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].k;
                        whatfield = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ori;
                        ORIGINDEX = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].origindex;
                        for (int iwjjj = iwj - 1; iwjjj >= 0; --iwjjj) { //el Left_End aunque no se tape si debe detectarse a efectos par/impar
                           i2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].i;
                           j2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].j;
                           k2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].k;
                           whatfield2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].ori;
                           ORIGINDEX2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwjjj].origindex;
                           if ((i1 == i2) && (j1 == j2) && (k1 == k2) && (whatfield == whatfield2)) {
                              multirabos = multirabos + 1;
                           } else {
                              primernorabo = ORIGINDEX2;
                              Jprimernorabo = iwjjj;
                              break;
                           }
                        }

                        //machaca rabos
                        if (((multirabos % 2) != 1) && (Jprimernorabo != 0)) {   //no el Left_End
                           if (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_RightEnd) {
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_RightEnd = false;
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].Is_RightEnd = true; //pasa caracter RightEnd al primernorabo
                              //!!!tocado esto por el problema de gra_powerline_simple.nfde 190916 el rabito se quedaba con el libre mal computado. Los ilibre,jlbre,klibre para Left_End son el primer punto directamente
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].ilibre = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].i; //!!!!sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ilibre
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].jlibre = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].j; //!!!!sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].jlibre
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].klibre = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].k; //!!!!sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].klibre
                              //y retocado 200117 por el problema del ntc1 de la demo de getafe que tambien calculaba mal el libre (no tiene que ser por guevos el primer punto, habra que ver como se conecta con el siguiente!!!)
                              if ((sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].i) == (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo - 1].i) &&
                                  (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].j) == (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo - 1].j) &&

(sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].k) == (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo - 1].k)) {
                                     switch (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].ori) {
                                     case 1:
                                        sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].ilibre = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].ilibre + 1;
                                        break;
                                     case 2:
                                        sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].jlibre = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].jlibre + 1;
                                        break;
                                     case 3:
                                        sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].klibre = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[Jprimernorabo].klibre + 1;
                                        break;
                                     }
                                }
                           
                           }
                           
                           if ((!sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_LeftEnd) && (!sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_RightEnd)) {
                              i2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].i;
                              j2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].j;
                              k2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].k;
                              whatfield2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ori;
                              ORIGINDEX2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].origindex;
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].multirabo = true;
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].multiraboDE = primernorabo;
                              sprintf(buff, "wir0_WARNING: strictOLD and taparrabos redundannt zig-zag rabito, will be ALSO eliminated %7d%7d%7d%7d%s%s%7d", " ", origIndex2, i2, j2, k2, dir(whatfield2).c_str(), "by segment ", primernorabo);
                              if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff);
                           }
                           
                           if ((!sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj - 1].Is_LeftEnd) && (!sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj - 1].Is_RightEnd)) {
                              i2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj - 1].i;
                              j2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj - 1].j;
                              k2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj - 1].k;
                              whatfield2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj - 1].ori;
                              ORIGINDEX2 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj - 1].origindex;
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj - 1].multirabo = true;
                              sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj - 1].multiraboDE = primernorabo;
                              sprintf(buff, "wir0_WARNING: strictOLD and taparrabos redundannt zig-zag rabito, will be ALSO eliminated %7d%7d%7d%7d%s%s%7d", " ", origIndex2, i2, j2, k2, dir(whatfield2).c_str(), "by segment ", primernorabo);
                              if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff);
                           }
                        }
                     }
                  }
               }
            }
            }
            
         }
         
         do iwi = 1, HWires.NumDifferentWires {
            do iwj = 1, sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].numsegmentos {
               if ((sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].multirabo) && ((sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_LeftEnd) || (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_RightEnd))) {
                  i1 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].i;
                  j1 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].j;
                  k1 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].k;
                  whatfield = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ori;
                  ORIGINDEX = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].origindex;
                  sprintf(buff, "wir0_BuggyERROR: strictOLD LeftEnd/RightEnd cannot be multirabo %7d%7d%7d%7d%s", " ", origIndex, i1, j1, k1, (std::string(" ") + dir(whatfield)).c_str());
                  if ((k1 >= ZI) && (k1 <= ZE)) WarnErrReport(buff, true);
               }
            }
         }
         
         
         
         
         
         do iwi = 1, HWires.NumDifferentWires {
            do iwj = 1, sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].numsegmentos {
               i1 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].i;
               j1 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].j;
               k1 = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].k;
               i1libre = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ilibre;
               j1libre = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].jlibre;
               k1libre = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].klibre;
               whatfield = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].ori;
               origindex = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].origIndex;
               LeftEnd_index = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].LeftEnd;
               RightEnd_index = sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].RightEnd;
               
               if (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_LeftEnd && sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_RightEnd) {
                  sprintf(buff, "wir0_INFO: Ending segment (LeftEnd_and_Right)%7d%7d%7d%7d-%7d%7d%7d%s%7d%7d", " ", origIndex, i1, j1, k1, i1libre, j1libre, k1libre, (std::string(" ") + dir(whatfield)).c_str(), LeftEnd_index, RightEnd_index);
                  if ((k1 >= ZI) && (k1 <= ZE) && control.verbose) WarnErrReport(buff);
               } else if (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_LeftEnd) {
                  sprintf(buff, "wir0_INFO: Ending segment (LeftEnd   )%7d%7d%7d%7d-%7d%7d%7d%s%7d", " ", origIndex, i1, j1, k1, i1libre, j1libre, k1libre, (std::string(" ") + dir(whatfield)).c_str(), LeftEnd_index);
                  if ((k1 >= ZI) && (k1 <= ZE) && control.verbose) WarnErrReport(buff);
               } else if (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].Is_RightEnd) {
                  sprintf(buff, "wir0_INFO: Ending segment (RightEnd   )%7d%7d%7d%7d-%7d%7d%7d%s%7d", " ", origIndex, i1, j1, k1, i1libre, j1libre, k1libre, (std::string(" ") + dir(whatfield)).c_str(), RightEnd_index);
                  if ((k1 >= ZI) && (k1 <= ZE) && control.verbose) WarnErrReport(buff);
               }
               if (sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].IsEnd_norLeft_norRight) {
                  if (control.connectendings) {
                     sprintf(buff, "wir0_INFO: Ending segment (other )%7d%7d%7d%7d-%7d%7d%7d%s", " ", origIndex, i1, j1, k1, i1libre, j1libre, k1libre, (std::string(" ") + dir(whatfield)).c_str());
                     if ((k1 >= ZI) && (k1 <= ZE) && control.verbose) WarnErrReport(buff);
                  } else {
                     sgg.Med(HWires.WireTipoMedio[iwi]).wire[1].segm[iwj].IsEnd_norLeft_norRight = false;

std::ostringstream buff_stream;
                     buff_stream << "wir0_WARNING: Resetting Ending segment (other ) to NON-ENDING" 
                                 << std::setw(7) << origIndex 
                                 << std::setw(7) << i1 
                                 << std::setw(7) << j1 
                                 << std::setw(7) << k1 
                                 << "-" 
                                 << std::setw(7) << i1libre 
                                 << std::setw(7) << j1libre 
                                 << std::setw(7) << k1libre 
                                 << dir(whatfield);
                     std::string buff = buff_stream.str();
                     if ((k1 > ZI) && (k1 <= ZE) && (whatfield != iEz)) WarnErrReport(buff);
                     if ((k1 >= ZI) && (k1 <= ZE) && (whatfield == iEz)) WarnErrReport(buff);
                  }
               }
            }
         }

         // Segment pre-counting including detection of duplicates

         for (int iwi = 1; iwi <= HWires.NumDifferentWires; ++iwi) {
            for (int iwj = 1; sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).numsegmentos; ++iwj) {
               sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwj].repetido = false;
            }
         }
         
         conta = 0;
         for (int iwi = 1; iwi <= HWires.NumDifferentWires; ++iwi) {
            for (int iwj = 1; iwj <= sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).numsegmentos; ++iwj) {
               if (((!sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwj].repetido) || control.strictOLD) && (!sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwj].multirabo)) {
                  i1 = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwj].i;
                  j1 = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwj].j;
                  k1 = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwj].k;
                  whatfield = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwj].ori;
                  
                  for (int iwjjj = iwj + 1; iwjjj <= sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).numsegmentos; ++iwjjj) {
                     i2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwjjj].i;
                     j2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwjjj].j;
                     k2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwjjj].k;
                     whatfield2 = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwjjj].ori;
                     repetido = (i1 == i2) && (j1 == j2) && (k1 == k2) && (whatfield == whatfield2);
                     if (repetido) {
                        if (!(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwjjj].Is_LeftEnd || 
                              sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwjjj].Is_RightEnd)) {

                           sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwjjj].repetido = repetido || 
                           sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwjjj].repetido;
                        } else if (!(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwj].Is_LeftEnd || 
                              sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwj].Is_RightEnd)) {

                           sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwj].repetido = repetido || 
                           sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwj].repetido;
                        } else {
                           // aviso but take a decision. md 260213 sometimes duplicates at start and end!!!!!!

                           if (((abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).Parallel_R_RightEnd) < 1.0e-12_RKIND_wires) && 
                           (abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).Series_R_RightEnd) < 1.0e-12_RKIND_wires) && 
                           (abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).Parallel_C_RightEnd) < 1.0e-12_RKIND_wires) && 
                           (abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).Series_C_RightEnd) > 1.0e7_RKIND_wires) && 
                           (abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).Parallel_L_RightEnd) < 1.0e-12_RKIND_wires) && 
                           (abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).Series_L_RightEnd) < 1.0e-12_RKIND_wires)) && 
                           //
                           ((abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).Parallel_R_LeftEnd) >= 1.0e-12_RKIND_wires) || 
                           (abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).Series_R_LeftEnd) >= 1.0e-12_RKIND_wires) || 
                           (abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).Parallel_C_LeftEnd) >= 1.0e-12_RKIND_wires) || 
                           (abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).Series_C_LeftEnd) <= 1.0e7_RKIND_wires) || 
                           (abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).Parallel_L_LeftEnd) >= 1.0e-12_RKIND_wires) || 
                           (abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).Series_L_LeftEnd) >= 1.0e-12_RKIND_wires))) {
                              if (sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwjjj].Is_RightEnd) {
                                 sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwjjj].repetido = true;
                                 if (control.strictOLD) {
                                    std::ostringstream buff_stream2;
                                    buff_stream2 << "wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Keeping both" 
                                                 << std::setw(7) << sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwj].origindex 
                                                 << std::setw(7) << sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwjjj].origindex 
                                                 << std::setw(7) << i1 
                                                 << std::setw(7) << j1 
                                                 << std::setw(7) << k1;
                                    std::string buff2 = buff_stream2.str();
                                    if ((k1 > ZI) && (k1 <= ZE) && (whatfield != iEz) && control.verbose) WarnErrReport(buff2);
                                    if ((k1 >= ZI) && (k1 <= ZE) && (whatfield == iEz) && control.verbose) WarnErrReport(buff2);
                                 } else {
                                    std::ostringstream buff_stream3;
                                    buff_stream3 << "wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Removing the second one (no RightEnd RLC)" 
                                                 << std::setw(7) << sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwj].origindex 
                                                 << std::setw(7) << sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwjjj].origindex 
                                                 << std::setw(7) << i1 
                                                 << std::setw(7) << j1 
                                                 << std::setw(7) << k1;
                                    std::string buff3 = buff_stream3.str();
                                    if ((k1 > ZI) && (k1 <= ZE) && (whatfield != iEz) && control.verbose) WarnErrReport(buff3);
                                    if ((k1 >= ZI) && (k1 <= ZE) && (whatfield == iEz) && control.verbose) WarnErrReport(buff3);
                                 }
                              } else if (sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwj].Is_RightEnd) {
                                 sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwj].repetido = true;
                                 if (control.strictOLD) {
                                    std::ostringstream buff_stream4;
                                    buff_stream4 << "wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Keeping both" 
                                                 << std::setw(7) << sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwj].origindex 
                                                 << std::setw(7) << sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwjjj].origindex 
                                                 << std::setw(7) << i1 
                                                 << std::setw(7) << j1 
                                                 << std::setw(7) << k1;
                                    std::string buff4 = buff_stream4.str();
                                    if ((k1 > ZI) && (k1 <= ZE) && (whatfield != iEz) && control.verbose) WarnErrReport(buff4);
                                    if ((k1 >= ZI) && (k1 <= ZE) && (whatfield == iEz) && control.verbose) WarnErrReport(buff4);
                                 } else {
                                    std::ostringstream buff_stream5;
                                    buff_stream5 << "wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Removing the first one (no RightEnd RLC)" 
                                                 << std::setw(7) << sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwj].origindex 
                                                 << std::setw(7) << sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm[iwjjj].origindex 
                                                 << std::setw(7) << i1 
                                                 << std::setw(7) << j1 
                                                 << std::setw(7) << k1;
                                    std::string buff5 = buff_stream5.str();
                                    if ((k1 > ZI) && (k1 <= ZE) && (whatfield != iEz) && control.verbose) WarnErrReport(buff5);
                                    if ((k1 >= ZI) && (k1 <= ZE) && (whatfield == iEz) && control.verbose) WarnErrReport(buff5);
                                 }
                              }
                           }

                           if (((abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).Parallel_R_RightEnd) >= 1.0e-12_RKIND_wires) || 
                           (abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).Series_R_RightEnd) >= 1.0e-12_RKIND_wires) || 
                           (abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).Parallel_C_RightEnd) >= 1.0e-12_RKIND_wires) || 
                           (abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).Series_C_RightEnd) <= 1.0e7_RKIND_wires) || 
                           (abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).Parallel_L_RightEnd) >= 1.0e-12_RKIND_wires) || 
                           (abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).Series_L_RightEnd) >= 1.0e-12_RKIND_wires)) ||

((std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Series_L_RightEnd) >= 1.0e-12)).and. &
                           //
                           ((std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Parallel_R_LeftEnd) < 1.0e-12).and. &
                           (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Series_R_LeftEnd) < 1.0e-12).and. &
                           (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Parallel_C_LeftEnd) < 1.0e-12).and. &
                           (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Series_C_LeftEnd) > 1.0e7).and. &
                           (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Parallel_L_LeftEnd) < 1.0e-12).and. &
                           (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Series_L_LeftEnd) < 1.0e-12)) ) then
                              if (sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].Is_LeftEnd) then
                                 sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].repetido = .true.
                                 if (control.strictOLD) then
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Keeping both', &
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].origindex, &
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].origindex,i1,j1,k1
                                 else
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Removing the second one (NO LeftEnd RLC, Only RightEnd RLC)', &
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].origindex, &
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].origindex,i1,j1,k1
                                 end if
                                 if ((k1 >  ZI).and.(k1 <= ZE).and.(whatfield /= iEz).and.control.verbose) call WarnErrReport(buff)
                                 if ((k1 >= ZI).and.(k1 <= ZE).and.(whatfield == iEz).and.control.verbose) call WarnErrReport(buff)
                              elseif (sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].Is_LeftEnd) then
                                 sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].repetido = .true.
                                 if (control.strictOLD) then
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Keeping both', &
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].origindex, &
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].origindex,i1,j1,k1
                                 else
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment. Removing the first one (NO LeftEnd RLC, Only RightEnd RLC)', &
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].origindex, &
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].origindex,i1,j1,k1
                                 end if
                                 if ((k1 >  ZI).and.(k1 <= ZE).and.(whatfield /= iEz).and.control.verbose) call WarnErrReport(buff)
                                 if ((k1 >= ZI).and.(k1 <= ZE).and.(whatfield == iEz).and.control.verbose) call WarnErrReport(buff)
                              end if
                           end if
                           //
                           if ( ((std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Parallel_R_LeftEnd) >= 1.0e-12).or. &
                           (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Series_R_LeftEnd) >= 1.0e-12).or. &
                           (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Parallel_C_LeftEnd) >= 1.0e-12).or. &
                           (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Series_C_LeftEnd) <= 1.0e7).or. &
                           (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Parallel_L_LeftEnd) >= 1.0e-12).or. &
                           (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Series_L_LeftEnd) >= 1.0e-12)).AND. &
                           //
                           ((std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Parallel_R_RightEnd) >= 1.0e-12).or. &
                           (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Series_R_RightEnd) >= 1.0e-12).or. &
                           (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Parallel_C_RightEnd) >= 1.0e-12).or. &
                           (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Series_C_RightEnd) <= 1.0e7).or. &
                           (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Parallel_L_RightEnd) >= 1.0e-12).or. &
                           (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Series_L_RightEnd) >= 1.0e-12)) ) then

                              sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].repetido = repetido.or. &
                              sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].repetido

                              if (  (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Parallel_R_RightEnd -    &
                              sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Parallel_R_LeftEnd) < 1.0e-12).and. &
                              (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Series_R_RightEnd -    &
                              sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Series_R_LeftEnd) < 1.0e-12).and. &
                              (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Parallel_C_RightEnd -    &
                              sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Parallel_C_LeftEnd) < 1.0e-12).and. &
                              (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Series_C_RightEnd -    &
                              sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Series_C_LeftEnd) < 1.0e-12).and. &
                              (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Parallel_L_RightEnd -    &
                              sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Parallel_L_LeftEnd) < 1.0e-12).and. &
                              (std::abs(sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Series_L_RightEnd -    &
                              sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].Series_L_LeftEnd) < 1.0e-12)) then
                                 if (control.strictOLD) then
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment with the same RLC. Keeping both', &
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].origindex, &
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].origindex,i1,j1,k1
                                 else
                                    write (buff,'(a,2i7,3i7)')  'wir0_INFO: Duplicate terminal LeftEnd and RightEnd Parallel segment with the same RLC. Will remove the second one', &
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwj].origindex, &
                                    sgg.Med(HWires.WireTipoMedio(iwi)).wire[1].segm[iwjjj].origindex,i1,j1,k1
                                 end if
                                 if ((k1 >  ZI).and.(k1 <= ZE).and.(whatfield /= iEz).and.control.verbose) call WarnErrReport(buff)

if ((k1 >= ZI) && (k1 <= ZE) && (whatfield == iEz) && control.verbose) WarnErrReport(buff);
                              else {
                                 std::ostringstream oss;
                                 oss << "wir0_ERROR: Duplicate terminal LeftEnd and RightEnd Parallel segment with non-null different RLC"
                                     << std::setw(7) << sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).origindex
                                     << std::setw(7) << sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwjjj).origindex
                                     << std::setw(7) << i1
                                     << std::setw(7) << j1
                                     << std::setw(7) << k1;
                                 std::string buff_str = oss.str();
                                 if ((k1 > ZI) && (k1 <= ZE) && (whatfield != iEz)) WarnErrReport(buff_str, true);
                                 if ((k1 >= ZI) && (k1 <= ZE) && (whatfield == iEz)) WarnErrReport(buff_str, true);
                              }
                           }
                        }
                     }
                  }
               }
            }
            //
            //clipping: in case of direct .nfde reading the PREPROCESSor has not clipped this data
            if ((i1 >= sgg.Alloc(whatfield).XI) &&
                (i1 <= sgg.Alloc(whatfield).XE) &&
                (j1 >= sgg.Alloc(whatfield).YI) &&
                (j1 <= sgg.Alloc(whatfield).YE) &&
                (k1 >= sgg.Alloc(whatfield).ZI) &&
                (k1 <= sgg.Alloc(whatfield).ZE)) {
               conta = conta + 1;
            }
         }
      }
   }
   if (conta == 0) therearewires = false;
   if (!therearewires) return;
   //Report duplicated segments
   for (iwi = 1; iwi <= HWires.NumDifferentWires; ++iwi) {
      for (iwj = 1; iwj <= sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).numsegmentos; ++iwj) {
         i1 = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).i;
         j1 = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).j;
         k1 = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).k;
         whatfield = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).ori;
         origindex = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).origIndex;
         if (sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).repetido && !sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).multirabo) {
            if (control.strictOLD) {
               std::ostringstream oss;
               oss << "wir0_WARNING: Keeping duplicate (Parallel) intra-WIRE segment"
                   << std::setw(7) << origindex
                   << std::setw(7) << i1
                   << std::setw(7) << j1
                   << std::setw(7) << k1
                   << " " << dir(whatfield);
               std::string buff_str = oss.str();
               if ((k1 >= ZI) && (k1 <= ZE) && control.verbose) WarnErrReport(buff_str);
            } else {
               std::ostringstream oss;
               oss << "wir0_WARNING: Removing duplicate (Parallel) intra-WIRE segment and voiding ASSOCIATED probes "
                   << std::setw(7) << origindex
                   << std::setw(7) << i1
                   << std::setw(7) << j1
                   << std::setw(7) << k1
                   << " " << dir(whatfield);
               std::string buff_str = oss.str();
               if ((k1 >= ZI) && (k1 <= ZE) && control.verbose) WarnErrReport(buff_str);
            }
         }
      }
   }

   HWires.NumCurrentSegments = conta;
   //inicializa ctes segmentos
   HWires.CurrentSegment.resize(HWires.NumCurrentSegments + 1);
   for (i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
      HWires.CurrentSegment[i1].ChargePlus = nullptr;
      HWires.CurrentSegment[i1].ChargeMinus = nullptr;
      HWires.CurrentSegment[i1].TipoWire = nullptr;
      HWires.CurrentSegment[i1].Efield_main2wire = nullptr;
      HWires.CurrentSegment[i1].Efield_wire2main = nullptr;
      //
      HWires.CurrentSegment[i1].R = 0.0;
      HWires.CurrentSegment[i1].Resist = 0.0;
      HWires.CurrentSegment[i1].Resist_devia = 0.0;
      HWires.CurrentSegment[i1].C = 0.0;
      HWires.CurrentSegment[i1].L = 0.0;
      HWires.CurrentSegment[i1].Lintrinsic = 0.0;
      HWires.CurrentSegment[i1].proc = false;
      HWires.CurrentSegment[i1].NumParallel = 1;
      HWires.CurrentSegment[i1].origindex = i1;
      HWires.CurrentSegment[i1].indexsegment = i1;
      HWires.CurrentSegment[i1].currentpast = 0.0;
      HWires.CurrentSegment[i1].current = 0.0;
      HWires.CurrentSegment[i1].inv_Lind_acum = 0.0;
      HWires.CurrentSegment[i1].Lind_acum = 0.0;
      HWires.CurrentSegment[i1].HEUR_safety = 0.0;
      HWires.CurrentSegment[i1].Lind = 0.0;
      HWires.CurrentSegment[i1].Lind_devia = 0.0;
      //!!! HWires.CurrentSegment[i1].logRoverR0 = 0.0;
      HWires.CurrentSegment[i1].delta = 0.0;
      HWires.CurrentSegment[i1].deltaTransv1 = 0.0;
      HWires.CurrentSegment[i1].deltaTransv2 = 0.0;
      HWires.CurrentSegment[i1].cte1 = 0.0;
      HWires.CurrentSegment[i1].cte2 = 0.0;
      HWires.CurrentSegment[i1].cte3 = 0.0;
      HWires.CurrentSegment[i1].cte1_for_devia = 0.0;
      HWires.CurrentSegment[i1].cte2_for_devia = 0.0;
      HWires.CurrentSegment[i1].cte3_for_devia = 0.0;
      HWires.CurrentSegment[i1].cte5 = 0.0;
      HWires.CurrentSegment[i1].ilibre = -1;
      HWires.CurrentSegment[i1].jlibre = -1;
      HWires.CurrentSegment[i1].klibre = -1;
      HWires.CurrentSegment[i1].i = -1;
      HWires.CurrentSegment[i1].j = -1;
      HWires.CurrentSegment[i1].k = -1;
      HWires.CurrentSegment[i1].tipofield = -1;
      HWires.CurrentSegment[i1].IsPMC = false;
      HWires.CurrentSegment[i1].orientadoalreves = false;
      HWires.CurrentSegment[i1].HasVsource = false;
      HWires.CurrentSegment[i1].IsShielded = false;
      HWires.CurrentSegment[i1].HasAbsorbing_RightEnd = false;
      HWires.CurrentSegment[i1].HasAbsorbing_LeftEnd = false;
      HWires.CurrentSegment[i1].HasParallel_RightEnd = false;
      HWires.CurrentSegment[i1].HasParallel_LeftEnd = false;
      HWires.CurrentSegment[i1].HasSeries_RightEnd = false;
      HWires.CurrentSegment[i1].HasSeries_LeftEnd = false;
      HWires.CurrentSegment[i1].IsEnd_norLeft_norRight = false;
      HWires.CurrentSegment[i1].Is_LeftEnd = false;
      HWires.CurrentSegment[i1].Is_RightEnd = false;
   }

   //assign segment info
   conta = 0;
   for (iwi = 1; iwi <= HWires.NumDifferentWires; ++iwi) {
      for (iwj = 1; iwj <= sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).numsegmentos; ++iwj) {
         i1libre = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).ilibre;
         j1libre = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).jlibre;
         k1libre = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).klibre;
         i1 = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).i;
         j1 = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).j;
         k1 = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).k;
         whatfield = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).ori;
         origindex = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).origIndex;
         IsEnd_norLeft_norRight = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).IsEnd_norLeft_norRight;
         Is_LeftEnd = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).Is_LeftEnd;
         Is_RightEnd = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).Is_RightEnd;
         if ((!sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).repetido || control.strictOLD) && !sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj).multirabo) {

// clipping: in case of direct .nfde reading the PREPROCESSor has not clipped this data
                  if ((i1 >= sgg.Alloc(whatfield).XI).and. &
                  (i1 <= sgg.Alloc(whatfield).XE).and. &
                  (j1 >= sgg.Alloc(whatfield).YI).and. &
                  (j1 <= sgg.Alloc(whatfield).YE).and. &
                  (k1 >= sgg.Alloc(whatfield).ZI).and. &
                  (k1 <= sgg.Alloc(whatfield).ZE)) {
                     conta = conta + 1;
                     HWires.CurrentSegment(conta).IsEnd_norLeft_norRight = IsEnd_norLeft_norRight;
                     HWires.CurrentSegment(conta).Is_LeftEnd = Is_LeftEnd;
                     HWires.CurrentSegment(conta).Is_RightEnd = Is_RightEnd;
                     HWires.CurrentSegment(conta).origindex = origindex;
                     HWires.CurrentSegment(conta).tipofield = whatfield;
                     HWires.CurrentSegment(conta).ilibre = i1libre;
                     HWires.CurrentSegment(conta).jlibre = j1libre;
                     HWires.CurrentSegment(conta).klibre = k1libre;
                     HWires.CurrentSegment(conta).i = i1;
                     HWires.CurrentSegment(conta).j = j1;
                     HWires.CurrentSegment(conta).k = k1;
                     HWires.CurrentSegment(conta).ie = i1;
                     HWires.CurrentSegment(conta).je = j1;
                     HWires.CurrentSegment(conta).ke = k1;
                     HWires.CurrentSegment(conta).indexmed = HWires.WireTipoMedio(iwi);
                     HWires.CurrentSegment(conta).TipoWire = sgg.Med(HWires.WireTipoMedio(iwi)).wire(1);
                     // only for the observation sign to match (not used in this routine)
                     if (iwj < sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).numsegmentos) {
                        if (!control.strictOLD) {
                           HWires.CurrentSegment(conta).orientadoalreves = &
                           (i1 > sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj+1).i).or. &
                           (j1 > sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj+1).j).or. &
                           (k1 > sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj+1).k);
                        } else {
                           HWires.CurrentSegment(conta).orientadoalreves = false; // later corrected
                        }
                     } else if (iwj > 1) {
                        // only for the observation sign to match (not used in this routine)
                        if (!control.strictOLD) {
                           HWires.CurrentSegment(conta).orientadoalreves = &
                           (i1 < sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj-1).i).or. &
                           (j1 < sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj-1).j).or. &
                           (k1 < sgg.Med(HWires.WireTipoMedio(iwi)).wire(1).segm(iwj-1).k);
                        } else {
                           HWires.CurrentSegment(conta).orientadoalreves = false; // later corrected
                        }
                     }
                     //
                     switch (HWires.CurrentSegment(conta).tipofield) {
                     case iEx:
                        // default
                        HWires.CurrentSegment(conta).Efield_wire2main = Ex(i1,j1,k1); 
                        HWires.CurrentSegment(conta).Efield_main2wire = Ex(i1,j1,k1); 
                        // al final se reapuntan los del thickness
                        HWires.CurrentSegment(conta).delta = 1.0_RKIND_wires / Idxe(i1);    // ojo esto de los delta habra que corregirlo para uniones
                        HWires.CurrentSegment(conta).deltaTransv1 = 1.0_RKIND_wires / Idyh(j1);
                        if (k1 <= sgg.ALLOC(iEz).ZE) { // esta corriente en el limite de los alloc nunca se precisa
                           // updateo por exceso para que haya una celda de solapamiento y la union topologica de
                           // hilos no se pierda por culpa de la puta particion MPI
                           HWires.CurrentSegment(conta).deltaTransv2 = 1.0_RKIND_wires / Idzh(k1);
                        } else { // no se comete error alguno
                           // es solo para que el indice no reviente 2012
                           HWires.CurrentSegment(conta).deltaTransv2 = 1.0_RKIND_wires / Idzh(k1-1);
                        }
                        // dama
                        HWires.CurrentSegment(conta).ie = i1+1;
                        HWires.CurrentSegment(conta).x = j1+0.25_RKIND_wires;
                        HWires.CurrentSegment(conta).y = k1+0.25_RKIND_wires; 
                        sggmiE = sggmiEx(i1,j1,k1); deembed_peclossyconformal_segments(sggmiE); sggmiEx(i1,j1,k1) = sggmiE; // por si se ha modificado !ojo agresivo 180220 !ojo cambiado aqui 170323 de sitio pq no se habian puesto los deltatrans y salia division por cero
                        // fin dama
                        break;
                     case iEy:    
                        // default
                        HWires.CurrentSegment(conta).Efield_wire2main = Ey(i1,j1,k1); 
                        HWires.CurrentSegment(conta).Efield_main2wire = Ey(i1,j1,k1);      
                        //
                        HWires.CurrentSegment(conta).delta = 1.0_RKIND_wires / Idye(j1);
                        if (k1 <= sgg.ALLOC(iEz).ZE) { // esta corriente en el limite de los alloc nunca se precisa
                           // updateo por exceso para que haya una celda de solapamiento y la union topologica de
                           // hilos no se pierda por culpa de la puta particion MPI
                           HWires.CurrentSegment(conta).deltaTransv1 = 1.0_RKIND_wires / Idzh(k1);
                        } else {
                           HWires.CurrentSegment(conta).deltaTransv1 = 1.0_RKIND_wires / Idzh(k1-1); // no se comete error alguno
                           // es solo para que el indice no reviente 2012
                        }
                        HWires.CurrentSegment(conta).deltaTransv2 = 1.0_RKIND_wires / Idxh(i1);
                        // dama
                        HWires.CurrentSegment(conta).je = j1+1;
                        HWires.CurrentSegment(conta).x = k1+0.25_RKIND_wires;
                        HWires.CurrentSegment(conta).y = i1+0.25_RKIND_wires;
                        sggmiE = sggmiEy(i1,j1,k1); deembed_peclossyconformal_segments(sggmiE); sggmiEy(i1,j1,k1) = sggmiE; // por si se ha modificado !ojo agresivo 180220 !ojo cambiado aqui 170323 de sitio pq no se habian puesto los deltatrans y salia division por cero
                        // fin dama
                        break;
                     case iEz:  
                        // default
                        HWires.CurrentSegment(conta).Efield_wire2main = Ez(i1,j1,k1); 
                        HWires.CurrentSegment(conta).Efield_main2wire = Ez(i1,j1,k1);  
                        //
                        HWires.CurrentSegment(conta).delta = 1.0_RKIND_wires / Idze(k1);
                        HWires.CurrentSegment(conta).deltaTransv1 = 1.0_RKIND_wires / Idxh(i1);
                        HWires.CurrentSegment(conta).deltaTransv2 = 1.0_RKIND_wires / Idyh(j1);
                        // dama
                        HWires.CurrentSegment(conta).ke = k1+1;
                        HWires.CurrentSegment(conta).x = i1+0.25_RKIND_wires;
                        HWires.CurrentSegment(conta).y = j1+0.25_RKIND_wires;
                        sggmiE = sggmiEz(i1,j1,k1); deembed_peclossyconformal_segments(sggmiE); sggmiEz(i1,j1,k1) = sggmiE; // por si se ha modificado !ojo agresivo 180220 !ojo cambiado aqui 170323 de sitio pq no se habian puesto los deltatrans y salia division por cero
                        // fin dama
                        break;
                     }
                  }
               } // del repetido
            }
         }
 
        // fin niapa 171216
      //
      HWires.NumCurrentSegments = conta;
      //

      // hacer agujeros

      for (i1 = 1; i1 <= HWires.NumCurrentSegments; i1++) {
         segmento = HWires.CurrentSegment(i1);
         i = segmento.i;
         j = segmento.j;
         k = segmento.k;

whatfield = segmento.tipofield;
            IsEnd_norLeft_norRight = segmento.IsEnd_norLeft_norRight;
            Is_LeftEnd = segmento.Is_LeftEnd;
            Is_RightEnd = segmento.Is_RightEnd;
            if ((i > SINPML_fullsize(whatfield).XI) &&
                (i < SINPML_fullsize(whatfield).XE) &&
                (j > SINPML_fullsize(whatfield).YI) &&
                (j < SINPML_fullsize(whatfield).YE) &&
                (k > SINPML_fullsize(whatfield).ZI) &&
                (k < SINPML_fullsize(whatfield).ZE)) {
                if (control.makeholes && (!IsEnd_norLeft_norRight) && (!Is_LeftEnd) && (!Is_RightEnd)) {
                    if (control.num_procs == 0) {
                        stoponerror(control.layoutnumber, control.num_procs, "Makeholes not available for MPI. Stoppping. ");
                    }
                    switch (whatfield) {
                        case iEx:
                            sggmiHx(i, j, k) = 1;
                            sggmiHx(i, j - 1, k) = 1;
                            sggmiHx(i, j, k - 1) = 1;
                            sggmiHx(i, j - 1, k - 1) = 1;
                            sggmiHx(i + 1, j, k) = 1;
                            sggmiHx(i + 1, j - 1, k) = 1;
                            sggmiHx(i + 1, j, k - 1) = 1;
                            sggmiHx(i + 1, j - 1, k - 1) = 1;
                            if (!sgg.med(sggmiEx(i, j, k)).Is.ThinWire) {
                                sggmiEx(i, j, k) = 1;
                                write(buff, "(a,3i7,a,i7,a)", "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ", i, j, k, " for WIRE-segment ", segmento.origIndex, " ", dir(whatfield));
                                WarnErrReport(buff);
                            }
                            if (!sgg.med(sggmiEy(i, j, k)).Is.ThinWire) {
                                sggmiEy(i, j, k) = 1;
                                write(buff, "(a,3i7,a,i7,a)", "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ", i, j, k, " for WIRE-segment ", segmento.origIndex, " ", dir(whatfield));
                                WarnErrReport(buff);
                            }
                            if (!sgg.med(sggmiEy(i, j - 1, k)).Is.ThinWire) {
                                sggmiEy(i, j - 1, k) = 1;
                                write(buff, "(a,3i7,a,i7,a)", "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ", i, j, k, " for WIRE-segment ", segmento.origIndex, " ", dir(whatfield));
                                WarnErrReport(buff);
                            }
                            if ((k <= sgg.alloc(iEz).ZE) && (k >= sgg.alloc(iEz).ZI)) {
                                if (!sgg.med(sggmiEz(i, j, k)).Is.ThinWire) {
                                    sggmiEz(i, j, k) = 1;
                                    write(buff, "(a,3i7,a,i7,a)", "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ", i, j, k, " for WIRE-segment ", segmento.origIndex, " ", dir(whatfield));
                                    WarnErrReport(buff);
                                }
                            }
                            if ((k - 1 <= sgg.alloc(iEz).ZE) && (k - 1 >= sgg.alloc(iEz).ZI)) {
                                if (!sgg.med(sggmiEz(i, j, k - 1)).Is.ThinWire) {
                                    sggmiEz(i, j, k - 1) = 1;
                                    write(buff, "(a,3i7,a,i7,a)", "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ", i, j, k, " for WIRE-segment ", segmento.origIndex, " ", dir(whatfield));
                                    WarnErrReport(buff);
                                }
                            }
                            if (!sgg.med(sggmiEy(i + 1, j, k)).Is.ThinWire) {
                                sggmiEy(i + 1, j, k) = 1;
                                write(buff, "(a,3i7,a,i7,a)", "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ", i, j, k, " for WIRE-segment ", segmento.origIndex, " ", dir(whatfield));
                                WarnErrReport(buff);
                            }
                            if (!sgg.med(sggmiEy(i + 1, j - 1, k)).Is.ThinWire) {
                                sggmiEy(i + 1, j - 1, k) = 1;
                                write(buff, "(a,3i7,a,i7,a)", "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ", i, j, k, " for WIRE-segment ", segmento.origIndex, " ", dir(whatfield));
                                WarnErrReport(buff);
                            }
                            if ((k <= sgg.alloc(iEz).ZE) && (k >= sgg.alloc(iEz).ZI)) {
                                if (!sgg.med(sggmiEz(i + 1, j, k)).Is.ThinWire) {
                                    sggmiEz(i + 1, j, k) = 1;
                                    write(buff, "(a,3i7,a,i7,a)", "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ", i, j, k, " for WIRE-segment ", segmento.origIndex, " ", dir(whatfield));
                                    WarnErrReport(buff);
                                }
                            }
                            if ((k - 1 <= sgg.alloc(iEz).ZE) && (k - 1 >= sgg.alloc(iEz).ZI)) {
                                if (!sgg.med(sggmiEz(i + 1, j, k - 1)).Is.ThinWire) {
                                    sggmiEz(i + 1, j, k - 1) = 1;
                                    write(buff, "(a,3i7,a,i7,a)", "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ", i, j, k, " for WIRE-segment ", segmento.origIndex, " ", dir(whatfield));
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
                                write(buff, "(a,3i7,a,i7,a)", "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ", i, j, k, " for WIRE-segment ", segmento.origIndex, " ", dir(whatfield));
                                WarnErrReport(buff);
                            }
                            if ((k <= sgg.alloc(iEz).ZE) && (k >= sgg.alloc(iEz).ZI)) {
                                if (!sgg.med(sggmiEz(i, j, k)).Is.ThinWire) {
                                    sggmiEz(i, j, k) = 1;
                                    write(buff, "(a,3i7,a,i7,a)", "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ", i, j, k, " for WIRE-segment ", segmento.origIndex, " ", dir(whatfield));
                                    WarnErrReport(buff);
                                }
                            }
                            if ((k - 1 <= sgg.alloc(iEz).ZE) && (k - 1 >= sgg.alloc(iEz).ZI)) {
                                if (!sgg.med(sggmiEz(i, j, k - 1)).Is.ThinWire) {
                                    sggmiEz(i, j, k - 1) = 1;
                                    write(buff, "(a,3i7,a,i7,a)", "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ", i, j, k, " for WIRE-segment ", segmento.origIndex, " ", dir(whatfield));
                                    WarnErrReport(buff);
                                }
                            }
                            if (!sgg.med(sggmiEx(i, j, k)).Is.ThinWire) {
                                sggmiEx(i, j, k) = 1;
                                write(buff, "(a,3i7,a,i7,a)", "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ", i, j, k, " for WIRE-segment ", segmento.origIndex, " ", dir(whatfield));
                                WarnErrReport(buff);
                            }
                            break;
                        default:
                            break;
                    }
                }
            }

// Note: The first line in the Fortran snippet is a continuation of a previous write statement.
            // It is assumed 'buff' is a std::string and 'dir' is a function or array access returning a string.
            // 'segmento' is likely a struct/class instance.
            // 'whatfield' is likely an integer or enum.
            // 'WarnErrReport' is a function.
            
            // Reconstructing the write statement from the previous context (implied):
            // write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
            // Since we only translate the provided chunk, we assume the previous write completed or this is a continuation.
            // However, the Fortran code shows:
            // i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
            // call WarnErrReport(buff)
            // This implies the previous write statement was split or this is a new write. 
            // Looking at the pattern, it seems like a write statement was started in a previous block (not shown) and continued here.
            // But wait, the line starts with 'i,j,k,...'. This is likely a continuation of a WRITE statement.
            // Let's assume the previous code handled the start of the string.
            // To be safe and preserve logic, I will assume 'buff' is already populated with the prefix.
            // Actually, looking at the indentation and structure, it's highly likely this is part of a WRITE statement.
            // Since I cannot see the start, I will assume the previous context filled 'buff' up to the point before 'i,j,k'.
            // Or, more likely, this is a standalone line in a WRITE statement.
            // Let's look at the next lines. They are full WRITE statements.
            // The first line is likely:
            // write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
            // But the Fortran code provided starts with 'i,j,k...'.
            // I will assume the previous block (not shown) did the first part.
            // However, to make this code compile and run correctly in isolation if needed, I should probably reconstruct the full write.
            // But the prompt says "Translate ONLY this part".
            // I will assume 'buff' is a std::string and append to it.
            
            // Re-reading the Fortran:
            // i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
            // call WarnErrReport(buff)
            // This is syntactically a continuation of a WRITE statement.
            // I will assume the previous code segment (not provided) started the WRITE.
            // To make it valid C++, I'll assume 'buff' is a std::string and we are appending.
            // But without the start, I can't know the format.
            // Let's look at the other WRITE statements. They use '(a,3i7,a,i7,a)'.
            // The arguments are: 'wir0_WARNING...', i, j, k, ' for WIRE-segment ', segmento%origIndex, ' '//dir(whatfield).
            // So the first line is likely the end of such a statement.
            // I will assume the previous code (not shown) did:
            // std::string buff = "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ";
            // And then formatted i, j, k.
            // Since I can't see it, I will just translate the visible part as an append to buff, assuming buff is already partially filled.
            // This is risky. 
            // Alternative: The Fortran code might be:
            // write (buff,'(a,3i7,a,i7,a)')  'wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ', i,j,k,' for WIRE-segment ',segmento%origIndex,' '//dir(whatfield)
            // And the line break is just for readability in the source.
            // If so, the first line 'i,j,k...' is just the continuation.
            // I will assume the previous context (not shown) handled the first part.
            // To be safe, I will write a comment and assume the user has the context.
            // However, I must output C++ code.
            // I will assume 'buff' is a std::string and append the rest.
            
            // Let's assume the previous code did:
            // buff = "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at ";
            // Then formatted i, j, k.
            // I will just append the rest.
            
            // Actually, looking at the indentation, it's inside an IF block.
            // I will assume 'buff' is a std::string.
            
            // Append the rest of the message
            buff += " for WIRE-segment " + std::to_string(segmento.origIndex) + dir(whatfield);
            WarnErrReport(buff);
            
            if (!sgg.med(sggmiEx(i-1,j,k))->Is.ThinWire) {
                sggmiEx(i-1,j,k) = 1;
                buff = "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at " + std::to_string(i) + std::to_string(j) + std::to_string(k) + " for WIRE-segment " + std::to_string(segmento.origIndex) + dir(whatfield);
                WarnErrReport(buff);
            }
            if ((k <= sgg.alloc(iEz).ZE) && (k >= sgg.alloc(iEz).ZI)) {
                if (!sgg.med(sggmiEz(i,j+1,k))->Is.ThinWire) {
                    sggmiEz(i,j+1,k) = 1;
                    buff = "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at " + std::to_string(i) + std::to_string(j) + std::to_string(k) + " for WIRE-segment " + std::to_string(segmento.origIndex) + dir(whatfield);
                    WarnErrReport(buff);
                }
            }
            if ((k-1 <= sgg.alloc(iEz).ZE) && (k-1 >= sgg.alloc(iEz).ZI)) {
                if (!sgg.med(sggmiEz(i,j+1,k-1))->Is.ThinWire) {
                    sggmiEz(i,j+1,k-1) = 1;
                    buff = "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at " + std::to_string(i) + std::to_string(j) + std::to_string(k) + " for WIRE-segment " + std::to_string(segmento.origIndex) + dir(whatfield);
                    WarnErrReport(buff);
                }
            }
            if (!sgg.med(sggmiEx(i,j+1,k))->Is.ThinWire) {
                sggmiEx(i,j+1,k) = 1;
                buff = "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at " + std::to_string(i) + std::to_string(j) + std::to_string(k) + " for WIRE-segment " + std::to_string(segmento.origIndex) + dir(whatfield);
                WarnErrReport(buff);
            }
            if (!sgg.med(sggmiEx(i-1,j+1,k))->Is.ThinWire) {
                sggmiEx(i-1,j+1,k) = 1;
                buff = "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at " + std::to_string(i) + std::to_string(j) + std::to_string(k) + " for WIRE-segment " + std::to_string(segmento.origIndex) + dir(whatfield);
                WarnErrReport(buff);
            }
            case (iEz):
                sggmiHz(i,j,k) = 1;
                sggmiHz(i-1,j,k) = 1;
                sggmiHz(i,j-1,k) = 1;
                sggmiHz(i-1,j-1,k) = 1;
                sggmiHz(i,j,k+1) = 1;
                sggmiHz(i-1,j,k+1) = 1;
                sggmiHz(i,j-1,k+1) = 1;
                sggmiHz(i-1,j-1,k+1) = 1;
                if (!sgg.med(sggmiEz(i,j,k))->Is.ThinWire) {
                    sggmiEz(i,j,k) = 1;
                    buff = "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at " + std::to_string(i) + std::to_string(j) + std::to_string(k) + " for WIRE-segment " + std::to_string(segmento.origIndex) + dir(whatfield);
                    WarnErrReport(buff);
                }
                if ((k <= sgg.alloc(iEx).ZE) && (k >= sgg.alloc(iEx).ZI)) {
                    if (!sgg.med(sggmiEx(i,j,k))->Is.ThinWire) {
                        sggmiEx(i,j,k) = 1;
                        buff = "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at " + std::to_string(i) + std::to_string(j) + std::to_string(k) + " for WIRE-segment " + std::to_string(segmento.origIndex) + dir(whatfield);
                        WarnErrReport(buff);
                    }
                    if (!sgg.med(sggmiEx(i-1,j,k))->Is.ThinWire) {
                        sggmiEx(i-1,j,k) = 1;
                        buff = "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at " + std::to_string(i) + std::to_string(j) + std::to_string(k) + " for WIRE-segment " + std::to_string(segmento.origIndex) + dir(whatfield);
                        WarnErrReport(buff);
                    }
                }
                if ((k <= sgg.alloc(iEy).ZE) && (k >= sgg.alloc(iEy).ZI)) {
                    if (!sgg.med(sggmiEy(i,j,k))->Is.ThinWire) {
                        sggmiEy(i,j,k) = 1;
                        buff = "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at " + std::to_string(i) + std::to_string(j) + std::to_string(k) + " for WIRE-segment " + std::to_string(segmento.origIndex) + dir(whatfield);
                        WarnErrReport(buff);
                    }
                    if (!sgg.med(sggmiEy(i,j-1,k))->Is.ThinWire) {
                        sggmiEy(i,j-1,k) = 1;
                        buff = "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at " + std::to_string(i) + std::to_string(j) + std::to_string(k) + " for WIRE-segment " + std::to_string(segmento.origIndex) + dir(whatfield);
                        WarnErrReport(buff);
                    }
                }
                if ((k+1 <= sgg.alloc(iEx).ZE) && (k+1 >= sgg.alloc(iEx).ZI)) {
                    if (!sgg.med(sggmiEx(i,j,k+1))->Is.ThinWire) {
                        sggmiEx(i,j,k+1) = 1;
                        buff = "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at " + std::to_string(i) + std::to_string(j) + std::to_string(k) + " for WIRE-segment " + std::to_string(segmento.origIndex) + dir(whatfield);
                        WarnErrReport(buff);
                    }
                    if (!sgg.med(sggmiEx(i-1,j,k+1))->Is.ThinWire) {
                        sggmiEx(i-1,j,k+1) = 1;
                        buff = "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at " + std::to_string(i) + std::to_string(j) + std::to_string(k) + " for WIRE-segment " + std::to_string(segmento.origIndex) + dir(whatfield);
                        WarnErrReport(buff);
                    }
                }
                if ((k+1 <= sgg.alloc(iEy).ZE) && (k+1 >= sgg.alloc(iEy).ZI)) {
                    if (!sgg.med(sggmiEy(i,j,k+1))->Is.ThinWire) {
                        sggmiEy(i,j,k+1) = 1;
                        buff = "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at " + std::to_string(i) + std::to_string(j) + std::to_string(k) + " for WIRE-segment " + std::to_string(segmento.origIndex) + dir(whatfield);
                        WarnErrReport(buff);
                    }
                    if (!sgg.med(sggmiEy(i,j-1,k+1))->Is.ThinWire) {
                        sggmiEy(i,j-1,k+1) = 1;
                        buff = "wir0_WARNING: Making a two-cell free-space thru-hole  (take care of possible open air leftovers) at " + std::to_string(i) + std::to_string(j) + std::to_string(k) + " for WIRE-segment " + std::to_string(segmento.origIndex) + dir(whatfield);
                        WarnErrReport(buff);
                    }
                }
                break;
            }
        }
    }
}

// !!!calculo y gestion autoinduccion
for (int i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
    jmed = HWires.CurrentSegment[i1].indexmed;
    desp = HWires.CurrentSegment[i1].delta;
    despT1 = HWires.CurrentSegment[i1].deltaTransv1;
    despT2 = HWires.CurrentSegment[i1].deltaTransv2;
    // Esto no debe ser preciso. la hipotesis es que vuelve promediando a una celda
}

if (control.wirethickness > 1) {
            despT1 = control.wirethickness * despT1;
            despT2 = control.wirethickness * despT2;
        }
        r0 = HWires.CurrentSegment[i1].TipoWire.Radius;

        if (r0 < 1e-30) {
            std::cout << "wir0_ERROR: ire radius cannot be null" << std::endl;
            if ((k1 >= ZI) && (k1 <= ZE)) {
                WarnErrReport(buff, true);
            }
        }

        if (r0 < 1e-9 * desp) {
            std::cout << "wir0_WARNING: WIRE radius too small " << r0 << std::endl;
            WarnErrReport(buff);
        }

        if ((r0 > 0.5 * despT1) || (r0 > 0.5 * despT2)) {
            std::cout << "wir0_WARNING: WIRE radius greater than half a space-step. Reduced to this limit " 
                      << std::min(0.5 * despT1, 0.5 * despT2) << std::endl;
            if ((HWires.CurrentSegment[i1].k >= ZI) && (HWires.CurrentSegment[i1].k <= ZE)) {
                WarnErrReport(buff);
            }
        }

        if (control.inductance_model == "berenger") {
            HWires.CurrentSegment[i1].Lind = 
                (1.0 / (4.0 * pi * InvMu[jmed])) * 
                (std::log((despT1 * despT1 + despT2 * despT2) / (4.0 * r0 * r0)) + 
                 despT1 / despT2 * std::atan(despT2 / despT1) + 
                 despT2 / despT1 * std::atan(despT1 / despT2) + 
                 pi * r0 * r0 / (despT2 * despT1) - 3.0);

            if ((r0 > 0.3 * despT1) || (r0 > 0.3 * despT2)) {
                HWires.CurrentSegment[i1].Lind = HWires.CurrentSegment[i1].Lind / 
                    (1.0 - pi * r0 * r0 / (despT1 * despT2));
            }
        } else if (control.inductance_model == "ledfelt") {
            HWires.CurrentSegment[i1].Lind = 
                (1.0 / (4.0 * pi * InvMu[jmed])) * 
                (std::log((despT1 * despT1 + despT2 * despT2) / (r0 * r0)) + 
                 despT1 / despT2 * std::atan(despT2 / despT1) + 
                 despT2 / despT1 * std::atan(despT1 / despT2) + 
                 pi * r0 * r0 / (16.0 * despT2 * despT1) - 3.0);

            if ((r0 > 0.3 * despT1) || (r0 > 0.3 * despT2)) {
                HWires.CurrentSegment[i1].Lind = HWires.CurrentSegment[i1].Lind / 
                    (1.0 - pi * r0 * r0 / (despT1 * despT2));
            }
        } else if (control.inductance_model == "boutayeb") {
            HWires.CurrentSegment[i1].Lind = 
                (1.0 / (4.0 * pi * InvMu[jmed])) * 
                (std::log((despT1 * despT1 + despT2 * despT2) / (4.0 * r0 * r0)) + 
                 despT1 / despT2 * std::atan(despT2 / despT1) + 
                 despT2 / despT1 * std::atan(despT1 / despT2) + 
                 pi * r0 * r0 / (despT2 * despT1) - 3.0);

            if ((r0 < 0.3 * despT1) || (r0 < 0.3 * despT2)) {
                HWires.CurrentSegment[i1].Lind = HWires.CurrentSegment[i1].Lind - 
                    0.57 / (4.0 * pi * InvMu[jmed]);
            }

            if ((r0 > 0.3 * despT1) || (r0 > 0.3 * despT2)) {
                HWires.CurrentSegment[i1].Lind = HWires.CurrentSegment[i1].Lind / 
                    (1.0 - pi * r0 * r0 / (despT1 * despT2));
            }
        } else {
            buff = "wir0_ERROR: Incorrect inductance model";
            WarnErrReport(buff, true);
        }

        if (HWires.CurrentSegment[i1].Lind < 0.0) {
            buff = "wir0_ERROR: Wrong self-inductance. ";
            WarnErrReport(buff, true);
        }

    } // end do i1

    LindProb.resize(HWires.NumCurrentSegments);
    for (i1 = 0; i1 < HWires.NumCurrentSegments; ++i1) {
        LindProb[i1] = true;
        HWires.CurrentSegment[i1].inv_Lind_acum = 1.0 / HWires.CurrentSegment[i1].Lind;
        HWires.CurrentSegment[i1].HEUR_safety = (sgg.dt * sgg.dt) / (eps0 * HWires.CurrentSegment[i1].deltaTransv1 * HWires.CurrentSegment[i1].deltaTransv2);
    }

    for (i1 = 0; i1 < HWires.NumCurrentSegments; ++i1) {
        if (LindProb[i1]) {
            org = &HWires.CurrentSegment[i1];
            org->NumParallel = 1;
            org->Lind_acum = org->Lind;
            
            for (j1 = i1 + 1; j1 < HWires.NumCurrentSegments; ++j1) {
                fin = &HWires.CurrentSegment[j1];
                if ((org->i == fin->i) && (org->j == fin->j) && (org->k == fin->k) && (org->tipofield == fin->tipofield)) {
                    org->NumParallel = org->NumParallel + 1;
                    if (control.stableradholland) {
                        org->Lind_acum = org->Lind_acum + fin->Lind;
                    }
                }
            }

            for (j1 = i1 + 1; j1 < HWires.NumCurrentSegments; ++j1) {
                fin = &HWires.CurrentSegment[j1];
                if ((org->i == fin->i) && (org->j == fin->j) && (org->k == fin->k) && (org->tipofield == fin->tipofield)) {
                    fin->NumParallel = org->NumParallel;
                    fin->Lind_acum = org->Lind_acum;
                    LindProb[j1] = false;
                }
            }
        }
    }

// !!!!!!!!!!!!!!!!!\E7\E7\E7\E7\E7\E7\E7\E7\E7\E7
        // !!!!!!!! hago 120715 mi criterior
        dtcritico = sgg.dt;
        for (int i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
            auto& dummy = HWires.CurrentSegment[i1];
            int jmed = dummy.indexmed;
            double desp = dummy.delta;
            double despT1 = dummy.deltaTransv1;
            double despT2 = dummy.deltaTransv2;
            double r0 = dummy.TipoWire->Radius;
            if (r0 < 1e-30) {
                std::ostringstream buff_stream;
                buff_stream << "wir0_ERROR: wire radius cannot be null";
                std::string buff = buff_stream.str();
                if ((k1 >= ZI) && (k1 <= ZE)) {
                    WarnErrReport(buff, true);
                }
            }
            // !!!!!!!!!!!!correccion bruta 13/07/15
            double deltadummy = dummy.inv_Lind_acum * dummy.HEUR_safety / 0.9; // EL 0.9 ES POR MI TRANQUILIDAD
            if (deltadummy > 1.0) {
                if (control.stableradholland) {
                    if (adjustl(control.inductance_model) == "boutayeb") {
                        double b = -4.0 * M_PI * InvMu(jmed) * (dummy.Lind * deltadummy) + 
                                   log(despT1 * despT1 + despT2 * despT2) + 
                                   (despT1 / despT2) * atan(despT2 / despT1) + 
                                   (despT2 / despT1) * atan(despT1 / despT2) - 3.0;
                        double a = M_PI / (despT2 * despT1);
                        if ((r0 < 0.3 * despT1) || (r0 < 0.3 * despT2)) {
                            b = b - 0.57;
                        }
                        double newr0 = sqrt(-Lambert(-a * exp(b) / 4.0) / a);
                        // !!!doublechecking
                        b = -4.0 * M_PI * InvMu(jmed) * (dummy.Lind) + 
                            log(despT1 * despT1 + despT2 * despT2) + 
                            (despT1 / despT2) * atan(despT2 / despT1) + 
                            (despT2 / despT1) * atan(despT1 / despT2) - 3.0;
                        a = M_PI / (despT2 * despT1);
                        if ((r0 < 0.3 * despT1) || (r0 < 0.3 * despT2)) {
                            b = b - 0.57;
                        }
                        double OLDR0 = sqrt(-Lambert(-a * exp(b) / 4.0) / a);
                        // !!!!!!!!!
                        std::ostringstream buff_stream2;
                        buff_stream2 << "wir0_WARNING: AUTOMATIC CORRECTION OF L/mu0=" << dummy.Lind / mu0 
                                     << " for r0=" << r0 << oldr0 << " " 
                                     << dummy.NumParallel << " wires at " 
                                     << dummy.i << dummy.j << dummy.k 
                                     << " to L/mu0=" << dummy.Lind * deltadummy / mu0 
                                     << " for newr0=" << newr0;
                        std::string buff = buff_stream2.str();
                    } else {
                        std::ostringstream buff_stream3;
                        buff_stream3 << "wir0_WARNING: AUTOMATIC CORRECTION OF L/mu0=" << dummy.Lind / mu0 
                                     << " " << dummy.NumParallel << " wires at " 
                                     << dummy.i << dummy.j << dummy.k 
                                     << " to L/mu0=" << dummy.Lind * deltadummy / mu0;
                        std::string buff = buff_stream3.str();
                    }
                    if ((dummy.k > ZI) && (dummy.k <= ZE)) {
                        WarnErrReport(buff);
                    }
                    dummy.Lind = dummy.Lind * deltadummy; // bajo repartiendo proporcialmente
                } else {
                    std::ostringstream buff_stream4;
                    buff_stream4 << "wir0_SEVEREWARNING: L/mu0=" << dummy.Lind / mu0 
                                 << " in " << dummy.NumParallel << " wires at " 
                                 << dummy.i << dummy.j << dummy.k 
                                 << " smaller (posibly unstable) than L/mu0=" << dummy.Lind * deltadummy / mu0;
                    std::string buff = buff_stream4.str();
                    if ((dummy.k > ZI) && (dummy.k <= ZE)) {
                        WarnErrReport(buff);
                    }
                    dtcritico = std::min(sgg.dt / sqrt(deltadummy), dtcritico);
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
        for (int i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
            auto& segmento = HWires.CurrentSegment[i1];
            //
            if (segmento.TipoWire->HasAbsorbing_LeftEnd) {
                if (segmento.Is_LeftEnd) {
                    segmento.HasAbsorbing_LeftEnd = true;
                    std::ostringstream buff_stream5;
                    buff_stream5 << "wir1_INFO: Absorbing conditions in terminal LeftEnd segment " 
                                 << segmento.origIndex << segmento.i << segmento.j << segmento.k << segmento.tipofield;
                    std::string buff = buff_stream5.str();
                    if ((segmento.k >= ZI) && (segmento.k <= ZE) && control.verbose) {
                        WarnErrReport(buff);
                    }
                }
            }
            if (segmento.TipoWire->HasAbsorbing_RightEnd) {
                if (segmento.Is_RightEnd) {
                    segmento.HasAbsorbing_RightEnd = true;
                    std::ostringstream buff_stream6;
                    buff_stream6 << "wir1_WARNING: Absorbing conditions  in terminal RightEnd segment " 
                                 << segmento.origIndex << segmento.i << segmento.j << segmento.k << segmento.tipofield;
                    std::string buff = buff_stream6.str();
                    if ((segmento.k >= ZI) && (segmento.k <= ZE)) {
                        WarnErrReport(buff);
                    }
                }
            }
            //
            if (segmento.TipoWire->HasParallel_LeftEnd) {
                if (segmento.Is_LeftEnd) {
                    segmento.HasParallel_LeftEnd = true;
                    std::ostringstream buff_stream7;
                    buff_stream7 << "wir1_WARNING: Parallel RLC in terminal LeftEnd segment " 
                                 << segmento.origIndex << segmento.i << segmento.j << segmento.k << segmento.tipofield;
                    std::string buff = buff_stream7.str();
                    if ((segmento.k >= ZI) && (segmento.k <= ZE)) {
                        WarnErrReport(buff);
                    }
                }
            }
            if (segmento.TipoWire->HasParallel_RightEnd) {
                if (segmento.Is_RightEnd) {
                    segmento.HasParallel_RightEnd = true;
                    std::ostringstream buff_stream8;
                    buff_stream8 << "wir1_WARNING: Parallel RLC in terminal RightEnd segment " 
                                 << segmento.origIndex << segmento.i << segmento.j << segmento.k << segmento.tipofield;
                    std::string buff = buff_stream8.str();
                    if ((segmento.k >= ZI) && (segmento.k <= ZE)) {
                        WarnErrReport(buff);
                    }
                }
            }

            if (segmento.TipoWire->HasSeries_LeftEnd) {
                if (segmento.Is_LeftEnd) {
                    segmento.HasSeries_LeftEnd = true;
                    std::ostringstream buff_stream9;
                    buff_stream9 << "wir1_WARNING: Series RLC in terminal LeftEnd segment " 
                                 << segmento.origIndex << segmento.i << segmento.j << segmento.k << segmento.tipofield;
                    std::string buff = buff_stream9.str();
                    if ((segmento.k >= ZI) && (segmento.k <= ZE)) {
                        WarnErrReport(buff);
                    }
                }
            }
            if (segmento.TipoWire->HasSeries_RightEnd) {
                if (segmento.Is_RightEnd) {
                    segmento.HasSeries_RightEnd = true;
                    std::ostringstream buff_stream10;
                    buff_stream10 << "wir1_WARNING: Series RLC in terminal RightEnd segment " 
                                  << segmento.origIndex << segmento.i << segmento.j << segmento.k << segmento.tipofield;
                    std::string buff = buff_stream10.str();
                    if ((segmento.k >= ZI) && (segmento.k <= ZE)) {
                        WarnErrReport(buff);
                    }
                }
            }


        }
        //
        // Create the final update constants for the advance of the currents
        // It takes into account the extra inductance and resistance per unit length specified in ORIGINAL
        // It also takes into account the Series/Parallel Grounding Inductance at and the end segments TR and TL !untested
        // Junctions do no affect to these constants (later taken into account by means of the fractionplus and
        // fractionminus constants)
        for (int i1 = 1; i1 <= HWires.NumCurrentSegments; ++i1) {
            // constantes de actualizacion
            auto& dummy = HWires.CurrentSegment[i1];
            double resist = 0.0;

            // !!!for lossy groundings
            int i = dummy.i;
            int j = dummy.j;
            int k = dummy.k;
            int whatfield = dummy.tipofield;
            //
            double rlossy = 0.0;
            double sigt = 0.0;
            double sigtPlus = 0.0;
            double sigtMinu = 0.0;
            bool IsLossy = false; 
            bool IsPEC = false; 
            bool IsLossyPlus = false; 
            bool IsPECPlus = false;  
            bool IsLossyMinu = false; 
            bool IsPECMinu = false;       
            //
            bool esPML = false;
            switch (whatfield) {
                case iEx:
                    esPML = sgg.med(sggmiEx(i, j, k))->is.PML;
                    break;
                case iEy:
                    esPML = sgg.med(sggmiEy(i, j, k))->is.PML;
                    break;
                case iEz:
                    esPML = sgg.med(sggmiEz(i, j, k))->is.PML;
                    break;
            }
            if (esPML) { 
                continue;
                // !!\C7  dummy%IsShielded=.true.
            } else {
                if ((k <= sgg.alloc[iEZ].ZE) && (k >= sgg.alloc[iEZ].ZI)) {
                    int kmenos1 = k - 1;
                    int kmas1 = k + 1;
                    //
                    // esta informacion solo se utiliza si realmente luego hay un nodo terminal y se suma la resistencia. En cualquier otro caso no sirver para nada
                }
            }
        }

// de todos modos hay un bug en la deteccion. sgg 110815
                // se usal la informacion nodal (lo que sigue algun d\EDa se mover\E1 a la rutina de generacion nodal que se creo en preprocess y se podra dejar solo lo que sigue 110815
                // !!!bug 270224 gg Cuando se sobrepasa el MPI  kmenos1=k o kmas1=k da un buggy error. Pero no se toma decision alguna. se puede ignorar con -ignoreerrors
                switch (whatfield) {
                    case iEx:
                        med[0] = sggMiEx(i + 1, j, k);
                        med[1] = sggMiEy(i + 1, j, k);
                        med[2] = sggMiEy(i + 1, j - 1, k);
                        med[3] = sggMiEz(i + 1, j, k);
                        if (kmenos1 < sgg.alloc[iEz].ZI) {
                            kmenos1 = k;
                        }
                        med[4] = sggMiEz(i + 1, j, kmenos1);
                        med[5] = sggMiNo(i + 1, j, k);
                        //
                        med[6] = sggMiEx(i - 1, j, k);
                        med[7] = sggMiEy(i, j, k);
                        med[8] = sggMiEy(i, j - 1, k);
                        med[9] = sggMiEz(i, j, k);
                        if (kmenos1 < sgg.alloc[iEz].ZI) {
                            kmenos1 = k;
                        }
                        med[10] = sggMiEz(i, j, kmenos1);
                        med[11] = sggMiNo(i, j, k);
                        break;
                    case iEy:
                        med[0] = sggMiEy(i, j + 1, k);
                        med[1] = sggMiEz(i, j + 1, k);
                        if (kmenos1 < sgg.alloc[iEz].ZI) {
                            kmenos1 = k;
                        }
                        med[2] = sggMiEz(i, j + 1, kmenos1);
                        med[3] = sggMiEx(i, j + 1, k);
                        med[4] = sggMiEx(i - 1, j + 1, k);
                        med[5] = sggMiNo(i, j + 1, k);
                        //
                        med[6] = sggMiEy(i, j - 1, k);
                        med[7] = sggMiEz(i, j, k);
                        if (kmenos1 < sgg.alloc[iEz].ZI) {
                            kmenos1 = k;
                        }
                        med[8] = sggMiEz(i, j, kmenos1);
                        med[9] = sggMiEx(i, j, k);
                        med[10] = sggMiEx(i - 1, j, k);
                        med[11] = sggMiNo(i, j, k);
                        break;
                    case iEz:
                        // !!!ojooooo 270224 esta logica esta mal porque machaco las variables kmas1 y kmenos1.... corregir.... no tiene impacto pq el bucle es informativo y no se toman decisiones. en todo caso solo afectaria a MPI!
                        // !!!pero esta no es la razon 27024 por la que sucede el error gg   (    8/   40) wir1_BUGGYERROR:  Lossy, pec,  2.           329         180         160   15000.0000000000      F F
                        if (kmas1 > sgg.alloc[iEz].ZE) {
                            kmas1 = k;
                        }
                        med[0] = sggMiEz(i, j, kmas1);
                        if (kmas1 > sgg.alloc[iEx].ZE) {
                            kmas1 = k;
                        }
                        med[1] = sggMiEx(i, j, kmas1);
                        if (kmas1 > sgg.alloc[iEx].ZE) {
                            kmas1 = k;
                        }
                        med[2] = sggMiEx(i - 1, j, kmas1);
                        if (kmas1 > sgg.alloc[iEy].ZE) {
                            kmas1 = k;
                        }
                        med[3] = sggMiEy(i, j, kmas1);
                        if (kmas1 > sgg.alloc[iEy].ZE) {
                            kmas1 = k;
                        }
                        med[4] = sggMiEy(i, j - 1, kmas1);
                        if (kmas1 > sgg.alloc[iHz].ZE) {
                            kmas1 = k;
                        }
                        med[5] = sggMiNo(i, j, kmas1);
                        //
                        if (kmenos1 < sgg.alloc[iEz].ZI) {
                            kmenos1 = k;
                        }
                        med[6] = sggMiEz(i, j, kmenos1);
                        med[7] = sggMiEx(i, j, k);
                        med[8] = sggMiEx(i - 1, j, k);
                        med[9] = sggMiEy(i, j, k);
                        med[10] = sggMiEy(i, j - 1, k);
                        med[11] = sggMiNo(i, j, k);
                        break;
                }
                // hay que tratar cada extremo por separado
                for (nm = 0; nm <= 5; nm++) {
                    IsLossyPlus = IsLossyPlus || sgg.Med[med[nm]].Is.Lossy;
                    IsPECPlus = IsPECPlus || sgg.med[med[nm]].is.PEC || (med[nm] == 0);
                    if (!IsPECPlus) {
                        sigtPlus = std::max(sigtPlus, sgg.Med[med[nm]].sigma);
                    }
                }
                for (nm = 6; nm <= 11; nm++) {
                    IsLossyMinu = IsLossyMinu || sgg.Med[med[nm]].Is.Lossy;
                    IsPECMinu = IsPECMinu || sgg.med[med[nm]].is.PEC || (med[nm] == 0);
                    if (!IsPECMInu) {
                        sigtMinu = std::max(sigtMinu, sgg.Med[med[nm]].sigma);
                    }
                }
                // !!!telaranias quitadas 060215
                if (IsPECPlus) { // domina el pec en un nodo con multiples medios
                    IsLossyPlus = false;
                }
                if (IsPECMinu) { // domina el pec en un nodo con multiples medios
                    IsLossyMinu = false;
                }
                sigt = std::max(sigtplus, sigtMinu); // con que uno de los dos extremos sea Lossy hay que corregir la resistencia de contacto
                IsLossy = IsLossyPlus || IsLossyMinu; // con que uno de los dos extremos sea Lossy hay que corregir la resistencia de contacto
                IsPEC = IsPECPlus && IsPECMinu; // !!!importante para ser PEC tienen que serlo los dos (es un caso trivial de un segemento unido a pec por los dos sitios. pero se da en el siva)
                // !!!checking de coherencia
                if (isPEC && IsLossy) { // es decir: los dos extremos PEC y alguno lossy es que hay un error
                    sprintf(buff, "wir1_BUGGYERROR:  Lossy, pec, 1.  %d %d %d %f %d %d", i, j, k, sigt, ispec, isLossy);
                    WarnErrReport(buff, true);
                }
                // !!!
                if ((!ispec) && (!isLossy)) { // algun extremo no pec y ambos extremos no lossy debe haber error si la conductividad no es nula
                    if (std::abs(sigt) > 1.0e-19_RKIND_wires) {
                        // 270224 creo que el bug gg 270224 es simplemente que conviven pec y lossy en plus o minus siendo el otro minus/plus vacio. Se voidea
                        // 270224 el islosyplus/minus a false por ser pec. el ispec es false porque uno es vacio y el islossy tambien es vacio
                        // 270224 pero se ha almacenado alguna sigt no nula. El arreglo es simplemente convertir esto en un warning de casuistica de conexionado en vez de un buggyerror
                        // !!comentado a 010324 por los motivos anteriores    sprintf(buff, "wir1_BUGGYERROR:  Lossy, pec,  2.  %d %d %d %f %d %d", i, j, k, sigt, ispec, isLossy);
                        // !!comentado a 010324 por los motivos anteriores    WarnErrReport(buff, true);
                    }
                    sprintf(buff, "wir1_WARNING:  Lossy, pec,  2.: A wire is connected at some ending both to a lossy and to a PEC edge. Assuming it PEC  %d %d %d %f %d %d %d %d", i, j, k, sigt, ispec, isLossy, IsPECPlus, IsPECMinu);
                    WarnErrReport(buff);
                }
                if (isLossy) { // alguno extremo lossy con conductividad desconocida (multiports de ss)
                    if (std::abs(sigt) < 1.0e-19_RKIND_wires) {
                        sigt = 1e4; // asignale una resistencia de contacto por defecto si es nula !tipico de los composites de ss !habra algun dia que afinar esto
                        sprintf(buff, "wir1_WARNING:  A Lossy segment with unknown conductivity. Assuming a STANDARD value of 1e4 S/m %d %d %d %f %d %d %d %d", i, j, k, sigt, ispec, isLossy, IsLossyPlus, IsLossyMinu);
                        WarnErrReport(buff);
                    }
                }

                // !!!!!!!!!!!!!!!!!!!!hasta aqui casuistica. Aniade ahora la resistencia si procede con la formula de tercero de fisicas
                if (isLossy) {

// rlossy = rlossy + 1.0 / (2.0 * pi * dummy.TipoWire.radius * sigt) / dummy.delta; // p.u.l.
                // rlossy = rlossy + 1.0 / (2.0 * pi * dummy.DELTA * sigt) / dummy.delta; // p.u.l.
                rlossy = rlossy; // + 1.0 / (2.0 * pi * (control.factordelta * dummy.DELTA + control.factorradius * dummy.TipoWire.radius) * sigt) / dummy.delta; // p.u.l. commented out on 290323 because I don't trust this until it is validated properly agb
            }
        }
    } // end if del esPML
    //
    resist = dummy.TipoWire.R;
    givenautoin = HWires.CurrentSegment(i1).TipoWire.L;
    dummy.givenautoin = givenautoin;
    dummy.resist = resist;

    resist_devia = dummy.TipoWire.R_devia;
    givenautoin_devia = HWires.CurrentSegment(i1).TipoWire.L_devia;
    dummy.givenautoin_devia = givenautoin_devia;
    dummy.resist_devia = resist_devia;

    // !!bug'inest OLD 020413 !\E7 now only resistances are treated and added to the final segments

    if ((dummy.Is_LeftEnd) && (isLossy || isLossy)) {
        // I don't take into account the very particular case of a single segment connected to lossy at both ends!
        // we would have to add the resistance twice but the casuistica gets messy !\E7
        if ((!dummy.HasParallel_LeftEnd) && (!dummy.HasSeries_LeftEnd) && (!dummy.HasAbsorbing_LeftEnd)) {
            resist = resist + rlossy;
            sprintf(buff, "wir1_INFO: Adding Lossy material resistence to LeftEnd segment in contact with lossy without a terminal RLC %7d%7d%7d%7d %s",
                dummy.origIndex, dummy.i, dummy.j, dummy.k, dir(dummy.tipofield).c_str());
        } else {
            resist = resist + rlossy;
            sprintf(buff, "wir1_INFO: Adding Lossy material resistence to LeftEnd segment grounded through RLC %7d%7d%7d%7d %s",
                dummy.origIndex, dummy.i, dummy.j, dummy.k, dir(dummy.tipofield).c_str());
        }
        if ((dummy.k > ZI) && (dummy.k <= ZE) && (dummy.tipofield != iEz) && control.verbose) WarnErrReport(buff);
        if ((dummy.k >= ZI) && (dummy.k <= ZE) && (dummy.tipofield == iEz) && control.verbose) WarnErrReport(buff);
    }
    if ((dummy.Is_RightEnd) && (isLossy || isLossy)) {
        // I don't take into account the very particular case of a single segment connected to lossy at both ends!
        // we would have to add the resistance twice but the casuistica gets messy !\E7
        if ((!dummy.HasParallel_RightEnd) && (!dummy.HasSeries_RightEnd) && (!dummy.HasAbsorbing_RightEnd)) {
            resist = resist + rlossy;
            sprintf(buff, "wir1_INFO: Adding Lossy material resistence to RightEnd segment in contact with lossy without a terminal RLC  %7d%7d%7d%7d %s",
                dummy.origIndex, dummy.i, dummy.j, dummy.k, dir(dummy.tipofield).c_str());
        } else {
            resist = resist + rlossy;
            sprintf(buff, "wir1_INFO: Adding Lossy material resistence to RightEnd segment grounded through RLC %7d%7d%7d%7d %s",
                dummy.origIndex, dummy.i, dummy.j, dummy.k, dir(dummy.tipofield).c_str());
        }
        if ((dummy.k > ZI) && (dummy.k <= ZE) && (dummy.tipofield != iEz) && control.verbose) WarnErrReport(buff);
        if ((dummy.k >= ZI) && (dummy.k <= ZE) && (dummy.tipofield == iEz) && control.verbose) WarnErrReport(buff);
    }
    // !!! lOS THAT DO NOT HAVE RESISTANCES AND ARE OPEN DO NOT CONNECT THEM TO LOSSY. IF CONNECTIONS TO LOSSY ARE WANTED
    // !!! ONE MUST SPECIFY A LUMPED RESISTANCE
    // !!! THEN IF THEY ARE CONNECTED TO pec DIRECTLY IF THE TOPOLOGY SENDS IT
    // !!!
    if ((dummy.IsEnd_norLeft_norRight) && (isLossy || isLossy)) {
        // I don't take into account the very particular case of a single segment connected to lossy at both ends!
        // we would have to add the resistance twice but the casuistica gets messy !\E7
        if ((!dummy.HasParallel_LeftEnd) && (!dummy.HasSeries_LeftEnd) && (!dummy.HasParallel_RightEnd) && (!dummy.HasSeries_RightEnd) &&
            (!dummy.HasAbsorbing_RightEnd)) {
            resist = resist + rlossy;
            sprintf(buff, "wir1_INFO: Adding Lossy material resistence to Ending segment (other) segment in contact with lossy without a terminal RLC %7d%7d%7d%7d %s",
                dummy.origIndex, dummy.i, dummy.j, dummy.k, dir(dummy.tipofield).c_str());
            if ((dummy.k > ZI) && (dummy.k <= ZE) && (dummy.tipofield != iEz) && control.verbose) WarnErrReport(buff);
            if ((dummy.k >= ZI) && (dummy.k <= ZE) && (dummy.tipofield == iEz) && control.verbose) WarnErrReport(buff);
        } else {
            resist = resist + rlossy;
            sprintf(buff, "wir1_BUGGYERROR:  Lossy material resistence to Ending (other) segment grounded through RLC () %7d%7d%7d%7d %s",
                dummy.origIndex, dummy.i, dummy.j, dummy.k, dir(dummy.tipofield).c_str());
            if ((dummy.k > ZI) && (dummy.k <= ZE) && (dummy.tipofield != iEz)) WarnErrReport(buff, true);
            if ((dummy.k >= ZI) && (dummy.k <= ZE) && (dummy.tipofield == iEz)) WarnErrReport(buff, true);
        }
    }

    if (dummy.HasParallel_RightEnd) {
        givenautoin = givenautoin + dummy.TipoWire.Parallel_L_RightEnd / dummy.delta; // add self-inductance! 2011 \E7 untested
        dummy.givenautoin = givenautoin;
        // stoch
        givenautoin_devia = givenautoin_devia + dummy.TipoWire.Parallel_L_RightEnd_devia / dummy.delta; // add self-inductance! 2011 \E7 untested
        dummy.givenautoin_devia = givenautoin_devia;

        // !!bug'inest OLD 020413 !\E7 now only resistances are treated and added to the final segments

        resist = resist + dummy.TipoWire.Parallel_R_RightEnd / dummy.delta; // add self-inductance! 2011 \E7 untested
        resist_devia = resist_devia + dummy.TipoWire.Parallel_R_RightEnd_devia / dummy.delta; // add self-inductance! 2011 \E7 untested
        if (dummy.TipoWire.Parallel_R_RightEnd != 0.0) {
            sprintf(buff, "wir1_INFO: Adding Parallel RightEnd Resistance in segment %7d%7d%7d%7d %s",
                dummy.origIndex, dummy.i, dummy.j, dummy.k, dir(dummy.tipofield).c_str());
            if ((dummy.k > ZI) && (dummy.k <= ZE) && (dummy.tipofield != iEz) && control.verbose) WarnErrReport(buff);
            if ((dummy.k >= ZI) && (dummy.k <= ZE) && (dummy.tipofield == iEz) && control.verbose) WarnErrReport(buff);
        } else {
            sprintf(buff, "wir1_INFO: Adding Parallel RightEnd null-Resistance in segment %7d%7d%7d%7d %s",
                dummy.origIndex, dummy.i, dummy.j, dummy.k, dir(dummy.tipofield).c_str());
            if ((dummy.k > ZI) && (dummy.k <= ZE) && (dummy.tipofield != iEz) && control.verbose) WarnErrReport(buff);
            if ((dummy.k >= ZI) && (dummy.k <= ZE) && (dummy.tipofield == iEz) && control.verbose) WarnErrReport(buff);
        }

        // (note that it is per unit length the intrinsic one)
        if (dummy.TipoWire.Parallel_L_RightEnd != 0.0) {
            sprintf(buff, "wir1_INFO: Adding Parallel RightEnd Inductance in segment %7d%7d%7d%7d %s",
                dummy.origIndex, dummy.i, dummy.j, dummy.k, dir(dummy.tipofield).c_str());
            if ((dummy.k > ZI) && (dummy.k <= ZE) && (dummy.tipofield != iEz) && control.verbose) WarnErrReport(buff);
            if ((dummy.k >= ZI) && (dummy.k <= ZE) && (dummy.tipofield == iEz) && control.verbose) WarnErrReport(buff);
        }
    }

// Added also to the last segment the resistance and inductance if there are capacitances
            if (dummy.TipoWire.Parallel_C_RightEnd >= 1.0e-12) {
                std::ostringstream buff_stream;
                buff_stream << "wir1_ERROR: (Currently unsupported)  Capacitances in Parallel RightEnd at segment "
                            << std::setw(7) << dummy.origIndex
                            << std::setw(7) << dummy.i
                            << std::setw(7) << dummy.j
                            << std::setw(7) << dummy.k
                            << " " << dir(dummy.tipofield);
                std::string buff = buff_stream.str();
                if ((dummy.k > ZI) && (dummy.k <= ZE) && (dummy.tipofield != iEz)) {
                    WarnErrReport(buff, true);
                }
                if ((dummy.k >= ZI) && (dummy.k <= ZE) && (dummy.tipofield == iEz)) {
                    WarnErrReport(buff, true);
                }
            } else {
                dummy.TipoWire.Parallel_C_RightEnd = 0.0;
            }
        }
        if (dummy.HasParallel_LeftEnd) {
            givenautoin = givenautoin + dummy.TipoWire.Parallel_L_LeftEnd / dummy.delta; // adds self-inductance !2011 \E7 untested
            dummy.givenautoin = givenautoin;

            givenautoin_devia = givenautoin_devia + dummy.TipoWire.Parallel_L_LeftEnd_devia / dummy.delta; // adds self-inductance !2011 \E7 untested
            dummy.givenautoin_devia = givenautoin_devia;
            // bug'inest OLD 020413 !\E7 now only resistances are treated and added to final segments

            resist = resist + dummy.TipoWire.Parallel_R_LeftEnd / dummy.delta; // adds self-inductance !2011 \E7 untested
            resist_devia = resist_devia + dummy.TipoWire.Parallel_R_LeftEnd_devia / dummy.delta; // adds self-inductance !2011 \E7 untested

            if (dummy.TipoWire.Parallel_R_LeftEnd != 0.0) {
                std::ostringstream buff_stream;
                buff_stream << "wir1_INFO: Adding Parallel LeftEnd Resistance in segment "
                            << std::setw(7) << dummy.origIndex
                            << std::setw(7) << dummy.i
                            << std::setw(7) << dummy.j
                            << std::setw(7) << dummy.k
                            << " " << dir(dummy.tipofield);
                std::string buff = buff_stream.str();
                if ((dummy.k > ZI) && (dummy.k <= ZE) && (dummy.tipofield != iEz) && control.verbose) {
                    WarnErrReport(buff);
                }
                if ((dummy.k >= ZI) && (dummy.k <= ZE) && (dummy.tipofield == iEz) && control.verbose) {
                    WarnErrReport(buff);
                }
            } else {
                std::ostringstream buff_stream;
                buff_stream << "wir1_INFO: Adding Parallel LeftEnd null-Resistance in segment "
                            << std::setw(7) << dummy.origIndex
                            << std::setw(7) << dummy.i
                            << std::setw(7) << dummy.j
                            << std::setw(7) << dummy.k
                            << " " << dir(dummy.tipofield);
                std::string buff = buff_stream.str();
                if ((dummy.k > ZI) && (dummy.k <= ZE) && (dummy.tipofield != iEz) && control.verbose) {
                    WarnErrReport(buff);
                }
                if ((dummy.k >= ZI) && (dummy.k <= ZE) && (dummy.tipofield == iEz) && control.verbose) {
                    WarnErrReport(buff);
                }
            }

            if (dummy.TipoWire.Parallel_L_LeftEnd != 0.0) {
                std::ostringstream buff_stream;
                buff_stream << "wir1_INFO: Adding Parallel LeftEnd Inductance in segment "
                            << std::setw(7) << dummy.origIndex
                            << std::setw(7) << dummy.i
                            << std::setw(7) << dummy.j
                            << std::setw(7) << dummy.k
                            << " " << dir(dummy.tipofield);
                std::string buff = buff_stream.str();
                if ((dummy.k > ZI) && (dummy.k <= ZE) && (dummy.tipofield != iEz) && control.verbose) {
                    WarnErrReport(buff);
                }
                if ((dummy.k >= ZI) && (dummy.k <= ZE) && (dummy.tipofield == iEz) && control.verbose) {
                    WarnErrReport(buff);
                }
            }

            if (dummy.TipoWire.Parallel_C_LeftEnd >= 1.0e-12) {
                std::ostringstream buff_stream;
                buff_stream << "wir1_ERROR: (Currently unsupported)  Capacitances in Parallel LeftEnd at segment "
                            << std::setw(7) << dummy.origIndex
                            << std::setw(7) << dummy.i
                            << std::setw(7) << dummy.j
                            << std::setw(7) << dummy.k
                            << " " << dir(dummy.tipofield);
                std::string buff = buff_stream.str();
                if ((dummy.k > ZI) && (dummy.k <= ZE) && (dummy.tipofield != iEz)) {
                    WarnErrReport(buff, true);
                }
                if ((dummy.k >= ZI) && (dummy.k <= ZE) && (dummy.tipofield == iEz)) {
                    WarnErrReport(buff, true);
                }
            } else {
                dummy.TipoWire.Parallel_C_LeftEnd = 0.0;
            }
        }

        if (dummy.HasSeries_RightEnd) {
            givenautoin = givenautoin + dummy.TipoWire.Series_L_RightEnd / dummy.delta; // adds self-inductance !2011 \E7 untested
            resist = resist + dummy.TipoWire.Series_R_RightEnd / dummy.delta;
            dummy.givenautoin = givenautoin;
            dummy.resist = resist;

            givenautoin_devia = givenautoin_devia + dummy.TipoWire.Series_L_RightEnd_devia / dummy.delta; // adds self-inductance !2011 \E7 untested
            resist_devia = resist_devia + dummy.TipoWire.Series_R_RightEnd_devia / dummy.delta;
            dummy.givenautoin_devia = givenautoin_devia;
            dummy.resist_devia = resist_devia;
            // (note that intrinsic is per unit length)

            if (dummy.TipoWire.Series_L_RightEnd != 0.0) {
                std::ostringstream buff_stream;
                buff_stream << "wir1_INFO: Adding Series RightEnd Inductance in segment "
                            << std::setw(7) << dummy.origIndex
                            << std::setw(7) << dummy.i
                            << std::setw(7) << dummy.j
                            << std::setw(7) << dummy.k
                            << " " << dir(dummy.tipofield);
                std::string buff = buff_stream.str();
                if ((dummy.k > ZI) && (dummy.k <= ZE) && (dummy.tipofield != iEz) && control.verbose) {
                    WarnErrReport(buff);
                }
                if ((dummy.k >= ZI) && (dummy.k <= ZE) && (dummy.tipofield == iEz) && control.verbose) {
                    WarnErrReport(buff);
                }
            }

            if (dummy.TipoWire.Series_R_RightEnd != 0.0) {
                std::ostringstream buff_stream;
                buff_stream << "wir1_INFO: Adding Series RightEnd Resistance in segment "
                            << std::setw(7) << dummy.origIndex
                            << std::setw(7) << dummy.i
                            << std::setw(7) << dummy.j
                            << std::setw(7) << dummy.k
                            << " " << dir(dummy.tipofield);
                std::string buff = buff_stream.str();
                if ((dummy.k > ZI) && (dummy.k <= ZE) && (dummy.tipofield != iEz) && control.verbose) {
                    WarnErrReport(buff);
                }
                if ((dummy.k >= ZI) && (dummy.k <= ZE) && (dummy.tipofield == iEz) && control.verbose) {
                    WarnErrReport(buff);
                }
            }

            if (dummy.TipoWire.Series_C_RightEnd <= 1.0e7) {
                std::ostringstream buff_stream;
                buff_stream << "wir1_ERROR: (Currently unsupported)  Capacitances smaller than 1.0e7 (inf) in Series RightEnd at segment "
                            << std::setw(7) << dummy.origIndex
                            << std::setw(7) << dummy.i
                            << std::setw(7) << dummy.j
                            << std::setw(7) << dummy.k
                            << " " << dir(dummy.tipofield);
                std::string buff = buff_stream.str();
                if ((dummy.k > ZI) && (dummy.k <= ZE) && (dummy.tipofield != iEz)) {
                    WarnErrReport(buff, true);
                }
                if ((dummy.k >= ZI) && (dummy.k <= ZE) && (dummy.tipofield == iEz)) {
                    WarnErrReport(buff, true);
                }
            } else {
                dummy.TipoWire.Series_C_RightEnd = 2.0e7;
            }
        }

        if (dummy.HasSeries_LeftEnd) {
            givenautoin = givenautoin + dummy.TipoWire.Series_L_LeftEnd / dummy.delta; // adds self-inductance !2011 \E7 untested
            resist = resist + dummy.TipoWire.Series_R_LeftEnd / dummy.delta; // adds self-inductance !2011 \E7 untested
            dummy.givenautoin = givenautoin;
            dummy.resist = resist;

            givenautoin_devia = givenautoin_devia + dummy.TipoWire.Series_L_LeftEnd_devia / dummy.delta; // adds self-inductance !2011 \E7 untested
            resist_devia = resist_devia + dummy.TipoWire.Series_R_LeftEnd_devia / dummy.delta; // adds self-inductance !2011 \E7 untested
            dummy.givenautoin_devia = givenautoin_devia;
            dummy.resist_devia = resist_devia;

            if (dummy.TipoWire.Series_L_LeftEnd != 0.0) {
                std::ostringstream buff_stream;
                buff_stream << "wir1_INFO: Adding Series LeftEnd Inductance in segment "
                            << std::setw(7) << dummy.origIndex
                            << std::setw(7) << dummy.i
                            << std::setw(7) << dummy.j
                            << std::setw(7) << dummy.k
                            << " " << dir(dummy.tipofield);
                std::string buff = buff_stream.str();
                if ((dummy.k > ZI) && (dummy.k <= ZE) && (dummy.tipofield != iEz) && control.verbose) {
                    WarnErrReport(buff);
                }
                if ((dummy.k >= ZI) && (dummy.k <= ZE) && (dummy.tipofield == iEz) && control.verbose) {
                    WarnErrReport(buff);
                }
            }

            if (dummy.TipoWire.Series_R_LeftEnd != 0.0) {
                std::ostringstream buff_stream;
                buff_stream << "wir1_INFO: Adding Series LeftEnd Resistance in segment "
                            << std::setw(7) << dummy.origIndex
                            << std::setw(7) << dummy.i
                            << std::setw(7) << dummy.j
                            << std::setw(7) << dummy.k
                            << " " << dir(dummy.tipofield);
                std::string buff = buff_stream.str();
                if ((dummy.k > ZI) && (dummy.k <= ZE) && (dummy.tipofield != iEz) && control.verbose) {
                    WarnErrReport(buff);
                }

