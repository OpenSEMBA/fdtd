#include <vector>
#include <memory>
#include <string>
#include <cstdint>
#include <functional>

// Assuming FDETYPES_m provides RKIND and RKIND_wires definitions.
// Typically these map to double or float. We assume double for scientific computing.
#ifndef RKIND
#define RKIND double
#endif
#ifndef RKIND_wires
#define RKIND_wires double
#endif

// Forward declarations for types defined in other modules or later in this file
struct CurrentSegments_t;
struct ChargeNodes_t;
struct wires_t;
struct source_t;
struct thick_t;
struct container_t;
struct TSegmentPtr_t;
struct TMultiline_t;

namespace wiresHolland_constants_m {

   constexpr int32_t MaxNumCurrentMinusPlus = 9;

   struct ChargeNodes_t {
      int32_t IndexNode;
      
      // Pointers to neighbours. Using raw pointers as in Fortran, but managed externally.
      CurrentSegments_t* CurrentPlus_1 = nullptr;
      CurrentSegments_t* CurrentMinus_1 = nullptr;
      CurrentSegments_t* CurrentPlus_2 = nullptr;
      CurrentSegments_t* CurrentMinus_2 = nullptr;
      CurrentSegments_t* CurrentPlus_3 = nullptr;
      CurrentSegments_t* CurrentMinus_3 = nullptr;
      CurrentSegments_t* CurrentPlus_4 = nullptr;
      CurrentSegments_t* CurrentMinus_4 = nullptr;
      CurrentSegments_t* CurrentPlus_5 = nullptr;
      CurrentSegments_t* CurrentMinus_5 = nullptr;
      CurrentSegments_t* CurrentPlus_6 = nullptr;
      CurrentSegments_t* CurrentMinus_6 = nullptr;
      CurrentSegments_t* CurrentPlus_7 = nullptr;
      CurrentSegments_t* CurrentMinus_7 = nullptr;
      CurrentSegments_t* CurrentPlus_8 = nullptr;
      CurrentSegments_t* CurrentMinus_8 = nullptr;
      CurrentSegments_t* CurrentPlus_9 = nullptr;
      CurrentSegments_t* CurrentMinus_9 = nullptr;

      bool IsMur = false;
      bool IsPeriodic = false;
      bool IsAttachedtoVoltage = false;
      bool IsPEC = false;
      bool HasIsource = false;
      bool Exists = false;
      bool Is_LeftEnd = false;
      bool Is_RightEnd = false;
      bool IsLossy = false;
      bool IsBackDownLeftMur = false;
      bool IsFrontUpRightMur = false;
      bool proc = false; // dama
      bool IsHeterogeneousJunction = false;
      bool IsInSingleRLCsegment = false;
      
      RKIND_wires cteMur = 0.0;
      RKIND_wires ctePlain = 0.0;
      RKIND_wires origctePlain = 0.0;
      RKIND_wires cteprop = 0.0;
      
      // to apply Mur. Needs extra storage everywhere but it is only 1D
      RKIND_wires ChargePresent = 0.0;
      RKIND_wires ChargePast = 0.0;
      
      ChargeNodes_t* NodeInside = nullptr;
      
      int32_t NumCurrentMinus = 0;
      int32_t NumCurrentPlus = 0;
      int32_t i = 0;
      int32_t j = 0;
      int32_t k = 0;
      
      source_t* Isource = nullptr;
      
      // 1-based indexing preserved via size 2*MaxNumCurrentMinusPlus + 1, index 0 unused
      std::vector<int32_t> YESsegment; 

      // Pointers initialized to null
      RKIND* already_YEEadvanced_byconformal_changedtoPECfield1 = nullptr;
      RKIND* already_YEEadvanced_byconformal_changedtoPECfield2 = nullptr;
      RKIND* already_YEEadvanced_byconformal_changedtoPECfield3 = nullptr;
      RKIND* already_YEEadvanced_byconformal_changedtoPECfield4 = nullptr;
      RKIND* already_YEEadvanced_byconformal_changedtoPECfield5 = nullptr;
      RKIND* already_YEEadvanced_byconformal_changedtoPECfield6 = nullptr;

#ifdef CompileWithMPI
      // For MPI purposes !only handled and initialized in MPIcomm
      CurrentSegments_t* MPISharedCurrent = nullptr;
#endif
   };

#ifdef CompileWithThickWires    
   struct container_t {
       RKIND* punt = nullptr;
       RKIND retardo = 0.0;      
       std::vector<RKIND> field_retard;  
   };
   
   struct thick_t {    
      int32_t Enumero = 0;
      int32_t Hnumero = 0;
      std::vector<container_t> Efield_wire2main;
      std::vector<container_t> Hfield_wire2main;
      std::vector<container_t> H_Efield_wire2main;
      std::vector<RKIND_wires> EArea;
      std::vector<RKIND_wires> rEArea;
      std::vector<RKIND_wires> HArea;
      std::vector<RKIND_wires> rHArea;
      std::vector<RKIND_wires> rEfractionArea;
      std::vector<RKIND_wires> Hsigno;
      std::vector<RKIND_wires> Hcte;  
      std::vector<int32_t> i;
      std::vector<int32_t> j;
      std::vector<int32_t> k;
      std::vector<int32_t> field;     
      bool Hplus = false;         
      std::vector<RKIND> Current_ret;    
      int32_t maxretardo = 0;
   };
#endif       

   struct CurrentSegments_t {
      int32_t IndexSegment = 0;
      int32_t NumParallel = 0;
      int32_t OrigIndex = 0;
      
      wires_t* TipoWire = nullptr;
      
      RKIND_wires Lind = 0.0;
      RKIND_wires inv_Lind_acum = 0.0;
      RKIND_wires HEUR_safety = 0.0;
      RKIND_wires Lind_acum = 0.0;
      
      RKIND_wires delta = 0.0;
      RKIND_wires deltaTransv1 = 0.0;
      RKIND_wires deltaTransv2 = 0.0;
      
      RKIND_wires givenautoin = 0.0;
      RKIND_wires resist = 0.0;
      
      RKIND_wires givenautoin_devia = 0.0;
      RKIND_wires resist_devia = 0.0;
      
      ChargeNodes_t* ChargePlus = nullptr;
      ChargeNodes_t* ChargeMinus = nullptr; // neighbours in the plus and Minus direction
      
      bool IsPMC = false;
      bool HasVsource = false;
      bool IsShielded = false;
      bool HasParallel_RightEnd = false;
      bool HasParallel_LeftEnd = false;
      bool HasSeries_RightEnd = false;
      bool HasSeries_LeftEnd = false;
      bool HasAbsorbing_RightEnd = false;
      bool HasAbsorbing_LeftEnd = false;
      
      bool Is_LeftEnd = false;
      bool Is_RightEnd = false;
      bool IsEnd_norLeft_norRight = false;
      bool proc = false;
      bool IsConformal = false;
      
      RKIND_wires cte1 = 0.0;
      RKIND_wires cte2 = 0.0;
      RKIND_wires cte3 = 0.0;
      RKIND_wires cte5 = 0.0;
      RKIND_wires FractionPlus = 0.0;
      RKIND_wires FractionMinus = 0.0;
      
      RKIND_wires Current = 0.0;
      RKIND_wires qplus_qminus = 0.0;
      
      RKIND_wires CurrentPast = 0.0; // added just for right observation
      // at the desired time step in observation.f90       
      
      RKIND* Efield_wire2main = nullptr;
      RKIND* Efield_main2wire = nullptr;
      
#ifdef CompileWithThickWires      
      thick_t thick;
#endif
//      real(kind=RKIND_wires)                           :: Efield_wire2main_past  !no sirve para nada 171216
      int32_t i = 0;
      int32_t j = 0;
      int32_t k = 0;
      int32_t indexmed = 0;
      int32_t ILIBRE = 0;
      int32_t JLIBRE = 0;
      int32_t KLIBRE = 0;
      // dama
      int32_t ie = 0;
      int32_t je = 0;
      int32_t ke = 0;
      RKIND_wires x = 0.0;
      RKIND_wires y = 0.0;
      RKIND_wires L = 0.0;
      RKIND_wires C = 0.0;
      RKIND_wires R = 0.0;
      RKIND_wires L_devia = 0.0;
      RKIND_wires C_devia = 0.0;
      RKIND_wires R_devia = 0.0;
      RKIND_wires cI = 0.0;
      RKIND_wires bI = 0.0;
      RKIND_wires Lintrinsic = 0.0;
      // fin dama
      int32_t tipofield = 0; // iEx,iEy o iEz
      bool orientadoalreves = false;
      
      source_t* Vsource = nullptr;
      
#ifdef CompileWithMPI
      // only required by the new MPI wires routines march'12 2012 bug multiwires MPI
      int32_t equivalentIndex = 0;
#endif

      // !!!crank-nicolson coefficients
      RKIND_wires upperdiag = 0.0;
      RKIND_wires diag = 0.0;
      RKIND_wires lowerdiag = 0.0;
      RKIND_wires rightCHminus = 0.0;
      RKIND_wires rightCHplus = 0.0;
      RKIND_wires rightCU = 0.0;
      RKIND_wires rightCUminus = 0.0;
      RKIND_wires rightCUplus = 0.0;
      // !!!!!!!!!!end crank-nicolson
      
      // !!!se aniade siempre aunque solo lo use stochastic
      RKIND_wires qplus_qminus_for_devia = 0.0;
      RKIND_wires current_for_devia = 0.0;
      RKIND_wires Efield_main2wire_for_devia = 0.0;
      RKIND_wires Lind_devia = 0.0;
      RKIND_wires cte1_for_devia = 0.0;
      RKIND_wires cte2_for_devia = 0.0;
      RKIND_wires cte3_for_devia = 0.0;  
   };

   // dama
   struct TSegmentPtr_t {
      CurrentSegments_t* ptr = nullptr;
   };

   struct TMultiline_t {
      int32_t NumParallel = 0;
      std::vector<TSegmentPtr_t*> Segments;
      std::vector<std::vector<RKIND_wires>> R;
      std::vector<std::vector<RKIND_wires>> L;
      std::vector<std::vector<RKIND_wires>> C;
      std::vector<std::vector<RKIND_wires>> b1I;
      std::vector<std::vector<RKIND_wires>> b2I;
      std::vector<std::vector<RKIND_wires>> b3I;
   };
   // !!!!!!!!!!!!fin dama

   struct ThinWires_t {
      int32_t NumMultilines = 0; // dama
      std::vector<TMultiline_t*> Multilines; // dama
      
      int32_t NumDifferentWires = 0;
      int32_t NumCurrentSegments = 0;
      int32_t NumChargeNodes = 0;
      
      std::vector<int32_t> WireTipoMedio;
      
      CurrentSegments_t NullSegment; // contiene informacion nula precisada por segmentos voided pero observados en la rutina de observacion 12/09/13
      ChargeNodes_t NullNode;
      
      std::vector<CurrentSegments_t*> CurrentSegment;
      std::vector<ChargeNodes_t*> ChargeNode;
      
#ifdef CompileWithMPI
      // For MPI purposes !only handled and initialized in MPIcomm
      std::vector<CurrentSegments_t*> MPIUpNeededCurrentSegment;
      std::vector<CurrentSegments_t*> MPIDownNeededCurrentSegment;
      int32_t NumNeededCurrentUpMPI = 0;
      int32_t NumNeededCurrentDownMPI = 0;
      std::vector<ChargeNodes_t*> MPIUpChargeNode;
      std::vector<ChargeNodes_t*> MPIDownChargeNode;
      // only required by the new MPI wires routines march'12 2012 bug multiwires MPI
      std::vector<CurrentSegments_t*> MPIUpSharedCurrentSegment;
      std::vector<CurrentSegments_t*> MPIDownSharedCurrentSegment;
      int32_t NumSharedCurrentUpMPI = 0;
      int32_t NumSharedCurrentDownMPI = 0;
#endif
      
      RKIND null_field = 0.0; // en los segmentos embeddeds y en los paralelos no hay acople entre thin-wire y medio
      RKIND_wires olddt = 0.0; // para permit scaling 141118
      // apunto a null_field el pointer field anterior en vez de al campo fdtd y lo obligo a ser cero
   };

   struct adyacc_t {
      bool Is = false;
      bool Parallel = false;
      bool IsHeterogeneousJunction = false;
      bool BothEndingsConnected = false;
      int32_t i = 0;
      int32_t j = 0;
      int32_t k = 0;
      // 1-based indexing preserved via size 3, index 0 unused
      std::vector<int32_t> YESsegment; 
   };

} // namespace wiresHolland_constants_m