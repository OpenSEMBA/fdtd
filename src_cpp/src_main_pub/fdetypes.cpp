```cpp
#include <vector>
#include <string>
#include <complex>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <algorithm>

// Note: MPI definitions would typically come from mpi.h if CompileWithMPI is defined
// Note: OpenMP definitions would typically come from omp.h if CompileWithOpenMP is defined

#ifdef CompileWithMPI
#include <mpi.h>
#endif

#ifdef CompileWithOpenMP
#include <omp.h>
#endif

// Handle Real/Int Kind definitions based on preprocessor flags
// Defaulting to standard sizes if not explicitly defined by build system

#ifndef CompileWithReal16
#ifndef CompileWithReal8
#ifndef CompileWithReal4
#define CompileWithReal4
#endif
#endif
#endif

#ifndef CompileWithInt4
#ifndef CompileWithInt2
#ifndef CompileWithInt1
#define CompileWithInt4
#endif
#endif
#endif

// Define basic types based on flags
#ifdef CompileWithReal16
    // Typically long double or a custom 128-bit float type
    using RealKind = long double;
    using RealKindWires = long double;
    using RealKindTiempo = long double;
    using CharKind = long double; // CKIND usually matches RKIND for complex
    #define RKIND_VAL 16
    #define DOUBLE_VAL 16
    #define SINGLE_VAL 4
    #define LONG_DOUBLE_VAL 16
#elif defined(CompileWithReal8)
    using RealKind = double;
    using RealKindWires = double;
    using RealKindTiempo = double;
    using CharKind = double;
    #define RKIND_VAL 8
    #define DOUBLE_VAL 8
    #define SINGLE_VAL 4
    #define LONG_DOUBLE_VAL 16
#else
    // Default to Real4 (float) but RKIND_wires and RKIND_tiempo default to double as per Fortran code logic
    using RealKind = float;
    using RealKindWires = double;
    using RealKindTiempo = double;
    using CharKind = double;
    #define RKIND_VAL 4
    #define DOUBLE_VAL 8
    #define SINGLE_VAL 4
    #define LONG_DOUBLE_VAL 16
#endif

#ifdef CompileWithInt1
    using IntKindMedia = int8_t;
    #ifdef CompileWithMPI
        // MPI_INTEGER1 is not standard in all MPI implementations, usually mapped to MPI_CHAR or similar
        // Assuming standard mapping for translation purposes
        const int INTEGERSIZE = MPI_CHAR; 
    #endif
    #define INTEGERSIZEOFMEDIAMATRICES 1
#elif defined(CompileWithInt2)
    using IntKindMedia = int16_t;
    #ifdef CompileWithMPI
        const int INTEGERSIZE = MPI_SHORT;
    #endif
    #define INTEGERSIZEOFMEDIAMATRICES 2
#else
    // Default Int4
    using IntKindMedia = int32_t;
    #ifdef CompileWithMPI
        const int INTEGERSIZE = MPI_INT;
    #endif
    #define INTEGERSIZEOFMEDIAMATRICES 4
#endif

// IKINDMTAG is always 4
const int IKINDMTAG = 4;

// Constants for RKIND values in integer form for parameters
const int2_t SINGLE_INT = 4;
const int2_t DOUBLE_INT = 8;
const int2_t LONG_DOUBLE_INT = 16;

// Helper for complex types
using ComplexKind = std::complex<RealKind>;
using ComplexKindWires = std::complex<RealKindWires>;

// MPI Size Constants (Approximations for translation, actual values depend on MPI implementation)
#ifdef CompileWithMPI
    const int REALSIZE = MPI_DOUBLE_PRECISION; // Defaulting to double precision for safety in translation unless specified
    const int REALSIZE_wires = MPI_DOUBLE_PRECISION;
    const int COMPLEXSIZE = MPI_DOUBLE_COMPLEX;
    const int REALSIZE_tiempo = MPI_DOUBLE_PRECISION;
#else
    const int REALSIZE = 8;
    const int REALSIZE_wires = 8;
    const int COMPLEXSIZE = 16;
    const int REALSIZE_tiempo = 8;
#endif

namespace FDETYPES_m {

    // Global Variables (Static to namespace to mimic module scope)
    int32_t quienmpi = 0;
    int32_t tamaniompi = 0;
    int32_t SUBCOMM_MPI = 0;
    int32_t SUBCOMM_MPI_conformal_probes = 0;
    int32_t MPI_conformal_probes_root = 0;

    // Parameters
    const int64_t maxmpibytes = (1LL << 27);
    const int32_t BuffObse = (1 << 10);
    const int64_t MaxMemoryProbes = (1LL << 37);
    const int32_t MaxProbes = 150000;
    const int topCPUtime = 10000000;
    const int BUFSIZE = 1024;
    const int BUFSIZE_LONG = 16384;

    // Priorities (Global variables)
    int32_t prior_BV = 0;
    int32_t prior_IB = 0;
    int32_t prior_pmlbody = 0;
    int32_t prior_AB = 0;
    int32_t prior_FDB = 0;
    int32_t prior_IS = 0;
    int32_t prior_AS = 0;
    int32_t prior_FDS = 0;
    int32_t prior_IL = 0;
    int32_t prior_AL = 0;
    int32_t prior_FDL = 0;
    int32_t prior_IP = 0;
    int32_t prior_AP = 0;
    int32_t prior_FDP = 0;
    int32_t prior_PEC = 0;
    int32_t prior_PMC = 0;
    int32_t prior_TG = 0;
    int32_t prior_CS = 0;
    int32_t prior_TW = 0;

    // Conformal flag
    bool input_conformal_flag = false;

    // Tunable Parameters
    const int2_t INTEGERSIZEOFMEDIAMATRICES_VAL = INTEGERSIZEOFMEDIAMATRICES;
    const int2_t IKINDMTAG_VAL = 4;

    const int2_t SINGLE_VAL = 4;
    const int2_t DOUBLE_VAL = 8;
    const int2_t LONG_DOUBLE_VAL = 16;

    const int2_t RKIND_VAL_INT = 
#ifdef CompileWithReal8
        DOUBLE_VAL;
#elif defined(CompileWithReal16)
        LONG_DOUBLE_VAL;
#else
        SINGLE_VAL;
#endif

    const int2_t RKIND_wires_VAL_INT = 
#ifdef CompileWithReal8
        DOUBLE_VAL;
#elif defined(CompileWithReal16)
        LONG_DOUBLE_VAL;
#else
        DOUBLE_VAL; // Default logic in Fortran
#endif

    const int2_t RKIND_tiempo_VAL_INT = 
#ifdef CompileWithReal8
        DOUBLE_VAL;
#elif defined(CompileWithReal16)
        LONG_DOUBLE_VAL;
#else
        DOUBLE_VAL;
#endif

    const int2_t CKIND_VAL_INT = 
#ifdef CompileWithReal8
        DOUBLE_VAL;
#elif defined(CompileWithReal16)
        LONG_DOUBLE_VAL;
#else
        DOUBLE_VAL;
#endif

    // RKIND constants
    const RealKind plusCPU_PML = 2.0; // Heuristic
    const RealKind heurCFL = 0.8;
    const RealKind pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148;
    const RealKind unmedio = 0.5;
    const ComplexKind mcPI2 = -ComplexKind(0.0, 1.0) * 2.0 * pi;

    const int32_t Down = 1, Up = 2, Left = 3, Right = 4, Back = 5, Front = 6;
    const int32_t iEx = 1, iEy = 2, iEz = 3, iHx = 4, iHy = 5, iHz = 6, centroide = 8, Nothing = 666;
    const int32_t iMEC = 51;
    const int32_t iMHC = 52;
    const int32_t iCur = 53;
    const int32_t iCurX = 54;
    const int32_t iCurY = 55;
    const int32_t iCurZ = 56;
    const int32_t mapvtk = 57;
    const int32_t iExC = 61;
    const int32_t iEyC = 62;
    const int32_t iEzC = 63;
    const int32_t iHxC = 64;
    const int32_t iHyC = 65;
    const int32_t iHzC = 66;
    const int32_t farfield = 67;
    const int32_t lineIntegral = 68;
    const int32_t iJx = 10 * iEx;
    const int32_t iJy = 10 * iEy;
    const int32_t iJz = 10 * iEz;
    const int32_t iQx = 10000 * iEx;
    const int32_t iQy = 10000 * iEy;
    const int32_t iQz = 10000 * iEz;
    const int32_t iVx = 1000 * iEx;
    const int32_t iVy = 1000 * iEy;
    const int32_t iVz = 1000 * iEz;
    const int32_t iBloqueJx = 100 * iEx;
    const int32_t iBloqueJy = 100 * iEy;
    const int32_t iBloqueJz = 100 * iEz;
    const int32_t iBloqueMx = 100 * iHx;
    const int32_t iBloqueMy = 100 * iHy;
    const int32_t iBloqueMz = 100 * iHz;

    const std::string SEPARADOR = "______________";
    const int32_t comi = 1, fine = 2, icoord = 1, jcoord = 2, kcoord = 3;

    const RealKind EPSILON_VACUUM = 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12;
    const RealKind MU_VACUUM = 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6;
    const RealKind c_vacuum = 1.0 / std::sqrt(EPSILON_VACUUM * MU_VACUUM);

    RealKind dt0; // Global variable for OLDrlo access

    const int32_t FACE_X = 1;
    const int32_t FACE_Y = 2;
    const int32_t FACE_Z = 3;
    const int32_t EDGE_X = 1;
    const int32_t EDGE_Y = 2;
    const int32_t EDGE_Z = 3;

    const std::string F_SOURCE_VOLTAGE = "VOLT";
    const std::string F_SOURCE_CURRENT = "CURR";

    const std::string fmt = 
#ifdef CompileWithReal4
        "(e27.17e3,11(e19.9e3))";
#elif defined(CompileWithReal8)
        "(12(e27.17e3))";
#elif defined(CompileWithReal16)
        "(12(e46.36e3))";
#else
        "(e27.17e3,11(e19.9e3))";
#endif

    // Derived Types

    struct tagtype_t {
        std::vector<std::string> tag;
        int32_t numertags = 0;
    };

    struct coorsxyz_t {
        std::vector<RealKind> x;
        std::vector<RealKind> y;
        std::vector<RealKind> z;
    };

    struct coorsxyzP_t {
        coorsxyz_t PhysCoor[6]; // 1-based index in Fortran, mapped to 0-5 here or kept as struct array
    };

    struct MedioExtra_t {
        int32_t pml_size = 0;
        int32_t index = 0;
        RealKind sigma = 0.0;
        RealKind sigmam = 0.0;
        bool exists = false;
    };

    struct logic_control_t {
        bool Wires = false;
        bool PMLbodies = false;
        bool MultiportS = false;
        bool AnisMultiportS = false;
        bool SGBCs = false;
        bool Lumpeds = false;
        bool EDispersives = false;
        bool MDispersives = false;
        bool PlaneWaveBoxes = false;
        bool Observation = false;
        bool FarFields = false;
        bool PMCBorders = false;
        bool PMLBorders = false;
        bool MurBorders = false;
        bool PECBorders = false;
        bool PeriodicBorders = false;
        bool Anisotropic = false;
        bool ThinSlot = false;
        bool NodalE = false;
        bool NodalH = false;
        bool MagneticMedia = false;
        bool PMLMagneticMedia = false;
        bool MTLNbundles = false;

        void reset() {
            Wires = false;
            PMLbodies = false;
            MultiportS = false;
            AnisMultiportS = false;
            SGBCs = false;
            Lumpeds = false;
            EDispersives = false;
            MDispersives = false;
            PlaneWaveBoxes = false;
            Observation = false;
            FarFields = false;
            PMCBorders = false;
            PMLBorders = false;
            MurBorders = false;
            PECBorders = false;
            Anisotropic = false;
            ThinSlot = false;
            NodalE = false;
            NodalH = false;
            PeriodicBorders = false;
            MagneticMedia = false;
            PMLMagneticMedia = false;
            MTLNbundles = false;
        }
    };

    struct Xlimit_t {
        int32_t XI = 0, XE = 0, NX = 0;
    };

    struct Ylimit_t {
        int32_t YI = 0, YE = 0, NY = 0;
    };

    struct Zlimit_t {
        int32_t ZI = 0, ZE = 0, NZ = 0;
    };

    struct limit_t {
        int32_t XI = 0, XE = 0, YI = 0, YE = 0, ZI = 0, ZE = 0, NX = 0, NY = 0, NZ = 0;
    };

    struct XYZlimit_t {
        int32_t XI = 0, XE = 0, YI = 0, YE = 0, ZI = 0, ZE = 0;
    };

    struct xyzlimit_scaled_t {
        int32_t XI = 0, XE = 0, YI = 0, YE = 0, ZI = 0, ZE = 0;
        RealKind xc = 0.0, yc = 0.0, zc = 0.0;
        int32_t Or = 0;
    };

    struct tagnumber_t {
        std::vector<std::vector<std::vector<IntKindMedia>>> x;
        std::vector<std::vector<std::vector<IntKindMedia>>> y;
        std::vector<std::vector<std::vector<IntKindMedia>>> z;
    };

    struct taglist_t {
        tagnumber_t edge;
        tagnumber_t face;

        IntKindMedia getFaceTag(int32_t field, int32_t i, int32_t j, int32_t k) const {
            switch (field) {
                case iHx: return face.x[i][j][k];
                case iHy: return face.y[i][j][k];
                case iHz: return face.z[i][j][k];
                default: return 0;
            }
        }

        IntKindMedia getEdgeTag(int32_t field, int32_t i, int32_t j, int32_t k) const {
            switch (field) {
                case iEx: return edge.x[i][j][k];
                case iEy: return edge.y[i][j][k];
                case iEz: return edge.z[i][j][k];
                default: return 0;
            }
        }
    };

    struct bounds_t {
        limit_t sggMiEx, sggMiEy, sggMiEz, sggMiHx, sggMiHy, sggMiHz;
        limit_t Ex, Ey, Ez, Hx, Hy, Hz;
        limit_t sweepEx, sweepEy, sweepEz, sweepHx, sweepHy, sweepHz;
        limit_t sweepSINPMLEx, sweepSINPMLEy, sweepSINPMLEz, sweepSINPMLHx, sweepSINPMLHy, sweepSINPMLHz;
        Xlimit_t dxe, dxh;
        Ylimit_t dye, dyh;
        Zlimit_t dze, dzh;
    };

    struct PML_t {
        RealKind orden[3][2] = {};
        RealKind CoeffReflPML[3][2] = {};
        int32_t NumLayers[3][2] = {};
    };

    struct fichevol_t {
        std::string Name;
        int32_t NumSamples = 0;
        RealKind DeltaSamples = 0.0;
        std::vector<RealKind> Samples;
    };

    struct fichevol_wires_t {
        std::string Name;
        int32_t NumSamples = 0;
        RealKindWires DeltaSamples = 0.0;
        std::vector<RealKindWires> Samples;
    };

    struct source_t {
        fichevol_wires_t Fichero;
        RealKindWires Resistance = 0.0;
        RealKindWires Multiplier = 0.0;
        int32_t i = 0, j = 0, k = 0;
    };

    struct NodalSource_t {
        fichevol_t Fichero;
        std::vector<xyzlimit_scaled_t> punto;
        int32_t numpuntos = 0;
        bool IsInitialValue = false;
        bool IsHard = false;
        bool IsElec = false;
    };

    struct WireDispersiveParams_t {
        int32_t numPoles = 0;
        std::vector<ComplexKind> res;
        std::vector<ComplexKind> p;
        ComplexKind d = ComplexKind(0.0, 0.0);
        ComplexKind e = ComplexKind(0.0, 0.0);
    };

    struct oriented_point_t {
        int32_t ori = 0;
        int32_t i = 0, j = 0, k = 0;
        int32_t origIndex = 0, ilibre = 0, jlibre = 0, klibre = 0;
        int32_t multiraboDE = 0;
        bool Is_LeftEnd = false;
        bool Is_RightEnd = false;
        bool IsEnd_norLeft_norRight = false;
        bool repetido = false;
        bool multirabo = false;
        bool orientadoalreves = false;
    };

#ifdef CompileWithMTLN
    struct Multiwires_t {};
#endif

    struct Wires_t {
        RealKindWires Radius = 0.0, R = 0.0, L = 0.0, C = 0.0;
        RealKindWires P_R = 0.0, P_L = 0.0, P_C = 0.0;
        RealKindWires Radius_devia = 0.0, R_devia = 0.0, L_devia = 0.0, C_devia = 0.0;
        std::vector<WireDispersiveParams_t> disp;
        int32_t numsegmentos = 0, NUMVOLTAGESOURCES = 0, NUMCURRENTSOURCES = 0;
        std::vector<oriented_point_t> segm;
        std::vector<source_t> Vsource;
        std::vector<source_t> Isource;
        bool VsourceExists = false;
        bool IsourceExists = false;
        bool HasParallel_LeftEnd = false;
        bool HasParallel_RightEnd = false;
        bool HasSeries_LeftEnd = false;
        bool HasSeries_RightEnd = false;
        bool HasAbsorbing_LeftEnd = false;
        bool HasAbsorbing_RightEnd = false;
        RealKindWires Parallel_R_RightEnd = 0.0, Parallel_R_LeftEnd = 0.0;
        RealKindWires Series_R_RightEnd = 0.0, Series_R_LeftEnd = 0.0;
        RealKindWires Parallel_L_RightEnd = 0.0, Parallel_L_LeftEnd = 0.0;
        RealKindWires Series_L_RightEnd = 0.0, Series_L_LeftEnd = 0.0;
        RealKindWires Parallel_C_RightEnd = 0.0, Parallel_C_LeftEnd = 0.0;
        RealKindWires Series_C_RightEnd = 0.0, Series_C_LeftEnd = 0.0;
        RealKindWires Parallel_R_RightEnd_devia = 0.0, Parallel_R_LeftEnd_devia = 0.0;
        RealKindWires Series_R_RightEnd_devia = 0.0, Series_R_LeftEnd_devia = 0.0;
        RealKindWires Parallel_L_RightEnd_devia = 0.0, Parallel_L_LeftEnd_devia = 0.0;
        RealKindWires Series_L_RightEnd_devia = 0.0, Series_L_LeftEnd_devia = 0.0;
        RealKindWires Parallel_C_RightEnd_devia = 0.0, Parallel_C_LeftEnd_devia = 0.0;
        RealKindWires Series_C_RightEnd_devia = 0.0, Series_C_LeftEnd_devia = 0.0;
        std::vector<WireDispersiveParams_t> disp_LeftEnd;
        std::vector<WireDispersiveParams_t> disp_RightEnd;
        int32_t LeftEnd = 0, RightEnd = 0;
    };

    struct SlantedNode_t {
        int32_t index = 0;
        RealKindWires x = 0.0, y = 0.0, z = 0.0;
        bool VsourceExists = false;
        bool IsourceExists = false;
        source_t* Vsource = nullptr;
        source_t* Isource = nullptr;
    };

    struct SlantedWires_t {
        RealKindWires radius = 0.0, R = 0.0, L = 0.0, C = 0.0;
        RealKindWires P_R = 0.0, P_L = 0.0, P_C = 0.0;
        std::vector<WireDispersiveParams_t> disp;
        int32_t LeftEnd = 0, RightEnd = 0;
        int32_t NumNodes = 0;
        std::vector<SlantedNode_t> nodes;
        bool HasParallel_LeftEnd = false;
        RealKindWires Parallel_R_LeftEnd = 0.0, Parallel_L_LeftEnd = 0.0, Parallel_C_LeftEnd = 0.0;
        bool HasParallel_RightEnd = false;
        RealKindWires Parallel_R_RightEnd = 0.0, Parallel_L_RightEnd = 0.0, Parallel_C_RightEnd = 0.0;
        bool HasSeries_LeftEnd = false;
        RealKindWires Series_R_LeftEnd = 0.0, Series_L_LeftEnd = 0.0, Series_C_LeftEnd = 0.0;
        bool HasSeries_RightEnd = false;
        RealKindWires Series_R_RightEnd = 0.0, Series_L_RightEnd = 0.0, Series_C_RightEnd = 0.0;
        std::vector<WireDispersiveParams_t> disp_LeftEnd;
        std::vector<WireDispersiveParams_t> disp_RightEnd;
    };

    struct Lumped_t {
        int32_t Orient = 0;
        RealKindWires R = 0.0, L = 0.0, C = 0.0;
        RealKindWires DiodB = 0.0, DiodIsat = 0.0;
        RealKindWires Rtime_on = 0.0, Rtime_off = 0.0;
        bool resistor = false;
        bool inductor = false;
        bool capacitor = false;
        bool diodo = false;
        RealKindWires R_devia = 0.0, L_devia = 0.0, C_devia = 0.0;
    };

    struct PMLbody_t {
        int32_t orient = 0;
    };

    struct Multiport_t {
        int32_t Multiportdir = 0;
        std::string multiportFileZ11, multiportFileZ22, multiportFileZ12, multiportFileZ21;
        std::vector<RealKind> epr, mur, sigma, sigmam, width;
        std::vector<RealKind> epr_devia, mur_devia, sigma_devia, sigmam_devia, width_devia;
        int32_t numcapas = 0;
    };

    struct AnisMultiport_t {
        int32_t Multiportdir = 0;
        std::string MultiportFileZ11, MultiportFileZ22, MultiportFileZ12, MultiportFileZ21;
        std::vector<RealKind> epr, mur, sigma, sigmam, width;
    };

    struct planeonde_t {
        RealKind INCERTMAX = 0.0;
        std::vector<RealKind> px, py, pz, ex, ey, ez, incert;
        int32_t esqx1 = 0, esqy1 = 0, esqz1 = 0, esqx2 = 0, esqy2 = 0, esqz2 = 0;
        fichevol_t Fichero;
        int32_t nummodes = 0;
        bool isRC = false;
    };

    struct Border_t {
        bool IsBackPEC = false, IsFrontPEC = false, IsLeftPEC = false, IsRightPEC = false, IsUpPEC = false, IsDownPEC = false;
        bool IsBackPMC = false, IsFrontPMC = false, IsLeftPMC = false, IsRightPMC = false, IsUpPMC = false, IsDownPMC = false;
        bool IsBackPML = false, IsFrontPML = false, IsLeftPML = false, IsRightPML = false, IsUpPML = false, IsDownPML = false;
        bool IsBackPeriodic = false, IsFrontPeriodic = false, IsLeftPeriodic = false, IsRightPeriodic = false, IsUpPeriodic = false, IsDownPeriodic = false;
        bool IsBackMUR = false, IsFrontMUR = false, IsLeftMUR = false, IsRightMUR = false, IsUpMUR = false, IsDownMUR = false;
    };

    struct direction_t {
        int32_t x = 0, y = 0, z = 0, orientation = 0;

        bool operator==(const direction_t& b) const {
            return x == b.x && y == b.y && z == b.z && orientation == b.orientation;
        }
    };

    struct observable_t {
        int32_t XI = 0, YI = 0, ZI = 0, XE = 0, YE = 0, ZE = 0, What = 0, Node = 0;
        int32_t Xtrancos = 0, Ytrancos = 0, Ztrancos = 0;
        std::vector<direction_t> line;
    };

    struct Obses_t {
        int32_t nP = 0;
        std::vector<observable_t> P;
        RealKind InitialTime = 0.0, FinalTime = 0.0, TimeStep = 0.0;
        RealKind InitialFreq = 0.0, FinalFreq = 0.0, FreqStep = 0.0;
        RealKind thetaStart = 0.0, thetaStop = 0.0, thetaStep = 0.0;
        RealKind phiStart = 0.0, phiStop = 0.0, phiStep = 0.0;
        std::string outputrequest;
        std::string FileNormalize;
        bool FreqDomain = false, TimeDomain = false, Saveall = false;
        bool TransFer = false, Volumic = false, Done = false, Begun = false, Flushed = false;
    };

    struct SharedElement_t {
        int32_t i = 0, j = 0, k = 0, field = 0, PropMed = 0, SharedMed = 0, times = 0;
    };

    struct Shared_t {
        int32_t Conta = 0, MaxConta = 10;
        std::vector<SharedElement_t> elem;
    };

    struct DispersiveParams_t {
        int32_t NumPolRes11 = 0, NumPolRes12 = 0, NumPolRes13 = 0, NumPolRes22 = 0, NumPolRes23 = 0, NumPolRes33 = 0;
        std::vector<ComplexKind> C11, A11, C12, A12, C13, A13, C22, A22, C23, A23, C33, A33;
        RealKind eps11 = 0.0, MU11 = 0.0, SIGMA11 = 0.0, SIGMAM11 = 0.0;
        RealKind eps12 = 0.0, MU12 = 0.0, SIGMA12 = 0.0, SIGMAM12 = 0.0;
        RealKind EPs13 = 0.0, MU13 = 0.0, SIGMA13 = 0.0, SIGMAM13 = 0.0;
        RealKind EPs22 = 0.0, MU22 = 0.0, SIGMA22 = 0.0, SIGMAM22 = 0.0;
        RealKind EPs23 = 0.0, MU23 = 0.0, SIGMA23 = 0.0, SIGMAM23 = 0.0;
        RealKind EPs33 = 0.0, MU33 = 0.0, SIGMA33 = 0.0, SIGMAM33 = 0.0;
    };

    struct Anisotropic_t {
        RealKind sigma[3][3] = {};
        RealKind epr[3][3] = {};
        RealKind mur[3][3] = {};
        RealKind sigmaM[3][3] = {};
    };

    struct Exists_t {
        bool PML = false, PEC = false, ConformalPEC = false, PMC = false, ThinWire = false, Multiwire = false, SlantedWire = false;
        bool EDispersive = false, MDispersive = false, EDispersiveAnis = false, MDispersiveAnis = false;
        bool ThinSlot = false, PMLbody = false, SGBC = false, SGBCDispersive = false, Lumped = false, Lossy = false;
        bool AnisMultiport = false, Multiport = false, MultiportPadding = false, Dielectric = false, Anisotropic = false;
        bool Volume = false, Line = false, Surface = false, Needed = false, Interfase = false;
        bool already_YEEadvanced_byconformal = false, split_and_useless = false;
    };

    struct MediaData_t {
        RealKind Priority = 0.0, Epr = 0.0, Sigma = 0.0, Mur = 0.0, SigmaM = 0.0;
        bool sigmareasignado = false;
        Exists_t Is;
        std::vector<Wires_t> Wire;
        std::vector<SlantedWires_t> SlantedWire;
        std::vector<PMLbody_t> PMLbody;
        std::vector<Multiport_t> Multiport;
        std::vector<AnisMultiport_t> AnisMultiport;
        std::vector<DispersiveParams_t> EDispersive;
        std::vector<DispersiveParams_t> MDispersive;
        std::vector<Anisotropic_t> Anisotropic;
        std::vector<Lumped_t> Lumped;
#ifdef CompileWithMTLN
        std::vector<Multiwires_t> Multiwire;
#endif
    };

    struct SGGFDTDINFO_t {
        std::vector<RealKindTiempo> tiempo;
        RealKindTiempo dt = 0.0;
        std::string extraswitches;