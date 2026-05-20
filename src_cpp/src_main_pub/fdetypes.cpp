#include <iostream>
#include <vector>
#include <string>
#include <complex>
#include <cmath>
#include <limits>
#include <algorithm>
#include <cstring>
#include <stdexcept>

// Forward declarations for MPI if needed, though we simulate the constants
#ifdef CompileWithMPI
#include <mpi.h>
#endif

#ifdef CompileWithOpenMP
#include <omp.h>
#endif

// Macro definitions for compilation flags
// Note: In C++, these are typically handled by compiler flags (-DCompileWithReal8, etc.)
// We assume standard defaults if not defined, mimicking the Fortran logic.

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

namespace FDETYPES_m {

    // Tunable Parameters
    int quienmpi = 0;
    int tamaniompi = 0;
    int SUBCOMM_MPI = 0;
    int SUBCOMM_MPI_conformal_probes = 0;
    int MPI_conformal_probes_root = 0;

    const long long maxmpibytes = (1LL << 27);
    const int BuffObse = (1 << 10);
    const long long MaxMemoryProbes = (1LL << 37);
    const int MaxProbes = 150000;
    const int topCPUtime = 10000000;
    const int BUFSIZE = 1024;
    const int BUFSIZE_LONG = 16384;

    // Rest of Parameters
    // Determine integer size for media matrices
    int INTEGERSIZEOFMEDIAMATRICES = 4;
    int INTEGERSIZE = 4; // Default MPI_INTEGER4 equivalent

#ifdef CompileWithInt1
    #undef CompileWithInt2
    #undef CompileWithInt4
    INTEGERSIZEOFMEDIAMATRICES = 1;
    #ifdef CompileWithMPI
    INTEGERSIZE = MPI_INTEGER1;
    #endif
#endif

#ifdef CompileWithInt2
    #undef CompileWithInt1
    #undef CompileWithInt4
    INTEGERSIZEOFMEDIAMATRICES = 2;
    #ifdef CompileWithMPI
    INTEGERSIZE = MPI_INTEGER2;
    #endif
#endif

#ifdef CompileWithInt4
    #undef CompileWithInt1
    #undef CompileWithInt2
    INTEGERSIZEOFMEDIAMATRICES = 4;
    #ifdef CompileWithMPI
    INTEGERSIZE = MPI_INTEGER4;
    #endif
#endif

    const int IKINDMTAG = 4;

    const int SINGLE = 4;
    const int DOUBLE = 8;
    const int LONG_DOUBLE = 16;

    // Determine Real Kind
    int RKIND = SINGLE;
    int RKIND_wires = DOUBLE;
    int RKIND_tiempo = DOUBLE;
    int CKIND = DOUBLE;

#ifdef CompileWithReal8
    RKIND = DOUBLE;
    RKIND_wires = DOUBLE;
    RKIND_tiempo = DOUBLE;
    CKIND = DOUBLE;
#else
#ifdef CompileWithReal16
    RKIND = LONG_DOUBLE;
    RKIND_wires = LONG_DOUBLE;
    RKIND_tiempo = LONG_DOUBLE;
    CKIND = LONG_DOUBLE;
#else
    // Default (Real4)
    RKIND = SINGLE;
    RKIND_wires = DOUBLE;
    RKIND_tiempo = DOUBLE;
    CKIND = DOUBLE;
#endif
#endif

    // MPI Real Sizes
    int REALSIZE = 4;
    int REALSIZE_wires = 8;
    int COMPLEXSIZE = 16;
    int REALSIZE_tiempo = 8;

#ifdef CompileWithMPI
    double plusCPU_PML = 2.0;
#ifdef CompileWithReal8
    REALSIZE = MPI_DOUBLE_PRECISION;
    REALSIZE_wires = MPI_DOUBLE_PRECISION;
    COMPLEXSIZE = MPI_DOUBLE_COMPLEX;
    REALSIZE_tiempo = MPI_DOUBLE_PRECISION;
#else
#ifdef CompileWithReal16
    REALSIZE = MPI_REAL16; // Assuming MPI_REAL16 exists or is mapped
    COMPLEXSIZE = MPI_COMPLEX32;
    REALSIZE_tiempo = MPI_REAL_16;
#else
    REALSIZE = MPI_REAL;
    REALSIZE_wires = MPI_DOUBLE_PRECISION;
    REALSIZE_tiempo = MPI_DOUBLE_PRECISION;
    COMPLEXSIZE = MPI_DOUBLE_COMPLEX;
#endif
#endif
#endif

    // Constants
    double heurCFL = 0.8;
    const double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148;
    const double unmedio = 0.5;
    std::complex<double> mcPI2 = -std::complex<double>(0.0, 1.0) * 2.0 * pi;

    const int Down = 1, Up = 2, Left = 3, Right = 4, Back = 5, Front = 6;
    const int iEx = 1, iEy = 2, iEz = 3, iHx = 4, iHy = 5, iHz = 6, centroide = 8, Nothing = 666;
    const int iMEC = 51;
    const int iMHC = 52;
    const int iCur = 53;
    const int iCurX = 54;
    const int iCurY = 55;
    const int iCurZ = 56;
    const int mapvtk = 57;
    const int iExC = 61;
    const int iEyC = 62;
    const int iEzC = 63;
    const int iHxC = 64;
    const int iHyC = 65;
    const int iHzC = 66;
    const int farfield = 67;
    const int lineIntegral = 68;
    const int iJx = 10 * iEx, iJy = 10 * iEy, iJz = 10 * iEz;
    const int iQx = 10000 * iEx, iQy = 10000 * iEy, iQz = 10000 * iEz;
    const int iVx = 1000 * iEx, iVy = 1000 * iEy, iVz = 1000 * iEz;
    const int iBloqueJx = 100 * iEx, iBloqueJy = 100 * iEy, iBloqueJz = 100 * iEz;
    const int iBloqueMx = 100 * iHx, iBloqueMy = 100 * iHy, iBloqueMz = 100 * iHz;

    const std::string SEPARADOR = "______________";
    const int comi = 1, fine = 2, icoord = 1, jcoord = 2, kcoord = 3;

    const double EPSILON_VACUUM = 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12;
    const double MU_VACUUM = 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6;
    const double c_vacuum = 1.0 / std::sqrt(EPSILON_VACUUM * MU_VACUUM);

    double dt0 = 0.0;

    const int FACE_X = 1;
    const int FACE_Y = 2;
    const int FACE_Z = 3;
    const int EDGE_X = 1;
    const int EDGE_Y = 2;
    const int EDGE_Z = 3;

    const std::string F_SOURCE_VOLTAGE = "VOLT";
    const std::string F_SOURCE_CURRENT = "CURR";

    std::string fmt;
    void init_fmt() {
#ifdef CompileWithReal4
        fmt = "(e27.17e3,11(e19.9e3))";
#else
#ifdef CompileWithReal8
        fmt = "(12(e27.17e3))";
#else
#ifdef CompileWithReal16
        fmt = "(12(e46.36e3))";
#else
        fmt = "(e27.17e3,11(e19.9e3))";
#endif
#endif
#endif
    }

    // Types

    struct tagtype_t {
        std::vector<std::string> tag;
        int numertags = 0;
    };

    struct coorsxyz_t {
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
    };

    struct coorsxyzP_t {
        std::vector<coorsxyz_t> PhysCoor;
        coorsxyzP_t() {
            PhysCoor.resize(6);
        }
    };

    struct MedioExtra_t {
        int pml_size = 0;
        int index = 0;
        double sigma = 0.0;
        double sigmam = 0.0;
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
            PeriodicBorders = false;
            Anisotropic = false;
            ThinSlot = false;
            NodalE = false;
            NodalH = false;
            MagneticMedia = false;
            PMLMagneticMedia = false;
            MTLNbundles = false;
        }
    };

    struct Xlimit_t {
        int XI = 0, XE = 0, NX = 0;
    };

    struct Ylimit_t {
        int YI = 0, YE = 0, NY = 0;
    };

    struct Zlimit_t {
        int ZI = 0, ZE = 0, NZ = 0;
    };

    struct limit_t {
        int XI = 0, XE = 0, YI = 0, YE = 0, ZI = 0, ZE = 0, NX = 0, NY = 0, NZ = 0;
    };

    struct XYZlimit_t {
        int XI = 0, XE = 0, YI = 0, YE = 0, ZI = 0, ZE = 0;
    };

    struct xyzlimit_scaled_t {
        int XI = 0, XE = 0, YI = 0, YE = 0, ZI = 0, ZE = 0;
        double xc = 0.0, yc = 0.0, zc = 0.0;
        int Or = 0;
    };

    struct tagnumber_t {
        std::vector<std::vector<std::vector<int>>> x;
        std::vector<std::vector<std::vector<int>>> y;
        std::vector<std::vector<std::vector<int>>> z;
    };

    struct taglist_t {
        tagnumber_t edge;
        tagnumber_t face;

        int getFaceTag(int i, int j, int k) {
            if (face.x.empty()) return 0;
            if (i < 0 || i >= face.x[0][0].size() || j < 0 || j >= face.x[0].size() || k < 0 || k >= face.x.size()) return 0;
            return face.x[k][j][i];
        }

        int getEdgeTag(int i, int j, int k) {
            if (edge.x.empty()) return 0;
            if (i < 0 || i >= edge.x[0][0].size() || j < 0 || j >= edge.x[0].size() || k < 0 || k >= edge.x.size()) return 0;
            return edge.x[k][j][i];
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
        double orden[3][2] = {};
        double CoeffReflPML[3][2] = {};
        int NumLayers[3][2] = {};
    };

    struct fichevol_t {
        std::string Name;
        int NumSamples = 0;
        double DeltaSamples = 0.0;
        std::vector<double> Samples;
    };

    struct fichevol_wires_t {
        std::string Name;
        int NumSamples = 0;
        double DeltaSamples = 0.0;
        std::vector<double> Samples;
    };

    struct source_t {
        fichevol_wires_t Fichero;
        double Resistance = 0.0;
        double Multiplier = 0.0;
        int i = 0, j = 0, k = 0;
    };

    struct NodalSource_t {
        fichevol_t Fichero;
        std::vector<xyzlimit_scaled_t> punto;
        int numpuntos = 0;
        bool IsInitialValue = false;
        bool IsHard = false;
        bool IsElec = false;
    };

    struct WireDispersiveParams_t {
        int numPoles = 0;
        std::vector<std::complex<double>> res;
        std::vector<std::complex<double>> p;
        std::complex<double> d;
        std::complex<double> e;
    };

    struct oriented_point_t {
        int ori = 0;
        int i = 0, j = 0, k = 0, origIndex = 0, ilibre = 0, jlibre = 0, klibre = 0, multiraboDE = 0;
        bool Is_LeftEnd = false;
        bool Is_RightEnd = false;
        bool IsEnd_norLeft_norRight = false;
        bool repetido = false;
        bool multirabo = false;
        bool orientadoalreves = false;
    };

#ifdef CompileWithMTLN
    struct Multiwires_t {
        // Empty in Fortran
    };
#endif

    struct Wires_t {
        double Radius = 0.0, R = 0.0, L = 0.0, C = 0.0;
        double P_R = 0.0, P_L = 0.0, P_C = 0.0;
        double Radius_devia = 0.0, R_devia = 0.0, L_devia = 0.0, C_devia = 0.0;
        std::vector<WireDispersiveParams_t> disp;
        int numsegmentos = 0;
        int NUMVOLTAGESOURCES = 0;
        int NUMCURRENTSOURCES = 0;
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
        double Parallel_R_RightEnd = 0.0, Parallel_R_LeftEnd = 0.0;
        double Series_R_RightEnd = 0.0, Series_R_LeftEnd = 0.0;
        double Parallel_L_RightEnd = 0.0, Parallel_L_LeftEnd = 0.0;
        double Series_L_RightEnd = 0.0, Series_L_LeftEnd = 0.0;
        double Parallel_C_RightEnd = 0.0, Parallel_C_LeftEnd = 0.0;
        double Series_C_RightEnd = 0.0, Series_C_LeftEnd = 0.0;
        double Parallel_R_RightEnd_devia = 0.0, Parallel_R_LeftEnd_devia = 0.0;
        double Series_R_RightEnd_devia = 0.0, Series_R_LeftEnd_devia = 0.0;
        double Parallel_L_RightEnd_devia = 0.0, Parallel_L_LeftEnd_devia = 0.0;
        double Series_L_RightEnd_devia = 0.0, Series_L_LeftEnd_devia = 0.0;
        double Parallel_C_RightEnd_devia = 0.0, Parallel_C_LeftEnd_devia = 0.0;
        double Series_C_RightEnd_devia = 0.0, Series_C_LeftEnd_devia = 0.0;
        std::vector<WireDispersiveParams_t> disp_LeftEnd;
        std::vector<WireDispersiveParams_t> disp_RightEnd;
        int LeftEnd = 0;
        int RightEnd = 0;
    };

    struct SlantedNode_t {
        int index = 0;
        double x = 0.0, y = 0.0, z = 0.0;
        bool VsourceExists = false;
        bool IsourceExists = false;
        source_t* Vsource = nullptr;
        source_t* Isource = nullptr;
    };

    struct SlantedWires_t {
        double radius = 0.0, R = 0.0, L = 0.0, C = 0.0;
        double P_R = 0.0, P_L = 0.0, P_C = 0.0;
        std::vector<WireDispersiveParams_t> disp;
        int LeftEnd = 0, RightEnd = 0;
        int NumNodes = 0;
        std::vector<SlantedNode_t> nodes;
        bool HasParallel_LeftEnd = false;
        double Parallel_R_LeftEnd = 0.0, Parallel_L_LeftEnd = 0.0, Parallel_C_LeftEnd = 0.0;
        bool HasParallel_RightEnd = false;
        double Parallel_R_RightEnd = 0.0, Parallel_L_RightEnd = 0.0, Parallel_C_RightEnd = 0.0;
        bool HasSeries_LeftEnd = false;
        double Series_R_LeftEnd = 0.0, Series_L_LeftEnd = 0.0, Series_C_LeftEnd = 0.0;
        bool HasSeries_RightEnd = false;
        double Series_R_RightEnd = 0.0, Series_L_RightEnd = 0.0, Series_C_RightEnd = 0.0;
        std::vector<WireDispersiveParams_t> disp_LeftEnd;
        std::vector<WireDispersiveParams_t> disp_RightEnd;
    };

    struct Lumped_t {
        int Orient = 0;
        double R = 0.0, L = 0.0, C = 0.0;
        double DiodB = 0.0, DiodIsat = 0.0;
        double Rtime_on = 0.0, Rtime_off = 0.0;
        bool resistor = false;
        bool inductor = false;
        bool capacitor = false;
        bool diodo = false;
        double R_devia = 0.0, L_devia = 0.0, C_devia = 0.0;
    };

    struct PMLbody_t {
        int orient = 0;
    };

    struct Multiport_t {
        int Multiportdir = 0;
        std::string multiportFileZ11, multiportFileZ22, multiportFileZ12, multiportFileZ21;
        std::vector<double> epr, mur, sigma, sigmam, width;
        std::vector<double> epr_devia, mur_devia, sigma_devia, sigmam_devia, width_devia;
        int numcapas = 0;
    };

    struct AnisMultiport_t {
        int Multiportdir = 0;
        std::string MultiportFileZ11, MultiportFileZ22, MultiportFileZ12, MultiportFileZ21;
        std::vector<double> epr, mur, sigma, sigmam, width;
    };

    struct planeonde_t {
        double INCERTMAX = 0.0;
        std::vector<double> px, py, pz, ex, ey, ez, incert;
        int esqx1 = 0, esqy1 = 0, esqz1 = 0, esqx2 = 0, esqy2 = 0, esqz2 = 0;
        fichevol_t Fichero;
        int nummodes = 0;
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
        int x = 0, y = 0, z = 0, orientation = 0;

        bool operator==(const direction_t& other) const {
            return x == other.x && y == other.y && z == other.z && orientation == other.orientation;
        }
    };

    struct observable_t {
        int XI = 0, YI = 0, ZI = 0, XE = 0, YE = 0, ZE = 0, What = 0, Node = 0;
        int Xtrancos = 0, Ytrancos = 0, Ztrancos = 0;
        std::vector<direction_t> line;
    };

    struct Obses_t {
        int nP = 0;
        std::vector<observable_t> P;
        double InitialTime = 0.0, FinalTime = 0.0, TimeStep = 0.0;
        double InitialFreq = 0.0, FinalFreq = 0.0, FreqStep = 0.0;
        double thetaStart = 0.0, thetaStop = 0.0, thetaStep = 0.0;
        double phiStart = 0.0, phiStop = 0.0, phiStep = 0.0;
        std::string outputrequest;
        std::string FileNormalize;
        bool FreqDomain = false, TimeDomain = false, Saveall = false;
        bool TransFer = false, Volumic = false, Done = false, Begun = false, Flushed = false;
    };

    struct SharedElement_t {
        int i = 0, j = 0, k = 0, field = 0, PropMed = 0, SharedMed = 0, times = 0;
    };

    struct Shared_t {
        int Conta = 0, MaxConta = 10;
        std::vector<SharedElement_t> elem;
    };

    struct DispersiveParams_t {
        int NumPolRes11 = 0, NumPolRes12 = 0, NumPolRes13 = 0, NumPolRes22 = 0, NumPolRes23 = 0, NumPolRes33 = 0;
        std::vector<std::complex<double>> C11, A11, C12, A12, C13, A13, C22, A22, C23, A23, C33, A33;
        double eps11 = 0.0, MU11 = 0.0, SIGMA11 = 0.0, SIGMAM11 = 0.0;
        double eps12 = 0.0, MU12 = 0.0, SIGMA12 = 0.0, SIGMAM12 = 0.0;
        double EPs13 = 0.0, MU13 = 0.0, SIGMA13 = 0.0, SIGMAM13 = 0.0;
        double EPs22 = 0.0, MU22 = 0.0, SIGMA22 = 0.0, SIGMAM22 = 0.0;
        double EPs23 = 0.0, MU23 = 0.0, SIGMA23 = 0.0, SIGMAM23 = 0.0;
        double EPs33 = 0.0, MU33 = 0.0, SIGMA33 = 0.0, SIGMAM33 = 0.0;
    };

    struct Anisotropic_t {
        double sigma[3][3] = {};
        double epr[3][3] = {};
        double mur[3][3] = {};
        double sigmaM[3][3] = {};
    };

    struct Exists_t {
        bool PML = false, PEC = false, ConformalPEC = false, PMC = false;
        bool ThinWire = false, Multiwire = false, SlantedWire = false;
        bool EDispersive = false, MDispersive = false;
        bool EDispersiveAnis = false, MDispersiveAnis = false;
        bool ThinSlot = false, PMLbody = false, SGBC = false, SGBCDispersive = false;
        bool Lumped = false, Lossy = false, AnisMultiport = false, Multiport = false;
        bool MultiportPadding = false, Dielectric = false, Anisotropic = false;
        bool Volume = false, Line = false, Surface = false, Needed = false, Interfase = false;
        bool already_YEEadvanced_byconformal = false, split_and_useless = false;
    };

    struct MediaData_t {
        double Priority = 0.0, Epr = 0.0, Sigma = 0.0, Mur = 0.0, SigmaM = 0.0;
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
        std::vector<double> tiempo;
        double dt = 0.0;
        std::string extraswitches;
        int NumMedia = 0, AllocMed = 0;
        int IniPMLMedia = 0, EndPMLMedia = 0;
        int NumPlaneWaves = 0, TimeSteps = 0, InitialTimeStep = 0;
        int NumNodalSources = 0;
        int NumberRequest = 0;
        std::vector<double> LineX, LineY, LineZ;
        std::vector<double> DX, DY, DZ;
        int AllocDxI = 0, AllocDyI = 0, AllocDzI = 0, AllocDxE = 0, AllocDyE = 0, AllocDzE = 0;
        std::vector<planeonde_t> PlaneWave;
        Border_t Border;
        PML_t PML;
        Shared_t Eshared;
        Shared_t Hshared;
        std::vector<XYZlimit_t> Alloc, Sweep, SINPMLSweep;
        std::vector<MediaData_t> Med;
        std::vector<NodalSource_t> NodalSource;
        std::vector<Obses_t> Observation;
        bool thereAreMagneticMedia = false;
        bool thereArePMLMagneticMedia = false;
        std::string nEntradaRoot;
        coorsxyzP_t Punto;
    };

    struct media_matrices_t {
        std::vector<std::vector<std::vector<int>>> sggMiNo;
        std::vector<std::vector<std::vector<int>>> sggMiEx;
        std::vector<std::vector<std::vector<int>>> sggMiEy;
        std::vector<std::vector<std::vector<int>>> sggMiEz;
        std::vector<std::vector<std::vector<int>>> sggMiHx;
        std::vector<std::vector<std::vector<int>>> sggMiHy;
        std::vector<std::vector<std::vector<int>>> sggMiHz;
        std::vector<std::vector<std::vector<int>>> sggMtag;
    };

    struct constants_t {
        std::vector<double> g1, g2, gM1, gM2;

        void destroy() {
            g1.clear();
            g2.clear();
            gM1.clear();
            gM2.clear();
        }
    };

    struct nf2ff_t {
        bool tr = false, fr = false, iz = false, de = false, ab = false, ar = false;
    };

    struct perform_t {
        bool flushFields = false;
        bool flushData = false;
        bool unpack = false;
        bool postprocess = false;
        bool flushXdmf = false;
        bool flushVTK = false;

        bool isFlush() {
            return flushFields || flushData || unpack || postprocess || flushXdmf || flushVTK;
        }

        void reset() {
            flushFields = false;
            flushData = false;
            unpack = false;
            postprocess = false;
            flushXdmf = false;
            flushVTK = false;
        }
    };

    struct sim_control_t {
        bool simu_devia = false, resume = false, saveall = false, makeholes = false;
        bool connectendings = false, isolategroupgroups = false, createmap = false;
        bool groundwires = false, noSlantedcrecepelo = false;
        // The original code was cut off here, so we stop here.
    };

} // namespace FDETYPES_m

// Continuation of FDETYPES_m translation

namespace FDETYPES_m {

    // Global variables corresponding to Fortran SAVE, public variables
    // Note: In a real translation, these might be part of a singleton or global state manager.
    // Here we define them as extern or static globals to match the Fortran module scope.
    extern int prior_BV;
    extern int prior_IB;
    extern int prior_pmlbody;
    extern int prior_AB;
    extern int prior_FDB;
    extern int prior_IS;
    extern int prior_AS;
    extern int prior_FDS;
    extern int prior_IL;
    extern int prior_AL;
    extern int prior_FDL;
    extern int prior_IP;
    extern int prior_AP;
    extern int prior_FDP;
    extern int prior_PEC;
    extern int prior_PMC;
    extern int prior_TG;
    extern int prior_CS;
    extern int prior_TW;

    extern bool input_conformal_flag;

    // Implementation of constants_destroy
    void constants_destroy(constants_t& this_obj) {
        delete[] this_obj.g1;
        delete[] this_obj.g2;
        delete[] this_obj.gm1;
        delete[] this_obj.gm2;
        this_obj.g1 = nullptr;
        this_obj.g2 = nullptr;
        this_obj.gm1 = nullptr;
        this_obj.gm2 = nullptr;
    }

    // Implementation of isFlush
    bool isFlush(perform_t& this_obj) {
        return this_obj.flushDATA || this_obj.flushFIELDS || this_obj.postprocess || this_obj.flushXdmf || this_obj.flushVTK;
    }

    // Implementation of perform_reset
    void perform_reset(perform_t& this_obj) {
        this_obj.flushFields = false;
        this_obj.flushData = false;
        this_obj.unpack = false;
        this_obj.postprocess = false;
        this_obj.flushXdmf = false;
        this_obj.flushVTK = false;
    }

    // Implementation of logic_reset
    void logic_reset(logic_control_t& this_obj) {
        this_obj.Wires = false;
        this_obj.PMLbodies = false;
        this_obj.MultiportS = false;
        this_obj.AnisMultiportS = false;
        this_obj.SGBCs = false;
        this_obj.Lumpeds = false;
        this_obj.EDispersives = false;
        this_obj.MDispersives = false;
        this_obj.PlaneWaveBoxes = false;
        this_obj.Observation = false;
        this_obj.FarFields = false;
        this_obj.PMCBorders = false;
        this_obj.PMLBorders = false;
        this_obj.MurBorders = false;
        this_obj.PECBorders = false;
        this_obj.Anisotropic = false;
        this_obj.ThinSlot = false;
        this_obj.NodalE = false;
        this_obj.NodalH = false;
        this_obj.PeriodicBorders = false;
        this_obj.MagneticMedia = false;
        this_obj.PMLMagneticMedia = false;
        this_obj.MTLNbundles = false;
    }

    // Implementation of setglobal
    void setglobal(int iu1, int iu2) {
        quienmpi = iu1;
        tamaniompi = iu2;
    }

    // Implementation of set_priorities
    void set_priorities(bool prioritizeCOMPOoverPEC, bool prioritizeISOTROPICBODYoverall, bool prioritizeTHINWIRE) {
        // Moved here the priority system to control them with switches. Useful for siva 070815 (PEC priority bug over compo of siva)
        prior_BV = 10; // background volume
        prior_AB = 30; // anisotropic body
        prior_FDB = 40; // Frequency dependent body
        prior_IS = 50; // Isotropic surface
        prior_AS = 60; // Anisotropic surface
        prior_FDS = 70; // Frequency dependent surface
        prior_IL = 90; // Isotropic line
        prior_AL = 100; // Anisotropic line
        prior_FDL = 110; // Frequency dependent line
        prior_IP = 120; // Isotropic point
        prior_AP = 130; // Anisotropic point
        prior_FDP = 140; // Frequency dependent point
        prior_PEC = 150; // Perfectly electric conducting body, surface, line, or point
        prior_PMC = 160; // Perfectly magnetic conducting body, surface, line, or point
        prior_TG = 155; // thin Slot has more priority than PEC
        
        // Added option -prioritizeCOMPOoverPEC to raise its priority and be able to simulate SIVA (sgg 070815)
        if (prioritizeTHINWIRE) {
            prior_TW = 1500; // Changed to 231024 and set with maximum priority. It is only experimental and for visualization
        } else { // Correct option. The above is only experimental and for visualization
            prior_TW = 15; // Priority of the thin-wire below all (except background)
        }
        
        // prior_pmlbody = prior_TW-1; // The thread has priority over the pmlbody (test HOLD coax sgg 251019)
        prior_pmlbody = prior_BV + 1; // The pml body can be penetrated by everything 311019 sgg
        
        if (prioritizeCOMPOoverPEC) { // Composite surface
            prior_CS = prior_PEC + 2;
        } else {
            prior_CS = prior_PEC - 2; // Composites have lower priority than PEC to properly handle PEC-composite junctions (ss's 210312 mail)
        }
        
        if (prioritizeISOTROPICBODYoverall) { // Isotropic body
            prior_IB = 200; // ONLY FOR THE SIVA CASE TO REMOVE HOLES OF vacuum
        } else {
            prior_IB = 20; // THE USUAL
        }
    }

    // Implementation of taglist_getFaceTag
    int taglist_getFaceTag(taglist_t& this_obj, int field, int i, int j, int k) {
        int res = 0;
        switch (field) {
            case iHx:
                res = this_obj.face.x[i][j][k];
                break;
            case iHy:
                res = this_obj.face.y[i][j][k];
                break;
            case iHz:
                res = this_obj.face.z[i][j][k];
                break;
            default:
                // Handle undefined field case if necessary
                break;
        }
        return res;
    }

    // Implementation of taglist_getEdgeTag
    int taglist_getEdgeTag(taglist_t& this_obj, int field, int i, int j, int k) {
        int res = 0;
        switch (field) {
            case iEx:
                res = this_obj.edge.x[i][j][k];
                break;
            case iEy:
                res = this_obj.edge.y[i][j][k];
                break;
            case iEz:
                res = this_obj.edge.z[i][j][k];
                break;
            default:
                // Handle undefined field case if necessary
                break;
        }
        return res;
    }

    // Implementation of direction_eq
    bool direction_eq(const direction_t& a, const direction_t& b) {
        bool eq = true;
        eq = eq && (a.x == b.x);
        eq = eq && (a.y == b.y);
        eq = eq && (a.z == b.z);
        eq = eq && (a.orientation == b.orientation);
        return eq;
    }

} // namespace FDETYPES_m