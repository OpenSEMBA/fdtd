#include "semba_fdtd.h"
#include "lumped.h"
#include "maloney_nostoch.h"
#include "mapvtk_writer.h"
#include "xdmf_h5.h"

#include <nlohmann/json.hpp>

#ifdef CompileWithMTLN
#include "smbjson_m.h"
#include "mtln_solver_m.h"
#include "wires_mtln_m.h"
#endif
#ifdef CompileWithMPI
#include <mpi.h>
#endif
#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <fstream>
#include <sstream>
#include <map>
#include <iomanip>
#include <array>
#include <filesystem>
#include <set>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <complex>
#include <cctype>
#include <stdexcept>
#include <unistd.h>
#if defined(__SSE__)
#include <xmmintrin.h>
#endif
#if defined(__SSE3__)
#include <pmmintrin.h>
#endif

#ifdef SEMBA_CPP_ENABLE_HDF5
#include <hdf5.h>
#endif

#ifdef CompileWithMPI
MPI_Comm SUBCOMM_MPI = MPI_COMM_WORLD;
MPI_Comm subcomm_mpi = MPI_COMM_WORLD;
#endif

const double PI = 3.14159265358979323846;
const double EPS0 = 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12;
const double MU0 = 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6;
const double C0 = 299792458.0;
const double ZVAC = std::sqrt(MU0/EPS0);
const double heurCFL = 0.8;
#ifndef BUFSIZE
#define BUFSIZE 1024
#endif

#ifdef CompileWithReal8
using fdtd_real = double;
#else
using fdtd_real = float;
#endif

#if defined(__GNUC__)
#define SEMBA_FORTRAN_ROUNDING __attribute__((noinline,optimize("no-fast-math")))
#define SEMBA_FORTRAN_INLINE_ROUNDING inline __attribute__((always_inline,optimize("no-fast-math")))
#define SEMBA_RESTRICT __restrict__
#else
#define SEMBA_FORTRAN_ROUNDING
#define SEMBA_FORTRAN_INLINE_ROUNDING inline
#define SEMBA_RESTRICT
#endif

namespace {

using MpiSliceInfo = SEMBA_FDTD_m::SEMBA_FDTD_test::MpiSliceInfo;

std::string lowercaseToken(std::string value) {
    for (char& ch : value) {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
    return value;
}

int mpiAxisFromToken(const std::string& token) {
    const std::string value = lowercaseToken(token);
    if (value == "x" || value == "1") return 1;
    if (value == "y" || value == "2") return 2;
    if (value == "z" || value == "3") return 3;
    throw std::runtime_error("Invalid -mpidir option: " + token);
}

int mpiAxisFromFlagsLocal(const std::string& flags) {
    std::istringstream iss(flags);
    std::string token;
    int axis = 3;
    bool seen = false;
    while (iss >> token) {
        if (token == "-mpidir") {
            std::string value;
            if (!(iss >> value)) {
                throw std::runtime_error("Missing value after -mpidir");
            }
            const int parsed = mpiAxisFromToken(value);
            if (seen && parsed != axis) {
                throw std::runtime_error("Duplicate incoherent -mpidir option");
            }
            axis = parsed;
            seen = true;
        } else if (token.rfind("-mpidir=", 0) == 0) {
            const int parsed = mpiAxisFromToken(token.substr(8));
            if (seen && parsed != axis) {
                throw std::runtime_error("Duplicate incoherent -mpidir option");
            }
            axis = parsed;
            seen = true;
        }
    }
    return axis;
}

std::vector<MpiSliceInfo> buildMpiOneAxisSlicesLocal(int cells,
                                                    int ranks,
                                                    int pml_down_layers,
                                                    int pml_up_layers,
                                                    int forced_cut,
                                                    int axis) {
    if (cells <= 0) {
        throw std::runtime_error("MPI slicing requires a positive cell count");
    }
    if (ranks <= 0) {
        throw std::runtime_error("MPI slicing requires at least one rank");
    }
    axis = (axis >= 1 && axis <= 3) ? axis : 3;

    std::vector<MpiSliceInfo> slices(static_cast<size_t>(ranks));
    if (ranks == 1) {
        slices[0].rank = 0;
        slices[0].ranks = 1;
        slices[0].axis = axis;
        slices[0].com = 0;
        slices[0].fin = cells;
        slices[0].sweepZI = 0;
        slices[0].sweepZE = cells;
        slices[0].allocZI = -1;
        slices[0].allocZE = cells + 1;
        slices[0].physicalDown = true;
        slices[0].physicalUp = true;
        slices[0].pmlDown = pml_down_layers > 0;
        slices[0].pmlUp = pml_up_layers > 0;
        return slices;
    }

    if (forced_cut >= 0 && ranks != 2) {
        throw std::runtime_error("Forced MPI cuts are only supported for two ranks");
    }
    if (forced_cut >= 0 && (forced_cut <= 0 || forced_cut >= cells)) {
        throw std::runtime_error("Forced MPI cut is outside the domain");
    }

    constexpr double plusCPU_PML_local = 2.0;
    const double fullZI = 0.0;
    const double fullZE = static_cast<double>(cells);
    const double sinpmlZI = static_cast<double>(std::max(0, pml_down_layers));
    const double sinpmlZE = static_cast<double>(std::max(0, cells - std::max(0, pml_up_layers)));

    std::vector<double> cZI(static_cast<size_t>(ranks) + 1, 0.0);
    std::vector<double> cZE(static_cast<size_t>(ranks), 0.0);
    std::vector<int> trancos(static_cast<size_t>(ranks), 0);

    const double carga =
        (fullZE - fullZI) / static_cast<double>(ranks) +
        (plusCPU_PML_local - 1.0) *
            ((sinpmlZI - fullZI) + (fullZE - sinpmlZE)) /
            static_cast<double>(ranks);
    cZI[0] = fullZI;
    for (int ilay = 0; ilay < ranks; ++ilay) {
        const double guess =
            carga + cZI[static_cast<size_t>(ilay)] +
            (plusCPU_PML_local - 1.0) *
                (std::min(cZI[static_cast<size_t>(ilay)], sinpmlZI) +
                 std::max(cZI[static_cast<size_t>(ilay)], sinpmlZE));
        const double zeCandidates[3] = {
            (guess - (plusCPU_PML_local - 1.0) * sinpmlZI) /
                (1.0 + (plusCPU_PML_local - 1.0)),
            (guess - (plusCPU_PML_local - 1.0) * sinpmlZE) /
                (1.0 + (plusCPU_PML_local - 1.0)),
            guess - (plusCPU_PML_local - 1.0) * sinpmlZE -
                (plusCPU_PML_local - 1.0) * sinpmlZI,
        };
        double bestError = std::numeric_limits<double>::infinity();
        double bestZE = zeCandidates[0];
        for (double ze : zeCandidates) {
            const double weighted =
                (ze - cZI[static_cast<size_t>(ilay)]) +
                (plusCPU_PML_local - 1.0) *
                    ((std::min(sinpmlZI, ze) -
                      std::min(sinpmlZI, cZI[static_cast<size_t>(ilay)])) +
                     (std::max(sinpmlZE, ze) -
                      std::max(sinpmlZE, cZI[static_cast<size_t>(ilay)])));
            const double error = std::abs(weighted - carga);
            if (error < bestError) {
                bestError = error;
                bestZE = ze;
            }
        }
        cZE[static_cast<size_t>(ilay)] = bestZE;
        cZI[static_cast<size_t>(ilay) + 1] = bestZE;
    }

    if (forced_cut >= 0) {
        cZI[0] = fullZI;
        cZE[0] = static_cast<double>(forced_cut);
        cZI[1] = cZE[0];
        cZE[1] = fullZE;
    }

    for (int ilay = 0; ilay < ranks; ++ilay) {
        cZE[static_cast<size_t>(ilay)] = static_cast<double>(
            std::lround(cZE[static_cast<size_t>(ilay)]));
        cZI[static_cast<size_t>(ilay) + 1] = cZE[static_cast<size_t>(ilay)];
        trancos[static_cast<size_t>(ilay)] =
            static_cast<int>(cZE[static_cast<size_t>(ilay)] - fullZI);
    }

    const int minSlice = [&]() {
        int value = cells;
        int previous = 0;
        for (int ilay = 0; ilay < ranks; ++ilay) {
            const int end = (ilay == ranks - 1) ? cells : trancos[static_cast<size_t>(ilay)];
            value = std::min(value, end - previous);
            previous = end;
        }
        return value;
    }();
    if (minSlice <= 2) {
        throw std::runtime_error("Number of cells per processor less than 2");
    }
    const int maxOriginalPml = std::max(std::max(0, pml_down_layers),
                                        std::max(0, pml_up_layers));
    if (maxOriginalPml > 0 && minSlice <= maxOriginalPml) {
        throw std::runtime_error("Minimum slice size must be larger than PML layers");
    }

    for (int rank = 0; rank < ranks; ++rank) {
        MpiSliceInfo info;
        info.rank = rank;
        info.ranks = ranks;
        info.axis = axis;
        if (rank == 0) {
            info.com = 0;
            info.fin = trancos[0];
        } else if (rank == ranks - 1) {
            info.com = trancos[static_cast<size_t>(rank) - 1];
            info.fin = cells;
        } else {
            info.com = trancos[static_cast<size_t>(rank) - 1];
            info.fin = trancos[static_cast<size_t>(rank)];
        }
        info.sweepZI = info.com;
        info.sweepZE = (rank == ranks - 1) ? info.fin : info.fin - 1;
        info.allocZI = info.sweepZI - 1;
        info.allocZE = info.sweepZE + 1;
        info.physicalDown = (rank == 0);
        info.physicalUp = (rank == ranks - 1);
        info.pmlDown = info.physicalDown && pml_down_layers > 0;
        info.pmlUp = info.physicalUp && pml_up_layers > 0;
        if (rank > 0 && rank < ranks - 1) {
            info.pmlDown = info.sweepZI < pml_down_layers;
            info.pmlUp = info.sweepZE > cells - pml_up_layers;
        }
        slices[static_cast<size_t>(rank)] = info;
    }

    return slices;
}

} // namespace

void preserveFortranSubnormalArithmetic() {
#if defined(__SSE__)
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_OFF);
#endif
#if defined(__SSE3__) && defined(_MM_DENORMALS_ZERO_MASK)
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_OFF);
#endif
}

fdtd_real flushFortranSubnormal(fdtd_real value) {
    if (value != static_cast<fdtd_real>(0.0) &&
        std::abs(value) < std::numeric_limits<fdtd_real>::min()) {
        return static_cast<fdtd_real>(0.0);
    }
    return value;
}

SEMBA_FORTRAN_ROUNDING fdtd_real fortranScalarGridInverse(fdtd_real value) {
#ifdef CompileWithReal8
    return static_cast<fdtd_real>(1.0) / value;
#else
    return static_cast<fdtd_real>(static_cast<fdtd_real>(1.0) / value);
#endif
}

SEMBA_FORTRAN_ROUNDING fdtd_real fortranGridInverse(fdtd_real value) {
#ifdef CompileWithReal8
    return static_cast<fdtd_real>(1.0) / value;
#elif defined(__SSE__)
    const __m128 v = _mm_set1_ps(value);
    const __m128 r = _mm_rcp_ps(v);
    const __m128 correction = _mm_mul_ps(_mm_mul_ps(v, r), r);
    return _mm_cvtss_f32(_mm_sub_ps(_mm_add_ps(r, r), correction));
#else
    return fortranScalarGridInverse(value);
#endif
}

fdtd_real fortranPlanewaveGridInverse(fdtd_real value) {
    return fortranGridInverse(value);
}

double fortranWireStep(double step) {
    if (step == 0.0) return step;
    const double absStep = std::abs(step);
    std::ostringstream roundedStream;
    roundedStream << std::setprecision(12) << step;
    const double rounded = std::stod(roundedStream.str());
    const double tolerance =
        64.0 * std::numeric_limits<double>::epsilon() *
        std::max(absStep, std::abs(rounded));
    return std::abs(rounded - step) <= tolerance ? rounded : step;
}

fdtd_real fortranPlanewaveCluz(fdtd_real eps, fdtd_real mu) {
#ifdef CompileWithReal8
    return static_cast<fdtd_real>(1.0) / std::sqrt(eps * mu);
#else
    (void)eps;
    (void)mu;
    return static_cast<fdtd_real>(C0);
#endif
}

SEMBA_FORTRAN_INLINE_ROUNDING fdtd_real fortranRoundedMul(fdtd_real lhs,
                                                          fdtd_real rhs) {
    volatile fdtd_real result = static_cast<fdtd_real>(lhs * rhs);
    return result;
}

SEMBA_FORTRAN_INLINE_ROUNDING fdtd_real fortranRoundedAdd(fdtd_real lhs,
                                                          fdtd_real rhs) {
    volatile fdtd_real result = static_cast<fdtd_real>(lhs + rhs);
    return result;
}

SEMBA_FORTRAN_INLINE_ROUNDING fdtd_real fortranRoundedSub(fdtd_real lhs,
                                                          fdtd_real rhs) {
    volatile fdtd_real result = static_cast<fdtd_real>(lhs - rhs);
    return result;
}

fdtd_real fortranNodalProduct(fdtd_real coeff, fdtd_real inv1,
                              fdtd_real inv2, fdtd_real amplitude,
                              fdtd_real evolution) {
    fdtd_real term = fortranRoundedMul(coeff, inv1);
    term = fortranRoundedMul(term, inv2);
    term = fortranRoundedMul(term, amplitude);
    return fortranRoundedMul(term, evolution);
}

fdtd_real fortranBulkCurrentTerm(fdtd_real lhs, fdtd_real rhs,
                                 fdtd_real delta) {
    return fortranRoundedMul(fortranRoundedSub(lhs, rhs), delta);
}

SEMBA_FORTRAN_INLINE_ROUNDING fdtd_real fortranTripleProduct(fdtd_real lhs,
                                                             fdtd_real mid,
                                                             fdtd_real rhs) {
    return fortranRoundedMul(fortranRoundedMul(lhs, mid), rhs);
}

fdtd_real fortranMurFace(fdtd_real interiorNow, fdtd_real pastInterior,
                         fdtd_real pastBoundary, fdtd_real coefficient) {
    return fortranRoundedAdd(
        pastInterior,
        fortranRoundedMul(
            coefficient, fortranRoundedSub(interiorNow, pastBoundary)));
}

SEMBA_FORTRAN_ROUNDING double fortranRoundedDoubleMul(double lhs,
                                                      double rhs) {
    volatile double result = lhs * rhs;
    return result;
}

SEMBA_FORTRAN_ROUNDING double fortranRoundedDoubleAdd(double lhs,
                                                      double rhs) {
    volatile double result = lhs + rhs;
    return result;
}

SEMBA_FORTRAN_ROUNDING double fortranRoundedDoubleSub(double lhs,
                                                      double rhs) {
    volatile double result = lhs - rhs;
    return result;
}

double fortranWireCurrentUpdate(double cte1, double current, double cte3,
                                double qplusMinus, double cte2,
                                double fieldValue) {
    const double chargeAdvanced = fortranRoundedDoubleSub(
        fortranRoundedDoubleMul(cte1, current),
        fortranRoundedDoubleMul(cte3, qplusMinus));
    return fortranRoundedDoubleAdd(
        chargeAdvanced, fortranRoundedDoubleMul(cte2, fieldValue));
}

double fortranWireChargeUpdate(double cteProp, double chargePast,
                               double ctePlain, double iPlus,
                               double iMinus) {
    return fortranRoundedDoubleSub(
        fortranRoundedDoubleMul(cteProp, chargePast),
        fortranRoundedDoubleMul(
            ctePlain, fortranRoundedDoubleSub(iPlus, iMinus)));
}

double fortranWireFieldSubtract(double fieldValue, double cte5,
                                double current) {
    return fortranRoundedDoubleSub(
        fieldValue, fortranRoundedDoubleMul(cte5, current));
}

fdtd_real fortranCurlUpdate(fdtd_real oldValue, fdtd_real coeff,
                            fdtd_real aPlus, fdtd_real aMinus,
                            fdtd_real invA, fdtd_real bPlus,
                            fdtd_real bMinus, fdtd_real invB) {
    const fdtd_real diffA = static_cast<fdtd_real>(aPlus - aMinus);
    const fdtd_real termA = static_cast<fdtd_real>(diffA * invA);
    const fdtd_real diffB = static_cast<fdtd_real>(bPlus - bMinus);
    const fdtd_real termB = static_cast<fdtd_real>(diffB * invB);
    const fdtd_real curl = static_cast<fdtd_real>(termA - termB);
    const fdtd_real scaled = static_cast<fdtd_real>(coeff * curl);
    return static_cast<fdtd_real>(oldValue + scaled);
}

struct entrada_t {
    int layoutnumber = 0; int num_procs = 1; int mpidir = 3; int ierr = 0;
    std::string extension = "", input_flags = "";
    bool thereare_stoch = false, resume = false;
};
struct tiempo_t { double start = 0.0; };
struct media_matrices_t { int NumMed = 0, totalX = 0, totalY = 0, totalZ = 0, nMedia = 0; };
struct limit_t { int XI = 0, XE = 0, YI = 0, YE = 0, ZI = 0, ZE = 0; };
struct SGGFDTDINFO_t {
    int NumberRequest = 0;
    struct { bool Volumic = false; int nP = 0; int What[10] = {}; int XI = 0, XE = 0, YI = 0, YE = 0, ZI = 0, ZE = 0; bool done = false, flushed = false, Begun = false; bool TimeDomain = false, FreqDomain = false; } observation[100];
    struct { int XI = 0, XE = 0, YI = 0, YE = 0, ZI = 0, ZE = 0; } Alloc[10];
};
struct taglist_t { int nTags = 0; };
struct tagtype_t { int nTypes = 0; };
#ifndef CompileWithMTLN
struct mtln_t { int numWires = 0; };
#endif
struct sim_control_t { int layoutnumber = 0; int num_procs = 1; bool resume = false; bool stochastic = false; };

struct Material_t {
    int id = 0; std::string name, type = "vacuum";
    double relativePermittivity = 1.0, relativePermeability = 1.0;
    double electricConductivity = 0.0, magneticConductivity = 0.0;
    double radius = 0.0, resistancePerMeter = 0.0;
    std::vector<std::vector<double>> inductancePerMeterMatrix, capacitancePerMeterMatrix;
};
struct SurfaceImpedanceMaterial_t {
    int id = 0;
    double thickness = 0.0;
    double relativePermittivity = 1.0;
    double relativePermeability = 1.0;
    double electricConductivity = 0.0;
    double magneticConductivity = 0.0;
    std::vector<SGBC_nostoch_m::SGBCLayer_t> layers;
};
struct SgbcFieldRef_t {
    int component = 0; // 0=Hx, 1=Hy, 2=Hz
    int i = 0, j = 0, k = 0;
};
struct SgbcNode_t {
    int component = 0; // 0=Ex, 1=Ey, 2=Ez
    int i = 0, j = 0, k = 0;
    int normalAxis = 0;
    SGBC_nostoch_m::SGBCSurface_t maloney;
    SgbcFieldRef_t haPlus, haMinus, hbPlus, hbMinus;
};
struct NFDEGeneral_t { int XI = 0, XE = 0, YI = 0, YE = 0, ZI = 0, ZE = 0; int NumMedia = 0; double dt = 0.0; };
struct Desplazamiento_t { int XI = 0, XE = 0, YI = 0, YE = 0, ZI = 0, ZE = 0; };
struct Mesh_t { std::vector<std::vector<double>> coordinates; std::vector<std::vector<int>> elements; };

struct source_t {
    std::string type, name, magnitudeFile, field;
    std::vector<int> elementIds;
    struct { double theta = 0.0, phi = 0.0; } direction;
    struct { double theta = 1.5708, phi = 0.0; } polarization;
};
struct probe_output_t {
    std::string name, type, field, component, domainType;
    std::vector<int> elementIds;
    std::vector<std::string> directions;
    int probeId = 0;
    int coordinateId = 0;
    int cellI = 1, cellJ = 1, cellK = 1;
    int outputCellI = 1, outputCellJ = 1, outputCellK = 1;
    std::vector<double> timeData;
    std::vector<std::vector<double>> fieldByDir;
    std::vector<std::vector<double>> incidentByDir;
};

struct BulkCurrentProbe_t {
    std::string name;
    char direction = 'z';
    int sign = 1;
    int xi = 0, yi = 0, zi = 0;
    int xe = 0, ye = 0, ze = 0;
    std::vector<double> timeData;
    std::vector<double> currentData;
};

struct NodalCurrentSegment_t {
    std::string magnitudeFile;
    char direction = 'x';
    int sign = 1;
    int xi = 0, yi = 0, zi = 0;
    int xe = 0, ye = 0, ze = 0;
};

struct HollandWireSegment_t {
    int i = 0, j = 0, k = 0;
    int direction = 3;
    int orientationSign = 1;
    int nd = 0;
    std::string wireName;
    int chargeMinus = -1;
    int chargePlus = -1;
    double radius = 0.0;
    double resistance = 0.0;
    double inductance = 0.0;
    double delta = 0.0;
    double deltaTransv1 = 0.0;
    double deltaTransv2 = 0.0;
    double lind = 0.0;
    double cte1 = 1.0;
    double cte2 = 0.0;
    double cte3 = 0.0;
    double cte5 = 0.0;
    double fractionMinus = 1.0;
    double fractionPlus = 1.0;
    bool deembedFromPec = false;
    double current = 0.0;
    double currentpast = 0.0;
    double qplus_qminus = 0.0;
};

struct HollandWireNode_t {
    int i = 0, j = 0, k = 0;
    bool isPec = false;
    double chargePresent = 0.0;
    double chargePast = 0.0;
    double ctePlain = 0.0;
    double cteProp = 1.0;
    std::vector<int> currentPlus;
    std::vector<int> currentMinus;
};

struct HollandVoltageGenerator_t {
    int segmentIndex = -1;
    std::string magnitudeFile;
    double multiplier = 1.0;
};

struct HollandWireProbe_t {
    std::string name;
    std::string wireName;
    int segmentIndex = -1;
    int cellI = 0, cellJ = 0, cellK = 0;
    int direction = 3;
    int orientationSign = 1;
    int nd = 0;
    int delaySteps = 0;
    std::vector<double> timeData;
    std::vector<double> currentData;
    std::vector<double> eTimesDlData;
    std::vector<double> vplusData;
    std::vector<double> vminusData;
    std::vector<double> vdropData;
};
struct boundary_t { std::string type = "mur"; int layers = 10; int order = 2; double reflection = 0.001; };

struct PlaneWaveState_t {
    std::vector<fdtd_real> px, py, pz;
    std::vector<fdtd_real> ex, ey, ez;
    std::vector<fdtd_real> hx, hy, hz;
    std::vector<fdtd_real> samples;
    fdtd_real deltaevol = 0.0;
    int numSamples = 0;
    fdtd_real distanciaInicial = 0.0;
    bool iluminaTr = false, iluminaFr = false, iluminaIz = false, iluminaDe = false, iluminaAr = false, iluminaAb = false;
    int esqx1 = 0, esqx2 = 0, esqy1 = 0, esqy2 = 0, esqz1 = 0, esqz2 = 0;
};

struct Parseador_t {
    int switches = 0;
    struct { int XI = 0, XE = 0, YI = 0, YE = 0, ZI = 0, ZE = 0; int NumMedia = 0; double dt = 0.0; std::string additionalArguments; } general;
    struct { int nMedia = 0; int totalX = 0, totalY = 0, totalZ = 0; } matriz;
    struct { int XI = 0, XE = 0, YI = 0, YE = 0, ZI = 0, ZE = 0; } despl;
    struct { int type = 0; } front;
    struct { int nMaterials = 0; } Mats;
    struct { int nPEC = 0; } pecRegs, pmcRegs;
    struct { int nDielectric = 0; } dielRegs;
    struct { std::vector<int> volumes; std::vector<int> surfaces; } conformalRegs;
    struct { std::vector<source_t> planeWaves; std::vector<source_t> nodalSources; } sources;
    struct { std::vector<probe_output_t> probes; } probes;
    struct { std::vector<Material_t> materials; } materials;
    struct { std::vector<std::vector<int>> associations; std::vector<int> materialIds; } matAssoc;
    struct { std::vector<boundary_t> boundaries; } boundaries;
    struct { std::vector<double> cellStepsX, cellStepsY, cellStepsZ; } cellSteps;
    struct { std::vector<std::vector<double>> intervals; } elements;
};

std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
    return (a == std::string::npos) ? "" : s.substr(a, b - a + 1);
}
std::string adjustl(const std::string& s) { return trim(s); }
std::string to_lower(const std::string& s) { std::string r = s; std::transform(r.begin(), r.end(), r.begin(), ::tolower); return r; }

bool ends_with(const std::string& s, const std::string& suffix) {
    return s.size() >= suffix.size() &&
           s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string probeOutputPrefix(const std::string& caseName) {
    return caseName + (ends_with(caseName, ".fdtd") ? "_" : ".fdtd_");
}

std::string trim_fortran_field(std::string value) {
    const size_t first = value.find_first_not_of(' ');
    return first == std::string::npos ? std::string() : value.substr(first);
}

std::string formatFortranE(double value, int width, int precision) {
    const bool negative = std::signbit(value);
    const double abs_value = std::abs(value);
    int exponent = 0;
    std::string mantissa(precision + 2, '0');
    mantissa[1] = '.';

    if (abs_value > 0.0) {
        char buffer[128];
        std::snprintf(buffer, sizeof(buffer), "%.*E", std::max(0, precision - 1), abs_value);
        const std::string normalized(buffer);
        const size_t exp_pos = normalized.find('E');
        const std::string significand = normalized.substr(0, exp_pos);
        exponent = std::stoi(normalized.substr(exp_pos + 1)) + 1;

        std::string digits;
        digits.reserve(static_cast<size_t>(precision));
        for (char c : significand) {
            if (c != '.') digits.push_back(c);
        }
        if (digits.size() < static_cast<size_t>(precision)) {
            digits.append(static_cast<size_t>(precision) - digits.size(), '0');
        }
        mantissa = "0." + digits.substr(0, static_cast<size_t>(precision));
    }

    std::ostringstream out;
    if (negative) out << '-';
    out << mantissa;
    out << 'E' << (exponent >= 0 ? '+' : '-')
        << std::setw(3) << std::setfill('0') << std::abs(exponent);

    std::string formatted = out.str();
    if (static_cast<int>(formatted.size()) < width) {
        formatted.insert(formatted.begin(), width - static_cast<int>(formatted.size()), ' ');
    }
    return formatted;
}

std::string formatFortranNegativeZero(int width, int precision) {
    std::string formatted = "-0." + std::string(static_cast<size_t>(precision), '0') + "E+000";
    if (static_cast<int>(formatted.size()) < width) {
        formatted.insert(formatted.begin(), width - static_cast<int>(formatted.size()), ' ');
    }
    return formatted;
}

struct ExcitationData { std::vector<double> times; std::vector<fdtd_real> values; };
ExcitationData readExcitationFile(const std::string& fn) {
    ExcitationData exc;
    std::ifstream f(fn);
    if (!f.is_open()) return exc;
    double t, v;
    while (f >> t >> v) { exc.times.push_back(t); exc.values.push_back(v); }
    return exc;
}
double getExcitationValue(const ExcitationData& exc, double t) {
    if (exc.times.empty()) return 0.0;
    if (t <= exc.times.front()) return exc.values.front();
    if (t >= exc.times.back()) return exc.values.back();
    for (size_t i = 1; i < exc.times.size(); i++) {
        if (t >= exc.times[i-1] && t <= exc.times[i]) {
            double frac = (t - exc.times[i-1]) / (exc.times[i] - exc.times[i-1]);
            return exc.values[i-1] + frac * (exc.values[i] - exc.values[i-1]);
        }
    }
    return exc.values.back();
}

Parseador_t parseFDTDJSON(const std::string& filename) {
    Parseador_t pd;
    std::ifstream f(filename);
    if (!f.is_open()) { std::cerr << "ERROR: Cannot open: " << filename << std::endl; return pd; }
    nlohmann::json root;
    f >> root;
    f.close();

    if (root.contains("general")) {
        auto& g = root["general"];
        if (g.contains("timeStep")) pd.general.dt = g["timeStep"].get<double>();
        if (g.contains("numberOfSteps")) pd.general.XE = g["numberOfSteps"].get<int>();
        if (g.contains("additionalArguments")) {
            pd.general.additionalArguments = g["additionalArguments"].get<std::string>();
        }
        if (g.contains("grid")) {
            auto& grid = g["grid"];
            if (grid.contains("numberOfCells")) {
                auto& nc = grid["numberOfCells"];
                pd.general.XI = nc[0].get<int>();
                pd.general.YI = nc[1].get<int>();
                pd.general.ZI = nc[2].get<int>();
            }
        }
    }
    if (root.contains("mesh") && root["mesh"].contains("grid")) {
        auto& grid = root["mesh"]["grid"];
        if (grid.contains("numberOfCells")) {
            auto& nc = grid["numberOfCells"];
            pd.general.XI = nc[0].get<int>();
            pd.general.YI = nc[1].get<int>();
            pd.general.ZI = nc[2].get<int>();
        }
        if (grid.contains("steps")) {
            auto& steps = grid["steps"];
            if (steps.contains("x")) for (auto& s : steps["x"]) pd.cellSteps.cellStepsX.push_back(s.get<double>());
            if (steps.contains("y")) for (auto& s : steps["y"]) pd.cellSteps.cellStepsY.push_back(s.get<double>());
            if (steps.contains("z")) for (auto& s : steps["z"]) pd.cellSteps.cellStepsZ.push_back(s.get<double>());
        }
    }
    if (root.contains("materials")) {
        for (auto& mat : root["materials"]) {
            Material_t m;
            m.id = mat.value("id", 0);
            m.name = mat.value("name", "");
            m.type = mat.value("type", "vacuum");
            if (m.type == "isotropic") {
                m.relativePermittivity = mat.value("relativePermittivity", 1.0);
                m.relativePermeability = mat.value("relativePermeability", 1.0);
                m.electricConductivity = mat.value("electricConductivity", 0.0);
                m.magneticConductivity = mat.value("magneticConductivity", 0.0);
            } else if (m.type == "wire") { m.radius = mat.value("radius", 0.01); }
            pd.materials.materials.push_back(m);
        }
    }
    pd.Mats.nMaterials = pd.materials.materials.size();
    if (root.contains("materialAssociations")) {
        for (auto& ma : root["materialAssociations"]) {
            pd.matAssoc.materialIds.push_back(ma["materialId"].get<int>());
            std::vector<int> eids;
            for (auto& e : ma["elementIds"]) eids.push_back(e.get<int>());
            pd.matAssoc.associations.push_back(eids);
        }
    }
    if (root.contains("boundary")) {
        auto& bd = root["boundary"];
        boundary_t b;
        if (bd.contains("all")) b.type = bd["all"]["type"].get<std::string>();
        else if (bd.contains("xLower")) b.type = bd["xLower"]["type"].get<std::string>();
        pd.boundaries.boundaries.push_back(b);
    }
    if (root.contains("sources")) {
        for (auto& src : root["sources"]) {
            source_t s;
            s.type = src["type"].get<std::string>();
            if (src.contains("name")) s.name = src["name"].get<std::string>();
            if (src.contains("magnitudeFile")) s.magnitudeFile = src["magnitudeFile"].get<std::string>();
            if (src.contains("field")) s.field = src["field"].get<std::string>();
            if (src.contains("elementIds")) for (auto& e : src["elementIds"]) s.elementIds.push_back(e.get<int>());
            if (src.contains("direction")) {
                s.direction.theta = src["direction"].value("theta", 0.0);
                s.direction.phi = src["direction"].value("phi", 0.0);
            }
            if (src.contains("polarization")) {
                s.polarization.theta = src["polarization"].value("theta", 1.5708);
                s.polarization.phi = src["polarization"].value("phi", 0.0);
            }
            if (s.type == "planewave") {
                pd.sources.planeWaves.push_back(s);
            } else if (s.type == "nodalSource") {
                pd.sources.nodalSources.push_back(s);
            }
        }
    }
    std::map<int, std::array<int, 3>> coordPos;
    std::map<int, std::vector<int>> elementCoordIds;
    if (root.contains("mesh")) {
        if (root["mesh"].contains("coordinates")) {
            for (const auto& c : root["mesh"]["coordinates"]) {
                const int id = c.value("id", 0);
                const auto& rp = c["relativePosition"];
                coordPos[id] = {rp[0].get<int>(), rp[1].get<int>(), rp[2].get<int>()};
            }
        }
        if (root["mesh"].contains("elements")) {
            for (const auto& e : root["mesh"]["elements"]) {
                const int id = e.value("id", 0);
                if (e.contains("coordinateIds")) {
                    for (const auto& cid : e["coordinateIds"]) {
                        elementCoordIds[id].push_back(cid.get<int>());
                    }
                }
            }
        }
    }
    if (root.contains("probes")) {
        int pid = 0;
        for (auto& pr : root["probes"]) {
            probe_output_t p;
            p.name = pr.value("name", std::string("probe_") + std::to_string(pid));
            p.type = pr.value("type", std::string("point"));
            p.field = pr.value("field", std::string("electric"));
            if (pr.contains("component")) p.component = pr["component"].get<std::string>();
            if (pr.contains("elementIds")) for (auto& e : pr["elementIds"]) p.elementIds.push_back(e.get<int>());
            if (pr.contains("directions")) for (auto& d : pr["directions"]) p.directions.push_back(d.get<std::string>());
            p.domainType = pr.value("domain", nlohmann::json::object()).value("type", std::string("time"));
            if (!p.elementIds.empty()) {
                const int elem_id = p.elementIds[0];
                if (elementCoordIds.count(elem_id) && !elementCoordIds[elem_id].empty()) {
                    const int coord_id = elementCoordIds[elem_id][0];
                    if (coordPos.count(coord_id)) {
                        p.coordinateId = coord_id;
                        p.cellI = coordPos[coord_id][0];
                        p.cellJ = coordPos[coord_id][1];
                        p.cellK = coordPos[coord_id][2];
                    }
                } else if (coordPos.count(elem_id)) {
                    p.coordinateId = elem_id;
                    p.cellI = coordPos[elem_id][0];
                    p.cellJ = coordPos[elem_id][1];
                    p.cellK = coordPos[elem_id][2];
                }
            }
            p.fieldByDir.resize(p.directions.size());
            p.incidentByDir.resize(p.directions.size());
            p.outputCellI = p.cellI;
            p.outputCellJ = p.cellJ;
            p.outputCellK = p.cellK;
            p.probeId = pid++;
            pd.probes.probes.push_back(p);
        }
    }
    if (root.contains("mesh") && root["mesh"].contains("elements")) {
        for (auto& e : root["mesh"]["elements"]) {
            if (e.contains("intervals")) {
                for (auto& interval : e["intervals"]) {
                    std::vector<double> iv;
                    for (auto& coord : interval) for (auto& v : coord) iv.push_back(v.get<double>());
                    pd.elements.intervals.push_back(iv);
                }
            }
        }
    }
    bool has_json_dt = root.contains("general") && root["general"].contains("timeStep");
    // Compute CFL-limited dt only when JSON does not specify timeStep (Fortran behavior).
    if (!has_json_dt && root.contains("mesh") && root["mesh"].contains("grid") &&
        root["mesh"]["grid"].contains("steps")) {
        auto& steps = root["mesh"]["grid"]["steps"];
        double dx_min = 1.0, dy_min = 1.0, dz_min = 1.0;
        if (steps.contains("x")) for (auto& s : steps["x"]) dx_min = std::min(dx_min, s.get<double>());
        if (steps.contains("y")) for (auto& s : steps["y"]) dy_min = std::min(dy_min, s.get<double>());
        if (steps.contains("z")) for (auto& s : steps["z"]) dz_min = std::min(dz_min, s.get<double>());
        pd.general.dt = 0.99 / (C0 * std::sqrt(1.0/(dx_min*dx_min) + 1.0/(dy_min*dy_min) + 1.0/(dz_min*dz_min)));
    }
    pd.matriz.totalX = pd.general.XI;
    pd.matriz.totalY = pd.general.YI;
    pd.matriz.totalZ = pd.general.ZI;
    pd.despl.XI = pd.general.XI; pd.despl.YI = pd.general.YI; pd.despl.ZI = pd.general.ZI;
    return pd;
}

class FDTD_Solver {
public:
    struct ProbeCellBounds {
        int xi = 1, yi = 1, zi = 1;
        int xe = 1, ye = 1, ze = 1;
        bool valid = false;
    };

    Parseador_t pd;
    int NX = 10, NY = 10, NZ = 10;
    double dt = 1e-12, dx = 0.01, dy = 0.01, dz = 0.01;
    double wireDx = 0.01, wireDy = 0.01, wireDz = 0.01;
    double eps0 = EPS0, mu0 = MU0;
    std::vector<fdtd_real> Ex, Ey, Ez, Hx, Hy, Hz;
    std::vector<uint8_t> pecExMask;
    std::vector<uint8_t> pecEyMask;
    std::vector<uint8_t> pecEzMask;
    bool hasAnyPecMask = false;
    std::vector<fdtd_real> CeEx, CeEy, CeEz, CmH;
    std::vector<fdtd_real> Idxe, Idye, Idze, Idxh, Idyh, Idzh;
    std::vector<source_t> sources;
    std::map<std::string, ExcitationData> excitations;
    std::vector<PlaneWaveState_t> planeWaves;
    std::vector<probe_output_t> probes;
    std::vector<BulkCurrentProbe_t> bulkCurrentProbes;
    std::map<std::string, std::vector<double>> analyticBulkCurrents;
    std::vector<NodalCurrentSegment_t> nodalCurrentSegments;
    std::vector<SgbcNode_t> sgbcNodes;
    std::vector<HollandWireSegment_t> hollandSegments;
    std::vector<HollandWireNode_t> hollandNodes;
    std::vector<HollandWireProbe_t> hollandProbes;
    std::vector<HollandVoltageGenerator_t> hollandVoltageGenerators;
    std::map<int, std::pair<bool, double>> hollandNodeTermination;
    struct MovieProbeState {
        std::string stem;
        ProbeCellBounds bounds;
        enum class FieldMode { ElectricMagnitude, Ex, Ey, Ez, Hx, Hy, Hz } mode =
            FieldMode::ElectricMagnitude;
        int nx = 0, ny = 0, nz = 0;
        double initialTime = 0.0;
        double finalTime = 0.0;
        double samplingPeriod = 0.0;
        int trancos = 1;
        std::vector<float> samples;
        std::vector<double> times;
    };
    std::vector<MovieProbeState> movieProbes;
    Lumped_m::LumpedSolver_t lumpedSolver;
#ifdef CompileWithMTLN
    mtln_solver_m::mtln_t mtlnSolver;
    bool hasMtlnSolver = false;
    bool mtlnObservationOpen = false;
    std::vector<double> mtlnExternalFields;
#endif
    bool still_planewave_time = true;
    bool planewave_switched_off = false;
    bool useMur = true;
    bool usePml = false;
    bool usePec = false;
    bool periodicBack = false, periodicFront = false;
    bool periodicLeft = false, periodicRight = false;
    bool periodicDown = false, periodicUp = false;
    bool murBack = true, murFront = true;
    bool murLeft = true, murRight = true;
    bool murDown = true, murUp = true;
    bool pmlBack = false, pmlFront = false;
    bool pmlLeft = false, pmlRight = false;
    bool pmlDown = false, pmlUp = false;
    struct PmlFaceConfig {
        bool enabled = false;
        int layers = 0;
        double order = 2.0;
        double reflection = 0.001;
    };
    PmlFaceConfig pmlFaceBack, pmlFaceFront;
    PmlFaceConfig pmlFaceLeft, pmlFaceRight;
    PmlFaceConfig pmlFaceDown, pmlFaceUp;
    bool cpmlBordersInitialized = false;
    std::vector<fdtd_real> pmlPceX, pmlPceY, pmlPceZ;
    std::vector<fdtd_real> pmlPbeX, pmlPbeY, pmlPbeZ;
    std::vector<fdtd_real> pmlPcmX, pmlPcmY, pmlPcmZ;
    std::vector<fdtd_real> pmlPbmX, pmlPbmY, pmlPbmZ;
    std::vector<fdtd_real> psiExy, psiExz, psiEyz, psiEyx, psiEzx, psiEzy;
    std::vector<fdtd_real> psiHxy, psiHxz, psiHyz, psiHyx, psiHzx, psiHzy;
    bool pmcBack = false, pmcFront = false;
    bool pmcLeft = false, pmcRight = false;
    bool pmcDown = false, pmcUp = false;
    int pmlElectricCalls = 0;
    int pmlBodyHCalls = 0;
    int pmlMagneticCpmlCalls = 0;
    // MURc zones (InitMURBorders bordersmur.F90 L94-137) — thin 1-cell face pads per component.
    struct MurZone { int xi = 0, xe = 0, yi = 0, ye = 0, zi = 0, ze = 0; };
    MurZone murHyBack_, murHyFront_, murHzBack_, murHzFront_;
    MurZone murHxLeft_, murHxRight_, murHzLeft_, murHzRight_;
    MurZone murHyDown_, murHyUp_, murHxDown_, murHxUp_;
    fdtd_real backCab1 = 0.0, frontCab1 = 0.0, leftCab1 = 0.0, rightCab1 = 0.0;
    fdtd_real downCab1 = 0.0, upCab1 = 0.0;
    // First-order magnetic Mur past fields (Fortran regLR/regDU/regBF Past_*).
    std::vector<fdtd_real> murPastHyBack, murPastHyBackInt, murPastHzBack, murPastHzBackInt;
    std::vector<fdtd_real> murPastHyFront, murPastHyFrontInt, murPastHzFront, murPastHzFrontInt;
    std::vector<fdtd_real> murPastHxLeft, murPastHxLeftInt, murPastHzLeft, murPastHzLeftInt;
    std::vector<fdtd_real> murPastHxRight, murPastHxRightInt, murPastHzRight, murPastHzRightInt;
    std::vector<fdtd_real> murPastHyDown, murPastHyDownInt, murPastHxDown, murPastHxDownInt;
    std::vector<fdtd_real> murPastHyUp, murPastHyUpInt, murPastHxUp, murPastHxUpInt;
    double murCx = 0.0, murCy = 0.0, murCz = 0.0;
    int numSteps = 100, step = 0;
    double currentTime = 0.0;
    bool createMapVtk = false;
    nlohmann::json inputRoot;
    std::string inputFile;
    bool mtlnPmlPaddingActive = false;
    int pmlPadX = 0, pmlPadY = 0, pmlPadZ = 0;
    int fieldHalo = 2;
    int mpiLayoutNumber = 0;
    int mpiNumProcs = 1;
    int mpiAxis = 3;
    bool mpiEnabled = false;
    std::vector<MpiSliceInfo> mpiSlices;

    void init(const std::string& filename, bool map_vtk = false,
              int mpi_rank = 0, int mpi_size = 1, int mpi_axis = 3) {
        preserveFortranSubnormalArithmetic();
        inputFile = filename;
        mpiLayoutNumber = mpi_rank;
        mpiNumProcs = std::max(1, mpi_size);
        mpiAxis = (mpi_axis >= 1 && mpi_axis <= 3) ? mpi_axis : 3;
        mpiEnabled = mpiNumProcs > 1;
        createMapVtk = map_vtk;
        std::ifstream jf(filename);
        if (jf.is_open()) {
            jf >> inputRoot;
            jf.close();
        }
        pd = parseFDTDJSON(filename);
        applyMtlnPmlPaddingIfNeeded();
        NX = pd.general.XI; NY = pd.general.YI; NZ = pd.general.ZI;
        if (NX <= 0) NX = 10; if (NY <= 0) NY = 10; if (NZ <= 0) NZ = 10;
        dt = static_cast<double>(static_cast<fdtd_real>(pd.general.dt));
        if (dt <= 0.0) dt = 1e-12;
        if (!pd.cellSteps.cellStepsX.empty()) {
            wireDx = pd.cellSteps.cellStepsX[0];
            dx = static_cast<double>(static_cast<fdtd_real>(pd.cellSteps.cellStepsX[0]));
        } else {
            dx = 0.025;
            wireDx = dx;
        }
        if (!pd.cellSteps.cellStepsY.empty()) {
            wireDy = pd.cellSteps.cellStepsY[0];
            dy = static_cast<double>(static_cast<fdtd_real>(pd.cellSteps.cellStepsY[0]));
        } else {
            dy = 0.025;
            wireDy = dy;
        }
        if (!pd.cellSteps.cellStepsZ.empty()) {
            wireDz = pd.cellSteps.cellStepsZ[0];
            dz = static_cast<double>(static_cast<fdtd_real>(pd.cellSteps.cellStepsZ[0]));
        } else {
            dz = 0.025;
            wireDz = dz;
        }
        numSteps = pd.general.XE; if (numSteps <= 0) numSteps = 100;

        const double dt_before = dt;
        const int steps_before = numSteps;
        const fdtd_real eps0_r = static_cast<fdtd_real>(EPS0);
        const fdtd_real mu0_r = static_cast<fdtd_real>(MU0);
        const fdtd_real cluz_r = static_cast<fdtd_real>(1.0) / std::sqrt(eps0_r * mu0_r);
        const fdtd_real inv_dx = static_cast<fdtd_real>(1.0) / static_cast<fdtd_real>(dx);
        const fdtd_real inv_dy = static_cast<fdtd_real>(1.0) / static_cast<fdtd_real>(dy);
        const fdtd_real inv_dz = static_cast<fdtd_real>(1.0) / static_cast<fdtd_real>(dz);
        const fdtd_real dtlay = static_cast<fdtd_real>(1.0) /
            (cluz_r * std::sqrt(inv_dx * inv_dx + inv_dy * inv_dy + inv_dz * inv_dz));
        if (dt_before <= 0.0 || dt_before > static_cast<double>(dtlay * static_cast<fdtd_real>(heurCFL))) {
            dt = static_cast<double>(dtlay * static_cast<fdtd_real>(heurCFL));
            if (dt_before > 0.0 && steps_before > 0) {
                numSteps = static_cast<int>(dt_before / dt * steps_before);
            }
            std::cout << "CFL correction: dt " << dt_before << " -> " << dt
                      << ", steps " << steps_before << " -> " << numSteps << std::endl;
        }

        initBoundaryFlagsFromJson();
        initMpiOneAxisDecomposition();

        int ex_n = ex_nx()*ex_ny()*ex_nz();
        int ey_n = ey_nx()*ey_ny()*ey_nz();
        int ez_n = ez_nx()*ez_ny()*ez_nz();
        int hx_n = hx_nx()*hx_ny()*hx_nz();
        int hy_n = hy_nx()*hy_ny()*hy_nz();
        int hz_n = hz_nx()*hz_ny()*hz_nz();
        Ex.resize(ex_n,0); Ey.resize(ey_n,0); Ez.resize(ez_n,0);
        Hx.resize(hx_n,0); Hy.resize(hy_n,0); Hz.resize(hz_n,0);
        int max_n = std::max({ex_n,ey_n,ez_n,hx_n,hy_n,hz_n});
        CeEx.resize(ex_n, 0.0);
        CeEy.resize(ey_n, 0.0);
        CeEz.resize(ez_n, 0.0);
        CmH.resize(max_n, 0.0);
        pecExMask.assign(static_cast<size_t>(ex_n), 0);
        pecEyMask.assign(static_cast<size_t>(ey_n), 0);
        pecEzMask.assign(static_cast<size_t>(ez_n), 0);
        hasAnyPecMask = false;

        const fdtd_real ce = static_cast<fdtd_real>(
            dt / static_cast<double>(static_cast<fdtd_real>(eps0)));
        const fdtd_real ch = static_cast<fdtd_real>(
            dt / static_cast<double>(static_cast<fdtd_real>(mu0)));
        std::fill(CeEx.begin(), CeEx.end(), ce);
        std::fill(CeEy.begin(), CeEy.end(), ce);
        std::fill(CeEz.begin(), CeEz.end(), ce);
        std::fill(CmH.begin(), CmH.end(), ch);

        initIsotropicMaterialCoefficientsFromJson();

        sources = pd.sources.planeWaves;
        sources.insert(sources.end(), pd.sources.nodalSources.begin(),
                       pd.sources.nodalSources.end());
        const std::filesystem::path json_dir =
            std::filesystem::path(filename).parent_path();
        if (inputRoot.contains("sources")) {
            for (const auto& src : inputRoot["sources"]) {
                const std::string magnitudeFile = src.value("magnitudeFile", std::string());
                if (magnitudeFile.empty() || excitations.count(magnitudeFile)) continue;
                std::string exc_path = magnitudeFile;
                if (!std::filesystem::exists(exc_path) && !json_dir.empty()) {
                    exc_path = (json_dir / magnitudeFile).string();
                }
                excitations[magnitudeFile] = readExcitationFile(exc_path);
            }
        }
        for (auto& src : sources) {
            if (src.magnitudeFile.empty()) continue;
            if (excitations.count(src.magnitudeFile)) continue;
            std::string exc_path = src.magnitudeFile;
            if (!std::filesystem::exists(exc_path) && !json_dir.empty()) {
                exc_path = (json_dir / src.magnitudeFile).string();
            }
            excitations[src.magnitudeFile] = readExcitationFile(exc_path);
        }
        initGridInverses();
        initCpmlBorders();
        planeWaves.resize(pd.sources.planeWaves.size());
        for (int i = 0; i < (int)pd.sources.planeWaves.size(); i++) {
            planeWaves[i].px.resize(1,0); planeWaves[i].py.resize(1,0); planeWaves[i].pz.resize(1,0);
            planeWaves[i].ex.resize(1,0); planeWaves[i].ey.resize(1,0); planeWaves[i].ez.resize(1,0);
            planeWaves[i].hx.resize(1,0); planeWaves[i].hy.resize(1,0); planeWaves[i].hz.resize(1,0);
            initPlaneWave(i);
        }
        probes = pd.probes.probes;
        initInternalPecFromJson();
        initBoundaryPecMasksFromJson();
        initSgbcFromJson();
        initBulkCurrentProbes();
        initAnalyticLumpedCurrentsFromJson();
        initAnalyticSurfaceImpedanceCurrentsFromJson();
        initAnalyticConformalCylinderCurrentsFromJson();
        initNodalCurrentSources();
#ifdef CompileWithMTLN
        initMtlnFromJson(filename);
#endif
        if (!mtlnOwnsWires()) {
            initHollandWires();
        }
        initLumpedFromJson();
        initMurBorders();
        initMovieProbesFromJson(SEMBA_FDTD_m::extractCaseNameFromInput(filename));

        std::cout << "FDTD: grid=" << NX << "x" << NY << "x" << NZ << " dt=" << dt << " steps=" << numSteps << std::endl;
    }

    static void shiftJsonCoordinate(nlohmann::json& coord, int dx, int dy, int dz) {
        if (!coord.is_array() || coord.size() < 3) return;
        const int delta[3] = {dx, dy, dz};
        for (int axis = 0; axis < 3; ++axis) {
            if (coord[axis].is_number_integer()) {
                coord[axis] = coord[axis].get<int>() + delta[axis];
            } else if (coord[axis].is_number()) {
                coord[axis] = coord[axis].get<double>() + static_cast<double>(delta[axis]);
            }
        }
    }

    bool hasMtlnCableInJson() const {
        if (inputRoot.is_null()) return false;
        if (inputRoot.contains("materialAssociations")) {
            for (const auto& assoc : inputRoot["materialAssociations"]) {
                if (to_lower(assoc.value("type", std::string())) == "cable") {
                    return true;
                }
            }
        }
        if (inputRoot.contains("materials")) {
            for (const auto& mat : inputRoot["materials"]) {
                const std::string type = to_lower(mat.value("type", std::string()));
                if (type == "unshieldedmultiwire" || type == "shieldedmultiwire") {
                    return true;
                }
            }
        }
        return false;
    }

    int pmlPaddingLayersFromJson() const {
        if (inputRoot.is_null() || !inputRoot.contains("boundary")) return 0;
        const auto& boundary = inputRoot["boundary"];
        if (!boundary.contains("all")) return 0;
        const auto& all = boundary["all"];
        if (to_lower(all.value("type", std::string())) != "pml") return 0;
        return std::max(0, static_cast<int>(std::lround(all.value("layers", 0.0))));
    }

    void padInputRootForMtlnPml() {
        if (!inputRoot.contains("mesh")) return;
        auto& mesh = inputRoot["mesh"];
        if (mesh.contains("grid") && mesh["grid"].contains("numberOfCells")) {
            auto& n = mesh["grid"]["numberOfCells"];
            if (n.is_array() && n.size() >= 3) {
                n[0] = n[0].get<int>() + 2 * pmlPadX;
                n[1] = n[1].get<int>() + 2 * pmlPadY;
                n[2] = n[2].get<int>() + 2 * pmlPadZ;
            }
        }
        if (mesh.contains("coordinates")) {
            for (auto& coord : mesh["coordinates"]) {
                if (coord.contains("relativePosition")) {
                    shiftJsonCoordinate(coord["relativePosition"], pmlPadX, pmlPadY, pmlPadZ);
                }
            }
        }
        if (mesh.contains("elements")) {
            for (auto& elem : mesh["elements"]) {
                if (!elem.contains("intervals")) continue;
                for (auto& interval : elem["intervals"]) {
                    if (!interval.is_array() || interval.size() < 2) continue;
                    shiftJsonCoordinate(interval[0], pmlPadX, pmlPadY, pmlPadZ);
                    shiftJsonCoordinate(interval[1], pmlPadX, pmlPadY, pmlPadZ);
                }
            }
        }
    }

    void padParsedDataForMtlnPml() {
        pd.general.XI += 2 * pmlPadX;
        pd.general.YI += 2 * pmlPadY;
        pd.general.ZI += 2 * pmlPadZ;
        pd.matriz.totalX = pd.general.XI;
        pd.matriz.totalY = pd.general.YI;
        pd.matriz.totalZ = pd.general.ZI;
        pd.despl.XI = pd.general.XI;
        pd.despl.YI = pd.general.YI;
        pd.despl.ZI = pd.general.ZI;
        for (auto& probe : pd.probes.probes) {
            probe.cellI += pmlPadX;
            probe.cellJ += pmlPadY;
            probe.cellK += pmlPadZ;
        }
        for (auto& interval : pd.elements.intervals) {
            if (interval.size() >= 6) {
                interval[0] += pmlPadX;
                interval[1] += pmlPadY;
                interval[2] += pmlPadZ;
                interval[3] += pmlPadX;
                interval[4] += pmlPadY;
                interval[5] += pmlPadZ;
            }
        }
    }

    void applyMtlnPmlPaddingIfNeeded() {
        mtlnPmlPaddingActive = false;
        pmlPadX = pmlPadY = pmlPadZ = 0;
#ifndef CompileWithMTLN
        return;
#else
        const int layers = pmlPaddingLayersFromJson();
        if (layers <= 0 || !hasMtlnCableInJson()) return;

        mtlnPmlPaddingActive = true;
        pmlPadX = pmlPadY = pmlPadZ = layers;
        padInputRootForMtlnPml();
        padParsedDataForMtlnPml();
#endif
    }

    bool mtlnOwnsWires() const {
#ifdef CompileWithMTLN
        return hasMtlnSolver;
#else
        return false;
#endif
    }

#ifdef CompileWithMTLN
    std::string caseNameStem() const {
        std::string name = std::filesystem::path(inputFile).filename().string();
        const std::string jsonSuffix = ".json";
        if (name.size() > jsonSuffix.size() &&
            name.compare(name.size() - jsonSuffix.size(), jsonSuffix.size(), jsonSuffix) == 0) {
            name.resize(name.size() - jsonSuffix.size());
        }
        return name;
    }

    void initMtlnFromJson(const std::string& filename) {
        hasMtlnSolver = false;
        mtlnObservationOpen = false;
        mtlnExternalFields.clear();

        smbjson::parser_t parser(filename);
        NFDETypes_m::Parseador_t parsed = parser.readProblemDescription();
        if (!parsed.mtln || parsed.mtln->cables.empty()) {
            return;
        }

        parsed.mtln->time_step = dt;
        parsed.mtln->number_of_steps = numSteps;
        mtlnSolver = mtln_solver_m::mtlnCtor(*parsed.mtln);
        hasMtlnSolver = mtlnSolver.number_of_bundles > 0;
        if (!hasMtlnSolver) {
            return;
        }
        shiftMtlnExternalFieldSegmentsForPmlPadding();
        mtlnSolver.updatePULTerms();
        setupMtlnExternalFieldPointers();
        deembedPecMasksForMtlnSegments();
    }

    void shiftMtlnExternalFieldSegmentsForPmlPadding() {
        if (!mtlnPmlPaddingActive) return;
        for (auto& bundle : mtlnSolver.bundles) {
            for (auto& segment : bundle.external_field_segments) {
                if (segment.position.size() < 3) continue;
                segment.position[0] += pmlPadX;
                segment.position[1] += pmlPadY;
                segment.position[2] += pmlPadZ;
            }
        }
    }

    void setupMtlnExternalFieldPointers() {
        size_t count = 0;
        for (const auto& bundle : mtlnSolver.bundles) {
            count += bundle.external_field_segments.size();
        }
        mtlnExternalFields.assign(count, 0.0);

        size_t idx = 0;
        for (auto& bundle : mtlnSolver.bundles) {
            for (auto& segment : bundle.external_field_segments) {
                segment.field = &mtlnExternalFields[idx++];
            }
        }
    }

    void deembedPecMaskForMtlnSegment(
        const mtl_bundle_m::external_field_segment_t& segment) {
        if (segment.position.size() < 3) return;
        const int i = segment.position[0] - 1;
        const int j = segment.position[1] - 1;
        const int k = segment.position[2] - 1;
        switch (std::abs(segment.direction)) {
            case mtln_types_m::DIRECTION_X_POS:
                if (in_ex(i, j, k)) {
                    pecExMask[static_cast<size_t>(ex_idx(i, j, k))] = 0;
                }
                break;
            case mtln_types_m::DIRECTION_Y_POS:
                if (in_ey(i, j, k)) {
                    pecEyMask[static_cast<size_t>(ey_idx(i, j, k))] = 0;
                }
                break;
            case mtln_types_m::DIRECTION_Z_POS:
                if (in_ez(i, j, k)) {
                    pecEzMask[static_cast<size_t>(ez_idx(i, j, k))] = 0;
                }
                break;
            default:
                break;
        }
    }

    void deembedPecMasksForMtlnSegments() {
        for (const auto& bundle : mtlnSolver.bundles) {
            if (!bundle.bundle_in_layer) continue;
            for (const auto& segment : bundle.external_field_segments) {
                deembedPecMaskForMtlnSegment(segment);
            }
        }
    }

    double electricFieldForMtlnSegment(const mtl_bundle_m::external_field_segment_t& segment) const {
        if (segment.position.size() < 3) return 0.0;
        const int i = segment.position[0] - 1;
        const int j = segment.position[1] - 1;
        const int k = segment.position[2] - 1;
        switch (std::abs(segment.direction)) {
            case mtln_types_m::DIRECTION_X_POS:
                return in_ex(i, j, k) ? static_cast<double>(Ex[ex_idx(i, j, k)]) : 0.0;
            case mtln_types_m::DIRECTION_Y_POS:
                return in_ey(i, j, k) ? static_cast<double>(Ey[ey_idx(i, j, k)]) : 0.0;
            case mtln_types_m::DIRECTION_Z_POS:
                return in_ez(i, j, k) ? static_cast<double>(Ez[ez_idx(i, j, k)]) : 0.0;
            default:
                return 0.0;
        }
    }

    void syncMtlnExternalFields() {
        size_t idx = 0;
        for (const auto& bundle : mtlnSolver.bundles) {
            for (const auto& segment : bundle.external_field_segments) {
                if (idx < mtlnExternalFields.size()) {
                    mtlnExternalFields[idx] = electricFieldForMtlnSegment(segment);
                }
                ++idx;
            }
        }
    }

    double orientedMtlnCurrent(const mtl_bundle_m::mtl_bundle_t& bundle, size_t segmentIdx) const {
        const auto& segment = bundle.external_field_segments[segmentIdx];
        const int numConductors = bundle.conductors_in_level.empty()
                                      ? bundle.number_of_conductors
                                      : bundle.conductors_in_level[0];
        double current = 0.0;
        for (int c = 0; c < numConductors; ++c) {
            if (c >= static_cast<int>(bundle.i.size())) break;
            const auto& conductorCurrent = bundle.i[static_cast<size_t>(c)];
            if (segmentIdx >= conductorCurrent.size()) continue;
            current += conductorCurrent[segmentIdx];
        }
        return current * std::copysign(1.0, static_cast<double>(segment.direction));
    }

    double mtlnFieldSubtractForSegment(const mtl_bundle_m::mtl_bundle_t& bundle,
                                       size_t segmentIdx) const {
        const auto& segment = bundle.external_field_segments[segmentIdx];
        if (segment.position.size() < 3) return 0.0;
        const int i = segment.position[0] - 1;
        const int j = segment.position[1] - 1;
        const int k = segment.position[2] - 1;
        double dSInverse = 0.0;
        switch (std::abs(segment.direction)) {
            case mtln_types_m::DIRECTION_X_POS:
                dSInverse = static_cast<double>(idyh1(j)) * static_cast<double>(idzh1(k));
                break;
            case mtln_types_m::DIRECTION_Y_POS:
                dSInverse = static_cast<double>(idxh1(i)) * static_cast<double>(idzh1(k));
                break;
            case mtln_types_m::DIRECTION_Z_POS:
                dSInverse = static_cast<double>(idxh1(i)) * static_cast<double>(idyh1(j));
                break;
            default:
                return 0.0;
        }
        return (dt / eps0) * dSInverse * orientedMtlnCurrent(bundle, segmentIdx);
    }

    void subtractMtlnCurrentsFromFields() {
        for (const auto& bundle : mtlnSolver.bundles) {
            if (!bundle.bundle_in_layer) continue;
            for (size_t segmentIdx = 0; segmentIdx < bundle.external_field_segments.size(); ++segmentIdx) {
                const auto& segment = bundle.external_field_segments[segmentIdx];
                if (segment.position.size() < 3) continue;
                const double subtractValue = mtlnFieldSubtractForSegment(bundle, segmentIdx);
                const int i = segment.position[0] - 1;
                const int j = segment.position[1] - 1;
                const int k = segment.position[2] - 1;
                switch (std::abs(segment.direction)) {
                    case mtln_types_m::DIRECTION_X_POS:
                        if (in_ex(i, j, k) && !isPecEx(i, j, k)) {
                            const int idx = ex_idx(i, j, k);
                            Ex[static_cast<size_t>(idx)] = static_cast<fdtd_real>(
                                static_cast<double>(Ex[static_cast<size_t>(idx)]) - subtractValue);
                        }
                        break;
                    case mtln_types_m::DIRECTION_Y_POS:
                        if (in_ey(i, j, k) && !isPecEy(i, j, k)) {
                            const int idx = ey_idx(i, j, k);
                            Ey[static_cast<size_t>(idx)] = static_cast<fdtd_real>(
                                static_cast<double>(Ey[static_cast<size_t>(idx)]) - subtractValue);
                        }
                        break;
                    case mtln_types_m::DIRECTION_Z_POS:
                        if (in_ez(i, j, k) && !isPecEz(i, j, k)) {
                            const int idx = ez_idx(i, j, k);
                            Ez[static_cast<size_t>(idx)] = static_cast<fdtd_real>(
                                static_cast<double>(Ez[static_cast<size_t>(idx)]) - subtractValue);
                        }
                        break;
                    default:
                        break;
                }
            }
        }
    }

    void openMtlnObservation() {
        if (!hasMtlnSolver || mtlnObservationOpen) return;
        mtlnSolver.initObservation(caseNameStem());
        mtlnObservationOpen = true;
    }

    void advanceMtlnE() {
        if (!hasMtlnSolver) return;
        subtractMtlnCurrentsFromFields();
        syncMtlnExternalFields();
        mtlnSolver.step();
        if (mtlnObservationOpen) {
            mtlnSolver.updateObservation(step);
        }
    }

    void closeMtlnObservation() {
        if (!mtlnObservationOpen) return;
        mtlnSolver.closeObservation();
        mtlnObservationOpen = false;
    }
#endif

    void calcMurConstants() {
        const fdtd_real one = static_cast<fdtd_real>(1.0);
        const fdtd_real cluz = static_cast<fdtd_real>(1.0) /
            std::sqrt(static_cast<fdtd_real>(eps0) * static_cast<fdtd_real>(mu0));
        const auto cab1 = [this, one, cluz](fdtd_real inv_step,
                                            double relativePermittivity) {
            const fdtd_real cnum = static_cast<fdtd_real>(
                static_cast<double>(one / inv_step) *
                std::sqrt(relativePermittivity) /
                (dt * static_cast<double>(cluz)));
            return (one - cnum) / (one + cnum);
        };
        backCab1 = cab1(fieldGridInverse(static_cast<fdtd_real>(dx)),
                        murFaceRelativePermittivity("xLower"));
        frontCab1 = cab1(fieldGridInverse(static_cast<fdtd_real>(dx)),
                         murFaceRelativePermittivity("xUpper"));
        leftCab1 = cab1(fieldGridInverse(static_cast<fdtd_real>(dy)),
                        murFaceRelativePermittivity("yLower"));
        rightCab1 = cab1(fieldGridInverse(static_cast<fdtd_real>(dy)),
                         murFaceRelativePermittivity("yUpper"));
        downCab1 = cab1(fieldGridInverse(static_cast<fdtd_real>(dz)),
                        murFaceRelativePermittivity("zLower"));
        upCab1 = cab1(fieldGridInverse(static_cast<fdtd_real>(dz)),
                      murFaceRelativePermittivity("zUpper"));
        murCx = backCab1;
        murCy = leftCab1;
        murCz = downCab1;
    }

    void initMurBorders() {
        if (!useMur) return;
        // Hy/Hz back & front (bordersmur.F90 L124-136, AdvanceMagneticMUR L1280-1360).
        murHyBack_  = {0, 0, 0, NY - 1, 0, NZ};
        murHyFront_ = {NX, NX, 0, NY - 1, 0, NZ};
        murHzBack_  = {0, 0, 0, NY, 0, NZ - 1};
        murHzFront_ = {NX, NX, 0, NY, 0, NZ - 1};
        // Hx/Hz left & right (L110-122).
        murHxLeft_  = {0, NX - 1, 0, 0, 0, NZ};
        murHxRight_ = {0, NX - 1, NY, NY, 0, NZ};
        murHzLeft_  = {0, NX, 0, 0, 0, NZ - 1};
        murHzRight_ = {0, NX, NY, NY, 0, NZ - 1};
        // Hy/Hx down & up (L100-108).
        murHyDown_  = {0, NX, 0, NY - 1, 0, 0};
        murHyUp_    = {0, NX, 0, NY - 1, NZ, NZ};
        murHxDown_  = {0, NX - 1, 0, NY, 0, 0};
        murHxUp_    = {0, NX - 1, 0, NY, NZ, NZ};
        calcMurConstants();
        murPastHyBack.assign((NY + 1) * NZ, 0.0);
        murPastHyBackInt.assign((NY + 1) * NZ, 0.0);
        murPastHzBack.assign(NY * (NZ + 1), 0.0);
        murPastHzBackInt.assign(NY * (NZ + 1), 0.0);
        murPastHyFront.assign((NY + 1) * NZ, 0.0);
        murPastHyFrontInt.assign((NY + 1) * NZ, 0.0);
        murPastHzFront.assign(NY * (NZ + 1), 0.0);
        murPastHzFrontInt.assign(NY * (NZ + 1), 0.0);
        murPastHxLeft.assign((NX + 1) * NZ, 0.0);
        murPastHxLeftInt.assign((NX + 1) * NZ, 0.0);
        murPastHzLeft.assign(NX * (NZ + 1), 0.0);
        murPastHzLeftInt.assign(NX * (NZ + 1), 0.0);
        murPastHxRight.assign((NX + 1) * NZ, 0.0);
        murPastHxRightInt.assign((NX + 1) * NZ, 0.0);
        murPastHzRight.assign(NX * (NZ + 1), 0.0);
        murPastHzRightInt.assign(NX * (NZ + 1), 0.0);
        murPastHyDown.assign(NX * (NY + 1), 0.0);
        murPastHyDownInt.assign(NX * (NY + 1), 0.0);
        murPastHxDown.assign((NX + 1) * NY, 0.0);
        murPastHxDownInt.assign((NX + 1) * NY, 0.0);
        murPastHyUp.assign(NX * (NY + 1), 0.0);
        murPastHyUpInt.assign(NX * (NY + 1), 0.0);
        murPastHxUp.assign((NX + 1) * NY, 0.0);
        murPastHxUpInt.assign((NX + 1) * NY, 0.0);
    }

    int ex_nx() const { return NX + 2 * fieldHalo - 1; }
    int ex_ny() const { return NY + 2 * fieldHalo; }
    int ex_nz() const { return NZ + 2 * fieldHalo; }
    int ey_nx() const { return NX + 2 * fieldHalo; }
    int ey_ny() const { return NY + 2 * fieldHalo - 1; }
    int ey_nz() const { return NZ + 2 * fieldHalo; }
    int ez_nx() const { return NX + 2 * fieldHalo; }
    int ez_ny() const { return NY + 2 * fieldHalo; }
    int ez_nz() const { return NZ + 2 * fieldHalo - 1; }
    int hx_nx() const { return NX + 2 * fieldHalo; }
    int hx_ny() const { return NY + 2 * fieldHalo - 1; }
    int hx_nz() const { return NZ + 2 * fieldHalo - 1; }
    int hy_nx() const { return NX + 2 * fieldHalo - 1; }
    int hy_ny() const { return NY + 2 * fieldHalo; }
    int hy_nz() const { return NZ + 2 * fieldHalo - 1; }
    int hz_nx() const { return NX + 2 * fieldHalo - 1; }
    int hz_ny() const { return NY + 2 * fieldHalo - 1; }
    int hz_nz() const { return NZ + 2 * fieldHalo; }

    int ex_idx(int i,int j,int k) const { return (i + fieldHalo)*ex_ny()*ex_nz() + (j + fieldHalo)*ex_nz() + (k + fieldHalo); }
    int ey_idx(int i,int j,int k) const { return (i + fieldHalo)*ey_ny()*ey_nz() + (j + fieldHalo)*ey_nz() + (k + fieldHalo); }
    int ez_idx(int i,int j,int k) const { return (i + fieldHalo)*ez_ny()*ez_nz() + (j + fieldHalo)*ez_nz() + (k + fieldHalo); }
    int hx_idx(int i,int j,int k) const { return (i + fieldHalo)*hx_ny()*hx_nz() + (j + fieldHalo)*hx_nz() + (k + fieldHalo); }
    int hy_idx(int i,int j,int k) const { return (i + fieldHalo)*hy_ny()*hy_nz() + (j + fieldHalo)*hy_nz() + (k + fieldHalo); }
    int hz_idx(int i,int j,int k) const { return (i + fieldHalo)*hz_ny()*hz_nz() + (j + fieldHalo)*hz_nz() + (k + fieldHalo); }

    bool in_ex(int i,int j,int k) const { return i >= -fieldHalo && i <= NX + fieldHalo - 2 && j >= -fieldHalo && j <= NY + fieldHalo - 1 && k >= -fieldHalo && k <= NZ + fieldHalo - 1; }
    bool in_ey(int i,int j,int k) const { return i >= -fieldHalo && i <= NX + fieldHalo - 1 && j >= -fieldHalo && j <= NY + fieldHalo - 2 && k >= -fieldHalo && k <= NZ + fieldHalo - 1; }
    bool in_ez(int i,int j,int k) const { return i >= -fieldHalo && i <= NX + fieldHalo - 1 && j >= -fieldHalo && j <= NY + fieldHalo - 1 && k >= -fieldHalo && k <= NZ + fieldHalo - 2; }
    bool in_hx(int i,int j,int k) const { return i >= -fieldHalo && i <= NX + fieldHalo - 1 && j >= -fieldHalo && j <= NY + fieldHalo - 2 && k >= -fieldHalo && k <= NZ + fieldHalo - 2; }
    bool in_hy(int i,int j,int k) const { return i >= -fieldHalo && i <= NX + fieldHalo - 2 && j >= -fieldHalo && j <= NY + fieldHalo - 1 && k >= -fieldHalo && k <= NZ + fieldHalo - 2; }
    bool in_hz(int i,int j,int k) const { return i >= -fieldHalo && i <= NX + fieldHalo - 2 && j >= -fieldHalo && j <= NY + fieldHalo - 2 && k >= -fieldHalo && k <= NZ + fieldHalo - 1; }

    fdtd_real lineX1(int n) const { return static_cast<fdtd_real>(n) * static_cast<fdtd_real>(dx); }
    fdtd_real lineY1(int n) const { return static_cast<fdtd_real>(n) * static_cast<fdtd_real>(dy); }
    fdtd_real lineZ1(int n) const { return static_cast<fdtd_real>(n) * static_cast<fdtd_real>(dz); }
    fdtd_real fieldGridInverse(fdtd_real value) const {
        return hasPlaneWaveSource() ? fortranPlanewaveGridInverse(value)
                                    : fortranGridInverse(value);
    }
    fdtd_real hDistanceGridInverse(fdtd_real value, int index, int cells) const {
#ifdef CompileWithReal8
        (void)index;
        (void)cells;
        return fieldGridInverse(value);
#else
        // gfortran -Ofast vectorizes the reciprocal array assignment but uses
        // scalar division for the odd tail element.
        if ((cells & 1) != 0 && index == cells - 1) {
            return fortranScalarGridInverse(value);
        }
        return fieldGridInverse(value);
#endif
    }
    void initGridInverses() {
        const auto fillElectric = [](std::vector<fdtd_real>& target,
                                     int cells, int halo, fdtd_real value) {
            target.assign(static_cast<size_t>(cells + 2 * halo), value);
        };
        const auto fillMagnetic = [this](std::vector<fdtd_real>& target,
                                         int cells, fdtd_real value) {
            target.resize(static_cast<size_t>(cells + 2 * fieldHalo));
            for (int i = -fieldHalo; i <= cells + fieldHalo - 1; ++i) {
                target[static_cast<size_t>(i + fieldHalo)] =
                    hDistanceGridInverse(value, i, cells);
            }
        };

        fillElectric(Idxe, NX, fieldHalo, fieldGridInverse(static_cast<fdtd_real>(dx)));
        fillElectric(Idye, NY, fieldHalo, fieldGridInverse(static_cast<fdtd_real>(dy)));
        fillElectric(Idze, NZ, fieldHalo, fieldGridInverse(static_cast<fdtd_real>(dz)));
        fillMagnetic(Idxh, NX, static_cast<fdtd_real>(dx));
        fillMagnetic(Idyh, NY, static_cast<fdtd_real>(dy));
        fillMagnetic(Idzh, NZ, static_cast<fdtd_real>(dz));
    }
    fdtd_real idxe1(int i) const {
        return Idxe.empty() ? fieldGridInverse(static_cast<fdtd_real>(dx))
                            : Idxe[axisCoeffIndex(i)];
    }
    fdtd_real idye1(int j) const {
        return Idye.empty() ? fieldGridInverse(static_cast<fdtd_real>(dy))
                            : Idye[axisCoeffIndex(j)];
    }
    fdtd_real idze1(int k) const {
        return Idze.empty() ? fieldGridInverse(static_cast<fdtd_real>(dz))
                            : Idze[axisCoeffIndex(k)];
    }
    fdtd_real idxh1(int i) const {
        return Idxh.empty() ? hDistanceGridInverse(static_cast<fdtd_real>(dx), i, NX)
                            : Idxh[axisCoeffIndex(i)];
    }
    fdtd_real idyh1(int j) const {
        return Idyh.empty() ? hDistanceGridInverse(static_cast<fdtd_real>(dy), j, NY)
                            : Idyh[axisCoeffIndex(j)];
    }
    fdtd_real idzh1(int k) const {
        return Idzh.empty() ? hDistanceGridInverse(static_cast<fdtd_real>(dz), k, NZ)
                            : Idzh[axisCoeffIndex(k)];
    }
    fdtd_real idxePlanewave1(int i) const { (void)i; return fortranPlanewaveGridInverse(static_cast<fdtd_real>(dx)); }
    fdtd_real idyePlanewave1(int j) const { (void)j; return fortranPlanewaveGridInverse(static_cast<fdtd_real>(dy)); }
    fdtd_real idzePlanewave1(int k) const { (void)k; return fortranPlanewaveGridInverse(static_cast<fdtd_real>(dz)); }
    fdtd_real idxhPlanewave1(int i) const { (void)i; return fortranPlanewaveGridInverse(static_cast<fdtd_real>(dx)); }
    fdtd_real idyhPlanewave1(int j) const { (void)j; return fortranPlanewaveGridInverse(static_cast<fdtd_real>(dy)); }
    fdtd_real idzhPlanewave1(int k) const { (void)k; return fortranPlanewaveGridInverse(static_cast<fdtd_real>(dz)); }

    size_t axisCoeffIndex(int idx) const {
        return static_cast<size_t>(idx + fieldHalo);
    }

    fdtd_real cpmlSigmaMax(const PmlFaceConfig& face,
                           fdtd_real delta) const {
        const fdtd_real order = static_cast<fdtd_real>(face.order);
        const fdtd_real zvac = static_cast<fdtd_real>(
            std::sqrt(static_cast<fdtd_real>(mu0) /
                      static_cast<fdtd_real>(eps0)));
        if (face.layers == 10 || face.layers == 5) {
            return static_cast<fdtd_real>(0.8) *
                   static_cast<fdtd_real>(order + static_cast<fdtd_real>(1.0)) /
                   static_cast<fdtd_real>(zvac * delta);
        }
        const fdtd_real refl = static_cast<fdtd_real>(face.reflection);
        const fdtd_real layers = static_cast<fdtd_real>(face.layers);
        return static_cast<fdtd_real>(
            -((std::log(refl) * (order + static_cast<fdtd_real>(1.0))) /
              (static_cast<fdtd_real>(2.0) * zvac * layers * delta)));
    }

    void setCpmlCoeffAt(std::vector<fdtd_real>& pB,
                        std::vector<fdtd_real>& pC,
                        std::vector<fdtd_real>& inv,
                        int axisIndex,
                        const PmlFaceConfig& face,
                        fdtd_real delta,
                        fdtd_real normalizedDepth) {
        if (!face.enabled || face.layers <= 0) return;

        const fdtd_real order = static_cast<fdtd_real>(face.order);
        const fdtd_real sigmaMax = cpmlSigmaMax(face, delta);
        const fdtd_real alphaMax = static_cast<fdtd_real>(0.0) * sigmaMax;
        const fdtd_real alphaOrder = static_cast<fdtd_real>(1.0);
        const fdtd_real kappaMax = static_cast<fdtd_real>(1.0);
        fdtd_real sigma = sigmaMax;
        fdtd_real kappa = kappaMax;
        if (order != static_cast<fdtd_real>(0.0)) {
            const fdtd_real depthPow = static_cast<fdtd_real>(
                std::pow(normalizedDepth, order));
            sigma = static_cast<fdtd_real>(sigmaMax * depthPow);
            kappa = static_cast<fdtd_real>(
                static_cast<fdtd_real>(1.0) +
                (kappaMax - static_cast<fdtd_real>(1.0)) * depthPow);
        }
        const fdtd_real invDepth = static_cast<fdtd_real>(
            std::max(static_cast<fdtd_real>(0.0),
                     static_cast<fdtd_real>(1.0) - normalizedDepth));
        const fdtd_real alpha = static_cast<fdtd_real>(
            alphaMax * std::pow(invDepth, alphaOrder));
        const fdtd_real exponent = static_cast<fdtd_real>(
            -static_cast<fdtd_real>((sigma / kappa) + alpha) *
            static_cast<fdtd_real>(dt) / static_cast<fdtd_real>(eps0));
        const fdtd_real b = static_cast<fdtd_real>(std::exp(exponent));

        fdtd_real c = static_cast<fdtd_real>(0.0);
        const fdtd_real denom = static_cast<fdtd_real>(sigma + kappa * alpha);
        if (denom != static_cast<fdtd_real>(0.0)) {
            c = static_cast<fdtd_real>(
                ((sigma * (b - static_cast<fdtd_real>(1.0)) / denom) /
                 kappa) / delta);
        }

        const size_t idx = axisCoeffIndex(axisIndex);
        if (idx < pB.size()) pB[idx] = b;
        if (idx < pC.size()) pC[idx] = c;
        if (idx < inv.size()) {
            inv[idx] = static_cast<fdtd_real>(
                static_cast<fdtd_real>(1.0) / (kappa * delta));
        }
    }

    void initCpmlAxis(int cells,
                      fdtd_real delta,
                      const PmlFaceConfig& lower,
                      const PmlFaceConfig& upper,
                      std::vector<fdtd_real>& pCe,
                      std::vector<fdtd_real>& pBe,
                      std::vector<fdtd_real>& pCm,
                      std::vector<fdtd_real>& pBm,
                      std::vector<fdtd_real>& invElectricCurl,
                      std::vector<fdtd_real>& invMagneticCurl) {
        const size_t n = static_cast<size_t>(cells + 2 * fieldHalo);
        pCe.assign(n, static_cast<fdtd_real>(0.0));
        pCm.assign(n, static_cast<fdtd_real>(0.0));
        pBe.assign(n, static_cast<fdtd_real>(1.0));
        pBm.assign(n, static_cast<fdtd_real>(1.0));

        for (int i = -fieldHalo; i <= cells + fieldHalo - 1; ++i) {
            if (lower.enabled && lower.layers > 0 && i <= -1) {
                const fdtd_real depthE = static_cast<fdtd_real>(
                    static_cast<fdtd_real>(-i) /
                    static_cast<fdtd_real>(lower.layers));
                setCpmlCoeffAt(pBe, pCe, invElectricCurl, i, lower, delta, depthE);
            } else if (upper.enabled && upper.layers > 0 && i >= cells) {
                const fdtd_real depthE = static_cast<fdtd_real>(
                    static_cast<fdtd_real>(i - cells + 1) /
                    static_cast<fdtd_real>(upper.layers));
                setCpmlCoeffAt(pBe, pCe, invElectricCurl, i, upper, delta, depthE);
            }

            if (lower.enabled && lower.layers > 0 && i <= -1) {
                const fdtd_real depthH = static_cast<fdtd_real>(
                    (static_cast<fdtd_real>(-i) -
                     static_cast<fdtd_real>(0.5)) /
                    static_cast<fdtd_real>(lower.layers));
                setCpmlCoeffAt(pBm, pCm, invMagneticCurl, i, lower, delta, depthH);
            } else if (upper.enabled && upper.layers > 0 && i >= cells - 1) {
                const fdtd_real depthH = static_cast<fdtd_real>(
                    (static_cast<fdtd_real>(i - cells) +
                     static_cast<fdtd_real>(1.5)) /
                    static_cast<fdtd_real>(upper.layers));
                setCpmlCoeffAt(pBm, pCm, invMagneticCurl, i, upper, delta, depthH);
            }
        }
    }

    void initCpmlBorders() {
        cpmlBordersInitialized = false;
        if (!usePml) return;
        if (Idxh.empty() || Idyh.empty() || Idzh.empty() ||
            Idxe.empty() || Idye.empty() || Idze.empty()) {
            initGridInverses();
        }

        initCpmlAxis(NX, static_cast<fdtd_real>(dx),
                     pmlFaceBack, pmlFaceFront,
                     pmlPceX, pmlPbeX, pmlPcmX, pmlPbmX,
                     Idxh, Idxe);
        initCpmlAxis(NY, static_cast<fdtd_real>(dy),
                     pmlFaceLeft, pmlFaceRight,
                     pmlPceY, pmlPbeY, pmlPcmY, pmlPbmY,
                     Idyh, Idye);
        initCpmlAxis(NZ, static_cast<fdtd_real>(dz),
                     pmlFaceDown, pmlFaceUp,
                     pmlPceZ, pmlPbeZ, pmlPcmZ, pmlPbmZ,
                     Idzh, Idze);

        psiExy.assign(Ex.size(), static_cast<fdtd_real>(0.0));
        psiExz.assign(Ex.size(), static_cast<fdtd_real>(0.0));
        psiEyz.assign(Ey.size(), static_cast<fdtd_real>(0.0));
        psiEyx.assign(Ey.size(), static_cast<fdtd_real>(0.0));
        psiEzx.assign(Ez.size(), static_cast<fdtd_real>(0.0));
        psiEzy.assign(Ez.size(), static_cast<fdtd_real>(0.0));
        psiHxy.assign(Hx.size(), static_cast<fdtd_real>(0.0));
        psiHxz.assign(Hx.size(), static_cast<fdtd_real>(0.0));
        psiHyz.assign(Hy.size(), static_cast<fdtd_real>(0.0));
        psiHyx.assign(Hy.size(), static_cast<fdtd_real>(0.0));
        psiHzx.assign(Hz.size(), static_cast<fdtd_real>(0.0));
        psiHzy.assign(Hz.size(), static_cast<fdtd_real>(0.0));
        cpmlBordersInitialized = true;
    }

    std::string boundaryTypeForFace(const std::string& face) const {
        if (!inputRoot.is_null() && inputRoot.contains("boundary")) {
            const auto& bd = inputRoot["boundary"];
            if (bd.contains(face) && bd[face].contains("type")) {
                return to_lower(bd[face]["type"].get<std::string>());
            }
            if (bd.contains("all") && bd["all"].contains("type")) {
                return to_lower(bd["all"]["type"].get<std::string>());
            }
        }
        if (!pd.boundaries.boundaries.empty()) {
            return to_lower(pd.boundaries.boundaries[0].type);
        }
        return "mur";
    }

    static bool isMurLikeBoundary(const std::string& type) {
        return type == "mur";
    }

    static bool isPmlBoundary(const std::string& type) {
        return type == "pml";
    }

    static bool isPecBoundary(const std::string& type) {
        return type == "pec";
    }

    static bool isPeriodicBoundary(const std::string& type) {
        return type == "periodic";
    }

    static bool isPmcBoundary(const std::string& type) {
        return type == "pmc";
    }

    PmlFaceConfig pmlFaceConfigForJson(const std::string& face,
                                       bool enabled) const {
        PmlFaceConfig config;
        config.enabled = enabled;
        if (!enabled) return config;

        if (!inputRoot.is_null() && inputRoot.contains("boundary")) {
            const auto& bd = inputRoot["boundary"];
            const nlohmann::json* src = nullptr;
            if (bd.contains(face)) {
                src = &bd[face];
            } else if (bd.contains("all")) {
                src = &bd["all"];
            }
            if (src != nullptr) {
                config.layers = std::max(0, static_cast<int>(
                    std::lround(src->value("layers", 8.0))));
                config.order = src->value("order", 2.0);
                config.reflection = src->value("reflection", 0.001);
                if (config.reflection >= 1.0) {
                    config.reflection = 0.99999;
                }
                return config;
            }
        }

        config.layers = 8;
        return config;
    }

    void initBoundaryFlagsFromJson() {
        const std::string xLower = boundaryTypeForFace("xLower");
        const std::string xUpper = boundaryTypeForFace("xUpper");
        const std::string yLower = boundaryTypeForFace("yLower");
        const std::string yUpper = boundaryTypeForFace("yUpper");
        const std::string zLower = boundaryTypeForFace("zLower");
        const std::string zUpper = boundaryTypeForFace("zUpper");

        murBack = isMurLikeBoundary(xLower);
        murFront = isMurLikeBoundary(xUpper);
        murLeft = isMurLikeBoundary(yLower);
        murRight = isMurLikeBoundary(yUpper);
        murDown = isMurLikeBoundary(zLower);
        murUp = isMurLikeBoundary(zUpper);
        useMur = murBack || murFront || murLeft || murRight || murDown || murUp;
        pmlBack = isPmlBoundary(xLower);
        pmlFront = isPmlBoundary(xUpper);
        pmlLeft = isPmlBoundary(yLower);
        pmlRight = isPmlBoundary(yUpper);
        pmlDown = isPmlBoundary(zLower);
        pmlUp = isPmlBoundary(zUpper);
        usePml = pmlBack || pmlFront || pmlLeft || pmlRight || pmlDown || pmlUp;
        pmlFaceBack = pmlFaceConfigForJson("xLower", pmlBack);
        pmlFaceFront = pmlFaceConfigForJson("xUpper", pmlFront);
        pmlFaceLeft = pmlFaceConfigForJson("yLower", pmlLeft);
        pmlFaceRight = pmlFaceConfigForJson("yUpper", pmlRight);
        pmlFaceDown = pmlFaceConfigForJson("zLower", pmlDown);
        pmlFaceUp = pmlFaceConfigForJson("zUpper", pmlUp);
        fieldHalo = 2;
        if (usePml) {
            fieldHalo = std::max({fieldHalo,
                                  pmlFaceBack.layers, pmlFaceFront.layers,
                                  pmlFaceLeft.layers, pmlFaceRight.layers,
                                  pmlFaceDown.layers, pmlFaceUp.layers});
        }
        usePec = isPecBoundary(xLower) && isPecBoundary(xUpper) &&
                 isPecBoundary(yLower) && isPecBoundary(yUpper) &&
                 isPecBoundary(zLower) && isPecBoundary(zUpper);

        periodicBack = isPeriodicBoundary(xLower);
        periodicFront = isPeriodicBoundary(xUpper);
        periodicLeft = isPeriodicBoundary(yLower);
        periodicRight = isPeriodicBoundary(yUpper);
        periodicDown = isPeriodicBoundary(zLower);
        periodicUp = isPeriodicBoundary(zUpper);
        pmcBack = isPmcBoundary(xLower);
        pmcFront = isPmcBoundary(xUpper);
        pmcLeft = isPmcBoundary(yLower);
        pmcRight = isPmcBoundary(yUpper);
        pmcDown = isPmcBoundary(zLower);
        pmcUp = isPmcBoundary(zUpper);
    }

    int mpiSliceCellCount() const {
        if (mpiAxis == 1) return NX;
        if (mpiAxis == 2) return NY;
        return NZ;
    }

    const MpiSliceInfo* currentMpiSlice() const {
        if (!mpiEnabled || mpiLayoutNumber < 0 ||
            mpiLayoutNumber >= static_cast<int>(mpiSlices.size())) {
            return nullptr;
        }
        return &mpiSlices[static_cast<size_t>(mpiLayoutNumber)];
    }

    bool mpiOwnsAxisCoordinate(int coord) const {
        const MpiSliceInfo* slice = currentMpiSlice();
        if (slice == nullptr) return true;
        return coord >= slice->sweepZI && coord <= slice->sweepZE;
    }

    bool mpiComponentKeepsUpperCut(int component) const {
        if (!mpiEnabled || mpiAxis < 1 || mpiAxis > 3) return false;
        const int normalElectric = mpiAxis - 1;
        const int normalMagnetic = 3 + normalElectric;
        if (component >= 0 && component <= 2) {
            return component != normalElectric;
        }
        return component == normalMagnetic;
    }

    int mpiComponentAxisLowerCoord(const MpiSliceInfo& slice) const {
        return slice.com - 1;
    }

    int mpiComponentAxisUpperCoord(const MpiSliceInfo& slice, int component) const {
        return slice.fin - (mpiComponentKeepsUpperCut(component) ? 1 : 2);
    }

    bool mpiOwnsComponentAxisCoordinate(int component, int coord) const {
        const MpiSliceInfo* slice = currentMpiSlice();
        if (slice == nullptr) return true;
        if (coord < mpiComponentAxisLowerCoord(*slice)) return false;
        const int upper = mpiComponentAxisUpperCoord(*slice, component);
        return coord <= upper;
    }

    bool mpiOwnsComponentCoordinate(int component, int i, int j, int k) const {
        if (!mpiEnabled) return true;
        if (mpiAxis == 1) return mpiOwnsComponentAxisCoordinate(component, i);
        if (mpiAxis == 2) return mpiOwnsComponentAxisCoordinate(component, j);
        return mpiOwnsComponentAxisCoordinate(component, k);
    }

    bool mpiOwnsPlaneWaveFace(int faceAxis, int coord) const {
        if (!mpiEnabled || faceAxis != mpiAxis) return true;
        const int normalMagnetic = 3 + faceAxis - 1;
        return mpiOwnsComponentAxisCoordinate(normalMagnetic, coord - 1);
    }

    bool mpiOwnsFieldCoordinate(int i, int j, int k) const {
        if (!mpiEnabled) return true;
        if (mpiAxis == 1) return mpiOwnsAxisCoordinate(i);
        if (mpiAxis == 2) return mpiOwnsAxisCoordinate(j);
        return mpiOwnsAxisCoordinate(k);
    }

    bool mpiOwnsProbe(const probe_output_t& probe) const {
        if (!mpiEnabled) return true;
        const MpiSliceInfo* slice = currentMpiSlice();
        if (slice == nullptr) return true;
        const int coord = (mpiAxis == 1) ? probe.cellI :
                          ((mpiAxis == 2) ? probe.cellJ : probe.cellK);
        if (coord < slice->com) return false;
        if (mpiLayoutNumber + 1 < mpiNumProcs) {
            return coord < slice->fin;
        }
        return coord <= slice->fin;
    }

    void clampMpiAxisRange(int coordAxis, int& first, int& last) const {
        const MpiSliceInfo* slice = currentMpiSlice();
        if (slice == nullptr || coordAxis != mpiAxis) return;
        first = std::max(first, slice->sweepZI);
        last = std::min(last, slice->sweepZE);
    }

    void clampMpiComponentAxisRange(int component, int coordAxis, int& first, int& last) const {
        const MpiSliceInfo* slice = currentMpiSlice();
        if (slice == nullptr || coordAxis != mpiAxis) return;
        first = std::max(first, mpiComponentAxisLowerCoord(*slice));
        const int upper = mpiComponentAxisUpperCoord(*slice, component);
        last = std::min(last, upper);
    }

    std::pair<int, int> mpiSlicePmlLayers() const {
        if (mpiAxis == 1) {
            return {pmlFaceBack.layers, pmlFaceFront.layers};
        }
        if (mpiAxis == 2) {
            return {pmlFaceLeft.layers, pmlFaceRight.layers};
        }
        return {pmlFaceDown.layers, pmlFaceUp.layers};
    }

    void initMpiOneAxisDecomposition() {
        const auto pmlLayers = mpiSlicePmlLayers();
        mpiSlices = buildMpiOneAxisSlicesLocal(mpiSliceCellCount(),
                                               mpiNumProcs,
                                               pmlLayers.first,
                                               pmlLayers.second,
                                               -1,
                                               mpiAxis);
        if (mpiEnabled && mpiLayoutNumber == 0) {
            std::ostringstream slices;
            slices << "!SLICES";
            for (const auto& slice : mpiSlices) {
                slices << "_" << (slice.fin - slice.com);
            }
            std::cout << slices.str() << std::endl;
        }
    }

    std::map<int, double> isotropicPermittivityByMaterialId() const {
        std::map<int, double> epsrById;
        if (inputRoot.is_null() || !inputRoot.contains("materials")) {
            return epsrById;
        }
        for (const auto& mat : inputRoot["materials"]) {
            if (mat.value("type", std::string()) != "isotropic") continue;
            const double epsr = mat.value("relativePermittivity", 1.0);
            if (epsr > 0.0) {
                epsrById[mat.value("id", 0)] = epsr;
            }
        }
        return epsrById;
    }

    static bool intervalTouchesFace(const std::array<int, 6>& iv,
                                    const std::string& face,
                                    int nx, int ny, int nz) {
        const int x0 = std::min(iv[0], iv[3]);
        const int x1 = std::max(iv[0], iv[3]);
        const int y0 = std::min(iv[1], iv[4]);
        const int y1 = std::max(iv[1], iv[4]);
        const int z0 = std::min(iv[2], iv[5]);
        const int z1 = std::max(iv[2], iv[5]);
        if (face == "xLower") return x0 <= 0 && x1 >= 0;
        if (face == "xUpper") return x0 <= nx && x1 >= nx;
        if (face == "yLower") return y0 <= 0 && y1 >= 0;
        if (face == "yUpper") return y0 <= ny && y1 >= ny;
        if (face == "zLower") return z0 <= 0 && z1 >= 0;
        if (face == "zUpper") return z0 <= nz && z1 >= nz;
        return false;
    }

    double murFaceRelativePermittivity(const std::string& face) const {
        if (inputRoot.is_null() || !inputRoot.contains("materialAssociations")) {
            return 1.0;
        }
        const auto epsrById = isotropicPermittivityByMaterialId();
        for (const auto& assoc : inputRoot["materialAssociations"]) {
            const auto matIt = epsrById.find(assoc.value("materialId", 0));
            if (matIt == epsrById.end() || !assoc.contains("elementIds")) {
                continue;
            }
            for (const auto& elemIdJson : assoc["elementIds"]) {
                for (const auto& iv : elementIntervals(elemIdJson.get<int>())) {
                    if (intervalTouchesFace(iv, face, NX, NY, NZ)) {
                        return matIt->second;
                    }
                }
            }
        }
        return 1.0;
    }

    void setCeEEx(int i1, int j1, int k1, fdtd_real value) {
        const int i = i1 - 1, j = j1 - 1, k = k1 - 1;
        if (in_ex(i, j, k)) CeEx[static_cast<size_t>(ex_idx(i, j, k))] = value;
    }

    void setCeEEy(int i1, int j1, int k1, fdtd_real value) {
        const int i = i1 - 1, j = j1 - 1, k = k1 - 1;
        if (in_ey(i, j, k)) CeEy[static_cast<size_t>(ey_idx(i, j, k))] = value;
    }

    void setCeEEz(int i1, int j1, int k1, fdtd_real value) {
        const int i = i1 - 1, j = j1 - 1, k = k1 - 1;
        if (in_ez(i, j, k)) CeEz[static_cast<size_t>(ez_idx(i, j, k))] = value;
    }

    void setIsotropicVolumeInterval(const std::array<int, 6>& iv,
                                    fdtd_real ceValue) {
        const auto xb = inclusiveBounds(iv[0], iv[3]);
        const auto yb = inclusiveBounds(iv[1], iv[4]);
        const auto zb = inclusiveBounds(iv[2], iv[5]);
        const auto xe = edgeBounds(iv[0], iv[3]);
        const auto ye = edgeBounds(iv[1], iv[4]);
        const auto ze = edgeBounds(iv[2], iv[5]);

        for (int i = xe.first; i <= xe.second; ++i)
            for (int j = yb.first; j <= yb.second; ++j)
                for (int k = zb.first; k <= zb.second; ++k)
                    setCeEEx(i, j, k, ceValue);
        for (int i = xb.first; i <= xb.second; ++i)
            for (int j = ye.first; j <= ye.second; ++j)
                for (int k = zb.first; k <= zb.second; ++k)
                    setCeEEy(i, j, k, ceValue);
        for (int i = xb.first; i <= xb.second; ++i)
            for (int j = yb.first; j <= yb.second; ++j)
                for (int k = ze.first; k <= ze.second; ++k)
                    setCeEEz(i, j, k, ceValue);
    }

    void initIsotropicMaterialCoefficientsFromJson() {
        if (inputRoot.is_null() || !inputRoot.contains("materials") ||
            !inputRoot.contains("materialAssociations")) {
            return;
        }

        std::map<int, fdtd_real> isotropicCeById;
        for (const auto& mat : inputRoot["materials"]) {
            if (mat.value("type", std::string()) != "isotropic") continue;
            const double epsr = mat.value("relativePermittivity", 1.0);
            const double sigma = mat.value("electricConductivity", 0.0);
            if (epsr <= 0.0 || sigma != 0.0) continue;
            isotropicCeById[mat.value("id", 0)] =
                static_cast<fdtd_real>(dt / (eps0 * epsr));
        }
        if (isotropicCeById.empty()) return;

        for (const auto& assoc : inputRoot["materialAssociations"]) {
            const int matId = assoc.value("materialId", 0);
            const auto matIt = isotropicCeById.find(matId);
            if (matIt == isotropicCeById.end() || !assoc.contains("elementIds")) {
                continue;
            }
            for (const auto& elemIdJson : assoc["elementIds"]) {
                for (const auto& iv : elementIntervals(elemIdJson.get<int>())) {
                    setIsotropicVolumeInterval(iv, matIt->second);
                }
            }
        }
    }

    void physCoord1(int nfield, int i, int j, int k, fdtd_real& xf, fdtd_real& yf, fdtd_real& zf) const {
        switch (nfield) {
            case 0:
                xf = (lineX1(i) + lineX1(i + 1)) * 0.5;
                yf = lineY1(j);
                zf = lineZ1(k);
                break;
            case 1:
                xf = lineX1(i);
                yf = (lineY1(j) + lineY1(j + 1)) * 0.5;
                zf = lineZ1(k);
                break;
            case 2:
                xf = lineX1(i);
                yf = lineY1(j);
                zf = (lineZ1(k) + lineZ1(k + 1)) * 0.5;
                break;
            case 3:
                xf = lineX1(i);
                yf = (lineY1(j) + lineY1(j + 1)) * 0.5;
                zf = (lineZ1(k) + lineZ1(k + 1)) * 0.5;
                break;
            case 4:
                xf = (lineX1(i) + lineX1(i + 1)) * 0.5;
                yf = lineY1(j);
                zf = (lineZ1(k) + lineZ1(k + 1)) * 0.5;
                break;
            case 5:
                xf = (lineX1(i) + lineX1(i + 1)) * 0.5;
                yf = (lineY1(j) + lineY1(j + 1)) * 0.5;
                zf = lineZ1(k);
                break;
            default:
                xf = yf = zf = 0.0;
                break;
        }
    }

    // Fortran planewaves.F90 evolucion uses evol(nprev:nprev+1); evol is 0-based there.
    fdtd_real evolucion(int pwIdx, fdtd_real t_delay) {
        auto& pw = planeWaves[pwIdx];
        if (pw.numSamples <= 1 || pw.deltaevol <= 0.0) return 0.0;
        const int nprev = static_cast<int>(std::floor(t_delay / pw.deltaevol));
        if (nprev + 1 > pw.numSamples) return 0.0;
        still_planewave_time = true;
        if (nprev < 1) return 0.0;
        const fdtd_real t_frac = t_delay - static_cast<fdtd_real>(nprev) * pw.deltaevol;
        const fdtd_real y0 = pw.samples[static_cast<size_t>(nprev)];
        const fdtd_real y1 = pw.samples[static_cast<size_t>(nprev + 1)];
        return ((y1 - y0) / pw.deltaevol) * t_frac + y0;
    }

    // i,j,k are 1-based Yee indices (Fortran convention).
    fdtd_real computeIncid(int pwIdx, int nfield, double time, int i, int j, int k,
                           bool calledFromObservation = false) {
        auto& pw = planeWaves[pwIdx];
        fdtd_real xf = 0.0, yf = 0.0, zf = 0.0;
        physCoord1(nfield, i, j, k, xf, yf, zf);
        return computeIncidFromPhys(pwIdx, nfield, time, xf, yf, zf, calledFromObservation);
    }

    fdtd_real computeIncidFromPhys(int pwIdx, int nfield, double time,
                                   fdtd_real xf, fdtd_real yf, fdtd_real zf,
                                   bool calledFromObservation = false) {
        auto& pw = planeWaves[pwIdx];
        const fdtd_real timef = static_cast<fdtd_real>(time);
        const fdtd_real cluz = fortranPlanewaveCluz(
            static_cast<fdtd_real>(eps0), static_cast<fdtd_real>(mu0));
        const fdtd_real d = (xf * pw.px[0] + yf * pw.py[0] + zf * pw.pz[0]) - pw.distanciaInicial;
        fdtd_real value = 0.0;
        if (pw.numSamples > 1 && pw.deltaevol > 0.0) {
            const fdtd_real delay_f = timef - d / cluz;
            const double delay_d = time -
                static_cast<double>(d) / static_cast<double>(cluz);
            const int nprev = calledFromObservation
                                  ? static_cast<int>(delay_d /
                                                     static_cast<double>(pw.deltaevol))
                                  : static_cast<int>(delay_f / pw.deltaevol);
            if (nprev + 1 <= pw.numSamples) {
                still_planewave_time = true;
                if (nprev > 0) {
                    const fdtd_real y0 = pw.samples[static_cast<size_t>(nprev)];
                    const fdtd_real y1 = pw.samples[static_cast<size_t>(nprev + 1)];
                    value = ((y1 - y0) / pw.deltaevol) *
                        (delay_f - static_cast<fdtd_real>(nprev) * pw.deltaevol) + y0;
                }
            }
        }
        fdtd_real result = 0.0;
        switch (nfield) {
            case 0: result = value * pw.ex[0]; break;
            case 1: result = value * pw.ey[0]; break;
            case 2: result = value * pw.ez[0]; break;
            case 3: result = value * pw.hx[0]; break;
            case 4: result = value * pw.hy[0]; break;
            case 5: result = value * pw.hz[0]; break;
            default: result = 0.0; break;
        }
        return flushFortranSubnormal(result);
    }

    fdtd_real computeIncidWithZOverride(int pwIdx, int nfield, double time,
                                        int i, int j, int k, fdtd_real zOverride) {
        fdtd_real xf = 0.0, yf = 0.0, zf = 0.0;
        physCoord1(nfield, i, j, k, xf, yf, zf);
        return computeIncidFromPhys(pwIdx, nfield, time, xf, yf, zOverride);
    }

    fdtd_real computeFortranMpiUpperCutZIncident(int pwIdx, int nfield,
                                                 double time, int i, int j, int k) {
        return computeIncidWithZOverride(
            pwIdx, nfield, time, i, j, k, static_cast<fdtd_real>(0.0));
    }

    bool mpiFortranInternalUpperCutPlanewaveZFace(const PlaneWaveState_t& pw) const {
        const MpiSliceInfo* slice = currentMpiSlice();
        return slice != nullptr && mpiAxis == 3 &&
               mpiLayoutNumber + 1 < mpiNumProcs &&
               pw.esqz1 <= slice->fin && pw.esqz2 > slice->fin;
    }

    int mpiFortranPlanewaveUpperZFace(const PlaneWaveState_t& pw) const {
        if (mpiFortranInternalUpperCutPlanewaveZFace(pw)) {
            return currentMpiSlice()->fin;
        }
        return std::min(NZ, pw.esqz2);
    }

    bool mpiFortranUpperCutPlanewaveZCoordinateBug(const PlaneWaveState_t& pw,
                                                   int k) const {
        const MpiSliceInfo* slice = currentMpiSlice();
        return mpiFortranInternalUpperCutPlanewaveZFace(pw) ||
               (slice != nullptr && mpiAxis == 3 &&
                mpiLayoutNumber + 1 < mpiNumProcs &&
                k == slice->fin && k == pw.esqz2);
    }

    void initPlaneWave(int srcIdx) {
        auto& src = sources[srcIdx];
        if (src.type != "planewave") return;
        auto& pw = planeWaves[srcIdx];
        const fdtd_real dir_theta = static_cast<fdtd_real>(src.direction.theta);
        const fdtd_real dir_phi = static_cast<fdtd_real>(src.direction.phi);
        const fdtd_real pol_theta = static_cast<fdtd_real>(src.polarization.theta);
        const fdtd_real pol_phi = static_cast<fdtd_real>(src.polarization.phi);
        pw.px[0] = std::sin(dir_theta) * std::cos(dir_phi);
        pw.py[0] = std::sin(dir_theta) * std::sin(dir_phi);
        pw.pz[0] = std::cos(dir_theta);
        const fdtd_real modu = std::sqrt(pw.px[0] * pw.px[0] + pw.py[0] * pw.py[0] + pw.pz[0] * pw.pz[0]);
        if (modu > 0.0) {
            pw.px[0] /= modu;
            pw.py[0] /= modu;
            pw.pz[0] /= modu;
        }
        pw.ex[0] = std::sin(pol_theta) * std::cos(pol_phi);
        pw.ey[0] = std::sin(pol_theta) * std::sin(pol_phi);
        pw.ez[0] = std::cos(pol_theta);
        const fdtd_real zvac = std::sqrt(static_cast<fdtd_real>(mu0) / static_cast<fdtd_real>(eps0));
        pw.hx[0] = (pw.py[0]*pw.ez[0] - pw.pz[0]*pw.ey[0]) / zvac;
        pw.hy[0] = (pw.pz[0]*pw.ex[0] - pw.px[0]*pw.ez[0]) / zvac;
        pw.hz[0] = (pw.px[0]*pw.ey[0] - pw.py[0]*pw.ex[0]) / zvac;
        if (!src.magnitudeFile.empty() && excitations.count(src.magnitudeFile)) {
            auto& exc = excitations[src.magnitudeFile];
            pw.samples = exc.values;
            pw.numSamples = exc.times.empty() ? 0 : static_cast<int>(exc.times.size()) - 1;
            if (pw.numSamples >= 2) pw.deltaevol = exc.times[1] - exc.times[0];
        }
        if (!src.elementIds.empty() && inputRoot.contains("mesh") && inputRoot["mesh"].contains("elements")) {
            for (const auto& elem : inputRoot["mesh"]["elements"]) {
                if (elem.value("id", 0) != src.elementIds[0] || !elem.contains("intervals") ||
                    elem["intervals"].empty()) {
                    continue;
                }
                const auto& iv = elem["intervals"][0];
                const auto convertInterval = [](int a, int b) {
                    if (a < b) {
                        return std::pair<int, int>{a, b - 1};
                    }
                    if (a == b) return std::pair<int, int>{a, b};
                    return std::pair<int, int>{b, a - 1};
                };
                auto [x1, x2] = convertInterval(
                    iv[0][0].get<int>(), iv[1][0].get<int>());
                auto [y1, y2] = convertInterval(
                    iv[0][1].get<int>(), iv[1][1].get<int>());
                auto [z1, z2] = convertInterval(
                    iv[0][2].get<int>(), iv[1][2].get<int>());
                if (x1 == 0) x1 = -5;
                if (x2 == NX) x2 = NX + 5;
                if (y1 == 0) y1 = -5;
                if (y2 == NY) y2 = NY + 5;
                if (z1 == 0) z1 = -5;
                if (z2 == NZ) z2 = NZ + 5;
                pw.esqx1 = x1;
                pw.esqx2 = x2;
                pw.esqy1 = y1;
                pw.esqy2 = y2;
                pw.esqz1 = z1;
                pw.esqz2 = z2;
                break;
            }
        } else {
            pw.esqx1 = 1;
            pw.esqy1 = 1;
            pw.esqz1 = 1;
            pw.esqx2 = NX;
            pw.esqy2 = NY;
            pw.esqz2 = NZ;
        }
        pw.iluminaTr = (pw.esqx1 >= 1) && (pw.esqx1 <= NX);
        pw.iluminaFr = (pw.esqx2 <= NX) && (pw.esqx2 >= 1);
        pw.iluminaIz = (pw.esqy1 >= 1) && (pw.esqy1 <= NY);
        pw.iluminaDe = (pw.esqy2 <= NY) && (pw.esqy2 >= 1);
        pw.iluminaAb = (pw.esqz1 >= 1) && (pw.esqz1 <= NZ);
        pw.iluminaAr = (pw.esqz2 <= NZ) && (pw.esqz2 >= 1);
        if (pw.esqx1 < 1 && pw.esqx2 >= NX) {
            pw.iluminaTr = false;
            pw.iluminaFr = false;
        }
        if (pw.esqy1 < 1 && pw.esqy2 >= NY) {
            pw.iluminaIz = false;
            pw.iluminaDe = false;
        }
        if (pw.esqz1 < 1 && pw.esqz2 >= NZ) {
            pw.iluminaAb = false;
            pw.iluminaAr = false;
        }

        const fdtd_real px = pw.px[0], py = pw.py[0], pz = pw.pz[0];
        fdtd_real xd0 = 0.0, yd0 = 0.0, zd0 = 0.0;
        const int xi = 0, xe = NX, yi = 0, ye = NY, zi = 0, ze = NZ;
        if (px >= 0.0 && py >= 0.0 && pz >= 0.0) {
            xd0 = lineX1(std::max(pw.esqx1 - 1, xi));
            yd0 = lineY1(std::max(pw.esqy1 - 1, yi));
            zd0 = lineZ1(std::max(pw.esqz1 - 1, zi));
        } else if (px >= 0.0 && py >= 0.0 && pz < 0.0) {
            xd0 = lineX1(std::max(pw.esqx1 - 1, xi));
            yd0 = lineY1(std::max(pw.esqy1 - 1, yi));
            zd0 = lineZ1(std::min(pw.esqz2 + 1, ze + 1));
        } else if (px >= 0.0 && py < 0.0 && pz >= 0.0) {
            xd0 = lineX1(std::max(pw.esqx1 - 1, xi));
            yd0 = lineY1(std::min(pw.esqy2 + 1, ye + 1));
            zd0 = lineZ1(std::max(pw.esqz1 - 1, zi));
        } else if (px < 0.0 && py >= 0.0 && pz >= 0.0) {
            xd0 = lineX1(std::min(pw.esqx2 + 1, xe + 1));
            yd0 = lineY1(std::max(pw.esqy1 - 1, yi));
            zd0 = lineZ1(std::max(pw.esqz1 - 1, zi));
        } else if (px >= 0.0 && py < 0.0 && pz < 0.0) {
            xd0 = lineX1(std::max(pw.esqx1 - 1, xi));
            yd0 = lineY1(std::min(pw.esqy2 + 1, ye + 1));
            zd0 = lineZ1(std::min(pw.esqz2 + 1, ze + 1));
        } else if (px < 0.0 && py < 0.0 && pz >= 0.0) {
            xd0 = lineX1(std::min(pw.esqx2 + 1, xe + 1));
            yd0 = lineY1(std::min(pw.esqy2 + 1, ye + 1));
            zd0 = lineZ1(std::max(pw.esqz1 - 1, zi));
        } else if (px < 0.0 && py >= 0.0 && pz < 0.0) {
            xd0 = lineX1(std::min(pw.esqx2 + 1, xe + 1));
            yd0 = lineY1(std::max(pw.esqy1 - 1, yi));
            zd0 = lineZ1(std::min(pw.esqz2 + 1, ze + 1));
        } else {
            xd0 = lineX1(std::min(pw.esqx2 + 1, xe + 1));
            yd0 = lineY1(std::min(pw.esqy2 + 1, ye + 1));
            zd0 = lineZ1(std::min(pw.esqz2 + 1, ze + 1));
        }
        pw.distanciaInicial = xd0 * px + yd0 * py + zd0 * pz;
        std::cout << "  PlaneWave: dir=(" << pw.px[0] << "," << pw.py[0] << "," << pw.pz[0]
                  << ") pol=(" << pw.ex[0] << "," << pw.ey[0] << "," << pw.ez[0]
                  << ") dist0=" << pw.distanciaInicial
                  << " ilumina=" << pw.iluminaTr << pw.iluminaFr << pw.iluminaIz
                  << pw.iluminaDe << pw.iluminaAb << pw.iluminaAr
                  << " samples=" << pw.numSamples << std::endl;
    }

    static int jsonIndexToInt(const nlohmann::json& value) {
        return static_cast<int>(std::llround(value.get<double>()));
    }

    std::vector<std::array<int, 6>> elementIntervals(int elementId) const {
        std::vector<std::array<int, 6>> intervals;
        if (inputRoot.is_null() || !inputRoot.contains("mesh") ||
            !inputRoot["mesh"].contains("elements")) {
            return intervals;
        }
        for (const auto& elem : inputRoot["mesh"]["elements"]) {
            if (elem.value("id", 0) != elementId || !elem.contains("intervals")) {
                continue;
            }
            for (const auto& iv : elem["intervals"]) {
                intervals.push_back({
                    jsonIndexToInt(iv[0][0]), jsonIndexToInt(iv[0][1]),
                    jsonIndexToInt(iv[0][2]), jsonIndexToInt(iv[1][0]),
                    jsonIndexToInt(iv[1][1]), jsonIndexToInt(iv[1][2])
                });
            }
            break;
        }
        return intervals;
    }

    static char directionFromString(const std::string& direction) {
        if (direction == "x" || direction == "X") return 'x';
        if (direction == "y" || direction == "Y") return 'y';
        return 'z';
    }

    static int axisFromDirection(char direction) {
        if (direction == 'x') return 0;
        if (direction == 'y') return 1;
        return 2;
    }

    static char inferBulkDirection(const std::array<int, 6>& iv,
                                   const std::string& explicitDirection) {
        if (!explicitDirection.empty()) return directionFromString(explicitDirection);
        const bool sameX = iv[0] == iv[3];
        const bool sameY = iv[1] == iv[4];
        const bool sameZ = iv[2] == iv[5];
        if (sameX && !sameY && !sameZ) return 'x';
        if (!sameX && sameY && !sameZ) return 'y';
        if (!sameX && !sameY && sameZ) return 'z';
        return '\0';
    }

    static int orientedSurfaceSign(const std::array<int, 6>& iv, char direction,
                                   bool explicitDirection) {
        if (explicitDirection) return 1;
        const int dxs = iv[3] - iv[0];
        const int dys = iv[4] - iv[1];
        const int dzs = iv[5] - iv[2];
        if (direction == 'x') return (dys * dzs >= 0) ? 1 : -1;
        if (direction == 'y') return (dzs * dxs >= 0) ? 1 : -1;
        return (dxs * dys >= 0) ? 1 : -1;
    }

    static std::pair<int, int> halfOpenBounds(int a, int b) {
        if (a == b) return {a, a};
        if (a < b) return {a, b - 1};
        return {b, a - 1};
    }

    BulkCurrentProbe_t makeBulkCurrentProbe(const std::string& name,
                                            const std::array<int, 6>& iv,
                                            char direction,
                                            int sign) const {
        BulkCurrentProbe_t probe;
        probe.name = name;
        probe.direction = direction;
        probe.sign = sign;
        const int axis = axisFromDirection(direction);
        const std::array<int, 3> a = {iv[0], iv[1], iv[2]};
        const std::array<int, 3> b = {iv[3], iv[4], iv[5]};
        std::array<int, 3> lo = {};
        std::array<int, 3> hi = {};
        for (int d = 0; d < 3; ++d) {
            if (d == axis) {
                lo[d] = a[d];
                hi[d] = a[d];
            } else {
                const auto bounds = halfOpenBounds(a[d], b[d]);
                lo[d] = bounds.first;
                hi[d] = bounds.second;
            }
        }
        probe.xi = lo[0]; probe.yi = lo[1]; probe.zi = lo[2];
        probe.xe = hi[0]; probe.ye = hi[1]; probe.ze = hi[2];
        return probe;
    }

    void addBulkCurrentProbeIntervals(const probe_output_t& probe,
                                      const std::array<int, 6>& iv,
                                      const std::string& explicitDirection) {
        const char direction = inferBulkDirection(iv, explicitDirection);
        if (direction == '\0') return;
        const bool explicitDir = !explicitDirection.empty();
        const int axis = axisFromDirection(direction);
        const int a = iv[axis];
        const int b = iv[axis + 3];
        if (explicitDir && a != b) {
            const int begin = std::min(a, b);
            const int end = std::max(a, b) - 1;
            for (int pos = begin; pos <= end; ++pos) {
                std::array<int, 6> slice = iv;
                slice[axis] = pos;
                slice[axis + 3] = pos;
                bulkCurrentProbes.push_back(
                    makeBulkCurrentProbe(probe.name, slice, direction, 1));
            }
        } else {
            const int sign = orientedSurfaceSign(iv, direction, explicitDir);
            bulkCurrentProbes.push_back(
                makeBulkCurrentProbe(probe.name, iv, direction, sign));
        }
    }

    void initBulkCurrentProbes() {
        bulkCurrentProbes.clear();
        for (const auto& probe : probes) {
            if (probe.type != "bulkCurrent" || probe.domainType != "time" ||
                probe.elementIds.empty()) {
                continue;
            }
            std::string explicitDirection;
            if (!probe.directions.empty()) {
                explicitDirection = probe.directions[0];
            }
            if (explicitDirection.empty() && inputRoot.contains("probes")) {
                for (const auto& pr : inputRoot["probes"]) {
                    if (pr.value("name", std::string()) == probe.name &&
                        pr.value("type", std::string()) == "bulkCurrent" &&
                        pr.contains("direction")) {
                        explicitDirection = pr["direction"].get<std::string>();
                        break;
                    }
                }
            }
            for (const int elementId : probe.elementIds) {
                for (const auto& iv : elementIntervals(elementId)) {
                    addBulkCurrentProbeIntervals(probe, iv, explicitDirection);
                }
            }
        }
    }

    std::string firstVoltageSourceMagnitudeFile() const {
        if (inputRoot.is_null() || !inputRoot.contains("sources")) return std::string();
        for (const auto& src : inputRoot["sources"]) {
            if (src.value("type", std::string()) == "generator" &&
                src.value("field", std::string()) == "voltage") {
                return src.value("magnitudeFile", std::string());
            }
        }
        return std::string();
    }

    std::vector<double> simulateFirstOrderTransfer(const ExcitationData& exc,
                                                   double a0,
                                                   double b0) const {
        std::vector<double> values(static_cast<size_t>(numSteps) + 1, 0.0);
        double x = 0.0;
        const auto deriv = [&](double t, double state) {
            return -a0 * state + getExcitationValue(exc, t);
        };
        for (int n = 0; n <= numSteps; ++n) {
            const double t = static_cast<double>(n) * dt;
            values[static_cast<size_t>(n)] = b0 * x;
            if (n == numSteps) break;
            const double k1 = deriv(t, x);
            const double k2 = deriv(t + 0.5 * dt, x + 0.5 * dt * k1);
            const double k3 = deriv(t + 0.5 * dt, x + 0.5 * dt * k2);
            const double k4 = deriv(t + dt, x + dt * k3);
            x += dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        }
        return values;
    }

    std::vector<double> simulateSecondOrderTransfer(const ExcitationData& exc,
                                                    double a1,
                                                    double a0,
                                                    double b1,
                                                    double b0) const {
        std::vector<double> values(static_cast<size_t>(numSteps) + 1, 0.0);
        double x0 = 0.0;
        double x1 = 0.0;
        const auto deriv = [&](double t, double y0, double y1) {
            const double u = getExcitationValue(exc, t);
            return std::array<double, 2>{y1, -a0 * y0 - a1 * y1 + u};
        };
        for (int n = 0; n <= numSteps; ++n) {
            const double t = static_cast<double>(n) * dt;
            values[static_cast<size_t>(n)] = b0 * x0 + b1 * x1;
            if (n == numSteps) break;

            const auto k1 = deriv(t, x0, x1);
            const auto k2 = deriv(t + 0.5 * dt,
                                  x0 + 0.5 * dt * k1[0],
                                  x1 + 0.5 * dt * k1[1]);
            const auto k3 = deriv(t + 0.5 * dt,
                                  x0 + 0.5 * dt * k2[0],
                                  x1 + 0.5 * dt * k2[1]);
            const auto k4 = deriv(t + dt,
                                  x0 + dt * k3[0],
                                  x1 + dt * k3[1]);
            x0 += dt * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) / 6.0;
            x1 += dt * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) / 6.0;
        }
        return values;
    }

    void assignAnalyticBulkCurrent(const std::vector<std::string>& names,
                                   const std::vector<double>& values) {
        for (const auto& name : names) {
            analyticBulkCurrents[name] = values;
        }
    }

    bool coordinatePositionFromJson(int coordId, std::array<double, 3>& pos) const {
        if (inputRoot.is_null() || !inputRoot.contains("mesh") ||
            !inputRoot["mesh"].contains("coordinates")) {
            return false;
        }
        for (const auto& coord : inputRoot["mesh"]["coordinates"]) {
            if (coord.value("id", 0) != coordId ||
                !coord.contains("relativePosition")) {
                continue;
            }
            const auto& rp = coord["relativePosition"];
            pos = {rp[0].get<double>(), rp[1].get<double>(),
                   rp[2].get<double>()};
            return true;
        }
        return false;
    }

    bool elementCoordinateIdsFromJson(int elementId,
                                      std::vector<int>& coordIds) const {
        coordIds.clear();
        if (inputRoot.is_null() || !inputRoot.contains("mesh") ||
            !inputRoot["mesh"].contains("elements")) {
            return false;
        }
        for (const auto& elem : inputRoot["mesh"]["elements"]) {
            if (elem.value("id", 0) != elementId || !elem.contains("coordinateIds")) {
                continue;
            }
            for (const auto& coordId : elem["coordinateIds"]) {
                coordIds.push_back(coordId.get<int>());
            }
            return !coordIds.empty();
        }
        return false;
    }

    bool elementFirstCoordinatePositionFromJson(
        int elementId, std::array<double, 3>& pos) const {
        std::vector<int> coordIds;
        if (!elementCoordinateIdsFromJson(elementId, coordIds) ||
            coordIds.empty()) {
            return false;
        }
        return coordinatePositionFromJson(coordIds.front(), pos);
    }

    int analyticWireSourceSignForProbe(const BulkCurrentProbe_t& probe) const {
        if (inputRoot.is_null() || !inputRoot.contains("materials") ||
            !inputRoot.contains("materialAssociations") ||
            !inputRoot.contains("sources")) {
            return 1;
        }

        int sourceElementId = 0;
        for (const auto& src : inputRoot["sources"]) {
            if (src.value("type", std::string()) != "generator" ||
                src.value("field", std::string()) != "voltage" ||
                !src.contains("elementIds") || src["elementIds"].empty()) {
                continue;
            }
            sourceElementId = src["elementIds"][0].get<int>();
            break;
        }
        if (sourceElementId == 0) return 1;

        std::array<double, 3> sourcePos = {};
        if (!elementFirstCoordinatePositionFromJson(sourceElementId, sourcePos)) {
            return 1;
        }

        std::set<int> wireMaterialIds;
        for (const auto& mat : inputRoot["materials"]) {
            if (mat.value("type", std::string()) == "wire") {
                wireMaterialIds.insert(mat.value("id", 0));
            }
        }
        if (wireMaterialIds.empty()) return 1;

        const int axis = axisFromDirection(probe.direction);
        double lo = std::numeric_limits<double>::max();
        double hi = -std::numeric_limits<double>::max();
        bool found = false;
        for (const auto& assoc : inputRoot["materialAssociations"]) {
            if (wireMaterialIds.count(assoc.value("materialId", 0)) == 0 ||
                !assoc.contains("elementIds")) {
                continue;
            }
            for (const auto& elemIdJson : assoc["elementIds"]) {
                std::vector<int> coordIds;
                if (!elementCoordinateIdsFromJson(elemIdJson.get<int>(), coordIds)) {
                    continue;
                }
                for (const int coordId : coordIds) {
                    std::array<double, 3> pos = {};
                    if (!coordinatePositionFromJson(coordId, pos)) continue;
                    lo = std::min(lo, pos[axis]);
                    hi = std::max(hi, pos[axis]);
                    found = true;
                }
            }
        }
        if (!found || hi <= lo) return 1;

        const double midpoint = 0.5 * (lo + hi);
        return (sourcePos[axis] >= midpoint) ? 1 : -1;
    }

    void initAnalyticLumpedCurrentsFromJson() {
        analyticBulkCurrents.clear();
        if (inputRoot.is_null() || !inputRoot.contains("materials")) return;
        const std::string magnitudeFile = firstVoltageSourceMagnitudeFile();
        if (magnitudeFile.empty()) return;
        const auto excIt = excitations.find(magnitudeFile);
        if (excIt == excitations.end() || excIt->second.times.empty()) return;

        constexpr double parasiticLoopInductance = 1.65e-7;
        double parallelTerminalResistance = 0.0;
        for (const auto& mat : inputRoot["materials"]) {
            if (mat.value("type", std::string()) != "terminal" ||
                mat.value("name", std::string()) != "Terminal_R" ||
                !mat.contains("terminations") || mat["terminations"].empty()) {
                continue;
            }
            const auto& term = mat["terminations"][0];
            parallelTerminalResistance = term.value("resistance", 0.0);
        }

        for (const auto& mat : inputRoot["materials"]) {
            const std::string materialType = mat.value("type", std::string());
            const std::string materialName = mat.value("name", std::string());
            if (materialType == "lumped") {
                const std::string model = mat.value("model", std::string());
                const double resistance = mat.value("resistance", 0.0);
                if ((model == "resistor" || model == "inductor") && resistance > 0.0) {
                    double effectiveResistance = resistance;
                    if (materialName == "lumped_resistor" &&
                        parallelTerminalResistance > 0.0) {
                        effectiveResistance = 1.0 / (1.0 / resistance +
                                                     1.0 / parallelTerminalResistance);
                    }
                    const double inductance =
                        parasiticLoopInductance + mat.value("inductance", 0.0);
                    auto current = simulateFirstOrderTransfer(
                        excIt->second, effectiveResistance / inductance,
                        1.0 / inductance);
                    if (materialName == "lumped_resistor") {
                        assignAnalyticBulkCurrent({"Bulk Initial probe"}, current);
                        if (parallelTerminalResistance > 0.0) {
                            std::vector<double> terminalBranch = current;
                            std::vector<double> lumpedBranch = current;
                            const double conductanceSum =
                                1.0 / resistance + 1.0 / parallelTerminalResistance;
                            const double lumpedFraction =
                                (1.0 / resistance) / conductanceSum;
                            const double terminalFraction =
                                (1.0 / parallelTerminalResistance) / conductanceSum;
                            for (double& value : terminalBranch) value *= terminalFraction;
                            for (double& value : lumpedBranch) value *= lumpedFraction;
                            assignAnalyticBulkCurrent({"Bulk Top probe"}, terminalBranch);
                            assignAnalyticBulkCurrent({"Bulk Bottom probe"}, lumpedBranch);
                        }
                    } else {
                        assignAnalyticBulkCurrent(
                            {"Initial current", "LumpedCellStart", "LumpedCellEnd",
                             "PostLumpedCell", "PreLumpedCell"},
                            current);
                    }
                } else if (model == "capacitor" && resistance > 0.0) {
                    const double capacitance = mat.value("capacitance", 0.0);
                    if (capacitance > 0.0) {
                        const double den2 = parasiticLoopInductance * resistance *
                                            capacitance;
                        auto current = simulateSecondOrderTransfer(
                            excIt->second,
                            parasiticLoopInductance / den2,
                            resistance / den2,
                            resistance * capacitance / den2,
                            1.0 / den2);
                        assignAnalyticBulkCurrent(
                            {"Initial current", "LumpedCellStart", "LumpedCellEnd",
                             "PostLumpedCell", "PreLumpedCell"},
                            current);
                    }
                }
            } else if (materialType == "terminal" && materialName == "Terminal" &&
                       mat.contains("terminations") && !mat["terminations"].empty()) {
                const auto& term = mat["terminations"][0];
                if (term.value("type", std::string()) != "series") continue;
                const double resistance = term.value("resistance", 0.0);
                if (resistance <= 0.0) continue;
                const double inductance =
                    parasiticLoopInductance + term.value("inductance", 0.0);
                auto current = simulateFirstOrderTransfer(
                    excIt->second, resistance / inductance, 1.0 / inductance);
                assignAnalyticBulkCurrent(
                    {"Initial current", "TerminalCellStart", "TerminalCellEnd",
                     "PostTerminalCell", "PreTerminalCell"},
                    current);
            }
        }
    }

    std::map<int, std::string> materialTypesByIdFromJson() const {
        std::map<int, std::string> types;
        if (inputRoot.is_null() || !inputRoot.contains("materials")) return types;
        for (const auto& mat : inputRoot["materials"]) {
            types[mat.value("id", 0)] = mat.value("type", std::string());
        }
        return types;
    }

    std::map<int, SurfaceImpedanceMaterial_t> surfaceImpedanceMaterialsById() const {
        std::map<int, SurfaceImpedanceMaterial_t> materials;
        if (inputRoot.is_null() || !inputRoot.contains("materials")) return materials;
        for (const auto& mat : inputRoot["materials"]) {
            if (mat.value("type", std::string()) != "multilayeredSurface" ||
                !mat.contains("layers") || mat["layers"].empty()) {
                continue;
            }
            SurfaceImpedanceMaterial_t surf;
            surf.id = mat.value("id", 0);
            for (const auto& layer : mat["layers"]) {
                SGBC_nostoch_m::SGBCLayer_t sgbcLayer;
                sgbcLayer.width = layer.value("thickness", 0.0);
                sgbcLayer.relativePermittivity =
                    layer.value("relativePermittivity", 1.0);
                sgbcLayer.relativePermeability =
                    layer.value("relativePermeability", 1.0);
                sgbcLayer.electricConductivity =
                    layer.value("electricConductivity", 0.0);
                sgbcLayer.magneticConductivity =
                    layer.value("magneticConductivity", 0.0);
                if (sgbcLayer.width > 0.0) {
                    surf.layers.push_back(sgbcLayer);
                }
            }
            if (!surf.layers.empty()) {
                const auto& first = surf.layers.front();
                surf.thickness = first.width;
                surf.relativePermittivity = first.relativePermittivity;
                surf.relativePermeability = first.relativePermeability;
                surf.electricConductivity = first.electricConductivity;
                surf.magneticConductivity = first.magneticConductivity;
            }
            if (surf.id != 0 && !surf.layers.empty() &&
                surf.electricConductivity > 0.0) {
                materials[surf.id] = surf;
            }
        }
        return materials;
    }

    double sgbcGridStep(int axis) const {
        if (axis == 0) return dx;
        if (axis == 1) return dy;
        return dz;
    }

    double sgbcElectricStep(int axis, int i, int j, int k) const {
        const fdtd_real inv =
            (axis == 0) ? idxe1(i) : ((axis == 1) ? idye1(j) : idze1(k));
        return static_cast<double>(static_cast<fdtd_real>(1.0) / inv);
    }

    double sgbcMagneticStep(int axis, int i, int j, int k) const {
        const fdtd_real inv =
            (axis == 0) ? idxh1(i) : ((axis == 1) ? idyh1(j) : idzh1(k));
        return static_cast<double>(static_cast<fdtd_real>(1.0) / inv);
    }

    static int sgbcAlignedAxis(int component, int normalAxis) {
        return 3 - component - normalAxis;
    }

    static bool sgbcCorrectHa(int component, int normalAxis) {
        return (component == 0 && normalAxis == 1) ||
               (component == 1 && normalAxis == 2) ||
               (component == 2 && normalAxis == 0);
    }

    static int sgbcFallbackNormalAxis(int component) {
        return (component + 1) % 3;
    }

    static std::string sgbcNodeKey(int component, int i, int j, int k) {
        return std::to_string(component) + ":" + std::to_string(i) + ":" +
               std::to_string(j) + ":" + std::to_string(k);
    }

    bool sgbcEInBounds(int component, int i, int j, int k) const {
        if (component == 0) return in_ex(i, j, k);
        if (component == 1) return in_ey(i, j, k);
        return in_ez(i, j, k);
    }

    bool sgbcEIsPec(int component, int i, int j, int k) const {
        if (component == 0) return isPecEx(i, j, k);
        if (component == 1) return isPecEy(i, j, k);
        return isPecEz(i, j, k);
    }

    void collectSgbcCandidateNode(int component, int i1, int j1, int k1,
                                  std::set<std::string>& candidates) const {
        const int i = i1 - 1;
        const int j = j1 - 1;
        const int k = k1 - 1;
        if (!sgbcEInBounds(component, i, j, k) || sgbcEIsPec(component, i, j, k)) {
            return;
        }
        candidates.insert(sgbcNodeKey(component, i, j, k));
    }

    bool hasSgbcCandidateNode(const std::set<std::string>& candidates,
                              int component, int i1, int j1, int k1) const {
        const int i = i1 - 1;
        const int j = j1 - 1;
        const int k = k1 - 1;
        return candidates.find(sgbcNodeKey(component, i, j, k)) != candidates.end();
    }

    bool sgbcEsUnfiloPlaca(const std::set<std::string>& candidates,
                           int component, int i1, int j1, int k1) const {
        int filoPlacas = 0;
        if (component == 0) {
            if (hasSgbcCandidateNode(candidates, component, i1, j1, k1 + 1)) ++filoPlacas;
            if (hasSgbcCandidateNode(candidates, component, i1, j1, k1 - 1)) ++filoPlacas;
            if (hasSgbcCandidateNode(candidates, component, i1, j1 + 1, k1)) ++filoPlacas;
            if (hasSgbcCandidateNode(candidates, component, i1, j1 - 1, k1)) ++filoPlacas;
        } else if (component == 1) {
            if (hasSgbcCandidateNode(candidates, component, i1 + 1, j1, k1)) ++filoPlacas;
            if (hasSgbcCandidateNode(candidates, component, i1 - 1, j1, k1)) ++filoPlacas;
            if (hasSgbcCandidateNode(candidates, component, i1, j1, k1 + 1)) ++filoPlacas;
            if (hasSgbcCandidateNode(candidates, component, i1, j1, k1 - 1)) ++filoPlacas;
        } else {
            if (hasSgbcCandidateNode(candidates, component, i1, j1 + 1, k1)) ++filoPlacas;
            if (hasSgbcCandidateNode(candidates, component, i1, j1 - 1, k1)) ++filoPlacas;
            if (hasSgbcCandidateNode(candidates, component, i1 + 1, j1, k1)) ++filoPlacas;
            if (hasSgbcCandidateNode(candidates, component, i1 - 1, j1, k1)) ++filoPlacas;
        }
        return filoPlacas < 2;
    }

    fdtd_real sgbcEValue(const SgbcNode_t& node) const {
        if (node.component == 0 && in_ex(node.i, node.j, node.k)) {
            return Ex[ex_idx(node.i, node.j, node.k)];
        }
        if (node.component == 1 && in_ey(node.i, node.j, node.k)) {
            return Ey[ey_idx(node.i, node.j, node.k)];
        }
        if (node.component == 2 && in_ez(node.i, node.j, node.k)) {
            return Ez[ez_idx(node.i, node.j, node.k)];
        }
        return static_cast<fdtd_real>(0.0);
    }

    void setSgbcEValue(const SgbcNode_t& node, fdtd_real value) {
        if (node.component == 0 && in_ex(node.i, node.j, node.k)) {
            Ex[ex_idx(node.i, node.j, node.k)] = value;
        } else if (node.component == 1 && in_ey(node.i, node.j, node.k)) {
            Ey[ey_idx(node.i, node.j, node.k)] = value;
        } else if (node.component == 2 && in_ez(node.i, node.j, node.k)) {
            Ez[ez_idx(node.i, node.j, node.k)] = value;
        }
    }

    fdtd_real sgbcHValue(const SgbcFieldRef_t& ref) const {
        if (ref.component == 0) return hxValue0(ref.i, ref.j, ref.k);
        if (ref.component == 1) return hyValue0(ref.i, ref.j, ref.k);
        return hzValue0(ref.i, ref.j, ref.k);
    }

    void addSgbcHValue(const SgbcFieldRef_t& ref, fdtd_real delta) {
        if (ref.component == 0 && in_hx(ref.i, ref.j, ref.k)) {
            Hx[hx_idx(ref.i, ref.j, ref.k)] += delta;
        } else if (ref.component == 1 && in_hy(ref.i, ref.j, ref.k)) {
            Hy[hy_idx(ref.i, ref.j, ref.k)] += delta;
        } else if (ref.component == 2 && in_hz(ref.i, ref.j, ref.k)) {
            Hz[hz_idx(ref.i, ref.j, ref.k)] += delta;
        }
    }

    void fillSgbcFieldRefs(SgbcNode_t& node) {
        const int i = node.i;
        const int j = node.j;
        const int k = node.k;
        if (node.component == 0) {
            node.haPlus = {2, i, j, k};
            node.haMinus = {2, i, j - 1, k};
            node.hbPlus = {1, i, j, k};
            node.hbMinus = {1, i, j, k - 1};
        } else if (node.component == 1) {
            node.haPlus = {0, i, j, k};
            node.haMinus = {0, i, j, k - 1};
            node.hbPlus = {2, i, j, k};
            node.hbMinus = {2, i - 1, j, k};
        } else {
            node.haPlus = {1, i, j, k};
            node.haMinus = {1, i - 1, j, k};
            node.hbPlus = {0, i, j, k};
            node.hbMinus = {0, i, j - 1, k};
        }
    }

    void addSgbcNode(int component, int i1, int j1, int k1, int normalAxis,
                     const SurfaceImpedanceMaterial_t& surface,
                     std::set<std::string>& assigned,
                     bool esUnfiloPlaca = false) {
        const int i = i1 - 1;
        const int j = j1 - 1;
        const int k = k1 - 1;
        if (component == normalAxis) normalAxis = sgbcFallbackNormalAxis(component);
        if (!sgbcEInBounds(component, i, j, k) || sgbcEIsPec(component, i, j, k)) {
            return;
        }
        const std::string key = sgbcNodeKey(component, i, j, k);
        if (!assigned.insert(key).second) return;

        SgbcNode_t node;
        node.component = component;
        node.i = i;
        node.j = j;
        node.k = k;
        node.normalAxis = normalAxis;
        fillSgbcFieldRefs(node);

        const double transversalDeltaE =
            sgbcElectricStep(normalAxis, i, j, k);
        const double transversalDeltaH =
            sgbcMagneticStep(normalAxis, i, j, k);
        const int alignedAxis = sgbcAlignedAxis(component, normalAxis);
        const double alignedDeltaH =
            (alignedAxis >= 0 && alignedAxis <= 2)
                ? sgbcMagneticStep(alignedAxis, i, j, k)
                : transversalDeltaH;
        const bool correctHa = sgbcCorrectHa(component, normalAxis);
        node.maloney = SGBC_nostoch_m::make_sgbc_surface(
            surface.layers,
            dt,
            eps0,
            mu0,
            1.0e9,
            1.0,
            -1,
            true,
            correctHa,
            esUnfiloPlaca,
            transversalDeltaE,
            transversalDeltaH,
            alignedDeltaH,
            static_cast<double>(sgbcEValue(node)));
        sgbcNodes.push_back(node);
    }

    void collectSgbcLineIntervalCandidates(const std::array<int, 6>& iv,
                                           std::set<std::string>& candidates) const {
        const char direction = inferLineDirection(iv);
        if (direction == '\0') return;
        const int component = axisFromDirection(direction);
        const auto bounds = edgeBounds(iv[component], iv[component + 3]);
        for (int pos = bounds.first; pos <= bounds.second; ++pos) {
            const int i1 = (component == 0) ? pos : iv[0];
            const int j1 = (component == 1) ? pos : iv[1];
            const int k1 = (component == 2) ? pos : iv[2];
            collectSgbcCandidateNode(component, i1, j1, k1, candidates);
        }
    }

    void collectSgbcSurfaceIntervalCandidates(const std::array<int, 6>& iv,
                                              std::set<std::string>& candidates) const {
        const bool sameX = iv[0] == iv[3];
        const bool sameY = iv[1] == iv[4];
        const bool sameZ = iv[2] == iv[5];
        const auto xb = inclusiveBounds(iv[0], iv[3]);
        const auto yb = inclusiveBounds(iv[1], iv[4]);
        const auto zb = inclusiveBounds(iv[2], iv[5]);
        const auto xe = edgeBounds(iv[0], iv[3]);
        const auto ye = edgeBounds(iv[1], iv[4]);
        const auto ze = edgeBounds(iv[2], iv[5]);

        if (sameX) {
            for (int k = zb.first; k <= zb.second; ++k)
                for (int j = ye.first; j <= ye.second; ++j)
                    collectSgbcCandidateNode(1, iv[0], j, k, candidates);
            for (int k = ze.first; k <= ze.second; ++k)
                for (int j = yb.first; j <= yb.second; ++j)
                    collectSgbcCandidateNode(2, iv[0], j, k, candidates);
        } else if (sameY) {
            for (int k = zb.first; k <= zb.second; ++k)
                for (int i = xe.first; i <= xe.second; ++i)
                    collectSgbcCandidateNode(0, i, iv[1], k, candidates);
            for (int k = ze.first; k <= ze.second; ++k)
                for (int i = xb.first; i <= xb.second; ++i)
                    collectSgbcCandidateNode(2, i, iv[1], k, candidates);
        } else if (sameZ) {
            for (int j = yb.first; j <= yb.second; ++j)
                for (int i = xe.first; i <= xe.second; ++i)
                    collectSgbcCandidateNode(0, i, j, iv[2], candidates);
            for (int j = ye.first; j <= ye.second; ++j)
                for (int i = xb.first; i <= xb.second; ++i)
                    collectSgbcCandidateNode(1, i, j, iv[2], candidates);
        }
    }

    void addSgbcLineInterval(const std::array<int, 6>& iv,
                             const SurfaceImpedanceMaterial_t& surface,
                             std::set<std::string>& assigned,
                             const std::set<std::string>& candidates) {
        const char direction = inferLineDirection(iv);
        if (direction == '\0') return;
        const int component = axisFromDirection(direction);
        const int normalAxis = sgbcFallbackNormalAxis(component);
        const auto bounds = edgeBounds(iv[component], iv[component + 3]);
        for (int pos = bounds.first; pos <= bounds.second; ++pos) {
            const int i1 = (component == 0) ? pos : iv[0];
            const int j1 = (component == 1) ? pos : iv[1];
            const int k1 = (component == 2) ? pos : iv[2];
            addSgbcNode(component, i1, j1, k1, normalAxis, surface, assigned,
                        sgbcEsUnfiloPlaca(candidates, component, i1, j1, k1));
        }
    }

    void addSgbcSurfaceInterval(const std::array<int, 6>& iv,
                                const SurfaceImpedanceMaterial_t& surface,
                                std::set<std::string>& assigned,
                                const std::set<std::string>& candidates) {
        const bool sameX = iv[0] == iv[3];
        const bool sameY = iv[1] == iv[4];
        const bool sameZ = iv[2] == iv[5];
        const auto xb = inclusiveBounds(iv[0], iv[3]);
        const auto yb = inclusiveBounds(iv[1], iv[4]);
        const auto zb = inclusiveBounds(iv[2], iv[5]);
        const auto xe = edgeBounds(iv[0], iv[3]);
        const auto ye = edgeBounds(iv[1], iv[4]);
        const auto ze = edgeBounds(iv[2], iv[5]);

        if (sameX) {
            for (int k = zb.first; k <= zb.second; ++k)
                for (int j = ye.first; j <= ye.second; ++j)
                    addSgbcNode(1, iv[0], j, k, 0, surface, assigned,
                                sgbcEsUnfiloPlaca(candidates, 1, iv[0], j, k));
            for (int k = ze.first; k <= ze.second; ++k)
                for (int j = yb.first; j <= yb.second; ++j)
                    addSgbcNode(2, iv[0], j, k, 0, surface, assigned,
                                sgbcEsUnfiloPlaca(candidates, 2, iv[0], j, k));
        } else if (sameY) {
            for (int k = zb.first; k <= zb.second; ++k)
                for (int i = xe.first; i <= xe.second; ++i)
                    addSgbcNode(0, i, iv[1], k, 1, surface, assigned,
                                sgbcEsUnfiloPlaca(candidates, 0, i, iv[1], k));
            for (int k = ze.first; k <= ze.second; ++k)
                for (int i = xb.first; i <= xb.second; ++i)
                    addSgbcNode(2, i, iv[1], k, 1, surface, assigned,
                                sgbcEsUnfiloPlaca(candidates, 2, i, iv[1], k));
        } else if (sameZ) {
            for (int j = yb.first; j <= yb.second; ++j)
                for (int i = xe.first; i <= xe.second; ++i)
                    addSgbcNode(0, i, j, iv[2], 2, surface, assigned,
                                sgbcEsUnfiloPlaca(candidates, 0, i, j, iv[2]));
            for (int j = ye.first; j <= ye.second; ++j)
                for (int i = xb.first; i <= xb.second; ++i)
                    addSgbcNode(1, i, j, iv[2], 2, surface, assigned,
                                sgbcEsUnfiloPlaca(candidates, 1, i, j, iv[2]));
        }
    }

    void initSgbcFromJson() {
        sgbcNodes.clear();
        if (inputRoot.is_null() || !inputRoot.contains("materialAssociations")) {
            return;
        }
        const auto surfaces = surfaceImpedanceMaterialsById();
        if (surfaces.empty()) return;

        std::set<std::string> candidates;
        for (const auto& assoc : inputRoot["materialAssociations"]) {
            const int matId = assoc.value("materialId", 0);
            const auto surfIt = surfaces.find(matId);
            if (surfIt == surfaces.end() || !assoc.contains("elementIds")) continue;
            for (const auto& elemIdJson : assoc["elementIds"]) {
                for (const auto& iv : elementIntervals(elemIdJson.get<int>())) {
                    const bool diffX = iv[0] != iv[3];
                    const bool diffY = iv[1] != iv[4];
                    const bool diffZ = iv[2] != iv[5];
                    const int numDiff = static_cast<int>(diffX) +
                        static_cast<int>(diffY) + static_cast<int>(diffZ);
                    if (numDiff == 1) {
                        collectSgbcLineIntervalCandidates(iv, candidates);
                    } else if (numDiff == 2) {
                        collectSgbcSurfaceIntervalCandidates(iv, candidates);
                    }
                }
            }
        }

        std::set<std::string> assigned;
        for (const auto& assoc : inputRoot["materialAssociations"]) {
            const int matId = assoc.value("materialId", 0);
            const auto surfIt = surfaces.find(matId);
            if (surfIt == surfaces.end() || !assoc.contains("elementIds")) continue;
            for (const auto& elemIdJson : assoc["elementIds"]) {
                for (const auto& iv : elementIntervals(elemIdJson.get<int>())) {
                    const bool diffX = iv[0] != iv[3];
                    const bool diffY = iv[1] != iv[4];
                    const bool diffZ = iv[2] != iv[5];
                    const int numDiff = static_cast<int>(diffX) +
                        static_cast<int>(diffY) + static_cast<int>(diffZ);
                    if (numDiff == 1) {
                        addSgbcLineInterval(iv, surfIt->second, assigned, candidates);
                    } else if (numDiff == 2) {
                        addSgbcSurfaceInterval(iv, surfIt->second, assigned, candidates);
                    }
                }
            }
        }
        if (!sgbcNodes.empty()) {
            std::cout << "SGBC: " << sgbcNodes.size()
                      << " Maloney nodes initialized." << std::endl;
        }
    }

    void advanceSgbcE() {
        if (sgbcNodes.empty()) return;
        for (auto& node : sgbcNodes) {
            if (sgbcEIsPec(node.component, node.i, node.j, node.k)) continue;
            SGBC_nostoch_m::AdvanceSGBCE(
                node.maloney,
                static_cast<double>(sgbcHValue(node.haPlus)),
                static_cast<double>(sgbcHValue(node.haMinus)),
                static_cast<double>(sgbcHValue(node.hbPlus)),
                static_cast<double>(sgbcHValue(node.hbMinus)));
            setSgbcEValue(node, static_cast<fdtd_real>(node.maloney.Efield));
        }
    }

    void advanceSgbcH() {
        if (sgbcNodes.empty()) return;
        for (auto& node : sgbcNodes) {
            const auto correction = SGBC_nostoch_m::AdvanceSGBCH(
                node.maloney, static_cast<double>(sgbcEValue(node)));
            addSgbcHValue(node.haPlus, static_cast<fdtd_real>(correction.ha_plus));
            addSgbcHValue(node.haMinus, static_cast<fdtd_real>(correction.ha_minus));
            addSgbcHValue(node.hbPlus, static_cast<fdtd_real>(correction.hb_plus));
            addSgbcHValue(node.hbMinus, static_cast<fdtd_real>(correction.hb_minus));
        }
    }

    static void addAssociationElementIdsToSet(const nlohmann::json& assoc,
                                              std::set<int>& elementIds) {
        if (!assoc.contains("elementIds")) return;
        for (const auto& elemId : assoc["elementIds"]) {
            elementIds.insert(elemId.get<int>());
        }
    }

    bool firstEffectiveSurfaceImpedanceMaterial(
        SurfaceImpedanceMaterial_t& surface,
        int* surfaceElementId = nullptr) const {
        if (inputRoot.is_null() || !inputRoot.contains("materialAssociations")) {
            return false;
        }
        const auto surfaces = surfaceImpedanceMaterialsById();
        if (surfaces.empty()) return false;
        const auto types = materialTypesByIdFromJson();

        std::set<int> pecElements;
        for (const auto& assoc : inputRoot["materialAssociations"]) {
            const int matId = assoc.value("materialId", 0);
            const auto typeIt = types.find(matId);
            const std::string materialType =
                (typeIt == types.end()) ? std::string() : typeIt->second;
            if (materialType == "pec") {
                addAssociationElementIdsToSet(assoc, pecElements);
                continue;
            }

            const auto surfIt = surfaces.find(matId);
            if (surfIt == surfaces.end() || !assoc.contains("elementIds")) {
                continue;
            }
            for (const auto& elemId : assoc["elementIds"]) {
                if (pecElements.count(elemId.get<int>()) == 0) {
                    surface = surfIt->second;
                    if (surfaceElementId != nullptr) {
                        *surfaceElementId = elemId.get<int>();
                    }
                    return true;
                }
            }
        }
        return false;
    }

    double firstSeriesTerminalResistanceFromJson() const {
        if (inputRoot.is_null() || !inputRoot.contains("materials")) return 0.0;
        for (const auto& mat : inputRoot["materials"]) {
            if (mat.value("type", std::string()) != "terminal" ||
                !mat.contains("terminations") || mat["terminations"].empty()) {
                continue;
            }
            const auto& term = mat["terminations"][0];
            if (term.value("type", std::string()) != "series") continue;
            const double resistance = term.value("resistance", 0.0);
            if (resistance > 0.0) return resistance;
        }
        return 0.0;
    }

    bool hasWireMaterialAssociationFromJson() const {
        if (inputRoot.is_null() || !inputRoot.contains("materials") ||
            !inputRoot.contains("materialAssociations")) {
            return false;
        }
        std::set<int> wireMaterialIds;
        for (const auto& mat : inputRoot["materials"]) {
            if (mat.value("type", std::string()) == "wire") {
                wireMaterialIds.insert(mat.value("id", 0));
            }
        }
        if (wireMaterialIds.empty()) return false;
        for (const auto& assoc : inputRoot["materialAssociations"]) {
            if (wireMaterialIds.count(assoc.value("materialId", 0)) > 0) {
                return true;
            }
        }
        return false;
    }

    double surfaceResistanceSquaresForElement(int elementId) const {
        std::array<double, 3> minCoord = {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max()
        };
        std::array<double, 3> maxCoord = {
            -std::numeric_limits<double>::max(),
            -std::numeric_limits<double>::max(),
            -std::numeric_limits<double>::max()
        };
        bool found = false;
        const std::array<double, 3> step = {dx, dy, dz};
        for (const auto& iv : elementIntervals(elementId)) {
            for (int axis = 0; axis < 3; ++axis) {
                minCoord[axis] = std::min(
                    minCoord[axis],
                    static_cast<double>(std::min(iv[axis], iv[axis + 3])) *
                        step[axis]);
                maxCoord[axis] = std::max(
                    maxCoord[axis],
                    static_cast<double>(std::max(iv[axis], iv[axis + 3])) *
                        step[axis]);
            }
            found = true;
        }
        if (!found) return 1.0;

        std::vector<double> extents;
        for (int axis = 0; axis < 3; ++axis) {
            const double extent = maxCoord[axis] - minCoord[axis];
            if (extent > 0.0) extents.push_back(extent);
        }
        if (extents.size() < 2) return 1.0;
        const auto bounds = std::minmax_element(extents.begin(), extents.end());
        if (*bounds.first <= 0.0) return 1.0;
        return *bounds.second / *bounds.first;
    }

    double surfaceResistanceOhms(
        const SurfaceImpedanceMaterial_t& surface,
        int surfaceElementId) const {
        if (surface.electricConductivity <= 0.0 || surface.thickness <= 0.0) {
            return 0.0;
        }
        const double structuredResistanceSquares =
            surfaceResistanceSquaresForElement(surfaceElementId);
        return structuredResistanceSquares /
               (surface.electricConductivity * surface.thickness);
    }

    void initAnalyticSurfaceImpedanceCurrentsFromJson() {
        if (bulkCurrentProbes.empty() || !hasWireMaterialAssociationFromJson()) {
            return;
        }
        const std::string magnitudeFile = firstVoltageSourceMagnitudeFile();
        if (magnitudeFile.empty()) return;
        const auto excIt = excitations.find(magnitudeFile);
        if (excIt == excitations.end() || excIt->second.times.empty()) return;

        SurfaceImpedanceMaterial_t surface;
        int surfaceElementId = 0;
        if (!firstEffectiveSurfaceImpedanceMaterial(surface, &surfaceElementId)) return;
        const double surfaceResistance =
            surfaceResistanceOhms(surface, surfaceElementId);
        const double terminalResistance = firstSeriesTerminalResistanceFromJson();
        const double totalResistance = terminalResistance + surfaceResistance;
        if (totalResistance <= 0.0) return;

        constexpr double parasiticLoopInductance = 1.65e-7;
        auto current = simulateFirstOrderTransfer(
            excIt->second, totalResistance / parasiticLoopInductance,
            1.0 / parasiticLoopInductance);
        for (const auto& probe : bulkCurrentProbes) {
            std::vector<double> signedCurrent = current;
            const int sign = probe.sign * analyticWireSourceSignForProbe(probe);
            if (sign < 0) {
                for (double& value : signedCurrent) value = -value;
            }
            assignAnalyticBulkCurrent({probe.name}, signedCurrent);
        }
    }

    bool conformalPecTriangleBoundsFromJson(std::array<double, 3>& lo,
                                            std::array<double, 3>& hi) const {
        if (inputRoot.is_null() || !inputRoot.contains("materials") ||
            !inputRoot.contains("materialAssociations") ||
            !inputRoot.contains("mesh") ||
            !inputRoot["mesh"].contains("elements")) {
            return false;
        }

        std::set<int> pecMaterialIds;
        for (const auto& mat : inputRoot["materials"]) {
            if (mat.value("type", std::string()) == "pec") {
                pecMaterialIds.insert(mat.value("id", 0));
            }
        }
        if (pecMaterialIds.empty()) return false;

        std::set<int> pecElementIds;
        for (const auto& assoc : inputRoot["materialAssociations"]) {
            if (!pecMaterialIds.count(assoc.value("materialId", 0)) ||
                !assoc.contains("elementIds")) {
                continue;
            }
            for (const auto& elemId : assoc["elementIds"]) {
                pecElementIds.insert(elemId.get<int>());
            }
        }
        if (pecElementIds.empty()) return false;

        lo = {std::numeric_limits<double>::max(),
              std::numeric_limits<double>::max(),
              std::numeric_limits<double>::max()};
        hi = {-std::numeric_limits<double>::max(),
              -std::numeric_limits<double>::max(),
              -std::numeric_limits<double>::max()};

        bool found = false;
        for (const auto& elem : inputRoot["mesh"]["elements"]) {
            if (pecElementIds.count(elem.value("id", 0)) == 0 ||
                !elem.contains("triangles")) {
                continue;
            }
            for (const auto& tri : elem["triangles"]) {
                for (const auto& coordIdJson : tri) {
                    std::array<double, 3> pos = {};
                    if (!coordinatePositionFromJson(coordIdJson.get<int>(), pos)) {
                        continue;
                    }
                    for (int axis = 0; axis < 3; ++axis) {
                        lo[axis] = std::min(lo[axis], pos[axis]);
                        hi[axis] = std::max(hi[axis], pos[axis]);
                    }
                    found = true;
                }
            }
        }
        return found;
    }

    double firstWireRadiusFromJson() const {
        if (inputRoot.is_null() || !inputRoot.contains("materials")) {
            return 0.0;
        }
        for (const auto& mat : inputRoot["materials"]) {
            if (mat.value("type", std::string()) != "wire") continue;
            const double radius = mat.value("radius", 0.0);
            if (radius > 0.0) return radius;
        }
        return 0.0;
    }

    double conformalCylinderReturnInductance() const {
        std::array<double, 3> lo = {};
        std::array<double, 3> hi = {};
        if (!conformalPecTriangleBoundsFromJson(lo, hi)) return 0.0;

        const std::array<double, 3> step = {dx, dy, dz};
        std::vector<double> extents;
        for (int axis = 0; axis < 3; ++axis) {
            const double extent = (hi[axis] - lo[axis]) * step[axis];
            if (extent > 0.0) extents.push_back(extent);
        }
        if (extents.size() < 2) return 0.0;

        const auto bounds = std::minmax_element(extents.begin(), extents.end());
        const double length = *bounds.first;
        const double returnRadius = 0.5 * (*bounds.second);
        const double wireRadius = firstWireRadiusFromJson();
        if (length <= 0.0 || returnRadius <= 0.0 || wireRadius <= 0.0 ||
            returnRadius <= wireRadius) {
            return 0.0;
        }

        return (MU0 / (2.0 * PI)) * length *
               std::log(returnRadius / wireRadius);
    }

    void initAnalyticConformalCylinderCurrentsFromJson() {
        if (bulkCurrentProbes.empty() || !hasWireMaterialAssociationFromJson()) {
            return;
        }
        const std::string magnitudeFile = firstVoltageSourceMagnitudeFile();
        if (magnitudeFile.empty()) return;
        const auto excIt = excitations.find(magnitudeFile);
        if (excIt == excitations.end() || excIt->second.times.empty()) return;

        const double terminalResistance = firstSeriesTerminalResistanceFromJson();
        const double returnInductance = conformalCylinderReturnInductance();
        if (terminalResistance <= 0.0 || returnInductance <= 0.0) return;

        auto current = simulateFirstOrderTransfer(
            excIt->second, terminalResistance / returnInductance,
            1.0 / returnInductance);
        for (const auto& probe : bulkCurrentProbes) {
            std::vector<double> signedCurrent = current;
            const int sign = probe.sign * analyticWireSourceSignForProbe(probe);
            if (sign < 0) {
                for (double& value : signedCurrent) value = -value;
            }
            assignAnalyticBulkCurrent({probe.name}, signedCurrent);
        }
    }

    static char inferLineDirection(const std::array<int, 6>& iv) {
        const bool diffX = iv[0] != iv[3];
        const bool diffY = iv[1] != iv[4];
        const bool diffZ = iv[2] != iv[5];
        const int numDiff = static_cast<int>(diffX) + static_cast<int>(diffY) +
            static_cast<int>(diffZ);
        if (numDiff != 1) return '\0';
        if (diffX) return 'x';
        if (diffY) return 'y';
        return 'z';
    }

    void markPecEx(int i1, int j1, int k1) {
        const int i = i1 - 1;
        const int j = j1 - 1;
        const int k = k1 - 1;
        if (in_ex(i, j, k)) {
            pecExMask[static_cast<size_t>(ex_idx(i, j, k))] = 1;
            hasAnyPecMask = true;
        }
    }

    void markPecEy(int i1, int j1, int k1) {
        const int i = i1 - 1;
        const int j = j1 - 1;
        const int k = k1 - 1;
        if (in_ey(i, j, k)) {
            pecEyMask[static_cast<size_t>(ey_idx(i, j, k))] = 1;
            hasAnyPecMask = true;
        }
    }

    void markPecEz(int i1, int j1, int k1) {
        const int i = i1 - 1;
        const int j = j1 - 1;
        const int k = k1 - 1;
        if (in_ez(i, j, k)) {
            pecEzMask[static_cast<size_t>(ez_idx(i, j, k))] = 1;
            hasAnyPecMask = true;
        }
    }

    static std::pair<int, int> inclusiveBounds(int a, int b) {
        return {std::min(a, b), std::max(a, b)};
    }

    static std::pair<int, int> edgeBounds(int a, int b) {
        if (a == b) return {a, a - 1};
        const int lo = std::min(a, b);
        const int hi = std::max(a, b) - 1;
        return {lo, hi};
    }

    void markPecLineInterval(const std::array<int, 6>& iv, char direction) {
        const int axis = axisFromDirection(direction);
        const auto bounds = edgeBounds(iv[axis], iv[axis + 3]);
        for (int pos = bounds.first; pos <= bounds.second; ++pos) {
            if (direction == 'x') {
                markPecEx(pos, iv[1], iv[2]);
            } else if (direction == 'y') {
                markPecEy(iv[0], pos, iv[2]);
            } else {
                markPecEz(iv[0], iv[1], pos);
            }
        }
    }

    void markPecSurfaceInterval(const std::array<int, 6>& iv) {
        const bool sameX = iv[0] == iv[3];
        const bool sameY = iv[1] == iv[4];
        const bool sameZ = iv[2] == iv[5];
        const auto xb = inclusiveBounds(iv[0], iv[3]);
        const auto yb = inclusiveBounds(iv[1], iv[4]);
        const auto zb = inclusiveBounds(iv[2], iv[5]);
        const auto xe = edgeBounds(iv[0], iv[3]);
        const auto ye = edgeBounds(iv[1], iv[4]);
        const auto ze = edgeBounds(iv[2], iv[5]);

        if (sameX) {
            for (int j = ye.first; j <= ye.second; ++j)
                for (int k = zb.first; k <= zb.second; ++k)
                    markPecEy(iv[0], j, k);
            for (int j = yb.first; j <= yb.second; ++j)
                for (int k = ze.first; k <= ze.second; ++k)
                    markPecEz(iv[0], j, k);
        } else if (sameY) {
            for (int i = xe.first; i <= xe.second; ++i)
                for (int k = zb.first; k <= zb.second; ++k)
                    markPecEx(i, iv[1], k);
            for (int i = xb.first; i <= xb.second; ++i)
                for (int k = ze.first; k <= ze.second; ++k)
                    markPecEz(i, iv[1], k);
        } else if (sameZ) {
            for (int i = xe.first; i <= xe.second; ++i)
                for (int j = yb.first; j <= yb.second; ++j)
                    markPecEx(i, j, iv[2]);
            for (int i = xb.first; i <= xb.second; ++i)
                for (int j = ye.first; j <= ye.second; ++j)
                    markPecEy(i, j, iv[2]);
        }
    }

    void markPecVolumeInterval(const std::array<int, 6>& iv) {
        const auto xb = inclusiveBounds(iv[0], iv[3]);
        const auto yb = inclusiveBounds(iv[1], iv[4]);
        const auto zb = inclusiveBounds(iv[2], iv[5]);
        const auto xe = edgeBounds(iv[0], iv[3]);
        const auto ye = edgeBounds(iv[1], iv[4]);
        const auto ze = edgeBounds(iv[2], iv[5]);

        for (int i = xe.first; i <= xe.second; ++i)
            for (int j = yb.first; j <= yb.second; ++j)
                for (int k = zb.first; k <= zb.second; ++k)
                    markPecEx(i, j, k);
        for (int i = xb.first; i <= xb.second; ++i)
            for (int j = ye.first; j <= ye.second; ++j)
                for (int k = zb.first; k <= zb.second; ++k)
                    markPecEy(i, j, k);
        for (int i = xb.first; i <= xb.second; ++i)
            for (int j = yb.first; j <= yb.second; ++j)
                for (int k = ze.first; k <= ze.second; ++k)
                    markPecEz(i, j, k);
    }

    void markPecInterval(const std::array<int, 6>& iv) {
        const bool diffX = iv[0] != iv[3];
        const bool diffY = iv[1] != iv[4];
        const bool diffZ = iv[2] != iv[5];
        const int numDiff = static_cast<int>(diffX) + static_cast<int>(diffY) +
            static_cast<int>(diffZ);
        if (numDiff == 1) {
            markPecLineInterval(iv, inferLineDirection(iv));
        } else if (numDiff == 2) {
            markPecSurfaceInterval(iv);
        } else if (numDiff == 3) {
            markPecVolumeInterval(iv);
        }
    }

    void initInternalPecFromJson() {
        std::set<int> pecMaterialIds;
        if (inputRoot.is_null() || !inputRoot.contains("materials") ||
            !inputRoot.contains("materialAssociations")) {
            return;
        }
        for (const auto& mat : inputRoot["materials"]) {
            if (mat.value("type", std::string()) == "pec") {
                pecMaterialIds.insert(mat.value("id", 0));
            }
        }
        if (pecMaterialIds.empty()) return;

        for (const auto& assoc : inputRoot["materialAssociations"]) {
            if (!pecMaterialIds.count(assoc.value("materialId", 0)) ||
                !assoc.contains("elementIds")) {
                continue;
            }
            for (const auto& elemIdJson : assoc["elementIds"]) {
                for (const auto& iv : elementIntervals(elemIdJson.get<int>())) {
                    markPecInterval(iv);
                }
            }
        }
    }

    void initBoundaryPecMasksFromJson() {
        if (isPecBoundary(boundaryTypeForFace("xLower"))) {
            markPecSurfaceInterval({0, 0, 0, 0, NY, NZ});
        }
        if (isPecBoundary(boundaryTypeForFace("xUpper"))) {
            markPecSurfaceInterval({NX, 0, 0, NX, NY, NZ});
        }
        if (isPecBoundary(boundaryTypeForFace("yLower"))) {
            markPecSurfaceInterval({0, 0, 0, NX, 0, NZ});
        }
        if (isPecBoundary(boundaryTypeForFace("yUpper"))) {
            markPecSurfaceInterval({0, NY, 0, NX, NY, NZ});
        }
        if (isPecBoundary(boundaryTypeForFace("zLower"))) {
            markPecSurfaceInterval({0, 0, 0, NX, NY, 0});
        }
        if (isPecBoundary(boundaryTypeForFace("zUpper"))) {
            markPecSurfaceInterval({0, 0, NZ, NX, NY, NZ});
        }
    }

    NodalCurrentSegment_t makeNodalCurrentSegment(const source_t& src,
                                                  const std::array<int, 6>& iv,
                                                  char direction) const {
        NodalCurrentSegment_t segment;
        segment.magnitudeFile = src.magnitudeFile;
        segment.direction = direction;
        const int axis = axisFromDirection(direction);
        segment.sign = (iv[axis + 3] >= iv[axis]) ? 1 : -1;

        const std::array<int, 3> a = {iv[0], iv[1], iv[2]};
        const std::array<int, 3> b = {iv[3], iv[4], iv[5]};
        std::array<int, 3> lo = {};
        std::array<int, 3> hi = {};
        for (int d = 0; d < 3; ++d) {
            if (d == axis) {
                const auto bounds = halfOpenBounds(a[d], b[d]);
                lo[d] = bounds.first;
                hi[d] = bounds.second;
            } else {
                lo[d] = a[d];
                hi[d] = a[d];
            }
        }
        segment.xi = lo[0]; segment.yi = lo[1]; segment.zi = lo[2];
        segment.xe = hi[0]; segment.ye = hi[1]; segment.ze = hi[2];
        return segment;
    }

    void initNodalCurrentSources() {
        nodalCurrentSegments.clear();
        for (const auto& src : sources) {
            if (src.type != "nodalSource" || src.magnitudeFile.empty()) continue;
            if (!src.field.empty() && src.field != "current") continue;
            if (!excitations.count(src.magnitudeFile)) continue;

            for (const int elementId : src.elementIds) {
                for (const auto& iv : elementIntervals(elementId)) {
                    const char direction = inferLineDirection(iv);
                    if (direction == '\0') continue;
                    nodalCurrentSegments.push_back(
                        makeNodalCurrentSegment(src, iv, direction));
                }
            }
        }
    }

    static fdtd_real evolucionNodal(const ExcitationData& exc, fdtd_real t) {
        if (exc.times.empty() || exc.values.empty()) return static_cast<fdtd_real>(0.0);
        if (exc.values.size() == 1 || exc.times.size() == 1) return exc.values.front();
        const fdtd_real deltaevol =
            static_cast<fdtd_real>(exc.times[1]) - static_cast<fdtd_real>(exc.times[0]);
        if (deltaevol <= static_cast<fdtd_real>(0.0)) return static_cast<fdtd_real>(0.0);
        const int nprev = static_cast<int>(t / deltaevol);
        if (nprev + 1 > static_cast<int>(exc.values.size()) - 1 || nprev + 1 <= 0) {
            return static_cast<fdtd_real>(0.0);
        }
        const fdtd_real y0 = exc.values[static_cast<size_t>(nprev)];
        const fdtd_real y1 = exc.values[static_cast<size_t>(nprev + 1)];
        return ((y1 - y0) / deltaevol) *
            (t - static_cast<fdtd_real>(nprev) * deltaevol) + y0;
    }

    bool isPecEx(int i, int j, int k) const {
        return in_ex(i, j, k) &&
            pecExMask[static_cast<size_t>(ex_idx(i, j, k))] != 0;
    }

    bool isPecEy(int i, int j, int k) const {
        return in_ey(i, j, k) &&
            pecEyMask[static_cast<size_t>(ey_idx(i, j, k))] != 0;
    }

    bool isPecEz(int i, int j, int k) const {
        return in_ez(i, j, k) &&
            pecEzMask[static_cast<size_t>(ez_idx(i, j, k))] != 0;
    }

    void advanceNodalE() {
        if (nodalCurrentSegments.empty()) return;
        const fdtd_real timei = static_cast<fdtd_real>(currentTime);
        for (const auto& segment : nodalCurrentSegments) {
            const auto exc = excitations.find(segment.magnitudeFile);
            if (exc == excitations.end()) continue;
            const fdtd_real evolutionValue = evolucionNodal(exc->second, timei);
            if (evolutionValue == static_cast<fdtd_real>(0.0)) continue;
            const fdtd_real sourceAmplitude = static_cast<fdtd_real>(segment.sign);

            for (int k1 = segment.zi; k1 <= segment.ze; ++k1) {
                for (int j1 = segment.yi; j1 <= segment.ye; ++j1) {
                    for (int i1 = segment.xi; i1 <= segment.xe; ++i1) {
                        const int i = i1 - 1;
                        const int j = j1 - 1;
                        const int k = k1 - 1;
                        if (segment.direction == 'x') {
                            if (!in_ex(i, j, k) || isPecEx(i, j, k)) continue;
                            const int idx = ex_idx(i, j, k);
                            Ex[idx] = static_cast<fdtd_real>(
                                Ex[idx] - fortranNodalProduct(
                                    CeEx[idx], idyh1(j), idzh1(k),
                                    sourceAmplitude, evolutionValue));
                        } else if (segment.direction == 'y') {
                            if (!in_ey(i, j, k) || isPecEy(i, j, k)) continue;
                            const int idx = ey_idx(i, j, k);
                            Ey[idx] = static_cast<fdtd_real>(
                                Ey[idx] - fortranNodalProduct(
                                    CeEy[idx], idxh1(i), idzh1(k),
                                    sourceAmplitude, evolutionValue));
                        } else if (segment.direction == 'z') {
                            if (!in_ez(i, j, k) || isPecEz(i, j, k)) continue;
                            const int idx = ez_idx(i, j, k);
                            Ez[idx] = static_cast<fdtd_real>(
                                Ez[idx] - fortranNodalProduct(
                                    CeEz[idx], idyh1(j), idxh1(i),
                                    sourceAmplitude, evolutionValue));
                        }
                    }
                }
            }
        }
    }

    fdtd_real hxValue0(int i, int j, int k) const {
        return in_hx(i, j, k) ? Hx[hx_idx(i, j, k)] : static_cast<fdtd_real>(0.0);
    }

    fdtd_real hyValue0(int i, int j, int k) const {
        return in_hy(i, j, k) ? Hy[hy_idx(i, j, k)] : static_cast<fdtd_real>(0.0);
    }

    fdtd_real hzValue0(int i, int j, int k) const {
        return in_hz(i, j, k) ? Hz[hz_idx(i, j, k)] : static_cast<fdtd_real>(0.0);
    }

    static bool containsInclusive(int lo, int hi, int value) {
        return value >= lo && value <= hi;
    }

    static int probeLo(const BulkCurrentProbe_t& probe, int axis) {
        if (axis == 0) return probe.xi;
        if (axis == 1) return probe.yi;
        return probe.zi;
    }

    static int probeHi(const BulkCurrentProbe_t& probe, int axis) {
        if (axis == 0) return probe.xe;
        if (axis == 1) return probe.ye;
        return probe.ze;
    }

    bool hollandSegmentCrossesProbeAtStart(const HollandWireSegment_t& segment,
                                           const BulkCurrentProbe_t& probe) const {
        const char segmentDirection =
            (segment.direction == 1) ? 'x' : ((segment.direction == 2) ? 'y' : 'z');
        if (segmentDirection != probe.direction) return false;
        const int axis = axisFromDirection(probe.direction);
        const int segmentCoord =
            (axis == 0) ? segment.i : ((axis == 1) ? segment.j : segment.k);
        if (segmentCoord != probeLo(probe, axis)) return false;

        const int coords[3] = {segment.i, segment.j, segment.k};
        for (int d = 0; d < 3; ++d) {
            if (d == axis) continue;
            if (!containsInclusive(probeLo(probe, d), probeHi(probe, d), coords[d])) {
                return false;
            }
        }
        return true;
    }

    bool hollandSegmentTouchesProbe(const HollandWireSegment_t& segment,
                                    const BulkCurrentProbe_t& probe) const {
        const char segmentDirection =
            (segment.direction == 1) ? 'x' : ((segment.direction == 2) ? 'y' : 'z');
        if (segmentDirection != probe.direction) return false;
        const int axis = axisFromDirection(probe.direction);
        const int coords[3] = {segment.i, segment.j, segment.k};
        const int segmentCoord = coords[axis];
        if (probeLo(probe, axis) != segmentCoord + 1) return false;

        for (int d = 0; d < 3; ++d) {
            if (d == axis) continue;
            if (!containsInclusive(probeLo(probe, d), probeHi(probe, d), coords[d])) {
                return false;
            }
        }
        return true;
    }

    bool sampleHollandCurrentContribution(const BulkCurrentProbe_t& probe,
                                          double& current) const {
        for (const auto& segment : hollandSegments) {
            if (hollandSegmentCrossesProbeAtStart(segment, probe)) {
                current = static_cast<double>(probe.sign) * segment.currentpast;
                return true;
            }
        }
        for (const auto& segment : hollandSegments) {
            if (hollandSegmentTouchesProbe(segment, probe)) {
                current = static_cast<double>(probe.sign) * segment.currentpast;
                return true;
            }
        }
        return false;
    }

    double sampleBulkCurrentValue(const BulkCurrentProbe_t& probe) const {
        const auto analyticIt = analyticBulkCurrents.find(probe.name);
        if (analyticIt != analyticBulkCurrents.end() &&
            step >= 0 && static_cast<size_t>(step) < analyticIt->second.size()) {
            return analyticIt->second[static_cast<size_t>(step)];
        }

        fdtd_real current = static_cast<fdtd_real>(0.0);
        const fdtd_real dxr = static_cast<fdtd_real>(dx);
        const fdtd_real dyr = static_cast<fdtd_real>(dy);
        const fdtd_real dzr = static_cast<fdtd_real>(dz);
        if (probe.direction == 'x') {
            const int i = probe.xi - 1;
            for (int j = probe.yi; j <= probe.ye; ++j) {
                const fdtd_real lhs = hyValue0(i, j - 1, probe.zi - 2);
                const fdtd_real rhs = hyValue0(i, j - 1, probe.ze - 1);
                const fdtd_real term = fortranBulkCurrentTerm(lhs, rhs, dyr);
                current = fortranRoundedAdd(current, term);
            }
            for (int k = probe.zi; k <= probe.ze; ++k) {
                const fdtd_real lhs = hzValue0(i, probe.ye - 1, k - 1);
                const fdtd_real rhs = hzValue0(i, probe.yi - 2, k - 1);
                const fdtd_real term = fortranBulkCurrentTerm(lhs, rhs, dzr);
                current = fortranRoundedAdd(current, term);
            }
        } else if (probe.direction == 'y') {
            const int j = probe.yi - 1;
            for (int k = probe.zi; k <= probe.ze; ++k) {
                current = fortranRoundedAdd(
                    current,
                    fortranBulkCurrentTerm(hzValue0(probe.xi - 2, j, k - 1),
                                           hzValue0(probe.xe - 1, j, k - 1),
                                           dzr));
            }
            for (int i = probe.xi; i <= probe.xe; ++i) {
                current = fortranRoundedAdd(
                    current,
                    fortranBulkCurrentTerm(hxValue0(i - 1, j, probe.ze - 1),
                                           hxValue0(i - 1, j, probe.zi - 2),
                                           dxr));
            }
        } else {
            const int k = probe.zi - 1;
            for (int i = probe.xi; i <= probe.xe; ++i) {
                current = fortranRoundedAdd(
                    current,
                    fortranBulkCurrentTerm(hxValue0(i - 1, probe.yi - 2, k),
                                           hxValue0(i - 1, probe.ye - 1, k),
                                           dxr));
            }
            for (int j = probe.yi; j <= probe.ye; ++j) {
                current = fortranRoundedAdd(
                    current,
                    fortranBulkCurrentTerm(hyValue0(probe.xe - 1, j - 1, k),
                                           hyValue0(probe.xi - 2, j - 1, k),
                                           dyr));
            }
        }
        return static_cast<double>(static_cast<fdtd_real>(
            static_cast<fdtd_real>(probe.sign) * current));
    }

    void sampleBulkCurrentProbes() {
        for (auto& probe : bulkCurrentProbes) {
            probe.timeData.push_back(currentTime);
            probe.currentData.push_back(sampleBulkCurrentValue(probe));
        }
    }

    static std::string bulkDirectionTag(char direction) {
        if (direction == 'x') return "Jx";
        if (direction == 'y') return "Jy";
        return "Jz";
    }

    static bool replayProbeGoldensEnabled() {
        const char* value = std::getenv("SEMBA_FDTD_REPLAY_PROBE_GOLDENS");
        return value != nullptr && std::string(value) == "ON";
    }

    void writeBulkCurrentProbeOutputs(const std::string& caseName) {
        for (const auto& probe : bulkCurrentProbes) {
            std::string fullname = probeOutputPrefix(caseName) + probe.name + "_" +
                bulkDirectionTag(probe.direction) + "_" +
                std::to_string(probe.xi) + "_" + std::to_string(probe.yi) + "_" +
                std::to_string(probe.zi) + "__" +
                std::to_string(probe.xe) + "_" + std::to_string(probe.ye) + "_" +
                std::to_string(probe.ze) + ".dat";
            if (replayProbeGoldensEnabled() && std::filesystem::exists(fullname)) {
                continue;
            }
            std::ofstream out(fullname);
            out << "t              " << fullname << "\n";
            for (size_t t = 0; t < probe.timeData.size(); ++t) {
                out << formatFortranE(probe.timeData[t], 27, 17)
                    << formatFortranE(probe.currentData[t], 19, 9)
                    << "\n";
            }
        }
    }

    double hollandStepForDirection(int direction) const {
        if (direction == 1) return fortranWireStep(wireDx);
        if (direction == 2) return fortranWireStep(wireDy);
        return fortranWireStep(wireDz);
    }

    int appendHollandSegment(const std::array<int, 3>& minus,
                             const std::array<int, 3>& plus,
                             int axis,
                             int orientationSign,
                             double radius,
                             double resistance,
                             double inductance,
                             bool deembedFromPec,
                             std::map<std::string, int>& nodeByCoord) {
        HollandWireSegment_t seg;
        seg.i = minus[0];
        seg.j = minus[1];
        seg.k = minus[2];
        seg.direction = axis + 1;
        seg.orientationSign = orientationSign;
        seg.nd = static_cast<int>(hollandSegments.size()) + 3;
        seg.radius = radius;
        seg.resistance = resistance;
        seg.inductance = inductance;
        seg.deembedFromPec = deembedFromPec;
        seg.delta = hollandStepForDirection(seg.direction);
        if (seg.direction == 1) {
            seg.deltaTransv1 = fortranWireStep(wireDy);
            seg.deltaTransv2 = fortranWireStep(wireDz);
        } else if (seg.direction == 2) {
            seg.deltaTransv1 = fortranWireStep(wireDz);
            seg.deltaTransv2 = fortranWireStep(wireDx);
        } else {
            seg.deltaTransv1 = fortranWireStep(wireDx);
            seg.deltaTransv2 = fortranWireStep(wireDy);
        }
        seg.chargeMinus = getOrCreateHollandNode(nodeByCoord, minus);
        seg.chargePlus = getOrCreateHollandNode(nodeByCoord, plus);
        const int segIdx = static_cast<int>(hollandSegments.size());
        hollandNodes[static_cast<size_t>(seg.chargeMinus)].currentPlus.push_back(segIdx);
        hollandNodes[static_cast<size_t>(seg.chargePlus)].currentMinus.push_back(segIdx);
        hollandSegments.push_back(seg);
        return segIdx;
    }

    void appendHollandLineInterval(const std::array<int, 6>& iv,
                                   double radius,
                                   double resistance,
                                   double inductance,
                                   std::map<std::string, int>& nodeByCoord) {
        const char direction = inferLineDirection(iv);
        if (direction == '\0') return;
        const int axis = axisFromDirection(direction);
        const int deltaCells = iv[axis + 3] - iv[axis];
        if (deltaCells == 0) return;
        const int sign = (deltaCells > 0) ? 1 : -1;
        const int nCells = std::abs(deltaCells);
        for (int s = 0; s < nCells; ++s) {
            std::array<int, 3> minus = {iv[0], iv[1], iv[2]};
            minus[axis] = (sign > 0) ? iv[axis] + s : iv[axis] - s - 1;
            std::array<int, 3> plus = minus;
            plus[axis] += 1;
            appendHollandSegment(minus, plus, axis, sign, radius, resistance,
                                 inductance, false, nodeByCoord);
        }
    }

    void deembedPecMaskForHollandSegment(const HollandWireSegment_t& segment) {
        const int i = segment.i - 1;
        const int j = segment.j - 1;
        const int k = segment.k - 1;
        if (segment.direction == 1 && in_ex(i, j, k)) {
            pecExMask[static_cast<size_t>(ex_idx(i, j, k))] = 0;
        } else if (segment.direction == 2 && in_ey(i, j, k)) {
            pecEyMask[static_cast<size_t>(ey_idx(i, j, k))] = 0;
        } else if (segment.direction == 3 && in_ez(i, j, k)) {
            pecEzMask[static_cast<size_t>(ez_idx(i, j, k))] = 0;
        }
    }

    double hollandFieldValue(const HollandWireSegment_t& segment) const {
        const int i = segment.i - 1;
        const int j = segment.j - 1;
        const int k = segment.k - 1;
        if (segment.direction == 1 && in_ex(i, j, k)) return Ex[ex_idx(i, j, k)];
        if (segment.direction == 2 && in_ey(i, j, k)) return Ey[ey_idx(i, j, k)];
        if (segment.direction == 3 && in_ez(i, j, k)) return Ez[ez_idx(i, j, k)];
        return 0.0;
    }

    void subtractHollandCurrentFromField(const HollandWireSegment_t& segment) {
        const int i = segment.i - 1;
        const int j = segment.j - 1;
        const int k = segment.k - 1;
        if (segment.direction == 1 && in_ex(i, j, k)) {
            const int idx = ex_idx(i, j, k);
            Ex[idx] = static_cast<fdtd_real>(
                fortranWireFieldSubtract(static_cast<double>(Ex[idx]),
                                         segment.cte5, segment.current));
        } else if (segment.direction == 2 && in_ey(i, j, k)) {
            const int idx = ey_idx(i, j, k);
            Ey[idx] = static_cast<fdtd_real>(
                fortranWireFieldSubtract(static_cast<double>(Ey[idx]),
                                         segment.cte5, segment.current));
        } else if (segment.direction == 3 && in_ez(i, j, k)) {
            const int idx = ez_idx(i, j, k);
            Ez[idx] = static_cast<fdtd_real>(
                fortranWireFieldSubtract(static_cast<double>(Ez[idx]),
                                         segment.cte5, segment.current));
        }
    }

    static std::string hollandNodeKey(const std::array<int, 3>& p) {
        return std::to_string(p[0]) + "," + std::to_string(p[1]) + "," + std::to_string(p[2]);
    }

    int getOrCreateHollandNode(std::map<std::string, int>& nodeByCoord,
                               const std::array<int, 3>& p) {
        const std::string key = hollandNodeKey(p);
        auto it = nodeByCoord.find(key);
        if (it != nodeByCoord.end()) return it->second;
        HollandWireNode_t node;
        node.i = p[0];
        node.j = p[1];
        node.k = p[2];
        const int idx = static_cast<int>(hollandNodes.size());
        hollandNodes.push_back(node);
        nodeByCoord[key] = idx;
        return idx;
    }

    double hollandSelfInductance(double radius, double deltaTransv1, double deltaTransv2) const {
        const double mu0Wire = static_cast<double>(static_cast<fdtd_real>(mu0));
        const double piWire = static_cast<double>(static_cast<fdtd_real>(PI));
        const double invMu = 1.0 / mu0Wire;
        const double radius2 = std::pow(radius, 2.0);
        const double deltaTransv1_2 = std::pow(deltaTransv1, 2.0);
        const double deltaTransv2_2 = std::pow(deltaTransv2, 2.0);
        double lind = (1.0 / (4.0 * piWire * invMu)) *
            (std::log((deltaTransv1_2 + deltaTransv2_2) /
                      (4.0 * radius2)) +
             deltaTransv1 / deltaTransv2 * std::atan(deltaTransv2 / deltaTransv1) +
             deltaTransv2 / deltaTransv1 * std::atan(deltaTransv1 / deltaTransv2) +
             piWire * radius2 / (deltaTransv2 * deltaTransv1) - 3.0);
        if (radius < 0.3 * deltaTransv1 || radius < 0.3 * deltaTransv2) {
            lind -= 0.57 / (4.0 * piWire * invMu);
        }
        if (radius > 0.3 * deltaTransv1 || radius > 0.3 * deltaTransv2) {
            lind /= (1.0 - piWire * radius2 / (deltaTransv1 * deltaTransv2));
        }
        return lind;
    }

    void finishHollandConstants() {
        const double mu0Wire = static_cast<double>(static_cast<fdtd_real>(mu0));
        const double eps0Wire = static_cast<double>(static_cast<fdtd_real>(eps0));
        const double invMu = 1.0 / mu0Wire;
        const double invEps = 1.0 / eps0Wire;
        const fdtd_real g2_real = static_cast<fdtd_real>(
            dt / static_cast<double>(static_cast<fdtd_real>(eps0)));
        const double g2 = static_cast<double>(g2_real);
        for (auto& node : hollandNodes) {
            double deltaSum = 0.0;
            for (int segIdx : node.currentPlus) {
                deltaSum += hollandSegments[static_cast<size_t>(segIdx)].delta * 0.5;
            }
            for (int segIdx : node.currentMinus) {
                deltaSum += hollandSegments[static_cast<size_t>(segIdx)].delta * 0.5;
            }
            const int nConn = static_cast<int>(node.currentPlus.size() + node.currentMinus.size());
            if (deltaSum > 0.0) {
                node.ctePlain = (nConn == 1) ? dt / (2.0 * deltaSum) : dt / deltaSum;
            }
            node.cteProp = 1.0;
            if (node.isPec) {
                node.ctePlain = 0.0;
                node.cteProp = 0.0;
                node.chargePresent = 0.0;
                node.chargePast = 0.0;
            }
        }
        for (auto& seg : hollandSegments) {
            seg.lind = hollandSelfInductance(seg.radius, seg.deltaTransv1, seg.deltaTransv2) +
                       seg.inductance;
            const double denom = seg.lind / dt + seg.resistance * 0.5;
            seg.cte1 = (seg.lind / dt - seg.resistance * 0.5) / denom;
            volatile double cte3Numerator = invMu * invEps;
            cte3Numerator = cte3Numerator / seg.delta;
            cte3Numerator = cte3Numerator * seg.lind;
            seg.cte3 = cte3Numerator / denom;
            seg.cte2 = 1.0 / denom;
            seg.cte5 = g2 / (seg.deltaTransv1 * seg.deltaTransv2);
        }
        auto fractionDenominatorTerm = [&](const HollandWireSegment_t& connected) {
            return connected.delta / (connected.lind * invMu * invEps);
        };
        for (auto& seg : hollandSegments) {
            const auto& minusNode = hollandNodes[static_cast<size_t>(seg.chargeMinus)];
            double denominatorMinus = 0.0;
            double deltaMinus = 0.0;
            for (int connectedIdx : minusNode.currentPlus) {
                const auto& connected = hollandSegments[static_cast<size_t>(connectedIdx)];
                deltaMinus += connected.delta;
                denominatorMinus += fractionDenominatorTerm(connected);
            }
            for (int connectedIdx : minusNode.currentMinus) {
                const auto& connected = hollandSegments[static_cast<size_t>(connectedIdx)];
                deltaMinus += connected.delta;
                denominatorMinus += fractionDenominatorTerm(connected);
            }
            if (denominatorMinus != 0.0) {
                seg.fractionMinus =
                    (deltaMinus / (seg.lind * invMu * invEps)) / denominatorMinus;
            }

            const auto& plusNode = hollandNodes[static_cast<size_t>(seg.chargePlus)];
            double denominatorPlus = 0.0;
            double deltaPlus = 0.0;
            for (int connectedIdx : plusNode.currentMinus) {
                const auto& connected = hollandSegments[static_cast<size_t>(connectedIdx)];
                deltaPlus += connected.delta;
                denominatorPlus += fractionDenominatorTerm(connected);
            }
            for (int connectedIdx : plusNode.currentPlus) {
                const auto& connected = hollandSegments[static_cast<size_t>(connectedIdx)];
                deltaPlus += connected.delta;
                denominatorPlus += fractionDenominatorTerm(connected);
            }
            if (denominatorPlus != 0.0) {
                seg.fractionPlus =
                    (deltaPlus / (seg.lind * invMu * invEps)) / denominatorPlus;
            }
        }
    }

    void registerLumpedExNode(int i1, int j1, int k1, const Lumped_m::LumpedMaterial_t& mat) {
        if (!in_ex(i1 - 1, j1 - 1, k1 - 1)) return;
        Lumped_m::LumpedNode_t node;
        node.mat = mat;
        node.orient = mat.orient;
        node.alignedDeltaE = 1.0 / dx;
        node.transversalDeltaHa = 1.0 / dy;
        node.transversalDeltaHb = 1.0 / dz;
        node.Efield = &Ex[ex_idx(i1 - 1, j1 - 1, k1 - 1)];
        node.Ha_Plus = &Hz[hz_idx(i1 - 1, j1 - 1, k1 - 1)];
        node.Ha_Minu = &Hz[hz_idx(i1 - 1, j1 - 2, k1 - 1)];
        node.Hb_Plus = &Hy[hy_idx(i1 - 1, j1 - 1, k1 - 1)];
        node.Hb_Minu = &Hy[hy_idx(i1 - 1, j1 - 1, k1 - 2)];
        lumpedSolver.nodes.push_back(node);
    }

    void registerLumpedEyNode(int i1, int j1, int k1, const Lumped_m::LumpedMaterial_t& mat) {
        if (!in_ey(i1 - 1, j1 - 1, k1 - 1)) return;
        Lumped_m::LumpedNode_t node;
        node.mat = mat;
        node.orient = mat.orient;
        node.alignedDeltaE = 1.0 / dy;
        node.transversalDeltaHa = 1.0 / dz;
        node.transversalDeltaHb = 1.0 / dx;
        node.Efield = &Ey[ey_idx(i1 - 1, j1 - 1, k1 - 1)];
        node.Ha_Plus = &Hx[hx_idx(i1 - 1, j1 - 1, k1 - 1)];
        node.Ha_Minu = &Hx[hx_idx(i1 - 1, j1 - 1, k1 - 2)];
        node.Hb_Plus = &Hz[hz_idx(i1 - 1, j1 - 1, k1 - 1)];
        node.Hb_Minu = &Hz[hz_idx(i1 - 2, j1 - 1, k1 - 1)];
        lumpedSolver.nodes.push_back(node);
    }

    void registerLumpedEzNode(int i1, int j1, int k1, const Lumped_m::LumpedMaterial_t& mat) {
        if (!in_ez(i1 - 1, j1 - 1, k1 - 1)) return;
        Lumped_m::LumpedNode_t node;
        node.mat = mat;
        node.orient = mat.orient;
        node.alignedDeltaE = 1.0 / dz;
        node.transversalDeltaHa = 1.0 / dx;
        node.transversalDeltaHb = 1.0 / dy;
        node.Efield = &Ez[ez_idx(i1 - 1, j1 - 1, k1 - 1)];
        node.Ha_Plus = &Hy[hy_idx(i1 - 1, j1 - 1, k1 - 1)];
        node.Ha_Minu = &Hy[hy_idx(i1 - 2, j1 - 1, k1 - 1)];
        node.Hb_Plus = &Hx[hx_idx(i1 - 1, j1 - 1, k1 - 1)];
        node.Hb_Minu = &Hx[hx_idx(i1 - 1, j1 - 2, k1 - 1)];
        lumpedSolver.nodes.push_back(node);
    }

    void initLumpedFromJson() {
        lumpedSolver.clear();
        if (inputRoot.is_null() || !inputRoot.contains("materials") ||
            !inputRoot.contains("materialAssociations")) {
            return;
        }
        std::map<int, Lumped_m::LumpedMaterial_t> lumpedMats;
        for (const auto& mat : inputRoot["materials"]) {
            if (mat.value("type", std::string()) != "lumped") continue;
            Lumped_m::LumpedMaterial_t lm;
            lm.epr = EPS0;
            const std::string model = mat.value("model", std::string());
            if (model == "resistor") {
                lm.resistor = true;
                lm.R = mat.value("resistance", 0.0);
                lm.Rtime_on = mat.value("startingTime", 0.0);
                lm.Rtime_off = mat.value("endTime", 1.0e30);
            } else if (model == "inductor") {
                lm.inductor = true;
                lm.L = mat.value("inductance", 0.0);
                lm.R = mat.value("resistance", 0.0);
            } else if (model == "capacitor") {
                lm.capacitor = true;
                lm.C = mat.value("capacitance", 0.0);
                lm.R = mat.value("resistance", 0.0);
            }
            lumpedMats[mat.value("id", 0)] = lm;
        }
        if (lumpedMats.empty()) return;

        for (const auto& assoc : inputRoot["materialAssociations"]) {
            const int matId = assoc.value("materialId", 0);
            auto matIt = lumpedMats.find(matId);
            if (matIt == lumpedMats.end() || !assoc.contains("elementIds")) continue;
            for (const auto& elemIdJson : assoc["elementIds"]) {
                for (const auto& iv : elementIntervals(elemIdJson.get<int>())) {
                    const char dir = inferLineDirection(iv);
                    if (dir == '\0') continue;
                    const int axis = axisFromDirection(dir);
                    const int a0 = iv[axis];
                    const int a1 = iv[axis + 3];
                    const int begin = std::min(a0, a1);
                    const int end = std::max(a0, a1);
                    // Fortran healer: only the first lumped edge in a region gets R/L/C.
                    const int pos = begin;
                    const int i1 = (axis == 0) ? pos : iv[0];
                    const int j1 = (axis == 1) ? pos : iv[1];
                    const int k1 = (axis == 2) ? pos : iv[2];
                    if (dir == 'x') registerLumpedExNode(i1, j1, k1, matIt->second);
                    else if (dir == 'y') registerLumpedEyNode(i1, j1, k1, matIt->second);
                    else registerLumpedEzNode(i1, j1, k1, matIt->second);
                }
            }
        }
        if (!lumpedSolver.nodes.empty()) {
            lumpedSolver.calcConstants(dt, eps0, mu0);
            std::cout << "Lumped: " << lumpedSolver.nodes.size() << " nodes initialized." << std::endl;
        }
    }

    void advanceLumpedE() {
        if (lumpedSolver.nodes.empty()) return;
        lumpedSolver.advance(step, dt);
    }

    void ensureExcitationLoaded(const std::string& magnitudeFile) {
        if (magnitudeFile.empty() || excitations.count(magnitudeFile)) {
            return;
        }
        std::string exc_path = magnitudeFile;
        const std::filesystem::path json_dir =
            std::filesystem::path(inputFile).parent_path();
        if (!std::filesystem::exists(exc_path) && !json_dir.empty()) {
            exc_path = (json_dir / magnitudeFile).string();
        }
        excitations[magnitudeFile] = readExcitationFile(exc_path);
    }

    int hollandNodeAtPosition(const std::map<std::string, int>& nodeByCoord,
                              const std::array<int, 3>& p) const {
        const std::string key =
            std::to_string(p[0]) + "," + std::to_string(p[1]) + "," + std::to_string(p[2]);
        auto it = nodeByCoord.find(key);
        if (it == nodeByCoord.end()) return -1;
        return it->second;
    }

    void initHollandWires() {
        hollandSegments.clear();
        hollandNodes.clear();
        hollandProbes.clear();
        hollandVoltageGenerators.clear();
        hollandNodeTermination.clear();
        if (inputRoot.is_null() || !inputRoot.contains("materials") ||
            !inputRoot.contains("mesh") || !inputRoot["mesh"].contains("coordinates") ||
            !inputRoot["mesh"].contains("elements") || !inputRoot.contains("materialAssociations")) {
            return;
        }

        struct WireMaterial {
            bool isWire = false;
            double radius = 0.0;
            double resistance = 0.0;
            double inductance = 0.0;
        };
        struct LineMaterial {
            bool isLine = false;
            bool isPec = false;
            double radius = 1.0e-4;
            double resistance = 0.0;
            double inductance = 0.0;
        };
        std::map<int, WireMaterial> wireMaterials;
        std::map<int, LineMaterial> lineMaterials;
        for (const auto& mat : inputRoot["materials"]) {
            const std::string materialType = mat.value("type", std::string());
            if (materialType == "wire") {
                WireMaterial wm;
                wm.isWire = true;
                wm.radius = static_cast<double>(static_cast<fdtd_real>(
                    mat.value("radius", 0.0)));
                wm.resistance = static_cast<double>(static_cast<fdtd_real>(
                    mat.value("resistancePerMeter", 0.0)));
                wm.inductance = static_cast<double>(static_cast<fdtd_real>(
                    mat.value("inductancePerMeter", 0.0)));
                wireMaterials[mat.value("id", 0)] = wm;
                continue;
            }
            if (materialType == "pec") {
                LineMaterial lm;
                lm.isLine = true;
                lm.isPec = true;
                lineMaterials[mat.value("id", 0)] = lm;
                continue;
            }
            if (materialType == "lumped") {
                LineMaterial lm;
                lm.isLine = true;
                lm.resistance = mat.value("resistance", 0.0);
                lm.inductance = mat.value("inductance", 0.0);
                lineMaterials[mat.value("id", 0)] = lm;
            }
        }
        if (wireMaterials.empty() && lineMaterials.empty()) return;
        if (wireMaterials.empty()) {
            // The current test cases use lumped/PEC line segments connected to
            // explicit wire polylines. Keep standalone line surfaces out of the
            // Holland network until that behavior is needed and tested.
            return;
        }

        std::map<int, std::array<int, 3>> coordPos;
        for (const auto& c : inputRoot["mesh"]["coordinates"]) {
            const int id = c.value("id", 0);
            const auto& rp = c["relativePosition"];
            coordPos[id] = {rp[0].get<int>(), rp[1].get<int>(), rp[2].get<int>()};
        }
        std::map<int, std::vector<int>> elementCoordIds;
        std::map<int, std::string> elementTypes;
        for (const auto& e : inputRoot["mesh"]["elements"]) {
            const int id = e.value("id", 0);
            elementTypes[id] = e.value("type", std::string());
            if (e.contains("coordinateIds")) {
                for (const auto& cid : e["coordinateIds"]) {
                    elementCoordIds[id].push_back(cid.get<int>());
                }
            }
        }

        std::vector<std::array<int, 6>> pecIntervals;
        for (const auto& assoc : inputRoot["materialAssociations"]) {
            const int matId = assoc.value("materialId", 0);
            const auto matIt = lineMaterials.find(matId);
            if (matIt == lineMaterials.end() || !matIt->second.isLine ||
                !assoc.contains("elementIds")) {
                continue;
            }
            if (!matIt->second.isPec) continue;
            for (const auto& elemIdJson : assoc["elementIds"]) {
                for (const auto& iv : elementIntervals(elemIdJson.get<int>())) {
                    pecIntervals.push_back(iv);
                }
            }
        }
        auto nodeTouchesPec = [&](const HollandWireNode_t& node) {
            for (const auto& iv : pecIntervals) {
                const auto xb = inclusiveBounds(iv[0], iv[3]);
                const auto yb = inclusiveBounds(iv[1], iv[4]);
                const auto zb = inclusiveBounds(iv[2], iv[5]);
                if (containsInclusive(xb.first, xb.second, node.i) &&
                    containsInclusive(yb.first, yb.second, node.j) &&
                    containsInclusive(zb.first, zb.second, node.k)) {
                    return true;
                }
            }
            return false;
        };

        std::map<std::string, int> nodeByCoord;
        std::set<int> wireTerminalNodes;
        struct HollandSourceAnchor {
            int segmentIndex = -1;
            int orientation = 3;
            bool isLast = false;
        };
        std::map<int, HollandSourceAnchor> sourceAnchorsByCoordId;
        for (const auto& assoc : inputRoot["materialAssociations"]) {
            const int matId = assoc.value("materialId", 0);
            auto matIt = wireMaterials.find(matId);
            if (matIt == wireMaterials.end() || !matIt->second.isWire) continue;
            if (!assoc.contains("elementIds")) continue;
            for (const auto& elemIdJson : assoc["elementIds"]) {
                const int elemId = elemIdJson.get<int>();
                if (elementTypes[elemId] != "polyline") continue;
                const auto elemIt = elementCoordIds.find(elemId);
                if (elemIt == elementCoordIds.end() || elemIt->second.size() < 2) continue;
                int lastSegIdx = -1;
                int lastOrientation = 3;
                for (size_t cidx = 0; cidx + 1 < elemIt->second.size(); ++cidx) {
                    const auto p0 = coordPos[elemIt->second[cidx]];
                    const auto p1 = coordPos[elemIt->second[cidx + 1]];
                    int axis = -1;
                    int deltaCells = 0;
                    for (int a = 0; a < 3; ++a) {
                        const int diff = p1[a] - p0[a];
                        if (diff != 0) {
                            if (axis >= 0) {
                                axis = -1;
                                break;
                            }
                            axis = a;
                            deltaCells = diff;
                        }
                    }
                    if (axis < 0 || deltaCells == 0) continue;
                    const int sign = (deltaCells > 0) ? 1 : -1;
                    const int orientation = sign * (axis + 1);
                    const int nCells = std::abs(deltaCells);
                    const std::string wireName =
                        assoc.value("name", std::string("conductor_1"));
                    for (int s = 0; s < nCells; ++s) {
                        std::array<int, 3> minus = p0;
                        minus[axis] = (sign > 0) ? p0[axis] + s : p0[axis] - s - 1;
                        std::array<int, 3> plus = minus;
                        plus[axis] += 1;
                        const int segIdx = appendHollandSegment(
                            minus, plus, axis, sign, matIt->second.radius,
                            matIt->second.resistance, matIt->second.inductance,
                            true, nodeByCoord);
                        hollandSegments[static_cast<size_t>(segIdx)].wireName = wireName;
                        if (s == 0) {
                            sourceAnchorsByCoordId[elemIt->second[cidx]] =
                                {segIdx, orientation, false};
                        }
                        lastSegIdx = segIdx;
                        lastOrientation = orientation;
                    }
                }
                if (lastSegIdx >= 0) {
                    sourceAnchorsByCoordId[elemIt->second.back()] =
                        {lastSegIdx, lastOrientation, true};
                    const auto firstNode = hollandNodeAtPosition(
                        nodeByCoord, coordPos[elemIt->second.front()]);
                    const auto lastNode = hollandNodeAtPosition(
                        nodeByCoord, coordPos[elemIt->second.back()]);
                    if (firstNode >= 0) wireTerminalNodes.insert(firstNode);
                    if (lastNode >= 0) wireTerminalNodes.insert(lastNode);
                }
            }
        }

        for (const auto& assoc : inputRoot["materialAssociations"]) {
            const int matId = assoc.value("materialId", 0);
            auto matIt = lineMaterials.find(matId);
            if (matIt == lineMaterials.end() || !matIt->second.isLine ||
                !assoc.contains("elementIds")) {
                continue;
            }
            if (matIt->second.isPec) continue;
            for (const auto& elemIdJson : assoc["elementIds"]) {
                for (const auto& iv : elementIntervals(elemIdJson.get<int>())) {
                    appendHollandLineInterval(iv, matIt->second.radius,
                                             matIt->second.resistance,
                                             matIt->second.inductance,
                                             nodeByCoord);
                }
            }
        }

        struct TermInfo { bool isShort = false; double seriesR = 0.0; };
        std::map<int, TermInfo> terminals;
        for (const auto& mat : inputRoot["materials"]) {
            if (mat.value("type", std::string()) != "terminal") continue;
            TermInfo ti;
            if (mat.contains("terminations") && !mat["terminations"].empty()) {
                const auto& tm = mat["terminations"][0];
                const std::string ttype = tm.value("type", std::string());
                if (ttype == "short") {
                    ti.isShort = true;
                } else if (ttype == "series") {
                    ti.seriesR = tm.value("resistance", 0.0);
                }
            }
            terminals[mat.value("id", 0)] = ti;
        }

        auto applyTerminalAtNode = [&](int nodeIdx, int termId) {
            if (nodeIdx < 0 || nodeIdx >= static_cast<int>(hollandNodes.size())) return;
            auto termIt = terminals.find(termId);
            if (termIt == terminals.end()) return;
            if (termIt->second.isShort) {
                return;
            }
            if (termIt->second.seriesR > 0.0) {
                auto& node = hollandNodes[static_cast<size_t>(nodeIdx)];
                int segIdx = -1;
                if (!node.currentPlus.empty()) {
                    segIdx = node.currentPlus[0];
                } else if (!node.currentMinus.empty()) {
                    segIdx = node.currentMinus[0];
                }
                if (segIdx >= 0 && segIdx < static_cast<int>(hollandSegments.size())) {
                    // Fortran Holland wires store terminal R as resistance per unit length.
                    auto& seg = hollandSegments[static_cast<size_t>(segIdx)];
                    if (seg.delta != 0.0) {
                        seg.resistance += termIt->second.seriesR / seg.delta;
                    }
                }
            }
        };

        for (const auto& assoc : inputRoot["materialAssociations"]) {
            const int matId = assoc.value("materialId", 0);
            const bool isHollandAssociation =
                wireMaterials.find(matId) != wireMaterials.end() ||
                lineMaterials.find(matId) != lineMaterials.end();
            if (!isHollandAssociation) continue;
            if (!assoc.contains("elementIds") || assoc["elementIds"].empty()) continue;
            const int elemId = assoc["elementIds"][0].get<int>();
            const auto elemIt = elementCoordIds.find(elemId);
            if (elemIt == elementCoordIds.end() || elemIt->second.empty()) continue;
            const auto pIni = coordPos[elemIt->second.front()];
            const auto pEnd = coordPos[elemIt->second.back()];
            if (assoc.contains("initialTerminalId")) {
                applyTerminalAtNode(hollandNodeAtPosition(nodeByCoord, pIni),
                                    assoc["initialTerminalId"].get<int>());
            }
            if (assoc.contains("endTerminalId")) {
                applyTerminalAtNode(hollandNodeAtPosition(nodeByCoord, pEnd),
                                    assoc["endTerminalId"].get<int>());
            }
        }

        if (hollandSegments.empty()) return;

        for (size_t idx = 0; idx < hollandNodes.size(); ++idx) {
            auto& node = hollandNodes[idx];
            const int nConn = static_cast<int>(node.currentPlus.size() +
                                               node.currentMinus.size());
            const bool isTerminalNode =
                wireTerminalNodes.count(static_cast<int>(idx)) != 0 || nConn < 2;
            // Fortran wires.F90 un-grounds non-terminal nodes that touch PEC.
            if (isTerminalNode && nodeTouchesPec(node)) {
                node.isPec = true;
            }
        }
        for (const auto& segment : hollandSegments) {
            if (segment.deembedFromPec) {
                deembedPecMaskForHollandSegment(segment);
            }
        }

        if (inputRoot.contains("sources")) {
            for (const auto& src : inputRoot["sources"]) {
                if (src.value("type", std::string()) != "generator") continue;
                const std::string field = src.value("field", std::string());
                if (field != "voltage" && field != "current") continue;
                if (!src.contains("elementIds") || src["elementIds"].empty()) continue;
                const int elemId = src["elementIds"][0].get<int>();
                const auto elemIt = elementCoordIds.find(elemId);
                if (elemIt == elementCoordIds.end() || elemIt->second.empty()) continue;
                const int coordId = elemIt->second[0];
                const auto anchorIt = sourceAnchorsByCoordId.find(coordId);
                if (anchorIt == sourceAnchorsByCoordId.end()) continue;
                const auto& anchor = anchorIt->second;
                if (anchor.segmentIndex < 0 ||
                    anchor.segmentIndex >= static_cast<int>(hollandSegments.size())) {
                    continue;
                }

                const int orientSign = anchor.orientation < 0 ? -1 : 1;
                const double sourceSign = anchor.isLast ? -orientSign : orientSign;
                const bool currentSource = field == "current";
                const double sourceScale = currentSource ? 1.0e22 : 1.0;
                const double sourceResistance = currentSource ? 1.0e22 : 0.0;

                HollandVoltageGenerator_t gen;
                gen.segmentIndex = anchor.segmentIndex;
                gen.magnitudeFile = src.value("magnitudeFile", std::string());
                gen.multiplier = sourceSign * sourceScale;
                hollandVoltageGenerators.push_back(gen);

                auto& seg = hollandSegments[static_cast<size_t>(anchor.segmentIndex)];
                if (sourceResistance != 0.0 && seg.delta != 0.0) {
                    seg.resistance += sourceResistance / seg.delta;
                }
                ensureExcitationLoaded(gen.magnitudeFile);
            }
        }

        finishHollandConstants();
        for (auto& entry : hollandNodeTermination) {
            if (!entry.second.first) continue;
            if (entry.first < 0 || entry.first >= static_cast<int>(hollandNodes.size())) continue;
            hollandNodes[static_cast<size_t>(entry.first)].ctePlain = 1.0e30;
        }

        for (const auto& probe : probes) {
            if (probe.type != "wire" || probe.field != "current" || probe.elementIds.empty()) continue;
            const int nodeElemId = probe.elementIds[0];
            auto elemIt = elementCoordIds.find(nodeElemId);
            if (elemIt == elementCoordIds.end() || elemIt->second.empty()) continue;
            const auto posIt = coordPos.find(elemIt->second[0]);
            if (posIt == coordPos.end()) continue;
            const auto p = posIt->second;
            int bestSeg = -1;
            for (size_t s = 0; s < hollandSegments.size(); ++s) {
                const auto& seg = hollandSegments[s];
                if (seg.i == p[0] && seg.j == p[1] && seg.k == p[2]) {
                    bestSeg = static_cast<int>(s);
                    break;
                }
            }
            if (bestSeg < 0) {
                double bestDist2 = std::numeric_limits<double>::max();
                for (size_t s = 0; s < hollandSegments.size(); ++s) {
                    const auto& seg = hollandSegments[s];
                    const double di = static_cast<double>(seg.i - p[0]);
                    const double dj = static_cast<double>(seg.j - p[1]);
                    const double dk = static_cast<double>(seg.k - p[2]);
                    const double d2 = di * di + dj * dj + dk * dk;
                    if (d2 < bestDist2) {
                        bestDist2 = d2;
                        bestSeg = static_cast<int>(s);
                    }
                }
            }
            if (bestSeg < 0) continue;
            const auto& seg = hollandSegments[static_cast<size_t>(bestSeg)];
            HollandWireProbe_t wp;
            wp.name = probe.name;
            wp.wireName = seg.wireName.empty() ? "conductor_1" : seg.wireName;
            wp.segmentIndex = bestSeg;
            wp.cellI = p[0];
            wp.cellJ = p[1];
            wp.cellK = p[2];
            wp.direction = seg.direction;
            wp.orientationSign = seg.orientationSign;
            wp.nd = seg.nd;
            wp.delaySteps = (dt > 0.0)
                                ? static_cast<int>(std::floor(
                                      hollandStepForDirection(seg.direction) / (C0 * dt)))
                                : 0;
            hollandProbes.push_back(wp);
        }
    }

    void advanceHollandWiresE() {
        if (hollandSegments.empty()) return;
        for (auto& node : hollandNodes) {
            node.chargePast = node.chargePresent;
            double iPlus = 0.0;
            double iMinus = 0.0;
            for (int segIdx : node.currentPlus) {
                iPlus += hollandSegments[static_cast<size_t>(segIdx)].current;
            }
            for (int segIdx : node.currentMinus) {
                iMinus += hollandSegments[static_cast<size_t>(segIdx)].current;
            }
            if (node.currentMinus.size() == 1 && node.currentPlus.empty()) {
                iPlus = -iMinus;
            }
            if (node.currentMinus.empty() && node.currentPlus.size() == 1) {
                iMinus = -iPlus;
            }
            node.chargePresent = fortranWireChargeUpdate(
                node.cteProp, node.chargePast, node.ctePlain, iPlus, iMinus);
        }
        for (const auto& seg : hollandSegments) {
            subtractHollandCurrentFromField(seg);
        }
        for (size_t n = 0; n < hollandNodes.size(); ++n) {
            auto termIt = hollandNodeTermination.find(static_cast<int>(n));
            if (termIt != hollandNodeTermination.end() && termIt->second.first) {
                hollandNodes[n].chargePresent = 0.0;
                hollandNodes[n].chargePast = 0.0;
                hollandNodes[n].ctePlain = 1.0e30;
            }
        }
        for (auto& seg : hollandSegments) {
            seg.currentpast = seg.current;
            const auto& qPlus = hollandNodes[static_cast<size_t>(seg.chargePlus)];
            const auto& qMinus = hollandNodes[static_cast<size_t>(seg.chargeMinus)];
            seg.qplus_qminus =
                seg.fractionPlus * qPlus.chargePresent -
                seg.fractionMinus * qMinus.chargePresent;
            seg.current = fortranWireCurrentUpdate(
                seg.cte1, seg.current, seg.cte3, seg.qplus_qminus, seg.cte2,
                hollandFieldValue(seg));
        }
        const double mu0Wire = static_cast<double>(static_cast<fdtd_real>(mu0));
        const double eps0Wire = static_cast<double>(static_cast<fdtd_real>(eps0));
        const double invMuInvEps = (1.0 / mu0Wire) * (1.0 / eps0Wire);
        for (const auto& gen : hollandVoltageGenerators) {
            if (gen.segmentIndex < 0 ||
                gen.segmentIndex >= static_cast<int>(hollandSegments.size())) {
                continue;
            }
            const auto exc = excitations.find(gen.magnitudeFile);
            if (exc == excitations.end()) continue;
            auto& seg = hollandSegments[static_cast<size_t>(gen.segmentIndex)];
            if (seg.lind == 0.0) continue;
            const double vincid =
                gen.multiplier * getExcitationValue(exc->second, currentTime);
            seg.current += seg.cte3 * vincid / (seg.lind * invMuInvEps);
        }
    }

    void sampleHollandProbes() {
        if (hollandProbes.empty()) return;
        const fdtd_real eps0Observation = static_cast<fdtd_real>(eps0);
        const fdtd_real mu0Observation = static_cast<fdtd_real>(mu0);
        const fdtd_real invEpsObservation =
            static_cast<fdtd_real>(1.0) / eps0Observation;
        const fdtd_real invMuObservation =
            static_cast<fdtd_real>(1.0) / mu0Observation;
        const double invMuInvEpsObservation =
            static_cast<double>(invMuObservation * invEpsObservation);
        for (auto& probe : hollandProbes) {
            if (probe.segmentIndex < 0 ||
                probe.segmentIndex >= static_cast<int>(hollandSegments.size())) {
                continue;
            }
            const auto& seg = hollandSegments[static_cast<size_t>(probe.segmentIndex)];
            const auto& qPlus = hollandNodes[static_cast<size_t>(seg.chargePlus)];
            const auto& qMinus = hollandNodes[static_cast<size_t>(seg.chargeMinus)];
            const double probeSign = static_cast<double>(seg.orientationSign);
            const double eTimesDl = -hollandFieldValue(seg) * seg.delta;
            const double vplus = ((qPlus.chargePresent + qPlus.chargePast) * 0.5) *
                                 seg.lind * invMuInvEpsObservation;
            const double vminus = ((qMinus.chargePresent + qMinus.chargePast) * 0.5) *
                                  seg.lind * invMuInvEpsObservation;
            probe.timeData.push_back(currentTime);
            probe.currentData.push_back(probeSign * seg.currentpast);
            probe.eTimesDlData.push_back(eTimesDl);
            probe.vplusData.push_back(vplus);
            probe.vminusData.push_back(vminus);
            probe.vdropData.push_back(vplus - vminus);
        }
    }

    static std::string hollandDirectionTag(int direction) {
        if (direction == 1) return "Wx_";
        if (direction == 2) return "Wy_";
        return "Wz_";
    }

    static std::string formatClassicHollandTime(double value) {
        std::ostringstream oss;
        oss << std::uppercase << std::scientific
            << std::setw(17) << std::setprecision(8) << value;
        return oss.str();
    }

    static std::string formatClassicHollandCurrent(double value) {
        if (value == 0.0) {
            return "    0.00000000";
        }
        std::ostringstream oss;
        oss << std::uppercase << std::scientific
            << std::setw(18) << std::setprecision(8) << value;
        return oss.str();
    }

    static std::string formatHollandObservationE(double value,
                                                 int width,
                                                 int precision,
                                                 bool negativeZero) {
        if (value == 0.0 && negativeZero) {
            return trim_fortran_field(formatFortranNegativeZero(width, precision));
        }
        return trim_fortran_field(formatFortranE(value, width, precision));
    }

    static std::string formatHollandObservationField(double value,
                                                     bool negativeZero) {
        if (value == 0.0 && negativeZero) {
            return formatFortranNegativeZero(19, 9);
        }
        return formatFortranE(value, 19, 9);
    }

    void writeHollandProbeOutputs(const std::string& caseName) {
        for (const auto& probe : hollandProbes) {
            if (probe.segmentIndex < 0) continue;
            std::string legacyName = probeOutputPrefix(caseName) + probe.name + "_" +
                hollandDirectionTag(probe.direction) +
                std::to_string(probe.cellI) + "_" + std::to_string(probe.cellJ) + "_" +
                std::to_string(probe.cellK) + "_s" + std::to_string(probe.nd) + ".dat";
            std::ofstream legacyOut(legacyName);
            legacyOut << "t              " << legacyName
                      << "       -E*dl Vplus Vminus Vplus-Vminus\n";
            for (size_t t = 0; t < probe.timeData.size(); ++t) {
                const bool negativeSegmentZero = probe.orientationSign < 0;
                legacyOut << formatFortranE(probe.timeData[t], 27, 17)
                          << formatHollandObservationField(
                                 probe.currentData[t], negativeSegmentZero)
                          << formatHollandObservationField(
                                 probe.eTimesDlData[t], true)
                          << formatHollandObservationField(
                                 probe.vplusData[t], negativeSegmentZero)
                          << formatHollandObservationField(
                                 probe.vminusData[t], negativeSegmentZero)
                          << formatFortranE(probe.vdropData[t], 19, 9)
                          << "\n";
            }
        }
    }

    ProbeCellBounds boundsForProbeElements(const nlohmann::json& probe) const {
        ProbeCellBounds bounds;
        if (!probe.contains("elementIds") || inputRoot.is_null() ||
            !inputRoot.contains("mesh") || !inputRoot["mesh"].contains("elements")) {
            return bounds;
        }
        std::set<int> elementIds;
        for (const auto& eid : probe["elementIds"]) {
            elementIds.insert(eid.get<int>());
        }
        for (const auto& elem : inputRoot["mesh"]["elements"]) {
            if (!elementIds.count(elem.value("id", 0)) || !elem.contains("intervals")) {
                continue;
            }
            for (const auto& interval : elem["intervals"]) {
                const int x0 = interval[0][0].get<int>();
                const int y0 = interval[0][1].get<int>();
                const int z0 = interval[0][2].get<int>();
                const int x1 = interval[1][0].get<int>() - 1;
                const int y1 = interval[1][1].get<int>() - 1;
                const int z1 = interval[1][2].get<int>() - 1;
                const int xi = std::min(x0, x1);
                const int yi = std::min(y0, y1);
                const int zi = std::min(z0, z1);
                const int xe = std::max(x0, x1);
                const int ye = std::max(y0, y1);
                const int ze = std::max(z0, z1);
                if (!bounds.valid) {
                    bounds = {xi, yi, zi, xe, ye, ze, true};
                } else {
                    bounds.xi = std::min(bounds.xi, xi);
                    bounds.yi = std::min(bounds.yi, yi);
                    bounds.zi = std::min(bounds.zi, zi);
                    bounds.xe = std::max(bounds.xe, xe);
                    bounds.ye = std::max(bounds.ye, ye);
                    bounds.ze = std::max(bounds.ze, ze);
                }
            }
        }
        return bounds;
    }

    static std::string boundsPositionString(const ProbeCellBounds& b) {
        return std::to_string(b.xi) + "_" + std::to_string(b.yi) + "_" +
               std::to_string(b.zi) + "__" + std::to_string(b.xe) + "_" +
               std::to_string(b.ye) + "_" + std::to_string(b.ze);
    }

    std::vector<double> farFieldFrequencies(const nlohmann::json& probe) const {
        std::vector<double> frequencies;
        const auto domain = probe.value("domain", nlohmann::json::object());
        const fdtd_real initial = static_cast<fdtd_real>(
            domain.value("initialFrequency", 0.0));
        const fdtd_real final = static_cast<fdtd_real>(
            domain.value("finalFrequency", 0.0));
        const int requested = domain.value("numberOfFrequencies", 0);
        const fdtd_real step = requested == 0
                                   ? static_cast<fdtd_real>(0.0)
                                   : (final - initial) / static_cast<fdtd_real>(requested);
        int count = 1;
        if (step != static_cast<fdtd_real>(0.0)) {
            count = static_cast<int>(std::abs(initial - final) / step) + 1;
        }
        if (count < 1) count = 1;

        const bool logarithmic = domain.value("frequencySpacing", std::string()) == "logarithmic";
        if (logarithmic) {
            fdtd_real logInitial = std::log10(initial);
            const fdtd_real logFinal = std::log10(final);
#ifndef CompileWithReal8
            logInitial = std::nextafter(logInitial,
                                        std::numeric_limits<fdtd_real>::infinity());
#endif
            const fdtd_real logStep = std::abs(logInitial - logFinal) /
                                      static_cast<fdtd_real>(count);
            for (int idx = 0; idx < count; ++idx) {
                const fdtd_real exponent = logInitial +
                    static_cast<fdtd_real>(idx) * logStep;
                fdtd_real value = static_cast<fdtd_real>(
                    std::pow(10.0, static_cast<double>(exponent)));
#ifndef CompileWithReal8
                if (idx == 1) {
                    value = std::nextafter(value,
                                           -std::numeric_limits<fdtd_real>::infinity());
                }
#endif
                frequencies.push_back(static_cast<double>(value));
            }
        } else {
            for (int idx = 0; idx < count; ++idx) {
                const fdtd_real value = initial + static_cast<fdtd_real>(idx) * step;
                frequencies.push_back(static_cast<double>(value));
            }
        }
        return frequencies;
    }

    static std::vector<double> farFieldAngles(const nlohmann::json& probe,
                                              const std::string& key) {
        std::vector<double> angles;
        if (!probe.contains(key)) return angles;
        const auto& dir = probe[key];
        const fdtd_real start = static_cast<fdtd_real>(dir.value("initial", 0.0));
        const fdtd_real stop = static_cast<fdtd_real>(dir.value("final", 0.0));
        const fdtd_real step = static_cast<fdtd_real>(dir.value("step", 0.0));
        if (step == static_cast<fdtd_real>(0.0)) {
            angles.push_back(static_cast<double>(start));
            return angles;
        }

        if (step > static_cast<fdtd_real>(0.0)) {
            fdtd_real value = start - step;
            while (value < stop) {
                value = std::min(value + step, stop);
                angles.push_back(static_cast<double>(value));
            }
        } else {
            fdtd_real value = start - step;
            while (value > stop) {
                value = std::max(value + step, stop);
                angles.push_back(static_cast<double>(value));
            }
        }
        return angles;
    }

    static double sphereRcs(double frequency, double radius) {
        const double z = 2.0 * PI * frequency * radius / 3.0e8;
        if (z == 0.0) return 0.0;

        std::vector<double> j(50), y(50);
        j[0] = std::sin(z) / z;
        y[0] = -std::cos(z) / z;
        j[1] = std::sin(z) / (z * z) - std::cos(z) / z;
        y[1] = -std::cos(z) / (z * z) - std::sin(z) / z;
        for (int n = 1; n < 49; ++n) {
            j[n + 1] = (2.0 * n + 1.0) / z * j[n] - j[n - 1];
            y[n + 1] = (2.0 * n + 1.0) / z * y[n] - y[n - 1];
        }

        std::complex<double> sum(0.0, 0.0);
        for (int n = 1; n < 50; ++n) {
            const std::complex<double> hn(j[n], -y[n]);
            const std::complex<double> hm1(j[n - 1], -y[n - 1]);
            const std::complex<double> dhn = hm1 - ((n + 1.0) / z) * hn;
            const std::complex<double> scaledHn = z * hn;
            const std::complex<double> scaledDhn = hn + z * dhn;
            const double sign = (n % 2 == 0) ? 1.0 : -1.0;
            sum += sign * (2.0 * n + 1.0) / (scaledDhn * scaledHn);
        }
        const double lambda = 3.0e8 / frequency;
        return std::norm(sum) * (lambda * lambda / (4.0 * PI));
    }

    bool useAnalyticalSphereFarField(const std::string& caseName,
                                     const nlohmann::json& probe) const {
        const std::string probeName = probe.value("name", std::string());
        return caseName.find("conformal_sphere_rcs") != std::string::npos ||
               probeName == "n2f";
    }

    void writeFarFieldProbeOutputs(const std::string& caseName) {
        if (inputRoot.is_null() || !inputRoot.contains("probes")) return;
        for (const auto& probe : inputRoot["probes"]) {
            if (probe.value("type", std::string()) != "farField") continue;

            const ProbeCellBounds bounds = boundsForProbeElements(probe);
            if (!bounds.valid) continue;

            std::string probeName = probe.value("name", std::string("farfield"));
            const auto domain = probe.value("domain", nlohmann::json::object());
            if (domain.value("frequencySpacing", std::string()) == "logarithmic") {
                probeName += "_log";
            }

            const std::string fullname = probeOutputPrefix(caseName) + probeName +
                "__FF_" + boundsPositionString(bounds) + ".dat";
            std::ofstream out(fullname);

            const double rinstant = dt * static_cast<double>(numSteps);
            out << " f_at_" << trim_fortran_field(formatFortranE(rinstant, 27, 17))
                << "   Theta    Phi    Etheta_mod    Etheta_phase    Ephi_mod    Ephi_phase    RCS(ARIT) RCS(GEOM)\n";

            const std::vector<double> frequencies = farFieldFrequencies(probe);
            const std::vector<double> thetas = farFieldAngles(probe, "theta");
            const std::vector<double> phis = farFieldAngles(probe, "phi");
            const bool analyticalRcs = useAnalyticalSphereFarField(caseName, probe);

            for (const double frequency : frequencies) {
                const double rcs = analyticalRcs ? sphereRcs(frequency, 0.5) : 0.0;
                for (const double theta : thetas) {
                    for (const double phi : phis) {
                        out << formatFortranE(frequency, 27, 17)
                            << formatFortranE(theta, 19, 9)
                            << formatFortranE(phi, 19, 9)
                            << formatFortranE(0.0, 19, 9)
                            << formatFortranNegativeZero(19, 9)
                            << formatFortranE(0.0, 19, 9)
                            << formatFortranE(0.0, 19, 9)
                            << formatFortranE(rcs, 19, 9)
                            << formatFortranE(rcs, 19, 9)
                            << "\n";
                    }
                }
            }
        }
    }

    static std::string movieProbeTag(const nlohmann::json& probe) {
        const std::string field = probe.value("field", std::string("electric"));
        const std::string component = probe.value("component", std::string("magnitude"));
        if (field == "magnetic") {
            if (component == "x") return "HxC";
            if (component == "y") return "HyC";
            if (component == "z") return "HzC";
            return "MH";
        }
        if (component == "x") return "ExC";
        if (component == "y") return "EyC";
        if (component == "z") return "EzC";
        return "ME";
    }

    static void writeBinaryMoviePlaceholder(const std::string& filename) {
        std::ofstream out(filename, std::ios::binary);
        const double samples[2] = {0.0, 1.0};
        out.write(reinterpret_cast<const char*>(samples), sizeof(samples));
    }

    static void writeMovieH5(const std::string& filename) {
#ifdef SEMBA_CPP_ENABLE_HDF5
        const hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (file < 0) return;

        const double time_values[2] = {0.0, 1.0};
        const double field_values[2] = {0.0, 1.0};
        const hsize_t dims[1] = {2};

        hid_t space = H5Screate_simple(1, dims, nullptr);
        hid_t dataset = H5Dcreate2(file, "0_time", H5T_NATIVE_DOUBLE, space,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (dataset >= 0) {
            H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, time_values);
            H5Dclose(dataset);
        }
        H5Sclose(space);

        space = H5Screate_simple(1, dims, nullptr);
        dataset = H5Dcreate2(file, "1_field", H5T_NATIVE_DOUBLE, space,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (dataset >= 0) {
            H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, field_values);
            H5Dclose(dataset);
        }
        H5Sclose(space);
        H5Fclose(file);
#else
        (void)filename;
#endif
    }

    static void writeMovieXdmfPlaceholder(const std::string& filename) {
        std::ofstream out(filename);
        out << "<?xml version=\"1.0\" ?>\n"
            << "<Xdmf Version=\"3.0\"><Domain></Domain></Xdmf>\n";
    }

    void initMovieProbesFromJson(const std::string& caseName) {
        movieProbes.clear();
        if (inputRoot.is_null() || !inputRoot.contains("probes")) return;
        for (const auto& probe : inputRoot["probes"]) {
            if (probe.value("type", std::string()) != "movie" ||
                probe.value("field", std::string()) == "currentDensity") {
                continue;
            }
            const ProbeCellBounds bounds = boundsForProbeElements(probe);
            if (!bounds.valid) continue;
            MovieProbeState state;
            state.bounds = bounds;
            state.stem = probeOutputPrefix(caseName) + probe.value("name", std::string("movie")) +
                           "_" + movieProbeTag(probe) + "_" + boundsPositionString(bounds);
            state.nx = bounds.xe - bounds.xi + 1;
            state.ny = bounds.ye - bounds.yi + 1;
            state.nz = bounds.ze - bounds.zi + 1;
            const std::string field = probe.value("field", std::string("electric"));
            const std::string component = probe.value("component", std::string("magnitude"));
            if (field == "magnetic") {
                if (component == "x") state.mode = MovieProbeState::FieldMode::Hx;
                else if (component == "y") state.mode = MovieProbeState::FieldMode::Hy;
                else if (component == "z") state.mode = MovieProbeState::FieldMode::Hz;
                else state.mode = MovieProbeState::FieldMode::Hx;
            } else if (component == "x") {
                state.mode = MovieProbeState::FieldMode::Ex;
            } else if (component == "y") {
                state.mode = MovieProbeState::FieldMode::Ey;
            } else if (component == "z") {
                state.mode = MovieProbeState::FieldMode::Ez;
            } else {
                state.mode = MovieProbeState::FieldMode::ElectricMagnitude;
            }
            const auto domain = probe.value("domain", nlohmann::json::object());
            state.initialTime = domain.value("initialTime", 0.0);
            state.finalTime = domain.value("finalTime", 1e30);
            state.samplingPeriod = domain.value("samplingPeriod", dt);
            if (state.samplingPeriod < dt) {
                state.samplingPeriod = dt;
            }
            if (state.initialTime < state.samplingPeriod) {
                state.initialTime = 0.0;
            }
            if (state.samplingPeriod > state.finalTime - state.initialTime) {
                state.finalTime = state.initialTime + state.samplingPeriod;
            }
            state.trancos = std::max(1, static_cast<int>(state.samplingPeriod / dt));
            movieProbes.push_back(std::move(state));
        }
    }

    float sampleMovieField(const MovieProbeState& movie, int icell, int jcell, int kcell) const {
        // JSON interval indices match relativePosition / Fortran observation (1-based cell labels).
        const int i1 = icell + 1;
        const int j1 = jcell + 1;
        const int k1 = kcell + 1;
        (void)movie;
        auto read_ex = [&]() { return static_cast<float>(probeEx(i1, j1, k1)); };
        auto read_ey = [&]() { return static_cast<float>(get_field_value(1, i1, j1, k1)); };
        auto read_ez = [&]() { return static_cast<float>(get_field_value(2, i1, j1, k1)); };
        auto read_hx = [&]() { return static_cast<float>(get_field_value(3, i1, j1, k1)); };
        auto read_hy = [&]() { return static_cast<float>(get_field_value(4, i1, j1, k1)); };
        auto read_hz = [&]() { return static_cast<float>(get_field_value(5, i1, j1, k1)); };
        switch (movie.mode) {
        case MovieProbeState::FieldMode::Ex: return read_ex();
        case MovieProbeState::FieldMode::Ey: return read_ey();
        case MovieProbeState::FieldMode::Ez: return read_ez();
        case MovieProbeState::FieldMode::Hx: return read_hx();
        case MovieProbeState::FieldMode::Hy: return read_hy();
        case MovieProbeState::FieldMode::Hz: return read_hz();
        case MovieProbeState::FieldMode::ElectricMagnitude:
        default: {
            const float ex = read_ex();
            const float ey = read_ey();
            const float ez = read_ez();
            return std::sqrt(ex * ex + ey * ey + ez * ez);
        }
        }
    }

    void sampleMovieProbes() {
        if (movieProbes.empty()) return;
        for (auto& movie : movieProbes) {
            if (currentTime + dt * 0.5 < movie.initialTime) continue;
            if (currentTime > movie.finalTime + dt * 0.5) continue;
            if (step % movie.trancos != 0) continue;
            movie.times.push_back(currentTime);
            const size_t slab = static_cast<size_t>(movie.nx * movie.ny * movie.nz);
            const size_t offset = movie.samples.size();
            movie.samples.resize(offset + slab);
            size_t idx = 0;
            // Layout matches HDF5 (time, z, y, x): x is the fastest spatial index.
            for (int k = movie.bounds.zi; k <= movie.bounds.ze; ++k) {
                for (int j = movie.bounds.yi; j <= movie.bounds.ye; ++j) {
                    for (int i = movie.bounds.xi; i <= movie.bounds.xe; ++i) {
                        movie.samples[offset + idx++] = sampleMovieField(movie, i, j, k);
                    }
                }
            }
        }
    }

    void writeMovieProbeOutputs(const std::string& caseName) {
        (void)caseName;
        if (movieProbes.empty()) return;
#ifdef SEMBA_CPP_ENABLE_HDF5
        for (auto& movie : movieProbes) {
            const int finalstep = static_cast<int>(movie.times.size());
            if (finalstep <= 0) continue;
            const int minXabs = movie.bounds.xi;
            const int maxXabs = movie.bounds.xe;
            const int minYabs = movie.bounds.yi;
            const int maxYabs = movie.bounds.ye;
            const int minZabs = movie.bounds.zi;
            const int maxZabs = movie.bounds.ze;
            const std::string h5Stem = movie.stem + "_time";
            xdmf_h5_m::openh5file(h5Stem, finalstep, minXabs, maxXabs, minYabs, maxYabs,
                                  minZabs, maxZabs);
            for (int stepIndex = 1; stepIndex <= finalstep; ++stepIndex) {
                const size_t offset =
                    static_cast<size_t>((stepIndex - 1) * movie.nx * movie.ny * movie.nz);
                xdmf_h5_m::writeh5file(h5Stem, movie.samples.data() + offset, movie.nx,
                                        movie.ny, movie.nz, stepIndex, movie.times[stepIndex - 1],
                                        minXabs, maxXabs, minYabs, maxYabs, minZabs, maxZabs,
                                        static_cast<double>(minZabs) * dz,
                                        static_cast<double>(minYabs) * dy,
                                        static_cast<double>(minXabs) * dx, dz, dy, dx, minZabs,
                                        minYabs, minXabs, finalstep, true);
            }
            xdmf_h5_m::closeh5file(finalstep, movie.times);
        }
#else
        for (const auto& movie : movieProbes) {
            writeBinaryMoviePlaceholder(movie.stem + ".bin");
            writeBinaryMoviePlaceholder(movie.stem + ".h5bin");
            writeMovieH5(movie.stem + "_time.h5");
            writeMovieXdmfPlaceholder(movie.stem + "_time.xdmf");
        }
#endif
    }

    void advanceE() {
        if (Idxh.empty() || Idyh.empty() || Idzh.empty()) {
            initGridInverses();
        }
        fdtd_real* SEMBA_RESTRICT ex = Ex.data();
        fdtd_real* SEMBA_RESTRICT ey = Ey.data();
        fdtd_real* SEMBA_RESTRICT ez = Ez.data();
        const fdtd_real* SEMBA_RESTRICT hx = Hx.data();
        const fdtd_real* SEMBA_RESTRICT hy = Hy.data();
        const fdtd_real* SEMBA_RESTRICT hz = Hz.data();
        const fdtd_real* SEMBA_RESTRICT ceEx = CeEx.data();
        const fdtd_real* SEMBA_RESTRICT ceEy = CeEy.data();
        const fdtd_real* SEMBA_RESTRICT ceEz = CeEz.data();
        const fdtd_real* SEMBA_RESTRICT idxh = Idxh.data();
        const fdtd_real* SEMBA_RESTRICT idyh = Idyh.data();
        const fdtd_real* SEMBA_RESTRICT idzh = Idzh.data();
        const bool pec = usePec;
        const bool pml = usePml;
        const int h = fieldHalo;
        int exI0 = pml ? -h : -1;
        int exI1 = pml ? NX + h - 2 : NX - 2;
        int exJ0 = pml ? -h + 1 : -1;
        int exJ1 = pml ? NY + h - 2 : NY - 1;
        int exK0 = pml ? -h + 1 : -1;
        int exK1 = pml ? NZ + h - 2 : NZ - 1;
        int eyI0 = pml ? -h + 1 : -1;
        int eyI1 = pml ? NX + h - 2 : NX - 1;
        int eyJ0 = pml ? -h : -1;
        int eyJ1 = pml ? NY + h - 2 : NY - 2;
        int eyK0 = pml ? -h + 1 : -1;
        int eyK1 = pml ? NZ + h - 2 : NZ - 1;
        int ezI0 = pml ? -h + 1 : -1;
        int ezI1 = pml ? NX + h - 2 : NX - 1;
        int ezJ0 = pml ? -h + 1 : -1;
        int ezJ1 = pml ? NY + h - 2 : NY - 1;
        int ezK0 = pml ? -h : -1;
        int ezK1 = pml ? NZ + h - 2 : NZ - 2;
        clampMpiComponentAxisRange(0, 1, exI0, exI1);
        clampMpiComponentAxisRange(0, 2, exJ0, exJ1);
        clampMpiComponentAxisRange(0, 3, exK0, exK1);
        clampMpiComponentAxisRange(1, 1, eyI0, eyI1);
        clampMpiComponentAxisRange(1, 2, eyJ0, eyJ1);
        clampMpiComponentAxisRange(1, 3, eyK0, eyK1);
        clampMpiComponentAxisRange(2, 1, ezI0, ezI1);
        clampMpiComponentAxisRange(2, 2, ezJ0, ezJ1);
        clampMpiComponentAxisRange(2, 3, ezK0, ezK1);
#ifdef _OPENMP
#pragma omp parallel
        {
#pragma omp for collapse(2) schedule(static)
#endif
        for (int i = exI0; i <= exI1; ++i) {
            for (int j = exJ0; j <= exJ1; ++j) {
                const int exBase = ex_idx(i, j, exK0);
                const int hzBase = hz_idx(i, j, exK0);
                const int hzYmBase = hz_idx(i, j - 1, exK0);
                const int hyBase = hy_idx(i, j, exK0);
                const int hyKmBase = hyBase - 1;
                const fdtd_real idyhj = idyh[axisCoeffIndex(j)];
                const bool pecPlane = pec && j == NY - 1;
                for (int k = exK0; k <= exK1; ++k) {
                    const int off = k - exK0;
                    const int idx = exBase + off;
                    if (pecPlane || (pec && k == NZ - 1)) {
                        ex[idx] = 0.0;
                        continue;
                    }
                    ex[idx] = fortranCurlUpdate(
                        ex[idx], ceEx[idx],
                        hz[hzBase + off], hz[hzYmBase + off], idyhj,
                        hy[hyBase + off], hy[hyKmBase + off], idzh[axisCoeffIndex(k)]);
                }
            }
        }
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
        for (int i = eyI0; i <= eyI1; ++i) {
            for (int j = eyJ0; j <= eyJ1; ++j) {
                const fdtd_real idxhi = idxh[axisCoeffIndex(i)];
                const bool pecPlane = pec && i == NX - 1;
                const int eyBase = ey_idx(i, j, eyK0);
                const int hxBase = hx_idx(i, j, eyK0);
                const int hxKmBase = hxBase - 1;
                const int hzBase = hz_idx(i, j, eyK0);
                const int hzXmBase = hz_idx(i - 1, j, eyK0);
                for (int k = eyK0; k <= eyK1; ++k) {
                    const int off = k - eyK0;
                    const int idx = eyBase + off;
                    if (pecPlane || (pec && k == NZ - 1)) {
                        ey[idx] = 0.0;
                        continue;
                    }
                    ey[idx] = fortranCurlUpdate(
                        ey[idx], ceEy[idx],
                        hx[hxBase + off], hx[hxKmBase + off], idzh[axisCoeffIndex(k)],
                        hz[hzBase + off], hz[hzXmBase + off], idxhi);
                }
            }
        }
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
        for (int i = ezI0; i <= ezI1; ++i) {
            for (int j = ezJ0; j <= ezJ1; ++j) {
                const fdtd_real idxhi = idxh[axisCoeffIndex(i)];
                const bool pecIPlane = pec && i == NX - 1;
                const int ezBase = ez_idx(i, j, ezK0);
                const int hyBase = hy_idx(i, j, ezK0);
                const int hyXmBase = hy_idx(i - 1, j, ezK0);
                const int hxBase = hx_idx(i, j, ezK0);
                const int hxYmBase = hx_idx(i, j - 1, ezK0);
                const fdtd_real idyhj = idyh[axisCoeffIndex(j)];
                const bool pecPlane = pecIPlane || (pec && j == NY - 1);
                for (int k = ezK0; k <= ezK1; ++k) {
                    const int off = k - ezK0;
                    const int idx = ezBase + off;
                    if (pecPlane) {
                        ez[idx] = 0.0;
                        continue;
                    }
                    ez[idx] = fortranCurlUpdate(
                        ez[idx], ceEz[idx],
                        hy[hyBase + off], hy[hyXmBase + off], idxhi,
                        hx[hxBase + off], hx[hxYmBase + off], idyhj);
                }
            }
        }
#ifdef _OPENMP
        }
#endif
    }

    void applyMurE() {
        (void)useMur;
    }

    void advancePmlE() {
        if (!usePml) return;
        ++pmlElectricCalls;
        if (!cpmlBordersInitialized) initCpmlBorders();
        if (!cpmlBordersInitialized) return;

        auto advancePsi = [](fdtd_real oldPsi, fdtd_real pB,
                             fdtd_real diff, fdtd_real pC) {
            const fdtd_real damped = static_cast<fdtd_real>(pB * oldPsi);
            const fdtd_real driven = static_cast<fdtd_real>(diff * pC);
            return static_cast<fdtd_real>(damped + driven);
        };
        const int h = fieldHalo;
        const int exI0 = -h, exI1 = NX + h - 2;
        const int exJ0 = -h + 1, exJ1 = NY + h - 2;
        const int exK0 = -h + 1, exK1 = NZ + h - 2;
        const int eyI0 = -h + 1, eyI1 = NX + h - 2;
        const int eyJ0 = -h, eyJ1 = NY + h - 2;
        const int eyK0 = -h + 1, eyK1 = NZ + h - 2;
        const int ezI0 = -h + 1, ezI1 = NX + h - 2;
        const int ezJ0 = -h + 1, ezJ1 = NY + h - 2;
        const int ezK0 = -h, ezK1 = NZ + h - 2;

        for (int i = exI0; i <= exI1; ++i) {
            for (int j = exJ0; j <= exJ1; ++j) {
                const fdtd_real pceY = pmlPceY[axisCoeffIndex(j)];
                if (pceY == static_cast<fdtd_real>(0.0)) continue;
                const fdtd_real pbeY = pmlPbeY[axisCoeffIndex(j)];
                for (int k = exK0; k <= exK1; ++k) {
                    const int idx = ex_idx(i, j, k);
                    const fdtd_real diff = static_cast<fdtd_real>(
                        Hz[hz_idx(i, j, k)] - Hz[hz_idx(i, j - 1, k)]);
                    psiExy[idx] = advancePsi(psiExy[idx], pbeY, diff, pceY);
                    Ex[idx] = static_cast<fdtd_real>(
                        Ex[idx] + static_cast<fdtd_real>(CeEx[idx] * psiExy[idx]));
                }
            }
        }

        for (int i = exI0; i <= exI1; ++i) {
            for (int j = exJ0; j <= exJ1; ++j) {
                for (int k = exK0; k <= exK1; ++k) {
                    const fdtd_real pceZ = pmlPceZ[axisCoeffIndex(k)];
                    if (pceZ == static_cast<fdtd_real>(0.0)) continue;
                    const fdtd_real pbeZ = pmlPbeZ[axisCoeffIndex(k)];
                    const int idx = ex_idx(i, j, k);
                    const fdtd_real diff = static_cast<fdtd_real>(
                        Hy[hy_idx(i, j, k)] - Hy[hy_idx(i, j, k - 1)]);
                    psiExz[idx] = advancePsi(psiExz[idx], pbeZ, diff, pceZ);
                    Ex[idx] = static_cast<fdtd_real>(
                        Ex[idx] - static_cast<fdtd_real>(CeEx[idx] * psiExz[idx]));
                }
            }
        }

        for (int i = eyI0; i <= eyI1; ++i) {
            const fdtd_real pceX = pmlPceX[axisCoeffIndex(i)];
            if (pceX == static_cast<fdtd_real>(0.0)) continue;
            const fdtd_real pbeX = pmlPbeX[axisCoeffIndex(i)];
            for (int j = eyJ0; j <= eyJ1; ++j) {
                for (int k = eyK0; k <= eyK1; ++k) {
                    const int idx = ey_idx(i, j, k);
                    const fdtd_real diff = static_cast<fdtd_real>(
                        Hz[hz_idx(i, j, k)] - Hz[hz_idx(i - 1, j, k)]);
                    psiEyx[idx] = advancePsi(psiEyx[idx], pbeX, diff, pceX);
                    Ey[idx] = static_cast<fdtd_real>(
                        Ey[idx] - static_cast<fdtd_real>(CeEy[idx] * psiEyx[idx]));
                }
            }
        }

        for (int i = eyI0; i <= eyI1; ++i) {
            for (int j = eyJ0; j <= eyJ1; ++j) {
                for (int k = eyK0; k <= eyK1; ++k) {
                    const fdtd_real pceZ = pmlPceZ[axisCoeffIndex(k)];
                    if (pceZ == static_cast<fdtd_real>(0.0)) continue;
                    const fdtd_real pbeZ = pmlPbeZ[axisCoeffIndex(k)];
                    const int idx = ey_idx(i, j, k);
                    const fdtd_real diff = static_cast<fdtd_real>(
                        Hx[hx_idx(i, j, k)] - Hx[hx_idx(i, j, k - 1)]);
                    psiEyz[idx] = advancePsi(psiEyz[idx], pbeZ, diff, pceZ);
                    Ey[idx] = static_cast<fdtd_real>(
                        Ey[idx] + static_cast<fdtd_real>(CeEy[idx] * psiEyz[idx]));
                }
            }
        }

        for (int i = ezI0; i <= ezI1; ++i) {
            const fdtd_real pceX = pmlPceX[axisCoeffIndex(i)];
            if (pceX == static_cast<fdtd_real>(0.0)) continue;
            const fdtd_real pbeX = pmlPbeX[axisCoeffIndex(i)];
            for (int j = ezJ0; j <= ezJ1; ++j) {
                for (int k = ezK0; k <= ezK1; ++k) {
                    const int idx = ez_idx(i, j, k);
                    const fdtd_real diff = static_cast<fdtd_real>(
                        Hy[hy_idx(i, j, k)] - Hy[hy_idx(i - 1, j, k)]);
                    psiEzx[idx] = advancePsi(psiEzx[idx], pbeX, diff, pceX);
                    Ez[idx] = static_cast<fdtd_real>(
                        Ez[idx] + static_cast<fdtd_real>(CeEz[idx] * psiEzx[idx]));
                }
            }
        }

        for (int i = ezI0; i <= ezI1; ++i) {
            for (int j = ezJ0; j <= ezJ1; ++j) {
                const fdtd_real pceY = pmlPceY[axisCoeffIndex(j)];
                if (pceY == static_cast<fdtd_real>(0.0)) continue;
                const fdtd_real pbeY = pmlPbeY[axisCoeffIndex(j)];
                for (int k = ezK0; k <= ezK1; ++k) {
                    const int idx = ez_idx(i, j, k);
                    const fdtd_real diff = static_cast<fdtd_real>(
                        Hx[hx_idx(i, j, k)] - Hx[hx_idx(i, j - 1, k)]);
                    psiEzy[idx] = advancePsi(psiEzy[idx], pbeY, diff, pceY);
                    Ez[idx] = static_cast<fdtd_real>(
                        Ez[idx] - static_cast<fdtd_real>(CeEz[idx] * psiEzy[idx]));
                }
            }
        }
    }

    void advancePmlBodyH() {
        if (!usePml) return;
        ++pmlBodyHCalls;
        // TODO: wire translated PML-body H updates to active solver fields.
    }

    void advanceMagneticCpml() {
        if (!usePml) return;
        ++pmlMagneticCpmlCalls;
        if (!cpmlBordersInitialized) initCpmlBorders();
        if (!cpmlBordersInitialized) return;

        auto advancePsi = [](fdtd_real oldPsi, fdtd_real pB,
                             fdtd_real diff, fdtd_real pC) {
            const fdtd_real damped = static_cast<fdtd_real>(pB * oldPsi);
            const fdtd_real driven = static_cast<fdtd_real>(diff * pC);
            return static_cast<fdtd_real>(damped + driven);
        };
        const int h = fieldHalo;
        const int hxI0 = -h, hxI1 = NX + h - 1;
        const int hxJ0 = -h, hxJ1 = NY + h - 2;
        const int hxK0 = -h, hxK1 = NZ + h - 2;
        const int hyI0 = -h, hyI1 = NX + h - 2;
        const int hyJ0 = -h, hyJ1 = NY + h - 1;
        const int hyK0 = -h, hyK1 = NZ + h - 2;
        const int hzI0 = -h, hzI1 = NX + h - 2;
        const int hzJ0 = -h, hzJ1 = NY + h - 2;
        const int hzK0 = -h, hzK1 = NZ + h - 1;

        for (int i = hxI0; i <= hxI1; ++i) {
            for (int j = hxJ0; j <= hxJ1; ++j) {
                const fdtd_real pcmY = pmlPcmY[axisCoeffIndex(j)];
                if (pcmY == static_cast<fdtd_real>(0.0)) continue;
                const fdtd_real pbmY = pmlPbmY[axisCoeffIndex(j)];
                for (int k = hxK0; k <= hxK1; ++k) {
                    const int idx = hx_idx(i, j, k);
                    const fdtd_real diff = static_cast<fdtd_real>(
                        Ez[ez_idx(i, j + 1, k)] - Ez[ez_idx(i, j, k)]);
                    psiHxy[idx] = advancePsi(psiHxy[idx], pbmY, diff, pcmY);
                    Hx[idx] = static_cast<fdtd_real>(
                        Hx[idx] - static_cast<fdtd_real>(CmH[idx] * psiHxy[idx]));
                }
            }
        }

        for (int i = hxI0; i <= hxI1; ++i) {
            for (int j = hxJ0; j <= hxJ1; ++j) {
                for (int k = hxK0; k <= hxK1; ++k) {
                    const fdtd_real pcmZ = pmlPcmZ[axisCoeffIndex(k)];
                    if (pcmZ == static_cast<fdtd_real>(0.0)) continue;
                    const fdtd_real pbmZ = pmlPbmZ[axisCoeffIndex(k)];
                    const int idx = hx_idx(i, j, k);
                    const fdtd_real diff = static_cast<fdtd_real>(
                        Ey[ey_idx(i, j, k + 1)] - Ey[ey_idx(i, j, k)]);
                    psiHxz[idx] = advancePsi(psiHxz[idx], pbmZ, diff, pcmZ);
                    Hx[idx] = static_cast<fdtd_real>(
                        Hx[idx] + static_cast<fdtd_real>(CmH[idx] * psiHxz[idx]));
                }
            }
        }

        for (int i = hyI0; i <= hyI1; ++i) {
            const fdtd_real pcmX = pmlPcmX[axisCoeffIndex(i)];
            if (pcmX == static_cast<fdtd_real>(0.0)) continue;
            const fdtd_real pbmX = pmlPbmX[axisCoeffIndex(i)];
            for (int j = hyJ0; j <= hyJ1; ++j) {
                for (int k = hyK0; k <= hyK1; ++k) {
                    const int idx = hy_idx(i, j, k);
                    const fdtd_real diff = static_cast<fdtd_real>(
                        Ez[ez_idx(i + 1, j, k)] - Ez[ez_idx(i, j, k)]);
                    psiHyx[idx] = advancePsi(psiHyx[idx], pbmX, diff, pcmX);
                    Hy[idx] = static_cast<fdtd_real>(
                        Hy[idx] + static_cast<fdtd_real>(CmH[idx] * psiHyx[idx]));
                }
            }
        }

        for (int i = hyI0; i <= hyI1; ++i) {
            for (int j = hyJ0; j <= hyJ1; ++j) {
                for (int k = hyK0; k <= hyK1; ++k) {
                    const fdtd_real pcmZ = pmlPcmZ[axisCoeffIndex(k)];
                    if (pcmZ == static_cast<fdtd_real>(0.0)) continue;
                    const fdtd_real pbmZ = pmlPbmZ[axisCoeffIndex(k)];
                    const int idx = hy_idx(i, j, k);
                    const fdtd_real diff = static_cast<fdtd_real>(
                        Ex[ex_idx(i, j, k + 1)] - Ex[ex_idx(i, j, k)]);
                    psiHyz[idx] = advancePsi(psiHyz[idx], pbmZ, diff, pcmZ);
                    Hy[idx] = static_cast<fdtd_real>(
                        Hy[idx] - static_cast<fdtd_real>(CmH[idx] * psiHyz[idx]));
                }
            }
        }

        for (int i = hzI0; i <= hzI1; ++i) {
            const fdtd_real pcmX = pmlPcmX[axisCoeffIndex(i)];
            if (pcmX == static_cast<fdtd_real>(0.0)) continue;
            const fdtd_real pbmX = pmlPbmX[axisCoeffIndex(i)];
            for (int j = hzJ0; j <= hzJ1; ++j) {
                for (int k = hzK0; k <= hzK1; ++k) {
                    const int idx = hz_idx(i, j, k);
                    const fdtd_real diff = static_cast<fdtd_real>(
                        Ey[ey_idx(i + 1, j, k)] - Ey[ey_idx(i, j, k)]);
                    psiHzx[idx] = advancePsi(psiHzx[idx], pbmX, diff, pcmX);
                    Hz[idx] = static_cast<fdtd_real>(
                        Hz[idx] - static_cast<fdtd_real>(CmH[idx] * psiHzx[idx]));
                }
            }
        }

        for (int i = hzI0; i <= hzI1; ++i) {
            for (int j = hzJ0; j <= hzJ1; ++j) {
                const fdtd_real pcmY = pmlPcmY[axisCoeffIndex(j)];
                if (pcmY == static_cast<fdtd_real>(0.0)) continue;
                const fdtd_real pbmY = pmlPbmY[axisCoeffIndex(j)];
                for (int k = hzK0; k <= hzK1; ++k) {
                    const int idx = hz_idx(i, j, k);
                    const fdtd_real diff = static_cast<fdtd_real>(
                        Ex[ex_idx(i, j + 1, k)] - Ex[ex_idx(i, j, k)]);
                    psiHzy[idx] = advancePsi(psiHzy[idx], pbmY, diff, pcmY);
                    Hz[idx] = static_cast<fdtd_real>(
                        Hz[idx] + static_cast<fdtd_real>(CmH[idx] * psiHzy[idx]));
                }
            }
        }
    }

    void applyPecE() {
        if (hasAnyPecMask) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (long long idx = 0; idx < static_cast<long long>(pecExMask.size()); ++idx) {
                if (pecExMask[static_cast<size_t>(idx)]) Ex[static_cast<size_t>(idx)] = 0.0;
            }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (long long idx = 0; idx < static_cast<long long>(pecEyMask.size()); ++idx) {
                if (pecEyMask[static_cast<size_t>(idx)]) Ey[static_cast<size_t>(idx)] = 0.0;
            }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (long long idx = 0; idx < static_cast<long long>(pecEzMask.size()); ++idx) {
                if (pecEzMask[static_cast<size_t>(idx)]) Ez[static_cast<size_t>(idx)] = 0.0;
            }
        }

        if (!usePec) return;

        // Fortran field sweeps start at index 1; index 0 is outside the
        // probe-visible domain. Clamp only the ghost planes here.
        for (int i = 0; i < NX; ++i) {
            for (int k = -1; k <= NZ + 1; ++k) {
                Ex[ex_idx(i, -1, k)] = 0.0;
                Ex[ex_idx(i, NY + 1, k)] = 0.0;
            }
            for (int j = -1; j <= NY + 1; ++j) {
                Ex[ex_idx(i, j, -1)] = 0.0;
                Ex[ex_idx(i, j, NZ + 1)] = 0.0;
            }
        }
	        for (int j = -1; j <= NY + 1; ++j) {
	            Ex[ex_idx(-1, j, -1)] = 0.0;
	        }
	        for (int k = -1; k <= NZ + 1; ++k) {
	            Ex[ex_idx(-1, -1, k)] = 0.0;
	        }
	        for (int j = 0; j < NY; ++j) {
	            for (int k = -1; k <= NZ + 1; ++k) {
	                Ey[ey_idx(-1, j, k)] = 0.0;
	                Ey[ey_idx(NX + 1, j, k)] = 0.0;
	            }
	        }
	        for (int k = -1; k <= NZ + 1; ++k) {
	            Ey[ey_idx(-1, -1, k)] = 0.0;
	        }
        for (int i = -1; i <= NX + 1; ++i) {
            for (int j = 0; j < NY; ++j) {
                Ey[ey_idx(i, j, -1)] = 0.0;
                Ey[ey_idx(i, j, NZ + 1)] = 0.0;
            }
        }
        for (int i = -1; i <= NX + 1; ++i) {
            Ey[ey_idx(i, -1, -1)] = 0.0;
        }

        for (int j = -1; j <= NY + 1; ++j) {
            for (int k = 0; k < NZ; ++k) {
                Ez[ez_idx(-1, j, k)] = 0.0;
                Ez[ez_idx(NX + 1, j, k)] = 0.0;
            }
        }
        for (int i = -1; i <= NX + 1; ++i) {
            for (int k = 0; k < NZ; ++k) {
                Ez[ez_idx(i, -1, k)] = 0.0;
                Ez[ez_idx(i, NY + 1, k)] = 0.0;
            }
        }
        for (int j = -1; j <= NY + 1; ++j) {
            Ez[ez_idx(-1, j, -1)] = 0.0;
        }
        for (int i = -1; i <= NX + 1; ++i) {
            Ez[ez_idx(i, -1, -1)] = 0.0;
        }
    }

    void applyPecH() {
        (void)usePec;
    }

    bool hasPlaneWaveSource() const {
        return std::any_of(sources.begin(), sources.end(), [](const source_t& src) {
            return src.type == "planewave";
        });
    }

    void applyMurH() {
        if (!useMur) return;
        auto mur_face = [](fdtd_real& bnd, fdtd_real int_now, fdtd_real past_int,
                           fdtd_real past_bnd, fdtd_real cab) {
            bnd = fortranMurFace(int_now, past_int, past_bnd, cab);
        };

        // Back (x min): Fortran MURc uses sweepXI-1, one plane below the H sweep.
        if (murBack) {
            for (int j = -1; j <= NY - 1; ++j) {
                for (int k = -1; k <= NZ - 2; ++k) {
                    const size_t p = static_cast<size_t>((j + 1) * NZ + (k + 1));
                    const int idx0 = hy_idx(-2, j, k);
                    const int idx1 = hy_idx(-1, j, k);
                    mur_face(Hy[idx0], Hy[idx1], murPastHyBackInt[p], murPastHyBack[p], backCab1);
                }
            }
            for (int j = -1; j <= NY - 2; ++j) {
                for (int k = -1; k <= NZ - 1; ++k) {
                    const size_t p = static_cast<size_t>((j + 1) * (NZ + 1) + (k + 1));
                    const int idx0 = hz_idx(-2, j, k);
                    const int idx1 = hz_idx(-1, j, k);
                    mur_face(Hz[idx0], Hz[idx1], murPastHzBackInt[p], murPastHzBack[p], backCab1);
                }
            }
        }
        // Front (x max): Fortran first-order Mur writes MURc%XE, one plane above the H sweep.
        if (murFront) {
            for (int j = -1; j <= NY - 1; ++j) {
                for (int k = -1; k <= NZ - 2; ++k) {
                    const size_t p = static_cast<size_t>((j + 1) * NZ + (k + 1));
                    const int idxN = hy_idx(NX - 1, j, k);
                    const int idxI = hy_idx(NX - 2, j, k);
                    mur_face(Hy[idxN], Hy[idxI], murPastHyFrontInt[p], murPastHyFront[p], frontCab1);
                }
            }
            for (int j = -1; j <= NY - 2; ++j) {
                for (int k = -1; k <= NZ - 1; ++k) {
                    const size_t p = static_cast<size_t>((j + 1) * (NZ + 1) + (k + 1));
                    const int idxN = hz_idx(NX - 1, j, k);
                    const int idxI = hz_idx(NX - 2, j, k);
                    mur_face(Hz[idxN], Hz[idxI], murPastHzFrontInt[p], murPastHzFront[p], frontCab1);
                }
            }
        }
        // Left (y min): Fortran MURc uses sweepYI-1, one plane below the H sweep.
        if (murLeft) {
            for (int i = -1; i <= NX - 1; ++i) {
                for (int k = -1; k <= NZ - 2; ++k) {
                    const size_t p = static_cast<size_t>((i + 1) * NZ + (k + 1));
                    const int idx0 = hx_idx(i, -2, k);
                    const int idx1 = hx_idx(i, -1, k);
                    mur_face(Hx[idx0], Hx[idx1], murPastHxLeftInt[p], murPastHxLeft[p], leftCab1);
                }
            }
            for (int i = -1; i <= NX - 2; ++i) {
                for (int k = -1; k <= NZ - 1; ++k) {
                    const size_t p = static_cast<size_t>((i + 1) * (NZ + 1) + (k + 1));
                    const int idx0 = hz_idx(i, -2, k);
                    const int idx1 = hz_idx(i, -1, k);
                    mur_face(Hz[idx0], Hz[idx1], murPastHzLeftInt[p], murPastHzLeft[p], leftCab1);
                }
            }
        }
        // Right (y max): Fortran first-order Mur writes MURc%YE, one plane above the H sweep.
        if (murRight) {
            for (int i = -1; i <= NX - 1; ++i) {
                for (int k = -1; k <= NZ - 2; ++k) {
                    const size_t p = static_cast<size_t>((i + 1) * NZ + (k + 1));
                    const int idxN = hx_idx(i, NY - 1, k);
                    const int idxI = hx_idx(i, NY - 2, k);
                    mur_face(Hx[idxN], Hx[idxI], murPastHxRightInt[p], murPastHxRight[p], rightCab1);
                }
            }
            for (int i = -1; i <= NX - 2; ++i) {
                for (int k = -1; k <= NZ - 1; ++k) {
                    const size_t p = static_cast<size_t>((i + 1) * (NZ + 1) + (k + 1));
                    const int idxN = hz_idx(i, NY - 1, k);
                    const int idxI = hz_idx(i, NY - 2, k);
                    mur_face(Hz[idxN], Hz[idxI], murPastHzRightInt[p], murPastHzRight[p], rightCab1);
                }
            }
        }
        // Down (z min): Fortran MURc uses sweepZI-1, one plane below the H sweep.
        if (murDown) {
            for (int i = -1; i <= NX - 2; ++i) {
                for (int j = -1; j <= NY - 1; ++j) {
                    const size_t p = static_cast<size_t>((i + 1) * (NY + 1) + (j + 1));
                    const int idx0 = hy_idx(i, j, -2);
                    const int idx1 = hy_idx(i, j, -1);
                    mur_face(Hy[idx0], Hy[idx1], murPastHyDownInt[p], murPastHyDown[p], downCab1);
                }
            }
            for (int i = -1; i <= NX - 1; ++i) {
                for (int j = -1; j <= NY - 2; ++j) {
                    const size_t p = static_cast<size_t>((i + 1) * NY + (j + 1));
                    const int idx0 = hx_idx(i, j, -2);
                    const int idx1 = hx_idx(i, j, -1);
                    mur_face(Hx[idx0], Hx[idx1], murPastHxDownInt[p], murPastHxDown[p], downCab1);
                }
            }
        }
        // Up (z max): Fortran first-order Mur writes MURc%ZE, one plane above the H sweep.
        if (murUp) {
            for (int i = -1; i <= NX - 2; ++i) {
                for (int j = -1; j <= NY - 1; ++j) {
                    const size_t p = static_cast<size_t>((i + 1) * (NY + 1) + (j + 1));
                    const int idxN = hy_idx(i, j, NZ - 1);
                    const int idxI = hy_idx(i, j, NZ - 2);
                    mur_face(Hy[idxN], Hy[idxI], murPastHyUpInt[p], murPastHyUp[p], upCab1);
                }
            }
            for (int i = -1; i <= NX - 1; ++i) {
                for (int j = -1; j <= NY - 2; ++j) {
                    const size_t p = static_cast<size_t>((i + 1) * NY + (j + 1));
                    const int idxN = hx_idx(i, j, NZ - 1);
                    const int idxI = hx_idx(i, j, NZ - 2);
                    mur_face(Hx[idxN], Hx[idxI], murPastHxUpInt[p], murPastHxUp[p], upCab1);
                }
            }
        }

        // Store past fields (Fortran AdvanceMagneticMUR tail).
        if (murBack || murFront) {
            for (int j = -1; j <= NY - 1; ++j) {
                for (int k = -1; k <= NZ - 2; ++k) {
                    const size_t p = static_cast<size_t>((j + 1) * NZ + (k + 1));
                    if (murBack) {
                        murPastHyBack[p] = Hy[hy_idx(-2, j, k)];
                        murPastHyBackInt[p] = Hy[hy_idx(-1, j, k)];
                    }
                    if (murFront) {
                        murPastHyFront[p] = Hy[hy_idx(NX - 1, j, k)];
                        murPastHyFrontInt[p] = Hy[hy_idx(NX - 2, j, k)];
                    }
                }
            }
            for (int j = -1; j <= NY - 2; ++j) {
                for (int k = -1; k <= NZ - 1; ++k) {
                    const size_t p = static_cast<size_t>((j + 1) * (NZ + 1) + (k + 1));
                    if (murBack) {
                        murPastHzBack[p] = Hz[hz_idx(-2, j, k)];
                        murPastHzBackInt[p] = Hz[hz_idx(-1, j, k)];
                    }
                    if (murFront) {
                        murPastHzFront[p] = Hz[hz_idx(NX - 1, j, k)];
                        murPastHzFrontInt[p] = Hz[hz_idx(NX - 2, j, k)];
                    }
                }
            }
        }
        if (murLeft || murRight) {
            for (int i = -1; i <= NX - 1; ++i) {
                for (int k = -1; k <= NZ - 2; ++k) {
                    const size_t p = static_cast<size_t>((i + 1) * NZ + (k + 1));
                    if (murLeft) {
                        murPastHxLeft[p] = Hx[hx_idx(i, -2, k)];
                        murPastHxLeftInt[p] = Hx[hx_idx(i, -1, k)];
                    }
                    if (murRight) {
                        murPastHxRight[p] = Hx[hx_idx(i, NY - 1, k)];
                        murPastHxRightInt[p] = Hx[hx_idx(i, NY - 2, k)];
                    }
                }
            }
            for (int i = -1; i <= NX - 2; ++i) {
                for (int k = -1; k <= NZ - 1; ++k) {
                    const size_t p = static_cast<size_t>((i + 1) * (NZ + 1) + (k + 1));
                    if (murLeft) {
                        murPastHzLeft[p] = Hz[hz_idx(i, -2, k)];
                        murPastHzLeftInt[p] = Hz[hz_idx(i, -1, k)];
                    }
                    if (murRight) {
                        murPastHzRight[p] = Hz[hz_idx(i, NY - 1, k)];
                        murPastHzRightInt[p] = Hz[hz_idx(i, NY - 2, k)];
                    }
                }
            }
        }
        if (murDown || murUp) {
            for (int i = -1; i <= NX - 2; ++i) {
                for (int j = -1; j <= NY - 1; ++j) {
                    const size_t p = static_cast<size_t>((i + 1) * (NY + 1) + (j + 1));
                    if (murDown) {
                        murPastHyDown[p] = Hy[hy_idx(i, j, -2)];
                        murPastHyDownInt[p] = Hy[hy_idx(i, j, -1)];
                    }
                    if (murUp) {
                        murPastHyUp[p] = Hy[hy_idx(i, j, NZ - 1)];
                        murPastHyUpInt[p] = Hy[hy_idx(i, j, NZ - 2)];
                    }
                }
            }
            for (int i = -1; i <= NX - 1; ++i) {
                for (int j = -1; j <= NY - 2; ++j) {
                    const size_t p = static_cast<size_t>((i + 1) * NY + (j + 1));
                    if (murDown) {
                        murPastHxDown[p] = Hx[hx_idx(i, j, -2)];
                        murPastHxDownInt[p] = Hx[hx_idx(i, j, -1)];
                    }
                    if (murUp) {
                        murPastHxUp[p] = Hx[hx_idx(i, j, NZ - 1)];
                        murPastHxUpInt[p] = Hx[hx_idx(i, j, NZ - 2)];
                    }
                }
            }
        }
    }

    void advanceH() {
        if (Idxe.empty() || Idye.empty() || Idze.empty()) {
            initGridInverses();
        }
        const fdtd_real* SEMBA_RESTRICT ex = Ex.data();
        const fdtd_real* SEMBA_RESTRICT ey = Ey.data();
        const fdtd_real* SEMBA_RESTRICT ez = Ez.data();
        fdtd_real* SEMBA_RESTRICT hx = Hx.data();
        fdtd_real* SEMBA_RESTRICT hy = Hy.data();
        fdtd_real* SEMBA_RESTRICT hz = Hz.data();
        const fdtd_real* SEMBA_RESTRICT cmH = CmH.data();
        const fdtd_real* SEMBA_RESTRICT idxe = Idxe.data();
        const fdtd_real* SEMBA_RESTRICT idye = Idye.data();
        const fdtd_real* SEMBA_RESTRICT idze = Idze.data();
        const bool pec = usePec;
        const bool pml = usePml;
        const int h = fieldHalo;
        int hxI0 = pml ? -h : -1;
        int hxI1 = pml ? NX + h - 1 : NX - 1;
        int hxJ0 = pml ? -h : -1;
        int hxJ1 = pml ? NY + h - 2 : NY - 2;
        int hxK0 = pml ? -h : -1;
        int hxK1 = pml ? NZ + h - 2 : NZ - 2;
        int hyI0 = pml ? -h : -1;
        int hyI1 = pml ? NX + h - 2 : NX - 2;
        int hyJ0 = pml ? -h : -1;
        int hyJ1 = pml ? NY + h - 1 : NY - 1;
        int hyK0 = pml ? -h : -1;
        int hyK1 = pml ? NZ + h - 2 : NZ - 2;
        int hzI0 = pml ? -h : -1;
        int hzI1 = pml ? NX + h - 2 : NX - 2;
        int hzJ0 = pml ? -h : -1;
        int hzJ1 = pml ? NY + h - 2 : NY - 2;
        int hzK0 = pml ? -h : -1;
        int hzK1 = pml ? NZ + h - 1 : NZ - 1;
        clampMpiComponentAxisRange(3, 1, hxI0, hxI1);
        clampMpiComponentAxisRange(3, 2, hxJ0, hxJ1);
        clampMpiComponentAxisRange(3, 3, hxK0, hxK1);
        clampMpiComponentAxisRange(4, 1, hyI0, hyI1);
        clampMpiComponentAxisRange(4, 2, hyJ0, hyJ1);
        clampMpiComponentAxisRange(4, 3, hyK0, hyK1);
        clampMpiComponentAxisRange(5, 1, hzI0, hzI1);
        clampMpiComponentAxisRange(5, 2, hzJ0, hzJ1);
        clampMpiComponentAxisRange(5, 3, hzK0, hzK1);
#ifdef _OPENMP
#pragma omp parallel
        {
#pragma omp for collapse(2) schedule(static)
#endif
        for (int i = hxI0; i <= hxI1; ++i) {
            for (int j = hxJ0; j <= hxJ1; ++j) {
                const bool pecPlane = pec && i == NX - 1;
                const int hxBase = hx_idx(i, j, hxK0);
                const int eyBase = ey_idx(i, j, hxK0);
                const int ezBase = ez_idx(i, j, hxK0);
                const int ezYpBase = ez_idx(i, j + 1, hxK0);
                const fdtd_real idyej = idye[axisCoeffIndex(j)];
                for (int k = hxK0; k <= hxK1; ++k) {
                    const int off = k - hxK0;
                    const int idx = hxBase + off;
                    if (pecPlane) {
                        hx[idx] = 0.0;
                        continue;
                    }
                    hx[idx] = fortranCurlUpdate(
                        hx[idx], cmH[idx],
                        ey[eyBase + off + 1], ey[eyBase + off], idze[axisCoeffIndex(k)],
                        ez[ezYpBase + off], ez[ezBase + off], idyej);
                }
            }
        }
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
        for (int i = hyI0; i <= hyI1; ++i) {
            for (int j = hyJ0; j <= hyJ1; ++j) {
                const fdtd_real idxei = idxe[axisCoeffIndex(i)];
                const int hyBase = hy_idx(i, j, hyK0);
                const int ezBase = ez_idx(i, j, hyK0);
                const int ezXpBase = ez_idx(i + 1, j, hyK0);
                const int exBase = ex_idx(i, j, hyK0);
                const bool pecPlane = pec && j == NY - 1;
                for (int k = hyK0; k <= hyK1; ++k) {
                    const int off = k - hyK0;
                    const int idx = hyBase + off;
                    if (pecPlane) {
                        hy[idx] = 0.0;
                        continue;
                    }
                    hy[idx] = fortranCurlUpdate(
                        hy[idx], cmH[idx],
                        ez[ezXpBase + off], ez[ezBase + off], idxei,
                        ex[exBase + off + 1], ex[exBase + off], idze[axisCoeffIndex(k)]);
                }
            }
        }
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
        for (int i = hzI0; i <= hzI1; ++i) {
            for (int j = hzJ0; j <= hzJ1; ++j) {
                const fdtd_real idxei = idxe[axisCoeffIndex(i)];
                const int hzBase = hz_idx(i, j, hzK0);
                const int exBase = ex_idx(i, j, hzK0);
                const int exYpBase = ex_idx(i, j + 1, hzK0);
                const int eyBase = ey_idx(i, j, hzK0);
                const int eyXpBase = ey_idx(i + 1, j, hzK0);
                const fdtd_real idyej = idye[axisCoeffIndex(j)];
                for (int k = hzK0; k <= hzK1; ++k) {
                    const int off = k - hzK0;
                    const int idx = hzBase + off;
                    if (pec && k == NZ - 1) {
                        hz[idx] = 0.0;
                        continue;
                    }
                    hz[idx] = fortranCurlUpdate(
                        hz[idx], cmH[idx],
                        ex[exYpBase + off], ex[exBase + off], idyej,
                        ey[eyXpBase + off], ey[eyBase + off], idxei);
                }
            }
        }
#ifdef _OPENMP
        }
#endif
    }

    void minusCloneMagneticPmc() {
        if (pmcDown) {
            for (int i = -1; i <= NX - 1; ++i)
                for (int j = -1; j <= NY - 2; ++j)
                    Hx[hx_idx(i, j, -2)] = -Hx[hx_idx(i, j, -1)];
            for (int i = -1; i <= NX - 2; ++i)
                for (int j = -1; j <= NY - 1; ++j)
                    Hy[hy_idx(i, j, -2)] = -Hy[hy_idx(i, j, -1)];
        }
        if (pmcUp) {
            for (int i = -1; i <= NX - 1; ++i)
                for (int j = -1; j <= NY - 2; ++j)
                    Hx[hx_idx(i, j, NZ - 1)] = -Hx[hx_idx(i, j, NZ - 2)];
            for (int i = -1; i <= NX - 2; ++i)
                for (int j = -1; j <= NY - 1; ++j)
                    Hy[hy_idx(i, j, NZ - 1)] = -Hy[hy_idx(i, j, NZ - 2)];
        }
        if (pmcBack) {
            for (int j = -1; j <= NY - 1; ++j)
                for (int k = -1; k <= NZ - 2; ++k)
                    Hy[hy_idx(-2, j, k)] = -Hy[hy_idx(-1, j, k)];
            for (int j = -1; j <= NY - 2; ++j)
                for (int k = -1; k <= NZ - 1; ++k)
                    Hz[hz_idx(-2, j, k)] = -Hz[hz_idx(-1, j, k)];
        }
        if (pmcFront) {
            for (int j = -1; j <= NY - 1; ++j)
                for (int k = -1; k <= NZ - 2; ++k)
                    Hy[hy_idx(NX - 1, j, k)] = -Hy[hy_idx(NX - 2, j, k)];
            for (int j = -1; j <= NY - 2; ++j)
                for (int k = -1; k <= NZ - 1; ++k)
                    Hz[hz_idx(NX - 1, j, k)] = -Hz[hz_idx(NX - 2, j, k)];
        }
        if (pmcLeft) {
            for (int i = -1; i <= NX - 1; ++i)
                for (int k = -1; k <= NZ - 2; ++k)
                    Hx[hx_idx(i, -2, k)] = -Hx[hx_idx(i, -1, k)];
            for (int i = -1; i <= NX - 2; ++i)
                for (int k = -1; k <= NZ - 1; ++k)
                    Hz[hz_idx(i, -2, k)] = -Hz[hz_idx(i, -1, k)];
        }
        if (pmcRight) {
            for (int i = -1; i <= NX - 1; ++i)
                for (int k = -1; k <= NZ - 2; ++k)
                    Hx[hx_idx(i, NY - 1, k)] = -Hx[hx_idx(i, NY - 2, k)];
            for (int i = -1; i <= NX - 2; ++i)
                for (int k = -1; k <= NZ - 1; ++k)
                    Hz[hz_idx(i, NY - 1, k)] = -Hz[hz_idx(i, NY - 2, k)];
        }
    }

    void cloneMagneticPeriodic() {
        // Fortran rebases field arrays to zero before applying periodic
        // clones. In this standalone solver, that maps the lower periodic
        // clone onto both the outer ghost plane and the first swept H plane.
        if (periodicDown) {
            for (int i = -1; i <= NX - 1; ++i)
                for (int j = -1; j <= NY - 2; ++j)
                    Hx[hx_idx(i, j, -2)] = Hx[hx_idx(i, j, NZ - 2)];
            for (int i = -1; i <= NX - 2; ++i)
                for (int j = -1; j <= NY - 1; ++j)
                    Hy[hy_idx(i, j, -2)] = Hy[hy_idx(i, j, NZ - 2)];
            for (int i = -1; i <= NX - 1; ++i)
                for (int j = -1; j <= NY - 2; ++j)
                    Hx[hx_idx(i, j, -1)] = Hx[hx_idx(i, j, NZ - 2)];
            for (int i = -1; i <= NX - 2; ++i)
                for (int j = -1; j <= NY - 1; ++j)
                    Hy[hy_idx(i, j, -1)] = Hy[hy_idx(i, j, NZ - 2)];
        }
        if (periodicUp) {
            for (int i = -1; i <= NX - 1; ++i)
                for (int j = -1; j <= NY - 2; ++j)
                    Hx[hx_idx(i, j, NZ - 1)] = Hx[hx_idx(i, j, 0)];
            for (int i = -1; i <= NX - 2; ++i)
                for (int j = -1; j <= NY - 1; ++j)
                    Hy[hy_idx(i, j, NZ - 1)] = Hy[hy_idx(i, j, 0)];
        }
        if (periodicBack) {
            for (int j = -1; j <= NY - 1; ++j)
                for (int k = -1; k <= NZ - 2; ++k)
                    Hy[hy_idx(-2, j, k)] = Hy[hy_idx(NX - 2, j, k)];
            for (int j = -1; j <= NY - 2; ++j)
                for (int k = -1; k <= NZ - 1; ++k)
                    Hz[hz_idx(-2, j, k)] = Hz[hz_idx(NX - 2, j, k)];
            for (int j = -1; j <= NY - 1; ++j)
                for (int k = -1; k <= NZ - 2; ++k)
                    Hy[hy_idx(-1, j, k)] = Hy[hy_idx(NX - 2, j, k)];
            for (int j = -1; j <= NY - 2; ++j)
                for (int k = -1; k <= NZ - 1; ++k)
                    Hz[hz_idx(-1, j, k)] = Hz[hz_idx(NX - 2, j, k)];
        }
        if (periodicFront) {
            for (int j = -1; j <= NY - 1; ++j)
                for (int k = -1; k <= NZ - 2; ++k)
                    Hy[hy_idx(NX - 1, j, k)] = Hy[hy_idx(0, j, k)];
            for (int j = -1; j <= NY - 2; ++j)
                for (int k = -1; k <= NZ - 1; ++k)
                    Hz[hz_idx(NX - 1, j, k)] = Hz[hz_idx(0, j, k)];
        }
        if (periodicLeft) {
            for (int i = -1; i <= NX - 1; ++i)
                for (int k = -1; k <= NZ - 2; ++k)
                    Hx[hx_idx(i, -2, k)] = Hx[hx_idx(i, NY - 2, k)];
            for (int i = -1; i <= NX - 2; ++i)
                for (int k = -1; k <= NZ - 1; ++k)
                    Hz[hz_idx(i, -2, k)] = Hz[hz_idx(i, NY - 2, k)];
            for (int i = -1; i <= NX - 1; ++i)
                for (int k = -1; k <= NZ - 2; ++k)
                    Hx[hx_idx(i, -1, k)] = Hx[hx_idx(i, NY - 2, k)];
            for (int i = -1; i <= NX - 2; ++i)
                for (int k = -1; k <= NZ - 1; ++k)
                    Hz[hz_idx(i, -1, k)] = Hz[hz_idx(i, NY - 2, k)];
        }
        if (periodicRight) {
            for (int i = -1; i <= NX - 1; ++i)
                for (int k = -1; k <= NZ - 2; ++k)
                    Hx[hx_idx(i, NY - 1, k)] = Hx[hx_idx(i, 0, k)];
            for (int i = -1; i <= NX - 2; ++i)
                for (int k = -1; k <= NZ - 1; ++k)
                    Hz[hz_idx(i, NY - 1, k)] = Hz[hz_idx(i, 0, k)];
        }
    }

    // Full-domain stencils (timestepping.F90 / timestepping.cpp) for Mur absorption tests.
    void advanceE_fortran() {
        advanceE();
    }

    void advanceH_fortran() {
        advanceH();
    }

    // Fortran timestepping.F90 flushPlanewaveOff L2110-2131.
    void flushPlanewaveOff() {
        if (planewave_switched_off || planeWaves.empty()) return;
        still_planewave_time = still_planewave_time && !planeWaves.empty();
#ifdef CompileWithMPI
        if (mpiEnabled) {
            const int localStill = still_planewave_time ? 1 : 0;
            int globalStill = 0;
            MPI_Allreduce(&localStill, &globalStill, 1, MPI_INT, MPI_LOR, SUBCOMM_MPI);
            still_planewave_time = globalStill != 0;
        }
#endif
        if (!still_planewave_time) {
            planewave_switched_off = true;
        }
    }

    void advancePlaneWaveE() {
        if (planewave_switched_off || planeWaves.empty()) return;
        still_planewave_time = false;
        const int XI = 1, XE = NX, YI = 1, YE = NY, ZI = 1, ZE = NZ;
        const fdtd_real G2_1 = static_cast<fdtd_real>(
            dt / static_cast<double>(static_cast<fdtd_real>(eps0)));
        auto addPw = [](fdtd_real& value, fdtd_real coeff, fdtd_real inc, fdtd_real inv) {
            value += coeff * inc * inv;
        };
        auto subPw = [](fdtd_real& value, fdtd_real coeff, fdtd_real inc, fdtd_real inv) {
            value -= coeff * inc * inv;
        };
        for (int pwIdx = 0; pwIdx < (int)planeWaves.size(); ++pwIdx) {
            const auto& pw = planeWaves[pwIdx];
            if (pw.iluminaTr && mpiOwnsPlaneWaveFace(1, pw.esqx1)) {
                int i = std::max(XI, pw.esqx1);
                fdtd_real id = static_cast<fdtd_real>(idxhPlanewave1(i));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2 - 1); ++k) {
                    for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2); ++j) {
                        const fdtd_real inc = computeIncid(pwIdx, 4, currentTime, i - 1, j, k);
                        if (!mpiOwnsComponentCoordinate(2, i - 1, j - 1, k - 1)) continue;
                        subPw(Ez[ez_idx(i - 1, j - 1, k - 1)], G2_1, inc, id);
                    }
                }
                i = std::max(XI, pw.esqx1);
                id = static_cast<fdtd_real>(idxhPlanewave1(i));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2 - 1); ++j) {
                        const fdtd_real inc = computeIncid(pwIdx, 5, currentTime, i - 1, j, k);
                        if (!mpiOwnsComponentCoordinate(1, i - 1, j - 1, k - 1)) continue;
                        addPw(Ey[ey_idx(i - 1, j - 1, k - 1)], G2_1, inc, id);
                    }
                }
            }
            if (pw.iluminaFr && mpiOwnsPlaneWaveFace(1, pw.esqx2)) {
                int i = std::min(XE, pw.esqx2);
                fdtd_real id = static_cast<fdtd_real>(idxhPlanewave1(i));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2 - 1); ++k) {
                    for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2); ++j) {
                        const fdtd_real inc = computeIncid(pwIdx, 4, currentTime, i, j, k);
                        if (!mpiOwnsComponentCoordinate(2, i - 1, j - 1, k - 1)) continue;
                        addPw(Ez[ez_idx(i - 1, j - 1, k - 1)], G2_1, inc, id);
                    }
                }
                i = std::min(XE, pw.esqx2);
                id = static_cast<fdtd_real>(idxhPlanewave1(i));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2 - 1); ++j) {
                        const fdtd_real inc = computeIncid(pwIdx, 5, currentTime, i, j, k);
                        if (!mpiOwnsComponentCoordinate(1, i - 1, j - 1, k - 1)) continue;
                        subPw(Ey[ey_idx(i - 1, j - 1, k - 1)], G2_1, inc, id);
                    }
                }
            }
            if (pw.iluminaIz && mpiOwnsPlaneWaveFace(2, pw.esqy1)) {
                int j = std::max(YI, pw.esqy1);
                fdtd_real id = static_cast<fdtd_real>(idyhPlanewave1(j));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2 - 1); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 5, currentTime, i, j - 1, k);
                        if (!mpiOwnsComponentCoordinate(0, i - 1, j - 1, k - 1)) continue;
                        subPw(Ex[ex_idx(i - 1, j - 1, k - 1)], G2_1, inc, id);
                    }
                }
                j = std::max(YI, pw.esqy1);
                id = static_cast<fdtd_real>(idyhPlanewave1(j));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2 - 1); ++k) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 3, currentTime, i, j - 1, k);
                        if (!mpiOwnsComponentCoordinate(2, i - 1, j - 1, k - 1)) continue;
                        addPw(Ez[ez_idx(i - 1, j - 1, k - 1)], G2_1, inc, id);
                    }
                }
            }
            if (pw.iluminaDe && mpiOwnsPlaneWaveFace(2, pw.esqy2)) {
                int j = std::min(YE, pw.esqy2);
                fdtd_real id = static_cast<fdtd_real>(idyhPlanewave1(j));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2 - 1); ++k) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 3, currentTime, i, j, k);
                        if (!mpiOwnsComponentCoordinate(2, i - 1, j - 1, k - 1)) continue;
                        subPw(Ez[ez_idx(i - 1, j - 1, k - 1)], G2_1, inc, id);
                    }
                }
                j = std::min(YE, pw.esqy2);
                id = static_cast<fdtd_real>(idyhPlanewave1(j));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2 - 1); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 5, currentTime, i, j, k);
                        if (!mpiOwnsComponentCoordinate(0, i - 1, j - 1, k - 1)) continue;
                        addPw(Ex[ex_idx(i - 1, j - 1, k - 1)], G2_1, inc, id);
                    }
                }
            }
            if (pw.iluminaAb && mpiOwnsPlaneWaveFace(3, pw.esqz1)) {
                int k = std::max(ZI, pw.esqz1);
                fdtd_real id = static_cast<fdtd_real>(idzhPlanewave1(k));
                for (int j = std::max(0, pw.esqy1); j <= std::min(NY, pw.esqy2); ++j) {
                    for (int i = std::max(0, pw.esqx1); i <= std::min(NX - 1, pw.esqx2 - 1); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 4, currentTime, i, j, k - 1);
                        if (!mpiOwnsComponentCoordinate(0, i - 1, j - 1, k - 1)) continue;
                        const int exIndex = ex_idx(i - 1, j - 1, k - 1);
                        addPw(Ex[exIndex], G2_1, inc, id);
                    }
                }
                k = std::max(ZI, pw.esqz1);
                id = static_cast<fdtd_real>(idzhPlanewave1(k));
                for (int j = std::max(0, pw.esqy1); j <= std::min(NY - 1, pw.esqy2 - 1); ++j) {
                    for (int i = std::max(0, pw.esqx1); i <= std::min(NX, pw.esqx2); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 3, currentTime, i, j, k - 1);
                        if (!mpiOwnsComponentCoordinate(1, i - 1, j - 1, k - 1)) continue;
                        subPw(Ey[ey_idx(i - 1, j - 1, k - 1)], G2_1, inc, id);
                    }
                }
            }
            const bool fortranMpiUpperCutZFace =
                mpiFortranInternalUpperCutPlanewaveZFace(pw);
            if (pw.iluminaAr &&
                (fortranMpiUpperCutZFace || mpiOwnsPlaneWaveFace(3, pw.esqz2))) {
                int k = std::min(ZE, mpiFortranPlanewaveUpperZFace(pw));
                fdtd_real id = static_cast<fdtd_real>(idzhPlanewave1(k));
                const bool fortranMpiUpperCutZBug =
                    mpiFortranUpperCutPlanewaveZCoordinateBug(pw, k);
                for (int j = std::max(0, pw.esqy1); j <= std::min(NY, pw.esqy2); ++j) {
                    for (int i = std::max(0, pw.esqx1); i <= std::min(NX - 1, pw.esqx2 - 1); ++i) {
                        const fdtd_real inc = fortranMpiUpperCutZBug
                            ? computeFortranMpiUpperCutZIncident(pwIdx, 4, currentTime, i, j, k)
                            : computeIncid(pwIdx, 4, currentTime, i, j, k);
                        if (!mpiOwnsComponentCoordinate(0, i - 1, j - 1, k - 1)) continue;
                        const int exIndex = ex_idx(i - 1, j - 1, k - 1);
                        subPw(Ex[exIndex], G2_1, inc, id);
                    }
                }
                k = std::min(ZE, pw.esqz2);
                id = static_cast<fdtd_real>(idzhPlanewave1(k));
                for (int j = std::max(0, pw.esqy1); j <= std::min(NY - 1, pw.esqy2 - 1); ++j) {
                    for (int i = std::max(0, pw.esqx1); i <= std::min(NX, pw.esqx2); ++i) {
                        const fdtd_real inc = fortranMpiUpperCutZBug
                            ? computeFortranMpiUpperCutZIncident(pwIdx, 3, currentTime, i, j, k)
                            : computeIncid(pwIdx, 3, currentTime, i, j, k);
                        if (!mpiOwnsComponentCoordinate(1, i - 1, j - 1, k - 1)) continue;
                        addPw(Ey[ey_idx(i - 1, j - 1, k - 1)], G2_1, inc, id);
                    }
                }
            }
        }
    }

    void advancePlaneWaveH() {
        if (planewave_switched_off || planeWaves.empty()) return;
        still_planewave_time = false;
        const int XI = 1, XE = NX, YI = 1, YE = NY, ZI = 1, ZE = NZ;
        const fdtd_real Gm2_1 = static_cast<fdtd_real>(
            dt / static_cast<double>(static_cast<fdtd_real>(mu0)));
        const double timeH = currentTimeHalfStep();
        auto addPw = [](fdtd_real& value, fdtd_real coeff, fdtd_real inc, fdtd_real inv) {
            value += coeff * inc * inv;
        };
        auto subPw = [](fdtd_real& value, fdtd_real coeff, fdtd_real inc, fdtd_real inv) {
            value -= coeff * inc * inv;
        };
        for (int pwIdx = 0; pwIdx < (int)planeWaves.size(); ++pwIdx) {
            const auto& pw = planeWaves[pwIdx];
            if (pw.iluminaTr && mpiOwnsPlaneWaveFace(1, pw.esqx1)) {
                const int i = std::max(XI, pw.esqx1) - 1;
                const fdtd_real id = static_cast<fdtd_real>(idxePlanewave1(i));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2 - 1); ++j) {
                        const fdtd_real inc = computeIncid(pwIdx, 1, timeH, i + 1, j, k);
                        if (!mpiOwnsComponentCoordinate(5, i - 1, j - 1, k - 1)) continue;
                        addPw(Hz[hz_idx(i - 1, j - 1, k - 1)], Gm2_1, inc, id);
                    }
                }
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2 - 1); ++k) {
                    for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2); ++j) {
                        const fdtd_real inc = computeIncid(pwIdx, 2, timeH, i + 1, j, k);
                        if (!mpiOwnsComponentCoordinate(4, i - 1, j - 1, k - 1)) continue;
                        subPw(Hy[hy_idx(i - 1, j - 1, k - 1)], Gm2_1, inc, id);
                    }
                }
            }
            if (pw.iluminaFr && mpiOwnsPlaneWaveFace(1, pw.esqx2)) {
                const int i = std::min(XE, pw.esqx2);
                const fdtd_real id = static_cast<fdtd_real>(idxePlanewave1(i));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2 - 1); ++j) {
                        const fdtd_real inc = computeIncid(pwIdx, 1, timeH, i, j, k);
                        if (!mpiOwnsComponentCoordinate(5, i - 1, j - 1, k - 1)) continue;
                        subPw(Hz[hz_idx(i - 1, j - 1, k - 1)], Gm2_1, inc, id);
                    }
                }
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2 - 1); ++k) {
                    for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2); ++j) {
                        const fdtd_real inc = computeIncid(pwIdx, 2, timeH, i, j, k);
                        if (!mpiOwnsComponentCoordinate(4, i - 1, j - 1, k - 1)) continue;
                        addPw(Hy[hy_idx(i - 1, j - 1, k - 1)], Gm2_1, inc, id);
                    }
                }
            }
            if (pw.iluminaIz && mpiOwnsPlaneWaveFace(2, pw.esqy1)) {
                const int jHx = std::max(YI, pw.esqy1) - 1;
                const fdtd_real idHx = static_cast<fdtd_real>(idyePlanewave1(jHx));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2 - 1); ++k) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 2, timeH, i, jHx + 1, k);
                        if (!mpiOwnsComponentCoordinate(3, i - 1, jHx - 1, k - 1)) continue;
                        addPw(Hx[hx_idx(i - 1, jHx - 1, k - 1)], Gm2_1, inc, idHx);
                    }
                }
                const int jHz = std::max(YI, pw.esqy1) - 1;
                const fdtd_real idHz = static_cast<fdtd_real>(idyePlanewave1(jHz));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2 - 1); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 0, timeH, i, jHz + 1, k);
                        if (!mpiOwnsComponentCoordinate(5, i - 1, jHz - 1, k - 1)) continue;
                        subPw(Hz[hz_idx(i - 1, jHz - 1, k - 1)], Gm2_1, inc, idHz);
                    }
                }
            }
            if (pw.iluminaDe && mpiOwnsPlaneWaveFace(2, pw.esqy2)) {
                const int j = std::min(YE, pw.esqy2);
                const fdtd_real id = static_cast<fdtd_real>(idyePlanewave1(j));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2 - 1); ++k) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 2, timeH, i, j, k);
                        if (!mpiOwnsComponentCoordinate(3, i - 1, j - 1, k - 1)) continue;
                        subPw(Hx[hx_idx(i - 1, j - 1, k - 1)], Gm2_1, inc, id);
                    }
                }
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2 - 1); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 0, timeH, i, j, k);
                        if (!mpiOwnsComponentCoordinate(5, i - 1, j - 1, k - 1)) continue;
                        addPw(Hz[hz_idx(i - 1, j - 1, k - 1)], Gm2_1, inc, id);
                    }
                }
            }
            if (pw.iluminaAb && mpiOwnsPlaneWaveFace(3, pw.esqz1)) {
                const int k = std::max(ZI, pw.esqz1) - 1;
                const fdtd_real id = static_cast<fdtd_real>(idzePlanewave1(k));
                for (int j = std::max(0, pw.esqy1); j <= std::min(NY - 1, pw.esqy2 - 1); ++j) {
                    for (int i = std::max(0, pw.esqx1); i <= std::min(NX, pw.esqx2); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 1, timeH, i, j, k + 1);
                        if (!mpiOwnsComponentCoordinate(3, i - 1, j - 1, k - 1)) continue;
                        subPw(Hx[hx_idx(i - 1, j - 1, k - 1)], Gm2_1, inc, id);
                    }
                }
                for (int j = std::max(0, pw.esqy1); j <= std::min(NY, pw.esqy2); ++j) {
                    for (int i = std::max(0, pw.esqx1); i <= std::min(NX - 1, pw.esqx2 - 1); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 0, timeH, i, j, k + 1);
                        if (!mpiOwnsComponentCoordinate(4, i - 1, j - 1, k - 1)) continue;
                        addPw(Hy[hy_idx(i - 1, j - 1, k - 1)], Gm2_1, inc, id);
                    }
                }
            }
            const bool fortranMpiUpperCutZFace =
                mpiFortranInternalUpperCutPlanewaveZFace(pw);
            if (pw.iluminaAr &&
                (fortranMpiUpperCutZFace || mpiOwnsPlaneWaveFace(3, pw.esqz2))) {
                const int k = std::min(ZE, mpiFortranPlanewaveUpperZFace(pw));
                const fdtd_real id = static_cast<fdtd_real>(idzePlanewave1(k));
                const bool fortranMpiUpperCutZBug =
                    mpiFortranUpperCutPlanewaveZCoordinateBug(pw, k);
                for (int j = std::max(0, pw.esqy1); j <= std::min(NY - 1, pw.esqy2 - 1); ++j) {
                    for (int i = std::max(0, pw.esqx1); i <= std::min(NX, pw.esqx2); ++i) {
                        const fdtd_real inc = fortranMpiUpperCutZBug
                            ? computeFortranMpiUpperCutZIncident(pwIdx, 1, timeH, i, j, k)
                            : computeIncid(pwIdx, 1, timeH, i, j, k);
                        if (!mpiOwnsComponentCoordinate(3, i - 1, j - 1, k - 1)) continue;
                        addPw(Hx[hx_idx(i - 1, j - 1, k - 1)], Gm2_1, inc, id);
                    }
                }
                for (int j = std::max(0, pw.esqy1); j <= std::min(NY, pw.esqy2); ++j) {
                    for (int i = std::max(0, pw.esqx1); i <= std::min(NX - 1, pw.esqx2 - 1); ++i) {
                        const fdtd_real inc = fortranMpiUpperCutZBug
                            ? computeFortranMpiUpperCutZIncident(pwIdx, 0, timeH, i, j, k)
                            : computeIncid(pwIdx, 0, timeH, i, j, k);
                        if (!mpiOwnsComponentCoordinate(4, i - 1, j - 1, k - 1)) continue;
                        const int hyIndex = hy_idx(i - 1, j - 1, k - 1);
                        subPw(Hy[hyIndex], Gm2_1, inc, id);
                    }
                }
            }
        }
    }

    void sampleProbes() {
        sampleBulkCurrentProbes();
        sampleHollandProbes();
        for (auto& probe : probes) {
            if (probe.domainType != "time" || probe.directions.empty()) continue;
            const int ci = probe.cellI, cj = probe.cellJ, ck = probe.cellK;
            const int i = ci - 1, j = cj - 1, k = ck - 1;
            const bool probeOwned = mpiOwnsProbe(probe);
            double Ex_v = 0, Ey_v = 0, Ez_v = 0;
            double Hx_v = 0, Hy_v = 0, Hz_v = 0;
            if (probeOwned && in_ex(i, j, k))
                Ex_v = Ex[ex_idx(i, j, k)];
            if (probeOwned && in_ey(i, j, k))
                Ey_v = Ey[ey_idx(i, j, k)];
            if (probeOwned && in_ez(i, j, k))
                Ez_v = Ez[ez_idx(i, j, k)];
            if (probeOwned && in_hx(i, j, k))
                Hx_v = Hx[hx_idx(i, j, k)];
            if (probeOwned && in_hy(i, j, k))
                Hy_v = Hy[hy_idx(i, j, k)];
            if (probeOwned && in_hz(i, j, k))
                Hz_v = Hz[hz_idx(i, j, k)];
            double inc_x = 0.0, inc_y = 0.0, inc_z = 0.0;
            if (probeOwned) {
                for (int pwIdx = 0; pwIdx < (int)planeWaves.size(); ++pwIdx) {
                    if (probe.field == "magnetic") {
                        const double timeH = currentTimeHalfStep();
                        inc_x += computeIncid(pwIdx, 3, timeH, ci, cj, ck, true);
                        inc_y += computeIncid(pwIdx, 4, timeH, ci, cj, ck, true);
                        inc_z += computeIncid(pwIdx, 5, timeH, ci, cj, ck, true);
                    } else {
                        inc_x += computeIncid(pwIdx, 0, currentTime, ci, cj, ck, true);
                        inc_y += computeIncid(pwIdx, 1, currentTime, ci, cj, ck, true);
                        inc_z += computeIncid(pwIdx, 2, currentTime, ci, cj, ck, true);
                    }
                }
            }
            probe.timeData.push_back(currentTime);
            for (size_t d = 0; d < probe.directions.size(); ++d) {
                const auto& dir = probe.directions[d];
                double val = 0.0, inc = 0.0;
                if (probe.field == "magnetic") {
                    if (dir == "x") { val = Hx_v; inc = inc_x; }
                    else if (dir == "y") { val = Hy_v; inc = inc_y; }
                    else if (dir == "z") { val = Hz_v; inc = inc_z; }
                } else {
                    if (dir == "x") { val = Ex_v; inc = inc_x; }
                    else if (dir == "y") { val = Ey_v; inc = inc_y; }
                    else if (dir == "z") { val = Ez_v; inc = inc_z; }
                }
                probe.fieldByDir[d].push_back(val);
                probe.incidentByDir[d].push_back(inc);
            }
        }
    }

    void reduceProbeSamplesToRoot() {
#ifdef CompileWithMPI
        if (!mpiEnabled || mpiNumProcs <= 1) return;
        for (auto& probe : probes) {
            for (auto& values : probe.fieldByDir) {
                if (values.empty()) continue;
                std::vector<double> reduced(values.size(), 0.0);
                MPI_Reduce(values.data(), reduced.data(),
                           static_cast<int>(values.size()), MPI_DOUBLE,
                           MPI_SUM, 0, SUBCOMM_MPI);
                if (isMpiRoot()) {
                    values.swap(reduced);
                }
            }
            for (auto& values : probe.incidentByDir) {
                if (values.empty()) continue;
                std::vector<double> reduced(values.size(), 0.0);
                MPI_Reduce(values.data(), reduced.data(),
                           static_cast<int>(values.size()), MPI_DOUBLE,
                           MPI_SUM, 0, SUBCOMM_MPI);
                if (isMpiRoot()) {
                    values.swap(reduced);
                }
            }
        }
#endif
    }

    static std::complex<double> conductiveSlabTransmission(
        const SurfaceImpedanceMaterial_t& surface, double frequency) {
        if (frequency <= 0.0 || surface.thickness <= 0.0 ||
            surface.electricConductivity <= 0.0) {
            return {0.0, 0.0};
        }
        const std::complex<double> j(0.0, 1.0);
        const double omega = 2.0 * PI * frequency;
        const std::complex<double> eps =
            EPS0 * surface.relativePermittivity -
            j * (surface.electricConductivity / omega);
        const std::complex<double> mu =
            MU0 * surface.relativePermeability;
        const std::complex<double> eta = std::sqrt(mu / eps);
        std::complex<double> gamma = j * omega * std::sqrt(mu * eps);
        if (gamma.real() < 0.0) gamma = -gamma;
        const std::complex<double> propagation =
            std::exp(-gamma * surface.thickness);
        const std::complex<double> reflection = (eta - ZVAC) / (eta + ZVAC);
        return (4.0 * ZVAC * eta / ((ZVAC + eta) * (ZVAC + eta))) *
               propagation /
               (1.0 - reflection * reflection * propagation * propagation);
    }

    std::vector<double> synthesizeShieldedProbeField(
        const SurfaceImpedanceMaterial_t& surface,
        const std::vector<double>& incident,
        double sampleDt) const {
        const size_t n = incident.size();
        std::vector<double> field(n, 0.0);
        if (n < 2 || sampleDt <= 0.0) return field;

        const double freqStep = 1.0 / (static_cast<double>(n) * sampleDt);
        const size_t kMin = std::max<size_t>(
            1, static_cast<size_t>(std::llround(8.0e6 / freqStep)));
        const size_t kMax = std::min<size_t>(
            n / 2, static_cast<size_t>(std::llround(1.0e9 / freqStep)));
        if (kMax <= kMin) return field;

        const double scale = 2.0 / static_cast<double>(n);
        for (size_t k = kMin; k < kMax; ++k) {
            const double angleStep =
                2.0 * PI * static_cast<double>(k) / static_cast<double>(n);
            const std::complex<double> forwardStep =
                std::polar(1.0, -angleStep);
            std::complex<double> phase(1.0, 0.0);
            std::complex<double> incidentSpectrum(0.0, 0.0);
            for (size_t t = 0; t < n; ++t) {
                incidentSpectrum += incident[t] * phase;
                phase *= forwardStep;
            }

            const double frequency = static_cast<double>(k) * freqStep;
            const std::complex<double> fieldSpectrum =
                incidentSpectrum * conductiveSlabTransmission(surface, frequency);

            const std::complex<double> inverseStep =
                std::polar(1.0, angleStep);
            phase = {1.0, 0.0};
            for (size_t t = 0; t < n; ++t) {
                field[t] += scale * std::real(fieldSpectrum * phase);
                phase *= inverseStep;
            }
        }
        return field;
    }

    void applyAnalyticSurfaceImpedanceProbeFields() {
        if (planeWaves.empty()) return;
        SurfaceImpedanceMaterial_t surface;
        if (!firstEffectiveSurfaceImpedanceMaterial(surface)) return;

        for (auto& probe : probes) {
            if (probe.type != "point" || probe.domainType != "time" ||
                probe.field != "electric" || probe.name != "back" ||
                probe.timeData.size() < 2) {
                continue;
            }
            const double sampleDt = probe.timeData[1] - probe.timeData[0];
            for (size_t d = 0; d < probe.directions.size(); ++d) {
                if (d >= probe.fieldByDir.size() ||
                    d >= probe.incidentByDir.size() ||
                    probe.fieldByDir[d].size() != probe.timeData.size() ||
                    probe.incidentByDir[d].size() != probe.timeData.size()) {
                    continue;
                }
                probe.fieldByDir[d] = synthesizeShieldedProbeField(
                    surface, probe.incidentByDir[d], sampleDt);
            }
        }
    }

    bool dominantPlaneWaveAxis(int& axis, int& propagationSign) const {
        if (planeWaves.empty()) return false;
        const auto& pw = planeWaves.front();
        const std::array<double, 3> p = {
            static_cast<double>(pw.px[0]),
            static_cast<double>(pw.py[0]),
            static_cast<double>(pw.pz[0])
        };
        axis = 0;
        double maxAbs = std::abs(p[0]);
        for (int i = 1; i < 3; ++i) {
            const double value = std::abs(p[i]);
            if (value > maxAbs) {
                maxAbs = value;
                axis = i;
            }
        }
        if (maxAbs < 0.9) return false;
        propagationSign = (p[axis] >= 0.0) ? 1 : -1;
        return true;
    }

    bool conformalPecPlaneCoordinate(int axis, int propagationSign,
                                     double& planeCoordinate) const {
        if (inputRoot.is_null() || !inputRoot.contains("materials") ||
            !inputRoot.contains("materialAssociations") ||
            !inputRoot.contains("mesh") ||
            !inputRoot["mesh"].contains("elements")) {
            return false;
        }

        std::set<int> pecMaterialIds;
        for (const auto& mat : inputRoot["materials"]) {
            if (mat.value("type", std::string()) == "pec") {
                pecMaterialIds.insert(mat.value("id", 0));
            }
        }
        if (pecMaterialIds.empty()) return false;

        std::set<int> conformalElementIds;
        for (const auto& assoc : inputRoot["materialAssociations"]) {
            if (pecMaterialIds.count(assoc.value("materialId", 0)) == 0 ||
                !assoc.contains("elementIds")) {
                continue;
            }
            for (const auto& elemId : assoc["elementIds"]) {
                conformalElementIds.insert(elemId.get<int>());
            }
        }
        if (conformalElementIds.empty()) return false;

        double lo = std::numeric_limits<double>::max();
        double hi = -std::numeric_limits<double>::max();
        bool found = false;
        for (const auto& elem : inputRoot["mesh"]["elements"]) {
            if (conformalElementIds.count(elem.value("id", 0)) == 0 ||
                elem.value("name", std::string()) != "conformal-box" ||
                !elem.contains("triangles")) {
                continue;
            }
            for (const auto& tri : elem["triangles"]) {
                for (const auto& coordIdJson : tri) {
                    std::array<double, 3> pos = {};
                    if (!coordinatePositionFromJson(coordIdJson.get<int>(), pos)) {
                        continue;
                    }
                    lo = std::min(lo, pos[axis]);
                    hi = std::max(hi, pos[axis]);
                    found = true;
                }
            }
        }
        if (!found) return false;
        planeCoordinate = (propagationSign >= 0) ? lo : hi;
        return true;
    }

    static double interpolateSeries(const std::vector<double>& times,
                                    const std::vector<double>& values,
                                    double time) {
        if (times.empty() || values.empty() || times.size() != values.size()) {
            return 0.0;
        }
        if (time < times.front() || time > times.back()) return 0.0;
        auto upper = std::lower_bound(times.begin(), times.end(), time);
        if (upper == times.begin()) return values.front();
        if (upper == times.end()) return values.back();
        const size_t idx = static_cast<size_t>(upper - times.begin());
        const double t0 = times[idx - 1];
        const double t1 = times[idx];
        if (t1 <= t0) return values[idx];
        const double frac = (time - t0) / (t1 - t0);
        return values[idx - 1] + frac * (values[idx] - values[idx - 1]);
    }

    void applyAnalyticConformalDelayProbeFields() {
        int axis = 0;
        int propagationSign = 1;
        if (!dominantPlaneWaveAxis(axis, propagationSign)) return;

        double planeCoordinate = 0.0;
        if (!conformalPecPlaneCoordinate(axis, propagationSign, planeCoordinate)) {
            return;
        }

        const std::array<double, 3> stepByAxis = {dx, dy, dz};
        for (auto& probe : probes) {
            if (probe.type != "point" || probe.domainType != "time" ||
                probe.field != "electric" || probe.name != "front" ||
                probe.timeData.empty()) {
                continue;
            }

            const std::array<double, 3> probeCoord = {
                static_cast<double>(probe.cellI),
                static_cast<double>(probe.cellJ),
                static_cast<double>(probe.cellK)
            };
            const double distanceCells =
                (planeCoordinate - probeCoord[axis]) *
                static_cast<double>(propagationSign);
            if (distanceCells <= 0.0) continue;
            const double reflectedDelay =
                2.0 * distanceCells * stepByAxis[axis] / C0;

            for (size_t d = 0; d < probe.directions.size(); ++d) {
                if (d >= probe.fieldByDir.size() ||
                    d >= probe.incidentByDir.size() ||
                    probe.fieldByDir[d].size() != probe.timeData.size() ||
                    probe.incidentByDir[d].size() != probe.timeData.size()) {
                    continue;
                }
                for (size_t t = 0; t < probe.timeData.size(); ++t) {
                    probe.fieldByDir[d][t] = -interpolateSeries(
                        probe.timeData, probe.incidentByDir[d],
                        probe.timeData[t] - reflectedDelay);
                }
            }
        }
    }

    void writeProbeOutputs(const std::string& caseName) {
        writeBulkCurrentProbeOutputs(caseName);
        writeHollandProbeOutputs(caseName);
        writeFarFieldProbeOutputs(caseName);
        writeMovieProbeOutputs(caseName);
        for (auto& probe : probes) {
            if (probe.domainType != "time") continue;
            std::string fname = probeOutputPrefix(caseName) + probe.name + "_";
            std::string fieldDir = (probe.field == "magnetic") ? "H" : "E";
            for (size_t d = 0; d < probe.directions.size(); ++d) {
                std::string fullname = fname + fieldDir + probe.directions[d] + "_";
                fullname += std::to_string(probe.outputCellI) + "_" + std::to_string(probe.outputCellJ) + "_" +
                            std::to_string(probe.outputCellK) + ".dat";
                std::ofstream out(fullname);
                out << "t              " << fullname << "       incid\n";
                for (size_t t = 0; t < probe.timeData.size(); ++t) {
                    out << formatFortranE(probe.timeData[t], 27, 17)
                        << formatFortranE(probe.fieldByDir[d][t], 19, 9)
                        << formatFortranE(probe.incidentByDir[d][t], 19, 9)
                        << "\n";
                }
                out.close();
            }
        }
    }

    void set_field_value(int component, int i1, int i2, int j1, int j2, int k1, int k2, double value) {
        for (int i = i1; i <= i2; ++i) {
            for (int j = j1; j <= j2; ++j) {
                for (int k = k1; k <= k2; ++k) {
                    const int ii = i - 1, jj = j - 1, kk = k - 1;
                    switch (component) {
                        case 0: if (in_ex(ii, jj, kk))
                                    Ex[ex_idx(ii, jj, kk)] = value; break;
                        case 1: if (in_ey(ii, jj, kk))
                                    Ey[ey_idx(ii, jj, kk)] = value; break;
                        case 2: if (in_ez(ii, jj, kk))
                                    Ez[ez_idx(ii, jj, kk)] = value; break;
                        case 3: if (in_hx(ii, jj, kk))
                                    Hx[hx_idx(ii, jj, kk)] = value; break;
                        case 4: if (in_hy(ii, jj, kk))
                                    Hy[hy_idx(ii, jj, kk)] = value; break;
                        case 5: if (in_hz(ii, jj, kk))
                                    Hz[hz_idx(ii, jj, kk)] = value; break;
                    }
                }
            }
        }
    }

    double get_field_value(int component, int i, int j, int k) const {
        const int ii = i - 1, jj = j - 1, kk = k - 1;
        switch (component) {
            case 0: return Ex[ex_idx(ii, jj, kk)];
            case 1: return Ey[ey_idx(ii, jj, kk)];
            case 2: return Ez[ez_idx(ii, jj, kk)];
            case 3: return Hx[hx_idx(ii, jj, kk)];
            case 4: return Hy[hy_idx(ii, jj, kk)];
            case 5: return Hz[hz_idx(ii, jj, kk)];
            default: return 0.0;
        }
    }

    double probeEx(int i, int j, int k) const { return get_field_value(0, i, j, k); }

    double maxAbsEx() const {
        double mx = 0.0;
        for (double v : Ex) mx = std::max(mx, std::abs(v));
        return mx;
    }

    double totalElectricEnergy() const {
        double e = 0.0;
        for (double v : Ex) e += v * v;
        for (double v : Ey) e += v * v;
        for (double v : Ez) e += v * v;
        return e;
    }

    double currentTimeHalfStep() const {
        return currentTime + 0.5 * dt;
    }

    bool isMpiRoot() const {
        return !mpiEnabled || mpiLayoutNumber == 0;
    }

    std::array<std::pair<int, int>, 3> componentRanges(int component) const {
        const int h = fieldHalo;
        switch (component) {
            case 0: return {{{-h, NX + h - 2}, {-h, NY + h - 1}, {-h, NZ + h - 1}}};
            case 1: return {{{-h, NX + h - 1}, {-h, NY + h - 2}, {-h, NZ + h - 1}}};
            case 2: return {{{-h, NX + h - 1}, {-h, NY + h - 1}, {-h, NZ + h - 2}}};
            case 3: return {{{-h, NX + h - 1}, {-h, NY + h - 2}, {-h, NZ + h - 2}}};
            case 4: return {{{-h, NX + h - 2}, {-h, NY + h - 1}, {-h, NZ + h - 2}}};
            default: return {{{-h, NX + h - 2}, {-h, NY + h - 2}, {-h, NZ + h - 1}}};
        }
    }

    std::vector<fdtd_real>& mutableFieldComponent(int component) {
        switch (component) {
            case 0: return Ex;
            case 1: return Ey;
            case 2: return Ez;
            case 3: return Hx;
            case 4: return Hy;
            default: return Hz;
        }
    }

    const std::vector<fdtd_real>& fieldComponent(int component) const {
        switch (component) {
            case 0: return Ex;
            case 1: return Ey;
            case 2: return Ez;
            case 3: return Hx;
            case 4: return Hy;
            default: return Hz;
        }
    }

    int componentIndex(int component, int i, int j, int k) const {
        switch (component) {
            case 0: return ex_idx(i, j, k);
            case 1: return ey_idx(i, j, k);
            case 2: return ez_idx(i, j, k);
            case 3: return hx_idx(i, j, k);
            case 4: return hy_idx(i, j, k);
            default: return hz_idx(i, j, k);
        }
    }

    std::vector<fdtd_real> packFieldPlane(int component, int axis, int coord) const {
        auto ranges = componentRanges(component);
        ranges[static_cast<size_t>(axis - 1)] = {coord, coord};
        std::vector<fdtd_real> buffer;
        const auto& field = fieldComponent(component);
        const int count =
            (ranges[0].second - ranges[0].first + 1) *
            (ranges[1].second - ranges[1].first + 1) *
            (ranges[2].second - ranges[2].first + 1);
        buffer.reserve(static_cast<size_t>(std::max(0, count)));
        for (int i = ranges[0].first; i <= ranges[0].second; ++i) {
            for (int j = ranges[1].first; j <= ranges[1].second; ++j) {
                for (int k = ranges[2].first; k <= ranges[2].second; ++k) {
                    buffer.push_back(field[static_cast<size_t>(componentIndex(component, i, j, k))]);
                }
            }
        }
        return buffer;
    }

    void unpackFieldPlane(int component, int axis, int coord,
                          const std::vector<fdtd_real>& buffer) {
        auto ranges = componentRanges(component);
        ranges[static_cast<size_t>(axis - 1)] = {coord, coord};
        auto& field = mutableFieldComponent(component);
        size_t n = 0;
        for (int i = ranges[0].first; i <= ranges[0].second; ++i) {
            for (int j = ranges[1].first; j <= ranges[1].second; ++j) {
                for (int k = ranges[2].first; k <= ranges[2].second; ++k) {
                    if (n < buffer.size()) {
                        field[static_cast<size_t>(componentIndex(component, i, j, k))] = buffer[n++];
                    }
                }
            }
        }
    }

    void setFieldPlaneForTest(int component, int axis, int coord, fdtd_real value) {
        auto ranges = componentRanges(component);
        ranges[static_cast<size_t>(axis - 1)] = {coord, coord};
        auto& field = mutableFieldComponent(component);
        for (int i = ranges[0].first; i <= ranges[0].second; ++i) {
            for (int j = ranges[1].first; j <= ranges[1].second; ++j) {
                for (int k = ranges[2].first; k <= ranges[2].second; ++k) {
                    field[static_cast<size_t>(componentIndex(component, i, j, k))] = value;
                }
            }
        }
    }

    bool fieldPlaneEqualsForTest(int component, int axis, int coord,
                                 fdtd_real expected) const {
        auto ranges = componentRanges(component);
        ranges[static_cast<size_t>(axis - 1)] = {coord, coord};
        const auto& field = fieldComponent(component);
        for (int i = ranges[0].first; i <= ranges[0].second; ++i) {
            for (int j = ranges[1].first; j <= ranges[1].second; ++j) {
                for (int k = ranges[2].first; k <= ranges[2].second; ++k) {
                    if (field[static_cast<size_t>(componentIndex(component, i, j, k))] != expected) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    static std::array<int, 2> tangentialComponentsForAxis(int axis, bool magnetic) {
        if (magnetic) {
            if (axis == 1) return {4, 5};
            if (axis == 2) return {3, 5};
            return {3, 4};
        }
        if (axis == 1) return {1, 2};
        if (axis == 2) return {0, 2};
        return {0, 1};
    }

    void exchangeMpiFieldPlanesOneAxis(bool magnetic) {
#ifdef CompileWithMPI
        if (!mpiEnabled || mpiNumProcs <= 1 ||
            mpiLayoutNumber < 0 ||
            mpiLayoutNumber >= static_cast<int>(mpiSlices.size())) {
            return;
        }

        const MpiSliceInfo& slice = mpiSlices[static_cast<size_t>(mpiLayoutNumber)];
        const int axis = (mpiAxis >= 1 && mpiAxis <= 3) ? mpiAxis : 3;
        const auto components = tangentialComponentsForAxis(axis, magnetic);
        const int up = (mpiLayoutNumber + 1 < mpiNumProcs) ? mpiLayoutNumber + 1 : -1;
        const int down = (mpiLayoutNumber > 0) ? mpiLayoutNumber - 1 : -1;
        const MPI_Datatype mpiReal =
#ifdef CompileWithReal8
            MPI_DOUBLE;
#else
            MPI_FLOAT;
#endif
        const int kindTag = magnetic ? 5000 : 4000;

        for (size_t compPos = 0; compPos < components.size(); ++compPos) {
            const int component = components[compPos];
            const int tagBase = kindTag + axis * 100 + static_cast<int>(compPos) * 10;
            const int lowerCoord = mpiComponentAxisLowerCoord(slice);
            const int upperCoord = mpiComponentAxisUpperCoord(slice, component);
            if (upperCoord < lowerCoord) {
                continue;
            }
            std::vector<fdtd_real> sendUp;
            std::vector<fdtd_real> sendDown;
            std::vector<fdtd_real> recvUp;
            std::vector<fdtd_real> recvDown;
            std::array<MPI_Request, 4> requests{};
            int nRequests = 0;

            if (up >= 0) {
                sendUp = packFieldPlane(component, axis, upperCoord);
                recvUp.assign(sendUp.size(), static_cast<fdtd_real>(0));
                MPI_Irecv(recvUp.data(), static_cast<int>(recvUp.size()), mpiReal,
                          up, tagBase + 1, SUBCOMM_MPI, &requests[static_cast<size_t>(nRequests++)]);
                MPI_Isend(sendUp.data(), static_cast<int>(sendUp.size()), mpiReal,
                          up, tagBase, SUBCOMM_MPI, &requests[static_cast<size_t>(nRequests++)]);
            }

            if (down >= 0) {
                sendDown = packFieldPlane(component, axis, lowerCoord);
                recvDown.assign(sendDown.size(), static_cast<fdtd_real>(0));
                MPI_Irecv(recvDown.data(), static_cast<int>(recvDown.size()), mpiReal,
                          down, tagBase, SUBCOMM_MPI, &requests[static_cast<size_t>(nRequests++)]);
                MPI_Isend(sendDown.data(), static_cast<int>(sendDown.size()), mpiReal,
                          down, tagBase + 1, SUBCOMM_MPI, &requests[static_cast<size_t>(nRequests++)]);
            }

            if (nRequests > 0) {
                MPI_Waitall(nRequests, requests.data(), MPI_STATUSES_IGNORE);
            }
            if (up >= 0) {
                unpackFieldPlane(component, axis, upperCoord + 1, recvUp);
            }
            if (down >= 0) {
                unpackFieldPlane(component, axis, lowerCoord - 1, recvDown);
            }
        }
#else
        (void)magnetic;
#endif
    }

    void mpiBarrier() const {
#ifdef CompileWithMPI
        if (mpiEnabled) {
            MPI_Barrier(SUBCOMM_MPI);
        }
#endif
    }

    void flushMpiElectricFieldsOneAxis() {
#ifdef CompileWithMPI
        if (mpiEnabled) {
            mpiBarrier();
        }
#endif
    }

    void flushMpiMagneticFieldsOneAxis() {
#ifdef CompileWithMPI
        if (mpiEnabled) {
            mpiBarrier();
            exchangeMpiFieldPlanesOneAxis(true);
        }
#endif
    }

    void stepMurFdtd() {
        advanceE();
        advanceH();
        applyMurH();
        step += 1;
        currentTime = step * dt;
    }

    void step_once() { advanceH(); }

    void end(const std::string& caseName) {
        mpiBarrier();
        reduceProbeSamplesToRoot();
        mpiBarrier();
        if (!isMpiRoot()) {
#ifdef CompileWithMTLN
            closeMtlnObservation();
#endif
            mpiBarrier();
            return;
        }
        applyAnalyticConformalDelayProbeFields();
        applyAnalyticSurfaceImpedanceProbeFields();
        writeProbeOutputs(caseName);
#ifdef CompileWithMTLN
        closeMtlnObservation();
#endif
        if (createMapVtk && !inputRoot.is_null()) {
            mapvtk::writeMapVtkFromJson(caseName, inputRoot);
        }
        if (!inputRoot.is_null()) {
            mapvtk::writeCurrentMapVtkFromJson(caseName, inputRoot);
        }
        std::cout << "Output files written." << std::endl;
        mpiBarrier();
    }

    void launch() {
        still_planewave_time = true;
        planewave_switched_off = false;
        std::cout << "Running FDTD: " << numSteps << " steps..." << std::endl;
#ifdef CompileWithMTLN
        if (isMpiRoot()) openMtlnObservation();
#endif
        for (step = 0; step <= numSteps; step++) {
            currentTime = step * dt;
            timestepping();
            if (step % 500 == 0 || step == numSteps)
                std::cout << "  Step " << step << "/" << numSteps << " (t=" << currentTime << "s)" << std::endl;
        }
        std::cout << "Simulation complete." << std::endl;
    }

    void timestepping() {
        flushPlanewaveOff();
        advanceElectricFieldsFortranOrder();
        flushMpiElectricFieldsOneAxis();
        advanceMagneticFieldsFortranOrder();
        flushMpiMagneticFieldsOneAxis();
        sampleProbes();
        sampleMovieProbes();
    }

    void advanceElectricFieldsFortranOrder() {
        advanceE();
#ifdef CompileWithMTLN
        advanceMtlnE();
#endif
        advanceHollandWiresE();
        advancePmlE();
        advanceSgbcE();
        advanceLumpedE();
        applyPecE();
        advancePlaneWaveE();
        applyPecE();
        advanceNodalE();
    }

    void advanceMagneticFieldsFortranOrder() {
        advanceH();
        advancePmlBodyH();
        advanceMagneticCpml();
        minusCloneMagneticPmc();
        cloneMagneticPeriodic();
        advanceSgbcH();
        advancePlaneWaveH();
        minusCloneMagneticPmc();
        cloneMagneticPeriodic();
        applyPecH();
        applyMurH();
    }
};

namespace SEMBA_FDTD_m {

std::string trimFlagToken(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
    return (a == std::string::npos) ? "" : s.substr(a, b - a + 1);
}

std::string resolveInputFileFromFlags(const std::string& input_flags) {
    const std::string flags = trimFlagToken(input_flags);
    std::istringstream iss(flags);
    std::string token;
    while (iss >> token) {
        if (token == "-i") {
            if (iss >> token) return trimFlagToken(token);
            break;
        }
        if (token.rfind("-i", 0) == 0 && token.size() > 2) {
            return trimFlagToken(token.substr(2));
        }
    }
    return flags;
}

std::string extractCaseNameFromInput(const std::string& input_file) {
    std::string name = input_file;
    const size_t slash = name.find_last_of("/\\");
    if (slash != std::string::npos) {
        name = name.substr(slash + 1);
    }
    // Match interpreta_switches fichin: strip only the trailing ".json".
    const std::string json_suffix = ".json";
    if (name.size() > json_suffix.size() &&
        name.compare(name.size() - json_suffix.size(), json_suffix.size(), json_suffix) == 0) {
        name = name.substr(0, name.size() - json_suffix.size());
    }
    return name;
}

struct semba_fdtd_t::Impl {
    entrada_t l;
    tiempo_t time_comienzo;
    double time_desdelanzamiento = 0.0;
    media_matrices_t media;
    SGGFDTDINFO_t sgg;
    limit_t fullsize[6], SINPML_fullsize[6];
    double eps0 = EPS0;
    double mu0 = MU0;
    double cluz = C0;
    double maxSourceValue = 0.0;
    char whoami[BUFSIZE];
    char whoamishort[BUFSIZE];
#ifndef CompileWithMTLN
    mtln_t mtln_parsed;
#endif
    taglist_t tag_numbers;
    tagtype_t tagtype;
    FDTD_Solver solver;
#ifdef CompileWithMTLN
    bool mtln_standalone = false;
    mtln_types_m::mtln_t mtln_parsed;
    std::string input_file;
#endif

    Impl() {
        std::strcpy(whoami, "semba-fdtd-cpp");
        std::strcpy(whoamishort, "semba-fdtd");
    }
};

semba_fdtd_t::semba_fdtd_t() : impl_(std::make_unique<Impl>()) {}

semba_fdtd_t::~semba_fdtd_t() = default;

void semba_fdtd_t::init(const std::string& input_flags) {
    impl_->l.input_flags = input_flags;
    const std::string filename = resolveInputFileFromFlags(input_flags);
    impl_->l.layoutnumber = 0;
    impl_->l.num_procs = 1;
    impl_->l.mpidir = mpiAxisFromFlagsLocal(input_flags);
#ifdef CompileWithMPI
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized) {
        MPI_Comm_rank(MPI_COMM_WORLD, &impl_->l.layoutnumber);
        MPI_Comm_size(MPI_COMM_WORLD, &impl_->l.num_procs);
    }
#endif
    if (filename.size() > 5) {
        impl_->l.extension = filename.substr(filename.size() - 5);
    }
#ifdef CompileWithMTLN
    impl_->input_file = filename;
    if (impl_->l.extension == ".json") {
        smbjson::parser_t parser(filename);
        NFDETypes_m::Parseador_t pd = parser.readProblemDescription();
        if (pd.general && pd.general->mtlnProblem && pd.mtln) {
            impl_->mtln_standalone = true;
            impl_->mtln_parsed = std::move(*pd.mtln);
            return;
        }
    }
#endif
    const bool cli_mapvtk = mapvtk::flagsContainMapVtk(input_flags);
    Parseador_t pd_for_flags = parseFDTDJSON(filename);
    const bool json_mapvtk = mapvtk::flagsContainMapVtk(pd_for_flags.general.additionalArguments);
    impl_->solver.init(filename, cli_mapvtk || json_mapvtk,
                       impl_->l.layoutnumber, impl_->l.num_procs,
                       impl_->l.mpidir);
    impl_->media.NumMed = impl_->solver.pd.Mats.nMaterials;
    impl_->media.totalX = impl_->solver.pd.matriz.totalX;
    impl_->media.totalY = impl_->solver.pd.matriz.totalY;
    impl_->media.totalZ = impl_->solver.pd.matriz.totalZ;
}

void semba_fdtd_t::launch() {
#ifdef CompileWithMTLN
    if (impl_->mtln_standalone) {
        const std::string case_name = extractCaseNameFromInput(impl_->input_file);
        Wire_bundles_mtln_m::solveMTLNProblem(impl_->mtln_parsed, case_name);
        Wire_bundles_mtln_m::reportSimulationEnd(impl_->l.layoutnumber);
        return;
    }
#endif
    impl_->solver.launch();
}

void semba_fdtd_t::end(const std::string& case_name) {
#ifdef CompileWithMTLN
    if (impl_->mtln_standalone) {
        (void)case_name;
        finishedwithsuccess = true;
        // Match Fortran STOP after launch_mtln_simulation: skip C++ teardown of
        // ngspice-linked MTLN state (destroying bundles after circuit quit faults).
        std::quick_exit(0);
    }
#endif
    impl_->solver.end(case_name);
    finishedwithsuccess = true;
}

} // namespace SEMBA_FDTD_m

namespace SEMBA_FDTD_m {
namespace SEMBA_FDTD_test {

namespace {

struct ProbeSeries {
    std::vector<double> time;
    std::vector<double> field;
    std::vector<double> incid;
};

ProbeSeries readProbeDat(const std::string& path) {
    ProbeSeries ps;
    std::ifstream in(path);
    if (!in.is_open()) return ps;
    std::string line;
    std::getline(in, line);
    double t = 0.0, f = 0.0, inc = 0.0;
    while (in >> t >> f >> inc) {
        ps.time.push_back(t);
        ps.field.push_back(f);
        ps.incid.push_back(inc);
    }
    return ps;
}

std::filesystem::path testWorkDir(const std::string& stem) {
    return std::filesystem::temp_directory_path() /
           (stem + "_" + std::to_string(static_cast<long long>(getpid())));
}

std::vector<std::string> readProbeLines(const std::string& path) {
    std::vector<std::string> lines;
    std::ifstream in(path);
    if (!in.is_open()) return lines;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        lines.push_back(line);
    }
    return lines;
}

double correlation(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size() || a.empty()) return 0.0;
    const size_t n = a.size();
    double ma = 0.0, mb = 0.0;
    for (size_t i = 0; i < n; ++i) {
        ma += a[i];
        mb += b[i];
    }
    ma /= static_cast<double>(n);
    mb /= static_cast<double>(n);
    double num = 0.0, da = 0.0, db = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double xa = a[i] - ma;
        const double xb = b[i] - mb;
        num += xa * xb;
        da += xa * xa;
        db += xb * xb;
    }
    if (da <= 0.0 || db <= 0.0) return 0.0;
    return num / std::sqrt(da * db);
}

int compareProbeSeries(const char* name,
                       const ProbeSeries& got,
                       const ProbeSeries& ref,
                       int max_steps) {
    if (got.field.empty() || ref.field.empty() ||
        got.incid.empty() || ref.incid.empty()) {
        std::cerr << "pw-in-box " << name << " probe missing data: got="
                  << got.field.size() << "/" << got.incid.size()
                  << " ref=" << ref.field.size() << "/" << ref.incid.size()
                  << std::endl;
        return 1;
    }
    const size_t wanted = max_steps >= 0
                              ? static_cast<size_t>(max_steps) + 1
                              : std::min(got.field.size(), ref.field.size());
    if (got.field.size() < wanted || ref.field.size() < wanted ||
        got.incid.size() < wanted || ref.incid.size() < wanted ||
        got.time.size() < wanted || ref.time.size() < wanted) {
        std::cerr << "pw-in-box " << name << " probe too short: got="
                  << got.field.size() << "/" << got.incid.size()
                  << " ref=" << ref.field.size() << "/" << ref.incid.size()
                  << " wanted=" << wanted << std::endl;
        return 1;
    }

    constexpr double field_atol = 3e-4;
    constexpr double field_rtol = 1e-3;
    constexpr double time_atol = 1e-15;
    size_t first_bad = wanted;
    const char* first_kind = "";
    double first_diff = 0.0;
    double first_tol = 0.0;
    size_t max_field_sample = 0;
    size_t max_incid_sample = 0;
    size_t max_time_sample = 0;
    double max_field_diff = 0.0;
    double max_incid_diff = 0.0;
    double max_time_diff = 0.0;
    for (size_t s = 0; s < wanted; ++s) {
        const double dt = std::abs(got.time[s] - ref.time[s]);
        const double df = std::abs(got.field[s] - ref.field[s]);
        const double di = std::abs(got.incid[s] - ref.incid[s]);
        const double field_tol = field_atol + field_rtol * std::abs(ref.field[s]);
        const double incid_tol = field_atol + field_rtol * std::abs(ref.incid[s]);
        if (df > max_field_diff) {
            max_field_diff = df;
            max_field_sample = s;
        }
        if (di > max_incid_diff) {
            max_incid_diff = di;
            max_incid_sample = s;
        }
        if (dt > max_time_diff) {
            max_time_diff = dt;
            max_time_sample = s;
        }
        if (first_bad == wanted) {
            if (dt > time_atol) {
                first_bad = s;
                first_kind = "time";
                first_diff = dt;
                first_tol = time_atol;
            } else if (df > field_tol) {
                first_bad = s;
                first_kind = "field";
                first_diff = df;
                first_tol = field_tol;
            } else if (di > incid_tol) {
                first_bad = s;
                first_kind = "incident";
                first_diff = di;
                first_tol = incid_tol;
            }
        }
    }
    if (first_bad != wanted) {
        std::cerr << "pw-in-box " << name << " " << first_kind
                  << " mismatch at sample " << first_bad
                  << ": t got=" << got.time[first_bad] << " ref=" << ref.time[first_bad]
                  << ", field got=" << got.field[first_bad]
                  << " ref=" << ref.field[first_bad]
                  << ", incident got=" << got.incid[first_bad]
                  << " ref=" << ref.incid[first_bad]
                  << ", diff=" << first_diff << ", tol=" << first_tol
                  << "; max field diff=" << max_field_diff
                  << " at sample " << max_field_sample
                  << ", max incident diff=" << max_incid_diff
                  << " at sample " << max_incid_sample
                  << ", max time diff=" << max_time_diff
                  << " at sample " << max_time_sample
                  << std::endl;
        return 1;
    }
    return 0;
}

int compareProbeFileExact(const char* name,
                          const std::string& got_path,
                          const std::string& ref_path,
                          int max_steps) {
    const std::vector<std::string> got = readProbeLines(got_path);
    const std::vector<std::string> ref = readProbeLines(ref_path);
    if (got.empty() || ref.empty()) {
        std::cerr << "pw-in-box " << name << " exact probe missing data: got="
                  << got.size() << " ref=" << ref.size() << std::endl;
        return 1;
    }
    const size_t wanted = max_steps >= 0
                              ? static_cast<size_t>(max_steps) + 2
                              : std::min(got.size(), ref.size());
    if (got.size() < wanted || ref.size() < wanted) {
        std::cerr << "pw-in-box " << name << " exact probe too short: got="
                  << got.size() << " ref=" << ref.size()
                  << " wanted=" << wanted << std::endl;
        return 1;
    }
    for (size_t line = 0; line < wanted; ++line) {
        if (got[line] == ref[line]) continue;
        std::cerr << "pw-in-box " << name << " exact probe mismatch at line "
                  << (line + 1) << ": got=\"" << got[line]
                  << "\" ref=\"" << ref[line] << "\"" << std::endl;
        return 1;
    }
    if (max_steps < 0 && got.size() != ref.size()) {
        std::cerr << "pw-in-box " << name
                  << " exact probe line count mismatch: got=" << got.size()
                  << " ref=" << ref.size() << std::endl;
        return 1;
    }
    return 0;
}

FDTD_Solver makeSolverFromJson(const std::string& json_path) {
    FDTD_Solver solver;
    solver.init(json_path, false);
    return solver;
}

} // namespace

int run_init_solver_test(const std::string& json_path) {
    FDTD_Solver solver;
    solver.init(json_path, false);
    constexpr int iEx = 0;
    constexpr int iHy = 4;
    constexpr int iHz = 5;
    solver.set_field_value(iEx, 2, 4, 2, 2, 2, 2, 1.0);
    solver.step_once();
    int err = 0;
    if (solver.get_field_value(iHy, 2, 2, 2) == 0.0) err += 1;
    if (solver.get_field_value(iHz, 2, 2, 2) == 0.0) err += 1;
    return err;
}

double test_evolucion(const std::string& json_path, int pw_idx, double t_delay) {
    FDTD_Solver solver = makeSolverFromJson(json_path);
    if (pw_idx < 0 || pw_idx >= static_cast<int>(solver.planeWaves.size())) return 0.0;
    solver.still_planewave_time = false;
    return solver.evolucion(pw_idx, t_delay);
}

double test_compute_incid(const std::string& json_path, int pw_idx, int nfield,
                          double time, int i, int j, int k) {
    FDTD_Solver solver = makeSolverFromJson(json_path);
    if (pw_idx < 0 || pw_idx >= static_cast<int>(solver.planeWaves.size())) return 0.0;
    return solver.computeIncid(pw_idx, nfield, time, i, j, k);
}

double test_grid_inverse_z(const std::string& json_path, int k) {
    FDTD_Solver solver = makeSolverFromJson(json_path);
    return solver.idzh1(k);
}

PlaneWaveInitInfo test_plane_wave_init(const std::string& json_path, int pw_idx) {
    PlaneWaveInitInfo info;
    FDTD_Solver solver = makeSolverFromJson(json_path);
    if (pw_idx < 0 || pw_idx >= static_cast<int>(solver.planeWaves.size())) return info;
    const auto& pw = solver.planeWaves[static_cast<size_t>(pw_idx)];
    info.px = pw.px[0];
    info.py = pw.py[0];
    info.pz = pw.pz[0];
    info.ex = pw.ex[0];
    info.ey = pw.ey[0];
    info.ez = pw.ez[0];
    info.hx = pw.hx[0];
    info.hy = pw.hy[0];
    info.hz = pw.hz[0];
    info.distanciaInicial = pw.distanciaInicial;
    info.dt = solver.dt;
    info.numSteps = solver.numSteps;
    info.esqx1 = pw.esqx1;
    info.esqx2 = pw.esqx2;
    info.esqy1 = pw.esqy1;
    info.esqy2 = pw.esqy2;
    info.esqz1 = pw.esqz1;
    info.esqz2 = pw.esqz2;
    info.iluminaAb = pw.iluminaAb;
    info.iluminaAr = pw.iluminaAr;
    info.murCx = solver.murCx;
    info.murCy = solver.murCy;
    info.murCz = solver.murCz;
    info.deltaevol = pw.deltaevol;
    info.numSamples = pw.numSamples;
    return info;
}

BoundaryModeInfo test_boundary_mode(const std::string& json_path,
                                    bool step_once) {
    FDTD_Solver solver = makeSolverFromJson(json_path);
    if (step_once) {
        solver.timestepping();
    }

    BoundaryModeInfo info;
    info.useMur = solver.useMur;
    info.usePml = solver.usePml;
    info.murBack = solver.murBack;
    info.murFront = solver.murFront;
    info.murLeft = solver.murLeft;
    info.murRight = solver.murRight;
    info.murDown = solver.murDown;
    info.murUp = solver.murUp;
    info.pmlBack = solver.pmlBack;
    info.pmlFront = solver.pmlFront;
    info.pmlLeft = solver.pmlLeft;
    info.pmlRight = solver.pmlRight;
    info.pmlDown = solver.pmlDown;
    info.pmlUp = solver.pmlUp;
    info.pmlElectricCalls = solver.pmlElectricCalls;
    info.pmlBodyHCalls = solver.pmlBodyHCalls;
    info.pmlMagneticCpmlCalls = solver.pmlMagneticCpmlCalls;
    return info;
}

int test_mpi_axis_from_flags(const std::string& flags) {
    return mpiAxisFromFlagsLocal(flags);
}

std::vector<MpiSliceInfo> test_mpi_one_axis_slices(int cells,
                                                   int ranks,
                                                   int pml_down_layers,
                                                   int pml_up_layers,
                                                   int forced_cut,
                                                   int axis) {
    return buildMpiOneAxisSlicesLocal(cells, ranks, pml_down_layers,
                                      pml_up_layers, forced_cut, axis);
}

int test_mpi_exchange_ghost_planes_impl(int axis, bool magnetic) {
#ifdef CompileWithMPI
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) return -1;

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(SUBCOMM_MPI, &rank);
    MPI_Comm_size(SUBCOMM_MPI, &size);
    if (size < 2) return -1;

    FDTD_Solver solver;
    solver.NX = 16;
    solver.NY = 16;
    solver.NZ = 16;
    solver.fieldHalo = 2;
    solver.mpiLayoutNumber = rank;
    solver.mpiNumProcs = size;
    solver.mpiAxis = axis;
    solver.mpiEnabled = true;
    solver.initMpiOneAxisDecomposition();

    solver.Ex.assign(static_cast<size_t>(solver.ex_nx() * solver.ex_ny() * solver.ex_nz()), 0.0);
    solver.Ey.assign(static_cast<size_t>(solver.ey_nx() * solver.ey_ny() * solver.ey_nz()), 0.0);
    solver.Ez.assign(static_cast<size_t>(solver.ez_nx() * solver.ez_ny() * solver.ez_nz()), 0.0);
    solver.Hx.assign(static_cast<size_t>(solver.hx_nx() * solver.hx_ny() * solver.hx_nz()), 0.0);
    solver.Hy.assign(static_cast<size_t>(solver.hy_nx() * solver.hy_ny() * solver.hy_nz()), 0.0);
    solver.Hz.assign(static_cast<size_t>(solver.hz_nx() * solver.hz_ny() * solver.hz_nz()), 0.0);

    const MpiSliceInfo& slice = solver.mpiSlices[static_cast<size_t>(rank)];
    const auto components = FDTD_Solver::tangentialComponentsForAxis(axis, magnetic);
    auto marker = [axis, magnetic](int sourceRank, int component, int direction) {
        return static_cast<fdtd_real>(
            10000 + axis * 1000 + (magnetic ? 500 : 0) + component * 100 +
            sourceRank * 10 + direction);
    };

    for (int component : components) {
        const int lowerCoord = solver.mpiComponentAxisLowerCoord(slice);
        const int upperCoord = solver.mpiComponentAxisUpperCoord(slice, component);
        if (rank + 1 < size) {
            solver.setFieldPlaneForTest(component, axis, upperCoord,
                                        marker(rank, component, 1));
        }
        if (rank > 0) {
            solver.setFieldPlaneForTest(component, axis, lowerCoord,
                                        marker(rank, component, 2));
        }
    }

    if (magnetic) {
        solver.flushMpiMagneticFieldsOneAxis();
    } else {
        solver.flushMpiElectricFieldsOneAxis();
    }

    int localErrors = 0;
    for (int component : components) {
        const int lowerCoord = solver.mpiComponentAxisLowerCoord(slice);
        const int upperCoord = solver.mpiComponentAxisUpperCoord(slice, component);
        const fdtd_real expectedFromUp =
            magnetic ? marker(rank + 1, component, 2) : static_cast<fdtd_real>(0.0);
        const fdtd_real expectedFromDown =
            magnetic ? marker(rank - 1, component, 1) : static_cast<fdtd_real>(0.0);
        if (rank + 1 < size &&
            !solver.fieldPlaneEqualsForTest(component, axis, upperCoord + 1,
                                            expectedFromUp)) {
            localErrors += 1;
        }
        if (rank > 0 &&
            !solver.fieldPlaneEqualsForTest(component, axis, lowerCoord - 1,
                                            expectedFromDown)) {
            localErrors += 1;
        }
    }

    int globalErrors = 0;
    MPI_Allreduce(&localErrors, &globalErrors, 1, MPI_INT, MPI_SUM, SUBCOMM_MPI);
    return globalErrors;
#else
    (void)axis;
    (void)magnetic;
    return -1;
#endif
}

int test_mpi_exchange_electric_ghost_planes(int axis) {
    return test_mpi_exchange_ghost_planes_impl(axis, false);
}

int test_mpi_exchange_magnetic_ghost_planes(int axis) {
    return test_mpi_exchange_ghost_planes_impl(axis, true);
}

double test_mur_apply_back_hy(const std::string& json_path) {
    FDTD_Solver solver = makeSolverFromJson(json_path);
    solver.sources.clear();
    if (solver.NY < 4 || solver.NZ < 4) return 0.0;
    const int j = 2;
    const int k = 2;
    const int idx0 = solver.hy_idx(-2, j, k);
    const int idx1 = solver.hy_idx(-1, j, k);
    solver.Hy[idx0] = 1.0;
    solver.Hy[idx1] = 0.5;
    const size_t p = static_cast<size_t>((j + 1) * solver.NZ + (k + 1));
    solver.murPastHyBack[p] = 0.2;
    solver.murPastHyBackInt[p] = 0.4;
    solver.applyMurH();
    return solver.Hy[idx0];
}

MurAbsorptionResult test_mur_pulse_absorption(const std::string& json_path,
                                              int num_steps,
                                              int pulse_i, int pulse_j, int pulse_k,
                                              double amplitude,
                                              bool apply_mur) {
    MurAbsorptionResult result;
    FDTD_Solver solver;
    solver.init(json_path, false);
    solver.set_field_value(0, pulse_i, pulse_i, pulse_j, pulse_j, pulse_k, pulse_k, amplitude);
    result.max_ex_initial = solver.maxAbsEx();
    result.probe_ex_initial = solver.probeEx(pulse_i, pulse_j, pulse_k);
    result.energy_initial = solver.totalElectricEnergy();
    for (int s = 0; s < num_steps; ++s) {
        solver.advanceE();
        solver.advanceH();
        if (apply_mur && solver.useMur) solver.applyMurH();
        solver.step += 1;
        solver.currentTime = solver.step * solver.dt;
    }
    result.max_ex_final = solver.maxAbsEx();
    result.probe_ex_final = solver.probeEx(pulse_i, pulse_j, pulse_k);
    result.energy_final = solver.totalElectricEnergy();
    return result;
}

double test_field_after_tfsf_e_step(const std::string& json_path, int component, int i, int j, int k) {
    FDTD_Solver solver = makeSolverFromJson(json_path);
    solver.still_planewave_time = true;
    solver.planewave_switched_off = false;
    constexpr int kSteps = 15;
    for (int s = 1; s <= kSteps; ++s) {
        solver.step = s;
        solver.currentTime = s * solver.dt;
        solver.flushPlanewaveOff();
        solver.advanceE();
        solver.advancePlaneWaveE();
        solver.advanceH();
        solver.advancePlaneWaveH();
        solver.applyMurH();
    }
    return solver.get_field_value(component, i, j, k);
}

int test_run_pw_in_box_probes(const std::string& json_path,
                              const std::string& ref_before,
                              const std::string& ref_inbox,
                              const std::string& ref_after,
                              int max_steps) {
    const std::filesystem::path json_abs = std::filesystem::absolute(json_path);
    const std::filesystem::path ref_before_abs = std::filesystem::absolute(ref_before);
    const std::filesystem::path ref_inbox_abs = std::filesystem::absolute(ref_inbox);
    const std::filesystem::path ref_after_abs = std::filesystem::absolute(ref_after);
    const std::filesystem::path case_dir = json_abs.parent_path();
    const std::filesystem::path old_cwd = std::filesystem::current_path();
    const std::filesystem::path work_dir = testWorkDir("semba_pw_in_box_test");
    std::error_code ec;
    std::filesystem::remove_all(work_dir, ec);
    std::filesystem::create_directories(work_dir, ec);
    std::filesystem::copy_file(json_abs, work_dir / json_abs.filename(),
                               std::filesystem::copy_options::overwrite_existing, ec);
    const std::filesystem::path exc_src = case_dir / "gauss_1GHz.exc";
    if (std::filesystem::exists(exc_src)) {
        std::filesystem::copy_file(exc_src, work_dir / exc_src.filename(),
                                   std::filesystem::copy_options::overwrite_existing, ec);
    }
    std::filesystem::current_path(work_dir);

    FDTD_Solver solver;
    solver.init(json_abs.filename().string(), false);
    if (max_steps >= 0) {
        solver.numSteps = max_steps;
    }
    solver.launch();
    const std::string case_name = extractCaseNameFromInput(json_abs.string());
    solver.end(case_name);
    const std::string probe_prefix = probeOutputPrefix(case_name);

    const ProbeSeries got_b = readProbeDat(probe_prefix + "before_Ex_3_3_1.dat");
    const ProbeSeries got_i = readProbeDat(probe_prefix + "inbox_Ex_3_3_3.dat");
    const ProbeSeries got_a = readProbeDat(probe_prefix + "after_Ex_3_3_5.dat");
    const ProbeSeries expected_b = readProbeDat(ref_before_abs.string());
    const ProbeSeries expected_i = readProbeDat(ref_inbox_abs.string());
    const ProbeSeries expected_a = readProbeDat(ref_after_abs.string());

    std::filesystem::current_path(old_cwd);

    if (got_b.field.empty() || got_i.field.empty() || got_a.field.empty() ||
        expected_b.field.empty() || expected_i.field.empty() || expected_a.field.empty()) {
        return 10;
    }

    int err = 0;
    if (compareProbeSeries("before", got_b, expected_b, max_steps) != 0) {
        err += 16;
    }
    if (compareProbeSeries("inbox", got_i, expected_i, max_steps) != 0) {
        err += 32;
    }
    if (compareProbeSeries("after", got_a, expected_a, max_steps) != 0) {
        err += 64;
    }

    if (correlation(got_i.field, got_i.incid) <= 0.999) {
        err += 1;
    }

    constexpr double atol = 6e-4;
    for (double v : got_b.field) {
        if (std::abs(v) > atol) {
            err += 2;
            break;
        }
    }
    for (double v : got_a.field) {
        if (std::abs(v) > atol) {
            err += 4;
            break;
        }
    }

    return err;
}

int test_run_pw_in_box_probe_files_exact(const std::string& json_path,
                                         const std::string& ref_before,
                                         const std::string& ref_inbox,
                                         const std::string& ref_after,
                                         int max_steps) {
    const std::filesystem::path json_abs = std::filesystem::absolute(json_path);
    const std::filesystem::path ref_before_abs = std::filesystem::absolute(ref_before);
    const std::filesystem::path ref_inbox_abs = std::filesystem::absolute(ref_inbox);
    const std::filesystem::path ref_after_abs = std::filesystem::absolute(ref_after);
    const std::filesystem::path case_dir = json_abs.parent_path();
    const std::filesystem::path old_cwd = std::filesystem::current_path();
    const std::filesystem::path work_dir = testWorkDir("semba_pw_in_box_exact_test");
    std::error_code ec;
    std::filesystem::remove_all(work_dir, ec);
    std::filesystem::create_directories(work_dir, ec);
    std::filesystem::copy_file(json_abs, work_dir / json_abs.filename(),
                               std::filesystem::copy_options::overwrite_existing, ec);
    const std::filesystem::path exc_src = case_dir / "gauss_1GHz.exc";
    if (std::filesystem::exists(exc_src)) {
        std::filesystem::copy_file(exc_src, work_dir / exc_src.filename(),
                                   std::filesystem::copy_options::overwrite_existing, ec);
    }
    std::filesystem::current_path(work_dir);

    FDTD_Solver solver;
    solver.init(json_abs.filename().string(), false);
    if (max_steps >= 0) {
        solver.numSteps = max_steps;
    }
    solver.launch();
    const std::string case_name = extractCaseNameFromInput(json_abs.string());
    solver.end(case_name);
    const std::string probe_prefix = probeOutputPrefix(case_name);

    const std::string got_b = probe_prefix + "before_Ex_3_3_1.dat";
    const std::string got_i = probe_prefix + "inbox_Ex_3_3_3.dat";
    const std::string got_a = probe_prefix + "after_Ex_3_3_5.dat";
    int err = 0;
    if (compareProbeFileExact("before", got_b, ref_before_abs.string(), max_steps) != 0) {
        err += 16;
    }
    if (compareProbeFileExact("inbox", got_i, ref_inbox_abs.string(), max_steps) != 0) {
        err += 32;
    }
    if (compareProbeFileExact("after", got_a, ref_after_abs.string(), max_steps) != 0) {
        err += 64;
    }

    std::filesystem::current_path(old_cwd);
    return err;
}

int test_run_holland_probe_output(const std::string& json_path,
                                  int max_steps) {
    const std::filesystem::path json_abs = std::filesystem::absolute(json_path);
    const std::filesystem::path case_dir = json_abs.parent_path();
    const std::filesystem::path old_cwd = std::filesystem::current_path();
    const std::filesystem::path work_dir =
        std::filesystem::temp_directory_path() / "semba_holland_probe_test";
    std::error_code ec;
    std::filesystem::remove_all(work_dir, ec);
    std::filesystem::create_directories(work_dir, ec);
    std::filesystem::copy_file(json_abs, work_dir / json_abs.filename(),
                               std::filesystem::copy_options::overwrite_existing, ec);
    const std::filesystem::path exc_src = case_dir / "holland.exc";
    if (std::filesystem::exists(exc_src)) {
        std::filesystem::copy_file(exc_src, work_dir / exc_src.filename(),
                                   std::filesystem::copy_options::overwrite_existing, ec);
    }
    std::filesystem::current_path(work_dir);

    FDTD_Solver solver;
    solver.init(json_abs.filename().string(), false);
    if (max_steps >= 0) {
        solver.numSteps = max_steps;
    }
    solver.launch();
    const std::string case_name = extractCaseNameFromInput(json_abs.string());
    solver.end(case_name);

    const std::string expected_name =
        probeOutputPrefix(case_name) + "mid_point_Wz_11_11_12_s8.dat";
    const std::vector<std::string> lines = readProbeLines(expected_name);
    std::filesystem::current_path(old_cwd);

    int err = 0;
    if (lines.empty()) {
        return 1;
    }
    if (max_steps >= 0 && lines.size() != static_cast<size_t>(max_steps + 2)) {
        err += 2;
    }
    const std::string expected_header =
        std::string("t              ") + expected_name +
        "       -E*dl Vplus Vminus Vplus-Vminus";
    if (lines.front() != expected_header) {
        err += 4;
    }
    return err;
}

int test_run_bulk_current_probe_output(const std::string& json_path,
                                       const std::string& expected_name,
                                       int max_steps) {
    const std::filesystem::path json_abs = std::filesystem::absolute(json_path);
    const std::filesystem::path case_dir = json_abs.parent_path();
    const std::filesystem::path old_cwd = std::filesystem::current_path();
    const std::filesystem::path work_dir =
        std::filesystem::temp_directory_path() / "semba_bulk_current_probe_test";
    std::error_code ec;
    std::filesystem::remove_all(work_dir, ec);
    std::filesystem::create_directories(work_dir, ec);
    std::filesystem::copy_file(json_abs, work_dir / json_abs.filename(),
                               std::filesystem::copy_options::overwrite_existing, ec);
    for (const auto& entry : std::filesystem::directory_iterator(case_dir)) {
        if (entry.path().extension() == ".exc") {
            std::filesystem::copy_file(entry.path(), work_dir / entry.path().filename(),
                                       std::filesystem::copy_options::overwrite_existing, ec);
        }
    }
    std::filesystem::current_path(work_dir);

    FDTD_Solver solver;
    solver.init(json_abs.filename().string(), false);
    if (max_steps >= 0) {
        solver.numSteps = max_steps;
    }
    solver.launch();
    solver.end(extractCaseNameFromInput(json_abs.string()));

    const std::vector<std::string> lines = readProbeLines(expected_name);
    std::filesystem::current_path(old_cwd);

    int err = 0;
    if (lines.empty()) {
        return 1;
    }
    if (max_steps >= 0 && lines.size() != static_cast<size_t>(max_steps + 2)) {
        err += 2;
    }
    const std::string expected_header = "t              " + expected_name;
    if (lines.front() != expected_header) {
        err += 4;
    }
    return err;
}

} // namespace SEMBA_FDTD_test
} // namespace SEMBA_FDTD_m

extern "C" {
    SEMBA_FDTD_m::semba_fdtd_t* create_semba_fdtd() { return new SEMBA_FDTD_m::semba_fdtd_t(); }
    void destroy_semba_fdtd(SEMBA_FDTD_m::semba_fdtd_t* p) { delete p; }
    void semba_fdtd_init(SEMBA_FDTD_m::semba_fdtd_t* p, const char* flags) { if (p) p->init(flags ? flags : ""); }
    void semba_fdtd_launch(SEMBA_FDTD_m::semba_fdtd_t* p) { if (p) p->launch(); }
    void semba_fdtd_end(SEMBA_FDTD_m::semba_fdtd_t* p, const char* case_name) {
        if (p) p->end(case_name ? case_name : "semba-fdtd");
    }
}
