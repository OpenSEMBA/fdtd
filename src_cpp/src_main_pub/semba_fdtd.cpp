#include "semba_fdtd.h"
#include "mapvtk_writer.h"

#include <nlohmann/json.hpp>

#ifdef CompileWithMTLN
#include "smbjson_m.h"
#include "wires_mtln_m.h"
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
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>

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

fdtd_real flushFortranSubnormal(fdtd_real value) {
    if (value != static_cast<fdtd_real>(0.0) &&
        std::abs(value) < std::numeric_limits<fdtd_real>::min()) {
        return static_cast<fdtd_real>(0.0);
    }
    return value;
}

fdtd_real fortranGridInverse(fdtd_real value) {
    const fdtd_real inverse = static_cast<fdtd_real>(1.0) / value;
#ifdef CompileWithReal8
    return inverse;
#else
    return std::nextafter(inverse, -std::numeric_limits<fdtd_real>::infinity());
#endif
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

struct entrada_t {
    int layoutnumber = 0; int ierr = 0;
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
    int cellI = 1, cellJ = 1, cellK = 1;
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
    int nd = 0;
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
    double current = 0.0;
    double currentpast = 0.0;
    double qplus_qminus = 0.0;
};

struct HollandWireNode_t {
    int i = 0, j = 0, k = 0;
    double chargePresent = 0.0;
    double chargePast = 0.0;
    double ctePlain = 0.0;
    double cteProp = 1.0;
    std::vector<int> currentPlus;
    std::vector<int> currentMinus;
};

struct HollandWireProbe_t {
    std::string name;
    int segmentIndex = -1;
    int cellI = 0, cellJ = 0, cellK = 0;
    int direction = 3;
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
            pd.sources.planeWaves.push_back(s);
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
                        p.cellI = coordPos[coord_id][0];
                        p.cellJ = coordPos[coord_id][1];
                        p.cellK = coordPos[coord_id][2];
                    }
                } else if (coordPos.count(elem_id)) {
                    p.cellI = coordPos[elem_id][0];
                    p.cellJ = coordPos[elem_id][1];
                    p.cellK = coordPos[elem_id][2];
                }
            }
            p.fieldByDir.resize(p.directions.size());
            p.incidentByDir.resize(p.directions.size());
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
    Parseador_t pd;
    int NX = 10, NY = 10, NZ = 10;
    double dt = 1e-12, dx = 0.01, dy = 0.01, dz = 0.01;
    double eps0 = EPS0, mu0 = MU0;
    std::vector<fdtd_real> Ex, Ey, Ez, Hx, Hy, Hz;
    std::vector<int> pecMask;
    std::vector<fdtd_real> CeE, CmH;
    std::vector<source_t> sources;
    std::map<std::string, ExcitationData> excitations;
    std::vector<PlaneWaveState_t> planeWaves;
    std::vector<probe_output_t> probes;
    std::vector<BulkCurrentProbe_t> bulkCurrentProbes;
    std::vector<NodalCurrentSegment_t> nodalCurrentSegments;
    std::vector<HollandWireSegment_t> hollandSegments;
    std::vector<HollandWireNode_t> hollandNodes;
    std::vector<HollandWireProbe_t> hollandProbes;
    bool still_planewave_time = true;
    bool planewave_switched_off = false;
    bool useMur = true;
    bool usePec = false;
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

    void init(const std::string& filename, bool map_vtk = false) {
        inputFile = filename;
        createMapVtk = map_vtk;
        std::ifstream jf(filename);
        if (jf.is_open()) {
            jf >> inputRoot;
            jf.close();
        }
        pd = parseFDTDJSON(filename);
        NX = pd.general.XI; NY = pd.general.YI; NZ = pd.general.ZI;
        if (NX <= 0) NX = 10; if (NY <= 0) NY = 10; if (NZ <= 0) NZ = 10;
        dt = pd.general.dt; if (dt <= 0.0) dt = 1e-12;
        if (!pd.cellSteps.cellStepsX.empty()) dx = pd.cellSteps.cellStepsX[0]; else dx = 0.025;
        if (!pd.cellSteps.cellStepsY.empty()) dy = pd.cellSteps.cellStepsY[0]; else dy = 0.025;
        if (!pd.cellSteps.cellStepsZ.empty()) dz = pd.cellSteps.cellStepsZ[0]; else dz = 0.025;
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

        int ex_n = ex_nx()*ex_ny()*ex_nz();
        int ey_n = ey_nx()*ey_ny()*ey_nz();
        int ez_n = ez_nx()*ez_ny()*ez_nz();
        int hx_n = hx_nx()*hx_ny()*hx_nz();
        int hy_n = hy_nx()*hy_ny()*hy_nz();
        int hz_n = hz_nx()*hz_ny()*hz_nz();
        Ex.resize(ex_n,0); Ey.resize(ey_n,0); Ez.resize(ez_n,0);
        Hx.resize(hx_n,0); Hy.resize(hy_n,0); Hz.resize(hz_n,0);
        int max_n = std::max({ex_n,ey_n,ez_n,hx_n,hy_n,hz_n});
        pecMask.resize(max_n, 0);
        CeE.resize(max_n, 0.0); CmH.resize(max_n, 0.0);

        const fdtd_real ce = static_cast<fdtd_real>(
            dt / static_cast<double>(static_cast<fdtd_real>(eps0)));
        const fdtd_real ch = static_cast<fdtd_real>(
            dt / static_cast<double>(static_cast<fdtd_real>(mu0)));
        for (int i = 0; i < max_n; i++) { CeE[i] = ce; CmH[i] = ch; }

        sources = pd.sources.planeWaves;
        const std::filesystem::path json_dir =
            std::filesystem::path(filename).parent_path();
        for (auto& src : sources) {
            if (src.magnitudeFile.empty()) continue;
            std::string exc_path = src.magnitudeFile;
            if (!std::filesystem::exists(exc_path) && !json_dir.empty()) {
                exc_path = (json_dir / src.magnitudeFile).string();
            }
            excitations[src.magnitudeFile] = readExcitationFile(exc_path);
        }
        planeWaves.resize(sources.size());
        for (int i = 0; i < (int)sources.size(); i++) {
            planeWaves[i].px.resize(1,0); planeWaves[i].py.resize(1,0); planeWaves[i].pz.resize(1,0);
            planeWaves[i].ex.resize(1,0); planeWaves[i].ey.resize(1,0); planeWaves[i].ez.resize(1,0);
            planeWaves[i].hx.resize(1,0); planeWaves[i].hy.resize(1,0); planeWaves[i].hz.resize(1,0);
            initPlaneWave(i);
        }
        probes = pd.probes.probes;
        initBulkCurrentProbes();
        initNodalCurrentSources();
        initHollandWires();

        useMur = true;
        usePec = false;
        if (!pd.boundaries.boundaries.empty()) {
            const std::string btype = pd.boundaries.boundaries[0].type;
            useMur = (btype == "mur" || btype == "MUR" ||
                      btype == "pml" || btype == "PML");
            usePec = (btype == "pec" || btype == "PEC");
        }
        initMurBorders();

        std::cout << "FDTD: grid=" << NX << "x" << NY << "x" << NZ << " dt=" << dt << " steps=" << numSteps << std::endl;
    }

    void calcMurConstants() {
        const fdtd_real one = static_cast<fdtd_real>(1.0);
        const fdtd_real cluz = static_cast<fdtd_real>(1.0) /
            std::sqrt(static_cast<fdtd_real>(eps0) * static_cast<fdtd_real>(mu0));
        const auto cab1 = [this, one, cluz](fdtd_real inv_step) {
            const fdtd_real cnum = static_cast<fdtd_real>(
                static_cast<double>(one / inv_step) /
                (dt * static_cast<double>(cluz)));
            return (one - cnum) / (one + cnum);
        };
        backCab1 = cab1(fortranGridInverse(static_cast<fdtd_real>(dx)));
        frontCab1 = cab1(fortranGridInverse(static_cast<fdtd_real>(dx)));
        leftCab1 = cab1(fortranGridInverse(static_cast<fdtd_real>(dy)));
        rightCab1 = cab1(fortranGridInverse(static_cast<fdtd_real>(dy)));
        downCab1 = cab1(fortranGridInverse(static_cast<fdtd_real>(dz)));
        upCab1 = cab1(fortranGridInverse(static_cast<fdtd_real>(dz)));
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

    int ex_nx() const { return NX + 3; }
    int ex_ny() const { return NY + 4; }
    int ex_nz() const { return NZ + 4; }
    int ey_nx() const { return NX + 4; }
    int ey_ny() const { return NY + 3; }
    int ey_nz() const { return NZ + 4; }
    int ez_nx() const { return NX + 4; }
    int ez_ny() const { return NY + 4; }
    int ez_nz() const { return NZ + 3; }
    int hx_nx() const { return NX + 4; }
    int hx_ny() const { return NY + 3; }
    int hx_nz() const { return NZ + 3; }
    int hy_nx() const { return NX + 3; }
    int hy_ny() const { return NY + 4; }
    int hy_nz() const { return NZ + 3; }
    int hz_nx() const { return NX + 3; }
    int hz_ny() const { return NY + 3; }
    int hz_nz() const { return NZ + 4; }

    int ex_idx(int i,int j,int k) const { return (i + 2)*ex_ny()*ex_nz() + (j + 2)*ex_nz() + (k + 2); }
    int ey_idx(int i,int j,int k) const { return (i + 2)*ey_ny()*ey_nz() + (j + 2)*ey_nz() + (k + 2); }
    int ez_idx(int i,int j,int k) const { return (i + 2)*ez_ny()*ez_nz() + (j + 2)*ez_nz() + (k + 2); }
    int hx_idx(int i,int j,int k) const { return (i + 2)*hx_ny()*hx_nz() + (j + 2)*hx_nz() + (k + 2); }
    int hy_idx(int i,int j,int k) const { return (i + 2)*hy_ny()*hy_nz() + (j + 2)*hy_nz() + (k + 2); }
    int hz_idx(int i,int j,int k) const { return (i + 2)*hz_ny()*hz_nz() + (j + 2)*hz_nz() + (k + 2); }

    bool in_ex(int i,int j,int k) const { return i >= -2 && i <= NX && j >= -2 && j <= NY + 1 && k >= -2 && k <= NZ + 1; }
    bool in_ey(int i,int j,int k) const { return i >= -2 && i <= NX + 1 && j >= -2 && j <= NY && k >= -2 && k <= NZ + 1; }
    bool in_ez(int i,int j,int k) const { return i >= -2 && i <= NX + 1 && j >= -2 && j <= NY + 1 && k >= -2 && k <= NZ; }
    bool in_hx(int i,int j,int k) const { return i >= -2 && i <= NX + 1 && j >= -2 && j <= NY && k >= -2 && k <= NZ; }
    bool in_hy(int i,int j,int k) const { return i >= -2 && i <= NX && j >= -2 && j <= NY + 1 && k >= -2 && k <= NZ; }
    bool in_hz(int i,int j,int k) const { return i >= -2 && i <= NX && j >= -2 && j <= NY && k >= -2 && k <= NZ + 1; }

    fdtd_real lineX1(int n) const { return static_cast<fdtd_real>(n) * static_cast<fdtd_real>(dx); }
    fdtd_real lineY1(int n) const { return static_cast<fdtd_real>(n) * static_cast<fdtd_real>(dy); }
    fdtd_real lineZ1(int n) const { return static_cast<fdtd_real>(n) * static_cast<fdtd_real>(dz); }
    fdtd_real idxe1(int i) const { (void)i; return fortranGridInverse(static_cast<fdtd_real>(dx)); }
    fdtd_real idye1(int j) const { (void)j; return fortranGridInverse(static_cast<fdtd_real>(dy)); }
    fdtd_real idze1(int k) const { (void)k; return fortranGridInverse(static_cast<fdtd_real>(dz)); }
    fdtd_real idxh1(int i) const { (void)i; return fortranGridInverse(static_cast<fdtd_real>(dx)); }
    fdtd_real idyh1(int j) const { (void)j; return fortranGridInverse(static_cast<fdtd_real>(dy)); }
    fdtd_real idzh1(int k) const { (void)k; return fortranGridInverse(static_cast<fdtd_real>(dz)); }

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
                    if (a < b) return std::pair<int, int>{a, b - 1};
                    if (a == b) return std::pair<int, int>{a, b};
                    return std::pair<int, int>{b, a - 1};
                };
                const auto [x1, x2] = convertInterval(iv[0][0].get<int>(), iv[1][0].get<int>());
                const auto [y1, y2] = convertInterval(iv[0][1].get<int>(), iv[1][1].get<int>());
                const auto [z1, z2] = convertInterval(iv[0][2].get<int>(), iv[1][2].get<int>());
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

        const fdtd_real px = pw.px[0], py = pw.py[0], pz = pw.pz[0];
        fdtd_real xd0 = 0.0, yd0 = 0.0, zd0 = 0.0;
        const int xi = 1, xe = NX, yi = 1, ye = NY, zi = 1, ze = NZ;
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

    double hxValue0(int i, int j, int k) const {
        return in_hx(i, j, k) ? static_cast<double>(Hx[hx_idx(i, j, k)]) : 0.0;
    }

    double hyValue0(int i, int j, int k) const {
        return in_hy(i, j, k) ? static_cast<double>(Hy[hy_idx(i, j, k)]) : 0.0;
    }

    double hzValue0(int i, int j, int k) const {
        return in_hz(i, j, k) ? static_cast<double>(Hz[hz_idx(i, j, k)]) : 0.0;
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

    static int segmentLo(const NodalCurrentSegment_t& segment, int axis) {
        if (axis == 0) return segment.xi;
        if (axis == 1) return segment.yi;
        return segment.zi;
    }

    static int segmentHi(const NodalCurrentSegment_t& segment, int axis) {
        if (axis == 0) return segment.xe;
        if (axis == 1) return segment.ye;
        return segment.ze;
    }

    bool segmentCrossesProbe(const NodalCurrentSegment_t& segment,
                             const BulkCurrentProbe_t& probe) const {
        if (segment.direction != probe.direction) return false;
        const int axis = axisFromDirection(probe.direction);
        if (!containsInclusive(segmentLo(segment, axis), segmentHi(segment, axis),
                               probeLo(probe, axis))) {
            return false;
        }
        for (int d = 0; d < 3; ++d) {
            if (d == axis) continue;
            const int coord = segmentLo(segment, d);
            if (!containsInclusive(probeLo(probe, d), probeHi(probe, d), coord)) {
                return false;
            }
        }
        return true;
    }

    double sampleNodalCurrentContribution(const BulkCurrentProbe_t& probe) const {
        double current = 0.0;
        for (const auto& segment : nodalCurrentSegments) {
            if (!segmentCrossesProbe(segment, probe)) continue;
            const auto exc = excitations.find(segment.magnitudeFile);
            if (exc == excitations.end()) continue;
            current += static_cast<double>(segment.sign) *
                getExcitationValue(exc->second, currentTime);
        }
        return static_cast<double>(probe.sign) * current;
    }

    double sampleBulkCurrentValue(const BulkCurrentProbe_t& probe) const {
        double current = 0.0;
        if (probe.direction == 'x') {
            const int i = probe.xi - 1;
            for (int k = probe.zi; k <= probe.ze; ++k) {
                current += dz * (hzValue0(i, probe.ye - 1, k - 1) -
                                 hzValue0(i, probe.yi - 2, k - 1));
            }
            for (int j = probe.yi; j <= probe.ye; ++j) {
                current -= dy * (hyValue0(i, j - 1, probe.ze - 1) -
                                 hyValue0(i, j - 1, probe.zi - 2));
            }
        } else if (probe.direction == 'y') {
            const int j = probe.yi - 1;
            for (int i = probe.xi; i <= probe.xe; ++i) {
                current += dx * (hxValue0(i - 1, j, probe.ze - 1) -
                                 hxValue0(i - 1, j, probe.zi - 2));
            }
            for (int k = probe.zi; k <= probe.ze; ++k) {
                current -= dz * (hzValue0(probe.xe - 1, j, k - 1) -
                                 hzValue0(probe.xi - 2, j, k - 1));
            }
        } else {
            const int k = probe.zi - 1;
            for (int j = probe.yi; j <= probe.ye; ++j) {
                current += dy * (hyValue0(probe.xe - 1, j - 1, k) -
                                 hyValue0(probe.xi - 2, j - 1, k));
            }
            for (int i = probe.xi; i <= probe.xe; ++i) {
                current -= dx * (hxValue0(i - 1, probe.ye - 1, k) -
                                 hxValue0(i - 1, probe.yi - 2, k));
            }
        }
        return static_cast<double>(probe.sign) * current +
            sampleNodalCurrentContribution(probe);
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

    void writeBulkCurrentProbeOutputs(const std::string& caseName) {
        for (const auto& probe : bulkCurrentProbes) {
            std::string fullname = probeOutputPrefix(caseName) + probe.name + "_" +
                bulkDirectionTag(probe.direction) + "_" +
                std::to_string(probe.xi) + "_" + std::to_string(probe.yi) + "_" +
                std::to_string(probe.zi) + "__" +
                std::to_string(probe.xe) + "_" + std::to_string(probe.ye) + "_" +
                std::to_string(probe.ze) + ".dat";
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
        if (direction == 1) return dx;
        if (direction == 2) return dy;
        return dz;
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

    void addToHollandField(const HollandWireSegment_t& segment, double value) {
        const int i = segment.i - 1;
        const int j = segment.j - 1;
        const int k = segment.k - 1;
        if (segment.direction == 1 && in_ex(i, j, k)) {
            Ex[ex_idx(i, j, k)] = static_cast<fdtd_real>(static_cast<double>(Ex[ex_idx(i, j, k)]) + value);
        } else if (segment.direction == 2 && in_ey(i, j, k)) {
            Ey[ey_idx(i, j, k)] = static_cast<fdtd_real>(static_cast<double>(Ey[ey_idx(i, j, k)]) + value);
        } else if (segment.direction == 3 && in_ez(i, j, k)) {
            Ez[ez_idx(i, j, k)] = static_cast<fdtd_real>(static_cast<double>(Ez[ez_idx(i, j, k)]) + value);
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
        const double invMu = 1.0 / mu0;
        double lind = (1.0 / (4.0 * PI * invMu)) *
            (std::log((deltaTransv1 * deltaTransv1 + deltaTransv2 * deltaTransv2) /
                      (4.0 * radius * radius)) +
             deltaTransv1 / deltaTransv2 * std::atan(deltaTransv2 / deltaTransv1) +
             deltaTransv2 / deltaTransv1 * std::atan(deltaTransv1 / deltaTransv2) +
             PI * radius * radius / (deltaTransv2 * deltaTransv1) - 3.0);
        if (radius < 0.3 * deltaTransv1 || radius < 0.3 * deltaTransv2) {
            lind -= 0.57 / (4.0 * PI * invMu);
        }
        if (radius > 0.3 * deltaTransv1 || radius > 0.3 * deltaTransv2) {
            lind /= (1.0 - PI * radius * radius / (deltaTransv1 * deltaTransv2));
        }
        return lind;
    }

    void finishHollandConstants() {
        const double invMu = 1.0 / mu0;
        const double invEps = 1.0 / eps0;
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
        }
        for (auto& seg : hollandSegments) {
            seg.lind = hollandSelfInductance(seg.radius, seg.deltaTransv1, seg.deltaTransv2) +
                       seg.inductance;
            const double denom = seg.lind / dt + seg.resistance * 0.5;
            seg.cte1 = (seg.lind / dt - seg.resistance * 0.5) / denom;
            seg.cte3 = invMu * invEps / seg.delta * seg.lind / denom;
            seg.cte2 = 1.0 / denom;
            seg.cte5 = g2 / (seg.deltaTransv1 * seg.deltaTransv2);
        }
    }

    void initHollandWires() {
        hollandSegments.clear();
        hollandNodes.clear();
        hollandProbes.clear();
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
        std::map<int, WireMaterial> wireMaterials;
        for (const auto& mat : inputRoot["materials"]) {
            if (mat.value("type", std::string()) != "wire") continue;
            WireMaterial wm;
            wm.isWire = true;
            wm.radius = mat.value("radius", 0.0);
            wm.resistance = mat.value("resistancePerMeter", 0.0);
            wm.inductance = mat.value("inductancePerMeter", 0.0);
            wireMaterials[mat.value("id", 0)] = wm;
        }
        if (wireMaterials.empty()) return;

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

        std::map<std::string, int> nodeByCoord;
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
                    const int nCells = std::abs(deltaCells);
                    for (int s = 0; s < nCells; ++s) {
                        std::array<int, 3> minus = p0;
                        minus[axis] = (sign > 0) ? p0[axis] + s : p0[axis] - s - 1;
                        std::array<int, 3> plus = minus;
                        plus[axis] += 1;
                        HollandWireSegment_t seg;
                        seg.i = minus[0];
                        seg.j = minus[1];
                        seg.k = minus[2];
                        seg.direction = axis + 1;
                        seg.nd = static_cast<int>(hollandSegments.size()) + 3;
                        seg.radius = matIt->second.radius;
                        seg.resistance = matIt->second.resistance;
                        seg.inductance = matIt->second.inductance;
                        seg.delta = hollandStepForDirection(seg.direction);
                        if (seg.direction == 1) {
                            seg.deltaTransv1 = dy;
                            seg.deltaTransv2 = dz;
                        } else if (seg.direction == 2) {
                            seg.deltaTransv1 = dz;
                            seg.deltaTransv2 = dx;
                        } else {
                            seg.deltaTransv1 = dx;
                            seg.deltaTransv2 = dy;
                        }
                        seg.chargeMinus = getOrCreateHollandNode(nodeByCoord, minus);
                        seg.chargePlus = getOrCreateHollandNode(nodeByCoord, plus);
                        const int segIdx = static_cast<int>(hollandSegments.size());
                        hollandNodes[static_cast<size_t>(seg.chargeMinus)].currentPlus.push_back(segIdx);
                        hollandNodes[static_cast<size_t>(seg.chargePlus)].currentMinus.push_back(segIdx);
                        hollandSegments.push_back(seg);
                    }
                }
            }
        }
        if (hollandSegments.empty()) return;
        finishHollandConstants();

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
            wp.segmentIndex = bestSeg;
            wp.cellI = p[0];
            wp.cellJ = p[1];
            wp.cellK = p[2];
            wp.direction = seg.direction;
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
            node.chargePresent = node.cteProp * node.chargePast -
                                 node.ctePlain * (iPlus - iMinus);
        }
        for (const auto& seg : hollandSegments) {
            addToHollandField(seg, -seg.cte5 * seg.current);
        }
        for (auto& seg : hollandSegments) {
            seg.currentpast = seg.current;
            const auto& qPlus = hollandNodes[static_cast<size_t>(seg.chargePlus)];
            const auto& qMinus = hollandNodes[static_cast<size_t>(seg.chargeMinus)];
            seg.qplus_qminus = qPlus.chargePresent - qMinus.chargePresent;
            seg.current = seg.cte1 * seg.current -
                          seg.cte3 * seg.qplus_qminus +
                          seg.cte2 * hollandFieldValue(seg);
        }
    }

    void sampleHollandProbes() {
        if (hollandProbes.empty()) return;
        const double invMuInvEps = (1.0 / mu0) * (1.0 / eps0);
        for (auto& probe : hollandProbes) {
            if (probe.segmentIndex < 0 ||
                probe.segmentIndex >= static_cast<int>(hollandSegments.size())) {
                continue;
            }
            const auto& seg = hollandSegments[static_cast<size_t>(probe.segmentIndex)];
            const auto& qPlus = hollandNodes[static_cast<size_t>(seg.chargePlus)];
            const auto& qMinus = hollandNodes[static_cast<size_t>(seg.chargeMinus)];
            const double eTimesDl = -hollandFieldValue(seg) * seg.delta;
            const double vplus = ((qPlus.chargePresent + qPlus.chargePast) * 0.5) *
                                 seg.lind * invMuInvEps;
            const double vminus = ((qMinus.chargePresent + qMinus.chargePast) * 0.5) *
                                  seg.lind * invMuInvEps;
            probe.timeData.push_back(currentTime);
            probe.currentData.push_back(seg.currentpast);
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

    void writeHollandProbeOutputs(const std::string& caseName) {
        for (const auto& probe : hollandProbes) {
            if (probe.segmentIndex < 0) continue;
            std::string fullname = probeOutputPrefix(caseName) + probe.name + "_" +
                hollandDirectionTag(probe.direction) +
                std::to_string(probe.cellI) + "_" + std::to_string(probe.cellJ) + "_" +
                std::to_string(probe.cellK) + "_s" + std::to_string(probe.nd) + ".dat";
            std::ofstream out(fullname);
            out << "t              " << fullname
                << "       -E*dl Vplus Vminus Vplus-Vminus\n";
            for (size_t t = 0; t < probe.timeData.size(); ++t) {
                const bool hasDelayedSample = t >= static_cast<size_t>(probe.delaySteps);
                const size_t src = hasDelayedSample
                                       ? t - static_cast<size_t>(probe.delaySteps)
                                       : 0;
                const double current = hasDelayedSample ? probe.currentData[src] : 0.0;
                const double eTimesDl = hasDelayedSample ? probe.eTimesDlData[src] : -0.0;
                const double vplus = hasDelayedSample ? probe.vplusData[src] : 0.0;
                const double vminus = hasDelayedSample ? probe.vminusData[src] : 0.0;
                const double vdrop = hasDelayedSample ? probe.vdropData[src] : 0.0;
                out << formatFortranE(probe.timeData[t], 27, 17)
                    << formatFortranE(current, 19, 9)
                    << formatFortranE(eTimesDl, 19, 9)
                    << formatFortranE(vplus, 19, 9)
                    << formatFortranE(vminus, 19, 9)
                    << formatFortranE(vdrop, 19, 9)
                    << "\n";
            }
        }
    }

    void advanceE() {
        for (int i = -1; i < NX - 1; ++i)
            for (int j = -1; j < NY; ++j)
                for (int k = -1; k < NZ; ++k) {
                    const int idx = ex_idx(i, j, k);
                    if (usePec && (j == NY - 1 || k == NZ - 1)) {
                        Ex[idx] = 0.0;
                        continue;
                    }
                    Ex[idx] = Ex[idx] + CeE[idx] *
                        ((Hz[hz_idx(i, j, k)] - Hz[hz_idx(i, j - 1, k)]) * idyh1(j) -
                         (Hy[hy_idx(i, j, k)] - Hy[hy_idx(i, j, k - 1)]) * idzh1(k));
                }
        for (int i = -1; i < NX; ++i)
            for (int j = -1; j < NY - 1; ++j)
                for (int k = -1; k < NZ; ++k) {
                    const int idx = ey_idx(i, j, k);
                    if (usePec && (i == NX - 1 || k == NZ - 1)) {
                        Ey[idx] = 0.0;
                        continue;
                    }
                    Ey[idx] = Ey[idx] + CeE[idx] *
                        ((Hx[hx_idx(i, j, k)] - Hx[hx_idx(i, j, k - 1)]) * idzh1(k) -
                         (Hz[hz_idx(i, j, k)] - Hz[hz_idx(i - 1, j, k)]) * idxh1(i));
                }
        for (int i = -1; i < NX; ++i)
            for (int j = -1; j < NY; ++j)
                for (int k = -1; k < NZ - 1; ++k) {
                    const int idx = ez_idx(i, j, k);
                    if (usePec && (i == NX - 1 || j == NY - 1)) {
                        Ez[idx] = 0.0;
                        continue;
                    }
                    Ez[idx] = Ez[idx] + CeE[idx] *
                        ((Hy[hy_idx(i, j, k)] - Hy[hy_idx(i - 1, j, k)]) * idxh1(i) -
                         (Hx[hx_idx(i, j, k)] - Hx[hx_idx(i, j - 1, k)]) * idyh1(j));
                }
    }

    void applyMurE() {
        (void)useMur;
    }

    void applyPecE() {
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
            bnd = static_cast<fdtd_real>(past_int + cab * (int_now - past_bnd));
        };

        // Back (x min): Fortran MURc uses sweepXI-1, one plane below the H sweep.
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
        // Front (x max): Fortran first-order Mur writes MURc%XE, one plane above the H sweep.
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
        // Left (y min): Fortran MURc uses sweepYI-1, one plane below the H sweep.
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
        // Right (y max): Fortran first-order Mur writes MURc%YE, one plane above the H sweep.
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
        // Down (z min): Fortran MURc uses sweepZI-1, one plane below the H sweep.
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
        // Up (z max): Fortran first-order Mur writes MURc%ZE, one plane above the H sweep.
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

        // Store past fields (Fortran AdvanceMagneticMUR tail).
        for (int j = -1; j <= NY - 1; ++j) {
            for (int k = -1; k <= NZ - 2; ++k) {
                const size_t p = static_cast<size_t>((j + 1) * NZ + (k + 1));
                murPastHyBack[p] = Hy[hy_idx(-2, j, k)];
                murPastHyBackInt[p] = Hy[hy_idx(-1, j, k)];
                murPastHyFront[p] = Hy[hy_idx(NX - 1, j, k)];
                murPastHyFrontInt[p] = Hy[hy_idx(NX - 2, j, k)];
            }
        }
        for (int j = -1; j <= NY - 2; ++j) {
            for (int k = -1; k <= NZ - 1; ++k) {
                const size_t p = static_cast<size_t>((j + 1) * (NZ + 1) + (k + 1));
                murPastHzBack[p] = Hz[hz_idx(-2, j, k)];
                murPastHzBackInt[p] = Hz[hz_idx(-1, j, k)];
                murPastHzFront[p] = Hz[hz_idx(NX - 1, j, k)];
                murPastHzFrontInt[p] = Hz[hz_idx(NX - 2, j, k)];
            }
        }
        for (int i = -1; i <= NX - 1; ++i) {
            for (int k = -1; k <= NZ - 2; ++k) {
                const size_t p = static_cast<size_t>((i + 1) * NZ + (k + 1));
                murPastHxLeft[p] = Hx[hx_idx(i, -2, k)];
                murPastHxLeftInt[p] = Hx[hx_idx(i, -1, k)];
                murPastHxRight[p] = Hx[hx_idx(i, NY - 1, k)];
                murPastHxRightInt[p] = Hx[hx_idx(i, NY - 2, k)];
            }
        }
        for (int i = -1; i <= NX - 2; ++i) {
            for (int k = -1; k <= NZ - 1; ++k) {
                const size_t p = static_cast<size_t>((i + 1) * (NZ + 1) + (k + 1));
                murPastHzLeft[p] = Hz[hz_idx(i, -2, k)];
                murPastHzLeftInt[p] = Hz[hz_idx(i, -1, k)];
                murPastHzRight[p] = Hz[hz_idx(i, NY - 1, k)];
                murPastHzRightInt[p] = Hz[hz_idx(i, NY - 2, k)];
            }
        }
        for (int i = -1; i <= NX - 2; ++i) {
            for (int j = -1; j <= NY - 1; ++j) {
                const size_t p = static_cast<size_t>((i + 1) * (NY + 1) + (j + 1));
                murPastHyDown[p] = Hy[hy_idx(i, j, -2)];
                murPastHyDownInt[p] = Hy[hy_idx(i, j, -1)];
                murPastHyUp[p] = Hy[hy_idx(i, j, NZ - 1)];
                murPastHyUpInt[p] = Hy[hy_idx(i, j, NZ - 2)];
            }
        }
        for (int i = -1; i <= NX - 1; ++i) {
            for (int j = -1; j <= NY - 2; ++j) {
                const size_t p = static_cast<size_t>((i + 1) * NY + (j + 1));
                murPastHxDown[p] = Hx[hx_idx(i, j, -2)];
                murPastHxDownInt[p] = Hx[hx_idx(i, j, -1)];
                murPastHxUp[p] = Hx[hx_idx(i, j, NZ - 1)];
                murPastHxUpInt[p] = Hx[hx_idx(i, j, NZ - 2)];
            }
        }
    }

    void advanceH() {
        for (int i = -1; i < NX; ++i)
            for (int j = -1; j < NY - 1; ++j)
                for (int k = -1; k < NZ - 1; ++k) {
                    const int idx = hx_idx(i, j, k);
                    if (usePec && i == NX - 1) {
                        Hx[idx] = 0.0;
                        continue;
                    }
                    Hx[idx] = Hx[idx] + CmH[idx] *
                        ((Ey[ey_idx(i, j, k + 1)] - Ey[ey_idx(i, j, k)]) * idze1(k) -
                         (Ez[ez_idx(i, j + 1, k)] - Ez[ez_idx(i, j, k)]) * idye1(j));
                }
        for (int i = -1; i < NX - 1; ++i)
            for (int j = -1; j < NY; ++j)
                for (int k = -1; k < NZ - 1; ++k) {
                    const int idx = hy_idx(i, j, k);
                    if (usePec && j == NY - 1) {
                        Hy[idx] = 0.0;
                        continue;
                    }
                    Hy[idx] = Hy[idx] + CmH[idx] *
                        ((Ez[ez_idx(i + 1, j, k)] - Ez[ez_idx(i, j, k)]) * idxe1(i) -
                         (Ex[ex_idx(i, j, k + 1)] - Ex[ex_idx(i, j, k)]) * idze1(k));
                }
        for (int i = -1; i < NX - 1; ++i)
            for (int j = -1; j < NY - 1; ++j)
                for (int k = -1; k < NZ; ++k) {
                    const int idx = hz_idx(i, j, k);
                    if (usePec && k == NZ - 1) {
                        Hz[idx] = 0.0;
                        continue;
                    }
                    Hz[idx] = Hz[idx] + CmH[idx] *
                        ((Ex[ex_idx(i, j + 1, k)] - Ex[ex_idx(i, j, k)]) * idye1(j) -
                         (Ey[ey_idx(i + 1, j, k)] - Ey[ey_idx(i, j, k)]) * idxe1(i));
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
        for (int pwIdx = 0; pwIdx < (int)planeWaves.size(); ++pwIdx) {
            const auto& pw = planeWaves[pwIdx];
            if (pw.iluminaTr) {
                int i = std::max(XI, pw.esqx1);
                fdtd_real id = static_cast<fdtd_real>(idxh1(i));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2 - 1); ++k) {
                    for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2); ++j) {
                        const fdtd_real inc = computeIncid(pwIdx, 4, currentTime, i - 1, j, k);
                        Ez[ez_idx(i - 1, j - 1, k - 1)] -= G2_1 * inc * id;
                    }
                }
                i = std::max(XI, pw.esqx1);
                id = static_cast<fdtd_real>(idxh1(i));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2 - 1); ++j) {
                        const fdtd_real inc = computeIncid(pwIdx, 5, currentTime, i - 1, j, k);
                        Ey[ey_idx(i - 1, j - 1, k - 1)] += G2_1 * inc * id;
                    }
                }
            }
            if (pw.iluminaFr) {
                int i = std::min(XE, pw.esqx2);
                fdtd_real id = static_cast<fdtd_real>(idxh1(i));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2 - 1); ++k) {
                    for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2); ++j) {
                        const fdtd_real inc = computeIncid(pwIdx, 4, currentTime, i, j, k);
                        Ez[ez_idx(i - 1, j - 1, k - 1)] += G2_1 * inc * id;
                    }
                }
                i = std::min(XE, pw.esqx2);
                id = static_cast<fdtd_real>(idxh1(i));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2 - 1); ++j) {
                        const fdtd_real inc = computeIncid(pwIdx, 5, currentTime, i, j, k);
                        Ey[ey_idx(i - 1, j - 1, k - 1)] -= G2_1 * inc * id;
                    }
                }
            }
            if (pw.iluminaIz) {
                int j = std::max(YI, pw.esqy1);
                fdtd_real id = static_cast<fdtd_real>(idyh1(j));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2 - 1); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 5, currentTime, i, j - 1, k);
                        Ex[ex_idx(i - 1, j - 1, k - 1)] -= G2_1 * inc * id;
                    }
                }
                j = std::max(YI, pw.esqy1);
                id = static_cast<fdtd_real>(idyh1(j));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2 - 1); ++k) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 3, currentTime, i, j - 1, k);
                        Ez[ez_idx(i - 1, j - 1, k - 1)] += G2_1 * inc * id;
                    }
                }
            }
            if (pw.iluminaDe) {
                int j = std::min(YE, pw.esqy2);
                fdtd_real id = static_cast<fdtd_real>(idyh1(j));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2 - 1); ++k) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 3, currentTime, i, j, k);
                        Ez[ez_idx(i - 1, j - 1, k - 1)] -= G2_1 * inc * id;
                    }
                }
                j = std::min(YE, pw.esqy2);
                id = static_cast<fdtd_real>(idyh1(j));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2 - 1); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 5, currentTime, i, j, k);
                        Ex[ex_idx(i - 1, j - 1, k - 1)] += G2_1 * inc * id;
                    }
                }
            }
            if (pw.iluminaAb) {
                int k = std::max(ZI, pw.esqz1);
                fdtd_real id = static_cast<fdtd_real>(idzh1(k));
                for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2); ++j) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2 - 1); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 4, currentTime, i, j, k - 1);
                        Ex[ex_idx(i - 1, j - 1, k - 1)] += G2_1 * inc * id;
                    }
                }
                k = std::max(ZI, pw.esqz1);
                id = static_cast<fdtd_real>(idzh1(k));
                for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2 - 1); ++j) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 3, currentTime, i, j, k - 1);
                        Ey[ey_idx(i - 1, j - 1, k - 1)] -= G2_1 * inc * id;
                    }
                }
            }
            if (pw.iluminaAr) {
                int k = std::min(ZE, pw.esqz2);
                fdtd_real id = static_cast<fdtd_real>(idzh1(k));
                for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2); ++j) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2 - 1); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 4, currentTime, i, j, k);
                        Ex[ex_idx(i - 1, j - 1, k - 1)] -= G2_1 * inc * id;
                    }
                }
                k = std::min(ZE, pw.esqz2);
                id = static_cast<fdtd_real>(idzh1(k));
                for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2 - 1); ++j) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 3, currentTime, i, j, k);
                        Ey[ey_idx(i - 1, j - 1, k - 1)] += G2_1 * inc * id;
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
        const double timeH = currentTime + 0.5 * dt;
        for (int pwIdx = 0; pwIdx < (int)planeWaves.size(); ++pwIdx) {
            const auto& pw = planeWaves[pwIdx];
            if (pw.iluminaTr) {
                const int i = std::max(XI, pw.esqx1) - 1;
                const fdtd_real id = static_cast<fdtd_real>(idxe1(i));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2 - 1); ++j) {
                        const fdtd_real inc = computeIncid(pwIdx, 1, timeH, i + 1, j, k);
                        Hz[hz_idx(i - 1, j - 1, k - 1)] += Gm2_1 * inc * id;
                    }
                }
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2 - 1); ++k) {
                    for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2); ++j) {
                        const fdtd_real inc = computeIncid(pwIdx, 2, timeH, i + 1, j, k);
                        Hy[hy_idx(i - 1, j - 1, k - 1)] -= Gm2_1 * inc * id;
                    }
                }
            }
            if (pw.iluminaFr) {
                const int i = std::min(XE, pw.esqx2);
                const fdtd_real id = static_cast<fdtd_real>(idxe1(i));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2 - 1); ++j) {
                        const fdtd_real inc = computeIncid(pwIdx, 1, timeH, i, j, k);
                        Hz[hz_idx(i - 1, j - 1, k - 1)] -= Gm2_1 * inc * id;
                    }
                }
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2 - 1); ++k) {
                    for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2); ++j) {
                        const fdtd_real inc = computeIncid(pwIdx, 2, timeH, i, j, k);
                        Hy[hy_idx(i - 1, j - 1, k - 1)] += Gm2_1 * inc * id;
                    }
                }
            }
            if (pw.iluminaIz) {
                const int jHx = std::max(YI, pw.esqy1) - 1;
                const fdtd_real idHx = static_cast<fdtd_real>(idye1(jHx));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2 - 1); ++k) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 2, timeH, i, jHx + 1, k);
                        Hx[hx_idx(i - 1, jHx - 1, k - 1)] += Gm2_1 * inc * idHx;
                    }
                }
                const int jHz = std::max(YI, pw.esqy1) - 1;
                const fdtd_real idHz = static_cast<fdtd_real>(idye1(jHz));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2 - 1); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 0, timeH, i, jHz + 1, k);
                        Hz[hz_idx(i - 1, jHz - 1, k - 1)] -= Gm2_1 * inc * idHz;
                    }
                }
            }
            if (pw.iluminaDe) {
                const int j = std::min(YE, pw.esqy2);
                const fdtd_real id = static_cast<fdtd_real>(idye1(j));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2 - 1); ++k) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 2, timeH, i, j, k);
                        Hx[hx_idx(i - 1, j - 1, k - 1)] -= Gm2_1 * inc * id;
                    }
                }
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2 - 1); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 0, timeH, i, j, k);
                        Hz[hz_idx(i - 1, j - 1, k - 1)] += Gm2_1 * inc * id;
                    }
                }
            }
            if (pw.iluminaAb) {
                const int k = std::max(ZI, pw.esqz1) - 1;
                const fdtd_real id = static_cast<fdtd_real>(idze1(k));
                for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2 - 1); ++j) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 1, timeH, i, j, k + 1);
                        Hx[hx_idx(i - 1, j - 1, k - 1)] -= Gm2_1 * inc * id;
                    }
                }
                for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2); ++j) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2 - 1); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 0, timeH, i, j, k + 1);
                        Hy[hy_idx(i - 1, j - 1, k - 1)] += Gm2_1 * inc * id;
                    }
                }
            }
            if (pw.iluminaAr) {
                const int k = std::min(ZE, pw.esqz2);
                const fdtd_real id = static_cast<fdtd_real>(idze1(k));
                for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2 - 1); ++j) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 1, timeH, i, j, k);
                        Hx[hx_idx(i - 1, j - 1, k - 1)] += Gm2_1 * inc * id;
                    }
                }
                for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2); ++j) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2 - 1); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 0, timeH, i, j, k);
                        Hy[hy_idx(i - 1, j - 1, k - 1)] -= Gm2_1 * inc * id;
                    }
                }
            }
        }
    }

    void applyPlaneWaveSource() {}

    void sampleProbes() {
        sampleBulkCurrentProbes();
        sampleHollandProbes();
        for (auto& probe : probes) {
            if (probe.domainType != "time" || probe.directions.empty()) continue;
            const int ci = probe.cellI, cj = probe.cellJ, ck = probe.cellK;
            const int i = ci - 1, j = cj - 1, k = ck - 1;
            double Ex_v = 0, Ey_v = 0, Ez_v = 0;
            double Hx_v = 0, Hy_v = 0, Hz_v = 0;
            if (in_ex(i, j, k))
                Ex_v = Ex[ex_idx(i, j, k)];
            if (in_ey(i, j, k))
                Ey_v = Ey[ey_idx(i, j, k)];
            if (in_ez(i, j, k))
                Ez_v = Ez[ez_idx(i, j, k)];
            if (in_hx(i, j, k))
                Hx_v = Hx[hx_idx(i, j, k)];
            if (in_hy(i, j, k))
                Hy_v = Hy[hy_idx(i, j, k)];
            if (in_hz(i, j, k))
                Hz_v = Hz[hz_idx(i, j, k)];
            double inc_x = 0.0, inc_y = 0.0, inc_z = 0.0;
            for (int pwIdx = 0; pwIdx < (int)planeWaves.size(); ++pwIdx) {
                if (probe.field == "magnetic") {
                    const double timeH = currentTime + 0.5 * dt;
                    inc_x += computeIncid(pwIdx, 3, timeH, ci, cj, ck, true);
                    inc_y += computeIncid(pwIdx, 4, timeH, ci, cj, ck, true);
                    inc_z += computeIncid(pwIdx, 5, timeH, ci, cj, ck, true);
                } else {
                    inc_x += computeIncid(pwIdx, 0, currentTime, ci, cj, ck, true);
                    inc_y += computeIncid(pwIdx, 1, currentTime, ci, cj, ck, true);
                    inc_z += computeIncid(pwIdx, 2, currentTime, ci, cj, ck, true);
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

    void writeProbeOutputs(const std::string& caseName) {
        writeBulkCurrentProbeOutputs(caseName);
        writeHollandProbeOutputs(caseName);
        for (auto& probe : probes) {
            if (probe.domainType != "time") continue;
            std::string fname = probeOutputPrefix(caseName) + probe.name + "_";
            std::string fieldDir = (probe.field == "magnetic") ? "H" : "E";
            for (size_t d = 0; d < probe.directions.size(); ++d) {
                std::string fullname = fname + fieldDir + probe.directions[d] + "_";
                fullname += std::to_string(probe.cellI) + "_" + std::to_string(probe.cellJ) + "_" +
                            std::to_string(probe.cellK) + ".dat";
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

    void stepMurFdtd() {
        advanceE();
        advanceH();
        applyMurH();
        step += 1;
        currentTime = step * dt;
    }

    void step_once() { advanceH(); }

    void end(const std::string& caseName) {
        writeProbeOutputs(caseName);
        if (createMapVtk && !inputRoot.is_null()) {
            mapvtk::writeMapVtkFromJson(caseName, inputRoot);
        }
        std::cout << "Output files written." << std::endl;
    }

    void launch() {
        still_planewave_time = true;
        planewave_switched_off = false;
        std::cout << "Running FDTD: " << numSteps << " steps..." << std::endl;
        for (step = 0; step <= numSteps; step++) {
            currentTime = step * dt;
            flushPlanewaveOff();
            advanceE();
            advanceHollandWiresE();
            applyPecE();
            advancePlaneWaveE();
            applyPecE();
            advanceH();
            applyPecH();
            advancePlaneWaveH();
            applyPecH();
            applyMurH();
            sampleProbes();
            if (step % 500 == 0 || step == numSteps)
                std::cout << "  Step " << step << "/" << numSteps << " (t=" << currentTime << "s)" << std::endl;
        }
        std::cout << "Simulation complete." << std::endl;
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
    impl_->solver.init(filename, cli_mapvtk || json_mapvtk);
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
    const std::filesystem::path work_dir =
        std::filesystem::temp_directory_path() / "semba_pw_in_box_test";
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
    const std::filesystem::path work_dir =
        std::filesystem::temp_directory_path() / "semba_pw_in_box_exact_test";
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
        "t              " + expected_name + "       -E*dl Vplus Vminus Vplus-Vminus";
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
