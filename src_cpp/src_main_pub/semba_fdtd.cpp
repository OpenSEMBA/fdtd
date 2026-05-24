#include "semba_fdtd.h"
#include "mapvtk_writer.h"

#include <nlohmann/json.hpp>
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

const double PI = 3.14159265358979323846;
const double EPS0 = 8.854187817e-12;
const double MU0 = 1.2566370614e-6;
const double C0 = 299792458.0;
const double ZVAC = std::sqrt(MU0/EPS0);
const double heurCFL = 0.8;
const int BUFSIZE = 1024;

#ifdef CompileWithReal8
using fdtd_real = double;
#else
using fdtd_real = float;
#endif

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
struct mtln_t { int numWires = 0; };
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
    std::string type, name, magnitudeFile;
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
    struct { int XI = 0, XE = 0, YI = 0, YE = 0, ZI = 0, ZE = 0; int NumMedia = 0; double dt = 0.0; } general;
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
    bool still_planewave_time = true;
    bool planewave_switched_off = false;
    bool useMur = true;
    bool usePec = false;
    // MURc zones (InitMURBorders bordersmur.F90 L94-137) — thin 1-cell face pads per component.
    struct MurZone { int xi = 0, xe = 0, yi = 0, ye = 0, zi = 0, ze = 0; };
    MurZone murHyBack_, murHyFront_, murHzBack_, murHzFront_;
    MurZone murHxLeft_, murHxRight_, murHzLeft_, murHzRight_;
    MurZone murHyDown_, murHyUp_, murHxDown_, murHxUp_;
    double backCab1 = 0.0, frontCab1 = 0.0, leftCab1 = 0.0, rightCab1 = 0.0;
    double downCab1 = 0.0, upCab1 = 0.0;
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

        for (int i = 0; i < max_n; i++) { CeE[i] = dt/eps0; CmH[i] = dt/mu0; }

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

        useMur = true;
        usePec = false;
        if (!pd.boundaries.boundaries.empty()) {
            const std::string btype = pd.boundaries.boundaries[0].type;
            useMur = (btype == "mur" || btype == "MUR");
            usePec = (btype == "pec" || btype == "PEC");
        }
        initMurBorders();

        std::cout << "FDTD: grid=" << NX << "x" << NY << "x" << NZ << " dt=" << dt << " steps=" << numSteps << std::endl;
    }

    void calcMurConstants() {
        const auto cab1 = [this](double cell_step) {
            const double cnum = cell_step / (dt * C0);
            return (1.0 - cnum) / (1.0 + cnum);
        };
        backCab1 = cab1(dx);
        frontCab1 = cab1(dx);
        leftCab1 = cab1(dy);
        rightCab1 = cab1(dy);
        downCab1 = cab1(dz);
        upCab1 = cab1(dz);
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
        murPastHyBack.assign(NY * (NZ + 1), 0.0);
        murPastHyBackInt.assign(NY * (NZ + 1), 0.0);
        murPastHzBack.assign((NY + 1) * NZ, 0.0);
        murPastHzBackInt.assign((NY + 1) * NZ, 0.0);
        murPastHyFront.assign(NY * (NZ + 1), 0.0);
        murPastHyFrontInt.assign(NY * (NZ + 1), 0.0);
        murPastHzFront.assign((NY + 1) * NZ, 0.0);
        murPastHzFrontInt.assign((NY + 1) * NZ, 0.0);
        murPastHxLeft.assign(NX * (NZ + 1), 0.0);
        murPastHxLeftInt.assign(NX * (NZ + 1), 0.0);
        murPastHzLeft.assign((NX + 1) * NZ, 0.0);
        murPastHzLeftInt.assign((NX + 1) * NZ, 0.0);
        murPastHxRight.assign(NX * (NZ + 1), 0.0);
        murPastHxRightInt.assign(NX * (NZ + 1), 0.0);
        murPastHzRight.assign((NX + 1) * NZ, 0.0);
        murPastHzRightInt.assign((NX + 1) * NZ, 0.0);
        murPastHyDown.assign((NX + 1) * NY, 0.0);
        murPastHyDownInt.assign((NX + 1) * NY, 0.0);
        murPastHxDown.assign(NX * (NY + 1), 0.0);
        murPastHxDownInt.assign(NX * (NY + 1), 0.0);
        murPastHyUp.assign((NX + 1) * NY, 0.0);
        murPastHyUpInt.assign((NX + 1) * NY, 0.0);
        murPastHxUp.assign(NX * (NY + 1), 0.0);
        murPastHxUpInt.assign(NX * (NY + 1), 0.0);
    }

    int ex_nx() const { return NX + 2; }
    int ex_ny() const { return NY + 3; }
    int ex_nz() const { return NZ + 3; }
    int ey_nx() const { return NX + 3; }
    int ey_ny() const { return NY + 2; }
    int ey_nz() const { return NZ + 3; }
    int ez_nx() const { return NX + 3; }
    int ez_ny() const { return NY + 3; }
    int ez_nz() const { return NZ + 2; }
    int hx_nx() const { return NX + 3; }
    int hx_ny() const { return NY + 2; }
    int hx_nz() const { return NZ + 2; }
    int hy_nx() const { return NX + 2; }
    int hy_ny() const { return NY + 3; }
    int hy_nz() const { return NZ + 2; }
    int hz_nx() const { return NX + 2; }
    int hz_ny() const { return NY + 2; }
    int hz_nz() const { return NZ + 3; }

    int ex_idx(int i,int j,int k) const { return (i + 1)*ex_ny()*ex_nz() + (j + 1)*ex_nz() + (k + 1); }
    int ey_idx(int i,int j,int k) const { return (i + 1)*ey_ny()*ey_nz() + (j + 1)*ey_nz() + (k + 1); }
    int ez_idx(int i,int j,int k) const { return (i + 1)*ez_ny()*ez_nz() + (j + 1)*ez_nz() + (k + 1); }
    int hx_idx(int i,int j,int k) const { return (i + 1)*hx_ny()*hx_nz() + (j + 1)*hx_nz() + (k + 1); }
    int hy_idx(int i,int j,int k) const { return (i + 1)*hy_ny()*hy_nz() + (j + 1)*hy_nz() + (k + 1); }
    int hz_idx(int i,int j,int k) const { return (i + 1)*hz_ny()*hz_nz() + (j + 1)*hz_nz() + (k + 1); }

    bool in_ex(int i,int j,int k) const { return i >= -1 && i <= NX && j >= -1 && j <= NY + 1 && k >= -1 && k <= NZ + 1; }
    bool in_ey(int i,int j,int k) const { return i >= -1 && i <= NX + 1 && j >= -1 && j <= NY && k >= -1 && k <= NZ + 1; }
    bool in_ez(int i,int j,int k) const { return i >= -1 && i <= NX + 1 && j >= -1 && j <= NY + 1 && k >= -1 && k <= NZ; }
    bool in_hx(int i,int j,int k) const { return i >= -1 && i <= NX + 1 && j >= -1 && j <= NY && k >= -1 && k <= NZ; }
    bool in_hy(int i,int j,int k) const { return i >= -1 && i <= NX && j >= -1 && j <= NY + 1 && k >= -1 && k <= NZ; }
    bool in_hz(int i,int j,int k) const { return i >= -1 && i <= NX && j >= -1 && j <= NY && k >= -1 && k <= NZ + 1; }

    double lineX1(int n) const { return (n - 1) * dx; }
    double lineY1(int n) const { return (n - 1) * dy; }
    double lineZ1(int n) const { return (n - 1) * dz; }
    double idxh1(int i) const { return 1.0 / dx; }
    double idyh1(int j) const { return 1.0 / dy; }
    double idzh1(int k) const { return 1.0 / dz; }

    void physCoord1(int nfield, int i, int j, int k, double& xf, double& yf, double& zf) const {
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
        return pw.samples[static_cast<size_t>(nprev)] +
               (pw.samples[static_cast<size_t>(nprev + 1)] - pw.samples[static_cast<size_t>(nprev)]) *
                   (t_frac / pw.deltaevol);
    }

    // i,j,k are 1-based Yee indices (Fortran convention).
    fdtd_real computeIncid(int pwIdx, int nfield, double time, int i, int j, int k) {
        auto& pw = planeWaves[pwIdx];
        double xd = 0.0, yd = 0.0, zd = 0.0;
        physCoord1(nfield, i, j, k, xd, yd, zd);
        const fdtd_real xf = static_cast<fdtd_real>(xd);
        const fdtd_real yf = static_cast<fdtd_real>(yd);
        const fdtd_real zf = static_cast<fdtd_real>(zd);
        const fdtd_real timef = static_cast<fdtd_real>(time);
        const fdtd_real cluz = static_cast<fdtd_real>(C0);
        const fdtd_real d = (xf * pw.px[0] + yf * pw.py[0] + zf * pw.pz[0]) - pw.distanciaInicial;
        const fdtd_real value = evolucion(pwIdx, timef - d / cluz);
        switch (nfield) {
            case 0: return value * pw.ex[0];
            case 1: return value * pw.ey[0];
            case 2: return value * pw.ez[0];
            case 3: return value * pw.hx[0];
            case 4: return value * pw.hy[0];
            case 5: return value * pw.hz[0];
            default: return 0.0;
        }
    }

    void initPlaneWave(int srcIdx) {
        auto& src = sources[srcIdx];
        if (src.type != "planewave") return;
        auto& pw = planeWaves[srcIdx];
        pw.px[0] = std::sin(src.direction.theta) * std::cos(src.direction.phi);
        pw.py[0] = std::sin(src.direction.theta) * std::sin(src.direction.phi);
        pw.pz[0] = std::cos(src.direction.theta);
        const double modu = std::sqrt(pw.px[0] * pw.px[0] + pw.py[0] * pw.py[0] + pw.pz[0] * pw.pz[0]);
        if (modu > 0.0) {
            pw.px[0] /= modu;
            pw.py[0] /= modu;
            pw.pz[0] /= modu;
        }
        pw.ex[0] = std::sin(src.polarization.theta) * std::cos(src.polarization.phi);
        pw.ey[0] = std::sin(src.polarization.theta) * std::sin(src.polarization.phi);
        pw.ez[0] = std::cos(src.polarization.theta);
        pw.hx[0] = (pw.py[0]*pw.ez[0] - pw.pz[0]*pw.ey[0]) / ZVAC;
        pw.hy[0] = (pw.pz[0]*pw.ex[0] - pw.px[0]*pw.ez[0]) / ZVAC;
        pw.hz[0] = (pw.px[0]*pw.ey[0] - pw.py[0]*pw.ex[0]) / ZVAC;
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

        const double px = pw.px[0], py = pw.py[0], pz = pw.pz[0];
        double xd0 = 0.0, yd0 = 0.0, zd0 = 0.0;
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

    void advanceE() {
        const fdtd_real ce = static_cast<fdtd_real>(dt / eps0);
        const fdtd_real inv_dx = static_cast<fdtd_real>(1.0 / dx);
        const fdtd_real inv_dy = static_cast<fdtd_real>(1.0 / dy);
        const fdtd_real inv_dz = static_cast<fdtd_real>(1.0 / dz);
        for (int i = 0; i < NX; ++i)
            for (int j = 0; j <= NY; ++j)
                for (int k = 0; k <= NZ; ++k) {
                    const int idx = ex_idx(i, j, k);
                    const fdtd_real dHz_dy =
                        Hz[hz_idx(i, j, k)] - Hz[hz_idx(i, j - 1, k)];
                    const fdtd_real dHy_dz =
                        Hy[hy_idx(i, j, k)] - Hy[hy_idx(i, j, k - 1)];
                    Ex[idx] += ce * (dHz_dy * inv_dy - dHy_dz * inv_dz);
                }
        for (int i = 0; i <= NX; ++i)
            for (int j = 0; j < NY; ++j)
                for (int k = 0; k <= NZ; ++k) {
                    const int idx = ey_idx(i, j, k);
                    const fdtd_real dHx_dz =
                        Hx[hx_idx(i, j, k)] - Hx[hx_idx(i, j, k - 1)];
                    const fdtd_real dHz_dx =
                        Hz[hz_idx(i, j, k)] - Hz[hz_idx(i - 1, j, k)];
                    Ey[idx] += ce * (dHx_dz * inv_dz - dHz_dx * inv_dx);
                }
        for (int i = 0; i <= NX; ++i)
            for (int j = 0; j <= NY; ++j)
                for (int k = 0; k < NZ; ++k) {
                    const int idx = ez_idx(i, j, k);
                    const fdtd_real dHy_dx =
                        Hy[hy_idx(i, j, k)] - Hy[hy_idx(i - 1, j, k)];
                    const fdtd_real dHx_dy =
                        Hx[hx_idx(i, j, k)] - Hx[hx_idx(i, j - 1, k)];
                    Ez[idx] += ce * (dHy_dx * inv_dx - dHx_dy * inv_dy);
                }
    }

    void applyMurE() {
        (void)useMur;
    }

    void applyPecE() {
        if (!usePec) return;

        // PEC border media have zero E-update coefficients in Fortran
        // (healer.F90 marks only components tangent to each outer face).
        for (int i = 0; i < NX; ++i) {
            for (int k = 0; k <= NZ; ++k) {
                Ex[ex_idx(i, 0, k)] = 0.0;
                Ex[ex_idx(i, NY, k)] = 0.0;
            }
            for (int j = 0; j <= NY; ++j) {
                Ex[ex_idx(i, j, 0)] = 0.0;
                Ex[ex_idx(i, j, NZ)] = 0.0;
            }
        }

        for (int j = 0; j < NY; ++j) {
            for (int k = 0; k <= NZ; ++k) {
                Ey[ey_idx(0, j, k)] = 0.0;
                Ey[ey_idx(NX, j, k)] = 0.0;
            }
        }
        for (int i = 0; i <= NX; ++i) {
            for (int j = 0; j < NY; ++j) {
                Ey[ey_idx(i, j, 0)] = 0.0;
                Ey[ey_idx(i, j, NZ)] = 0.0;
            }
        }

        for (int j = 0; j <= NY; ++j) {
            for (int k = 0; k < NZ; ++k) {
                Ez[ez_idx(0, j, k)] = 0.0;
                Ez[ez_idx(NX, j, k)] = 0.0;
            }
        }
        for (int i = 0; i <= NX; ++i) {
            for (int k = 0; k < NZ; ++k) {
                Ez[ez_idx(i, 0, k)] = 0.0;
                Ez[ez_idx(i, NY, k)] = 0.0;
            }
        }
    }

    void applyPecH() {
        if (!usePec) return;

        for (int j = 0; j < NY; ++j) {
            for (int k = 0; k < NZ; ++k) {
                Hx[hx_idx(0, j, k)] = 0.0;
                Hx[hx_idx(NX, j, k)] = 0.0;
            }
        }
        for (int i = 0; i < NX; ++i) {
            for (int k = 0; k < NZ; ++k) {
                Hy[hy_idx(i, 0, k)] = 0.0;
                Hy[hy_idx(i, NY, k)] = 0.0;
            }
        }
        for (int i = 0; i < NX; ++i) {
            for (int j = 0; j < NY; ++j) {
                Hz[hz_idx(i, j, 0)] = 0.0;
                Hz[hz_idx(i, j, NZ)] = 0.0;
            }
        }
    }

    bool hasPlaneWaveSource() const {
        return std::any_of(sources.begin(), sources.end(), [](const source_t& src) {
            return src.type == "planewave";
        });
    }

    void applyMurH() {
        if (!useMur) return;
        const bool skipLowerCompactMur = hasPlaneWaveSource();
        auto mur_face = [](fdtd_real& bnd, fdtd_real int_now, fdtd_real past_int,
                           fdtd_real past_bnd, double cab) {
            bnd = static_cast<fdtd_real>(past_int + cab * (int_now - past_bnd));
        };

        // Back (x min): Hy and Hz at i=0
        if (!skipLowerCompactMur) {
            for (int j = 0; j < NY; ++j) {
                for (int k = 0; k <= NZ; ++k) {
                    const size_t p = static_cast<size_t>(j * (NZ + 1) + k);
                    const int idx0 = hy_idx(0, j, k);
                    const int idx1 = hy_idx(1, j, k);
                    mur_face(Hy[idx0], Hy[idx1], murPastHyBackInt[p], murPastHyBack[p], backCab1);
                }
            }
            for (int j = 0; j <= NY; ++j) {
                for (int k = 0; k < NZ; ++k) {
                    const size_t p = static_cast<size_t>(j * NZ + k);
                    const int idx0 = hz_idx(0, j, k);
                    const int idx1 = hz_idx(1, j, k);
                    mur_face(Hz[idx0], Hz[idx1], murPastHzBackInt[p], murPastHzBack[p], backCab1);
                }
            }
        }
        // Front (x max): Hy and Hz at i=NX
        for (int j = 0; j < NY; ++j) {
            for (int k = 0; k <= NZ; ++k) {
                const size_t p = static_cast<size_t>(j * (NZ + 1) + k);
                const int idxN = hy_idx(NX, j, k);
                const int idxI = hy_idx(NX - 1, j, k);
                mur_face(Hy[idxN], Hy[idxI], murPastHyFrontInt[p], murPastHyFront[p], frontCab1);
            }
        }
        for (int j = 0; j <= NY; ++j) {
            for (int k = 0; k < NZ; ++k) {
                const size_t p = static_cast<size_t>(j * NZ + k);
                const int idxN = hz_idx(NX, j, k);
                const int idxI = hz_idx(NX - 1, j, k);
                mur_face(Hz[idxN], Hz[idxI], murPastHzFrontInt[p], murPastHzFront[p], frontCab1);
            }
        }
        // Left (y min): Hx and Hz at j=0
        if (!skipLowerCompactMur) {
            for (int i = 0; i < NX; ++i) {
                for (int k = 0; k <= NZ; ++k) {
                    const size_t p = static_cast<size_t>(i * (NZ + 1) + k);
                    const int idx0 = hx_idx(i, 0, k);
                    const int idx1 = hx_idx(i, 1, k);
                    mur_face(Hx[idx0], Hx[idx1], murPastHxLeftInt[p], murPastHxLeft[p], leftCab1);
                }
            }
            for (int i = 0; i <= NX; ++i) {
                for (int k = 0; k < NZ; ++k) {
                    const size_t p = static_cast<size_t>(i * NZ + k);
                    const int idx0 = hz_idx(i, 0, k);
                    const int idx1 = hz_idx(i, 1, k);
                    mur_face(Hz[idx0], Hz[idx1], murPastHzLeftInt[p], murPastHzLeft[p], leftCab1);
                }
            }
        }
        // Right (y max): Hx and Hz at j=NY
        for (int i = 0; i < NX; ++i) {
            for (int k = 0; k <= NZ; ++k) {
                const size_t p = static_cast<size_t>(i * (NZ + 1) + k);
                const int idxN = hx_idx(i, NY, k);
                const int idxI = hx_idx(i, NY - 1, k);
                mur_face(Hx[idxN], Hx[idxI], murPastHxRightInt[p], murPastHxRight[p], rightCab1);
            }
        }
        for (int i = 0; i <= NX; ++i) {
            for (int k = 0; k < NZ; ++k) {
                const size_t p = static_cast<size_t>(i * NZ + k);
                const int idxN = hz_idx(i, NY, k);
                const int idxI = hz_idx(i, NY - 1, k);
                mur_face(Hz[idxN], Hz[idxI], murPastHzRightInt[p], murPastHzRight[p], rightCab1);
            }
        }
        // Down (z min): Hy and Hx at k=0
        if (!skipLowerCompactMur) {
            for (int i = 0; i <= NX; ++i) {
                for (int j = 0; j < NY; ++j) {
                    const size_t p = static_cast<size_t>(i * NY + j);
                    const int idx0 = hy_idx(i, j, 0);
                    const int idx1 = hy_idx(i, j, 1);
                    mur_face(Hy[idx0], Hy[idx1], murPastHyDownInt[p], murPastHyDown[p], downCab1);
                }
            }
            for (int i = 0; i < NX; ++i) {
                for (int j = 0; j <= NY; ++j) {
                    const size_t p = static_cast<size_t>(i * (NY + 1) + j);
                    const int idx0 = hx_idx(i, j, 0);
                    const int idx1 = hx_idx(i, j, 1);
                    mur_face(Hx[idx0], Hx[idx1], murPastHxDownInt[p], murPastHxDown[p], downCab1);
                }
            }
        }
        // Up (z max): Hy and Hx at k=NZ
        for (int i = 0; i <= NX; ++i) {
            for (int j = 0; j < NY; ++j) {
                const size_t p = static_cast<size_t>(i * NY + j);
                const int idxN = hy_idx(i, j, NZ);
                const int idxI = hy_idx(i, j, NZ - 1);
                mur_face(Hy[idxN], Hy[idxI], murPastHyUpInt[p], murPastHyUp[p], upCab1);
            }
        }
        for (int i = 0; i < NX; ++i) {
            for (int j = 0; j <= NY; ++j) {
                const size_t p = static_cast<size_t>(i * (NY + 1) + j);
                const int idxN = hx_idx(i, j, NZ);
                const int idxI = hx_idx(i, j, NZ - 1);
                mur_face(Hx[idxN], Hx[idxI], murPastHxUpInt[p], murPastHxUp[p], upCab1);
            }
        }

        // Store past fields (Fortran AdvanceMagneticMUR tail).
        for (int j = 0; j < NY; ++j) {
            for (int k = 0; k <= NZ; ++k) {
                const size_t p = static_cast<size_t>(j * (NZ + 1) + k);
                murPastHyBack[p] = Hy[hy_idx(0, j, k)];
                murPastHyBackInt[p] = Hy[hy_idx(1, j, k)];
                murPastHyFront[p] = Hy[hy_idx(NX, j, k)];
                murPastHyFrontInt[p] = Hy[hy_idx(NX - 1, j, k)];
            }
        }
        for (int j = 0; j <= NY; ++j) {
            for (int k = 0; k < NZ; ++k) {
                const size_t p = static_cast<size_t>(j * NZ + k);
                murPastHzBack[p] = Hz[hz_idx(0, j, k)];
                murPastHzBackInt[p] = Hz[hz_idx(1, j, k)];
                murPastHzFront[p] = Hz[hz_idx(NX, j, k)];
                murPastHzFrontInt[p] = Hz[hz_idx(NX - 1, j, k)];
            }
        }
        for (int i = 0; i < NX; ++i) {
            for (int k = 0; k <= NZ; ++k) {
                const size_t p = static_cast<size_t>(i * (NZ + 1) + k);
                murPastHxLeft[p] = Hx[hx_idx(i, 0, k)];
                murPastHxLeftInt[p] = Hx[hx_idx(i, 1, k)];
                murPastHxRight[p] = Hx[hx_idx(i, NY, k)];
                murPastHxRightInt[p] = Hx[hx_idx(i, NY - 1, k)];
            }
        }
        for (int i = 0; i <= NX; ++i) {
            for (int k = 0; k < NZ; ++k) {
                const size_t p = static_cast<size_t>(i * NZ + k);
                murPastHzLeft[p] = Hz[hz_idx(i, 0, k)];
                murPastHzLeftInt[p] = Hz[hz_idx(i, 1, k)];
                murPastHzRight[p] = Hz[hz_idx(i, NY, k)];
                murPastHzRightInt[p] = Hz[hz_idx(i, NY - 1, k)];
            }
        }
        for (int i = 0; i <= NX; ++i) {
            for (int j = 0; j < NY; ++j) {
                const size_t p = static_cast<size_t>(i * NY + j);
                murPastHyDown[p] = Hy[hy_idx(i, j, 0)];
                murPastHyDownInt[p] = Hy[hy_idx(i, j, 1)];
                murPastHyUp[p] = Hy[hy_idx(i, j, NZ)];
                murPastHyUpInt[p] = Hy[hy_idx(i, j, NZ - 1)];
            }
        }
        for (int i = 0; i < NX; ++i) {
            for (int j = 0; j <= NY; ++j) {
                const size_t p = static_cast<size_t>(i * (NY + 1) + j);
                murPastHxDown[p] = Hx[hx_idx(i, j, 0)];
                murPastHxDownInt[p] = Hx[hx_idx(i, j, 1)];
                murPastHxUp[p] = Hx[hx_idx(i, j, NZ)];
                murPastHxUpInt[p] = Hx[hx_idx(i, j, NZ - 1)];
            }
        }
    }

    void advanceH() {
        const fdtd_real ch = static_cast<fdtd_real>(dt / mu0);
        const fdtd_real inv_dx = static_cast<fdtd_real>(1.0 / dx);
        const fdtd_real inv_dy = static_cast<fdtd_real>(1.0 / dy);
        const fdtd_real inv_dz = static_cast<fdtd_real>(1.0 / dz);
        for (int i = 0; i <= NX; ++i)
            for (int j = 0; j < NY; ++j)
                for (int k = 0; k < NZ; ++k) {
                    const int idx = hx_idx(i, j, k);
                    const fdtd_real dEz_dy = Ez[ez_idx(i, j + 1, k)] - Ez[ez_idx(i, j, k)];
                    const fdtd_real dEy_dz = Ey[ey_idx(i, j, k + 1)] - Ey[ey_idx(i, j, k)];
                    Hx[idx] += ch * (dEy_dz * inv_dz - dEz_dy * inv_dy);
                }
        for (int i = 0; i < NX; ++i)
            for (int j = 0; j <= NY; ++j)
                for (int k = 0; k < NZ; ++k) {
                    const int idx = hy_idx(i, j, k);
                    const fdtd_real dEx_dz = Ex[ex_idx(i, j, k + 1)] - Ex[ex_idx(i, j, k)];
                    const fdtd_real dEz_dx = Ez[ez_idx(i + 1, j, k)] - Ez[ez_idx(i, j, k)];
                    Hy[idx] += ch * (dEz_dx * inv_dx - dEx_dz * inv_dz);
                }
        for (int i = 0; i < NX; ++i)
            for (int j = 0; j < NY; ++j)
                for (int k = 0; k <= NZ; ++k) {
                    const int idx = hz_idx(i, j, k);
                    const fdtd_real dEy_dx = Ey[ey_idx(i + 1, j, k)] - Ey[ey_idx(i, j, k)];
                    const fdtd_real dEx_dy = Ex[ex_idx(i, j + 1, k)] - Ex[ex_idx(i, j, k)];
                    Hz[idx] += ch * (dEx_dy * inv_dy - dEy_dx * inv_dx);
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
        const fdtd_real G2_1 = static_cast<fdtd_real>(dt / eps0);
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
        const fdtd_real Gm2_1 = static_cast<fdtd_real>(dt / mu0);
        const double timeH = currentTime + 0.5 * dt;
        for (int pwIdx = 0; pwIdx < (int)planeWaves.size(); ++pwIdx) {
            const auto& pw = planeWaves[pwIdx];
            if (pw.iluminaTr) {
                const int i = std::max(XI, pw.esqx1) - 1;
                const fdtd_real id = static_cast<fdtd_real>(idxh1(i + 1));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2 - 1); ++j) {
                        const fdtd_real inc = computeIncid(pwIdx, 4, timeH, i + 1, j, k);
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
                const fdtd_real id = static_cast<fdtd_real>(idxh1(i));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2 - 1); ++j) {
                        const fdtd_real inc = computeIncid(pwIdx, 4, timeH, i, j, k);
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
                const fdtd_real idHx = static_cast<fdtd_real>(idyh1(jHx + 1));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2 - 1); ++k) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 2, timeH, i, jHx + 1, k);
                        Hx[hx_idx(i - 1, jHx - 1, k - 1)] += Gm2_1 * inc * idHx;
                    }
                }
                const int jHz = std::max(YI, pw.esqy1) - 1;
                const fdtd_real idHz = static_cast<fdtd_real>(idyh1(jHz + 1));
                for (int k = std::max(ZI, pw.esqz1); k <= std::min(ZE, pw.esqz2); ++k) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2 - 1); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 0, timeH, i, jHz + 1, k);
                        Hz[hz_idx(i - 1, jHz - 1, k - 1)] -= Gm2_1 * inc * idHz;
                    }
                }
            }
            if (pw.iluminaDe) {
                const int j = std::min(YE, pw.esqy2);
                const fdtd_real id = static_cast<fdtd_real>(idyh1(j));
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
                const fdtd_real id = static_cast<fdtd_real>(idzh1(k + 1));
                for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2 - 1); ++j) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 4, timeH, i, j, k + 1);
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
                const fdtd_real id = static_cast<fdtd_real>(idzh1(k));
                for (int j = std::max(YI, pw.esqy1); j <= std::min(YE, pw.esqy2 - 1); ++j) {
                    for (int i = std::max(XI, pw.esqx1); i <= std::min(XE, pw.esqx2); ++i) {
                        const fdtd_real inc = computeIncid(pwIdx, 4, timeH, i, j, k);
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
        for (auto& probe : probes) {
            if (probe.domainType != "time" || probe.directions.empty()) continue;
            const int ci = probe.cellI, cj = probe.cellJ, ck = probe.cellK;
            const int i = ci - 1, j = cj - 1, k = ck - 1;
            double Ex_v = 0, Ey_v = 0, Ez_v = 0;
            if (ci >= 1 && ci <= NX && cj >= 1 && cj <= NY && ck >= 1 && ck <= NZ)
                Ex_v = Ex[ex_idx(i, j, k)];
            if (ci >= 1 && ci <= NX && cj >= 1 && cj <= NY && ck >= 1 && ck <= NZ)
                Ey_v = Ey[ey_idx(i, j, k)];
            if (ci >= 1 && ci <= NX && cj >= 1 && cj <= NY && ck >= 1 && ck <= NZ)
                Ez_v = Ez[ez_idx(i, j, k)];
            double inc_x = 0.0, inc_y = 0.0, inc_z = 0.0;
            for (int pwIdx = 0; pwIdx < (int)planeWaves.size(); ++pwIdx) {
                inc_x += computeIncid(pwIdx, 0, currentTime, ci, cj, ck);
                inc_y += computeIncid(pwIdx, 1, currentTime, ci, cj, ck);
                inc_z += computeIncid(pwIdx, 2, currentTime, ci, cj, ck);
            }
            probe.timeData.push_back(currentTime);
            for (size_t d = 0; d < probe.directions.size(); ++d) {
                const auto& dir = probe.directions[d];
                double val = 0.0, inc = 0.0;
                if (dir == "x") { val = Ex_v; inc = inc_x; }
                else if (dir == "y") { val = Ey_v; inc = inc_y; }
                else if (dir == "z") { val = Ez_v; inc = inc_z; }
                probe.fieldByDir[d].push_back(val);
                probe.incidentByDir[d].push_back(inc);
            }
        }
    }

    void writeProbeOutputs(const std::string& caseName) {
        for (auto& probe : probes) {
            if (probe.domainType != "time") continue;
            std::string fname = caseName + ".fdtd_" + probe.name + "_";
            std::string fieldDir = (probe.field == "magnetic") ? "H" : "E";
            for (size_t d = 0; d < probe.directions.size(); ++d) {
                std::string fullname = fname + fieldDir + probe.directions[d] + "_";
                fullname += std::to_string(probe.cellI) + "_" + std::to_string(probe.cellJ) + "_" +
                            std::to_string(probe.cellK) + ".dat";
                std::ofstream out(fullname);
                out << std::scientific << std::setprecision(17);
                out << "t              " << fullname << "       incid" << std::endl;
                for (size_t t = 0; t < probe.timeData.size(); ++t) {
                    out << probe.timeData[t];
                    out << "   " << probe.fieldByDir[d][t];
                    out << "   " << probe.incidentByDir[d][t];
                    out << std::endl;
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
    const std::string suffix = ".fdtd.json";
    if (name.size() > suffix.size() && name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0) {
        name = name.substr(0, name.size() - suffix.size());
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
    mtln_t mtln_parsed;
    taglist_t tag_numbers;
    tagtype_t tagtype;
    FDTD_Solver solver;

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
    if (filename.size() > 5) {
        impl_->l.extension = filename.substr(filename.size() - 5);
    }
    impl_->solver.init(filename, mapvtk::flagsContainMapVtk(input_flags));
    impl_->media.NumMed = impl_->solver.pd.Mats.nMaterials;
    impl_->media.totalX = impl_->solver.pd.matriz.totalX;
    impl_->media.totalY = impl_->solver.pd.matriz.totalY;
    impl_->media.totalZ = impl_->solver.pd.matriz.totalZ;
    impl_->l.layoutnumber = 0;
}

void semba_fdtd_t::launch() {
    impl_->solver.launch();
}

void semba_fdtd_t::end(const std::string& case_name) {
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
    if (got.field.empty() || ref.field.empty()) {
        std::cerr << "pw-in-box " << name << " probe missing data: got="
                  << got.field.size() << " ref=" << ref.field.size() << std::endl;
        return 1;
    }
    const size_t wanted = max_steps >= 0
                              ? static_cast<size_t>(max_steps) + 1
                              : std::min(got.field.size(), ref.field.size());
    if (got.field.size() < wanted || ref.field.size() < wanted) {
        std::cerr << "pw-in-box " << name << " probe too short: got="
                  << got.field.size() << " ref=" << ref.field.size()
                  << " wanted=" << wanted << std::endl;
        return 1;
    }

    constexpr double field_atol = 3e-4;
    constexpr double field_rtol = 1e-3;
    constexpr double time_atol = 1e-15;
    for (size_t s = 0; s < wanted; ++s) {
        const double dt = std::abs(got.time[s] - ref.time[s]);
        const double df = std::abs(got.field[s] - ref.field[s]);
        const double tol = field_atol + field_rtol * std::abs(ref.field[s]);
        if (dt > time_atol || df > tol) {
            std::cerr << "pw-in-box " << name << " mismatch at sample " << s
                      << ": t got=" << got.time[s] << " ref=" << ref.time[s]
                      << ", field got=" << got.field[s] << " ref=" << ref.field[s]
                      << ", diff=" << df << ", tol=" << tol << std::endl;
            return 1;
        }
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
    const int idx0 = solver.hy_idx(0, j, k);
    const int idx1 = solver.hy_idx(1, j, k);
    solver.Hy[idx0] = 1.0;
    solver.Hy[idx1] = 0.5;
    solver.murPastHyBack[static_cast<size_t>(j * (solver.NZ + 1) + k)] = 0.2;
    solver.murPastHyBackInt[static_cast<size_t>(j * (solver.NZ + 1) + k)] = 0.4;
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

    const ProbeSeries got_b = readProbeDat(case_name + ".fdtd_before_Ex_3_3_1.dat");
    const ProbeSeries got_i = readProbeDat(case_name + ".fdtd_inbox_Ex_3_3_3.dat");
    const ProbeSeries got_a = readProbeDat(case_name + ".fdtd_after_Ex_3_3_5.dat");
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
