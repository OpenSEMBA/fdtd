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
#include <functional>

const double PI = 3.14159265358979323846;
const double EPS0 = 8.854187817e-12;
const double MU0 = 1.2566370614e-6;
const double C0 = 299792458.0;
const double ZVAC = std::sqrt(MU0/EPS0);
const int BUFSIZE = 1024;

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
    std::vector<double> timeData, valueData;
};
struct boundary_t { std::string type = "mur"; int layers = 10; int order = 2; double reflection = 0.001; };

struct PlaneWaveState_t {
    std::vector<double> px, py, pz;
    std::vector<double> ex, ey, ez;
    std::vector<double> hx, hy, hz;
    std::vector<double> samples;
    double deltaevol = 0.0;
    int numSamples = 0;
    double distanciaInicial = 0.0;
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

struct ExcitationData { std::vector<double> times, values; };
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
    // Always compute CFL-limited dt from mesh steps if steps are provided
    if (root.contains("mesh") && root["mesh"].contains("grid") && root["mesh"]["grid"].contains("steps")) {
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
    std::vector<double> Ex, Ey, Ez, Hx, Hy, Hz;
    std::vector<int> pecMask;
    std::vector<double> CeE, CmH;
    std::vector<source_t> sources;
    std::map<std::string, ExcitationData> excitations;
    std::vector<PlaneWaveState_t> planeWaves;
    std::vector<probe_output_t> probes;
    bool still_planewave_time = false;
    int numSteps = 100, step = 0;
    double currentTime = 0.0;

    void init(const std::string& filename) {
        pd = parseFDTDJSON(filename);
        NX = pd.general.XI; NY = pd.general.YI; NZ = pd.general.ZI;
        if (NX <= 0) NX = 10; if (NY <= 0) NY = 10; if (NZ <= 0) NZ = 10;
        dt = pd.general.dt; if (dt <= 0.0) dt = 1e-12;
        if (!pd.cellSteps.cellStepsX.empty()) dx = pd.cellSteps.cellStepsX[0]; else dx = 0.025;
        if (!pd.cellSteps.cellStepsY.empty()) dy = pd.cellSteps.cellStepsY[0]; else dy = 0.025;
        if (!pd.cellSteps.cellStepsZ.empty()) dz = pd.cellSteps.cellStepsZ[0]; else dz = 0.025;
        numSteps = pd.general.XE; if (numSteps <= 0) numSteps = 100;

        int ex_n = (NX+1)*NY*NZ, ey_n = NX*(NY+1)*NZ, ez_n = NX*NY*(NZ+1);
        int hx_n = NX*(NY+1)*(NZ+1), hy_n = (NX+1)*NY*(NZ+1), hz_n = (NX+1)*(NY+1)*NZ;
        Ex.resize(ex_n,0); Ey.resize(ey_n,0); Ez.resize(ez_n,0);
        Hx.resize(hx_n,0); Hy.resize(hy_n,0); Hz.resize(hz_n,0);
        int max_n = std::max({ex_n,ey_n,ez_n,hx_n,hy_n,hz_n});
        pecMask.resize(max_n, 0);
        CeE.resize(max_n, 0.0); CmH.resize(max_n, 0.0);

        for (int i = 0; i < max_n; i++) { CeE[i] = dt/eps0; CmH[i] = dt/mu0; }

        sources = pd.sources.planeWaves;
        for (auto& src : sources) {
            if (!src.magnitudeFile.empty()) excitations[src.magnitudeFile] = readExcitationFile(src.magnitudeFile);
        }
        planeWaves.resize(sources.size());
        for (int i = 0; i < (int)sources.size(); i++) {
            planeWaves[i].px.resize(1,0); planeWaves[i].py.resize(1,0); planeWaves[i].pz.resize(1,0);
            planeWaves[i].ex.resize(1,0); planeWaves[i].ey.resize(1,0); planeWaves[i].ez.resize(1,0);
            planeWaves[i].hx.resize(1,0); planeWaves[i].hy.resize(1,0); planeWaves[i].hz.resize(1,0);
            initPlaneWave(i);
        }
        probes = pd.probes.probes;

        std::cout << "FDTD: grid=" << NX << "x" << NY << "x" << NZ << " dt=" << dt << " steps=" << numSteps << std::endl;
    }

    int ex_idx(int i,int j,int k) { return i*NY*NZ + j*NZ + k; }
    int ey_idx(int i,int j,int k) { return i*(NY+1)*NZ + j*NZ + k; }
    int ez_idx(int i,int j,int k) { return i*NY*(NZ+1) + j*(NZ+1) + k; }
    int hx_idx(int i,int j,int k) { return i*(NY+1)*(NZ+1) + j*(NZ+1) + k; }
    int hy_idx(int i,int j,int k) { return i*NY*(NZ+1) + j*(NZ+1) + k; }
    int hz_idx(int i,int j,int k) { return i*(NY+1)*NZ + j*NZ + k; }

    double evolucion(int pwIdx, double t_delay) {
        auto& pw = planeWaves[pwIdx];
        if (pw.numSamples <= 1) return 0.0;
        if (pw.deltaevol <= 0.0) return 0.0;
        double t_min = 0.0, t_max = (pw.numSamples - 1) * pw.deltaevol;
        if (t_delay < t_min || t_delay > t_max) return 0.0;
        int nprev = (int)std::floor(t_delay / pw.deltaevol);
        if (nprev < 0) nprev = 0;
        if (nprev + 1 >= pw.numSamples) {
            if (nprev < pw.numSamples) { still_planewave_time = true; return pw.samples[nprev]; }
            return 0.0;
        }
        still_planewave_time = true;
        double t_frac = t_delay - nprev * pw.deltaevol;
        double val = pw.samples[nprev] + (pw.samples[nprev+1] - pw.samples[nprev]) * (t_frac / pw.deltaevol);
        if (std::abs(val) > 1e6) return 0.0;
        return val;
    }

    double computeIncid(int pwIdx, int nfield, double time, int i, int j, int k) {
        auto& pw = planeWaves[pwIdx];
        double xf=0, yf=0, zf=0;
        switch(nfield) {
            case 0: xf=(i>=0&&i<NX)?(i*dx+dx/2):0; yf=(j>=0&&j<NY)?(j*dy):0; zf=(k>=0&&k<NZ)?(k*dz):0; break;
            case 1: xf=(i>=0&&i<NX)?(i*dx):0; yf=(j>=0&&j<=NY)?(j*dy+dy/2):0; zf=(k>=0&&k<NZ)?(k*dz):0; break;
            case 2: xf=(i>=0&&i<NX)?(i*dx):0; yf=(j>=0&&j<NY)?(j*dy):0; zf=(k>=0&&k<=NZ)?(k*dz+dz/2):0; break;
            case 3: xf=(i>=0&&i<NX)?(i*dx):0; yf=(j>=0&&j<=NY)?(j*dy+dy/2):0; zf=(k>=0&&k<=NZ)?(k*dz+dz/2):0; break;
            case 4: xf=(i>=0&&i<NX)?(i*dx+dx/2):0; yf=(j>=0&&j<NY)?(j*dy):0; zf=(k>=0&&k<=NZ)?(k*dz+dz/2):0; break;
            case 5: xf=(i>=0&&i<NX)?(i*dx+dx/2):0; yf=(j>=0&&j<=NY)?(j*dy+dy/2):0; zf=(k>=0&&k<NZ)?(k*dz):0; break;
            default: return 0.0;
        }
        double d = (xf*pw.px[0] + yf*pw.py[0] + zf*pw.pz[0]) - pw.distanciaInicial;
        double t_arrival = time - d / C0;
        double value = evolucion(pwIdx, t_arrival);
        switch(nfield) {
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
        pw.px[0] = std::sin(src.direction.theta * PI / 180.0) * std::cos(src.direction.phi * PI / 180.0);
        pw.py[0] = std::sin(src.direction.theta * PI / 180.0) * std::sin(src.direction.phi * PI / 180.0);
        pw.pz[0] = std::cos(src.direction.theta * PI / 180.0);
        pw.ex[0] = std::sin(src.polarization.theta * PI / 180.0) * std::cos(src.polarization.phi * PI / 180.0);
        pw.ey[0] = std::sin(src.polarization.theta * PI / 180.0) * std::sin(src.polarization.phi * PI / 180.0);
        pw.ez[0] = std::cos(src.polarization.theta * PI / 180.0);
        pw.hx[0] = (pw.py[0]*pw.ez[0] - pw.pz[0]*pw.ey[0]) / ZVAC;
        pw.hy[0] = (pw.pz[0]*pw.ex[0] - pw.px[0]*pw.ez[0]) / ZVAC;
        pw.hz[0] = (pw.px[0]*pw.ey[0] - pw.py[0]*pw.ex[0]) / ZVAC;
        if (!src.magnitudeFile.empty() && excitations.count(src.magnitudeFile)) {
            auto& exc = excitations[src.magnitudeFile];
            pw.samples = exc.values;
            pw.numSamples = exc.times.size();
            if (pw.numSamples >= 2) pw.deltaevol = exc.times[1] - exc.times[0];
        }
        pw.esqx1 = 0; pw.esqx2 = NX;
        pw.esqy1 = 0; pw.esqy2 = NY;
        pw.esqz1 = 0; pw.esqz2 = NZ;
        double max_d = -1e30;
        double coords[8][3] = {
            {0,0,0}, {NX*dx,0,0}, {0,NY*dy,0}, {NX*dx,NY*dy,0},
            {0,0,NZ*dz}, {NX*dx,0,NZ*dz}, {0,NY*dy,NZ*dz}, {NX*dx,NY*dy,NZ*dz}
        };
        for (int c = 0; c < 8; c++) {
            double d = coords[c][0]*pw.px[0] + coords[c][1]*pw.py[0] + coords[c][2]*pw.pz[0];
            if (d > max_d) max_d = d;
        }
        pw.distanciaInicial = max_d;
        pw.iluminaTr = (pw.px[0] > 0);
        pw.iluminaFr = (pw.px[0] < 0);
        pw.iluminaIz = (pw.py[0] > 0);
        pw.iluminaDe = (pw.py[0] < 0);
        pw.iluminaAb = (pw.pz[0] > 0);
        pw.iluminaAr = (pw.pz[0] < 0);
        std::cout << "  PlaneWave: dir=(" << pw.px[0] << "," << pw.py[0] << "," << pw.pz[0]
                  << ") pol=(" << pw.ex[0] << "," << pw.ey[0] << "," << pw.ez[0]
                  << ") samples=" << pw.numSamples << std::endl;
    }

    void advanceE() {
        int ex_n=(NX+1)*NY*NZ, ey_n=NX*(NY+1)*NZ, ez_n=NX*NY*(NZ+1);
        for (int i=0; i<=NX; i++)
            for (int j=0; j<NY; j++)
                for (int k=0; k<NZ; k++) {
                    int idx = ex_idx(i,j,k);
                    double dHz_dy = (j+1<NY) ? (Hz[hz_idx(i,j+1,k)] - Hz[hz_idx(i,j,k)]) : (-Hz[hz_idx(i,j,k)]);
                    double dHy_dz = (k+1<NZ) ? (Hy[hy_idx(i,j,k+1)] - Hy[hy_idx(i,j,k)]) : (-Hy[hy_idx(i,j,k)]);
                    Ex[idx] += CeE[idx] * (dHz_dy/dz - dHy_dz/dy);
                }
        for (int i=0; i<NX; i++)
            for (int j=0; j<=NY; j++)
                for (int k=0; k<NZ; k++) {
                    int idx = ey_idx(i,j,k);
                    double dHz_dx = (i+1<NX) ? (Hz[hz_idx(i+1,j,k)] - Hz[hz_idx(i,j,k)]) : (-Hz[hz_idx(i,j,k)]);
                    double dHx_dz = (k+1<NZ) ? (Hx[hx_idx(i,j,k+1)] - Hx[hx_idx(i,j,k)]) : (-Hx[hx_idx(i,j,k)]);
                    Ey[idx] += CeE[ex_n+idx] * (dHz_dx/dx - dHx_dz/dz);
                }
        for (int i=0; i<NX; i++)
            for (int j=0; j<NY; j++)
                for (int k=0; k<=NZ; k++) {
                    int idx = ez_idx(i,j,k);
                    double dHx_dy = (j+1<NY) ? (Hx[hx_idx(i,j+1,k)] - Hx[hx_idx(i,j,k)]) : (-Hx[hx_idx(i,j,k)]);
                    double dHy_dx = (i+1<NX) ? (Hy[hy_idx(i+1,j,k)] - Hy[hy_idx(i,j,k)]) : (-Hy[hy_idx(i,j,k)]);
                    Ez[idx] += CeE[ex_n+ey_n+idx] * (dHx_dy/dx - dHy_dx/dy);
                }
    }

    void advanceH() {
        int hx_n=NX*(NY+1)*(NZ+1), hy_n=(NX+1)*NY*(NZ+1), hz_n=(NX+1)*(NY+1)*NZ;
        for (int i=0; i<NX; i++)
            for (int j=0; j<=NY; j++)
                for (int k=0; k<=NZ; k++) {
                    int idx = hx_idx(i,j,k);
                    double dEz_dy = (j+1<=NY) ? (Ez[ez_idx(i,j+1,k)] - Ez[ez_idx(i,j,k)]) : (-Ez[ez_idx(i,j,k)]);
                    double dEy_dz = (k+1<=NZ) ? (Ey[ey_idx(i,j,k+1)] - Ey[ey_idx(i,j,k)]) : (-Ey[ey_idx(i,j,k)]);
                    Hx[idx] += CmH[idx] * (dEz_dy/dy - dEy_dz/dz);
                }
        for (int i=0; i<=NX; i++)
            for (int j=0; j<NY; j++)
                for (int k=0; k<=NZ; k++) {
                    int idx = hy_idx(i,j,k);
                    double dEx_dz = (k+1<=NZ) ? (Ex[ex_idx(i,j,k+1)] - Ex[ex_idx(i,j,k)]) : (-Ex[ex_idx(i,j,k)]);
                    double dEz_dx = (i+1<=NX) ? (Ez[ez_idx(i+1,j,k)] - Ez[ez_idx(i,j,k)]) : (-Ez[ez_idx(i,j,k)]);
                    Hy[idx] += CmH[hx_n+idx] * (dEx_dz/dz - dEz_dx/dx);
                }
        for (int i=0; i<=NX; i++)
            for (int j=0; j<=NY; j++)
                for (int k=0; k<NZ; k++) {
                    int idx = hz_idx(i,j,k);
                    double dEy_dx = (i+1<=NX) ? (Ey[ey_idx(i+1,j,k)] - Ey[ey_idx(i,j,k)]) : (-Ey[ey_idx(i,j,k)]);
                    double dEx_dy = (j+1<=NY) ? (Ex[ex_idx(i,j+1,k)] - Ex[ex_idx(i,j,k)]) : (-Ex[ex_idx(i,j,k)]);
                    Hz[idx] += CmH[hx_n+hy_n+idx] * (dEy_dx/dx - dEx_dy/dy);
                }
    }

    void advancePlaneWaveE() {
        for (int pwIdx = 0; pwIdx < (int)planeWaves.size(); pwIdx++) {
            auto& pw = planeWaves[pwIdx];
            double G2_1 = dt / eps0;
            if (pw.iluminaTr) {
                int i = 0; double Id = dz;
                for (int k=0;k<NZ;k++) for (int j=0;j<NY;j++)
                    Ez[ez_idx(i,j,k)] -= G2_1 * computeIncid(pwIdx,4,currentTime,i-1,j,k) * Id;
                for (int k=0;k<NZ;k++) for (int j=0;j<NY;j++)
                    Ey[ey_idx(i,j,k)] += G2_1 * computeIncid(pwIdx,5,currentTime,i-1,j,k) * Id;
            }
            if (pw.iluminaFr) {
                int i = NX; double Id = dz;
                for (int k=0;k<NZ;k++) for (int j=0;j<NY;j++)
                    Ez[ez_idx(i,j,k)] += G2_1 * computeIncid(pwIdx,4,currentTime,i,j,k) * Id;
                for (int k=0;k<NZ;k++) for (int j=0;j<NY;j++)
                    Ey[ey_idx(i,j,k)] -= G2_1 * computeIncid(pwIdx,5,currentTime,i,j,k) * Id;
            }
            if (pw.iluminaIz) {
                int j = 0; double Id = dz;
                for (int k=0;k<NZ;k++) for (int i=0;i<=NX;i++)
                    Ex[ex_idx(i,j,k)] -= G2_1 * computeIncid(pwIdx,5,currentTime,i,j-1,k) * Id;
                for (int k=0;k<NZ;k++) for (int i=0;i<NX;i++)
                    Ez[ez_idx(i,j,k)] += G2_1 * computeIncid(pwIdx,3,currentTime,i,j-1,k) * Id;
            }
            if (pw.iluminaDe) {
                int j = NY; double Id = dz;
                for (int k=0;k<NZ;k++) for (int i=0;i<NX;i++)
                    Ez[ez_idx(i,j,k)] -= G2_1 * computeIncid(pwIdx,3,currentTime,i,j,k) * Id;
                for (int k=0;k<NZ;k++) for (int i=0;i<=NX;i++)
                    Ex[ex_idx(i,j,k)] += G2_1 * computeIncid(pwIdx,5,currentTime,i,j,k) * Id;
            }
            if (pw.iluminaAb) {
                int k = 0; double Id = dy;
                for (int j=0;j<NY;j++) for (int i=0;i<=NX;i++)
                    Ex[ex_idx(i,j,k)] += G2_1 * computeIncid(pwIdx,4,currentTime,i,j,k-1) * Id;
                for (int j=0;j<=NY;j++) for (int i=0;i<NX;i++)
                    Ey[ey_idx(i,j,k)] -= G2_1 * computeIncid(pwIdx,3,currentTime,i,j,k-1) * Id;
            }
            if (pw.iluminaAr) {
                int k = NZ; double Id = dy;
                for (int j=0;j<NY;j++) for (int i=0;i<=NX;i++)
                    Ex[ex_idx(i,j,k)] -= G2_1 * computeIncid(pwIdx,4,currentTime,i,j,k) * Id;
                for (int j=0;j<=NY;j++) for (int i=0;i<NX;i++)
                    Ey[ey_idx(i,j,k)] += G2_1 * computeIncid(pwIdx,3,currentTime,i,j,k) * Id;
            }
        }
    }

    void advancePlaneWaveH() {
        for (int pwIdx = 0; pwIdx < (int)planeWaves.size(); pwIdx++) {
            auto& pw = planeWaves[pwIdx];
            double Gm2_1 = dt / mu0;
            double timeH = currentTime + 0.5 * dt;
            if (pw.iluminaTr) {
                int i = 0; double Id = dy;
                for (int k=0;k<=NZ;k++) for (int j=0;j<=NY;j++)
                    Hz[hz_idx(i,j,k)] += Gm2_1 * computeIncid(pwIdx,1,timeH,i+1,j,k) * Id;
                for (int k=0;k<=NZ;k++) for (int j=0;j<NY;j++)
                    Hy[hy_idx(i,j,k)] -= Gm2_1 * computeIncid(pwIdx,2,timeH,i+1,j,k) * Id;
            }
            if (pw.iluminaFr) {
                int i = NX; double Id = dy;
                for (int k=0;k<=NZ;k++) for (int j=0;j<=NY;j++)
                    Hz[hz_idx(i,j,k)] -= Gm2_1 * computeIncid(pwIdx,1,timeH,i,j,k) * Id;
                for (int k=0;k<=NZ;k++) for (int j=0;j<NY;j++)
                    Hy[hy_idx(i,j,k)] += Gm2_1 * computeIncid(pwIdx,2,timeH,i,j,k) * Id;
            }
            if (pw.iluminaIz) {
                int j = 0; double Id = dx;
                for (int k=0;k<=NZ;k++) for (int i=0;i<NX;i++)
                    Hx[hx_idx(i,j,k)] += Gm2_1 * computeIncid(pwIdx,2,timeH,i,j+1,k) * Id;
                for (int k=0;k<=NZ;k++) for (int i=0;i<=NX;i++)
                    Hz[hz_idx(i,j,k)] -= Gm2_1 * computeIncid(pwIdx,0,timeH,i,j+1,k) * Id;
            }
            if (pw.iluminaDe) {
                int j = NY; double Id = dx;
                for (int k=0;k<=NZ;k++) for (int i=0;i<NX;i++)
                    Hx[hx_idx(i,j,k)] -= Gm2_1 * computeIncid(pwIdx,2,timeH,i,j,k) * Id;
                for (int k=0;k<=NZ;k++) for (int i=0;i<=NX;i++)
                    Hz[hz_idx(i,j,k)] += Gm2_1 * computeIncid(pwIdx,0,timeH,i,j,k) * Id;
            }
            if (pw.iluminaAb) {
                int k = 0; double Id = dx;
                for (int j=0;j<=NY;j++) for (int i=0;i<NX;i++)
                    Hx[hx_idx(i,j,k)] -= Gm2_1 * computeIncid(pwIdx,1,timeH,i,j,k+1) * Id;
                for (int j=0;j<NY;j++) for (int i=0;i<=NX;i++)
                    Hy[hy_idx(i,j,k)] += Gm2_1 * computeIncid(pwIdx,0,timeH,i,j,k+1) * Id;
            }
            if (pw.iluminaAr) {
                int k = NZ; double Id = dx;
                for (int j=0;j<=NY;j++) for (int i=0;i<NX;i++)
                    Hx[hx_idx(i,j,k)] += Gm2_1 * computeIncid(pwIdx,1,timeH,i,j,k) * Id;
                for (int j=0;j<NY;j++) for (int i=0;i<=NX;i++)
                    Hy[hy_idx(i,j,k)] -= Gm2_1 * computeIncid(pwIdx,0,timeH,i,j,k) * Id;
            }
        }
    }

    void applyPlaneWaveSource() {}

    void sampleProbes() {
        for (auto& probe : probes) {
            if (probe.domainType != "time" || probe.elementIds.empty()) continue;
            int ci = probe.elementIds[0], cj = (probe.elementIds.size()>1)?probe.elementIds[1]:1, ck = (probe.elementIds.size()>2)?probe.elementIds[2]:1;
            double Ex_v=0,Ey_v=0,Ez_v=0;
            if (ci<=NX && cj<NY && ck<NZ) Ex_v = Ex[ex_idx(ci,cj,ck)];
            if (ci<NX && cj<=NY && ck<NZ) Ey_v = Ey[ey_idx(ci,cj,ck)];
            if (ci<NX && cj<NY && ck<=NZ) Ez_v = Ez[ez_idx(ci,cj,ck)];
            probe.timeData.push_back(currentTime);
            for (auto& dir : probe.directions) {
                if (dir=="x") probe.valueData.push_back(Ex_v);
                else if (dir=="y") probe.valueData.push_back(Ey_v);
                else if (dir=="z") probe.valueData.push_back(Ez_v);
                else probe.valueData.push_back(0.0);
            }
        }
    }

    void writeProbeOutputs(const std::string& caseName) {
        for (auto& probe : probes) {
            if (probe.domainType != "time") continue;
            std::string fname = caseName + ".fdtd_" + probe.name + "_";
            std::string fieldDir = (probe.field=="magnetic") ? "H" : "E";
            for (size_t d=0; d<probe.directions.size(); d++) {
                std::string fullname = fname + fieldDir + probe.directions[d] + "_";
                int ci=1,cj=1,ck=1;
                if (!probe.elementIds.empty()) {
                    ci=probe.elementIds[0];
                    if (probe.elementIds.size()>1) cj=probe.elementIds[1];
                    if (probe.elementIds.size()>2) ck=probe.elementIds[2];
                }
                fullname += std::to_string(ci)+"_"+std::to_string(cj)+"_"+std::to_string(ck)+".dat";
                std::ofstream out(fullname);
                out << std::scientific << std::setprecision(17);
                out << "t              " << fullname << "       incid" << std::endl;
                for (size_t t=0; t<probe.timeData.size(); t++) {
                    out << probe.timeData[t];
                    out << "   " << probe.valueData[t];
                    out << "   0.000000000E+000";
                    out << std::endl;
                }
                out.close();
            }
        }
    }

    void launch() {
        std::cout << "Running FDTD: " << numSteps << " steps..." << std::endl;
        for (step = 0; step <= numSteps; step++) {
            currentTime = step * dt;
            if (step == 0) { sampleProbes(); }
            else { advanceE(); advancePlaneWaveE(); advanceH(); advancePlaneWaveH(); sampleProbes(); }
            if (step % 500 == 0 || step == numSteps)
                std::cout << "  Step " << step << "/" << numSteps << " (t=" << currentTime << "s)" << std::endl;
        }
        std::cout << "Simulation complete." << std::endl;
    }

    void end(const std::string& caseName) {
        writeProbeOutputs(caseName);
        std::cout << "Output files written." << std::endl;
    }
};

namespace SEMBA_FDTD_m {
struct semba_fdtd_t {
    entrada_t l; tiempo_t time_comienzo; double time_desdelanzamiento=0.0;
    media_matrices_t media; SGGFDTDINFO_t sgg; limit_t fullsize[6], SINPML_fullsize[6];
    double eps0, mu0, cluz; double maxSourceValue=0.0;
    char whoami[BUFSIZE], whoamishort[BUFSIZE];
    mtln_t mtln_parsed; taglist_t tag_numbers; tagtype_t tagtype;
    bool finishedwithsuccess = false;
    FDTD_Solver solver;
    semba_fdtd_t() : eps0(EPS0), mu0(MU0), cluz(C0) {
        strcpy(whoami,"semba-fdtd-cpp"); strcpy(whoamishort,"semba-fdtd");
    }
    void init(const std::string& input_flags="") {
        l.input_flags = input_flags;
        std::string filename = input_flags;
        if (filename.size()>5) l.extension = filename.substr(filename.size()-5);
        solver.init(filename);
        media.NumMed = solver.pd.Mats.nMaterials;
        media.totalX = solver.pd.matriz.totalX;
        media.totalY = solver.pd.matriz.totalY;
        media.totalZ = solver.pd.matriz.totalZ;
        l.layoutnumber = 1;
    }
    void launch() { solver.launch(); }
    void end(const std::string& cn) { solver.end(cn); finishedwithsuccess = true; }
};
}

extern "C" {
    SEMBA_FDTD_m::semba_fdtd_t* create_semba_fdtd() { return new SEMBA_FDTD_m::semba_fdtd_t(); }
    void destroy_semba_fdtd(SEMBA_FDTD_m::semba_fdtd_t* p) { delete p; }
    void semba_fdtd_init(SEMBA_FDTD_m::semba_fdtd_t* p, const char* flags) { if(p) p->init(flags?flags:""); }
    void semba_fdtd_launch(SEMBA_FDTD_m::semba_fdtd_t* p) { if(p) p->launch(); }
    void semba_fdtd_end(SEMBA_FDTD_m::semba_fdtd_t* p) { if(p) p->end("semba-fdtd"); }
}
