#include <string>
#include <nlohmann/json.hpp>
#include <vector>
#include <memory>
#include <iostream>
#include <algorithm>
#include <cstring>
struct mesh_t { std::vector<std::vector<double>> coordinates; std::vector<std::vector<int>> elements; };
struct Mesh_t { std::vector<std::vector<double>> coordinates; std::vector<std::vector<int>> elements; };
struct IdChildTable_t { IdChildTable_t() {} IdChildTable_t(const nlohmann::json&, const nlohmann::json&, const std::string&) {} bool hasKey(const std::string&) const { return false; } int size() const { return 0; } nlohmann::json operator[](const std::string&) const { return nlohmann::json(); } };
struct NFDEGeneral_t { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; int NumMedia=0; double dt=0.0; };
struct Desplazamiento_t { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; };
struct MTLN_t { int numWires=0; };
struct SGGFDTDINFO_t { int NumberRequest=0; struct { bool Volumic=false; int nP=0; int What[10]={}; int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; bool done=false,flushed=false,Begun=false,TimeDomain=false,FreqDomain=false; } observation[100]; struct { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; } Alloc[10]; };
struct limit_t { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; };
struct sim_control_t { int layoutnumber=0; int num_procs=1; bool resume=false; bool stochastic=false; };
struct fhash_tbl_t { void insert(int, const std::vector<int>&) {} void insert(int, const std::vector<std::vector<int>>&) {} void insert(int, const std::vector<double>&) {} void insert(int, const std::vector<std::vector<double>>&) {} };
struct MatrizMedios_t { int nMedia=0; };
struct Materials_t { int nMaterials=0; };
struct PECRegions_t { int nPEC=0; };
struct DielectricRegions_t { int nDielectric=0; };
struct LossyThinSurfaces_t { int nSurfaces=0; };
struct Frontera_t { int type=0; };
struct Planewave_t { int nPW=0; };
struct NodalSource_t { int nSrc=0; };
namespace smbjson_m {
struct parser_t {
    nlohmann::json root;
    parser_t() : root(nlohmann::json::object()) {}
    NFDEGeneral_t readGeneral() { return NFDEGeneral_t{}; }
    Mesh_t readMesh() { return Mesh_t{}; }
    void readMTLN(MTLN_t&) {}
    Desplazamiento_t readGrid() { return Desplazamiento_t{}; }
    std::vector<std::string> readMedia() { return {}; }
    std::vector<std::string> readMaterials() { return {}; }
    PECRegions_t readPECRegions() { return PECRegions_t{}; }
    DielectricRegions_t readDielectricRegions() { return DielectricRegions_t{}; }
    LossyThinSurfaces_t readLossyThinSurfaces() { return LossyThinSurfaces_t{}; }
    std::vector<Frontera_t> readBoundaries() { return {}; }
    std::vector<Planewave_t> readPlaneWaves() { return {}; }
    std::vector<NodalSource_t> readNodalSources() { return {}; }
    bool parse(const std::string&) { return true; }
};
}
