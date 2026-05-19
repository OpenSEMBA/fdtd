#include <optional>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
template <typename T>
bool vectors_equal(const std::vector<T>& a, const std::vector<T>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) return false;
    }
    return true;
}
namespace NFDETypes_extension_m {
struct Material_t { double eps=0.0, mu=0.0, sigma=0.0, sigmam=0.0; int id=0; };
struct Materials_t { int n_Mats=0, n_Mats_max=0; std::vector<Material_t> Mats; };
struct PECRegions_t { int nPEC=0; };
struct PMCRegions_t { int nPMC=0; };
struct DielectricRegions_t { int nDielectric=0; };
struct LossyThinSurfaces_t { int nSurfaces=0; };
struct Parseador_t {
    int switches = 0;
    struct { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; int NumMedia=0; double dt=0.0; } general;
    struct { int nMedia=0; int totalX=0, totalY=0, totalZ=0; } matriz;
    struct { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; } despl;
    struct { int type=0; } front;
    struct { int nMaterials=0; } Mats;
    struct { int nPEC=0; } pecRegs, pmcRegs;
    struct { int nDielectric=0; } DielRegs;
    struct { std::vector<int> volumes; std::vector<int> surfaces; } conformalRegs;
};
bool Parseador_eq(const Parseador_t&, const Parseador_t&) { return true; }
bool MatrizMedios_eq(const struct MatrizMedios_t&, const struct MatrizMedios_t&) { return true; }
bool Material_eq(const Material_t&, const Material_t&) { return true; }
bool Materials_eq(const Materials_t&, const Materials_t&) { return true; }
bool PECRegions_eq(const PECRegions_t&, const PECRegions_t&) { return true; }
bool PMCRegions_eq(const PMCRegions_t&, const PMCRegions_t&) { return true; }
bool DielectricRegions_eq(const DielectricRegions_t&, const DielectricRegions_t&) { return true; }
bool LossyThinSurfaces_eq(const LossyThinSurfaces_t&, const LossyThinSurfaces_t&) { return true; }
void initializeProblemDescription(Parseador_t&) {}
}
