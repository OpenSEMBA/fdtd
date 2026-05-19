#include <vector>
#include <string>
#include <iostream>
struct Parseador_t { int dummy=0; };
struct mesh_t { std::vector<std::vector<double>> coordinates; std::vector<std::vector<int>> elements; };
struct SGGFDTDINFO_t { int NumberRequest=0; struct { bool Volumic=false; int nP=0; int What[10]={}; int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; bool done=false,flushed=false,Begun=false,TimeDomain=false,FreqDomain=false; } observation[100]; struct { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; } Alloc[10]; };
struct limit_t { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; };
struct conformal_geometry_t {};
struct PECElements_t {};
void addCoordinates(mesh_t&, const Parseador_t&) {}
void addElements(mesh_t&, const Parseador_t&) {}
void buildMesh(mesh_t&, const Parseador_t&, const SGGFDTDINFO_t&, const limit_t&, conformal_geometry_t&, PECElements_t&) {}
void buildGeometry(mesh_t&, const Parseador_t&, const SGGFDTDINFO_t&, const limit_t&, conformal_geometry_t&, PECElements_t&) {}
