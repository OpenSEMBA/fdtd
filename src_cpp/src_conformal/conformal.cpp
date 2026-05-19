#include <vector>
#include <memory>
#include <typeinfo>
#include <algorithm>
#include <iostream>
namespace conformal_m {
struct conformal_geometry_t { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; };
struct PECElements_t { int nPEC=0; };
void buildCellMap(void*, const void*) {}
void buildConformalGeometry(conformal_geometry_t*, const void*, const void*, const void*, const void*, const void*, bool) {}
void buildConformalPECElements(PECElements_t*, const void*, const void*, bool) {}
void buildConformalPECElements(PECElements_t*, const void*, const void*, const void*, bool) {}
}
