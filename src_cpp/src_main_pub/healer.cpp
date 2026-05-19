#include <vector>
#include <string>
#include <iostream>
struct SGGFDTDINFO_t { int NumberRequest=0; struct { bool Volumic=false; int nP=0; int What[10]={}; int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; bool done=false,flushed=false,Begun=false,TimeDomain=false,FreqDomain=false; } observation[100]; struct { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; } Alloc[10]; };
struct limit_t { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; };
struct media_matrices_t { int NumMed=0; };
struct conformal_geometry_t {};
struct PECElements_t {};
void fillBoundaryFaceIfAllEdgesPEC(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) {}
void fillEdgesInsideVolumeX(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) {}
void fillEdgesInsideVolumeY(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) {}
void fillEdgesInsideVolumeZ(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) {}
void healPEC(const SGGFDTDINFO_t&, const limit_t&, media_matrices_t&, conformal_geometry_t&, PECElements_t&) {}
void healPMC(const SGGFDTDINFO_t&, const limit_t&, media_matrices_t&, conformal_geometry_t&, PECElements_t&) {}
void healDielectric(const SGGFDTDINFO_t&, const limit_t&, media_matrices_t&, conformal_geometry_t&, PECElements_t&) {}
void healThinSlot(const SGGFDTDINFO_t&, const limit_t&, media_matrices_t&, conformal_geometry_t&, PECElements_t&) {}
void healMTLNWires(const SGGFDTDINFO_t&, const limit_t&, media_matrices_t&, conformal_geometry_t&, PECElements_t&) {}
