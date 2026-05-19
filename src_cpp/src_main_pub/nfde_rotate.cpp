#include <vector>
#include <string>
#include <iostream>
struct SGGFDTDINFO_t { int NumberRequest=0; struct { bool Volumic=false; int nP=0; int What[10]={}; int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; bool done=false,flushed=false,Begun=false,TimeDomain=false,FreqDomain=false; } observation[100]; struct { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; } Alloc[10]; };
struct limit_t { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; };
struct media_matrices_t { int NumMed=0; };
struct mtln_t { int numWires=0; };
void rotate_generateSpaceSteps(const SGGFDTDINFO_t&, const limit_t&, media_matrices_t&, int) {}
void rotate_generateCurrent_Field_Sources(const SGGFDTDINFO_t&, const limit_t&, media_matrices_t&, int) {}
void rotate_generatePlaneWaves(const SGGFDTDINFO_t&, const limit_t&, media_matrices_t&, int) {}
void rotate_generateBoxSources(const SGGFDTDINFO_t&, const limit_t&, media_matrices_t&, int) {}
void rotate_generateFronteras(const SGGFDTDINFO_t&, const limit_t&, media_matrices_t&, int) {}
