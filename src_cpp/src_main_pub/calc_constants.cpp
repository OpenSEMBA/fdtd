#include <vector>
#include <cmath>
struct SGGFDTDINFO_t { int NumberRequest=0; struct { bool Volumic=false; int nP=0; int What[10]={}; int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; bool done=false,flushed=false,Begun=false,TimeDomain=false,FreqDomain=false; } observation[100]; struct { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; } Alloc[10]; int NumMedia=0; double dt=0.0; };
struct media_matrices_t { int NumMed=0; };
double getEps0(const SGGFDTDINFO_t&) { return 8.854187817e-12; }
double getMu0(const SGGFDTDINFO_t&) { return 1.25663706212e-6; }
double getC(const SGGFDTDINFO_t&) { return 299792458.0; }
double getCellSize(const SGGFDTDINFO_t&, int) { return 1.0; }
double getDx(const SGGFDTDINFO_t&) { return 1.0; }
double getDy(const SGGFDTDINFO_t&) { return 1.0; }
double getDz(const SGGFDTDINFO_t&) { return 1.0; }
