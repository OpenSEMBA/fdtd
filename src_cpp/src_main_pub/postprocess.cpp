#include <vector>
#include <string>
#include <iostream>
#include <cstring>
namespace OBSERVATION_m {
struct Observation_t { struct { int P[100]; } P; bool TimeDomain=false; bool FreqDomain=false; };
}
struct SGGFDTDINFO_t { int NumberRequest=0; struct { bool Volumic=false; int nP=0; int What[10]={}; int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; bool done=false,flushed=false,Begun=false,TimeDomain=false,FreqDomain=false; } observation[100]; struct { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; } Alloc[10]; };
const int BUFSIZE = 1024;
void chninstant(const std::string&, const std::vector<double>&, const double*, int, int, int, int, int, int, int, int, int, int) {}
void postprocess(const SGGFDTDINFO_t&, const std::string&, int, int) {}
void writeProbeData(const std::string&, const std::vector<double>&, int) {}
