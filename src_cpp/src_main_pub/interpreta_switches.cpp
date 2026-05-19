#include <vector>
#include <string>
#include <iostream>
struct SGGFDTDINFO_t { int NumberRequest=0; struct { bool Volumic=false; int nP=0; int What[10]={}; int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; bool done=false,flushed=false,Begun=false,TimeDomain=false,FreqDomain=false; } observation[100]; struct { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; } Alloc[10]; };
struct sim_control_t { int layoutnumber=0; int num_procs=1; bool resume=false; bool stochastic=false; };
struct media_matrices_t { int NumMed=0; };
struct limit_t { int XI=0,XE=0,YI=0,YE=0,ZI=0,ZE=0; };
struct entrada_t { int layout=0; };
void interpreta_switches(SGGFDTDINFO_t&, sim_control_t&, media_matrices_t&, limit_t&, int, char**) {}
void readSwitches(int, char**, sim_control_t&, SGGFDTDINFO_t&) {}
