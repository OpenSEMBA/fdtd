// Minimal observation stubs for semba-outputs until observation.cpp is fully translated.
#include <string>
#include <vector>

struct SGGFDTDINFO_t;
struct media_matrices_t;
struct taglist_t;
struct limit_t;
struct bounds_t;
struct sim_control_t;
struct nf2ff_t;

using RKIND = double;
using RKIND_tiempo = double;

void InitObservation(SGGFDTDINFO_t& sgg, media_matrices_t& media, taglist_t& tag_numbers,
                     bool& ThereAreObservation, bool& ThereAreWires, bool& ThereAreFarFields,
                     int& initialtimestep, RKIND_tiempo& lastexecutedtime,
                     const std::vector<limit_t>& SINPML_fullsize, RKIND eps00, RKIND mu00,
                     bounds_t& bounds, sim_control_t& control) {
    (void)sgg;
    (void)media;
    (void)tag_numbers;
    (void)SINPML_fullsize;
    (void)eps00;
    (void)mu00;
    (void)bounds;
    (void)control;
    ThereAreObservation = false;
    ThereAreWires = false;
    ThereAreFarFields = false;
    initialtimestep = 0;
    lastexecutedtime = 0.0;
}

void UpdateObservation(SGGFDTDINFO_t& sgg, media_matrices_t& media, taglist_t& tag_numbers,
                       int timestep, int ini_save,
                       const std::vector<std::vector<std::vector<RKIND>>>& Ex,
                       const std::vector<std::vector<std::vector<RKIND>>>& Ey,
                       const std::vector<std::vector<std::vector<RKIND>>>& Ez,
                       const std::vector<std::vector<std::vector<RKIND>>>& Hx,
                       const std::vector<std::vector<std::vector<RKIND>>>& Hy,
                       const std::vector<std::vector<std::vector<RKIND>>>& Hz,
                       const std::vector<RKIND>& dxe, const std::vector<RKIND>& dye,
                       const std::vector<RKIND>& dze, const std::vector<RKIND>& dxh,
                       const std::vector<RKIND>& dyh, const std::vector<RKIND>& dzh,
                       const std::string& wiresflavor, const std::vector<limit_t>& SINPML_fullsize,
                       bool wirecrank, bool noconformalmapvtk, bounds_t& bounds) {
    (void)sgg;
    (void)media;
    (void)tag_numbers;
    (void)timestep;
    (void)ini_save;
    (void)Ex;
    (void)Ey;
    (void)Ez;
    (void)Hx;
    (void)Hy;
    (void)Hz;
    (void)dxe;
    (void)dye;
    (void)dze;
    (void)dxh;
    (void)dyh;
    (void)dzh;
    (void)wiresflavor;
    (void)SINPML_fullsize;
    (void)wirecrank;
    (void)noconformalmapvtk;
    (void)bounds;
}

void FlushObservationFiles(SGGFDTDINFO_t& sgg, int nInit, int FinalInstant, int layoutnumber,
                           int num_procs, const std::vector<RKIND>& dxe, const std::vector<RKIND>& dye,
                           const std::vector<RKIND>& dze, const std::vector<RKIND>& dxh,
                           const std::vector<RKIND>& dyh, const std::vector<RKIND>& dzh,
                           bounds_t& b, bool singlefilewrite, const nf2ff_t& facesNF2FF, bool flushff) {
    (void)sgg;
    (void)nInit;
    (void)FinalInstant;
    (void)layoutnumber;
    (void)num_procs;
    (void)dxe;
    (void)dye;
    (void)dze;
    (void)dxh;
    (void)dyh;
    (void)dzh;
    (void)b;
    (void)singlefilewrite;
    (void)facesNF2FF;
    (void)flushff;
}

void CloseObservationFiles(SGGFDTDINFO_t& sgg, int layoutnumber, int num_procs, bool singlefilewrite,
                           int initialtimestep, RKIND_tiempo lastexecutedtime, bool resume) {
    (void)sgg;
    (void)layoutnumber;
    (void)num_procs;
    (void)singlefilewrite;
    (void)initialtimestep;
    (void)lastexecutedtime;
    (void)resume;
}

void DestroyObservation(SGGFDTDINFO_t& sgg) { (void)sgg; }

void unpacksinglefiles() {}

#ifdef CompileWithMTLN
void InitObservationMTLN(const std::string& nEntradaRoot) { (void)nEntradaRoot; }
void UpdateObservationMTLN(int step) { (void)step; }
void CloseObservationFilesMTLN() {}
#endif
