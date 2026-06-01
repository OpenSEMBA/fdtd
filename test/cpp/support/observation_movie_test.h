#ifndef OBSERVATION_MOVIE_TEST_H
#define OBSERVATION_MOVIE_TEST_H

#include <cstdint>
#include <string>
#include <vector>

namespace Observa_m {

using RKIND = double;
using RKIND_tiempo = double;

struct XYZlimit_t {
    int32_t XI = 0;
    int32_t XE = 0;
    int32_t YI = 0;
    int32_t YE = 0;
    int32_t ZI = 0;
    int32_t ZE = 0;
};

struct limit_t : XYZlimit_t {
    int32_t NX = 0;
    int32_t NY = 0;
    int32_t NZ = 0;
};

struct observable_t {
    int32_t XI = 0;
    int32_t YI = 0;
    int32_t ZI = 0;
    int32_t XE = 0;
    int32_t YE = 0;
    int32_t ZE = 0;
    int32_t Xtrancos = 1;
    int32_t Ytrancos = 1;
    int32_t Ztrancos = 1;
    int32_t What = 0;
};

struct Obses_movie_t {
    RKIND_tiempo TimeStep = 0.0;
    RKIND_tiempo InitialTime = 0.0;
    RKIND_tiempo FinalTime = 0.0;
    int32_t nP = 0;
    bool Volumic = false;
    bool TimeDomain = false;
    bool FreqDomain = false;
    bool Saveall = false;
    bool TransFer = false;
    bool Done = false;
    bool Begun = false;
    bool Flushed = false;
    RKIND InitialFreq = 0.0;
    RKIND FinalFreq = 0.0;
    RKIND FreqStep = 0.0;
    std::string outputrequest;
    std::vector<observable_t> P;
};

struct SGGFDTDINFO_t {
    int32_t NumMedia = 0;
    int32_t NumberRequest = 0;
    RKIND_tiempo dt = 0.0;
    std::vector<RKIND_tiempo> tiempo;
    std::vector<XYZlimit_t> Sweep;
    std::vector<XYZlimit_t> SINPMLSweep;
    std::vector<XYZlimit_t> alloc;
    std::vector<Obses_movie_t> Observation;
};

struct media_matrices_t {};
struct taglist_t {};
struct bounds_t {};
struct nf2ff_t {
    bool tr = false;
    bool fr = false;
    bool iz = false;
    bool de = false;
    bool ab = false;
    bool ar = false;
};

struct sim_control_t {
    int32_t layoutnumber = 0;
    int32_t num_procs = 0;
    int32_t mpidir = 3;
    int32_t finalTimeStep = 0;
    std::string nEntradaRoot;
    std::string wiresflavor;
    bool resume = false;
    bool saveall = false;
    bool NF2FFDecim = false;
    bool simu_devia = false;
    bool singleFileWrite = false;
    nf2ff_t facesNF2FF;
};

struct output_item_t {
    int unit = 0;
    std::string path;
};

struct output_movie_t {
    int timeswritten = 0;
    std::vector<output_item_t> item;
};

std::string get_temp_dir();
XYZlimit_t create_xyz_limit(int xi, int yi, int zi, int xe, int ye, int ze);
limit_t create_limit(int xi, int xe, int yi, int ye, int zi, int ze, int nx, int ny, int nz);
SGGFDTDINFO_t create_base_sgg();
Obses_movie_t define_time_movie_observation();
media_matrices_t create_media(const std::vector<XYZlimit_t>& alloc);
taglist_t create_tag_list(const std::vector<XYZlimit_t>& alloc);
nf2ff_t create_faces_nf2ff(bool tr, bool fr, bool iz, bool de, bool ab, bool ar);
sim_control_t create_control_flags(int layoutnumber, int num_procs, int mpidir, int finaltimestep,
                                   const std::string& nEntradaRoot, const std::string& wiresflavor,
                                   bool resume, bool saveall, bool nf2ffDecim, bool simu_devia,
                                   bool singlefilewrite, const nf2ff_t& facesNF2FF);

void InitObservation(SGGFDTDINFO_t& sgg, media_matrices_t& media, taglist_t& tag_numbers,
                     bool& ThereAreObservation, bool& ThereAreWires, bool& ThereAreFarFields,
                     int& initialtimestep, RKIND_tiempo& lastexecutedtime,
                     const std::vector<limit_t>& SINPML_fullsize, RKIND eps, RKIND mu,
                     bounds_t& bounds, sim_control_t& control);

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
                       bool wirecrank, bool noconformalmapvtk, bounds_t& bounds);

std::vector<output_movie_t>& GetOutput();

struct dummyFields_t {
    std::vector<std::vector<std::vector<RKIND>>> Ex;
    std::vector<std::vector<std::vector<RKIND>>> Ey;
    std::vector<std::vector<std::vector<RKIND>>> Ez;
    std::vector<std::vector<std::vector<RKIND>>> Hx;
    std::vector<std::vector<std::vector<RKIND>>> Hy;
    std::vector<std::vector<std::vector<RKIND>>> Hz;
    std::vector<RKIND> dxe;
    std::vector<RKIND> dye;
    std::vector<RKIND> dze;
    std::vector<RKIND> dxh;
    std::vector<RKIND> dyh;
    std::vector<RKIND> dzh;

    void createDummyFields(int lower, int upper, RKIND delta);
};

} // namespace Observa_m

#endif
