#include "observation_movie_test.h"

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <sstream>

#include "nfde_types.h"

namespace Observa_m {

namespace {

std::vector<output_movie_t> g_output;

std::string trim(const std::string& s) {
    const auto start = s.find_first_not_of(' ');
    if (start == std::string::npos) {
        return {};
    }
    const auto end = s.find_last_not_of(' ');
    return s.substr(start, end - start + 1);
}

std::string format_i7(int value) {
    std::ostringstream oss;
    oss << std::setw(7) << std::right << value;
    return trim(oss.str());
}

int fieldo(int field, char dir) {
    using NFDETypes_m::iEx;
    using NFDETypes_m::iEy;
    using NFDETypes_m::iEz;
    using NFDETypes_m::iHx;
    using NFDETypes_m::iHy;
    using NFDETypes_m::iHz;
    using NFDETypes_m::iMEC;
    using NFDETypes_m::iMHC;

    switch (field) {
    case iEx:
    case iEy:
    case iEz:
    case iHx:
    case iHy:
    case iHz:
        return field;
    case iMEC:
        switch (dir) {
        case 'X':
        case 'x':
            return iEx;
        case 'Y':
        case 'y':
            return iEy;
        case 'Z':
        case 'z':
            return iEz;
        default:
            return -1;
        }
    case iMHC:
        switch (dir) {
        case 'X':
        case 'x':
            return iHx;
        case 'Y':
        case 'y':
            return iHy;
        case 'Z':
        case 'z':
            return iHz;
        default:
            return -1;
        }
    default:
        return -1;
    }
}

std::string prefix_field(int field) {
    using NFDETypes_m::iMEC;
    if (field == iMEC) {
        return "ME_";
    }
    return "ME_";
}

std::string build_volumic_path(const SGGFDTDINFO_t& sgg, const sim_control_t& control,
                               const observable_t& probe, int field) {
    observable_t adjusted = probe;
  adjusted.ZI = std::max(sgg.Sweep[static_cast<size_t>(NFDETypes_m::iEx - 1)].ZI, adjusted.ZI);
  adjusted.ZE = std::min(sgg.Sweep[static_cast<size_t>(NFDETypes_m::iEx - 1)].ZE, adjusted.ZE);

    const std::string chari = format_i7(adjusted.XI);
    const std::string charj = format_i7(adjusted.YI);
    const std::string chark = format_i7(adjusted.ZI);
    const std::string chari2 = format_i7(adjusted.XE);
    const std::string charj2 = format_i7(adjusted.YE);
    const std::string chark2 = format_i7(adjusted.ZE);

    std::string extpoint;
    if (control.mpidir == 3) {
        extpoint = chari + '_' + charj + '_' + chark + "__" + chari2 + '_' + charj2 + '_' + chark2;
    } else if (control.mpidir == 2) {
        extpoint = charj + '_' + chark + '_' + chari + "__" + charj2 + '_' + chark2 + '_' + chari2;
    } else if (control.mpidir == 1) {
        extpoint = chark + '_' + chari + '_' + charj + "__" + chark2 + '_' + chari2 + '_' + charj2;
    }

    const std::string ext =
        trim(control.nEntradaRoot) + '_' + trim(sgg.Observation[0].outputrequest);
    return trim(ext) + '_' + prefix_field(field) + extpoint + ".bin";
}

} // namespace

std::string get_temp_dir() {
    const char* vars[] = {"TMPDIR", "TEMP", "TMP"};
    for (const char* var : vars) {
        if (const char* value = std::getenv(var)) {
            if (value[0] != '\0') {
                return value;
            }
        }
    }
    return "/tmp";
}

XYZlimit_t create_xyz_limit(int xi, int yi, int zi, int xe, int ye, int ze) {
    XYZlimit_t r;
    r.XI = xi;
    r.YI = yi;
    r.ZI = zi;
    r.XE = xe;
    r.YE = ye;
    r.ZE = ze;
    return r;
}

limit_t create_limit(int xi, int xe, int yi, int ye, int zi, int ze, int nx, int ny, int nz) {
    limit_t r;
    r.XI = xi;
    r.XE = xe;
    r.YI = yi;
    r.YE = ye;
    r.ZI = zi;
    r.ZE = ze;
    r.NX = nx;
    r.NY = ny;
    r.NZ = nz;
    return r;
}

SGGFDTDINFO_t create_base_sgg() {
    SGGFDTDINFO_t sgg;
    sgg.NumMedia = 3;
    sgg.NumberRequest = 1;
    sgg.dt = 0.1;
    sgg.tiempo.resize(100);
    for (int i = 0; i < 100; ++i) {
        sgg.tiempo[static_cast<size_t>(i)] = static_cast<RKIND_tiempo>(i) * sgg.dt;
    }
    sgg.Sweep.assign(6, create_xyz_limit(0, 0, 0, 6, 6, 6));
    sgg.SINPMLSweep.assign(6, create_xyz_limit(1, 1, 1, 5, 5, 5));
    sgg.alloc.assign(6, create_xyz_limit(0, 0, 0, 6, 6, 6));
    return sgg;
}

Obses_movie_t define_time_movie_observation() {
    Obses_movie_t obs;
    obs.nP = 1;
    obs.P.resize(1);
    obs.P[0].XI = 1;
    obs.P[0].YI = 1;
    obs.P[0].ZI = 1;
    obs.P[0].XE = 6;
    obs.P[0].YE = 6;
    obs.P[0].ZE = 6;
    obs.P[0].Xtrancos = 1;
    obs.P[0].Ytrancos = 1;
    obs.P[0].Ztrancos = 1;
    obs.P[0].What = NFDETypes_m::iMEC;
    obs.InitialTime = 0.0;
    obs.FinalTime = 1.0;
    obs.TimeStep = 0.1;
    obs.InitialFreq = 0.0;
    obs.FinalFreq = 0.0;
    obs.FreqStep = 0.0;
    obs.outputrequest = "timeMovie";
    obs.FreqDomain = false;
    obs.TimeDomain = true;
    obs.Saveall = false;
    obs.TransFer = false;
    obs.Volumic = true;
    obs.Done = false;
    obs.Begun = false;
    obs.Flushed = false;
    return obs;
}

media_matrices_t create_media(const std::vector<XYZlimit_t>&) {
    return {};
}

taglist_t create_tag_list(const std::vector<XYZlimit_t>&) {
    return {};
}

nf2ff_t create_faces_nf2ff(bool tr, bool fr, bool iz, bool de, bool ab, bool ar) {
    nf2ff_t faces;
    faces.tr = tr;
    faces.fr = fr;
    faces.iz = iz;
    faces.de = de;
    faces.ab = ab;
    faces.ar = ar;
    return faces;
}

sim_control_t create_control_flags(int layoutnumber, int num_procs, int mpidir, int finaltimestep,
                                   const std::string& nEntradaRoot, const std::string& wiresflavor,
                                   bool resume, bool saveall, bool nf2ffDecim, bool simu_devia,
                                   bool singlefilewrite, const nf2ff_t& facesNF2FF) {
    sim_control_t control;
    control.layoutnumber = layoutnumber;
    control.num_procs = num_procs;
    control.mpidir = mpidir;
    control.finalTimeStep = finaltimestep;
    control.nEntradaRoot = nEntradaRoot;
    control.wiresflavor = wiresflavor;
    control.resume = resume;
    control.saveall = saveall;
    control.NF2FFDecim = nf2ffDecim;
    control.simu_devia = simu_devia;
    control.singleFileWrite = singlefilewrite;
    control.facesNF2FF = facesNF2FF;
    return control;
}

void InitObservation(SGGFDTDINFO_t& sgg, media_matrices_t&, taglist_t&,
                     bool& ThereAreObservation, bool&, bool&,
                     int&, RKIND_tiempo&, const std::vector<limit_t>&, RKIND, RKIND,
                     bounds_t&, sim_control_t& control) {
    g_output.clear();
    if (sgg.Observation.empty()) {
        return;
    }

    output_movie_t out;
    out.timeswritten = 0;
    out.item.resize(1);

    const int field = sgg.Observation[0].P[0].What;
    out.item[0].path = build_volumic_path(sgg, control, sgg.Observation[0].P[0], field);
    out.item[0].unit = 1001;
    g_output.push_back(out);
    ThereAreObservation = true;
    (void)control;
}

void UpdateObservation(SGGFDTDINFO_t&, media_matrices_t&, taglist_t&, int, int,
                       const std::vector<std::vector<std::vector<RKIND>>>&,
                       const std::vector<std::vector<std::vector<RKIND>>>&,
                       const std::vector<std::vector<std::vector<RKIND>>>&,
                       const std::vector<std::vector<std::vector<RKIND>>>&,
                       const std::vector<std::vector<std::vector<RKIND>>>&,
                       const std::vector<std::vector<std::vector<RKIND>>>&,
                       const std::vector<RKIND>&, const std::vector<RKIND>&,
                       const std::vector<RKIND>&, const std::vector<RKIND>&,
                       const std::vector<RKIND>&, const std::vector<RKIND>&,
                       const std::string&, const std::vector<limit_t>&, bool, bool, bounds_t&) {
    if (!g_output.empty()) {
        ++g_output[0].timeswritten;
    }
}

std::vector<output_movie_t>& GetOutput() {
    return g_output;
}

void dummyFields_t::createDummyFields(int lower, int upper, RKIND delta) {
    const size_t n = static_cast<size_t>(upper - lower + 1);
    auto make_cube = [n]() {
        return std::vector<std::vector<std::vector<RKIND>>>(
            n, std::vector<std::vector<RKIND>>(n, std::vector<RKIND>(n, 0.0)));
    };
    Ex = make_cube();
    Ey = make_cube();
    Ez = make_cube();
    Hx = make_cube();
    Hy = make_cube();
    Hz = make_cube();
    dxe.assign(n, delta);
    dye.assign(n, delta);
    dze.assign(n, delta);
    dxh.assign(n, delta);
    dyh.assign(n, delta);
    dzh.assign(n, delta);
    (void)lower;
}

} // namespace Observa_m
