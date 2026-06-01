#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "wires_mtln_m.h"

#include "fdetypes_mtln.h"
#include "mtln_solver_m.h"
#include "mtln_types.h"

#ifdef CompileWithMPI
#include <mpi.h>
extern MPI_Comm subcomm_mpi;
#endif

using mtln_solver_m::mtln_t;
using mtln_types_m::DIRECTION_X_POS;
using mtln_types_m::DIRECTION_Y_POS;
using mtln_types_m::DIRECTION_Z_POS;

namespace FDETYPES_m {
struct MediaData_t {
    struct {
        bool pec = false;
        bool lossy = false;
    } is;
};
struct SGGFDTDINFO_t {
    std::vector<XYZlimit_t> Alloc;
    std::vector<MediaData_t> Med;
    double dt = 0.0;
};
} // namespace FDETYPES_m

using FDETYPES_m::SGGFDTDINFO_t;
using mtl_bundle_m::external_field_segment_t;

void print11(int, const std::string& msg) {
    std::cout << msg << std::endl;
}

namespace Wire_bundles_mtln_m {

#ifdef CompileWithMTLN

double eps0 = 0.0;
double mu0 = 0.0;
mtln_t mtln_solver;
std::vector<std::vector<int>> indexMap;

inline bool isEmbeddedInPECorLossy(int media, const SGGFDTDINFO_t& sgg) {
    return media == 0 || sgg.Med[static_cast<size_t>(media)].is.pec ||
           sgg.Med[static_cast<size_t>(media)].is.lossy;
}

inline void readGridIndices(int& i, int& j, int& k, const external_field_segment_t& field_segment) {
    i = field_segment.position[0];
    j = field_segment.position[1];
    k = field_segment.position[2];
}

void InitWires_mtln(
    const SGGFDTDINFO_t& sgg,
    std::vector<std::vector<std::vector<double>>>& Ex,
    std::vector<std::vector<std::vector<double>>>& Ey,
    std::vector<std::vector<std::vector<double>>>& Ez,
    const std::vector<std::vector<std::vector<int>>>& sggMiEx,
    const std::vector<std::vector<std::vector<int>>>& sggMiEy,
    const std::vector<std::vector<std::vector<int>>>& sggMiEz,
    const std::vector<std::vector<std::vector<int>>>&,
    const std::vector<std::vector<std::vector<int>>>&,
    const std::vector<std::vector<std::vector<int>>>&,
    double eps00,
    double mu00,
    const mtln_types_m::mtln_t& mtln_parsed,
    bool& thereAreMTLNbundles,
    double& dtcritico) {
    eps0 = eps00;
    mu0 = mu00;

    std::array<FDETYPES_m::XYZlimit_t, 6> alloc{};
    if (sgg.Alloc.size() >= 6) {
        for (int i = 0; i < 6; ++i) {
            alloc[static_cast<size_t>(i)] = sgg.Alloc[static_cast<size_t>(i)];
        }
    }

#ifdef CompileWithMPI
    mtln_solver = mtln_solver_m::mtlnCtor(mtln_parsed, alloc);
    MPI_Barrier(subcomm_mpi);
#else
    (void)alloc;
    mtln_solver = mtln_solver_m::mtlnCtor(mtln_parsed);
#endif

    if (mtln_solver.number_of_bundles >= 1) {
        thereAreMTLNbundles = true;
    } else {
        thereAreMTLNbundles = false;
        return;
    }

    if (mtln_solver.dt < dtcritico) {
        dtcritico = mtln_solver.dt;
    }

    for (int m = 0; m < mtln_solver.number_of_bundles; ++m) {
        auto& bundle = mtln_solver.bundles[static_cast<size_t>(m)];
        if (!bundle.bundle_in_layer) {
            continue;
        }
        for (size_t n = 0; n < bundle.external_field_segments.size(); ++n) {
            int i = 0;
            int j = 0;
            int k = 0;
            readGridIndices(i, j, k, bundle.external_field_segments[n]);
            const int dir = std::abs(bundle.external_field_segments[n].direction);
            if (dir == DIRECTION_X_POS) {
                if (isEmbeddedInPECorLossy(sggMiEx[static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(k)], sgg)) {
                    bundle.external_field_segments[n].field = &mtln_solver.null_field;
                } else {
                    bundle.external_field_segments[n].field = &Ex[static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(k)];
                }
            } else if (dir == DIRECTION_Y_POS) {
                if (isEmbeddedInPECorLossy(sggMiEy[static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(k)], sgg)) {
                    bundle.external_field_segments[n].field = &mtln_solver.null_field;
                } else {
                    bundle.external_field_segments[n].field = &Ey[static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(k)];
                }
            } else if (dir == DIRECTION_Z_POS) {
                if (isEmbeddedInPECorLossy(sggMiEz[static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(k)], sgg)) {
                    bundle.external_field_segments[n].field = &mtln_solver.null_field;
                } else {
                    bundle.external_field_segments[n].field = &Ez[static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(k)];
                }
            }
        }
    }

    mtln_solver.updatePULTerms();
}

double computeFieldFromCurrent(int m, int n, const SGGFDTDINFO_t& sgg,
                               const std::vector<double>& Idxh, const std::vector<double>& Idyh,
                               const std::vector<double>& Idzh);

double getOrientedCurrent(int m, int n) {
    const auto& bundle = mtln_solver.bundles[static_cast<size_t>(m)];
    const int direction = bundle.external_field_segments[static_cast<size_t>(n)].direction;
    double curr = 0.0;
    const int num_conductors =
        bundle.conductors_in_level.empty() ? bundle.number_of_conductors : bundle.conductors_in_level[0];
    for (int i = 0; i < num_conductors; ++i) {
        curr += bundle.i[static_cast<size_t>(i)][static_cast<size_t>(n)];
    }
    return curr * std::copysign(1.0, static_cast<double>(direction));
}

double computeFieldFromCurrent(int m, int n, const SGGFDTDINFO_t& sgg, const std::vector<double>& Idxh,
                               const std::vector<double>& Idyh, const std::vector<double>& Idzh) {
    int i = 0;
    int j = 0;
    int k = 0;
    readGridIndices(i, j, k, mtln_solver.bundles[static_cast<size_t>(m)].external_field_segments[static_cast<size_t>(n)]);
    const int direction =
        mtln_solver.bundles[static_cast<size_t>(m)].external_field_segments[static_cast<size_t>(n)].direction;
    double dS_inverse = 0.0;
    const int abs_dir = std::abs(direction);
    if (abs_dir == 1) {
        dS_inverse = Idyh[static_cast<size_t>(j)] * Idzh[static_cast<size_t>(k)];
    } else if (abs_dir == 2) {
        dS_inverse = Idxh[static_cast<size_t>(i)] * Idzh[static_cast<size_t>(k)];
    } else if (abs_dir == 3) {
        dS_inverse = Idxh[static_cast<size_t>(i)] * Idyh[static_cast<size_t>(j)];
    }
    const double factor = (sgg.dt / eps0) * dS_inverse;
    return factor * getOrientedCurrent(m, n);
}

void AdvanceWiresE_mtln(const SGGFDTDINFO_t& sgg, const std::vector<double>& Idxh,
                        const std::vector<double>& Idyh, const std::vector<double>& Idzh, double eps00,
                        double mu00) {
    eps0 = eps00;
    mu0 = mu00;
    mtln_solver.null_field = 0.0;

    for (int m = 0; m < mtln_solver.number_of_bundles; ++m) {
        auto& bundle = mtln_solver.bundles[static_cast<size_t>(m)];
        if (!bundle.bundle_in_layer) {
            continue;
        }
        for (size_t n = 0; n < bundle.external_field_segments.size(); ++n) {
            double* punt = bundle.external_field_segments[n].field;
            if (punt != nullptr) {
                *punt -= computeFieldFromCurrent(m, static_cast<int>(n), sgg, Idxh, Idyh, Idzh);
            }
        }
    }
    mtln_solver.null_field = 0.0;
    mtln_solver.step();
}

mtln_t* GetSolverPtr() {
    return &mtln_solver;
}

void solveMTLNProblem(const mtln_types_m::mtln_t& mtln_parsed, const std::string& nEntradaRoot) {
    if (mtln_solver.number_of_bundles > 0) {
        mtln_solver.network_manager.circuit.quit();
    }
    mtln_solver = mtln_solver_m::mtlnCtor(mtln_parsed);
    mtln_solver.updatePULTerms();
    mtln_solver.initObservation(nEntradaRoot);
    mtln_solver.run();
    mtln_solver.closeObservation();
    if (mtln_solver.number_of_bundles > 0) {
        mtln_solver.network_manager.circuit.quit();
    }
}

void reportSimulationEnd(int layoutnumber) {
    print11(layoutnumber, "MTLN simulation finished.");
}

#endif

} // namespace Wire_bundles_mtln_m

#ifdef CompileWithMTLN
void InitWires_mtln(const SGGFDTDINFO_t& sgg, std::vector<std::vector<std::vector<double>>>& Ex,
                    std::vector<std::vector<std::vector<double>>>& Ey,
                    std::vector<std::vector<std::vector<double>>>& Ez,
                    const std::vector<std::vector<std::vector<int>>>& sggMiEx,
                    const std::vector<std::vector<std::vector<int>>>& sggMiEy,
                    const std::vector<std::vector<std::vector<int>>>& sggMiEz,
                    const std::vector<std::vector<std::vector<int>>>& sggMiHx,
                    const std::vector<std::vector<std::vector<int>>>& sggMiHy,
                    const std::vector<std::vector<std::vector<int>>>& sggMiHz, double eps00, double mu00,
                    const mtln_types_m::mtln_t& mtln_parsed, bool& thereAreMTLNbundles, double& dtcritico) {
    Wire_bundles_mtln_m::InitWires_mtln(sgg, Ex, Ey, Ez, sggMiEx, sggMiEy, sggMiEz, sggMiHx, sggMiHy, sggMiHz,
                                        eps00, mu00, mtln_parsed, thereAreMTLNbundles, dtcritico);
}

void AdvanceWiresE_mtln(const SGGFDTDINFO_t& sgg, const std::vector<double>& Idxh,
                        const std::vector<double>& Idyh, const std::vector<double>& Idzh, double eps00,
                        double mu00) {
    Wire_bundles_mtln_m::AdvanceWiresE_mtln(sgg, Idxh, Idyh, Idzh, eps00, mu00);
}
#endif
