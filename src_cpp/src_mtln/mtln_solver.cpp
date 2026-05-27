#include "mtln_solver_m.h"

#include <cmath>
#include <fstream>
#include <iostream>

#include "utils_m.h"

#ifdef CompileWithMPI
#include <mpi.h>
extern MPI_Comm SUBCOMM_MPI;
#endif

namespace mtln_solver_m {

using mtln_types_m::TERMINAL_NODE_SIDE_END;
using mtln_types_m::TERMINAL_NODE_SIDE_INI;

mtln_t mtlnCtor(const parsed_mtln_t& parsed) {
    mtln_t res;
    preprocess_t pre = mtln_preprocess_m::preprocess(parsed);
    if (pre.bundles.empty()) {
        res.number_of_bundles = 0;
        return res;
    }
    res.dt = pre.dt;
    res.time = 0.0;
    res.final_time = pre.final_time;
    res.bundles = std::move(pre.bundles);
    res.number_of_bundles = static_cast<int>(res.bundles.size());
    res.network_manager = pre.network_manager;
    res.number_of_steps = parsed.number_of_steps;
    res.null_field = 0.0;
    res.updateBundlesTimeStep(res.dt);
    res.initNodes();
    return res;
}

mtln_t mtlnCtor(const parsed_mtln_t& parsed, const std::array<FDETYPES_m::XYZlimit_t, 6>& alloc) {
    mtln_t res;
#ifdef CompileWithMPI
    MPI_Barrier(SUBCOMM_MPI);
#endif
    preprocess_t pre = mtln_preprocess_m::preprocess(parsed, alloc);
    if (pre.bundles.empty()) {
        res.number_of_bundles = 0;
        return res;
    }
    res.dt = pre.dt;
    res.time = 0.0;
    res.final_time = pre.final_time;
    res.bundles = std::move(pre.bundles);
    res.number_of_bundles = static_cast<int>(res.bundles.size());
    res.network_manager = pre.network_manager;
    res.number_of_steps = parsed.number_of_steps;
    res.null_field = 0.0;
    res.updateBundlesTimeStep(res.dt);
    res.initNodes();
    return res;
}

void mtln_t::initNodes() {
    for (auto& network : network_manager.networks) {
        for (auto& node : network.nodes) {
            node.v = 0.0;
            node.i = 0.0;
        }
    }
}

void mtln_t::step() {
    setExternalLongitudinalField();
    advanceBundlesVoltage();
    advanceNWVoltage();
    advanceBundlesCurrent();
    updateProbes();
    advanceTime();
}

void mtln_t::step_alone() {
    advanceBundlesVoltage();
    advanceNWVoltage();
    advanceBundlesCurrent();
    updateProbes();
    advanceTime();
}

void mtln_t::setExternalLongitudinalField() {
#ifdef CompileWithMPI
    MPI_Barrier(SUBCOMM_MPI);
#endif
    for (int i = 0; i < number_of_bundles; ++i) {
        if (bundles[static_cast<size_t>(i)].bundle_in_layer) {
            bundles[static_cast<size_t>(i)].setExternalLongitudinalField();
        }
    }
}

void mtln_t::advanceBundlesVoltage() {
    for (int i = 0; i < number_of_bundles; ++i) {
        if (bundles[static_cast<size_t>(i)].bundle_in_layer) {
            bundles[static_cast<size_t>(i)].updateGenerators(time, dt);
            bundles[static_cast<size_t>(i)].advanceVoltage();
        }
    }
}

void mtln_t::advanceNWVoltage() {
    if (number_of_bundles == 0) {
        return;
    }
    for (auto& network : network_manager.networks) {
        for (auto& node : network.nodes) {
            const int b = node.bundle_number - 1;
            const int c = node.conductor_number - 1;
            const int i_idx = node.i_index;
            if (b >= 0 && b < number_of_bundles && bundles[static_cast<size_t>(b)].bundle_in_layer &&
                i_idx >= 0 &&
                i_idx < static_cast<int>(bundles[static_cast<size_t>(b)].i[static_cast<size_t>(c)].size())) {
                node.i = bundles[static_cast<size_t>(b)].i[static_cast<size_t>(c)][static_cast<size_t>(i_idx)];
            }
        }
    }

    network_manager.advanceVoltage();

    for (auto& network : network_manager.networks) {
        for (auto& node : network.nodes) {
            const int b = node.bundle_number - 1;
            const int c = node.conductor_number - 1;
            const int v_idx = node.v_index;
            if (b < 0 || b >= number_of_bundles) {
                continue;
            }
            if (!node.open) {
                if (bundles[static_cast<size_t>(b)].bundle_in_layer && v_idx >= 0 &&
                    v_idx < static_cast<int>(bundles[static_cast<size_t>(b)].v[static_cast<size_t>(c)].size())) {
                    bundles[static_cast<size_t>(b)].v[static_cast<size_t>(c)][static_cast<size_t>(v_idx)] = node.v;
                }
            } else if (node.side == TERMINAL_NODE_SIDE_INI) {
                std::vector<double> idiff_row(
                    bundles[static_cast<size_t>(b)].i_diff[0][static_cast<size_t>(c)].begin(),
                    bundles[static_cast<size_t>(b)].i_diff[0][static_cast<size_t>(c)].end());
                std::vector<double> i_row;
                for (int k = 0; k < bundles[static_cast<size_t>(b)].number_of_conductors; ++k) {
                    i_row.push_back(bundles[static_cast<size_t>(b)].i[static_cast<size_t>(k)][0]);
                }
                bundles[static_cast<size_t>(b)].v[static_cast<size_t>(c)][0] -=
                    2.0 * utils_m::dot_product(idiff_row, i_row);
            } else if (node.side == TERMINAL_NODE_SIDE_END) {
                const int n = bundles[static_cast<size_t>(b)].number_of_divisions;
                std::vector<double> idiff_row(
                    bundles[static_cast<size_t>(b)].i_diff[static_cast<size_t>(n)][static_cast<size_t>(c)].begin(),
                    bundles[static_cast<size_t>(b)].i_diff[static_cast<size_t>(n)][static_cast<size_t>(c)].end());
                std::vector<double> i_row;
                for (int k = 0; k < bundles[static_cast<size_t>(b)].number_of_conductors; ++k) {
                    i_row.push_back(bundles[static_cast<size_t>(b)].i[static_cast<size_t>(k)][static_cast<size_t>(n - 1)]);
                }
                bundles[static_cast<size_t>(b)].v[static_cast<size_t>(c)][static_cast<size_t>(n)] +=
                    2.0 * utils_m::dot_product(idiff_row, i_row);
            }
        }
    }
}

void mtln_t::advanceBundlesCurrent() {
#ifdef CompileWithMPI
    MPI_Barrier(SUBCOMM_MPI);
#endif
    for (int i = 0; i < number_of_bundles; ++i) {
        if (bundles[static_cast<size_t>(i)].bundle_in_layer) {
            bundles[static_cast<size_t>(i)].advanceCurrent();
        }
    }
}

void mtln_t::advanceTime() {
    time += dt;
}

void mtln_t::updateProbes() {
    for (int i = 0; i < number_of_bundles; ++i) {
        auto& bundle = bundles[static_cast<size_t>(i)];
        if (!bundle.probes.empty() && bundle.bundle_in_layer) {
            for (auto& probe : bundle.probes) {
                if (probe.in_layer) {
                    probe.update(time, bundle.v, bundle.i);
                }
            }
        }
    }
}

int mtln_t::getTimeRange() const {
    return static_cast<int>(std::floor(final_time / dt));
}

int mtln_t::getTimeRange(double time_in) const {
    return static_cast<int>(std::floor(time_in / dt));
}

void mtln_t::updateBundlesTimeStep(double dt_in) {
    for (auto& bundle : bundles) {
        bundle.dt = dt_in;
    }
}

void mtln_t::updatePULTerms() {
    for (auto& bundle : bundles) {
        if (bundle.bundle_in_layer) {
            bundle.updateLRTerms();
            bundle.updateCGTerms();
            for (auto& probe : bundle.probes) {
                probe.resizeFrames(getTimeRange(), bundle.number_of_conductors);
            }
        }
    }
}

void mtln_t::runUntil(double final_time_in) {
    const int max_steps = getTimeRange(final_time_in);
    for (int step_idx = 0; step_idx <= max_steps; ++step_idx) {
        advanceBundlesVoltage();
        advanceNWVoltage();
        advanceBundlesCurrent();
        updateProbes();
        advanceTime();
        updateObservation(step_idx);
    }
}

void mtln_t::run() {
    const int max_steps = getTimeRange();
    for (int step_idx = 0; step_idx <= max_steps; ++step_idx) {
        advanceBundlesVoltage();
        advanceNWVoltage();
        advanceBundlesCurrent();
        updateProbes();
        advanceTime();
        updateObservation(step_idx);
    }
}

void mtln_t::initObservation(const std::string& nEntradaRoot) {
    if (bundles.empty()) {
        return;
    }
    int unit = 2000;
    for (auto& bundle : bundles) {
        for (auto& probe : bundle.probes) {
#ifdef CompileWithMPI
            if (!probe.in_layer) {
                continue;
            }
#endif
            probe.unit = unit++;
            probe.output_path = nEntradaRoot + "_" + probe.name + ".dat";
            std::ofstream out(probe.output_path);
            if (!out.is_open()) {
                continue;
            }
            std::string header = "time";
            const int n_conductors =
                probe.val.empty() ? bundle.number_of_conductors : static_cast<int>(probe.val[0].size());
            for (int k = 0; k < n_conductors; ++k) {
                header += " conductor_" + std::to_string(k + 1);
            }
            out << header << '\n';
        }
    }
}

void mtln_t::updateObservation(int) {
    for (auto& bundle : bundles) {
        for (auto& probe : bundle.probes) {
#ifdef CompileWithMPI
            if (!probe.in_layer) {
                continue;
            }
#endif
            if (probe.output_path.empty()) {
                continue;
            }
            std::ofstream out(probe.output_path, std::ios::app);
            if (!out.is_open()) {
                continue;
            }
            out << probe.t[0];
            for (const auto& val : probe.val[0]) {
                out << ' ' << val;
            }
            out << '\n';
        }
    }
}

void mtln_t::closeObservation() {
    for (auto& bundle : bundles) {
        for (auto& probe : bundle.probes) {
#ifdef CompileWithMPI
            if (!probe.in_layer) {
                continue;
            }
#endif
            probe.output_path.clear();
        }
    }
}

} // namespace mtln_solver_m
