#ifndef MTLN_SOLVER_M_H
#define MTLN_SOLVER_M_H

#include <array>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

#include "fdetypes_mtln.h"
#include "mtl_bundle_m.h"
#include "mtln_preprocess_m.h"
#include "network_manager_m.h"

namespace mtln_solver_m {

using mtl_bundle_m::mtl_bundle_t;
using mtln_preprocess_m::parsed_mtln_t;
using mtln_preprocess_m::preprocess_t;
using network_manager_m::network_manager_t;

struct mtln_t {
    double time = 0.0;
    double dt = 0.0;
    double final_time = 0.0;
    std::vector<mtl_bundle_t> bundles;
    network_manager_t network_manager;
    int number_of_bundles = 0;
    int number_of_steps = 0;
    double null_field = 0.0;

    void updateBundlesTimeStep(double dt_in);
    void updatePULTerms();
    void initNodes();
    int getTimeRange() const;
    int getTimeRange(double time) const;
    void updateProbes();
    void advanceNWVoltage();
    void advanceBundlesVoltage();
    void advanceBundlesCurrent();
    void advanceTime();
    void step();
    void step_alone();
    void setExternalLongitudinalField();
    void runUntil(double final_time_in);
    void run();
    void initObservation(const std::string& nEntradaRoot);
    void updateObservation(int step);
    void closeObservation();
};

mtln_t mtlnCtor(const parsed_mtln_t& parsed);
mtln_t mtlnCtor(const parsed_mtln_t& parsed, const std::array<FDETYPES_m::XYZlimit_t, 6>& alloc);

} // namespace mtln_solver_m

#endif
