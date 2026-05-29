#ifndef MTL_M_H
#define MTL_M_H

#include <optional>
#include <string>
#include <vector>

#include "dispersive_m.h"
#include "mtln_types.h"

namespace mtl_m {

using mtln_types_m::multipolar_expansion_t;
using mtln_types_m::segment_t;
using mtln_types_m::transfer_impedance_per_meter_t;

#ifdef CompileWithMPI
struct communicator_t {
    int field_index = -1;
    int v_index = -1;
    int comm_task = 0;
    int comm_type = 0;
    int delta_rank = 0;
};

struct comm_t {
    std::vector<communicator_t> comms;
    int rank = 0;
};
#endif

struct mtl_t {
    std::string name;
    int number_of_conductors = 0;
    std::vector<std::vector<std::vector<double>>> lpul;
    std::vector<std::vector<std::vector<double>>> cpul;
    std::vector<std::vector<std::vector<double>>> rpul;
    std::vector<std::vector<std::vector<double>>> gpul;
    std::vector<double> step_size;
    std::vector<std::vector<std::vector<double>>> du;
    dispersive_m::lumped_t lumped_elements;
    double time = 0.0;
    double dt = 0.0;
    std::string parent_name;
    int conductor_in_parent = 0;
    transfer_impedance_per_meter_t transfer_impedance;
    std::vector<transfer_impedance_per_meter_t> initial_connector_transfer_impedances;
    std::vector<transfer_impedance_per_meter_t> end_connector_transfer_impedances;
    std::vector<segment_t> segments;

#ifdef CompileWithMPI
    comm_t mpi_comm;
    std::vector<std::vector<int>> layer_indices;
    bool bundle_in_layer = true;
#endif

    void setTimeStep(int numberOfSteps, double finalTime);
    void checkTimeStep(bool getMax, std::optional<double> dt = std::nullopt);
    void allocatePULMatrices();
    void computeLCParameters(const multipolar_expansion_t& multipolar_expansion);
    void computeLCParametersFromRadius(double rad);
    void initLC(const std::vector<std::vector<double>>& lpul, const std::vector<std::vector<double>>& cpul);
    void initRG(const std::vector<std::vector<double>>& rpul, const std::vector<std::vector<double>>& gpul);
    void initDirections();
    double getMaxTimeStep();
    std::vector<std::vector<double>> getPhaseVelocities();

#ifdef CompileWithMPI
    void initCommunicators(const std::vector<int>& alloc_z);
    void initStepSizeAndFieldSegments(const std::vector<double>& step_size,
                                      const std::vector<segment_t>& segments,
                                      const std::vector<std::vector<int>>& layer_indices);
#endif
};

struct transmission_line_level_t {
    std::vector<mtl_t> lines;
};

struct transmission_line_bundle_t {
    std::vector<transmission_line_level_t> levels;
};

mtl_t mtl_shielded(const std::vector<std::vector<double>>& lpul,
                   const std::vector<std::vector<double>>& cpul,
                   const std::vector<std::vector<double>>& rpul,
                   const std::vector<std::vector<double>>& gpul,
                   const std::vector<double>& step_size,
                   const std::string& name,
                   const std::vector<segment_t>& segments,
                   double dt,
                   const std::string& parent_name,
                   int conductor_in_parent,
                   const transfer_impedance_per_meter_t& transfer_impedance,
                   std::optional<std::vector<std::vector<int>>> layer_indices = std::nullopt,
                   std::optional<bool> bundle_in_layer = std::nullopt,
                   std::optional<std::vector<int>> alloc_z = std::nullopt);

mtl_t mtl_unshielded(const std::vector<std::vector<double>>& lpul,
                     const std::vector<std::vector<double>>& cpul,
                     const std::vector<std::vector<double>>& rpul,
                     const std::vector<std::vector<double>>& gpul,
                     const std::vector<double>& step_size,
                     const std::string& name,
                     const std::vector<segment_t>& segments,
                     double dt,
                     const std::vector<multipolar_expansion_t>& multipolar_expansion,
                     double radius,
                     std::optional<std::vector<std::vector<int>>> layer_indices = std::nullopt,
                     std::optional<bool> bundle_in_layer = std::nullopt,
                     std::optional<std::vector<int>> alloc_z = std::nullopt);

void checkPULDimensions(const std::vector<std::vector<double>>& lpul,
                        const std::vector<std::vector<double>>& cpul,
                        const std::vector<std::vector<double>>& rpul,
                        const std::vector<std::vector<double>>& gpul);

} // namespace mtl_m

#endif
