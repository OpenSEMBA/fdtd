#ifndef MTL_BUNDLE_M_H
#define MTL_BUNDLE_M_H

#include <string>
#include <vector>

#include "dispersive_m.h"
#include "mtl_m.h"

namespace mtl_bundle_m {

struct external_field_segment_t {
    std::vector<int> position = {0, 0, 0};
    int direction = 0;
};

struct mtl_bundle_t {
    std::string name;
    std::vector<std::vector<std::vector<double>>> lpul;
    std::vector<std::vector<std::vector<double>>> cpul;
    std::vector<std::vector<std::vector<double>>> rpul;
    std::vector<std::vector<std::vector<double>>> gpul;
    int number_of_conductors = 0;
    int number_of_divisions = 0;
    std::vector<double> step_size;
    std::vector<std::vector<double>> v;
    std::vector<std::vector<double>> i;
    std::vector<std::vector<double>> i_prev;
    std::vector<std::vector<double>> v_source;
    std::vector<std::vector<double>> i_source;
    std::vector<std::vector<double>> e_L;
    std::vector<std::vector<std::vector<double>>> du;
    double time = 0.0;
    double dt = 0.0;
    dispersive_m::transfer_impedance_t transfer_impedance;

    void initialAllocation();
    void mergePULMatrices(const std::vector<mtl_m::transmission_line_level_t>& levels);
    void mergeDispersiveMatrices(const std::vector<mtl_m::transmission_line_level_t>& levels);
    void addTransferImpedance(int conductor_out, const std::vector<int>& range_in,
                              const mtln_types_m::transfer_impedance_per_meter_t& transfer_impedance);
};

mtl_bundle_t mtl_bundle_ctor(const std::vector<mtl_m::transmission_line_level_t>& levels,
                             const std::string& name = {});

int countNumberOfConductors(const std::vector<mtl_m::transmission_line_level_t>& levels);
std::vector<external_field_segment_t> buildExternalFieldSegments(
    const std::vector<mtl_m::transmission_line_level_t>& levels);

} // namespace mtl_bundle_m

#endif
