#include "mtl_bundle_m.h"

namespace mtl_bundle_m {

namespace {

void copyBlock(std::vector<std::vector<std::vector<double>>>& dest, int destOffset,
               const std::vector<std::vector<std::vector<double>>>& src, int srcN) {
    const int nDiv = static_cast<int>(dest.size());
    for (int seg = 0; seg < nDiv; ++seg) {
        for (int i = 0; i < srcN; ++i) {
            for (int j = 0; j < srcN; ++j) {
                dest[static_cast<size_t>(seg)][static_cast<size_t>(destOffset + i)][static_cast<size_t>(destOffset + j)] =
                    src[static_cast<size_t>(seg)][static_cast<size_t>(i)][static_cast<size_t>(j)];
            }
        }
    }
}

} // namespace

int countNumberOfConductors(const std::vector<mtl_m::transmission_line_level_t>& levels) {
    int res = 0;
    for (const auto& level : levels) {
        for (const auto& line : level.lines) {
            res += line.number_of_conductors;
        }
    }
    return res;
}

std::vector<external_field_segment_t> buildExternalFieldSegments(
    const std::vector<mtl_m::transmission_line_level_t>& levels) {
    std::vector<external_field_segment_t> res;
    if (levels.empty() || levels[0].lines.empty()) {
        return res;
    }
    const auto& segments = levels[0].lines[0].segments;
    res.resize(segments.size());
    for (size_t i = 0; i < segments.size(); ++i) {
        res[i].position = {segments[i].x, segments[i].y, segments[i].z};
        res[i].direction = segments[i].orientation;
    }
    return res;
}

void mtl_bundle_t::initialAllocation() {
    lpul.assign(static_cast<size_t>(number_of_divisions),
                std::vector<std::vector<double>>(static_cast<size_t>(number_of_conductors),
                                                 std::vector<double>(static_cast<size_t>(number_of_conductors), 0.0)));
    cpul.assign(static_cast<size_t>(number_of_divisions + 1),
                std::vector<std::vector<double>>(static_cast<size_t>(number_of_conductors),
                                                 std::vector<double>(static_cast<size_t>(number_of_conductors), 0.0)));
    gpul = cpul;
    rpul = lpul;
    du = lpul;
    v.assign(static_cast<size_t>(number_of_conductors),
             std::vector<double>(static_cast<size_t>(number_of_divisions + 1), 0.0));
    i.assign(static_cast<size_t>(number_of_conductors),
             std::vector<double>(static_cast<size_t>(number_of_divisions), 0.0));
    i_prev = i;
    e_L = i;
    v_source = v;
    i_source = i;
}

void mtl_bundle_t::mergePULMatrices(const std::vector<mtl_m::transmission_line_level_t>& levels) {
    int n_sum = 0;
    for (const auto& level : levels) {
        for (const auto& line : level.lines) {
            const int n = line.number_of_conductors;
            copyBlock(lpul, n_sum, line.lpul, n);
            copyBlock(cpul, n_sum, line.cpul, n);
            copyBlock(rpul, n_sum, line.rpul, n);
            copyBlock(gpul, n_sum, line.gpul, n);
            copyBlock(du, n_sum, line.du, n);
            n_sum += n;
        }
    }
}

void mtl_bundle_t::mergeDispersiveMatrices(const std::vector<mtl_m::transmission_line_level_t>& levels) {
    int number_of_poles = 0;
    for (const auto& level : levels) {
        for (const auto& line : level.lines) {
            number_of_poles = std::max(number_of_poles, line.lumped_elements.number_of_poles);
        }
    }
    transfer_impedance =
        dispersive_m::transfer_impedance_t(number_of_conductors, number_of_poles, number_of_divisions, dt);

    int n_sum = 0;
    for (const auto& level : levels) {
        for (const auto& line : level.lines) {
            const int n = line.number_of_conductors;
            const auto& lumped = line.lumped_elements;
            for (int div = 0; div < number_of_divisions; ++div) {
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        transfer_impedance.q1[static_cast<size_t>(div)][static_cast<size_t>(n_sum + i)]
                            [static_cast<size_t>(n_sum + j)] = lumped.q1[static_cast<size_t>(div)][static_cast<size_t>(i)][static_cast<size_t>(j)];
                        transfer_impedance.q2[static_cast<size_t>(div)][static_cast<size_t>(n_sum + i)]
                            [static_cast<size_t>(n_sum + j)] = lumped.q2[static_cast<size_t>(div)][static_cast<size_t>(i)][static_cast<size_t>(j)];
                        transfer_impedance.q3[static_cast<size_t>(div)][static_cast<size_t>(n_sum + i)]
                            [static_cast<size_t>(n_sum + j)] = lumped.q3[static_cast<size_t>(div)][static_cast<size_t>(i)][static_cast<size_t>(j)];
                        transfer_impedance.q1_sum[static_cast<size_t>(div)][static_cast<size_t>(n_sum + i)]
                            [static_cast<size_t>(n_sum + j)] = lumped.q1_sum[static_cast<size_t>(div)][static_cast<size_t>(i)][static_cast<size_t>(j)];
                        transfer_impedance.q2_sum[static_cast<size_t>(div)][static_cast<size_t>(n_sum + i)]
                            [static_cast<size_t>(n_sum + j)] = lumped.q2_sum[static_cast<size_t>(div)][static_cast<size_t>(i)][static_cast<size_t>(j)];
                        transfer_impedance.d[static_cast<size_t>(div)][static_cast<size_t>(n_sum + i)]
                            [static_cast<size_t>(n_sum + j)] = lumped.d[static_cast<size_t>(div)][static_cast<size_t>(i)][static_cast<size_t>(j)];
                        transfer_impedance.e[static_cast<size_t>(div)][static_cast<size_t>(n_sum + i)]
                            [static_cast<size_t>(n_sum + j)] = lumped.e[static_cast<size_t>(div)][static_cast<size_t>(i)][static_cast<size_t>(j)];
                    }
                    transfer_impedance.q3_phi[static_cast<size_t>(div)][static_cast<size_t>(n_sum + i)] =
                        lumped.q3_phi[static_cast<size_t>(div)][static_cast<size_t>(i)];
                    transfer_impedance.phi[static_cast<size_t>(div)][static_cast<size_t>(n_sum + i)] =
                        lumped.phi[static_cast<size_t>(div)][static_cast<size_t>(i)];
                }
            }
            n_sum += n;
        }
    }
}

void mtl_bundle_t::addTransferImpedance(int conductor_out, const std::vector<int>& range_in,
                                        const mtln_types_m::transfer_impedance_per_meter_t& zt) {
    transfer_impedance.addTransferImpedance(conductor_out, range_in, zt);
}

mtl_bundle_t mtl_bundle_ctor(const std::vector<mtl_m::transmission_line_level_t>& levels,
                             const std::string& name) {
    mtl_bundle_t res;
    res.name = name;
    res.number_of_conductors = countNumberOfConductors(levels);
    if (!levels.empty() && !levels[0].lines.empty()) {
        res.dt = levels[0].lines[0].dt;
        res.step_size = levels[0].lines[0].step_size;
        res.number_of_divisions = static_cast<int>(res.step_size.size());
    }
    res.initialAllocation();
    res.mergePULMatrices(levels);
    res.mergeDispersiveMatrices(levels);
    return res;
}

} // namespace mtl_bundle_m
