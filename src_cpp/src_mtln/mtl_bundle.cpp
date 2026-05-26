#include "mtl_bundle_m.h"

#include "utils_m.h"

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
    res.external_field_segments = buildExternalFieldSegments(levels);
    res.initialAllocation();
    res.mergePULMatrices(levels);
    res.mergeDispersiveMatrices(levels);
#ifdef CompileWithMPI
    if (!levels.empty() && !levels[0].lines.empty()) {
        res.bundle_in_layer = levels[0].lines[0].bundle_in_layer;
        res.layer_indices = levels[0].lines[0].layer_indices;
        res.mpi_comm = levels[0].lines[0].mpi_comm;
    }
#endif
    return res;
}

namespace {

std::vector<std::vector<double>> matrix_at(
    const std::vector<std::vector<std::vector<double>>>& arr, int seg) {
    return arr[static_cast<size_t>(seg)];
}

} // namespace

void mtl_bundle_t::updateLRTerms() {
    const int nd = number_of_divisions;
    const int nc = number_of_conductors;
    i_term.assign(static_cast<size_t>(nd),
                  std::vector<std::vector<double>>(static_cast<size_t>(nc),
                                                   std::vector<double>(static_cast<size_t>(nc), 0.0)));
    v_diff = i_term;

    for (int seg = 0; seg < nd; ++seg) {
        const auto& du_seg = matrix_at(du, seg);
        const auto& lpul_seg = matrix_at(lpul, seg);
        const auto& rpul_seg = matrix_at(rpul, seg);
        std::vector<std::vector<double>> b_plus(nc, std::vector<double>(nc, 0.0));
        std::vector<std::vector<double>> b_minus(nc, std::vector<double>(nc, 0.0));
        for (int i = 0; i < nc; ++i) {
            for (int j = 0; j < nc; ++j) {
                const double d_ij = transfer_impedance.d[static_cast<size_t>(seg)][static_cast<size_t>(i)][static_cast<size_t>(j)];
                const double e_ij = transfer_impedance.e[static_cast<size_t>(seg)][static_cast<size_t>(i)][static_cast<size_t>(j)];
                const double q1_ij = std::real(
                    transfer_impedance.q1_sum[static_cast<size_t>(seg)][static_cast<size_t>(i)][static_cast<size_t>(j)]);
                b_plus[static_cast<size_t>(i)][static_cast<size_t>(j)] =
                    lpul_seg[static_cast<size_t>(i)][static_cast<size_t>(j)] / dt + 0.5 * d_ij + e_ij / dt +
                    0.5 * rpul_seg[static_cast<size_t>(i)][static_cast<size_t>(j)] + q1_ij;
                b_minus[static_cast<size_t>(i)][static_cast<size_t>(j)] =
                    lpul_seg[static_cast<size_t>(i)][static_cast<size_t>(j)] / dt - 0.5 * d_ij + e_ij / dt -
                    0.5 * rpul_seg[static_cast<size_t>(i)][static_cast<size_t>(j)] - q1_ij;
            }
        }
        const auto f1 = utils_m::matmul(du_seg, b_plus);
        const auto f2 = utils_m::matmul(du_seg, b_minus);
        v_diff[static_cast<size_t>(seg)] = utils_m::inv(f1);
        i_term[static_cast<size_t>(seg)] = utils_m::matmul(v_diff[static_cast<size_t>(seg)], f2);
    }
}

void mtl_bundle_t::updateCGTerms() {
    const int nd = number_of_divisions;
    const int nc = number_of_conductors;
    const int ncp1 = nd + 1;
    v_term.assign(static_cast<size_t>(ncp1),
                  std::vector<std::vector<double>>(static_cast<size_t>(nc),
                                                   std::vector<double>(static_cast<size_t>(nc), 0.0)));
    i_diff.assign(static_cast<size_t>(ncp1),
                  std::vector<std::vector<double>>(static_cast<size_t>(nc),
                                                   std::vector<double>(static_cast<size_t>(nc), 0.0)));

    for (int seg = 0; seg < ncp1; ++seg) {
        std::vector<std::vector<double>> ext_du(nc, std::vector<double>(nc, 0.0));
        if (seg == 0) {
            ext_du = matrix_at(du, 0);
        } else if (seg == nd) {
            ext_du = matrix_at(du, nd - 1);
        } else {
            const auto& du_a = matrix_at(du, seg - 1);
            const auto& du_b = matrix_at(du, seg);
            for (int i = 0; i < nc; ++i) {
                for (int j = 0; j < nc; ++j) {
                    ext_du[static_cast<size_t>(i)][static_cast<size_t>(j)] =
                        0.5 * (du_a[static_cast<size_t>(i)][static_cast<size_t>(j)] +
                               du_b[static_cast<size_t>(i)][static_cast<size_t>(j)]);
                }
            }
        }
        std::vector<std::vector<double>> lhs(nc, std::vector<double>(nc, 0.0));
        std::vector<std::vector<double>> rhs(nc, std::vector<double>(nc, 0.0));
        const auto& cpul_seg = matrix_at(cpul, seg);
        const auto& gpul_seg = matrix_at(gpul, seg);
        for (int i = 0; i < nc; ++i) {
            for (int j = 0; j < nc; ++j) {
                double sum_lhs = 0.0;
                double sum_rhs = 0.0;
                for (int k = 0; k < nc; ++k) {
                    sum_lhs += ext_du[static_cast<size_t>(i)][static_cast<size_t>(k)] *
                               (cpul_seg[static_cast<size_t>(k)][static_cast<size_t>(j)] / dt);
                    sum_rhs += ext_du[static_cast<size_t>(i)][static_cast<size_t>(k)] *
                               (cpul_seg[static_cast<size_t>(k)][static_cast<size_t>(j)] / dt);
                }
                // Fortran: matmul(extended_du, cpul/dt) + 0.5*gpul (matrix add, not inside matmul).
                lhs[static_cast<size_t>(i)][static_cast<size_t>(j)] =
                    sum_lhs + 0.5 * gpul_seg[static_cast<size_t>(i)][static_cast<size_t>(j)];
                rhs[static_cast<size_t>(i)][static_cast<size_t>(j)] =
                    sum_rhs - 0.5 * gpul_seg[static_cast<size_t>(i)][static_cast<size_t>(j)];
            }
        }
        i_diff[static_cast<size_t>(seg)] = utils_m::inv(lhs);
        v_term[static_cast<size_t>(seg)] = utils_m::matmul(i_diff[static_cast<size_t>(seg)], rhs);
    }
}

void mtl_bundle_t::updateGenerators(double time_in, double dt_in) {
    for (auto& gen : generators) {
        if (!gen.in_layer) {
            continue;
        }
        const int c = gen.conductor - 1;
        const int idx = gen.index - 1;
        const double du_cc = du[static_cast<size_t>(idx)][static_cast<size_t>(c)][static_cast<size_t>(c)];
        const double val = 0.5 * (gen.interpolate(time_in + dt_in) + gen.interpolate(time_in));
        if (gen.source_type == mtln_types_m::SOURCE_TYPE_VOLTAGE) {
            v_source[static_cast<size_t>(c)][static_cast<size_t>(idx)] = val / du_cc;
        } else if (gen.source_type == mtln_types_m::SOURCE_TYPE_CURRENT) {
            v_source[static_cast<size_t>(c)][static_cast<size_t>(idx)] = val * gen.resistance / du_cc;
        }
    }
}

void mtl_bundle_t::advanceVoltage() {
    const int nc = number_of_conductors;
    for (int seg = 1; seg < number_of_divisions; ++seg) {
        std::vector<double> v_col(static_cast<size_t>(nc));
        std::vector<double> di(static_cast<size_t>(nc));
        std::vector<double> du_is(nc, 0.0);
        for (int c = 0; c < nc; ++c) {
            v_col[static_cast<size_t>(c)] = v[static_cast<size_t>(c)][static_cast<size_t>(seg)];
            di[static_cast<size_t>(c)] =
                i[static_cast<size_t>(c)][static_cast<size_t>(seg)] - i[static_cast<size_t>(c)][static_cast<size_t>(seg - 1)];
        }
        for (int c = 0; c < nc; ++c) {
            for (int k = 0; k < nc; ++k) {
                du_is[static_cast<size_t>(c)] +=
                    du[static_cast<size_t>(seg)][static_cast<size_t>(c)][static_cast<size_t>(k)] *
                    i_source[static_cast<size_t>(k)][static_cast<size_t>(seg)];
            }
        }
        const auto v_new = utils_m::matmul_vec(v_term[static_cast<size_t>(seg)], v_col);
        for (int c = 0; c < nc; ++c) {
            double corr = 0.0;
            for (int k = 0; k < nc; ++k) {
                corr += i_diff[static_cast<size_t>(seg)][static_cast<size_t>(c)][static_cast<size_t>(k)] *
                        (di[static_cast<size_t>(k)] + du_is[static_cast<size_t>(k)]);
            }
            v[static_cast<size_t>(c)][static_cast<size_t>(seg)] = v_new[static_cast<size_t>(c)] - corr;
        }
    }
}

void mtl_bundle_t::advanceCurrent() {
    transfer_impedance.updateQ3Phi();
    i_prev = i;
    const int nc = number_of_conductors;
    for (int seg = 0; seg < number_of_divisions; ++seg) {
        std::vector<double> i_col(static_cast<size_t>(nc));
        std::vector<double> dv(static_cast<size_t>(nc));
        std::vector<double> vsrc(static_cast<size_t>(nc));
        for (int c = 0; c < nc; ++c) {
            i_col[static_cast<size_t>(c)] = i[static_cast<size_t>(c)][static_cast<size_t>(seg)];
            dv[static_cast<size_t>(c)] =
                v[static_cast<size_t>(c)][static_cast<size_t>(seg + 1)] - v[static_cast<size_t>(c)][static_cast<size_t>(seg)] -
                e_L[static_cast<size_t>(c)][static_cast<size_t>(seg)] * step_size[static_cast<size_t>(seg)];
            double src = 0.0;
            for (int k = 0; k < nc; ++k) {
                src += du[static_cast<size_t>(seg)][static_cast<size_t>(c)][static_cast<size_t>(k)] *
                       v_source[static_cast<size_t>(k)][static_cast<size_t>(seg)];
            }
            vsrc[static_cast<size_t>(c)] = src;
        }
        std::vector<double> du_q3(nc, 0.0);
        for (int j = 0; j < nc; ++j) {
            for (int k = 0; k < nc; ++k) {
                du_q3[static_cast<size_t>(j)] +=
                    du[static_cast<size_t>(seg)][static_cast<size_t>(j)][static_cast<size_t>(k)] *
                    std::real(transfer_impedance.q3_phi[static_cast<size_t>(seg)][static_cast<size_t>(k)]);
            }
        }
        std::vector<double> q3_term(nc, 0.0);
        for (int c = 0; c < nc; ++c) {
            for (int j = 0; j < nc; ++j) {
                q3_term[static_cast<size_t>(c)] +=
                    v_diff[static_cast<size_t>(seg)][static_cast<size_t>(c)][static_cast<size_t>(j)] *
                    du_q3[static_cast<size_t>(j)];
            }
        }
        for (int c = 0; c < nc; ++c) {
            double term1 = utils_m::matmul_vec(i_term[static_cast<size_t>(seg)], i_col)[static_cast<size_t>(c)];
            double term2 = 0.0;
            for (int k = 0; k < nc; ++k) {
                term2 += v_diff[static_cast<size_t>(seg)][static_cast<size_t>(c)][static_cast<size_t>(k)] *
                         (dv[static_cast<size_t>(k)] - vsrc[static_cast<size_t>(k)]);
            }
            i[static_cast<size_t>(c)][static_cast<size_t>(seg)] =
                term1 - term2 - q3_term[static_cast<size_t>(c)];
        }
    }
    transfer_impedance.updatePhi(i_prev, i);
}

void mtl_bundle_t::setExternalLongitudinalField() {
    if (conductors_in_level.empty()) {
        return;
    }
    const int ncond = conductors_in_level[0];
    for (int c = 0; c < ncond; ++c) {
        for (int seg = 0; seg < number_of_divisions; ++seg) {
            const auto& efs = external_field_segments[static_cast<size_t>(seg)];
            const double sign = static_cast<double>(efs.direction) / std::abs(static_cast<double>(efs.direction));
            e_L[static_cast<size_t>(c)][static_cast<size_t>(seg)] =
                (efs.field != nullptr) ? (*efs.field) * sign : 0.0;
        }
    }
}

void mtl_bundle_t::addProbe(int index, int probe_type, const std::string& name,
                            const std::vector<double>& position,
                            const std::vector<std::vector<int>>& layer_indices) {
    probes.push_back(probes_m::probeCtor(index, probe_type, dt, position, name, layer_indices));
}

void mtl_bundle_t::addGenerator(int index, int conductor, int gen_type, double resistance,
                                const std::string& path,
                                const std::vector<std::vector<int>>& layer_indices) {
    generators.push_back(
        generators_m::generatorCtor(index, conductor, gen_type, resistance, path, layer_indices));
    if (index < 1 || conductor < 1) {
        return;
    }
    const size_t seg = static_cast<size_t>(index - 1);
    const size_t cond = static_cast<size_t>(conductor - 1);
    if (seg < du.size() && cond < du[seg].size() && cond < du[seg][cond].size()) {
        const double du_cc = du[seg][cond][cond];
        if (du_cc != 0.0) {
            if (gen_type == mtln_types_m::SOURCE_TYPE_VOLTAGE) {
                rpul[seg][cond][cond] += resistance / du_cc;
            } else if (gen_type == mtln_types_m::SOURCE_TYPE_CURRENT) {
                rpul[seg][cond][cond] += resistance / du_cc;
            }
        }
    }
}

void mtl_bundle_t::setConnectorTransferImpedance(int, int, const std::vector<int>&,
                                               const mtln_types_m::transfer_impedance_per_meter_t&) {
    // Connector transfer impedance at bundle level is handled during preprocess setup.
}

} // namespace mtl_bundle_m
