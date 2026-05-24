#include <stdexcept>

#include "dispersive_m.h"
#include "rational_approximation_m.h"

namespace dispersive_m {

namespace {

using mtln_types_m::TRANSFER_IMPEDANCE_DIRECTION_BOTH;
using mtln_types_m::TRANSFER_IMPEDANCE_DIRECTION_INWARDS;
using mtln_types_m::TRANSFER_IMPEDANCE_DIRECTION_OUTWARDS;

Complex zeroComplex() {
    return {0.0, 0.0};
}

bool isCouplingInwards(int direction) {
    return direction == TRANSFER_IMPEDANCE_DIRECTION_INWARDS ||
           direction == TRANSFER_IMPEDANCE_DIRECTION_BOTH;
}

bool isCouplingOutwards(int direction) {
    return direction == TRANSFER_IMPEDANCE_DIRECTION_OUTWARDS ||
           direction == TRANSFER_IMPEDANCE_DIRECTION_BOTH;
}

std::vector<std::vector<std::vector<Complex>>> sumQComponents(
    const std::vector<std::vector<std::vector<std::vector<Complex>>>>& a) {
    if (a.empty() || a[0].empty() || a[0][0].empty() || a[0][0][0].empty()) {
        return {};
    }
    const int nd = static_cast<int>(a.size());
    const int nc = static_cast<int>(a[0].size());
    std::vector<std::vector<std::vector<Complex>>> res(
        static_cast<size_t>(nd),
        std::vector<std::vector<Complex>>(static_cast<size_t>(nc), std::vector<Complex>(static_cast<size_t>(nc), zeroComplex())));
    for (int i = 0; i < nd; ++i) {
        for (int j = 0; j < nc; ++j) {
            for (int k = 0; k < nc; ++k) {
                for (int p = 0; p < static_cast<int>(a[i][j][k].size()); ++p) {
                    res[static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(k)] +=
                        a[static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(k)][static_cast<size_t>(p)];
                }
            }
        }
    }
    return res;
}

void allocateDispersive(dispersive_t& res, int number_of_conductors, int number_of_poles,
                        int number_of_divisions) {
    res.number_of_conductors = number_of_conductors;
    res.number_of_poles = number_of_poles;
    res.number_of_divisions = number_of_divisions;

    const Complex z = zeroComplex();
    res.phi.assign(static_cast<size_t>(number_of_divisions),
                   std::vector<std::vector<Complex>>(
                       static_cast<size_t>(number_of_conductors),
                       std::vector<Complex>(static_cast<size_t>(number_of_poles), z)));
    auto make4 = [&](int np) {
        return std::vector<std::vector<std::vector<std::vector<Complex>>>>(
            static_cast<size_t>(number_of_divisions),
            std::vector<std::vector<std::vector<Complex>>>(
                static_cast<size_t>(number_of_conductors),
                std::vector<std::vector<Complex>>(static_cast<size_t>(number_of_conductors),
                                                  std::vector<Complex>(static_cast<size_t>(np), z))));
    };
    res.q1 = make4(number_of_poles);
    res.q2 = make4(number_of_poles);
    res.q3 = make4(number_of_poles);
    res.d.assign(static_cast<size_t>(number_of_divisions),
                 std::vector<std::vector<RKIND>>(static_cast<size_t>(number_of_conductors),
                                                 std::vector<RKIND>(static_cast<size_t>(number_of_conductors), 0.0)));
    res.e = res.d;
    res.q1_sum.assign(static_cast<size_t>(number_of_divisions),
                      std::vector<std::vector<Complex>>(static_cast<size_t>(number_of_conductors),
                                                        std::vector<Complex>(static_cast<size_t>(number_of_conductors), z)));
    res.q2_sum = res.q1_sum;
    res.q3_phi.assign(static_cast<size_t>(number_of_divisions),
                      std::vector<Complex>(static_cast<size_t>(number_of_conductors), z));
}

} // namespace

dispersive_t::dispersive_t(int number_of_conductors, int number_of_poles, int number_of_divisions,
                           RKIND_TIEMPO dt_in) {
    dt = dt_in;
    allocateDispersive(*this, number_of_conductors, number_of_poles, number_of_divisions);
}

void dispersive_t::increaseOrder(int number_of_poles_new) {
    dispersive_t newer(number_of_conductors, number_of_poles_new, number_of_divisions, dt);
    for (int i = 0; i < number_of_divisions; ++i) {
        for (int j = 0; j < number_of_conductors; ++j) {
            for (int k = 0; k < number_of_conductors; ++k) {
                for (int p = 0; p < number_of_poles; ++p) {
                    newer.q1[static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(k)][static_cast<size_t>(p)] =
                        q1[static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(k)][static_cast<size_t>(p)];
                    newer.q2[static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(k)][static_cast<size_t>(p)] =
                        q2[static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(k)][static_cast<size_t>(p)];
                    newer.q3[static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(k)][static_cast<size_t>(p)] =
                        q3[static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(k)][static_cast<size_t>(p)];
                }
            }
            for (int p = 0; p < number_of_poles; ++p) {
                newer.phi[static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(p)] =
                    phi[static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(p)];
            }
        }
    }
    q1 = std::move(newer.q1);
    q2 = std::move(newer.q2);
    q3 = std::move(newer.q3);
    phi = std::move(newer.phi);
    number_of_poles = number_of_poles_new;
}

void dispersive_t::updatePhi(const std::vector<std::vector<RKIND>>& i_prev,
                             const std::vector<std::vector<RKIND>>& i_now) {
    for (int k = 0; k < number_of_poles; ++k) {
        for (int i_div = 0; i_div < number_of_divisions; ++i_div) {
            for (int i = 0; i < number_of_conductors; ++i) {
                Complex sum = zeroComplex();
                for (int j = 0; j < number_of_conductors; ++j) {
                    Complex term1 = zeroComplex();
                    Complex term2 = zeroComplex();
                    Complex term3 = zeroComplex();
                    for (int p = 0; p < number_of_poles; ++p) {
                        if (p == k) {
                            term3 += q3[static_cast<size_t>(i_div)][static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(p)] *
                                     phi[static_cast<size_t>(i_div)][static_cast<size_t>(j)][static_cast<size_t>(p)];
                        }
                    }
                    term1 += q1[static_cast<size_t>(i_div)][static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(k)] *
                             Complex(i_now[static_cast<size_t>(j)][static_cast<size_t>(i_div)], 0.0);
                    term2 += q2[static_cast<size_t>(i_div)][static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(k)] *
                             Complex(i_prev[static_cast<size_t>(j)][static_cast<size_t>(i_div)], 0.0);
                    sum += term1 + term2 + term3;
                }
                phi[static_cast<size_t>(i_div)][static_cast<size_t>(i)][static_cast<size_t>(k)] = sum;
            }
        }
    }
}

void dispersive_t::updateQ3Phi() {
    for (auto& row : q3_phi) {
        std::fill(row.begin(), row.end(), zeroComplex());
    }
    for (int i_div = 0; i_div < number_of_divisions; ++i_div) {
        for (int i = 0; i < number_of_conductors; ++i) {
            for (int j = 0; j < number_of_conductors; ++j) {
                Complex dot = zeroComplex();
                for (int p = 0; p < number_of_poles; ++p) {
                    dot += q3[static_cast<size_t>(i_div)][static_cast<size_t>(i)][static_cast<size_t>(j)][static_cast<size_t>(p)] *
                           phi[static_cast<size_t>(i_div)][static_cast<size_t>(j)][static_cast<size_t>(p)];
                }
                q3_phi[static_cast<size_t>(i_div)][static_cast<size_t>(i)] += dot;
            }
        }
    }
}

lumped_t::lumped_t(int number_of_conductors, int number_of_poles, int number_of_divisions,
                   RKIND_TIEMPO dt_in)
    : dispersive_t(number_of_conductors, number_of_poles, number_of_divisions, dt_in) {}

bool lumped_t::positionIsEmpty(int index, int conductor) const {
    const int i = index - 1;
    const int c = conductor - 1;
    if (d[static_cast<size_t>(i)][static_cast<size_t>(c)][static_cast<size_t>(c)] != 0.0 ||
        e[static_cast<size_t>(i)][static_cast<size_t>(c)][static_cast<size_t>(c)] != 0.0) {
        return false;
    }
    for (int p = 0; p < number_of_poles; ++p) {
        if (q1[static_cast<size_t>(i)][static_cast<size_t>(c)][static_cast<size_t>(c)][static_cast<size_t>(p)] != zeroComplex() ||
            q2[static_cast<size_t>(i)][static_cast<size_t>(c)][static_cast<size_t>(c)][static_cast<size_t>(p)] != zeroComplex() ||
            q3[static_cast<size_t>(i)][static_cast<size_t>(c)][static_cast<size_t>(c)][static_cast<size_t>(p)] != zeroComplex()) {
            return false;
        }
    }
    return true;
}

void lumped_t::addDispersiveLumpedInConductor(int index, int conductor,
                                              const rational_approximation_m::pol_res_t& connector) {
    const int i = index - 1;
    const int c = conductor - 1;
    d[static_cast<size_t>(i)][static_cast<size_t>(c)][static_cast<size_t>(c)] += connector.r;
    e[static_cast<size_t>(i)][static_cast<size_t>(c)][static_cast<size_t>(c)] += connector.l;
    if (connector.number_of_poles != 0) {
        for (int p = 0; p < connector.number_of_poles; ++p) {
            q1[static_cast<size_t>(i)][static_cast<size_t>(c)][static_cast<size_t>(c)][static_cast<size_t>(p)] -= connector.q1[static_cast<size_t>(p)];
            q2[static_cast<size_t>(i)][static_cast<size_t>(c)][static_cast<size_t>(c)][static_cast<size_t>(p)] -= connector.q2[static_cast<size_t>(p)];
            q3[static_cast<size_t>(i)][static_cast<size_t>(c)][static_cast<size_t>(c)][static_cast<size_t>(p)] -= connector.q3[static_cast<size_t>(p)];
        }
    }
}

void lumped_t::addDispersiveLumped(int index, int conductor, const transfer_impedance_per_meter_t& model) {
    if (!positionIsEmpty(index, conductor)) {
        throw std::runtime_error("Dispersive connector already in conductor at position");
    }
    auto connector = rational_approximation_m::pol_resCtor(model, dt);
    if (connector.number_of_poles > number_of_poles) {
        increaseOrder(connector.number_of_poles);
    }
    addDispersiveLumpedInConductor(index, conductor, connector);
    q1_sum = sumQComponents(q1);
    q2_sum = sumQComponents(q2);
}

transfer_impedance_t::transfer_impedance_t(int number_of_conductors, int number_of_poles,
                                           int number_of_divisions, RKIND_TIEMPO dt_in)
    : dispersive_t(number_of_conductors, number_of_poles, number_of_divisions, dt_in) {}

void transfer_impedance_t::addTransferImpedanceInConductors(
    int conductor_1, int conductor_2, const rational_approximation_m::pol_res_t& connector) {
    const int c1 = conductor_1 - 1;
    const int c2 = conductor_2 - 1;
    for (int i = 0; i < number_of_divisions; ++i) {
        d[static_cast<size_t>(i)][static_cast<size_t>(c1)][static_cast<size_t>(c2)] -= connector.r;
        e[static_cast<size_t>(i)][static_cast<size_t>(c1)][static_cast<size_t>(c2)] -= connector.l;
        if (connector.number_of_poles != 0) {
            for (int p = 0; p < connector.number_of_poles; ++p) {
                q1[static_cast<size_t>(i)][static_cast<size_t>(c1)][static_cast<size_t>(c2)][static_cast<size_t>(p)] +=
                    connector.q1[static_cast<size_t>(p)];
                q2[static_cast<size_t>(i)][static_cast<size_t>(c1)][static_cast<size_t>(c2)][static_cast<size_t>(p)] +=
                    connector.q2[static_cast<size_t>(p)];
                q3[static_cast<size_t>(i)][static_cast<size_t>(c1)][static_cast<size_t>(c2)][static_cast<size_t>(p)] +=
                    connector.q3[static_cast<size_t>(p)];
            }
        }
    }
}

void transfer_impedance_t::setTransferImpedanceInConductors(
    int index, int conductor_1, int conductor_2, const rational_approximation_m::pol_res_t& connector) {
    const int i = index - 1;
    const int c1 = conductor_1 - 1;
    const int c2 = conductor_2 - 1;
    d[static_cast<size_t>(i)][static_cast<size_t>(c1)][static_cast<size_t>(c2)] = -connector.r;
    e[static_cast<size_t>(i)][static_cast<size_t>(c1)][static_cast<size_t>(c2)] = -connector.l;
    if (connector.number_of_poles != 0) {
        for (int p = 0; p < connector.number_of_poles; ++p) {
            q1[static_cast<size_t>(i)][static_cast<size_t>(c1)][static_cast<size_t>(c2)][static_cast<size_t>(p)] =
                connector.q1[static_cast<size_t>(p)];
            q2[static_cast<size_t>(i)][static_cast<size_t>(c1)][static_cast<size_t>(c2)][static_cast<size_t>(p)] =
                connector.q2[static_cast<size_t>(p)];
            q3[static_cast<size_t>(i)][static_cast<size_t>(c1)][static_cast<size_t>(c2)][static_cast<size_t>(p)] =
                connector.q3[static_cast<size_t>(p)];
        }
    }
}

void transfer_impedance_t::addTransferImpedance(int conductor_out, const std::vector<int>& range_in,
                                                const transfer_impedance_per_meter_t& model) {
    auto connector = rational_approximation_m::pol_resCtor(model, dt);
    if (connector.number_of_poles > number_of_poles) {
        increaseOrder(connector.number_of_poles);
    }
    for (int in : range_in) {
        if (isCouplingInwards(connector.direction)) {
            addTransferImpedanceInConductors(in, conductor_out, connector);
        }
        if (isCouplingOutwards(connector.direction)) {
            addTransferImpedanceInConductors(conductor_out, in, connector);
        }
    }
    q1_sum = sumQComponents(q1);
    q2_sum = sumQComponents(q2);
}

void transfer_impedance_t::setTransferImpedance(int index, int conductor_out,
                                                const std::vector<int>& range_in,
                                                const transfer_impedance_per_meter_t& model) {
    auto connector = rational_approximation_m::pol_resCtor(model, dt);
    if (connector.number_of_poles > number_of_poles) {
        increaseOrder(connector.number_of_poles);
    }
    for (int in : range_in) {
        if (isCouplingInwards(connector.direction)) {
            setTransferImpedanceInConductors(index, in, conductor_out, connector);
        }
        if (isCouplingOutwards(connector.direction)) {
            setTransferImpedanceInConductors(index, conductor_out, in, connector);
        }
    }
    q1_sum = sumQComponents(q1);
    q2_sum = sumQComponents(q2);
}

std::size_t flatSize4(const std::vector<std::vector<std::vector<std::vector<Complex>>>>& a) {
    if (a.empty()) return 0;
    return a.size() * a[0].size() * a[0][0].size() * a[0][0][0].size();
}

std::size_t flatSize3(const std::vector<std::vector<std::vector<Complex>>>& a) {
    if (a.empty()) return 0;
    return a.size() * a[0].size() * a[0][0].size();
}

std::size_t flatSize3Real(const std::vector<std::vector<std::vector<RKIND>>>& a) {
    if (a.empty()) return 0;
    return a.size() * a[0].size() * a[0][0].size();
}

std::size_t flatSize2(const std::vector<std::vector<Complex>>& a) {
    if (a.empty()) return 0;
    return a.size() * a[0].size();
}

} // namespace dispersive_m
