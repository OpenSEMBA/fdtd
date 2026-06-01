#ifndef RATIONAL_APPROXIMATION_M_H
#define RATIONAL_APPROXIMATION_M_H

#include <complex>
#include <cmath>
#include <vector>

#include "mtln_types.h"

namespace rational_approximation_m {

using mtln_types_m::RKIND;
using mtln_types_m::RKIND_TIEMPO;
using mtln_types_m::transfer_impedance_per_meter_t;

struct pol_res_t {
    std::vector<std::complex<RKIND>> q1;
    std::vector<std::complex<RKIND>> q2;
    std::vector<std::complex<RKIND>> q3;
    RKIND r = 0.0;
    RKIND l = 0.0;
    int number_of_poles = 0;
    int direction = 0;
};

inline pol_res_t pol_resCtor(const transfer_impedance_per_meter_t& model, RKIND_TIEMPO dt) {
    pol_res_t res;
    res.r = model.resistive_term;
    res.l = model.inductive_term;
    res.number_of_poles = static_cast<int>(model.poles.size());
    res.q1.resize(static_cast<size_t>(res.number_of_poles));
    res.q2.resize(static_cast<size_t>(res.number_of_poles));
    res.q3.resize(static_cast<size_t>(res.number_of_poles));

    if (res.number_of_poles != 0) {
        for (int i = 0; i < res.number_of_poles; ++i) {
            const std::complex<RKIND> alpha = model.residues[static_cast<size_t>(i)] / model.poles[static_cast<size_t>(i)];
            const std::complex<RKIND> beta = model.poles[static_cast<size_t>(i)] * dt;
            res.q1[static_cast<size_t>(i)] =
                -(alpha / beta) * (std::exp(beta) - beta - RKIND(1));
            res.q2[static_cast<size_t>(i)] =
                -(alpha / beta) * (RKIND(1) + std::exp(beta) * (beta - RKIND(1)));
            res.q3[static_cast<size_t>(i)] = -std::exp(beta);
        }
    }
    res.direction = model.direction;
    return res;
}

} // namespace rational_approximation_m

#endif
