#ifndef DISPERSIVE_M_H
#define DISPERSIVE_M_H

#include <complex>
#include <cstddef>
#include <vector>

#include "mtln_types.h"
#include "rational_approximation_m.h"

namespace dispersive_m {

using mtln_types_m::RKIND;
using mtln_types_m::RKIND_TIEMPO;
using mtln_types_m::transfer_impedance_per_meter_t;

using Complex = std::complex<RKIND>;

struct dispersive_t {
    RKIND_TIEMPO dt = 0.0;
    int number_of_divisions = 0;
    int number_of_conductors = 0;
    int number_of_poles = 0;

    std::vector<std::vector<std::vector<Complex>>> phi;
    std::vector<std::vector<std::vector<std::vector<Complex>>>> q1;
    std::vector<std::vector<std::vector<std::vector<Complex>>>> q2;
    std::vector<std::vector<std::vector<std::vector<Complex>>>> q3;
    std::vector<std::vector<std::vector<RKIND>>> d;
    std::vector<std::vector<std::vector<RKIND>>> e;
    std::vector<std::vector<std::vector<Complex>>> q1_sum;
    std::vector<std::vector<std::vector<Complex>>> q2_sum;
    std::vector<std::vector<Complex>> q3_phi;

    dispersive_t() = default;
    dispersive_t(int number_of_conductors, int number_of_poles, int number_of_divisions, RKIND_TIEMPO dt);

    void increaseOrder(int number_of_poles_new);
    void updateQ3Phi();
};

struct lumped_t : dispersive_t {
    lumped_t() = default;
    lumped_t(int number_of_conductors, int number_of_poles, int number_of_divisions, RKIND_TIEMPO dt);

    bool positionIsEmpty(int index, int conductor) const;
    void addDispersiveLumped(int index, int conductor, const transfer_impedance_per_meter_t& model);

private:
    void addDispersiveLumpedInConductor(int index, int conductor,
                                        const rational_approximation_m::pol_res_t& connector);
};

struct transfer_impedance_t : dispersive_t {
    transfer_impedance_t() = default;
    transfer_impedance_t(int number_of_conductors, int number_of_poles, int number_of_divisions,
                         RKIND_TIEMPO dt);

    void addTransferImpedance(int conductor_out, const std::vector<int>& range_in,
                              const transfer_impedance_per_meter_t& model);
    void setTransferImpedance(int index, int conductor_out, const std::vector<int>& range_in,
                              const transfer_impedance_per_meter_t& model);

private:
    void addTransferImpedanceInConductors(int conductor_1, int conductor_2,
                                          const rational_approximation_m::pol_res_t& connector);
    void setTransferImpedanceInConductors(int index, int conductor_1, int conductor_2,
                                         const rational_approximation_m::pol_res_t& connector);
};

std::size_t flatSize4(const std::vector<std::vector<std::vector<std::vector<Complex>>>>& a);
std::size_t flatSize3(const std::vector<std::vector<std::vector<Complex>>>& a);
std::size_t flatSize3Real(const std::vector<std::vector<std::vector<RKIND>>>& a);
std::size_t flatSize2(const std::vector<std::vector<Complex>>& a);

} // namespace dispersive_m

#endif
