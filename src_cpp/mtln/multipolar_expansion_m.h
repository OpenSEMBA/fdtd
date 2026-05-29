#ifndef MULTIPOLAR_EXPANSION_M_H
#define MULTIPOLAR_EXPANSION_M_H

#include <vector>

#include "mtln_types.h"

namespace multipolar_expansion_m {

using mtln_types_m::RKIND;
using mtln_types_m::box_2d_t;
using mtln_types_m::field_reconstruction_t;
using mtln_types_m::multipolar_coefficient_t;
using mtln_types_m::multipolar_expansion_t;

RKIND multipolarExpansion2D(const std::vector<RKIND>& position,
                            const std::vector<multipolar_coefficient_t>& ab,
                            const std::vector<RKIND>& expansionCenter);

RKIND getAveragePotential(const field_reconstruction_t& potential,
                          const box_2d_t& innerBox,
                          const box_2d_t& outerBox);

std::vector<std::vector<RKIND>> getCellCapacitanceOnBox(
    const multipolar_expansion_t& multipolarExpansionParameters,
    const box_2d_t& cellBox);

std::vector<std::vector<RKIND>> getCellInductanceOnBox(
    const multipolar_expansion_t& multipolarExpansionParameters,
    const box_2d_t& cellBox);

} // namespace multipolar_expansion_m

#endif
