#ifndef GENERATORS_M_H
#define GENERATORS_M_H

#include <string>
#include <vector>

#include "mtln_types.h"

namespace generators_m {

using mtln_types_m::RKIND;
using mtln_types_m::RKIND_TIEMPO;
using mtln_types_m::SOURCE_TYPE_CURRENT;
using mtln_types_m::SOURCE_TYPE_UNDEFINED;
using mtln_types_m::SOURCE_TYPE_VOLTAGE;

struct generator_t {
    int index = 0;
    int conductor = 0;
    std::vector<RKIND> value;
    std::vector<RKIND_TIEMPO> time;
    RKIND resistance = 0.0;
    int source_type = SOURCE_TYPE_UNDEFINED;
    bool in_layer = true;

    void initGenerator(const std::string& path);
    RKIND interpolate(RKIND_TIEMPO t) const;
};

generator_t generatorCtor(int index,
                          int conductor,
                          int gen_type,
                          RKIND resistance,
                          const std::string& path,
                          const std::vector<std::vector<int>>& layer_indices = {});

} // namespace generators_m

#endif
