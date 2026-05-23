#ifndef OBSERVATION_PREPROCESS_H
#define OBSERVATION_PREPROCESS_H

#include <vector>
#include <cstdint>
#include "observation_types.h"

namespace Observa_m {

struct probe_t {
    int32_t what = 0;
};

struct Obses_t_full : Obses_t {
    bool done = false;
    bool flushed = false;
    bool begun = false;
    bool Saveall = false;
    std::vector<probe_t> P;
};

void preprocess_observation_full(Obses_t_full& observation, output_t& privateOutput,
                                 const std::vector<RKIND_tiempo>& time,
                                 int finaltimestep, RKIND_tiempo dt, bool saveall);

} // namespace Observa_m

#endif
