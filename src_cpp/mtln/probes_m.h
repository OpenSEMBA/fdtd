#ifndef PROBES_M_H
#define PROBES_M_H

#include <cstdint>
#include <string>
#include <vector>

#include "mtln_types.h"

namespace probes_m {

using mtln_types_m::PROBE_TYPE_CURRENT;
using mtln_types_m::PROBE_TYPE_VOLTAGE;
using mtln_types_m::RKIND;
using mtln_types_m::RKIND_TIEMPO;

struct probe_t {
    int type = 0;
    std::vector<RKIND> t;
    std::vector<std::vector<RKIND>> val;
    RKIND_TIEMPO dt = 0.0;
    int index = 0;
    int current_frame = 0;
    int unit = 0;
    std::string name;
    bool in_layer = true;
    std::string output_path;

    void resizeFrames(int num_frames, int number_of_conductors);
    void update(RKIND_TIEMPO time,
                const std::vector<std::vector<RKIND>>& v,
                const std::vector<std::vector<RKIND>>& i);
    void saveFrame(RKIND_TIEMPO time, const std::vector<RKIND>& values);
};

probe_t probeCtor(int index,
                  int probe_type,
                  RKIND_TIEMPO dt,
                  const std::vector<RKIND>& position,
                  const std::string& name,
                  const std::vector<std::vector<int32_t>>& layer_indices = {});

} // namespace probes_m

#endif
