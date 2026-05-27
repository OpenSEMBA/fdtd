#include "probes_m.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#ifdef CompileWithMPI
#include <mpi.h>
extern MPI_Comm SUBCOMM_MPI;
#endif

namespace probes_m {

probe_t probeCtor(int index, int probe_type, RKIND_TIEMPO dt_in,
                  const std::vector<RKIND>& position, const std::string& name,
                  const std::vector<std::vector<int32_t>>& layer_indices) {
    probe_t res;
    res.type = probe_type;
    res.index = index;
    res.dt = dt_in;
    res.current_frame = 1;

#ifdef CompileWithMPI
    if (!layer_indices.empty()) {
        int sizeof_comm = 1;
        MPI_Comm_size(SUBCOMM_MPI, &sizeof_comm);
        if (sizeof_comm > 1) {
            res.in_layer = false;
            int found_slice_idx = -1;
            for (int k = 0; k < static_cast<int>(layer_indices.size()); ++k) {
                if (index >= layer_indices[static_cast<size_t>(k)][0] &&
                    index <= layer_indices[static_cast<size_t>(k)][1] + 1) {
                    res.in_layer = true;
                    found_slice_idx = k;
                    break;
                }
            }
            int layer_index = 0;
            if (res.in_layer && found_slice_idx >= 0) {
                for (int k = 0; k < found_slice_idx; ++k) {
                    layer_index += layer_indices[static_cast<size_t>(k)][1] + 1 -
                                   (layer_indices[static_cast<size_t>(k)][0] - 1);
                }
                layer_index += index - layer_indices[static_cast<size_t>(found_slice_idx)][0] + 1;
            }
            res.index = layer_index;
        }
    }
#endif

    res.name = name + "_";
    if (probe_type == PROBE_TYPE_VOLTAGE) {
        res.name += "V";
    } else if (probe_type == PROBE_TYPE_CURRENT) {
        res.name += "I";
    } else {
        throw std::runtime_error("Undefined probe");
    }

    auto trim = [](std::string s) {
        s.erase(0, s.find_first_not_of(" \t"));
        s.erase(s.find_last_not_of(" \t") + 1);
        return s;
    };
    res.name += "_" + trim(std::to_string(static_cast<int>(position[0]))) + "_" +
                trim(std::to_string(static_cast<int>(position[1]))) + "_" +
                trim(std::to_string(static_cast<int>(position[2])));
    return res;
}

void probe_t::resizeFrames(int, int number_of_conductors) {
    t.resize(1);
    val.resize(1, std::vector<RKIND>(static_cast<size_t>(number_of_conductors), 0.0));
    t[0] = 0.0;
}

void probe_t::update(RKIND_TIEMPO time_in,
                     const std::vector<std::vector<RKIND>>& v,
                     const std::vector<std::vector<RKIND>>& i) {
    if (type == PROBE_TYPE_VOLTAGE) {
        std::vector<RKIND> col_v;
        if (!v.empty() && index > 0 && index <= static_cast<int>(v[0].size())) {
            for (const auto& row : v) {
                col_v.push_back(row[static_cast<size_t>(index - 1)]);
            }
        }
        saveFrame(time_in, col_v);
    } else if (type == PROBE_TYPE_CURRENT) {
        if (index == static_cast<int>(i[0].size()) + 1) {
            std::vector<RKIND> col_i;
            if (!i.empty() && index - 1 > 0 && index - 1 <= static_cast<int>(i[0].size())) {
                for (const auto& row : i) {
                    col_i.push_back(row[static_cast<size_t>(index - 2)]);
                }
            }
            saveFrame(time_in + 0.5 * dt, col_i);
        } else {
            std::vector<RKIND> col_i;
            if (!i.empty() && index > 0 && index <= static_cast<int>(i[0].size())) {
                for (const auto& row : i) {
                    col_i.push_back(row[static_cast<size_t>(index - 1)]);
                }
            }
            saveFrame(time_in + 0.5 * dt, col_i);
        }
    }
}

void probe_t::saveFrame(RKIND_TIEMPO time_in, const std::vector<RKIND>& values) {
    if (t.size() < 1) {
        t.resize(1);
    }
    if (val.size() < 1) {
        val.resize(1);
    }
    t[0] = time_in;
    if (val[0].size() != values.size()) {
        val[0].resize(values.size());
    }
    for (size_t k = 0; k < values.size(); ++k) {
        val[0][k] = values[k];
    }
}

} // namespace probes_m
