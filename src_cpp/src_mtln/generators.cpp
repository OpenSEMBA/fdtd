#include "generators_m.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

#ifdef CompileWithMPI
#include <mpi.h>
extern MPI_Comm SUBCOMM_MPI;
#endif

namespace generators_m {

void generator_t::initGenerator(const std::string& path) {
    if (path.empty()) {
        time.clear();
        value.clear();
        return;
    }
    int line_count = 0;
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << path << std::endl;
        return;
    }
    RKIND_TIEMPO t;
    RKIND val;
    while (file >> t >> val) {
        ++line_count;
    }
    file.close();

    time.resize(static_cast<size_t>(line_count));
    value.resize(static_cast<size_t>(line_count));
    file.open(path);
    for (int i = 0; i < line_count; ++i) {
        if (!(file >> time[static_cast<size_t>(i)] >> value[static_cast<size_t>(i)])) {
            throw std::runtime_error("Error reading excitation file");
        }
    }
}

RKIND generator_t::interpolate(RKIND_TIEMPO t) const {
    if (time.empty()) {
        return 0.0;
    }
    int index = 0;
    RKIND_TIEMPO max_timediff = -std::numeric_limits<RKIND_TIEMPO>::infinity();
    for (int i = 0; i < static_cast<int>(time.size()); ++i) {
        const RKIND_TIEMPO timediff = time[static_cast<size_t>(i)] - t;
        if (timediff <= 0.0 && timediff >= max_timediff) {
            max_timediff = timediff;
            index = i + 1;
        }
    }
    if (index == 0) {
        index = 1;
    }
    if (index >= static_cast<int>(time.size())) {
        index = static_cast<int>(time.size()) - 1;
    }
    const int idx0 = index - 1;
    const RKIND x1 = time[static_cast<size_t>(idx0)];
    const RKIND y1 = value[static_cast<size_t>(idx0)];
    RKIND x2 = x1;
    RKIND y2 = y1;
    if (idx0 + 1 < static_cast<int>(time.size())) {
        x2 = time[static_cast<size_t>(idx0 + 1)];
        y2 = value[static_cast<size_t>(idx0 + 1)];
    }
    if (x2 == x1) {
        return y1;
    }
    return (t * (y2 - y1) + x2 * y1 - x1 * y2) / (x2 - x1);
}

generator_t generatorCtor(int index, int conductor, int gen_type, RKIND resistance_in,
                          const std::string& path,
                          const std::vector<std::vector<int>>& layer_indices_opt) {
    generator_t res;
    res.index = index;
    res.conductor = conductor;
    res.resistance = resistance_in;
    res.source_type = gen_type;

#ifdef CompileWithMPI
    if (!layer_indices_opt.empty()) {
        int sizeof_comm = 0;
        MPI_Comm_size(SUBCOMM_MPI, &sizeof_comm);
        if (sizeof_comm > 1) {
            res.in_layer = false;
            int slice = 0;
            for (int i = 0; i < static_cast<int>(layer_indices_opt.size()); ++i) {
                const int lower = layer_indices_opt[static_cast<size_t>(i)][0];
                const int upper = layer_indices_opt[static_cast<size_t>(i)][1];
                if (index >= lower && index <= upper + 1) {
                    res.in_layer = true;
                    slice = i + 1;
                }
            }
            int layer_index = 0;
            if (res.in_layer) {
                for (int i = 0; i < slice - 1; ++i) {
                    const int lower = layer_indices_opt[static_cast<size_t>(i)][0];
                    const int upper = layer_indices_opt[static_cast<size_t>(i)][1];
                    layer_index += upper + 1 - (lower - 1);
                }
                const int last_slice_idx = slice - 1;
                layer_index += res.index - layer_indices_opt[static_cast<size_t>(last_slice_idx)][0] + 1;
            }
            res.index = layer_index;
        }
    }
#endif

    if (res.in_layer) {
        res.initGenerator(path);
    }
    return res;
}

} // namespace generators_m
