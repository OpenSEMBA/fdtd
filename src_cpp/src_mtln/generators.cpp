#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <optional>
#include <cstdint>

// Assuming these types are defined in mtln_types_m
enum SourceType {
    SOURCE_TYPE_UNDEFINED = 0,
    SOURCE_TYPE_VOLTAGE = 1,
    SOURCE_TYPE_CURRENT = 2
};

// Assuming RKIND and RKIND_tiempo are double based on typical Fortran real(8) usage
using RKIND = double;
using RKIND_tiempo = double;

// MPI stubs if CompileWithMPI is defined
#ifdef CompileWithMPI
#include <mpi.h>
extern MPI_Comm SUBCOMM_MPI;
#endif

namespace generators_m {

    struct generator_t {
        int index;
        int conductor;
        std::vector<RKIND> value;
        std::vector<RKIND_tiempo> time;
        RKIND resistance;
        SourceType source_type = SOURCE_TYPE_UNDEFINED;
        bool in_layer = true;

        void initGenerator(const std::string& path) {
            RKIND val;
            RKIND_tiempo t;
            int io;
            int line_count = 0;
            int i;

            if (path.empty()) {
                // error
                time.resize(0);
                value.resize(0);
                return;
            }

            // First pass: count the number of lines
            line_count = 0;
            std::ifstream file(path);
            if (!file.is_open()) {
                std::cerr << "Error opening file: " << path << std::endl;
                return;
            }
            while (file >> t >> val) {
                line_count++;
            }
            file.close();

            // Allocate arrays with the exact size needed
            time.resize(line_count);
            value.resize(line_count);

            // Second pass: fill the arrays
            file.open(path);
            if (!file.is_open()) {
                std::cerr << "Error reopening file: " << path << std::endl;
                return;
            }
            for (i = 0; i < line_count; ++i) {
                if (!(file >> time[i] >> value[i])) {
                    file.close();
                    throw std::runtime_error("Error reading excitation file");
                }
            }
            file.close();
        }

        RKIND interpolate(RKIND_tiempo t) const {
            RKIND res;
            RKIND x1, x2;
            RKIND y1, y2;
            int index;
            
            // Calculate time differences
            std::vector<RKIND_tiempo> timediff(time.size());
            for (size_t i = 0; i < time.size(); ++i) {
                timediff[i] = time[i] - t;
            }

            // Find maxloc where timediff <= 0
            // Fortran maxloc returns the index of the first maximum element.
            // We want the largest index i such that time[i] <= t (since timediff <= 0).
            // However, maxloc on (timediff <= 0) finds the first true value if we treat booleans as 1/0?
            // No, maxloc on logical array returns index of first .true.
            // But we want the last time point <= t for interpolation usually? 
            // Let's look at the logic:
            // timediff = this%time - t.
            // If time[i] <= t, timediff[i] <= 0.
            // maxloc(timediff, 1, (timediff) <= 0) finds the index of the first element where timediff <= 0.
            // Since time is likely sorted ascending, the first element <= t is the start of the interval containing t.
            
            index = -1;
            for (int i = 0; i < static_cast<int>(timediff.size()); ++i) {
                if (timediff[i] <= 0.0) {
                    index = i + 1; // Fortran is 1-based
                    break;
                }
            }

            if (index == 0) index = 1;
            
            // Adjust to 0-based for C++ vector access
            int idx0 = index - 1;
            
            x1 = time[idx0];
            y1 = value[idx0];
            
            if (idx0 + 1 >= static_cast<int>(time.size())) {
                x2 = x1;
                y2 = y1;
            } else {
                x2 = time[idx0 + 1];
                y2 = value[idx0 + 1];
            }
            
            if (x2 == x1) {
                res = y1;
            } else {
                res = (t * (y2 - y1) + x2 * y1 - x1 * y2) / (x2 - x1);
            }
            
            return res;
        }
    };

    inline generator_t generatorCtor(int index, int conductor, int gen_type, RKIND resistance, const std::string& path, const std::vector<std::vector<int>>& layer_indices_opt) {
        generator_t res;
        res.index = index;
        res.conductor = conductor;
        res.resistance = resistance;
        res.source_type = static_cast<SourceType>(gen_type);
        
#ifdef CompileWithMPI
        if (!layer_indices_opt.empty()) {
            int sizeof_comm = 0;
            int ierr;
            MPI_Comm_size(SUBCOMM_MPI, &sizeof_comm, &ierr);
            
            if (sizeof_comm > 1) {
                res.in_layer = false;
                int slice = 0;
                int num_layers = static_cast<int>(layer_indices_opt.size());
                for (int i = 0; i < num_layers; ++i) {
                    // layer_indices is 1-based in Fortran, but passed as vector of vectors.
                    // Assuming the input vector is 0-based but contains 1-based indices from Fortran logic?
                    // Or is the input vector already adjusted?
                    // The Fortran code uses layer_indices(i, 1) and layer_indices(i, 2).
                    // Let's assume the passed vector is 0-indexed in C++ but holds 1-based indices.
                    int lower = layer_indices_opt[i][0];
                    int upper = layer_indices_opt[i][1];
                    
                    if (index >= lower && index <= upper + 1) {
                        res.in_layer = true;
                        slice = i + 1; // 1-based slice index
                    }
                }

                int layer_index = 0;
                if (res.in_layer) {
                    for (int i = 0; i < slice - 1; ++i) {
                        int lower = layer_indices_opt[i][0];
                        int upper = layer_indices_opt[i][1];
                        layer_index = layer_index + upper + 1 - (lower - 1);
                    }
                    // Note: 'i' in the do loop above is the loop variable. 
                    // In Fortran, after the loop, 'i' retains its last value? No, usually not guaranteed or scope-dependent.
                    // But here 'slice' is used. The loop runs from 1 to slice-1.
                    // The last 'i' in the loop is slice-1.
                    // So we need to access layer_indices(slice-1, 1).
                    // In C++, slice is 1-based, so index is slice-1.
                    int last_slice_idx = slice - 1;
                    int lower = layer_indices_opt[last_slice_idx][0];
                    layer_index = layer_index + res.index - lower + 1;
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

}