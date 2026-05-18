#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdint>

// Assuming FDETYPES_m provides these types
// In a real translation, these would come from an include file
#ifndef RKIND
#define RKIND double
#endif

#ifndef RKIND_TIEMPO
#define RKIND_TIEMPO double
#endif

#ifndef PROBE_TYPE_CURRENT
#define PROBE_TYPE_CURRENT 1
#endif

#ifndef PROBE_TYPE_VOLTAGE
#define PROBE_TYPE_VOLTAGE 2
#endif

// MPI stubs if not compiled with MPI
#ifdef CompileWithMPI
#include <mpi.h>
extern MPI_Comm SUBCOMM_MPI;
#else
// Dummy MPI functions for non-MPI build
#define MPI_COMM_SIZE(comm, size, ierr) do { *size = 1; *ierr = 0; } while(0)
#endif

namespace probes_m {

    struct probe_t {
        int type;
        std::vector<RKIND> t;
        std::vector<std::vector<RKIND>> val;
        RKIND_TIEMPO dt;
        int index;
        int current_frame;
        int unit;
        std::string name;
        bool in_layer = true;

        void resizeFrames(int num_frames, int number_of_conductors);
        void update(RKIND_TIEMPO t, const std::vector<std::vector<RKIND>>& v, const std::vector<std::vector<RKIND>>& i);
        void saveFrame(RKIND_TIEMPO time, const std::vector<RKIND>& values);
    };

    inline probe_t probeCtor(int index, int probe_type, RKIND_TIEMPO dt, 
                             const std::vector<RKIND>& position, 
                             const std::string& name, 
                             const std::vector<std::vector<int32_t>>& layer_indices = {}) {
        probe_t res;
        res.type = probe_type;
        res.index = index;
        res.dt = dt;
        res.current_frame = 1;

#ifdef CompileWithMPI
        if (!layer_indices.empty()) {
            int sizeof_comm = 1;
            int ierr = 0;
            MPI_COMM_SIZE(SUBCOMM_MPI, &sizeof_comm, &ierr);
            
            if (sizeof_comm > 1) {
                res.in_layer = false;
                int slice = 0;
                for (int i = 0; i < static_cast<int>(layer_indices.size()); ++i) {
                    if (index >= layer_indices[i][0] && index <= layer_indices[i][1] + 1) {
                        res.in_layer = true;
                        slice = i + 1; // 1-based index for logic below
                    }
                }

                int layer_index = 0;
                if (res.in_layer) {
                    for (int i = 0; i < slice - 1; ++i) {
                        layer_index = layer_index + layer_indices[i][1] + 1 - (layer_indices[i][0] - 1);
                    }
                    // Note: In Fortran, 'i' in the last line refers to 'slice' from the previous loop context
                    // which is the index where the condition was met. 
                    // However, the loop variable 'i' goes out of scope or retains last value.
                    // In C++, we need to track the specific index.
                    // Let's refactor slightly to match logic exactly.
                    
                    // Re-evaluating the Fortran logic:
                    // do i = 1, size(layer_indices,1)
                    //    if (...) then res%in_layer=.true.; slice = i; end if
                    // end do
                    // ...
                    // do i = 1, slice - 1 ...
                    // ...
                    // layer_index = layer_index + res%index - layer_indices(i,1) + 1
                    // Here 'i' is the value of 'slice' from the previous loop.
                    
                    // To be safe and clear:
                    int found_slice_idx = -1;
                    for (int k = 0; k < static_cast<int>(layer_indices.size()); ++k) {
                        if (index >= layer_indices[k][0] && index <= layer_indices[k][1] + 1) {
                            found_slice_idx = k;
                            break;
                        }
                    }
                    
                    if (found_slice_idx != -1) {
                        for (int k = 0; k < found_slice_idx; ++k) {
                            layer_index = layer_index + layer_indices[k][1] + 1 - (layer_indices[k][0] - 1);
                        }
                        layer_index = layer_index + index - layer_indices[found_slice_idx][0] + 1;
                    }
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

        // Format position
        std::string a = std::to_string(static_cast<int>(position[0]));
        std::string b = std::to_string(static_cast<int>(position[1]));
        std::string c = std::to_string(static_cast<int>(position[2]));
        
        // trim and adjustl equivalent (std::to_string doesn't have leading spaces usually, but just in case)
        auto trim = [](std::string& s) {
            s.erase(0, s.find_first_not_of(" \t"));
            s.erase(s.find_last_not_of(" \t") + 1);
        };
        trim(a);
        trim(b);
        trim(c);

        res.name += "_" + a + "_" + b + "_" + c;

        return res;
    }

    inline void probe_t::resizeFrames(int num_frames, int number_of_conductors) {
        // A single frame suffices: data is written to disk each step right after
        // it is stored, so there is no need to buffer the full time history.
        t.resize(1);
        val.resize(1, std::vector<RKIND>(number_of_conductors));
        t[0] = 0.0;
        for (auto& row : val) {
            for (auto& val_elem : row) {
                val_elem = 0.0;
            }
        }
    }

    inline void probe_t::update(RKIND_TIEMPO t, const std::vector<std::vector<RKIND>>& v, const std::vector<std::vector<RKIND>>& i) {
        if (type == PROBE_TYPE_VOLTAGE) {
            // v(:,this%index) -> v[index-1] if 1-based, but vector is 0-based.
            // Fortran: v(:, index). If index is 1-based, it's v[index-1] in 0-based C++.
            // However, the vector passed is already sliced.
            // In Fortran: v(:, this%index) means all rows, column this%index.
            // Assuming v is (conductors, probes) or similar.
            // Let's assume the signature matches the intent: pass the column vector.
            // We need to extract the column.
            std::vector<RKIND> col_v;
            if (!v.empty() && this->index > 0 && this->index <= static_cast<int>(v[0].size())) {
                for (const auto& row : v) {
                    col_v.push_back(row[this->index - 1]);
                }
            }
            saveFrame(t, col_v);
        } else if (type == PROBE_TYPE_CURRENT) {
            if (this->index == static_cast<int>(i[0].size()) + 1) {
                // i(:, this%index - 1)
                std::vector<RKIND> col_i;
                if (!i.empty() && this->index - 1 > 0 && this->index - 1 <= static_cast<int>(i[0].size())) {
                    for (const auto& row : i) {
                        col_i.push_back(row[this->index - 2]);
                    }
                }
                saveFrame(t + 0.5 * this->dt, col_i);
            } else {
                // i(:, this%index)
                std::vector<RKIND> col_i;
                if (!i.empty() && this->index > 0 && this->index <= static_cast<int>(i[0].size())) {
                    for (const auto& row : i) {
                        col_i.push_back(row[this->index - 1]);
                    }
                }
                saveFrame(t + 0.5 * this->dt, col_i);
            }
        }
    }

    inline void probe_t::saveFrame(RKIND_TIEMPO time, const std::vector<RKIND>& values) {
        // Always overwrite slot 1; the caller flushes to disk each step.
        if (t.size() < 1) t.resize(1);
        if (val.size() < 1) val.resize(1);
        
        t[0] = time;
        if (val[0].size() != values.size()) {
            val[0].resize(values.size());
        }
        for (size_t k = 0; k < values.size(); ++k) {
            val[0][k] = values[k];
        }
    }

} // namespace probes_m