#include <vector>
#include <string>
#include <optional>
#include <iostream>
#include <stdexcept>

// Forward declarations and stubs for external types/modules to make the code compile conceptually
// In a real translation, these would be defined in their respective header files.

namespace FDETYPES_m {
    constexpr int RKIND_TIEMPO = 8; // Assuming double precision
}

namespace mtln_types_m {
    struct probe_t {};
    struct mtln_t {
        double time_step;
        int number_of_steps;
        std::vector<probe_t> probes;
        std::vector<void*> wireGenerators; // Stub for generators
        std::vector<void*> networks; // Stub for networks
        std::vector<void*> cables; // Stub for cables
    };
}

namespace mtl_bundle_m {
    struct cable_abstract_t {};
    struct cable_level_t {
        std::vector<cable_abstract_t> cables;
    };
    struct cable_bundle_t {
        std::vector<cable_level_t> levels;
    };
    struct transmission_line_level_t {
        struct line_t {
            std::string name;
            int number_of_conductors;
            int conductor_in_parent;
            double transfer_impedance;
            std::vector<void*> initial_connector_transfer_impedances; // Stub
            std::vector<void*> end_connector_transfer_impedances; // Stub
        };
        std::vector<line_t> lines;
    };
    struct transmission_line_bundle_t {
        std::vector<transmission_line_level_t> levels;
    };
    struct transfer_impedance_per_meter_t {
        bool has_transfer_impedance() const { return false; }
    };
    struct mtl_t {
        std::string parent_name;
        std::string name;
        int number_of_conductors;
        int conductor_in_parent;
        std::vector<std::vector<double>> du; // Stub for 2D vector access
        std::vector<int> conductors_in_level;
        
        void addTransferImpedance(int outer, const std::vector<int>& inner, double z) {}
        void setConnectorTransferImpedance(int pos1, int pos2, const std::vector<int>& inner, const transfer_impedance_per_meter_t& zt) {}
    };
    struct mtl_bundle_t {
        std::vector<std::vector<double>> du; // Stub
        std::vector<int> conductors_in_level;
        
        void addTransferImpedance(int outer, const std::vector<int>& inner, double z) {}
        void setConnectorTransferImpedance(int pos1, int pos2, const std::vector<int>& inner, const transfer_impedance_per_meter_t& zt) {}
    };
}

namespace network_manager_m {
    struct network_manager_t {};
}

namespace mtl_m {
    struct XYZlimit_t {
        double min, max;
    };
}

namespace Report_m {
    inline void WarnErrReport(const std::string& msg) {
        std::cerr << "Warning/Error: " << msg << std::endl;
    }
}

namespace fhash {
    struct fhash_key_t {
        std::string key;
        bool operator<(const fhash_key_t& other) const { return key < other.key; }
    };
    using key = fhash_key_t;
    
    struct fhash_tbl_t {
        std::vector<fhash_key_t> keys;
        std::vector<int> values;

        void insert(const fhash_key_t& k, int v) {
            keys.push_back(k);
            values.push_back(v);
        }
        
        bool find(const fhash_key_t& k, int& v) const {
            for(size_t i=0; i<keys.size(); ++i) {
                if(keys[i].key == k.key) {
                    v = values[i];
                    return true;
                }
            }
            return false;
        }
    };
}

// MPI Stub
#ifdef CompileWithMPI
#include <mpi.h>
extern MPI_Comm SUBCOMM_MPI;
extern MPI_Comm subcomm_mpi;
#endif

namespace mtln_preprocess_m {

    constexpr int XPOS = 1;
    constexpr int XNEG = -1;
    constexpr int YPOS = 2;
    constexpr int YNEG = -2;
    constexpr int ZPOS = 3;
    constexpr int ZNEG = -3;

    struct preprocess_t {
        std::vector<mtl_bundle_m::mtl_bundle_t> bundles;
        network_manager_m::network_manager_t network_manager;
        fhash::fhash_tbl_t conductors_before_cable;
        fhash::fhash_tbl_t cable_name_to_bundle_id;
        double final_time;
        double dt;

        // Methods
        std::vector<mtl_bundle_m::mtl_bundle_t> buildMTLBundles(const std::vector<mtl_bundle_m::transmission_line_bundle_t>& line_bundles) {
            std::vector<mtl_bundle_m::mtl_bundle_t> res;
            // Stub implementation
            for (const auto& lb : line_bundles) {
                mtl_bundle_m::mtl_bundle_t bundle;
                // Logic to populate bundle would go here
                res.push_back(bundle);
            }
            return res;
        }

        network_manager_m::network_manager_t buildNetworkManager(const std::vector<void*>& networks) {
            network_manager_m::network_manager_t res;
            // Stub implementation
            return res;
        }

        void buildNetwork() {
            // Stub
        }

        void connectNodeToGround() {
            // Stub
        }

        void connectNodes() {
            // Stub
        }

        void connectNodesToNetworkCircuit() {
            // Stub
        }

        void addNodeWithId() {
            // Stub
        }

        void addProbesWithId(const std::vector<mtln_types_m::probe_t>& probes) {
            // Stub
        }

        void addGenerators(const std::vector<void*>& generators) {
            // Stub
        }
    };

    // Helper functions declared outside the class for the interface
    std::vector<mtl_bundle_m::cable_bundle_t> buildCableBundles(const std::vector<void*>& cables);
    std::vector<mtl_bundle_m::transmission_line_bundle_t> buildLineBundles(const std::vector<mtl_bundle_m::cable_bundle_t>& cable_bundles, double dt, const std::vector<mtl_m::XYZlimit_t>& alloc);
    std::vector<mtl_bundle_m::transmission_line_bundle_t> buildLineBundles(const std::vector<mtl_bundle_m::cable_bundle_t>& cable_bundles, double dt);
    fhash::fhash_tbl_t mapCablesToBundlesId(const std::vector<mtl_bundle_m::transmission_line_bundle_t>& line_bundles, const std::vector<mtl_bundle_m::mtl_bundle_t>& bundles);

    preprocess_t preprocess(const mtln_types_m::mtln_t& parsed, const std::optional<std::vector<mtl_m::XYZlimit_t>>& alloc = std::nullopt) {
        preprocess_t res;
        fhash::fhash_tbl_t cable_name_to_bundle_id;
        std::vector<mtl_bundle_m::transmission_line_bundle_t> line_bundles;
        std::vector<mtl_bundle_m::cable_bundle_t> cable_bundles;

#ifdef CompileWithMPI
        int ierr, rank;
        MPI_COMM_RANK(SUBCOMM_MPI, &rank, &ierr);
#endif

        res.final_time = parsed.time_step * parsed.number_of_steps;
        res.dt = parsed.time_step;

#ifdef CompileWithMPI
        MPI_Barrier(subcomm_mpi, &ierr);
#endif

        // Group cables into bundles
        cable_bundles = buildCableBundles(parsed.cables);

#ifdef CompileWithMPI
        MPI_Barrier(subcomm_mpi, &ierr);
#endif

        // Create mtl objects from cables
        if (alloc.has_value()) {
            line_bundles = buildLineBundles(cable_bundles, res.dt, alloc.value());
        } else {
            line_bundles = buildLineBundles(cable_bundles, res.dt);
        }

#ifdef CompileWithMPI
        MPI_Barrier(subcomm_mpi, &ierr);
#endif

        // create mtl_bundles from mtl objects
        res.bundles = res.buildMTLBundles(line_bundles);
        if (res.bundles.empty()) {
            return res;
        }

        res.cable_name_to_bundle_id = mapCablesToBundlesId(line_bundles, res.bundles);

        // Probes
        res.addProbesWithId(parsed.probes);

        // Generators
        res.addGenerators(parsed.wireGenerators);

        res.network_manager = res.buildNetworkManager(parsed.networks);

        return res;
    }

    std::vector<int> conductorsInLevel(const mtl_bundle_m::transmission_line_bundle_t& line) {
        std::vector<int> res(line.levels.size(), 0);
        for (size_t i = 0; i < line.levels.size(); ++i) {
            for (const auto& l : line.levels[i].lines) {
                res[i] += l.number_of_conductors;
            }
        }
        return res;
    }

    int findConductorsBeforeCable(const std::string& name, const mtl_bundle_m::transmission_line_level_t& level) {
        int res = 0;
        for (const auto& l : level.lines) {
            if (l.name != name) {
                res += l.number_of_conductors;
            } else {
                return res;
            }
        }
        return res;
    }

    int findOuterConductorNumber(const mtl_bundle_m::mtl_t& line, const mtl_bundle_m::transmission_line_level_t& level, int conductors_in_level) {
        return findConductorsBeforeCable(line.parent_name, level) +
               conductors_in_level +
               line.conductor_in_parent;
    }

    std::vector<int> findInnerConductorRange(const mtl_bundle_m::mtl_t& line, const mtl_bundle_m::transmission_line_level_t& level, int conductors_in_level) {
        int start = findConductorsBeforeCable(line.name, level) + conductors_in_level;
        std::vector<int> res(line.number_of_conductors);
        for (int k = 0; k < line.number_of_conductors; ++k) {
            res[k] = start + k + 1; // 1-based indexing implied by Fortran logic usually, but vector is 0-based. 
            // Fortran: [(k, k = 1, line%number_of_conductors)] generates 1, 2, ..., N
            // If we want to preserve exact values, we generate 1..N.
            // However, the previous calculation `findConductorsBeforeCable` returns a count.
            // Let's assume the range is relative to the start of the cable's conductors in the bundle.
            // The Fortran code: `findConductorsBeforeCable(...) + conductors_in_level + [(k, k = 1, N)]`
            // This implies the result is a vector of indices.
            res[k] = start + k; 
        }
        return res;
    }

    void setBundleTransferImpedance(mtl_bundle_m::mtl_bundle_t& bundle, const mtl_bundle_m::transmission_line_bundle_t& line) {
        std::vector<int> conductors_in_level = conductorsInLevel(line);
        bundle.conductors_in_level = conductors_in_level;

        for (size_t i = 1; i < line.levels.size(); ++i) {
            for (const auto& l : line.levels[i].lines) {
                int conductor_out = findOuterConductorNumber(l, line.levels[i-1], 0); // Sum of previous levels handled inside findOuter? No, passed as arg.
                // Wait, the Fortran call: `findOuterConductorNumber(line%levels(i)%lines(j), line%levels(i-1), sum(conductors_in_level(1:i-2)))`
                // Note: Fortran arrays are 1-based. `1:i-2` when i=2 is empty sum = 0.
                // In C++, `conductors_in_level` is 0-based. `1:i-2` corresponds to indices `0` to `i-3`.
                // Let's recalculate the sum properly.
                
                int sum_prev = 0;
                for (size_t k = 0; k < i - 1; ++k) { // 1 to i-2 in 1-based is 0 to i-3 in 0-based? 
                    // If i=2 (3rd level in 1-based? No, i starts at 2 in Fortran loop `do i = 2`).
                    // i=2 (Fortran) -> index 1 (C++). `1:i-2` -> `1:0` -> empty. Sum=0.
                    // i=3 (Fortran) -> index 2 (C++). `1:1` -> index 0. Sum=conductors_in_level[0].
                    // So we sum indices 0 to i-2 (exclusive of i-1).
                    // Actually, let's just map the logic directly.
                    // Fortran `sum(conductors_in_level(1:i-2))`
                    // If i=2, range is 1 to 0 -> 0.
                    // If i=3, range is 1 to 1 -> index 0.
                    // If i=4, range is 1 to 2 -> indices 0, 1.
                    // So we sum `conductors_in_level[0]` to `conductors_in_level[i-3]`.
                    if (k < i - 1) { // This logic is tricky. Let's stick to the loop variable `i` from Fortran.
                         // In C++, `i` is index. Fortran `i` is index+1.
                         // Let's use 0-based `idx` for levels.
                    }
                }
                
                // Let's rewrite the loop using 0-based indexing for clarity
                // Original: do i = 2, size(line%levels)
                // C++: for (size_t i = 1; i < line.levels.size(); ++i)
                // sum(conductors_in_level(1:i-2)) -> sum of first i-1 elements?
                // i=2 (Fortran) -> sum(1:0) -> 0.
                // i=3 (Fortran) -> sum(1:1) -> element 1.
                // So for C++ index `i`, we sum elements `0` to `i-2`.
                
                int sum_prev_levels = 0;
                for (size_t k = 0; k < i - 1; ++k) {
                    sum_prev_levels += conductors_in_level[k];
                }

                int conductor_out = findOuterConductorNumber(l, line.levels[i-1], sum_prev_levels);
                
                // findInnerConductorRange(line%levels(i)%lines(j), line%levels(i), sum(conductors_in_level(1:i-1)))
                // sum(1:i-1) -> sum of first i elements (indices 0 to i-1).
                int sum_curr_levels = 0;
                for (size_t k = 0; k < i; ++k) {
                    sum_curr_levels += conductors_in_level[k];
                }
                
                std::vector<int> range_in = findInnerConductorRange(l, line.levels[i], sum_curr_levels);
                
                bundle.addTransferImpedance(conductor_out, range_in, l.transfer_impedance);
            }
        }

        if (line.levels.size() > 1) {
            for (const auto& l : line.levels[1].lines) { // Level 2 in Fortran is index 1 in C++
                // findOuterConductorNumber(line%levels(2)%lines(j), line%levels(1), sum(conductors_in_level(1:0)))
                // sum(1:0) is 0.
                int conductor_out = findOuterConductorNumber(l, line.levels[0], 0);
                
                // findInnerConductorRange(line%levels(2)%lines(j), line%levels(2), sum(conductors_in_level(1:1)))
                // sum(1:1) is element 1 (index 0).
                int sum_curr = conductors_in_level[0];
                std::vector<int> range_in = findInnerConductorRange(l, line.levels[1], sum_curr);
                
                int conductor_in_parent = l.conductor_in_parent;
                
                if (!l.initial_connector_transfer_impedances.empty()) {
                    // Assuming index `conductor_in_parent` is valid. 
                    // Note: Fortran might be 1-based, C++ 0-based. Adjust if necessary.
                    // Assuming the stub handles indexing or it's 0-based.
                    if (conductor_in_parent < l.initial_connector_transfer_impedances.size()) {
                        // Stub access
                        // const auto& zt = l.initial_connector_transfer_impedances[conductor_in_parent];
                        // if (zt.has_transfer_impedance()) {
                        //     bundle.setConnectorTransferImpedance(1, conductor_out, range_in, zt);
                        // }
                    }
                }
                
                if (!l.end_connector_transfer_impedances.empty()) {
                    if (conductor_in_parent < l.end_connector_transfer_impedances.size()) {
                        // const auto& zt = l.end_connector_transfer_impedances[conductor_in_parent];
                        // if (zt.has_transfer_impedance()) {
                        //     bundle.setConnectorTransferImpedance(bundle.du.size(), conductor_out, range_in, zt);
                        // }
                    }
                }
            }
        }
    }

    void mapConductorsBeforeCable(fhash::fhash_tbl_t& conductors_before_cable, const mtl_bundle_m::transmission_line_bundle_t& line) {
        std::vector<int> conductors_in_level = conductorsInLevel(line);
        // Implementation continues...
    }

}

conductors_before_cable.set(key(line.levels[0].lines[0].name), 0);
        for (i = 1; i < line.levels.size(); i++) {
            for (j = 0; j < line.levels[i].lines.size(); j++) {
                range_in = findInnerConductorRange(line.levels[i].lines[j], line.levels[i], std::accumulate(conductors_in_level.begin(), conductors_in_level.begin() + i, 0));
                if (range_in.size() != 0) {
                    conductors_before_cable.set(key(line.levels[i].lines[j].name), range_in[0] - 1);
                } else {
                    throw std::runtime_error("range in cannot be empty");
                }
            }
        }
    }

    std::vector<mtl_bundle_t> buildMTLBundles(preprocess_t& this_obj, const std::vector<transmission_line_bundle_t>& lines) {
        mtl_bundle_t res[lines.size()];
        fhash_tbl_t conductors_before_cable;
        int i;
#ifdef CompileWithMPI
        int ierr;
#endif

#ifdef CompileWithMPI
        mpi_barrier(subcomm_mpi, ierr);
#endif
        res.resize(lines.size());
        for (i = 0; i < lines.size(); i++) {
            res[i] = mtldCtor(lines[i].levels, lines[i].levels[0].lines[0].name);
            if (res[i].dt < this_obj.dt) {
                this_obj.dt = res[i].dt;
            }
            setBundleTransferImpedance(res[i], lines[i]);
            mapConductorsBeforeCable(conductors_before_cable, lines[i]);
        }
        this_obj.conductors_before_cable = conductors_before_cable;
        return res;
    }

    mtl_t buildLineFromCable(cable_t& cable, double dt, const std::vector<std::vector<int>>& layer_indices, bool bundle_in_layer, const std::vector<int>& alloc_z) {
        mtl_t res;
        
        int conductor_in_parent = 0;
        std::string parent_name;

        if (auto* shielded = dynamic_cast<shielded_multiwire_t*>(&cable)) {
            if (shielded->parent_cable) {
                parent_name = shielded->parent_cable->name;
                conductor_in_parent = shielded->conductor_in_parent;
            } else {
                parent_name = "unassigned_parent";
                conductor_in_parent = -1;
            }
            res = mtl_shielded(
                lpul = shielded->inductance_per_meter,
                cpul = shielded->capacitance_per_meter,
                rpul = shielded->resistance_per_meter,
                gpul = shielded->conductance_per_meter,
                step_size = shielded->step_size,
                name = shielded->name,
                segments = shielded->segments,
                dt = dt,
                parent_name = parent_name,
                conductor_in_parent = conductor_in_parent,
                transfer_impedance = shielded->transfer_impedance
#ifdef CompileWithMPI
                , layer_indices = layer_indices, bundle_in_layer = bundle_in_layer, alloc_z = alloc_z
#endif
            );
        } else if (auto* unshielded = dynamic_cast<unshielded_multiwire_t*>(&cable)) {
            res = mtl_unshielded(
                lpul = unshielded->cell_inductance_per_meter,
                cpul = unshielded->cell_capacitance_per_meter,
                rpul = unshielded->resistance_per_meter,
                gpul = unshielded->conductance_per_meter,
                step_size = unshielded->step_size,
                name = unshielded->name,
                segments = unshielded->segments,
                dt = dt,
                multipolar_expansion = unshielded->multipolar_expansion,
                radius = unshielded->radius
#ifdef CompileWithMPI
                , layer_indices = layer_indices, bundle_in_layer = bundle_in_layer, alloc_z = alloc_z
#endif
            );
        }
        if (cable.initial_connector) {
            addInitialConnector(res, *cable.initial_connector);
        } else {
            res.initial_connector_transfer_impedances.resize(0);
        }
        if (cable.end_connector) {
            addEndConnector(res, *cable.end_connector);
        } else {
            res.end_connector_transfer_impedances.resize(0);
        }
        
        return res;
    }

    std::vector<transmission_line_bundle_t> buildLineBundles(const std::vector<cable_bundle_t>& cable_bundles, double dt, const std::vector<XYZlimit_t>& alloc) {
        transmission_line_bundle_t res[cable_bundles.size()];
        int i, j, k;
        int nb, nl, nc;
        std::vector<std::vector<int>> layer_indices;
        bool bundle_in_layer = false;
        std::vector<int> alloc_z(2);
        if (alloc.size() > 0) {
            alloc_z[0] = alloc[2].zi;
            alloc_z[1] = alloc[2].ze;
        }
        nb = cable_bundles.size();
        res.resize(nb);
        for (i = 0; i < nb; i++) {
            if (alloc.size() > 0) {
                if (layer_indices.size() > 0) layer_indices.clear();
                bundle_in_layer = isBundleInLayer(cable_bundles[i].levels[0].cables[0].ptr, alloc_z);
                if (bundle_in_layer) {
                    layer_indices = findIndicesInLayer(cable_bundles[i].levels[0].cables[0].ptr, alloc_z);
                } else {
                    layer_indices.resize(0, 2);
                    for (auto& row : layer_indices) {
                        row[0] = 0;
                        row[1] = 0;
                    }
                }
            }
            nl = cable_bundles[i].levels.size();
            res[i].levels.resize(nl);
            for (j = 0; j < nl; j++) {
                nc = cable_bundles[i].levels[j].cables.size();
                res[i].levels[j].lines.resize(nc);
                for (k = 0; k < nc; k++) {
                    if (alloc.size() > 0) {
                        res[i].levels[j].lines[k] = buildLineFromCable(*cable_bundles[i].levels[j].cables[k].ptr, dt, layer_indices, bundle_in_layer, alloc_z);
                    } else {
                        res[i].levels[j].lines[k] = buildLineFromCable(*cable_bundles[i].levels[j].cables[k].ptr, dt);
                    }
                }
            }
        }
        return res;
    }

    bool isBundleInLayer(cable_t& cable, const std::vector<int>& alloc_z) {
        int n = 0;
        bool in_layer = false;
        for (int i = 0; i < cable.segments.size(); i++) {
            if (isSegmentWithinAllocBox(cable.segments, i, alloc_z)) {
                if (!in_layer) {
                    in_layer = true;
                }
            } else {
                if (in_layer) {
                    in_layer = false;
                    n++;
                }
            }
        }
        if (in_layer) n++;
        return (n != 0);
    }

    std::vector<std::vector<int>> findIndicesInLayer(cable_t& cable, const std::vector<int>& alloc_z) {
        // Implementation details omitted as they are not fully provided in the snippet
        return {};
    }

// Note: The following code assumes the existence of the class definitions
        // and helper functions mentioned in the Fortran code (e.g., cable_t, segment_t, etc.)
        // and translates the logic into C++ functions within a namespace.

        // Helper function declaration (assumed to be defined elsewhere or in the same namespace)
        // bool isSegmentWithinAllocBox(const std::vector<segment_t>& segs, int i, const std::vector<int>& z);

        std::vector<std::vector<int>> getSegmentRanges(const cable_t& cable, const std::vector<int>& alloc_z) {
            int n = 0;
            // Precount
            for (int i = 0; i < static_cast<int>(cable.segments.size()); ++i) {
                if (isSegmentWithinAllocBox(cable.segments, i, alloc_z)) {
                    // Logic handled in the loop below to count transitions
                }
            }
            
            // Reset and count properly
            n = 0;
            bool in_layer = false;
            for (int i = 0; i < static_cast<int>(cable.segments.size()); ++i) {
                if (isSegmentWithinAllocBox(cable.segments, i, alloc_z)) {
                    if (!in_layer) {
                        in_layer = true;
                    }
                } else {
                    if (in_layer) {
                        in_layer = false;
                        n++;
                    }
                }
            }
            if (in_layer) n++;

            std::vector<std::vector<int>> res(n, std::vector<int>(2));
            int current_n = 0;
            in_layer = false;
            for (int i = 0; i < static_cast<int>(cable.segments.size()); ++i) {
                if (isSegmentWithinAllocBox(cable.segments, i, alloc_z)) {
                    if (!in_layer) {
                        res[current_n][0] = i + 1; // Fortran is 1-based, assuming output should match or be adjusted. 
                                                   // Fortran code: res(n,1) = i. If i is 1-based index, it stays 1-based.
                        in_layer = true;
                    }
                } else {
                    if (in_layer) {
                        res[current_n][1] = i; // Fortran code: res(n, 2) = i-1. 
                                               // If i is the index of the segment NOT in layer, the previous one was i-1.
                        in_layer = false;
                        current_n++;
                    }
                }
            }
            if (in_layer) {
                res[current_n][1] = static_cast<int>(cable.segments.size()); // Fortran: i - 1 where i is loop limit + 1? 
                                                                            // Fortran loop: do i = 1, size. If in_layer at end, res(n,2) = i - 1.
                                                                            // If loop finishes, i becomes size + 1. So res(n,2) = size.
            }
            
            return res;
        }

        bool isSegmentWithinAllocBox(const std::vector<segment_t>& segs, int i, const std::vector<int>& z) {
            // Fortran: segs(i)%z >= z(1) .and. segs(i)%z <= z(2)
            // Assuming i is 1-based index in Fortran, convert to 0-based for C++ vector
            if (i < 1 || i > static_cast<int>(segs.size())) return false;
            return (segs[i - 1].z >= z[0]) && (segs[i - 1].z <= z[1]);
        }

        cable_bundle_t buildCableBundleFromParent(const cable_abstract_t& parent, const std::vector<cable_abstract_t>& cables) {
            cable_bundle_t res;
            cable_level_t level;

            res.levels.resize(1);
            res.levels[0].cables.resize(1);

            level.cables.resize(1);
            level.cables[0].ptr = parent.ptr;

            res.levels[0] = level;

            while (findNextLevel(level, cables) != 0) {
                appendLevel(res.levels, level);
            }

            return res;
        }

        void appendLevel(std::vector<cable_level_t>& levels, const cable_level_t& newLevel) {
            std::vector<cable_level_t> oldLevels = levels;
            levels.resize(oldLevels.size() + 1);
            for (size_t i = 0; i < oldLevels.size(); ++i) {
                levels[i] = oldLevels[i];
            }
            levels[oldLevels.size()] = newLevel;
        }

        int findNextLevel(cable_level_t& curr_level, const std::vector<cable_abstract_t>& cs) {
            cable_level_t next_level;
            int next_level_size = 0;

            // Count size
            for (size_t i = 0; i < curr_level.cables.size(); ++i) {
                for (size_t j = 0; j < cs.size(); ++j) {
                    cable_t* ptr = cs[j].ptr;
                    if (shielded_multiwire_t* shielded = dynamic_cast<shielded_multiwire_t*>(ptr)) {
                        if (shielded->parent_cable == curr_level.cables[i].ptr) {
                            next_level_size++;
                        }
                    }
                }
            }

            next_level.cables.resize(next_level_size);
            int n = 0;
            for (size_t i = 0; i < curr_level.cables.size(); ++i) {
                for (size_t j = 0; j < cs.size(); ++j) {
                    cable_t* ptr = cs[j].ptr;
                    if (shielded_multiwire_t* shielded = dynamic_cast<shielded_multiwire_t*>(ptr)) {
                        if (shielded->parent_cable == curr_level.cables[i].ptr) {
                            n++;
                            next_level.cables[n - 1].ptr = cs[j].ptr;
                        }
                    }
                }
            }

            curr_level = next_level;
            return static_cast<int>(curr_level.cables.size());
        }

        std::vector<cable_abstract_t> findParentCables(const std::vector<cable_abstract_t>& cables) {
            std::vector<int> parent_ids;
            
            // MPI barrier omitted as per instruction to convert #ifdef to C++ logic 
            // (assuming single thread or external synchronization if needed, but here just skipping MPI call)

            for (size_t i = 0; i < cables.size(); ++i) {
                cable_t* ptr = cables[i].ptr;
                if (unshielded_multiwire_t* unshielded = dynamic_cast<unshielded_multiwire_t*>(ptr)) {
                    parent_ids.push_back(static_cast<int>(i + 1)); // 1-based index stored
                } else if (shielded_multiwire_t* shielded = dynamic_cast<shielded_multiwire_t*>(ptr)) {
                    if (shielded->parent_cable == nullptr) {
                        parent_ids.push_back(static_cast<int>(i + 1));
                    }
                }
            }

            std::vector<cable_abstract_t> res(parent_ids.size());
            for (size_t i = 0; i < parent_ids.size(); ++i) {
                // parent_ids are 1-based, convert to 0-based for access
                res[i].ptr = cables[parent_ids[i] - 1].ptr;
            }

            return res;
        }

        std::vector<cable_abstract_t> buildCableBundles(const std::vector<cable_abstract_t>& cables) {
            std::vector<cable_abstract_t> parents = findParentCables(cables);
            std::vector<cable_bundle_t> cable_bundles(parents.size());
            
            for (size_t i = 0; i < parents.size(); ++i) {
                cable_bundles[i] = buildCableBundleFromParent(parents[i], cables);
            }

            // Note: The Fortran function returns cable_bundles but the variable name in result is cable_bundles.
            // The C++ function returns std::vector<cable_abstract_t> in the signature above, which is incorrect based on Fortran.
            // Fortran: function buildCableBundles(cables) result(cable_bundles)
            //          type(cable_abstract_t), dimension(:), allocatable :: cable_bundles
            // Wait, looking at the code:
            // allocate(cable_bundles(size(parents)))
            // cable_bundles(i) = buildCableBundleFromParent(...)
            // buildCableBundleFromParent returns type(cable_bundle_t)
            // So cable_bundles should be type(cable_bundle_t).
            // The Fortran declaration says: type(cable_abstract_t), dimension(:), allocatable :: cable_bundles
            // This looks like a typo in the Fortran code provided or my interpretation.
            // Let's look closer:
            // function buildCableBundles(cables) result(cable_bundles)
            // ...
            // type(cable_abstract_t), dimension(:), allocatable :: cable_bundles
            // ...
            // cable_bundles(i) = buildCableBundleFromParent(...)
            // buildCableBundleFromParent returns type(cable_bundle_t).
            // Assigning cable_bundle_t to cable_abstract_t? This implies polymorphism or a bug in the source.
            // Assuming the return type of buildCableBundles is std::vector<cable_bundle_t> based on usage.
            
            // Re-reading Fortran:
            // type(cable_abstract_t), dimension(:), allocatable :: cable_bundles
            // This is the RESULT variable declaration.
            // But it is assigned a cable_bundle_t.
            // This is likely a mistake in the provided Fortran snippet or cable_bundle_t inherits from cable_abstract_t.
            // Given "class (cable_t)" and "type(cable_abstract_t)", it's possible cable_bundle_t is a cable_abstract_t.
            // However, usually bundles are containers.
            // Let's assume the return type is std::vector<cable_bundle_t> as that matches the logic.
            // If the Fortran strictly requires returning cable_abstract_t, we'd need to cast, but that's unsafe.
            // I will return std::vector<cable_bundle_t> as it's the logical type.
            
            // Correction: The prompt asks to translate ONLY this part. I will output the function returning the correct type.
            // If the user insists on strict translation, I might need to return std::vector<cable_abstract_t> and cast, 
            // but that's bad C++. I'll stick to the logical type.
            
            // Actually, let's look at the variable name in the result: `cable_bundles`.
            // The declaration is `type(cable_abstract_t), dimension(:), allocatable :: cable_bundles`.
            // This is very suspicious. I will assume it should be `type(cable_bundle_t)`.
            
            return cable_bundles; // This line is unreachable in the snippet provided because the return type mismatch prevents compilation in C++ if we follow strict typing.
                                  // I will define the function to return std::vector<cable_bundle_t>.
        }
        
        // Redefining buildCableBundles with correct return type
        std::vector<cable_bundle_t> buildCableBundlesCorrected(const std::vector<cable_abstract_t>& cables) {
            std::vector<cable_abstract_t> parents = findParentCables(cables);
            std::vector<cable_bundle_t> cable_bundles(parents.size());
            
            for (size_t i = 0; i < parents.size(); ++i) {
                cable_bundles[i] = buildCableBundleFromParent(parents[i], cables);
            }

            return cable_bundles;
        }

        fhash_tbl_t mapCablesToBundlesId(const std::vector<transmission_line_bundle_t>& lines, const std::vector<mtl_bundle_t>& bundles) {
            fhash_tbl_t res;
            for (size_t i = 0; i < lines.size(); ++i) {
                for (size_t j = 0; j < lines[i].levels.size(); ++j) {
                    for (size_t k = 0; k < lines[i].levels[j].lines.size(); ++k) {
                        res.set(lines[i].levels[j].lines[k].name, static_cast<int>(i + 1)); // 1-based index
                    }
                }
            }
            return res;
        }

        fhash_tbl_t mapCablesToBundles(const std::vector<transmission_line_bundle_t>& lines, const std::vector<mtl_bundle_t>& bundles) {
            fhash_tbl_t res;
            for (size_t i = 0; i < lines.size(); ++i) {
                for (size_t j = 0; j < lines[i].levels.size(); ++j) {
                    for (size_t k = 0; k < lines[i].levels[j].lines.size(); ++k) {
                        res.set(lines[i].levels[j].lines[k].name, bundles[i]);
                    }
                }
            }
            return res;
        }

        std::vector<std::string> writeParallelRLCnode(const nw_node_t& node, const termination_t& termination, const std::string& end_node) {
            std::vector<std::string> res;
            std::string buff;
            std::string termination_r, termination_l, termination_c, line_c, line_g;

            termination_c = std::to_string(termination.capacitance);
            termination_r = std::to_string(termination.resistance);
            termination_l = std::to_string(termination.inductance);
            line_c = std::to_string(node.line_c_per_meter * node.step / 2);
            
            // The function returns an empty vector in the Fortran snippet provided
            return res;
        }

        std::vector<std::string> writeSeriesRLCnode(const nw_node_t& node, const termination_t& termination, const std::string& end_node) {
            std::vector<std::string> res;
            std::string buff;
            // Function body is empty in the provided snippet
            return res;
        }

std::string termination_r, termination_l, termination_c, line_c, line_g, generator_r;

        termination_c = std::to_string(termination.capacitance);
        termination_r = std::to_string(termination.resistance);
        termination_l = std::to_string(termination.inductance);
        line_c = std::to_string(node.line_c_per_meter * node.step / 2);
        res.resize(0);
        if (termination.source.path_to_excitation != "") {
            buff = "R" + node.name + " " + node.name + " " + node.name + "_S " + termination_r;
            appendToStringArray(res, buff);
            buff = "L" + node.name + " " + node.name + " " + node.name + "_S " + termination_l;
            appendToStringArray(res, buff);
            buff = "C" + node.name + " " + node.name + " " + node.name + "_S " + termination_c;
            appendToStringArray(res, buff);

            generator_r = std::to_string(termination.source.resistance);
            if (termination.source.source_type == SOURCE_TYPE_VOLTAGE) {
                buff = "V" + node.name + "_S " + node.name + "_S " + node.name + "_genR dc 0";
                appendToStringArray(res, buff);
                buff = "R" + node.name + "_S " + node.name + "_genR " + end_node + " " + generator_r;
                appendToStringArray(res, buff);
            } else if (termination.source.source_type == SOURCE_TYPE_CURRENT) {
                buff = "I" + node.name + "_S " + end_node + " " + node.name + "_S dc 0";
                appendToStringArray(res, buff);
                if (termination.source.resistance != 1.0e22_rkind) {
                    buff = "R" + node.name + "_S " + end_node + " " + node.name + "_S " + generator_r;
                    appendToStringArray(res, buff);
                }
            }
        } else {
            buff = "R" + node.name + " " + node.name + " " + end_node + termination_r;
            appendToStringArray(res, buff);
            buff = "L" + node.name + " " + node.name + " " + end_node + termination_l;
            appendToStringArray(res, buff);
            buff = "C" + node.name + " " + node.name + " " + end_node + termination_c;
            appendToStringArray(res, buff);
        }
        buff = "I" + node.name + " " + node.name + " 0 dc 0";
        appendToStringArray(res, buff);
        buff = "CL" + node.name + " " + node.name + " 0 " + line_c;
        appendToStringArray(res, buff);

        if (node.line_g_per_meter != 0) {
            line_g = std::to_string(1.0 / (node.line_g_per_meter * node.step / 2));
            buff = "GL" + node.name + " " + node.name + " 0 " + line_g;
            appendToStringArray(res, buff);
        }

    }

    std::vector<std::string> writeNetwork_circuitNode(const nw_node_t& node, const termination_t& termination, const std::string& end_node) {
        std::vector<std::string> res;
        std::string buff;
        std::string line_c, line_g, short_r, generator_r;
        short_r = std::to_string(1e-10);
        line_c = std::to_string(node.line_c_per_meter * node.step / 2);
        res.resize(0);

        if (termination.source.path_to_excitation != "") {
            generator_r = std::to_string(termination.source.resistance);
            buff = "R" + node.name + " " + node.name + " " + node.name + "_S " + short_r;
            appendToStringArray(res, buff);
            if (termination.source.source_type == SOURCE_TYPE_VOLTAGE) {
                buff = "V" + node.name + "_S " + node.name + "_S " + node.name + "_genR dc 0";
                appendToStringArray(res, buff);
                buff = "R" + node.name + "_S " + node.name + "_genR " + end_node + " " + generator_r;
                appendToStringArray(res, buff);
            } else if (termination.source.source_type == SOURCE_TYPE_CURRENT) {
                buff = "I" + node.name + "_S " + end_node + " " + node.name + "_S dc 0";
                appendToStringArray(res, buff);
                if (termination.source.resistance != 1.0e22_rkind) {
                    buff = "R" + node.name + "_S " + end_node + " " + node.name + "_S " + generator_r;
                    appendToStringArray(res, buff);
                }
            }
        } else {
            buff = "R" + node.name + " " + node.name + " " + end_node + " " + short_r;
            appendToStringArray(res, buff);
        }

        buff = "I" + node.name + " " + node.name + " 0 dc 0";
        appendToStringArray(res, buff);
        buff = "CL" + node.name + " " + node.name + " 0 " + line_c;
        appendToStringArray(res, buff);

        if (node.line_g_per_meter != 0) {
            line_g = std::to_string(1.0 / (node.line_g_per_meter * node.step / 2));
            buff = "GL" + node.name + " " + node.name + " 0 " + line_g;
            appendToStringArray(res, buff);
        }

        return res;
    }

    std::vector<std::string> writeModelNode(const nw_node_t& node, const termination_t& termination, const std::string& end_node) {
        std::vector<std::string> res;
        std::string buff;
        std::string model_name, model_file;
        std::string line_c, line_g, generator_r;
        line_c = std::to_string(node.line_c_per_meter * node.step / 2);
        res.resize(0);

        model_name = termination.model.name;
        model_file = termination.model.file;

        buff = ".include " + model_file;
        appendToStringArray(res, buff);

        if (termination.source.path_to_excitation != "") {
            generator_r = std::to_string(termination.source.resistance);
            if (termination.source.source_type == SOURCE_TYPE_VOLTAGE) {
                buff = "V" + node.name + "_S " + node.name + " " + node.name + "_genR dc 0";
                appendToStringArray(res, buff);
                buff = "R" + node.name + "_S " + node.name + "_genR " + node.name + "_S " + generator_r;
                appendToStringArray(res, buff);
            } else if (termination.source.source_type == SOURCE_TYPE_CURRENT) {
                buff = "I" + node.name + "_S " + node.name + "_S " + node.name + " dc 0";
                appendToStringArray(res, buff);
                if (termination.source.resistance != 1.0e22_rkind) {
                    buff = "R" + node.name + "_S " + node.name + "_S " + node.name + " " + generator_r;
                    appendToStringArray(res, buff);
                }
            }
            buff = "x" + node.name + " " + node.name + "_S " + end_node + " " + model_name;
            appendToStringArray(res, buff);
        } else {
            buff = "x" + node.name + " " + node.name + " " + end_node + " " + model_name;
            appendToStringArray(res, buff);
        }

        buff = "I" + node.name + " " + node.name + " 0 dc 0";
        appendToStringArray(res, buff);
        buff = "CL" + node.name + " " + node.name + " 0 " + line_c;
        appendToStringArray(res, buff);

        if (node.line_g_per_meter != 0) {
            line_g = std::to_string(1.0 / (node.line_g_per_meter * node.step / 2));

buff = "GL" + node.name + " " + node.name + " 0 " + line_g;
        appendToStringArray(res, buff);
    }


}

std::vector<std::string> writeSeriesRLnode(const nw_node_t& node, const termination_t& termination, const std::string& end_node) {
    std::vector<std::string> res;
    std::string buff;
    std::string termination_r, termination_l, line_c, line_g, generator_r;

    termination_r = std::to_string(termination.resistance);
    termination_l = std::to_string(termination.inductance);
    line_c = std::to_string(node.line_c_per_meter * node.step / 2);
    res.reserve(1);

    buff = "R" + node.name + " " + node.name + "_R " + node.name + " " + termination_r;
    appendToStringArray(res, buff);
    if (termination.source.path_to_excitation != "") {
        generator_r = std::to_string(termination.source.resistance);
        buff = "L" + node.name + " " + node.name + "_R " + node.name + "_S" + " " + termination_l;
        appendToStringArray(res, buff);
        if (termination.source.source_type == SOURCE_TYPE_VOLTAGE) {
            buff = "V" + node.name + "_S " + node.name + "_S " + node.name + "_genR dc 0";
            appendToStringArray(res, buff);
            buff = "R" + node.name + "_S " + node.name + "_genR " + end_node + " " + generator_r;
            appendToStringArray(res, buff);
        } else if (termination.source.source_type == SOURCE_TYPE_CURRENT) {
            buff = "I" + node.name + "_S " + end_node + " " + node.name + "_S dc 0";
            appendToStringArray(res, buff);
            if (termination.source.resistance != 1.0e22) {
                buff = "R" + node.name + "_S " + end_node + " " + node.name + "_S " + generator_r;
                appendToStringArray(res, buff);
            }
        }
    } else {
        buff = "L" + node.name + " " + node.name + "_R " + end_node + " " + termination_l;
        appendToStringArray(res, buff);
    }
    buff = "I" + node.name + " " + node.name + " 0 dc 0";
    appendToStringArray(res, buff);
    buff = "CL" + node.name + " " + node.name + " 0 " + line_c;
    appendToStringArray(res, buff);

    if (node.line_g_per_meter != 0) {
        line_g = std::to_string(1.0 / (node.line_g_per_meter * node.step / 2));
        buff = "GL" + node.name + " " + node.name + " 0 " + line_g;
        appendToStringArray(res, buff);
    }

    return res;
}


std::vector<std::string> writeSeriesNode(const nw_node_t& node, const termination_t& termination, const std::string& end_node) {
    std::vector<std::string> res;

    if (termination.capacitance >= 1e22) {
        res = writeSeriesRLnode(node, termination, end_node);
    } else {
        res = writeSeriesRLCnode(node, termination, end_node);
    }

    return res;
}

void appendToStringArray(std::vector<std::string>& arr, const std::string& str) {
    // This has been implemented because there seems to be a bug in gfortran: 
    // https://fortran-lang.discourse.group/t/read-data-and-append-it-to-array-best-practice/1915
    // and arr = [ arr, str ] can't be used.
    std::vector<std::string> old_arr = arr;
    arr.clear();
    arr.resize(old_arr.size() + 1);
    for (size_t i = 0; i < old_arr.size(); ++i) {
        arr[i] = old_arr[i];
    }
    arr[arr.size() - 1] = str;
}

std::vector<std::string> writeShortNode(const nw_node_t& node, const termination_t& termination, const std::string& end_node) {
    std::vector<std::string> res;
    std::string buff;
    std::string short_R, line_c, line_g, generator_r;

    short_R = std::to_string(1e-10);
    line_c = std::to_string(node.line_c_per_meter * node.step / 2);

    res.reserve(1);
    if (termination.source.path_to_excitation != "") {
        generator_r = std::to_string(termination.source.resistance);
        buff = "R" + node.name + " " + node.name + " " + node.name + "_S" + " " + short_R;
        appendToStringArray(res, buff);
        if (termination.source.source_type == SOURCE_TYPE_VOLTAGE) {
            buff = "V" + node.name + "_S " + node.name + "_S " + node.name + "_genR dc 0";
            appendToStringArray(res, buff);
            buff = "R" + node.name + "_S " + node.name + "_genR " + end_node + " " + generator_r;
            appendToStringArray(res, buff);
        } else if (termination.source.source_type == SOURCE_TYPE_CURRENT) {
            buff = "I" + node.name + "_S " + end_node + " " + node.name + "_S dc 0";
            appendToStringArray(res, buff);
            if (termination.source.resistance != 1.0e22) {
                buff = "R" + node.name + "_S " + end_node + " " + node.name + "_S " + generator_r;
                appendToStringArray(res, buff);
            }
        }
    } else {
        buff = "R" + node.name + " " + node.name + " " + end_node + " " + short_R;
        appendToStringArray(res, buff);
    }
    buff = "I" + node.name + " " + node.name + " 0 dc 0";
    appendToStringArray(res, buff);
    buff = "CL" + node.name + " " + node.name + " 0 " + line_c;
    appendToStringArray(res, buff);
    
    if (node.line_g_per_meter != 0) {
        line_g = std::to_string(1.0 / (node.line_g_per_meter * node.step / 2));
        buff = "GL" + node.name + " " + node.name + " 0 " + line_g;
        appendToStringArray(res, buff);
    }

    return res;
}

std::vector<std::string> writeOpenNode(const nw_node_t& node, const termination_t& termination, const std::string& end_node) {
    std::vector<std::string> res;
    std::string buff;
    std::string line_c, line_g;

    line_c = std::to_string(node.line_c_per_meter * node.step / 2);

    res.reserve(1);
    buff = "R" + node.name + " " + node.name + " " + end_node + " 1e22";
    appendToStringArray(res, buff);
    buff = "I" + node.name + " " + node.name + " 0 dc 0";
    appendToStringArray(res, buff);
    buff = "CL" + node.name + " " + node.name + " 0 " + line_c;
    appendToStringArray(res, buff);
    
    if (node.line_g_per_meter != 0) {
        line_g = std::to_string(1.0 / (node.line_g_per_meter * node.step / 2));
        buff = "GL" + node.name + " " + node.name + " 0 " + line_g;
        appendToStringArray(res, buff);
    }

    return res;
}

std::vector<std::string> writeNodeDescription(const nw_node_t& node, const termination_t& termination, const std::string& end_node) {
    std::vector<std::string> res;
    if (termination.termination_type == TERMINATION_SERIES) {
        res = writeSeriesNode(node, termination, end_node);
    } else if (termination.termination_type == TERMINATION_PARALLEL) {
        res = writeParallelRLCNode(node, termination, end_node);

else if (termination.termination_type == TERMINATION_RsLCp) {
            res = writeXsYZpNode(node, termination, end_node, "RLC");
        }
        else if (termination.termination_type == TERMINATION_LsRCp) {
            res = writeXsYZpNode(node, termination, end_node, "LRC");
        }
        else if (termination.termination_type == TERMINATION_CsLRp) {
            res = writeXsYZpNode(node, termination, end_node, "CLR");
        }
        else if (termination.termination_type == TERMINATION_RLsCp) {
            res = writeXYsZpNode(node, termination, end_node, "RLC");
        }
        else if (termination.termination_type == TERMINATION_RCsLp) {
            res = writeXYsZpNode(node, termination, end_node, "RCL");
        }
        else if (termination.termination_type == TERMINATION_LCsRp) {
            res = writeXYsZpNode(node, termination, end_node, "LCR");
        }
        else if (termination.termination_type == TERMINATION_SHORT) {
            res = writeShortNode(node, termination, end_node);
        }
        else if (termination.termination_type == TERMINATION_OPEN) {
            res = writeOpenNode(node, termination, end_node);
        }
        else if (termination.termination_type == TERMINATION_CIRCUIT) {
            res = writeModelNode(node, termination, end_node);
        }
        else if (termination.termination_type == TERMINATION_NETWORK) {
            res = writeNetwork_circuitNode(node, termination, end_node);
        }
        else if (termination.termination_type == TERMINATION_UNDEFINED) {
            WarnErrReport("writeNodeDescription: undefined termination at " + node.name, true);
        }
    }

    std::vector<std::string> writeXYsZpNode(const nw_node_t& node, const termination_t& termination, const std::string& end_node, const std::string& XYZ) {
        std::vector<std::string> res;
        std::string buff;
        std::string node_name = node.name;
        std::string termination_x, termination_y, termination_z, line_c, line_g, generator_r;

        if (XYZ == "RLC" || XYZ == "LRC") {
            termination_x = std::to_string(termination.resistance);
            termination_y = std::to_string(termination.inductance);
            termination_z = std::to_string(termination.capacitance);
        }
        else if (XYZ == "LCR" || XYZ == "CLR") {
            termination_x = std::to_string(termination.inductance);
            termination_y = std::to_string(termination.capacitance);
            termination_z = std::to_string(termination.resistance);
        }
        else if (XYZ == "CRL" || XYZ == "RCL") {
            termination_x = std::to_string(termination.capacitance);
            termination_y = std::to_string(termination.resistance);
            termination_z = std::to_string(termination.inductance);
        }

        line_c = std::to_string(node.line_c_per_meter * node.step / 2.0);

        buff = XYZ.substr(0, 1) + node_name + " " + node_name + " " + node_name + "_X " + termination_x;
        appendToStringArray(res, buff);

        if (termination.source.path_to_excitation != "") {
            generator_r = std::to_string(termination.source.resistance);
            buff = XYZ.substr(1, 1) + node_name + " " + node_name + "_X " + node_name + "_S " + termination_y;
            appendToStringArray(res, buff);
            buff = XYZ.substr(2, 1) + node_name + " " + node_name + " " + node_name + "_S " + termination_z;
            appendToStringArray(res, buff);

            if (termination.source.source_type == SOURCE_TYPE_VOLTAGE) {
                buff = "V" + node_name + "_S " + node_name + "_S " + node_name + "_genR dc 0";
                appendToStringArray(res, buff);
                buff = "R" + node_name + "_S " + node_name + "_genR " + end_node + " " + generator_r;
                appendToStringArray(res, buff);
            }
            else if (termination.source.source_type == SOURCE_TYPE_CURRENT) {
                buff = "I" + node_name + "_S " + end_node + " " + node_name + "_S dc 0";
                appendToStringArray(res, buff);
                if (termination.source.resistance != 1.0e22) {
                    buff = "R" + node_name + "_S " + end_node + " " + node_name + "_S " + generator_r;
                    appendToStringArray(res, buff);
                }
            }
        }
        else {
            buff = XYZ.substr(1, 1) + node_name + " " + node_name + "_X " + end_node + " " + termination_y;
            appendToStringArray(res, buff);
            buff = XYZ.substr(2, 1) + node_name + " " + node_name + " " + end_node + " " + termination_z;
            appendToStringArray(res, buff);
        }

        buff = "I" + node_name + " " + node_name + " 0 dc 0";
        appendToStringArray(res, buff);
        buff = "CL" + node_name + " " + node_name + " 0 " + line_c;
        appendToStringArray(res, buff);

        if (node.line_g_per_meter != 0) {
            line_g = std::to_string(1.0 / (node.line_g_per_meter * node.step / 2.0));
            buff = "GL" + node_name + " " + node_name + " 0 " + line_g;
            appendToStringArray(res, buff);
        }

        return res;
    }

    std::vector<std::string> writeXsYZpNode(const nw_node_t& node, const termination_t& termination, const std::string& end_node, const std::string& XYZ) {
        std::vector<std::string> res;
        std::string buff;
        std::string node_name = node.name;
        std::string termination_x, termination_y, termination_z, line_c, line_g, generator_r;

        if (XYZ == "RLC" || XYZ == "RCL") {
            termination_x = std::to_string(termination.resistance);
            termination_y = std::to_string(termination.inductance);
            termination_z = std::to_string(termination.capacitance);
        }
        else if (XYZ == "LRC" || XYZ == "LCR") {
            termination_x = std::to_string(termination.inductance);
            termination_y = std::to_string(termination.resistance);
            termination_z = std::to_string(termination.capacitance);
        }
        else if (XYZ == "CLR" || XYZ == "CRL") {
            termination_x = std::to_string(termination.capacitance);
            termination_y = std::to_string(termination.resistance);
            termination_z = std::to_string(termination.inductance);
        }

        line_c = std::to_string(node.line_c_per_meter * node.step / 2.0);

        res.push_back(XYZ.substr(0, 1) + node_name + " " + node_name + " " + node_name + "_p " + termination_x);

        if (termination.source.path_to_excitation != "") {
            generator_r = std::to_string(termination.source.resistance);
            buff = XYZ.substr(1, 1) + node_name + " " + node_name + "_p " + node_name + "_S " + termination_y;
            appendToStringArray(res, buff);
            buff = XYZ.substr(2, 1) + node_name + " " + node_name + "_p " + node_name + "_S " + termination_z;
            appendToStringArray(res, buff);

            if (termination.source.source_type == SOURCE_TYPE_VOLTAGE) {
                buff = "V" + node_name + "_S " + node_name + "_S " + node_name + "_genR dc 0";
                appendToStringArray(res, buff);
                buff = "R" + node_name + "_S " + node_name + "_genR " + end_node + " " + generator_r;
                appendToStringArray(res, buff);
            }
            else if (termination.source.source_type == SOURCE_TYPE_CURRENT) {
                buff = "I" + node_name + "_S " + end_node + " " + node_name + "_S dc 0";
                appendToStringArray(res, buff);
                if (termination.source.resistance != 1.0e22) {
                    buff = "R" + node_name + "_S " + end_node + " " + node_name + "_S " + generator_r;
                    appendToStringArray(res, buff);
                }
            }
        }
        
        // Note: The Fortran code snippet ends abruptly here. 
        // Assuming the rest of the logic follows the pattern of writeXYsZpNode or similar structure,
        // but since the input is cut off, we stop here or complete the minimal valid syntax.
        // However, to be a valid function, it needs a return. 
        // Based on the cut-off, I will assume the rest of the function body is missing in the prompt.
        // I will provide the code up to the cut-off point wrapped in a valid function structure.
        // Since the prompt asks to translate ONLY this part, and the part is incomplete,
        // I will complete the function logically based on the previous function's structure if implied,
        // but strictly speaking, I should only translate what is there.
        // The Fortran code ends with an incomplete if statement. 
        // I will assume the user wants the translation of the provided text.
        // The provided text ends inside an if block. 
        // I will close the braces to make it valid C++ for the provided snippet.
        
        return res;
    }

appendToStringArray(res, buff);
                }
            }
        } else {
            buff = trim(XYZ.substr(1, 1) + node.name + " " + node.name + "_p " + end_node + " " + termination_y);
            appendToStringArray(res, buff);
            buff = trim(XYZ.substr(2, 1) + node.name + " " + node.name + "_p " + end_node + " " + termination_z);
            appendToStringArray(res, buff);
        }
        buff = trim("I" + node.name + " " + node.name + " 0 " + " dc 0");
        appendToStringArray(res, buff);
        buff = trim("CL" + node.name + " " + node.name + " 0 " + line_c);
        appendToStringArray(res, buff);

        if (node.line_g_per_meter != 0) {
            std::ostringstream line_g_stream;
            line_g_stream << 1.0 / (node.line_g_per_meter * node.step / 2);
            std::string line_g_str = line_g_stream.str();
            buff = trim(trim("GL" + node.name) + " " + trim(node.name) + " 0 " + trim(line_g_str));
            appendToStringArray(res, buff);
        }

    }

    nw_node_t preprocess_t::addNodeWithId(const terminal_node_t& node) const {
        int stat;
        int d;
        nw_node_t res;
        std::string sConductor;
        int conductor_number;

        this->conductors_before_cable.get(key(node.belongs_to_cable.name), conductor_number);
        conductor_number = conductor_number + node.conductor_in_cable;

        this->cable_name_to_bundle_id.get(key(node.belongs_to_cable.name), d, stat);

        if (stat != 0) return res;
        
        std::ostringstream sConductorStream;
        sConductorStream << node.conductor_in_cable;
        sConductor = sConductorStream.str();
        
        res.name = trim(node.belongs_to_cable.name) + "_" + trim(sConductor) + "_" + nodeSideToString(node.side);
        res.v = 0.0;
        res.i = 0.0;
        res.bundle_number = d;
        res.conductor_number = conductor_number;

        {
            int v_index, i_index;
            double line_c_per_meter, line_g_per_meter, step;
            if (node.side == TERMINAL_NODE_SIDE_INI) {
                v_index = this->bundles[d].v.second.min();
                i_index = this->bundles[d].i.second.min();
                line_c_per_meter = this->bundles[d].cpul[this->bundles[d].cpul.first.min(), conductor_number, conductor_number];
                line_g_per_meter = this->bundles[d].gpul[this->bundles[d].gpul.first.min(), conductor_number, conductor_number];
                step = this->bundles[d].du[this->bundles[d].du.first.min(), conductor_number, conductor_number];
                res.side = TERMINAL_NODE_SIDE_INI;

            } else if (node.side == TERMINAL_NODE_SIDE_END) {
                v_index = this->bundles[d].v.second.max();
                i_index = this->bundles[d].i.second.max();
                line_c_per_meter = this->bundles[d].cpul[this->bundles[d].cpul.first.max(), conductor_number, conductor_number];
                line_g_per_meter = this->bundles[d].gpul[this->bundles[d].gpul.first.max(), conductor_number, conductor_number];
                step = this->bundles[d].du[this->bundles[d].du.first.max(), conductor_number, conductor_number];
                res.side = TERMINAL_NODE_SIDE_END;
            }
            res.v_index = v_index;
            res.i_index = i_index;
            res.line_c_per_meter = line_c_per_meter;
            res.line_g_per_meter = line_g_per_meter;
            res.step = step;
        }

        if (node.termination.termination_type == TERMINATION_open) res.open = true;
        res.source = node.termination.source;

        return res;
    }

    std::string preprocess_t::nodeSideToString(int side) const {
        if (side == TERMINAL_NODE_SIDE_INI) {
            return "initial";
        } else if (side == TERMINAL_NODE_SIDE_END) {
            return "end";
        }
        return "";
    }

    void preprocess_t::connectNodeToGround(const terminal_connection_t& terminal_connection, std::vector<nw_node_t>& nodes, std::vector<std::string>& description) {
        std::vector<nw_node_t> aux_nodes = nodes;
        nodes.clear();
        nodes.resize(aux_nodes.size() + 1);

        nw_node_t new_node = this->addNodeWithId(terminal_connection.nodes[0]);
        nodes[aux_nodes.size()] = new_node;
        for (size_t i = 0; i < aux_nodes.size(); ++i) {
            nodes[i] = aux_nodes[i];
        }

        std::vector<std::string> node_description = writeNodeDescription(new_node, terminal_connection.nodes[0].termination, "0");
        std::vector<std::string> old_description = description;
        description.clear();
        description.resize(old_description.size() + node_description.size());
        for (size_t i = 0; i < old_description.size(); ++i) {
            description[i] = old_description[i];
        }
        for (size_t i = 0; i < node_description.size(); ++i) {
            description[old_description.size() + i] = node_description[i];
        }
    }

    void preprocess_t::connectNodesToNetworkCircuit(const terminal_connection_t& terminal_connection, std::vector<nw_node_t>& nodes, std::vector<std::string>& description) {
        std::vector<nw_node_t> aux_nodes = nodes;
        nodes.clear();
        nodes.resize(aux_nodes.size() + terminal_connection.nodes.size());

        for (int i = 0; i < terminal_connection.nodes.size(); ++i) {
            nw_node_t new_node = this->addNodeWithId(terminal_connection.nodes[i]);
            nodes[aux_nodes.size() + i] = new_node;

            std::ostringstream str_term_stream;
            str_term_stream << terminal_connection.nodes[i].termination.networkCircuitNode;
            std::string str_term = str_term_stream.str();
            std::string network_circuit_node = trim(terminal_connection.network_circuit.circuit_name) + "_" + trim(str_term);

            std::vector<std::string> node_description = writeNodeDescription(new_node, terminal_connection.nodes[i].termination, trim(network_circuit_node));

            std::vector<std::string> old_description;
            if (description.size() > 0) {
                old_description = description;
            } else {
                old_description.clear();
            }
            
            description.clear();
            description.resize(old_description.size() + node_description.size());
            for (size_t j = 0; j < old_description.size(); ++j) {
                description[j] = old_description[j];
            }
            for (size_t j = 0; j < node_description.size(); ++j) {
                description[old_description.size() + j] = node_description[j];
            }
        }
        for (size_t i = 0; i < aux_nodes.size(); ++i) {
            nodes[i] = aux_nodes[i];
        }
    }

    void preprocess_t::connectNodes(const terminal_connection_t& terminal_connection, std::vector<nw_node_t>& nodes, std::vector<std::string>& description) {
        std::vector<nw_node_t> aux_nodes = nodes;
        nodes.clear();
        nodes.resize(aux_nodes.size() + terminal_connection.nodes.size());

        std::string interior_node = trim(terminal_connection.nodes[0].belongs_to_cable.name) + "_" +
                                    trim(terminal_connection.nodes[1].belongs_to_cable.name) + "_inter";

        for (int i = 0; i < terminal_connection.nodes.size(); ++i) {
            nw_node_t new_node = this->addNodeWithId(terminal_connection.nodes[i]);
            nodes[aux_nodes.size() + i] = new_node;
        }
    }

node_description = writeNodeDescription(new_node, terminal_connection.nodes(i).termination, interior_node);

        if (old_description.size() > 0) { 
            old_description.clear();
            old_description.resize(description.size());
        }
        old_description = description;

        description.clear();
        description.resize(old_description.size() + node_description.size());
        for (size_t k = 0; k < old_description.size(); ++k) {
            description[k] = old_description[k];
        }
        for (size_t k = 0; k < node_description.size(); ++k) {
            description[old_description.size() + k] = node_description[k];
        }
    
        }
        for (size_t i = 0; i < aux_nodes.size(); ++i) {
            nodes[i] = aux_nodes[i];
        }
    }


    network_t buildNetwork(preprocess_t& this_obj, const terminal_network_t& terminal_network) {
        std::vector<nw_node_t> nodes;
        std::vector<std::string> description;
        std::vector<std::string> listOfModels;
        network_t res;
        int i;
        std::vector<terminal_connection_t> network_circuit_connections;
        std::vector<terminal_connection_t> node2node_connections;

        filterConnections(terminal_network.connections, network_circuit_connections, node2node_connections);

        listOfModels.clear();
        description.clear();
        nodes.clear();
        
        for (i = 0; i < static_cast<int>(node2node_connections.size()); ++i) {
            if (node2node_connections[i].nodes.size() == 1) { 
                this_obj.connectNodeToGround(node2node_connections[i], nodes, description);
            } else if (node2node_connections[i].nodes.size() > 1) { 
                this_obj.connectNodes(node2node_connections[i], nodes, description);
            }
        }

        for (i = 0; i < static_cast<int>(network_circuit_connections.size()); ++i) { 
            addCircuitModel(description, network_circuit_connections[i].network_circuit, listOfModels);
            addCircuitInstance(description, network_circuit_connections[i].network_circuit);
        }

        for (i = 0; i < static_cast<int>(network_circuit_connections.size()); ++i) { 
            this_obj.connectNodesToNetworkCircuit(network_circuit_connections[i], nodes, description);
        }

        res = networkCtor(nodes, description);

        return res;
    }

    bool isModelIncluded(const std::vector<std::string>& listOfModels, const std::string& model) {
        if (listOfModels.size() == 0) { 
            return false;
        }
        for (size_t i = 0; i < listOfModels.size(); ++i) {
            if (model == listOfModels[i]) {
                return true;
            }
        }
        return false;
    }

    void addCircuitInstance(std::vector<std::string>& description, const network_circuit_t& network_circuit) {
        std::string buff;

        std::string ports = " ";
        for (int i = 0; i < network_circuit.number_of_nodes; ++i) {
            std::stringstream str_term;
            str_term << (i + 1);
            ports = ports + trim(network_circuit.circuit_name) + "_" + trim(str_term.str()) + " ";
        }

        buff = trim("x" + trim(network_circuit.circuit_name) + " " + trim(ports) + " " + trim(network_circuit.model_name));
        appendToStringArray(description, buff);    

    }

    void addCircuitModel(std::vector<std::string>& description, const network_circuit_t& network_circuit, std::vector<std::string>& listOfModels) {
        std::string buff;

        std::string ports = " ";
        std::stringstream str_term;
        int i;

        buff = trim(network_circuit.model_file);
        if (isModelIncluded(buff, listOfModels)) return;

        appendToStringArray(listOfModels, buff);    

        buff = trim(".include " + network_circuit.model_file);
        appendToStringArray(description, buff);    

    }

    void filterConnections(const std::vector<terminal_connection_t>& all_conn, std::vector<terminal_connection_t>& subckt_conn, std::vector<terminal_connection_t>& node_conn) {
        int i, j, subckt_size, node_size, numberOfNodes, numberOfCktNodes;
        bool is_ckt;

        subckt_size = 0;
        node_size = 0;

        for (i = 0; i < static_cast<int>(all_conn.size()); ++i) {
            if (all_conn[i].network_circuit.number_of_nodes != -1) { 
                subckt_size = subckt_size + 1;
            } else {
                node_size = node_size + 1;
            }
        }


        subckt_conn.resize(subckt_size);
        node_conn.resize(node_size);
        subckt_size = 1;
        node_size = 1;

        is_ckt = true;

        for (i = 0; i < static_cast<int>(all_conn.size()); ++i) {
            if (all_conn[i].network_circuit.number_of_nodes != -1) { 
                subckt_conn[subckt_size - 1] = all_conn[i];
                subckt_size = subckt_size + 1;
            } else { 
                node_conn[node_size - 1] = all_conn[i];
                node_size = node_size + 1;
            }
        }
    }

    void endDescription(std::vector<std::string>& description) {
        std::string buff;

        buff = ".end";
        appendToStringArray(description, buff);
        
        buff = "NULL";
        appendToStringArray(description, buff);
        
    }

    void addNetworksDescription(std::vector<std::string>& description, const std::vector<network_t>& networks) {
        int i, j; 
        std::string buff;
        for (i = 0; i < static_cast<int>(networks.size()); ++i) {
            for (j = 0; j < static_cast<int>(networks[i].description.size()); ++j) {
                buff = networks[i].description[j];
                appendToStringArray(description, buff);
            }
        }
    }

    void addAnalysis(std::vector<std::string>& description, double final_time, double dt, int print_step) {
        std::string buff;
        std::string sTime, sdt, sDelta, sPrint;
        int i;

        std::stringstream ssTime, ssdt, ssDelta, ssPrint;
        ssTime << std::scientific << std::setprecision(2) << final_time;
        ssdt << std::scientific << std::setprecision(2) << dt;
        ssDelta << std::scientific << std::setprecision(2) << dt/200;
        ssPrint << std::scientific << std::setprecision(2) << final_time/print_step;
        
        sTime = ssTime.str();
        sdt = ssdt.str();
        sDelta = ssDelta.str();
        sPrint = ssPrint.str();

        buff = trim(".option reltol = 0.005 gmin=1e-50");
        appendToStringArray(description, buff);       
        buff = trim(".tran " + sdt + " " + sTime + " 0 " + sDelta);
        appendToStringArray(description, buff);       
    }

    void addSavedNodes(std::vector<std::string>& description, const std::vector<network_t>& networks) {
        std::string buff;
        std::string saved_nodes;
        int i, j;
        for (j = 0; j < static_cast<int>(networks.size()); ++j) {
            for (i = 0; i < static_cast<int>(networks[j].nodes.size()); ++i) {
                saved_nodes = ".save  V1" + trim(networks[j].nodes[i].name) + "#branch ";
                saved_nodes = saved_nodes + trim(networks[j].nodes[i].name) + " ";
                buff = trim(saved_nodes);
                appendToStringArray(description, buff);
            }
        }
        

    }


    network_manager_t buildNetworkManager(preprocess_t& this_obj, const std::vector<terminal_network_t>& terminal_networks) {
        std::vector<network_t> networks;
        network_manager_t res;
        std::vector<std::string> description;
        std::string buff;

int i, n;
        std::vector<bool> network_in_MPIslice;

#ifdef CompileWithMPI
        int j, k, d, stat;
#endif
        network_in_MPIslice.resize(terminal_networks.size(), true);
        n = terminal_networks.size();
#ifdef CompileWithMPI
        
        for (i = 0; i < terminal_networks.size(); i++) {
            for (j = 0; j < terminal_networks[i].connections.size(); j++) {
                for (k = 0; k < terminal_networks[i].connections[j].nodes.size(); k++) {
                    if (terminal_networks[i].connections[j].nodes[k].belongs_to_cable != nullptr) {
                        this->cable_name_to_bundle_id.get(
                            key(terminal_networks[i].connections[j].nodes[k].belongs_to_cable->name),
                            d,
                            stat);
                        if (stat != 0) network_in_MPIslice[i] = network_in_MPIslice[i] && false;
                        if (!this->bundles[d].bundle_in_layer) network_in_MPIslice[i] = network_in_MPIslice[i] && false;
                    } else {
                        network_in_MPIslice[i] = network_in_MPIslice[i] && false;
                    }
                }
            }
        }
        n = std::count(network_in_MPIslice.begin(), network_in_MPIslice.end(), true);
#endif
        
        networks.resize(n);
        n = 0;
        for (i = 0; i < network_in_MPIslice.size(); i++) {
            if (network_in_MPIslice[i]) { 
                n++;
                networks[n - 1] = this->buildNetwork(terminal_networks[i]);
            }
        }
        
        description.resize(0);
        buff = "* network description message";
        appendToStringArray(description, buff);
        // description = ["* network description message"]
        addNetworksDescription(description, networks);
        addAnalysis(description, this->final_time, this->dt, 100);
        addSavedNodes(description, networks);
        endDescription(description);        
        res = network_managerCtor(networks, description, this->final_time, this->dt);

    }

    void addGenerators(preprocess_t& this, std::vector<parsed_generator_t>& parsed_generators) {
        int i, d, stat, n;

        for (i = 0; i < parsed_generators.size(); i++) {
            this->cable_name_to_bundle_id.get(key(parsed_generators[i].attached_to_cable->name),
                                               d,
                                               stat);
            if (stat != 0) return;
            this->conductors_before_cable.get(parsed_generators[i].attached_to_cable->name, n);
            this->bundles[d].addGenerator(index = parsed_generators[i].index,
                                              conductor = n + parsed_generators[i].conductor,
                                              gen_type = parsed_generators[i].generator_type,
                                              resistance = parsed_generators[i].resistance,
                                              path = parsed_generators[i].path_to_excitation
#ifdef CompileWithMPI                                               
                                             , layer_indices = this->bundles[d].layer_indices
#endif                        
                                          );

        }
    }


    void addProbesWithId(preprocess_t& this, std::vector<parsed_probe_t>& parsed_probes) {
        int i, d, stat;
        mtl_bundle_t tbundle;
        std::string probe_name;

        for (i = 0; i < parsed_probes.size(); i++) {
            this->cable_name_to_bundle_id.get(key(parsed_probes[i].attached_to_cable->name),
                                               d,
                                               stat);

            if (stat != 0) return;
            probe_name = parsed_probes[i].probe_name + "_" + this->bundles[d].name;
            

            this->bundles[d].addProbe(index = parsed_probes[i].index,
                                          probe_type = parsed_probes[i].probe_type,
                                          name = probe_name,
                                          position = parsed_probes[i].probe_position
#ifdef CompileWithMPI                                               
                                          , layer_indices = this->bundles[d].layer_indices
#endif                        
                                          );
        }
    }


}

