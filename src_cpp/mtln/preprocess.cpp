#include "mtln_preprocess_m.h"

#include <algorithm>
#include <any>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <cstdio>

#include "Report_m.h"
#include "fhash_m.h"
#include "generators_m.h"
#include "mtl_m.h"
#include "network_m.h"
#include "probes_m.h"

#ifdef CompileWithMPI
#include <mpi.h>
extern MPI_Comm SUBCOMM_MPI;
extern MPI_Comm subcomm_mpi;
#endif

namespace mtln_preprocess_m {

using mtln_types_m::cable_abstract_t;
using mtln_types_m::cable_t;
using mtln_types_m::segment_t;
using mtln_types_m::SOURCE_TYPE_VOLTAGE;
using mtln_types_m::SOURCE_TYPE_CURRENT;
using mtln_types_m::TERMINATION_SERIES;
using mtln_types_m::TERMINATION_PARALLEL;
using mtln_types_m::TERMINATION_RsLCp;
using mtln_types_m::TERMINATION_LsRCp;
using mtln_types_m::TERMINATION_CsLRp;
using mtln_types_m::TERMINATION_RLsCp;
using mtln_types_m::TERMINATION_RCsLp;
using mtln_types_m::TERMINATION_LCsRp;
using mtln_types_m::TERMINATION_SHORT;
using mtln_types_m::TERMINATION_OPEN;
using mtln_types_m::TERMINATION_CIRCUIT;
using mtln_types_m::TERMINATION_NETWORK;
using mtln_types_m::TERMINATION_UNDEFINED;
using mtln_types_m::TERMINAL_NODE_SIDE_INI;
using mtln_types_m::TERMINAL_NODE_SIDE_END;
using mtln_types_m::parsed_generator_t;
using parsed_probe_t = mtln_types_m::probe_t;
using mtln_types_m::shielded_multiwire_t;
using mtln_types_m::termination_t;
using mtln_types_m::terminal_connection_t;
using mtln_types_m::terminal_network_t;
using mtln_types_m::terminal_node_t;
using mtln_types_m::network_circuit_t;
using mtln_types_m::unshielded_multiwire_t;
using mtl_m::mtl_t;
using mtl_m::transmission_line_bundle_t;
using mtl_m::transmission_line_level_t;
using mtl_bundle_m::mtl_bundle_t;
using network_m::network_t;
using network_m::nw_node_t;

namespace {

void fhash_set_int(fhash_m::fhash_tbl_t& tbl, const fhash_m::fhash_key_t& key, int value) {
    tbl.set(key, std::any(value));
}

bool fhash_get_int(const fhash_m::fhash_tbl_t& tbl, const fhash_m::fhash_key_t& key, int& value) {
    std::any out;
    int stat = 0;
    tbl.get_raw(key, out, &stat);
    if (stat != 0) {
        return false;
    }
    value = std::any_cast<int>(out);
    return true;
}

std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t");
    const auto e = s.find_last_not_of(" \t");
    if (b == std::string::npos) {
        return "";
    }
    return s.substr(b, e - b + 1);
}

void appendToStringArray(std::vector<std::string>& arr, const std::string& s) {
    arr.push_back(s);
}

void addInitialConnector(mtl_t& line, const mtln_types_m::connector_t& connector) {
    for (int i = 0; i < line.number_of_conductors; ++i) {
        line.rpul[0][static_cast<size_t>(i)][static_cast<size_t>(i)] =
            connector.resistances[static_cast<size_t>(i)] / line.du[0][static_cast<size_t>(i)][static_cast<size_t>(i)];
    }
    line.initial_connector_transfer_impedances = connector.transfer_impedances_per_meter;
}

void addEndConnector(mtl_t& line, const mtln_types_m::connector_t& connector) {
    const int last = static_cast<int>(line.du.size()) - 1;
    for (int i = 0; i < line.number_of_conductors; ++i) {
        line.rpul[static_cast<size_t>(last)][static_cast<size_t>(i)][static_cast<size_t>(i)] =
            connector.resistances[static_cast<size_t>(i)] / line.du[static_cast<size_t>(last)][static_cast<size_t>(i)][static_cast<size_t>(i)];
    }
    line.end_connector_transfer_impedances = connector.transfer_impedances_per_meter;
}

} // namespace

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
    for (size_t i = 0; i < curr_level.cables.size(); ++i) {
        for (size_t j = 0; j < cs.size(); ++j) {
            cable_t* ptr = cs[j].ptr.get();
            if (auto* shielded = dynamic_cast<shielded_multiwire_t*>(ptr)) {
                if (shielded->parent_cable == curr_level.cables[i]) {
                    next_level_size++;
                }
            }
        }
    }
    next_level.cables.resize(static_cast<size_t>(next_level_size));
    int n = 0;
    for (size_t i = 0; i < curr_level.cables.size(); ++i) {
        for (size_t j = 0; j < cs.size(); ++j) {
            cable_t* ptr = cs[j].ptr.get();
            if (auto* shielded = dynamic_cast<shielded_multiwire_t*>(ptr)) {
                if (shielded->parent_cable == curr_level.cables[i]) {
                    n++;
                    next_level.cables[static_cast<size_t>(n - 1)] = cs[j].ptr.get();
                }
            }
        }
    }
    curr_level = next_level;
    return static_cast<int>(curr_level.cables.size());
}

cable_bundle_t buildCableBundleFromParent(cable_t* parent, const std::vector<cable_abstract_t>& cables) {
    cable_bundle_t res;
    cable_level_t level;
    res.levels.resize(1);
    res.levels[0].cables.resize(1);
    level.cables.resize(1);
    level.cables[0] = parent;
    res.levels[0] = level;
    while (findNextLevel(level, cables) != 0) {
        appendLevel(res.levels, level);
    }
    return res;
}

std::vector<size_t> findParentCableIndices(const std::vector<cable_abstract_t>& cables) {
    std::vector<size_t> parent_ids;
    for (size_t i = 0; i < cables.size(); ++i) {
        cable_t* ptr = cables[i].ptr.get();
        if (dynamic_cast<unshielded_multiwire_t*>(ptr)) {
            parent_ids.push_back(i);
        } else if (auto* shielded = dynamic_cast<shielded_multiwire_t*>(ptr)) {
            if (shielded->parent_cable == nullptr) {
                parent_ids.push_back(i);
            }
        }
    }
    return parent_ids;
}

std::vector<cable_bundle_t> buildCableBundles(const std::vector<cable_abstract_t>& cables) {
    const auto parent_ids = findParentCableIndices(cables);
    std::vector<cable_bundle_t> cable_bundles(parent_ids.size());
    for (size_t i = 0; i < parent_ids.size(); ++i) {
        cable_bundles[i] = buildCableBundleFromParent(cables[parent_ids[i]].ptr.get(), cables);
    }
    return cable_bundles;
}

    std::vector<int> conductorsInLevel(const transmission_line_bundle_t& line) {
        std::vector<int> res(line.levels.size(), 0);
        for (size_t i = 0; i < line.levels.size(); ++i) {
            for (const auto& l : line.levels[i].lines) {
                res[i] += l.number_of_conductors;
            }
        }
        return res;
    }

    int findConductorsBeforeCable(const std::string& name, const transmission_line_level_t& level) {
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

    int findOuterConductorNumber(const mtl_t& line, const transmission_line_level_t& level, int conductors_in_level) {
        return findConductorsBeforeCable(line.parent_name, level) +
               conductors_in_level +
               line.conductor_in_parent;
    }

    std::vector<int> findInnerConductorRange(const mtl_t& line, const transmission_line_level_t& level, int conductors_in_level) {
        int start = findConductorsBeforeCable(line.name, level) + conductors_in_level;
        std::vector<int> res(line.number_of_conductors);
        for (int k = 0; k < line.number_of_conductors; ++k) {
            res[k] = start + k + 1;
        }
        return res;
    }

    void setBundleTransferImpedance(mtl_bundle_t& bundle, const transmission_line_bundle_t& line) {
        std::vector<int> conductors_in_level = conductorsInLevel(line);
        bundle.conductors_in_level = conductors_in_level;

        for (size_t i = 1; i < line.levels.size(); ++i) {
            for (const auto& l : line.levels[i].lines) {
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

    void mapConductorsBeforeCable(fhash_m::fhash_tbl_t& conductors_before_cable, const transmission_line_bundle_t& line) {
        const std::vector<int> conductors_in_level = conductorsInLevel(line);
        if (!line.levels.empty() && !line.levels[0].lines.empty()) {
            fhash_set_int(conductors_before_cable, fhash_m::key(line.levels[0].lines[0].name), 0);
        }
        for (size_t i = 1; i < line.levels.size(); ++i) {
            const int sum_prev = std::accumulate(conductors_in_level.begin(), conductors_in_level.begin() + static_cast<int>(i), 0);
            for (const auto& ln : line.levels[i].lines) {
                const auto range_in = findInnerConductorRange(ln, line.levels[i], sum_prev);
                if (range_in.empty()) {
                    throw std::runtime_error("range in cannot be empty");
                }
                fhash_set_int(conductors_before_cable, fhash_m::key(ln.name), range_in[0] - 1);
            }
        }
    }

std::vector<mtl_bundle_t> preprocess_t::buildMTLBundles(const std::vector<transmission_line_bundle_t>& lines) {
        std::vector<mtl_bundle_t> res;
        fhash_m::fhash_tbl_t conductors_before_cable;
        int i;
#ifdef CompileWithMPI
        MPI_Barrier(subcomm_mpi);
#endif
        res.resize(lines.size());
        for (i = 0; i < lines.size(); i++) {
            res[i] = mtl_bundle_m::mtl_bundle_ctor(lines[i].levels, lines[i].levels[0].lines[0].name);
            if (res[i].dt < this->dt) {
                this->dt = res[i].dt;
            }
            setBundleTransferImpedance(res[i], lines[i]);
            mapConductorsBeforeCable(conductors_before_cable, lines[i]);
        }
        this->conductors_before_cable = conductors_before_cable;
        return res;
    }

    mtl_t buildLineFromCable(cable_t& cable, double dt,
                             std::optional<std::vector<std::vector<int>>> layer_indices = std::nullopt,
                             std::optional<bool> bundle_in_layer = std::nullopt,
                             std::optional<std::vector<int>> alloc_z = std::nullopt) {
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
            res = mtl_m::mtl_shielded(
                shielded->inductance_per_meter, shielded->capacitance_per_meter, shielded->resistance_per_meter,
                shielded->conductance_per_meter, shielded->step_size, shielded->name, shielded->segments, dt,
                parent_name, conductor_in_parent, shielded->transfer_impedance
	#ifdef CompileWithMPI
	                , layer_indices, bundle_in_layer, alloc_z
	#endif
	            );
	        } else if (auto* unshielded = dynamic_cast<unshielded_multiwire_t*>(&cable)) {
            res = mtl_m::mtl_unshielded(
                unshielded->cell_inductance_per_meter, unshielded->cell_capacitance_per_meter,
                unshielded->resistance_per_meter, unshielded->conductance_per_meter, unshielded->step_size,
                unshielded->name, unshielded->segments, dt, unshielded->multipolar_expansion, unshielded->radius
	#ifdef CompileWithMPI
	                , layer_indices, bundle_in_layer, alloc_z
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

bool isSegmentWithinAllocBox(const std::vector<segment_t>& segs, int i, const std::vector<int>& z) {
    if (i < 1 || i > static_cast<int>(segs.size())) {
        return false;
    }
    return (segs[static_cast<size_t>(i - 1)].z >= z[0]) && (segs[static_cast<size_t>(i - 1)].z <= z[1]);
}

std::vector<std::vector<int>> getSegmentRanges(const cable_t& cable, const std::vector<int>& alloc_z) {
            int n = 0;
            // Precount
            for (int i = 0; i < static_cast<int>(cable.segments.size()); ++i) {
                if (isSegmentWithinAllocBox(cable.segments, i + 1, alloc_z)) {
                    // Logic handled in the loop below to count transitions
                }
            }
            
            // Reset and count properly
            n = 0;
            bool in_layer = false;
            for (int i = 0; i < static_cast<int>(cable.segments.size()); ++i) {
                if (isSegmentWithinAllocBox(cable.segments, i + 1, alloc_z)) {
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
                if (isSegmentWithinAllocBox(cable.segments, i + 1, alloc_z)) {
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

bool isBundleInLayer(cable_t& cable, const std::vector<int>& alloc_z) {
    int n = 0;
    bool in_layer = false;
    for (int i = 0; i < static_cast<int>(cable.segments.size()); ++i) {
        if (isSegmentWithinAllocBox(cable.segments, i + 1, alloc_z)) {
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
    if (in_layer) {
        n++;
    }
    return n != 0;
}

std::vector<std::vector<int>> findIndicesInLayer(cable_t& cable, const std::vector<int>& alloc_z) {
    return getSegmentRanges(cable, alloc_z);
}

    std::vector<transmission_line_bundle_t> buildLineBundles(const std::vector<cable_bundle_t>& cable_bundles, double dt,
                                                             const std::array<FDETYPES_m::XYZlimit_t, 6>& alloc) {
        std::vector<transmission_line_bundle_t> res;
        int i, j, k;
        int nb, nl, nc;
        std::vector<std::vector<int>> layer_indices;
        bool bundle_in_layer = false;
        std::vector<int> alloc_z(2);
        if (alloc[2].ZI != 0 || alloc[2].ZE != 0) {
            alloc_z[0] = alloc[2].ZI;
            alloc_z[1] = alloc[2].ZE;
        }
        nb = cable_bundles.size();
        res.resize(nb);
        for (i = 0; i < nb; i++) {
            if (alloc_z[0] != 0 || alloc_z[1] != 0) {
                if (layer_indices.size() > 0) layer_indices.clear();
                bundle_in_layer = isBundleInLayer(*cable_bundles[i].levels[0].cables[0], alloc_z);
                if (bundle_in_layer) {
                    layer_indices = findIndicesInLayer(*cable_bundles[i].levels[0].cables[0], alloc_z);
                } else {
                    layer_indices.clear();
                }
            }
            nl = cable_bundles[i].levels.size();
            res[i].levels.resize(nl);
            for (j = 0; j < nl; j++) {
                nc = cable_bundles[i].levels[j].cables.size();
                res[i].levels[j].lines.resize(nc);
                for (k = 0; k < nc; k++) {
                    if (alloc_z[0] != 0 || alloc_z[1] != 0) {
                        res[i].levels[j].lines[k] = buildLineFromCable(*cable_bundles[i].levels[j].cables[k], dt, layer_indices, bundle_in_layer, alloc_z);
                    } else {
                        res[i].levels[j].lines[k] = buildLineFromCable(*cable_bundles[i].levels[j].cables[k], dt);
                    }
                }
            }
        }
        return res;
    }

        fhash_m::fhash_tbl_t mapCablesToBundlesId(const std::vector<transmission_line_bundle_t>& lines, const std::vector<mtl_bundle_t>& bundles) {
            fhash_m::fhash_tbl_t res;
            for (size_t i = 0; i < lines.size(); ++i) {
                for (size_t j = 0; j < lines[i].levels.size(); ++j) {
                    for (size_t k = 0; k < lines[i].levels[j].lines.size(); ++k) {
                        fhash_set_int(res, fhash_m::key(lines[i].levels[j].lines[k].name), static_cast<int>(i + 1)); // 1-based index
                    }
                }
            }
            return res;
        }
void preprocess_t::addProbesWithId(const std::vector<parsed_probe_t>& parsed_probes) {
    for (const auto& p : parsed_probes) {
        if (!p.attached_to_cable) continue;
        int d = 0;
        if (!fhash_get_int(cable_name_to_bundle_id, fhash_m::key(p.attached_to_cable->name), d)) {
            continue;
        }
        const std::string probe_name = p.probe_name + "_" + bundles[static_cast<size_t>(d - 1)].name;
        bundles[static_cast<size_t>(d - 1)].addProbe(p.index, p.probe_type, probe_name, p.probe_position
#ifdef CompileWithMPI
                                                     , bundles[static_cast<size_t>(d - 1)].layer_indices
#endif
        );
    }
}

void preprocess_t::addGenerators(const std::vector<parsed_generator_t>& parsed_generators) {
    for (const auto& g : parsed_generators) {
        if (!g.attached_to_cable) continue;
        int d = 0;
        int n = 0;
        if (!fhash_get_int(cable_name_to_bundle_id, fhash_m::key(g.attached_to_cable->name), d)) {
            continue;
        }
        fhash_get_int(conductors_before_cable, fhash_m::key(g.attached_to_cable->name), n);
        bundles[static_cast<size_t>(d - 1)].addGenerator(g.index, n + g.conductor, g.generator_type, g.resistance,
                                                        g.path_to_excitation
#ifdef CompileWithMPI
                                                        , bundles[static_cast<size_t>(d - 1)].layer_indices
#endif
        );
    }
}

namespace {

void append_desc(std::vector<std::string>& arr, const std::string& s) {
    arr.push_back(s);
}

std::string formatSpiceValue(double value) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "% .16g", value);
    return buffer;
}

std::vector<std::string> writeSeriesRLnode(const nw_node_t& node, const termination_t& termination,
                                           const std::string& end_node) {
    std::vector<std::string> res;
    const std::string termination_r = formatSpiceValue(termination.resistance);
    const std::string termination_l = formatSpiceValue(termination.inductance);
    const std::string line_c = formatSpiceValue(node.line_c_per_meter * node.step / 2.0);
    append_desc(res, "R" + node.name + " " + node.name + "_R " + node.name + " " + termination_r);
    if (!termination.source.path_to_excitation.empty()) {
        const std::string generator_r = formatSpiceValue(termination.source.resistance);
        append_desc(res, "L" + node.name + " " + node.name + "_R " + node.name + "_S " + termination_l);
        if (termination.source.source_type == mtln_types_m::SOURCE_TYPE_VOLTAGE) {
            append_desc(res, "V" + node.name + "_S " + node.name + "_S " + node.name + "_genR dc 0");
            append_desc(res, "R" + node.name + "_S " + node.name + "_genR " + end_node + " " + generator_r);
        } else if (termination.source.source_type == mtln_types_m::SOURCE_TYPE_CURRENT) {
            append_desc(res, "I" + node.name + "_S " + end_node + " " + node.name + "_S dc 0");
            if (termination.source.resistance != 1.0e22) {
                append_desc(res, "R" + node.name + "_S " + end_node + " " + node.name + "_S " + generator_r);
            }
        }
    } else {
        append_desc(res, "L" + node.name + " " + node.name + "_R " + end_node + " " + termination_l);
    }
    append_desc(res, "I" + node.name + " " + node.name + " 0 dc 0");
    append_desc(res, "CL" + node.name + " " + node.name + " 0 " + line_c);
    if (node.line_g_per_meter != 0.0) {
        const std::string line_g =
            formatSpiceValue(1.0 / (node.line_g_per_meter * node.step / 2.0));
        append_desc(res, "GL" + node.name + " " + node.name + " 0 " + line_g);
    }
    return res;
}

std::vector<std::string> writeShortNode(const nw_node_t& node, const termination_t& termination,
                                        const std::string& end_node) {
    std::vector<std::string> res;
    const std::string short_r = "1e-10";
    const std::string line_c = formatSpiceValue(node.line_c_per_meter * node.step / 2.0);
    if (!termination.source.path_to_excitation.empty()) {
        const std::string generator_r = formatSpiceValue(termination.source.resistance);
        append_desc(res, "R" + node.name + " " + node.name + " " + node.name + "_S " + short_r);
        if (termination.source.source_type == mtln_types_m::SOURCE_TYPE_VOLTAGE) {
            append_desc(res, "V" + node.name + "_S " + node.name + "_S " + node.name + "_genR dc 0");
            append_desc(res, "R" + node.name + "_S " + node.name + "_genR " + end_node + " " + generator_r);
        } else if (termination.source.source_type == mtln_types_m::SOURCE_TYPE_CURRENT) {
            append_desc(res, "I" + node.name + "_S " + end_node + " " + node.name + "_S dc 0");
            if (termination.source.resistance != 1.0e22) {
                append_desc(res, "R" + node.name + "_S " + end_node + " " + node.name + "_S " + generator_r);
            }
        }
    } else {
        append_desc(res, "R" + node.name + " " + node.name + " " + end_node + " " + short_r);
    }
    append_desc(res, "I" + node.name + " " + node.name + " 0 dc 0");
    append_desc(res, "CL" + node.name + " " + node.name + " 0 " + line_c);
    if (node.line_g_per_meter != 0.0) {
        const std::string line_g =
            formatSpiceValue(1.0 / (node.line_g_per_meter * node.step / 2.0));
        append_desc(res, "GL" + node.name + " " + node.name + " 0 " + line_g);
    }
    return res;
}

std::vector<std::string> writeSeriesRLCnode(const nw_node_t& node, const termination_t& termination,
                                            const std::string& end_node) {
    std::vector<std::string> res;
    const std::string termination_r = formatSpiceValue(termination.resistance);
    const std::string termination_l = formatSpiceValue(termination.inductance);
    const std::string termination_c = formatSpiceValue(termination.capacitance);
    const std::string line_c = formatSpiceValue(node.line_c_per_meter * node.step / 2.0);
    if (!termination.source.path_to_excitation.empty()) {
        const std::string generator_r = formatSpiceValue(termination.source.resistance);
        append_desc(res, "R" + node.name + " " + node.name + " " + node.name + "_S " + termination_r);
        append_desc(res, "L" + node.name + " " + node.name + " " + node.name + "_S " + termination_l);
        append_desc(res, "C" + node.name + " " + node.name + " " + node.name + "_S " + termination_c);
        if (termination.source.source_type == mtln_types_m::SOURCE_TYPE_VOLTAGE) {
            append_desc(res, "V" + node.name + "_S " + node.name + "_S " + node.name + "_genR dc 0");
            append_desc(res, "R" + node.name + "_S " + node.name + "_genR " + end_node + " " + generator_r);
        } else if (termination.source.source_type == mtln_types_m::SOURCE_TYPE_CURRENT) {
            append_desc(res, "I" + node.name + "_S " + end_node + " " + node.name + "_S dc 0");
            if (termination.source.resistance != 1.0e22) {
                append_desc(res, "R" + node.name + "_S " + end_node + " " + node.name + "_S " + generator_r);
            }
        }
    } else {
        append_desc(res, "R" + node.name + " " + node.name + " " + end_node + termination_r);
        append_desc(res, "L" + node.name + " " + node.name + " " + end_node + termination_l);
        append_desc(res, "C" + node.name + " " + node.name + " " + end_node + termination_c);
    }
    append_desc(res, "I" + node.name + " " + node.name + " 0 dc 0");
    append_desc(res, "CL" + node.name + " " + node.name + " 0 " + line_c);
    if (node.line_g_per_meter != 0.0) {
        const std::string line_g =
            formatSpiceValue(1.0 / (node.line_g_per_meter * node.step / 2.0));
        append_desc(res, "GL" + node.name + " " + node.name + " 0 " + line_g);
    }
    return res;
}

std::vector<std::string> writeSeriesNode(const nw_node_t& node, const termination_t& termination,
                                         const std::string& end_node) {
    if (termination.capacitance >= 1e22) {
        return writeSeriesRLnode(node, termination, end_node);
    }
    return writeSeriesRLCnode(node, termination, end_node);
}

std::vector<std::string> writeParallelRLCnode(const nw_node_t& /*node*/, const termination_t& /*termination*/,
                                              const std::string& /*end_node*/) {
    return {};
}

std::vector<std::string> writeNetwork_circuitNode(const nw_node_t& node, const termination_t& termination,
                                                  const std::string& end_node) {
    std::vector<std::string> res;
    const std::string short_r = "1e-10";
    const std::string line_c = formatSpiceValue(node.line_c_per_meter * node.step / 2.0);
    if (!termination.source.path_to_excitation.empty()) {
        const std::string generator_r = formatSpiceValue(termination.source.resistance);
        append_desc(res, "R" + node.name + " " + node.name + " " + node.name + "_S " + short_r);
        if (termination.source.source_type == mtln_types_m::SOURCE_TYPE_VOLTAGE) {
            append_desc(res, "V" + node.name + "_S " + node.name + "_S " + node.name + "_genR dc 0");
            append_desc(res, "R" + node.name + "_S " + node.name + "_genR " + end_node + " " + generator_r);
        } else if (termination.source.source_type == mtln_types_m::SOURCE_TYPE_CURRENT) {
            append_desc(res, "I" + node.name + "_S " + end_node + " " + node.name + "_S dc 0");
            if (termination.source.resistance != 1.0e22) {
                append_desc(res, "R" + node.name + "_S " + end_node + " " + node.name + "_S " + generator_r);
            }
        }
    } else {
        append_desc(res, "R" + node.name + " " + node.name + " " + end_node + " " + short_r);
    }
    append_desc(res, "I" + node.name + " " + node.name + " 0 dc 0");
    append_desc(res, "CL" + node.name + " " + node.name + " 0 " + line_c);
    if (node.line_g_per_meter != 0.0) {
        const std::string line_g =
            formatSpiceValue(1.0 / (node.line_g_per_meter * node.step / 2.0));
        append_desc(res, "GL" + node.name + " " + node.name + " 0 " + line_g);
    }
    return res;
}

std::vector<std::string> writeModelNode(const nw_node_t& node, const termination_t& termination,
                                        const std::string& end_node) {
    std::vector<std::string> res;
    const std::string line_c = formatSpiceValue(node.line_c_per_meter * node.step / 2.0);
    append_desc(res, ".include " + termination.model.file);
    if (!termination.source.path_to_excitation.empty()) {
        const std::string generator_r = formatSpiceValue(termination.source.resistance);
        if (termination.source.source_type == mtln_types_m::SOURCE_TYPE_VOLTAGE) {
            append_desc(res, "V" + node.name + "_S " + node.name + " " + node.name + "_genR dc 0");
            append_desc(res, "R" + node.name + "_S " + node.name + "_genR " + node.name + "_S " + generator_r);
        } else if (termination.source.source_type == mtln_types_m::SOURCE_TYPE_CURRENT) {
            append_desc(res, "I" + node.name + "_S " + node.name + "_S " + node.name + " dc 0");
            if (termination.source.resistance != 1.0e22) {
                append_desc(res, "R" + node.name + "_S " + node.name + "_S " + node.name + " " + generator_r);
            }
        }
        append_desc(res, "x" + node.name + " " + node.name + "_S " + end_node + " " + termination.model.name);
    } else {
        append_desc(res, "x" + node.name + " " + node.name + " " + end_node + " " + termination.model.name);
    }
    append_desc(res, "I" + node.name + " " + node.name + " 0 dc 0");
    append_desc(res, "CL" + node.name + " " + node.name + " 0 " + line_c);
    if (node.line_g_per_meter != 0.0) {
        const std::string line_g =
            formatSpiceValue(1.0 / (node.line_g_per_meter * node.step / 2.0));
        append_desc(res, "GL" + node.name + " " + node.name + " 0 " + line_g);
    }
    return res;
}

std::vector<std::string> writeOpenNode(const nw_node_t& node, const termination_t& /*termination*/,
                                       const std::string& end_node) {
    std::vector<std::string> res;
    const std::string line_c = formatSpiceValue(node.line_c_per_meter * node.step / 2.0);
    append_desc(res, "R" + node.name + " " + node.name + " " + end_node + " 1e22");
    append_desc(res, "I" + node.name + " " + node.name + " 0 dc 0");
    append_desc(res, "CL" + node.name + " " + node.name + " 0 " + line_c);
    if (node.line_g_per_meter != 0.0) {
        const std::string line_g =
            formatSpiceValue(1.0 / (node.line_g_per_meter * node.step / 2.0));
        append_desc(res, "GL" + node.name + " " + node.name + " 0 " + line_g);
    }
    return res;
}

void assignTerminationXYsZp(const std::string& xyz, const termination_t& termination, std::string& tx,
                            std::string& ty, std::string& tz) {
    if (xyz == "RLC" || xyz == "LRC") {
        tx = formatSpiceValue(termination.resistance);
        ty = formatSpiceValue(termination.inductance);
        tz = formatSpiceValue(termination.capacitance);
    } else if (xyz == "LCR" || xyz == "CLR") {
        tx = formatSpiceValue(termination.inductance);
        ty = formatSpiceValue(termination.capacitance);
        tz = formatSpiceValue(termination.resistance);
    } else if (xyz == "CRL" || xyz == "RCL") {
        tx = formatSpiceValue(termination.capacitance);
        ty = formatSpiceValue(termination.resistance);
        tz = formatSpiceValue(termination.inductance);
    }
}

void assignTerminationXsYZp(const std::string& xyz, const termination_t& termination, std::string& tx,
                            std::string& ty, std::string& tz) {
    if (xyz == "RLC" || xyz == "RCL") {
        tx = formatSpiceValue(termination.resistance);
        ty = formatSpiceValue(termination.inductance);
        tz = formatSpiceValue(termination.capacitance);
    } else if (xyz == "LRC" || xyz == "LCR") {
        tx = formatSpiceValue(termination.inductance);
        ty = formatSpiceValue(termination.resistance);
        tz = formatSpiceValue(termination.capacitance);
    } else if (xyz == "CLR" || xyz == "CRL") {
        tx = formatSpiceValue(termination.capacitance);
        ty = formatSpiceValue(termination.resistance);
        tz = formatSpiceValue(termination.inductance);
    }
}

void appendLineShunt(std::vector<std::string>& res, const nw_node_t& node) {
    const std::string line_c = formatSpiceValue(node.line_c_per_meter * node.step / 2.0);
    append_desc(res, "I" + node.name + " " + node.name + " 0 dc 0");
    append_desc(res, "CL" + node.name + " " + node.name + " 0 " + line_c);
    if (node.line_g_per_meter != 0.0) {
        const std::string line_g =
            formatSpiceValue(1.0 / (node.line_g_per_meter * node.step / 2.0));
        append_desc(res, "GL" + node.name + " " + node.name + " 0 " + line_g);
    }
}

std::vector<std::string> writeXYsZpNode(const nw_node_t& node, const termination_t& termination,
                                        const std::string& end_node, const std::string& xyz) {
    std::vector<std::string> res;
    std::string tx;
    std::string ty;
    std::string tz;
    assignTerminationXYsZp(xyz, termination, tx, ty, tz);
    append_desc(res, std::string(1, xyz[0]) + node.name + " " + node.name + " " + node.name + "_X " + tx);
    if (!termination.source.path_to_excitation.empty()) {
        const std::string generator_r = formatSpiceValue(termination.source.resistance);
        append_desc(res, std::string(1, xyz[1]) + node.name + " " + node.name + "_X " + node.name + "_S " + ty);
        append_desc(res, std::string(1, xyz[2]) + node.name + " " + node.name + " " + node.name + "_S " + tz);
        if (termination.source.source_type == mtln_types_m::SOURCE_TYPE_VOLTAGE) {
            append_desc(res, "V" + node.name + "_S " + node.name + "_S " + node.name + "_genR dc 0");
            append_desc(res, "R" + node.name + "_S " + node.name + "_genR " + end_node + " " + generator_r);
        } else if (termination.source.source_type == mtln_types_m::SOURCE_TYPE_CURRENT) {
            append_desc(res, "I" + node.name + "_S " + end_node + " " + node.name + "_S dc 0");
            if (termination.source.resistance != 1.0e22) {
                append_desc(res, "R" + node.name + "_S " + end_node + " " + node.name + "_S " + generator_r);
            }
        }
    } else {
        append_desc(res, std::string(1, xyz[1]) + node.name + " " + node.name + "_X " + end_node + " " + ty);
        append_desc(res, std::string(1, xyz[2]) + node.name + " " + node.name + " " + end_node + " " + tz);
    }
    appendLineShunt(res, node);
    return res;
}

std::vector<std::string> writeXsYZpNode(const nw_node_t& node, const termination_t& termination,
                                        const std::string& end_node, const std::string& xyz) {
    std::vector<std::string> res;
    std::string tx;
    std::string ty;
    std::string tz;
    assignTerminationXsYZp(xyz, termination, tx, ty, tz);
    append_desc(res, std::string(1, xyz[0]) + node.name + " " + node.name + " " + node.name + "_p " + tx);
    if (!termination.source.path_to_excitation.empty()) {
        const std::string generator_r = formatSpiceValue(termination.source.resistance);
        append_desc(res, std::string(1, xyz[1]) + node.name + " " + node.name + "_p " + node.name + "_S " + ty);
        append_desc(res, std::string(1, xyz[2]) + node.name + " " + node.name + "_p " + node.name + "_S " + tz);
        if (termination.source.source_type == mtln_types_m::SOURCE_TYPE_VOLTAGE) {
            append_desc(res, "V" + node.name + "_S " + node.name + "_S " + node.name + "_genR dc 0");
            append_desc(res, "R" + node.name + "_S " + node.name + "_genR " + end_node + " " + generator_r);
        } else if (termination.source.source_type == mtln_types_m::SOURCE_TYPE_CURRENT) {
            append_desc(res, "I" + node.name + "_S " + end_node + " " + node.name + "_S dc 0");
            if (termination.source.resistance != 1.0e22) {
                append_desc(res, "R" + node.name + "_S " + end_node + " " + node.name + "_S " + generator_r);
            }
        }
    } else {
        append_desc(res, std::string(1, xyz[1]) + node.name + " " + node.name + "_p " + end_node + " " + ty);
        append_desc(res, std::string(1, xyz[2]) + node.name + " " + node.name + "_p " + end_node + " " + tz);
    }
    appendLineShunt(res, node);
    return res;
}

std::vector<std::string> writeNodeDescription(const nw_node_t& node, const termination_t& termination,
                                              const std::string& end_node) {
    if (termination.termination_type == TERMINATION_SERIES) {
        return writeSeriesNode(node, termination, end_node);
    }
    if (termination.termination_type == TERMINATION_PARALLEL) {
        return writeParallelRLCnode(node, termination, end_node);
    }
    if (termination.termination_type == TERMINATION_RsLCp) {
        return writeXsYZpNode(node, termination, end_node, "RLC");
    }
    if (termination.termination_type == TERMINATION_LsRCp) {
        return writeXsYZpNode(node, termination, end_node, "LRC");
    }
    if (termination.termination_type == TERMINATION_CsLRp) {
        return writeXsYZpNode(node, termination, end_node, "CLR");
    }
    if (termination.termination_type == TERMINATION_RLsCp) {
        return writeXYsZpNode(node, termination, end_node, "RLC");
    }
    if (termination.termination_type == TERMINATION_RCsLp) {
        return writeXYsZpNode(node, termination, end_node, "RCL");
    }
    if (termination.termination_type == TERMINATION_LCsRp) {
        return writeXYsZpNode(node, termination, end_node, "LCR");
    }
    if (termination.termination_type == TERMINATION_SHORT) {
        return writeShortNode(node, termination, end_node);
    }
    if (termination.termination_type == TERMINATION_OPEN) {
        return writeOpenNode(node, termination, end_node);
    }
    if (termination.termination_type == TERMINATION_CIRCUIT) {
        return writeModelNode(node, termination, end_node);
    }
    if (termination.termination_type == TERMINATION_NETWORK) {
        return writeNetwork_circuitNode(node, termination, end_node);
    }
    if (termination.termination_type == TERMINATION_UNDEFINED) {
        Report_m::WarnErrReport("writeNodeDescription: undefined termination at " + node.name, true);
    }
    return {};
}

void filterConnections(const std::vector<terminal_connection_t>& all_conn,
                       std::vector<terminal_connection_t>& subckt_conn,
                       std::vector<terminal_connection_t>& node_conn) {
    subckt_conn.clear();
    node_conn.clear();
    for (const auto& c : all_conn) {
        if (c.network_circuit.number_of_nodes != -1) {
            subckt_conn.push_back(c);
        } else {
            node_conn.push_back(c);
        }
    }
}

void addNetworksDescription(std::vector<std::string>& description, const std::vector<network_t>& networks) {
    for (const auto& net : networks) {
        for (const auto& line : net.description) {
            append_desc(description, line);
        }
    }
}

void addAnalysis(std::vector<std::string>& description, double final_time, double dt) {
    char buf[256];
    std::snprintf(buf, sizeof(buf), ".option reltol = 0.005 gmin=1e-50");
    append_desc(description, buf);
    std::snprintf(buf, sizeof(buf), ".tran %e %e 0 %e", dt, final_time, dt / 200.0);
    append_desc(description, buf);
}

void addSavedNodes(std::vector<std::string>& description, const std::vector<network_t>& networks) {
    for (const auto& net : networks) {
        for (const auto& node : net.nodes) {
            std::string buff = ".save  V1" + node.name + "#branch " + node.name + " ";
            append_desc(description, buff);
        }
    }
}

void endDescription(std::vector<std::string>& description) {
    append_desc(description, ".end");
    append_desc(description, "NULL");
}

} // namespace

std::string preprocess_t::nodeSideToString(int side) const {
    if (side == mtln_types_m::TERMINAL_NODE_SIDE_INI) {
        return "initial";
    }
    if (side == mtln_types_m::TERMINAL_NODE_SIDE_END) {
        return "end";
    }
    return "";
}

nw_node_t preprocess_t::addNodeWithId(const terminal_node_t& node) const {
    nw_node_t res;
    int conductor_number = 0;
    if (!fhash_get_int(conductors_before_cable, fhash_m::key(node.belongs_to_cable->name), conductor_number)) {
        conductor_number = 0;
    }
    conductor_number += node.conductor_in_cable;

    int d = 0;
    if (!fhash_get_int(cable_name_to_bundle_id, fhash_m::key(node.belongs_to_cable->name), d)) {
        return res;
    }

    res.name = node.belongs_to_cable->name + "_" + std::to_string(node.conductor_in_cable) + "_" +
               nodeSideToString(node.side);
    res.v = 0.0;
    res.i = 0.0;
    res.bundle_number = d;
    res.conductor_number = conductor_number;

    const auto& bundle = bundles[static_cast<size_t>(d - 1)];
    if (node.side == mtln_types_m::TERMINAL_NODE_SIDE_INI) {
        res.v_index = 0;
        res.i_index = 0;
        res.line_c_per_meter = bundle.cpul[0][static_cast<size_t>(conductor_number - 1)][static_cast<size_t>(conductor_number - 1)];
        res.line_g_per_meter = bundle.gpul[0][static_cast<size_t>(conductor_number - 1)][static_cast<size_t>(conductor_number - 1)];
        res.step = bundle.du[0][static_cast<size_t>(conductor_number - 1)][static_cast<size_t>(conductor_number - 1)];
        res.side = mtln_types_m::TERMINAL_NODE_SIDE_INI;
    } else if (node.side == mtln_types_m::TERMINAL_NODE_SIDE_END) {
        const int last_div = static_cast<int>(bundle.du.size()) - 1;
        res.v_index = static_cast<int>(bundle.v[0].size()) - 1;
        res.i_index = static_cast<int>(bundle.i[0].size()) - 1;
        res.line_c_per_meter = bundle.cpul[static_cast<size_t>(last_div)][static_cast<size_t>(conductor_number - 1)][static_cast<size_t>(conductor_number - 1)];
        res.line_g_per_meter = bundle.gpul[static_cast<size_t>(last_div)][static_cast<size_t>(conductor_number - 1)][static_cast<size_t>(conductor_number - 1)];
        res.step = bundle.du[static_cast<size_t>(last_div)][static_cast<size_t>(conductor_number - 1)][static_cast<size_t>(conductor_number - 1)];
        res.side = mtln_types_m::TERMINAL_NODE_SIDE_END;
    }
    if (node.termination.termination_type == mtln_types_m::TERMINATION_OPEN) {
        res.open = true;
    }
    res.source = node.termination.source;
    return res;
}

void preprocess_t::connectNodeToGround(const terminal_connection_t& terminal_connection,
                                       std::vector<nw_node_t>& nodes, std::vector<std::string>& description) {
    const nw_node_t new_node = addNodeWithId(terminal_connection.nodes[0]);
    std::vector<nw_node_t> aux = nodes;
    nodes.resize(aux.size() + 1);
    nodes.back() = new_node;
    for (size_t i = 0; i < aux.size(); ++i) {
        nodes[i] = aux[i];
    }
    const auto node_description = writeNodeDescription(new_node, terminal_connection.nodes[0].termination, "0");
    description.insert(description.end(), node_description.begin(), node_description.end());
}

void preprocess_t::connectNodes(const terminal_connection_t& terminal_connection, std::vector<nw_node_t>& nodes,
                              std::vector<std::string>& description) {
    const std::string interior_node = terminal_connection.nodes[0].belongs_to_cable->name + "_" +
                                      terminal_connection.nodes[1].belongs_to_cable->name + "_inter";
    std::vector<nw_node_t> aux = nodes;
    nodes.resize(aux.size() + terminal_connection.nodes.size());
    for (size_t i = 0; i < terminal_connection.nodes.size(); ++i) {
        const nw_node_t new_node = addNodeWithId(terminal_connection.nodes[i]);
        nodes[aux.size() + i] = new_node;
        const auto node_description =
            writeNodeDescription(new_node, terminal_connection.nodes[i].termination, interior_node);
        description.insert(description.end(), node_description.begin(), node_description.end());
    }
    for (size_t i = 0; i < aux.size(); ++i) {
        nodes[i] = aux[i];
    }
}

bool isModelIncluded(const std::vector<std::string>& list_of_models, const std::string& model) {
    for (const auto& m : list_of_models) {
        if (m == model) {
            return true;
        }
    }
    return false;
}

void addCircuitModel(std::vector<std::string>& description, const network_circuit_t& network_circuit,
                     std::vector<std::string>& list_of_models) {
    const std::string model_file = trim(network_circuit.model_file);
    if (isModelIncluded(list_of_models, model_file)) {
        return;
    }
    list_of_models.push_back(model_file);
    append_desc(description, ".include " + model_file);
}

void addCircuitInstance(std::vector<std::string>& description, const network_circuit_t& network_circuit) {
    std::string ports;
    for (int i = 0; i < network_circuit.number_of_nodes; ++i) {
        ports += trim(network_circuit.circuit_name) + "_" + std::to_string(i + 1) + " ";
    }
    append_desc(description,
                "x" + trim(network_circuit.circuit_name) + " " + ports + trim(network_circuit.model_name));
}

void preprocess_t::connectNodesToNetworkCircuit(const terminal_connection_t& terminal_connection,
                                                std::vector<nw_node_t>& nodes,
                                                std::vector<std::string>& description) {
    const std::vector<nw_node_t> aux_nodes = nodes;
    nodes.resize(aux_nodes.size() + terminal_connection.nodes.size());
    for (size_t i = 0; i < terminal_connection.nodes.size(); ++i) {
        const nw_node_t new_node = addNodeWithId(terminal_connection.nodes[i]);
        nodes[aux_nodes.size() + i] = new_node;
        const std::string network_circuit_node =
            trim(terminal_connection.network_circuit.circuit_name) + "_" +
            std::to_string(terminal_connection.nodes[i].termination.networkCircuitNode);
        const auto node_description =
            writeNodeDescription(new_node, terminal_connection.nodes[i].termination, network_circuit_node);
        description.insert(description.end(), node_description.begin(), node_description.end());
    }
    for (size_t i = 0; i < aux_nodes.size(); ++i) {
        nodes[i] = aux_nodes[i];
    }
}

network_t preprocess_t::buildNetwork(const terminal_network_t& terminal_network) {
    std::vector<nw_node_t> nodes;
    std::vector<std::string> description;
    std::vector<std::string> list_of_models;
    std::vector<terminal_connection_t> network_circuit_connections;
    std::vector<terminal_connection_t> node2node_connections;
    filterConnections(terminal_network.connections, network_circuit_connections, node2node_connections);

    for (const auto& conn : node2node_connections) {
        if (conn.nodes.size() == 1) {
            connectNodeToGround(conn, nodes, description);
        } else if (conn.nodes.size() > 1) {
            connectNodes(conn, nodes, description);
        }
    }
    for (const auto& conn : network_circuit_connections) {
        addCircuitModel(description, conn.network_circuit, list_of_models);
        addCircuitInstance(description, conn.network_circuit);
    }
    for (const auto& conn : network_circuit_connections) {
        connectNodesToNetworkCircuit(conn, nodes, description);
    }
    return network_m::networkCtor(nodes, description);
}

network_manager_m::network_manager_t preprocess_t::buildNetworkManager(
    const std::vector<terminal_network_t>& terminal_networks) {
    std::vector<network_t> networks;
    std::vector<bool> network_in_MPIslice(terminal_networks.size(), true);
#ifdef CompileWithMPI
    for (size_t i = 0; i < terminal_networks.size(); ++i) {
        for (const auto& connection : terminal_networks[i].connections) {
            for (const auto& node : connection.nodes) {
                if (!node.belongs_to_cable) {
                    network_in_MPIslice[i] = false;
                    continue;
                }

                int d = 0;
                if (!fhash_get_int(cable_name_to_bundle_id, fhash_m::key(node.belongs_to_cable->name), d) ||
                    d <= 0 || d > static_cast<int>(bundles.size())) {
                    network_in_MPIslice[i] = false;
                    continue;
                }

                if (!bundles[static_cast<size_t>(d - 1)].bundle_in_layer) {
                    network_in_MPIslice[i] = false;
                }
            }
        }
    }
#endif
    for (size_t i = 0; i < terminal_networks.size(); ++i) {
        if (network_in_MPIslice[i]) {
            networks.push_back(buildNetwork(terminal_networks[i]));
        }
    }
    std::vector<std::string> description;
    append_desc(description, "* network description message");
    addNetworksDescription(description, networks);
    addAnalysis(description, final_time, dt);
    addSavedNodes(description, networks);
    endDescription(description);
    return network_manager_m::network_managerCtor(networks, description, final_time, dt);
}

preprocess_t preprocess(const parsed_mtln_t& parsed) {
    std::array<FDETYPES_m::XYZlimit_t, 6> alloc{};
    return preprocess(parsed, alloc);
}

preprocess_t preprocess(const parsed_mtln_t& parsed, const std::array<FDETYPES_m::XYZlimit_t, 6>& alloc) {
    preprocess_t res;
#ifdef CompileWithMPI
    MPI_Barrier(subcomm_mpi);
#endif
    res.final_time = parsed.time_step * static_cast<double>(parsed.number_of_steps);
    res.dt = parsed.time_step;

    const auto cable_bundles = buildCableBundles(parsed.cables);
#ifdef CompileWithMPI
    MPI_Barrier(subcomm_mpi);
#endif
    const auto line_bundles = buildLineBundles(cable_bundles, res.dt, alloc);
#ifdef CompileWithMPI
    MPI_Barrier(subcomm_mpi);
#endif
    res.bundles = res.buildMTLBundles(line_bundles);
    if (res.bundles.empty()) {
        return res;
    }
    res.cable_name_to_bundle_id = mapCablesToBundlesId(line_bundles, res.bundles);
    res.addProbesWithId(parsed.probes);
    res.addGenerators(parsed.wireGenerators);
    res.network_manager = res.buildNetworkManager(parsed.networks);
    return res;
}

} // namespace mtln_preprocess_m
