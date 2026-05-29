#ifndef MTLN_PREPROCESS_M_H
#define MTLN_PREPROCESS_M_H

#include <optional>
#include <string>
#include <vector>

#include "fhash_m.h"
#include "fdetypes_mtln.h"
#include "mtl_bundle_m.h"
#include "mtl_m.h"
#include "mtln_types.h"
#include "network_manager_m.h"
#include "network_m.h"

namespace mtln_preprocess_m {

using parsed_mtln_t = mtln_types_m::mtln_t;
using parsed_probe_t = mtln_types_m::probe_t;
using parsed_generator_t = mtln_types_m::parsed_generator_t;

struct cable_level_t {
    std::vector<mtln_types_m::cable_t*> cables;
};

struct cable_bundle_t {
    std::vector<cable_level_t> levels;
};

struct preprocess_t {
    std::vector<mtl_bundle_m::mtl_bundle_t> bundles;
    network_manager_m::network_manager_t network_manager;
    fhash_m::fhash_tbl_t conductors_before_cable;
    fhash_m::fhash_tbl_t cable_name_to_bundle_id;
    double final_time = 0.0;
    double dt = 0.0;

    std::vector<mtl_bundle_m::mtl_bundle_t> buildMTLBundles(
        const std::vector<mtl_m::transmission_line_bundle_t>& line_bundles);
    network_manager_m::network_manager_t buildNetworkManager(
        const std::vector<mtln_types_m::terminal_network_t>& terminal_networks);
    network_m::network_t buildNetwork(const mtln_types_m::terminal_network_t& terminal_network);
    void connectNodeToGround(const mtln_types_m::terminal_connection_t& terminal_connection,
                             std::vector<network_m::nw_node_t>& nodes,
                             std::vector<std::string>& description);
    void connectNodes(const mtln_types_m::terminal_connection_t& terminal_connection,
                      std::vector<network_m::nw_node_t>& nodes,
                      std::vector<std::string>& description);
    void connectNodesToNetworkCircuit(const mtln_types_m::terminal_connection_t& terminal_connection,
                                      std::vector<network_m::nw_node_t>& nodes,
                                      std::vector<std::string>& description);
    network_m::nw_node_t addNodeWithId(const mtln_types_m::terminal_node_t& node) const;
    void addProbesWithId(const std::vector<parsed_probe_t>& probes);
    void addGenerators(const std::vector<parsed_generator_t>& generators);

private:
    std::string nodeSideToString(int side) const;
};

preprocess_t preprocess(const parsed_mtln_t& parsed);
preprocess_t preprocess(const parsed_mtln_t& parsed, const std::array<FDETYPES_m::XYZlimit_t, 6>& alloc);

std::vector<int> conductorsInLevel(const mtl_m::transmission_line_bundle_t& line);
int findConductorsBeforeCable(const std::string& name, const mtl_m::transmission_line_level_t& level);

} // namespace mtln_preprocess_m

#endif
