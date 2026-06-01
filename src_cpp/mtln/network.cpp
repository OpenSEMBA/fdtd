#include "network_m.h"

namespace network_m {

int countNodes(const std::vector<terminal_connection_t>& connections) {
    int count = 0;
    for (const auto& connection : connections) {
        count += static_cast<int>(connection.nodes.size());
    }
    return count;
}

network_t networkCtor(const std::vector<nw_node_t>& nodes, const std::vector<std::string>& description) {
    network_t res;
    res.nodes = nodes;
    res.description = description;
    res.number_of_nodes = static_cast<int>(nodes.size());
    return res;
}

} // namespace network_m
