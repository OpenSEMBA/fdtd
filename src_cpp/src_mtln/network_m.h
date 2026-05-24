#ifndef NETWORK_M_H
#define NETWORK_M_H

#include <string>
#include <vector>

#include "mtln_types.h"

namespace network_m {

using mtln_types_m::node_source_t;
using mtln_types_m::terminal_connection_t;

struct nw_node_t {
    std::string name;
    node_source_t source;
    int source_type = 0;
    double line_c_per_meter = 0.0;
    double line_g_per_meter = 0.0;
    double step = 0.0;
    double v = 0.0;
    double i = 0.0;
    int bundle_number = 0;
    int conductor_number = 0;
    int v_index = 0;
    int i_index = 0;
    int side = 0;
    bool open = false;
};

struct network_t {
    int number_of_nodes = 0;
    std::vector<nw_node_t> nodes;
    std::vector<std::string> description;
};

int countNodes(const std::vector<terminal_connection_t>& connections);
network_t networkCtor(const std::vector<nw_node_t>& nodes, const std::vector<std::string>& description);

} // namespace network_m

#endif
