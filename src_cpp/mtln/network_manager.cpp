#include "network_manager_m.h"

#include <cstdlib>

namespace network_manager_m {

namespace {

void appendToString_tArray(std::vector<circuit_m::string_t>& arr, const circuit_m::string_t& str) {
    const auto old_arr = arr;
    arr.resize(old_arr.size() + 1);
    for (size_t i = 0; i < old_arr.size(); ++i) {
        arr[i] = old_arr[i];
    }
    arr.back() = str;
}

std::vector<mtln_types_m::node_source_t> copy_sources(const std::vector<network_t>& networks) {
    std::vector<mtln_types_m::node_source_t> res;
    for (const auto& network : networks) {
        for (const auto& node : network.nodes) {
            res.push_back(node.source);
        }
    }
    return res;
}

std::vector<circuit_m::string_t> copy_node_names(const std::vector<network_t>& networks) {
    std::vector<circuit_m::string_t> res;
    for (const auto& network : networks) {
        for (const auto& node : network.nodes) {
            res.emplace_back(node.name, static_cast<int>(node.name.size()));
        }
    }
    res.emplace_back("time", 4);
    return res;
}

} // namespace

network_manager_t network_managerCtor(const std::vector<network_t>& networks,
                                      const std::vector<std::string>& description,
                                      RKIND_TIEMPO final_time,
                                      RKIND_TIEMPO dt) {
    (void)final_time;
    network_manager_t res;
    const bool printInput = (std::getenv("MTLN_PRINT_NETLIST") != nullptr);
    res.dt = dt;
    res.time = 0.0;
    res.networks = networks;
    res.circuit.init(copy_node_names(networks), copy_sources(networks));
    res.circuit.dt = dt;
    res.circuit.readInput(description, printInput);
    res.circuit.setModStopTimes(dt);
    return res;
}

void network_manager_t::updateNetworkVoltages() {
    for (auto& network : networks) {
        for (auto& node : network.nodes) {
            node.v = circuit.getNodeVoltage(node.name);
        }
    }
}

void network_manager_t::updateCircuitCurrentsFromNetwork() {
    for (const auto& network : networks) {
        for (const auto& node : network.nodes) {
            circuit.updateNodeCurrent(node.name, node.i);
        }
    }
}

void network_manager_t::advanceVoltage() {
    updateCircuitCurrentsFromNetwork();
    circuit.step();
    circuit.time = circuit.time + circuit.dt;
    updateNetworkVoltages();
}

} // namespace network_manager_m
