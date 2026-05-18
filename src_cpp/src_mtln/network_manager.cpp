#include <vector>
#include <string>
#include <memory>

// Forward declarations for types defined in other modules
// These would typically be in network_m.hpp, circuit_m.hpp, mtln_types_m.hpp, FDETYPES_m.hpp

namespace network_manager_m {

    // Assuming these types exist from included modules
    // struct network_t;
    // struct circuit_t;
    // struct node_source_t;
    // struct string_t;
    
    // Constants from FDETYPES_m
    using RKIND = double;
    using RKIND_TIEMPO = double;

    struct network_manager_t {
        std::vector<std::shared_ptr<struct network_t>> networks;
        struct circuit_t circuit;
        RKIND time;
        RKIND dt;

        // Methods
        void advanceVoltage();
        void updateCircuitCurrentsFromNetwork();
        void updateNetworkVoltages();
    };

    // Helper function declarations
    void appendToString_tArray(std::vector<std::shared_ptr<struct string_t>>& arr, const std::shared_ptr<struct string_t>& str);
    std::vector<std::shared_ptr<struct node_source_t>> copy_sources(const std::vector<std::shared_ptr<struct network_t>>& networks);
    std::vector<std::shared_ptr<struct string_t>> copy_node_names(const std::vector<std::shared_ptr<struct network_t>>& networks);
    network_manager_t network_managerCtor(const std::vector<std::shared_ptr<struct network_t>>& networks, 
                                          const std::vector<std::string>& description, 
                                          RKIND_TIEMPO final_time, 
                                          RKIND_TIEMPO dt);

    // Implementation of helper functions

    void appendToString_tArray(std::vector<std::shared_ptr<struct string_t>>& arr, const std::shared_ptr<struct string_t>& str) {
        // This has been implemented because there seems to be a bug in gfortran: 
        // https://fortran-lang.discourse.group/t/read-data-and-append-it-to-array-best-practice/1915
        // and arr = [ arr, str ] can't be used.
        auto old_arr = arr;
        arr.clear();
        arr.resize(old_arr.size() + 1);
        for (size_t i = 0; i < old_arr.size(); ++i) {
            arr[i] = old_arr[i];
        }
        arr[arr.size() - 1] = str;
    }

    std::vector<std::shared_ptr<struct node_source_t>> copy_sources(const std::vector<std::shared_ptr<struct network_t>>& networks) {
        std::vector<std::shared_ptr<struct node_source_t>> res;
        int n = 0;
        for (size_t i = 0; i < networks.size(); ++i) {
            for (size_t j = 0; j < networks[i]->nodes.size(); ++j) {
                n++;
            }
        }
        res.resize(n);
        n = 0;
        for (size_t i = 0; i < networks.size(); ++i) {
            for (size_t j = 0; j < networks[i]->nodes.size(); ++j) {
                res[n]->path_to_excitation = networks[i]->nodes[j]->source->path_to_excitation;
                res[n]->source_type = networks[i]->nodes[j]->source->source_type;
                res[n]->resistance = networks[i]->nodes[j]->source->resistance;
                n++;
            }
        }
        return res;
    }

    std::vector<std::shared_ptr<struct string_t>> copy_node_names(const std::vector<std::shared_ptr<struct network_t>>& networks) {
        std::vector<std::shared_ptr<struct string_t>> res;
        // Allocate with 0 size initially
        res.clear();
        for (size_t i = 0; i < networks.size(); ++i) {
            for (size_t j = 0; j < networks[i]->nodes.size(); ++j) {
                std::string trimmed_name = networks[i]->nodes[j]->name;
                // Assuming trim function exists or is handled by string_t constructor
                // In C++, std::string doesn't have leading/trailing whitespace trimming by default
                // We assume string_t handles this or the data is already clean
                auto temp = std::make_shared<struct string_t>(trimmed_name);
                appendToString_tArray(res, temp);
            }
        }
        auto time_str = std::make_shared<struct string_t>("time", 4);
        appendToString_tArray(res, time_str);
        return res;
    }

    network_manager_t network_managerCtor(const std::vector<std::shared_ptr<struct network_t>>& networks, 
                                          const std::vector<std::string>& description, 
                                          RKIND_TIEMPO final_time, 
                                          RKIND_TIEMPO dt) {
        network_manager_t res;
        bool printInput = true;
        res.dt = dt;
        res.time = 0.0;
        res.networks = networks;
        res.circuit.init(copy_node_names(networks), copy_sources(networks));
        res.circuit.dt = dt;
#ifdef CompileWithRelease
        printInput = false;
#endif        
        res.circuit.readInput(description, printInput);
        res.circuit.setModStopTimes(dt);
        return res;
    }

    void network_manager_t::updateNetworkVoltages() {
        for (size_t i = 0; i < networks.size(); ++i) {
            for (size_t j = 0; j < networks[i]->number_of_nodes; ++j) {
                networks[i]->nodes[j]->v = circuit.getNodeVoltage(networks[i]->nodes[j]->name);
            }
        }
    }

    void network_manager_t::updateCircuitCurrentsFromNetwork() {
        for (size_t i = 0; i < networks.size(); ++i) {
            for (size_t j = 0; j < networks[i]->number_of_nodes; ++j) {
                circuit.updateNodeCurrent(networks[i]->nodes[j]->name, networks[i]->nodes[j]->i);
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