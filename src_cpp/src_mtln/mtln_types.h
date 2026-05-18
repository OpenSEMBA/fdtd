#include <vector>
#include <string>
#include <complex>
#include <memory>
#include <algorithm>
#include <iostream>

// Forward declarations for types defined in other modules
// Assuming FDETYPES_m provides direction_t, BUFSIZE, RKIND, RKIND_TIEMPO
// Assuming cable_t is defined elsewhere or in a separate header
// For this translation, we assume cable_t is a forward declared class or defined in another header.
// If cable_t is not defined here, it must be included.
// Since the prompt implies this is a single file translation request but uses external types,
// we will assume necessary headers are included or types are forward declared.

// Placeholder for external types if not provided in the snippet
// In a real scenario, these would be in "FDETYPES_m.h" or similar.
#ifndef FDETYPES_M_H
#define FDETYPES_M_H
#include <cstdint>

// Assuming RKIND is double and RKIND_TIEMPO is double based on typical usage
// BUFSIZE is likely an integer constant
constexpr int BUFSIZE = 256; 
using RKIND = double;
using RKIND_TIEMPO = double;

// direction_t is likely an enum or base class
enum class direction_t {
    // Values depend on FDETYPES_m, assuming standard directions or empty base
    // Since segment_t extends direction_t, direction_t is likely a base class or enum
    // If it's an enum, inheritance in C++ structs isn't direct. 
    // Given segment_t has extra fields, direction_t is likely a base struct.
    // We will define it as an empty base struct for inheritance purposes if needed,
    // or assume it contains common direction fields. 
    // Without FDETYPES_m content, we assume it's a base struct.
};

#endif

// Forward declaration of cable_t as it is used as a pointer base class
class cable_t;

namespace mtln_types_m {

    // Constants
    constexpr int TERMINATION_UNDEFINED = -1;
    constexpr int TERMINATION_SHORT = 1;
    constexpr int TERMINATION_OPEN = 2;
    constexpr int TERMINATION_SERIES = 3;
    constexpr int TERMINATION_PARALLEL = 4;
    constexpr int TERMINATION_RsLCp = 5;
    constexpr int TERMINATION_RLsCp = 6;
    constexpr int TERMINATION_LsRCp = 7;
    constexpr int TERMINATION_CsLRp = 8;
    constexpr int TERMINATION_RCsLp = 9;
    constexpr int TERMINATION_LCsRp = 10;
    constexpr int TERMINATION_CIRCUIT = 11;
    constexpr int TERMINATION_NETWORK = 12

    constexpr int TERMINAL_NODE_SIDE_UNDEFINED = -1;
    constexpr int TERMINAL_NODE_SIDE_INI = 1;
    constexpr int TERMINAL_NODE_SIDE_END = 2;

    constexpr int TRANSFER_IMPEDANCE_DIRECTION_UNDEFINED = -1;
    constexpr int TRANSFER_IMPEDANCE_DIRECTION_INWARDS = 1;
    constexpr int TRANSFER_IMPEDANCE_DIRECTION_OUTWARDS = 2;
    constexpr int TRANSFER_IMPEDANCE_DIRECTION_BOTH = 3;

    constexpr int PROBE_TYPE_UNDEFINED = -1;
    constexpr int PROBE_TYPE_VOLTAGE = 1;
    constexpr int PROBE_TYPE_CURRENT = 2;

    constexpr int SOURCE_TYPE_UNDEFINED = -1;
    constexpr int SOURCE_TYPE_VOLTAGE = 1;
    constexpr int SOURCE_TYPE_CURRENT = 2;

    constexpr int DIRECTION_X_POS = 1;
    constexpr int DIRECTION_X_NEG = -1;
    constexpr int DIRECTION_Y_POS = 2;
    constexpr int DIRECTION_Y_NEG = -2;
    constexpr int DIRECTION_Z_POS = 3;
    constexpr int DIRECTION_Z_NEG = -3;

    // Derived Types

    struct parsed_generator_t {
        std::string path_to_excitation = "";
        int generator_type = SOURCE_TYPE_UNDEFINED;
        cable_t* attached_to_cable = nullptr;
        RKIND resistance = 0.0;
        int index = -1;
        int conductor = -1;

        bool operator==(const parsed_generator_t& other) const {
            return wire_source_eq(*this, other);
        }
    };

    struct node_source_t {
        std::string path_to_excitation = "";
        int source_type = SOURCE_TYPE_UNDEFINED;
        RKIND resistance = 0.0;
    };

    struct terminal_circuit_t {
        std::string file = "";
        std::string name = "";
    };

    struct termination_t {
        int termination_type = TERMINATION_UNDEFINED;
        RKIND resistance = 0.0;
        RKIND inductance = 0.0;
        RKIND capacitance = 1e22;
        node_source_t source;
        terminal_circuit_t model;
        int networkCircuitNode = -1;

        bool operator==(const termination_t& other) const {
            return termination_eq(*this, other);
        }
    };

    struct terminal_node_t {
        cable_t* belongs_to_cable = nullptr;
        int conductor_in_cable = 0;
        int side = TERMINAL_NODE_SIDE_UNDEFINED;
        termination_t termination;

        bool operator==(const terminal_node_t& other) const {
            return terminal_node_eq(*this, other);
        }
    };

    struct network_circuit_t {
        std::string model_file = "";
        std::string model_name = "";
        std::string circuit_name = "";
        int number_of_nodes = -1;
        int nodeId = -1;
    };

    struct terminal_connection_t {
        std::vector<terminal_node_t> nodes;
        network_circuit_t network_circuit;

        bool operator==(const terminal_connection_t& other) const {
            return terminal_connection_eq(*this, other);
        }

        void add_node(const terminal_node_t& node) {
            terminal_connection_add_node(*this, node);
        }
    };

    struct terminal_network_t {
        std::vector<terminal_connection_t> connections;

        bool operator==(const terminal_network_t& other) const {
            return terminal_network_eq(*this, other);
        }

        void add_connection(const terminal_connection_t& connection) {
            terminal_network_add_connection(*this, connection);
        }
    };

    struct transfer_impedance_per_meter_t {
        RKIND inductive_term = 0.0;
        RKIND resistive_term = 0.0;
        std::vector<std::complex<RKIND>> poles;
        std::vector<std::complex<RKIND>> residues;
        int direction = TRANSFER_IMPEDANCE_DIRECTION_UNDEFINED;

        bool operator==(const transfer_impedance_per_meter_t& other) const {
            return transfer_impedance_per_meter_eq(*this, other);
        }

        bool has_transfer_impedance() const {
            return ::mtln_types_m::has_transfer_impedance(*this);
        }
    };

    struct connector_t {
        int id = 0;
        std::vector<RKIND> resistances;
        std::vector<transfer_impedance_per_meter_t> transfer_impedances_per_meter;

        bool operator==(const connector_t& other) const {
            return connector_eq(*this, other);
        }
    };

    struct multipolar_coefficient_t {
        RKIND a = 0.0;
        RKIND b = 0.0;

        bool operator==(const multipolar_coefficient_t& other) const {
            return multipolar_coefficient_eq(*this, other);
        }
    };

    struct field_reconstruction_t {
        RKIND inner_region_average_potential = 0.0;
        std::vector<RKIND> expansion_center = {0.0, 0.0};
        std::vector<multipolar_coefficient_t> ab;
        std::vector<RKIND> conductor_potentials;

        bool operator==(const field_reconstruction_t& other) const {
            return field_reconstruction_eq(*this, other);
        }
    };

    struct box_2d_t {
        std::vector<RKIND> min = {0.0, 0.0};
        std::vector<RKIND> max = {0.0, 0.0};

        bool operator==(const box_2d_t& other) const {
            return box_2d_eq(*this, other);
        }
    };

    struct multipolar_expansion_t {
        box_2d_t inner_region;
        std::vector<field_reconstruction_t> electric;
        std::vector<field_reconstruction_t> magnetic;

        bool operator==(const multipolar_expansion_t& other) const {
            return multipolar_expansion_eq(*this, other);
        }
    };

    // segment_t extends direction_t. Since direction_t is likely a base struct,
    // we inherit from it. If direction_t is just an enum, this mapping is tricky.
    // Assuming direction_t is a struct with common direction fields.
    struct segment_t : public direction_t {
        box_2d_t dualBox;
        RKIND d1 = 0.0;
        RKIND d2 = 0.0;
    };

    // Forward declaration needed for cable_t methods
    class cable_t;

    struct cable_t {
        std::string name;
        std::vector<RKIND> step_size;
        std::vector<segment_t> segments;
        connector_t* initial_connector = nullptr;
        connector_t* end_connector = nullptr;
        std::string tag; // Fixed size buffer handled by std::string or char array
        int n_segments = 0;

        bool operator==(const cable_t& other) const {
            return cable_eq(*this, other);
        }
    };

    struct unshielded_multiwire_t : public cable_t {
        std::vector<std::vector<RKIND>> cell_inductance_per_meter;
        std::vector<std::vector<RKIND>> cell_capacitance_per_meter;
        std::vector<std::vector<RKIND>> resistance_per_meter;
        std::vector<std::vector<RKIND>> conductance_per_meter;
        RKIND radius = 0.0;
        std::vector<multipolar_expansion_t> multipolar_expansion;
    };

    struct shielded_multiwire_t : public cable_t {
        std::vector<std::vector<RKIND>> resistance_per_meter;
        std::vector<std::vector<RKIND>> conductance_per_meter;
        std::vector<std::vector<RKIND>> inductance_per_meter;
        std::vector<std::vector<RKIND>> capacitance_per_meter;
        transfer_impedance_per_meter_t transfer_impedance;
        cable_t* parent_cable = nullptr;
        int conductor_in_parent = -1;
    };

    struct probe_t {
        cable_t* attached_to_cable = nullptr;
        int index = 0;
        int probe_type = PROBE_TYPE_UNDEFINED;
        std::string probe_name;
        std::vector<RKIND> probe_position = {0.0, 0.0, 0.0};

        bool operator==(const probe_t& other) const {
            return probe_eq(*this, other);
        }
    };

    struct cable_abstract_t {
        cable_t* ptr = nullptr;
    };

    struct mtln_t {
        std::vector<cable_abstract_t> cables;
        std::vector<terminal_network_t> networks;
        std::vector<probe_t> probes;
        std::vector<parsed_generator_t> wireGenerators;
        std::vector<connector_t>* connectors = nullptr; // Pointer to vector as in Fortran
        RKIND_TIEMPO time_step = 0.0;
        int number_of_steps = 0;
        int n_sh = 0;
        int n_unsh = 0;

        bool operator==(const mtln_t& other) const {
            return mtln_eq(*this, other);
        }
    };

    // Function Implementations

    inline bool mtln_eq(const mtln_t& a, const mtln_t& b) {
        if (a.time_step != b.time_step) return false;
        if (a.number_of_steps != b.number_of_steps) return false;

        if (a.cables.size() != b.cables.size()) return false;
        for (size_t i = 0; i < a.cables.size(); ++i) {
            // Check if pointers are equal. Note: Fortran == on pointers checks association.
            // In C++, we check if both are null or both point to the same object.
            if (a.cables[i].ptr != b.cables[i].ptr) return false;
        }

        if (a.probes.size() != b.probes.size()) return false;
        for (size_t i = 0; i < a.probes.size(); ++i) {
            if (!(a.probes[i] == b.probes[i])) return false;
        }

        if (a.networks.size() != b.networks.size()) return false;
        for (size_t i = 0; i < a.networks.size(); ++i) {
            if (!(a.networks[i] == b.networks[i])) return false;
        }

        return true;
    }

    inline bool transfer_impedance_per_meter_eq(const transfer_impedance_per_meter_t& a, const transfer_impedance_per_meter_t& b) {
        return (a.inductive_term == b.inductive_term) &&
               (a.resistive_term == b.resistive_term) &&
               (a.poles == b.poles) &&
               (a.residues == b.residues) &&
               (a.direction == b.direction);
    }

    inline bool multipolar_coefficient_eq(const multipolar_coefficient_t& a, const multipolar_coefficient_t& b) {
        return (a.a == b.a) && (a.b == b.b);
    }

    inline bool field_reconstruction_eq(const field_reconstruction_t& lhs, const field_reconstruction_t& rhs) {
        bool res = true;
        res = res && (lhs.inner_region_average_potential == rhs.inner_region_average_potential);
        res = res && (lhs.expansion_center == rhs.expansion_center);
        res = res && (lhs.ab.size() > 0 && rhs.ab.size() > 0); // allocated check
        res = res && (lhs.ab == rhs.ab);
        res = res && (lhs.conductor_potentials.size() > 0 && rhs.conductor_potentials.size() > 0); // allocated check
        res = res && (lhs.conductor_potentials == rhs.conductor_potentials);
        return res;
    }

    inline bool box_2d_eq(const box_2d_t& a, const box_2d_t& b) {
        return (a.min == b.min) && (a.max == b.max);
    }

    inline bool multipolar_expansion_eq(const multipolar_expansion_t& a, const multipolar_expansion_t& b) {
        bool res = true;
        res = res && (a.inner_region == b.inner_region);
        res = res && (a.electric.size() > 0 && b.electric.size() > 0);
        res = res && (a.electric == b.electric);
        res = res && (a.magnetic.size() > 0 && b.magnetic.size() > 0);
        res = res && (a.magnetic == b.magnetic);
        return res;
    }

    inline bool cable_eq(const cable_t& a, const cable_t& b) {
        bool res = true;
        res = res && (a.name == b.name);
        res = res && (a.step_size == b.step_size);
        res = res && (a.segments.size() == b.segments.size());
        for (size_t i = 0; i < a.segments.size(); ++i) {
            res = res && (a.segments[i] == b.segments[i]);
        }

        // Check initial_connector
        if (a.initial_connector && b.initial_connector) {
            res = res && (a.initial_connector == b.initial_connector);
        } else if (!a.initial_connector && !b.initial_connector) {
            res = res && true;
        } else {
            res = res && false;
        }

        // Check end_connector
        if (a.end_connector && b.end_connector) {
            res = res && (a.end_connector == b.end_connector);
        } else if (!a.end_connector && !b.end_connector) {
            res = res && true;
        } else {
            res = res && false;
        }

        // Type-specific checks using dynamic_cast
        const shielded_multiwire_t* shielded_a = dynamic_cast<const shielded_multiwire_t*>(&a);
        const shielded_multiwire_t* shielded_b = dynamic_cast<const shielded_multiwire_t*>(&b);
        const unshielded_multiwire_t* unshielded_a = dynamic_cast<const unshielded_multiwire_t*>(&a);
        const unshielded_multiwire_t* unshielded_b = dynamic_cast<const unshielded_multiwire_t*>(&b);

        if (shielded_a) {
            if (shielded_b) {
                res = res && (shielded_a->inductance_per_meter == shielded_b->inductance_per_meter);
                res = res && (shielded_a->capacitance_per_meter == shielded_b->capacitance_per_meter);
                res = res && (shielded_a->resistance_per_meter == shielded_b->resistance_per_meter);
                res = res && (shielded_a->conductance_per_meter == shielded_b->conductance_per_meter);
                res = res && (shielded_a->transfer_impedance == shielded_b->transfer_impedance);
                res = res && (shielded_a->conductor_in_parent == shielded_b->conductor_in_parent);
                
                if (shielded_a->parent_cable && shielded_b->parent_cable) {
                    res = res && (shielded_a->parent_cable == shielded_b->parent_cable);
                } else if (!shielded_a->parent_cable && !shielded_b->parent_cable) {
                    res = res && true;
                } else {
                    res = res && false;
                }
            } else {
                res = false;
            }
        } else if (unshielded_a) {
            if (unshielded_b) {
                res = res && (unshielded_a->multipolar_expansion == unshielded_b->multipolar_expansion);
                res = res && (unshielded_a->cell_inductance_per_meter == unshielded_b->cell_inductance_per_meter);
                res = res && (unshielded_a->cell_capacitance_per_meter == unshielded_b->cell_capacitance_per_meter);
            } else {
                res = false;
            }
        }

        return res;
    }

    inline bool connector_eq(const connector_t& a, const connector_t& b) {
        return (a.id == b.id) &&
               (a.resistances == b.resistances) &&
               (a.transfer_impedances_per_meter == b.transfer_impedances_per_meter);
    }

    inline bool wire_source_eq(const parsed_generator_t& a, const parsed_generator_t& b) {
        bool res = (a.path_to_excitation == b.path_to_excitation) &&
                   (a.generator_type == b.generator_type) &&
                   (a.resistance == b.resistance) &&
                   (a.index == b.index);
        
        if (!a.attached_to_cable || !b.attached_to_cable) {
            res = res && false;
        } else {
            res = res && (a.attached_to_cable == b.attached_to_cable);
        }
        return res;
    }

    inline bool termination_eq(const termination_t& a, const termination_t& b) {
        return (a.termination_type == b.termination_type) &&
               (a.resistance == b.resistance) &&
               (a.inductance == b.inductance) &&
               (a.capacitance == b.capacitance) &&
               (a.source.path_to_excitation == b.source.path_to_excitation) &&
               (a.source.source_type == b.source.source_type);
    }

    inline bool probe_eq(const probe_t& a, const probe_t& b) {
        bool res = (a.index == b.index) &&
                   (a.probe_type == b.probe_type) &&
                   (a.probe_name == b.probe_name) &&
                   (a.probe_position == b.probe_position);

        if (!a.attached_to_cable && !b.attached_to_cable) {
            res = res && true;
        } else if ((a.attached_to_cable && !b.attached_to_cable) || (!a.attached_to_cable && b.attached_to_cable)) {
            res = res && false;
        } else {
            res = res && (a.attached_to_cable == b.attached_to_cable);
        }
        return res;
    }

    inline bool terminal_node_eq(const terminal_node_t& a, const terminal_node_t& b) {
        bool res = (a.conductor_in_cable == b.conductor_in_cable) &&
                   (a.side == b.side) &&
                   (a.termination == b.termination);

        if (!a.belongs_to_cable && !b.belongs_to_cable) {
            res = res && true;
        } else if ((a.belongs_to_cable && !b.belongs_to_cable) || (!a.belongs_to_cable && b.belongs_to_cable)) {
            res = res && false;
        } else {
            res = res && (a.belongs_to_cable == b.belongs_to_cable);
        }
        return res;
    }

    inline bool terminal_connection_eq(const terminal_connection_t& a, const terminal_connection_t& b) {
        if (a.nodes.size() != b.nodes.size()) return false;
        for (size_t i = 0; i < a.nodes.size(); ++i) {
            if (!(a.nodes[i] == b.nodes[i])) return false;
        }
        return true;
    }

    inline bool terminal_network_eq(const terminal_network_t& a, const terminal_network_t& b) {
        if (a.connections.size() != b.connections.size()) return false;
        for (size_t i = 0; i < a.connections.size(); ++i) {
            if (!(a.connections[i] == b.connections[i])) return false;
        }
        return true;
    }

    inline void terminal_connection_add_node(terminal_connection_t& this_obj, const terminal_node_t& node) {
        this_obj.nodes.push_back(node);
    }

    inline bool has_transfer_impedance(const transfer_impedance_per_meter_t& this_obj) {
        return (this_obj.resistive_term != 0) || (this_obj.inductive_term != 0) ||
               (this_obj.poles.size() != 0) || (this_obj.residues.size() != 0);
    }

    inline void terminal_network_add_connection(terminal_network_t& this_obj, const terminal_connection_t& connection) {
        this_obj.connections.push_back(connection);
    }

} // namespace mtln_types_m