#ifndef TEST_MTLN_TYPES_H
#define TEST_MTLN_TYPES_H

#include <gtest/gtest.h>

#include "mtln_types.h"

using namespace mtln_types_m;

TEST(mtln, mtln_types) {
    int err = 0;

    node_source_t node_source;
    node_source.path_to_excitation = "path";
    node_source.source_type = SOURCE_TYPE_VOLTAGE;

    terminal_node_t t;
    t.termination = termination_t{};
    t.termination.source = node_source;
    t.termination.termination_type = TERMINATION_SERIES;
    t.termination.resistance = 150;
    t.termination.inductance = 0.0;
    t.termination.capacitance = 1e22;

    if (t.termination.source.path_to_excitation != "path") {
        err += 1;
    }
    EXPECT_EQ(err, 0);
}

TEST(mtln, mtln_derived_types) {
    int err = 0;
    const std::string square_excitation = "coaxial_line_paul_8_6_0.25_square.exc";

    terminal_node_t node;
    node.conductor_in_cable = 1;
    node.side = TERMINAL_NODE_SIDE_INI;

    node_source_t node_source;
    node_source.path_to_excitation = square_excitation;
    node_source.source_type = SOURCE_TYPE_VOLTAGE;

    node.termination.termination_type = TERMINATION_SERIES;
    node.termination.resistance = 150;
    node.termination.inductance = 0.0;
    node.termination.capacitance = 1e22;
    node.termination.source = node_source;

    terminal_connection_t connection;
    connection.add_node(node);

    terminal_network_t network;
    network.add_connection(connection);

    if (network.connections[0].nodes[0].termination.source.path_to_excitation !=
        square_excitation) {
        err += 1;
    }
    EXPECT_EQ(err, 0);
}

#endif
