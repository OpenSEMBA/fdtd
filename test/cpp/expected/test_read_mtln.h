#ifndef TEST_READ_MTLN_H
#define TEST_READ_MTLN_H

#ifdef CompileWithMTLN

#include "test_smbjson_helpers.h"

using namespace mtln_types_m;

namespace {

inline std::vector<std::vector<double>> matrix2x2(double a, double b, double c, double d) {
    return {{a, b}, {c, d}};
}

inline std::vector<std::vector<double>> matrix8x8Blocks() {
    auto block = matrix2x2(2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07);
    std::vector<std::vector<double>> m(8, std::vector<double>(8, 0.0));
    for (int b = 0; b < 4; ++b) {
        int o = b * 2;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j) {
                m[o + i][o + j] = block[i][j];
            }
    }
    return m;
}

inline std::vector<std::vector<double>> capacitance8x8Blocks() {
    auto cap = matrix2x2(105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12);
    std::vector<std::vector<double>> m(8, std::vector<double>(8, 0.0));
    for (int b = 0; b < 4; ++b) {
        int o = b * 2;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j) {
                m[o + i][o + j] = cap[i][j];
            }
    }
    return m;
}

inline void setXSegments(std::vector<segment_t>& segments, int y, int z, int count, int xStart) {
    segments.resize(count);
    for (int i = 0; i < count; ++i) {
        segments[i].x = xStart + i + 1;
        segments[i].y = y;
        segments[i].z = z;
        segments[i].orientation = DIRECTION_X_POS;
    }
}

inline void setYNegSegments(std::vector<segment_t>& segments, int x, int z, int count) {
    segments.resize(count);
    for (int i = 0; i < count; ++i) {
        segments[i].x = x;
        segments[i].y = 7 - (i + 1);
        segments[i].z = z;
        segments[i].orientation = DIRECTION_Y_NEG;
    }
}

} // namespace

inline Parseador_t expectedReadMtln() {
    Parseador_t expected = buildExpected();

    expected.general->dt = 1e-12;
    expected.general->nmax = 1000;

    expected.matriz->totalX = 101;
    expected.matriz->totalY = 8;
    expected.matriz->totalZ = 3;

    expected.despl->nX = 1;
    expected.despl->nY = 1;
    expected.despl->nZ = 1;
    expected.despl->desX.resize(1);
    expected.despl->desY.resize(1);
    expected.despl->desZ.resize(1);
    expected.despl->desX[0] = 0.1;
    expected.despl->desY[0] = 0.1;
    expected.despl->desZ[0] = 0.1;
    expected.despl->mx1 = 0;
    expected.despl->mx2 = 100;
    expected.despl->my1 = 0;
    expected.despl->my2 = 7;
    expected.despl->mz1 = 0;
    expected.despl->mz2 = 2;

    for (int i = 0; i < 6; ++i)
        expected.front->tipoFrontera[i] = F_MUR;

    expected.nodSrc->NodalSource.resize(1);
    expected.nodSrc->n_C1P_max = 1;
    expected.nodSrc->n_C2P_max = 1;
    expected.nodSrc->n_nodSrc = 1;
    expected.nodSrc->n_nodSrc_max = 1;
    expected.nodSrc->NodalSource[0].nombre = "gauss.exc";
    expected.nodSrc->NodalSource[0].isElec = true;
    expected.nodSrc->NodalSource[0].isHard = false;
    expected.nodSrc->NodalSource[0].isInitialValue = false;
    expected.nodSrc->NodalSource[0].n_C2P = 1;
    expected.nodSrc->NodalSource[0].c2P.resize(1);
    auto& nodC2P = expected.nodSrc->NodalSource[0].c2P[0];
    nodC2P.Xi = 2;
    nodC2P.Xe = 3;
    nodC2P.Yi = 7;
    nodC2P.Ye = 7;
    nodC2P.Zi = 1;
    nodC2P.Ze = 1;
    nodC2P.xc = 1;
    nodC2P.yc = 0;
    nodC2P.zc = 0;
    nodC2P.*(&NFDETypes_m::coords_scaled_t::Or) = 1;
    nodC2P.tag = "";

    expected.mtln->time_step = 1e-12;
    expected.mtln->number_of_steps = 1000;

    expected.mtln->connectors.resize(4);
    expected.mtln->connectors[0].id = 24;
    expected.mtln->connectors[0].resistances = {100.0e-3};
    expected.mtln->connectors[0].transfer_impedances_per_meter.clear();
    expected.mtln->connectors[1].id = 25;
    expected.mtln->connectors[1].resistances = {19.0};
    expected.mtln->connectors[1].transfer_impedances_per_meter.clear();
    expected.mtln->connectors[2].id = 204;
    expected.mtln->connectors[2].resistances = {100.0e-3};
    expected.mtln->connectors[2].transfer_impedances_per_meter.resize(1);
    expected.mtln->connectors[2].transfer_impedances_per_meter[0].direction =
        TRANSFER_IMPEDANCE_DIRECTION_INWARDS;
    expected.mtln->connectors[2].transfer_impedances_per_meter[0].resistive_term = 3.33;
    expected.mtln->connectors[2].transfer_impedances_per_meter[0].inductive_term = 2.6e-9;
    expected.mtln->connectors[3].id = 205;
    expected.mtln->connectors[3].resistances = {19.0};
    expected.mtln->connectors[3].transfer_impedances_per_meter.resize(1);
    expected.mtln->connectors[3].transfer_impedances_per_meter[0].direction =
        TRANSFER_IMPEDANCE_DIRECTION_INWARDS;
    expected.mtln->connectors[3].transfer_impedances_per_meter[0].resistive_term = 609.3;
    expected.mtln->connectors[3].transfer_impedances_per_meter[0].inductive_term = 2.6e-9;

    auto& conn = expected.mtln->connectors;

    expected.mtln->cables.resize(9);

    expected.mtln->cables[0].ptr = std::make_unique<unshielded_multiwire_t>();
    auto* c1 = static_cast<unshielded_multiwire_t*>(expected.mtln->cables[0].ptr.get());
    initializeCablePULParameters(c1, 1);
    c1->name = "line_0_0";
    c1->cell_inductance_per_meter = {{5.481553487168089e-07}};
    c1->cell_capacitance_per_meter = {{2.0270004E-11}};
    c1->resistance_per_meter = {{22.9e-3}};
    c1->multipolar_expansion.clear();
    c1->step_size.resize(9, 0.1);
    setXSegments(c1->segments, 7, 1, 9, 0);
    c1->n_segments = 9;
    c1->initial_connector = &conn[0];
    c1->end_connector = nullptr;

    expected.mtln->cables[1].ptr = std::make_unique<shielded_multiwire_t>();
    auto* c2 = static_cast<shielded_multiwire_t*>(expected.mtln->cables[1].ptr.get());
    initializeCablePULParameters(c2, 1);
    c2->name = "line_1_0";
    c2->inductance_per_meter = {{8.802075200000001e-08}};
    c2->capacitance_per_meter = {{5.5840010E-10}};
    c2->resistance_per_meter = {{3.9e-3}};
    c2->step_size.resize(9, 0.1);
    setXSegments(c2->segments, 7, 1, 9, 0);
    c2->transfer_impedance.direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS;
    c2->transfer_impedance.resistive_term = 0.0;
    c2->transfer_impedance.inductive_term = 8.9e-9;
    c2->parent_cable = c1;
    c2->conductor_in_parent = 1;
    c2->n_segments = 9;
    c2->initial_connector = &conn[2];
    c2->end_connector = nullptr;

    expected.mtln->cables[2].ptr = std::make_unique<shielded_multiwire_t>();
    auto* c3 = static_cast<shielded_multiwire_t*>(expected.mtln->cables[2].ptr.get());
    initializeCablePULParameters(c3, 8);
    c3->name = "line_2_0";
    c3->inductance_per_meter = matrix8x8Blocks();
    c3->capacitance_per_meter = capacitance8x8Blocks();
    for (int i = 0; i < 8; ++i)
        c3->resistance_per_meter[i][i] = 62.0e-3;
    c3->step_size.resize(9, 0.1);
    setXSegments(c3->segments, 7, 1, 9, 0);
    c3->transfer_impedance.direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS;
    c3->transfer_impedance.resistive_term = 0.0;
    c3->transfer_impedance.inductive_term = 4.2e-9;
    c3->parent_cable = c2;
    c3->conductor_in_parent = 1;
    c3->n_segments = 9;

    expected.mtln->cables[3].ptr = std::make_unique<unshielded_multiwire_t>();
    auto* c4 = static_cast<unshielded_multiwire_t*>(expected.mtln->cables[3].ptr.get());
    initializeCablePULParameters(c4, 1);
    c4->name = "line_0_1";
    c4->cell_inductance_per_meter = {{6.482560773828984e-07}};
    c4->cell_capacitance_per_meter = {{1.7140003E-11}};
    c4->resistance_per_meter = {{11.8e-3}};
    c4->multipolar_expansion.clear();
    c4->step_size.resize(8, 0.1);
    setXSegments(c4->segments, 7, 1, 8, 9);
    c4->n_segments = 8;
    c4->initial_connector = &conn[1];
    c4->end_connector = nullptr;

    expected.mtln->cables[4].ptr = std::make_unique<shielded_multiwire_t>();
    auto* c5 = static_cast<shielded_multiwire_t*>(expected.mtln->cables[4].ptr.get());
    initializeCablePULParameters(c5, 1);
    c5->name = "line_1_1";
    c5->inductance_per_meter = {{1.37228e-07}};
    c5->capacitance_per_meter = {{3.2310005E-10}};
    c5->resistance_per_meter = {{12.2e-3}};
    c5->conductance_per_meter = {{0.0}};
    c5->step_size.resize(8, 0.1);
    setXSegments(c5->segments, 7, 1, 8, 9);
    c5->transfer_impedance.direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS;
    c5->transfer_impedance.resistive_term = 0.0;
    c5->transfer_impedance.inductive_term = 7.4e-9;
    c5->parent_cable = c4;
    c5->conductor_in_parent = 1;
    c5->n_segments = 8;
    c5->initial_connector = &conn[3];

    expected.mtln->cables[5].ptr = std::make_unique<shielded_multiwire_t>();
    auto* c6 = static_cast<shielded_multiwire_t*>(expected.mtln->cables[5].ptr.get());
    initializeCablePULParameters(c6, 2);
    c6->name = "line_2_4";
    c6->inductance_per_meter = matrix2x2(2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07);
    c6->capacitance_per_meter = matrix2x2(105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12);
    for (int i = 0; i < 2; ++i)
        c6->resistance_per_meter[i][i] = 62.0e-3;
    c6->step_size.resize(8, 0.1);
    setXSegments(c6->segments, 7, 1, 8, 9);
    c6->transfer_impedance.direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS;
    c6->transfer_impedance.resistive_term = 0.0;
    c6->transfer_impedance.inductive_term = 4.2e-9;
    c6->parent_cable = c5;
    c6->conductor_in_parent = 1;
    c6->n_segments = 8;

    expected.mtln->cables[6].ptr = std::make_unique<unshielded_multiwire_t>();
    auto* c7 = static_cast<unshielded_multiwire_t*>(expected.mtln->cables[6].ptr.get());
    initializeCablePULParameters(c7, 1);
    c7->name = "line_0_2";
    c7->cell_inductance_per_meter = {{5.802145885361537e-07}};
    c7->cell_capacitance_per_meter = {{1.9150003E-11}};
    c7->resistance_per_meter = {{17.3e-3}};
    c7->multipolar_expansion.clear();
    c7->step_size.resize(7, 0.1);
    setYNegSegments(c7->segments, 10, 1, 7);
    c7->n_segments = 7;

    expected.mtln->cables[7].ptr = std::make_unique<shielded_multiwire_t>();
    auto* c8 = static_cast<shielded_multiwire_t*>(expected.mtln->cables[7].ptr.get());
    initializeCablePULParameters(c8, 1);
    c8->name = "line_1_2";
    c8->inductance_per_meter = {{9.1890502e-08}};
    c8->capacitance_per_meter = {{4.7190007E-10}};
    c8->resistance_per_meter = {{6.5e-3}};
    c8->conductance_per_meter = {{0.0}};
    c8->step_size.resize(7, 0.1);
    setYNegSegments(c8->segments, 10, 1, 7);
    c8->transfer_impedance.direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS;
    c8->transfer_impedance.resistive_term = 0.0;
    c8->transfer_impedance.inductive_term = 3.0e-9;
    c8->parent_cable = c7;
    c8->conductor_in_parent = 1;
    c8->n_segments = 7;

    expected.mtln->cables[8].ptr = std::make_unique<shielded_multiwire_t>();
    auto* c9 = static_cast<shielded_multiwire_t*>(expected.mtln->cables[8].ptr.get());
    initializeCablePULParameters(c9, 6);
    c9->name = "line_2_5";
    auto blockL = matrix2x2(2.4382084E-07, 4.7377505E-08, 4.7377508E-08, 2.4382081E-07);
    auto blockC = matrix2x2(105.5e-12, -20.5e-12, -20.5e-12, 105.5e-12);
    c9->inductance_per_meter.assign(6, std::vector<double>(6, 0.0));
    c9->capacitance_per_meter.assign(6, std::vector<double>(6, 0.0));
    for (int b = 0; b < 3; ++b) {
        int o = b * 2;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j) {
                c9->inductance_per_meter[o + i][o + j] = blockL[i][j];
                c9->capacitance_per_meter[o + i][o + j] = blockC[i][j];
            }
    }
    for (int i = 0; i < 6; ++i)
        c9->resistance_per_meter[i][i] = 62.0e-3;
    c9->step_size.resize(7, 0.1);
    setYNegSegments(c9->segments, 10, 1, 7);
    c9->transfer_impedance.direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS;
    c9->transfer_impedance.resistive_term = 0.0;
    c9->transfer_impedance.inductive_term = 4.2e-9;
    c9->parent_cable = c8;
    c9->conductor_in_parent = 1;
    c9->n_segments = 7;

    expected.mtln->probes.resize(7);
    expected.mtln->probes[0].attached_to_cable = c1;
    expected.mtln->probes[0].index = 1;
    expected.mtln->probes[0].probe_type = PROBE_TYPE_VOLTAGE;
    expected.mtln->probes[0].probe_name = "b1_terminal_voltage";
    expected.mtln->probes[0].probe_position = {1, 7, 1};
    expected.mtln->probes[1].attached_to_cable = c1;
    expected.mtln->probes[1].index = 1;
    expected.mtln->probes[1].probe_type = PROBE_TYPE_CURRENT;
    expected.mtln->probes[1].probe_name = "b1_terminal_current";
    expected.mtln->probes[1].probe_position = {1, 7, 1};
    expected.mtln->probes[2].attached_to_cable = c1;
    expected.mtln->probes[2].index = 10;
    expected.mtln->probes[2].probe_type = PROBE_TYPE_CURRENT;
    expected.mtln->probes[2].probe_name = "junction_current";
    expected.mtln->probes[2].probe_position = {10, 7, 1};
    expected.mtln->probes[3].attached_to_cable = c4;
    expected.mtln->probes[3].index = 1;
    expected.mtln->probes[3].probe_type = PROBE_TYPE_CURRENT;
    expected.mtln->probes[3].probe_name = "junction_current";
    expected.mtln->probes[3].probe_position = {10, 7, 1};
    expected.mtln->probes[4].attached_to_cable = c7;
    expected.mtln->probes[4].index = 1;
    expected.mtln->probes[4].probe_type = PROBE_TYPE_CURRENT;
    expected.mtln->probes[4].probe_name = "junction_current";
    expected.mtln->probes[4].probe_position = {10, 7, 1};
    expected.mtln->probes[5].attached_to_cable = c4;
    expected.mtln->probes[5].index = 9;
    expected.mtln->probes[5].probe_type = PROBE_TYPE_CURRENT;
    expected.mtln->probes[5].probe_name = "b2_terminal_current";
    expected.mtln->probes[5].probe_position = {18, 7, 1};
    expected.mtln->probes[6].attached_to_cable = c7;
    expected.mtln->probes[6].index = 8;
    expected.mtln->probes[6].probe_type = PROBE_TYPE_CURRENT;
    expected.mtln->probes[6].probe_name = "b3_terminal_current";
    expected.mtln->probes[6].probe_position = {10, 0, 1};

    expected.mtln->networks.resize(4);

    // NETWORK 1
    expected.mtln->networks[0].connections.resize(10);
    for (int i = 0; i < 10; ++i)
        expected.mtln->networks[0].connections[i].nodes.resize(1);

    auto& n1 = expected.mtln->networks[0].connections;
    n1[0].nodes[0].conductor_in_cable = 1;
    n1[0].nodes[0].side = TERMINAL_NODE_SIDE_INI;
    n1[0].nodes[0].belongs_to_cable = c1;
    n1[0].nodes[0].termination.termination_type = TERMINATION_SERIES;
    n1[0].nodes[0].termination.resistance = 0.7e-3;

    n1[1].nodes[0].conductor_in_cable = 1;
    n1[1].nodes[0].side = TERMINAL_NODE_SIDE_INI;
    n1[1].nodes[0].belongs_to_cable = c2;
    n1[1].nodes[0].termination.termination_type = TERMINATION_SERIES;
    n1[1].nodes[0].termination.resistance = 1e-6;

    for (int i = 2; i < 10; ++i) {
        n1[i].nodes[0].side = TERMINAL_NODE_SIDE_INI;
        n1[i].nodes[0].belongs_to_cable = c3;
        n1[i].nodes[0].conductor_in_cable = i - 1;
    }
    n1[2].nodes[0].termination.termination_type = TERMINATION_RsLCp;
    n1[2].nodes[0].termination.resistance = 50.0;
    n1[2].nodes[0].termination.inductance = 30e-12;
    n1[2].nodes[0].termination.capacitance = 60e-9;
    n1[3].nodes[0].termination.termination_type = TERMINATION_SERIES;
    n1[3].nodes[0].termination.resistance = 1e10;
    n1[4].nodes[0].termination.termination_type = TERMINATION_RsLCp;
    n1[4].nodes[0].termination.resistance = 50.0;
    n1[4].nodes[0].termination.inductance = 30e-12;
    n1[4].nodes[0].termination.capacitance = 60e-9;
    n1[5].nodes[0].termination.termination_type = TERMINATION_SERIES;
    n1[5].nodes[0].termination.resistance = 1e10;
    n1[6].nodes[0].termination.termination_type = TERMINATION_RsLCp;
    n1[6].nodes[0].termination.resistance = 50.0;
    n1[6].nodes[0].termination.inductance = 30e-12;
    n1[6].nodes[0].termination.capacitance = 60e-9;
    n1[7].nodes[0].termination.termination_type = TERMINATION_SERIES;
    n1[7].nodes[0].termination.resistance = 1e10;
    n1[8].nodes[0].termination.termination_type = TERMINATION_RsLCp;
    n1[8].nodes[0].termination.resistance = 50.0;
    n1[8].nodes[0].termination.inductance = 30e-12;
    n1[8].nodes[0].termination.capacitance = 60e-9;
    n1[9].nodes[0].termination.termination_type = TERMINATION_SERIES;
    n1[9].nodes[0].termination.resistance = 1e10;

    // NETWORK 2
    expected.mtln->networks[1].connections.resize(10);
    expected.mtln->networks[1].connections[0].nodes.resize(3);
    expected.mtln->networks[1].connections[1].nodes.resize(3);
    for (int i = 2; i < 10; ++i)
        expected.mtln->networks[1].connections[i].nodes.resize(2);

    auto& n2 = expected.mtln->networks[1].connections;
    n2[0].nodes[0].conductor_in_cable = 1;
    n2[0].nodes[0].side = TERMINAL_NODE_SIDE_END;
    n2[0].nodes[0].belongs_to_cable = c1;
    n2[0].nodes[0].termination.termination_type = TERMINATION_SERIES;
    n2[0].nodes[0].termination.resistance = 1e-6;
    n2[0].nodes[1].conductor_in_cable = 1;
    n2[0].nodes[1].side = TERMINAL_NODE_SIDE_INI;
    n2[0].nodes[1].belongs_to_cable = c4;
    n2[0].nodes[1].termination.termination_type = TERMINATION_SERIES;
    n2[0].nodes[1].termination.resistance = 1e-6;
    n2[0].nodes[2].conductor_in_cable = 1;
    n2[0].nodes[2].side = TERMINAL_NODE_SIDE_INI;
    n2[0].nodes[2].belongs_to_cable = c7;
    n2[0].nodes[2].termination.termination_type = TERMINATION_SERIES;
    n2[0].nodes[2].termination.resistance = 1e-6;
    n2[1].nodes[0].conductor_in_cable = 1;
    n2[1].nodes[0].side = TERMINAL_NODE_SIDE_END;
    n2[1].nodes[0].belongs_to_cable = c2;
    n2[1].nodes[0].termination.termination_type = TERMINATION_SERIES;
    n2[1].nodes[0].termination.resistance = 1e-6;
    n2[1].nodes[1].conductor_in_cable = 1;
    n2[1].nodes[1].side = TERMINAL_NODE_SIDE_INI;
    n2[1].nodes[1].belongs_to_cable = c5;
    n2[1].nodes[1].termination.termination_type = TERMINATION_SERIES;
    n2[1].nodes[1].termination.resistance = 1e-6;
    n2[1].nodes[2].conductor_in_cable = 1;
    n2[1].nodes[2].side = TERMINAL_NODE_SIDE_INI;
    n2[1].nodes[2].belongs_to_cable = c8;
    n2[1].nodes[2].termination.termination_type = TERMINATION_SERIES;
    n2[1].nodes[2].termination.resistance = 1e-6;

    for (int i = 2; i < 8; ++i) {
        int fi = i + 1;
        n2[i].nodes[0].conductor_in_cable = fi - 2;
        n2[i].nodes[0].side = TERMINAL_NODE_SIDE_END;
        n2[i].nodes[0].belongs_to_cable = c3;
        n2[i].nodes[0].termination.termination_type = TERMINATION_SERIES;
        n2[i].nodes[0].termination.resistance = 1e-6;
        n2[i].nodes[1].conductor_in_cable = fi - 2;
        n2[i].nodes[1].side = TERMINAL_NODE_SIDE_INI;
        n2[i].nodes[1].belongs_to_cable = c9;
        n2[i].nodes[1].termination.termination_type = TERMINATION_SERIES;
        n2[i].nodes[1].termination.resistance = 1e-6;
    }
    for (int i = 8; i < 10; ++i) {
        int fi = i + 1;
        n2[i].nodes[0].conductor_in_cable = fi - 2;
        n2[i].nodes[0].side = TERMINAL_NODE_SIDE_END;
        n2[i].nodes[0].belongs_to_cable = c3;
        n2[i].nodes[0].termination.termination_type = TERMINATION_SERIES;
        n2[i].nodes[0].termination.resistance = 1e-6;
        n2[i].nodes[1].conductor_in_cable = fi - 8;
        n2[i].nodes[1].side = TERMINAL_NODE_SIDE_INI;
        n2[i].nodes[1].belongs_to_cable = c6;
        n2[i].nodes[1].termination.termination_type = TERMINATION_SERIES;
        n2[i].nodes[1].termination.resistance = 1e-6;
    }

    // NETWORK 3
    expected.mtln->networks[2].connections.resize(4);
    for (int i = 0; i < 4; ++i)
        expected.mtln->networks[2].connections[i].nodes.resize(1);
    auto& n3 = expected.mtln->networks[2].connections;
    n3[0].nodes[0].conductor_in_cable = 1;
    n3[0].nodes[0].side = TERMINAL_NODE_SIDE_END;
    n3[0].nodes[0].belongs_to_cable = c4;
    n3[0].nodes[0].termination.termination_type = TERMINATION_SERIES;
    n3[0].nodes[0].termination.resistance = 1;
    n3[1].nodes[0].conductor_in_cable = 1;
    n3[1].nodes[0].side = TERMINAL_NODE_SIDE_END;
    n3[1].nodes[0].belongs_to_cable = c5;
    n3[1].nodes[0].termination.termination_type = TERMINATION_SERIES;
    n3[1].nodes[0].termination.resistance = 1e-6;
    for (int i = 2; i < 4; ++i) {
        n3[i].nodes[0].side = TERMINAL_NODE_SIDE_END;
        n3[i].nodes[0].belongs_to_cable = c6;
        n3[i].nodes[0].conductor_in_cable = i - 1;
    }
    n3[2].nodes[0].termination.termination_type = TERMINATION_SERIES;
    n3[2].nodes[0].termination.resistance = 50;
    n3[3].nodes[0].termination.termination_type = TERMINATION_SERIES;
    n3[3].nodes[0].termination.resistance = 50;

    // NETWORK 4
    expected.mtln->networks[3].connections.resize(8);
    for (int i = 0; i < 8; ++i)
        expected.mtln->networks[3].connections[i].nodes.resize(1);
    auto& n4 = expected.mtln->networks[3].connections;
    n4[0].nodes[0].conductor_in_cable = 1;
    n4[0].nodes[0].side = TERMINAL_NODE_SIDE_END;
    n4[0].nodes[0].belongs_to_cable = c7;
    n4[0].nodes[0].termination.termination_type = TERMINATION_SERIES;
    n4[0].nodes[0].termination.resistance = 0.7e-3;
    n4[1].nodes[0].conductor_in_cable = 1;
    n4[1].nodes[0].side = TERMINAL_NODE_SIDE_END;
    n4[1].nodes[0].belongs_to_cable = c8;
    n4[1].nodes[0].termination.termination_type = TERMINATION_SERIES;
    n4[1].nodes[0].termination.resistance = 1e-6;
    for (int i = 2; i < 8; ++i) {
        n4[i].nodes[0].side = TERMINAL_NODE_SIDE_END;
        n4[i].nodes[0].belongs_to_cable = c9;
        n4[i].nodes[0].conductor_in_cable = i - 1;
        n4[i].nodes[0].termination.termination_type = TERMINATION_SERIES;
        n4[i].nodes[0].termination.resistance = 50;
    }

    return expected;
}

#endif // CompileWithMTLN

#endif // TEST_READ_MTLN_H
