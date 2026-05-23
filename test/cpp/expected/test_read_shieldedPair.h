#ifndef TEST_READ_SHIELDEDPAIR_H
#define TEST_READ_SHIELDEDPAIR_H

#ifdef CompileWithMTLN

#include "test_smbjson_helpers.h"

using namespace mtln_types_m;

inline Parseador_t expectedReadShieldedPair() {
    Parseador_t expected = buildExpected();

    expected.general->dt = 0.43e-10;
    expected.general->nmax = 700;

    expected.matriz->totalX = 151;
    expected.matriz->totalY = 151;
    expected.matriz->totalZ = 151;

    expected.despl->nX = 1;
    expected.despl->nY = 1;
    expected.despl->nZ = 1;
    expected.despl->desX.resize(1);
    expected.despl->desY.resize(1);
    expected.despl->desZ.resize(1);
    expected.despl->desX[0] = 0.180;
    expected.despl->desY[0] = 0.180;
    expected.despl->desZ[0] = 0.0504;
    expected.despl->mx1 = 0;
    expected.despl->mx2 = 150;
    expected.despl->my1 = 0;
    expected.despl->my2 = 150;
    expected.despl->mz1 = 0;
    expected.despl->mz2 = 150;

    for (int i = 0; i < 6; ++i) {
        expected.front->tipoFrontera[i] = F_PML;
        expected.front->propiedadesPML[i].numCapas = 6;
        expected.front->propiedadesPML[i].orden = 2.0;
        expected.front->propiedadesPML[i].refl = 0.0001;
    }

    expected.pecRegs->nLins = 0;
    expected.pecRegs->nLins_max = 0;
    expected.pecRegs->nSurfs = 1;
    expected.pecRegs->nSurfs_max = 1;
    expected.pecRegs->nVols = 0;
    expected.pecRegs->nVols_max = 0;
    expected.pecRegs->Surfs.resize(1);
    expected.pecRegs->Surfs[0].Xi = 20;
    expected.pecRegs->Surfs[0].Xe = 129;
    expected.pecRegs->Surfs[0].Yi = 20;
    expected.pecRegs->Surfs[0].Ye = 129;
    expected.pecRegs->Surfs[0].Zi = 74;
    expected.pecRegs->Surfs[0].Ze = 74;
    expected.pecRegs->Surfs[0].Xtrancos = 1;
    expected.pecRegs->Surfs[0].Ytrancos = 1;
    expected.pecRegs->Surfs[0].Ztrancos = 1;
    expected.pecRegs->Surfs[0].Or = 3;
    expected.pecRegs->Surfs[0].tag = "material5@layer5";

    expected.plnSrc->collection.resize(1);
    expected.plnSrc->collection[0].nombre_fichero = "shielded_pair.exc";
    expected.plnSrc->collection[0].atributo = "LOCKED";
    expected.plnSrc->collection[0].coor1[0] = 10;
    expected.plnSrc->collection[0].coor1[1] = 10;
    expected.plnSrc->collection[0].coor1[2] = 10;
    expected.plnSrc->collection[0].coor2[0] = 139;
    expected.plnSrc->collection[0].coor2[1] = 139;
    expected.plnSrc->collection[0].coor2[2] = 139;
    expected.plnSrc->collection[0].theta = 3.1416;
    expected.plnSrc->collection[0].phi = 0.0;
    expected.plnSrc->collection[0].alpha = 1.5708;
    expected.plnSrc->collection[0].beta = -1.5708;
    expected.plnSrc->collection[0].isRC = false;
    expected.plnSrc->collection[0].numModes = 1;
    expected.plnSrc->collection[0].INCERTMAX = 0.0;
    expected.plnSrc->nc = 1;
    expected.plnSrc->nC_max = 1;

    expected.mtln->time_step = 0.43e-10;
    expected.mtln->number_of_steps = 700;

    expected.mtln->cables.resize(2);

    expected.mtln->cables[0].ptr = std::make_unique<unshielded_multiwire_t>();
    auto* line0 = static_cast<unshielded_multiwire_t*>(expected.mtln->cables[0].ptr.get());
    initializeCablePULParameters(line0, 1);
    line0->name = "line_0";
    line0->resistance_per_meter = {{22.9e-3}};
    line0->step_size = {0.0504, 0.180, 0.180, 0.180, 0.0504};
    line0->segments.resize(5);
    line0->segments[0].x = 75;
    line0->segments[0].y = 71;
    line0->segments[0].z = 74;
    line0->segments[0].orientation = DIRECTION_Z_POS;
    for (int i = 1; i < 4; ++i) {
        line0->segments[i].x = 75;
        line0->segments[i].y = 69 + i + 1;
        line0->segments[i].z = 75;
        line0->segments[i].orientation = DIRECTION_Y_POS;
    }
    line0->segments[4].x = 75;
    line0->segments[4].y = 74;
    line0->segments[4].z = 74;
    line0->segments[4].orientation = DIRECTION_Z_NEG;
    line0->initial_connector = nullptr;
    line0->end_connector = nullptr;

    expected.mtln->cables[1].ptr = std::make_unique<shielded_multiwire_t>();
    auto* line1 = static_cast<shielded_multiwire_t*>(expected.mtln->cables[1].ptr.get());
    initializeCablePULParameters(line1, 2);
    line1->name = "line_1";
    line1->inductance_per_meter = {
        {3.13182309e-07, 7.45674981e-08},
        {7.45674981e-08, 3.13182309e-07}};
    line1->capacitance_per_meter = {
        {85.0e-12, -20.5e-12},
        {-20.5e-12, 85.0e-12}};
    line1->step_size = {0.0504, 0.180, 0.180, 0.180, 0.0504};
    line1->segments.resize(5);
    line1->segments[0].x = 75;
    line1->segments[0].y = 71;
    line1->segments[0].z = 74;
    line1->segments[0].orientation = DIRECTION_Z_POS;
    for (int i = 1; i < 4; ++i) {
        line1->segments[i].x = 75;
        line1->segments[i].y = 69 + i + 1;
        line1->segments[i].z = 75;
        line1->segments[i].orientation = DIRECTION_Y_POS;
    }
    line1->segments[4].x = 75;
    line1->segments[4].y = 74;
    line1->segments[4].z = 74;
    line1->segments[4].orientation = DIRECTION_Z_NEG;
    line1->transfer_impedance.direction = TRANSFER_IMPEDANCE_DIRECTION_BOTH;
    line1->transfer_impedance.resistive_term = 0.0;
    line1->transfer_impedance.inductive_term = 4.0e-9;
    line1->transfer_impedance.poles.clear();
    line1->transfer_impedance.residues.clear();
    line1->parent_cable = line0;
    line1->conductor_in_parent = 1;
    line1->initial_connector = nullptr;
    line1->end_connector = nullptr;

    expected.mtln->probes.resize(4);
    expected.mtln->probes[0].attached_to_cable = line0;
    expected.mtln->probes[0].index = 1;
    expected.mtln->probes[0].probe_type = PROBE_TYPE_CURRENT;
    expected.mtln->probes[0].probe_name = "wire_end";
    expected.mtln->probes[0].probe_position = {75, 71, 74};

    expected.mtln->probes[1].attached_to_cable = line0;
    expected.mtln->probes[1].index = 1;
    expected.mtln->probes[1].probe_type = PROBE_TYPE_VOLTAGE;
    expected.mtln->probes[1].probe_name = "wire_end";
    expected.mtln->probes[1].probe_position = {75, 71, 74};

    expected.mtln->probes[2].attached_to_cable = line0;
    expected.mtln->probes[2].index = 6;
    expected.mtln->probes[2].probe_type = PROBE_TYPE_CURRENT;
    expected.mtln->probes[2].probe_name = "wire_start";
    expected.mtln->probes[2].probe_position = {75, 74, 74};

    expected.mtln->probes[3].attached_to_cable = line0;
    expected.mtln->probes[3].index = 6;
    expected.mtln->probes[3].probe_type = PROBE_TYPE_VOLTAGE;
    expected.mtln->probes[3].probe_name = "wire_start";
    expected.mtln->probes[3].probe_position = {75, 74, 74};

    expected.mtln->networks.resize(2);

    expected.mtln->networks[0].connections.resize(3);
    for (int c = 0; c < 3; ++c) {
        expected.mtln->networks[0].connections[c].nodes.resize(1);
        auto& node = expected.mtln->networks[0].connections[c].nodes[0];
        node.side = TERMINAL_NODE_SIDE_INI;
        node.termination.termination_type = TERMINATION_SERIES;
        node.termination.resistance = 50.0;
    }
    expected.mtln->networks[0].connections[0].nodes[0].conductor_in_cable = 1;
    expected.mtln->networks[0].connections[0].nodes[0].belongs_to_cable = line0;
    expected.mtln->networks[0].connections[1].nodes[0].conductor_in_cable = 1;
    expected.mtln->networks[0].connections[1].nodes[0].belongs_to_cable = line1;
    expected.mtln->networks[0].connections[2].nodes[0].conductor_in_cable = 2;
    expected.mtln->networks[0].connections[2].nodes[0].belongs_to_cable = line1;

    expected.mtln->networks[1].connections.resize(3);
    for (int c = 0; c < 3; ++c) {
        expected.mtln->networks[1].connections[c].nodes.resize(1);
        auto& node = expected.mtln->networks[1].connections[c].nodes[0];
        node.side = TERMINAL_NODE_SIDE_END;
        node.termination.termination_type = TERMINATION_SERIES;
        node.termination.resistance = 50.0;
    }
    expected.mtln->networks[1].connections[0].nodes[0].conductor_in_cable = 1;
    expected.mtln->networks[1].connections[0].nodes[0].belongs_to_cable = line0;
    expected.mtln->networks[1].connections[1].nodes[0].conductor_in_cable = 1;
    expected.mtln->networks[1].connections[1].nodes[0].belongs_to_cable = line1;
    expected.mtln->networks[1].connections[2].nodes[0].conductor_in_cable = 2;
    expected.mtln->networks[1].connections[2].nodes[0].belongs_to_cable = line1;

    expected.mtln->connectors.clear();

    return expected;
}

#endif // CompileWithMTLN

#endif // TEST_READ_SHIELDEDPAIR_H
