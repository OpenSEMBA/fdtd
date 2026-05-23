#ifndef TEST_READ_CONNECTEDWIRES_H
#define TEST_READ_CONNECTEDWIRES_H

#ifdef CompileWithMTLN

#include "test_smbjson_helpers.h"

using namespace mtln_types_m;

inline Parseador_t expectedReadConnectedWires() {
    Parseador_t expected = buildExpected();

    expected.general->dt = 1e-12;
    expected.general->nmax = 1000;

    expected.matriz->totalX = 61;
    expected.matriz->totalY = 61;
    expected.matriz->totalZ = 61;

    expected.despl->nX = 1;
    expected.despl->nY = 1;
    expected.despl->nZ = 1;
    expected.despl->desX.resize(1);
    expected.despl->desY.resize(1);
    expected.despl->desZ.resize(1);
    expected.despl->desX[0] = 0.01;
    expected.despl->desY[0] = 0.01;
    expected.despl->desZ[0] = 0.01;
    expected.despl->mx1 = 0;
    expected.despl->mx2 = 60;
    expected.despl->my1 = 0;
    expected.despl->my2 = 60;
    expected.despl->mz1 = 0;
    expected.despl->mz2 = 60;

    for (int i = 0; i < 6; ++i) {
        expected.front->tipoFrontera[i] = F_PML;
        expected.front->propiedadesPML[i].numCapas = 6;
        expected.front->propiedadesPML[i].orden = 2.0;
        expected.front->propiedadesPML[i].refl = 0.001;
    }

    expected.pecRegs->nLins = 0;
    expected.pecRegs->nLins_max = 0;
    expected.pecRegs->nSurfs = 1;
    expected.pecRegs->nSurfs_max = 1;
    expected.pecRegs->nVols = 0;
    expected.pecRegs->nVols_max = 0;
    expected.pecRegs->Surfs.resize(1);
    expected.pecRegs->Surfs[0].Xi = 25;
    expected.pecRegs->Surfs[0].Xe = 44;
    expected.pecRegs->Surfs[0].Yi = 20;
    expected.pecRegs->Surfs[0].Ye = 29;
    expected.pecRegs->Surfs[0].Zi = 30;
    expected.pecRegs->Surfs[0].Ze = 30;
    expected.pecRegs->Surfs[0].Xtrancos = 1;
    expected.pecRegs->Surfs[0].Ytrancos = 1;
    expected.pecRegs->Surfs[0].Ztrancos = 1;
    expected.pecRegs->Surfs[0].Or = 3;
    expected.pecRegs->Surfs[0].tag = "aluminum@ground_plane";

    expected.mtln->time_step = 1e-12;
    expected.mtln->number_of_steps = 1000;

    expected.mtln->cables.resize(2);
    expected.mtln->cables[0].ptr = std::make_unique<unshielded_multiwire_t>();
    expected.mtln->cables[1].ptr = std::make_unique<unshielded_multiwire_t>();
    auto* cable1 = static_cast<unshielded_multiwire_t*>(expected.mtln->cables[0].ptr.get());
    auto* cable2 = static_cast<unshielded_multiwire_t*>(expected.mtln->cables[1].ptr.get());

    cable1->name = "cable1";
    initializeCablePULParameters(cable1);
    cable1->step_size.resize(10, 0.01);
    cable1->segments.resize(10);
    for (int i = 0; i < 2; ++i) {
        cable1->segments[i].x = 27;
        cable1->segments[i].y = 25;
        cable1->segments[i].z = 29 + i + 1;
        cable1->segments[i].orientation = DIRECTION_Z_POS;
    }
    for (int i = 2; i < 10; ++i) {
        cable1->segments[i].x = 24 + i + 1;
        cable1->segments[i].y = 25;
        cable1->segments[i].z = 32;
        cable1->segments[i].orientation = DIRECTION_X_POS;
    }
    cable1->initial_connector = nullptr;
    cable1->end_connector = nullptr;

    cable2->name = "cable2";
    initializeCablePULParameters(cable2);
    cable2->step_size.resize(10, 0.01);
    cable2->segments.resize(10);
    for (int i = 0; i < 8; ++i) {
        cable2->segments[i].x = 34 + i + 1;
        cable2->segments[i].y = 25;
        cable2->segments[i].z = 32;
        cable2->segments[i].orientation = DIRECTION_X_POS;
    }
    for (int i = 8; i < 10; ++i) {
        cable2->segments[i].x = 43;
        cable2->segments[i].y = 25;
        cable2->segments[i].z = 40 - (i + 1);
        cable2->segments[i].orientation = DIRECTION_Z_NEG;
    }
    cable2->initial_connector = nullptr;
    cable2->end_connector = nullptr;

    expected.mtln->probes.resize(2);
    expected.mtln->probes[0].attached_to_cable = cable1;
    expected.mtln->probes[0].index = 1;
    expected.mtln->probes[0].probe_type = PROBE_TYPE_CURRENT;
    expected.mtln->probes[0].probe_name = "wire_start";
    expected.mtln->probes[0].probe_position = {27, 25, 30};

    expected.mtln->probes[1].attached_to_cable = cable2;
    expected.mtln->probes[1].index = 11;
    expected.mtln->probes[1].probe_type = PROBE_TYPE_CURRENT;
    expected.mtln->probes[1].probe_name = "wire_end";
    expected.mtln->probes[1].probe_position = {43, 25, 30};

    expected.mtln->networks.resize(3);

    expected.mtln->networks[0].connections.resize(1);
    expected.mtln->networks[0].connections[0].nodes.resize(1);
    auto& n00 = expected.mtln->networks[0].connections[0].nodes[0];
    n00.conductor_in_cable = 1;
    n00.side = TERMINAL_NODE_SIDE_INI;
    n00.belongs_to_cable = cable1;
    n00.termination.termination_type = TERMINATION_SHORT;
    n00.termination.source.path_to_excitation = "ramp.exc";
    n00.termination.source.source_type = SOURCE_TYPE_VOLTAGE;

    expected.mtln->networks[1].connections.resize(1);
    expected.mtln->networks[1].connections[0].nodes.resize(2);
    auto& n10 = expected.mtln->networks[1].connections[0].nodes[0];
    n10.conductor_in_cable = 1;
    n10.side = TERMINAL_NODE_SIDE_END;
    n10.belongs_to_cable = cable1;
    n10.termination.termination_type = TERMINATION_SERIES;
    n10.termination.resistance = 50;
    auto& n11 = expected.mtln->networks[1].connections[0].nodes[1];
    n11.conductor_in_cable = 1;
    n11.side = TERMINAL_NODE_SIDE_INI;
    n11.belongs_to_cable = cable2;
    n11.termination.termination_type = TERMINATION_SHORT;

    expected.mtln->networks[2].connections.resize(1);
    expected.mtln->networks[2].connections[0].nodes.resize(1);
    auto& n20 = expected.mtln->networks[2].connections[0].nodes[0];
    n20.conductor_in_cable = 1;
    n20.side = TERMINAL_NODE_SIDE_END;
    n20.belongs_to_cable = cable2;
    n20.termination.termination_type = TERMINATION_SHORT;

    return expected;
}

#endif // CompileWithMTLN

#endif // TEST_READ_CONNECTEDWIRES_H
