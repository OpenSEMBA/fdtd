#ifndef TEST_READ_TOWELHANGER_H
#define TEST_READ_TOWELHANGER_H

#ifdef CompileWithMTLN

#include "test_smbjson_helpers.h"

using namespace mtln_types_m;

inline Parseador_t expectedReadTowelHanger() {
    Parseador_t expected = buildExpected();

    expected.general->dt = 1e-12;
    expected.general->nmax = 2000;

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
    expected.pecRegs->Surfs[0].tag = "copper@ground_plane";

    expected.mtln->time_step = 1e-12;
    expected.mtln->number_of_steps = 2000;

    expected.mtln->cables.resize(1);
    expected.mtln->cables[0].ptr = std::make_unique<unshielded_multiwire_t>();
    auto* wire = static_cast<unshielded_multiwire_t*>(expected.mtln->cables[0].ptr.get());
    wire->name = "wire";
    initializeCablePULParameters(wire);
    wire->step_size.resize(20, 0.01);
    wire->segments.resize(20);

    for (int i = 0; i < 2; ++i) {
        wire->segments[i].x = 27;
        wire->segments[i].y = 25;
        wire->segments[i].z = 29 + i + 1;
        wire->segments[i].orientation = DIRECTION_Z_POS;
    }
    for (int i = 2; i < 10; ++i) {
        wire->segments[i].x = 24 + i + 1;
        wire->segments[i].y = 25;
        wire->segments[i].z = 32;
        wire->segments[i].orientation = DIRECTION_X_POS;
    }
    for (int i = 0; i < 8; ++i) {
        wire->segments[10 + i].x = 34 + i + 1;
        wire->segments[10 + i].y = 25;
        wire->segments[10 + i].z = 32;
        wire->segments[10 + i].orientation = DIRECTION_X_POS;
    }
    for (int i = 8; i < 10; ++i) {
        wire->segments[10 + i].x = 43;
        wire->segments[10 + i].y = 25;
        wire->segments[10 + i].z = 40 - (i + 1);
        wire->segments[10 + i].orientation = DIRECTION_Z_NEG;
    }
    wire->initial_connector = nullptr;
    wire->end_connector = nullptr;

    expected.mtln->probes.resize(3);
    expected.mtln->probes[0].attached_to_cable = wire;
    expected.mtln->probes[0].index = 1;
    expected.mtln->probes[0].probe_type = PROBE_TYPE_CURRENT;
    expected.mtln->probes[0].probe_name = "wire_start";
    expected.mtln->probes[0].probe_position = {27, 25, 30};

    expected.mtln->probes[1].attached_to_cable = wire;
    expected.mtln->probes[1].index = 21;
    expected.mtln->probes[1].probe_type = PROBE_TYPE_CURRENT;
    expected.mtln->probes[1].probe_name = "wire_end";
    expected.mtln->probes[1].probe_position = {43, 25, 30};

    expected.mtln->probes[2].attached_to_cable = wire;
    expected.mtln->probes[2].index = 11;
    expected.mtln->probes[2].probe_type = PROBE_TYPE_CURRENT;
    expected.mtln->probes[2].probe_name = "wire_mid";
    expected.mtln->probes[2].probe_position = {35, 25, 32};

    expected.mtln->networks.resize(2);

    expected.mtln->networks[0].connections.resize(1);
    expected.mtln->networks[0].connections[0].nodes.resize(1);
    auto& n0 = expected.mtln->networks[0].connections[0].nodes[0];
    n0.conductor_in_cable = 1;
    n0.side = TERMINAL_NODE_SIDE_INI;
    n0.belongs_to_cable = wire;
    n0.termination.termination_type = TERMINATION_SERIES;
    n0.termination.resistance = 50.0;
    n0.termination.source.path_to_excitation = "towelHanger.exc";
    n0.termination.source.source_type = SOURCE_TYPE_VOLTAGE;

    expected.mtln->networks[1].connections.resize(1);
    expected.mtln->networks[1].connections[0].nodes.resize(1);
    auto& n1 = expected.mtln->networks[1].connections[0].nodes[0];
    n1.conductor_in_cable = 1;
    n1.side = TERMINAL_NODE_SIDE_END;
    n1.belongs_to_cable = wire;
    n1.termination.termination_type = TERMINATION_SHORT;

    return expected;
}

#endif // CompileWithMTLN

#endif // TEST_READ_TOWELHANGER_H
