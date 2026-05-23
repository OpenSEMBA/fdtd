#ifndef TEST_READ_HOLLAND1981_UNSHIELDED_H
#define TEST_READ_HOLLAND1981_UNSHIELDED_H

#ifdef CompileWithMTLN

#include "test_smbjson_helpers.h"

using namespace mtln_types_m;

inline Parseador_t expectedReadHolland1981Unshielded() {
    Parseador_t expected = buildExpected();

    expected.general->dt = 30e-12;
    expected.general->nmax = 1000;

    expected.matriz->totalX = 21;
    expected.matriz->totalY = 21;
    expected.matriz->totalZ = 23;

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
    expected.despl->mx2 = 20;
    expected.despl->my1 = 0;
    expected.despl->my2 = 20;
    expected.despl->mz1 = 0;
    expected.despl->mz2 = 22;

    for (int i = 0; i < 6; ++i) {
        expected.front->tipoFrontera[i] = F_PML;
        expected.front->propiedadesPML[i].numCapas = 6;
        expected.front->propiedadesPML[i].orden = 2.0;
        expected.front->propiedadesPML[i].refl = 0.001;
    }

    expected.plnSrc->collection.resize(1);
    expected.plnSrc->collection[0].nombre_fichero = "holland.exc";
    expected.plnSrc->collection[0].atributo = "LOCKED";
    expected.plnSrc->collection[0].coor1[0] = 1;
    expected.plnSrc->collection[0].coor1[1] = 1;
    expected.plnSrc->collection[0].coor1[2] = 1;
    expected.plnSrc->collection[0].coor2[0] = 18;
    expected.plnSrc->collection[0].coor2[1] = 18;
    expected.plnSrc->collection[0].coor2[2] = 20;
    expected.plnSrc->collection[0].theta = 1.5708;
    expected.plnSrc->collection[0].phi = 0.0;
    expected.plnSrc->collection[0].alpha = 0.0;
    expected.plnSrc->collection[0].beta = 0.0;
    expected.plnSrc->collection[0].isRC = false;
    expected.plnSrc->collection[0].numModes = 1;
    expected.plnSrc->collection[0].INCERTMAX = 0.0;
    expected.plnSrc->nc = 1;
    expected.plnSrc->nC_max = 1;

    expected.mtln->time_step = 30e-12;
    expected.mtln->number_of_steps = 1000;

    expected.mtln->cables.resize(1);
    expected.mtln->cables[0].ptr = std::make_unique<unshielded_multiwire_t>();
    auto* cable = static_cast<unshielded_multiwire_t*>(expected.mtln->cables[0].ptr.get());
    cable->name = "single_unshielded_multiwire";
    initializeCablePULParameters(cable);
    cable->step_size.resize(10, 0.1);
    cable->segments.resize(10);
    for (int i = 0; i < 10; ++i) {
        cable->segments[i].x = 11;
        cable->segments[i].y = 11;
        cable->segments[i].z = i + 7;
        cable->segments[i].orientation = DIRECTION_Z_POS;
    }
    cable->initial_connector = nullptr;
    cable->end_connector = nullptr;

    expected.mtln->probes.resize(1);
    expected.mtln->probes[0].attached_to_cable = cable;
    expected.mtln->probes[0].index = 6;
    expected.mtln->probes[0].probe_type = PROBE_TYPE_CURRENT;
    expected.mtln->probes[0].probe_name = "mid_point";
    expected.mtln->probes[0].probe_position = {11, 11, 12};

    expected.mtln->networks.resize(2);

    expected.mtln->networks[0].connections.resize(1);
    expected.mtln->networks[0].connections[0].nodes.resize(1);
    auto& n0 = expected.mtln->networks[0].connections[0].nodes[0];
    n0.conductor_in_cable = 1;
    n0.side = TERMINAL_NODE_SIDE_INI;
    n0.belongs_to_cable = cable;
    n0.termination.termination_type = TERMINATION_OPEN;

    expected.mtln->networks[1].connections.resize(1);
    expected.mtln->networks[1].connections[0].nodes.resize(1);
    auto& n1 = expected.mtln->networks[1].connections[0].nodes[0];
    n1.conductor_in_cable = 1;
    n1.side = TERMINAL_NODE_SIDE_END;
    n1.belongs_to_cable = cable;
    n1.termination.termination_type = TERMINATION_OPEN;

    return expected;
}

#endif // CompileWithMTLN

#endif // TEST_READ_HOLLAND1981_UNSHIELDED_H
