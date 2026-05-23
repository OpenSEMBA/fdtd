#ifndef TEST_READ_UNSHIELDED_MULTIWIRES_MULTIPOLAR_EXPANSION_H
#define TEST_READ_UNSHIELDED_MULTIWIRES_MULTIPOLAR_EXPANSION_H

#ifdef CompileWithMTLN

#include "test_smbjson_helpers.h"

using namespace mtln_types_m;

inline Parseador_t expectedReadUnshieldedMultiwiresMultipolarExpansion() {
    Parseador_t expected = buildExpected();

    expected.general->dt = 3e-11;
    expected.general->nmax = 1100;

    expected.matriz->totalX = 31;
    expected.matriz->totalY = 31;
    expected.matriz->totalZ = 31;

    expected.despl->nX = 1;
    expected.despl->nY = 1;
    expected.despl->nZ = 1;
    expected.despl->desX.resize(1);
    expected.despl->desY.resize(1);
    expected.despl->desZ.resize(1);
    expected.despl->desX[0] = 0.2;
    expected.despl->desY[0] = 0.2;
    expected.despl->desZ[0] = 0.2;
    expected.despl->mx1 = 0;
    expected.despl->mx2 = 30;
    expected.despl->my1 = 0;
    expected.despl->my2 = 30;
    expected.despl->mz1 = 0;
    expected.despl->mz2 = 30;

    for (int i = 0; i < 6; ++i) {
        expected.front->tipoFrontera[i] = F_PML;
        expected.front->propiedadesPML[i].numCapas = 8;
        expected.front->propiedadesPML[i].orden = 2;
        expected.front->propiedadesPML[i].refl = 0.001;
    }

    expected.pecRegs->nSurfs = 1;
    expected.pecRegs->nSurfs_max = 1;
    expected.pecRegs->Surfs.resize(1);

    expected.plnSrc->collection.resize(1);
    expected.plnSrc->collection[0].nombre_fichero = "unshielded_50ns.exc";
    expected.plnSrc->collection[0].atributo = "LOCKED";
    expected.plnSrc->collection[0].coor1[0] = 1;
    expected.plnSrc->collection[0].coor1[1] = 1;
    expected.plnSrc->collection[0].coor1[2] = 1;
    expected.plnSrc->collection[0].coor2[0] = 28;
    expected.plnSrc->collection[0].coor2[1] = 28;
    expected.plnSrc->collection[0].coor2[2] = 28;
    expected.plnSrc->collection[0].theta = 1.5708;
    expected.plnSrc->collection[0].phi = 0.0;
    expected.plnSrc->collection[0].alpha = 0.0;
    expected.plnSrc->collection[0].beta = 0.0;
    expected.plnSrc->collection[0].isRC = false;
    expected.plnSrc->collection[0].numModes = 1;
    expected.plnSrc->collection[0].INCERTMAX = 0.0;
    expected.plnSrc->nc = 1;
    expected.plnSrc->nC_max = 1;

    expected.mtln->time_step = 3e-11;
    expected.mtln->number_of_steps = 1100;

    expected.mtln->cables.resize(1);
    expected.mtln->cables[0].ptr = std::make_unique<unshielded_multiwire_t>();
    auto* cable = static_cast<unshielded_multiwire_t*>(expected.mtln->cables[0].ptr.get());
    initializeCablePULParameters(cable, 2);
    cable->name = "unshielded_pair";

    cable->multipolar_expansion.resize(1);
    auto& me = cable->multipolar_expansion[0];
    me.inner_region.min = {-0.0265000002, -0.0310000002};
    me.inner_region.max = {0.03550000020000001, 0.0310000002};

    me.electric.resize(2);
    me.electric[0].conductor_potentials = {1.0, 0.5909272203987278};
    me.electric[0].expansion_center = {-0.004970886788455953, 6.610694092023349e-07};
    me.electric[0].inner_region_average_potential = 0.5608636261599323;
    me.electric[0].ab.resize(1);
    me.electric[0].ab[0].a = 0.9488836986256424;
    me.electric[0].ab[0].b = 0.0;

    me.electric[1].conductor_potentials = {0.8497110567446987, 1.0};
    me.electric[1].expansion_center = {0.009920513440028656, 6.949869591535922e-07};
    me.electric[1].inner_region_average_potential = 0.8070848243572611;
    me.electric[1].ab.resize(1);
    me.electric[1].ab[0].a = 1.3644011168458479;
    me.electric[1].ab[0].b = 0.0;

    me.magnetic = me.electric;

    cable->step_size.resize(15, 0.2);
    cable->segments.resize(15);
    for (int i = 0; i < 15; ++i) {
        cable->segments[i].x = 2;
        cable->segments[i].y = 11;
        cable->segments[i].z = 6 + i + 1;
        cable->segments[i].orientation = DIRECTION_Z_POS;
    }
    cable->initial_connector = nullptr;
    cable->end_connector = nullptr;

    expected.mtln->probes.resize(1);
    expected.mtln->probes[0].attached_to_cable = cable;
    expected.mtln->probes[0].index = 8;
    expected.mtln->probes[0].probe_type = PROBE_TYPE_CURRENT;
    expected.mtln->probes[0].probe_name = "test";
    expected.mtln->probes[0].probe_position = {2, 11, 14};

    expected.mtln->networks.resize(2);

    expected.mtln->networks[0].connections.resize(2);
    for (int c = 0; c < 2; ++c) {
        expected.mtln->networks[0].connections[c].nodes.resize(1);
        expected.mtln->networks[0].connections[c].nodes[0].side = TERMINAL_NODE_SIDE_INI;
        expected.mtln->networks[0].connections[c].nodes[0].belongs_to_cable = cable;
        expected.mtln->networks[0].connections[c].nodes[0].termination.termination_type =
            TERMINATION_OPEN;
    }
    expected.mtln->networks[0].connections[0].nodes[0].conductor_in_cable = 1;
    expected.mtln->networks[0].connections[1].nodes[0].conductor_in_cable = 2;

    expected.mtln->networks[1].connections.resize(2);
    for (int c = 0; c < 2; ++c) {
        expected.mtln->networks[1].connections[c].nodes.resize(1);
        expected.mtln->networks[1].connections[c].nodes[0].side = TERMINAL_NODE_SIDE_END;
        expected.mtln->networks[1].connections[c].nodes[0].belongs_to_cable = cable;
        expected.mtln->networks[1].connections[c].nodes[0].termination.termination_type =
            TERMINATION_OPEN;
    }
    expected.mtln->networks[1].connections[0].nodes[0].conductor_in_cable = 1;
    expected.mtln->networks[1].connections[1].nodes[0].conductor_in_cable = 2;

    return expected;
}

#endif // CompileWithMTLN

#endif // TEST_READ_UNSHIELDED_MULTIWIRES_MULTIPOLAR_EXPANSION_H
