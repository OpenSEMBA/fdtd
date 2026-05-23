#ifndef TEST_READ_AIRPLANE_H
#define TEST_READ_AIRPLANE_H

#include "test_smbjson_helpers.h"

inline Parseador_t expectedReadAirplane() {
    Parseador_t expected = buildExpected();

    // Expected general info.
    expected.general->dt = 1.2971876e-09;
    expected.general->nmax = 23;

    // Expected media matrix.
    expected.matriz->totalX = 51;
    expected.matriz->totalY = 51;
    expected.matriz->totalZ = 51;

    // Expected grid.
    expected.despl->nX = 1;
    expected.despl->nY = 1;
    expected.despl->nZ = 1;
    expected.despl->desX.resize(1);
    expected.despl->desY.resize(1);
    expected.despl->desZ.resize(1);
    expected.despl->desX[0] = 0.32419496007084553;
    expected.despl->desY[0] = 0.12839303226115248;
    expected.despl->desZ[0] = 0.3621442456908099;
    expected.despl->originx = 1.0;
    expected.despl->originy = 2.0;
    expected.despl->originz = 3.0;
    expected.despl->mx1 = 0;
    expected.despl->mx2 = 50;
    expected.despl->my1 = 0;
    expected.despl->my2 = 50;
    expected.despl->mz1 = 0;
    expected.despl->mz2 = 50;

    // Expected boundaries.
    for (int i = 0; i < 6; ++i) {
        expected.front->tipoFrontera[i] = F_PML;
        expected.front->propiedadesPML[i].numCapas = 10;
        expected.front->propiedadesPML[i].orden = 2.0;
        expected.front->propiedadesPML[i].refl = 0.001;
    }

    // Expected sources.
    expected.nodSrc->n_nodSrc = 1;
    expected.nodSrc->n_nodSrc_max = 1;
    expected.nodSrc->n_C1P_max = 1;
    expected.nodSrc->n_C2P_max = 1;
    expected.nodSrc->NodalSource.resize(1);
    expected.nodSrc->NodalSource[0].nombre = "gauss.exc";
    expected.nodSrc->NodalSource[0].isElec = true;
    expected.nodSrc->NodalSource[0].isHard = false;
    expected.nodSrc->NodalSource[0].isInitialValue = false;
    expected.nodSrc->NodalSource[0].c2P.resize(1);
    expected.nodSrc->NodalSource[0].n_C2P = 1;
    expected.nodSrc->NodalSource[0].c2P[0].Or = iEz;
    expected.nodSrc->NodalSource[0].c2P[0].Xi = 5;
    expected.nodSrc->NodalSource[0].c2P[0].Xe = 5;
    expected.nodSrc->NodalSource[0].c2P[0].Yi = 30;
    expected.nodSrc->NodalSource[0].c2P[0].Ye = 30;
    expected.nodSrc->NodalSource[0].c2P[0].Zi = 39;
    expected.nodSrc->NodalSource[0].c2P[0].Ze = 46;
    expected.nodSrc->NodalSource[0].c2P[0].tag = "";
    expected.nodSrc->NodalSource[0].c2P[0].xc = 0.0;
    expected.nodSrc->NodalSource[0].c2P[0].yc = 0.0;
    expected.nodSrc->NodalSource[0].c2P[0].zc = 1.0;

    return expected;
}

#endif // TEST_READ_AIRPLANE_H
