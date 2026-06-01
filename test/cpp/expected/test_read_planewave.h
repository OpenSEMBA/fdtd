#ifndef TEST_READ_PLANEWAVE_H
#define TEST_READ_PLANEWAVE_H

#include "test_smbjson_helpers.h"

inline Parseador_t expectedReadPlanewave() {
    Parseador_t expected = buildExpected();

    // Expected general info.
    expected.general->dt = 10e-12;
    expected.general->nmax = 2000;

    // Expected media matrix.
    expected.matriz->totalX = 11;
    expected.matriz->totalY = 11;
    expected.matriz->totalZ = 11;

    // Expected grid.
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
    expected.despl->mx2 = 10;
    expected.despl->my1 = 0;
    expected.despl->my2 = 10;
    expected.despl->mz1 = 0;
    expected.despl->mz2 = 10;

    // Expected boundaries.
    for (int i = 0; i < 6; ++i)
        expected.front->tipoFrontera[i] = F_MUR;

    // Expected sources.
    expected.plnSrc->collection.resize(1);
    expected.plnSrc->collection[0].nombre_fichero = "gauss.exc";
    expected.plnSrc->collection[0].atributo = "LOCKED";
    expected.plnSrc->collection[0].coor1[0] = 1;
    expected.plnSrc->collection[0].coor1[1] = 1;
    expected.plnSrc->collection[0].coor1[2] = 1;
    expected.plnSrc->collection[0].coor2[0] = 8;
    expected.plnSrc->collection[0].coor2[1] = 8;
    expected.plnSrc->collection[0].coor2[2] = 8;
    expected.plnSrc->collection[0].theta = 0.0;
    expected.plnSrc->collection[0].phi = 0.0;
    expected.plnSrc->collection[0].alpha = 1.5708;
    expected.plnSrc->collection[0].beta = 0.0;
    expected.plnSrc->collection[0].isRC = false;
    expected.plnSrc->collection[0].numModes = 1;
    expected.plnSrc->collection[0].INCERTMAX = 0.0;
    expected.plnSrc->nc = 1;
    expected.plnSrc->nC_max = 1;

    // Expected probes - sonda
    expected.Sonda->len_cor_max = 3;
    expected.Sonda->length = 2;
    expected.Sonda->length_max = 2;
    expected.Sonda->collection.resize(2);

    expected.Sonda->collection[0].outputrequest = "electric_field_point_probe";
    expected.Sonda->collection[0].type1 = NP_T1_PLAIN;
    expected.Sonda->collection[0].type2 = NP_T2_TIME;
    expected.Sonda->collection[0].filename = " ";
    expected.Sonda->collection[0].tstart = 0.0;
    expected.Sonda->collection[0].tstop = 0.0;
    expected.Sonda->collection[0].tstep = 0.0;
    expected.Sonda->collection[0].fstart = 0.0;
    expected.Sonda->collection[0].fstop = 0.0;
    expected.Sonda->collection[0].fstep = 0.0;
    expected.Sonda->collection[0].cordinates.resize(3);
    expected.Sonda->collection[0].len_cor = 3;
    for (int c = 0; c < 3; ++c) {
        expected.Sonda->collection[0].cordinates[c].tag = "electric_field_point_probe";
        expected.Sonda->collection[0].cordinates[c].Xi = 4;
        expected.Sonda->collection[0].cordinates[c].Yi = 4;
        expected.Sonda->collection[0].cordinates[c].Zi = 4;
    }
    expected.Sonda->collection[0].cordinates[0].Or = NP_COR_EX;
    expected.Sonda->collection[0].cordinates[1].Or = NP_COR_EY;
    expected.Sonda->collection[0].cordinates[2].Or = NP_COR_EZ;

    expected.Sonda->collection[1].outputrequest = "magnetic_field_point_probe";
    expected.Sonda->collection[1].type1 = NP_T1_PLAIN;
    expected.Sonda->collection[1].type2 = NP_T2_TIME;
    expected.Sonda->collection[1].filename = " ";
    expected.Sonda->collection[1].tstart = 0.0;
    expected.Sonda->collection[1].tstop = 0.0;
    expected.Sonda->collection[1].tstep = 0.0;
    expected.Sonda->collection[1].fstart = 0.0;
    expected.Sonda->collection[1].fstop = 0.0;
    expected.Sonda->collection[1].fstep = 0.0;
    expected.Sonda->collection[1].cordinates.resize(1);
    expected.Sonda->collection[1].len_cor = 1;
    expected.Sonda->collection[1].cordinates[0].tag = "magnetic_field_point_probe";
    expected.Sonda->collection[1].cordinates[0].Xi = 6;
    expected.Sonda->collection[1].cordinates[0].Yi = 6;
    expected.Sonda->collection[1].cordinates[0].Zi = 6;
    expected.Sonda->collection[1].cordinates[0].Or = NP_COR_HX;

    return expected;
}

#endif // TEST_READ_PLANEWAVE_H
