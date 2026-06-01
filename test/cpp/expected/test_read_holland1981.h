#ifndef TEST_READ_HOLLAND1981_H
#define TEST_READ_HOLLAND1981_H

#include "test_smbjson_helpers.h"

inline Parseador_t expectedReadHolland1981() {
    Parseador_t expected = buildExpected();

    // Expected general info.
    expected.general->dt = 30e-12;
    expected.general->nmax = 1000;

    // Expected media matrix.
    expected.matriz->totalX = 21;
    expected.matriz->totalY = 21;
    expected.matriz->totalZ = 23;

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
    expected.despl->mx2 = 20;
    expected.despl->my1 = 0;
    expected.despl->my2 = 20;
    expected.despl->mz1 = 0;
    expected.despl->mz2 = 22;

    // Expected boundaries.
    for (int i = 0; i < 6; ++i) {
        expected.front->tipoFrontera[i] = F_PML;
        expected.front->propiedadesPML[i].numCapas = 6;
        expected.front->propiedadesPML[i].orden = 2.0;
        expected.front->propiedadesPML[i].refl = 0.001;
    }

    // Expected sources.
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

    // Expected probes - sonda
    expected.Sonda->length = 1;
    expected.Sonda->length_max = 1;
    expected.Sonda->len_cor_max = 1;
    expected.Sonda->collection.resize(1);
    expected.Sonda->collection[0].outputrequest = "mid_point";
    expected.Sonda->collection[0].type1 = NP_T1_PLAIN;
    expected.Sonda->collection[0].type2 = NP_T2_TIME;
    expected.Sonda->collection[0].filename = " ";
    expected.Sonda->collection[0].tstart = 0.0;
    expected.Sonda->collection[0].tstop = 0.0;
    expected.Sonda->collection[0].tstep = 0.0;
    expected.Sonda->collection[0].fstart = 0.0;
    expected.Sonda->collection[0].fstop = 0.0;
    expected.Sonda->collection[0].fstep = 0.0;
    expected.Sonda->collection[0].cordinates.resize(1);
    expected.Sonda->collection[0].len_cor = 1;
    expected.Sonda->collection[0].cordinates[0].tag = "mid_point";
    expected.Sonda->collection[0].cordinates[0].Xi = 8;
    expected.Sonda->collection[0].cordinates[0].Yi = 0;
    expected.Sonda->collection[0].cordinates[0].Zi = 0;
    expected.Sonda->collection[0].cordinates[0].Or = NP_COR_WIRECURRENT;

    // Expected thin wires
    expected.tWires->tw.resize(1);
    expected.tWires->tw[0].rad = 0.02;
    expected.tWires->tw[0].dispfile = "";
    expected.tWires->tw[0].dispfile_LeftEnd = "";
    expected.tWires->tw[0].dispfile_RightEnd = "";
    expected.tWires->tw[0].n_twc = 10;
    expected.tWires->tw[0].n_twc_max = 10;
    expected.tWires->tw[0].twc.resize(10);
    for (int i = 0; i < 10; ++i) {
        expected.tWires->tw[0].twc[i].srcfile = "None";
        expected.tWires->tw[0].twc[i].srctype = "None";
        expected.tWires->tw[0].twc[i].i = 11;
        expected.tWires->tw[0].twc[i].j = 11;
        expected.tWires->tw[0].twc[i].K = 7 + i;
        expected.tWires->tw[0].twc[i].d = 3; // DIR_Z from cells module
        expected.tWires->tw[0].twc[i].nd = 3 + i;
        expected.tWires->tw[0].twc[i].tag = "material1@layer2";
    }
    expected.tWires->tw[0].tl = MATERIAL_CONS;
    expected.tWires->tw[0].tr = MATERIAL_CONS;
    expected.tWires->tw[0].LeftEnd = 1;
    expected.tWires->tw[0].RightEnd = 2;
    expected.tWires->n_tw = 1;
    expected.tWires->n_tw_max = 1;

    return expected;
}

#endif // TEST_READ_HOLLAND1981_H
