#ifndef TEST_READ_THINSLOT_H
#define TEST_READ_THINSLOT_H

#include "test_smbjson_helpers.h"

inline Parseador_t expectedReadThinSlot() {
    Parseador_t expected = buildExpected();

    // Expected general info.
    expected.general->dt = 10e-12;
    expected.general->nmax = 2000;

    // Expected media matrix.
    expected.matriz->totalX = 5;
    expected.matriz->totalY = 5;
    expected.matriz->totalZ = 51;

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
    expected.despl->mx2 = 4;
    expected.despl->my1 = 0;
    expected.despl->my2 = 4;
    expected.despl->mz1 = 0;
    expected.despl->mz2 = 50;

    // Expected boundaries.
    expected.front->tipoFrontera[F_XL - 1] = F_PER;
    expected.front->tipoFrontera[F_XU - 1] = F_PER;
    expected.front->tipoFrontera[F_YL - 1] = F_PER;
    expected.front->tipoFrontera[F_YU - 1] = F_PER;
    expected.front->tipoFrontera[F_ZL - 1] = F_MUR;
    expected.front->tipoFrontera[F_ZU - 1] = F_MUR;

    // Expected sources.
    expected.plnSrc->collection.resize(1);
    expected.plnSrc->collection[0].nombre_fichero = "gauss.exc";
    expected.plnSrc->collection[0].atributo = "LOCKED";
    expected.plnSrc->collection[0].coor1[0] = 0;
    expected.plnSrc->collection[0].coor1[1] = 0;
    expected.plnSrc->collection[0].coor1[2] = 2;
    expected.plnSrc->collection[0].coor2[0] = 3;
    expected.plnSrc->collection[0].coor2[1] = 3;
    expected.plnSrc->collection[0].coor2[2] = 47;
    expected.plnSrc->collection[0].theta = 0.0;
    expected.plnSrc->collection[0].phi = 0.0;
    expected.plnSrc->collection[0].alpha = 1.5708;
    expected.plnSrc->collection[0].beta = 0.0;
    expected.plnSrc->collection[0].isRC = false;
    expected.plnSrc->collection[0].numModes = 1;
    expected.plnSrc->collection[0].INCERTMAX = 0.0;
    expected.plnSrc->nc = 1;
    expected.plnSrc->nC_max = 1;

    // Materials - PEC square
    expected.pecRegs->nVols = 0;
    expected.pecRegs->nSurfs = 1;
    expected.pecRegs->nLins = 0;
    expected.pecRegs->nVols_max = 0;
    expected.pecRegs->nSurfs_max = 1;
    expected.pecRegs->nLins_max = 0;
    expected.pecRegs->Vols.resize(0);
    expected.pecRegs->Surfs.resize(1);
    expected.pecRegs->Lins.resize(0);
    expected.pecRegs->Surfs[0].Or = +iEz;
    expected.pecRegs->Surfs[0].Xi = 0;
    expected.pecRegs->Surfs[0].Xe = 3;
    expected.pecRegs->Surfs[0].Yi = 0;
    expected.pecRegs->Surfs[0].Ye = 3;
    expected.pecRegs->Surfs[0].Zi = 25;
    expected.pecRegs->Surfs[0].Ze = 25;
    expected.pecRegs->Surfs[0].tag = "copper@square";

    // Thin slot
    expected.tSlots->n_tg = 1;
    expected.tSlots->tg.resize(1);
    expected.tSlots->tg[0].width = 3e-3;
    expected.tSlots->tg[0].n_tgc = 2;
    expected.tSlots->tg[0].tgc.resize(2);
    expected.tSlots->tg[0].tgc[0].i = 1;
    expected.tSlots->tg[0].tgc[0].j = 2;
    expected.tSlots->tg[0].tgc[0].K = 25;
    expected.tSlots->tg[0].tgc[0].node = 0;
    expected.tSlots->tg[0].tgc[0].dir = iEx;
    expected.tSlots->tg[0].tgc[0].Or = -1;
    expected.tSlots->tg[0].tgc[0].tag = "3mm-gap@slot";
    expected.tSlots->tg[0].tgc[1] = expected.tSlots->tg[0].tgc[0];
    expected.tSlots->tg[0].tgc[1].i = 2;

    // Expected probes - sonda
    expected.Sonda->len_cor_max = 3;
    expected.Sonda->length = 2;
    expected.Sonda->length_max = 2;
    expected.Sonda->collection.resize(2);
    for (int i = 0; i < 2; ++i) {
        expected.Sonda->collection[i].type1 = NP_T1_PLAIN;
        expected.Sonda->collection[i].type2 = NP_T2_TIME;
        expected.Sonda->collection[i].filename = " ";
        expected.Sonda->collection[i].tstart = 0.0;
        expected.Sonda->collection[i].tstop = 0.0;
        expected.Sonda->collection[i].tstep = 0.0;
        expected.Sonda->collection[i].fstart = 0.0;
        expected.Sonda->collection[i].fstop = 0.0;
        expected.Sonda->collection[i].fstep = 0.0;
        expected.Sonda->collection[i].cordinates.resize(3);
        expected.Sonda->collection[i].cordinates[0].Or = NP_COR_EX;
        expected.Sonda->collection[i].cordinates[1].Or = NP_COR_EY;
        expected.Sonda->collection[i].cordinates[2].Or = NP_COR_EZ;
        expected.Sonda->collection[i].len_cor = 3;
    }
    // Point probe at front
    expected.Sonda->collection[0].outputrequest = "front";
    for (int c = 0; c < 3; ++c) {
        expected.Sonda->collection[0].cordinates[c].tag = "front";
        expected.Sonda->collection[0].cordinates[c].Xi = 2;
        expected.Sonda->collection[0].cordinates[c].Yi = 2;
        expected.Sonda->collection[0].cordinates[c].Zi = 10;
    }
    // Point probe at back
    expected.Sonda->collection[1].outputrequest = "back";
    for (int c = 0; c < 3; ++c) {
        expected.Sonda->collection[1].cordinates[c].tag = "back";
        expected.Sonda->collection[1].cordinates[c].Xi = 2;
        expected.Sonda->collection[1].cordinates[c].Yi = 2;
        expected.Sonda->collection[1].cordinates[c].Zi = 40;
    }

    return expected;
}

#endif // TEST_READ_THINSLOT_H
