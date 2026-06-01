#ifndef TEST_READ_DIELECTRICSLAB_H
#define TEST_READ_DIELECTRICSLAB_H

#include "test_smbjson_helpers.h"

inline Parseador_t expectedReadDielectricSlab() {
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
    expected.front->tipoFrontera[F_XL - 1] = F_PEC;
    expected.front->tipoFrontera[F_XU - 1] = F_PEC;
    expected.front->tipoFrontera[F_YL - 1] = F_PMC;
    expected.front->tipoFrontera[F_YU - 1] = F_PMC;
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

    // Dielectric slab region
    expected.DielRegs->nVols = 1;
    expected.DielRegs->nSurfs = 0;
    expected.DielRegs->nLins = 0;
    expected.DielRegs->nVols_max = 1;
    expected.DielRegs->nSurfs_max = 0;
    expected.DielRegs->Vols.resize(1);
    expected.DielRegs->Surfs.resize(0);
    expected.DielRegs->Lins.resize(0);
    expected.DielRegs->Vols[0].c1P.resize(0);
    expected.DielRegs->Vols[0].c2P.resize(1);
    expected.DielRegs->Vols[0].n_C1P = 0;
    expected.DielRegs->Vols[0].n_C2P = 1;
    expected.DielRegs->Vols[0].sigma = 0.0;
    expected.DielRegs->Vols[0].eps = 2.1 * EPSILON_VACUUM;
    expected.DielRegs->Vols[0].mu = MU_VACUUM;
    expected.DielRegs->Vols[0].sigmam = 0.0;
    expected.DielRegs->Vols[0].c2P[0].Or = 0;
    expected.DielRegs->Vols[0].c2P[0].Xi = 0;
    expected.DielRegs->Vols[0].c2P[0].Xe = 3;
    expected.DielRegs->Vols[0].c2P[0].Yi = 0;
    expected.DielRegs->Vols[0].c2P[0].Ye = 3;
    expected.DielRegs->Vols[0].c2P[0].Zi = 20;
    expected.DielRegs->Vols[0].c2P[0].Ze = 29;
    expected.DielRegs->Vols[0].c2P[0].tag = "teflon@slab";

    // Expected probes - sonda
    expected.Sonda->len_cor_max = 3;
    expected.Sonda->length = 3;
    expected.Sonda->length_max = 3;
    expected.Sonda->collection.resize(3);
    for (int i = 0; i < 3; ++i) {
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
    // Point probe in dielectric slab
    expected.Sonda->collection[1].outputrequest = "inner";
    for (int c = 0; c < 3; ++c) {
        expected.Sonda->collection[1].cordinates[c].tag = "inner";
        expected.Sonda->collection[1].cordinates[c].Xi = 2;
        expected.Sonda->collection[1].cordinates[c].Yi = 2;
        expected.Sonda->collection[1].cordinates[c].Zi = 25;
    }
    // Point probe at back
    expected.Sonda->collection[2].outputrequest = "back";
    for (int c = 0; c < 3; ++c) {
        expected.Sonda->collection[2].cordinates[c].tag = "back";
        expected.Sonda->collection[2].cordinates[c].Xi = 2;
        expected.Sonda->collection[2].cordinates[c].Yi = 2;
        expected.Sonda->collection[2].cordinates[c].Zi = 40;
    }

    return expected;
}

#endif // TEST_READ_DIELECTRICSLAB_H
