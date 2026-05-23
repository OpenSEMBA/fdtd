#ifndef TEST_READ_SPHERE_H
#define TEST_READ_SPHERE_H

#include "test_smbjson_helpers.h"

inline Parseador_t expectedReadSphere() {
    Parseador_t expected = buildExpected();

    // Expected general info.
    expected.general->dt = 3.85167e-11;
    expected.general->nmax = 100;

    // Expected media matrix.
    expected.matriz->totalX = 81;
    expected.matriz->totalY = 81;
    expected.matriz->totalZ = 81;

    // Expected grid.
    expected.despl->nX = 1;
    expected.despl->nY = 1;
    expected.despl->nZ = 1;
    expected.despl->desX.resize(1);
    expected.despl->desY.resize(1);
    expected.despl->desZ.resize(1);
    expected.despl->desX[0] = 0.025;
    expected.despl->desY[0] = 0.025;
    expected.despl->desZ[0] = 0.025;
    expected.despl->mx1 = 0;
    expected.despl->mx2 = 80;
    expected.despl->my1 = 0;
    expected.despl->my2 = 80;
    expected.despl->mz1 = 0;
    expected.despl->mz2 = 80;

    // Expected boundaries.
    for (int i = 0; i < 6; ++i) {
        expected.front->tipoFrontera[i] = F_PML;
        expected.front->propiedadesPML[i].numCapas = 10;
        expected.front->propiedadesPML[i].orden = 2;
        expected.front->propiedadesPML[i].refl = 0.001;
    }

    // Expected material regions.
    expected.pecRegs->nSurfs = 1;
    expected.pecRegs->nSurfs_max = 1;
    expected.pecRegs->Surfs.resize(1);
    // -- specific surfs not included do NOT use comparison --

    // Expected sources.
    expected.plnSrc->collection.resize(1);
    expected.plnSrc->collection[0].nombre_fichero = "gauss.exc";
    expected.plnSrc->collection[0].atributo = "LOCKED";
    expected.plnSrc->collection[0].coor1[0] = 0;
    expected.plnSrc->collection[0].coor1[1] = 0;
    expected.plnSrc->collection[0].coor1[2] = 0;
    expected.plnSrc->collection[0].coor2[0] = 79;
    expected.plnSrc->collection[0].coor2[1] = 79;
    expected.plnSrc->collection[0].coor2[2] = 79;
    expected.plnSrc->collection[0].theta = 1.5707963268;
    expected.plnSrc->collection[0].phi = 0.0;
    expected.plnSrc->collection[0].alpha = 1.5707963268;
    expected.plnSrc->collection[0].beta = 4.7123889804;
    expected.plnSrc->collection[0].isRC = false;
    expected.plnSrc->collection[0].numModes = 1;
    expected.plnSrc->collection[0].INCERTMAX = 0.0;
    expected.plnSrc->nc = 1;
    expected.plnSrc->nC_max = 1;

    // Expected probes - old far-field sonda
    expected.oldSONDA->n_probes = 1;
    expected.oldSONDA->n_probes_max = 1;
    expected.oldSONDA->probes.resize(1);
    expected.oldSONDA->probes[0].n_FarField = 1;
    expected.oldSONDA->probes[0].n_FarField_max = 1;
    expected.oldSONDA->probes[0].FarField.resize(1);
    auto& ff = expected.oldSONDA->probes[0].FarField[0].probe;
    ff.grname = " ";
    ff.outputrequest = "FarField_log_";
    ff.tstart = 0.0;
    ff.tstop = 0.0;
    ff.tstep = 0.0;
    ff.fstart = 1e6;
    ff.fstop = 1e9;
    ff.fstep = (1e9 - 1e6) / 5;
    ff.FileNormalize = "gauss.exc";
    ff.i = {2, 77};
    ff.j = {2, 77};
    ff.K = {2, 77};
    ff.node.resize(0);
    ff.n_cord = 2;
    ff.n_cord_max = 2;
    ff.thetastart = 0.0;
    ff.thetastop = 180.0;
    ff.thetastep = 90.0;
    ff.phistart = 0.0;
    ff.phistop = 360.0;
    ff.phistep = 90.0;

    expected.VolPrb->length = 1;
    expected.VolPrb->length_max = 1;
    expected.VolPrb->len_cor_max = 2;
    expected.VolPrb->collection.resize(1);
    expected.VolPrb->collection[0].cordinates.resize(1);
    expected.VolPrb->collection[0].len_cor = 1;
    expected.VolPrb->collection[0].cordinates[0].Xi = 2;
    expected.VolPrb->collection[0].cordinates[0].Xe = 77;
    expected.VolPrb->collection[0].cordinates[0].Yi = 2;
    expected.VolPrb->collection[0].cordinates[0].Ye = 77;
    expected.VolPrb->collection[0].cordinates[0].Zi = 2;
    expected.VolPrb->collection[0].cordinates[0].Ze = 77;
    expected.VolPrb->collection[0].cordinates[0].Or = iExC;
    expected.VolPrb->collection[0].cordinates[0].Xtrancos = 1;
    expected.VolPrb->collection[0].cordinates[0].Ytrancos = 1;
    expected.VolPrb->collection[0].cordinates[0].Ztrancos = 1;
    expected.VolPrb->collection[0].cordinates[0].tag = "";
    expected.VolPrb->collection[0].tstart = 0.0;
    expected.VolPrb->collection[0].tstop = 0.0;
    expected.VolPrb->collection[0].tstep = 1e-9;
    expected.VolPrb->collection[0].fstart = 0.0;
    expected.VolPrb->collection[0].fstop = 0.0;
    expected.VolPrb->collection[0].fstep = 0.0;
    expected.VolPrb->collection[0].outputrequest = "electric_field_movie";
    expected.VolPrb->collection[0].filename = " ";
    expected.VolPrb->collection[0].type2 = NP_T2_TIME;

    return expected;
}

#endif // TEST_READ_SPHERE_H
