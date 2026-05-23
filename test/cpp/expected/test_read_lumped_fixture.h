#ifndef TEST_READ_LUMPED_FIXTURE_H
#define TEST_READ_LUMPED_FIXTURE_H

#include "test_smbjson_helpers.h"

inline Parseador_t expectedReadLumpedFixture() {
    Parseador_t expected = buildExpected();

    // Expected general info
    expected.general->dt = 7.7033e-12;
    expected.general->nmax = 389;

    // Expected media matrix
    expected.matriz->totalX = 21;
    expected.matriz->totalY = 9;
    expected.matriz->totalZ = 10;

    // Expected grid
    expected.despl->nX = 1;
    expected.despl->nY = 1;
    expected.despl->nZ = 1;
    expected.despl->desX.resize(1);
    expected.despl->desY.resize(1);
    expected.despl->desZ.resize(1);
    expected.despl->desX[0] = 0.005;
    expected.despl->desY[0] = 0.005;
    expected.despl->desZ[0] = 0.005;
    expected.despl->mx1 = 0;
    expected.despl->my1 = 0;
    expected.despl->mz1 = 0;
    expected.despl->mx2 = 20;
    expected.despl->my2 = 8;
    expected.despl->mz2 = 9;

    // Expected boundaries
    for (int i = 0; i < 6; ++i)
        expected.front->tipoFrontera[i] = F_MUR;

    // Expected material regions
    expected.pecRegs->nVols = 0;
    expected.pecRegs->nSurfs = 6;
    expected.pecRegs->nLins = 0;
    expected.pecRegs->nVols_max = 0;
    expected.pecRegs->nSurfs_max = 6;
    expected.pecRegs->nLins_max = 0;
    expected.pecRegs->Vols.resize(0);
    expected.pecRegs->Surfs.resize(6);
    expected.pecRegs->Lins.resize(0);

    // Left side - Surface 1
    expected.pecRegs->Surfs[0].Or = +iEx;
    expected.pecRegs->Surfs[0].Xi = 2;
    expected.pecRegs->Surfs[0].Xe = 2;
    expected.pecRegs->Surfs[0].Yi = 2;
    expected.pecRegs->Surfs[0].Ye = 5;
    expected.pecRegs->Surfs[0].Zi = 2;
    expected.pecRegs->Surfs[0].Ze = 6;
    expected.pecRegs->Surfs[0].tag = "pec@left_side";

    // Left side - Surface 2
    expected.pecRegs->Surfs[1].Or = +iEz;
    expected.pecRegs->Surfs[1].Xi = 2;
    expected.pecRegs->Surfs[1].Xe = 8;
    expected.pecRegs->Surfs[1].Yi = 2;
    expected.pecRegs->Surfs[1].Ye = 5;
    expected.pecRegs->Surfs[1].Zi = 2;
    expected.pecRegs->Surfs[1].Ze = 2;
    expected.pecRegs->Surfs[1].tag = "pec@left_side";

    // Left side - Surface 3
    expected.pecRegs->Surfs[2].Or = +iEz;
    expected.pecRegs->Surfs[2].Xi = 2;
    expected.pecRegs->Surfs[2].Xe = 8;
    expected.pecRegs->Surfs[2].Yi = 2;
    expected.pecRegs->Surfs[2].Ye = 5;
    expected.pecRegs->Surfs[2].Zi = 7;
    expected.pecRegs->Surfs[2].Ze = 7;
    expected.pecRegs->Surfs[2].tag = "pec@left_side";

    // Right side - Surface 1
    expected.pecRegs->Surfs[3].Or = +iEz;
    expected.pecRegs->Surfs[3].Xi = 11;
    expected.pecRegs->Surfs[3].Xe = 17;
    expected.pecRegs->Surfs[3].Yi = 2;
    expected.pecRegs->Surfs[3].Ye = 5;
    expected.pecRegs->Surfs[3].Zi = 2;
    expected.pecRegs->Surfs[3].Ze = 2;
    expected.pecRegs->Surfs[3].tag = "pec@right_side";

    // Right side - Surface 2
    expected.pecRegs->Surfs[4].Or = +iEz;
    expected.pecRegs->Surfs[4].Xi = 11;
    expected.pecRegs->Surfs[4].Xe = 17;
    expected.pecRegs->Surfs[4].Yi = 2;
    expected.pecRegs->Surfs[4].Ye = 5;
    expected.pecRegs->Surfs[4].Zi = 7;
    expected.pecRegs->Surfs[4].Ze = 7;
    expected.pecRegs->Surfs[4].tag = "pec@right_side";

    // Right side - Surface 3
    expected.pecRegs->Surfs[5].Or = +iEx;
    expected.pecRegs->Surfs[5].Xi = 18;
    expected.pecRegs->Surfs[5].Xe = 18;
    expected.pecRegs->Surfs[5].Yi = 2;
    expected.pecRegs->Surfs[5].Ye = 5;
    expected.pecRegs->Surfs[5].Zi = 2;
    expected.pecRegs->Surfs[5].Ze = 6;
    expected.pecRegs->Surfs[5].tag = "pec@right_side";

    // Expected dielectric regions (including lumped resistor)
    expected.DielRegs->nVols = 0;
    expected.DielRegs->nSurfs = 0;
    expected.DielRegs->nLins = 1;
    expected.DielRegs->nVols_max = 0;
    expected.DielRegs->nSurfs_max = 0;
    expected.DielRegs->nLins_max = 1;
    expected.DielRegs->Vols.resize(0);
    expected.DielRegs->Surfs.resize(0);
    expected.DielRegs->Lins.resize(1);

    // Lumped resistor line
    expected.DielRegs->Lins[0].c1P.resize(0);
    expected.DielRegs->Lins[0].c2P.resize(1);
    expected.DielRegs->Lins[0].n_C1P = 0;
    expected.DielRegs->Lins[0].n_C2P = 1;
    expected.DielRegs->Lins[0].c2P[0].Or = iEx;
    expected.DielRegs->Lins[0].c2P[0].Xi = 9;
    expected.DielRegs->Lins[0].c2P[0].Xe = 10;
    expected.DielRegs->Lins[0].c2P[0].Yi = 4;
    expected.DielRegs->Lins[0].c2P[0].Ye = 4;
    expected.DielRegs->Lins[0].c2P[0].Zi = 7;
    expected.DielRegs->Lins[0].c2P[0].Ze = 7;
    expected.DielRegs->Lins[0].c2P[0].tag = "100ohm_resistor@lumped_line";
    expected.DielRegs->Lins[0].sigma = 0.0;
    expected.DielRegs->Lins[0].eps = EPSILON_VACUUM;
    expected.DielRegs->Lins[0].mu = MU_VACUUM;
    expected.DielRegs->Lins[0].sigmam = 0.0;
    expected.DielRegs->Lins[0].R = 100.0;
    expected.DielRegs->Lins[0].Rtime_on = 0.0;
    expected.DielRegs->Lins[0].Rtime_off = 1.0;
    expected.DielRegs->Lins[0].resistor = true;
    expected.DielRegs->Lins[0].orient = 1;
    expected.DielRegs->Lins[0].DiodOri = 1;

    // Expected sources
    expected.nodSrc->n_nodSrc = 1;
    expected.nodSrc->n_nodSrc_max = 1;
    expected.nodSrc->n_C1P_max = 1;
    expected.nodSrc->n_C2P_max = 1;
    expected.nodSrc->NodalSource.resize(1);
    expected.nodSrc->NodalSource[0].nombre = "predefinedExcitation.1.exc";
    expected.nodSrc->NodalSource[0].isElec = true;
    expected.nodSrc->NodalSource[0].isHard = false;
    expected.nodSrc->NodalSource[0].isInitialValue = false;
    expected.nodSrc->NodalSource[0].c2P.resize(1);
    expected.nodSrc->NodalSource[0].n_C2P = 1;
    expected.nodSrc->NodalSource[0].c2P[0].Or = iEx;
    expected.nodSrc->NodalSource[0].c2P[0].Xi = 9;
    expected.nodSrc->NodalSource[0].c2P[0].Xe = 10;
    expected.nodSrc->NodalSource[0].c2P[0].Yi = 4;
    expected.nodSrc->NodalSource[0].c2P[0].Ye = 4;
    expected.nodSrc->NodalSource[0].c2P[0].Zi = 2;
    expected.nodSrc->NodalSource[0].c2P[0].Ze = 2;
    expected.nodSrc->NodalSource[0].c2P[0].tag = "";
    expected.nodSrc->NodalSource[0].c2P[0].xc = 1.0;
    expected.nodSrc->NodalSource[0].c2P[0].yc = 0.0;
    expected.nodSrc->NodalSource[0].c2P[0].zc = 0.0;

    // Expected probes - electric field point probe
    expected.Sonda->length = 1;
    expected.Sonda->length_max = 1;
    expected.Sonda->len_cor_max = 3;
    expected.Sonda->collection.resize(1);
    expected.Sonda->collection[0].outputrequest = "e_probe";
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
        expected.Sonda->collection[0].cordinates[c].Xi = 10;
        expected.Sonda->collection[0].cordinates[c].Yi = 3;
        expected.Sonda->collection[0].cordinates[c].Zi = 7;
        expected.Sonda->collection[0].cordinates[c].tag = "e_probe";
    }
    expected.Sonda->collection[0].cordinates[0].Or = NP_COR_EX;
    expected.Sonda->collection[0].cordinates[1].Or = NP_COR_EY;
    expected.Sonda->collection[0].cordinates[2].Or = NP_COR_EZ;

    // Bulk current probe
    expected.BloquePrb->n_bp = 1;
    expected.BloquePrb->n_bp_max = 1;
    expected.BloquePrb->bp.resize(1);
    expected.BloquePrb->bp[0].outputrequest = "Bulk probe";
    expected.BloquePrb->bp[0].FileNormalize = " ";
    expected.BloquePrb->bp[0].type2 = NP_T2_TIME;
    expected.BloquePrb->bp[0].tstart = 0.0;
    expected.BloquePrb->bp[0].tstop = 0.0;
    expected.BloquePrb->bp[0].tstep = 0.0;
    expected.BloquePrb->bp[0].fstart = 0.0;
    expected.BloquePrb->bp[0].fstop = 0.0;
    expected.BloquePrb->bp[0].fstep = 0.0;
    expected.BloquePrb->bp[0].i1 = 6;
    expected.BloquePrb->bp[0].i2 = 6;
    expected.BloquePrb->bp[0].j1 = 1;
    expected.BloquePrb->bp[0].j2 = 6;
    expected.BloquePrb->bp[0].k1 = 6;
    expected.BloquePrb->bp[0].k2 = 7;
    expected.BloquePrb->bp[0].skip = 1;
    expected.BloquePrb->bp[0].nml = iEx;
    expected.BloquePrb->bp[0].t = BcELECT;
    expected.BloquePrb->bp[0].tag = "Bulk probe";

    return expected;
}

#endif // TEST_READ_LUMPED_FIXTURE_H
