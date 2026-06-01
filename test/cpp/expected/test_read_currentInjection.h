#ifndef TEST_READ_CURRENTINJECTION_H
#define TEST_READ_CURRENTINJECTION_H

#include "test_smbjson_helpers.h"

inline Parseador_t expectedReadCurrentInjection() {
    Parseador_t expected = buildExpected();

    // Expected general info.
    expected.general->dt = 1e-12;
    expected.general->nmax = 1000;

    // Expected media matrix.
    expected.matriz->totalX = 21;
    expected.matriz->totalY = 21;
    expected.matriz->totalZ = 21;

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
    expected.despl->my1 = 0;
    expected.despl->mz1 = 0;
    expected.despl->mx2 = 20;
    expected.despl->my2 = 20;
    expected.despl->mz2 = 20;

    // Expected boundaries.
    for (int i = 0; i < 6; ++i)
        expected.front->tipoFrontera[i] = F_MUR;

    // Expected material regions.
    expected.pecRegs->nVols = 0;
    expected.pecRegs->nSurfs = 1;
    expected.pecRegs->nLins = 1;
    expected.pecRegs->nVols_max = 0;
    expected.pecRegs->nSurfs_max = 1;
    expected.pecRegs->nLins_max = 1;
    expected.pecRegs->Vols.resize(0);
    expected.pecRegs->Surfs.resize(1);
    expected.pecRegs->Lins.resize(1);

    // Body
    expected.pecRegs->Surfs[0].Or = +iEz;
    expected.pecRegs->Surfs[0].Xi = 5;
    expected.pecRegs->Surfs[0].Xe = 14;
    expected.pecRegs->Surfs[0].Yi = 5;
    expected.pecRegs->Surfs[0].Ye = 14;
    expected.pecRegs->Surfs[0].Zi = 10;
    expected.pecRegs->Surfs[0].Ze = 10;
    expected.pecRegs->Surfs[0].tag = "aluminum@body";

    // Exit line
    expected.pecRegs->Lins[0].Or = +iEy;
    expected.pecRegs->Lins[0].Xi = 10;
    expected.pecRegs->Lins[0].Xe = 10;
    expected.pecRegs->Lins[0].Yi = 15;
    expected.pecRegs->Lins[0].Ye = 19;
    expected.pecRegs->Lins[0].Zi = 10;
    expected.pecRegs->Lins[0].Ze = 10;
    expected.pecRegs->Lins[0].tag = "aluminum@exit";

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
    expected.nodSrc->NodalSource[0].c2P[0].Or = iEy;
    expected.nodSrc->NodalSource[0].c2P[0].Xi = 10;
    expected.nodSrc->NodalSource[0].c2P[0].Xe = 10;
    expected.nodSrc->NodalSource[0].c2P[0].Yi = 0;
    expected.nodSrc->NodalSource[0].c2P[0].Ye = 4;
    expected.nodSrc->NodalSource[0].c2P[0].Zi = 10;
    expected.nodSrc->NodalSource[0].c2P[0].Ze = 10;
    expected.nodSrc->NodalSource[0].c2P[0].tag = "";
    expected.nodSrc->NodalSource[0].c2P[0].xc = 0.0;
    expected.nodSrc->NodalSource[0].c2P[0].yc = 1.0;
    expected.nodSrc->NodalSource[0].c2P[0].zc = 0.0;

    // Expected probes - bloqueprobes
    expected.BloquePrb->n_bp = 2;
    expected.BloquePrb->n_bp_max = 2;
    expected.BloquePrb->bp.resize(2);

    expected.BloquePrb->bp[0].outputrequest = "bulk_current_at_entry";
    expected.BloquePrb->bp[0].FileNormalize = " ";
    expected.BloquePrb->bp[0].type2 = NP_T2_TIME;
    expected.BloquePrb->bp[0].tstart = 0.0;
    expected.BloquePrb->bp[0].tstop = 0.0;
    expected.BloquePrb->bp[0].tstep = 0.0;
    expected.BloquePrb->bp[0].fstart = 0.0;
    expected.BloquePrb->bp[0].fstop = 0.0;
    expected.BloquePrb->bp[0].fstep = 0.0;
    expected.BloquePrb->bp[0].i1 = 9;
    expected.BloquePrb->bp[0].i2 = 10;
    expected.BloquePrb->bp[0].j1 = 2;
    expected.BloquePrb->bp[0].j2 = 2;
    expected.BloquePrb->bp[0].k1 = 9;
    expected.BloquePrb->bp[0].k2 = 10;
    expected.BloquePrb->bp[0].skip = 1;
    expected.BloquePrb->bp[0].nml = iEy;
    expected.BloquePrb->bp[0].t = BcELECT;
    expected.BloquePrb->bp[0].tag = "bulk_current_at_entry";

    expected.BloquePrb->bp[1].outputrequest = "bulk_current_at_exit";
    expected.BloquePrb->bp[1].FileNormalize = " ";
    expected.BloquePrb->bp[1].type2 = NP_T2_TIME;
    expected.BloquePrb->bp[1].tstart = 0.0;
    expected.BloquePrb->bp[1].tstop = 0.0;
    expected.BloquePrb->bp[1].tstep = 0.0;
    expected.BloquePrb->bp[1].fstart = 0.0;
    expected.BloquePrb->bp[1].fstop = 0.0;
    expected.BloquePrb->bp[1].fstep = 0.0;
    expected.BloquePrb->bp[1].i1 = 9;
    expected.BloquePrb->bp[1].i2 = 10;
    expected.BloquePrb->bp[1].j1 = 17;
    expected.BloquePrb->bp[1].j2 = 17;
    expected.BloquePrb->bp[1].k1 = 9;
    expected.BloquePrb->bp[1].k2 = 10;
    expected.BloquePrb->bp[1].skip = 1;
    expected.BloquePrb->bp[1].nml = iEy;
    expected.BloquePrb->bp[1].t = BcELECT;
    expected.BloquePrb->bp[1].tag = "bulk_current_at_exit";

    return expected;
}

#endif // TEST_READ_CURRENTINJECTION_H
