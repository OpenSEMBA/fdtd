#ifndef TEST_READ_SGBC_H
#define TEST_READ_SGBC_H

#include "test_smbjson_helpers.h"

inline Parseador_t expectedReadSgbc() {
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

    // Expected materials - PECs
    expected.pecRegs->nSurfs = 1;
    expected.pecRegs->nLins = 0;
    expected.pecRegs->nVols_max = 0;
    expected.pecRegs->nSurfs_max = 1;
    expected.pecRegs->nLins_max = 0;
    expected.pecRegs->Vols.resize(0);
    expected.pecRegs->Surfs.resize(1);

    // 2x2 PEC square
    expected.pecRegs->Surfs[0].Or = +iEz;
    expected.pecRegs->Surfs[0].Xi = 3;
    expected.pecRegs->Surfs[0].Xe = 4;
    expected.pecRegs->Surfs[0].Yi = 3;
    expected.pecRegs->Surfs[0].Ye = 4;
    expected.pecRegs->Surfs[0].Zi = 3;
    expected.pecRegs->Surfs[0].Ze = 3;
    expected.pecRegs->Surfs[0].tag = "material1@layer1";

    // Composites
    expected.LossyThinSurfs->cs.resize(2);
    expected.LossyThinSurfs->length = 2;
    expected.LossyThinSurfs->length_max = 2;
    expected.LossyThinSurfs->nC_max = 1;

    // 2-layer composite
    expected.LossyThinSurfs->cs[0].c.resize(1);
    expected.LossyThinSurfs->cs[0].nc = 1;
    expected.LossyThinSurfs->cs[0].files = "2-layers-composite";
    expected.LossyThinSurfs->cs[0].c[0].tag = "2-layers-composite@layer2";
    expected.LossyThinSurfs->cs[0].c[0].Or = +iEy;
    expected.LossyThinSurfs->cs[0].c[0].Xi = 3;
    expected.LossyThinSurfs->cs[0].c[0].Xe = 4;
    expected.LossyThinSurfs->cs[0].c[0].Yi = 3;
    expected.LossyThinSurfs->cs[0].c[0].Ye = 3;
    expected.LossyThinSurfs->cs[0].c[0].Zi = 3;
    expected.LossyThinSurfs->cs[0].c[0].Ze = 4;
    expected.LossyThinSurfs->cs[0].numcapas = 2;
    expected.LossyThinSurfs->cs[0].thk = {1e-3, 5e-3};
    expected.LossyThinSurfs->cs[0].sigma = {2e-4, 0.0};
    expected.LossyThinSurfs->cs[0].eps = {1.3 * EPSILON_VACUUM, 1.3 * EPSILON_VACUUM};
    expected.LossyThinSurfs->cs[0].mu = {MU_VACUUM, MU_VACUUM};
    expected.LossyThinSurfs->cs[0].sigmam = {0.0, 0.0};
    expected.LossyThinSurfs->cs[0].thk_devia = {0.0, 0.0};
    expected.LossyThinSurfs->cs[0].sigma_devia = {0.0, 0.0};
    expected.LossyThinSurfs->cs[0].eps_devia = {0.0, 0.0};
    expected.LossyThinSurfs->cs[0].mu_devia = {0.0, 0.0};
    expected.LossyThinSurfs->cs[0].sigmam_devia = {0.0, 0.0};

    // 3-layer composite
    expected.LossyThinSurfs->cs[1].c.resize(1);
    expected.LossyThinSurfs->cs[1].nc = 1;
    expected.LossyThinSurfs->cs[1].files = "3-layers-composite";
    expected.LossyThinSurfs->cs[1].c[0].tag = "3-layers-composite@layer3";
    expected.LossyThinSurfs->cs[1].c[0].Or = +iEx;
    expected.LossyThinSurfs->cs[1].c[0].Xi = 3;
    expected.LossyThinSurfs->cs[1].c[0].Xe = 3;
    expected.LossyThinSurfs->cs[1].c[0].Yi = 3;
    expected.LossyThinSurfs->cs[1].c[0].Ye = 4;
    expected.LossyThinSurfs->cs[1].c[0].Zi = 3;
    expected.LossyThinSurfs->cs[1].c[0].Ze = 4;
    expected.LossyThinSurfs->cs[1].numcapas = 3;
    expected.LossyThinSurfs->cs[1].thk = {1e-3, 5e-3, 1e-3};
    expected.LossyThinSurfs->cs[1].sigma = {2e-4, 0.0, 0.0};
    expected.LossyThinSurfs->cs[1].eps = {EPSILON_VACUUM, EPSILON_VACUUM, EPSILON_VACUUM};
    expected.LossyThinSurfs->cs[1].mu = {MU_VACUUM, 1.3 * MU_VACUUM, MU_VACUUM};
    expected.LossyThinSurfs->cs[1].sigmam = {0.0, 0.0, 1e-4};
    expected.LossyThinSurfs->cs[1].thk_devia = {0.0, 0.0, 0.0};
    expected.LossyThinSurfs->cs[1].sigma_devia = {0.0, 0.0, 0.0};
    expected.LossyThinSurfs->cs[1].eps_devia = {0.0, 0.0, 0.0};
    expected.LossyThinSurfs->cs[1].mu_devia = {0.0, 0.0, 0.0};
    expected.LossyThinSurfs->cs[1].sigmam_devia = {0.0, 0.0, 0.0};

    return expected;
}

#endif // TEST_READ_SGBC_H
