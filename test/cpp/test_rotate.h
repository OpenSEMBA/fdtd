#ifndef TEST_ROTATE_H
#define TEST_ROTATE_H

#include <cmath>
#include <complex>
#include <gtest/gtest.h>

#include "nfde_rotate_m.h"
#include "test_rotate_helpers.h"

using namespace NFDETypes_m;
using rotate_test::kRotTol;

TEST(rotate, rotate_spaceSteps_test) {
    Parseador_t this_obj = rotate_test::makeSpaceStepsProblem();
    int mpidir = 2;

    this_obj.matriz->totalX = 10;
    this_obj.matriz->totalY = 20;
    this_obj.matriz->totalZ = 30;
    this_obj.despl->nX = 5;
    this_obj.despl->nY = 15;
    this_obj.despl->nZ = 25;
    this_obj.despl->mx1 = 1;
    this_obj.despl->my1 = 11;
    this_obj.despl->mz1 = 21;
    this_obj.despl->mx2 = 6;
    this_obj.despl->my2 = 16;
    this_obj.despl->mz2 = 26;

    nfde_rotate_m::rotate_generateSpaceSteps(this_obj, mpidir);

    EXPECT_EQ(this_obj.matriz->totalX, 30);
    EXPECT_EQ(this_obj.matriz->totalY, 10);
    EXPECT_EQ(this_obj.matriz->totalZ, 20);
    EXPECT_EQ(this_obj.despl->nX, 25);
    EXPECT_EQ(this_obj.despl->nY, 5);
    EXPECT_EQ(this_obj.despl->nZ, 15);
    EXPECT_EQ(this_obj.despl->mx1, 21);
    EXPECT_EQ(this_obj.despl->my1, 1);
    EXPECT_EQ(this_obj.despl->mz1, 11);
    EXPECT_EQ(this_obj.despl->mx2, 26);
    EXPECT_EQ(this_obj.despl->my2, 6);
    EXPECT_EQ(this_obj.despl->mz2, 16);

    mpidir = 1;
    this_obj.matriz->totalX = 10;
    this_obj.matriz->totalY = 20;
    this_obj.matriz->totalZ = 30;
    this_obj.despl->nX = 5;
    this_obj.despl->nY = 15;
    this_obj.despl->nZ = 25;
    this_obj.despl->mx1 = 1;
    this_obj.despl->my1 = 11;
    this_obj.despl->mz1 = 21;
    this_obj.despl->mx2 = 6;
    this_obj.despl->my2 = 16;
    this_obj.despl->mz2 = 26;

    nfde_rotate_m::rotate_generateSpaceSteps(this_obj, mpidir);

    EXPECT_EQ(this_obj.matriz->totalX, 20);
    EXPECT_EQ(this_obj.matriz->totalY, 30);
    EXPECT_EQ(this_obj.matriz->totalZ, 10);
    EXPECT_EQ(this_obj.despl->nX, 15);
    EXPECT_EQ(this_obj.despl->nY, 25);
    EXPECT_EQ(this_obj.despl->nZ, 5);
    EXPECT_EQ(this_obj.despl->mx1, 11);
    EXPECT_EQ(this_obj.despl->my1, 21);
    EXPECT_EQ(this_obj.despl->mz1, 1);
    EXPECT_EQ(this_obj.despl->mx2, 16);
    EXPECT_EQ(this_obj.despl->my2, 26);
    EXPECT_EQ(this_obj.despl->mz2, 6);

    rotate_test::freeSpaceStepsProblem(this_obj);
}

TEST(rotate, rotate_currentFieldSources_test) {
    Parseador_t this_obj;
    rotate_test::setupCurrentFieldSources(this_obj);

    nfde_rotate_m::rotate_generateCurrent_Field_Sources(this_obj, 2);

    EXPECT_EQ(this_obj.nodSrc->NodalSource[0].c1P[0].Xi, 3);
    EXPECT_EQ(this_obj.nodSrc->NodalSource[0].c1P[0].Yi, 1);
    EXPECT_EQ(this_obj.nodSrc->NodalSource[0].c1P[0].Zi, 2);
    EXPECT_EQ(this_obj.nodSrc->NodalSource[0].c1P[0].Or, iEy);
    EXPECT_EQ(this_obj.nodSrc->NodalSource[0].c2P[0].Xi, 6);
    EXPECT_EQ(this_obj.nodSrc->NodalSource[0].c2P[0].Yi, 4);
    EXPECT_EQ(this_obj.nodSrc->NodalSource[0].c2P[0].Zi, 5);
    EXPECT_EQ(this_obj.nodSrc->NodalSource[0].c2P[0].Or, iEz);

    rotate_test::freeCurrentFieldSources(this_obj);
}

TEST(rotate, rotate_planeWaves_test) {
    Parseador_t this_obj;
    rotate_test::setupPlaneWaves(this_obj);
    const double theta = this_obj.plnSrc->collection[0].theta;
    const double phi = this_obj.plnSrc->collection[0].phi;
    const double alpha = this_obj.plnSrc->collection[0].alpha;
    const double beta = this_obj.plnSrc->collection[0].beta;

    nfde_rotate_m::rotate_generatePlaneWaves(this_obj, 2);

    EXPECT_EQ(this_obj.plnSrc->collection[0].coor1[0], 3);
    EXPECT_EQ(this_obj.plnSrc->collection[0].coor1[1], 1);
    EXPECT_EQ(this_obj.plnSrc->collection[0].coor1[2], 2);
    EXPECT_EQ(this_obj.plnSrc->collection[0].coor2[0], 6);
    EXPECT_EQ(this_obj.plnSrc->collection[0].coor2[1], 4);
    EXPECT_EQ(this_obj.plnSrc->collection[0].coor2[2], 5);
    EXPECT_NEAR(this_obj.plnSrc->collection[0].theta,
              std::atan2(std::sqrt(std::cos(theta) * std::cos(theta) +
                                   std::cos(phi) * std::cos(phi) * std::sin(theta) * std::sin(theta)),
                         std::sin(phi) * std::sin(theta)),
              kRotTol);
    EXPECT_NEAR(this_obj.plnSrc->collection[0].phi,
              std::atan2(std::cos(phi) * std::sin(theta), std::cos(theta)), kRotTol);
    EXPECT_NEAR(this_obj.plnSrc->collection[0].alpha,
              std::atan2(std::sqrt(std::cos(alpha) * std::cos(alpha) +
                                   std::cos(beta) * std::cos(beta) * std::sin(alpha) * std::sin(alpha)),
                         std::sin(beta) * std::sin(alpha)),
              kRotTol);
    EXPECT_NEAR(this_obj.plnSrc->collection[0].beta,
              std::atan2(std::cos(beta) * std::sin(alpha), std::cos(alpha)), kRotTol);

    rotate_test::freePlaneWaves(this_obj);
    rotate_test::setupPlaneWaves(this_obj);

    nfde_rotate_m::rotate_generatePlaneWaves(this_obj, 1);

    EXPECT_EQ(this_obj.plnSrc->collection[0].coor1[0], 2);
    EXPECT_EQ(this_obj.plnSrc->collection[0].coor1[1], 3);
    EXPECT_EQ(this_obj.plnSrc->collection[0].coor1[2], 1);
    EXPECT_EQ(this_obj.plnSrc->collection[0].coor2[0], 5);
    EXPECT_EQ(this_obj.plnSrc->collection[0].coor2[1], 6);
    EXPECT_EQ(this_obj.plnSrc->collection[0].coor2[2], 4);
    EXPECT_NEAR(this_obj.plnSrc->collection[0].theta,
              std::atan2(std::sqrt(std::cos(theta) * std::cos(theta) +
                                   std::sin(phi) * std::sin(phi) * std::sin(theta) * std::sin(theta)),
                         std::cos(phi) * std::sin(theta)),
              kRotTol);
    EXPECT_NEAR(this_obj.plnSrc->collection[0].phi,
              std::atan2(std::cos(theta), std::sin(phi) * std::sin(theta)), kRotTol);
    EXPECT_NEAR(this_obj.plnSrc->collection[0].alpha,
              std::atan2(std::sqrt(std::cos(alpha) * std::cos(alpha) +
                                   std::sin(beta) * std::sin(beta) * std::sin(alpha) * std::sin(alpha)),
                         std::cos(beta) * std::sin(alpha)),
              kRotTol);
    EXPECT_NEAR(this_obj.plnSrc->collection[0].beta,
              std::atan2(std::cos(alpha), std::sin(beta) * std::sin(alpha)), kRotTol);

    rotate_test::freePlaneWaves(this_obj);
}

TEST(rotate, rotate_boxSources_test) {
    Parseador_t this_obj;
    rotate_test::setupBoxSources(this_obj);

    nfde_rotate_m::rotate_generateBoxSources(this_obj, 2);
    EXPECT_EQ(this_obj.boxSrc->Vols[0].coor1[0], 3);
    EXPECT_EQ(this_obj.boxSrc->Vols[0].coor1[1], 1);
    EXPECT_EQ(this_obj.boxSrc->Vols[0].coor1[2], 2);
    EXPECT_EQ(this_obj.boxSrc->Vols[0].coor2[0], 6);
    EXPECT_EQ(this_obj.boxSrc->Vols[0].coor2[1], 4);
    EXPECT_EQ(this_obj.boxSrc->Vols[0].coor2[2], 5);

    rotate_test::freeBoxSources(this_obj);
    rotate_test::setupBoxSources(this_obj);

    nfde_rotate_m::rotate_generateBoxSources(this_obj, 1);
    EXPECT_EQ(this_obj.boxSrc->Vols[0].coor1[0], 2);
    EXPECT_EQ(this_obj.boxSrc->Vols[0].coor1[1], 3);
    EXPECT_EQ(this_obj.boxSrc->Vols[0].coor1[2], 1);
    EXPECT_EQ(this_obj.boxSrc->Vols[0].coor2[0], 5);
    EXPECT_EQ(this_obj.boxSrc->Vols[0].coor2[1], 6);
    EXPECT_EQ(this_obj.boxSrc->Vols[0].coor2[2], 4);

    rotate_test::freeBoxSources(this_obj);
}

TEST(rotate, rotate_fronteras_test) {
    Parseador_t this_obj;
    rotate_test::setupFronteras(this_obj);

    nfde_rotate_m::rotate_generateFronteras(this_obj, 2);
    EXPECT_EQ(this_obj.front->tipoFrontera[0], 5);
    EXPECT_EQ(this_obj.front->tipoFrontera[1], 6);
    EXPECT_EQ(this_obj.front->tipoFrontera[2], 1);
    EXPECT_EQ(this_obj.front->tipoFrontera[3], 2);
    EXPECT_EQ(this_obj.front->tipoFrontera[4], 3);
    EXPECT_EQ(this_obj.front->tipoFrontera[5], 4);
    EXPECT_DOUBLE_EQ(this_obj.front->propiedadesPML[0].orden, 5.0);
    EXPECT_DOUBLE_EQ(this_obj.front->propiedadesPML[1].orden, 6.0);
    EXPECT_DOUBLE_EQ(this_obj.front->propiedadesPML[2].orden, 1.0);
    EXPECT_DOUBLE_EQ(this_obj.front->propiedadesPML[3].orden, 2.0);
    EXPECT_DOUBLE_EQ(this_obj.front->propiedadesPML[4].orden, 3.0);
    EXPECT_DOUBLE_EQ(this_obj.front->propiedadesPML[5].orden, 4.0);
    EXPECT_DOUBLE_EQ(this_obj.front->propiedadesPML[0].refl, 0.5);
    EXPECT_DOUBLE_EQ(this_obj.front->propiedadesPML[1].refl, 0.6);
    EXPECT_EQ(this_obj.front->propiedadesPML[0].numCapas, 50);
    EXPECT_EQ(this_obj.front->propiedadesPML[5].numCapas, 40);

    rotate_test::freeFronteras(this_obj);
    rotate_test::setupFronteras(this_obj);

    nfde_rotate_m::rotate_generateFronteras(this_obj, 1);
    EXPECT_EQ(this_obj.front->tipoFrontera[0], 3);
    EXPECT_EQ(this_obj.front->tipoFrontera[5], 2);
    EXPECT_DOUBLE_EQ(this_obj.front->propiedadesPML[0].orden, 3.0);
    EXPECT_DOUBLE_EQ(this_obj.front->propiedadesPML[5].orden, 2.0);
    EXPECT_EQ(this_obj.front->propiedadesPML[0].numCapas, 30);
    EXPECT_EQ(this_obj.front->propiedadesPML[5].numCapas, 20);

    rotate_test::freeFronteras(this_obj);
}

TEST(rotate, rotate_pecs_test) {
    Parseador_t this_obj;
    rotate_test::setupPECs(this_obj);

    nfde_rotate_m::rotate_generatePECs(this_obj, 2);
    EXPECT_EQ(this_obj.pecRegs->Vols[0].Xi, 3);
    EXPECT_EQ(this_obj.pecRegs->Vols[0].Yi, 1);
    EXPECT_EQ(this_obj.pecRegs->Vols[0].Zi, 2);
    EXPECT_EQ(this_obj.pecRegs->Vols[0].Or, iEy);
    EXPECT_EQ(this_obj.pecRegs->Surfs[0].Xi, 9);
    EXPECT_EQ(this_obj.pecRegs->Surfs[0].Or, iEz);
    EXPECT_EQ(this_obj.pecRegs->Lins[0].Xi, 15);
    EXPECT_EQ(this_obj.pecRegs->Lins[0].Or, iEx);

    rotate_test::freePECs(this_obj);
    rotate_test::setupPECs(this_obj);

    nfde_rotate_m::rotate_generatePECs(this_obj, 1);
    EXPECT_EQ(this_obj.pecRegs->Vols[0].Xi, 2);
    EXPECT_EQ(this_obj.pecRegs->Vols[0].Or, iEz);
    EXPECT_EQ(this_obj.pecRegs->Surfs[0].Or, iEx);
    EXPECT_EQ(this_obj.pecRegs->Lins[0].Or, iEy);

    rotate_test::freePECs(this_obj);
}

TEST(rotate, rotate_nonMetals_test) {
    Parseador_t this_obj;
    rotate_test::setupNONMetals(this_obj);

    nfde_rotate_m::rotate_generateNONMetals(this_obj, 2);
    EXPECT_EQ(this_obj.DielRegs->Vols[0].c1P[0].Xi, 3);
    EXPECT_EQ(this_obj.DielRegs->Vols[0].c1P[0].Or, iEy);
    EXPECT_EQ(this_obj.DielRegs->Vols[0].c2P[1].Or, iEy);
    EXPECT_EQ(this_obj.DielRegs->Lins[0].c2P[1].Or, iEx);

    rotate_test::freeNONMetals(this_obj);
    rotate_test::setupNONMetals(this_obj);

    nfde_rotate_m::rotate_generateNONMetals(this_obj, 1);
    EXPECT_EQ(this_obj.DielRegs->Vols[0].c1P[0].Xi, 2);
    EXPECT_EQ(this_obj.DielRegs->Vols[0].c1P[0].Or, iEz);
    EXPECT_EQ(this_obj.DielRegs->Surfs[0].c2P[1].Or, iEx);
    EXPECT_EQ(this_obj.DielRegs->Lins[0].c2P[1].Or, iEy);

    rotate_test::freeNONMetals(this_obj);
}

TEST(rotate, rotate_thinWires_test) {
    Parseador_t this_obj;
    rotate_test::setupThinWires(this_obj);

    nfde_rotate_m::rotate_generateThinWires(this_obj, 2);
    EXPECT_EQ(this_obj.tWires->tw[0].twc[0].i, 3);
    EXPECT_EQ(this_obj.tWires->tw[0].twc[0].j, 1);
    EXPECT_EQ(this_obj.tWires->tw[0].twc[0].K, 2);
    EXPECT_EQ(this_obj.tWires->tw[0].twc[0].d, iEy);
    EXPECT_EQ(this_obj.tWires->tw[0].twc[1].i, 7);
    EXPECT_EQ(this_obj.tWires->tw[0].twc[1].d, iEx);

    rotate_test::freeThinWires(this_obj);
    rotate_test::setupThinWires(this_obj);

    nfde_rotate_m::rotate_generateThinWires(this_obj, 1);
    EXPECT_EQ(this_obj.tWires->tw[0].twc[0].i, 2);
    EXPECT_EQ(this_obj.tWires->tw[0].twc[0].K, 1);
    EXPECT_EQ(this_obj.tWires->tw[0].twc[0].d, iEz);
    EXPECT_EQ(this_obj.tWires->tw[0].twc[1].d, iEy);

    rotate_test::freeThinWires(this_obj);
}

TEST(rotate, rotate_slantedWires_test) {
    Parseador_t this_obj;
    rotate_test::setupSlantedWires(this_obj);

    nfde_rotate_m::rotate_generateSlantedWires(this_obj, 2);
    EXPECT_DOUBLE_EQ(this_obj.sWires->sw[0].swc[0].x, 3.0);
    EXPECT_DOUBLE_EQ(this_obj.sWires->sw[0].swc[0].y, 1.0);
    EXPECT_DOUBLE_EQ(this_obj.sWires->sw[0].swc[0].z, 2.0);
    EXPECT_DOUBLE_EQ(this_obj.sWires->sw[0].swc[1].x, 7.0);

    rotate_test::freeSlantedWires(this_obj);
    rotate_test::setupSlantedWires(this_obj);

    nfde_rotate_m::rotate_generateSlantedWires(this_obj, 1);
    EXPECT_DOUBLE_EQ(this_obj.sWires->sw[0].swc[0].x, 2.0);
    EXPECT_DOUBLE_EQ(this_obj.sWires->sw[0].swc[0].y, 3.0);
    EXPECT_DOUBLE_EQ(this_obj.sWires->sw[0].swc[0].z, 1.0);

    rotate_test::freeSlantedWires(this_obj);
}

TEST(rotate, rotate_thinSlots_test) {
    Parseador_t this_obj;
    rotate_test::setupThinSlots(this_obj);

    nfde_rotate_m::rotate_generateThinSlots(this_obj, 2);
    EXPECT_EQ(this_obj.tSlots->tg[0].tgc[0].i, 3);
    EXPECT_EQ(this_obj.tSlots->tg[0].tgc[0].dir, iEy);
    EXPECT_EQ(this_obj.tSlots->tg[0].tgc[1].dir, iEz);

    rotate_test::freeThinSlots(this_obj);
    rotate_test::setupThinSlots(this_obj);

    nfde_rotate_m::rotate_generateThinSlots(this_obj, 1);
    EXPECT_EQ(this_obj.tSlots->tg[0].tgc[0].i, 2);
    EXPECT_EQ(this_obj.tSlots->tg[0].tgc[0].dir, iEz);
    EXPECT_EQ(this_obj.tSlots->tg[0].tgc[1].dir, iEx);

    rotate_test::freeThinSlots(this_obj);
}

TEST(rotate, rotate_lossyThinSurface_test) {
    Parseador_t this_obj;
    rotate_test::setupLossyThinSurface(this_obj);

    nfde_rotate_m::rotate_generateLossyThinSurface(this_obj, 2);
    EXPECT_EQ(this_obj.LossyThinSurfs->cs[0].c[0].Xi, 3);
    EXPECT_EQ(this_obj.LossyThinSurfs->cs[0].c[0].Or, iEy);
    EXPECT_EQ(this_obj.LossyThinSurfs->cs[0].c[1].Or, iEz);

    rotate_test::freeLossyThinSurface(this_obj);
    rotate_test::setupLossyThinSurface(this_obj);

    nfde_rotate_m::rotate_generateLossyThinSurface(this_obj, 1);
    EXPECT_EQ(this_obj.LossyThinSurfs->cs[0].c[0].Xi, 2);
    EXPECT_EQ(this_obj.LossyThinSurfs->cs[0].c[0].Or, iEz);
    EXPECT_EQ(this_obj.LossyThinSurfs->cs[0].c[1].Or, iEx);

    rotate_test::freeLossyThinSurface(this_obj);
}

TEST(rotate, rotate_freqDependMaterials_test) {
    Parseador_t this_obj;
    rotate_test::setupFDMs(this_obj);

    nfde_rotate_m::rotate_generateFDMs(this_obj, 2);
    EXPECT_EQ(this_obj.frqDepMats->Vols[0].a11[0], std::complex<double>(33, 0));
    EXPECT_EQ(this_obj.frqDepMats->Vols[0].a33[0], std::complex<double>(22, 0));
    EXPECT_DOUBLE_EQ(this_obj.frqDepMats->Vols[0].eps11, 33.0);
    EXPECT_DOUBLE_EQ(this_obj.frqDepMats->Vols[0].eps33, 22.0);
    EXPECT_EQ(this_obj.frqDepMats->Vols[0].K11, 33);
    EXPECT_EQ(this_obj.frqDepMats->Vols[0].Km33, 22);

    rotate_test::freeFDMs(this_obj);
    rotate_test::setupFDMs(this_obj);

    nfde_rotate_m::rotate_generateFDMs(this_obj, 1);
    EXPECT_EQ(this_obj.frqDepMats->Vols[0].a11[0], std::complex<double>(22, 0));
    EXPECT_EQ(this_obj.frqDepMats->Vols[0].a33[0], std::complex<double>(11, 0));
    EXPECT_DOUBLE_EQ(this_obj.frqDepMats->Vols[0].mu11, 22.0);
    EXPECT_EQ(this_obj.frqDepMats->Vols[0].Km11, 22);

    rotate_test::freeFDMs(this_obj);
}

TEST(rotate, rotate_sondas_test) {
    Parseador_t this_obj;
    rotate_test::setupSONDAs(this_obj);
    const double thetaStart = 1.5;
    const double phiStart = 1.5;
    const double thetaStop = 2.0;
    const double phiStop = 2.0;

    nfde_rotate_m::rotate_generateSONDAs(this_obj, 2);
    EXPECT_EQ(this_obj.oldSONDA->probes[0].Electric[0].probe.i[0], 3);
    EXPECT_EQ(this_obj.oldSONDA->probes[0].Electric[0].probe.j[0], 1);
    EXPECT_EQ(this_obj.oldSONDA->probes[0].Electric[0].probe.K[0], 2);
    EXPECT_NEAR(this_obj.oldSONDA->probes[0].FarField[0].probe.thetastart,
                std::atan2(std::sqrt(std::cos(thetaStart) * std::cos(thetaStart) +
                                     std::cos(phiStart) * std::cos(phiStart) *
                                         std::sin(thetaStart) * std::sin(thetaStart)),
                           std::sin(phiStart) * std::sin(thetaStart)),
                kRotTol);

    rotate_test::freeSONDAs(this_obj);
    rotate_test::setupSONDAs(this_obj);

    nfde_rotate_m::rotate_generateSONDAs(this_obj, 1);
    EXPECT_EQ(this_obj.oldSONDA->probes[0].Magnetic[0].probe.i[0], 2);
    EXPECT_EQ(this_obj.oldSONDA->probes[0].Magnetic[0].probe.j[0], 3);
    EXPECT_EQ(this_obj.oldSONDA->probes[0].Magnetic[0].probe.K[0], 1);
    EXPECT_NEAR(this_obj.oldSONDA->probes[0].FarField[0].probe.phistart,
                std::atan2(std::cos(thetaStart), std::sin(phiStart) * std::sin(thetaStart)),
                kRotTol);

    rotate_test::freeSONDAs(this_obj);
}

TEST(rotate, rotate_masSondas_test) {
    Parseador_t this_obj;
    rotate_test::setupMasSondas(this_obj, 2);

    nfde_rotate_m::rotate_generateMasSondas(this_obj, 2);
    EXPECT_EQ(this_obj.Sonda->collection[0].cordinates[0].Xi, 3);
    EXPECT_EQ(this_obj.Sonda->collection[0].cordinates[0].Or, NP_COR_EZ);
    EXPECT_EQ(this_obj.Sonda->collection[1].cordinates[0].Or, NP_COR_EX);

    rotate_test::freeMasSondas(this_obj);
    rotate_test::setupMasSondas(this_obj, 2);

    nfde_rotate_m::rotate_generateMasSondas(this_obj, 1);
    EXPECT_EQ(this_obj.Sonda->collection[0].cordinates[0].Or, NP_COR_EX);
    EXPECT_EQ(this_obj.Sonda->collection[1].cordinates[0].Or, NP_COR_EY);

    rotate_test::freeMasSondas(this_obj);
}

TEST(rotate, rotate_bloqueProbes_test) {
    Parseador_t this_obj;
    rotate_test::setupBloqueProbes(this_obj, 2);

    nfde_rotate_m::rotate_generateBloqueProbes(this_obj, 2);
    EXPECT_EQ(this_obj.BloquePrb->bp[0].i1, 3);
    EXPECT_EQ(this_obj.BloquePrb->bp[0].j1, 1);
    EXPECT_EQ(this_obj.BloquePrb->bp[1].i2, 10);

    rotate_test::freeBloqueProbes(this_obj);
    rotate_test::setupBloqueProbes(this_obj, 2);

    nfde_rotate_m::rotate_generateBloqueProbes(this_obj, 1);
    EXPECT_EQ(this_obj.BloquePrb->bp[0].i1, 2);
    EXPECT_EQ(this_obj.BloquePrb->bp[0].j1, 3);
    EXPECT_EQ(this_obj.BloquePrb->bp[1].k1, 5);

    rotate_test::freeBloqueProbes(this_obj);
}

TEST(rotate, rotate_volumicProbes_test) {
    Parseador_t this_obj;
    rotate_test::setupVolumicProbes(this_obj, 2);

    nfde_rotate_m::rotate_generateVolumicProbes(this_obj, 2);
    EXPECT_EQ(this_obj.VolPrb->collection[0].cordinates[0].Xi, 1);
    EXPECT_EQ(this_obj.VolPrb->collection[0].cordinates[0].Or, 1);
    EXPECT_EQ(this_obj.VolPrb->collection[1].cordinates[0].Ze, 10);

    rotate_test::freeVolumicProbes(this_obj);
    rotate_test::setupVolumicProbes(this_obj, 2);

    nfde_rotate_m::rotate_generateVolumicProbes(this_obj, 1);
    EXPECT_EQ(this_obj.VolPrb->collection[0].cordinates[0].Yi, 2);
    EXPECT_EQ(this_obj.VolPrb->collection[1].cordinates[0].Or, 2);

    rotate_test::freeVolumicProbes(this_obj);
}

#ifdef CompileWithMTLN
TEST(rotate, rotate_mtln_test) {
    Parseador_t this_obj;
    rotate_test::setupMtln(this_obj);

    nfde_rotate_m::rotate_mtln(this_obj, 2);
    EXPECT_EQ(this_obj.mtln->cables[0].ptr->segments[0].x, 3);
    EXPECT_EQ(this_obj.mtln->cables[0].ptr->segments[0].y, 1);
    EXPECT_EQ(this_obj.mtln->cables[0].ptr->segments[0].z, 2);
    EXPECT_EQ(this_obj.mtln->cables[0].ptr->segments[0].orientation, 2);

    rotate_test::resetMtlnSegment(this_obj);
    nfde_rotate_m::rotate_mtln(this_obj, 1);
    EXPECT_EQ(this_obj.mtln->cables[0].ptr->segments[0].x, 2);
    EXPECT_EQ(this_obj.mtln->cables[0].ptr->segments[0].y, 3);
    EXPECT_EQ(this_obj.mtln->cables[0].ptr->segments[0].z, 1);
    EXPECT_EQ(this_obj.mtln->cables[0].ptr->segments[0].orientation, 3);

    rotate_test::freeMtln(this_obj);
}
#endif

#endif
