#ifndef TEST_ROTATE_HELPERS_H
#define TEST_ROTATE_HELPERS_H

#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include "nfde_types.h"

namespace rotate_test {

constexpr double kRotTol = 1e-10;

using namespace NFDETypes_m;

inline Parseador_t makeSpaceStepsProblem() {
    Parseador_t p;
    p.despl = new Desplazamiento_t();
    p.matriz = new MatrizMedios_t();
    return p;
}

inline void freeSpaceStepsProblem(Parseador_t& p) {
    delete p.despl;
    delete p.matriz;
    p.despl = nullptr;
    p.matriz = nullptr;
}

inline void setupCurrentFieldSources(Parseador_t& p) {
    p.nodSrc = new NodSource_t();
    p.nodSrc->n_nodSrc = 1;
    p.nodSrc->NodalSource.resize(1);
    p.nodSrc->NodalSource[0].n_C1P = 1;
    p.nodSrc->NodalSource[0].n_C2P = 1;
    p.nodSrc->NodalSource[0].c1P.resize(1);
    p.nodSrc->NodalSource[0].c2P.resize(1);
    p.nodSrc->NodalSource[0].c1P[0].Xi = 1;
    p.nodSrc->NodalSource[0].c1P[0].Yi = 2;
    p.nodSrc->NodalSource[0].c1P[0].Zi = 3;
    p.nodSrc->NodalSource[0].c1P[0].Or = iEx;
    p.nodSrc->NodalSource[0].c2P[0].Xi = 4;
    p.nodSrc->NodalSource[0].c2P[0].Yi = 5;
    p.nodSrc->NodalSource[0].c2P[0].Zi = 6;
    p.nodSrc->NodalSource[0].c2P[0].Or = iEy;
}

inline void freeCurrentFieldSources(Parseador_t& p) {
    delete p.nodSrc;
    p.nodSrc = nullptr;
}

inline void setupPlaneWaves(Parseador_t& p) {
    p.plnSrc = new PlaneWaves_t();
    p.plnSrc->nc = 1;
    p.plnSrc->collection.resize(1);
    p.plnSrc->collection[0].coor1[0] = 1;
    p.plnSrc->collection[0].coor1[1] = 2;
    p.plnSrc->collection[0].coor1[2] = 3;
    p.plnSrc->collection[0].coor2[0] = 4;
    p.plnSrc->collection[0].coor2[1] = 5;
    p.plnSrc->collection[0].coor2[2] = 6;
    p.plnSrc->collection[0].theta = 0.5;
    p.plnSrc->collection[0].phi = 0.3;
    p.plnSrc->collection[0].alpha = 0.7;
    p.plnSrc->collection[0].beta = 0.2;
}

inline void freePlaneWaves(Parseador_t& p) {
    delete p.plnSrc;
    p.plnSrc = nullptr;
}

inline void setupBoxSources(Parseador_t& p) {
    p.boxSrc = new Boxes_t();
    p.boxSrc->nVols = 1;
    p.boxSrc->Vols.resize(1);
    p.boxSrc->Vols[0].coor1[0] = 1;
    p.boxSrc->Vols[0].coor1[1] = 2;
    p.boxSrc->Vols[0].coor1[2] = 3;
    p.boxSrc->Vols[0].coor2[0] = 4;
    p.boxSrc->Vols[0].coor2[1] = 5;
    p.boxSrc->Vols[0].coor2[2] = 6;
}

inline void freeBoxSources(Parseador_t& p) {
    delete p.boxSrc;
    p.boxSrc = nullptr;
}

inline void setupFronteras(Parseador_t& p) {
    p.front = new Frontera_t();
    for (int i = 0; i < 6; ++i) {
        p.front->tipoFrontera[i] = i + 1;
        p.front->propiedadesPML[i].orden = static_cast<double>(i + 1);
        p.front->propiedadesPML[i].refl = 0.1 * (i + 1);
        p.front->propiedadesPML[i].numCapas = (i + 1) * 10;
    }
}

inline void freeFronteras(Parseador_t& p) {
    delete p.front;
    p.front = nullptr;
}

inline void setupPECs(Parseador_t& p) {
    p.pecRegs = new PECRegions_t();
    p.pecRegs->nVols = 1;
    p.pecRegs->nSurfs = 1;
    p.pecRegs->nLins = 1;
    p.pecRegs->Vols.resize(1);
    p.pecRegs->Surfs.resize(1);
    p.pecRegs->Lins.resize(1);
    p.pecRegs->Vols[0] = {1, 4, 2, 5, 3, 6, 1, 1, 1, iEx, ""};
    p.pecRegs->Surfs[0] = {7, 10, 8, 11, 9, 12, 1, 1, 1, iEy, ""};
    p.pecRegs->Lins[0] = {13, 16, 14, 17, 15, 18, 1, 1, 1, iEz, ""};
}

inline void freePECs(Parseador_t& p) {
    delete p.pecRegs;
    p.pecRegs = nullptr;
}

inline coords_t makeCoord(int xi, int yi, int zi, int or_val) {
    coords_t c;
    c.Xi = xi;
    c.Yi = yi;
    c.Zi = zi;
    c.Or = or_val;
    return c;
}

inline void setupNONMetals(Parseador_t& p) {
    p.DielRegs = new DielectricRegions_t();
    p.DielRegs->nVols = 1;
    p.DielRegs->nSurfs = 1;
    p.DielRegs->nLins = 1;
    p.DielRegs->Vols.resize(1);
    p.DielRegs->Surfs.resize(1);
    p.DielRegs->Lins.resize(1);

    auto& vol = p.DielRegs->Vols[0];
    vol.n_C1P = 2;
    vol.n_C2P = 2;
    vol.c1P = {makeCoord(1, 2, 3, iEx), makeCoord(4, 5, 6, iEy)};
    vol.c2P = {makeCoord(7, 8, 9, iEz), makeCoord(10, 11, 12, iEx)};

    auto& surf = p.DielRegs->Surfs[0];
    surf.n_C1P = 2;
    surf.n_C2P = 2;
    surf.c1P = {makeCoord(13, 14, 15, iEy), makeCoord(16, 17, 18, iEz)};
    surf.c2P = {makeCoord(19, 20, 21, iEx), makeCoord(22, 23, 24, iEy)};

    auto& lin = p.DielRegs->Lins[0];
    lin.n_C1P = 2;
    lin.n_C2P = 2;
    lin.c1P = {makeCoord(25, 26, 27, iEz), makeCoord(28, 29, 30, iEx)};
    lin.c2P = {makeCoord(31, 32, 33, iEy), makeCoord(34, 35, 36, iEz)};
}

inline void freeNONMetals(Parseador_t& p) {
    delete p.DielRegs;
    p.DielRegs = nullptr;
}

inline void setupThinWires(Parseador_t& p) {
    p.tWires = new ThinWires_t();
    p.tWires->n_tw = 1;
    p.tWires->tw.resize(1);
    p.tWires->tw[0].n_twc = 2;
    p.tWires->tw[0].twc.resize(2);
    p.tWires->tw[0].twc[0] = {"type1", "source.dat", 1, 2, 3, 1, iEx, 0.5, "wire1"};
    p.tWires->tw[0].twc[1] = {"type1", "source.dat", 5, 6, 7, 2, iEz, 0.8, "wire2"};
}

inline void freeThinWires(Parseador_t& p) {
    delete p.tWires;
    p.tWires = nullptr;
}

inline void setupSlantedWires(Parseador_t& p) {
    p.sWires = new SlantedWiresInfo_t();
    p.sWires->n_sw = 1;
    p.sWires->sw.resize(1);
    p.sWires->sw[0].n_swc = 2;
    p.sWires->sw[0].swc.resize(2);
    p.sWires->sw[0].swc[0] = {"type1", "source.dat", 1.0, 2.0, 3.0, 1, 0.5, "wire1"};
    p.sWires->sw[0].swc[1] = {"type1", "source.dat", 5.0, 6.0, 7.0, 2, 0.8, "wire2"};
}

inline void freeSlantedWires(Parseador_t& p) {
    delete p.sWires;
    p.sWires = nullptr;
}

inline void setupThinSlots(Parseador_t& p) {
    p.tSlots = new ThinSlots_t();
    p.tSlots->n_tg = 1;
    p.tSlots->tg.resize(1);
    p.tSlots->tg[0].n_tgc = 2;
    p.tSlots->tg[0].tgc.resize(2);
    p.tSlots->tg[0].tgc[0] = {1, 2, 3, 1, iEx, 1, "slot1"};
    p.tSlots->tg[0].tgc[1] = {4, 5, 6, 2, iEy, 1, "slot1"};
}

inline void freeThinSlots(Parseador_t& p) {
    delete p.tSlots;
    p.tSlots = nullptr;
}

inline void setupLossyThinSurface(Parseador_t& p) {
    p.LossyThinSurfs = new LossyThinSurfaces_t();
    p.LossyThinSurfs->length = 1;
    p.LossyThinSurfs->cs.resize(1);
    p.LossyThinSurfs->cs[0].nc = 2;
    p.LossyThinSurfs->cs[0].c.resize(2);
    p.LossyThinSurfs->cs[0].c[0] = {1, 4, 2, 5, 3, 6, 1, 2, 3, iEx, "tag1"};
    p.LossyThinSurfs->cs[0].c[1] = {7, 10, 8, 11, 9, 12, 4, 5, 6, iEy, "tag2"};
}

inline void freeLossyThinSurface(Parseador_t& p) {
    delete p.LossyThinSurfs;
    p.LossyThinSurfs = nullptr;
}

inline void setupFDMs(Parseador_t& p) {
    p.frqDepMats = new FreqDepenMaterials_t();
    p.frqDepMats->nVols = 1;
    p.frqDepMats->Vols.resize(1);
    auto& v = p.frqDepMats->Vols[0];
    v.a11 = {std::complex<double>(11, 0)};
    v.a12 = {std::complex<double>(12, 0)};
    v.a13 = {std::complex<double>(13, 0)};
    v.a22 = {std::complex<double>(22, 0)};
    v.a23 = {std::complex<double>(23, 0)};
    v.a33 = {std::complex<double>(33, 0)};
    v.am11 = {std::complex<double>(11, 0)};
    v.am12 = {std::complex<double>(12, 0)};
    v.am13 = {std::complex<double>(13, 0)};
    v.am22 = {std::complex<double>(22, 0)};
    v.am23 = {std::complex<double>(23, 0)};
    v.am33 = {std::complex<double>(33, 0)};
    v.b11 = {std::complex<double>(11, 0)};
    v.b12 = {std::complex<double>(12, 0)};
    v.b13 = {std::complex<double>(13, 0)};
    v.b22 = {std::complex<double>(22, 0)};
    v.b23 = {std::complex<double>(23, 0)};
    v.b33 = {std::complex<double>(33, 0)};
    v.bm11 = {std::complex<double>(11, 0)};
    v.bm12 = {std::complex<double>(12, 0)};
    v.bm13 = {std::complex<double>(13, 0)};
    v.bm22 = {std::complex<double>(22, 0)};
    v.bm23 = {std::complex<double>(23, 0)};
    v.bm33 = {std::complex<double>(33, 0)};
    v.eps11 = 11;
    v.eps12 = 12;
    v.eps13 = 13;
    v.eps22 = 22;
    v.eps23 = 23;
    v.eps33 = 33;
    v.mu11 = 11;
    v.mu12 = 12;
    v.mu13 = 13;
    v.mu22 = 22;
    v.mu23 = 23;
    v.mu33 = 33;
    v.sigma11 = 11;
    v.sigma12 = 12;
    v.sigma13 = 13;
    v.sigma22 = 22;
    v.sigma23 = 23;
    v.sigma33 = 33;
    v.sigmam11 = 11;
    v.sigmam12 = 12;
    v.sigmam13 = 13;
    v.sigmam22 = 22;
    v.sigmam23 = 23;
    v.sigmam33 = 33;
    v.K11 = 11;
    v.K12 = 12;
    v.K13 = 13;
    v.K22 = 22;
    v.K23 = 23;
    v.K33 = 33;
    v.Km11 = 11;
    v.Km12 = 12;
    v.Km13 = 13;
    v.Km22 = 22;
    v.Km23 = 23;
    v.Km33 = 33;
}

inline void freeFDMs(Parseador_t& p) {
    delete p.frqDepMats;
    p.frqDepMats = nullptr;
}

inline void setupSONDAs(Parseador_t& p) {
    p.oldSONDA = new Sondas_t();
    p.oldSONDA->n_probes = 1;
    p.oldSONDA->probes.resize(1);
    p.oldSONDA->probes[0].n_Electric = 1;
    p.oldSONDA->probes[0].n_Magnetic = 1;
    p.oldSONDA->probes[0].n_FarField = 1;
    p.oldSONDA->probes[0].Electric.resize(1);
    p.oldSONDA->probes[0].Magnetic.resize(1);
    p.oldSONDA->probes[0].FarField.resize(1);
    p.oldSONDA->probes[0].Electric[0].probe.i = {1};
    p.oldSONDA->probes[0].Electric[0].probe.j = {2};
    p.oldSONDA->probes[0].Electric[0].probe.K = {3};
    p.oldSONDA->probes[0].Electric[0].probe.n_cord = 1;
    p.oldSONDA->probes[0].Magnetic[0].probe.i = {1};
    p.oldSONDA->probes[0].Magnetic[0].probe.j = {2};
    p.oldSONDA->probes[0].Magnetic[0].probe.K = {3};
    p.oldSONDA->probes[0].Magnetic[0].probe.n_cord = 1;
    p.oldSONDA->probes[0].FarField[0].probe.thetastart = 1.5;
    p.oldSONDA->probes[0].FarField[0].probe.phistart = 1.5;
    p.oldSONDA->probes[0].FarField[0].probe.thetastop = 2.0;
    p.oldSONDA->probes[0].FarField[0].probe.phistop = 2.0;
}

inline void freeSONDAs(Parseador_t& p) {
    delete p.oldSONDA;
    p.oldSONDA = nullptr;
}

inline void initMasSondaCoords(MasSonda_t& probe, int xi, int xe, int yi, int ye, int zi, int ze,
                               int xtr, int ytr, int ztr, int or_val, int type1, int type2) {
    probe.len_cor = 1;
    probe.cordinates.resize(1);
    probe.cordinates[0].Xi = xi;
    probe.cordinates[0].Xe = xe;
    probe.cordinates[0].Yi = yi;
    probe.cordinates[0].Ye = ye;
    probe.cordinates[0].Zi = zi;
    probe.cordinates[0].Ze = ze;
    probe.cordinates[0].Xtrancos = xtr;
    probe.cordinates[0].Ytrancos = ytr;
    probe.cordinates[0].Ztrancos = ztr;
    probe.cordinates[0].Or = or_val;
    probe.type1 = type1;
    probe.type2 = type2;
}

inline void setupMasSondas(Parseador_t& p, int n_probes) {
    p.Sonda = new MasSondas_t();
    p.Sonda->length = n_probes;
    p.Sonda->collection.resize(n_probes);
    for (int i = 0; i < n_probes; ++i) {
        p.Sonda->collection[i].filename = "probe.dat";
        p.Sonda->collection[i].outputrequest = "output.txt";
    }
    initMasSondaCoords(p.Sonda->collection[0], 1, 4, 2, 5, 3, 6, 1, 2, 3, NP_COR_EY, 1, 1);
    if (n_probes > 1) {
        initMasSondaCoords(p.Sonda->collection[1], 5, 8, 6, 9, 7, 10, 1, 2, 3, NP_COR_EZ, 2, 2);
    }
}

inline void freeMasSondas(Parseador_t& p) {
    delete p.Sonda;
    p.Sonda = nullptr;
}

inline void setupBloqueProbes(Parseador_t& p, int n_probes) {
    p.BloquePrb = new BloqueProbes_t();
    p.BloquePrb->n_bp = n_probes;
    p.BloquePrb->bp.resize(n_probes);
    p.BloquePrb->bp[0] = {0, 0, 0, 0, 0, 0, "", 1, 1, 4, 2, 5, 3, 6, 1, iEx, true, "output.txt", "probe1"};
    if (n_probes > 1) {
        p.BloquePrb->bp[1] = {0, 0, 0, 0, 0, 0, "", 2, 5, 8, 6, 9, 7, 10, 2, iEx, true, "output.txt", "probe2"};
    }
}

inline void freeBloqueProbes(Parseador_t& p) {
    delete p.BloquePrb;
    p.BloquePrb = nullptr;
}

inline void setupVolumicProbes(Parseador_t& p, int n_probes) {
    p.VolPrb = new VolProbes_t();
    p.VolPrb->length = n_probes;
    p.VolPrb->collection.resize(n_probes);
    for (int i = 0; i < n_probes; ++i) {
        p.VolPrb->collection[i].len_cor = 1;
        p.VolPrb->collection[i].cordinates.resize(1);
    }
    p.VolPrb->collection[0].cordinates[0] = {1, 4, 2, 5, 3, 6, 1, 2, 3, 1, "probe1"};
    if (n_probes > 1) {
        p.VolPrb->collection[1].cordinates[0] = {5, 8, 6, 9, 7, 10, 1, 2, 3, 2, "probe2"};
    }
}

inline void freeVolumicProbes(Parseador_t& p) {
    delete p.VolPrb;
    p.VolPrb = nullptr;
}

#ifdef CompileWithMTLN
inline void setupMtln(Parseador_t& p) {
    p.mtln = new mtln_t();
    p.mtln->cables.resize(1);
    p.mtln->cables[0].ptr = std::make_unique<cable_t>();
    p.mtln->cables[0].ptr->segments.resize(1);
    p.mtln->cables[0].ptr->segments[0].x = 1;
    p.mtln->cables[0].ptr->segments[0].y = 2;
    p.mtln->cables[0].ptr->segments[0].z = 3;
    p.mtln->cables[0].ptr->segments[0].orientation = 1;
}

inline void resetMtlnSegment(Parseador_t& p) {
    p.mtln->cables[0].ptr->segments[0].x = 1;
    p.mtln->cables[0].ptr->segments[0].y = 2;
    p.mtln->cables[0].ptr->segments[0].z = 3;
    p.mtln->cables[0].ptr->segments[0].orientation = 1;
}

inline void freeMtln(Parseador_t& p) {
    delete p.mtln;
    p.mtln = nullptr;
}
#endif

} // namespace rotate_test

#endif
