#ifndef NFDE_TYPES_EXTENSION_M_H
#define NFDE_TYPES_EXTENSION_M_H

#include "../src_main_pub/nfde_types.h"

namespace NFDETypes_extension_m {

using namespace NFDETypes_m;

// ============================================================================
// initializeProblemDescription
// ============================================================================

inline void initializeProblemDescription(Parseador_t& pD) {
    pD.general = new NFDEGeneral_t();
    pD.matriz = new MatrizMedios_t();
    pD.despl = new Desplazamiento_t();
    pD.front = new Frontera_t();

    pD.Mats = new Materials_t();
    pD.Mats->n_Mats = 3;
    pD.Mats->n_Mats_max = 3;
    pD.Mats->Mats.resize(3);

    pD.Mats->Mats[0].id = 1;
    pD.Mats->Mats[0].eps = EPSILON_VACUUM;
    pD.Mats->Mats[0].mu = MU_VACUUM;
    pD.Mats->Mats[0].sigma = 0.0;
    pD.Mats->Mats[0].sigmam = 0.0;

    pD.Mats->Mats[1].id = 2;
    pD.Mats->Mats[1].eps = EPSILON_VACUUM;
    pD.Mats->Mats[1].mu = MU_VACUUM;
    pD.Mats->Mats[1].sigma = SIGMA_PEC;
    pD.Mats->Mats[1].sigmam = 0.0;

    pD.Mats->Mats[2].id = 3;
    pD.Mats->Mats[2].eps = EPSILON_VACUUM;
    pD.Mats->Mats[2].mu = MU_VACUUM;
    pD.Mats->Mats[2].sigma = 0.0;
    pD.Mats->Mats[2].sigmam = SIGMA_PMC;

    pD.pecRegs = new PECRegions_t();
    pD.pmcRegs = new PECRegions_t();
    pD.DielRegs = new DielectricRegions_t();
    pD.LossyThinSurfs = new LossyThinSurfaces_t();
    pD.plnSrc = new PlaneWaves_t();
    pD.nodSrc = new NodSource_t();
    pD.oldSONDA = new Sondas_t();
    pD.Sonda = new MasSondas_t();
    pD.BloquePrb = new BloqueProbes_t();
    pD.VolPrb = new VolProbes_t();
    pD.tWires = new ThinWires_t();
    pD.tSlots = new ThinSlots_t();
    pD.conformalRegs = new ConformalPECRegions_t();
}

// ============================================================================
// Top-level Parseador comparison (pointer-based, must be free function)
// ============================================================================

inline bool parseador_eq(const Parseador_t& a, const Parseador_t& b, bool ignoreRegions = false) {
    if (a.switches != b.switches) return false;
    if (!(*a.general == *b.general)) return false;
    if (!(*a.matriz == *b.matriz)) return false;
    if (!(*a.despl == *b.despl)) return false;
    if (!(*a.front == *b.front)) return false;
    if (!(*a.Mats == *b.Mats)) return false;

    if (!ignoreRegions) {
        if (!(*a.pecRegs == *b.pecRegs)) return false;
        if (!(*a.pmcRegs == *b.pmcRegs)) return false;
        if (!(*a.DielRegs == *b.DielRegs)) return false;
        if (!(*a.LossyThinSurfs == *b.LossyThinSurfs)) return false;
    }

    if (!(*a.plnSrc == *b.plnSrc)) return false;
    if (!(*a.nodSrc == *b.nodSrc)) return false;
    if (!(*a.oldSONDA == *b.oldSONDA)) return false;
    if (!(*a.Sonda == *b.Sonda)) return false;
    if (!(*a.BloquePrb == *b.BloquePrb)) return false;
    if (!(*a.VolPrb == *b.VolPrb)) return false;
    if (!(*a.tWires == *b.tWires)) return false;
    if (!(*a.tSlots == *b.tSlots)) return false;
    return true;
}

} // namespace NFDETypes_extension_m

#endif // NFDE_TYPES_EXTENSION_M_H
