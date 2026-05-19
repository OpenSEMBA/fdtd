// This file corresponds to NFDETypes_extension_m.f90
// Includes for types defined in NFDETypes_m should be added here.
// #include "NFDETypes.h"
#include <optional>
#include <memory>
#include <string>
#include <vector>

// Forward declarations or includes for types used in this module
// Assuming NFDETypes.h contains definitions for Parseador_t, MatrizMedios_t, etc.

namespace NFDETypes_extension_m {

    // Placeholder for constants if they are not defined in NFDETypes_m
    // In a real scenario, these would come from the included header or be defined here.
    // extern const double EPSILON_VACUUM;
    // extern const double MU_VACUUM;
    // extern const double SIGMA_PEC;
    // extern const double SIGMA_PMC;

    // Helper function to check equality of two vectors
    template <typename T>
    bool vectors_equal(const std::vector<T>& a, const std::vector<T>& b) {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] != b[i]) return false;
        }
        return true;
    }

    // Helper function to check equality of two vectors of pointers/objects with operator==
    template <typename T>
    bool vector_objects_equal(const std::vector<T>& a, const std::vector<T>& b) {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); ++i) {
            if (!(a[i] == b[i])) return false;
        }
        return true;
    }

    // Helper for associated logic simulation
    // In C++, std::vector handles memory. "Associated" roughly maps to "not empty" or explicit pointer checks.
    // For std::vector, we assume if size > 0, it's "associated".
    // The Fortran logic:
    // if (associated(a%Lins) .eqv. associated(b%Lins)) then
    //    if (associated(a%Lins)) then ... compare ...
    // else if (associated(a%Lins)) then
    //    if (size(a%Lins) /= 0) return .false.
    // else
    //    if (size(b%Lins) /= 0) return .false.
    // end if
    // This simplifies to: if sizes differ, return false. If both empty, true. If both non-empty, compare.
    // However, the Fortran code explicitly checks `associated` which implies pointer validity.
    // In C++ std::vector, we just check size. If the Fortran code relies on null pointers for "not associated",
    // and C++ uses vectors, we map "associated" to "size() > 0".

    inline bool compare_associated_logic(const std::vector<int>& a, const std::vector<int>& b) {
        bool a_assoc = !a.empty();
        bool b_assoc = !b.empty();
        
        if (a_assoc == b_assoc) {
            if (a_assoc) {
                return vectors_equal(a, b);
            }
            return true;
        } else {
            if (a_assoc) return false;
            else return false;
        }
    }
    
    // Generic version for any type with operator==
    template <typename T>
    bool compare_associated_logic_vec(const std::vector<T>& a, const std::vector<T>& b) {
        bool a_assoc = !a.empty();
        bool b_assoc = !b.empty();
        
        if (a_assoc == b_assoc) {
            if (a_assoc) {
                return vector_objects_equal(a, b);
            }
            return true;
        } else {
            if (a_assoc) return false;
            else return false;
        }
    }

    // Note: The following functions are marked 'elemental' in Fortran, meaning they can be applied element-wise.
    // In C++, we just implement them as regular functions.

    bool Parseador_eq(const struct Parseador_t& a, const struct Parseador_t& b);
    bool NFDEGeneral_eq(const struct NFDEGeneral_t& a, const struct NFDEGeneral_t& b);
    bool desplazamiento_eq(const struct desplazamiento_t& a, const struct desplazamiento_t& b); // Assuming this type exists
    bool MatrizMedios_eq(const struct MatrizMedios_t& a, const struct MatrizMedios_t& b);
    bool Material_eq(const struct Material_t& a, const struct Material_t& b);
    bool Materials_eq(const struct Materials_t& a, const struct Materials_t& b);
    bool pecregions_eq(const struct PECRegions_t& a, const struct PECRegions_t& b);
    bool dielectric_eq(const struct dielectric_t& a, const struct dielectric_t& b);
    bool dielectricregions_eq(const struct DielectricRegions_t& a, const struct DielectricRegions_t& b);
    bool freqdepenmaterial_eq(const struct FreqDepenMaterial_t& a, const struct FreqDepenMaterial_t& b);
    bool freqdepenmaterials_eq(const struct FreqDepenMaterials_t& a, const struct FreqDepenMaterials_t& b);
    bool anisotropicbody_eq(const struct ANISOTROPICbody_t& a, const struct ANISOTROPICbody_t& b);
    bool anisotropicelements_eq(const struct ANISOTROPICelements_t& a, const struct ANISOTROPICelements_t& b);
    bool LossyThinSurface_eq(const struct LossyThinSurface_t& a, const struct LossyThinSurface_t& b);
    bool LossyThinSurfaces_eq(const struct LossyThinSurfaces_t& a, const struct LossyThinSurfaces_t& b);
    bool frontera_eq(const struct frontera_t& a, const struct frontera_t& b); // Assuming this type exists
    bool fronteraPML_eq(const struct fronteraPML_t& a, const struct fronteraPML_t& b); // Assuming this type exists
    bool box_eq(const struct Box_t& a, const struct Box_t& b);
    bool boxes_eq(const struct Boxes_t& a, const struct Boxes_t& b);
    bool planewaves_eq(const struct planewaves_t& a, const struct planewaves_t& b); // Assuming this type exists
    bool planewave_eq(const struct PlaneWave_t& a, const struct PlaneWave_t& b);
    bool curr_field_src_eq(const struct curr_field_src_t& a, const struct curr_field_src_t& b); // Assuming this type exists
    bool nodsource_eq(const struct nodsource_t& a, const struct nodsource_t& b); // Assuming this type exists
    bool coords_eq(const struct coords_t& a, const struct coords_t& b); // Assuming this type exists
    bool coords_scaled_eq(const struct coords_scaled_t& a, const struct coords_scaled_t& b); // Assuming this type exists
    bool abstractSonda_eq(const struct abstractSonda_t& a, const struct abstractSonda_t& b); // Assuming this type exists
    bool sonda_eq(const struct sonda_t& a, const struct sonda_t& b); // Assuming this type exists
    bool sondas_eq(const struct sondas_t& a, const struct sondas_t& b); // Assuming this type exists
    bool massonda_eq(const struct massonda_t& a, const struct massonda_t& b); // Assuming this type exists
    bool massondas_eq(const struct massondas_t& a, const struct massondas_t& b); // Assuming this type exists
    bool bloqueprobe_eq(const struct bloqueprobe_t& a, const struct bloqueprobe_t& b); // Assuming this type exists
    bool bloqueprobes_eq(const struct bloqueprobes_t& a, const struct bloqueprobes_t& b); // Assuming this type exists
    bool volprobe_eq(const struct volprobe_t& a, const struct volprobe_t& b); // Assuming this type exists
    bool volprobes_eq(const struct volprobes_t& a, const struct volprobes_t& b); // Assuming this type exists
    bool FarField_Sonda_eq(const struct FarField_Sonda_t& a, const struct FarField_Sonda_t& b); // Assuming this type exists
    bool Electric_Sonda_eq(const struct Electric_Sonda_t& a, const struct Electric_Sonda_t& b); // Assuming this type exists
    bool Magnetic_Sonda_eq(const struct Magnetic_Sonda_t& a, const struct Magnetic_Sonda_t& b); // Assuming this type exists
    bool NormalElectric_Sonda_eq(const struct NormalElectric_Sonda_t& a, const struct NormalElectric_Sonda_t& b); // Assuming this type exists
    bool NormalMagnetic_Sonda_eq(const struct NormalMagnetic_Sonda_t& a, const struct NormalMagnetic_Sonda_t& b); // Assuming this type exists
    bool SurfaceElectricCurrent_Sonda_eq(const struct SurfaceElectricCurrent_Sonda_t& a, const struct SurfaceElectricCurrent_Sonda_t& b); // Assuming this type exists
    bool SurfaceMagneticCurrent_Sonda_eq(const struct SurfaceMagneticCurrent_Sonda_t& a, const struct SurfaceMagneticCurrent_Sonda_t& b); // Assuming this type exists
    bool ThinWireComp_eq(const struct ThinWireComp_t& a, const struct ThinWireComp_t& b);
    bool ThinWire_eq(const struct ThinWire_t& a, const struct ThinWire_t& b);
    bool ThinWires_eq(const struct ThinWires_t& a, const struct ThinWires_t& b);
    bool SlantedWireComp_eq(const struct SlantedWireComp_t& a, const struct SlantedWireComp_t& b);
    bool SlantedWire_eq(const struct SlantedWire_t& a, const struct SlantedWire_t& b);
    bool SlantedWires_eq(const struct SlantedWiresInfo_t& a, const struct SlantedWiresInfo_t& b);
    bool ThinSlotComp_eq(const struct ThinSlotComp_t& a, const struct ThinSlotComp_t& b);
    bool ThinSlot_eq(const struct ThinSlot_t& a, const struct ThinSlot_t& b);
    bool ThinSlots_eq(const struct ThinSlots_t& a, const struct ThinSlots_t& b);

    void initializeProblemDescription(struct Parseador_t& pD) {
        // Assuming Parseador_t has public members corresponding to the fields
        // allocate(pD%general) -> Assuming general is a pointer or unique_ptr in C++ struct
        // Since we don't have the definition of Parseador_t, we assume standard C++ equivalents.
        // If Parseador_t is defined in NFDETypes.h, we assume it has smart pointers or raw pointers.
        // For the sake of translation, we assume the C++ struct has been adapted to use std::unique_ptr or similar for allocated members.
        
        // Note: Without the actual struct definitions, this is a best-effort translation of the logic.
        // The user must ensure Parseador_t and its members are correctly defined in C++.
        
        // pD.general.reset(new NFDEGeneral_t()); // Example
        // pD.matriz.reset(new MatrizMedios_t());
        // pD.despl.reset(new desplazamiento_t());
        // pD.front.reset(new frontera_t());
        
        // pD.mats.n_Mats = 3;
        // pD.mats.n_Mats_max = 3;
        // pD.mats.mats.resize(3);
        
        // pD.mats.mats[0].id = 1;
        // pD.mats.mats[0].eps = EPSILON_VACUUM;
        // pD.mats.mats[0].mu = MU_VACUUM;
        // pD.mats.mats[0].sigma = 0.0;
        // pD.mats.mats[0].sigmam = 0.0;
        
        // pD.mats.mats[1].id = 2;
        // pD.mats.mats[1].eps = EPSILON_VACUUM;
        // pD.mats.mats[1].mu = MU_VACUUM;
        // pD.mats.mats[1].sigma = SIGMA_PEC;
        // pD.mats.mats[1].sigmam = 0.0;
        
        // pD.mats.mats[2].id = 3;
        // pD.mats.mats[2].eps = EPSILON_VACUUM;
        // pD.mats.mats[2].mu = MU_VACUUM;
        // pD.mats.mats[2].sigma = 0.0;
        // pD.mats.mats[2].sigmam = SIGMA_PMC;
        
        // pD.pecRegs.lins.resize(0);
        // pD.pecRegs.surfs.resize(0);
        // pD.pecRegs.vols.resize(0);
        
        // pD.pmcRegs.lins.resize(0);
        // pD.pmcRegs.surfs.resize(0);
        // pD.pmcRegs.vols.resize(0);
        
        // pD.DielRegs.lins.resize(0);
        // pD.DielRegs.surfs.resize(0);
        // pD.DielRegs.vols.resize(0);
        
        // pD.LossyThinSurfs.cs.resize(0);
        
        // pD.frqDepMats = ... // Depends on type
        // pD.aniMats = ...
        
        // pD.plnSrc.collection.resize(0);
        // pD.nodSrc.NodalSource.resize(0);
        
        // pD.Sonda.len_cor_max = 0;
        // pD.Sonda.length = 0;
        // pD.Sonda.length_max = 0;
        // pD.Sonda.collection.resize(0);
        
        // pD.oldSONDA.probes.resize(0);
        // pD.oldSONDA.n_probes = 0;
        // pD.oldSONDA.n_probes_max = 0;
        
        // pD.BloquePrb.bp.resize(0);
        
        // pD.VolPrb.collection.resize(0);
        
        // pD.tWires.tw.resize(0);
        // pD.tWires.n_tw = 0;
        // pD.tWires.n_tw = 0; // Duplicate assignment in Fortran
        
        // pD.sWires = ...
        // pD.tSlots.tg.resize(0);
        
        #ifdef CompileWithMTLN
            // pD.mtln.cables.resize(0);
            // pD.mtln.probes.resize(0);
            // pD.mtln.networks.resize(0);
            // pD.mtln.connectors.resize(0);
        #endif
        
        // pD.conformalRegs.volumes.resize(0);
        // pD.conformalRegs.surfaces.resize(0);
    }

    bool Parseador_eq(const struct Parseador_t& a, const struct Parseador_t& b) {
        if (!(a.switches == b.switches)) return false;
        if (!(a.general == b.general)) return false;
        if (!(a.matriz == b.matriz)) return false;
        if (!(a.despl == b.despl)) return false;
        if (!(a.front == b.front)) return false;
        if (!(a.Mats == b.Mats)) return false;
        if (!(a.pecRegs == b.pecRegs)) return false;
        if (!(a.pmcRegs == b.pmcRegs)) return false;
        if (!(a.DielRegs == b.DielRegs)) return false;
        if (!(a.LossyThinSurfs == b.LossyThinSurfs)) return false;
        if (!(a.frqDepMats == b.frqDepMats)) return false;
        if (!(a.aniMats == b.aniMats)) return false;
        if (!(a.boxSrc == b.boxSrc)) return false;
        if (!(a.plnSrc == b.plnSrc)) return false;
        if (!(a.nodSrc == b.nodSrc)) return false;
        if (!(a.oldSONDA == b.oldSONDA)) return false;
        if (!(a.Sonda == b.Sonda)) return false;
        if (!(a.BloquePrb == b.BloquePrb)) return false;
        if (!(a.VolPrb == b.VolPrb)) return false;
        if (!(a.tWires == b.tWires)) return false;
        if (!(a.sWires == b.sWires)) return false;
        if (!(a.tSlots == b.tSlots)) return false;
        return true;
    }

    bool MatrizMedios_eq(const struct MatrizMedios_t& a, const struct MatrizMedios_t& b) {
        if (a.totalX != b.totalX) return false;
        if (a.totalY != b.totalY) return false;
        if (a.totalZ != b.totalZ) return false;
        return true;
    }

    bool Material_eq(const struct Material_t& a, const struct Material_t& b) {
        if (a.eps != b.eps) return false;
        if (a.mu != b.mu) return false;
        if (a.sigma != b.sigma) return false;
        if (a.sigmam != b.sigmam) return false;
        if (a.id != b.id) return false;
        return true;
    }

    bool Materials_eq(const struct Materials_t& a, const struct Materials_t& b) {
        if (a.n_Mats != b.n_Mats) return false;
        if (a.n_Mats_max != b.n_Mats_max) return false;
        if (!vectors_equal(a.Mats, b.Mats)) return false; // Assuming Mats is a vector of Material_t or similar
        return true;
    }

    bool pecregions_eq(const struct PECRegions_t& a, const struct PECRegions_t& b) {
        if (a.nLins != b.nLins) return false;
        if (a.nLins_max != b.nLins_max) return false;
        if (a.nSurfs != b.nSurfs) return false;
        if (a.nSurfs_max != b.nSurfs_max) return false;
        if (a.nVols != b.nVols) return false;
        if (a.nVols_max != b.nVols_max) return false;
        
        if (!compare_associated_logic_vec(a.Lins, b.Lins)) return false;
        if (!compare_associated_logic_vec(a.Surfs, b.Surfs)) return false;
        if (!compare_associated_logic_vec(a.Vols, b.Vols)) return false;
        
        return true;
    }

    bool dielectric_eq(const struct dielectric_t& a, const struct dielectric_t& b) {
        if (a.n_C1P != b.n_C1P) return false;
        if (a.n_C2P != b.n_C2P) return false;
        if (a.sigma != b.sigma) return false;
        if (a.eps != b.eps) return false;
        if (a.mu != b.mu) return false;
        if (a.sigmam != b.sigmam) return false;
        if (a.Rtime_on != b.Rtime_on) return false;
        if (a.Rtime_off != b.Rtime_off) return false;
        if (a.R != b.R) return false;
        if (a.L != b.L) return false;
        if (a.C != b.C) return false;
        if (a.R_devia != b.R_devia) return false;
        if (a.L_devia != b.L_devia) return false;
        if (a.C_devia != b.C_devia) return false;
        if (a.DiodB != b.DiodB) return false;
        if (a.DiodIsat != b.DiodIsat) return false;
        if (a.DiodOri != b.DiodOri) return false;
        if (a.orient != b.orient) return false;
        if (a.resistor != b.resistor) return false;
        if (a.inductor != b.inductor) return false;
        if (a.capacitor != b.capacitor) return false;
        if (a.diodo != b.diodo) return false;
        if (a.plain != b.plain) return false;
        if (a.PMLbody != b.PMLbody) return false;
        
        // Assuming C1P and C2P are vectors or pointers
        if (a.C1P.empty()) return false;
        if (b.C1P.empty()) return false;
        if (!vectors_equal(a.C1P, b.C1P)) return false;
        
        if (a.C2P.empty()) return false;
        if (b.C2P.empty()) return false;
        if (!vectors_equal(a.C2P, b.C2P)) return false;
        
        return true;
    }

    bool freqdepenmaterial_eq(const struct FreqDepenMaterial_t& a, const struct FreqDepenMaterial_t& b) {
        if (!vectors_equal(a.a11, b.a11)) return false;
        if (!vectors_equal(a.b11, b.b11)) return false;
        if (!vectors_equal(a.am11, b.am11)) return false;
        if (!vectors_equal(a.bm11, b.bm11)) return false;
        if (!vectors_equal(a.a12, b.a12)) return false;
        if (!vectors_equal(a.b12, b.b12)) return false;
        if (!vectors_equal(a.am12, b.am12)) return false;
        if (!vectors_equal(a.bm12, b.bm12)) return false;
        if (!vectors_equal(a.a13, b.a13)) return false;
        if (!vectors_equal(a.b13, b.b13)) return false;
        if (!vectors_equal(a.am13, b.am13)) return false;
        if (!vectors_equal(a.bm13, b.bm13)) return false;
        if (!vectors_equal(a.a22, b.a22)) return false;
        if (!vectors_equal(a.b22, b.b22)) return false;
        if (!vectors_equal(a.am22, b.am22)) return false;
        if (!vectors_equal(a.bm22, b.bm22)) return false;
        if (!vectors_equal(a.a23, b.a23)) return false;
        if (!vectors_equal(a.b23, b.b23)) return false;
        if (!vectors_equal(a.am23, b.am23)) return false;
        if (!vectors_equal(a.bm23, b.bm23)) return false;
        if (!vectors_equal(a.a33, b.a33)) return false;
        if (!vectors_equal(a.b33, b.b33)) return false;
        if (!vectors_equal(a.am33, b.am33)) return false;
        if (!vectors_equal(a.bm33, b.bm33)) return false;
        if (!vectors_equal(a.alpha, b.alpha)) return false;
        if (!vectors_equal(a.beta, b.beta)) return false;
        if (!vectors_equal(a.gamma, b.gamma)) return false;
        if (!vectors_equal(a.alpham, b.alpham)) return false;
        if (!vectors_equal(a.betam, b.betam)) return false;
        if (!vectors_equal(a.gammam, b.gammam)) return false;
        
        // coords_eq(a%c(1:a%n_c), b%c(1:b%n_c))
        // Assuming c is a vector of coords_t
        if (a.n_c != b.n_c) return false; // Lengths must match for element-wise comparison
        if (!vector_objects_equal(a.c, b.c)) return false;
        
        if (a.eps11 != b.eps11) return false;
        if (a.eps12 != b.eps12) return false;
        if (a.eps13 != b.eps13) return false;
        if (a.eps22 != b.eps22) return false;
        if (a.eps23 != b.eps23) return false;
        if (a.eps33 != b.eps33) return false;
        if (a.mu11 != b.mu11) return false;
        if (a.mu12 != b.mu12) return false;
        if (a.mu13 != b.mu13) return false;
        if (a.mu22 != b.mu22) return false;
        if (a.mu23 != b.mu23) return false;
        if (a.mu33 != b.mu33) return false;
        if (a.sigma11 != b.sigma11) return false;
        if (a.sigma12 != b.sigma12) return false;
        if (a.sigma13 != b.sigma13) return false;
        if (a.sigma22 != b.sigma22) return false;
        if (a.sigma23 != b.sigma23) return false;
        if (a.sigma33 != b.sigma33) return false;
        if (a.sigmam11 != b.sigmam11) return false;
        if (a.sigmam12 != b.sigmam12) return false;
        if (a.sigmam13 != b.sigmam13) return false;
        if (a.sigmam22 != b.sigmam22) return false;
        if (a.sigmam23 != b.sigmam23) return false;
        if (a.sigmam33 != b.sigmam33) return false;
        if (a.K11 != b.K11) return false;
        if (a.Km11 != b.Km11) return false;
        if (a.K12 != b.K12) return false;
        if (a.Km12 != b.Km12) return false;
        if (a.K13 != b.K13) return false;
        if (a.Km13 != b.Km13) return false;
        if (a.K22 != b.K22) return false;
        if (a.Km22 != b.Km22) return false;
        if (a.K23 != b.K23) return false;
        if (a.Km23 != b.Km23) return false;
        if (a.K33 != b.K33) return false;
        if (a.Km33 != b.Km33) return false;
        if (a.L != b.L) return false;
        if (a.Lm != b.Lm) return false;
        if (a.n_c != b.n_c) return false;
        if (a.files != b.files) return false;
        
        return true;
    }

    bool freqdepenmaterials_eq(const struct FreqDepenMaterials_t& a, const struct FreqDepenMaterials_t& b) {
        if (a.nVols != b.nVols) return false;
        if (a.nSurfs != b.nSurfs) return false;
        if (a.nLins != b.nLins) return false;
        if (a.nVols_max != b.nVols_max) return false;
        if (a.nSurfs_max != b.nSurfs_max) return false;
        if (a.nLins_max != b.nLins_max) return false;
        if (a.n_c_max != b.n_c_max) return false;
        if (!vector_objects_equal(a.Vols, b.Vols)) return false;
        if (!vector_objects_equal(a.Surfs, b.Surfs)) return false;
        if (!vector_objects_equal(a.Lins, b.Lins)) return false;
        return true;
    }

    bool anisotropicbody_eq(const struct ANISOTROPICbody_t& a, const struct ANISOTROPICbody_t& b) {
        if (!vectors_equal(a.c1P, b.c1P)) return false;
        if (!vectors_equal(a.c2P, b.c2P)) return false;
        if (!vectors_equal(a.sigma, b.sigma)) return false;
        if (!vectors_equal(a.eps, b.eps)) return false;
        if (!vectors_equal(a.mu, b.mu)) return false;
        if (!vectors_equal(a.sigmam, b.sigmam)) return false;
        if (a.n_C1P != b.n_C1P) return false;
        if (a.n_C2P != b.n_C2P) return false;
        return true;
    }

    bool anisotropicelements_eq(const struct ANISOTROPICelements_t& a, const struct ANISOTROPICelements_t& b) {
        if (a.nVols != b.nVols) return false;
        if (a.nSurfs != b.nSurfs) return false;
        if (a.nLins != b.nLins) return false;
        if (a.nVols_max != b.nVols_max) return false;
        if (a.nSurfs_max != b.nSurfs_max) return false;
        if (a.nLins_max != b.nLins_max) return false;
        if (a.n_C1P_max != b.n_C1P_max) return false;
        if (a.n_C2P_max != b.n_C2P_max) return false;
        if (!vector_objects_equal(a.Vols, b.Vols)) return false;
        if (!vector_objects_equal(a.Surfs, b.Surfs)) return false;
        if (!vector_objects_equal(a.Lins, b.Lins)) return false;
        return true;
    }

    bool dielectricregions_eq(const struct DielectricRegions_t& a, const struct DielectricRegions_t& b) {
        if (a.nVols != b.nVols) return false;
        if (a.nSurfs != b.nSurfs) return false;
        if (a.nLins != b.nLins) return false;
        if (a.nVols_max != b.nVols_max) return false;
        if (a.nSurfs_max != b.nSurfs_max) return false;
        if (a.nLins_max != b.nLins_max) return false;
        if (a.n_C1P_max != b.n_C1P_max) return false;
        if (a.n_C2P_max != b.n_C2P_max) return false;
        
        if (!compare_associated_logic_vec(a.Lins, b.Lins)) return false;
        if (!compare_associated_logic_vec(a.Surfs, b.Surfs)) return false;
        if (!compare_associated_logic_vec(a.Vols, b.Vols)) return false;
        
        return true;
    }

    bool LossyThinSurface_eq(const struct LossyThinSurface_t& a, const struct LossyThinSurface_t& b) {
        if (a.nc != b.nc) return false;
        if (a.files != b.files) return false;
        if (a.numcapas != b.numcapas) return false;
        if (!vectors_equal(a.c, b.c)) return false;
        if (!vectors_equal(a.sigma, b.sigma)) return false;
        if (!vectors_equal(a.eps, b.eps)) return false;
        if (!vectors_equal(a.mu, b.mu)) return false;
        if (!vectors_equal(a.sigmam, b.sigmam)) return false;
        if (!vectors_equal(a.thk, b.thk)) return false;
        if (!vectors_equal(a.sigma_devia, b.sigma_devia)) return false;
        if (!vectors_equal(a.eps_devia, b.eps_devia)) return false;
        if (!vectors_equal(a.mu_devia, b.mu_devia)) return false;
        if (!vectors_equal(a.sigmam_devia, b.sigmam_devia)) return false;
        if (!vectors_equal(a.thk_devia, b.thk_devia)) return false;
        return true;
    }

    bool LossyThinSurfaces_eq(const struct LossyThinSurfaces_t& a, const struct LossyThinSurfaces_t& b) {
        if (!vector_objects_equal(a.cs, b.cs)) return false;
        if (a.length != b.length) return false;
        if (a.length_max != b.length_max) return false;
        if (a.nC_max != b.nC_max) return false;
        return true;
    }

    bool ThinWireComp_eq(const struct ThinWireComp_t& a, const struct ThinWireComp_t& b) {
        if (a.srctype != b.srctype) return false;
        if (a.srcfile != b.srcfile) return false;
        if (a.i != b.i) return false;
        if (a.j != b.j) return false;
        if (a.K != b.K) return false;
        if (a.nd != b.nd) return false;
        if (a.d != b.d) return false;
        if (a.m != b.m) return false;
        if (a.tag != b.tag) return false;
        return true;
    }

    bool ThinWire_eq(const struct ThinWire_t& a, const struct ThinWire_t& b) {
        if (a.rad != b.rad) return false;
        if (a.disp != b.disp) return false;
        if (a.dispfile != b.dispfile) return false;
        if (a.res != b.res) return false;
        if (a.ind != b.ind) return false;
        if (a.cap != b.cap) return false;
        if (a.P_res != b.P_res) return false;
        if (a.P_ind != b.P_ind) return false;
        if (a.P_cap != b.P_cap) return false;
        if (a.dispfile_LeftEnd != b.dispfile_LeftEnd) return false;
        if (a.R_LeftEnd != b.R_LeftEnd) return false;
        if (a.L_LeftEnd != b.L_LeftEnd) return false;
        if (a.C_LeftEnd != b.C_LeftEnd) return false;
        if (a.dispfile_RightEnd != b.dispfile_RightEnd) return false;
        if (a.R_RightEnd != b.R_RightEnd) return false;
        if (a.L_RightEnd != b.L_RightEnd) return false;
        if (a.C_RightEnd != b.C_RightEnd) return false;
        if (a.LeftEnd != b.LeftEnd) return false;
        if (a.RightEnd != b.RightEnd) return false;
        if (a.tl != b.tl) return false;
        if (a.tr != b.tr) return false;
        if (!vector_objects_equal(a.twc, b.twc)) return false;
        if (a.n_twc != b.n_twc) return false;
        if (a.n_twc_max != b.n_twc_max) return false;
        return true;
    }

    bool ThinWires_eq(const struct ThinWires_t& a, const struct ThinWires_t& b) {
        if (a.n_tw != b.n_tw) return false;
        if (a.n_tw_max != b.n_tw_max) return false;
        if (!compare_associated_logic_vec(a.tw, b.tw)) return false;
        return true;
    }

    bool SlantedWireComp_eq(const struct SlantedWireComp_t& a, const struct SlantedWireComp_t& b) {
        if (a.srctype != b.srctype) return false;
        if (a.srcfile != b.srcfile) return false;
        if (a.x != b.x) return false;
        if (a.y != b.y) return false;
        if (a.z != b.z) return false;
        if (a.nd != b.nd) return false;
        if (a.m != b.m) return false;
        if (a.tag != b.tag) return false;
        return true;
    }

    bool SlantedWire_eq(const struct SlantedWire_t& a, const struct SlantedWire_t& b) {
        if (a.rad != b.rad) return false;
        if (a.disp != b.disp) return false;
        if (a.dispfile != b.dispfile) return false;
        if (a.res != b.res) return false;
        if (a.ind != b.ind) return false;
        if (a.cap != b.cap) return false;
        if (a.P_res != b.P_res) return false;
        if (a.P_ind != b.P_ind) return false;
        if (a.P_cap != b.P_cap) return false;
        if (a.dispfile_LeftEnd != b.dispfile_LeftEnd) return false;
        if (a.R_LeftEnd != b.R_LeftEnd) return false;
        if (a.L_LeftEnd != b.L_LeftEnd) return false;
        if (a.C_LeftEnd != b.C_LeftEnd) return false;
        if (a.dispfile_RightEnd != b.dispfile_RightEnd) return false;
        if (a.R_RightEnd != b.R_RightEnd) return false;
        if (a.L_RightEnd != b.L_RightEnd) return false;
        if (a.C_RightEnd != b.C_RightEnd) return false;
        if (a.LeftEnd != b.LeftEnd) return false;
        if (a.RightEnd != b.RightEnd) return false;
        if (a.tl != b.tl) return false;
        if (a.tr != b.tr) return false;
        if (a.n_swc != b.n_swc) return false;
        if (a.n_swc_max != b.n_swc_max) return false;
        return true;
    }

    bool SlantedWires_eq(const struct SlantedWiresInfo_t& a, const struct SlantedWiresInfo_t& b) {
        if (a.sw.empty()) return false;
        if (b.sw.empty()) return false;
        if (!vector_objects_equal(a.sw, b.sw)) return false;
        if (a.n_sw != b.n_sw) return false;
        if (a.n_sw_max != b.n_sw_max) return false;
        return true;
    }

    bool ThinSlotComp_eq(const struct ThinSlotComp_t& a, const struct ThinSlotComp_t& b) {
        if (a.i != b.i) return false;
        if (a.j != b.j) return false;
        if (a.K != b.K) return false;
        if (a.node != b.node) return false;
        if (a.dir != b.dir) return false;
        if (a.Or != b.Or) return false;
        if (a.tag != b.tag) return false;
        return true;
    }

    bool ThinSlot_eq(const struct ThinSlot_t& a, const struct ThinSlot_t& b) {
        if (a.width != b.width) return false;
        if (a.n_tgc != b.n_tgc) return false;
        if (a.n_tgc_max != b.n_tgc_max) return false;
        if (!compare_associated_logic_vec(a.tgc, b.tgc)) return false;
        return true;
    }

    bool ThinSlots_eq(const struct ThinSlots_t& a, const struct ThinSlots_t& b) {
        if (a.n_tg != b.n_tg) return false;
        if (a.n_tg_max != b.n_tg_max) return false;
        if (!compare_associated_logic_vec(a.tg, b.tg)) return false;
        return true;
    }

    bool NFDEGeneral_eq(const struct NFDEGeneral_t& a, const struct NFDEGeneral_t& b) {
        if (a.dt != b.dt) return false;
        if (a.nmax != b.nmax) return false;
        return true;
    }

    bool box_eq(const struct Box_t& a, const struct Box_t& b) {
        if (a.nombre_fichero != b.nombre_fichero) return false;
        if (!vectors_equal(a.coor1, b.coor1)) return false;
        if (!vectors_equal(a.coor2, b.coor2)) return false;
        return true;
    }

    bool boxes_eq(const struct Boxes_t& a, const struct Boxes_t& b) {
        if (a.nVols != b.nVols) return false;
        if (a.nVols_max != b.nVols_max) return false;
        if (!vector_objects_equal(a.Vols, b.Vols)) return false;
        return true;
    }

    bool planewave_eq(const struct PlaneWave_t& a, const struct PlaneWave_t& b) {
        if (a.nombre_fichero != b.nombre_fichero) return false;
        if (a.atributo != b.atributo) return false;
        if (a.coor1.size() != b.coor1.size()) return false;
        for (size_t i = 0; i < a.coor1.size(); ++i) {
            if (a.coor1[i] != b.coor1[i]) return false;
        }
        if (a.coor2.size() != b.coor2.size()) return false;
        for (size_t i = 0; i < a.coor2.size(); ++i) {
            if (a.coor2[i] != b.coor2[i]) return false;
        }
        if (a.theta != b.theta) return false;
        if (a.phi != b.phi) return false;
        if (a.alpha != b.alpha) return false;
        if (a.beta != b.beta) return false;
        if (a.isRC != b.isRC) return false;
        return true;
    }

    // Other interface procedures are declared but not defined in the provided snippet.
    // They would follow similar patterns.

} // namespace NFDETypes_extension_m

// This file is a continuation of the translation.
// It contains equality comparison functions for various derived types.
// Assumes that the necessary headers (types.h, etc.) are included in the main translation unit.
// Assumes that std::vector, std::string, and std::optional (if used for pointers) are available.
// Assumes that 'associated' logic is handled by checking if a pointer is non-null or if a vector is non-empty,
// depending on how the Fortran allocatable/pointer arrays were translated.
// For this translation, we assume:
// - Fortran pointers/allocatables are translated to std::vector<T>* or std::optional<std::vector<T>>.
// - Given the pattern `associated(a%x)` and `size(a%x)`, it strongly suggests pointers to arrays.
// - We will assume the types have a member that is a pointer to a vector or a vector itself.
// - However, looking at `associated(a%collection)`, it's likely a pointer.
// - Let's assume the C++ types have members like `std::vector<T>* collection` or similar.
// - To be safe and generic, I will assume the C++ types have a method or property that indicates "associated" (non-null)
//   and a way to get the size and content.
// - Actually, a common pattern in Fortran-to-C++ translation for allocatables/pointers is to use `std::vector` directly
//   and check `!vec.empty()` for `associated`. But `associated` in Fortran specifically checks if the pointer is allocated/associated.
//   If the translation used `std::vector*`, then `associated` becomes `ptr != nullptr`.
//   Let's assume the C++ types have pointers to vectors for the "associated" checks.
//   Example: `type(PlaneWaves_t)` has `std::vector<int>* collection`.

// Helper lambda or function to check equality of two vectors pointed to by pointers
// This is a bit tricky without knowing the exact C++ type definitions.
// I will assume the C++ types have members that are pointers to std::vector<T>.
// And I will assume there is a way to compare the vectors.

// Since I cannot see the type definitions, I will write the code assuming:
// 1. The types have members that are pointers to std::vector<T>.
// 2. The pointers are named exactly as in Fortran.
// 3. The vectors support `==` operator or we use a loop.
// 4. `associated(a%x)` translates to `a.x != nullptr`.
// 5. `size(a%x)` translates to `a.x->size()`.
// 6. `a%x == b%x` translates to `*a.x == *b.x` or element-wise comparison.

// Note: The `#define CHECK_PTR` macro in Fortran is tricky in C++.
// I will inline the logic or use a helper function.

// Let's assume the following structure for a type with a pointer member `ptr`:
// struct Type_t {
//     std::vector<T>* ptr;
//     // ... other members
// };

// If the translation used `std::vector<T>` directly (no pointer), then `associated` is `!empty()`.
// But `associated` in Fortran is about pointer association.
// I will stick to the pointer assumption as it's more faithful to `associated`.

// However, to make the code compile, I need to know the exact types.
// Since I don't, I will use a generic approach assuming the members are pointers to vectors.

// If the previous chunks defined the types, this chunk just uses them.

// Let's write the functions.

namespace my_namespace { // Replace with actual namespace if known

// Helper to check if two pointers to vectors are equal
// This is a template helper.
template <typename T>
bool vectors_equal(const std::vector<T>* a, const std::vector<T>* b) {
    if (!a && !b) return true;
    if (!a || !b) return false;
    if (a->size() != b->size()) return false;
    for (size_t i = 0; i < a->size(); ++i) {
        if ((*a)[i] != (*b)[i]) return false;
    }
    return true;
}

// Specialization for std::string if needed, but vector<string> == works.
// Specialization for double/float if needed.

// For logical arrays, we assume they are stored as std::vector<bool> or std::vector<int>.
// The `==` operator for std::vector<bool> is tricky, so we might need a custom comparator.
// But let's assume `std::vector<bool>` works or they are stored as `std::vector<int>`.

// Let's assume the types are defined in a header that is included.

// Function: planewaves_eq
inline bool planewaves_eq(const PlaneWaves_t& a, const PlaneWaves_t& b) {
    if (a.nc != b.nc) return false;
    if (a.nC_max != b.nC_max) return false;
    
    // Check collection
    if (a.collection != b.collection) {
        // One is null, the other is not
        if (a.collection) {
            if (!a.collection->empty()) return false;
        } else {
            if (!b.collection->empty()) return false;
        }
    } else {
        // Both are null or both are non-null
        if (a.collection) {
            if (!vectors_equal(a.collection, b.collection)) return false;
        }
    }
    return true;
}

// Function: curr_field_src_eq
inline bool curr_field_src_eq(const Curr_Field_Src_t& a, const Curr_Field_Src_t& b) {
    if (a.n_C1P != b.n_C1P) return false;
    if (a.n_C2P != b.n_C2P) return false;
    
    // Check c1P
    if (a.c1P != b.c1P) {
        if (a.c1P) {
            if (!a.c1P->empty()) return false;
        } else {
            if (!b.c1P->empty()) return false;
        }
    } else {
        if (a.c1P) {
            if (!vectors_equal(a.c1P, b.c1P)) return false;
        }
    }
    
    // Check c2P
    if (a.c2P != b.c2P) {
        if (a.c2P) {
            if (!a.c2P->empty()) return false;
        } else {
            if (!b.c2P->empty()) return false;
        }
    } else {
        if (a.c2P) {
            if (!vectors_equal(a.c2P, b.c2P)) return false;
        }
    }
    
    if (a.nombre != b.nombre) return false;
    if (a.isElec != b.isElec) return false;
    if (a.isHard != b.isHard) return false;
    if (a.isInitialValue != b.isInitialValue) return false;
    
    return true;
}

// Function: nodsource_eq
inline bool nodsource_eq(const NodSource_t& a, const NodSource_t& b) {
    if (a.n_nodSrc != b.n_nodSrc) return false;
    if (a.n_nodSrc_max != b.n_nodSrc_max) return false;
    if (a.n_C1P_max != b.n_C1P_max) return false;
    if (a.n_C2P_max != b.n_C2P_max) return false;
    
    if (a.NodalSource != b.NodalSource) {
        if (a.NodalSource) {
            if (!a.NodalSource->empty()) return false;
        } else {
            if (!b.NodalSource->empty()) return false;
        }
    } else {
        if (a.NodalSource) {
            if (!vectors_equal(a.NodalSource, b.NodalSource)) return false;
        }
    }
    
    return true;
}

// Function: fronteraPML_eq
inline bool fronteraPML_eq(const FronteraPML_t& a, const FronteraPML_t& b) {
    if (a.orden != b.orden) return false;
    if (a.refl != b.refl) return false;
    if (a.numCapas != b.numCapas) return false;
    return true;
}

// Function: frontera_eq
inline bool frontera_eq(const Frontera_t& a, const Frontera_t& b) {
    // Assuming tipoFrontera is an array/vector
    if (a.tipoFrontera.size() != b.tipoFrontera.size()) return false;
    for (size_t i = 0; i < a.tipoFrontera.size(); ++i) {
        if (a.tipoFrontera[i] != b.tipoFrontera[i]) return false;
    }
    
    // Assuming propiedadesPML is an array/vector
    if (a.propiedadesPML.size() != b.propiedadesPML.size()) return false;
    for (size_t i = 0; i < a.propiedadesPML.size(); ++i) {
        if (a.propiedadesPML[i] != b.propiedadesPML[i]) return false;
    }
    
    return true;
}

// Function: desplazamiento_eq
inline bool desplazamiento_eq(const Desplazamiento_t& a, const Desplazamiento_t& b) {
    if (!a.desX || !a.desY || !a.desZ) return false;
    if (!b.desX || !b.desY || !b.desZ) return false;
    
    if (a.desX->size() != b.desX->size()) return false;
    if (a.desY->size() != b.desY->size()) return false;
    if (a.desZ->size() != b.desZ->size()) return false;
    
    for (size_t i = 0; i < a.desX->size(); ++i) {
        if ((*a.desX)[i] != (*b.desX)[i]) return false;
    }
    for (size_t i = 0; i < a.desY->size(); ++i) {
        if ((*a.desY)[i] != (*b.desY)[i]) return false;
    }
    for (size_t i = 0; i < a.desZ->size(); ++i) {
        if ((*a.desZ)[i] != (*b.desZ)[i]) return false;
    }
    
    if (a.nX != b.nX) return false;
    if (a.nY != b.nY) return false;
    if (a.nZ != b.nZ) return false;
    
    if (a.mx1 != b.mx1) return false;
    if (a.my1 != b.my1) return false;
    if (a.mz1 != b.mz1) return false;
    
    if (a.mx2 != b.mx2) return false;
    if (a.my2 != b.my2) return false;
    if (a.mz2 != b.mz2) return false;
    
    if (a.originx != b.originx) return false;
    if (a.originy != b.originy) return false;
    if (a.originz != b.originz) return false;
    
    return true;
}

// Function: coords_eq
inline bool coords_eq(const coords_t& a, const coords_t& b) {
    if (a.xi != b.xi) return false;
    if (a.xe != b.xe) return false;
    if (a.yi != b.yi) return false;
    if (a.ye != b.ye) return false;
    if (a.zi != b.zi) return false;
    if (a.ze != b.ze) return false;
    if (a.xtrancos != b.xtrancos) return false;
    if (a.ytrancos != b.ytrancos) return false;
    if (a.ztrancos != b.ztrancos) return false;
    if (a.or != b.or) return false;
    if (a.tag != b.tag) return false;
    return true;
}

// Function: coords_scaled_eq
inline bool coords_scaled_eq(const coords_scaled_t& a, const coords_scaled_t& b) {
    if (a.Xi != b.Xi) return false;
    if (a.Xe != b.Xe) return false;
    if (a.Yi != b.Yi) return false;
    if (a.Ye != b.Ye) return false;
    if (a.Zi != b.Zi) return false;
    if (a.Ze != b.Ze) return false;
    if (a.xc != b.xc) return false;
    if (a.yc != b.yc) return false;
    if (a.zc != b.zc) return false;
    if (a.Or != b.Or) return false;
    if (a.tag != b.tag) return false;
    return true;
}

// Function: FarField_Sonda_eq
inline bool FarField_Sonda_eq(const FarField_Sonda_t& a, const FarField_Sonda_t& b) {
    if (a.probe != b.probe) return false;
    return true;
}

// Function: Electric_Sonda_eq
inline bool Electric_Sonda_eq(const Electric_Sonda_t& a, const Electric_Sonda_t& b) {
    if (a.probe != b.probe) return false;
    return true;
}

// Function: Magnetic_Sonda_eq
inline bool Magnetic_Sonda_eq(const Magnetic_Sonda_t& a, const Magnetic_Sonda_t& b) {
    if (a.probe != b.probe) return false;
    return true;
}

// Function: NormalElectric_Sonda_eq
inline bool NormalElectric_Sonda_eq(const NormalElectric_Sonda_t& a, const NormalElectric_Sonda_t& b) {
    if (a.probe != b.probe) return false;
    if (a.n_nml != b.n_nml) return false;
    if (a.n_nml_max != b.n_nml_max) return false;
    
    if (a.nml != b.nml) {
        if (a.nml) {
            if (!a.nml->empty()) return false;
        } else {
            if (!b.nml->empty()) return false;
        }
    } else {
        if (a.nml) {
            if (!vectors_equal(a.nml, b.nml)) return false;
        }
    }
    
    return true;
}

// Function: NormalMagnetic_Sonda_eq
inline bool NormalMagnetic_Sonda_eq(const NormalMagnetic_Sonda_t& a, const NormalMagnetic_Sonda_t& b) {
    if (a.probe != b.probe) return false;
    if (a.n_nml != b.n_nml) return false;
    if (a.n_nml_max != b.n_nml_max) return false;
    
    if (a.nml != b.nml) {
        if (a.nml) {
            if (!a.nml->empty()) return false;
        } else {
            if (!b.nml->empty()) return false;
        }
    } else {
        if (a.nml) {
            if (!vectors_equal(a.nml, b.nml)) return false;
        }
    }
    
    return true;
}

// Function: SurfaceElectricCurrent_Sonda_eq
inline bool SurfaceElectricCurrent_Sonda_eq(const SurfaceElectricCurrent_Sonda_t& a, const SurfaceElectricCurrent_Sonda_t& b) {
    if (a.probe != b.probe) return false;
    if (a.n_nml != b.n_nml) return false;
    if (a.n_nml_max != b.n_nml_max) return false;
    
    if (a.nml != b.nml) {
        if (a.nml) {
            if (!a.nml->empty()) return false;
        } else {
            if (!b.nml->empty()) return false;
        }
    } else {
        if (a.nml) {
            if (!vectors_equal(a.nml, b.nml)) return false;
        }
    }
    
    return true;
}

// Function: SurfaceMagneticCurrent_Sonda_eq
inline bool SurfaceMagneticCurrent_Sonda_eq(const SurfaceMagneticCurrent_Sonda_t& a, const SurfaceMagneticCurrent_Sonda_t& b) {
    if (a.probe != b.probe) return false;
    if (a.n_nml != b.n_nml) return false;
    if (a.n_nml_max != b.n_nml_max) return false;
    
    if (a.nml != b.nml) {
        if (a.nml) {
            if (!a.nml->empty()) return false;
        } else {
            if (!b.nml->empty()) return false;
        }
    } else {
        if (a.nml) {
            if (!vectors_equal(a.nml, b.nml)) return false;
        }
    }
    
    return true;
}

// Function: abstractSonda_eq
inline bool abstractSonda_eq(const abstractSonda_t& a, const abstractSonda_t& b) {
    if (a.n_FarField != b.n_FarField) return false;
    if (a.n_Electric != b.n_Electric) return false;
    if (a.n_Magnetic != b.n_Magnetic) return false;
    if (a.n_NormalElectric != b.n_NormalElectric) return false;
    if (a.n_NormalMagnetic != b.n_NormalMagnetic) return false;
    if (a.n_SurfaceElectricCurrent != b.n_SurfaceElectricCurrent) return false;
    if (a.n_SurfaceMagneticCurrent != b.n_SurfaceMagneticCurrent) return false;
    if (a.n_FarField_max != b.n_FarField_max) return false;
    if (a.n_Electric_max != b.n_Electric_max) return false;
    if (a.n_Magnetic_max != b.n_Magnetic_max) return false;
    if (a.n_NormalElectric_max != b.n_NormalElectric_max) return false;
    if (a.n_NormalMagnetic_max != b.n_NormalMagnetic_max) return false;
    if (a.n_SurfaceElectricCurrent_max != b.n_SurfaceElectricCurrent_max) return false;
    if (a.n_SurfaceMagneticCurrent_max != b.n_SurfaceMagneticCurrent_max) return false;

    // Helper lambda for checking pointer members
    auto check_ptr = [](const auto* a_ptr, const auto* b_ptr) -> bool {
        if (a_ptr != b_ptr) {
            if (a_ptr) {
                if (!a_ptr->empty()) return false;
            } else {
                if (!b_ptr->empty()) return false;
            }
        } else {
            if (a_ptr) {
                if (!vectors_equal(a_ptr, b_ptr)) return false;
            }
        }
        return true;
    };

    if (!check_ptr(a.FarField, b.FarField)) return false;
    if (!check_ptr(a.Electric, b.Electric)) return false;
    if (!check_ptr(a.Magnetic, b.Magnetic)) return false;
    if (!check_ptr(a.NormalElectric, b.NormalElectric)) return false;
    if (!check_ptr(a.NormalMagnetic, b.NormalMagnetic)) return false;
    if (!check_ptr(a.SurfaceElectricCurrent, b.SurfaceElectricCurrent)) return false;
    if (!check_ptr(a.SurfaceMagneticCurrent, b.SurfaceMagneticCurrent)) return false;

    return true;
}

// Function: sondas_eq
inline bool sondas_eq(const Sondas_t& a, const Sondas_t& b) {
    if (a.n_probes != b.n_probes) return false;
    if (a.n_probes_max != b.n_probes_max) return false;
    
    if (a.probes != b.probes) {
        if (a.probes) {
            if (!a.probes->empty()) return false;
        } else {
            if (!b.probes->empty()) return false;
        }
    } else {
        if (a.probes) {
            if (!vectors_equal(a.probes, b.probes)) return false;
        }
    }
    
    return true;
}

// Function: sonda_eq
inline bool sonda_eq(const Sonda_t& a, const Sonda_t& b) {
    if (a.grname != b.grname) return false;
    
    // Check i
    if (a.i != b.i) {
        if (a.i) {
            if (!a.i->empty()) return false;
        } else {
            if (!b.i->empty()) return false;
        }
    } else {
        if (a.i) {
            if (!vectors_equal(a.i, b.i)) return false;
        }
    }
    
    // Check j
    if (a.j != b.j) {
        if (a.j) {
            if (!a.j->empty()) return false;
        } else {
            if (!b.j->empty()) return false;
        }
    } else {
        if (a.j) {
            if (!vectors_equal(a.j, b.j)) return false;
        }
    }
    
    // Check k
    if (a.k != b.k) {
        if (a.k) {
            if (!a.k->empty()) return false;
        } else {
            if (!b.k->empty()) return false;
        }
    } else {
        if (a.k) {
            if (!vectors_equal(a.k, b.k)) return false;
        }
    }
    
    // Check node
    if (a.node != b.node) {
        if (a.node) {
            if (!a.node->empty()) return false;
        } else {
            if (!b.node->empty()) return false;
        }
    } else {
        if (a.node) {
            if (!vectors_equal(a.node, b.node)) return false;
        }
    }
    
    if (a.n_cord != b.n_cord) return false;
    if (a.n_cord_max != b.n_cord_max) return false;
    if (a.tstart != b.tstart) return false;
    if (a.tstop != b.tstop) return false;
    if (a.tstep != b.tstep) return false;
    if (a.outputrequest != b.outputrequest) return false;
    if (a.fstart != b.fstart) return false;
    if (a.fstop != b.fstop) return false;
    if (a.fstep != b.fstep) return false;
    if (a.phistart != b.phistart) return false;
    if (a.phistop != b.phistop) return false;
    if (a.phistep != b.phistep) return false;
    if (a.thetastart != b.thetastart) return false;
    if (a.thetastop != b.thetastop) return false;
    if (a.thetastep != b.thetastep) return false;
    if (a.FileNormalize != b.FileNormalize) return false;
    
    return true;
}

// Function: masSonda_eq
inline bool masSonda_eq(const MasSonda_t& a, const MasSonda_t& b) {
    if (a.filename != b.filename) return false;
    if (a.type1 != b.type1) return false;
    if (a.type2 != b.type2) return false;
    if (a.outputrequest != b.outputrequest) return false;
    if (a.len_cor != b.len_cor) return false;
    if (a.tstart != b.tstart) return false;
    if (a.tstop != b.tstop) return false;
    if (a.tstep != b.tstep) return false;
    if (a.fstart != b.fstart) return false;
    if (a.fstop != b.fstop) return false;
    if (a.fstep != b.fstep) return false;
    
    if (a.cordinates != b.cordinates) {
        if (a.cordinates) {
            if (!a.cordinates->empty()) return false;
        } else {
            if (!b.cordinates->empty()) return false;
        }
    } else {
        if (a.cordinates) {
            if (!vectors_equal(a.cordinates, b.cordinates)) return false;
        }
    }
    
    return true;
}

// Function: MasSondas_eq
inline bool MasSondas_eq(const MasSondas_t& a, const MasSondas_t& b) {
    if (a.length != b.length) return false;
    if (a.length_max != b.length_max) return false;
    if (a.len_cor_max != b.len_cor_max) return false;
    
    if (a.collection != b.collection) {
        if (a.collection) {
            if (!a.collection->empty()) return false;
        } else {
            if (!b.collection->empty()) return false;
        }
    } else {
        if (a.collection) {
            // This is a vector of objects, not a vector of primitives.
            // We need to compare each element using the appropriate == operator.
            if (a.collection->size() != b.collection->size()) return false;
            for (size_t i = 0; i < a.collection->size(); ++i) {
                if (!(*a.collection)[i] == (*b.collection)[i]) return false;
            }
        }
    }
    
    return true;
}

// Function: bloqueprobe_eq
inline bool bloqueprobe_eq(const BloqueProbe_t& a, const BloqueProbe_t& b) {
    if (a.tstart != b.tstart) return false;
    if (a.tstop != b.tstop) return false;
    if (a.tstep != b.tstep) return false;
    if (a.fstart != b.fstart) return false;
    if (a.fstop != b.fstop) return false;
    if (a.fstep != b.fstep) return false;
    if (a.FileNormalize != b.FileNormalize) return false; // Assuming trim is handled by string comparison
    if (a.type2 != b.type2) return false;
    if (a.i1 != b.i1) return false;
    if (a.i2 != b.i2) return false;
    if (a.j1 != b.j1) return false;
    if (a.j2 != b.j2) return false;
    if (a.k1 != b.k1) return false;
    if (a.k2 != b.k2) return false;
    if (a.skip != b.skip) return false;
    if (a.nml != b.nml) return false;
    if (a.t != b.t) return false;
    if (a.outputrequest != b.outputrequest) return false;
    if (a.tag != b.tag) return false;
    
    return true;
}

// Function: bloqueprobes_eq
inline bool bloqueprobes_eq(const BloqueProbes_t& a, const BloqueProbes_t& b) {
    if (a.n_bp != b.n_bp) return false;
    if (a.n_bp_max != b.n_bp_max) return false;
    
    if (a.bp != b.bp) {
        if (a.bp) {
            if (!a.bp->empty()) return false;
        } else {
            if (!b.bp->empty()) return false;
        }
    } else {
        if (a.bp) {
            if (a.bp->size() != b.bp->size()) return false;
            for (size_t i = 0; i < a.bp->size(); ++i) {
                if (!(*a.bp)[i] == (*b.bp)[i]) return false;
            }
        }
    }
    
    return true;
}

// Function: volprobe_eq
inline bool volprobe_eq(const VolProbe_t& a, const VolProbe_t& b) {
    if (a.tstart != b.tstart) return false;
    if (a.tstop != b.tstop) return false;
    if (a.tstep != b.tstep) return false;
    if (a.outputrequest != b.outputrequest) return false;
    if (a.len_cor != b.len_cor) return false;
    if (a.fstart != b.fstart) return false;
    if (a.fstop != b.fstop) return false;
    if (a.fstep != b.fstep) return false;
    if (a.type2 != b.type2) return false;
    if (a.filename != b.filename) return false;
    
    if (a.cordinates != b.cordinates) {
        if (a.cordinates) {
            if (!a.cordinates->empty()) return false;
        } else {
            if (!b.cordinates->empty()) return false;
        }
    } else {
        if (a.cordinates) {
            if (a.cordinates->size() != b.cordinates->size()) return false;
            for (size_t i = 0; i < a.cordinates->size(); ++i) {
                if (!(*a.cordinates)[i] == (*b.cordinates)[i]) return false;
            }
        }
    }
    
    return true;
}

// Function: volprobes_eq
inline bool volprobes_eq(const VolProbes_t& a, const VolProbes_t& b) {
    if (a.length != b.length) return false;
    if (a.length_max != b.length_max) return false;
    if (a.len_cor_max != b.len_cor_max) return false;
    
    if (a.collection != b.collection) {
        if (a.collection) {
            if (!a.collection->empty()) return false;
        } else {
            if (!b.collection->empty()) return false;
        }
    } else {
        if (a.collection) {
            if (a.collection->size() != b.collection->size()) return false;
            for (size_t i = 0; i < a.collection->size(); ++i) {
                if (!(*a.collection)[i] == (*b.collection)[i]) return false;
            }
        }
    }
    
    return true;
}

} // namespace