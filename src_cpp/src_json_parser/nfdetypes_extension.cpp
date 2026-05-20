#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <memory>

// Forward declarations and includes for types defined in NFDETypes_m
// Assuming these types are defined in a header "NFDETypes.h"
// #include "NFDETypes.h"

// Constants assumed to be defined elsewhere or here
#ifndef EPSILON_VACUUM
#define EPSILON_VACUUM 8.854187817e-12
#endif
#ifndef MU_VACUUM
#define MU_VACUUM 1.2566370614e-6
#endif
#ifndef SIGMA_PEC
#define SIGMA_PEC 1.0e10 // Approximation for PEC
#endif
#ifndef SIGMA_PMC
#define SIGMA_PMC 1.0e10 // Approximation for PMC
#endif

// Placeholder types to make the code compile structurally. 
// In a real scenario, these would be fully implemented classes/structs 
// matching the Fortran derived types from NFDETypes_m.

struct Switches_t {
    // Placeholder fields
    int dummy = 0;
    bool operator==(const Switches_t& other) const { return dummy == other.dummy; }
    bool operator!=(const Switches_t& other) const { return !(*this == other); }
};

struct Parseador_t; // Forward declaration

struct NFDEGeneral_t {
    double dt = 0.0;
    int nmax = 0;
    bool operator==(const NFDEGeneral_t& other) const {
        return (dt == other.dt) && (nmax == other.nmax);
    }
    bool operator!=(const NFDEGeneral_t& other) const { return !(*this == other); }
};

struct MatrizMedios_t {
    int totalX = 0;
    int totalY = 0;
    int totalZ = 0;
    bool operator==(const MatrizMedios_t& other) const {
        return (totalX == other.totalX) && (totalY == other.totalY) && (totalZ == other.totalZ);
    }
    bool operator!=(const MatrizMedios_t& other) const { return !(*this == other); }
};

struct Material_t {
    double eps = 0.0;
    double mu = 0.0;
    double sigma = 0.0;
    double sigmam = 0.0;
    int id = 0;
    bool operator==(const Material_t& other) const {
        return (eps == other.eps) && (mu == other.mu) && (sigma == other.sigma) && 
               (sigmam == other.sigmam) && (id == other.id);
    }
    bool operator!=(const Material_t& other) const { return !(*this == other); }
};

struct Materials_t {
    int n_Mats = 0;
    int n_Mats_max = 0;
    std::vector<Material_t> mats;
    bool operator==(const Materials_t& other) const {
        if (n_Mats != other.n_Mats) return false;
        if (n_Mats_max != other.n_Mats_max) return false;
        if (mats.size() != other.mats.size()) return false;
        for (size_t i = 0; i < mats.size(); ++i) {
            if (mats[i] != other.mats[i]) return false;
        }
        return true;
    }
    bool operator!=(const Materials_t& other) const { return !(*this == other); }
};

struct PECRegions_t {
    int nLins = 0;
    int nLins_max = 0;
    int nSurfs = 0;
    int nSurfs_max = 0;
    int nVols = 0;
    int nVols_max = 0;
    std::vector<int> Lins; // Placeholder type
    std::vector<int> Surfs; // Placeholder type
    std::vector<int> Vols; // Placeholder type
    
    bool operator==(const PECRegions_t& other) const {
        if (nLins != other.nLins) return false;
        if (nLins_max != other.nLins_max) return false;
        if (nSurfs != other.nSurfs) return false;
        if (nSurfs_max != other.nSurfs_max) return false;
        if (nVols != other.nVols) return false;
        if (nVols_max != other.nVols_max) return false;
        
        // Check pointers/associations logic
        bool aLinsAssoc = !Lins.empty();
        bool bLinsAssoc = !other.Lins.empty();
        if (aLinsAssoc != bLinsAssoc) return false;
        if (aLinsAssoc) {
            if (Lins != other.Lins) return false;
        }

        bool aSurfsAssoc = !Surfs.empty();
        bool bSurfsAssoc = !other.Surfs.empty();
        if (aSurfsAssoc != bSurfsAssoc) return false;
        if (aSurfsAssoc) {
            if (Surfs != other.Surfs) return false;
        }

        bool aVolsAssoc = !Vols.empty();
        bool bVolsAssoc = !other.Vols.empty();
        if (aVolsAssoc != bVolsAssoc) return false;
        if (aVolsAssoc) {
            if (Vols != other.Vols) return false;
        }

        return true;
    }
    bool operator!=(const PECRegions_t& other) const { return !(*this == other); }
};

struct PMCRegions_t {
    int nLins = 0;
    int nLins_max = 0;
    int nSurfs = 0;
    int nSurfs_max = 0;
    int nVols = 0;
    int nVols_max = 0;
    std::vector<int> Lins;
    std::vector<int> Surfs;
    std::vector<int> Vols;

    bool operator==(const PMCRegions_t& other) const {
        if (nLins != other.nLins) return false;
        if (nLins_max != other.nLins_max) return false;
        if (nSurfs != other.nSurfs) return false;
        if (nSurfs_max != other.nSurfs_max) return false;
        if (nVols != other.nVols) return false;
        if (nVols_max != other.nVols_max) return false;
        
        bool aLinsAssoc = !Lins.empty();
        bool bLinsAssoc = !other.Lins.empty();
        if (aLinsAssoc != bLinsAssoc) return false;
        if (aLinsAssoc) {
            if (Lins != other.Lins) return false;
        }

        bool aSurfsAssoc = !Surfs.empty();
        bool bSurfsAssoc = !other.Surfs.empty();
        if (aSurfsAssoc != bSurfsAssoc) return false;
        if (aSurfsAssoc) {
            if (Surfs != other.Surfs) return false;
        }

        bool aVolsAssoc = !Vols.empty();
        bool bVolsAssoc = !other.Vols.empty();
        if (aVolsAssoc != bVolsAssoc) return false;
        if (aVolsAssoc) {
            if (Vols != other.Vols) return false;
        }

        return true;
    }
    bool operator!=(const PMCRegions_t& other) const { return !(*this == other); }
};

struct DielectricRegions_t {
    int nVols = 0;
    int nSurfs = 0;
    int nLins = 0;
    int nVols_max = 0;
    int nSurfs_max = 0;
    int nLins_max = 0;
    int n_C1P_max = 0;
    int n_C2P_max = 0;
    std::vector<int> Lins;
    std::vector<int> Surfs;
    std::vector<int> Vols;

    bool operator==(const DielectricRegions_t& other) const {
        if (nVols != other.nVols) return false;
        if (nSurfs != other.nSurfs) return false;
        if (nLins != other.nLins) return false;
        if (nVols_max != other.nVols_max) return false;
        if (nSurfs_max != other.nSurfs_max) return false;
        if (nLins_max != other.nLins_max) return false;
        if (n_C1P_max != other.n_C1P_max) return false;
        if (n_C2P_max != other.n_C2P_max) return false;

        bool aLinsAssoc = !Lins.empty();
        bool bLinsAssoc = !other.Lins.empty();
        if (aLinsAssoc != bLinsAssoc) return false;
        if (aLinsAssoc) {
            if (Lins != other.Lins) return false;
        }

        bool aSurfsAssoc = !Surfs.empty();
        bool bSurfsAssoc = !other.Surfs.empty();
        if (aSurfsAssoc != bSurfsAssoc) return false;
        if (aSurfsAssoc) {
            if (Surfs != other.Surfs) return false;
        }

        bool aVolsAssoc = !Vols.empty();
        bool bVolsAssoc = !other.Vols.empty();
        if (aVolsAssoc != bVolsAssoc) return false;
        if (aVolsAssoc) {
            if (Vols != other.Vols) return false;
        }

        return true;
    }
    bool operator!=(const DielectricRegions_t& other) const { return !(*this == other); }
};

struct Dielectric_t {
    int n_C1P = 0;
    int n_C2P = 0;
    double sigma = 0.0;
    double eps = 0.0;
    double mu = 0.0;
    double sigmam = 0.0;
    bool Rtime_on = false;
    bool Rtime_off = false;
    double R = 0.0;
    double L = 0.0;
    double C = 0.0;
    double R_devia = 0.0;
    double L_devia = 0.0;
    double C_devia = 0.0;
    double DiodB = 0.0;
    double DiodIsat = 0.0;
    double DiodOri = 0.0;
    int orient = 0;
    bool resistor = false;
    bool inductor = false;
    bool capacitor = false;
    bool diodo = false;
    bool plain = false;
    bool PMLbody = false;
    std::vector<int> C1P;
    std::vector<int> C2P;

    bool operator==(const Dielectric_t& other) const {
        if (n_C1P != other.n_C1P) return false;
        if (n_C2P != other.n_C2P) return false;
        if (sigma != other.sigma) return false;
        if (eps != other.eps) return false;
        if (mu != other.mu) return false;
        if (sigmam != other.sigmam) return false;
        if (Rtime_on != other.Rtime_on) return false;
        if (Rtime_off != other.Rtime_off) return false;
        if (R != other.R) return false;
        if (L != other.L) return false;
        if (C != other.C) return false;
        if (R_devia != other.R_devia) return false;
        if (L_devia != other.L_devia) return false;
        if (C_devia != other.C_devia) return false;
        if (DiodB != other.DiodB) return false;
        if (DiodIsat != other.DiodIsat) return false;
        if (DiodOri != other.DiodOri) return false;
        if (orient != other.orientation) return false; // Assuming 'orient' matches 'orientation' or similar
        if (resistor != other.resistor) return false;
        if (inductor != other.inductor) return false;
        if (capacitor != other.capacitor) return false;
        if (diodo != other.diodo) return false;
        if (plain != other.plain) return false;
        if (PMLbody != other.PMLbody) return false;
        if (C1P != other.C1P) return false;
        if (C2P != other.C2P) return false;
        return true;
    }
    bool operator!=(const Dielectric_t& other) const { return !(*this == other); }
};

struct Coords_t {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    bool operator==(const Coords_t& other) const {
        return (x == other.x) && (y == other.y) && (z == other.z);
    }
    bool operator!=(const Coords_t& other) const { return !(*this == other); }
};

struct FreqDepenMaterial_t {
    std::vector<double> a11, b11, am11, bm11;
    std::vector<double> a12, b12, am12, bm12;
    std::vector<double> a13, b13, am13, bm13;
    std::vector<double> a22, b22, am22, bm22;
    std::vector<double> a23, b23, am23, bm23;
    std::vector<double> a33, b33, am33, bm33;
    std::vector<double> alpha, beta, gamma;
    std::vector<double> alpham, betam, gammam;
    std::vector<Coords_t> c;
    double eps11 = 0.0, eps12 = 0.0, eps13 = 0.0;
    double eps22 = 0.0, eps23 = 0.0, eps33 = 0.0;
    double mu11 = 0.0, mu12 = 0.0, mu13 = 0.0;
    double mu22 = 0.0, mu23 = 0.0, mu33 = 0.0;
    double sigma11 = 0.0, sigma12 = 0.0, sigma13 = 0.0;
    double sigma22 = 0.0, sigma23 = 0.0, sigma33 = 0.0;
    double sigmam11 = 0.0, sigmam12 = 0.0, sigmam13 = 0.0;
    double sigmam22 = 0.0, sigmam23 = 0.0, sigmam33 = 0.0;
    double K11 = 0.0, Km11 = 0.0;
    double K12 = 0.0, Km12 = 0.0;
    double K13 = 0.0, Km13 = 0.0;
    double K22 = 0.0, Km22 = 0.0;
    double K23 = 0.0, Km23 = 0.0;
    double K33 = 0.0, Km33 = 0.0;
    double L = 0.0, Lm = 0.0;
    int n_c = 0;
    std::string files;

    bool operator==(const FreqDepenMaterial_t& other) const {
        if (a11 != other.a11) return false;
        if (b11 != other.b11) return false;
        if (am11 != other.am11) return false;
        if (bm11 != other.bm11) return false;
        if (a12 != other.a12) return false;
        if (b12 != other.b12) return false;
        if (am12 != other.am12) return false;
        if (bm12 != other.bm12) return false;
        if (a13 != other.a13) return false;
        if (b13 != other.b13) return false;
        if (am13 != other.am13) return false;
        if (bm13 != other.bm13) return false;
        if (a22 != other.a22) return false;
        if (b22 != other.b22) return false;
        if (am22 != other.am22) return false;
        if (bm22 != other.bm22) return false;
        if (a23 != other.a23) return false;
        if (b23 != other.b23) return false;
        if (am23 != other.am23) return false;
        if (bm23 != other.bm23) return false;
        if (a33 != other.a33) return false;
        if (b33 != other.b33) return false;
        if (am33 != other.am33) return false;
        if (bm33 != other.bm33) return false;
        if (alpha != other.alpha) return false;
        if (beta != other.beta) return false;
        if (gamma != other.gamma) return false;
        if (alpham != other.alpham) return false;
        if (betam != other.betam) return false;
        if (gammam != other.gammam) return false;
        
        // Compare coords arrays
        if (c.size() != other.c.size()) return false;
        for (size_t i = 0; i < c.size(); ++i) {
            if (c[i] != other.c[i]) return false;
        }

        if (eps11 != other.eps11) return false;
        if (eps12 != other.eps12) return false;
        if (eps13 != other.eps13) return false;
        if (eps22 != other.eps22) return false;
        if (eps23 != other.eps23) return false;
        if (eps33 != other.eps33) return false;
        if (mu11 != other.mu11) return false;
        if (mu12 != other.mu12) return false;
        if (mu13 != other.mu13) return false;
        if (mu22 != other.mu22) return false;
        if (mu23 != other.mu23) return false;
        if (mu33 != other.mu33) return false;
        if (sigma11 != other.sigma11) return false;
        if (sigma12 != other.sigma12) return false;
        if (sigma13 != other.sigma13) return false;
        if (sigma22 != other.sigma22) return false;
        if (sigma23 != other.sigma23) return false;
        if (sigma33 != other.sigma33) return false;
        if (sigmam11 != other.sigmam11) return false;
        if (sigmam12 != other.sigmam12) return false;
        if (sigmam13 != other.sigmam13) return false;
        if (sigmam22 != other.sigmam22) return false;
        if (sigmam23 != other.sigmam23) return false;
        if (sigmam33 != other.sigmam33) return false;
        if (K11 != other.K11) return false;
        if (Km11 != other.Km11) return false;
        if (K12 != other.K12) return false;
        if (Km12 != other.Km12) return false;
        if (K13 != other.K13) return false;
        if (Km13 != other.Km13) return false;
        if (K22 != other.K22) return false;
        if (Km22 != other.Km22) return false;
        if (K23 != other.K23) return false;
        if (Km23 != other.Km23) return false;
        if (K33 != other.K33) return false;
        if (Km33 != other.Km33) return false;
        if (L != other.L) return false;
        if (Lm != other.Lm) return false;
        if (n_c != other.n_c) return false;
        if (files != other.files) return false;
        
        return true;
    }
    bool operator!=(const FreqDepenMaterial_t& other) const { return !(*this == other); }
};

struct FreqDepenMaterials_t {
    int nVols = 0;
    int nSurfs = 0;
    int nLins = 0;
    int nVols_max = 0;
    int nSurfs_max = 0;
    int nLins_max = 0;
    int n_c_max = 0;
    std::vector<FreqDepenMaterial_t> Vols;
    std::vector<FreqDepenMaterial_t> Surfs;
    std::vector<FreqDepenMaterial_t> Lins;

    bool operator==(const FreqDepenMaterials_t& other) const {
        if (nVols != other.nVols) return false;
        if (nSurfs != other.nSurfs) return false;
        if (nLins != other.nLins) return false;
        if (nVols_max != other.nVols_max) return false;
        if (nSurfs_max != other.nSurfs_max) return false;
        if (nLins_max != other.nLins_max) return false;
        if (n_c_max != other.n_c_max) return false;
        if (Vols != other.Vols) return false;
        if (Surfs != other.Surfs) return false;
        if (Lins != other.Lins) return false;
        return true;
    }
    bool operator!=(const FreqDepenMaterials_t& other) const { return !(*this == other); }
};

struct ANISOTROPICbody_t {
    std::vector<int> c1P;
    std::vector<int> c2P;
    std::vector<double> sigma;
    std::vector<double> eps;
    std::vector<double> mu;
    std::vector<double> sigmam;
    int n_C1P = 0;
    int n_C2P = 0;

    bool operator==(const ANISOTROPICbody_t& other) const {
        if (c1P != other.c1P) return false;
        if (c2P != other.c2P) return false;
        if (sigma != other.sigma) return false;
        if (eps != other.eps) return false;
        if (mu != other.mu) return false;
        if (sigmam != other.sigmam) return false;
        if (n_C1P != other.n_C1P) return false;
        if (n_C2P != other.n_C2P) return false;
        return true;
    }
    bool operator!=(const ANISOTROPICbody_t& other) const { return !(*this == other); }
};

struct ANISOTROPICelements_t {
    int nVols = 0;
    int nSurfs = 0;
    int nLins = 0;
    int nVols_max = 0;
    int nSurfs_max = 0;
    int nLins_max = 0;
    int n_C1P_max = 0;
    int n_C2P_max = 0;
    std::vector<ANISOTROPICbody_t> Vols;
    std::vector<ANISOTROPICbody_t> Surfs;
    std::vector<ANISOTROPICbody_t> Lins;

    bool operator==(const ANISOTROPICelements_t& other) const {
        if (nVols != other.nVols) return false;
        if (nSurfs != other.nSurfs) return false;
        if (nLins != other.nLins) return false;
        if (nVols_max != other.nVols_max) return false;
        if (nSurfs_max != other.nSurfs_max) return false;
        if (nLins_max != other.nLins_max) return false;
        if (n_C1P_max != other.n_C1P_max) return false;
        if (n_C2P_max != other.n_C2P_max) return false;
        if (Vols != other.Vols) return false;
        if (Surfs != other.Surfs) return false;
        if (Lins != other.Lins) return false;
        return true;
    }
    bool operator!=(const ANISOTROPICelements_t& other) const { return !(*this == other); }
};

struct LossyThinSurface_t {
    int nc = 0;
    std::string files;
    int numcapas = 0;
    std::vector<Coords_t> c;
    std::vector<double> sigma;
    std::vector<double> eps;
    std::vector<double> mu;
    std::vector<double> sigmam;
    std::vector<double> thk;
    std::vector<double> sigma_devia;
    std::vector<double> eps_devia;
    std::vector<double> mu_devia;
    std::vector<double> sigmam_devia;
    std::vector<double> thk_devia;

    bool operator==(const LossyThinSurface_t& other) const {
        if (nc != other.nc) return false;
        if (files != other.files) return false;
        if (numcapas != other.numcapas) return false;
        if (c != other.c) return false;
        if (sigma != other.sigma) return false;
        if (eps != other.eps) return false;
        if (mu != other.mu) return false;
        if (sigmam != other.sigmam) return false;
        if (thk != other.thk) return false;
        if (sigma_devia != other.sigma_devia) return false;
        if (eps_devia != other.eps_devia) return false;
        if (mu_devia != other.mu_devia) return false;
        if (sigmam_devia != other.sigmam_devia) return false;
        if (thk_devia != other.thk_devia) return false;
        return true;
    }
    bool operator!=(const LossyThinSurface_t& other) const { return !(*this == other); }
};

struct LossyThinSurfaces_t {
    std::vector<LossyThinSurface_t> cs;
    int length = 0;
    int length_max = 0;
    int nC_max = 0;

    bool operator==(const LossyThinSurfaces_t& other) const {
        if (cs != other.cs) return false;
        if (length != other.length) return false;
        if (length_max != other.length_max) return false;
        if (nC_max != other.nC_max) return false;
        return true;
    }
    bool operator!=(const LossyThinSurfaces_t& other) const { return !(*this == other); }
};

struct ThinWireComp_t {
    int srctype = 0;
    std::string srcfile;
    int i = 0;
    int j = 0;
    int K = 0;
    int nd = 0;
    double d = 0.0;
    int m = 0;
    std::string tag;

    bool operator==(const ThinWireComp_t& other) const {
        if (srctype != other.srctype) return false;
        if (srcfile != other.srcfile) return false;
        if (i != other.i) return false;
        if (j != other.j) return false;
        if (K != other.K) return false;
        if (nd != other.nd) return false;
        if (d != other.d) return false;
        if (m != other.m) return false;
        if (tag != other.tag) return false;
        return true;
    }
    bool operator!=(const ThinWireComp_t& other) const { return !(*this == other); }
};

struct ThinWire_t {
    double rad = 0.0;
    bool disp = false;
    std::string dispfile;
    double res = 0.0;
    double ind = 0.0;
    double cap = 0.0;
    double P_res = 0.0;
    double P_ind = 0.0;
    double P_cap = 0.0;
    std::string dispfile_LeftEnd;
    double R_LeftEnd = 0.0;
    double L_LeftEnd = 0.0;
    double C_LeftEnd = 0.0;
    std::string dispfile_RightEnd;
    double R_RightEnd = 0.0;
    double L_RightEnd = 0.0;
    double C_RightEnd = 0.0;
    Coords_t LeftEnd;
    Coords_t RightEnd;
    double tl = 0.0;
    double tr = 0.0;
    std::vector<ThinWireComp_t> twc;
    int n_twc = 0;
    int n_twc_max = 0;

    bool operator==(const ThinWire_t& other) const {
        if (rad != other.rad) return false;
        if (disp != other.disp) return false;
        if (dispfile != other.dispfile) return false;
        if (res != other.res) return false;
        if (ind != other.ind) return false;
        if (cap != other.cap) return false;
        if (P_res != other.P_res) return false;
        if (P_ind != other.P_ind) return false;
        if (P_cap != other.P_cap) return false;
        if (dispfile_LeftEnd != other.dispfile_LeftEnd) return false;
        if (R_LeftEnd != other.R_LeftEnd) return false;
        if (L_LeftEnd != other.L_LeftEnd) return false;
        if (C_LeftEnd != other.C_LeftEnd) return false;
        if (dispfile_RightEnd != other.dispfile_RightEnd) return false;
        if (R_RightEnd != other.R_RightEnd) return false;
        if (L_RightEnd != other.L_RightEnd) return false;
        if (C_RightEnd != other.C_RightEnd) return false;
        if (LeftEnd != other.LeftEnd) return false;
        if (RightEnd != other.RightEnd) return false;
        if (tl != other.tl) return false;
        if (tr != other.tr) return false;
        if (twc != other.twc) return false;
        if (n_twc != other.n_twc) return false;
        if (n_twc_max != other.n_twc_max) return false;
        return true;
    }
    bool operator!=(const ThinWire_t& other) const { return !(*this == other); }
};

struct ThinWires_t {
    std::vector<ThinWire_t> tw;
    int n_tw = 0;
    int n_tw_max = 0;

    bool operator==(const ThinWires_t& other) const {
        if (n_tw != other.n_tw) return false;
        if (n_tw_max != other.n_tw_max) return false;
        
        bool aTwAssoc = !tw.empty();
        bool bTwAssoc = !other.tw.empty();
        if (aTwAssoc != bTwAssoc) return false;
        if (aTwAssoc) {
            if (tw.size() != other.tw.size()) return false;
            for (size_t i = 0; i < tw.size(); ++i) {
                if (tw[i] != other.tw[i]) return false;
            }
        }
        return true;
    }
    bool operator!=(const ThinWires_t& other) const { return !(*this == other); }
};

struct SlantedWireComp_t {
    int srctype = 0;
    std::string srcfile;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    int nd = 0;
    int m = 0;
    std::string tag;

    bool operator==(const SlantedWireComp_t& other) const {
        if (srctype != other.srctype) return false;
        if (srcfile != other.srcfile) return false;
        if (x != other.x) return false;
        if (y != other.y) return false;
        if (z != other.z) return false;
        if (nd != other.nd) return false;
        if (m != other.m) return false;
        if (tag != other.tag) return false;
        return true;
    }
    bool operator!=(const SlantedWireComp_t& other) const { return !(*this == other); }
};

struct SlantedWire_t {
    double rad = 0.0;
    bool disp = false;
    std::string dispfile;
    double res = 0.0;
    double ind = 0.0;
    double cap = 0.0;
    double P_res = 0.0;
    double P_ind = 0.0;
    double P_cap = 0.0;
    std::string dispfile_LeftEnd;
    double R_LeftEnd = 0.0;
    double L_LeftEnd = 0.0;
    double C_LeftEnd = 0.0;
    std::string dispfile_RightEnd;
    double R_RightEnd = 0.0;
    double L_RightEnd = 0.0;
    double C_RightEnd = 0.0;
    Coords_t LeftEnd;
    Coords_t RightEnd;
    double tl = 0.0;
    double tr = 0.0;
    int n_swc = 0;
    int n_swc_max = 0;

    bool operator==(const SlantedWire_t& other) const {
        if (rad != other.rad) return false;
        if (disp != other.disp) return false;
        if (dispfile != other.dispfile) return false;
        if (res != other.res) return false;
        if (ind != other.ind) return false;
        if (cap != other.cap) return false;
        if (P_res != other.P_res) return false;
        if (P_ind != other.P_ind) return false;
        if (P_cap != other.P_cap) return false;
        if (dispfile_LeftEnd != other.dispfile_LeftEnd) return false;
        if (R_LeftEnd != other.R_LeftEnd) return false;
        if (L_LeftEnd != other.L_LeftEnd) return false;
        if (C_LeftEnd != other.C_LeftEnd) return false;
        if (dispfile_RightEnd != other.dispfile_RightEnd) return false;
        if (R_RightEnd != other.R_RightEnd) return false;
        if (L_RightEnd != other.L_RightEnd) return false;
        if (C_RightEnd != other.C_RightEnd) return false;
        if (LeftEnd != other.LeftEnd) return false;
        if (RightEnd != other.RightEnd) return false;
        if (tl != other.tl) return false;
        if (tr != other.tr) return false;
        if (n_swc != other.n_swc) return false;
        if (n_swc_max != other.n_swc_max) return false;
        return true;
    }
    bool operator!=(const SlantedWire_t& other) const { return !(*this == other); }
};

struct SlantedWiresInfo_t {
    std::vector<SlantedWire_t> sw;
    int n_sw = 0;
    int n_sw_max = 0;

    bool operator==(const SlantedWiresInfo_t& other) const {
        if (sw.size() != other.sw.size()) return false;
        for (size_t i = 0; i < sw.size(); ++i) {
            if (sw[i] != other.sw[i]) return false;
        }
        if (n_sw != other.n_sw) return false;
        if (n_sw_max != other.n_sw_max) return false;
        return true;
    }
    bool operator!=(const SlantedWiresInfo_t& other) const { return !(*this == other); }
};

struct ThinSlotComp_t {
    int i = 0;
    int j = 0;
    int K = 0;
    int node = 0;
    int dir = 0;
    int Or = 0;
    std::string tag;

    bool operator==(const ThinSlotComp_t& other) const {
        if (i != other.i) return false;
        if (j != other.j) return false;
        if (K != other.K) return false;
        if (node != other.node) return false;
        if (dir != other.dir) return false;
        if (Or != other.Or) return false;
        if (tag != other.tag) return false;
        return true;
    }
    bool operator!=(const ThinSlotComp_t& other) const { return !(*this == other); }
};

struct ThinSlot_t {
    double width = 0.0;
    int n_tgc = 0;
    int n_tgc_max = 0;
    std::vector<ThinSlotComp_t> tgc;

    bool operator==(const ThinSlot_t& other) const {
        if (width != other.width) return false;
        if (n_tgc != other.n_tgc) return false;
        if (n_tgc_max != other.n_tgc_max) return false;
        
        bool aTgcAssoc = !tgc.empty();
        bool bTgcAssoc = !other.tgc.empty();
        if (aTgcAssoc != bTgcAssoc) return false;
        if (aTgcAssoc) {
            if (tgc != other.tgc) return false;
        }
        return true;
    }
    bool operator!=(const ThinSlot_t& other) const { return !(*this == other); }
};

struct ThinSlots_t {
    int n_tg = 0;
    int n_tg_max = 0;
    std::vector<ThinSlot_t> tg;

    bool operator==(const ThinSlots_t& other) const {
        if (n_tg != other.n_tg) return false;
        if (n_tg_max != other.n_tg_max) return false;
        
        bool aTgAssoc = !tg.empty();
        bool bTgAssoc = !other.tg.empty();
        if (aTgAssoc != bTgAssoc) return false;
        if (aTgAssoc) {
            if (tg != other.tg) return false;
        }
        return true;
    }
    bool operator!=(const ThinSlots_t& other) const { return !(*this == other); }
};

struct Box_t {
    std::string nombre_fichero;
    std::vector<double> coor1;
    std::vector<double> coor2;

    bool operator==(const Box_t& other) const {
        if (nombre_fichero != other.nombre_fichero) return false;
        if (coor1 != other.coor1) return false;
        if (coor2 != other.coor2) return false;
        return true;
    }
    bool operator!=(const Box_t& other) const { return !(*this == other); }
};

struct Boxes_t {
    int nVols = 0;
    int nVols_max = 0;
    std::vector<Box_t> Vols;

    bool operator==(const Boxes_t& other) const {
        if (nVols != other.nVols) return false;
        if (nVols_max != other.nVols_max) return false;
        if (Vols != other.Vols) return false;
        return true;
    }
    bool operator!=(const Boxes_t& other) const { return !(*this == other); }
};

struct PlaneWave_t {
    std::string nombre_fichero;
    int atributo = 0;
    std::vector<double> coor1;
    std::vector<double> coor2;
    double theta = 0.0;
    double phi = 0.0;
    double alpha = 0.0;
    double beta = 0.0;
    bool isRC = false;

    bool operator==(const PlaneWave_t& other) const {
        if (nombre_fichero != other.nombre_fichero) return false;
        if (atributo != other.atributo) return false;
        if (coor1 != other.coor1) return false;
        if (coor2 != other.coor2) return false;
        if (theta != other.theta) return false;
        if (phi != other.phi) return false;
        if (alpha != other.alpha) return false;
        if (beta != other.beta) return false;
        if (isRC != other.isRC) return false;
        return true;
    }
    bool operator!=(const PlaneWave_t& other) const { return !(*this == other); }
};

struct PlaneWaves_t {
    std::vector<PlaneWave_t> collection;
    // Other fields might exist but are not in the snippet
    bool operator==(const PlaneWaves_t& other) const {
        if (collection != other.collection) return false;
        return true;
    }
    bool operator!=(const PlaneWaves_t& other) const { return !(*this == other); }
};

// Additional placeholder types for the rest of the Parseador_t members
struct Desplazamiento_t { bool operator==(const Desplazamiento_t&) const { return true; } bool operator!=(const Desplazamiento_t&) const { return false; } };
struct Frontera_t { bool operator==(const Frontera_t&) const { return true; } bool operator!=(const Frontera_t&) const { return false; } };
struct FronteraPML_t { bool operator==(const FronteraPML_t&) const { return true; } bool operator!=(const FronteraPML_t&) const { return false; } };
struct CurrFieldSrc_t { bool operator==(const CurrFieldSrc_t&) const { return true; } bool operator!=(const CurrFieldSrc_t&) const { return false; } };
struct NodSource_t { bool operator==(const NodSource_t&) const { return true; } bool operator!=(const NodSource_t&) const { return false; } };
struct CoordsScaled_t { bool operator==(const CoordsScaled_t&) const { return true; } bool operator!=(const CoordsScaled_t&) const { return false; } };
struct AbstractSonda_t { bool operator==(const AbstractSonda_t&) const { return true; } bool operator!=(const AbstractSonda_t&) const { return false; } };
struct Sonda_t { bool operator==(const Sonda_t&) const { return true; } bool operator!=(const Sonda_t&) const { return false; } };
struct Sondas_t { bool operator==(const Sondas_t&) const { return true; } bool operator!=(const Sondas_t&) const { return false; } };
struct Massonda_t { bool operator==(const Massonda_t&) const { return true; } bool operator!=(const Massonda_t&) const { return false; } };
struct Massondas_t { bool operator==(const Massondas_t&) const { return true; } bool operator!=(const Massondas_t&) const { return false; } };
struct Bloqueprobe_t { bool operator==(const Bloqueprobe_t&) const { return true; } bool operator!=(const Bloqueprobe_t&) const { return false; } };
struct Bloqueprobes_t { bool operator==(const Bloqueprobes_t&) const { return true; } bool operator!=(const Bloqueprobes_t&) const { return false; } };
struct Volprobe_t { bool operator==(const Volprobe_t&) const { return true; } bool operator!=(const Volprobe_t&) const { return false; } };
struct Volprobes_t { bool operator==(const Volprobes_t&) const { return true; } bool operator!=(const Volprobes_t&) const { return false; } };
struct FarField_Sonda_t { bool operator==(const FarField_Sonda_t&) const { return true; } bool operator!=(const FarField_Sonda_t&) const { return false; } };
struct Electric_Sonda_t { bool operator==(const Electric_Sonda_t&) const { return true; } bool operator!=(const Electric_Sonda_t&) const { return false; } };
struct Magnetic_Sonda_t { bool operator==(const Magnetic_Sonda_t&) const { return true; } bool operator!=(const Magnetic_Sonda_t&) const { return false; } };
struct NormalElectric_Sonda_t { bool operator==(const NormalElectric_Sonda_t&) const { return true; } bool operator!=(const NormalElectric_Sonda_t&) const { return false; } };
struct NormalMagnetic_Sonda_t { bool operator==(const NormalMagnetic_Sonda_t&) const { return true; } bool operator!=(const NormalMagnetic_Sonda_t&) const { return false; } };
struct SurfaceElectricCurrent_Sonda_t { bool operator==(const SurfaceElectricCurrent_Sonda_t&) const { return true; } bool operator!=(const SurfaceElectricCurrent_Sonda_t&) const { return false; } };
struct SurfaceMagneticCurrent_Sonda_t { bool operator==(const SurfaceMagneticCurrent_Sonda_t&) const { return true; } bool operator!=(const SurfaceMagneticCurrent_Sonda_t&) const { return false; } };
struct OldSonda_t { 
    std::vector<Bloqueprobe_t> probes;
    int n_probes = 0;
    int n_probes_max = 0;
    bool operator==(const OldSonda_t& other) const {
        if (probes != other.probes) return false;
        if (n_probes != other.n_probes) return false;
        if (n_probes_max != other.n_probes_max) return false;
        return true;
    }
    bool operator!=(const OldSonda_t& other) const { return !(*this == other); }
};
struct ConformalRegs_t {
    std::vector<int> volumes;
    std::vector<int> surfaces;
    bool operator==(const ConformalRegs_t& other) const {
        if (volumes != other.volumes) return false;
        if (surfaces != other.surfaces) return false;
        return true;
    }
    bool operator!=(const ConformalRegs_t& other) const { return !(*this == other); }
};
struct MTLN_t {
    std::vector<int> cables;
    std::vector<int> probes;
    std::vector<int> networks;
    std::vector<int> connectors;
    bool operator==(const MTLN_t& other) const {
        if (cables != other.cables) return false;
        if (probes != other.probes) return false;
        if (networks != other.networks) return false;
        if (connectors != other.connectors) return false;
        return true;
    }
    bool operator!=(const MTLN_t& other) const { return !(*this == other); }
};

struct Parseador_t {
    Switches_t switches;
    NFDEGeneral_t general;
    MatrizMedios_t matriz;
    Desplazamiento_t despl;
    Frontera_t front;
    FronteraPML_t frontPML; // Assuming this exists based on interface
    Materials_t Mats;
    PECRegions_t pecRegs;
    PMCRegions_t pmcRegs;
    DielectricRegions_t DielRegs;
    LossyThinSurfaces_t LossyThinSurfs;
    FreqDepenMaterials_t frqDepMats;
    ANISOTROPICelements_t aniMats;
    Boxes_t boxSrc;
    PlaneWaves_t plnSrc;
    NodSource_t nodSrc;
    OldSonda_t oldSONDA;
    Sondas_t Sonda;
    Bloqueprobes_t BloquePrb;
    Volprobes_t VolPrb;
    ThinWires_t tWires;
    SlantedWiresInfo_t sWires;
    ThinSlots_t tSlots;
    ConformalRegs_t conformalRegs;
    
#ifdef CompileWithMTLN
    MTLN_t mtln;
#endif

    bool operator==(const Parseador_t& other) const {
        if (switches != other.switches) return false;
        if (general != other.general) return false;
        if (matriz != other.matriz) return false;
        if (despl != other.despl) return false;
        if (front != other.front) return false;
        if (Mats != other.Mats) return false;
        if (pecRegs != other.pecRegs) return false;
        if (pmcRegs != other.pmcRegs) return false;
        if (DielRegs != other.DielRegs) return false;
        if (LossyThinSurfs != other.LossyThinSurfs) return false;
        if (frqDepMats != other.frqDepMats) return false;
        if (aniMats != other.aniMats) return false;
        if (boxSrc != other.boxSrc) return false;
        if (plnSrc != other.plnSrc) return false;
        if (nodSrc != other.nodSrc) return false;
        if (oldSONDA != other.oldSONDA) return false;
        if (Sonda != other.Sonda) return false;
        if (BloquePrb != other.BloquePrb) return false;
        if (VolPrb != other.VolPrb) return false;
        if (tWires != other.tWires) return false;
        if (sWires != other.sWires) return false;
        if (tSlots != other.tSlots) return false;
        if (conformalRegs != other.conformalRegs) return false;
#ifdef CompileWithMTLN
        if (mtln != other.mtln) return false;
#endif
        return true;
    }
    bool operator!=(const Parseador_t& other) const { return !(*this == other); }
};

namespace NFDETypes_extension_m {

    void initializeProblemDescription(Parseador_t& pD) {
        pD.general = NFDEGeneral_t();
        pD.matriz = MatrizMedios_t();
        pD.despl = Desplazamiento_t();
        pD.front = Frontera_t();

        pD.Mats.n_Mats = 3;
        pD.Mats.n_Mats_max = 3;
        pD.Mats.mats.resize(3);

        pD.Mats.mats[0].id = 1;
        pD.Mats.mats[0].eps = EPSILON_VACUUM;
        pD.Mats.mats[0].mu = MU_VACUUM;
        pD.Mats.mats[0].sigma = 0.0;
        pD.Mats.mats[0].sigmam = 0.0;

        pD.Mats.mats[1].id = 2;
        pD.Mats.mats[1].eps = EPSILON_VACUUM;
        pD.Mats.mats[1].mu = MU_VACUUM;
        pD.Mats.mats[1].sigma = SIGMA_PEC;
        pD.Mats.mats[1].sigmam = 0.0;

        pD.Mats.mats[2].id = 3;
        pD.Mats.mats[2].eps = EPSILON_VACUUM;
        pD.Mats.mats[2].mu = MU_VACUUM;
        pD.Mats.mats[2].sigma = 0.0;
        pD.Mats.mats[2].sigmam = SIGMA_PMC;

        pD.pecRegs = PECRegions_t();
        pD.pecRegs.Lins.clear();
        pD.pecRegs.Surfs.clear();
        pD.pecRegs.Vols.clear();

        pD.pmcRegs = PMCRegions_t();
        pD.pmcRegs.Lins.clear();
        pD.pmcRegs.Surfs.clear();
        pD.pmcRegs.Vols.clear();

        pD.DielRegs = DielectricRegions_t();
        pD.DielRegs.Lins.clear();
        pD.DielRegs.Surfs.clear();
        pD.DielRegs.Vols.clear();
      
        pD.LossyThinSurfs = LossyThinSurfaces_t();
        pD.LossyThinSurfs.cs.clear();

        pD.frqDepMats = FreqDepenMaterials_t();
        pD.aniMats = ANISOTROPICelements_t();
        
        pD.plnSrc = PlaneWaves_t();
        pD.plnSrc.collection.clear();
        pD.nodSrc = NodSource_t();
        pD.nodSrc.collection.clear(); // Assuming collection field for NodSource
        pD.boxSrc = Boxes_t();
        
        pD.Sonda = Sondas_t();
        pD.Sonda.len_cor_max = 0; // Placeholder
        pD.Sonda.length = 0; // Placeholder
        pD.Sonda.length_max = 0; // Placeholder
        pD.Sonda.collection.clear(); // Placeholder
        
        pD.oldSONDA = OldSonda_t();
        pD.oldSONDA.probes.clear();
        pD.oldSONDA.n_probes = 0;
        pD.oldSONDA.n_probes_max = 0;

        pD.BloquePrb = Bloqueprobes_t();
        pD.BloquePrb.bp.clear(); // Placeholder
        
        pD.VolPrb = Volprobes_t();
        pD.VolPrb.collection.clear(); // Placeholder
        
        pD.tWires = ThinWires_t();
        pD.tWires.tw.clear();
        pD.tWires.n_tw = 0;
        pD.tWires.n_tw_max = 0; // Fixed typo from original Fortran

        pD.sWires = SlantedWiresInfo_t();
        pD.tSlots = ThinSlots_t();
        pD.tSlots.tg.clear();

#ifdef CompileWithMTLN
        pD.mtln = MTLN_t();
        pD.mtln.cables.clear();
        pD.mtln.probes.clear();
        pD.mtln.networks.clear();
        pD.mtln.connectors.clear();
#endif

        pD.conformalRegs = ConformalRegs_t();
        pD.conformalRegs.volumes.clear();
        pD.conformalRegs.surfaces.clear();
    }

    bool Parseador_eq(const Parseador_t& a, const Parseador_t& b) {
        return a == b;
    }

    bool NFDEGeneral_eq(const NFDEGeneral_t& a, const NFDEGeneral_t& b) {
        return a == b;
    }

    bool desplazamiento_eq(const Desplazamiento_t& a, const Desplazamiento_t& b) {
        return a == b;
    }

    bool MatrizMedios_eq(const MatrizMedios_t& a, const MatrizMedios_t& b) {
        return a == b;
    }

    bool Material_eq(const Material_t& a, const Material_t& b) {
        return a == b;
    }

    bool Materials_eq(const Materials_t& a, const Materials_t& b) {
        return a == b;
    }

    bool pecregions_eq(const PECRegions_t& a, const PECRegions_t& b) {
        return a == b;
    }

    bool dielectric_eq(const Dielectric_t& a, const Dielectric_t& b) {
        return a == b;
    }

    bool dielectricregions_eq(const DielectricRegions_t& a, const DielectricRegions_t& b) {
        return a == b;
    }

    bool freqdepenmaterial_eq(const FreqDepenMaterial_t& a, const FreqDepenMaterial_t& b) {
        return a == b;
    }

    bool freqdepenmaterials_eq(const FreqDepenMaterials_t& a, const FreqDepenMaterials_t& b) {
        return a == b;
    }

    bool anisotropicbody_eq(const ANISOTROPICbody_t& a, const ANISOTROPICbody_t& b) {
        return a == b;
    }

    bool anisotropicelements_eq(const ANISOTROPICelements_t& a, const ANISOTROPICelements_t& b) {
        return a == b;
    }

    bool LossyThinSurface_eq(const LossyThinSurface_t& a, const LossyThinSurface_t& b) {
        return a == b;
    }

    bool LossyThinSurfaces_eq(const LossyThinSurfaces_t& a, const LossyThinSurfaces_t& b) {
        return a == b;
    }

    bool ThinWireComp_eq(const ThinWireComp_t& a, const ThinWireComp_t& b) {
        return a == b;
    }

    bool ThinWire_eq(const ThinWire_t& a, const ThinWire_t& b) {
        return a == b;
    }

    bool ThinWires_eq(const ThinWires_t& a, const ThinWires_t& b) {
        return a == b;
    }

    bool SlantedWireComp_eq(const SlantedWireComp_t& a, const SlantedWireComp_t& b) {
        return a == b;
    }

    bool SlantedWire_eq(const SlantedWire_t& a, const SlantedWire_t& b) {
        return a == b;
    }

    bool SlantedWires_eq(const SlantedWiresInfo_t& a, const SlantedWiresInfo_t& b) {
        return a == b;
    }

    bool ThinSlotComp_eq(const ThinSlotComp_t& a, const ThinSlotComp_t& b) {
        return a == b;
    }

    bool ThinSlot_eq(const ThinSlot_t& a, const ThinSlot_t& b) {
        return a == b;
    }

    bool ThinSlots_eq(const ThinSlots_t& a, const ThinSlots_t& b) {
        return a == b;
    }

    bool box_eq(const Box_t& a, const Box_t& b) {
        return a == b;
    }

    bool boxes_eq(const Boxes_t& a, const Boxes_t& b) {
        return a == b;
    }

    bool planewave_eq(const PlaneWave_t& a, const PlaneWave_t& b) {
        return a == b;
    }

    bool planewaves_eq(const PlaneWaves_t& a, const PlaneWaves_t& b) {
        return a == b;
    }

    bool curr_field_src_eq(const CurrFieldSrc_t& a, const CurrFieldSrc_t& b) {
        return a == b;
    }

    bool nodsource_eq(const NodSource_t& a, const NodSource_t& b) {
        return a == b;
    }

    bool coords_eq(const Coords_t& a, const Coords_t& b) {
        return a == b;
    }

    bool coords_scaled_eq(const CoordsScaled_t& a, const CoordsScaled_t& b) {
        return a == b;
    }

    bool abstractSonda_eq(const AbstractSonda_t& a, const AbstractSonda_t& b) {
        return a == b;
    }

    bool sonda_eq(const Sonda_t& a, const Sonda_t& b) {
        return a == b;
    }

    bool sondas_eq(const Sondas_t& a, const Sondas_t& b) {
        return a == b;
    }

    bool massonda_eq(const Massonda_t& a, const Massonda_t& b) {
        return a == b;
    }

    bool massondas_eq(const Massondas_t& a, const Massondas_t& b) {
        return a == b;
    }

    bool bloqueprobe_eq(const Bloqueprobe_t& a, const Bloqueprobe_t& b) {
        return a == b;
    }

    bool bloqueprobes_eq(const Bloqueprobes_t& a, const Bloqueprobes_t& b) {
        return a == b;
    }

    bool volprobe_eq(const Volprobe_t& a, const Volprobe_t& b) {
        return a == b;
    }

    bool volprobes_eq(const Volprobes_t& a, const Volprobes_t& b) {
        return a == b;
    }

    bool FarField_Sonda_eq(const FarField_Sonda_t& a, const FarField_Sonda_t& b) {
        return a == b;
    }

    bool Electric_Sonda_eq(const Electric_Sonda_t& a, const Electric_Sonda_t& b) {
        return a == b;
    }

    bool Magnetic_Sonda_eq(const Magnetic_Sonda_t& a, const Magnetic_Sonda_t& b) {
        return a == b;
    }

    bool NormalElectric_Sonda_eq(const NormalElectric_Sonda_t& a, const NormalElectric_Sonda_t& b) {
        return a == b;
    }

    bool NormalMagnetic_Sonda_eq(const NormalMagnetic_Sonda_t& a, const NormalMagnetic_Sonda_t& b) {
        return a == b;
    }

    bool SurfaceElectricCurrent_Sonda_eq(const SurfaceElectricCurrent_Sonda_t& a, const SurfaceElectricCurrent_Sonda_t& b) {
        return a == b;
    }

    bool SurfaceMagneticCurrent_Sonda_eq(const SurfaceMagneticCurrent_Sonda_t& a, const SurfaceMagneticCurrent_Sonda_t& b) {
        return a == b;
    }

    bool frontera_eq(const Frontera_t& a, const Frontera_t& b) {
        return a == b;
    }

    bool fronteraPML_eq(const FronteraPML_t& a, const FronteraPML_t& b) {
        return a == b;
    }

} // namespace NFDETypes_extension_m

if (a.INCERTMAX != b.INCERTMAX) return false;
   if (a.numModes != b.numModes) return false;
   return true;
}

bool planewaves_eq(const PlaneWaves_t& a, const PlaneWaves_t& b) {
   if (a.nc != b.nc) return false;
   if (a.nC_max != b.nC_max) return false;
   if (a.collection.empty() && b.collection.empty()) {
      return true;
   } else if (a.collection.empty() != b.collection.empty()) {
      return false;
   } else {
      if (a.collection.size() != b.collection.size()) return false;
      for (size_t i = 0; i < a.collection.size(); ++i) {
         if (a.collection[i] != b.collection[i]) return false;
      }
      return true;
   }
}

bool curr_field_src_eq(const Curr_Field_Src_t& a, const Curr_Field_Src_t& b) {
   if (a.n_C1P != b.n_C1P) return false;
   if (a.n_C2P != b.n_C2P) return false;
   
   // Check c1P
   if (a.c1P.empty() && b.c1P.empty()) {
      // OK
   } else if (a.c1P.empty() != b.c1P.empty()) {
      return false;
   } else {
      if (a.c1P.size() != b.c1P.size()) return false;
      for (size_t i = 0; i < a.c1P.size(); ++i) {
         if (a.c1P[i] != b.c1P[i]) return false;
      }
   }

   // Check c2P
   if (a.c2P.empty() && b.c2P.empty()) {
      // OK
   } else if (a.c2P.empty() != b.c2P.empty()) {
      return false;
   } else {
      if (a.c2P.size() != b.c2P.size()) return false;
      for (size_t i = 0; i < a.c2P.size(); ++i) {
         if (a.c2P[i] != b.c2P[i]) return false;
      }
   }

   if (a.nombre != b.nombre) return false;
   if (a.isElec != b.isElec) return false;
   if (a.isHard != b.isHard) return false;
   if (a.isInitialValue != b.isInitialValue) return false;
   return true;
}

bool nodsource_eq(const NodSource_t& a, const NodSource_t& b) {
   if (a.n_nodSrc != b.n_nodSrc) return false;
   if (a.n_nodSrc_max != b.n_nodSrc_max) return false;
   if (a.n_C1P_max != b.n_C1P_max) return false;
   if (a.n_C2P_max != b.n_C2P_max) return false;

   if (a.NodalSource.empty() && b.NodalSource.empty()) {
      return true;
   } else if (a.NodalSource.empty() != b.NodalSource.empty()) {
      return false;
   } else {
      if (a.NodalSource.size() != b.NodalSource.size()) return false;
      for (size_t i = 0; i < a.NodalSource.size(); ++i) {
         if (a.NodalSource[i] != b.NodalSource[i]) return false;
      }
      return true;
   }
}

bool fronteraPML_eq(const FronteraPML_t& a, const FronteraPML_t& b) {
   if (a.orden != b.orden) return false;
   if (a.refl != b.refl) return false;
   if (a.numCapas != b.numCapas) return false;
   return true;
}

bool frontera_eq(const Frontera_t& a, const Frontera_t& b) {
   if (a.tipoFrontera != b.tipoFrontera) return false;
   if (a.propiedadesPML != b.propiedadesPML) return false;
   return true;
}

bool desplazamiento_eq(const Desplazamiento_t& a, const Desplazamiento_t& b) {
   if (a.desX.empty() || a.desY.empty() || a.desZ.empty() ||
       b.desX.empty() || b.desY.empty() || b.desZ.empty()) {
      return false;
   }

   if (a.desX.size() != b.desX.size()) return false;
   if (a.desY.size() != b.desY.size()) return false;
   if (a.desZ.size() != b.desZ.size()) return false;

   for (size_t i = 0; i < a.desX.size(); ++i) {
      if (a.desX[i] != b.desX[i]) return false;
   }
   for (size_t i = 0; i < a.desY.size(); ++i) {
      if (a.desY[i] != b.desY[i]) return false;
   }
   for (size_t i = 0; i < a.desZ.size(); ++i) {
      if (a.desZ[i] != b.desZ[i]) return false;
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

bool coords_eq(const coords_t& a, const coords_t& b) {
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

bool coords_scaled_eq(const coords_scaled_t& a, const coords_scaled_t& b) {
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

bool FarField_Sonda_eq(const FarField_Sonda_t& a, const FarField_Sonda_t& b) {
   if (a.probe != b.probe) return false;
   return true;
}

bool Electric_Sonda_eq(const Electric_Sonda_t& a, const Electric_Sonda_t& b) {
   if (a.probe != b.probe) return false;
   return true;
}

bool Magnetic_Sonda_eq(const Magnetic_Sonda_t& a, const Magnetic_Sonda_t& b) {
   if (a.probe != b.probe) return false;
   return true;
}

bool NormalElectric_Sonda_eq(const NormalElectric_Sonda_t& a, const NormalElectric_Sonda_t& b) {
   if (a.probe != b.probe) return false;
   if (a.n_nml != b.n_nml) return false;
   if (a.n_nml_max != b.n_nml_max) return false;

   if (a.nml.empty() && b.nml.empty()) {
      return true;
   } else if (a.nml.empty() != b.nml.empty()) {
      return false;
   } else {
      if (a.nml.size() != b.nml.size()) return false;
      for (size_t i = 0; i < a.nml.size(); ++i) {
         if (a.nml[i] != b.nml[i]) return false;
      }
      return true;
   }
}

bool NormalMagnetic_Sonda_eq(const NormalMagnetic_Sonda_t& a, const NormalMagnetic_Sonda_t& b) {
   if (a.probe != b.probe) return false;
   if (a.n_nml != b.n_nml) return false;
   if (a.n_nml_max != b.n_nml_max) return false;

   if (a.nml.empty() && b.nml.empty()) {
      return true;
   } else if (a.nml.empty() != b.nml.empty()) {
      return false;
   } else {
      if (a.nml.size() != b.nml.size()) return false;
      for (size_t i = 0; i < a.nml.size(); ++i) {
         if (a.nml[i] != b.nml[i]) return false;
      }
      return true;
   }
}

bool SurfaceElectricCurrent_Sonda_eq(const SurfaceElectricCurrent_Sonda_t& a, const SurfaceElectricCurrent_Sonda_t& b) {
   if (a.probe != b.probe) return false;
   if (a.n_nml != b.n_nml) return false;
   if (a.n_nml_max != b.n_nml_max) return false;

   if (a.nml.empty() && b.nml.empty()) {
      return true;
   } else if (a.nml.empty() != b.nml.empty()) {
      return false;
   } else {
      if (a.nml.size() != b.nml.size()) return false;
      for (size_t i = 0; i < a.nml.size(); ++i) {
         if (a.nml[i] != b.nml[i]) return false;
      }
      return true;
   }
}

bool SurfaceMagneticCurrent_Sonda_eq(const SurfaceMagneticCurrent_Sonda_t& a, const SurfaceMagneticCurrent_Sonda_t& b) {
   if (a.probe != b.probe) return false;
   if (a.n_nml != b.n_nml) return false;
   if (a.n_nml_max != b.n_nml_max) return false;

   if (a.nml.empty() && b.nml.empty()) {
      return true;
   } else if (a.nml.empty() != b.nml.empty()) {
      return false;
   } else {
      if (a.nml.size() != b.nml.size()) return false;
      for (size_t i = 0; i < a.nml.size(); ++i) {
         if (a.nml[i] != b.nml[i]) return false;
      }
      return true;
   }
}

bool abstractSonda_eq(const abstractSonda_t& a, const abstractSonda_t& b) {
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

   auto check_ptr = [](const std::vector<SomeType>& a_field, const std::vector<SomeType>& b_field) -> bool {
       if (a_field.empty() && b_field.empty()) return true;
       if (a_field.empty() != b_field.empty()) return false;
       if (a_field.size() != b_field.size()) return false;
       for (size_t i = 0; i < a_field.size(); ++i) {
           if (a_field[i] != b_field[i]) return false;
       }
       return true;
   };

   // Note: The specific types for FarField, Electric, etc. inside abstractSonda_t
   // are assumed to be comparable vectors. Adjust types as necessary based on definition.
   if (!check_ptr(a.FarField, b.FarField)) return false;
   if (!check_ptr(a.Electric, b.Electric)) return false;
   if (!check_ptr(a.Magnetic, b.Magnetic)) return false;
   if (!check_ptr(a.NormalElectric, b.NormalElectric)) return false;
   if (!check_ptr(a.NormalMagnetic, b.NormalMagnetic)) return false;
   if (!check_ptr(a.SurfaceElectricCurrent, b.SurfaceElectricCurrent)) return false;
   if (!check_ptr(a.SurfaceMagneticCurrent, b.SurfaceMagneticCurrent)) return false;

   return true;
}

bool sondas_eq(const Sondas_t& a, const Sondas_t& b) {
   if (a.n_probes != b.n_probes) return false;
   if (a.n_probes_max != b.n_probes_max) return false;

   if (a.probes.empty() && b.probes.empty()) {
      return true;
   } else if (a.probes.empty() != b.probes.empty()) {
      return false;
   } else {
      if (a.probes.size() != b.probes.size()) return false;
      for (size_t i = 0; i < a.probes.size(); ++i) {
         if (a.probes[i] != b.probes[i]) return false;
      }
      return true;
   }
}

bool sonda_eq(const Sonda_t& a, const Sonda_t& b) {
   if (a.grname != b.grname) return false;

   auto check_int_vec = [](const std::vector<int>& a_vec, const std::vector<int>& b_vec) -> bool {
       if (a_vec.empty() && b_vec.empty()) return true;
       if (a_vec.empty() != b_vec.empty()) return false;
       if (a_vec.size() != b_vec.size()) return false;
       for (size_t i = 0; i < a_vec.size(); ++i) {
           if (a_vec[i] != b_vec[i]) return false;
       }
       return true;
   };

   if (!check_int_vec(a.i, b.i)) return false;
   if (!check_int_vec(a.j, b.j)) return false;
   if (!check_int_vec(a.k, b.k)) return false;
   if (!check_int_vec(a.node, b.node)) return false;

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

bool masSonda_eq(const MasSonda_t& a, const MasSonda_t& b) {
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

   if (a.cordinates.empty() && b.cordinates.empty()) {
      return true;
   } else if (a.cordinates.empty() != b.cordinates.empty()) {
      return false;
   } else {
      if (a.cordinates.size() != b.cordinates.size()) return false;
      for (size_t i = 0; i < a.cordinates.size(); ++i) {
         if (a.cordinates[i] != b.cordinates[i]) return false;
      }
      return true;
   }
}

bool MasSondas_eq(const MasSondas_t& a, const MasSondas_t& b) {
   if (a.length != b.length) return false;
   if (a.length_max != b.length_max) return false;
   if (a.len_cor_max != b.len_cor_max) return false;

   if (a.collection.empty() && b.collection.empty()) {
      return true;
   } else if (a.collection.empty() != b.collection.empty()) {
      return false;
   } else {
      if (a.collection.size() != b.collection.size()) return false;
      for (size_t i = 0; i < a.collection.size(); ++i) {
         if (a.collection[i] != b.collection[i]) return false;
      }
      return true;
   }
}

bool bloqueprobe_eq(const BloqueProbe_t& a, const BloqueProbe_t& b) {
   if (a.tstart != b.tstart) return false;
   if (a.tstop != b.tstop) return false;
   if (a.tstep != b.tstep) return false;
   if (a.fstart != b.fstart) return false;
   if (a.fstop != b.fstop) return false;
   if (a.fstep != b.fstep) return false;
   if (a.FileNormalize != b.FileNormalize) return false;
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

bool bloqueprobes_eq(const BloqueProbes_t& a, const BloqueProbes_t& b) {
   if (a.n_bp != b.n_bp) return false;
   if (a.n_bp_max != b.n_bp_max) return false;

   if (a.bp.empty() && b.bp.empty()) {
      return true;
   } else if (a.bp.empty() != b.bp.empty()) {
      return false;
   } else {
      if (a.bp.size() != b.bp.size()) return false;
      for (size_t i = 0; i < a.bp.size(); ++i) {
         if (a.bp[i] != b.bp[i]) return false;
      }
      return true;
   }
}

bool volprobe_eq(const VolProbe_t& a, const VolProbe_t& b) {
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

   if (a.cordinates.empty() && b.cordinates.empty()) {
      return true;
   } else if (a.cordinates.empty() != b.cordinates.empty()) {
      return false;
   } else {
      if (a.cordinates.size() != b.cordinates.size()) return false;
      for (size_t i = 0; i < a.cordinates.size(); ++i) {
         if (a.cordinates[i] != b.cordinates[i]) return false;
      }
      return true;
   }
}

bool volprobes_eq(const VolProbes_t& a, const VolProbes_t& b) {
   if (a.length != b.length) return false;
   if (a.length_max != b.length_max) return false;
   if (a.len_cor_max != b.len_cor_max) return false;

   if (a.collection.empty() && b.collection.empty()) {
      return true;
   } else if (a.collection.empty() != b.collection.empty()) {
      return false;
   } else {
      if (a.collection.size() != b.collection.size()) return false;
      for (size_t i = 0; i < a.collection.size(); ++i) {
         if (a.collection[i] != b.collection[i]) return false;
      }
      return true;
   }
}