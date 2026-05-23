#ifndef NFDE_TYPES_H
#define NFDE_TYPES_H

#include <vector>
#include <string>
#include <complex>
#include <cstdint>
#include <utility>

// Assuming FDETYPES_m provides RKIND
// #include "FDETYPES_m.hpp" 

// Assuming conformal_types_m provides triangle_t, interval_t, bufsize
// #include "conformal_types_m.hpp"

// Assuming mtln_types_m is used if CompileWithMTLN is defined
// #ifdef CompileWithMTLN
// #include "mtln_types_m.hpp"
// #endif


namespace NFDETypes_m {

    using RK = double;

    // CONSTANTS FOR THE PARSER
    constexpr double SIGMA_PEC = 1e19;
    constexpr double SIGMA_PMC = 1e19;
    constexpr double EPSILON_VACUUM = 8.854187817e-12;
    constexpr double MU_VACUUM = 1.2566370614e-6;

    // PROBES
    constexpr int32_t NP_T1_PLAIN = 0;
    constexpr int32_t NP_T1_AMBOS = 2;
    constexpr int32_t NP_T2_TIME = 0;
    constexpr int32_t NP_T2_FREQ = 1;
    constexpr int32_t NP_T2_TRANSFER = 2;
    constexpr int32_t NP_T2_TIMEFREQ = 3;
    constexpr int32_t NP_T2_TIMETRANSF = 4;
    constexpr int32_t NP_T2_FREQTRANSF = 5;
    constexpr int32_t NP_T2_TIMEFRECTRANSF = 6;
    constexpr int32_t NP_COR_EX = 0;
    constexpr int32_t NP_COR_EY = 1;
    constexpr int32_t NP_COR_EZ = 2;
    constexpr int32_t NP_COR_HX = 3;
    constexpr int32_t NP_COR_HY = 4;
    constexpr int32_t NP_COR_HZ = 5;
    constexpr int32_t NP_COR_WIRECURRENT = 6;
    constexpr int32_t NP_COR_DDP = 7;
    constexpr int32_t NP_COR_LINE = 8;
    constexpr int32_t NP_COR_CHARGE = 9;
    // Probe field type constants
    constexpr int32_t iEx = 1, iEy = 2, iEz = 3;
    constexpr int32_t iHx = 4, iHy = 5, iHz = 6;
    constexpr int32_t iMEC = 51, iMHC = 52;
    constexpr int32_t iCur = 53, iCurX = 54, iCurY = 55, iCurZ = 56;
    constexpr int32_t mapvtk = 57;
    constexpr int32_t iExC = 61, iEyC = 62, iEzC = 63;
    constexpr int32_t iHxC = 64, iHyC = 65, iHzC = 66;
    constexpr int32_t farfield = 67, lineIntegral = 68;
    constexpr int32_t centroide = 8, Nothing = 666;
    constexpr int32_t iJx = 10*iEx, iJy = 10*iEy, iJz = 10*iEz;
    constexpr int32_t iBloqueJx = 100*iEx, iBloqueJy = 100*iEy, iBloqueJz = 100*iEz;
    constexpr int32_t iBloqueMx = 100*iHx, iBloqueMy = 100*iHy, iBloqueMz = 100*iHz;
    constexpr bool BcELECT = true;
    constexpr bool BcMAGNE = false;

    // THIN WIRES
    constexpr int32_t MATERIAL_CONS = 0;
    constexpr int32_t MATERIAL_absorbing = 100;
    constexpr int32_t Parallel_CONS = 1;
    constexpr int32_t SERIES_CONS = 2;
    constexpr int32_t DISPERSIVE_CONS = 3;

    // BORDERS
    constexpr int32_t F_PEC = 1;
    constexpr int32_t F_PMC = 2;
    constexpr int32_t F_PER = 4;
    constexpr int32_t F_MUR = 7;
    constexpr int32_t F_PML = 9;
    constexpr int32_t F_XL = 1;
    constexpr int32_t F_XU = 2;
    constexpr int32_t F_YL = 3;
    constexpr int32_t F_YU = 4;
    constexpr int32_t F_ZL = 5;
    constexpr int32_t F_ZU = 6;
    constexpr int32_t F_TIMEFRECTRANSF = 0;

    // rlc y diodos
    constexpr int32_t inductor = 20;
    constexpr int32_t capacitor = 21;
    constexpr int32_t resistor = 22;
    constexpr int32_t diodo = 23;
    constexpr int32_t Dielectric = 24;
    constexpr int32_t PMLbody = 25;

    // TYPES

    // Basic cordinate type for two points and orientation
    struct coords_t {
        int32_t Xi = -1;
        int32_t Xe = -1;
        int32_t Yi = -1;
        int32_t Ye = -1;
        int32_t Zi = -1;
        int32_t Ze = -1;
        int32_t Xtrancos = 1;
        int32_t Ytrancos = 1;
        int32_t Ztrancos = 1;
        int32_t Or = 0; // field orientation
        std::string tag; // Assuming BUFSIZE is handled by std::string or fixed char array. 
                         // std::string is safer and more C++ idiomatic.
    bool operator==(const coords_t& other) const {
        if (!(Xi == other.Xi)) return false;
        if (!(Xe == other.Xe)) return false;
        if (!(Yi == other.Yi)) return false;
        if (!(Ye == other.Ye)) return false;
        if (!(Zi == other.Zi)) return false;
        if (!(Ze == other.Ze)) return false;
        if (!(Xtrancos == other.Xtrancos)) return false;
        if (!(Ytrancos == other.Ytrancos)) return false;
        if (!(Ztrancos == other.Ztrancos)) return false;
        if (!(Or == other.Or)) return false;
        if (!(tag == other.tag)) return false;
        return true;
    }

    };

    struct coords_scaled_t {
        int32_t Xi = -1;
        int32_t Xe = -1;
        int32_t Yi = -1;
        int32_t Ye = -1;
        int32_t Zi = -1;
        int32_t Ze = -1;
        double xc = 0.0;
        double yc = 0.0;
        double zc = 0.0;
        int32_t Or = 0; // field orientation nuevo 2015
        std::string tag;
    bool operator==(const coords_scaled_t& other) const {
        if (!(Xi == other.Xi)) return false;
        if (!(Xe == other.Xe)) return false;
        if (!(Yi == other.Yi)) return false;
        if (!(Ye == other.Ye)) return false;
        if (!(Zi == other.Zi)) return false;
        if (!(Ze == other.Ze)) return false;
        if (!(xc == other.xc)) return false;
        if (!(yc == other.yc)) return false;
        if (!(zc == other.zc)) return false;
        if (!(Or == other.Or)) return false;
        if (!(tag == other.tag)) return false;
        return true;
    }

    };

    // Basic constants for materials
    struct Material_t {
        double eps = 0.0;
        double mu = 0.0;
        double sigma = 0.0;
        double sigmam = 0.0;
        int32_t id = 0;
    bool operator==(const Material_t& other) const {
        if (!(eps == other.eps)) return false;
        if (!(mu == other.mu)) return false;
        if (!(sigma == other.sigma)) return false;
        if (!(sigmam == other.sigmam)) return false;
        if (!(id == other.id)) return false;
        return true;
    }

    };

    // New Class which is a collection of different materials
    struct Materials_t {
        int32_t n_Mats = 0;
        int32_t n_Mats_max = 0;
        std::vector<Material_t> Mats; // Converted from pointer to vector for safety and ease
                                      // Note: Original was pointer, managed manually. 
                                      // Vector handles memory. If manual management is needed, 
                                      // use std::unique_ptr<std::vector<Material_t>> or raw pointer.
                                      // Given the complexity of Fortran pointer management, 
                                      // a vector is the standard C++ replacement for allocatable arrays.
                                      // However, to strictly mimic "pointer" behavior if needed:
                                      // std::vector<Material_t>* Mats = nullptr;
                                      // But std::vector is preferred.
    bool operator==(const Materials_t& other) const {
        if (!(n_Mats == other.n_Mats)) return false;
        if (!(n_Mats_max == other.n_Mats_max)) return false;
        if (!(Mats == other.Mats)) return false;
        return true;
    }

    };

    // Identifies conformal PEC "media"
    struct ConformalPECElements_t {
        std::vector<std::vector<double>> triangles; // triangle_t = std::vector<double>
        std::vector<std::pair<double, double>> intervals; // interval_t = std::pair<double, double>
        std::string tag;
    };

    struct ConformalPECRegions_t {
        std::vector<ConformalPECElements_t> volumes;
        std::vector<ConformalPECElements_t> surfaces;
    };

    struct edge_t {
        int32_t cell[3];
        int32_t direction = -1;
        double ratio = -1;
        double material_coords[2];
    };

    struct face_t {
        int32_t cell[3];
        int32_t direction = -1;
        double ratio = -1;
    };

    struct conformal_edge_media_t {
        std::vector<edge_t> edges;
        double ratio;
        int32_t n_elements;
    };

    struct conformal_face_media_t {
        std::vector<face_t> faces;
        double ratio;
        int32_t n_elements;
    };

    struct ConformalMedia_t {
        int32_t n_edges_media = 0;
        int32_t n_faces_media = 0;
        std::vector<conformal_face_media_t> face_media;
        std::vector<conformal_edge_media_t> edge_media;
        double time_step_scale_factor = 1.0;
        std::string tag;
    };

    // Locates all the different PEC media found
    struct PECRegions_t {
        int32_t nVols = 0;
        int32_t nSurfs = 0;
        int32_t nLins = 0;
        int32_t nVols_max = 0;
        int32_t nSurfs_max = 0;
        int32_t nLins_max = 0;
        std::vector<coords_t> Vols;
        std::vector<coords_t> Surfs;
        std::vector<coords_t> Lins;
    bool operator==(const PECRegions_t& other) const {
        if (!(nVols == other.nVols)) return false;
        if (!(nSurfs == other.nSurfs)) return false;
        if (!(nLins == other.nLins)) return false;
        if (!(nVols_max == other.nVols_max)) return false;
        if (!(nSurfs_max == other.nSurfs_max)) return false;
        if (!(nLins_max == other.nLins_max)) return false;
        if (!(Vols == other.Vols)) return false;
        if (!(Surfs == other.Surfs)) return false;
        if (!(Lins == other.Lins)) return false;
        return true;
    }

    };

    // Defines a Non Metal Body
    struct Dielectric_t {
        std::vector<coords_t> c1P;
        std::vector<coords_t> c2P;
        double sigma = 0.0;
        double eps = 0.0;
        double mu = 0.0;
        double sigmam = 0.0;
        int32_t n_C1P = 0;
        int32_t n_C2P = 0;
        // Lumped vars.
        double Rtime_on = 0.0;
        double Rtime_off = 0.0;
        double R = 0.0;
        double L = 0.0;
        double C = 0.0;
        // Stoch vars.
        double R_devia = 0.0;
        double L_devia = 0.0;
        double C_devia = 0.0;
        //
        double DiodB = 0.0;
        double DiodIsat = 0.0;
        int32_t DiodOri = 0;
        // Berenger's waveports
        int32_t orient = 0;
        bool resistor = false;
        bool inductor = false;
        bool capacitor = false;
        bool diodo = false;
        bool plain = false;
        bool PMLbody = false;
    bool operator==(const Dielectric_t& other) const {
        if (!(c1P == other.c1P)) return false;
        if (!(c2P == other.c2P)) return false;
        if (!(sigma == other.sigma)) return false;
        if (!(eps == other.eps)) return false;
        if (!(mu == other.mu)) return false;
        if (!(sigmam == other.sigmam)) return false;
        if (!(n_C1P == other.n_C1P)) return false;
        if (!(n_C2P == other.n_C2P)) return false;
        if (!(Rtime_on == other.Rtime_on)) return false;
        if (!(Rtime_off == other.Rtime_off)) return false;
        if (!(R == other.R)) return false;
        if (!(L == other.L)) return false;
        if (!(C == other.C)) return false;
        if (!(R_devia == other.R_devia)) return false;
        if (!(L_devia == other.L_devia)) return false;
        if (!(C_devia == other.C_devia)) return false;
        if (!(DiodB == other.DiodB)) return false;
        if (!(DiodIsat == other.DiodIsat)) return false;
        if (!(DiodOri == other.DiodOri)) return false;
        if (!(orient == other.orient)) return false;
        if (!(resistor == other.resistor)) return false;
        if (!(inductor == other.inductor)) return false;
        if (!(capacitor == other.capacitor)) return false;
        if (!(diodo == other.diodo)) return false;
        if (!(plain == other.plain)) return false;
        if (!(PMLbody == other.PMLbody)) return false;
        return true;
    }

    };

    // Locates all the different Non Metal Media found
    struct DielectricRegions_t {
        std::vector<Dielectric_t> Vols;
        std::vector<Dielectric_t> Surfs;
        std::vector<Dielectric_t> Lins;
        int32_t nVols = 0;
        int32_t nSurfs = 0;
        int32_t nLins = 0;
        int32_t nVols_max = 0;
        int32_t nSurfs_max = 0;
        int32_t nLins_max = 0;
        int32_t n_C1P_max = 0;
        int32_t n_C2P_max = 0;
    bool operator==(const DielectricRegions_t& other) const {
        if (!(Vols == other.Vols)) return false;
        if (!(Surfs == other.Surfs)) return false;
        if (!(Lins == other.Lins)) return false;
        if (!(nVols == other.nVols)) return false;
        if (!(nSurfs == other.nSurfs)) return false;
        if (!(nLins == other.nLins)) return false;
        if (!(nVols_max == other.nVols_max)) return false;
        if (!(nSurfs_max == other.nSurfs_max)) return false;
        if (!(nLins_max == other.nLins_max)) return false;
        if (!(n_C1P_max == other.n_C1P_max)) return false;
        if (!(n_C2P_max == other.n_C2P_max)) return false;
        return true;
    }

    };

    // type that defines the information of a frequency dependent material
    struct FreqDepenMaterial_t {
        std::vector<std::complex<double>> a11;
        std::vector<std::complex<double>> b11;
        std::vector<std::complex<double>> am11;
        std::vector<std::complex<double>> bm11;
        std::vector<std::complex<double>> a12;
        std::vector<std::complex<double>> b12;
        std::vector<std::complex<double>> am12;
        std::vector<std::complex<double>> bm12;
        std::vector<std::complex<double>> a13;
        std::vector<std::complex<double>> b13;
        std::vector<std::complex<double>> am13;
        std::vector<std::complex<double>> bm13;
        std::vector<std::complex<double>> a22;
        std::vector<std::complex<double>> b22;
        std::vector<std::complex<double>> am22;
        std::vector<std::complex<double>> bm22;
        std::vector<std::complex<double>> a23;
        std::vector<std::complex<double>> b23;
        std::vector<std::complex<double>> am23;
        std::vector<std::complex<double>> bm23;
        std::vector<std::complex<double>> a33;
        std::vector<std::complex<double>> b33;
        std::vector<std::complex<double>> am33;
        std::vector<std::complex<double>> bm33;
        std::vector<double> alpha;
        std::vector<double> beta;
        std::vector<double> gamma;
        std::vector<double> alpham;
        std::vector<double> betam;
        std::vector<double> gammam;
        std::vector<coords_t> c;
        double eps11 = 0.0;
        double eps12 = 0.0;
        double eps13 = 0.0;
        double eps22 = 0.0;
        double eps23 = 0.0;
        double eps33 = 0.0;
        double mu11 = 0.0;
        double mu12 = 0.0;
        double mu13 = 0.0;
        double mu22 = 0.0;
        double mu23 = 0.0;
        double mu33 = 0.0;
        double sigma11 = 0.0;
        double sigma12 = 0.0;
        double sigma13 = 0.0;
        double sigma22 = 0.0;
        double sigma23 = 0.0;
        double sigma33 = 0.0;
        double sigmam11 = 0.0;
        double sigmam12 = 0.0;
        double sigmam13 = 0.0;
        double sigmam22 = 0.0;
        double sigmam23 = 0.0;
        double sigmam33 = 0.0;
        int32_t K11 = 0;
        int32_t Km11 = 0;
        int32_t K12 = 0;
        int32_t Km12 = 0;
        int32_t K13 = 0;
        int32_t Km13 = 0;
        int32_t K22 = 0;
        int32_t Km22 = 0;
        int32_t K23 = 0;
        int32_t Km23 = 0;
        int32_t K33 = 0;
        int32_t Km33 = 0;
        int32_t L = 0;
        int32_t Lm = 0;
        int32_t n_c = 0;
        std::string files = " ";
    };

    // type that defines the list of frequency dependent materials
    struct FreqDepenMaterials_t {
        std::vector<FreqDepenMaterial_t> Vols;
        std::vector<FreqDepenMaterial_t> Surfs;
        std::vector<FreqDepenMaterial_t> Lins;
        int32_t nVols = 0;
        int32_t nSurfs = 0;
        int32_t nLins = 0;
        int32_t nVols_max = 0;
        int32_t nSurfs_max = 0;
        int32_t nLins_max = 0;
        int32_t n_c_max = 0;
    };

    // Type for the ANISOTROPIC body, surface and lines
    struct ANISOTROPICbody_t {
        std::vector<coords_t> c1P;
        std::vector<coords_t> c2P;
        double sigma[3][3];
        double eps[3][3];
        double mu[3][3];
        double sigmam[3][3];
        int32_t n_C1P = 0;
        int32_t n_C2P = 0;
    };

    // Type that contains the elements found in the nfde File
    struct ANISOTROPICelements_t {
        std::vector<ANISOTROPICbody_t> Vols;
        std::vector<ANISOTROPICbody_t> Surfs;
        std::vector<ANISOTROPICbody_t> Lins;
        int32_t nVols = 0;
        int32_t nSurfs = 0;
        int32_t nLins = 0;
        int32_t nVols_max = 0;
        int32_t nSurfs_max = 0;
        int32_t nLins_max = 0;
        int32_t n_C1P_max = 0;
        int32_t n_C2P_max = 0;
    };

    // Defines a Comp Surface
    struct LossyThinSurface_t {
        std::vector<coords_t> c;
        std::vector<double> sigma;
        std::vector<double> eps;
        std::vector<double> mu;
        std::vector<double> sigmam;
        std::vector<double> thk;
        // for_devia
        std::vector<double> sigma_devia;
        std::vector<double> eps_devia;
        std::vector<double> mu_devia;
        std::vector<double> sigmam_devia;
        std::vector<double> thk_devia;
        //
        int32_t nc = 0;
        std::string files = " ";
        int32_t numcapas;
    bool operator==(const LossyThinSurface_t& other) const {
        if (!(c == other.c)) return false;
        if (!(sigma == other.sigma)) return false;
        if (!(eps == other.eps)) return false;
        if (!(mu == other.mu)) return false;
        if (!(sigmam == other.sigmam)) return false;
        if (!(thk == other.thk)) return false;
        if (!(sigma_devia == other.sigma_devia)) return false;
        if (!(eps_devia == other.eps_devia)) return false;
        if (!(mu_devia == other.mu_devia)) return false;
        if (!(sigmam_devia == other.sigmam_devia)) return false;
        if (!(thk_devia == other.thk_devia)) return false;
        if (!(nc == other.nc)) return false;
        if (!(files == other.files)) return false;
        if (!(numcapas == other.numcapas)) return false;
        return true;
    }

    };

    // Locates all the different Comp media found
    struct LossyThinSurfaces_t {
        std::vector<LossyThinSurface_t> cs;
        int32_t length = 0;
        int32_t length_max = 0;
        int32_t nC_max = 0;
    bool operator==(const LossyThinSurfaces_t& other) const {
        if (!(cs == other.cs)) return false;
        if (!(length == other.length)) return false;
        if (!(length_max == other.length_max)) return false;
        if (!(nC_max == other.nC_max)) return false;
        return true;
    }

    };

    // Component for Thin Wires
    struct ThinWireComp_t {
        std::string srctype;
        std::string srcfile;
        int32_t i = -1;
        int32_t j = -1;
        int32_t K = -1;
        int32_t nd = -1;
        int32_t d = -1;
        double m = 0.0;
        std::string tag;
    bool operator==(const ThinWireComp_t& other) const {
        if (!(srctype == other.srctype)) return false;
        if (!(srcfile == other.srcfile)) return false;
        if (!(i == other.i)) return false;
        if (!(j == other.j)) return false;
        if (!(K == other.K)) return false;
        if (!(nd == other.nd)) return false;
        if (!(d == other.d)) return false;
        if (!(m == other.m)) return false;
        if (!(tag == other.tag)) return false;
        return true;
    }

    };

    // ThinWire component that defines the overall properties of the definition of ThinWires
    struct ThinWire_t {
        std::vector<ThinWireComp_t> twc;
        double rad = 0;
        double rad_devia = 0;
        bool disp = false;
        std::string dispfile;
        double res = 0;
        double res_devia = 0;
        double ind = 0;
        double ind_devia = 0;
        double cap = 0;
        double cap_devia = 0;
        double P_res = 0;
        double P_ind = 0;
        double P_cap = 0;
        std::string dispfile_LeftEnd;
        double R_LeftEnd = 0;
        double R_LeftEnd_devia = 0;
        double L_LeftEnd = 0;
        double L_LeftEnd_devia = 0;
        double C_LeftEnd = 0;
        double C_LeftEnd_devia = 0;
        std::string dispfile_RightEnd;
        double R_RightEnd = 0;
        double R_RightEnd_devia = 0;
        double L_RightEnd = 0;
        double L_RightEnd_devia = 0;
        double C_RightEnd = 0;
        double C_RightEnd_devia = 0;
        int32_t LeftEnd = 0;
        int32_t RightEnd = 0;
        // Components
        int32_t tl = 0;
        int32_t tr = 0;
        int32_t n_twc = 0;
        int32_t n_twc_max = 0;
    bool operator==(const ThinWire_t& other) const {
        if (!(twc == other.twc)) return false;
        if (!(rad == other.rad)) return false;
        if (!(rad_devia == other.rad_devia)) return false;
        if (!(disp == other.disp)) return false;
        if (!(dispfile == other.dispfile)) return false;
        if (!(res == other.res)) return false;
        if (!(res_devia == other.res_devia)) return false;
        if (!(ind == other.ind)) return false;
        if (!(ind_devia == other.ind_devia)) return false;
        if (!(cap == other.cap)) return false;
        if (!(cap_devia == other.cap_devia)) return false;
        if (!(P_res == other.P_res)) return false;
        if (!(P_ind == other.P_ind)) return false;
        if (!(P_cap == other.P_cap)) return false;
        if (!(dispfile_LeftEnd == other.dispfile_LeftEnd)) return false;
        if (!(R_LeftEnd == other.R_LeftEnd)) return false;
        if (!(R_LeftEnd_devia == other.R_LeftEnd_devia)) return false;
        if (!(L_LeftEnd == other.L_LeftEnd)) return false;
        if (!(L_LeftEnd_devia == other.L_LeftEnd_devia)) return false;
        if (!(C_LeftEnd == other.C_LeftEnd)) return false;
        if (!(C_LeftEnd_devia == other.C_LeftEnd_devia)) return false;
        if (!(dispfile_RightEnd == other.dispfile_RightEnd)) return false;
        if (!(R_RightEnd == other.R_RightEnd)) return false;
        if (!(R_RightEnd_devia == other.R_RightEnd_devia)) return false;
        if (!(L_RightEnd == other.L_RightEnd)) return false;
        if (!(L_RightEnd_devia == other.L_RightEnd_devia)) return false;
        if (!(C_RightEnd == other.C_RightEnd)) return false;
        if (!(C_RightEnd_devia == other.C_RightEnd_devia)) return false;
        if (!(LeftEnd == other.LeftEnd)) return false;
        if (!(RightEnd == other.RightEnd)) return false;
        if (!(tl == other.tl)) return false;
        if (!(tr == other.tr)) return false;
        if (!(n_twc == other.n_twc)) return false;
        if (!(n_twc_max == other.n_twc_max)) return false;
        return true;
    }

    };

    // List of the different thin wires that were found in the file
    struct ThinWires_t {
        std::vector<ThinWire_t> tw;
        int32_t n_tw = 0;
        int32_t n_tw_max = 0;
    bool operator==(const ThinWires_t& other) const {
        if (!(tw == other.tw)) return false;
        if (!(n_tw == other.n_tw)) return false;
        if (!(n_tw_max == other.n_tw_max)) return false;
        return true;
    }

    };

    // Component for Slanted Wires
    struct SlantedWireComp_t {
        std::string srctype;
        std::string srcfile;
        double x = -1.0;
        double y = -1.0;
        double z = -1.0;
        int32_t nd = -1;
        double m = 0.0;
        std::string tag;
    bool operator==(const SlantedWireComp_t& other) const {
        if (!(srctype == other.srctype)) return false;
        if (!(srcfile == other.srcfile)) return false;
        if (!(x == other.x)) return false;
        if (!(y == other.y)) return false;
        if (!(z == other.z)) return false;
        if (!(nd == other.nd)) return false;
        if (!(m == other.m)) return false;
        if (!(tag == other.tag)) return false;
        return true;
    }

    };

    // ThinWire component that defines the overall properties of the definition of ThinWires (Slanted)
    struct SlantedWire_t {
        std::vector<SlantedWireComp_t> swc;
        double rad = 0;
        bool disp = false;
        std::string dispfile;
        double res = 0;
        double ind = 0;
        double cap = 0;
        double P_res = 0;
        double P_ind = 0;
        double P_cap = 0;
        std::string dispfile_LeftEnd;
        double R_LeftEnd = 0;
        double L_LeftEnd = 0;
        double C_LeftEnd = 0;
        std::string dispfile_RightEnd;
        double R_RightEnd = 0;
        double L_RightEnd = 0;
        double C_RightEnd = 0;
        int32_t LeftEnd = 0;
        int32_t RightEnd = 0;
        // Components
        int32_t tl = 0;
        int32_t tr = 0;
        int32_t n_swc = 0;
        int32_t n_swc_max = 0;
    bool operator==(const SlantedWire_t& other) const {
        if (!(swc == other.swc)) return false;
        if (!(rad == other.rad)) return false;
        if (!(disp == other.disp)) return false;
        if (!(dispfile == other.dispfile)) return false;
        if (!(res == other.res)) return false;
        if (!(ind == other.ind)) return false;
        if (!(cap == other.cap)) return false;
        if (!(P_res == other.P_res)) return false;
        if (!(P_ind == other.P_ind)) return false;
        if (!(P_cap == other.P_cap)) return false;
        if (!(dispfile_LeftEnd == other.dispfile_LeftEnd)) return false;
        if (!(R_LeftEnd == other.R_LeftEnd)) return false;
        if (!(L_LeftEnd == other.L_LeftEnd)) return false;
        if (!(C_LeftEnd == other.C_LeftEnd)) return false;
        if (!(dispfile_RightEnd == other.dispfile_RightEnd)) return false;
        if (!(R_RightEnd == other.R_RightEnd)) return false;
        if (!(L_RightEnd == other.L_RightEnd)) return false;
        if (!(C_RightEnd == other.C_RightEnd)) return false;
        if (!(LeftEnd == other.LeftEnd)) return false;
        if (!(RightEnd == other.RightEnd)) return false;
        if (!(tl == other.tl)) return false;
        if (!(tr == other.tr)) return false;
        if (!(n_swc == other.n_swc)) return false;
        if (!(n_swc_max == other.n_swc_max)) return false;
        return true;
    }

    };

    // List of the different thin wires that were found in the file (Slanted)
    struct SlantedWiresInfo_t {
        std::vector<SlantedWire_t> sw;
        int32_t n_sw = 0;
        int32_t n_sw_max = 0;
    bool operator==(const SlantedWiresInfo_t& other) const {
        if (!(sw == other.sw)) return false;
        if (!(n_sw == other.n_sw)) return false;
        if (!(n_sw_max == other.n_sw_max)) return false;
        return true;
    }

    };

    // Component for Thin Slots
    struct ThinSlotComp_t {
        int32_t i = 0;
        int32_t j = 0;
        int32_t K = 0;
        int32_t node = 0;
        int32_t dir = -1;
        int32_t Or = -1;
        std::string tag;
    bool operator==(const ThinSlotComp_t& other) const {
        if (!(i == other.i)) return false;
        if (!(j == other.j)) return false;
        if (!(K == other.K)) return false;
        if (!(node == other.node)) return false;
        if (!(dir == other.dir)) return false;
        if (!(Or == other.Or)) return false;
        if (!(tag == other.tag)) return false;
        return true;
    }

    };

    // ThinSlot component that defines the overall properties of the definition of ThinSlots in ORIGINAL
    struct ThinSlot_t {
        std::vector<ThinSlotComp_t> tgc;
        double width = 0;
        int32_t n_tgc = 0;
        int32_t n_tgc_max = 0;
    bool operator==(const ThinSlot_t& other) const {
        if (!(tgc == other.tgc)) return false;
        if (!(width == other.width)) return false;
        if (!(n_tgc == other.n_tgc)) return false;
        if (!(n_tgc_max == other.n_tgc_max)) return false;
        return true;
    }

    };

    // List of the different thin Slots that were found in the file
    struct ThinSlots_t {
        std::vector<ThinSlot_t> tg;
        int32_t n_tg = 0;
        int32_t n_tg_max = 0;
    bool operator==(const ThinSlots_t& other) const {
        if (!(tg == other.tg)) return false;
        if (!(n_tg == other.n_tg)) return false;
        if (!(n_tg_max == other.n_tg_max)) return false;
        return true;
    }

    };

    // PML Border Type
    struct FronteraPML_t {
        double orden = 2.0;
        double refl = 1e-3;
        int32_t numCapas = 8;
    bool operator==(const FronteraPML_t& other) const {
        if (!(orden == other.orden)) return false;
        if (!(refl == other.refl)) return false;
        if (!(numCapas == other.numCapas)) return false;
        return true;
    }

    };

    // Tipo de la frontera
    struct Frontera_t {
        int32_t tipoFrontera[6];
        FronteraPML_t propiedadesPML[6];
    bool operator==(const Frontera_t& other) const {
        for (int _i = 0; _i < 6; ++_i) {
            if (tipoFrontera[_i] != other.tipoFrontera[_i]) return false;
            if (!(propiedadesPML[_i] == other.propiedadesPML[_i])) return false;
        }
        return true;
    }

    };

    // type to define the new probe object
    struct MasSonda_t {
        std::string filename;
        std::vector<coords_t> cordinates;
        double tstart;
        double tstop;
        double tstep;
        double fstart;
        double fstop;
        double fstep;
        int32_t type1;
        int32_t type2;
        int32_t len_cor = 0;
        std::string outputrequest;
    bool operator==(const MasSonda_t& other) const {
        if (!(filename == other.filename)) return false;
        if (!(cordinates == other.cordinates)) return false;
        if (!(tstart == other.tstart)) return false;
        if (!(tstop == other.tstop)) return false;
        if (!(tstep == other.tstep)) return false;
        if (!(fstart == other.fstart)) return false;
        if (!(fstop == other.fstop)) return false;
        if (!(fstep == other.fstep)) return false;
        if (!(type1 == other.type1)) return false;
        if (!(type2 == other.type2)) return false;
        if (!(len_cor == other.len_cor)) return false;
        if (!(outputrequest == other.outputrequest)) return false;
        return true;
    }

    };

    // type that defines a list of probes to be appended and accessed
    struct MasSondas_t {
        std::vector<MasSonda_t> collection;
        int32_t length = 0;
        int32_t length_max = 0;
        int32_t len_cor_max = 0;
    bool operator==(const MasSondas_t& other) const {
        if (!(collection == other.collection)) return false;
        if (!(length == other.length)) return false;
        if (!(length_max == other.length_max)) return false;
        if (!(len_cor_max == other.len_cor_max)) return false;
        return true;
    }

    };

    // This type contains the basic information in nearly all the different PROBES
    struct Sonda_t {
        std::string grname;
        std::vector<int32_t> i;
        std::vector<int32_t> j;
        std::vector<int32_t> K;
        std::vector<int32_t> node;
        int32_t n_cord = 0;
        int32_t n_cord_max = 0;
        RK tstart = 0.0;
        RK tstop = 0.0;
        RK tstep = 0.0;
        std::string outputrequest;
        RK fstart = 0.0;
        RK fstop = 0.0;
        RK fstep = 0.0;
        RK phistart = 0.0;
        RK phistop = 0.0;
        RK phistep = 0.0;
        RK thetastart = 0.0;
        RK thetastop = 0.0;
        RK thetastep = 0.0;
        std::string FileNormalize;
    bool operator==(const Sonda_t& other) const {
        if (!(grname == other.grname)) return false;
        if (!(i == other.i)) return false;
        if (!(j == other.j)) return false;
        if (!(K == other.K)) return false;
        if (!(node == other.node)) return false;
        if (!(n_cord == other.n_cord)) return false;
        if (!(n_cord_max == other.n_cord_max)) return false;
        if (!(tstart == other.tstart)) return false;
        if (!(tstop == other.tstop)) return false;
        if (!(tstep == other.tstep)) return false;
        if (!(outputrequest == other.outputrequest)) return false;
        if (!(fstart == other.fstart)) return false;
        if (!(fstop == other.fstop)) return false;
        if (!(fstep == other.fstep)) return false;
        if (!(phistart == other.phistart)) return false;
        if (!(phistop == other.phistop)) return false;
        if (!(phistep == other.phistep)) return false;
        if (!(thetastart == other.thetastart)) return false;
        if (!(thetastop == other.thetastop)) return false;
        if (!(thetastep == other.thetastep)) return false;
        if (!(FileNormalize == other.FileNormalize)) return false;
        return true;
    }

    };

template<typename T>
struct FortranPointer {
    T* ptr = nullptr;
    FortranPointer() : ptr(nullptr) {}
    FortranPointer(T* p) : ptr(p) {}
    operator T*() const { return ptr; }
    T* operator->() const { return ptr; }
    T& operator*() const { return *ptr; }
};

// Type: FarField_Sonda_t
struct FarField_Sonda_t {
    Sonda_t probe;
    bool operator==(const FarField_Sonda_t& other) const {
        if (!(probe == other.probe)) return false;
        return true;
    }

};

// Type: Electric_Sonda_t
struct Electric_Sonda_t {
    Sonda_t probe;
    bool operator==(const Electric_Sonda_t& other) const {
        if (!(probe == other.probe)) return false;
        return true;
    }

};

// Type: Magnetic_Sonda_t
struct Magnetic_Sonda_t {
    Sonda_t probe;
    bool operator==(const Magnetic_Sonda_t& other) const {
        if (!(probe == other.probe)) return false;
        return true;
    }

};

// Type: NormalElectric_Sonda_t
struct NormalElectric_Sonda_t {
    Sonda_t probe;
    std::vector<int32_t> nml; // dimension(:), pointer
    int32_t n_nml = 0;
    int32_t n_nml_max = 0;
    bool operator==(const NormalElectric_Sonda_t& other) const {
        if (!(probe == other.probe)) return false;
        if (!(nml == other.nml)) return false;
        if (!(n_nml == other.n_nml)) return false;
        if (!(n_nml_max == other.n_nml_max)) return false;
        return true;
    }

};

// Type: NormalMagnetic_Sonda_t
struct NormalMagnetic_Sonda_t {
    Sonda_t probe;
    std::vector<int32_t> nml; // dimension(:), pointer
    int32_t n_nml = 0;
    int32_t n_nml_max = 0;
    bool operator==(const NormalMagnetic_Sonda_t& other) const {
        if (!(probe == other.probe)) return false;
        if (!(nml == other.nml)) return false;
        if (!(n_nml == other.n_nml)) return false;
        if (!(n_nml_max == other.n_nml_max)) return false;
        return true;
    }

};

// Type: SurfaceElectricCurrent_Sonda_t
struct SurfaceElectricCurrent_Sonda_t {
    Sonda_t probe;
    std::vector<int32_t> nml; // dimension(:), pointer
    int32_t n_nml = 0;
    int32_t n_nml_max = 0;
    bool operator==(const SurfaceElectricCurrent_Sonda_t& other) const {
        if (!(probe == other.probe)) return false;
        if (!(nml == other.nml)) return false;
        if (!(n_nml == other.n_nml)) return false;
        if (!(n_nml_max == other.n_nml_max)) return false;
        return true;
    }

};

// Type: SurfaceMagneticCurrent_Sonda_t
struct SurfaceMagneticCurrent_Sonda_t {
    Sonda_t probe;
    std::vector<int32_t> nml; // dimension(:), pointer
    int32_t n_nml = 0;
    int32_t n_nml_max = 0;
    bool operator==(const SurfaceMagneticCurrent_Sonda_t& other) const {
        if (!(probe == other.probe)) return false;
        if (!(nml == other.nml)) return false;
        if (!(n_nml == other.n_nml)) return false;
        if (!(n_nml_max == other.n_nml_max)) return false;
        return true;
    }

};

// Type: abstractSonda_t
struct abstractSonda_t {
    int32_t n_FarField = 0;
    int32_t n_Electric = 0;
    int32_t n_Magnetic = 0;
    int32_t n_NormalElectric = 0;
    int32_t n_NormalMagnetic = 0;
    int32_t n_SurfaceElectricCurrent = 0;
    int32_t n_SurfaceMagneticCurrent = 0;
    
    int32_t n_FarField_max = 0;
    int32_t n_Electric_max = 0;
    int32_t n_Magnetic_max = 0;
    int32_t n_NormalElectric_max = 0;
    int32_t n_NormalMagnetic_max = 0;
    int32_t n_SurfaceElectricCurrent_max = 0;
    int32_t n_SurfaceMagneticCurrent_max = 0;
    
    std::vector<FarField_Sonda_t> FarField; // dimension(:), pointer
    std::vector<Electric_Sonda_t> Electric; // dimension(:), pointer
    std::vector<Magnetic_Sonda_t> Magnetic; // dimension(:), pointer
    std::vector<NormalElectric_Sonda_t> NormalElectric; // dimension(:), pointer
    std::vector<NormalMagnetic_Sonda_t> NormalMagnetic; // dimension(:), pointer
    std::vector<SurfaceElectricCurrent_Sonda_t> SurfaceElectricCurrent; // dimension(:), pointer
    std::vector<SurfaceMagneticCurrent_Sonda_t> SurfaceMagneticCurrent; // dimension(:), pointer
    bool operator==(const abstractSonda_t& other) const {
        if (!(n_FarField == other.n_FarField)) return false;
        if (!(n_Electric == other.n_Electric)) return false;
        if (!(n_Magnetic == other.n_Magnetic)) return false;
        if (!(n_NormalElectric == other.n_NormalElectric)) return false;
        if (!(n_NormalMagnetic == other.n_NormalMagnetic)) return false;
        if (!(n_SurfaceElectricCurrent == other.n_SurfaceElectricCurrent)) return false;
        if (!(n_SurfaceMagneticCurrent == other.n_SurfaceMagneticCurrent)) return false;
        if (!(n_FarField_max == other.n_FarField_max)) return false;
        if (!(n_Electric_max == other.n_Electric_max)) return false;
        if (!(n_Magnetic_max == other.n_Magnetic_max)) return false;
        if (!(n_NormalElectric_max == other.n_NormalElectric_max)) return false;
        if (!(n_NormalMagnetic_max == other.n_NormalMagnetic_max)) return false;
        if (!(n_SurfaceElectricCurrent_max == other.n_SurfaceElectricCurrent_max)) return false;
        if (!(n_SurfaceMagneticCurrent_max == other.n_SurfaceMagneticCurrent_max)) return false;
        if (!(FarField == other.FarField)) return false;
        if (!(Electric == other.Electric)) return false;
        if (!(Magnetic == other.Magnetic)) return false;
        if (!(NormalElectric == other.NormalElectric)) return false;
        if (!(NormalMagnetic == other.NormalMagnetic)) return false;
        if (!(SurfaceElectricCurrent == other.SurfaceElectricCurrent)) return false;
        if (!(SurfaceMagneticCurrent == other.SurfaceMagneticCurrent)) return false;
        return true;
    }

};

// Type: Sondas_t
struct Sondas_t {
    std::vector<abstractSonda_t> probes; // dimension(:), pointer
    int32_t n_probes = 0;
    int32_t n_probes_max = 0;
    bool operator==(const Sondas_t& other) const {
        if (!(probes == other.probes)) return false;
        if (!(n_probes == other.n_probes)) return false;
        if (!(n_probes_max == other.n_probes_max)) return false;
        return true;
    }

};

// Type: BloqueProbe_t
struct BloqueProbe_t {
    RK tstart = 0.0;
    RK tstop = 0.0;
    RK tstep = 0.0;
    RK fstart = 0.0;
    RK fstop = 0.0;
    RK fstep = 0.0;
    std::string FileNormalize; // len=BUFSIZE
    int32_t type2 = 0;
    int32_t i1 = 0;
    int32_t i2 = 0;
    int32_t j1 = 0;
    int32_t j2 = 0;
    int32_t k1 = 0;
    int32_t k2 = 0;
    int32_t skip = 0;
    int32_t nml = 0;
    bool t = false;
    std::string outputrequest; // len=BUFSIZE
    std::string tag; // len=BUFSIZE
    bool operator==(const BloqueProbe_t& other) const {
        if (!(tstart == other.tstart)) return false;
        if (!(tstop == other.tstop)) return false;
        if (!(tstep == other.tstep)) return false;
        if (!(fstart == other.fstart)) return false;
        if (!(fstop == other.fstop)) return false;
        if (!(fstep == other.fstep)) return false;
        if (!(FileNormalize == other.FileNormalize)) return false;
        if (!(type2 == other.type2)) return false;
        if (!(i1 == other.i1)) return false;
        if (!(i2 == other.i2)) return false;
        if (!(j1 == other.j1)) return false;
        if (!(j2 == other.j2)) return false;
        if (!(k1 == other.k1)) return false;
        if (!(k2 == other.k2)) return false;
        if (!(skip == other.skip)) return false;
        if (!(nml == other.nml)) return false;
        if (!(t == other.t)) return false;
        if (!(outputrequest == other.outputrequest)) return false;
        if (!(tag == other.tag)) return false;
        return true;
    }

};

// Type: BloqueProbes_t
struct BloqueProbes_t {
    std::vector<BloqueProbe_t> bp; // dimension(:), pointer
    int32_t n_bp = 0;
    int32_t n_bp_max = 0;
    bool operator==(const BloqueProbes_t& other) const {
        if (!(bp == other.bp)) return false;
        if (!(n_bp == other.n_bp)) return false;
        if (!(n_bp_max == other.n_bp_max)) return false;
        return true;
    }

};

// Type: VolProbe_t
struct VolProbe_t {
    std::vector<NFDETypes_m::coords_t> cordinates; // dimension(:), pointer
    RK tstart = 0.0;
    RK tstop = 0.0;
    RK tstep = 0.0;
    std::string outputrequest; // len=BUFSIZE
    int32_t len_cor = 0;
    RK fstart = 0.0;
    RK fstop = 0.0;
    RK fstep = 0.0;
    int32_t type2 = 0;
    std::string filename; // len=BUFSIZE
    bool operator==(const VolProbe_t& other) const {
        if (!(cordinates == other.cordinates)) return false;
        if (!(tstart == other.tstart)) return false;
        if (!(tstop == other.tstop)) return false;
        if (!(tstep == other.tstep)) return false;
        if (!(outputrequest == other.outputrequest)) return false;
        if (!(len_cor == other.len_cor)) return false;
        if (!(fstart == other.fstart)) return false;
        if (!(fstop == other.fstop)) return false;
        if (!(fstep == other.fstep)) return false;
        if (!(type2 == other.type2)) return false;
        if (!(filename == other.filename)) return false;
        return true;
    }

};

// Type: VolProbes_t
struct VolProbes_t {
    std::vector<VolProbe_t> collection; // dimension(:), pointer
    int32_t length = 0;
    int32_t length_max = 0;
    int32_t len_cor_max = 0;
    bool operator==(const VolProbes_t& other) const {
        if (!(collection == other.collection)) return false;
        if (!(length == other.length)) return false;
        if (!(length_max == other.length_max)) return false;
        if (!(len_cor_max == other.len_cor_max)) return false;
        return true;
    }

};

// Type: Box_t
struct Box_t {
    std::string nombre_fichero; // len=BUFSIZE
    int32_t coor1[3];
    int32_t coor2[3];
    bool operator==(const Box_t& other) const {
        if (nombre_fichero != other.nombre_fichero) return false;
        for (int _i = 0; _i < 3; ++_i) {
            if (coor1[_i] != other.coor1[_i]) return false;
            if (coor2[_i] != other.coor2[_i]) return false;
        }
        return true;
    }

};

// Type: Boxes_t
struct Boxes_t {
    std::vector<Box_t> Vols; // dimension(:), pointer
    int32_t nVols = 0;
    int32_t nVols_max = 0;
    bool operator==(const Boxes_t& other) const {
        if (!(Vols == other.Vols)) return false;
        if (!(nVols == other.nVols)) return false;
        if (!(nVols_max == other.nVols_max)) return false;
        return true;
    }

};

// Type: PlaneWave_t
struct PlaneWave_t {
    std::string nombre_fichero; // len=BUFSIZE
    std::string atributo; // len=BUFSIZE
    int32_t coor1[3];
    int32_t coor2[3];
    RK theta = 0.0;
    RK phi = 0.0;
    RK alpha = 0.0;
    RK beta = 0.0;
    bool isRC = false;
    RK INCERTMAX = 0.0;
    int32_t numModes = 0;
    bool operator==(const PlaneWave_t& other) const {
        if (nombre_fichero != other.nombre_fichero) return false;
        if (atributo != other.atributo) return false;
        for (int _i = 0; _i < 3; ++_i) {
            if (coor1[_i] != other.coor1[_i]) return false;
            if (coor2[_i] != other.coor2[_i]) return false;
        }
        if (!(theta == other.theta)) return false;
        if (!(phi == other.phi)) return false;
        if (!(alpha == other.alpha)) return false;
        if (!(beta == other.beta)) return false;
        if (!(isRC == other.isRC)) return false;
        if (!(INCERTMAX == other.INCERTMAX)) return false;
        if (!(numModes == other.numModes)) return false;
        return true;
    }

};

// Type: PlaneWaves_t
struct PlaneWaves_t {
    std::vector<PlaneWave_t> collection; // dimension(:), pointer
    int32_t nc = 0;
    int32_t nC_max = 0;
    bool operator==(const PlaneWaves_t& other) const {
        if (!(collection == other.collection)) return false;
        if (!(nc == other.nc)) return false;
        if (!(nC_max == other.nC_max)) return false;
        return true;
    }

};

// Type: Curr_Field_Src_t
struct Curr_Field_Src_t {
    std::vector<NFDETypes_m::coords_scaled_t> c1P; // dimension(:), pointer
    std::vector<NFDETypes_m::coords_scaled_t> c2P; // dimension(:), pointer
    std::string nombre; // len=BUFSIZE
    int32_t n_C1P = 0;
    int32_t n_C2P = 0;
    bool isElec = false;
    bool isHard = false;
    bool isInitialValue = false;
    bool operator==(const Curr_Field_Src_t& other) const {
        if (!(c1P == other.c1P)) return false;
        if (!(c2P == other.c2P)) return false;
        if (!(nombre == other.nombre)) return false;
        if (!(n_C1P == other.n_C1P)) return false;
        if (!(n_C2P == other.n_C2P)) return false;
        if (!(isElec == other.isElec)) return false;
        if (!(isHard == other.isHard)) return false;
        if (!(isInitialValue == other.isInitialValue)) return false;
        return true;
    }

};

// Type: NodSource_t
struct NodSource_t {
    std::vector<Curr_Field_Src_t> NodalSource; // dimension(:), pointer
    int32_t n_nodSrc = 0;
    int32_t n_nodSrc_max = 0;
    int32_t n_C1P_max = 0;
    int32_t n_C2P_max = 0;
    bool operator==(const NodSource_t& other) const {
        if (!(NodalSource == other.NodalSource)) return false;
        if (!(n_nodSrc == other.n_nodSrc)) return false;
        if (!(n_nodSrc_max == other.n_nodSrc_max)) return false;
        if (!(n_C1P_max == other.n_C1P_max)) return false;
        if (!(n_C2P_max == other.n_C2P_max)) return false;
        return true;
    }

};

// Type: MatrizMedios_t
struct MatrizMedios_t {
    int32_t totalX = 0;
    int32_t totalY = 0;
    int32_t totalZ = 0;
    bool operator==(const MatrizMedios_t& other) const {
        if (!(totalX == other.totalX)) return false;
        if (!(totalY == other.totalY)) return false;
        if (!(totalZ == other.totalZ)) return false;
        return true;
    }

};

// Type: NFDEGeneral_t
struct NFDEGeneral_t {
    RK dt = 0.0;
    int32_t nmax = 0;
    bool mtlnProblem = false;
    bool operator==(const NFDEGeneral_t& other) const {
        if (!(dt == other.dt)) return false;
        if (!(nmax == other.nmax)) return false;
        if (!(mtlnProblem == other.mtlnProblem)) return false;
        return true;
    }

};

// Type: Desplazamiento_t
struct Desplazamiento_t {
    std::vector<RK> desX; // dimension(:), pointer
    std::vector<RK> desY; // dimension(:), pointer
    std::vector<RK> desZ; // dimension(:), pointer
    int32_t nX = 0;
    int32_t mx1 = 0;
    int32_t mx2 = 0;
    int32_t nY = 0;
    int32_t my1 = 0;
    int32_t my2 = 0;
    int32_t nZ = 0;
    int32_t mz1 = 0;
    int32_t mz2 = 0;
    RK originx = 0.0;
    RK originy = 0.0;
    RK originz = 0.0;
    bool operator==(const Desplazamiento_t& other) const {
        if (!(desX == other.desX)) return false;
        if (!(desY == other.desY)) return false;
        if (!(desZ == other.desZ)) return false;
        if (!(nX == other.nX)) return false;
        if (!(mx1 == other.mx1)) return false;
        if (!(mx2 == other.mx2)) return false;
        if (!(nY == other.nY)) return false;
        if (!(my1 == other.my1)) return false;
        if (!(my2 == other.my2)) return false;
        if (!(nZ == other.nZ)) return false;
        if (!(mz1 == other.mz1)) return false;
        if (!(mz2 == other.mz2)) return false;
        if (!(originx == other.originx)) return false;
        if (!(originy == other.originy)) return false;
        if (!(originz == other.originz)) return false;
        return true;
    }

};

// Type: Parseador_t
struct Parseador_t {
    std::string switches = " "; // len=BUFSIZE
    NFDEGeneral_t* general = nullptr;
    MatrizMedios_t* matriz = nullptr;
    Desplazamiento_t* despl = nullptr;
    Frontera_t* front = nullptr;
    Materials_t* Mats = nullptr;
    PECRegions_t* pecRegs = nullptr;
    PECRegions_t* pmcRegs = nullptr;
    DielectricRegions_t* DielRegs = nullptr;
    LossyThinSurfaces_t* LossyThinSurfs = nullptr;
    FreqDepenMaterials_t* frqDepMats = nullptr;
    ANISOTROPICelements_t* aniMats = nullptr;
    Boxes_t* boxSrc = nullptr;
    PlaneWaves_t* plnSrc = nullptr;
    NodSource_t* nodSrc = nullptr;
    Sondas_t* oldSONDA = nullptr;
    MasSondas_t* Sonda = nullptr;
    BloqueProbes_t* BloquePrb = nullptr;
    VolProbes_t* VolPrb = nullptr;
    ThinWires_t* tWires = nullptr;
    SlantedWiresInfo_t* sWires = nullptr;
    ThinSlots_t* tSlots = nullptr;
    ConformalPECRegions_t* conformalRegs = nullptr;
#ifdef CompileWithMTLN
    mtln_t* mtln = nullptr;
#endif
};

// Type: t_linea_t
struct t_linea_t {
    int32_t LEN = 0;
    std::string dato; // len=BUFSIZE
    bool operator==(const t_linea_t& other) const {
        if (!(LEN == other.LEN)) return false;
        if (!(dato == other.dato)) return false;
        return true;
    }

};

// Type: t_NFDE_FILE_t
struct t_NFDE_FILE_t {
    int64_t targ = 0;
    int64_t numero = 0;
    std::vector<t_linea_t> lineas; // dimension(:), pointer
    bool thereare_stoch = false;
    bool operator==(const t_NFDE_FILE_t& other) const {
        if (!(targ == other.targ)) return false;
        if (!(numero == other.numero)) return false;
        if (!(lineas == other.lineas)) return false;
        if (!(thereare_stoch == other.thereare_stoch)) return false;
        return true;
    }

};

} // namespace NFDETypes_m

#endif // NFDE_TYPES_H