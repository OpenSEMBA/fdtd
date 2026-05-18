#include <vector>
#include <string>
#include <complex>
#include <cstdint>

// Assuming FDETYPES_m provides RKIND
// #include "FDETYPES_m.hpp" 

// Assuming conformal_types_m provides triangle_t, interval_t, bufsize
// #include "conformal_types_m.hpp"

// Assuming mtln_types_m is used if CompileWithMTLN is defined
// #ifdef CompileWithMTLN
// #include "mtln_types_m.hpp"
// #endif

namespace NFDETypes_m {

    // Constants from FDETYPES_m or similar, assumed to be available
    // In a real translation, these would be defined in the included headers
    // extern const int RKIND; 
    // For the purpose of this translation, we assume RKIND is defined elsewhere
    // or we define a placeholder. Since the prompt asks to preserve names,
    // we rely on the included headers to provide RKIND.
    
    // If RKIND is not found in headers, we might need to define it locally 
    // or assume it's 8 (double) based on typical usage, but strict adherence 
    // to "preserve names" means we shouldn't invent definitions. 
    // However, C++ requires definitions. We will assume the headers provide it.
    
    // Placeholder for RKIND if not provided by headers, typically 8 for double
    // Note: In a real project, this should come from FDETYPES_m
    // const int RKIND = 8; 

    using RK = double; // Assuming RKIND=8 maps to double. If RKIND=4, use float.
                       // Standard Fortran real(kind=8) is double.

    // CONSTANTS FOR THE PARSER
    // global variable stochastic (commented out in Fortran, so ignored)

    // MATERIALS
    constexpr double SIGMA_PEC = 1e19;
    constexpr double SIGMA_PMC = 1e19;

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
    };

    // Basic constants for materials
    struct Material_t {
        double eps = 0.0;
        double mu = 0.0;
        double sigma = 0.0;
        double sigmam = 0.0;
        int32_t id = 0;
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
    };

    // Identifies conformal PEC "media"
    struct ConformalPECElements_t {
        std::vector<triangle_t> triangles; // Assuming triangle_t is defined in conformal_types_m
        std::vector<interval_t> intervals; // Assuming interval_t is defined in conformal_types_m
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
    };

    // Locates all the different Comp media found
    struct LossyThinSurfaces_t {
        std::vector<LossyThinSurface_t> cs;
        int32_t length = 0;
        int32_t length_max = 0;
        int32_t nC_max = 0;
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
    };

    // List of the different thin wires that were found in the file
    struct ThinWires_t {
        std::vector<ThinWire_t> tw;
        int32_t n_tw = 0;
        int32_t n_tw_max = 0;
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
    };

    // List of the different thin wires that were found in the file (Slanted)
    struct SlantedWiresInfo_t {
        std::vector<SlantedWire_t> sw;
        int32_t n_sw = 0;
        int32_t n_sw_max = 0;
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
    };

    // ThinSlot component that defines the overall properties of the definition of ThinSlots in ORIGINAL
    struct ThinSlot_t {
        std::vector<ThinSlotComp_t> tgc;
        double width = 0;
        int32_t n_tgc = 0;
        int32_t n_tgc_max = 0;
    };

    // List of the different thin Slots that were found in the file
    struct ThinSlots_t {
        std::vector<ThinSlot_t> tg;
        int32_t n_tg = 0;
        int32_t n_tg_max = 0;
    };

    // PML Border Type
    struct FronteraPML_t {
        double orden = 2.0;
        double refl = 1e-3;
        int32_t numCapas = 8;
    };

    // Tipo de la frontera
    struct Frontera_t {
        int32_t tipoFrontera[6];
        FronteraPML_t propiedadesPML[6];
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
    };

    // type that defines a list of probes to be appended and accessed
    struct MasSondas_t {
        std::vector<MasSonda_t> collection;
        int32_t length = 0;
        int32_t length_max = 0;
        int32_t len_cor_max = 0;
    };

    // This type contains the basic information in nearly all the different PROBES
    struct Sonda_t {
        std::string grname;
        std::vector<int32_t> i;
        std::vector<int32_t> j;
        std::vector<int32_t> K;
    };

} // namespace NFDETypes_m

#include <vector>
#include <string>
#include <memory>
#include <cstdint>

// Assuming these types are defined in previous chunks or headers
// #include "coords_t.h"
// #include "coords_scaled_t.h"
// #include "Frontera_t.h"
// #include "Materials_t.h"
// #include "PECRegions_t.h"
// #include "DielectricRegions_t.h"
// #include "LossyThinSurfaces_t.h"
// #include "FreqDepenMaterials_t.h"
// #include "ANISOTROPICelements_t.h"
// #include "ThinWires_t.h"
// #include "SlantedWiresInfo_t.h"
// #include "ThinSlots_t.h"
// #include "ConformalPECRegions_t.h"
// #include "mtln_t.h"
// #include "MasSondas_t.h"

// Constants and Types assumed from context
#ifndef BUFSIZE
#define BUFSIZE 256
#endif

#ifndef RK
#define RK 8
#endif

#ifndef RKIND
#define RKIND 8
#endif

using RK = double;

// Forward declarations for types referenced in pointers
struct NFDEGeneral_t;
struct MatrizMedios_t;
struct Desplazamiento_t;
struct Frontera_t;
struct Materials_t;
struct PECRegions_t;
struct DielectricRegions_t;
struct LossyThinSurfaces_t;
struct FreqDepenMaterials_t;
struct ANISOTROPICelements_t;
struct Boxes_t;
struct PlaneWaves_t;
struct NodSource_t;
struct Sondas_t;
struct MasSondas_t;
struct BloqueProbes_t;
struct VolProbes_t;
struct ThinWires_t;
struct SlantedWiresInfo_t;
struct ThinSlots_t;
struct ConformalPECRegions_t;
struct mtln_t;

// Helper to simulate Fortran NULL pointer for derived types
template <typename T>
struct FortranPointer {
    T* ptr = nullptr;
    FortranPointer() : ptr(nullptr) {}
    FortranPointer(T* p) : ptr(p) {}
    operator T*() const { return ptr; }
    T* operator->() const { return ptr; }
    T& operator*() const { return *ptr; }
};

// Type: Sonda_t
struct Sonda_t {
    std::vector<int32_t> node; // dimension(:), pointer
    int32_t n_cord = 0;
    int32_t n_cord_max = 0;
    RK tstart = 0.0;
    RK tstop = 0.0;
    RK tstep = 0.0;
    std::string outputrequest; // len=BUFSIZE
    RK fstart = 0.0;
    RK fstop = 0.0;
    RK fstep = 0.0;
    RK phistart = 0.0;
    RK phistop = 0.0;
    RK phistep = 0.0;
    RK thetastart = 0.0;
    RK thetastop = 0.0;
    RK thetastep = 0.0;
    std::string FileNormalize; // len=BUFSIZE
};

// Type: FarField_Sonda_t
struct FarField_Sonda_t {
    Sonda_t probe;
};

// Type: Electric_Sonda_t
struct Electric_Sonda_t {
    Sonda_t probe;
};

// Type: Magnetic_Sonda_t
struct Magnetic_Sonda_t {
    Sonda_t probe;
};

// Type: NormalElectric_Sonda_t
struct NormalElectric_Sonda_t {
    Sonda_t probe;
    std::vector<int32_t> nml; // dimension(:), pointer
    int32_t n_nml = 0;
    int32_t n_nml_max = 0;
};

// Type: NormalMagnetic_Sonda_t
struct NormalMagnetic_Sonda_t {
    Sonda_t probe;
    std::vector<int32_t> nml; // dimension(:), pointer
    int32_t n_nml = 0;
    int32_t n_nml_max = 0;
};

// Type: SurfaceElectricCurrent_Sonda_t
struct SurfaceElectricCurrent_Sonda_t {
    Sonda_t probe;
    std::vector<int32_t> nml; // dimension(:), pointer
    int32_t n_nml = 0;
    int32_t n_nml_max = 0;
};

// Type: SurfaceMagneticCurrent_Sonda_t
struct SurfaceMagneticCurrent_Sonda_t {
    Sonda_t probe;
    std::vector<int32_t> nml; // dimension(:), pointer
    int32_t n_nml = 0;
    int32_t n_nml_max = 0;
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
};

// Type: Sondas_t
struct Sondas_t {
    std::vector<abstractSonda_t> probes; // dimension(:), pointer
    int32_t n_probes = 0;
    int32_t n_probes_max = 0;
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
};

// Type: BloqueProbes_t
struct BloqueProbes_t {
    std::vector<BloqueProbe_t> bp; // dimension(:), pointer
    int32_t n_bp = 0;
    int32_t n_bp_max = 0;
};

// Type: VolProbe_t
struct VolProbe_t {
    std::vector<coords_t> cordinates; // dimension(:), pointer
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
};

// Type: VolProbes_t
struct VolProbes_t {
    std::vector<VolProbe_t> collection; // dimension(:), pointer
    int32_t length = 0;
    int32_t length_max = 0;
    int32_t len_cor_max = 0;
};

// Type: Box_t
struct Box_t {
    std::string nombre_fichero; // len=BUFSIZE
    int32_t coor1[3];
    int32_t coor2[3];
};

// Type: Boxes_t
struct Boxes_t {
    std::vector<Box_t> Vols; // dimension(:), pointer
    int32_t nVols = 0;
    int32_t nVols_max = 0;
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
};

// Type: PlaneWaves_t
struct PlaneWaves_t {
    std::vector<PlaneWave_t> collection; // dimension(:), pointer
    int32_t nc = 0;
    int32_t nC_max = 0;
};

// Type: Curr_Field_Src_t
struct Curr_Field_Src_t {
    std::vector<coords_scaled_t> c1P; // dimension(:), pointer
    std::vector<coords_scaled_t> c2P; // dimension(:), pointer
    std::string nombre; // len=BUFSIZE
    int32_t n_C1P = 0;
    int32_t n_C2P = 0;
    bool isElec = false;
    bool isHard = false;
    bool isInitialValue = false;
};

// Type: NodSource_t
struct NodSource_t {
    std::vector<Curr_Field_Src_t> NodalSource; // dimension(:), pointer
    int32_t n_nodSrc = 0;
    int32_t n_nodSrc_max = 0;
    int32_t n_C1P_max = 0;
    int32_t n_C2P_max = 0;
};

// Type: MatrizMedios_t
struct MatrizMedios_t {
    int32_t totalX = 0;
    int32_t totalY = 0;
    int32_t totalZ = 0;
};

// Type: NFDEGeneral_t
struct NFDEGeneral_t {
    RK dt = 0.0;
    int32_t nmax = 0;
    bool mtlnProblem = false;
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
};

// Type: t_NFDE_FILE_t
struct t_NFDE_FILE_t {
    int64_t targ = 0;
    int64_t numero = 0;
    std::vector<t_linea_t> lineas; // dimension(:), pointer
    bool thereare_stoch = false;
};